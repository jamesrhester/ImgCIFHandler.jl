# Intelligently extract something from a zip file

# Assumes that the server supports ranges

struct ZipArchive <: CompressedArchive
    local_cache::AbstractString
    original_url::URI
    file_directory::Dict{String, Tuple{Int64, Int64}}
    insecure::Bool  #flag that original_url should be adjusted
end

"""
    ZipArchive(local_cache, original_url; insecure = false

Create a ZipArchive instance. If `insecure` is true, allow
downloading by ftp/http even if ftps/https was requested.
"""
ZipArchive(local_cache, original_url; insecure = false) = begin

    if make_insecure(original_url) == original_url
        insecure = true
    end
    
    hs = handles_secure(original_url)

    if !hs && insecure
        real_url = make_insecure(original_url)
        use_insecure = true
    elseif !hs
        throw(error("Cannot securely contact server for $original_url. Consider
using `insecure` keyword to allow non-SSL connections"))
    else
        real_url = original_url
        use_insecure = false
    end

    @debug "Real url $real_url"
    
    if supports_range_header(real_url)

        # Get total size of file

        full_size = get_file_size(request(real_url))

        # Some simple servers can't handle empty start bytes
        
#       if full_size <= 0
            range_str = "bytes=-65000" #hope this one understands
#        else
#            range_str = "$(full_size-65000)-$(full_size - 1)"
#        end

        @debug "Using $range_str for download"
        
        # We generate the index

        io = IOBuffer()
        req = request(real_url, headers = ["Range" => range_str], output = io)

        @debug "Downloaded 65000 last bytes from $real_url" req

        cdir, cdir_len, num_entries = process_eocd(io, maxoffset = 100)

        @debug "Have EOCD" cdir, cdir_len, num_entries
    
        # Now download the full end of the zip file

        seek(io, 0)
        req = request(real_url, headers=["Range" => "bytes=$cdir-"], output = io)

        @debug "Result of central directory request" req
        
        file_dir = interpret_cdfh(io, 0, num_entries)
        ZipArchive(local_cache, URI(original_url), Dict(file_dir), use_insecure)
    else

        @debug "No range headers supported"
        
        ZipArchive(local_cache, URI(original_url), Dict(), insecure)
    end
end

get_real_url(a::ZipArchive) = a.insecure ? make_insecure(a.original_url) : a.original_url

download_images_os(a::ZipArchive, ext_info) = begin

    loc = get_local_dir(a)
    need_to_get = filter( x-> !has_local_version(a, x), ext_info; view = true)
    @debug "For $ext_info need to get $need_to_get"
    if size(need_to_get, 1) > 0
        arch_paths = need_to_get.archive_path
        for ntg in eachrow(need_to_get)
            ap = ntg.archive_path
            start, comp_len = get_file_offset(a, ap)
            if start >= 0
                io = IOBuffer()
                download("$(get_real_url(a))", io, headers = ["Range" => "bytes=$start-$(start + comp_len)"])
                zip_deflate(io, 0, comp_len, local_equivalent(a, ntg))
            else
                @warn "Unable to download $ntg"
            end
        end
    end
    
    @debug "Expect files at" local_equivalent.(Ref(a), eachrow(need_to_get))
    @assert all( x-> has_local_version(a, x), eachrow(need_to_get))

end

get_file_offset(a::ZipArchive, fname) = begin
    get(a.file_directory, fname, (-1,0))
end

peek_image(a::ZipArchive, cif_block::CifContainer; entry_no=0, max_down = 2.0e7) = begin

    # Find something already present if possible
    
    p = get_any_local(a, cif_block)
    if !isnothing(p) return p end

    ei = first(get_loop(cif_block, "_array_data_external_data.id"))
    uri = getproperty(ei,"_array_data_external_data.uri")
    path_in_arch = getproperty(ei, "_array_data_external_data.archive_path")

    if length(a.file_directory) == 0

        # Get entire Zip file

        if !ispath(local_equivalent(a))

            # Determine file size
        
            fs = get_file_size(request(uri))
        
            if fs > max_down
                throw(error("File size $fs too large (max $max_down). Consider using -l or -s options after unpacking locally, or increasing maximum download using -m"))
            elseif fs == 0
                throw(error("Unable to determine file size for $uri. Download and unpack locally then use -l or -s options"))
            end
    
            @debug "Now downloading $uri"
        
            Downloads.download(uri, local_equivalent(a))

        end

        # Extract an image

        @debug "Extracting $(local_equivalent(a, ei)) from $uri at $(local_equivalent(a))"

        full_archive = open(local_equivalent(a),"r")
        zr = ZipReader(read(full_archive))
        if !(path_in_arch in zip_names(zr))
            throw(error("Bad archive $uri: missing $path_in_arch"))
        end

        in_mem = zip_readentry(zr, path_in_arch)
        fo = open(local_equivalent(a, ei), "w")
        write(fo, in_mem)
        close(fo)
    else

        # Download the specific file

        start, file_len = get_file_offset(a, path_in_arch)

        if start < 0
            throw(error("File $path_in_arch not found in $uri"))
        end

        io = IOBuffer()
        resp = request("$(get_real_url(a))", headers = ["Range" => "bytes=$start-$(start+file_len)"],
                       output = io)

        @debug "Obtained compressed chunk" resp
        
        zip_deflate(io, 0, file_len, local_equivalent(a, ei))
    end
        
    return local_equivalent(a, ei)
    
end

"""
    Extract the nth file from a Zip archive    
"""
get_nth_zip_file(zip_archive, n::Int, target_dir) = begin

    all_files = get_zip_dir(zip_archive)

    @debug "File information" all_files
    
    fname, (file_offset, comp_len) = all_files[n]
    zf = open(zip_archive, "r")
    zip_deflate(zf, file_offset, comp_len, joinpath(target_dir, fname))
    
end

access_remote_zip(zip_url, start, finish) = begin

    if !supports_range_header(zip_url)
        throw(error("Server doesn't support range: $zip_url"))
    end

    hdrs = request(zip_url)   #just the headers
    file_size = get_file_size(hdrs)
    if finish > file_size
        throw(error("End pos $finish too large, file is only $file_size"))
    end

    io = IOBuffer(sizehint = finish - start)
    result = download(zip_url, output = io, headers=["Range"=>"bytes=$start-$finish"])
    return io
end

"""
    Return a list of files in the zip file
"""
get_zip_dir(zipfile::String) = begin
    
    bytes = open(zipfile,"r")

    get_zip_dir(bytes)
end

get_zip_dir(io::IO) = begin

    cdir, cdir_len, num_entries = process_eocd(io, maxoffset = 100)

    @debug "Have EOCD" cdir, cdir_len, num_entries
    
    interpret_cdfh(io, cdir, num_entries)
end

"""
   Determine values and offsets from end record
"""
process_eocd(io::IO; maxoffset = 64000) = begin

    seekend(io)    #start at the end

    end_pos = position(io)
    
    if end_pos < 22
        throw(error("Not enough bytes, need at least 20 to find Zip signature"))
    end
    
    # find start of eocd relative to end

    step_back = 20
    seek(io, end_pos - step_back)  #minimum value is 20

    while true

        sig = readle(io, UInt32)

        if sig == 0x06054b50
            break
        end
        if position(io) < 2 || end_pos - position(io) > maxoffset
            throw(error("Zip EOCD not found, supply more bytes"))
        end

        step_back += 2
        seek(io, end_pos - step_back)
    end

    # Little endian according to spec
    skip(io, 6)
    num_entries = readle(io, UInt16)
    cdir_len = readle(io, UInt32)
    cdir_offset = readle(io, UInt32)

    return Int64(cdir_offset), Int64(cdir_len), Int64(num_entries) 
end

"""
    Extract information from the central directory. The file length
returned includes header information.
"""
interpret_cdfh(io::IO, offset, num_entries) = begin

    file_list = []
    while length(file_list) < num_entries

        @debug "Processing file" length(file_list) offset

        fname, file_loc, comp_len, offset = get_zip_file_entry(io, offset)

        @debug "New offset" fname offset
        
        push!(file_list, (fname, (file_loc, comp_len)))
    end

    return file_list
end

"""
    Return information for one file
"""
get_zip_file_entry(io::IO, offset) = begin

    seek(io, offset)
    sig = readle(io, UInt32)
    if sig != 0x02014b50

        @debug "Can't find sig" sig
        
        throw(error("Incorrect offset for file entry: $offset $sig"))
    end

    
    seek(io, offset + 20)
    comp_len = readle(io, UInt32)
    seek(io, offset + 28)
    nlen = readle(io, UInt16)
    mlen = readle(io, UInt16)
    klen = readle(io, UInt16)
    seek(io, offset + 42)
    file_loc = readle(io, UInt32)

    fname = transcode(String, read(io, nlen))
    next_entry = offset + 46 + nlen + mlen + klen
    chunk_len = 30 + nlen + mlen + comp_len
    return fname, file_loc, chunk_len, next_entry

end

"""
    Deflate a stream into the supplied file
"""
zip_deflate(io::IO, offset, chunk_len, dest_file) = begin

    @debug "Deflating" offset dest_file
    
    seek(io, offset)

    sig = readle(io, UInt32)
    if sig != 0x04034b50

        @debug "Can't find sig" sig
        
        throw(error("Incorrect offset for file entry: $offset $sig"))
    end

    # Get general purpose bit flag
    seek(io, offset + 6)
    bitf = readle(io, UInt16)

    @debug "Bit flags" bitf
    
    if bitf & 0x1 << 3 > 0
        @warn "Zip: no compressed size available in local file header"
    end

    seek(io, offset + 8)
    decomp_type = readle(io, UInt16)
    seek(io, offset + 18)
    comp_size = readle(io, UInt32)
    decomp_size = readle(io, UInt32)
    seek(io, offset + 26)
    nlen = readle(io, UInt16)
    mlen = readle(io, UInt16)

    fname = transcode(String, read(io, nlen))

    if comp_size == 0    #didn't know at time of writing
        comp_size = chunk_len -30 - nlen - mlen
    end
    
    @debug "Decomp $fname starting at $(offset + 30 + nlen + mlen)" decomp_type Int64(decomp_size) Int64(comp_size) Int64(nlen) Int64(mlen)
    
    seek(io, offset + 30 + nlen + mlen)

    df = open(dest_file, "w")
    if decomp_type == 0
        write(df, read(io, comp_size))
        return
    end

    clipped = IOBuffer(read(io, comp_size))

    if decomp_type == 8    #ordinary default
        g = CodecZlib.DeflateDecompressorStream(clipped)
    elseif decomp_type == 9
        g = CodecInflate64.Deflate64DecompressorStream(clipped)
    else
        throw(error("Unrecognised zip file compression $decomp_type"))
    end

    write(df, read(g))
    close(df)
end

# Copied from ZipArchive.jl

readle(io::IO, ::Type{UInt64}) = UInt64(readle(io, UInt32)) | UInt64(readle(io, UInt32))<<32
readle(io::IO, ::Type{UInt32}) = UInt32(readle(io, UInt16)) | UInt32(readle(io, UInt16))<<16
readle(io::IO, ::Type{UInt16}) = UInt16(read(io, UInt8)) | UInt16(read(io, UInt8))<<8
readle(io::IO, ::Type{UInt8}) = read(io, UInt8)
