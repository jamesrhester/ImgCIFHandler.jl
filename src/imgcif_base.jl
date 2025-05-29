
get_image_ids(c::CifContainer) = begin
    return c["_array_data.binary_id"]
end

# Archive types
"""
In order to access images, an archive must be supplied
to the `imgload` routine.

An archive is either remote, or local. A local archive may
have a remote equivalent, such that it is a cached version.

All archives must support the `local_equivalent` method,
to supply the cached local location of the referenced
frame, and the `download_images_os` method, to return
a specified frame.
"""
abstract type ImageArchive end

"""
A `BulkArchive` consists of a single file containing
multiple image frames.
"""
abstract type BulkArchive <: ImageArchive end

abstract type CompressedArchive <: BulkArchive end
abstract type CompressedTarArchive <: CompressedArchive end

"""
In an `AddressableArchive` each frame can be referenced
by an individual, unique URL
""" 
abstract type AddressableArchive <: ImageArchive end

"""
A tar archive has been created by tar and then optionally
compressed. The type parameter T refers to the compression used in
the tar archive. To add a new compression type add it
to the `decomp_option` method.
"""
struct TarArchive{T} <: CompressedTarArchive
    local_cache::AbstractString
    original_url::URI
end

# Zip archives are defined separately

"""
When frames are not bundled into a separate archive format, we have
a directory archive.
"""
struct DirectoryArchive <: AddressableArchive
    local_cache::AbstractString
    original_url::URI
end

DirectoryArchive(local_dir, all_urls::Vector{URI}) = begin

    DirectoryArchive(local_dir, get_constant_part(all_urls))
        
end

create_archives(c::Cif; kwargs...) = create_archive(first(c).second; kwargs...)

"""
    create_archives(c::CifContainer; subs = Dict())

Create a series of ImageArchive objects, one for each distinct location
in `c`. Note that a distinct location is not necessarily a distinct URL,
as rsync: and file: URLs may be unique for each frame but correspond
to a single directory. `subs` is an optional dictionary keyed by
URL, where the value is an existing local directory containing the
contents of that URL, uncompressed and with the same directory
structure as the archive, so that archive_path can be used to access
the image frame.
"""
create_archives(c::CifContainer; subs = Dict()) = begin

    cat_name = "_array_data_external_data"
    arch_type = haskey(c, "$cat_name.archive_format") ? unique(c["$cat_name.archive_format"]) : nothing
    
    dl_info = unique(c["$cat_name.uri"])
    test_uris = URI.(dl_info)

    @debug "Unique URIs before resolution" test_uris

    # catch relative URIs and make them absolute

    my_dir = URI(scheme="file", path = "$(c.original_file)")
    test_uris = map(x -> resolvereference(my_dir, x), test_uris)

    @debug "URIs are" test_uris

    schemes = unique([ x.scheme for x in test_uris])

    # Check for sanity

    if length(schemes) > 1
        throw(error("More than one URI scheme present in file: $schemes"))
    end
    schemes = schemes[]

    if !isnothing(arch_type) && length(arch_type) > 1
        throw(error("More than one archive type present in file: $arch_type"))
    elseif !isnothing(arch_type)
        arch_type = arch_type[]
    end
        
    create_archives(test_uris; arch_type = arch_type, subs = subs, root_dir = "$(dirname(c.original_file))")
    
end

create_archives(u::Vector{URI}; arch_type = nothing, subs = Dict(), root_dir="") = begin

    arch_list = ImageArchive[]

    uq_u = unique(u)

    base_dirs = map(uq_u) do one_uri
        if one_uri in keys(subs)
            @debug "Found local sub" subs[one_uri]
            subs[one_uri]
        elseif "Universal" in keys(subs)
            subs["Universal"]
        else
            if !isempty(subs)
                @warn "No local substitute for $one_uri found"
            end
            mktempdir()
        end
    end

    if arch_type in ("TBZ","TGZ","TXZ","TAR")

        for (one_uri, bd) in zip(uq_u, base_dirs)
            push!(arch_list, TarArchive{Val(Symbol(arch_type))}(bd, one_uri))
        end
        return arch_list

    end

    if arch_type == nothing

        if first(u).scheme == "file" || first(u).scheme == ""

            # No need to cache at all
            
            return [LocalArchive(root_dir)]

        else

            const_part = "$(get_constant_part(u))"
            if const_part in keys(subs)
                @debug "Found local sub:" subs[const_part]
                base_dir = subs[const_part]
            else
                @debug "$const_part not in subs list" keys(subs)
                base_dir = base_dirs[1]
            end

            return [DirectoryArchive(base_dir, u)]
    
        end
    end
    
    if arch_type == "ZIP"
        for (one_uri, bd) in zip(uq_u, base_dirs)
            push!(arch_list, ZipArchive(bd, one_uri))
        end
        return arch_list
    end
        
    throw(error("Unrecognised archive type $arch_type"))
end

create_archives(u::URI; kwargs...) = create_archives([u]; kwargs...)

decomp_option(::TarArchive{T}) where {T} = begin
    @debug "Type parameter is" T
    if T == Val(:TGZ) return "-z" end
    if T == Val(:TBZ) return "-j" end
    if T == Val(:TXZ) return "-J" end
    return nothing
end

get_local_dir(a::BulkArchive) = a.local_cache

struct LocalArchive <: AddressableArchive
    root_dir::String    #For relative URIs
end

# Getting an image from an archive

download_images_os(a::ImageArchive, ext_info) = begin
    @warn "No image download defined for $(typeof(a))"
end

download_images_os(a::TarArchive, ext_info) = begin

    cmd_list = Cmd[]
    loc = get_local_dir(a)
    need_to_get = filter( x-> !has_local_version(a, x), ext_info; view=true)

    @debug "For $ext_info need to get $need_to_get"
    if size(need_to_get, 1) > 0
        arch_paths = need_to_get.archive_path
        push!(cmd_list, Cmd(`curl -s --show-error $(ext_info.uri[1])`,ignorestatus=true))
        j = decomp_option(a)
        if !isnothing(j)
            push!(cmd_list, `tar -C $loc -x $j -f - --occurrence $arch_paths`)
        else
            push!(cmd_list, `tar -C $loc -x -f - --occurrence $arch_paths`)
        end
                  
        @debug "Commands to extract" cmd_list
        
        try
            run(pipeline(cmd_list...))
        catch exc
            @debug "Finished downloading" exc
        end

        # decompress if necessary

        decompress_frames(a, need_to_get)

        # verify

        @debug "Expect files at" local_equivalent.(Ref(a), eachrow(need_to_get))
        @assert all( x-> has_local_version(a, x), eachrow(need_to_get))
    end

end

"""
    download_images_os(a::DirectoryArchive, ext_info)

This will not check for local existence of the file, so check
using `has_local_version` before calling this.
"""
download_images_os(a::DirectoryArchive, ext_info; max_down=2e7) = begin

    c = get_prefix(ext_info)

    ext_scheme = a.original_url.scheme
    
    aurl = "$(a.original_url)"
    for r in eachrow(ext_info)

        full_uri = "$(getproperty(r,"$(c)uri"))"
        if !startswith("$full_uri", aurl)
            throw(error("Asked to download from $full_uri but archive is for $aurl"))
        end
        
        if ext_scheme == "rsync"

            # SBGRID does not download the full directory hierarchy
            # from the URL, therefore we insert a "." so that we can
            # control it

            download_uri = joinpath(aurl,".",full_uri[(length(aurl)+2):end])

            @debug "Rsync download address" download_uri
            run(`rsync -avR $download_uri $(a.local_directory)`)

        elseif ext_scheme in ("https", "http")

            r = request(full_uri)
            fs = get_file_size(r)
        
            if fs == 0 && max_down > 0

                throw(error("Cannot determine file size for $uri. Download manually and use -l or -s options"))
                
            elseif fs > max_down

                throw(error("$uri size $fs exceeds maximum download limit of $max_down, please use -m option or download manually and use -l / -s options"))
            end
        
            @debug "Downloading $(local_equivalent(a, ei)) from $uri"
            
            Downloads.download(full_uri, output = local_equivalent(a, r))
        else
            throw(error("Unrecognised URI scheme $ext_scheme"))
        end

        # Sanity check
        @assert ispath(local_equivalent(a, r))
            
    end

end

download_images_os(a::LocalArchive, ext_info) = begin
    for r in eachrow(ext_info)
        if !ispath(local_equivalent(a, r))
            throw(error("File $local_equivalent does not exist"))
        end
    end
end

"""
    get_prefix(ext_info)

Return the category name prefix for ext_info
"""
get_prefix(ext_info) = begin
    n = first(names(ext_info))
    test_prefix = "_array_data_external_data."
    if occursin(test_prefix, n)
        return test_prefix
    end

    return ""
end

"""
    check_match_uri

Make sure the URL of the compressed archive matches the URL
provided in `ext_info`
"""
check_match_uri(a::CompressedArchive, ext_info) = begin
    
    c = get_prefix(ext_info)
    test_uri = URI(getproperty(ext_info, "$(c)uri"))

    test_base = URI(scheme="file", path=pwd())
    if resolvereference(test_base, test_uri) == test_uri
        return test_uri == a.original_url
    else
        # is extremely unportable relative URL, don't try
        return true
    end
end

"""
    get_constant_part(v::Vector{URI})

Return a URI formed from the constant part of the URIs listed in v
"""
get_constant_part(v::Vector{URI}) = begin
    
    # find unchanging portion of path

    new_path = v[1].path
    
    path_elements = map(x -> splitpath(x.path), v)
    for pe in 1:length(path_elements[1])
        all_entries = [x[pe] for x in path_elements]
        if length(unique(all_entries)) == 1
            continue
        else
            new_path = joinpath(splitpath(v[1].path)[1:pe-1])
        end
    end

    b = v[1]
    URI(scheme=b.scheme, host = b.host,
        path = new_path, fragment = b.fragment,
        port = b.port, query = b.query,
        userinfo = b.userinfo)

end

local_equivalent(a::CompressedArchive, ext_info) = begin
    
    c = get_prefix(ext_info)

    # Check that URIs match

    if !check_match_uri(a, ext_info)
        @debug "Non-matching URI" getproperty(ext_info, "$(c)uri") a.original_url
        return nothing
    end
    
    b = joinpath(a.local_cache, getproperty(ext_info, "$(c)archive_path"))
    if "file_compression" in names(ext_info)
        return b * "_final"
    end
    @debug "Local equivalent:" ext_info b
    return b
end

"""
    local_equivalent(a::LocalArchive, ext_info)

The local file name corresponding to the information provided in `ext_info`
"""
local_equivalent(a::LocalArchive, ext_info) = begin
    c = get_prefix(ext_info)
    base_part = URI(getproperty(ext_info,"$(c)uri")).path
    if !isabspath(base_part)
        base_part = joinpath(a.root_dir, base_part)
    end
    if "$(c)file_compression" in names(ext_info)
        return base_part * "_final"
    else
        return base_part
    end
end

"""
    local_equivalent(a::DirectoryArchive, ext_info)

The local file name corresponding to the information provided in `ext_info`.
"""

local_equivalent(a::DirectoryArchive, ext_info) = begin

    c = get_prefix(ext_info)
    u = URI(getproperty(ext_info,"$(c)uri"))
    base_part = "$(a.original_url)"
    unique_part = "$u"[length(base_part)+2:end]
    base = joinpath(a.local_cache, unique_part)
    if "$(c)file_compression" in names(ext_info)
        return base * "_final"
    else
        return base
    end

end

"""
    local_equivalent(a::CompressedArchive)

The local file name for the archive itself, if downloaded
in full.
"""
local_equivalent(a::CompressedArchive) = begin
    lname = basename(a.original_url.path)
    return joinpath(get_local_dir(a), lname)
end

has_local_version(a::ImageArchive, ext_info) = begin
    l = local_equivalent(a, ext_info)
    if isnothing(l)
        return true
    end
    @debug "Local path is" l
    ispath(l)
end

#===== Image loading =====#

"""
    imgload(c::Block,array_ids, a::ImageArchive)

Return the image referenced in CIF Block `c` corresponding to the specified raw array identifiers,
using `ImageArchive` a
"""
imgload(c::CifContainer,bin_ids, a::ImageArchive) = begin

    dl_info = external_specs_from_bin_ids(bin_ids,c)
    imgload(a, dl_info)
end

imgload(c::CifContainer, bin_ids, a::Vector{<:ImageArchive}) = begin

    dl_info = external_specs_from_bin_ids(bin_ids, c)
    for la in a
        i = imgload(la, dl_info)
        if !isnothing(i)
            return i
        end
    end
    
end

"""
    imgload(a::ImageArchive, ext_info::DataFrame)

Return an image from `a` obtained by summing a series of images corresponding to the 
entries in `ext_info`.
"""
imgload(a::ImageArchive, ext_info::DataFrame) = begin

    # Fix any leftover long column names

    rename!(x-> replace(x,"_array_data_external_data." => ""),ext_info)
    
    temp_locals = local_equivalent.(Ref(a), eachrow(ext_info))
    if nothing in temp_locals
        @debug "Archive not relevant" a ext_info
        return nothing
    end

    if !all(x -> has_local_version(a, x), eachrow(ext_info))
        download_images_os(a, ext_info)
    end
    
    # Now accumulate the image values

    path = "path" in names(ext_info) ? ext_info.path[1] : nothing
    frame = "frame" in names(ext_info) ? parse(Int64,ext_info.frame[1]) : nothing

    # Do a precheck

    can_load = check_format(temp_locals[1],Val(Symbol(ext_info.format[1]));path=path,frame=frame)
    if !can_load[1]
        return can_load[2]
    end
    final_image = imgload(temp_locals[1],Val(Symbol(ext_info.format[1]));path=path,frame=frame)
    for fr_no in 2:size(ext_info,1)
        path = "path" in names(ext_info) ? ext_info.path[fr_no] : nothing
        frame = "frame" in names(ext_info) ? parse(Int64,ext_info.frame[fr_no]) : nothing
        @debug "Accumulating frame $fr_no" path frame
        
        final_image += imgload(temp_locals[fr_no],Val(Symbol(ext_info.format[1]));path=path,frame=frame)
    end
    return final_image
end

imgload(c::CifContainer, scan_id, frame_no::Int, a::Vector{ImageArchive}) = begin
    bin_ids = bin_id_from_scan_frame(c, scan_id, frame_no)
    imgload(c, bin_ids, a)
end

imgload(c::CifContainer, scan_id, frame_no::Int, a::ImageArchive) = imgload(c, scan_id, frame_no, [a])

"""
    imgload(c::CIF)

Return the image referenced by the first encountered `_array_data.binary_id` in the
first block of CIF file `c`. To avoid re-downloading when multiple images are
referenced, first create an archive using `create_archive(c)` and pass this as
the second argument.
"""
imgload(c::Cif) = begin
    cc = first(c).second
    a = create_archives(cc)
    imgload(cc, a)
end

"""
    imgload(c::CifContainer, a::ImageArchive)

Return the image referenced by the first encountered `_array_data.binary_id` in
CIF block `c`. `a` is an archive created using `create_archive`.
"""
imgload(c::CifContainer, a::Vector{<:ImageArchive}) = begin
    f_id = c["_array_data.binary_id"][1]
    imgload(c, f_id, a)
end

imgload(c::CifContainer, a::ImageArchive) = imgload(c, [a])

imgload(s::AbstractString) = begin
    imgload(Cif(s))
end

"""
decompress_frames
"""
decompress_frames(ext_info::DataFrame, loc, file_list) = begin
    fc = ext_info.file_compression[1]
    map(file_list) do tl
        final_file = open(joinpath(loc,tl*"final"),"w")
        if fc == "GZ"
            run(pipeline(`zcat $tl`,final_file))
        elseif fc == "BZ2"
            run(pipeline(`bzcat $tl`,final_file))
        end
        close(final_file)
        tl*"_final"
    end
end

decompress_frames(a::ImageArchive, ext_info) = begin

    if "file_compression" in names(ext_info)
        fc = ext_info.file_compression
        final_name = local_equivalent(a, ext_info)
        compressed_name = local_equivalent(a, ext_info)[1:end-6]
        final_file = open(final_name,"w")
        if fc == "GZ"
            run(pipeline(`zcat $compressed_name`,final_file))
        elseif fc == "BZ2"
            run(pipeline(`bzcat $compressed_name`,final_file))
        end
        close(final_file)
    end
end

make_absolute_uri(c::CifContainer,u::AbstractString) = begin
    resolvereference(URI(c.original_file),URI(u))
end

# Format checking

#==
We define a generic method to pre-check image files before attempting
to load them. This should catch syntax errors, unsupported features,
and so forth. Particular formats can define this method. It returns
a tuple `(bool, [message])` where `bool` is true for success, and
`messages` contains a description of the error on failure.
==#

"""
    check_format(loc,fmt;kwargs...)

Check that the image in local file `loc` with format `fmt` is
valid. Returns `(false,"message")` on failure and `(true,"")` on
success.  
"""
check_format(loc,fmt;kwargs...) = (true,"")

"""
    peek_image(a::LocalArchive, cif_block::CifContainer)

Find the name of an image in `a`.
"""
peek_image(a::LocalArchive, cif_block::CifContainer; kwargs...) = begin

   get_any_local(a, cif_block)
    
end

"""
    peek_image(a::DirectoryArchive, cif_block::CifContainer; entry_no=0, check_name=true, extract=true)

Find the name of an image in `a`.
"""
peek_image(a::DirectoryArchive, cif_block::CifContainer; kwargs...) = begin

    c = "_array_data_external_data"
    if haskey(cif_block,"$c.id")
        r = first(get_loop(cif_block, "$c.id"))
        if !has_local_version(a, r)
            download_images_os(a, r)
        end
        return local_equivalent(a, r)
    end
    
end

peek_image(a::TarArchive, cif_block::CifContainer; entry_no=0, check_name=true, extract=true, kwargs...) = begin

    # Find something already present if possible
    
    p = get_any_local(a, cif_block)
    if !isnothing(p) return p end
    
    cmd_list = Cmd[]
    uri = a.original_url
    
    push!(cmd_list, Cmd(`curl -s $uri`,ignorestatus=true))
    j = decomp_option(a)
    if !isnothing(j)
        push!(cmd_list, Cmd(`tar -t -v $j -f -`,ignorestatus=true))
    else
        push!(cmd_list, Cmd(`tar -t -v -f -`, ignorestatus=true))
    end
    
    awkstr1 =  "\$3 > 0 { print \$NF }"
    awkstr2 =  "\$3 > 0 && FNR >= $entry_no { exit }"
    push!(cmd_list, `awk -e $awkstr1 -e $awkstr2`)

    @debug "Peeking into $uri for $(typeof(a)) starting at $entry_no"
    @debug "Command list is $cmd_list"
    fname = nothing
    try
        fname = readchomp(pipeline(cmd_list...))
    catch exc
        @debug "Finished downloading" exc
    end

    # Just the last file is the one we want
    fname = split(fname,"\n")[end]
    
    @debug fname

    if fname != nothing && extract
        pop!(cmd_list)    #no awk statement
        if !isnothing(j)
            cmd_list[2] = Cmd(`tar -C $(a.local_cache) -x -v $j -f - --occurrence $fname`)
        else
            cmd_list[2] = Cmd(`tar -C $(a.local_cache) -x -v -f - --occurrence $fname`)
        end

        @debug "Downloading" cmd_list
        
        try
            run(pipeline(cmd_list...))
        catch exc
            @debug "Finished downloading" exc
        end
    end

    if fname != nothing && check_name
        c = "_array_data_external_data"
        my_dir = URI(scheme="file", path = "$(cif_block.original_file)")
        if haskey(cif_block,"$c.archive_path")
            pos = indexin([fname],cif_block["$c.archive_path"])[]
            if pos != nothing
                ref_uri = resolvereference(my_dir, cif_block["$c.uri"][pos])
                if ref_uri == uri
                    return joinpath(a.local_cache,fname)
                else
                    @debug "Couldn't match $fname" uri ref_uri
                end
            end
        end
        return nothing
    end

    return joinpath(a.local_cache,fname)

end

get_file_size(r) = 0

get_file_size(r::Response) = begin
    h = Dict(r.headers)
    if "content-length" in keys(h)
        return parse(Int64, h["content-length"])
    end
    return 0
end

"""
    ping_archive(u::ImageArchive)

Check that an archive is accessible
"""
ping_archive(a::BulkArchive, cif_block::CifContainer) = begin

    c = "_array_data_external_data"
    ei = get_loop(cif_block, "$c.id")
    all_uris = unique(getproperty(ei, "$c.uri"))
    bad = filter(all_uris) do au
        !success(pipeline(`curl -s $au -r 0-0 --fail`, devnull))
    end
    if length(bad) > 0
        @debug "Following URLs not available" bad
        return nothing
    end

    return all_uris
end

ping_archive(a::AddressableArchive, cif_block::CifContainer) = begin

    peek_image(a, cif_block)
    
end

"""
    get_any_local(a::BulkArchive, cif_block::CifContainer)

Return a filename that is locally cached
"""
get_any_local(a::ImageArchive, cif_block::CifContainer) = begin

    catname = "_array_data_external_data"
    if haskey(cif_block,"$catname.id")
        for r in eachrow(get_loop(cif_block, "$catname.id"))
            if has_local_version(a, r)
                return local_equivalent(a, r)
            end
        end
    end

end

#=  Determine if a server supports the "range" header. Issues here

- The server may not provide an accept-ranges header unless you actually
do a GET (nginx on Zenodo)
- HTTPS negotiation might fail.
=#

"""
    Handles secure connections if requested
"""
handles_secure(server_url::AbstractString) = begin
    req = request(server_url, headers = ["Range" => "bytes=0-10"],
                  throw = false)
    return !is_ssl_error(req)
end

handles_secure(server_url::URIs.URI) = handles_secure("$server_url")

supports_range_header(server_url::AbstractString) = begin

    downloader = Downloader()

    req = request(server_url, headers = ["Range"=>"bytes=0-10"],
                  throw = false)

    if has_accept_ranges(req)
        return true
    end
    
    # Some servers appear not to return the partial header unless
    # you actually do a POST. HEAD is not enough. So our strategy
    # is to do a 10-byte download and deliberately abort, then
    # study the headers.

    @debug "Trying a proper download to check range header" server_url
    
    easy_hook = (easy, info) -> begin
        Curl.setopt(easy, Curl.CURLOPT_MAXFILESIZE, 10)
    end

    downloader.easy_hook = easy_hook

    req = request(server_url, headers=["Range"=>"bytes=0-100"], throw = false,
                  output = mktemp()[2], downloader = downloader)

    return is_partial_download(req)
    
end

supports_range_header(server_url::URIs.URI) = supports_range_header("$server_url")

is_ssl_error(r::RequestError) = r.code == 35

is_ssl_error(r::Response) = false

is_partial_download(r::RequestError) = r.response.status == 206

is_partial_download(r::Response) = r.status == 206

has_accept_ranges(r::RequestError) = get(Dict(r.response.headers), "accept-ranges", "none") == "bytes"

has_accept_ranges(r::Response) = get(Dict(r.headers), "accept-ranges", "none") == "bytes"

make_insecure(s::String) = "$(make_insecure(URI(s)))"
make_insecure(u::URI) = begin

    if u.scheme[end] != 's' return u end

    URI(scheme = u.scheme[1:end-1], host = u.host, port = u.port,
        fragment = u.fragment, path = u.path, userinfo = u.userinfo,
        query = u.query)
    
end
