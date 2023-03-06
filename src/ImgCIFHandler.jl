# Provide methods for accessing imgCIF data

#==

# Background

imgCIF is a text format for describing raw data images collected on
crystallographic instruments. As such it is not suited for actually
containing the raw data, but instead provides tags pointing to the
data storage location.

This library provides methods for bringing that image data into
Julia, including handling the variety of formats provided, and for
general manipulation of the information found in the imgCIF.

==#
"""
Handle items in an imgCIF file. See method `imgload`.
"""
module ImgCIFHandler

using CrystalInfoFramework
using FilePaths
using DataFrames
using InvertedIndices
import Tar
using CodecBzip2
using CodecZlib
using Downloads

# See format-specific includes for more using statements

using TranscodingStreams
using URIs
using SimpleBufferStream
using LinearAlgebra
using ImageBinarization
using ImageFiltering
using Statistics
using Rotations

export create_archives  #register an archive
export imgload         #Load raw data
export peek_image      #Find first image in archive
export make_absolute_uri #Use Cif block contents to make absolute URI
export get_detector_axis_settings #Get axis settings for particular frame
export get_beam_centre
export get_gonio_axes
export get_pixel_coordinates
export get_detector_distance   #Get distance of flat detector from crystal
export get_surface_axes #get the axes used to locate pixels on the detector
export get_id_sequence #get a list of sequential binary ids from the same scan
export get_axis_vector #get the vector for an axis_id
export find_peaks #Find some peaks in an image
export bin_id_from_scan_frame #Convert scan/frame_no into binary id
export scan_frame_from_bin_id #Convert bin_id into scan/frame
export peak_to_frames #Calculate all peak appearances
export get_dependency_chain #Get the dependent axes of an axis

export Peak
export intensity, coords, frame, scan, dist #Working with peaks

export ImageArchive
export has_local_version #for ImageArchive types
export get_constant_part #for RsyncArchive work

include("hdf_image.jl")
include("cbf_image.jl")
include("adsc_image.jl")
include("imgcif.jl")
include("recip.jl")

get_image_ids(c::CifContainer) = begin
    return c["_array_data.binary_id"]
end

# Archive types
"""
In order to access images, an archive must be specified,
and then supplied to the `imgload` routine. To avoid
re-downloading on each call, this archive should be
created separately.

An archive is either remote, or local. A local archive may
have a remote equivalent, such that it is a cached version.
We provide routines to convert from remote to local, which
involves downloading selected files.
"""
abstract type ImageArchive end

"""
A `BulkArchive` consists of a single file containing all image
frames.
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
A tar archive has been created by tar and then optionall
compressed. The type parameter T refers to the compression used in
the tar archive. To add a new compression type add it
to the `decomp_option` method.
"""
struct TarArchive{T} <: CompressedTarArchive
    local_cache::AbstractString
    original_url::URI
end

struct ZipArchive <: CompressedArchive
    local_cache::AbstractString
    original_url::URI
end

create_archives(c::Cif; kwargs...) = create_archive(first(c).second; kwargs...)

"""
    create_archives(c::CifContainer; subs = Dict())

Create a series of ImageArchive objects, one for each distinct location
in `c`. Note that a distinct location is not necessarily a distinct URL,
as rsync: and file: URLs may be unique for each frame but correspond
to a single directory.
"""
create_archives(c::CifContainer; subs = Dict()) = begin

    cat_name = "_array_data_external_data"
    arch_type = haskey(c, "$cat_name.archive_format") ? unique(c["$cat_name.archive_format"]) : nothing
    
    dl_info = unique(c["$cat_name.uri"])
    test_uris = URI.(dl_info)

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

    arch_list = []

    uq_u = unique(u)

    base_dirs = map(uq_u) do one_uri
        if one_uri in keys(subs)
            @debug "Found local sub" subs[one_uri]
            subs[one_uri]
        else
            mktempdir()
        end
    end

    if arch_type in ("TBZ","TGZ","TXZ","TAR")
        for (one_uri, bd) in zip(uq_u, base_dirs)
            push!(arch_list, TarArchive{Val(Symbol(arch_type))}(bd, one_uri))
        end
        return arch_list
    end

    if arch_type == "ZIP"
        for (one_uri, bd) in zip(uq_u, base_dirs)
            push!(arch_list, ZipArchive(bd, one_uri))
        end
        return arch_list
    end
    
    # "file" schemes pointing to compressed archives are covered by
    # the types above. Do not mix schemes.

    if first(u).scheme == "file" || first(u).scheme == ""
        return [LocalArchive(root_dir)]
    end

    if first(u).scheme == "rsync"
        const_part = "$(get_constant_part(u))"
        if const_part in keys(subs)
            @debug "Found local sub:" subs[const_part]
            base_dir = subs[const_part]
        else
            @debug "$const_part not in subs list" keys(subs)
            base_dir = base_dirs[1]
        end
        
        return [RsyncArchive(base_dir, u)]
    end
    
    throw(error("Unrecognised archive type or URI scheme $arch_type, $u"))
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

"""
An RsyncArchive corresponds to a dataset for which each
frame is addressable through a unique rsync: URL. 
`original_url` in this case corresponds to the longest
unique part of the URL, and this is replaced by
`local_directory` in the local filesystem.
"""
struct RsyncArchive <: AddressableArchive
    local_directory::AbstractString
    original_url::URI
end

RsyncArchive(local_dir, all_urls::Vector{URI}) = begin

    RsyncArchive(local_dir, get_constant_part(all_urls))
        
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
        push!(cmd_list, Cmd(`curl -s $(ext_info.uri[1])`,ignorestatus=true))
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
    download_images_os(a::RsyncArchive, ext_info)

As rsync takes care of checking if the local version is equivalent,
we do not need to keep track ourselves. But we do to save rsync
calls.
"""
download_images_os(a::RsyncArchive, ext_info) = begin

    c = get_prefix(ext_info)

    # SBGRID does not download the full directory hierarchy
    # from the URL, therefore we insert a "." so that we can
    # control it

    aurl = "$(a.original_url)"
    for r in eachrow(ext_info)
        full_uri = getproperty(r,"$(c)uri")
        if !startswith(full_uri, aurl)
            throw(error("Asked to download from $full_uri but archive is for $aurl"))
        end
        download_uri = joinpath(aurl,".",full_uri[(length(aurl)+2):end])
        @debug "Rsync download address" download_uri
        run(`rsync -avR $download_uri $(a.local_directory)`)

        #verify
        @assert ispath(local_equivalent(a, r))
    end

end

download_images_os(a::LocalArchive, ext_info) = begin
    for r in eachrow(ext_info)
        @assert ispath(local_equivalent(a, r))
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
    URI(scheme="rsync", host = b.host,
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
    local_equivalent(a::RsyncArchive, ext_info)

The local file name corresponding to the information provided in `ext_info`.
Note we use rsync -avR to preserve the directory hierarchy.
"""
local_equivalent(a::RsyncArchive, ext_info) = begin

    c = get_prefix(ext_info)
    u = URI(getproperty(ext_info,"$(c)uri"))
    base_part = "$(a.original_url)"
    unique_part = "$u"[length(base_part)+2:end]
    base = joinpath(a.local_directory, unique_part)
    if "$(c)file_compression" in names(ext_info)
        return base * "_final"
    else
        return base
    end
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

imgload(c::CifContainer, bin_ids, a::Vector) = begin

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

imgload(c::CifContainer, scan_id, frame_no::Int, a::ImageArchive) = begin
    bin_ids = bin_id_from_scan_frame(c, scan_id, frame_no)
    imgload(c, bin_ids, a)
end

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
    imgload(cc, first(a))
end

"""
    imgload(c::CifContainer, a::ImageArchive)

Return the image referenced by the first encountered `_array_data.binary_id` in
CIF block `c`. `a` is an archive created using `create_archive`.
"""
imgload(c::CifContainer, a::ImageArchive) = begin
    f_id = c["_array_data.binary_id"][1]
    imgload(c, f_id, a)
end

imgload(p::AbstractPath) = begin
    imgload(Cif(p,native=true))
end

imgload(s::AbstractString) = begin
    imgload(Path(s))
end
     
download_images_os(uri, ext_info, local_copy, arch_paths)  = begin

    # Use OS pipelines to download efficiently

    arch_type = "archive_format" in names(ext_info) ? ext_info.archive_format[1] : nothing

    cmd_list = Cmd[]
    loc = mktempdir()
    decomp_option = "-v"

    if arch_type == "TGZ" decomp_option = "-z" end
    if arch_type == "TBZ" decomp_option = "-j" end
    if arch_type == "TXZ" decomp_option = "-J" end
    if arch_type in ("TGZ","TBZ","TAR","TXZ")
        if local_copy == nothing
            push!(cmd_list, Cmd(`curl -s $uri`,ignorestatus=true))
        else
            push!(cmd_list, `cat $local_copy`)
        end
        push!(cmd_list, `tar -C $loc -x $decomp_option -f - --occurrence $arch_paths`)
        temp_locals = joinpath.(Ref(loc),arch_paths)
    else
        if local_copy == nothing
            temp_locals = [joinpath(loc,"temp_download")]
            push!(cmd_list, `curl $uri -o $(temp_locals[1])`)
        else
            temp_locals = [local_copy]
        end
    end
    @debug "Command list" cmd_list
    if length(cmd_list) > 0
        try
            run(pipeline(cmd_list...))
        catch exc
            @debug "Finished downloading" exc
        end
    end

    # Now the final files are listed in $temp_locals

    if arch_type == "ZIP"   #has been downloaded to local storage
        run(`unzip $(temp_locals[1]) $arch_paths -d $loc`)
        temp_locals = joinpath.(Ref(loc),arch_paths)
    end

    if "file_compression" in names(ext_info)
        temp_locals = decompress_frames(ext_info, loc, temp_locals)
    else
        temp_locals = map(temp_locals) do tl
            if isnothing(local_copy)
                tl
            elseif tl != local_copy && isfile(local_copy) #ok to move
                mv(tl,joinpath(loc,tl*"_final"))
                tl*"_final"
            elseif isfile(local_copy)
                local_copy
            else
                tl
            end
        end
    end

    @debug "Final files for reading are" temp_locals

    return temp_locals
end

"""
UNTESTED. Wasn't working properly last time it was checked.
"""
download_images_native(uri, ext_info, local_copy, arch_paths) = begin
    # Set up input stream
    stream = BufferStream()
    have_file = false   #changes to true when file found
    loc = mktempdir() #Where the final file is found
    # Parse the URI to catch local files
    u = URI(uri)
    arch_type = "archive_format" in names(ext_info) ? ext_info.archive_format[1] : nothing


    # Begin asynchronous section. Thanks to Julia Discourse for the technique!
    @sync begin
        @async try
            Downloads.download("$uri",stream;verbose=true)
        catch exc
            if !have_file
                @error "Problem downloading $uri" exc
            end
        finally
            close(stream)
        end

        decomp = stream
        if !(arch_type in (nothing,"TAR"))
            if arch_type == "TGZ"
                decomp = GzipDecompressorStream(stream)
            elseif arch_type == "TBZ"
                decomp = Bzip2DecompressorStream(stream)
            end
        end
        # Now handle having an internal directory structure
        if arch_path != nothing
            full_path = joinpath(loc,arch_path)
            if arch_type in ("TGZ","TBZ","TAR")
                # callback to abort after reading
                abort_callback(x) = begin
                    @info("Found: $have_file")
                    if have_file == true
                        if ispath(full_path) && stat(full_path).size > 0
                            @info "$(stat(full_path))"
                            cp(full_path,full_path*"_1")
                            throw(error("Extracted one file to $full_path"))
                        else
                            @info("Can't yet see file at path $full_path or size=0")
                        end     
                    end
                    if x.path == arch_path
                        @info("Extracting", x)
                        have_file = true
                        return true
                    else
                        @info("Ignoring", x)
                        return false
                    end
                end

                @async try
                    Tar.extract(abort_callback,decomp,loc)
                catch exc
                    if !have_file || !ispath(full_path) || isopen(full_path)
                        @info("File at $full_path is $(stat(full_path))")
                        @error "Untar problem" exc
                    end
                finally
                    loc = joinpath(loc,arch_path)
                    close(decomp)
                end
                
            elseif arch_type == "ZIP"
                @async begin
                    w = ZipFile.Reader(stream)
                    loc,final_file = mktemp()
                    for f in w.files
                        if f == arch_path
                            write(final_file,read(f))
                            close(final_file)
                            break
                        end
                    end
                end
            end
        else
            loc,final_file = mktemp()
            @async begin
                count = write(final_file,read(stream))
                println("$count bytes read")
                close(final_file)
            end
        end
    end   #all @asyncs should finish before proceeding
    #
    println("Extracted file to $loc")
    # Apply any final decompression
    endloc = loc
    fdecomp = nothing
    if file_compression == "GZ"
        fdecomp = GzipDecompressor()
    elseif file_compression == "BZ2"
        fdecomp = Bzip2Decompressor()
    end
    if fdecomp != nothing
        endloc,unc_file = mktemp()
        out_str = TranscodingStream(fdecomp,open(loc,"r"))
        write(unc_file,out_str)
        close(unc_file)
        println("Decompressed file is $endloc")
    end
    return [endloc]
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

## Some utility files

"""
    list_archive(uri;n=5)

Given an archive `uri`, list the first `n` members. A negative
number for `n` lists all members. Uses `SimpleBufferStream`.
"""
list_archive(u::URI;n=5,compressed=nothing) = begin
    counter = 1

    # Set up chain of streams

    dldstream = BufferStream()
    unzipstream = BufferStream()
    untarstream = BufferStream()

    # Our header information

    hdrs = nothing

    function task_chain(in_str,out_str)
        @async begin
            while !eof(in_str)
                write(out_str, readavailable(in_str))
            end
            close(out_str)
        end
    end
    
    function do_dld(out_stream)
        @async begin
            try
                Downloads.download("$u",outstream;verbose=true)
            catch exc
                if counter < n
                    @error "Problem downloading $uri" exc
                end
            finally
                close(outstream)
            end
        end
    end

    if !(compressed in (nothing,"TAR"))
        if compressed == "TGZ"
            decomp = GzipDecompressorStream(in_stream)
        elseif compressed == "TBZ"
            decomp = Bzip2DecompressorStream(in_stream)
        end
    end        

    # Now handle having an internal directory structure
    
    if compressed in ("TGZ","TBZ","TAR")

            # callback to abort after listing

            abort_callback(x) = begin
                counter = counter + 1
                if if n > 0 && counter > n
                    throw(error("Made it to $n"))
                end
                @info(x)
                return true
            end

            @async try
                hdrs = Tar.list(abort_callback,decomp)
            catch exc
                if counter < n
                    @error "Untar problem" exc
                end
            finally
                close(decomp)
            end
        end
    end
    return hdrs
end

"""
    peek_image(a::ImageArchive, cif_block::CifContainer; entry_no=0, check_name=true)

Find the name of an image in `a`, searching for `entry_no`, and
check if this image is listed in `cif_block` if `check_name`
is true.
"""
peek_image(a::LocalArchive, cif_block::CifContainer; entry_no=0, check_name=true) = begin

    # For a local archive, some files are automatically present locally; all we need to do is
    # find the first one in cif_block that matches a local file

    catname = "_array_data_external_data"
    if haskey(cif_block,"$catname.id")
        for r in eachrow(get_loop(cif_block, "$catname.id"))
            if has_local_version(a, r)
                return local_equivalent(a, r)
            end
        end
    end
    
end

"""
    peek_image(a::RsyncArchive, cif_block::CifContainer; entry_no=0, check_name=true, extract=true)

Find the name of an image in `a`, searching for `entry_no`, and
check if this image is listed in `cif_block` if `check_name`
is true. If `extract` the image is also retrieved to the local
cache.
"""
peek_image(a::RsyncArchive, cif_block::CifContainer; entry_no=0, check_name=true, extract=true) = begin

    # For an Rsync archive, some files are automatically present
    # locally; all we need to do is find the first one in cif_block
    # that matches a local file

    c = "_array_data_external_data"
    if haskey(cif_block,"$c.id")
        r = first(get_loop(cif_block, "$c.id"))
        if !has_local_version(a, r)
            aurl = "$(a.original_url)"
            full_uri = getproperty(r, "$(c).uri")
            download_uri = joinpath(aurl, ".", full_uri[(length(aurl)+2):end])
            rscmd = Cmd(`rsync -avR $download_uri $(a.local_directory)`)
            @debug "Obtaining file..." rscmd
            run(rscmd)
        end
        return local_equivalent(a, r)
    end
    
end

peek_image(a::TarArchive, cif_block::CifContainer; entry_no=0, check_name=true, extract=true) = begin

    #TODO: detect already-present images

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

"""
    peek_image(URI,archive_type,cif_block::CifContainer;entry_no=0, check_name=true)

Find the name of an image in archive of type `archive_type` at `URL`, searching
from entry number `entry_no`,and check that this image is available in `cif_block`
if `check_name` is true.
"""
peek_image(uri::URI,arch_type,cif_block::CifContainer;entry_no=0,check_name=true) = begin

    cmd_list = Cmd[]
    
    if arch_type == "ZIP"
        throw(error("Peeking into file not supported for ZIP"))
    end
    
    decomp_option = ""
    if arch_type == "TGZ" decomp_option = "-z"
    elseif arch_type == "TBZ" decomp_option = "-j"
    elseif arch_type == "TXZ" decomp_option = "-J"
    else throw(error("Unrecognised archive type $arch_type"))
    end

    push!(cmd_list, Cmd(`curl -s $uri`,ignorestatus=true))
    push!(cmd_list, Cmd(`tar -t -v $decomp_option -f -`,ignorestatus=true))
    awkstr1 =  "\$3 > 0 { print \$NF }"
    awkstr2 =  "\$3 > 0 && FNR >= $entry_no { exit }"
    push!(cmd_list, `awk -e $awkstr1 -e $awkstr2`)

    @debug "Peeking into $uri for $arch_type starting at $entry_no"
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
    if fname != nothing && check_name
        if haskey(cif_block,"_array_data_external_data.archive_path")
            pos = indexin([fname],cif_block["_array_data_external_data.archive_path"])[]
            if pos != nothing && cif_block["_array_data_external_data.uri"][pos] == "$uri"
                return fname
            end
        end
        return nothing
    end
    return fname
end

peek_image(u::URI,arch_type) = peek_image(u,arch_type,Block{String}())  #for testing

end
