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
using ZipArchives

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
export ping_archive    #Check that URL exists
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
include("kcd_image.jl")
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
A tar archive has been created by tar and then optionally
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

struct HDFArchive <: BulkArchive
    local_cache::AbstractString
    original_url::URI
end

HDFArchive(local_dir, all_urls::Vector{URI}) = begin

    HDFArchive(local_dir, get_constant_part(all_urls))
        
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
    if arch_type === nothing
        arch_type = unique(c["$cat_name.format"])
    end
    
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

    arch_list = []

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

    if arch_type == "HDF5"
        const_part = "$(get_constant_part(u))"
        if const_part in keys(subs)
            @debug "Found local sub:" subs[const_part]
            base_dir = subs[const_part]
        else
            @debug "$const_part not in subs list" keys(subs)
            base_dir = base_dirs[1]
        end

        return [HDFArchive(base_dir, u)]
    
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
        if !startswith("$full_uri", aurl)
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

local_equivalent(a::HDFArchive, ext_info) = begin

    c = get_prefix(ext_info)
    u = URI(getproperty(ext_info,"$(c)uri"))
    base_part = "$(a.original_url)"
    unique_part = "$u"[length(base_part)+2:end]
    base = joinpath(a.local_cache, unique_part)
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
    peek_image(a::RsyncArchive, cif_block::CifContainer; entry_no=0, check_name=true, extract=true)

Find the name of an image in `a`.
"""
peek_image(a::RsyncArchive, cif_block::CifContainer; kwargs...) = begin

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

peek_image(a::HDFArchive, cif_block::CifContainer; max_down = 2e7, kwargs...) = begin

    ei = first(get_loop(cif_block, "_array_data_external_data.id"))
    
    if !has_local_version(a, ei)
        
        cmd_list = Cmd[]
        uri = getproperty(ei,"_array_data_external_data.uri")

        r = request(uri)
        fs = get_file_size(r)
        
        if fs == 0

            throw(error("Cannot determine file size for $uri. Download manually and use -l or -s options"))
        elseif fs > max_down

            throw(error("$uri size $fs exceeds maximum download limit of $max_down, please use -m option or download manually and use -l / -s options"))
        end
        
        @debug "Downloading $(local_equivalent(a, ei)) from $uri"
        
        run(Cmd(`curl -s $uri --output $(local_equivalent(a, ei))`))
    end
    
    return local_equivalent(a, ei)
end

peek_image(a::ZipArchive, cif_block::CifContainer; entry_no=0, max_down = 2.0e7) = begin

    # Find something already present if possible
    
    p = get_any_local(a, cif_block)
    if !isnothing(p) return p end
    
    # Get entire Zip file
    
    ei = first(get_loop(cif_block, "_array_data_external_data.id"))
    uri = getproperty(ei,"_array_data_external_data.uri")

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
    path_in_arch = getproperty(ei, "_array_data_external_data.archive_path")
    if !(path_in_arch in zip_names(zr))
        throw(error("Bad archive $uri: missing $path_in_arch"))
    end

    in_mem = zip_readentry(zr, path_in_arch)
    fo = open(local_equivalent(a, ei), "w")
    write(fo, in_mem)
    close(fo)

    return local_equivalent(a, ei)
    
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

end
