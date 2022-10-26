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

export imgload         #Load raw data
export peek_image      #Find first image in archive
export make_absolute_uri #Use Cif block contents to make absolute URI
export get_detector_axis_settings #Get axis settings for particular frame
export get_beam_centre
export get_gonio_axes
export get_pixel_coordinates
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

include("hdf_image.jl")
include("cbf_image.jl")
include("adsc_image.jl")
include("imgcif.jl")
include("recip.jl")

get_image_ids(c::CifContainer) = begin
    return c["_array_data.binary_id"]
end
    
"""
    imgload(c::Block,array_ids;local_version=Dict())

Return the image referenced in CIF Block `c` corresponding to the specified raw array identifiers.
`local_version` gives local copies for URLs listed in `c`. `cached` is a list of files that are
already available locally as a dictionary `{URI, path}=>local_path`. On failure returns 
a string containing an error message.
"""
imgload(c::CifContainer,bin_ids;local_version=Dict(), cached = Dict()) = begin

    dl_info = external_specs_from_bin_ids(bin_ids,c)
    
    local_copy =  get(local_version,"$(dl_info[1,:uri])", nothing)
    @debug "Loading image from" local_copy
    
    imgload(dl_info,local_copy = local_copy, cached = cached)
end

"""
    imgload(ext_info::DataFrame;local_copy = nothing, cached=cached)

Return an image obtained by summing a series of images corresponding to the entries in `ext_info`.
If local_copy is a file, it is used as the archive source instead of the requested
URL.
"""
imgload(ext_info::DataFrame;local_copy = nothing, cached= Dict()) = begin

    # Practical restrictions: only one uri and archive type

    uri = unique(ext_info.full_uri)
    if length(uri) > 1
        throw(error("Multiple images must be at the same URI"))
    end

    uri = uri[]
    arch_paths = "archive_path" in names(ext_info) ? ext_info.archive_path : nothing

    # Check cache

    not_in_cache = []
    if !isnothing(arch_paths)
        not_in_cache = filter( x-> !haskey(cached, ("$uri", x)) || cached[("$uri",x)] == "" , arch_paths)
    else
        not_in_cache = !haskey(cached,("$uri",nothing))
    end

    @debug "Not found in cache" not_in_cache uri
    
    if length(not_in_cache) > 0

        @debug "Not all requested files cached, downloading"
        temp_locals = download_images_os(uri, ext_info, local_copy, arch_paths)
    else
        @debug "Using cached files" arch_paths
        if arch_paths == nothing
            temp_locals = [cached[("$uri",nothing)]]
        else
            temp_locals = map(ap -> cached[("$uri", ap)], arch_paths)
        end
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

imgload(c::CifContainer,scan_id,frame_no::Int;kwargs...) = begin
    bin_ids = bin_id_from_scan_frame(c,scan_id,frame_no)
    imgload(c,bin_ids;kwargs...)
end

"""
    imgload(c::CIF)

Return the image referenced by the first encountered `_array_data.binary_id` in the
first block of CIF file `c`.
"""
imgload(c::Cif) = begin
    b = first(c).second
    f_id = b["_array_data.binary_id"][1]
    imgload(b,f_id)
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
    if arch_type in ("TGZ","TBZ","TAR")
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
    peek_image(URI,archive_type,cif_block::CifContainer;entry_no=0)

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
