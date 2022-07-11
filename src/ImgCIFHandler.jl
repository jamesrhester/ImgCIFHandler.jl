# Provide methods for accessing imgCIF data

#==

# Background

imgCIF is a text format for describing raw data images collected on
crystallographic instruments. As such it is not suited for actually
containing the raw data, but instead provides tags pointing to the
data storage location.

This library provides methods for bringing that image data into
Julia, including handling the variety of formats provided. The only
function exported is `imgload`.

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

export imgload         #Load raw data
export peek_image      #Find first image in archive
export make_absolute_uri #Use Cif block contents to make absolute URI
export get_detector_axis_settings #Get axis settings for particular frame
export get_beam_centre
export get_gonio_axes
export get_pixel_coordinates
export get_surface_axes #get the axes used to locate pixels on the detector
export get_id_sequence #get a list of sequential binary ids from the same scan

include("hdf_image.jl")
include("cbf_image.jl")
include("adsc_image.jl")
include("imgcif.jl")

get_image_ids(c::CifContainer) = begin
    return c["_array_data.binary_id"]
end
    
"""
    imgload(c::Block,array_ids)

Return the image referenced in CIF Block `c` corresponding to the specified raw array identifiers.
`local_version` gives local copies for URLs listed in `c`. On failure returns a string containing
an error message.
"""
imgload(c::CifContainer,bin_ids;local_version=Dict()) = begin

    dl_info = external_specs_from_bin_ids(bin_ids,c)
    
    local_copy =  get(local_version,"$(dl_info[1,:uri])", nothing)
    @debug "Loading image from" local_copy
    
    imgload(dl_info,local_copy = local_copy)
end

"""
    imgload(img_info::DataFrame,local_copy=nothing)

Return the raw 2D data specified by the information in
`img_info`. `local_copy` is a local copy of the URI referenced
in `img_info`, if present.
"""
imgload(ext_info::DataFrame;kwargs...) = begin
    # May switch later to native if we can get Tar to terminate early
    imgload_os(ext_info;kwargs...)
end

imgload_native(uri::URI,format::Val;arch_type=nothing,arch_path=nothing,file_compression=nothing,kwargs...)= begin
    # Set up input stream
    stream = BufferStream()
    have_file = false   #changes to true when file found
    loc = mktempdir() #Where the final file is found
    # Parse the URI to catch local files
    u = URI(uri)

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
    imgload(endloc,format;kwargs...)
end

"""
    imgload_os(ext_info::DataFrame;local_copy = nothing)

Return an image obtained by summing a series of images corresponding to the entries in `ext_info`
"""
imgload_os(ext_info::DataFrame;local_copy = nothing) = begin

    # Practical restrictions: only one uri and archive type

    if length(unique(ext_info.full_uri)) > 1
        throw(error("Multiple images must be at the same URI"))
    end

    uri = ext_info.full_uri[1]

    # Use OS pipelines to download efficiently

    cmd_list = Cmd[]
    loc = mktempdir()
    decomp_option = "-v"

    arch_type = "archive_format" in names(ext_info) ? ext_info.archive_format[1] : nothing
    arch_paths = "archive_path" in names(ext_info) ? ext_info.archive_path : nothing
    
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
        fc = ext_info.file_compression[1]
    else
        fc = nothing
    end
    
    temp_locals = map(temp_locals) do tl
        if !isnothing(fc)
            final_file = open(joinpath(loc,tl*"final"),"w")
            if fc == "GZ"
                run(pipeline(`zcat $tl`,final_file))
            elseif fc == "BZ2"
                run(pipeline(`bzcat $tl`,final_file))
            end
            close(final_file)
            tl*"_final"
        else
            if tl != local_copy #ok to move
                mv(tl,joinpath(loc,tl*"_final"))
                tl*"_final"
            else
                local_copy
            end
        end
    end

    @debug "Final files for reading are" temp_locals

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

imgload(c::CifContainer,frame::Int;scan=nothing,diffrn=nothing) = begin
    println("Not implemented yet")
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

frame_from_frame_id(c::CifContainer,frame::String,scan,diffrn) = begin
    # The frame number is provided in _diffrn_scan_frame.frame_number
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
