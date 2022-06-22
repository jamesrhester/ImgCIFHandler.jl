#
#  Auto-install all required items
#
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ImgCIFHandler#main
using ImageInTerminal, Colors,ImageContrastAdjustment
using ArgParse
using CrystalInfoFramework,FilePaths,URIs, Tar, Images

# Tests for imgCIF files
const test_list = []
const test_list_with_img = []
const test_full_list = []
const dictionary_checks = []

# These macros will become more sophisticated
# with time to make a nice printout. For now
# they simply collect the tests into three
# lists.
macro noimgcheck(description, check)
    # extract name
    #func_name = expr.args[1].args[1] #:=->:call
    #push!(test_list,func_name)
    quote
        x = $(esc(check))
        push!(test_list,($(esc(description)),x))
    end
end

macro imgcheck(description, check)
    quote
        x = $(esc(check))
        push!(test_list_with_img,($(esc(description)),x))
    end
end

macro fullcheck(description, check)
    quote
        x = $(esc(check))
        push!(test_full_list,($(esc(description)), x))
    end
end

#macro plaincheck(testname,testexpr,message)
#end

# Checks that do not require an image

@noimgcheck "Required items" required_items(incif) = begin
    messages = []
    info_check = ("_array_structure_list.dimension",
                  "_array_structure_list.index",
                  "_array_structure_list.precedence",
                  "_array_structure_list.direction",
                  "_array_structure_list_axis.displacement_increment",
                  "_array_structure_list_axis.displacement",
                  "_diffrn_scan_axis.axis_id",
                  "_diffrn_scan_axis.angle_start",
                  "_diffrn_scan_axis.displacement_start",
                  "_diffrn_scan_axis.angle_range",
                  "_diffrn_scan_axis.displacement_increment",
                  "_diffrn_scan_axis.displacement_range",
                  "_diffrn_scan_axis.angle_increment",
                  "_diffrn_scan.frames",
                  "_axis.depends_on",
                  "_axis.type",
                  "_axis.vector[1]",
                  "_axis.vector[2]",
                  "_axis.vector[3]",
                  "_axis.offset[1]",
                  "_axis.offset[2]",
                  "_axis.offset[3]",
                  "_diffrn_radiation_wavelength.id",
                  "_diffrn_radiation_wavelength.value",
                  "_diffrn_radiation.type",
                  "_array_data.array_id",
                  "_array_data.binary_id")
    for ic in info_check
        if !haskey(incif,ic)
            push!(messages,(false,"Required item $ic is missing"))
        end
    end
    return messages
end

# Make sure there is a data specification
@noimgcheck "Data source" data_source(incif) = begin
    messages = []
    if !haskey(incif,"_array_data.data") && !haskey(incif,"_array_data_external_data.id")
        push!(messages,(false,"No source of image data specified"))
        return messages
    end
    if haskey(incif,"_array_data.data")
        push!(messages,(true,"WARNING:raw data included in file, processing will be slow"))
    end
    p = URI(incif["_array_data_external_data.uri"][1])
    if p.scheme == "file" || p.scheme == nothing
        push!(messages,(true,"WARNING: external data stored in local file system, this is not portable"))
    end
    if resolvereference("file:///dummy",p) != p
        push!(messages,(true,"WARNING: relative URI $p provided, this is not portable"))
    end
    return messages
end

@noimgcheck "Axes defined" axes_defined(incif) = begin
    all_axes = vcat(incif["_axis.id"],[nothing])
    messages = []
    test_values = ("_axis.depends_on","_diffrn_scan_frame_axis.axis_id",
                   "_diffrn_scan_axis.axis_id",
                   "_diffrn_measurement_axis.axis_id",
                   "_diffrn_detector_axis.axis_id",
                   "_array_structure_list_axis.axis_id")
    for tv in test_values
        if !haskey(incif,tv) continue end
        unk = setdiff(skipmissing(incif[tv]),all_axes)
        if length(unk) > 0
            push!(messages,(false,"Undefined axes in $tv: $unk"))
        end
    end
    return messages
end

@noimgcheck "Our limitations" our_limitations(incif) = begin
    messages = []
    if haskey(incif,"_array_data.array_id") && length(unique(incif["_array_data.array_id"])) > 1
        push!(messages,(false,"WARNING: cannot currently correctly check files with more than one data array structure"))
    end
    if haskey(incif,"_diffrn_detector.id") && length(unique(incif["_diffrn_detector.id"])) > 1
        push!(messages,(false,"WARNING: cannot currently correctly check files with more than one detector"))
    end
    if haskey(incif,"_diffrn_detector_element.id") && length(unique(incif["_diffrn_detector_element.id"])) > 1
        push!(messages(false,"WARNING: cannot currently correctly check files with more than one detector element"))
    end
    if haskey(incif,"_diffrn.id") && length(incif["_diffrn.id"]) > 1
        push!(messages(false,"WARNING: cannot currently correctly check files with more than one set of diffraction conditions (more than one diffrn.id)"))
    end
    return messages
end

# Check that the detector translation is as expected
@noimgcheck "Detector translation" trans_is_neg(incif) = begin
    messages = []
    # find the detector translation based on the idea that it will
    # be a translation axis directly dependent on 2 theta
    axes = get_loop(incif,"_axis.id")
    tt = filter(axes) do r
        getproperty(r,"_axis.equipment")==    "detector" && getproperty(r,"_axis.type")==     "rotation" &&  getproperty(r,"_axis.depends_on")==nothing
    end
    if size(tt,1) != 1
        return [(true,"Warning: can't identify two theta axis $tt")]
    end
    axname = getproperty(tt,"_axis.id")[]
    det = filter(axes) do r
            getproperty(r,"_axis.equipment") == "detector" && getproperty(r,"_axis.type") == "translation" && getproperty(r,"_axis.depends_on") == axname
    end
    if size(det,1) != 1
        detnames = getproperty(det,"_axis.id")
        return [(true,"Warning: can't identify detector translation axis, have $detnames")]
    end
    det = first(det)
    detname = getproperty(det,"_axis.id")
    # check that translation is negative Z
    signv = sign(parse(Float64,getproperty(det,"_axis.vector[3]")))
    signo = signv*parse(Float64,getproperty(det,"_axis.offset[3]"))
    if signo == 1
        push!(messages,(false,"Detector translation $detname is positive"))
    end
    if signv == 1
        push!(messages,(false,"Detector translation axis $detname points towards the source"))
    end
    av1 = parse(Float64,getproperty(det,"_axis.vector[1]"))
    av2 = parse(Float64,getproperty(det,"_axis.vector[2]"))
    if av1 != 0 || av2 != 0
        push!(messages,(false,"Detector translation is not parallel to beam"))
    end
    return messages
end

@noimgcheck "Scan range" scan_range(incif) = begin
    # Check that scan ranges are correct
    # The scan range is the total range from start of first step
    # to the end of the final step.
    
    messages = []

    cn = "_diffrn_scan_axis."  #for brevity

    # account for unlooped case
    
    scan_loop = get_loop(incif,cn*"axis_id")
    if size(scan_loop,1) == 0 && haskey(incif,cn*"axis_id")
        loop_names = [cn*"axis_id",cn*"displacement_range",cn*"angle_range",cn*"angle_start",
                      cn*"angle_increment",cn*"displacement_increment",cn*"displacement_start"]
        for n in loop_names
            if !haskey(incif,n) incif[n] = [missing] end
        end
        create_loop!(incif,loop_names)
        scan_loop = get_loop(incif,cn*"axis_id")
    end

    actual = filter(scan_loop) do r
        !ismissing(getproperty(r,cn*"displacement_range")) || !ismissing(getproperty(r,cn*"angle_range"))
    end
    
    if size(actual,2) == 0
        return [(false,"No scan range has been specified")]
    end

    for or in eachrow(actual)
        if ismissing(getproperty(or,cn*"displacement_range"))
            stem = "angle_"
        else
            stem = "displacement_"
        end
        
        start = parse(Float64,getproperty(or,cn*stem*"start"))
        stepsize = parse(Float64,getproperty(or,cn*stem*"increment"))
        range = parse(Float64,getproperty(or,cn*stem*"range"))
        if stepsize == 0
            if range == 0
                continue
            else
                push!(messages,(false,"Non-zero range with zero stepsize for $or"))
                continue
            end
        end
        
        numsteps = range/stepsize + 1

        # Look up how many steps should be there

        if cn*"scan_id" in names(scan_loop)
            scan_id = getproperty(or,cn*"scan_id")
            scan_info = incif[Dict("_diffrn_scan.id"=>scan_id)]
            nosteps = parse(Float64,scan_info[1,"_diffrn_scan.frames"])
        else
            scan_id = "Scan01"
            nosteps = parse(Float64,incif["_diffrn_scan.frames"][])
        end

        # And now check
        if abs(nosteps - numsteps) > 0.3
            push!(messages,(false,"Range/increment do not match number of steps $nosteps for scan $scan_id, expected $numsteps"))
        else
            push!(messages,(true,"Range/increment match number of steps $nosteps for scan $scan_id (expected $numsteps)"))
        end
    end
    return messages
end

@noimgcheck "All frames present" check_frames(incif) = begin
    # Check that all frames in a scan have source data
    messages = []
    basename = "_diffrn_scan."

    # Get expected number of frames per scan

    if !haskey(incif,basename*"frames")
        return messages
    end
    
    numsteps = parse.(Int64,incif[basename*"frames"])
    if length(numsteps) > 1 || haskey(incif,basename*"id")
        scan_names = incif[basename*"id"]
    else
        scan_names = ["SCAN01"]
    end
    
    # Check that they are all there

    bname = "_diffrn_scan_frame."
    if length(numsteps) == 1 && !haskey(incif,bname*"scan_id")
        incif[bname*"scan_id"] = ["SCAN01"]
        create_loop!(incif,[bname*"scan_id",bname*"frame_number"])
    end

    scan_info = get_loop(incif,"_diffrn_scan_frame.scan_id")
    
    for (one_scan,nsteps) in zip(scan_names,numsteps)
        os = filter(scan_info) do r
            getproperty(r,"_diffrn_scan_frame.scan_id") == one_scan
        end
        fnumbers = parse.(Int64,os[!,"_diffrn_scan_frame.frame_number"])
        ufn = unique(fnumbers)
        if minimum(ufn) != 1
            push!(messages,(false,"First frame not specified for scan $one_scan"))
            continue
        end
        if maximum(ufn) != nsteps
            push!(messages,(false,"Last frame $nsteps not specified for scan $one_scan (last is $(maximum(ufn)))"))
            continue
        end
        if length(ufn) != nsteps
            push!(messages,(false,"Only $(length(ufn)) frames specified for scan $one_scan instead of $nsteps"))
        else
            push!(messages,(true,"All frames present and correct for $one_scan"))
        end
    end
    return messages
end

@noimgcheck "Detector surface axes used properly" check_surface_axes(incif) = begin
    messages = []

    surf_axes = unique(incif["_array_structure_list_axis.axis_id"])

    # these axes should not be listed in axis setting lists

    if haskey(incif,"_diffrn_scan_axis.axis_id")
        if !isdisjoint(surf_axes,incif["_diffrn_scan_axis.axis_id"])
            push!(messages,(false,"Detector surface axes $surf_axes should not be listed in _diffrn_scan_axis"))
        end
    end
    
    if haskey(incif,"_diffrn_scan_frame_axis.axis_id") &&
        !isdisjoint(surf_axes,incif["_diffrn_scan_frame_axis.axis_id"])
        push!(messages,(false,"Detector surface axes $surf_axes should not be listed in _diffrn_scan_frame_axis"))
    end
    return messages
end

@noimgcheck "Pixel size and origin described correctly" check_pixel_coords(incif) = begin
    messages = []

    origin = parse.(Float64,incif["_array_structure_list_axis.displacement"]) #in mm
    origin_set = incif["_array_structure_list_axis.axis_set_id"]
    if 0 in origin
        push!(messages,(false,"Origin of array coordinates should be in centre of pixel but at least one value of _array_structure_list_axis.displacement is 0"))
    end

    # Check that pixels are displaced by half their size
    
    if haskey(incif,"_array_element_size.size") # in metres!!
        elsize = parse.(Float64,incif["_array_element_size.size"])*10^3
        elind = parse.(Int64,incif["_array_element_size.index"])
        for oneind in elind
            asind = indexin([oneind],parse.(Int64,incif["_array_structure_list.index"]))[]
            setind = incif["_array_structure_list.axis_set_id"][asind]
            disp_ind = indexin([setind],origin_set)[]
            disp = origin[disp_ind]
            @debug disp elsize[oneind]
            if !isapprox(abs(disp), elsize[oneind]/2.0,atol=0.000001)
                push!(messages,(false,"Pixel axis displacement $disp mm is not at centre of pixel of corresponding dimension $(elsize[oneind]) mm"))
            end
        end
    end

    return messages
end

@noimgcheck "Check calculated beam centre" check_beam_centre(incif) = begin
    messages = []

    centre1,centre2,index1,index2 = get_beam_centre(incif)
    if index1 < 0 || index2 < 0
        push!(messages,(false,"Beam centre in pixels $index1, $index2 for zero axis settings has at least one negative value, which implies a position off the detector"))
    end

    dimensions = parse.(Int64,incif["_array_structure_list.dimension"])
    precs = parse.(Int64,incif["_array_structure_list.precedence"])

    # index1 is lowest precedence, index2 is highest

    indexify = indexin([2,1],precs)
    for (pix_cent, prec) in zip((index1,index2),indexify)
        if pix_cent > dimensions[prec]
            push!(messages,(false,"Beam centre coordinate in pixels $pix_cent for zero axis settings is greater than maximum dimension $(dimensions[prec]), which implies a position off the detector"))
        end
    end

    return messages
end

@noimgcheck "Check principal axis is aligned with X" check_principal_axis(incif) = begin
    messages = []

    ga,gt = get_gonio_axes(incif)
    axloop = get_loop(incif,"_axis.id")
    filt_axis = filter(row -> row["_axis.equipment"] == "goniometer" &&
                       row["_axis.type"] == "rotation" &&
                       row["_axis.depends_on"] == nothing, axloop, view=true)
    if size(filt_axis,1) != 1
        push!(messages,(false,"Incorrectly specified goniometer axes: more than one base goniometer axis (ie _axis.depends_on is '.')"))
    else
        row = filt_axis[1,:]
        if parse(Float64,row["_axis.vector[1]"]) != 1 || parse(Float64,row["_axis.vector[2]"]) != 0 || parse(Float64,row["_axis.vector[3]"]) != 0
            push!(messages,(false,"Principal axis $(row["_axis.id"]) is not the same as the X axis"))
        end
    end
    return messages
end

# Check that the image is described correctly
# We assume a single array structure for all data
const img_types = Dict(UInt8 =>"unsigned 8-bit integer",
                       UInt16=>"unsigned 16-bit integer",
                       UInt32=>"unsigned 32-bit integer",
                       Int8  =>"signed 8-bit integer",
                       Int16 =>"signed 16-bit integer",
                       Int32 =>"signed 32-bit integer",
                       Float32 =>"signed 32-bit real IEEE",
                       Float64 =>"signed 64-bit real IEEE",
                       ComplexF32 =>"signed 32-bit complex IEEE"
                       )

@imgcheck "Image type and dimensions" img_dims(incif,img,img_id) = begin
    messages = []
    fast_pixels = size(img)[1]
    slow_pixels = size(img)[2]
    fast_pos,slow_pos = indexin(["1","2"],incif["_array_structure_list.precedence"])
    dims = parse.(Int32,incif["_array_structure_list.dimension"])
    if dims[fast_pos] != fast_pixels
        push!(messages,(false,"Stated fast dimension $fast_pos does not match actual image dimension $fast_pixels"))
    end
    if dims[slow_pos] != slow_pixels
        push!(messages,(false,"Stated slow dimension $slow_pos does not match actual image dimension $slow_pixels"))
    end
    if haskey(incif,"_array_structure.encoding_type") && img_types[eltype(img)] != incif["_array_structure.encoding_type"][1]
        push!(messages,(false,"Stated encoding $(incif["_array_structure.encoding_type"][1]) does not match array element type $(eltype(img))"))
    end
    if haskey(incif,"_array_structure.byte_order") && haskey(incif,"_array_data.external_data_id")
        push!(messages,(true,"WARNING: byte order provided in file containing external data pointers"))
    end
    if haskey(incif,"_array_structure.compression") && incif["_array_structure.compression"][1] != "none" && haskey(incif,"_array_data.external_data_id")
        push!(messages,(true,"Externally-provided data by definition is uncompressed but compression is specified as $(incif["array_structure.compression"][1])"))
    end
    return messages
end

#######################################################################################
#
#             Check the full archive
#
#######################################################################################
"""
    Download all archives listed in `incif`, returning a dictionary {url -> local file}.
    If `get_full` is false, just start the download and abort as soon as a file is found.
    `subs` is a dictionary of local file equivalents to urls.
"""
download_archives(incif;get_full=false,pick=1,subs=Dict()) = begin
    urls = unique(incif["_array_data_external_data.uri"])
    if !isempty(subs)
        for u in urls
            if !(u in keys(subs))
                throw(error("No local file substitution found for $u; have $(keys(subs)) only"))
            end
        end
        if get_full return subs end
    end
    result_dict = Dict{String,String}()
    for u in urls
        if !isempty(subs)
            loc = URI("file://"*subs[u])
        else
            loc = URI(make_absolute_uri(incif,u))
        end
        pos = indexin([u],incif["_array_data_external_data.uri"])[]
        arch_type = nothing
        if haskey(incif,"_array_data_external_data.archive_format")
            arch_type = incif["_array_data_external_data.archive_format"][pos]
        end
        x = ""
        if !get_full && arch_type in ("TBZ","TGZ","TAR")
            try
                x = peek_image(loc,arch_type,incif,check_name=false,entry_no=pick)
            catch exn
                @debug "Peeking into $u gave error:" exn
            end
            @debug "Got $x for $u"
            result_dict[u] = x == nothing ? "" : x
        elseif get_full
            result_dict[u] = Downloads.download(u)
        else
            @warn "Partial downloading not available for $u: try full downloading with option -f"
        end
    end
    return result_dict     
end

test_archive_present(local_archives) = begin
    messages = []
    for (k,v) in local_archives
        if v == ""
            push!(messages, (false, "Unable to access $k"))
        end
    end
    return messages
end

@fullcheck "All members present" member_check(incif,all_archives) = begin
    messages = []
end

##############
#
#   Utility routines
#
##############

"""
    Find an image ID that is available in the archive. `archive_list` is
    a dictionary of URI -> item pairs where `item` is an archive member
    if have_full is false, or else it is an entire archive otherwise. 
"""
find_load_id(incif,archive_list,have_full) = begin

    known_paths = incif["_array_data_external_data.archive_path"]
    
    # First see if our frame is in the file

    if !have_full
        @debug archive_list
        for (k,v) in archive_list
            if v in known_paths
                pos = indexin([v],known_paths)[]
                # a level of indirection
                ext_data_id = incif["_array_data_external_data.id"][pos]
                vv = incif["_array_data.external_data_id"]
                pos = indexin([ext_data_id],vv)[]
                return incif["_array_data.binary_id"][pos]
            end
        end

        throw(error("Sample image not found, try full archive (option -f)"))
    end

    # Get a list of all archive members
    
end

verdict(msg_list) = begin
    ok = reduce((x,y)-> x[1] & y[1],msg_list;init=true)
    println(ok ? "PASS" : "FAIL")
    for (isok,message) in msg_list
        println("   "*message)
    end
    return ok
end

"""
Display an image in the terminal
"""
display_check_image(im;logscale=true,cut_ratio=1000) = begin
    #alg = Equalization(nbins=256,maxval = floor(maximum(im)/10))
    if maximum(im) > 1.0
        im = im/maximum(im)
    end
    clamp_low,clamp_high = find_best_cutoff(im,cut_ratio=cut_ratio)
    alg = LinearStretching(src_maxval = clamp_high)
    im_new = adjust_histogram(im,alg)
    println("Image for checking")
    #println("Max, min for adjusted image: $(maximum(im_new)), $(minimum(im_new))\n")
    imshow(Gray.(im_new))
    println("\n")
    return im_new
end

"""
Save an image
"""
save_check_image(im,imfn;logscale=true,cut_ratio=1000) = begin
    
    if maximum(im) > 1.0
        im = im/maximum(im)
    end
    clamp_low,clamp_high = find_best_cutoff(im,cut_ratio=cut_ratio)
    alg = LinearStretching(src_maxval = clamp_high)
    im_new = adjust_histogram(im,alg)
    save(imfn, im_new)
    return im_new
end

"""
Find the best value for displaying a diffraction image, calculated as the first
intensity bin that has 1000 times less points than the highest. Based on the
logic that the highest value will be a "typical" background point. Skip any
that correspond to negative values as these are likely to have a different
meaning.
"""
find_best_cutoff(im;cut_ratio=1000) = begin
    edges,bins = build_histogram(im)
    # Find largest number of points > 0
    maxpts = 0
    maxpos = 0
    for i in 1:length(bins)
        if first(edges) + (i-1)*step(edges) < 0 continue end
        if bins[i] > maxpts
            maxpts = bins[i]
            maxpos = i
        end
        if bins[i] < maxpts break end
    end
    cutoff = maxpts/cut_ratio
    maxbin = maxpos
    for i in maxpos+1:length(bins)
        if bins[i] < cutoff
            maxbin = i
            break
        end
    end
    maxval = first(edges) + maxbin*step(edges)
    minval = first(edges) + maxpos*step(edges)
    #println("Min,max $minval($maxpos),$maxval($maxbin)")
    return minval,maxval
end

#==

     End of routines for individual checks

==#

run_img_checks(incif;images=false,always=false,full=false,connected=false,pick=1,subs=Dict(),savepng=false) = begin
    ok = true
    println("Running checks (no image download)")
    println("="^40*"\n")
    for (desc,one_test) in test_list
        print("\nTesting: $desc: ")
        ok = ok & verdict(one_test(incif))
    end
    testimage = [[]]  # for consistency

    if !connected return (ok, testimage) end
    
    # Test archive access

    println("Testing presence of archive:")
    
    all_archives = download_archives(incif;get_full=full,pick=pick,subs=subs)

    print("\nTesting: All archives are accessible: ")
    
    ok = ok & verdict(test_archive_present(all_archives))
    
    # Test with an image
    
    if length(test_list_with_img) > 0 && ((ok && images) || always)
        testimage = [[]]
        println("\nRunning checks with downloaded images")
        println("="^40*"\n")

        # Choose image to load

        load_id = find_load_id(incif,all_archives,full)
        
        try
            if full
                testimage = imgload(incif,load_id;local_version=all_archives)
            else
                testimage = imgload(incif,load_id;local_version=subs)
            end
        catch e
            @debug e
            verdict([(false,"Unable to access image $load_id: $e")])
            rethrow()
        end
        
        if savepng
            imgfn=string(incif.original_file,".png") 
            save_check_image(testimage, imgfn, logscale=false)
        else
            display_check_image(testimage,logscale=false)
        end
        
        for (desc,one_test) in test_list_with_img
            print("\nTesting image $load_id: $desc: ")
            ok = ok & verdict(one_test(incif,testimage,load_id))
        end
    end
    
    # Tests requiring fully-downloaded archives

    if full
        for (desc,one_test) in test_full_list
            print("\nTesting full archive: $desc:")
            ok = ok & verdict(one_test(incif,all_archives))
        end
    end
    return (ok,testimage)
end

parse_cmdline(d) = begin
    s = ArgParseSettings(d)
    @add_arg_table! s begin
        "-i", "--check-images"
        help = "Also perform checks on the images"
        nargs = 0
        "-j", "--always-check-images"
        help = "Check images even if non-image checks fail"
        nargs = 0
        "-o", "--output-png"
        help = "Save downloaded image as <filename>.png"
        nargs = 0
        "-d", "--dictionary"
        help = "Check conformance to <dictionary>"
        nargs = 1
        default = [""]
        arg_type = String
        "--dictionary-only"
        help = "Check dictionary conformance only (-d option must be provided)"
        nargs = 0
        "-f", "--full-download"
        help = "Fully download archive for image and archive checks (required for ZIP)"
        nargs = 0
        "-n", "--no-internet"
        help = "Do not check raw image archive existence or contents"
        nargs = 0
        "-p", "--pick"
        nargs = 1
        default = [1]
        arg_type = Int64
        help = "Use this entry number in the archive for checking"
        "-s", "--sub"
        nargs = 2
        metavar = ["original","local"]
        action = "append_arg"
        help = "Use <local> file in place of URI <original> (for testing). Use -f to access whole file. All URIs in <filename> must be provided, option may be used multiple times"
        "filename"
        help = "Name of imgCIF data file to check"
        required = true
        "blockname"
        help = "Block name to check. If missing, the first block is checked"
        required = false
    end
    parse_args(s)
end

if abspath(PROGRAM_FILE) == @__FILE__
    parsed_args = parse_cmdline("Check contents of imgCIF files")
    #println("$parsed_args")
    incif = Cif(Path(parsed_args["filename"]))
    if isnothing(parsed_args["blockname"])
        blockname = first(incif).first
    else
        blockname = parsed_args["blockname"]
    end
    subs = Dict(parsed_args["sub"])
    println("\n ImgCIF checker version 0.3\n")
    println("Checking block $blockname in $(incif.original_file)\n")
    if parsed_args["dictionary"] != [""]
    end
    result,img = run_img_checks(incif[blockname],
                                images=parsed_args["check-images"],
                                always=parsed_args["always-check-images"],
                                full = parsed_args["full-download"],
                                connected = !parsed_args["no-internet"],
                                pick = parsed_args["pick"][],
                                subs = subs,
                                savepng = parsed_args["output-png"]
                                )
    println("\n====End of Checks====")
    if result exit(0) else exit(1) end
end

    
