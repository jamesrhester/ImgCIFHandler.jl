#
#  Auto-install all required items
#
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ImgCIFHandler#main
using ImageInTerminal, Colors,ImageContrastAdjustment
using ImageTransformations
using ArgParse
using CrystalInfoFramework,FilePaths,URIs, Tar, Images
using Luxor

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
        if any(x->ismissing.(x),incif[ic])
            push!(messages,(false,"At least one value for required item $ic is missing"))
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

    # Check that axes have a length

    vs = getindex.(Ref(incif),["_axis.id","_axis.vector[1]","_axis.vector[2]","_axis.vector[3]"])
    for (nm,v1,v2,v3) in zip(vs...)
        v = parse.(Float64,[v1,v2,v3])
        if v[1]^2 + v[2]^2 + v[3]^2 <= 0
            push!(messages,(false,"Axis vector for $nm has length <= 0"))
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

    actual = filter(scan_loop) do r
        !ismissing(getproperty(r,cn*"displacement_range")) || !ismissing(getproperty(r,cn*"angle_range"))
    end
    
    if size(actual,2) == 0
        return [(false,"No scan range has been specified")]
    end

    for or in eachrow(actual)
        dr = getproperty(or,cn*"displacement_range")
        if isnothing(dr) || ismissing(dr)  #really should not allow missing
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

#             Check the full archive

"""
    download_uris(incif;subs=Dict)

    Download uris from `incif` to local files, where `subs` contains
    local file equivalents to urls {url -> local file}
"""
download_uris(incif;subs=Dict())= begin

    result_dict = Dict{String,String}()

    # Get all referenced URIs
    
    urls = unique(incif["_array_data_external_data.uri"])

    # If we have local substitutes, make sure they cover all uris found
    
    for u in urls
        if u in keys(subs)
            result_dict[u] = URI("file://"*subs[u])
            continue
        end
        loc = URI(make_absolute_uri(incif,u))
        result_dict[u] = Downloads.download(loc)
    end

    return result_dict
end

"""
    get_archive_member_name(incif;pick=1,subs=Dict())

    For each distinct URI in `incif`, return the name of a member file that is listed
    in `incif`. `subs` is a dictionary of local file equivalents to urls. `pick`
    selects the nth member of the archive instead of the first.
"""
get_archive_member_name(incif;pick=1,subs=Dict()) = begin

    @debug "Looking for member in" subs
    
    # Get all referenced URIs
    
    urls = unique(incif["_array_data_external_data.uri"])

    result_dict = Dict{String,String}()

    # Cycle through URIs, looking either for everything or one thing
    
    for u in urls

        # Construct actual URI
        
        if haskey(subs,u)
            loc = URI("file://"*subs[u])
        else
            if !isempty(subs)
                @warn "Warning: no substitute for $u"
            end
            loc = URI(make_absolute_uri(incif,u))
        end

        # Find archive type
        
        pos = indexin([u],incif["_array_data_external_data.uri"])[]
        arch_type = nothing
        if haskey(incif,"_array_data_external_data.archive_format")
            arch_type = incif["_array_data_external_data.archive_format"][pos]
        end

        # Extract a single file if we have an archive
        
        x = ""
        if arch_type in ("TBZ","TGZ","TAR")
            try
                x = peek_image(loc,arch_type,incif,check_name=false,entry_no=pick)
            catch exn
                @debug "Peeking into $u gave error:" exn
            end
            @debug "Got $x for $u"
            result_dict[u] = x == nothing ? "" : x
        elseif arch_type === nothing  #No unpacking possible eg HDF5/single frame
            #TODO: actually check the file exists using http request
            result_dict[u] = "$loc"
        else
            @warn "Partial downloading not available for $u: try full downloading with option -f"
        end
    end

    # result_dict links uri to a filename, which is contained in that uri

    @debug "Final result" result_dict
    
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


#   Utility routines

"""
    Find an image ID that is obtainable. `archive_list` is
    a dictionary of URI -> item pairs where `item` is an archive member
    or the URI itself if the URI is a complete frame.
"""
find_load_id(incif,archive_list) = begin

    external_info = get_loop(incif,"_array_data_external_data.id")
    known_uris = unique(incif["_array_data_external_data.uri"])
    if haskey(incif,"_array_data_external_data.archive_path")
        known_paths = incif["_array_data_external_data.archive_path"]
    else
        known_paths = known_uris
    end
    
    # First see if our frame is in the file

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
    create_check_image(incif,im;logscale=true,cut_ratio=1000,gravity=true)

Create a contrast enhanced image from the matrix `im`. If `gravity` is true, 
information in `incif` is used, if available, to calculate the rotation
of the image so that down on the image 
becomes down in real space. Returns the unrotated image and the rotation.
"""
create_check_image(incif,im;logscale=true,cut_ratio=1000,gravity=true) = begin
    
    # First try to improve contrast
    
    #alg = Equalization(nbins=256,maxval = floor(maximum(im)/10))
    if maximum(im) > 1.0
        im = im/maximum(im)
    end
    clamp_low,clamp_high = find_best_cutoff(im,cut_ratio=cut_ratio)
    alg = LinearStretching(src_maxval = clamp_high)
    im_new = adjust_histogram(im,alg)
    @debug "Max, min for adjusted image:" maximum(im_new) minimum(im_new)

    # Adjust geometry if we know gravity

    if !gravity || !("gravity" in incif["_axis.equipment"]) return im_new,0 end

    gravity = indexin(["gravity"],incif["_axis.equipment"])[]
    grav_vec = parse.(Float64,
                      [
                          incif["_axis.vector[1]"][gravity],
                          incif["_axis.vector[2]"][gravity],
                          incif["_axis.vector[3]"][gravity],
                      ]
                      )
    norm_grav = sqrt(grav_vec[1]^2+grav_vec[2]^2)
    grav_vec = grav_vec/norm_grav
    
    corner_loc = get_pixel_coordinates(incif,0,0)
    fast_dir = get_pixel_coordinates(incif,0,1) - corner_loc

    norm_fast = sqrt(fast_dir[1]^2+fast_dir[2]^2)
    fast_dir = fast_dir/norm_fast

    # Make "down" (PNG fast direction) in direction of gravity.
    
    # Assume square pixels. We want rotation matrix
    # R[theta][fast_x,fast_y] = [grav_x,grav_y]
    # where one of grav_x,grav_y = 0 and we require one of the
    # components of fast to be basically zero
    #
    # So fast_x * cos theta  -fast_y*sin theta  = grav_x
    # sin theta * fast_x +  fast_y * cos_theta = grav_y
    #

    # Simplest is just to try the 4 possible rotations

    rot = nothing
    for d in 0:3
        r=d*90
        try_x = fast_dir[1]*cosd(r) - fast_dir[2]*sind(r)
        try_y = fast_dir[1]*sind(r) + fast_dir[2]*cosd(r)
        if isapprox(try_x,grav_vec[1],atol=0.1) && isapprox(try_y,grav_vec[2],atol=0.1)
            rot = d
            @debug "Image should be rotated by " rot*90
            break
        end
    end

    if rot == nothing
        println("Warning: unable to rotate image to match gravity")
        rot = 0
    end

    #im_new = rotl90(im_new,rot)
    
    return im_new,rot*90
end

"""
    annotate_check_image(image,rotation,beam_centre,names,filename; border=30)

Add axes and beam centre to `image`, saving the result in
`filename`. `names` contains the names of the fast and slow
axes, in that order. `border` is the size of the border around
the image.
"""
annotate_check_image(im, rot, incif;border=30,scan_id=nothing,frame_no=nothing) = begin

    width_orig = size(im,2)
    height_orig = size(im,1)   #in display, height is fast direction
    
    # Rotate the image matrix

    im_new = rotl90(im,div(rot,90))
    
    # Get the beam centre

    if isnothing(scan_id)
        _,_,slow_c,fast_c = get_beam_centre(incif)
    else
        _,_,slow_c,fast_c = get_beam_centre(incif,scan_id,frame_no)
    end
    fast_n,slow_n = get_surface_axes(incif)

    # Transform beam centre coordinates

    slow_c = slow_c - width_orig/2
    fast_c = fast_c - height_orig/2
    new_slow_c = cosd(rot)* slow_c + sind(rot)* fast_c
    new_fast_c = -sind(rot)* slow_c + cosd(rot)* fast_c
    new_slow_c = new_slow_c + width_orig/2
    new_fast_c = new_fast_c + height_orig/2
                             
    # TODO: scale to useful size
    
    width_new = size(im_new,2)   #slow direction
    height_new = size(im_new,1)  #fast direction

    @debug "Old, new height, width" height_orig width_orig height_new width_new

    # Start the drawing
    
    Luxor.Drawing(width_new+width_orig+border*4,maximum((height_new,height_orig))+border*2+20, :rec)
    background("white")

    # Place rotated and original images

    placeimage(Gray.(im_new),Point(border,border))
    placeimage(Gray.(im),Point(border*3 + width_new,border))
    setcolor("black")
    label("Laboratory frame",0,Point(border+width_new/2,border*2+height_new+10))
    label("Original",0,Point(border*3+width_new+width_orig/2,border*2+height_orig+10))
    
    # draw a square for the beam centre
    
    setcolor("red")
    box(Point(new_slow_c+border,new_fast_c+border),4,4,:fill)
    draw_axes(height_new,width_new,rot,(fast_n,slow_n))
    snapshot(fname="$(incif.original_file)"*".png")
end

draw_axes(height,width,angle,names;border=30) = begin
    fast,slow = names
    @debug "Fast, slow axis names" fast slow
    @debug "Rotate a/c wise by " angle
    
    # rotate and origin to corner
    origin(width/2+border,height/2+border)
    rotate(deg2rad(-1*angle))
    translate(Point(-(width+border)/2,-(height+border)/2))
    # draw some lines

    setcolor("black")
    
    # The X coordinate is across (i.e. the slow png direction)

    arrow(O,Point(width/2,0))

    # And Y is down (i.e. the fast png direction)

    arrow(O,Point(0,(height/2.0)))

    # annotate the lines

    translate(Point(width/4.0,0))

    # switch the label on the bottom to get the text upright
    if angle == 180 # slow axis on bottom
	rotate(pi)
	label(slow,:N)
	rotate(pi)
    else
	label(slow,:S)
    end
    translate(Point(-width/4.0,height/4.0))
    if angle == 90  # fast on bottom
        @debug "Rotating $fast"
	rotate(pi/2)
	label(fast,:N)
	rotate(-pi/2)
    else
	rotate(-pi/2)
	label(fast,:S)
    end
end

show_check_image(im::AbstractArray,rot) = begin
    println("Image for checking")
    imshow(Gray.(rotl90(im,div(rot,90))))
    println("\n")
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

"""
    fix_loops!(cif)

Make all items that are unlooped but we need to be looped, into
loops
"""
fix_loops!(incif) = begin

    loop_groups = (("_diffrn_scan_axis",["axis_id","displacement_range","angle_range",
                                         "angle_start","angle_increment",
                                         "displacement_increment","displacement_start"]),
                   ("_diffrn_scan",["id","frames"])
                   )

    for (cn,lnames) in loop_groups
        test_loop = get_loop(incif,cn*"."*lnames[1])
        if size(test_loop,1) == 0 && haskey(incif,cn*"."*lnames[1])
            for n in lnames
                if !haskey(incif,cn*"."*n) incif[cn*"."*n] = [missing] end
            end
            create_loop!(incif,map(x->cn*"."*x,lnames))
        end
    end
end

"""
    run_img_checks(incif;images=false,always=false,full=false,connected=false,pick=1,subs=Dict(),savepng=false,accum=1)

Test `incif` for conformance to imgCIF consistency requirements. Meaning of arguments:
`images`: run checks on downloaded image
`always`: always run checks on images, even if previous checks fail
`full`: download all referenced archives in full
`connected`: perform checks that require an internet connection
`pick`: archive member to download for image checks
`subs`: dictionary of uri -> local file correspondences to avoid downloading
`savepng`: output annotated check image to a file
`accum`: accumulate this many frames to make the check image
"""
run_img_checks(incif;images=false,always=false,full=false,connected=false,pick=1,subs=Dict(),savepng=false,accum=1) = begin
    ok = true
    println("Running checks (no image download)")
    println("="^40*"\n")
    for (desc,one_test) in test_list
        print("\nTesting: $desc: ")
        res = []
        try
            res = one_test(incif)
        catch e
            @debug "Error during test" e
            res = [(false,"Unable to carry out test, assume missing or bad value")]
        end
        ok = ok & verdict(res)
    end
    testimage = [[]]  # for consistency

    if !connected return (ok, testimage) end
    
    # Test archive access

    println("Testing presence of archive:")

    if full
        subs = download_uris(incif,subs)
    end
    
    all_archives = get_archive_member_name(incif;pick=pick,subs=subs)

    print("\nTesting: All archives are accessible: ")
    
    ok = ok & verdict(test_archive_present(all_archives))
    
    # Test with an image
    
    if length(test_list_with_img) > 0 && ((ok && images) || always)
        testimage = [[]]
        println("\nRunning checks with downloaded images")
        println("="^40*"\n")

        # Choose image(s) to load

        load_id = find_load_id(incif,all_archives)

        if accum > 1
            load_ids = get_id_sequence(incif,load_id,accum)
        else
            load_ids = [load_id]
        end
        
        try
            testimage = imgload(incif,load_ids;local_version=subs)
        catch e
            @debug e
            verdict([(false,"Unable to access image $load_id: $e")])
            rethrow()
        end

        # Output an image
        
        new_im,rot = create_check_image(incif,testimage,logscale=false)
        imgfn = nothing
        if savepng
            scan_id,frame_no = ImgCIFHandler.scan_frame_from_bin_id(load_id,incif)
            annotate_check_image(new_im,rot,incif,scan_id=scan_id,frame_no=frame_no)
        else
            show_check_image(new_im,rot)
        end
        # Run the image checks
        
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
        "-a", "--accumulate"
        help = "Accumulate multiple images for check"
        nargs = 1
        default = [1]
        metavar = ["number"]
        arg_type = Int64
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
    incif = Cif(FilePaths.Path(parsed_args["filename"]),native=true)
    if isnothing(parsed_args["blockname"])
        blockname = first(incif).first
    else
        blockname = parsed_args["blockname"]
    end
    subs = Dict(parsed_args["sub"])
    println("\n ImgCIF checker version 2022-06-29\n")
    println("Checking block $blockname in $(incif.original_file)\n")
    if parsed_args["dictionary"] != [""]
    end

    # Fix loops

    fix_loops!(incif[blockname])
    
    result,img = run_img_checks(incif[blockname],
                                images=parsed_args["check-images"],
                                always=parsed_args["always-check-images"],
                                full = parsed_args["full-download"],
                                connected = !parsed_args["no-internet"],
                                pick = parsed_args["pick"][],
                                accum = parsed_args["accumulate"][],
                                subs = subs,
                                savepng = parsed_args["output-png"]
                                )
    println("\n====End of Checks====")
    if result exit(0) else exit(1) end
end

    
