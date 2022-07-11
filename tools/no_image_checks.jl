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
