# Utility routines for performing non-image imgCIF operations
"""
    get_detector_axes(imgcif::CifContainer;check=false)

Return a list of detector axes, together with type (rotation/translation).
If `check` is true, raise an error if any axes have not been assigned to a detector.
Do not include axes describing the pixel array.
"""
get_detector_axes(imgcif::CifContainer;check=false) = begin
    axis_loop = get_loop(imgcif,"_axis.id")
    pix_axes = String.(imgcif["_array_structure_list_axis.axis_id"])

    filt_axis = filter(row -> row["_axis.equipment"] == "detector" && !(row["_axis.id"] in pix_axes),axis_loop, view=true)
    det_axes = String.(filt_axis[!,"_axis.id"])

    if check
        other_det_axes = get_loop(imgcif,"_diffrn_detector_axis.axis_id")
        if !issetequal(String.(other_det_axes[!,"_diffrn_detector_axis.axis_id"]),det_axes)
            throw(error("Detector axes in _diffrn_detector_axis do not match axis category"))
        end
    end
    return det_axes,String.(filt_axis[!,"_axis.type"])
end

"""
    get_detector_axis_settings(imgcif,scanid,frameno)

Return the settings of all detector axes for the specified frame together
with axis names and types in order names,types,settings. Values are
computed from _diffrn_scan_axis information if _diffrn_scan_frame_axis
is missing. Set `scanid` to `nothing` or `missing` if no scans are
explicitly defined in the file.
"""
get_detector_axis_settings(imgcif::CifContainer,scanid,frameno) = begin
    axis_names, axis_types = get_detector_axes(imgcif)
    if haskey(imgcif,"_diffrn_scan_frame_axis.axis_id")
        return get_explicit_frame(imgcif,scanid,frameno)
    end

    # Calculate settings

    cat = "_diffrn_scan_axis"
    if !haskey(imgcif,"$cat.axis_id")

        # Assume axis settings are default values
        return zip(axis_names, axis_types, fill("0",length(axis_names)))
    end

    scan_frames = get_loop(imgcif,"$cat.axis_id")
    our_scan = filter(row->row["$cat.scan_id"] == scanid, scan_frames,view=true)
    positions = map(axis_names,axis_types) do one_axis,one_type
        our_axis = filter(row->row["$cat.axis_id"] == one_axis,our_scan,view=true)
        if size(our_axis)[1] == 0
            0
        else
            our_axis = our_axis[1,:]
            isrot = one_type == "rotation"
            if isrot
                start = our_axis["$cat.angle_start"]
                increment = our_axis["$cat.angle_increment"]
            else
                start = our_axis["$cat.displacement_start"]
                increment = our_axis["$cat.displacement_increment"]
            end
            start = !ismissing(start) && start != nothing ? parse(Float64,start) : 0
            increment = !ismissing(increment) && increment != nothing ? parse(Float64,increment) : 0
            start + (frameno-1)*increment
        end
    end
    return axis_names,axis_types,positions
end

# Version given a filename
get_detector_axis_settings(filename::AbstractString,args...) = begin
    get_detector_axis_settings(first(Cif(Path(filename))).second,args...)
end

# Version if no scanid or frameid
get_detector_axis_settings(imgcif::CifContainer) = begin
    axis_names,axis_types = get_detector_axes(imgcif)
    return axis_names,axis_types,fill(0.0,length(axis_names))
end

"""
    get_gonio_axes(imgcif::CifContainer;check=false)

Return a list of goniometer axes, together with type (rotation/translation).
"""
get_gonio_axes(imgcif::CifContainer) = begin
    axis_loop = get_loop(imgcif,"_axis.id")

    filt_axis = filter(row -> row["_axis.equipment"] == "goniometer",axis_loop, view=true)
    gonio_axes = String.(filt_axis[!,"_axis.id"])

    return gonio_axes,String.(filt_axis[!,"_axis.type"])
end

"""
    get_surface_axes(imgcif::CifContainer)

Return the axes describing the detector surface together with the origin
in X,Y coordinates
"""
get_surface_axes(imgcif::CifContainer) = begin
    
end
