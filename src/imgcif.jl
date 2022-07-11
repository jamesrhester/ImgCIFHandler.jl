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
    get_detector_axis_settings(first(Cif(Path(filename),native=true)).second,args...)
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

Return the axis names describing the detector surface and return in order
fast,slow
"""
get_surface_axes(incif::CifContainer) = begin
    names = incif["_array_structure_list_axis.axis_id"]
    # assume no fancy multi-axis detector pixel directions
    axis_sets = incif["_array_structure_list_axis.axis_set_id"]
    precedence = parse.(Int64,incif["_array_structure_list.precedence"])
    prec_axis_set = incif["_array_structure_list.axis_set_id"]
    speeds = indexin([1,2],precedence)
    axis_id = prec_axis_set[speeds]
    axis_name_ind = indexin(axis_id,axis_sets)
    fast,slow = String.(names[axis_name_ind])
end

"""
    scan_frame_from_img_name(u,name,cif_block)

Return the scan_id, frame_id for the given `name` at URL `u` according to contents
of `cif_block`. If `name` is `nothing`, `u` must point to a single frame. If `cif_block`
contains no scan or frame information, (nothing,1) is returned.
"""
scan_frame_from_img_name(u,name,cif_block) = begin
    ext_loop = get_loop(cif_block,"_array_data_external_data.id")
    fel = filter(row->row["_array_data_external_data.uri"] == u,ext_loop,view=true)
    if size(fel,1) == 0 throw(error("$u not found in cif block")) end
    if name != nothing
        fel = filter(row->row["_array_data_external_data.archive_path"]==name,fel,view=true)
    end
    if size(fel,1) != 1
        throw(error("$u $name does not correspond to a single frame (got $(size(fel,1)) answers)"))
    end

    # now find the frame number external_id -> binary_id -> frame_id -> frame_no

    ext_id = fel[!,"_array_data_external_data.id"][]
    array_loop = filter(row->row["_array_data.external_data_id"]==ext_id,
                        get_loop(cif_block,"_array_data.binary_id"),view=true)
    bin_id = array_loop[!,"_array_data.binary_id"][]
    return scan_frame_from_bin_id(bin_id,cif_block)
end

"""
    scan_frame_from_bin_id(binary_id,cif_block)

Return the scan_id, frame_id for the given `binary_id` according to contents
of `cif_block`. If `cif_block` contains no scan or frame information, 
(nothing,1) is returned.
"""
scan_frame_from_bin_id(bin_id,cif_block) = begin

    # assume that frame_id is globally unique, not per scan
    
    frame_loop = filter(row->row["_diffrn_data_frame.binary_id"]==bin_id,
                      get_loop(cif_block,"_diffrn_data_frame.id"),view=true)
    frame_id = frame_loop[!,"_diffrn_data_frame.id"][]

    # get scan id
    
    scan_loop = filter(row->row["_diffrn_scan_frame.frame_id"]==frame_id,
                       get_loop(cif_block,"_diffrn_scan_frame.frame_id"), view=true)
    if !haskey(cif_block,"_diffrn_scan.id")
        scan_id = nothing
    elseif length(cif_block["_diffrn_scan.id"])==1
        scan_id = cif_block["_diffrn_scan.id"][]
    else
        scan_id = scan_loop[!,"_diffrn_scan_frame.scan_id"][]
    end
    frame_no = scan_loop[!,"_diffrn_scan_frame.frame_number"][]
    return scan_id,parse(Int64,frame_no)
end

"""
    external_specs_from_bin_ids(binary_ids::Vector,incif)

Get the download specifications for each of the provided binary ids, returning
a DataFrame. 
"""
external_specs_from_bin_ids(bin_ids::Vector,c) = begin
    ext_cat = "_array_data_external_data"   #for convenience
    cat = "_array_data"
    
    img_loop = get_loop(c,"$cat.binary_id")
    if !("$cat.external_data_id" in names(img_loop))
        @error "$(c.original_file) does not contain external data pointers" img_loop
    end
    info = filter(row -> row["$cat.binary_id"] in bin_ids, img_loop,view=true)
    if size(info)[1] == 0
        throw(error("No array data with ids $bin_ids found"))
    end

    # One level of indirection

    ext_ids = info[!,"$cat.external_data_id"]
    ext_loop = get_loop(c,"$ext_cat.id")
    info = filter(row -> row["$ext_cat.id"] in ext_ids, ext_loop)

    rename!(x-> replace(x,"$ext_cat." => ""),info)

    transform!(info, "uri" => (x->make_absolute_uri.(Ref(c),x)) => :full_uri)
    
    @debug "Images for $bin_ids: $info"

    #==
    return_dict["ext_loc"] = "$ext_cat.path" in all_cols ? info["$ext_cat.path"] : nothing
    return_dict["ext_format"] = "$ext_cat.format" in all_cols ? Val(Symbol(info["$ext_cat.format"])) : nothing
    return_dict["ext_comp"] = "$ext_cat.archive_format" in all_cols ? info["$ext_cat.archive_format"] : nothing
    return_dict["ext_ap"] = "$ext_cat.archive_path" in all_cols ? info["$ext_cat.archive_path"] : nothing
    return_dict["ext_frame"] = "$ext_cat.frame" in all_cols ? info["$ext_cat.frame"] : nothing
    return return_dict
    ==#

    return info
end

external_specs_from_bin_ids(bin_id,c) = begin
    external_specs_from_bin_ids([bin_id],c)
end

"""
    get_id_sequence(incif,binary_id,total)

Return the binary ids of a sequence of frames starting from the supplied `binary_id`
and containing no more than `total` frames.
"""
get_id_sequence(incif,binary_id,total) = begin

    # Get scan name and frame number of start frame

    sc,fr = scan_frame_from_bin_id(binary_id,incif)

    @debug "Start frame is" sc fr
    
    # Get total number of frames for this scan
    
    if haskey(incif,"_diffrn_scan.frames")

        scan_loop = get_loop(incif,"_diffrn_scan.frames")

        if isnothing(sc)
            all_frames = parse(Int64,incif["_diffrn_scan.frames"][])
        else
            our_scan = filter(row->row["_diffrn_scan.id"]==sc,scan_loop,view=true)
            all_frames = parse(Int64,our_scan[!,"_diffrn_scan.frames"][])
        end
    else
        all_frames = length(incif[!,"_diffrn_scan_frame.frame_id"])
    end

    # Calculate the final frame number
    
    final_fr = minimum((all_frames,fr+total-1))

    # Collect the binary ids of each frame

    f_loop = get_loop(incif,"_diffrn_scan_frame.frame_id")
    if haskey(incif,"_diffrn_scan_frame.scan_id") && sc != nothing
        f_loop = filter(row->row["_diffrn_scan_frame.scan_id"]==sc,f_loop,view=true)
    end

    ddf_loop = get_loop(incif,"_diffrn_data_frame.id")

    bin_ids = map(fr:final_fr) do f
        f_id = filter(row->row["_diffrn_scan_frame.frame_number"] == "$f",f_loop)
        f_id = f_id[!,"_diffrn_scan_frame.frame_id"][]
        bin_id = filter(row->row["_diffrn_data_frame.id"]==f_id,ddf_loop,view=true)
        bin_id[!,"_diffrn_data_frame.binary_id"][]
    end

    @assert bin_ids[1] == binary_id   #sanity check
    @debug "$total frames from $binary_id" bin_ids
    return bin_ids
end

