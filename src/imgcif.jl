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
get_detector_axis_settings(imgcif::CifContainer, args...) = begin
    axis_names, _ = get_detector_axes(imgcif)
    get_axis_settings(imgcif, axis_names, args...)
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

get_gonio_axis_settings(filename::AbstractString,args...) = begin
    get_gonio_axis_settings(first(Cif(Path(filename),native=true)).second,args...)
end

get_gonio_axis_settings(imgcif::CifContainer, args...) = begin
    axis_names, _ = get_gonio_axes(imgcif)
    get_axis_settings(imgcif, axis_names, args...)
end

get_measurement_axis(filename,scanid) = begin
    get_measurement_axis(first(Cif(Path(filename),native=true)).second,scanid)
end

get_measurement_axis(imgcif::CifContainer,scanid) = begin

    c = "_diffrn_scan_axis"
    dsa_loop = get_loop(imgcif,"$c.angle_range")
    dsaf = filter(r -> r["$c.scan_id"]==scanid,dsa_loop)
    filter!(dsaf) do r
        r["$c.angle_range"]!=nothing && parse(Float64,r["$c.angle_range"])!= 0
    end
    if size(dsaf,1) > 1
        throw(error("Cannot handle more than one scan axis"))
    end
    if size(dsaf,1) == 0
        throw(error("No scan carried out"))
    end
    return dsaf[!,"$c.axis_id"][]
end

get_axis_settings(imgcif::CifContainer,axis_names, args...) = begin

    axis_loop = get_loop(imgcif,"_axis.id")
    axis_types = map(axis_names) do an
        id = indexin([an],imgcif["_axis.id"])[]
        imgcif["_axis.type"][id]
    end

    if haskey(imgcif,"_diffrn_scan_frame_axis.axis_id") && length(args) > 0
        return get_explicit_frame(imgcif, args...)
    end

    # Calculate settings

    cat = "_diffrn_scan_axis"
    if !haskey(imgcif,"$cat.axis_id") || length(args) == 0

        # Assume axis settings are default values
        return axis_names, axis_types, fill(0.0,length(axis_names))
    end

    scanid, frameno = args
    
    scan_frames = get_loop(imgcif,"$cat.axis_id")

    # If we have a choice of scans filter the one we want
    
    if haskey(imgcif,"_diffrn_scan.id") && length(imgcif["_diffrn_scan.id"]) > 1 
        our_scan = filter(row->row["$cat.scan_id"] == scanid, scan_frames,view=true)
    else
        our_scan = scan_frames
    end
    
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
    return axis_names, axis_types, positions
end

"""
    get_gonio_axes(imgcif::CifContainer;check=false)

Return a list of goniometer axes, together with type (rotation/translation). Axes are
in order of dependency chain, with base axis first. Typically omega-chi-phi.
"""
get_gonio_axes(imgcif::CifContainer) = begin
    axis_loop = get_loop(imgcif,"_axis.id")

    filt_axis = filter(row -> row["_axis.equipment"] == "goniometer",axis_loop, view=true)
    gonio_axes = String.(filt_axis[!,"_axis.id"])

    return sort_axes(imgcif,gonio_axes),String.(filt_axis[!,"_axis.type"])
end

sort_axes(imgcif::CifContainer,axis_list) = begin
    
    axis_dep_order = []
    axis_loop = get_loop(imgcif,"_axis.id")
    axis_loop = filter(r->r["_axis.id"] in axis_list,axis_loop)
    test_list = vcat(axis_list,nothing)
    ax_name = nothing

    while true
        ax_dep = filter(r->r["_axis.depends_on"]==ax_name && ax_name in test_list,axis_loop)
        if size(ax_dep,1) == 0 break end
        if size(ax_dep,1) > 1
            throw(error("More than one axis depends on $ax_name: $ax_dep"))
        end
        new_axis = ax_dep[!,"_axis.id"][]
        if new_axis in axis_dep_order
            throw(error("Dependency loop in axis.depends_on!"))
        end
        push!(axis_dep_order,new_axis)
        ax_name = new_axis
    end

    if length(setdiff(axis_dep_order,axis_list))!= 0
        @warn "Not all provided axes are in the dependency chain" axis_dep_order axis_list
    end

    # @debug "Axes in dependency order" axis_dep_order
    
    return axis_dep_order
end

"""
    get_dependency_chain(imgcif::CifContainer,axis_id)

Return the list of dependent axes for `axis_id`, with `axis_id` the first
entry in the returned list.
"""
get_dependency_chain(imgcif::CifContainer,axis_id) = begin
    axis_loop = get_loop(imgcif,"_axis.id")
    dep_chain = [axis_id]
    current_axis = axis_id
    while true
        new_ax = filter(r->r["_axis.id"]==current_axis,axis_loop)
        @assert size(new_ax,1) == 1
        dep_axis = new_ax[!,"_axis.depends_on"][]
        if dep_axis == nothing break end
        if dep_axis in dep_chain
            @error "Dependency loop for axes" dep_axis dep_chain ex=ErrorException
            throw(error("Dependency loop: $dep_axis is already in $dep_chain"))
        end
        push!(dep_chain,dep_axis)
        current_axis = dep_axis
    end
    return dep_chain
end

"""
    get_surface_axes(imgcif::CifContainer)

Return the axis names describing the detector surface and return in order
fast,slow.
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
    get_axis_vector(incif::CifContainer,axis_id)

Return the direction of `axis_id`
"""
get_axis_vector(incif::CifContainer,axis_id) = begin
    axis_ind = indexin([axis_id],incif["_axis.id"])[]
    s = "_axis.vector"
    map([1,2,3]) do i
        parse(Float64,incif["$s[$i]"][axis_ind])
    end
end

get_axis_offset(incif::CifContainer,axis_id) = begin
    axis_ind = indexin([axis_id],incif["_axis.id"])[]
    s = "_axis.offset"
    map([1,2,3]) do i
        parse(Float64,incif["$s[$i]"][axis_ind])
    end
end

"""
    get_axis_vector(incif::CifContainer,axis_id,scan_id,frame_no)

Return the direction of `axis_id` for `frame_no` of `scan_id`. The
axis vector of the nth axis up the stack is the result of applying
all underlying rotations to its starting vector.
"""
get_axis_vector(incif::CifContainer,axis_id,scan_id,frame_no) = begin

    get_axis_poise(incif,axis_id,scan_id,frame_no)[1]

end

"""
    get_axis_vector(incif::CifContainer,axis_id,axis_vals)

Return the direction of `axis_id` when the axes take values
given by `axis_vals`, which is a dictionary of name=>setting.
The axis vector of the nth axis up the stack is the result of applying
all underlying rotations to its starting vector.
"""
get_axis_vector(incif::CifContainer,axis_id,axis_vals) = begin

    start_vector = get_axis_vector(incif, axis_id)  #starting vector
    dep_chain = get_dependency_chain(incif, axis_id)

    if length(dep_chain) == 1 return start_vector end
    
    ax_poise = get_axis_vector.(Ref(incif), dep_chain[2:end], Ref(axis_vals))

    # Apply rotations in order

    current = start_vector
    for (axname,axpz) in zip(reverse(dep_chain[2:end]),reverse(ax_poise))
        axpos = axis_vals[axname]
        @debug "Now rotating along $axpz ($axname) by $axpos"
        current = AngleAxis(axpos*pi/180,axpz...)*current
    end

    return current
end

"""
    get_axis_poise(incif::CifContainer,axis_id,scan_id,frame_no)

Return the direction and offset of `axis_id` for `frame_no` of `scan_id`. The
axis vector of the nth axis up the stack is the result of applying
all underlying rotations to its starting vector.
"""
get_axis_poise(incif::CifContainer, axis_id, args...) = begin
    
    start_vector = get_axis_vector(incif, axis_id)  #starting vector
    start_offset = get_axis_offset(incif, axis_id)  #starting offset

    dep_chain = get_dependency_chain(incif, axis_id)

    if length(dep_chain) == 1 return start_vector, start_offset end

    # Get starting positions
    
    ax_vec = get_axis_vector.(Ref(incif), dep_chain[2:end])
    ax_offset = get_axis_offset.(Ref(incif), dep_chain[2:end])
    ax_nm,ax_types,ax_pos = get_axis_settings(incif,dep_chain[2:end], args...)

    # Apply rotations/translations in order

    current = start_vector
    current_shift = start_offset

    @debug "Shift at beginning" current_shift
    
    for (n, v, o, t, p) in zip(ax_nm, ax_vec, ax_offset, ax_types, ax_pos)
        if t == "rotation"
            @debug "Now rotating $current by $p" n
            rotmat = AngleAxis(p*pi/180,v...)
            current = rotmat * current
            current_shift = rotmat * current_shift + o
            @debug "Shift" current_shift rotmat o
        else
            @debug "Now translating $current_shift by $p along $v" n
            current_shift = current_shift + o + p * v
            @debug "After" current_shift
        end
        
    end

    return current, current_shift
end

"""
    unrotate(incif::CifContainer,pt,scan_id,frame_no)

Calculate the position of reciprocal lattice pt observed in
`scan_id` at `frame_no` when all axes are at their reference positions.
"""
unrotate(incif::CifContainer, pt, args...) = begin
    ga = reverse(get_gonio_axes(incif)[1])
    gv = get_axis_vector.(Ref(incif), ga, args...)
    gs = get_axis_settings(incif, ga, args...)[end]

    @debug "Axis vectors and settings:" gv gs
    
    # And unrotate
    
    current = pt

    @debug "Point now" current
    for (ax,pos) in zip(gv,gs)
        rot_mat = AngleAxis(-pos*pi/180,ax...)
        current = rot_mat*current
        @debug "Point now" current rot_mat
    end
    return current
end

"""
    rotate_gonio(incif::CifContainer,pt,scan_id,frame_no)

Rotate reciprocal space `pt` according to the goniometer
settings for `scan_id` at `frame_no`.
"""
rotate_gonio(incif::CifContainer,pt, args...) = begin
    ga = reverse(get_gonio_axes(incif)[1])  #base is first
    gs = get_axis_settings(incif, ga, args...)[end]
    gv = get_axis_vector.(Ref(incif),ga)    #starting values

    # Move through the axes

    current_pt = pt
    for i in 1:length(ga)
        pos = gs[i]
        rot_mat = AngleAxis(pos*pi/180,gv[i]...)
        current_pt = rot_mat*current_pt
    end

    @debug "Final axis vectors" gv

    return current_pt
end

"""
    rotate_gonio(incif::CifContainer,pt,scan_id)

Rotate reciprocal space `pt` according to the goniometer
settings for `scan_id`, setting the scan axis to zero.
Useful for comparison with cbflib.
"""
rotate_gonio(incif::CifContainer,pt,scan_id) = begin
    
    ma = get_measurement_axis(incif,scan_id)
    ga = reverse(get_gonio_axes(incif)[1])  #base is first
    gs = get_axis_settings(incif, ga, scan_id,1)[end]
    gv = get_axis_vector.(Ref(incif),ga)  #starting values

    # Move through the axes

    current_pt = pt
    for i in 1:length(ga)

        # Do not touch the measurement axis
        
        if ga[i] == ma continue end

        pos = gs[i]
        rot_mat = AngleAxis(pos*pi/180,gv[i]...)

        current_pt = rot_mat*current_pt
    end

    return current_pt
end

get_recip_point(incif::CifContainer,slow,fast,scan_id,frame_no) = begin

    pixel_coord =  get_pixel_coordinates(incif,slow,fast,scan_id,frame_no)

    # Transform to reciprocal space

    length = norm(pixel_coord)
    lambda = parse(Float64,incif["_diffrn_radiation_wavelength.value"][])
    pixel_coord = pixel_coord / (length*lambda) + [0,0,1/lambda]

    # Unrotate

    @debug "Recip point before unrotating" pixel_coord
    
    unrotate(incif,pixel_coord,scan_id,frame_no)

end

# Detector calculations

"""
    calculate_position(incif::CifContainer, axis_id, pos, scan_id, frame_no)

Return laboratory coordinates for `pos` of `axis_id` with everything positioned
for `frame_no` of `scan_id`. Only makes sense for translation axes. If
`scan_id` and `frame_no` missing, will use reference positions.
"""
calculate_position(incif::CifContainer, axis_id, pos, args...) = begin

    vector, offset = get_axis_poise(incif, axis_id, args...)
    dist = pixel_to_length(incif, axis_id, pos)
    @debug "Axis $axis_id pos $pos" dist offset
    disp = offset + dist*vector
    return disp, offset
end

"""
    pixel_to_length(incif::CifContainer, axis_id, pix;pix_origin=false)

Convert a pixel to a distance from the origin of detector coordinates.
If `pixel_origin` is true, the distance is from the centre of the origin
pixel, not the origin of the detector axes.
"""
pixel_to_length(incif::CifContainer, axis_id, pix;pix_origin=false) = begin

    c = "_array_structure_list_axis"
    
    if !(axis_id in incif["$c.axis_id"])
        throw(error("$axis_id is not a detector pixel axis"))
    end

    row = indexin([axis_id],incif["$c.axis_id"])[]
    start = pix_origin ? 0.0 : parse(Float64,incif["$c.displacement"][row])

    return start + parse(Float64,incif["$c.displacement_increment"][row])*pix
end

get_pixel_coordinates(incif::CifContainer,slow,fast,args...) = begin

    fast_ax, slow_ax = get_surface_axes(incif)
    fast_disp, f_orig = calculate_position(incif, fast_ax, fast, args...)
    slow_disp, s_orig = calculate_position(incif, slow_ax, slow, args...)

    if f_orig != s_orig
        throw(error("Detector surface axes $slow_ax, $fast_ax have different origins: $s_orig, $f_orig. Suggest making one depend on the other in axis description"))
    end

    fast_disp + (slow_disp - s_orig)
end

"""
    get_beam_centre(incif::CifContainer,scan_id,frame_no)

Return the beam centre coordinates for `frame_no` of `scan_id`,
taking into account all axis positions. Return pixel coordinates
slow, fast, and mm coordinates slow, fast as 4 values. The mm
coordinates are from the centre of the origin pixel, not the
origin of the detector coordinates.
"""
get_beam_centre(incif::CifContainer,args...) = begin

    # Follow the algorithm of cbflib

    o_coords = get_pixel_coordinates(incif,0,0,args...)
    s_step = get_pixel_coordinates(incif,1,0,args...) - o_coords
    f_step = get_pixel_coordinates(incif,0,1,args...) - o_coords
    
    # Check are linearly independent

    det = s_step[1]*f_step[2] - s_step[2]*f_step[1]
    if det == 0.0
        throw(error("Detector pixel axes are parallel"))
    end

    # Calculate distance to x=y=0 in pixel coordinates

    index = [-f_step[2]*o_coords[1] + f_step[1]*o_coords[2],
               s_step[2]*o_coords[1] - s_step[1]*o_coords[2]]/det

    fast_ax, slow_ax = get_surface_axes(incif)

    centre = pixel_to_length.(Ref(incif), [slow_ax, fast_ax], index, pix_origin=true)
    return centre, index
end

"""
    get_detector_normal(incif::CifContainer, scan_id, frame_no)

Analog of cbf routine.  Note that curved detectors must raise an error.
"""
get_detector_normal(incif::CifContainer, scan_id, frame_no) = begin

    # Follow the algorithm of cbflib

    # TODO no answer for curved detector

    get_pixel_normal(incif, 0, 0, scan_id, frame_no)
end

"""
    get_pixel_normal(incif::CifContainer, slow, fast, scan_id, frame_no)

Calculate the normal to the detector a pixel `slow, fast` when the detector
is positioned for frame `frame_no` of scan `scan_id`.
"""
get_pixel_normal(incif::CifContainer, slow, fast, scan_id, frame_no) = begin

    # As for cbflib

    pixel00 = get_pixel_coordinates(incif, slow - 0.5, fast - 0.5, scan_id, frame_no)
    pixel01 = get_pixel_coordinates(incif, slow - 0.5, fast + 0.5, scan_id, frame_no)
    pixel10 = get_pixel_coordinates(incif, slow + 0.5, fast - 0.5, scan_id, frame_no)
    pixel01 = pixel01 - pixel00
    pixel10 = pixel10 - pixel00

    normal = cross(pixel01, pixel10)

    if dot(normal, normal) <= 0
        throw(error("Cannot calculate normal to pixel $slow, $fast"))
    end

    normal = LinearAlgebra.normalize(normal)

    return normal
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
(nothing,1) is returned. The reverse operation is `bin_id_from_scan_frame`
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
    bin_id_from_scan_frame(cc::CifContainer,scan_id,frame_no)

Return the binary identifier(s) for the given `frame_no` of scan `scan_id`.
`scan_id` is ignored if only one scan is present. Return is an array of
`bin_id`.
"""
bin_id_from_scan_frame(cc,scan_id,frame_no) = begin

    # First turn a frame_no into a frame_id

    c = "_diffrn_scan_frame"
    dsf_loop = get_loop(cc,"$c.frame_number")
    if haskey(cc,"_diffrn_scan.id") && length(cc["_diffrn_scan.id"])>1
        dsff = filter(r->r["$c.scan_id"] == scan_id,dsf_loop)
    else
        dsff = dsf_loop
    end

    dsfn = filter(r->r["$c.frame_number"] == "$frame_no",dsff)

    if size(dsfn,1) > 1
        throw(error("More than one frame with same number: $dsfn"))
    end

    if size(dsfn,1) == 0
        return []
    end

    frame_id = dsfn[!,"$c.frame_id"][]

    # Then turn a scan_id / frame_id into a bin_id
    
    c = "_diffrn_data_frame"
    all_frames = get_loop(cc,"$c.binary_id")
    aff = filter(all_frames) do r
        r["$c.id"] == frame_id
    end

    # Then the binary ids are the answer!

    return aff[:,"$c.binary_id"]

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

"""
    get_scan_axis(cc,scan_id)

Return the rotation axis that moves during scan `scan_id` together with the
start, finish and increment values.
"""
get_scan_axis(cc,scan_id) = begin
    c = "_diffrn_scan_axis"
    dsl = get_loop(cc,"$c.scan_id")
    dslf = filter(row->row["$c.scan_id"] == scan_id,dsl)
    filter!(dslf) do r
        s = r["$c.angle_range"]
        !isnothing(s) && parse(Float64,s) != 0
    end
    if size(dslf,1) > 1
        throw(error("Unable to process joint axis movements during scan $scan_id"))
    elseif size(dslf,1) == 0
        throw(error("No axis rotated during scan $scan_id"))
    end
    start = parse(Float64,dslf[!,"$c.angle_start"][])
    finish = start + parse(Float64,dslf[!,"$c.angle_range"][])
    increment = parse(Float64,dslf[!,"$c.angle_increment"][])
    return dslf[!,"$c.axis_id"][],start,finish,increment
end

struct Peak
    scan_id::String
    frame_no::Int64
    slow::Float64
    fast::Float64
    intensity::Float64
end

Peak(scan_id::String, frame_no, peak_info) = begin
    x, intensity = peak_info
    Peak(scan_id, frame_no, x[2], x[1], intensity)
end

Peak(scan_id::String, frame_no::Int64, slow, fast) = Peak(scan_id, frame_no, slow, fast, 0)

Peak(scan_id::String, frame_no::Float64, slow, fast, i) = begin
    Peak(scan_id, round(Int64,frame_no) ,slow, fast, i)
end
         
intensity(p::Peak) = p.intensity
coords(p::Peak) = p.slow, p.fast
frame(p::Peak) = p.frame_no
scan(p::Peak) = p.scan_id

dist(p1::Peak, p2::Peak) = begin
    sqrt((p1.slow - p2.slow)^2 + (p1.fast - p2.fast)^2)
end

Base.show(io::IO, p::Peak) = begin
    write(io, "$(p.scan_id) $(p.frame_no) $(p.slow) $(p.fast) $(p.intensity)")
end

"""
    find_peaks(im; mindist = 10)

A primitive peak-finding algorithm. No attempt is made to find all
peaks. All pixels with values less than 1/2000 of the maximum are
set to zero, then a Niblack binarisation algorithm is applied to
produce a peak mask, which is reapplied to the original image
and local maxima found. `mindist` specifies the minimum distance
in pixels that peaks must be separated in order to be kept. The
routine returns fast, slow, max intensity
"""
find_peaks(im; mindist=10) = begin

    # Find the maximum peak (not spurion) intensity
    
    maxint = maximum(im)
    c = argmax(im)
    avg = ceil(mean(im)) #proxy for background
    @debug "Assume background is" avg
    i = 0
    while true
        bound_x_lower = max(c[1]-5,1)
        bound_x_upper = min(c[1]+5,size(im,1))
        bound_y_lower = max(c[2]-5,1)
        bound_y_upper = min(c[2]+5,size(im,2))
        view_area = im[bound_x_lower:bound_x_upper,bound_y_lower:bound_y_upper]
        @debug "Maximum for $c" maximum(view_area)
        view_area[argmax(view_area)] = 0
        if maximum(view_area) > avg
            @debug "Max nearby" maximum(view_area)
            break
        end
        im[c] = 0
        maxint = maximum(im)
        c = argmax(im)
        i+=1
        if i>100 break end   #sanity check
    end

    threshold = maxint/20
    if threshold < 10  #weak image or no spots?
        threshold = max(maxint / 2, 10)
    end
    t_im = map(x-> x > threshold ? x : eltype(im)(0), im)

    @debug "Maximum, threshold" maxint, threshold
    # binarize

    binarize!(t_im,Niblack())

    # mask original

    m_im = .*(t_im,im)

    # find local maxima and remove spurious ones

    candidates = findlocalmaxima(m_im)

    maxvals = getindex.(Ref(im),candidates)
    thresh = maximum(maxvals)/20

    @debug "Peak reject threshold" thresh
    
    c_with_int = map(candidates) do c
        bound_x_lower = max(c[1]-10,1)
        bound_x_upper = min(c[1]+10,size(m_im,1))
        bound_y_lower = max(c[2]-10,1)
        bound_y_upper = min(c[2]+10,size(m_im,2))
        view_area = m_im[bound_x_lower:bound_x_upper,bound_y_lower:bound_y_upper]
        @debug "Maximum for $c" maximum(view_area)

        # detect and remove single-pixel peaks

        view_area[findlocalmaxima(view_area)[1]] = 0
        
        c, maximum(view_area)
    end

    # Filter peaks that are too close each other or not intense enough

    keepers = Base.filter( c_with_int ) do c
        
        for (x, intensity) in c_with_int

            if intensity < thresh return false end
            
            if x == c[1] continue end

            d = sqrt((x[1] - c[1][1])^2 + (x[2] - c[1][2])^2)

            if d < mindist
                @debug "Too close" d x c
                return false
            end

        end
        true

    end

    return keepers
end

"""
    peak_to_frames(p::Peak, cc::CifContainer;single=false)

Given peak `p`, return an array of peaks that the pixel should appear,
if at all. If only peaks from the same scan should be checked,
`single` should be true.  
"""
peak_to_frames(p::Peak, cc;single=false) = begin

    # Get useful constants
    
    lambda = parse(Float64,cc["_diffrn_radiation_wavelength.value"][])
    c = "_array_structure_list"
    asl = get_loop(cc,"$c.axis_set_id")
    fast_num = filter(row->row["$c.precedence"]=="1",asl)
    fast_num = parse(Int64,fast_num[!,"$c.dimension"][])
    slow_num = filter(row->row["$c.precedence"]=="2",asl)
    slow_num = parse(Int64,slow_num[!,"$c.dimension"][])

    @debug "Peak to frames" lambda slow_num fast_num
    
    # Get zero-rotation reciprocal coordinates
    
    filename = "$(cc.original_file)"
    slow,fast = coords(p)
    scan_id = scan(p)
    frame_no = frame(p)
    recip_coords = get_recip_point(cc, slow, fast, scan_id, frame_no)

    # Loop over scans looking for intersections

    found_list = Peak[]
    all_scans = single ? [scan_id] : cc["_diffrn_scan.id"]
    for one_scan in all_scans

        # Get number of frames for later

        scan_ind = indexin([one_scan],cc["_diffrn_scan.id"])[]
        no_frames = parse(Int64,cc["_diffrn_scan.frames"][scan_ind])
        
        # Get reciprocal coordinates at start of scan
        
        begin_pt = get_scan_start(cc,recip_coords,one_scan)

        @debug "At zero scan axis:" rotate_gonio(cc,recip_coords,one_scan)

        @debug "Reference, start of scan" recip_coords begin_pt

        # Calculate rotation to intersection
        
        scan_axis,scan_begin,finish, incr = get_scan_axis(cc,one_scan)
        rot_vec = get_axis_vector(cc,scan_axis,one_scan,frame_no)
        hits = ewald_intersect(lambda,rot_vec,begin_pt)

        if hits === nothing

            @debug "No intersections for $one_scan"
            continue

        end

        rot_ang = rot_angle.(Ref(begin_pt),hits,Ref(rot_vec))

        # Check if we need to rotate in the opposite direction
        
        fn = map(rot_ang) do x
            if sign(x)==sign(incr) || isapprox(x,0,atol=0.01)
                1 + x/incr
            else
                1 + sign(incr)*(360 - abs(x))/incr
            end
        end

        @debug "Intersections for $pixel_coords in $one_scan" rot_ang fn

        # Determine which ones we will see
        
        for (h,n,r) in zip(hits,fn,rot_ang)
            if n <= no_frames
                x,y = get_detector_coords(cc,one_scan,n,h)

                @debug "Det coords" x y
                if x > 0 && y > 0 && x <= slow_num && y <= fast_num
                    @debug "Found!"
                    push!(found_list,Peak(one_scan,n,x,y, intensity(p)))
                end
            else
                @debug "Rejecting, frame no $n > $no_frames"
            end
        end
    end

    @debug "Found: $found_list" pixel_coords scan_id frame_no
    
    return found_list
end

"""
    get_scan_start(cc,pt,scan_id)

Calculate the reciprocal lattice location of `pt` at the 
start of `scan_id`
"""
get_scan_start(cc,pt,scan_id) = begin
    rotate_gonio(cc,pt,scan_id,1)
end

"""
    get_detector_coords(cc,scan_id,frame_no,r_coord)

Find the pixel coordinates on the detector positioned for `scan_id` and
`frame_no` for reciprocal space coordinate `r_coord`.
"""
get_detector_coords(cc::CifContainer,scan_id,frame_no,r_coord) = begin

    # If the ray is (0,0,1/lambda) + d((0,0,1/lambda)-r_coord)
    # And the plane is described by a normal and point on the plane

    lambda = parse(Float64,cc["_diffrn_radiation_wavelength.value"][])
    det_point = get_pixel_coordinates(cc, 0, 0, scan_id, frame_no)
    normal = get_detector_normal(cc, scan_id, frame_no)
    det_intersect = detector_intersect(r_coord, lambda, normal, det_point)
    det_coords = lab_to_det(cc, det_intersect, scan_id, frame_no)
end

"""
    lab_to_det(cc,det_point)

Calculate the detector coordinates of `det_point`, which is a point on
the detector expressed in laboratory coordinates. The coordinates are
obtained by projecting the vector from the detector origin to the
point onto the two detector axes.
"""
lab_to_det(cc,det_point,scan_id,frame_no) = begin

    # Get the origin

    det_origin = get_pixel_coordinates(cc,0,0,scan_id,frame_no)
    point_vec = det_point .- det_origin
    slow_size = get_pixel_coordinates(cc,1,0,scan_id,frame_no) .- det_origin
    fast_size = get_pixel_coordinates(cc,0,1,scan_id,frame_no) .- det_origin
    slow_dir = LinearAlgebra.normalize(slow_size)
    fast_dir = LinearAlgebra.normalize(fast_size)
    slow_coord = dot(point_vec,slow_dir)/norm(slow_size)
    fast_coord = dot(point_vec,fast_dir)/norm(fast_size)
    return slow_coord,fast_coord
end
