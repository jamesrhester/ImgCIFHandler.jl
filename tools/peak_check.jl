# Routines for checking peak locations
"""
    accumulate_peaks(incif; peakvals = [], maxframes=10, maxpeaks=10, cached = Vector{ImageArchive})

Work through up to `maxframes` images from incif until up to
`maxpeaks` peaks have been found that have predicted positions in at
least one of the other frames. Use `cached` as a source for images. If
`peakvals` is not empty, use those as a source of peaks. Each entry in
`peakvals` is of the form `scan, frameno, fast, slow` .  """
accumulate_peaks(incif; peakvals = [], scan_id = nothing, maxframes =
10, maxpeaks = 10, cached = Vector{ImageArchive}) = begin

    good_peaks = Vector{Peak}[]

    if length(peakvals) == 0
        
        # Choose the scan and frames (first 10 basically)

        if isnothing(scan_id)
            scan_id = incif["_diffrn_scan.id"][1]
        end
        
        no_frames = parse(Int64,incif["_diffrn_scan.frames"][1])
        no_frames = min(maxframes, no_frames)

        # Find identifiers

        frame_ids = bin_id_from_scan_frame.(Ref(incif),scan_id, 1:no_frames)
        filter!(x->length(x) == 1, frame_ids)  #can't handle multi-image frames for now

        if length(frame_ids) == 0
            throw(error("Frame $scan_id/$frame_no not found or consists of multiple images which is not yet supported"))
        end

        frame_ids = map(x->x[],frame_ids)

        # Now get the peaks
        
        for fi in 1:no_frames

            full_img = imgload(incif, frame_ids[fi], cached)
            apply_mask!(incif, full_img)
            peaks = Peak.(scan_id, fi, find_peaks(full_img))
            for p in peaks
                merge_peak!(good_peaks, p, incif)
            end

            if length(good_peaks) > maxpeaks
                break
            end
        end
        
    else    # we have been provided with peaks
        
        for p in peakvals

            # Refine target frame
            
            _, new_p = create_peak_area(incif, p;skip=1, range=1,
                                                cached = cached)
            
            @debug "Old, new peaks:" p new_p

            merge_peak!(good_peaks, new_p, incif)
        end

    end
            
    @debug "Found $(length(good_peaks)) multi peaks" good_peaks

    # tabulate the peaks

    println("Scan Frame slow fast")
    for each_row in good_peaks
        println("$(first(each_row))")
        for each_equiv in each_row[2:end]
            println("   $each_equiv")
        end
    end
    
    return good_peaks
end

"""
    If `new_peak` is not a previously seen peak from another frame, and it
    has predicted positions in other scans, add it in.
"""
merge_peak!(good_peaks::Vector{Vector{Peak}}, new_peak::Peak, incif) = begin

    @debug "Is it unique?" new_peak
    
    for g in 1:length(good_peaks)
        for gg in 1:length(good_peaks[g])
            orig_peak = good_peaks[g][gg]
            if abs(frame(orig_peak)-frame(new_peak)) < 5
                if dist(new_peak, orig_peak) < 5
                    # too close
                    if intensity(new_peak) > intensity(orig_peak)
                        # replace with this one
                        @debug "Found better peak" new_peak orig_peak
                        predicted = peak_to_frames(new_peak, incif)
                        if length(predicted) > 1
                            good_peaks[g] = predicted
                        else
                            @warn "Not replacing peak" orig_peak new_peak
                        end
                    end
                    @debug "Duplicate peak" new_peak orig_peak
                    return
                end
            end
        end
    end
    predicted = peak_to_frames(new_peak, incif)
    if length(predicted) > 1
        push!(good_peaks, peak_to_frames(new_peak, incif))
    end
end

"""
    create_peak_image(incif,peaklist::Vector{Vector{Peak}};skip=1,range=5)

Peak list contains a list of observed peaks and equivalents in form
for each independent peak. This
routine presents an excerpt of the relevant scans near the predicted
positions. The excerpt is created by adding all images over a range of
`range` frames each side of the predicted frame, using every `skip`
image. So the default values will use 5 images, where each image is
separated by 1 frames.  
"""
create_peak_image(incif, predicted::Vector{Vector{Vector{Peak}}}; redo=true, skip=1,range=5,kwargs...) = begin

    if skip < 0 || range < 0
        skip, range = guess_skip_range(incif)
    end
    
    # Create the image excerpts

    all_images = []
    
    for one_grid_point in predicted

        group_images = []

        for one_pred_set in one_grid_point
            image_set =  []
            for p in one_pred_set
                @debug "Processing:" p
                one_img, max_peak = create_peak_area(incif, p;skip=skip,range=range,kwargs...)
                @debug "max peak" max_peak
                if frame(max_peak) != frame(p)
                    @warn "Max peak intensity $(intensity(max_peak)) in frame $(frame(max_peak)), not $(frame(p))"
                end
                push!(image_set,one_img)
            end
            push!(group_images,image_set)
        end
        push!(all_images, group_images)
    end

    # Lay them out nicely

    # Find the beam centre coords
    
    detx_pos = indexin(["detx", "ele1_x"], incif["_axis.id"])
    if isnothing(detx_pos[1])
        detx_pos = detx_pos[2]
    else
        detx_pos = detx_pos[1]
    end
    
    cent_coords = incif["_axis.offset[1]"][detx_pos], incif["_axis.offset[2]"][detx_pos]
    
    caption = "$(basename(incif.original_file)) Centre $(cent_coords[1]) $(cent_coords[2])"
    display_peak_table(all_images,"$(incif.original_file)"*"_peaks"*".png", caption = caption)
end

create_peak_area(incif,p::Peak;skip=1,range=5,window=10, cached = Vector{ImageArchive}) = begin

    # find out our maximum frame number

    scan_id = scan(p)
    if haskey(incif,"_diffrn_scan.id") && length(incif["_diffrn_scan.id"]) > 1
        flid = indexin([scan_id],incif["_diffrn_scan.id"])[]
        max_frame = incif["_diffrn_scan.frames"][flid]
    else
        max_frame = incif["_diffrn_scan.frames"][]
    end

    max_frame = parse(Int64,max_frame)
    
    frame_no = frame(p)
    frame_nos = collect(max(1,frame_no-range*skip):skip:min(max_frame,frame_no+range*skip))
    frame_ids = bin_id_from_scan_frame.(Ref(incif), scan_id, frame_nos)
    filter!(x->length(x) == 1, frame_ids)  #can't handle multi-image frames for now

    if length(frame_ids) == 0
        throw(error("Frame $scan_id/$frame_no not found or consists of multiple images which is not yet supported"))
    end

    frame_ids = map(x->x[],frame_ids)
    
    @debug "Frames for peak:" frame_ids

    slow, fast = coords(p)
    slow = Int(round(slow))
    fast = Int(round(fast))
    full_img = imgload(incif, frame_ids[1], cached)
    fast_m, slow_m = size(full_img)

    # accumulate all the peaks

    frame_max = -1
    max_val = -1
    max_coords = fast, slow
    corners = (max(1,slow-window),max(1,fast-window),min(slow_m,slow+window),min(fast_m,fast+window))
    acc_img = zeros(corners[4]-corners[2]+1, corners[3]-corners[1]+1)

    @debug "Corners" corners
    
    # find an overload value

    if haskey(incif,"_array_intensities.overload")
        overload = parse(Float64, incif["_array_intensities.overload"][])
    else
        overload = typemax(eltype(full_img))
    end
    
    for fid in frame_ids
        full_img = imgload(incif, fid, cached)

    # Window down and remove negative/overload/mask values

        full_img = full_img[corners[2]:corners[4],corners[1]:corners[3]]
        for i in eachindex(full_img)
            full_img[i] = (full_img[i] < 0 || full_img[i] >= overload ) ? 0 : full_img[i]
        end

        max_img = maximum(full_img)
        if max_img > max_val
            frame_max = fid
            max_val = max_img
            max_coords = Tuple(argmax(full_img))
            @debug "New peak pos" fid max_coords
        end

        acc_img += full_img
    end

    @debug "Range, max" frame_ids frame_max

    frame_max = indexin([frame_max],frame_ids)[]
    
    max_peak = Peak(scan_id, frame_nos[frame_max], max_coords[2]+corners[1], max_coords[1] + corners[2], max_val)
    acc_img/maximum(acc_img), max_peak
end

"""
    display_peak_table(image_list;border=30)

`image_list` contains a list of image lists. Each list should be laid out
horizontally. All images are assumed to be the same dimensions.
"""
display_peak_table(image_list, fname; caption = "") = begin

    if length(image_list) == 0
        @warn "No peaks found, use --peakval to manually indicate peak positions"
        return
    end
    
    # Calculate layout parameters based on image size

    # Get the longest and tallest entries
    # Outer dimension is grid step, next is row, inner is across
    
    max_down = maximum(length.(image_list))
    max_across = maximum( map( x-> maximum(length.(x)), image_list))
    side = sqrt(length(image_list))

    @debug "Max width, height for $side x $side grid" max_across max_down
    
    one_ht = one_wid = 0
    test_list = first(image_list)
    for row in test_list
        for cell in row
            h, w = size(cell)
            one_ht = h > one_ht ? h : one_ht
            one_wid = w > one_wid ? w : one_wid
        end
    end
    
    border = Int(round(one_wid))

    @debug "For peak table" one_ht one_wid border

    # Drawing is an n x n grid of peak displays, each max_across wide
    # and max_down high. Add 20 at end for caption.
    
    Luxor.Drawing((one_wid+border)*max_across*side + border,
                  (one_ht+border)*max_down*side + border + 20,
                  fname)
    
    background("white")

    @debug "Picture dimensions" (one_wid + border)*max_across*side + border (one_ht + border)*max_down*side + border

    sethue("red")
    fontsize(12)

    for (n, one_group) in enumerate(image_list)
        
        grid_x, grid_y = divrem(n-1, side)
        grid_x *= (one_wid+border)*max_across + border
        grid_y *= (one_ht + border)*max_down + border

        overall_score = 0
        for (i,row) in enumerate(one_group)
            row_pos = 2*border + (i-1)*one_ht + grid_y
            row_score = 0.0
            for (j,col) in enumerate(row)
                col_pos = 2*border+(j-1)*one_wid + grid_x
                placeimage(Gray.(col),Point(col_pos,row_pos))

                # score it

                score = score_image(col)
                text("$score", Point(col_pos, row_pos))
                row_score += score
            end

            # output the average score
            row_av = round(row_score/length(row), digits=2)
            fontsize(14)
            text("$(row_av)", Point(grid_x, row_pos + one_ht/2))
            fontsize(12)
            overall_score += row_av
        end

        full_score = round(overall_score/length(one_group), digits=2)
        fontsize(14)
        sethue("blue")
        text("$full_score", Point(grid_x, grid_y + 2*border + length(one_group)*one_ht))
        sethue("red")
        fontsize(12)
    end

    # add caption
    
    sethue("black")
    fontsize(16)
    text(caption, Point(border, (one_ht + border)*max_down*side + border + 10))

    # output

    finish()

    preview()
    # Luxor.shapshot(fname=fname)
end

"""
    guess_skip_range(incif::CifContainer)

Guess values for the total range over which to accumulate images, and how
many to skip to cover this range.
"""
guess_skip_range(incif) = begin

    # We assume a range of +/- 5 degrees is sufficient

    all_steps = incif["_diffrn_scan_axis.angle_increment"]
    all_steps = filter(x -> !ismissing(x) && !isnothing(x), all_steps)
    all_steps = [ abs(parse.(Float64, a)) for a in all_steps ]
    step = maximum(all_steps)
    range = ceil(Int64, 1.0/step)

    # Now we want no more than 5 images either side but assume a skip
    # of more than 3 might miss a peak

    @debug "Initial range" range
    
    if range > 5
        skip = ceil(Int64, range / 5.0)
        if skip > 3
            skip = 3
        end
        range = round(Int64, range / skip)
    else
        skip = 1
    end
    
    return skip, range
end

"""
Assign a number to an image, scoring how well it shows a peak. If
the peak is too weak, zero. If it is too close to the edges, 0.5.
Otherwise 1.
"""
score_image(img)  = begin
    md = median(img)
    mx = maximum(img)
    if mx < 2*md
        return 0.0
    end
    
    if mx > 2*md &&  mx < 3*md
        return 0.5
    end

    amax = argmax(img)
    margin = 0.25 .* size(img)  #outer quarter
    for idx in (1,2)

        @debug amax margin size(img)
        if amax[idx] < margin[idx] || amax[idx] > size(img)[idx] - margin[idx]
            return 0.5
        end
    end

    return 1 
end

"""
Experimental routine for finding peaks
"""
peak_find(img) = begin

    # Segment

    segs = fast_scanning(img, median(img))
    
end
