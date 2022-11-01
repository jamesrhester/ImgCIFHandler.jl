# Routines for checking peak locations
"""
    accumulate_peaks(incif; peakvals = [], maxframes=10, maxpeaks=10, local_version = Dict(), cached = Dict())

Work through up to `maxframes` images from incif until
up to `maxpeaks` peaks have been found that have predicted
positions in at least one of the other frames. Use `local_version`
and `cached` as a source for images. If `peakvals` is not empty,
use those as a source of peaks. Each entry in `peakvals` is
of the form `scan, frameno, fast, slow` .
"""
accumulate_peaks(incif; peakvals = [], scan_id = nothing, maxframes = 10, maxpeaks = 10, local_version = Dict(), cached = Dict()) = begin

    
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

            full_img = imgload(incif, frame_ids[fi], local_version=local_version, cached = cached)
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
            merge_peak!(good_peaks, p, incif)
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
    create_peak_image(incif,peaklist::Vector{Vector{Peak}};skip=3,range=2)

Peak list contains a list of observed peaks and equivalents in form
for each independent peak. This
routine presents an excerpt of the relevant scans near the predicted
positions. The excerpt is created by adding all images over a range of
`range` frames each side of the predicted frame, using every `skip`
image. So the default values will use 5 images, where each image is
separated by 3 frames, for a total "frame range" of 13.  
"""
create_peak_image(incif,predicted::Vector{Vector{Peak}};skip=-1,range=-1,kwargs...) = begin

    if skip < 0 || range < 0
        skip, range = guess_skip_range(incif)
    end
    
    # Create the image excerpts

    all_images = []
    for one_pred_set in predicted
        image_set =  []
        for p in one_pred_set
            @debug "Processing:" p
            one_img = create_peak_area(incif,p;skip=skip,range=range,kwargs...)
            push!(image_set,one_img)
        end
        push!(all_images,image_set)
    end

    # Lay them out nicely

    display_peak_table(all_images,"$(incif.original_file)"*"_peaks"*".png")
end

create_peak_area(incif,p::Peak;skip=3,range=2,window=10, local_version=Dict(), cached=Dict()) = begin

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

    # accumulate all the peaks

    full_img = imgload(incif,frame_ids,local_version=local_version, cached = cached)
    fast_m,slow_m = size(full_img)

    # Window down and remove negative values

    slow, fast = coords(p)
    slow = Int(round(slow))
    fast = Int(round(fast))
    
    corners = (max(1,slow-window),max(1,fast-window),min(slow_m,slow+window),min(fast_m,fast+window))

    full_img = full_img[corners[2]:corners[4],corners[1]:corners[3]]
    for i in eachindex(full_img)
        full_img[i] = full_img[i] < 0 ? 0 : full_img[i]
    end
                     
    full_img/maximum(full_img)
end

"""
    display_peak_table(image_list;border=30)

`image_list` contains a list of image lists. Each list should be laid out
horizontally. All images must be the same
dimensions.
"""
display_peak_table(image_list,fname) = begin

    # Calculate layout parameters based on image size

    one_ht = one_wid = 0
    for row in image_list
        for cell in row
            h, w = size(cell)
            one_ht = h > one_ht ? h : one_ht
            one_wid = w > one_wid ? w : one_wid
        end
    end
    
    border = Int(round(one_wid/2))

    max_len = maximum(length.(image_list))

    @debug "For peak table" one_ht one_wid max_len
    Luxor.Drawing((one_wid+border)*max_len + border,
                  (one_ht+border)*length(image_list) + border,
                  fname)
    background("white")

    for (i,row) in enumerate(image_list)
        for (j,col) in enumerate(row)
            placeimage(Gray.(col),Point(2*border+(j-1)*one_wid,2*border+(i-1)*one_ht))
        end
    end

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
    range = ceil(Int64, 5.0/step)

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
