# Routines for checking peak locations

"""
    create_peak_image(incif,peaklist;skip=3,range=2)

Peak list contains a list of observed peaks in form
[scan,frameno,slow,fast] where `slow` and `fast` are the slow and fast
directions on the detector.  This routine calculates the predicted
location of each peak, then presents an excerpt of the relevant scans
near that point. The excerpt is created by adding all images over a
range of `range` frames each side of the predicted frame, using every
`skip` image. So the default values will use 5 images, where
each image is separated by 3 frames, for a total "frame range" of
13.  
"""
create_peak_image(incif,peaklist;skip=3,range=2,kwargs...) = begin

    # Calculated the predicted peaks
    
    predicted = map(peaklist) do p
        scan_id,frame_no,slow,fast = p
        peak_to_frames((slow,fast),scan_id,frame_no,incif)
    end

    @debug "Predicted peaks:" peaklist predicted
    
    # Create the image excerpts

    all_images = []
    for one_pred_set in predicted
        image_set =  []
        for one_pred in one_pred_set
            scan_id,frame_no,slow,fast = one_pred
            @debug "Processing:" scan_id frame_no slow fast
            one_img = create_peak_area(incif,scan_id,frame_no,slow,fast;skip=skip,range=range,kwargs...)
            push!(image_set,one_img)
        end
        push!(all_images,image_set)
    end

    # Lay them out nicely

    display_peak_table(all_images,"$(incif.original_file)"*"_peaks"*".png")
end

create_peak_area(incif,scan_id,frame_no,slow,fast;skip=3,range=2,window=10,local_version=Dict()) = begin

    # find out our maximum frame number

    if haskey(incif,"_diffrn_scan.id") && length(incif["_diffrn_scan.id"]) > 1
        flid = indexin([scan_id],incif["_diffrn_scan.id"])[]
        max_frame = incif["_diffrn_scan.frames"][flid]
    else
        max_frame = incif["_diffrn_scan.frames"][]
    end

    max_frame = parse(Int64,max_frame)
    frame_no = Int(round(frame_no))
    frame_nos = collect(max(1,frame_no-range*skip):skip:min(max_frame,frame_no+range*skip))
    frame_ids = bin_id_from_scan_frame.(Ref(incif),scan_id,frame_nos)
    filter!(x->length(x) == 1, frame_ids)  #can't handle multi-image frames for now

    if length(frame_ids) == 0
        throw(error("Frame $scan_id/$frame_no not found or consists of multiple images which is not yet supported"))
    end

    frame_ids = map(x->x[],frame_ids)
    
    @debug "Frames for peak:" frame_ids

    # accumulate all the peaks

    full_img = imgload(incif,frame_ids,local_version=local_version)
    fast_m,slow_m = size(full_img)

    # Window down

    slow = Int(round(slow))
    fast = Int(round(fast))
    
    corners = (max(1,slow-window),max(1,fast-window),min(slow_m,slow+window),min(fast_m,fast+window))

    full_img = full_img[corners[2]:corners[4],corners[1]:corners[3]]
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

    one_ht,one_wid = size(image_list[1][1])
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
