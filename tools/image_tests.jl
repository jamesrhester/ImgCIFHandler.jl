#
#  If there are issues with missing packages, please run
#  "julia install_image_tests.jl" to get them installed
#
import Pkg

# Uncomment following for release version, comment to use
# local development versions of packages

Pkg.activate(@__DIR__)

using ImgCIFHandler
using ImageInTerminal, Colors,ImageContrastAdjustment
using ImageTransformations
using ArgParse
using CrystalInfoFramework,URIs, Tar
using Luxor
using Statistics
using Dates

# Tests for imgCIF files
const test_list = []
const test_list_with_img = []
const test_full_list = []
const dictionary_checks = []

include("check_macros.jl")
include("no_image_checks.jl")
include("image_checks.jl")
include("peak_check.jl")

is_archive(u) = begin
    long_ext = u[end-6:end]
    short_ext = u[end-2:end]
    long_ext in [".tar.gz", "tar.bz2"] || short_ext in ["tgz","zip"]
end

"""
    get_archive_member_name(incif, local_archive:ImageArchive;pick=1)

    Return a local equivalent for an entry in `incif`, using `local_archive`.

"""
get_archive_member_name(incif, local_archive::ImageArchive; pick=1) = begin

    peek_image(a, incif)
    
end

test_archive_present(incif, local_archives) = begin

    messages = []
    for la in local_archives
        if ping_archive(la, incif) == nothing
            push!(messages, (false, "Unable to access $la"))
        end
    end

    return messages
    
end

@fullcheck "All members present" member_check(incif,all_archives) = begin
    messages = []
end


#   Utility routines

"""
    Find an image ID that is present locally
"""
find_load_id(incif, local_archive) = begin

    external_info = get_loop(incif,"_array_data_external_data.id")

    for ext_row in eachrow(external_info)
        for la in local_archive
            if has_local_version(la, ext_row)
                ext_data_id = getproperty(ext_row, "_array_data_external_data.id")
                vv = incif["_array_data.external_data_id"]
                pos = indexin([ext_data_id],vv)[]
                return incif["_array_data.binary_id"][pos]
            end
        end
    end

    @error "No recognised frame files found" local_archive
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
    apply_mask!(incif,im)

Detect masked pixels and change their value to zero. A masked pixel is
any pixel less than the _array_intensities.underload, or greater than
or equal to _array_intensities.overload. Assumes only a single array
type in `incif`.
"""
apply_mask!(incif,im::Array{T,2}) where T = begin
    overload = typemax(T)+1
    if haskey(incif,"_array_intensities.overload")
        overload = parse(Float64,incif["_array_intensities.overload"][])
    end
    underload = -1
    if haskey(incif,"_array_intensities.underload")
        underload = parse(Float64,incif["_array_intensities.underload"][])
    end

    @debug "Pixels >= $underload and < $overload are valid"

    # debugging information

    masked = count(x-> x >= overload || x < underload, im)
    @debug "Number of masked pixels:" masked

    if masked > 0
        for i in eachindex(im)
            im[i] = im[i] >= overload || im[i] < underload ? 0 : im[i]
        end
    end
    return im
end

"""
    improve_contrast(im;cut_ratio=1000)

Return a new image that has contrast-enhanced `im`. See
`find_best_cutoff` for meaning of `cut_ratio`
"""
improve_contrast(im;cut_ratio=1000) = begin
    clamp_low,clamp_high = find_best_cutoff(im,cut_ratio=cut_ratio)
    alg = LinearStretching(src_maxval = clamp_high)
    im_new = adjust_histogram(im,alg)
    #alg = Equalization(nbins=256,maxval = maximum(im_new))
    #im_new = adjust_histogram(im_new,alg)
    
    @debug "Max, min for adjusted image:" maximum(im_new) minimum(im_new)

    return im_new
end

"""
    create_check_image(incif,im;logscale=true,cut_ratio=1000,gravity=true,
                                peaks=false)

Create a contrast enhanced image from the matrix `im`. If `gravity` is true, 
information in `incif` is used, if available, to calculate the rotation
of the image so that down on the image 
becomes down in real space. Returns the unrotated image and the rotation.
If `peaks` is true, peaks are searched for.
"""
create_check_image(incif,im;logscale=true,cut_ratio=1000,gravity=true,do_peaks=false) = begin

    # First detect mask and set to zero

    apply_mask!(incif,im)

    # Get a list of peaks for later

    peaks = []
    
    if do_peaks
        peaks = find_peaks(im)
    end

    @debug "Peaks found at" peaks
    
    # Then try to improve contrast
    if maximum(im) > 1.0
        im = im/maximum(im)
    end

    im_new = improve_contrast(im)

    # Make sure the image is that seen from the front of the
    # detector. The image is presented by image display tools
    # with the origin at top left, with fast direction down.
    # So if slow is across, fast down, the cross product is
    # the view direction. View direction should be [0,0,-1]
    # for away from source in imgCIF coordinates.

    fast_v,slow_v = get_axis_vector.(Ref(incif),get_surface_axes(incif))
    view_d = slow_v[1]*fast_v[2] - slow_v[2]*fast_v[1] #z comp of cross product
    do_transpose = view_d > 0   # need to transpose
    
    # Adjust geometry if we know gravity

    if !gravity || !("gravity" in incif["_axis.equipment"]) return im_new,do_transpose,0,peaks end

    gravity = indexin(["gravity"],incif["_axis.equipment"])[]
    grav_vec = get_axis_vector(incif,incif["_axis.id"][gravity])
    norm_grav = sqrt(grav_vec[1]^2+grav_vec[2]^2)
    grav_vec = grav_vec/norm_grav
    
    corner_loc = get_pixel_coordinates(incif,0,0)
    if !do_transpose
        fast_dir = get_pixel_coordinates(incif,0,1) - corner_loc
    else
        fast_dir = get_pixel_coordinates(incif,1,0) - corner_loc
    end
    
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

    return im_new,do_transpose,rot*90,peaks
end

"""
    annotate_check_image(image,transpose,rotation,beam_centre,names,filename; border=30, hsize=512,peaks=[],scan_id=nothing,frame_no=nothing)

Add axes and beam centre to `image`, saving the result in `filename`. Metadata 
relating to the image is read from `incif`, and the final image is rotated by `rot`
and (prior to this) transposed if `transp` is true.
`border` is the size of the border around
the image. `hsize` is the target maximum horizontal dimension of the check image,
in pixels. `peaks` is a list of the most intense peaks found in the image,
which should have circles drawn around them. If `scan_id` and `frame_no` are
not `nothing`, they are used to calculate the beam centre for the corresponding
detector position, otherwise the beam centre for zero detector motion is used.
"""
annotate_check_image(im, transp, rot, incif;border=30,hsize=512,scan_id=nothing,frame_no=nothing,peaks=[]) = begin

    # In the following there are three images: the original image `im`,
    # the transposed rotated scaled image `im_new` and the reference image
    # next to it `im_ref`
    
    width_orig = size(im,2)
    height_orig = size(im,1)
    
    # Transpose and rotate the image matrix

    if transp
        im_new = rotl90(permutedims(im),div(rot,90))
    else
        im_new = rotl90(im,div(rot,90))
    end
    
    # Scale to nearest to requested size that is multiple of 2
    # But not smaller

    im_ref = im
    scale_factor = hsize/size(im_new,2)
    if scale_factor < 1
        scale_factor = floor(-1*log2(scale_factor))
        @debug "Scaling images by 2^$scale_factor"
        for i in 1:Int(scale_factor)
            im_new = restrict(im_new)
            im_ref = restrict(im_ref)
        end
        scale_factor = 2.0^(-1*scale_factor)
    else
        scale_factor = 1
    end
    
    width_new = size(im_new,2)   #slow direction
    height_new = size(im_new,1)  #fast direction

    width_ref = size(im_ref,2)
    height_ref = size(im_ref,1)   #in display, height is fast direction
    
    # Get the beam centre

    if isnothing(scan_id)
        _,cent_pix = get_beam_centre(incif)
    else
        _,cent_pix = get_beam_centre(incif,scan_id,frame_no)
    end

    # Transform beam centre coordinates; scale and rotation

    new_slow_c,new_fast_c = calc_new_coords(cent_pix,transp,rot,scale_factor,width_orig,height_orig)

    cent_pix = cent_pix * scale_factor

    @debug "Old, new height, width" height_orig width_orig height_new width_new

    # Start the drawing
    
    Luxor.Drawing(width_new+width_ref+border*4,maximum((height_new,height_ref))+border*2+20, :rec)
    background("white")

    # Place rotated and original images

    placeimage(Gray.(im_new),Point(border,border))
    placeimage(Gray.(im_ref),Point(border*3 + width_new,border))
    setcolor("black")
    label("Laboratory frame",0,Point(border+width_new/2,border*2+height_new+10))
    label("Original",0,Point(border*3+width_new+width_ref/2,border*2+height_ref+10))
    
    # draw a red beam centre
    
    setcolor("red")
    box(Point(new_slow_c+border,new_fast_c+border),4,4,:fill)
    box(Point(border*3+width_new+cent_pix[1],border+cent_pix[2]),4,4,:fill)

    # draw the axis labels
    
    ax_labels = get_surface_axes(incif)
    if transp ax_labels = reverse(ax_labels) end
    draw_axes(height_new,width_new,rot,ax_labels)

    peaks = [x[1] for x in peaks]   #ignore intensity
    draw_peaks(width_orig,height_orig,scale_factor,transp,rot,peaks)

    # output
    
    snapshot(fname="$(incif.original_file)"*".png")
    finish()
end

""" 
    calc_new_coords(old,transp,rot,scale,width,height)

Calculate the coordinates on the canvas resulting from a possible
transpose (boolean `transp`) followed by rotation in
degrees of `rot` and scaling of `scale` where the rotation is around
the centre of an image of `height` and `width`. Original coordinates
in `old` as (slow,fast) where the slow direction is conventionally
horizontal and the origin is at pixel (1,1).
"""
calc_new_coords(old,transp,rot,scale,width,height) = begin

    @debug "Rotation" rot isodd(div(rot,90))
    if isodd(div(rot,90)) && !transp || transp && iseven(div(rot,90))
        new_width = height*scale
        new_height = width*scale
    else
        new_width = width*scale
        new_height = height*scale
    end

    trans_mat = transp ? [0 1; 1 0] : [1 0; 0 1]
    slow,fast = old
    slow = slow - width/2
    fast = fast - height/2

    rot_mat = [cosd(rot) sind(rot); -sind(rot) cosd(rot)]
    new_slow,new_fast = scale*rot_mat*trans_mat*[slow,fast] + [new_width/2,new_height/2] 

    @debug "Transpose" transpose
    @debug "Rotating $old by $rot" rot_mat
    @debug "New coords" new_slow new_fast
    @debug "New width, height" new_width new_height
    
    return (new_slow,new_fast)
end

draw_axes(height,width,angle,names;border=30) = begin
    fast,slow = names
    @debug "Fast, slow axis names" fast slow
    @debug "Rotate a/c wise by " angle
    
    # rotate and origin to corner

    origin(width/2+border,height/2+border)

    # show image centre
    
    setcolor("green")
    cropmarks(Point(0,0),4,4)

    # move to axis origin

    rotate(deg2rad(-1*angle))

    width_new = abs(cosd(angle)* width + sind(angle)* height)
    height_new = abs(-sind(angle)* width + cosd(angle)* height)

    translate(Point(-(width_new+border)/2,-(height_new+border)/2))
    
    # draw some lines

    setcolor("black")
    
    # The X coordinate is across (i.e. the slow png direction)

    arrow(O,Point(width_new/2,0))

    # And Y is down (i.e. the fast png direction)

    arrow(O,Point(0,(height_new/2.0)))

    # annotate the lines

    translate(Point(width_new/4.0,0))

    # switch the label on the bottom to get the text upright
    if angle == 180 # slow axis on bottom
	rotate(pi)
	label(slow,:N)
	rotate(pi)
    else
	label(slow,:S)
    end
    translate(Point(-width_new/4.0,height_new/4.0))
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

draw_peaks(width,height,scale,transp,angle,peaks;border=30) = begin

    @debug "Drawing $(length(peaks)) peaks"
    
    setcolor("blue")
    origin(Point(border,border))   #top left
    for op in peaks
        slow,fast = calc_new_coords((op[2],op[1]),transp,angle,scale,width,height)
        circle(Point(slow,fast),5,:stroke)
    end 
end

show_check_image(im::AbstractArray,transp,rot) = begin
    println("Image for checking")
    if transp im = permutedims(im) end
    display(Gray.(rotl90(im,div(rot,90))))
    println("\n")
end

"""
Find the best value for displaying a diffraction image, calculated as the first
intensity bin that has 1000 times less points than the highest. Based on the
logic that the highest value will be a "typical" background point. Skip any
that correspond to negative values as these are likely to have a different
meaning.  Drop any points that are separated by 40 bins from the next
lowest point.
"""
find_best_cutoff(im;cut_ratio=1000) = begin
    edges,bins = ImageContrastAdjustment.build_histogram(im)
    dropped = 0
    while sum(bins[end-10:end-1]) == 0 && dropped < 10
        @debug "Dropping high outlier" maximum(im)
        im[argmax(im)] = median(im)
        dropped += 1
        edges,bins = ImageContrastAdjustment.build_histogram(im)
    end

    @debug "Finished removing high outliers" edges bins
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
    run_img_checks(incif;images=false,always=false,full=false,max_down=1,connected=false,pick=1,subs=Dict(),savepng=false,accum=1,skip=false,peak_check=false,peakvals=[], grid=["detx","dety",1.0,3])

Test `incif` for conformance to imgCIF consistency requirements. Meaning of arguments:
`images`: run checks on downloaded image
`always`: always run checks on images, even if previous checks fail
`full`: download all referenced archives in full
`max_down`: maximum automatic download size
`connected`: perform checks that require an internet connection
`pick`: archive member to download for image checks
`subs`: dictionary of uri -> local file correspondences to avoid downloading
`savepng`: output annotated check image to a file
`accum`: accumulate this many frames to make the check image
`skip`: skip non-image checks
`peak_check`: generate an image for checking peak locations
`peakvals`: list of peaks to include in check
`grid`: generate a grid of peak checks
"""
run_img_checks(incif;images=false,always=false,full=false,max_down=1,connected=false,pick=1,subs=Dict(),savepng=false,accum=1,skip=false,peak_check=false,peakvals=[], grid=[1, 0]) = begin
    ok = true

    if !skip
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
    end
    
    # Test archive access

    println("Testing presence of archive:")

    local_archive = create_archives(incif, subs = subs)

    @debug "Created local archive" local_archive

    if connected
        print("\nTesting: All archives are accessible: ")
        
        ok = ok & verdict(test_archive_present(incif, local_archive))

    end
    
    # Test with an image

    testimage = [[]]

    if (ok && images) || always
        println("\nRunning checks with downloaded images")
        println("="^40*"\n")

        # Choose image(s) to load
        
        peek_image(first(local_archive), incif, max_down = max_down)

        load_id = find_load_id(incif, local_archive)

        if accum > 1
            load_ids = get_id_sequence(incif,load_id,accum)
        else
            load_ids = [load_id]
        end
        
        try
            testimage = imgload(incif, load_ids, local_archive)
        catch e
            @debug e
            verdict([(false,"Unable to access image $load_id: $e")])
            rethrow()
        end

        @debug "Loaded image " load_id
        
        if testimage isa String   #we have failed
            verdict([(false,testimage)])
            return (false,[[]])
        end

        # Output an image
        
        new_im,transp,rot,peaks = create_check_image(incif,testimage,logscale=false, do_peaks = false)
        scan_id,frame_no = ImgCIFHandler.scan_frame_from_bin_id(load_id,incif)
        
        imgfn = nothing
        if savepng
            annotate_check_image(new_im,transp,rot,incif,scan_id=scan_id,frame_no=frame_no, peaks = peaks, hsize = 1024)
        else
            show_check_image(new_im,transp,rot)
        end
        
        # Run the image checks
        
        for (desc,one_test) in test_list_with_img
            print("\nTesting image $load_id: $desc: ")
            ok = ok & verdict(one_test(incif,testimage,load_id))
        end

        # Generate the peak checking images

        if peak_check

            # Don't do it unless we have local files

            generate_peak_check(incif, peakvals = peakvals, cached = local_archive, grid = grid)

         end
    end
    
    # Tests requiring fully-downloaded archives

    if full
        for (desc,one_test) in test_full_list
            print("\nTesting full archive: $desc:")
            ok = ok & verdict(one_test(incif, local_archive))
        end
    end
    return (ok,testimage)
end

"""
 Generate a grid of peak check images (or just one) showing how well the
 predicted peak positions match. Note our hard-coded axis names.
"""
generate_peak_check(incif; peakvals = [], cached = [], grid = [1, 0]) = begin

    fast_name, slow_name = get_surface_axes(incif)

    detx_row = indexin([fast_name], incif["_axis.id"])
    if isnothing(detx_row[1])
        throw(error("No axis named $fast_name found in axis list but a detector axis
has this name"))
    end

    detx_row = detx_row[]
    
    x_offset = parse(Float64, incif["_axis.offset[1]"][detx_row])
    y_offset = parse(Float64, incif["_axis.offset[2]"][detx_row])

    predicted = Vector{Vector{Peak}}[]
    
    for grid_x_pos in -grid[2]:grid[2]

        incif["_axis.offset[1]"][detx_row] = "$(x_offset + grid_x_pos * grid[1])"
        
        @debug "Grid x now" incif["_axis.offset[1]"][detx_row]
        
        for grid_y_pos in -grid[2]:grid[2]
            
            incif["_axis.offset[2]"][detx_row] = "$(y_offset + grid_y_pos * grid[1])"

            @debug "Grid y now" incif["_axis.offset[2]"][detx_row]
            
            new_predicted = accumulate_peaks(incif, peakvals=peakvals, cached=cached)

            @debug "Predicted peaks" new_predicted

            push!(predicted, new_predicted)
            
        end

    end

    incif["_axis.offset[1]"][detx_row] = "$(x_offset)"
    incif["_axis.offset[2]"][detx_row] = "$(y_offset)"

    # Now plot them all in a nice grid

    create_peak_image(incif, predicted, cached=cached, skip=-1)

end

read_peak_file(filename) = begin
    lines = readlines(filename)
    filter!(x -> length(x) > 0 && x[1] != '#', lines)
    vals = []

    for l in lines
        append!(vals, split(l))
    end

    peaks = Peak[]
    for v in 1:4:length(vals)
        push!(peaks, Peak(String(vals[v]), parse(Int64,vals[v+1]),
                          parse(Float64, vals[v+3]), parse(Float64, vals[v+2])
                          )
              )
    end
    return peaks
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
        "-m", "--max-download"
        help = "Maximum archive size to download automatically, in megabytes"
        nargs = 1
        default = [1]
        arg_type = Int64
        "-n", "--no-internet"
        help = "Do not check raw image archive existence or contents"
        nargs = 0
        "-p", "--pick"
        nargs = 1
        default = [1]
        arg_type = Int64
        help = "Use this entry number in the archive for checking"
        "--peaks"
        nargs = 0
        help = "Generate a peak check image. Implies -f (full download). Suggest also using -s if a fully downloaded archive is locally available. See also option --peakval"
        "--peakval"
        nargs = 4
        metavar = ["scan","frame","fast","slow"]
        action = "append_arg"
        help = "Coordinates of a peak to include in peak check. Implies --peaks. <fast> is the fast direction on the detector, typically horizontal. If only one scan, <scan> is
ignored (but must be provided)."
        "--peak-file"
        nargs = 1
        help = "Name of file containing peak information. Peaks information is given as scanid frameno fast slow. All values are whitespace separated. Lines beginning with `#` are ignored."
        "-g", "--grid"
        nargs = 2
        metavar = ["spacing", "number"]
        default = ["1.0", "0"]
        help = "Create a grid of peak checks, varying the x,y offsets of detx by <spacing>, creating <number> steps in each direction"
        "-s", "--sub"
        nargs = 2
        metavar = ["original","local"]
        action = "append_arg"
        help = "Use <local> file in place of archive URI <original> (for testing). Use -f to access whole file. All URIs in <filename> must be provided, option may be used multiple times. If <local> is a directory, its contents are assumed to be the unpacked archive."

        "-l", "--local"
        nargs = 1
        help = "Use the directory given in <local> in place of all URIs found in the file. This is the equivalent of -s when there is one URI. Useful if all URIs have been unpacked into directory <local>. For non-archive URLs, this directory is mapped to the common part of all URLs listed in the imgCIF file."
        
        "--skip"
        nargs = 0
        help = "Only check actual images"
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
    incif = Cif(parsed_args["filename"])
    if isnothing(parsed_args["blockname"])
        blockname = first(incif).first
    else
        blockname = parsed_args["blockname"]
    end

    subs = Dict(parsed_args["sub"])

    if !isempty(parsed_args["local"])
        subs["Universal"] = parsed_args["local"][]
    end

    # Convert local directories to absolute paths

    for (k,v) in subs
        subs[k] = abspath(expanduser(v))
    end
    
    println("\n ImgCIF checker version 2025-03-07\n")
    println("Checking block $blockname in $(incif.original_file)\n")
    println(now())
    if parsed_args["dictionary"] != [""]
    end

    @debug "Parsed arguments" parsed_args subs
    
    # Fix loops

    fix_loops!(incif[blockname])

    # Fix logic

    if parsed_args["output-png"] parsed_args["check-images"] = true end
    if length(parsed_args["peakval"])>0
        parsed_args["peaks"] = true
    end

    # Get peaks
    
    peakvals = map(parsed_args["peakval"]) do p
        Peak(p[1], parse(Int64,p[2]) , parse(Float64, p[4]),
             parse(Float64, p[3]))
    end

    if parsed_args["peak-file"] != []
        append!(peakvals, read_peak_file(parsed_args["peak-file"][]))
    end

    @debug "Peaks" peakvals
    
    if parsed_args["peaks"] || parsed_args["peak-file"] != []
        parsed_args["output-png"] = true
        parsed_args["full-download"] = true
    end

    grid = parse.(Float64, parsed_args["grid"])
    
    result,img = run_img_checks(incif[blockname],
                                images=parsed_args["check-images"],
                                always=parsed_args["always-check-images"],
                                full = parsed_args["full-download"],
                                max_down = parsed_args["max-download"][] * 1e6,
                                connected = !parsed_args["no-internet"],
                                pick = parsed_args["pick"][],
                                accum = parsed_args["accumulate"][],
                                subs = subs,
                                savepng = parsed_args["output-png"],
                                skip = parsed_args["skip"],
                                peak_check = parsed_args["peaks"],
                                peakvals = peakvals,
                                grid = grid
                                )
    println("\n====End of Checks====")
    if result exit(0) else exit(1) end
end

    
