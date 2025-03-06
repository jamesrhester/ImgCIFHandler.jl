"""
    imgload(handle,::Val{:KCD})

Read an image from an Kappa CCD file.
"""
imgload(filename::AbstractString,::Val{:KCD};path=nothing,frame=nothing) = begin
    loc = open(filename,"r")
    header = read_kcd_header(loc)
    dim1 = header["x dimension"]
    dim2 = header["y dimension"]
    nbreadout = header["number of readouts"]
    expect_size = 2 * dim1 * dim2 * nbreadout

    # Calculate offset to find start of image

    fsize = filesize(filename)
    seekstart(loc)

    @debug "Skipping" basename(filename) fsize expect_size fsize-expect_size
    
    skip(loc, fsize - expect_size)

    @debug "Next byte is" peek(loc)
    
    binary = read(loc, 2 * dim1 * dim2)  #Sequence of UInt16

    @debug "Final bytes are" binary[end-1] binary[end] 
    
    better = ltoh.(reinterpret(UInt16, binary))

    i = nbreadout - 1
    while i > 0
        binary_img = read(loc, 2 * dim1 * dim2)

        @debug "Next readout, final bytes" binary_img[end-1] binary_img[end]
        
        better += ltoh.(reinterpret(UInt16, binary_img))
        i -= 1
    end

    data = reshape(better,(dim1, dim2))
    return data
end

read_kcd_header(loc) = begin
    seekstart(loc)
    full_text = readuntil(loc, " "^20)
    lines = lowercase.(split(full_text,"\n"))[2:end] #first line is ignored
    
    filter!( x-> length(x) > 0, lines)

    info = map(lines) do l
        words = strip.(split(l, "="))
        if length(words) > 1
            words[1] => words[2:end]
        else
            words[1] => words[1]
        end
    end

    d = Dict(info)

    @debug "Header info" d
    
    for n in ("x dimension", "y dimension", "number of readouts")
        d[n] = parse(UInt64, d[n][])
    end

    return d
end
