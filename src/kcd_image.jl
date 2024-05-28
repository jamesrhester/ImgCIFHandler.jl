"""
    imgload(handle,::Val{:KCD})

Read an image from an Kappa CCD file.
"""
imgload(filename::AbstractString,::Val{:KCD};path=nothing,frame=nothing) = begin
    loc = open(filename,"r")
    header = read_kcd_header(loc)
    dim1 = header["X dimension"]
    dim2 = header["Y dimension"]
    nbreadout = header["Number of readouts"]
    expect_size = 2 * dim1 * dim2 * nbreadout

    # Seek from end of file to find start of image

    seekend(loc)
    skip(-1*expect_size)
    
    binary = read(loc, 2 * dim1 * dim2)  #Sequence of UInt16
    better = reinterpret(UInt16, binary)

    i = nbreadout - 1
    while i > 0
        binary = read(loc, 2 * dim1 * dim2)
        better += reinterpret(UInt16, binary)
        i -= 1
    end
    
    # get the endianness right
    better = ltoh.(better)
    data = reshape(better,(dim1, dim2))
    return data
end

read_kcd_header(loc) = begin
    seekstart(loc)
    full_text = readuntil(loc, "          ")
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
    for n in ("X dimension", "Y dimension", "Number of readouts")
        d[n] = parse(UInt64, d[n])
    end

    return d
end
