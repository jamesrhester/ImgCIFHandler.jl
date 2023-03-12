"""
    imgload(handle,::Val{:SMV})

Read an image from an SMV-formatted file. Only 2-byte integers currently
recognised, should be expanded to cover all known options.
"""
imgload(filename::AbstractString,::Val{:SMV};path=nothing,frame=nothing) = begin
    loc = open(filename,"r")
    header = read_adsc_header(loc)
    dim1 = header["size1"]
    dim2 = header["size2"]
    seek(loc,header["header_bytes"])
    binary = read(loc,2*dim1*dim2)  #Sequence of UInt8
    if haskey(header,"bitmapsize")
        println("Bitmap present, ignored")
        #bmap = read_bitmap(loc,header)
        #data = apply_bitmap(better,bmap)
    end
    better = reinterpret(UInt16,binary)
    # get the endianness right
    if header["byte_order"]=="little_endian" 
        better = ltoh.(better)
    else
        better = ntoh.(better)
    end
    data = reshape(better,(dim1,dim2))
    return data
end

read_adsc_header(loc) = begin
    seekstart(loc)
    line = readline(loc)
    header = Dict{String,Any}()
    while !occursin("}",line)
        if occursin("=",line)
            (key,val) = split(String(line),"=")
            header[lowercase(strip(key))] = strip(val,[' ',';','\n','\r']) #an array
        end
        line = readline(loc)
    end
    # parse things that we need
    header["size1"] = parse(UInt32,header["size1"])
    header["size2"] = parse(UInt32,header["size2"])
    header["header_bytes"] = parse(UInt32,header["header_bytes"])
    if haskey(header,"bitmapsize")
        header["bitmapsize"] = parse(UInt32,header["bitmapsize"])
    end
    @debug "Header, image info" header["header_bytes"] (header["size1"]*header["size2"]*2)
    return header
end

# Read any bitmap in io. The first characters are `BRLE`.
read_bitmap(io,header) = begin
    bmap_len = header["bitmapsize"]
    image_len = header["size1"]*header["size2"]*2
    seek(io,image_len + header["header_bytes"])
    marker = read(io,4)
    # Should be 'BRLE'
    marker = convert.(Char,marker)
    if !(marker[1] == 'B' && marker[2] == 'R' && marker[3] == 'L' && marker[4] == 'E')
        throw(error("RLE Bitmap marker not found: found $marker"))
    end
    raw = read(io,bmap_len-4)
    return reinterpret(UInt16,raw)
end

# Bitmap is a run-length encoded bitmap. It describes a mask, where the
# value for the mask (1/0) is in the MSB and the remainder of the integer is
# the number of entries containing that value.
apply_bitmap(raw,bmap) = begin
    pos = 1
    while pos < length(raw)
        for one_section in raw
        end
    end
end

# Recreate the bitmap
expand_bitmap(bmap,dim1,dim2) = begin
    tbf = BitVector()
    for one_section in 1:length(bmap)
        rlen = bmap[one_section] & UInt32(32767)
        val =  bmap[one_section] & UInt32(32768) != 0
        append!(tbf,fill(val,rlen))
        println("$one_section:$(length(tbf)) after adding $rlen x $val")
    end
    return tbf
end

            
