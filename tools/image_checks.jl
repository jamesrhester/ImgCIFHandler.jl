# Check that the image is described correctly
# We assume a single array structure for all data
const img_types = Dict(UInt8 =>"unsigned 8-bit integer",
                       UInt16=>"unsigned 16-bit integer",
                       UInt32=>"unsigned 32-bit integer",
                       Int8  =>"signed 8-bit integer",
                       Int16 =>"signed 16-bit integer",
                       Int32 =>"signed 32-bit integer",
                       Float32 =>"signed 32-bit real IEEE",
                       Float64 =>"signed 64-bit real IEEE",
                       ComplexF32 =>"signed 32-bit complex IEEE"
                       )

@imgcheck "Image type and dimensions" img_dims(incif,img,img_id) = begin
    messages = []
    fast_pixels = size(img)[1]
    slow_pixels = size(img)[2]
    fast_pos,slow_pos = indexin(["1","2"],incif["_array_structure_list.precedence"])
    dims = parse.(Int32,incif["_array_structure_list.dimension"])
    if dims[fast_pos] != fast_pixels
        push!(messages,(false,"Stated fast dimension $fast_pos does not match actual image dimension $fast_pixels"))
    end
    if dims[slow_pos] != slow_pixels
        push!(messages,(false,"Stated slow dimension $slow_pos does not match actual image dimension $slow_pixels"))
    end
    if haskey(incif,"_array_structure.encoding_type") && img_types[eltype(img)] != incif["_array_structure.encoding_type"][1]
        push!(messages,(false,"Stated encoding $(incif["_array_structure.encoding_type"][1]) does not match array element type $(eltype(img))"))
    end
    if haskey(incif,"_array_structure.byte_order") && haskey(incif,"_array_data.external_data_id")
        push!(messages,(true,"WARNING: byte order provided in file containing external data pointers"))
    end
    if haskey(incif,"_array_structure.compression") && incif["_array_structure.compression"][1] != "none" && haskey(incif,"_array_data.external_data_id")
        push!(messages,(true,"Externally-provided data by definition is uncompressed but compression is specified as $(incif["array_structure.compression"][1])"))
    end
    return messages
end

@imgcheck "Overloaded values present" mask_check(incif,img,img_id) = begin
    messages = []

    # Assume that more values of typemax than of preceding 10 values is likely a
    # sign of overload

    maxval = typemax(eltype(img))

    num_max = count(x->x == maxval,img)
    num_lower = count(x-> x > maxval - 10 && x < maxval,img)
    if num_max > num_lower
        if !haskey(incif,"_array_intensities.overload")
            push!(messages,(false,"$num_max apparently overloaded intensities present (values of $maxval), but _array_intensities.overload is missing"))
        end
    end
    if haskey(incif,"_array_intensities.overload")
        ovld = parse(Float64,incif["_array_intensities.overload"][])
        if ovld > maxval
            push!(messages(false,"_array_intensities.overload value of $ovld is greater than maximum possible value for this image element type $maxval"))
        end
    end
    return messages
end
