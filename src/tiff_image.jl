# HDF5 support with filters

using TiffImages

"""
    imgload(handle, ::Val{:TIFF}; frame=1)

Read from file TIFF file `handle`, returning only frame number
`frame`, counting from 1.
"""
imgload(loc::AbstractString, ::Val{:TIFF}; path = nothing, frame=1) = begin
    f = TiffImages.load(loc, lazyio = true)
    @debug "Opening TIFF file at" loc frame
    if frame >0 && frame <= size(f)[3]
        result = float.(f[:,:,frame])
    else
        result = nothing
    end
    @debug "$(size(result))"
    return result
end
