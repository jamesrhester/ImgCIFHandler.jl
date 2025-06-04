# HDF5 support with filters

using HDF5
using H5Zblosc,H5Zbzip2,H5Zlz4,H5Zzstd

# Uncomment when H5Zbitshuffle is officially registered

using H5Zbitshuffle 

"""
    imgload(handle,::Val{:HDF};path="/";frame=1)

Read a 2D HDF5 image from file `handle` located
at HDF5 location `path`, returning only frame number
`frame`, counting from 1
"""
imgload(loc::AbstractString,::Val{:HDF5};path="/",frame=1) = begin
    f = h5open(loc)
    @debug "Opening HDF5 file at" loc path frame
    if !ismissing(frame) && !isnothing(frame)
        result = f[path][:,:,frame]
    else
        result = read(f[path])
    end
    @debug "$(size(result))"
    return result
end

"""
Acceptable filters: the bundled filters of the
HDF5 filter plugin package, see 
https://portal.hdfgroup.org/display/support/HDF5+Filter+Plugins

BZIP2, LZF, BLOSC, MAFISC, LZ4, Bitshuffle, and ZFP
"""
const allowed_filters = (32000,32001,32002,32004,32008,32013)

"""
    check_format(loc::AbstractString, Val{:HDF5};path="/",frame=1)

Check that the given data set uses allowed filters
"""
check_format(loc::AbstractString, ::Val{:HDF5};path="/",frame=1) = begin
    messages = ""
    f = h5open(loc)
    ds = f[path]
    filt_pipeline = HDF5.get_create_properties(ds).filters
    for filt in filt_pipeline
        if filt isa HDF5.Filters.ExternalFilter
            filter_id = HDF5.Filters.filterid(filt)
        else
            filter_id = HDF5.Filters.filterid(typeof(filt))
        end

        # IDs 0 - 255 are maintained by HDF5 group as part of core HDF5 library
        
        if filter_id > 255 && !(filter_id in allowed_filters)
            messages *= "\nFilter ID $(filter_id) used in $loc/$path is not allowed\n"
        end
    end
    close(f)
    if messages != "" return (false,messages) end
    return (true,"")
end

