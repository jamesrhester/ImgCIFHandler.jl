# Load in the image part of a CBF file, ignoring all else
import Base.Libc:FILE

using CBFlib_small_jll
# Comment out the above line and
# uncomment below to use your locally installed libcbf instead
# of that provided by CBFlib_small_jll

#const libcbf="libcbf"  #

# libcbf routines that we need
#
# libcbf expects the address of a cbf_handle_struct *
# typedefed to cbf_handle. The address doesn't have
# to point to an actual cbf_handle_struct but must
# be something that libcbf can update.

mutable struct CBF_Handle_Struct end

# cbflib expects us to allocate the memory for CBF Detectors
mutable struct CBF_Detector_Struct
    cbf_positioner::Ptr{Cvoid}
    displacement::NTuple{2,Float64}
    increment::NTuple{2,Float64}
    axes::Csize_t
    index::NTuple{2,Csize_t}
    cbf_handle::Ptr{Cvoid}
    element::Cint
end

mutable struct CBF_Handle
    handle::Ptr{CBF_Handle_Struct}
end

mutable struct Det_Handle
    handle::Ptr{CBF_Detector_Struct}
end

cbf_make_handle() = begin
    handle = CBF_Handle(0)
    finalizer(cbf_free_handle!,handle)
    err_no = ccall((:cbf_make_handle,libcbf),Cint,(Ref{CBF_Handle},),handle)
    cbf_error(err_no)
    #println("Our handle is $(handle.handle)")
    return handle
end

cbf_free_handle!(handle::CBF_Handle) = begin
    q = time_ns()
    #error_string = "$q: Finalizing CBF Handle $(handle.handle)"
    #t = @task println(error_string)
    #schedule(t)
    err_no = ccall((:cbf_free_handle,libcbf),Cint,(Ptr{CBF_Handle_Struct},),handle.handle)
    cbf_error(err_no, extra = "while finalising")
    return 0
end

cbf_read_file(filename) = begin
    handle = cbf_make_handle()
    f = open(filename,"r")
    fptr = Base.Libc.FILE(f)
    #println("Our file pointer is $fptr")
    err_no = ccall((:cbf_read_widefile,libcbf),Cint,(Ptr{CBF_Handle_Struct},FILE,Cint),handle.handle,fptr,0)
    cbf_error(err_no, extra = "while trying to read $filename")
    return handle
end

cbf_write_file(handle,filename) = begin
    f = open(filename,"w")
    fptr = Base.Libc.FILE(f)
    enc = 0x02 # CR at end of line
    err_no = ccall((:cbf_write_widefile,libcbf),Cint,(Ptr{CBF_Handle_Struct},FILE,Cint,Cint,Cint,Cint),handle.handle,fptr,0,1,0,enc)
    cbf_error(err_no, extra = "while writing $filename")
    close(f)
end

#
# Low-level CBF functions
#

cbf_require_category(handle::CBF_Handle,catname::AbstractString) = begin
    errno = ccall((:cbf_require_category,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Cstring),handle.handle,catname)
    cbf_error(errno,extra = "while requiring category $catname")
end

cbf_find_category(handle::CBF_Handle,catname::AbstractString) = begin
    errno = ccall((:cbf_find_category,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Cstring),handle.handle,catname)
    if cbf_error_dict[errno] == :CBF_NOTFOUND return false end
    cbf_error(errno,extra = "while finding category $catname")
    return true
end

cbf_require_column(handle::CBF_Handle,colname::AbstractString) = begin
    errno = ccall((:cbf_require_column,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Cstring),handle.handle,colname)
    cbf_error(errno,extra = "while requiring column $colname")
end

cbf_find_column(handle::CBF_Handle,colname::AbstractString) = begin
    errno = ccall((:cbf_find_column,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Cstring),handle.handle,colname)
    if cbf_error_dict[errno] == :CBF_NOTFOUND return false end
    cbf_error(errno,extra = "while finding column $colname")
    return true
end

cbf_find_row(handle::CBF_Handle,rowval::AbstractString) = begin
    errno = ccall((:cbf_find_row,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Cstring),handle.handle,rowval)
    if errno == 0 return true end
    if cbf_error_dict[errno] == :CBF_NOTFOUND return false end
    cbf_error(errno, extra = "while finding row value $rowval")
end

cbf_set_value(handle::CBF_Handle,newval) = begin
    errno = ccall((:cbf_set_value,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Cstring),handle.handle,"$newval")
    cbf_error(errno,extra = "while setting value to $newval")
end

cbf_next_row(handle::CBF_Handle) = begin
    errno = ccall((:cbf_next_row,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},),handle.handle)
    if cbf_error_dict[errno] == :CBF_NOTFOUND return false end
    cbf_error(errno, extra = "while iterating rows")
    return true
end

cbf_rewind_row(handle::CBF_Handle) = begin
    errno = ccall((:cbf_rewind_row,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},),handle.handle)
    cbf_error(errno, extra = "when rewinding to first row")
end

cbf_new_row(handle::CBF_Handle) = begin
    errno = ccall((:cbf_new_row,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},),handle.handle)
    cbf_error(errno, extra = "when making a new row")
end

cbf_count_rows(handle::CBF_Handle) = begin
    nrows = Ref{Cuint}(0)
    errno = ccall((:cbf_count_rows,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Ref{Cuint}),handle.handle,nrows)
    cbf_error(errno, extra = "when counting rows")
    return nrows[]
end

cbf_get_value(handle::CBF_Handle) = begin
    newval = Ref{Float64}(0)
    errno = ccall((:cbf_get_doublevalue,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Ref{Float64}),handle.handle,newval)
    cbf_error(errno,extra = "while getting value ")
    return newval[]
end

cbf_get_string_value(handle::CBF_Handle) = begin
    newval = Ref{Ptr{UInt8}}(0)
    errno = ccall((:cbf_get_value,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Ptr{Ptr{UInt8}}),handle.handle,newval)
    cbf_error(errno,extra = "while setting value to $newval")
    return unsafe_string(newval[])
end

cbf_require_diffrn_id(handle,newid) = begin
    current = Ref{Ptr{UInt8}}(0)
    errno = ccall((:cbf_require_diffrn_id,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Ptr{Ptr{UInt8}},Cstring),
                  handle.handle,current,newid)
    cbf_error(errno,extra = "while requiring diffrn_id")

    # Recover the string

    return unsafe_string(current[])

end

cbf_set_diffrn_id(handle,newid) = begin
    errno = ccall((:cbf_set_diffrn_id,libcbf),Cint,
                  (Ptr{CBF_Handle_Struct},Cstring),
                  handle.handle,newid)
    cbf_error(errno,extra = "while setting diffrn_id")
end

#
# Mid-level functions
#

"""
    cbf_set_axis_positions(handle::CBF_Handle,names,types,positions)

Set all of the axes in `names` with `types` of `rotation` or
`translation` to the respective `positions`, if those axes are
missing from the `diffrn_scan_frame_axis` and `diffrn_scan_axis` loops
"""
cbf_set_axis_positions(handle::CBF_Handle,names,types,positions) = begin

    @debug names,types,positions
    for (an,at,ap) in zip(names,types,positions)
        cbf_require_category(handle,"diffrn_scan_frame_axis")
        cbf_require_column(handle,"axis_id")
        if !cbf_find_row(handle,an)
            cbf_new_row(handle)
            @debug "Setting value to" an
            cbf_set_value(handle,an)
            if at=="rotation"
                cbf_require_column(handle,"angle")
                @debug "Setting angular value to " ap
                cbf_set_value(handle,"$ap")
            else
                cbf_require_column(handle,"displacement")
                @debug "Setting translation to " ap
                cbf_set_value(handle,"$ap")
            end
        end

        # And diffrn_scan_axis

        cbf_require_category(handle,"diffrn_scan_axis")
        cbf_require_column(handle,"axis_id")
        if !cbf_find_row(handle,an)
            cbf_new_row(handle)
            @debug "diffrn_scan_axis" an
            cbf_set_value(handle,an)
            if at == "rotation"
                cbf_require_column(handle,"angle_increment")
                cbf_set_value(handle,"0.0")
            else
                cbf_require_column(handle,"displacement_increment")
                cbf_set_value(handle,"0.0")
            end
        end 
    end
end

"""
    cbf_fix_detector_axes(handle,axis_names)

If detector axis names have been left out of _diffrn_detector_axis,
this will fix it, as long as there is only one diffrn.id
"""
cbf_fix_detector_axes(handle,axis_names) = begin

    # Get the name of the single detector
    
    haveit = cbf_find_category(handle,"diffrn_detector")
    if haveit && cbf_count_rows(handle) > 1 return end
    cbf_find_column(handle,"id")
    detname = cbf_get_string_value(handle)
    @debug detname
    
    # And for consistency number of axes should be true
    # Assume two axes on detector surface

    cbf_require_column(handle,"number_of_axes")
    cbf_set_value(handle,"$(length(axis_names))-2")
    
    # Add necessary rows to diffrn_detector_axis

    cbf_require_category(handle,"diffrn_detector_axis")
    for an in axis_names
        cbf_require_column(handle,"axis_id")
        if !cbf_find_row(handle,an)
            cbf_new_row(handle)
            cbf_set_value(handle,an)
            cbf_require_column(handle,"detector_id")
            cbf_set_value(handle,detname)
        end
    end
end

"""
    add_default_items(handle,itemlist)

Add all data names in itemlist to `handle`, if missing. `itemlist`
has format `[(dataname,origin_val),(,),...]` where `origin_val`
is the value to copy for dataname if present. If it has a period
character, it is a data name whose value is used, otherwise it
is the value itself.
"""
add_default_items(handle::CBF_Handle,itemlist) = begin
    for (one_item,origin_item) in itemlist
        cat,obj = split(one_item,".")
        cat = cat[2:end]

        if '.' in origin_item
            origin_cat,origin_obj = split(origin_item,".")
            origin_cat = origin_cat[2:end]
            
        # See if we have a value available

            if cbf_find_category(handle,origin_cat)
                if cbf_count_rows(handle) == 1 && cbf_find_column(handle,origin_obj)
                    def_val = cbf_get_string_value(handle)
                    @debug "Found $def_val for $origin_item"
                end
            else
                @debug "No category for $origin_item"
            end
        else
            def_val = origin_item
        end

        cbf_require_category(handle,cat)

        nrows = cbf_count_rows(handle)

        if !cbf_find_column(handle,obj)
            cbf_require_column(handle,obj)
            @debug "Adding $one_item to $nrows rows"
            if nrows == 0
                cbf_set_value(handle,def_val)
            else
                for i in 1:nrows
                    cbf_set_value(handle,def_val)
                    cbf_next_row(handle)
                end
            end
        else
            @debug "$one_item already present"
        end
    end

end

"""
    make_nice_for_cbf(handle)

Add any missing items so that cbflib doesn't complain
"""
make_nice_for_cbf(handle::CBF_Handle) = begin
    
    actual = cbf_require_diffrn_id(handle,"DEFAULT")

    # we have a special default value for _array_structure.id
 
    need_items = [("_array_structure.id","1"),
                  ("_diffrn_detector.diffrn_id","_diffrn.id"),
                  ("_diffrn_detector_element.detector_id","_diffrn_detector.id"),
                  ("_diffrn_detector_element.id","ELEMENT"),
                  ("_diffrn_data_frame.array_id","_array_structure.id"),
                  ("_diffrn_data_frame.detector_element_id","_diffrn_detector_element.id")]

    add_default_items(handle,need_items)
    
end

cbf_get_arraysize(handle) = begin
    compression = Ref{Cuint}(0)
    bid =  Ref{Cint}(0)
    elsize = Ref{Csize_t}(0)
    elsigned = Ref{Cint}(0)
    elunsigned = Ref{Cint}(0)
    elements = Ref{Csize_t}(0)
    minelem = Ref{Cint}(0)
    maxelem = Ref{Cint}(0)
    isreal = Ref{Cint}(0)
    byteorder = Ref{Ptr{UInt8}}(0)
    fast = Ref{Csize_t}(0)
    mid = Ref{Csize_t}(0)
    slow = Ref{Csize_t}(0)
    padding = Ref{Csize_t}(0)
    err_no = ccall((:cbf_get_arrayparameters_wdims,libcbf),Cint,
              (Ptr{CBF_Handle_Struct},
               Ref{Cuint}, #compression
               Ref{Cint}, #binary_id
               Ref{Csize_t}, #elsize
               Ref{Cint},
               Ref{Cint},
               Ref{Csize_t}, #elements
               Ref{Cint},
               Ref{Cint},
               Ref{Cint}, #is real
               Ref{Ptr{UInt8}}, #byteorder
               Ref{Csize_t}, #fast
               Ref{Csize_t}, #mid
               Ref{Csize_t}, #slow
               Ref{Csize_t}, #padding
               )
              , handle.handle,compression, bid, elsize,elsigned,elunsigned,
              elements,minelem, maxelem, isreal, byteorder,fast,mid,slow,padding)
    cbf_error(err_no, extra = "while trying to get array size")
    #println("Compression type $(compression[]) for binary id $(bid[])")
    #println("elements = $(elements[]) of size $(elsize[])")
    #println("Dims $(fast[]) x $(mid[]) x $(slow[])")
    #println("Is real? $(isreal[])")
    bo = @GC.preserve unsafe_string(byteorder[])
    #println("Check: byte order is $bo")
    if isreal[] != 0    # real numbers
        if elsize[] == 4
            elt = Float32
        elseif elsize[] == 8
            elt = Float64
        else
            throw(error("No real type with size $(elsize[])"))
        end
    else   # integers
        if elsize[] == 4
            elt = Int32
        elseif elsize[] == 8
            elt = Int64
        elseif elsize[] == 2
            elt = Int16
        else throw(error("No integer type with size $(elsize[])"))
        end
    end
    return elements[],elt,fast[],mid[],slow[]
end

cbf_get_realarray(handle,data_array) = begin
    bid = Ref{Cint}(0)
    num_read = Ref{Csize_t}(0)
    fast,mid = size(data_array)
    elsize = sizeof(eltype(data_array))
    err_no = ccall((:cbf_get_realarray,libcbf),Cint,
                   (Ptr{CBF_Handle_Struct},
                    Ref{Cint}, #binary id
                    Ptr{Float64},
                    Csize_t,
                    Csize_t,
                    Ref{Csize_t}
                    ),handle.handle,bid,data_array,elsize,fast*mid,num_read)
    cbf_error(err_no, extra = "while reading in array")
    #println("Read in $(num_read[]) values")
    return num_read[]
end

cbf_get_integerarray(handle,data_array) = begin
    bid = Ref{Cint}(0)
    num_read = Ref{Csize_t}(0)
    fast,mid = size(data_array)
    elt = eltype(data_array)
    elsize = sizeof(elt)
    is_signed = Ref{Cint}(signed(elt) == elt ? 1 : 0)
    err_no = ccall((:cbf_get_integerarray,libcbf),Cint,
                   (Ptr{CBF_Handle_Struct},
                    Ref{Cint}, #binary id
                    Ptr{Float64},
                    Csize_t,
                    Ref{Cint},
                    Csize_t,
                    Ref{Csize_t}
                    ),handle.handle,bid,data_array,elsize,is_signed,fast*mid,num_read)
    cbf_error(err_no, extra = "while reading in array")
    #println("Read in $(num_read[]) values")
    return num_read[]
end

#
#  Extra routines just because they are useful
#

cbf_get_beam_center(handle) = begin
    index1 = Ref{Float64}(0.0)
    index2 = Ref{Float64}(0.0)
    centre1 = Ref{Float64}(0.0)
    centre2 = Ref{Float64}(0.0)
    @debug "Det info is: $(handle.handle)"
    err_no = ccall((:cbf_get_beam_center,libcbf),Cint,
                   (Ptr{CBF_Detector_Struct},
                    Ref{Float64},
                    Ref{Float64},
                    Ref{Float64},
                    Ref{Float64}),
                   handle.handle,index1,index2,centre1,centre2)
    cbf_error(err_no, extra = "while getting beam centre")
    return centre1[],centre2[],index1[],index2[]
end

cbf_get_pixel_coordinates(handle,slowcoord,fastcoord) = begin
    coordx = Ref{Float64}(0.0)
    coordy = Ref{Float64}(0.0)
    coordz = Ref{Float64}(0.0)
    err_no = ccall((:cbf_get_pixel_coordinates,libcbf),Cint,
                   (Ptr{CBF_Detector_Struct},
                   Cdouble,
                   Cdouble,
                   Ref{Cdouble},
                   Ref{Cdouble},
                   Ref{Cdouble}),
                   handle.handle,slowcoord,fastcoord,coordx,coordy,coordz)
    cbf_error(err_no, extra = "while getting pixel coordinate $slowcoord,$fastcoord")
    return coordx[],coordy[],coordz[]
end

"""
    imgload(filename,::Val{:CBF})

Return a single image from `filename`, which should be in (mini)CBF format and contain a single
frame.
"""
imgload(filename::AbstractString,::Val{:CBF};path=nothing,frame=nothing) = begin
    handle = cbf_read_file(filename)
    if !cbf_find_category(handle,"array_data")
        throw(error("category array_data not present in $filename"))
    end
    if !cbf_find_column(handle,"data")
        throw(error("item _array_data.data missing in $filename"))
    end
    err_no = ccall((:cbf_rewind_row,libcbf),Cint,(Ptr{CBF_Handle_Struct},),handle.handle)
    cbf_error(err_no)
    # Find the dimensions of the array
    total_size,elt,fast,mid,slow = cbf_get_arraysize(handle)
    # Create an array for this information
    data_array = Array{elt,2}(undef,fast,mid)
    if elt <: AbstractFloat
        num_read = cbf_get_realarray(handle,data_array)
    else
        num_read = cbf_get_integerarray(handle,data_array)
    end
    if num_read != total_size
        throw(error("Read failure: expected $total_size, read $num_read"))
    end
    return data_array
end

"""
    prepare_detector(filename,scan,frame)

Return a CBF detector object (`det_handle`) and detector information `det_data`
 ready for calculations positioned according
to `frame` of `scan`. `scan` and `frame` may be omitted in which case the
zero settings are used.

 `destruct_detector` must be called when the detector is no longer used.
The `det_data` item returned must be preserved from garbage collection until
the returned handle is no longer needed.
"""
prepare_detector(filename,args...) = begin
    # determine axis settings for scan

    axis_names,types,postns = get_detector_axis_settings(filename,args...)

    handle = cbf_read_file(filename)

    cbf_set_axis_positions(handle,axis_names,types,postns)
    cbf_fix_detector_axes(handle,axis_names)

    make_nice_for_cbf(handle)

    @debug "Set axes to" axis_names types postns
    
    # now go and find them

    det_data = CBF_Detector_Struct(Ptr{Cvoid}(),
                                     (0.0,0.0),
                                     (0.0,0.0),
                                     0,(0,0),handle.handle,
                                   0)
    det_handle = Det_Handle(pointer_from_objref(det_data))
    @debug "Det handle: $det_handle, $(pointer_from_objref(det_handle))"
    @debug "det data: $(pointer_from_objref(det_data))"
    err_no = ccall((:cbf_construct_detector,libcbf),Cint,
                   (Ptr{CBF_Handle_Struct},
                    Ref{Det_Handle},
                    Cuint),
                   handle.handle,det_handle,0)
    cbf_error(err_no, extra = "while constructing detector")
    @debug "Det info is: $det_data"
    return det_handle,det_data
end

destruct_detector(det_handle) = begin
    err_no = ccall((:cbf_free_detector,libcbf),Cint,
                                   (Ptr{CBF_Detector_Struct},),det_handle.handle)
    cbf_error(err_no, extra = "while destructing detector")
end

"""
    get_beam_centre(filename,scan,frame)

Return the beam centre for `frame` in `scan`, taking into account all
axis positions. We return `det_data` as must be preserved until the
handle is destroyed.
"""
get_beam_centre(filename::AbstractString,args...) = begin

    det_handle, det_data = prepare_detector(filename,args...)    
    centre1,centre2,index1,index2 = cbf_get_beam_center(det_handle)
    GC.@preserve det_data destruct_detector(det_handle)
    
    return centre1,centre2,index1,index2

end

get_beam_centre(incif::CifContainer,args...) = begin
    get_beam_centre("$(incif.original_file)",args...)
end

"""
   get_pixel_coordinates(incif::CifContainer,fast_coord,slow_coord,scan,frame)

Return the coordinates in lab space of the nominated pixel, with axes oriented
as for `frame` in `scan`. `scan` and `frame` are optional.
"""
get_pixel_coordinates(filename::AbstractString,slow_coord,fast_coord,args...) = begin
    dh, dd = prepare_detector(filename,args...)
    x,y,z = cbf_get_pixel_coordinates(dh,slow_coord,fast_coord)
    GC.@preserve dd destruct_detector(dh)
    return [x,y,z]
end

get_pixel_coordinates(incif::CifContainer,args...) = begin
    filename = incif.original_file
    get_pixel_coordinates(filename,args...)
end

cbf_error(val;extra="") = begin
    if val == 0 return end
    throw(error("CBF Error: $(cbf_error_dict[val]) $extra"))
    # TODO: actually do an AND for multiple errors
end

const cbf_error_dict = Dict(
          0 => :CBF_SUCCESS        ,
 0x00000001 => :CBF_FORMAT         , 
 0x00000002 => :CBF_ALLOC          , 
 0x00000004 => :CBF_ARGUMENT       , 
 0x00000008 => :CBF_ASCII          , 
 0x00000010 => :CBF_BINARY         , 
 0x00000020 => :CBF_BITCOUNT        ,
 0x00000040 => :CBF_ENDOFDATA       ,
 0x00000080 => :CBF_FILECLOSE       ,
 0x00000100 => :CBF_FILEOPEN        ,
 0x00000200 => :CBF_FILEREAD        ,
 0x00000400 => :CBF_FILESEEK        ,
 0x00000800 => :CBF_FILETELL        ,
 0x00001000 => :CBF_FILEWRITE       ,
 0x00002000 => :CBF_IDENTICAL       ,
 0x00004000 => :CBF_NOTFOUND        ,
 0x00008000 => :CBF_OVERFLOW        ,
 0x00010000 => :CBF_UNDEFINED       ,
 0x00020000 => :CBF_NOTIMPLEMENTED  ,
 0x00040000 => :CBF_NOCOMPRESSION   ,  
 0x00080000 => :CBF_H5ERROR         ,
 0x00100000 => :CBF_H5DIFFERENT     ,
 0x00200000 => :CBF_SIZE            )
