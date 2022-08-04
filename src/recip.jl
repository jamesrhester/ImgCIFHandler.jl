# Reciprocal space calculations
"""
    ewald_intersect(lambda,rot_axis,pt)

Find the points in reciprocal space where `pt`
goes through the Ewald sphere when rotated about
`rot_axis`. See gamedev.stackexchange.com/75756.

In our coordinate system the beam direction is
+ve Z. Therefore, the Ewald sphere is centred
at 0,0,1/lambda.
"""
ewald_intersect(lambda,rot,pt) = begin

    # normalise rot_axis just in case

    rot_axis = LinearAlgebra.normalize(rot)

    # centre of rotation is projection of
    # point onto the rotation axis * normal

    rot_c = dot(pt,rot_axis)*rot_axis
    rot_rad = norm(cross(pt,rot_axis))

    @debug "Point $pt lies $rot_rad from rotation centre $rot_c"
    
    ewald_c = [0,0,1.0/lambda]
    ewald_r = 1.0/lambda

    # Distance of rotation plane from sphere centre
    
    sep = dot(rot_axis,(rot_c - ewald_c))

    if abs(sep) > ewald_r
        @warn "No intersections of $pt for wavelength $lambda"
        return []
    end

    # Calculate circle of intersection plane-sphere

    plane_c = ewald_c + sep*rot_axis
    plane_r = sqrt(ewald_r^2 - sep^2)
    dist = norm(plane_c - rot_c)

    @debug "Plane of rotation $sep from Ewald sphere centre"
    @debug "Plane intersection centre $plane_c, radius $plane_r"
    
    # Now calculate intersection of the two circles

    h = 0.5 + (rot_rad^2 - plane_r^2)/(2dist^2)

    @debug "Circle intersection $h from Ewald centre"
    
    c_int = rot_c + h*(plane_c - rot_c)
    x = rot_rad^2 - h^2*dist^2
    if x < 0   # no intersections
        return nothing
    end
    r_int = sqrt(x)

    @debug "Intersection point is $r_int from Ewald centre" c_int
    
    t = LinearAlgebra.normalize(cross((plane_c - rot_c),rot_axis))
    p0 = c_int + r_int*t
    p1 = c_int - r_int*t

    @debug "Intersections" p0 p1

    return p0,p1
end

"""
    rot_angle(start,finish,cor,axis)

Calculate the angle of rotation from `start` to `finish`
around centre of rotation `cor` lying on rotation axis
`axis`.
"""
rot_angle(start,finish,axis) = begin

    # Calculate centre of rotation

    rot_c = dot(start,axis)*axis

    # Get the start and finish vector relative to centre
    
    fn = LinearAlgebra.normalize(finish - rot_c)
    sn = LinearAlgebra.normalize(start - rot_c)

    q = rotation_between(sn,fn)

    # Rotation angle will be 0 to 180
    
    ra = rotation_angle(q)*180/pi
    short_axis = rotation_axis(q)
    if dot(short_axis,axis) < 0

        # We have to go in the other direction

        ra = -1 * ra
    end
    @debug "Rot from $start to $finish around $axis is $ra"
    return ra
end

"""
    detector_intersect(ray,lambda,normal,plane_point)

Calculate the coordinates of the intersection of `ray` with a plane
described by the normal to the plane and a point on the plane. `ray`
is a point on a line stretching from (0,0,1/lambda).
"""
detector_intersect(ray,lambda,normal,plane_point) = begin

    # Intersection of ray with plane is l_0 + l*d where
    # d is given by (p_0 - l_0).n/l.n. l_0 point on line
    # p_0 point on plane. We use the real space normal
    # and origin (0,0,0) at the crystal to simplify 

    ewald_origin = [0,0,1.0/lambda]

    l = ray .- ewald_origin
    d = dot(plane_point,normal)/dot(l,normal)
    @debug "Intersection is $d from $ewald_origin"
    return d*l 
end

test_intersections() = begin

    # Some test cases

    @assert ewald_intersect(0.5,[0,1,0],[0,5,0]) == []

    ei = ewald_intersect(0.5,[0,1,0],[1,0,0])
    println(ei)
    ra = rot_angle([1,0,0],ei[1])
    ra2 = rot_angle([1,0,0],ei[2])

    println("Rotation angles for $ei are $ra, $ra2")
end
