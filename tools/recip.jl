# Reciprocal space calculations
using LinearAlgebra

"""
    ewald_intersect(lambda,rot_axis,pt)

Find the points in reciprocal space where `pt`
goes through the Ewald sphere when rotated about
`rot_axis`. See gamedev.stackexchange.com/75756.
"""
ewald_intersect(lambda,rot_axis,pt) = begin

    # normalise rot_axis just in case

    normalize!(rot_axis)

    # centre of rotation is projection of
    # point onto the rotation axis * normal

    rot_c = dot(pt,rot_axis)*rot_axis
    rot_rad = norm(cross(pt,rot_axis))

    @debug "Point $pt lies $rot_rad from rotation centre $rot_c"
    
    ewald_c = [0,0,-1.0/lambda]
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
    r_int = sqrt(rot_rad^2 - h^2*dist^2)

    @debug "Intersection point is $r_int from Ewald centre" c_int
    
    t = normalize(cross((plane_c - rot_c),rot_axis))
    p0 = c_int + r_int*t
    p1 = c_int - r_int*t

    @debug "Intersections" p0 p1

    return p0,p1
end

"""
    rot_angle(start,finish,axis)

Calculate the angle of rotation from `start` to `finish`
around `axis`, which is perpendicular to `start` and `finish`.
"""
rot_angle(start,finish,axis) = begin
    ra = acosd(dot(normalize!(finish),normalize!(start)))
    cr = cross(start,finish)
    if dot(cr,axis) < 0 ra = -1*ra end
    return ra,cr
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
