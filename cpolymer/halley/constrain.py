# constraints - a geometrical constraints solver

'''
http://www.halley.cc/code/python/constraints.py

A geometrical constraints solver for finding intersection and proximity.

ABSTRACT

A constraint is a description of a geometric rule, such as what positions
the tip of a clock's hand may be found, or where the nearest point on a
plane may be found.

Constraints are like mathematical sets in that the combination or
intersection (multiplication) of two constraints returns a simpler
constraint.  For example, if two infinite lines are not co-linear, then
they will intersect at one point (a Unity) or at no points (Nowhere).

Supported constraint types currently include Nowhere, Unity, Duality,
Linear, Planar, Circular, Spherical, and Anywhere.

SYNOPSIS

    >> import vectors ; from vectors import *
    >> import constraints ; from constraints import *

    >> c = Circular(point=V.O, normal=V.X, radius=1)
    >> l = Linear(pointA=V.O, pointB=V.Y)

    >> c * l
       Duality( V(0.0, 1.0, 0.0), V(0.0, -1.0, 0.0) )

    >> c.closest( V(0.0, 1.0, 1.0) )
       V(0.0, 0.7071067811, 0.7071067811)

    >> c.accepts( V(0.0, -1.0, 0.0) )
       True

AUTHOR

    Ed Halley (ed@halley.cc) 15 May 2005

'''

__all__ = [ 'Anywhere', 'Nowhere',
            'Unity', 'Duality',
            'Linear', 'Circular', 'Spherical' ]

from . import vectors ; from .vectors import V, zero, equal, EPSILON
import math ; from math import sqrt
import random
#----------------------------------------------------------------------------

class Constraint:
    '''The base class of all simple geometrical constraints.'''

    def __init__(self): pass
    def combine(self, other): raise NotImplemented
    def closest(self, point): raise NotImplemented
    def accepts(self, point): raise NotImplemented
    def __mul__(self, other): return self.combine(other)
    def __rmul__(self, other): return other.combine(self)
    def __str__(self): return self.__repr__()
    def __bool__(self): return True

#----------------------------------------------------------------------------

class Anywhere (Constraint):
    '''If a point is completely unconstrained, it can be Anywhere.'''

    def __repr__(self): return "Anywhere()"

    def combine(self, other):
        "Anywhere combined with any other constraint yields the other."
        return other

    def closest(self, point):
        "The closest point inside Anywhere from a given point is that point."
        return point

    def accepts(self, point):
        "Any point is accepted."
        return True

#----------------------------------------------------------------------------

class Nowhere (Constraint):
    "If constraints on a point are inconsistent, it can be Nowhere."

    def __repr__(self): return "Nowhere()"
    def __bool__(self): return False

    def combine(self, other):
        "Nowhere combined with any other constraint yields Nowhere."
        return self

    def closest(self, point):
        "There is no closest point to Nowhere, but we return the original."
        return point

    def accepts(self, point):
        "No point is accepted."
        return False

#----------------------------------------------------------------------------

class Unity (Constraint):
    "If a point can only be in a single place, this constraint is a Unity."

    def __init__(self, point):
        Constraint.__init__(self)
        self.point = V(point)

    def __repr__(self): return "Unity(" + repr(self.point) + ")"

    def combine(self, other):
        "Unity combined with any other constraint returns Unity or Nowhere."
        if other.accepts(self.point): return self
        return Nowhere()

    def closest(self, point):
        "The closest point within a Unity is the one described."
        return self.point

    def accepts(self, point):
        "One point is accepted."
        return self.point == point

#----------------------------------------------------------------------------

class Duality (Constraint):
    "If a point can only be in one of two places, this is a Duality."

    def __init__(self, pointA, pointB):
        Constraint.__init__(self)
        if pointA == pointB:
            raise ValueError('Duality must refer to two different points.')
        self.pointA = V(pointA)
        self.pointB = V(pointB)

    def __repr__(self):
        return "Duality" + repr((self.pointA, self.pointB))

    def combine(self, other):
        "Duality versus another constraint is Duality, Unity or Nowhere."
        a = other.accepts(self.pointA)
        b = other.accepts(self.pointB)
        if a and b: return self
        if a: return Unity(self.pointA)
        if b: return Unity(self.pointB)
        return Nowhere()

    def closest(self, point):
        "The closer of the two points is returned."
        dA = self.pointA.dsquared(point)
        dB = self.pointA.dsquared(point)
        if dA <= dB: return self.pointA
        return self.pointB

    def accepts(self, point):
        "Either of our points is accepted."
        if self.pointA == point: return True
        if self.pointB == point: return True
        return False

#----------------------------------------------------------------------------

def is_point_on_line(la, dir, p):
    if not zero(dir[0]):
        alpha = (p[0] - la[0]) / dir[0]
        y = dir[1] * alpha + la[1]
        z = dir[2] * alpha + la[2]
        return equal(y, p[1]) and equal(z, p[2])
    elif not zero(dir[1]):
        alpha = (p[1] - la[1]) / dir[1]
        z = dir[2] * alpha + la[2]
        return equal(la[0], p[0]) and equal(z, p[2])
    elif not zero(dir[2]): # else:
        return equal(la[0], p[0]) and equal(la[1], p[1])
    return False

def closest_point_on_line(la, dir, p):
    off = p - la
    det = off.dot(dir) / dir.dot(dir)
    return la + dir*det

class Linear (Constraint):
    "A Linear constraint allows a point anywhere along an infinite line."
    
    def __init__(self, pointA, pointB):
        Constraint.__init__(self)
        if pointA == pointB:
            raise ValueError('A Linear must refer to two different points.')
        self.pointA = V(pointA)
        self.pointB = V(pointB)
        self.dir = self.pointB - self.pointA

    def __repr__(self):
        return "Linear" + repr((self.pointA, self.pointB))

    def combine(self, other):
        "Linear versus something else is Linear, Duality, Unity or Nowhere."
        if isinstance(other, Anywhere): return self
        if isinstance(other, Nowhere): return other
        if isinstance(other, Unity): return other.combine(self)
        if isinstance(other, Duality): return other.combine(self)
        if isinstance(other, Linear): return self._linear_linear(other)
        if isinstance(other, Planar): return other.combine(self)
        if isinstance(other, Circular): return other.combine(self)
        if isinstance(other, Spherical): return other.combine(self)
        return Nowhere()

    def _linear_linear(self, other):
        # linear vs linear
        # are we colinear (self) or parallel (Nowhere)?
        dirS = self.pointB - self.pointA
        dirO = other.pointB - other.pointA
        crossSO = dirO * dirS
        if crossSO.zero():
            if self.accepts(other.pointA): return self
            return Nowhere()
        # make an arbitrary line joining our two lines; get their crosses
        # if the crosses are parallel, the three form a coplanar triangle
        # if the crosses are not parallel, our lines don't meet (Nowhere)
        span = self.pointA - other.pointA
        crossO = span * dirO
        if crossO.zero(): return Unity(self.pointA)
        crossS = span * dirS
        if crossS.zero(): return Unity(other.pointA)
        crossXSO = crossO * crossS
        if not crossXSO.zero(): return Nowhere()
        # find the intersection (Unity)
        dist = 0.0
        if not zero(crossSO[2]):
            dist = crossO[2] / crossSO[2]
        elif not zero(crossAB[1]):
            dist = crossO[1] / crossSO[1]
        elif not zero(crossAB[0]): # else:
            dist = crossO[0] / crossSO[0]
        return Unity(self.pointA + (dirS*dist))

    def closest(self, point):
        "A point on the line is returned."
        return closest_point_on_line(self.pointA, self.dir, point)

    def accepts(self, point):
        "True if the point is on the line."
        return is_point_on_line(self.pointA, self.dir, point)

def _linear_linear_test():
    from testing import __ok__
    # simple intersect at origin, commutative checks
    u = V(-1, 0, 0) ; v = -u ; c1 = Linear(u, v)
    r = V(0, -1, 0) ; s = -r ; c2 = Linear(r, s)
    __ok__(str(c1 * c2) == "Unity(V(0.0, 0.0, 0.0))")
    __ok__(str(c2 * c1) == "Unity(V(0.0, 0.0, 0.0))")
    c1 = Linear(v, u)
    __ok__(str(c1 * c2) == "Unity(V(0.0, 0.0, 0.0))")
    __ok__(str(c2 * c1) == "Unity(V(0.0, 0.0, 0.0))")
    # intersect away from origin
    o = V(0.125, 0.25, 0.5)
    u = V(-2, 3, -4) ; v = -u ; u += o ; v += o ; c1 = Linear(u, v)
    r = V(-3, 8, 21) ; s = -r ; r += o ; s += o ; c2 = Linear(r, s)
    __ok__(str(c2 * c1) == str(Unity(o)))
    __ok__(str(c1 * c2) == str(Unity(o)))
    # non-intersecting
    u = V(-1, 0, 0) ; v = -u
    u += V(0.0, 0.0, 0.5)
    v += V(0.0, 0.0, 0.5)
    c1 = Linear(u, v)
    r = V(0, -1, 0) ; s = -r
    c2 = Linear(r, s)
    __ok__(str(c1 * c2) == "Nowhere()")
    # intersecting
    r += 0.5
    s += 0.5
    c2 = Linear(r, s)
    __ok__(str(c1 * c2) == "Unity(V(0.5, 0.0, 0.5))")
    # colinear
    u = s + V(0, 0.5, 0)
    v = r + V(0, 0.5, 0)
    c1 = Linear(u, v)
    __ok__(str(c1 * c2) == str(c1))

#----------------------------------------------------------------------------

def distance_from_plane(ppoint, normal, p):
    # negative if below the plane
    v = p - ppoint
    return v.dot(normal)

def is_point_on_plane(ppoint, normal, p):
    return zero(distance_from_plane(ppoint, normal, p))

def closest_point_on_plane(ppoint, normal, p):
    dist = distance_from_plane(ppoint, normal, p)
    return p - normal*dist

def intersect_plane_and_line(ppoint, normal, p, dir):
    # only handles actual intersecting cases (dot != 0)
    dot = normal.dot(dir)
    dots = normal.dot(ppoint) - normal.dot(p)
    alpha = dots / dot
    intersect = p + dir*alpha
    return intersect

class Planar (Constraint):
    "A Planar constraint allows a point anywhere on an infinite plane."
    
    def __init__(self, point, normal):
        Constraint.__init__(self)
        self.point = V(point)
        self.normal = V(normal).normalize()
        self.d = self.normal.dot(self.point)
        #d = -distance_from_plane(self.point, self.normal, V())
        #__ok__(equal(d, self.d))

    def __repr__(self):
        return "Planar" + repr((self.point, self.normal))

    def combine(self, other):
        "Planar versus another constraint."
        if isinstance(other, Anywhere): return self
        if isinstance(other, Nowhere): return other
        if isinstance(other, Unity): return other.combine(self)
        if isinstance(other, Duality): return other.combine(self)
        if isinstance(other, Linear): return self._planar_linear(other)
        if isinstance(other, Planar): return self._planar_planar(other)
        if isinstance(other, Circular): return other.combine(self)
        if isinstance(other, Spherical): return other.combine(self)
        return Nowhere()

    def _planar_linear(self, other):
        dir = other.dir
        dot = self.normal.dot(dir)
        # are we coplanar (Linear) or parallel (Nowhere)?
        if zero(dot):
            if self.accepts(other.pointA): return other
            return Nowhere()
        # find the intersection (Unity)
        intersect = intersect_plane_and_line(self.point, self.normal,
                                             other.pointA, other.dir)
        return Unity(intersect)

    def _planar_planar(self, other):
        crossSO = self.normal * other.normal
        # are we coplanar (Planar) or parallel (Nowhere)?
        if crossSO.zero():
            if self.accepts(other.point): return self
            return Nowhere()
        # find the intersection (Linear)
        dir = crossSO
        point = V()
        if not zero(crossSO[2]):
            point[0] = ( other.normal[1] * self.d -
                         self.normal[1] * other.d ) / crossSO[2]
            point[1] = ( self.normal[0] * other.d -
                         other.normal[0] * self.d ) / crossSO[2]
            point[2] = 0.0
        elif not zero(crossSO[1]):
            point[0] = ( self.normal[2] * other.d -
                         other.normal[2] * self.d ) / crossSO[1]
            point[1] = 0.0
            point[2] = ( other.normal[0] * self.d -
                         self.normal[0] * other.d ) / crossSO[1]
        elif not zero(crossSO[0]): # else:
            point[0] = 0.0
            point[1] = ( other.normal[1] * self.d -
                         self.normal[1] * other.d ) / crossSO[0]
            point[2] = ( self.normal[0] * other.d -
                         other.normal[0] * self.d ) / crossSO[0]
        return Linear(point, point + dir)

    def closest(self, point):
        "The closer of the two points is returned."
        return closest_point_on_plane(self.pointA, self.pointB, point)

    def accepts(self, point):
        "True if the point is on the plane."
        return is_point_on_plane(self.point, self.normal, point)

def _planar_linear_test():
    from testing import __ok__
    p = Planar(V(3,0,-1), V(0,3,0))
    __ok__(str(p) == "Planar(V(3.0, 0.0, -1.0), V(0.0, 1.0, 0.0))")
    l = Linear(V(2,1,1), V(0,-1,1))
    u = p*l
    __ok__(str(u) == "Unity(V(1.0, 0.0, 1.0))")
    p = Planar(V(3,-7,-1), V(1,1,0))
    __ok__(equal(sqrt(2)/2, p.normal[1]))

def _planar_planar_test():
    from testing import __ok__
    p = Planar(V(3,0,-1), V(0,3,0))
    q = Planar(V(5,0,-5), V(1,1,0))
    pq = p*q
    __ok__(isinstance(pq, Linear))
    __ok__(equal(pq.pointB[1], pq.pointA[1]))
    q = Planar(V(5,0,-5), V(0,1,0))
    pq = p*q
    __ok__(pq is p)
    q = Planar(V(5,1,-5), V(0,1,0))
    pq = p*q
    __ok__(isinstance(pq, Nowhere))

#----------------------------------------------------------------------------

def is_point_on_circle(center, normal, radius, p):
    return (is_point_on_plane(center, normal, p) and
            is_point_on_sphere(center, radius, p))

def closest_point_on_circle(center, normal, radius, p):
    planar = closest_point_on_plane(center, normal, p)
    spoke = planar - center
    if not spoke: spoke = V(1, 0, 0) * normal
    if not spoke: spoke = V(0, 1, 0) * normal
    spoke = spoke.magnitude(radius)
    return center + spoke

class Circular (Constraint):
    "A Circular constraint allows a point anywhere along a circle's edge."
    
    def __init__(self, point, normal, radius):
        Constraint.__init__(self)
        if zero(radius):
            raise ValueError('A circle must have a non-zero radius.')
        self.point = V(point)
        self.normal = V(normal).normalize()
        self.radius = abs(float(radius))
        self.d = self.normal.dot(self.point)

    def __repr__(self):
        return "Circular" + repr((self.point, self.normal, self.radius))
    
    def get_random(self):
        if self.normal[0] != 0 or self.normal[1] != 0:
            x0 = V(-self.normal[1],self.normal[0],0)
        else:
            x0 = V(-self.normal[2],0,self.normal[0])
        x0 = x0.normalize()
        self.normal = self.normal.normalize()
        x1 = self.normal.cross(x0)
        alpha = 2*math.pi*random.random()
        return self.point + self.radius * ( x0 * math.cos(alpha) + x1 * math.sin(alpha))        
    def generate_N(self,N):
        if self.normal[0] != 0 or self.normal[1] != 0:
            x0 = V(-self.normal[1],self.normal[0],0)
        else:
            x0 = V(-self.normal[2],0,self.normal[0])
        x0 = x0.normalize()
        self.normal = self.normal.normalize()
        x1 = self.normal.cross(x0)
        alphas = [ (2*3.1415*i)/N for i in range(N) ]
        return [self.point + self.radius * ( x0 * math.cos(alpha) + x1 * math.sin(alpha)) for alpha in  alphas]
    def combine(self, other):
        "Circular versus another constraint gives Unity, Duality or Nowhere."
        if isinstance(other, Anywhere): return self
        if isinstance(other, Nowhere): return other
        if isinstance(other, Unity): return other.combine(self)
        if isinstance(other, Duality): return other.combine(self)
        if isinstance(other, Linear): return self._circular_linear(other)
        if isinstance(other, Planar): return self._circular_planar(other)
        if isinstance(other, Circular): return self._circular_circular(other)
        if isinstance(other, Spherical): return other.combine(self)
        return Nowhere()

    def _circular_linear(self, other):
        dir = other.pointB - other.pointA
        dot = self.normal.dot(dir)
        # are we coplanar (Linear) or parallel (Nowhere)?
        if zero(dot):
            return self._circular_coplanar_linear(other)
        # find the intersection (Unity)
        intersect = intersect_plane_and_line(self.point, self.normal,
                                             other.pointA, other.dir)
        if self.accepts(intersect):
            return Unity(intersect)
        return Nowhere()

    def _circular_coplanar_linear(self, other):
        # might find unity, might find duality, might find nowhere
        mark = other.closest(self.point)
        aa = mark.dsquared(self.point)
        cc = self.radius * self.radius
        if aa > cc + EPSILON:
            return Nowhere()
        if aa < cc - EPSILON:
            chord = V(other.dir)
            b = sqrt(cc - aa)
            chord = chord.magnitude(b)
            return Duality(mark + chord, mark - chord)
        return Unity(mark)

    def _circular_planar(self, other):
        # first planar planar, which may provide a coplanar line
        planar = Planar(self.point, self.normal)
        intersect = planar.combine(other)
        if isinstance(intersect, Nowhere): return intersect
        if isinstance(intersect, Planar): return self
        return self._circular_coplanar_linear(intersect)

    def _circular_circular(self, other):
        # decide between coplanar and offplanar solutions
        planarS = Planar(self.point, self.normal)
        planarO = Planar(other.point, other.normal)
        intersect = planarS.combine(planarO)
        if isinstance(intersect, Nowhere): return intersect
        if isinstance(intersect, Planar):
            return self._circular_coplanar_circular(other)
        # the line intersection may cut through zero, one or both circles
        reducedS = self.combine(intersect)
        reducedO = other.combine(intersect)
        return reducedS.combine(reducedO)

    def _circular_coplanar_circular(self, other):
        cc = other.point - self.point
        d = cc.magnitude()
        rS = self.radius
        rO = other.radius
        # too far apart?
        if (rS + rO) < d - EPSILON:
            return Nowhere()
        # equal or circumscribed?
        if zero(d):
            if equal(rS, rO):
                return self
            return Nowhere()
        # tangent?
        if equal(d, rS+rO):
            cc = cc.magnitude(rS)
            return Unity(self.point + cc)
        # find the midpoint where the "radical line" (chord) crosses cc
        term = d*d - rS*rS + rO*rO
        xS = term / 2*d
        xO = d - xS
        # find the length of the radical line
        leg = 4*d*d*rS*rS - term*term
        yy = 1/d * sqrt(leg)
        y = yy/2
        # find the intersection points
        cross = cc * self.normal
        chord = cc * cross
        chord = chord.magnitude(y)
        mark = cc
        mark = mark.magnitude(xS)
        mark += self.point
        return Duality(mark + chord, mark - chord)

    def closest(self, point):
        "A point on the circle is returned."
        return closest_point_on_circle(self.point, self.normal,
                                       self.radius, point)

    def accepts(self, point):
        "True if the point is on the circle."
        return is_point_on_circle(self.point, self.normal,
                                  self.radius, point)

def _circular_linear_test():
    from testing import __ok__
    c = Circular(V(), V.Y, 1)
    l = Linear(V.Y, V(2,-1,0))
    __ok__(str(c*l) == "Unity(V(1.0, 0.0, 0.0))")
    l = Linear(V.Y, -V.Y)
    __ok__(str(c*l) == "Nowhere()")
    l = Linear(V(0.5, 0, 0.5), V(0.5, 0, -0.5))
    cl = c*l
    __ok__(isinstance(cl, Duality))
    __ok__(cl.pointA[0] == cl.pointB[0])

def _circular_planar_test():
    from testing import __ok__
    # coplanar
    c = Circular(V(-1,0,0),V.Y,1)
    p = Planar(V( 1,0,0),V.Y)
    cp = c * p
    __ok__(cp is c)
    # parallel planar
    c = Circular(V(-1,0,0),V.Y,1)
    p = Planar(V( 1,1,0),V.Y)
    cp = c * p
    __ok__(isinstance(cp, Nowhere))
    # grab the tangent
    c = Circular(V(0,0,0),V.Y,1)
    p = Planar(V(0,1,0),V(1,1,0))
    cp = c * p
    __ok__(isinstance(cp, Unity))
    __ok__(equal(cp.point[0], 1))
    # grab a duality
    c = Circular(V(0,0,0),V.Y,1)
    p = Planar(V(0,.5,0),V(1,1,0))
    cp = c * p
    __ok__(isinstance(cp, Duality))
    __ok__(equal(abs(cp.pointA[2]), 0.5*sqrt(3)))

def _circular_circular_test():
    from testing import __ok__
    # tangent at origin
    c1 = Circular(V(-1,0,0),V.Y,1)
    c2 = Circular(V( 1,0,0),V.Y,1)
    c1c2 = c1 * c2
    __ok__(isinstance(c1c2, Unity))
    __ok__(c1c2.point == V())
    # cross off origin
    c1 = Circular(V.O,V.Y,1)
    c2 = Circular(V.X,V.Y,1)
    c1c2 = c1 * c2
    __ok__(isinstance(c1c2, Duality))
    __ok__(c1c2.pointA[0] == c1c2.pointB[0])
    __ok__(equal(c1c2.pointA[1], -c1c2.pointB[1]))
    __ok__(equal(abs(c1c2.pointA[1]), 0.5*sqrt(3)))

#----------------------------------------------------------------------------

def is_point_on_sphere(center, r, p):
    return equal(r, center.distance(p))

def closest_point_on_sphere(center, r, p):
    dir = p - center
    if dir.zero():
        return center + V(r, 0, 0)
    dir = dir.magnitude(r)
    return center + dir

class Spherical (Constraint):
    "A Spherical constraint allows a point anywhere on a sphere's surface."
    
    def __init__(self, point, radius):
        Constraint.__init__(self)
        if zero(radius):
            raise ValueError('A sphere must have a non-zero radius.')
        self.point = V(point)
        self.radius = abs(float(radius))

    def __repr__(self):
        return "Spherical" + repr((self.point, self.radius))

    def get_random(self):
        
        phi = 2*math.pi*random.random()
        theta = math.pi*random.random()
        return self.point + self.radius * ( V( [math.sin(theta)*math.cos(phi),
                                                math.sin(theta)*math.sin(phi),
                                                math.cos(theta)]))
    def combine(self, other):
        "Spherical versus another constraint."
        if isinstance(other, Anywhere): return self
        if isinstance(other, Nowhere): return other
        if isinstance(other, Unity): return other.combine(self)
        if isinstance(other, Duality): return other.combine(self)
        if isinstance(other, Linear): return self._spherical_linear(other)
        if isinstance(other, Planar): return self._spherical_planar(other)
        if isinstance(other, Circular): return self._spherical_circular(other)
        if isinstance(other, Spherical): return self._spherical_spherical(other)
        return Nowhere()

    def _spherical_linear(self, other):
        mark = other.closest(self.point)
        aa = mark.dsquared(self.point)
        cc = self.radius * self.radius
        if aa > cc + EPSILON:
            return Nowhere()
        if aa < cc - EPSILON:
            chord = V(other.dir)
            b = sqrt(cc - aa)
            chord = chord.magnitude(b)
            return Duality(mark + chord, mark - chord)
        return Unity(mark)

    def _spherical_planar(self, other):
        mark = other.closest(self.point)
        aa = mark.dsquared(self.point)
        cc = self.radius * self.radius
        if aa > cc + EPSILON:
            return Nowhere()
        if aa < cc - EPSILON:
            b = sqrt(cc - aa)
            return Circular(mark, mark - point, b)
        return Unity(mark)

    def _spherical_circular(self, other):
        # find our circle in their plane, then circular vs circular
        planar = Planar(other.point, other.normal)
        intersect = self._spherical_planar(planar)
        if not isinstance(intersect, Circular): return intersect
        return intersect.combine(other)

    def _spherical_spherical(self, other):
        cc = other.point - self.point
        d = cc.magnitude()
        rS = self.radius
        rO = other.radius
        # too far apart?
        if (rS + rO) < d - EPSILON:
            return Nowhere()
        # equal or circumscribed?
        if zero(d):
            if equal(rS, rO):
                return self
            return Nowhere()
        # tangent?
        if equal(d, rS+rO):
            cc = cc.magnitude(rS)
            return Unity(self.point + cc)
        if max(rS,rO) - (d+min(rO,rS)) > 0:
            #One sphere included in another
            return Nowhere()
        # find the midpoint where the "radical line" (chord) crosses cc
        #print d,rS,rO
        #print rS - (d+rO) 
        term = d*d - rS*rS + rO*rO
        #print term
        #term =  rS*rS - rO*rO
        xS = term / (2*d)
        #print xS
        xO = d - xS
        # find the length of the radical line
        #leg = 4*d*d*rS*rS - term*term
        #yy = 1/d * sqrt(leg)
        #y = yy/2
        y = sqrt(rS*rS-xO*xO)
        # find the intersection points
        
        cc.normalize()
        #print cc,xS
        mark = cc
        mark = mark.magnitude(xO)
        mark += self.point
        return Circular(mark, cc, y)

    def closest(self, point):
        "A point on the sphere is returned."
        return closest_point_on_sphere(self.point, self.radius, point)

    def accepts(self, point):
        "True if the point is on the sphere."
        return is_point_on_sphere(self.point, self.radius, point)

def _spherical_linear_test():
    from testing import __ok__
    # no intersection is Nowhere
    # tangent intersection is Unity
    # chord intersection is Duality
    pass

def _spherical_planar_test():
    from testing import __ok__
    # no intersection is Nowhere
    # tangent intersection is Unity
    # chord intersection is Circular
    pass

def _spherical_circular_test():
    from testing import __ok__
    pass

def _spherical_spherical_test():
    from testing import __ok__
    pass

#----------------------------------------------------------------------------

def __test__():
    from testing import __ok__, __report__
    print('Testing geometrical constraint classes...')

    # trivial constraints
    a1 = Anywhere()
    n1 = Nowhere()
    c1 = Circular(V.O, V.X, 1)
    __ok__(isinstance(a1 * a1, Anywhere))
    __ok__(isinstance(a1 * n1, Nowhere))
    __ok__(isinstance(n1 * n1, Nowhere))
    __ok__(isinstance(a1 * c1, Circular))
    __ok__((a1 * c1) is c1)
    __ok__(isinstance(n1 * c1, Nowhere))

    # parametric constraints
    _linear_linear_test()
    _planar_linear_test()
    _planar_planar_test()
    _circular_linear_test()
    _circular_planar_test()
    _circular_circular_test()
    _spherical_linear_test()
    _spherical_planar_test()
    _spherical_circular_test()
    _spherical_spherical_test()

    __report__()

if __name__ == '__main__':
    raise Exception('This module is not a stand-alone script.  Import it in a program.')
