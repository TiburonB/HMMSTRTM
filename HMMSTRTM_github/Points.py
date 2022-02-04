# Points.py 
# TLB 08/23/2020
# This file contains the Point object and associated functions
# A Point constitutes a mean (x,y,z) value on a cartesian grid and stderror values:
        # (xerr,yerr,zerr)
import numpy as np

class Point:
    def __init__ (self, xu, yu, zu, xerr = 0, yerr = 0, zerr = 0):
        self.x = xu
        self.y = yu
        self.z = zu
        self.xerr = xerr
        self.yerr = yerr
        self.zerr = zerr
        return


    @classmethod
    def from_array(cls, arr):
        return cls(arr[0], arr[1], arr[2])

    def no_err(self):
        if self.xerr == 0 and self.yerr == 0 and self.zerr == 0:
            return True
        return False

    def mean_distance(self, point2):
        dx2 = (self.x-point2.x)**2 
        dy2 = (self.y-point2.y)**2 
        dz2 = (self.z-point2.z)**2 
        ud  = (dx2 + dy2 + dz2)**0.5
        return ud

    def distance(self, point2):
        return self.mean_distance(point2)

    def equals(self, point2):
        
        if self.no_err() and point2.no_err(): 
            # if there are no stderrs associated with either point then 
            # return true when the mean distance == 0 
            if self.mean_distance(point2) == 0: 
                return True
            else:
                return False
        else:
            print("Point stderrs not accounted for yet, returning False.")
            return False

    def __str__(self):
        return "(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")"

    def __add__ (self, p2):
        return Point(self.x + p2.x, self.y + p2.y, self.z + p2.z)
    def __sub__ (self, p2):
        return Point(self.x - p2.x, self.y - p2.y, self.z - p2.z)
    def __mul__(self, p2):
        return self.x * p2.x + self.y * p2.y + self.z * p2.z
    def __pow__(self, other):
        """Return VectorxVector (cross product) or Vectorxscalar."""
        if isinstance(other, Point):
            a, b, c = self.x, self.y, self.z
            d, e, f = other.x, other.y, other.z
            c1 = numpy.linalg.det(numpy.array(((b, c), (e, f))))
            c2 = -numpy.linalg.det(numpy.array(((a, c), (d, f))))
            c3 = numpy.linalg.det(numpy.array(((a, b), (d, e))))
            return Point(c1, c2, c3)
        else:
            return Point(self.x * other, self.y*other, self.z*other)

    def calc_dihedral(v1, v2, v3, v4):
        ab = v1 - v2
        cb = v3 - v2
        db = v4 - v3
        u = ab ** cb
        v = db ** cb
        w = u ** v
        angle = u.angle(v)
        # Determine sign of angle
        try:
            if cb.angle(w) > 0.001:
                angle = -angle
        except ZeroDivisionError:
            return math.pi
        return angle





class PointSet:
    def __init__(self, points):
        self.points = points
        self.n = len(self.points)
        self.sdev = 0

    @classmethod
    def from_mat(cls, poss):
        points = []
        for pos in poss:
             points.append(Point.from_array(pos))
        return cls(points)
    
    
    @classmethod
    def fromcoords(cls, xs, ys, zs):
        points = []
        for i in range(0,len(xs)):
             points.append(Point(xs[i],ys[i],zs[i]))
        return cls(points)

    def setsdev(self, val):
        self.sdev = val
        for point in self.points:
            point.xerr = val
            point.yerr = val
            point.zerr = val
        return

    def getRMSD(self, set2):
        sum_distance = 0
        for i in range(0,self.n):
            sum_distance += (self.points[i].distance(set2.points[i]))**2

        sum_avg_squared_dev = 0
        for i in range(0,self.n):
            sum_avg_squared_dev += self.points[i].xerr**2
            sum_avg_squared_dev += set2.points[i].xerr**2
            sum_avg_squared_dev += self.points[i].yerr**2
            sum_avg_squared_dev += set2.points[i].yerr**2
            sum_avg_squared_dev += self.points[i].zerr**2
            sum_avg_squared_dev += set2.points[i].zerr**2
        
        RMSD = sum_distance + sum_avg_squared_dev 
        RMSD = RMSD / self.n
        return RMSD


