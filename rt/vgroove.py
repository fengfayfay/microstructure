from math import sqrt
from vec3 import Vec3, dot
from ray import Ray, Hit, reflect, TMAX
from plane import Plane

#
# VGroove with axis along y
#
class VGroove:
    def __init__(self, cosh):
        nz = cosh
        nx = 1-sqrt(nz*nz)
        height = 1/nx
        self.lside = Plane( Vec3(0.0, 0.0, -height), Vec3( nx, 0.0, nz) )
        self.rside = Plane( Vec3(0.0, 0.0, -height), Vec3(-nx, 0.0, nz) )

    def intersect(self, ray):
        t1 = None
        t2 = None
        if dot(ray.D, self.lside.N) < 0: # front facing
            t1 = self.lside.intersect(ray)
            #print('intersect left face', t1)
        if dot(ray.D, self.rside.N) < 0: # front facing
            t2 = self.rside.intersect(ray)
            #print('intersect right face', t2)

        tmin = None
        if t1 is None or t2 is None:
            if t1 is None:
               tmin = t2
            if t2 is None:
               tmin = t1
        else:
            tmin = min(t1, t2)

        hit = None
        if tmin is not None:
            p = ray(tmin)
            if p.z <= 0.:
                if   tmin == t1:
                    #print('hit left side')
                    hit = Hit(p, self.lside.N)
                elif tmin == t2:
                    #print('hit right side')
                    hit = Hit(p, self.rside.N)
                #print(hit)
        return hit

    def trace(self, ray):
        hits = []
        hit = self.intersect(ray)
        r = ray
        bounce = 1
        while hit:
            r = Ray(hit.P, reflect(hit.N, r.D))
            hit.r = r
            hit.bounce = bounce
            bounce += 1
            hits.append(hit)
            hit = self.intersect(r)
        return hits
       
        
