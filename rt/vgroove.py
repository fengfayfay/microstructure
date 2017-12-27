from math import sqrt, cos, radians
from vec3 import Vec3, dot, halfvector
from ray import Ray, Hit, reflect, TMAX
from plane import Plane

#
# VGroove with axis along y
#
class VGroove:
    def __init__(self, h):
        nz = cos(radians(h))
        nx = sqrt(1-nz*nz)
        height = nz/nx
        #print(nx, nz, height)
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
                    hit = Hit(p, self.lside.N, 'left')
                elif tmin == t2:
                    #print('hit right side')
                    hit = Hit(p, self.rside.N, 'right')
                #print(hit)
        return hit

    def trace(self, ray):
        hits = []
        hit = self.intersect(ray)
        bounce = 1
        while hit:
            hit.bounce = bounce
            hit.i = ray.D
            hit.r = reflect(hit.N, ray.D)
            hits.append(hit)
            ray = Ray(hit.P, hit.r)
            bounce += 1
            hit = self.intersect(ray)
        return hits
       
def G1(N,H,I):
    HdotI = dot(H,I)
    if HdotI <= 0.:
        return 0.
    NdotH = dot(N,H)
    NdotI = dot(N,I)
    return 2. * NdotH * NdotI / HdotI

def G2(N,H,I,O):
    return min(1,G1(N,H,I),G1(N,H,O))

