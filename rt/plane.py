from vec3 import dot
from ray import TMAX

class Plane:
    def __init__(self, P, N):
        self.P = P
        self.N = N.norm()

    def intersect(self, ray):
        O = ray.O
        D = ray.D
        P = self.P
        N = self.N

        # Return the distance from O to the intersection of the ray (O, D) with the 
        # plane (P, N), or +inf if there is no intersection.
        # O and P are 3D points, D and N (normal) are normalized vectors.
        denom = dot(D, N)
        if abs(denom) < 1e-6:
            return TMAX
        d = dot(P - O, N) / denom
        if d < 0:
            return TMAX
        return d
