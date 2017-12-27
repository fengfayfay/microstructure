from math import cos, sin, sqrt

class Vec3():
    def __init__(self, x, y, z):
        (self.x, self.y, self.z) = (x, y, z)
    def __str__(self):
        return 'Vec3({}, {}, {})'.format(self.x, self.y, self.z)
    def __mul__(self, other):
        return Vec3(self.x * other, self.y * other, self.z * other)
    def __add__(self, other):
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)
    def __sub__(self, other):
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)
    def __neg__(self):
        return Vec3(-self.x, -self.y, -self.z)
    def dot(self, other):
        return (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
    def __abs__(self):
        return self.dot(self)
    def norm(self):
        mag = sqrt(abs(self))
        if mag == 0.:
            return self
        return self * (1.0 / mag)

X = Vec3(1.,0.,0.)
Y = Vec3(0.,1.,0.)
Z = Vec3(0.,0.,1.)


def FromAngle(theta):
    return Vec3(sin(theta), 0., cos(theta))

def FromCosAngle(costheta):
    sintheta = sqrt(1.-costheta*costheta)
    return Vec3(sintheta, 0., costheta)


def dot(v, w):
    return v.dot(w)

#
# assume N and I are unit vectors
# assume I points away from the surface
#
def reflect(N, I):
    #print('reflect({},{}'.format(N,I))
    assert dot(N,I) >= 0
    return -I + N * (2*dot(N,I))

def halfvector(I,O):
    return (I+O).norm()

