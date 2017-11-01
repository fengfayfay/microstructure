from math import sqrt

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

def dot(v, w):
    return v.dot(w)
