import math;
import vec3;


def makeVector(theta, phi):
    sinTheta = math.sin(theta)
    cosTheta = math.cos(theta)
    sinPhi = math.sin(phi)
    cosPhi = math.cos(phi)
    return vec3.Vec3(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta)

def multiply(rotation, w):
    x = vec3.dot(rotation[0], w)
    y = vec3.dot(rotation[1], w)
    z = vec3.dot(rotation[2], w)
    return vec3.Vec3(x, y, z) 


class Rotate:

    def __init__(self, wh):
        self.wh = wh
        self.computeRotation()
    
    def computeRotation(self):
        self.phi = math.atan2(self.wh.y, self.wh.x)
        #print(math.degrees(self.phi))
        
        gamma = -self.phi
        sinGamma = math.sin(gamma)
        cosGamma = math.cos(gamma)
        self.rotation = []
        self.rotation.append(vec3.Vec3(cosGamma, -sinGamma, 0))
        self.rotation.append(vec3.Vec3(sinGamma, cosGamma, 0))
        self.rotation.append(vec3.Vec3(0, 0, 1))
        self.inverse = []
        self.inverse.append(vec3.Vec3(cosGamma, sinGamma, 0))
        self.inverse.append(vec3.Vec3(-sinGamma, cosGamma, 0))
        self.inverse.append(vec3.Vec3(0, 0, 1))

    def rotate(self, w):
        wr = multiply(self.rotation, w)
        return wr

    def inverse_rotate(self, w):
        wr = multiply(self.inverse, w)
        return wr


               
        
