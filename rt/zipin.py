from math import cos, sin, radians, degrees, ceil, pi
from vec3 import Vec3, Z
from vgroove import G1

EPSILON = 1e-8

def clamp(val, min, max):
    if val < min: val = min
    if val > max: val = max
    return val

def far(theta, phi):

    y_max =  sin(phi+theta)
    y_min =  sin(phi-theta)
    area = y_max - y_min
    y_min = max(0, y_min)

    #print('far',y_min, y_max)

    hits = []
    n = 1
    psi_min = theta + phi
    while True:
        psi_max = psi_min + 2*theta
        y1 = sin(psi_min)
        y2 = sin(psi_max)
        y1 = clamp(y1, y_min, y_max)
        y2 = clamp(y2, y_min, y_max)
        a = abs(y2-y1)
        #print(degrees(psi_min),y1,degrees(psi_max),y2,a)
        if a > EPSILON:
            sign = -1 if n % 2 == 0 else 1
            #xi_min = (-1 if n % 2 == 0 else 1) * (pi - phi - 2*n*theta)
            xi_min = sign * (pi - psi_min - theta)
            #print(degrees(psi_min+theta),degrees(xi_min))
            hits += [(degrees(xi_min), a/area, n, 'left')]
        n += 1
        psi_min = psi_max
        if psi_min > pi:
            break

    return hits

def near(theta, phi):
    if theta < phi:
        return []

    y_max =  sin(phi+theta)
    y_min =  sin(phi-theta)
    area = y_max - y_min
    y_max = min(0, y_max)
    #print('near',y_min, y_max)

    hits = []
    n = 1
    psi_min = phi - theta
    while True:
        psi_max = psi_min - 2*theta
        y1 = sin(psi_min)
        y2 = sin(psi_max)
        y1 = clamp(y1, y_min, y_max)
        y2 = clamp(y2, y_min, y_max)
        a = abs(y2-y1)
        #print(degrees(psi_min),y1,degrees(psi_max),y2,a)
        if a > EPSILON:
            sign = -1 if n % 2 == 1 else 1
            #xi_min = sign * (pi - phi - 2*n*theta)
            xi_min = sign * (pi + psi_min - theta)
            #print(y1,y2,degrees(psi_min-theta),degrees(xi_min))
            hits += [(degrees(xi_min), a/area, n, 'right')]
        n += 1
        psi_min = psi_max
        if psi_min < -pi:
            break

    return hits


def zipin(theta, phi):
    assert theta > 0.
    return far(theta,phi)+near(theta, phi)

def surface(phi, theta):
    H = Vec3(cos(theta), 0.0, sin(theta))
    I = Vec3(sin(phi), 0.0, cos(phi))
    return I, H

def Gzipin(theta, phi):
    N = Z
    for hit in zipin(theta,phi):
        angle, area, bounce, side = hit
        if side == 'left' and bounce == 1:
            H, I = surface(phi, theta)
            return area*G1(N,H,I)
    return 0

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print('python3 zipin.py theta phi')
        sys.exit(0)

    theta = radians(float(sys.argv[1]))
    phi = radians(float(sys.argv[2]))

    hits = zipin(theta, phi)
    print(hits)

