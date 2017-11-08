from math import cos, sin, radians
import sys

def clamp(x, xmax):
    x /= xmax
    if x < -1.: return -1.
    if x >  1.: return  1
    return x

def back(theta, psi, phi):
    costheta = cos(radians(theta))
    sintheta = sin(radians(theta))
    cospsi = cos(radians(psi))
    sinpsi = sin(radians(psi))
    cosphi = cos(radians(phi))
    sinphi = sin(radians(phi))
    #print(sinphi,cosphi)
    #print(sinpsi,cospsi)
    #print(-sinpsi, costheta-cospsi, (costheta-cospsi)*sinphi, -sinpsi+(costheta-cospsi)*sinphi)
    return clamp(-sinpsi+(costheta-cospsi)*sinphi/cosphi, sintheta)


def far(theta, phi):
    hits = []
    bounce = 1
    sign = 1.
    psi = theta
    xi = 180-phi-2*theta
    while True:
        xmin = back(theta, psi, phi)
        xmax = back(theta, min(psi+2*theta, 180-phi), phi)
        if xmin < xmax:
            hits.append((sign*xi,(xmax-xmin)/2,bounce,'left'))
            #print('f', bounce, sign*xi)
        psi += 2*theta
        if psi > 180-phi:
            break
        xi -= 2*theta
        bounce += 1
        sign *= -1
    return hits

def near(theta, phi):
    hits = []
    if phi <= theta: 
        bounce = 1
        sign = -1.
        psi = -theta
        xi = 180+phi-2*theta
        while True:
            xmax = back(theta, psi, phi)
            xmin = back(theta, max(psi-2*theta, -180-phi), phi)
            if xmin < xmax:
                hits.append((sign*xi,(xmax-xmin)/2,bounce,'right'))
            psi -= 2*theta
            if psi < -180-phi:
                break
            xi -= 2*theta
            bounce += 1
            sign *= -1
    return hits

def zipin(theta, phi):
    return near(theta, phi) + far(theta, phi)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print('python3 zipin.py theta phi')
        sys.exit(0)

    theta = float(sys.argv[1])
    phi = float(sys.argv[2])

    print(zipin(theta, phi))

