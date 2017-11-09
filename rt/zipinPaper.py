import math as mt
import sys as system

def clamp(x, xmax):
    x /= xmax
    if x < -1.: return -1.
    if x >  1.: return  1
    return x


def isNinety(phi):
    return abs(phi - mt.pi * .5) < .000001

def far(theta, phi):

    gAngle = 2.0 * theta 
    psi_min = mt.pi - 2.0 *(theta + phi)
    if isNinety(phi):
        y_h_max = 0
    else:
        y_h_max = mt.tan(theta)/mt.tan(phi)

    s = (1.0-y_h_max) * mt.sin(theta+phi)
    offset = mt.asin(s)
    psi_max = mt.pi - (theta + phi) - offset

    n_min= mt.ceil(psi_min/gAngle)
    n_max = mt.ceil(psi_max/gAngle)
    xi_min = mt.pi - phi - n_min * gAngle
    xi_max = mt.pi - phi - n_max * gAngle

    if n_min%2:
        xi_min *= -1 
    if n_max%2:
        xi_max *= -1 
   
    hits = []
    hits.append((n_min, xi_min))
    hits.append((n_max, xi_max)) 
    return hits

def near(theta, phi):
    gAngle = 2.0 * theta 
    psi_min = mt.pi - 2.0 * phi
    if isNinety(phi):
        y_h_max = 0
    else:
        y_h_max = mt.tan(theta)/mt.tan(phi)

    s = (1.0-y_h_max) * mt.sin(theta-phi)
    offset = mt.asin(s)
    psi_max = mt.pi - (theta + phi) + offset

    n_min= mt.ceil(psi_min/gAngle)
    n_max = mt.ceil(psi_max/gAngle)
    xi_min = mt.pi + phi - n_min * gAngle
    xi_max = mt.pi + phi - n_max * gAngle

    if n_min%2 == 0:
        xi_min *= -1 
    if n_max%2 == 0:
        xi_max *= -1 
   
    hits = []
    hits.append((n_min, xi_min))
    hits.append((n_max, xi_max)) 
    return hits

def zipinPaper(theta, phi):
    thetaR = mt.radians(theta)
    phiR = mt.radians(phi)

    if phiR > thetaR:
        return far(thetaR, phiR)
    else:
        return near(thetaR, phiR)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print('python3 zipin.py theta phi')
        sys.exit(0)

    theta = float(sys.argv[1])
    phi = float(sys.argv[2])

    print(zipinPaper(theta, phi))

