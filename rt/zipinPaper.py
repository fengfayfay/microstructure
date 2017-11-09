import math as mt
import sys as system

def clamp(x, xmax):
    x /= xmax
    if x < -1.: return -1.
    if x >  1.: return  1
    return x


def isNinety(phi):
    return abs(phi - mt.pi * .5) < .000001

def XIntersect(psi, theta, phi):
    angleA = mt.pi - (psi+theta+phi)
    thetaPlusPhi = theta + phi
    ycOverH =  mt.sin(angleA)/mt.sin(theta+phi)
    h = mt.cos(theta)
    y = h * (1.0 - ycOverH)
    x = y * mt.tan(theta) + y * mt.tan(phi) - mt.sin(theta) 
    return x


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

    x_max_intersect = XIntersect(psi_max, theta, phi)
    x_min_intersect = XIntersect(psi_min, theta, phi)

    if n_min%2:
        xi_min *= -1 
    if n_max%2:
        xi_max *= -1 
   
    hits = []
    hits.append((mt.degrees(xi_min), x_min_intersect, n_min, 'left'))
    hits.append((mt.degrees(xi_max), x_max_intersect, n_max, 'left')) 
    return hits

def near(theta, phi):
    gAngle = 2.0 * theta 
    psi_max = mt.pi - 2.0 * phi
    psi_min = mt.pi - theta -phi 

    n_min= mt.ceil(psi_min/gAngle)
    n_max = mt.ceil(psi_max/gAngle)
    xi_min = mt.pi + phi - n_min * gAngle
    xi_max = mt.pi + phi - n_max * gAngle

    if n_min%2 == 0:
        xi_min *= -1 
    if n_max%2 == 0:
        xi_max *= -1 
   
    x_min_intersect = 0
    x_max_intersect = 0
    hits = []
    hits.append((mt.degrees(xi_min), x_min_intersect, n_min, 'right'))
    hits.append((mt.degrees(xi_max), x_max_intersect, n_max, 'right')) 
    return hits

def zipinPaper(theta, phi):
    thetaR = mt.radians(theta)
    phiR = mt.radians(phi)

    xmax = mt.sin(thetaR)
    print((-xmax, xmax))

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

