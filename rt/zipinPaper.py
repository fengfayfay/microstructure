import math as mt
import sys as system

MAXBOUNCE = 10 

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


def far(theta, phi, maxX, hasMaxX, range):

    gAngle = 2.0 * theta 
    psi_min = mt.pi - 2.0 *(theta + phi)

    if isNinety(phi):
        y_h_max = 0
    else:
        if hasMaxX:
            y_h_max = 1
            offset = 0
        else:
            y_h_max =  2.0 * mt.tan(theta)  / (mt.tan(theta) + mt.tan(phi)) 
            s = (1.0-y_h_max) * mt.sin(theta+phi)
            offset = mt.asin(s)

    psi_max = mt.pi - (theta + phi) - offset

    n_min= mt.ceil(psi_min/gAngle)
    n_max = mt.ceil(psi_max/gAngle)
    xi_min = mt.pi - phi - n_min * gAngle
    xi_max = mt.pi - phi - n_max * gAngle

    x_max_intersect = XIntersect(psi_max, theta, phi)
    x_min_intersect = XIntersect(psi_min, theta, phi)

    x_critical_intersect = XIntersect(n_min * gAngle, theta, phi);
    print('critical value of x: ' + repr(x_critical_intersect))
    print('ratio of min x: ' +  repr((x_critical_intersect-x_min_intersect)/(2.0*mt.sin(theta))))

    if n_min%2:
        xi_min *= -1 
    if n_max%2:
        xi_max *= -1 
   
    hits = []
    min_ratio = (x_critical_intersect - x_min_intersect)/range
    max_ratio = (x_max_intersect - x_critical_intersect)/range
    hits.append(((mt.degrees(xi_min), n_min), min_ratio, (x_min_intersect, x_critical_intersect), 'left'))
    hits.append(((mt.degrees(xi_max), n_max), max_ratio, (x_critical_intersect, x_max_intersect), 'left')) 
    return hits

def near(theta, phi, range):
    gAngle = 2.0 * theta 
    psi_max = mt.pi - 2.0 * phi
    psi_min = mt.pi - theta -phi 

    x_near_min = mt.cos(theta) * mt.tan(phi)
    #x_near_min = mt.cos(theta) * mt.tan(theta-phi)
    print('minimum x value for near hit: ' + repr(x_near_min))

    n_min= mt.ceil(psi_min/gAngle)
    n_max = mt.ceil(psi_max/gAngle)
    xi_min = mt.pi + phi - n_min * gAngle
    xi_max = mt.pi + phi - n_max * gAngle

    if n_min%2 == 0:
        xi_min *= -1 
    if n_max%2 == 0:
        xi_max *= -1 
    
    x_min_intersect = x_near_min
    x_max_intersect = mt.sin(theta)
    hits = []
    if n_min == n_max:
        ratio = (x_max_intersect-x_min_intersect)/range
        hits.append(((mt.degrees(xi_min), n_min), ratio, (x_min_intersect, x_max_intersect), 'right'))
    else:
        x_critical_intersect = (x_max_intersect + x_min_intersect) * 0.5
        if phi != theta:
            yc = mt.cos(theta) * (1.0 - mt.sin(2.0 * n_min * theta + theta - phi)/ mt.sin(theta-phi))
            x_critical_intersect = mt.sin(theta) - yc * mt.tan(theta)    
        min_ratio = (x_critical_intersect - x_min_intersect)/range
        max_ratio = (x_max_intersect - x_critical_intersect)/range
        hits.append(((mt.degrees(xi_min), n_min), min_ratio, (x_min_intersect, x_critical_intersect), 'right'))
        hits.append(((mt.degrees(xi_max), n_max), max_ratio, (x_critical_intersect, x_max_intersect), 'right')) 
 
    return (hits, x_near_min)

def zipinPaper(theta, phi, isDegrees = True):

    if isDegrees:
        thetaR = mt.radians(theta)
        phiR = mt.radians(phi)
    else:
        thetaR = theta
        phiR = phi

    xmax = mt.sin(thetaR)
    print((-xmax, xmax))

    if phiR > thetaR:
        return far(thetaR, phiR, mt.sin(thetaR), False, xmax*2)
    else:
        (hits, x_near_min) = near(thetaR, phiR, xmax*2) 
        hits += far(thetaR, phiR, x_near_min, True, xmax*2)
        return hits


def zipinBrdf(phi, chi):
    n = 1
    while n < MAXBOUNCE:
        theta = (180 - phi - chi) * .5 / n
        hits = zipinPaper(theta, phi)
        print(hits)
        n += 1
    
    

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print('python3 zipin.py theta phi')
        sys.exit(0)

    #theta = float(sys.argv[1])
    phi = float(sys.argv[1])
    chi = float(sys.argv[2])

    #print(zipinPaper(theta, phi))
    print(zipinBrdf(phi, chi))

