import math as mt
import sys as system

MAXBOUNCE = 10 

def clamp(x, xmax):
    x /= xmax
    if x < -1.: return -1.
    if x >  1.: return  1
    return x


def addHit(hits, xi, bcount,  side, z, interval, range):
    ratio = (interval[1] - interval[0])/range
    
    if ratio > 1e-6:
        xi = mt.degrees(xi)
        hits.append((xi, ratio, bcount, side, z))
    

def far(theta, phi, maxX, hasNear, hitrange):

    gAngle = 2.0 * theta 
    psi_min = mt.pi - 2.0 *(theta + phi)

    sinThetaPhi = mt.sin(theta+phi)

    psi_max = mt.pi - (theta+phi)

    if hasNear == False and sinThetaPhi > 0:
        #this is exactly Gi because when all the rays go to the left side
        #the max left side visible is Gi
        zMax = hitrange * mt.cos(phi)/ sinThetaPhi 
        psi_max -= mt.asin((1.0-zMax) * sinThetaPhi)
    else:
        zMax = 1

    #print(mt.degrees(psi_min), mt.degrees(psi_max))

    n_min= mt.ceil(psi_min/gAngle)
    n_max = mt.ceil(psi_max/gAngle)
    '''
    if psi_min == gAngle * n_min:
        n_min += 1
    if psi_max == gAngle * n_max:
        n_max += 1
    '''

    xi_min = mt.pi - phi - n_min * gAngle
    xi_max = mt.pi - phi - n_max * gAngle

    
    x_min_intersect = -hitrange * .5
    x_max_intersect = maxX
    xrange = maxX - x_min_intersect
    
    if n_min%2:
        xi_min *= -1 
    if n_max%2:
        xi_max *= -1 

    hits = []
    if n_max > n_min:
        #compute z critical intersect, length from circle center
        k = n_min
        criticalAngle = mt.pi - (theta+phi) - k* gAngle
        if criticalAngle <= 0:
            z = 0
        else:
            z = mt.sin(criticalAngle)/sinThetaPhi 
    
        x_critical_intersect = x_min_intersect + (1-z) / zMax  * xrange

        addHit(hits, xi_min, n_min, 'left', 1-z, (x_min_intersect, x_critical_intersect), hitrange)
        addHit(hits, xi_max, n_max, 'left', zMax - (1-z), (x_critical_intersect, x_max_intersect), hitrange)
    else:
        addHit(hits, xi_min, n_min, 'left', zMax, (x_min_intersect, x_max_intersect), hitrange)
   
    return hits

def near(theta, phi, hitrange):
    gAngle = 2.0 * theta 
   
    #Feng's correction to near psi computation from Paper 
    psi_max = mt.pi -(theta-phi)
    psi_min = mt.pi - 2 * (theta-phi)
    

    ##Zipin Paper near psi computation 
    #psi_max = mt.pi - 2.0 * phi
    #psi_min = mt.pi - theta - phi

    x_near_min = mt.cos(theta) * mt.tan(phi)

    n_min= mt.ceil(psi_min/gAngle)
    n_max = mt.ceil(psi_max/gAngle)

    #print(n_min, n_max)
    '''
    if psi_min == gAngle * n_min:
        n_min += 1
    if psi_max == gAngle * n_max:
        n_max += 1
    '''

    xi_min = mt.pi + phi - n_min * gAngle
    xi_max = mt.pi + phi - n_max * gAngle
    #print(xi_min, xi_max);
  
    ###PROBLEM WITH ZIPIN#MUST#FIGUREOUTWHYTHISMAKESENSE### 
    '''
    if (xi_min > mt.pi * .5) or (xi_max > mt.pi * .5):
        print(n_min, n_max)
        print (mt.degrees(xi_min), mt.degrees(xi_max))
        print (mt.degrees(theta), mt.degrees(phi))
    '''
    

    if n_min%2 == 0:
        xi_min *= -1 
    if n_max%2 == 0:
        xi_max *= -1 

    #print (n_min, n_max)
    
    #xi_min = -xi_min
    #xi_max = -xi_max
    
    x_min_intersect = x_near_min
    x_max_intersect = mt.sin(theta)
    #print (x_min_intersect, x_max_intersect, hitrange)

    hits = []
    if n_min == n_max or theta - phi < .000001:
        addHit(hits, xi_min, n_min, 'right', 1.0, (x_min_intersect, x_max_intersect), hitrange)
    else:
        k = n_min
        criticalAngle = mt.pi - (theta-phi) - k * gAngle
        if criticalAngle <= 0:
            z = 0
        else:
            z = mt.sin(criticalAngle) /mt.sin(theta-phi)
        x_critical_intersect = x_min_intersect + (x_max_intersect - x_min_intersect) * z
        addHit(hits, xi_min, n_min, 'right', 1.0-z, (x_critical_intersect, x_max_intersect), hitrange)
        addHit(hits, xi_max, n_max, 'right', z, (x_min_intersect, x_critical_intersect), hitrange)
    return (hits, x_near_min)

def zipin(theta, phi):

    thetaR = mt.radians(theta)
    phiR = mt.radians(phi)

    xmax = mt.sin(thetaR)

    if phiR > thetaR:
        return far(thetaR, phiR, mt.sin(thetaR), False, xmax*2)
    else:
        (hits, x_near_min) = near(thetaR, phiR, xmax*2) 
        hits += far(thetaR, phiR, x_near_min, True,  xmax*2)
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

