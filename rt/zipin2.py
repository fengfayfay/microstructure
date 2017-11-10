from math import cos, sin, radians, ceil

def clamp(val, min, max):
    if val < min: val = min
    if val > max: val = max
    return val

def far(theta, phi):

    delta = theta + phi
    y_min =  0
    y_max =  sin(radians(delta))
    #print(y_min, y_max, y_max-y_min)
    psi_min = 180 -   delta
    psi_max = 180 - 2*delta
    n_min = int(ceil(psi_min/(2*theta)))
    n_max = int(ceil(psi_max/(2*theta)))

    hits = []

    y1 = sin(radians(2*(n_min-1)*theta+theta+phi))
    y2 = sin(radians(2*n_min*theta+theta+phi))
    #print(n_min, y1, y2)
    y1 = clamp( y1, y_min, y_max)
    y2 = clamp( y2, y_min, y_max)
    #print(n_min, y1, y2)
    xi_min = (-1 if n_min % 2 == 0 else 1) * (180 - phi - 2*n_min*theta)
    hits += [(xi_min, abs(y2-y1), n_min, 'left')]

    if n_min != n_max:
        y1 = sin(radians(2*(n_max-1)*theta+theta+phi))
        y2 = sin(radians(2*n_max*theta+theta+phi))
        y1 = clamp( y1, y_min, y_max)
        y2 = clamp( y2, y_min, y_max)
        #print(n_max, y1, y2)
        xi_max = (-1 if n_max % 2 == 0 else 1) * (180 - phi - 2*n_max*theta)
        hits += [(xi_max, abs(y2-y1), n_max, 'left')]

    return hits

def near(theta, phi):
    if theta < phi:
        return []

    delta = theta - phi
    y_max =  sin(radians(delta))
    y_min =  0
    #print(y_min, y_max, y_max-y_min)
    psi_max = 180 - 2*delta
    psi_min = 180 -   delta
    n_max = int(ceil(psi_max/(2*theta)))
    n_min = int(ceil(psi_min/(2*theta)))

    hits = []

    #psi1 = 2*(n_min-1)*theta+theta-phi
    #psi2 = 2*n_min*theta+theta-phi
    #print('min', n_min, psi1, psi2)
    y1 = sin(radians(2*(n_min-1)*theta+theta-phi))
    y2 = sin(radians(2*n_min*theta+theta-phi))
    y1 = clamp( y1, y_min, y_max)
    y2 = clamp( y2, y_min, y_max)
    #print('min', n_min, y1, y2)
    xi_min = (-1 if n_min % 2 == 1 else 1) * (180 + phi - 2*n_min*theta)
    hits += [(xi_min, abs(y2-y1), n_min, 'right')]

    if n_max != n_min:
        y1 = sin(radians(2*(n_max-1)*theta+theta-phi))
        y2 = sin(radians(2*n_max*theta+theta-phi))
        y1 = clamp( y1, y_min, y_max)
        y2 = clamp( y2, y_min, y_max)
        #print('max', n_max, y1, y2)
        xi_max = (-1 if n_max % 2 == 1 else 1) * (180 + phi - 2*n_max*theta)
        hits += [(xi_max, abs(y2-y1), n_max, 'right')]

    return hits

def zipin(theta, phi):
    return far(theta,phi)+near(theta, phi)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print('python3 zipin.py theta phi')
        sys.exit(0)

    theta = float(sys.argv[1])
    phi = float(sys.argv[2])

    print(zipin(theta, phi))
    ymax = sin(radians(phi+theta))
    ymin = sin(radians(phi-theta))
    print(ymin, ymax)

