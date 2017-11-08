from math import asin, atan2, degrees
import numpy as np

def histogram(hits):
    NI = 180
    NBOUNCES = 10
    hist = np.zeros((NI,NBOUNCES), dtype=np.int)
    for hit in hits:
        h = hit[-1]
        bounce = h.bounce
        deg = int(degrees(atan2(h.r.x,h.r.z))+0.5)
        hist[deg+90,bounce] += 1
    return hist / len(hits)

