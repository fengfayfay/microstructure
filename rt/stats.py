import numpy as np

def histogram(hits):
    NCOSI = 20
    NBOUNCES = 10
    hist = np.zeros((NCOSI,NBOUNCES), dtype=np.int)
    for hit in hits:
        bounce = hit.bounce
        cosr = int(NCOSI*hit.r.D.z)
        hist[cosr,bounce-1] += 1
    return hist / len(hits)


        
