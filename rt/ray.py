from vec3 import dot
TMAX = 1e36

class Ray():
    def __init__(self, O, D):
        self.O = O
        self.D = D # D.norm()
        self.hits = []

    def __str__(self):
        return 'Ray({}, {})'.format(self.O, self.D)

    def __call__(self, t):
        return self.O + self.D * t

class Hit():
    def __init__(self, P, N, side=None, bounce=None):
        self.P = P
        self.N = N
        self.side = side
        self.bounce = bounce

    def __str__(self):
        return 'Hit({}, {})'.format(self.P, self.N)


# assume I heads towards the surface
def reflect(N, I):
    #print('reflect({},{}'.format(N,I))
    return I - N*(2*dot(N,I))
