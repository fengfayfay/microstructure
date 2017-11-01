from vec3 import Vec3
from ray import Ray
from vgroove import VGroove
from mc import uniform

v = VGroove(1.)

lside = 0
rside = 0
for i in range(1000):
    r = Ray( Vec3(uniform(-1.0, 1.0), 0.0, 0.0), Vec3(0.0, 0.0, -1.0))
    #print(r)
    t = v.intersect(r)
    p = r(t)
    if p.x < 0.:
        lside += 1
    else:
        rside += 1
print(lside, rside)
