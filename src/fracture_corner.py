import porepy as pp
import numpy as np
from src.rotate import rotate
import math
def corner_fractur(b,fi,L,ntr):
    fracture = []
    x = L/2
    y = L/2

    x1 = 0 - 100*b
    y1 = L/2

    x3 = L + 100*b
    y3 = L/2

    origin = (x, y)
    point1 = (x1, y1)
    point3 = (x3, y3)

    point_nov1 = rotate(origin, point1, math.radians(fi))
    point_nov3 = rotate(origin, point3, math.radians(fi))

    frac = pp.LineFracture(np.array([[point_nov1[0], point_nov3[0]], [point_nov1[1], point_nov3[1]]]))
    fracture.append(frac)

    for i in range(ntr):
        x1 = point_nov1[0] + i*b*math.cos( math.radians(90) + math.radians(fi))
        y1 = point_nov1[1] + i*b*math.sin( math.radians(90) + math.radians(fi))

        x3 = point_nov3[0] + i*b*math.cos( math.radians(90) + math.radians(fi))
        y3 = point_nov3[1] + i*b*math.sin( math.radians(90) + math.radians(fi))

        x11 = point_nov1[0] - i*b*math.cos( math.radians(90) + math.radians(fi))
        y11 = point_nov1[1] - i*b*math.sin( math.radians(90) + math.radians(fi))

        x31 = point_nov3[0] - i*b*math.cos( math.radians(90) + math.radians(fi))
        y31 = point_nov3[1] - i*b*math.sin( math.radians(90) + math.radians(fi))


        frac = pp.LineFracture(np.array([[x1, x3], [y1, y3]]))
        fracture.append(frac)
        frac = pp.LineFracture(np.array([[x11, x31], [y11, y31]]))
        fracture.append(frac)

    return fracture

