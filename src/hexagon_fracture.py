import porepy as pp
import numpy as np
import math


def hexagon_fractur(L,b):
    ax = b*math.sqrt(3)/2
    ay = b*1/2
    fracture = []
    Nx = int(L/(2*ax)) +2
    Ny = int(L/(2 * b + 2 * ay))+2
    for i in range(Nx):
        for j in range(Ny):
            frac1 = pp.LineFracture(np.array([[0 + 2 * ax * i, ax + 2 * ax * i], [ay + (2*b+2*ay)*j, 0+ (2*b+2*ay)*j]]))
            fracture.append(frac1)

            if i < Nx-1:
                frac2 = pp.LineFracture(np.array([[ax + 2 * ax * i, 2*ax + 2 * ax * i],
                                                  [0 + (2 * b + 2 * ay) * j, ay + (2 * b + 2 * ay) * j]]))
                fracture.append(frac2)

            if j < Ny-1:
                frac3 = pp.LineFracture(np.array([[0 + 2 * ax * i, 0 + 2 * ax * i],
                                                  [ay + (2 * b + 2 * ay) * j, ay + b + (2 * b + 2 * ay) * j]]))
                fracture.append(frac3)

            if i < Nx and j < Ny-1:
                frac4 = pp.LineFracture(np.array([[0 + 2 * ax * i, ax + 2 * ax * i],
                                                  [ay + b + (2*b+2*ay)*j, b+2*ay + (2*b+2*ay)*j]]))
                fracture.append(frac4)

            if i < Nx-1 and j < Ny - 1:
                frac5 = pp.LineFracture(np.array([[ax + 2 * ax * i, 2*ax + 2 * ax * i],
                                                  [0+b+2*ay + (2*b+2*ay)*j, ay + b + (2*b+2*ay)*j]]))
                fracture.append(frac5)
            #
            if i < Nx and j < Ny -1:
                frac6 = pp.LineFracture(np.array([[ax + 2 * ax * i, ax + 2 * ax * i],
                                                  [2*ay+b + (2 * b + 2 * ay) * j, 2*ay + 2*b + (2 * b + 2 * ay) * j]]))
                fracture.append(frac6)
    return fracture
