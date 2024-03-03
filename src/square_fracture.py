import porepy as pp
import numpy as np
def square_fractur(L,b):
    x = int(L / b)
    X = []
    z = 0
    for i in range(x):
        z = z + b / 2
        X.append(z)
        z = z + b / 2
    fracture = []
    N = len(X)
    for j in range(N):
        frac1 = pp.LineFracture(np.array([[X[j], X[j]], [0, L]]))
        fracture.append(frac1)
        frac2 = pp.LineFracture(np.array([[0, L], [X[j], X[j]]]))
        fracture.append(frac2)
    return fracture
