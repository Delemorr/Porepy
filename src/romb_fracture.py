import porepy as pp
import numpy as np
def romb_fractur(L,b):
    x = int(L / b)
    X = []
    z = 0
    for i in range(x):
        z = z + b / 2
        X.append(z)
        z = z + b / 2
    fracture = []
    N = len(X)+1
    for j in range(N):

        if j >= 1:
            frac1 = pp.LineFracture(np.array([[X[j-1],0], [0, X[j-1]]]))
            fracture.append(frac1)

            frac2 = pp.LineFracture(np.array([[0,X[j-1]], [X[N-j-1], L]]))
            fracture.append(frac2)

            frac3 = pp.LineFracture(np.array([[X[j-1],L], [0, X[N-j-1]]]))
            fracture.append(frac3)

            frac4 = pp.LineFracture(np.array([[L,X[j-1]], [X[j-1], L]]))
            fracture.append(frac4)

    return fracture
