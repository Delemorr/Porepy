import porepy as pp
import numpy as np
def triangle_fractur(L,b):
    ax = b/2
    ay = b/2*(3)**(1/2)
    Nx = int(L/b) +1
    Ny = int(L/(2*ay)) + 2
    X = []
    X2 =[]
    Y = []
    fracture = []
    for i in range (Nx):
        x = ax + 2*ax*i
        X.append(x)

    for i in range (Nx):
        x = b*i
        X2.append(x)

    for j in range (Ny):
        y = ay + 2*ay*j
        Y.append(y)

    g=(Y[-2] + ay - L)/2

    Y = []
    for j in range (Ny):
        y = ay + 2*ay*j - g
        Y.append(y)

    for j in range(1,Nx+1):

        if j <= Ny:
            frac1 = pp.LineFracture(np.array([[X[j-1],0],
                                              [-g, Y[j-1]]]))
            fracture.append(frac1)

        if j > Ny:
            frac1 = pp.LineFracture(np.array([[X[j-1],X2[j-Ny]],
                                              [-g, Y[-1]]]))
            fracture.append(frac1)

        if j <= Ny -1:
            frac2 = pp.LineFracture(np.array([[X2[j],0],
                                              [Y[-1], Y[Ny-j-1]]]))
            fracture.append(frac2)

        if j > Ny:
            frac2 = pp.LineFracture(np.array([[X2[j-1],X[j-Ny-1]],
                                              [Y[-1], -g]]))
            fracture.append(frac2)

        if j <= Ny-1:
            frac3 = pp.LineFracture(np.array([[X[Nx-j-1],X[-1]],
                                              [-g, Y[j]-ay]]))
            fracture.append(frac3)

        if j <= Ny-1:
            frac4 = pp.LineFracture(np.array([[X2[-j],X[-1]],
                                              [Y[-1], Y[-j]-ay]]))
            fracture.append(frac4)



    for j in range(Ny*2):
        frac5 = pp.LineFracture(np.array([[0,X[-1]],
                                          [-g+j*ay, -g+j*ay]]))
        fracture.append(frac5)

    return fracture
