import numpy as np
def ekv(L,ro,mu,P1,P2,model):
    "Матрица:"
    sd = model.mdg.subdomains()[0]
    data = model.mdg.subdomain_data(sd)
    df = data['parameters']['mobility']['darcy_flux']
    x = sd.face_centers[0]
    y = sd.face_normals[0]
    A = sd.face_areas
    indx1 = x == 0
    df = df * y
    print(A[indx1])
    f1 = df[indx1] * ro / (mu*A[indx1])
    # print(f1)
    # print(df[indx1])
    Q1 = sum(abs(f1))

    "Трещины:"
    sd = model.mdg.subdomains()
    n = int(len(sd))
    Qtr = 0
    j = 0
    for i in range(1,n):
        x = np.round(sd[i].face_centers[0],  decimals=2)
        indx = x == 0
        if True in indx:
            j = j + 1
            data_tr = model.mdg.subdomain_data(sd[i])
            ff = data_tr['parameters']['mobility']['darcy_flux'] * ro / (mu)
            if indx[0] == True:
                Qtr = Qtr + abs(ff[0])
            else:
                Qtr = Qtr + abs(ff[-1])

    print("количество трещин на границе:" + str(j))
    k = (Q1+Qtr) * mu*L/(L*ro * (P1 - P2))
    return k
