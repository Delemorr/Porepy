
def ekv(z, L,ro,mu,P1,P2,model):
    sd = model.mdg.subdomains()[0]
    data = model.mdg.subdomain_data(sd)
    df = data['parameters']['mobility']['darcy_flux']

    x = sd.face_centers[0]
    y = sd.face_normals[0]

    indx1 = x == 0
    indx2 = x == L

    df = df * y

    f1 = df[indx1]
    f2 = df[indx2]
    Q1 = sum(f1) * ro / (mu*z)
    Q2 = sum(f2) * ro / (mu*z)
    e = (Q1-Q2)/(Q1+Q2 / 2)
    if e > 0.1:
        raise ValueError("Не совпадает расход на границах Q1 = " + str(Q1) + "Q2 = " + str(Q2))
    k = Q1 * mu / (ro * (P1 - P2))
    return k
