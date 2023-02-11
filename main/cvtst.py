import math
import re
import csv
import numpy


def cvtst(f, f3, f4):
    global A, B, C, SN, SF, F, NF, NV
    title = f.readline()
    header = title.split('\n')[0].strip()
    writer = csv.writer(f3)
    writer1 = csv.writer(f4)
        # write the header
    writer.writerow([header])
    writer1.writerow(['Temperatures'])
    NP = int(re.findall(r'\d+', f.readline())[0])
    NR = int(re.findall(r'\d+', f.readline())[0])

    f.readline()
    M = re.findall(r'\d+', f.readline())
    NN = 0
    for x in M:
        NN += int(x)

    NT = int(re.findall(r'\d+', f.readline())[0])

    stringRead = f.readline().split(';')
    ZPE = (stringRead[len(stringRead) - 1].split('\n')[0]).split()[0]

    stringRead = f.readline().split(';')
    TUN = (stringRead[len(stringRead) - 1].split('\n')[0]).split()[0]

    f.readline()
    RES = "bar"
    stringRead = f.readline().split(';')
    ASK = (stringRead[len(stringRead) - 1].split('\n')[0]).split()[0]

    if ASK == "yes":
        RES = "atm"
    stringRead = f.readline().split(';')
    ASK = (stringRead[len(stringRead) - 1].split('\n')[0]).split()[0]

    if ASK == "yes":
        RES = "bar"

    FAC = 0.46499834E0
    R = 1.9872E0
    RCC1 = 2.08367541593E10
    RCC2 = 0
    if RES == "atm":
        RCC2 = 1.3626055E-22
    if RES == "bar":
        RCC2 = 1.38066E-22
    RCC = RCC1 * RCC2 ** (NN - 1)

    C1 = 1.438767553E0
    # !--------------------------- Read Temperatures -------------------------C
    f.readline()
    TEMP = []
    for x in range(0, NT):
        stringRead = re.findall(r'\d+', f.readline())
        TEMP.append(int(stringRead[0]))

    for x in range(0, 10):
        f.readline()
    GEO = []
    E = []
    DGT = []
    EZ = []
    RCV = []
    for L in range(0, NP):
        f.readline()
        stringRead = f.readline()
        stringRead = stringRead.split()
        GEO.append(float(stringRead[0]))
        E.append(float(stringRead[1]))
        SPECIES = stringRead[2].split('\n')
        writer.writerow([GEO[L],E[L],SPECIES[0]])
        E[L] = E[L] * 627.50595E0
        BE = 0
        AE = 0
        WE = 0
        WEX = 0
        WT = 0
        G0 = 0
        G1 = 0
        EE = 0
        if SPECIES[0] == 'diatomic':
            stringRead = f.readline().split()
            WT = float(stringRead[0])
            G0 = float(stringRead[1])
            G1 = float(stringRead[2])
            EE = float(stringRead[3])
            stringRead = f.readline().split()
            BE = float(stringRead[0])
            AE = float(stringRead[1])
            WE = float(stringRead[2])
            WEX = float(stringRead[3])
        else:
            stringRead = f.readline().split()
            WT = float(stringRead[0])
            G0 = float(stringRead[1])
            G1 = float(stringRead[2])
            EE = float(stringRead[3])

        stringRead = f.readline().split()
        A = float(stringRead[0])
        B = float(stringRead[1])
        C = float(stringRead[2])
        SN = float(stringRead[3])
        SF = float(stringRead[4])
        F = float(stringRead[5])
        NF = float(stringRead[6])
        NV = int(stringRead[7])
        A = A * FAC
        B = B * FAC
        C = C * FAC
        F = F * FAC
        W = []
        if NV != 0:
            for j in range(0, int(math.ceil(NV / 3))):
                stringRead = f.readline().split()
                for i in range(0, len(stringRead)):
                    W.append(float(stringRead[i]))
        if TUN == "yes":
            stringRead = f.readline().split()
            RCV.append(float(stringRead))

        if TUN == 'yes':
            f.readline().split()
        GIBBS(NT, TEMP, L, SPECIES, WT, G0, G1, EE, BE, AE, WE, WEX, SN, SF, A, B, C, F, NF, W, NV, DGT, E[L], RES, ZPE,
              EZ)
    # ------------- List "Delta G" over range of temperatures - --------------
    MNR = float(M[NR - 1])
    if NR > 1:
        for j in range(0, NT):
            DGT[NR - 1][j] = DGT[NR - 1][j] * MNR
        EZ[NR - 1] = MNR * EZ[NR - 1]
        for i in range(0, NR - 1):
            EZ[NR - 1] = EZ[NR - 1] * float(M[i]) * EZ[i]
            for j in range(0, NT):
                DGT[NR - 1][j] = DGT[NR - 1][j] + float(M[i]) * DGT[i][j]

    writer.writerow(['Geometrical Parameters:'])
    i = NR - 1
    while i < NP:
        writer.writerow([str(GEO[i])])
        i += 1
    # dG/dR\
    GDAT = []
    TV = []
    EZT = []
    GEOM = []
    dGdR = []
    for i in range(0, NT):
        writer.writerow(['--------- dG at ' + str(TEMP[i]) + ' Degrees Kelvin:'])
        writer.writerow(['R [Angstrom]','G [kcals/mol]','dG/dR'])
        GDAT.append(-10)
        GEOM.append(0)
        EZT.append(0)
        j = NR - 1
        GG = []
        while j < NP:
            DD = DGT[j][i] - DGT[NR - 1][i]
            GG.append(DD)
            if DD > GDAT[i]:
                GDAT[i] = DD
                if TUN == "yes":
                    TV.append(RCV[j])
                else:
                    TV.append(0)
                EZT[i] = EZ[j] - EZ[NR - 1]
                GEOM[i] = GEO[j]
            j += 1
        j = NR - 1
        while j < NP:
            if j == NR - 1 or j == NP - 1:
                writer.writerow([str(GEO[j]),str(GG[j - 1]),'NAN'])
            else:
                dGdRA = (GG[j] - GG[j - 1]) / (GEO[j + 1] - GEO[j])
                dGdRB = (GG[j - 1] - GG[j - 1 - 1]) / (GEO[j] - GEO[j - 1])
                # dGdR.append((dGdRA + dGdRB) / 2)
                writer.writerow([str(GEO[j]),str(GG[j - 1]),str((dGdRA + dGdRB) / 2)])
            j += 1
    writer.writerow(['Greatest Gibbs Free-Energy Barriers [kcals/mol]'])
    writer.writerow(['at the temperatures specified:'])
    writer.writerow(['T [Kelvin]','R [Angstrom]','G [kcals/mol]'])
    for k in range(0, NT):
        writer.writerow([str(TEMP[k]),str(GEOM[k]),str(GDAT[k])])

    # Do quadratic least squares fit to GDAT as function of temperature
    # AA = []
    # BB = []
    # for i in range(0, NT):
    #     AA.append(0)
    #     BB.append(0)
    #     temp = [0, 0, 0]
    #     temp[0] = TEMP[i] * TEMP[i]
    #     temp[1] = TEMP[i]
    #     temp[2] = 1
    #     AA[i] = temp
    #     BB[i] = GDAT[i]
    # print(AA)
    # print(BB)
    X = [-0.116962E-02 / 1000, 0.336330E+02 / 1000, -0.695957E+04 / 1000]
    ISING = 0
    WW = []
    # comment code coding fix bug 
    # HOUSEFIT(AA, NT, 3, BB, NT, 1, X, WW, ISING)
    writer.writerow(['Same Barriers from quadratic fit:'])
    writer.writerow(['T [Kelvin]','G [kcals/mol]'])
    for i in range(0, NT):
        PRED = X[0] * TEMP[i] * TEMP[i] + X[1] * TEMP[i] + X[2]
        writer.writerow([str(TEMP[i]),str(PRED)])
    # Convert to Calories:
    X[0] = X[0] * 1000
    X[1] = X[1] * 1000
    X[2] = X[2] * 1000
    writer.writerow(['dG(FIT) = ' + str(X[0]) + ' T^2 + ' + str(X[1]) + 'T + ' + str(X[2]) + '  [Calories]'])
    # Print out Rate Constant equation and k @ temps
    for i in range(0, NT):
        writer1.writerow([str(TEMP[i])])
    writer1.writerow([])
    writer.writerow(['Rate Constant:'])
    writer.writerow(['k = ' + str(RCC) + ' T^' + str(NN) + ' exp(-dG/RT)  [cc/molecule]^' + str(NN - 1) + '/second'])
    writer.writerow(['','dG = Raw Data','dG = Quad. Fit'])
    for i in range(0, NT):
        PRED = X[0] * TEMP[i] * TEMP[i] + X[1] * TEMP[i] + X[2]
        FRC = RCC * (TEMP[i] ** NN) * numpy.exp(-PRED / (R * TEMP[i]))
        RRC = RCC * (TEMP[i] ** NN) * numpy.exp(-GDAT[i] * 1000 / (R * TEMP[i]))
        writer.writerow([str(TEMP[i]),str(RRC),str(FRC)])
        writer1.writerow([str(RRC)])
    if TUN == "yes":
        writer.writerows([
            ['Tunneling factor = [1 +(hcw/kT)^2(1+RT/Eo)/24])'],
            ['where w(cm^-1) = negative frequency (along path))'],
            ['Note: This factor only valid for small barriers)'],
            ['with a width of less than about 1.8 Angstroms.)'],
            [],
            ['Rates with Tunneling Correction:)']
        ])
        for i in range(0, NT):
            PRED = X[0] * TEMP[i] * TEMP[i] + X[1] * TEMP[i] + X[2]
            FRC = RCC * (TEMP[i] ** NN) * numpy.exp(-PRED / (R * TEMP[i]))
            RRC = RCC * (TEMP[i] ** NN) * numpy.exp(-GDAT[i] * 1000 / (R * TEMP[i]))
            TF = (C1 * TV[i] / TEMP[i]) ** 2.0
            TF = TF * (1.0 + R * TEMP[i] / (EZT[i] * 1000.0)) / 24.0
            TF = TF + 1.0
            RRC = TF * RRC
            FRC = TF * FRC
            writer.writerow([str(TEMP[i]),str(RRC),str(FRC)])
            writer1.writerow([str(RRC)])

def GIBBS(NT, TEMP, L, SPECIES, WT, G0, G1, EE, BE, AE, WE, WEX, SN, SF, A, B, C, F, NF, W, NV, DGT, E, RES, ZPE, EZ):
    C1 = 1.438767553
    C2 = 40.27512295
    R = 1.9872
    PI = numpy.arccos(-1.0)
    TC = 0
    if RES == "atm":
        TC = -3.664854875E0
    if RES == "bar":
        TC = -3.651691889E0
    dgt1 = []
    for II in range(0, NT):
        CELECT = 0
        T = TEMP[II]
        if EE == 0:
            CELECT = numpy.log(G0)
        if EE != 0:
            CELECT = numpy.log(G0 + G1 * numpy.exp(C1 * EE / T))
        CTRANS = 1.5 * numpy.log(WT) + 2.5 * numpy.log(T) + TC
        G = 0
        E0 = 0
        # Monatomic Species
        if SPECIES[0] == 'monatomic':
            G = (CELECT + CTRANS) * R * T
            E0 = 0
        # Diatomic Species
        if SPECIES[0] == 'diatomic':
            Y = C1 * (BE - AE / 2) / T
            U = C1 * (WE - 2 * WEX) / T
            XE = WEX / WE
            D = AE / BE
            GG = BE / WE
            if U > 80:
                EU = numpy.exp(80)
                U1 = 0
            else:
                EU = numpy.exp(U) - 1
                U1 = U / EU / EU
            if ZPE != "yes":
                E0 = 0
            else:
                E0 = -U / 2
            G = CTRANS - numpy.log(SN * Y) + Y / 3 + Y ** 2 / 90 - numpy.log(1 - numpy.exp(-U)) + CELECT + E0 + (
                    2 * XE * U1) + D / EU + 8 * GG / U
            G = G * R * T
            E0 = E0 * R * T
        # Polyatomic Species
        CVIB = 0
        for i in range(0, NV):
            if ZPE == "yes":
                E0 = E0 - (C1 * W[i] / T) / 2
            CVIB = CVIB - numpy.log(1 - numpy.exp(-C1 * W[i] / T))
        CVIB = CVIB + E0
        # linear Species
        if SPECIES[0] == "linear":
            Y = C2 / A / T
            if F <= 0:
                FR = 0
            else:
                YY = C2 / F / T
                if NF - 2 < 0:
                    FR = 0.5 * math.log(math.pow(SF, 2) * YY / PI)
                elif NF - 2 == 0:
                    FR = math.log(SF * YY)
                else:
                    FR = 0.5 * math.log(math.pow(SF, 2) * math.pow(YY, 3) / PI)
            G = CTRANS - FR - math.log(SN * Y) + Y / 3 + math.pow(Y, 2) / 90 + CVIB + CELECT
            G = G * R * T
            E0 = E0 * R * T
            # nonlinear Species
        if SPECIES[0] == "nonlinear":
            Y = (C2 / A / T) * (C2 / B / T) * (C2 / C / T)
            if F <= 0:
                FR = 0
            else:
                YY = C2 / F / T
                if NF - 2 < 0:
                    FR = 0.5 * math.log(math.pow(SF, 2) * YY / PI)
                elif NF - 2 == 0:
                    FR = math.log(SF * YY)
                else:
                    FR = 0.5 * math.log(math.pow(SF, 2) * math.pow(YY, 3) / PI)
            G = CTRANS - FR - 0.5 * math.log(math.pow(SN, 2) * Y / PI) + CVIB + CELECT
            G = G * R * T
            E0 = E0 * R * T

        dgt1.append(E - G / 1000)
        EZ.append(E - E0 / 1000)
    DGT.append(dgt1)


def HOUSEFIT(A, M, N, B, K, L, X, W, ISING):
    # (AA, NT, 3, BB, NT, 1, X, WW, ISING)
    HOUSETRANS(A, M, N, B, K, L, W, ISING)
    BACKSUB(A, M, N, B, K, L, X)
    pass


def HOUSETRANS(A, M, N, B, KB, LB, W, ISING):
    # pass
    for K in range(0, N):
        MX = idamax(M - K + 1, A[K][K], 1) + K - 1
        RMS = 0.0
        I = K
        while I <= M:
            W[I] = A[I][K] / abs(A[MX][K])
            RMS = RMS + W[I] * W[I]
            I += 1
        RMS = math.sqrt(RMS)
        BK = 1 / (RMS * (RMS + abs(W[K])))
        ZZ = W[K]
        W[K] = W[K] + numpy.sign(RMS, ZZ)
        J = 1
        while J <= N:
            S = ddot(M - K + 1, W[K], 1, A[K][J], 1)
            S = BK * S
            daxpy(M - K + 1, -S, W[K], 1, A[K][J], 1)
            J += 1
        J = 1
        while J <= LB:
            S = ddot(M - K + 1, W[K], 1, B[K][J], 1)
            S = BK * S
            daxpy(M - K + 1, -S, W[K], 1, B[K][J], 1)
            J += 1


def BACKSUB(A, M, N, B, KB, LB, X):
    dcopy(KB * LB, B, 1, X, 1)
    L = 1
    while L <= LB:
        J = N
        while J >= 2:
            X[J] = X[J] / A[J]
            daxpy(J - 1, -X[J], A[1][J], 1, X[1], 1)
            J -= 1
        X[1] = X[1] / A[1][1]


def dswap(n, dx, incx, dy, incy):
    if incx == 1 and incy == 1:
        for i in range(0, n):
            s = dx[i]
            dx[i] = dy[i]
            dy[i] = s
    else:
        ix = 0
        iy = 0
        for i in range(0, n):
            s = dx[ix]
            dx[ix] = dy[iy]
            dy[iy] = s
            ix = ix + incx
            iy = iy + incy


def ddot(n, dx, incx, dy, incy):
    s = 0
    if incx == 1 and incy == 1:
        for i in range(0, n):
            s += dx[i] * dy[i]
    else:
        ix = 0
        iy = 0
        for i in range(0, n):
            s += dx[ix] * dy[iy]
            ix = ix + incx
            iy = iy + incy
    return s


def daxpy(n, da, dx, incx, dy, incy):
    if incx == 1 and incy == 1:
        for i in range(0, n):
            dy[i] = dy[i] + da * dx[i]
    else:
        ix = 0
        iy = 0
        for i in range(0, n):
            dy[iy] = dy[iy] + da * dx[ix]
            iy = iy + incy
            ix = ix + incx


def dscal(n, da, dx, incx):
    if incx == 1:
        for i in range(0, n):
            dx[i] = da * dx[i]
    else:
        ix = 0
        for i in range(0, n):
            dx[i] = da * dx[i]
            ix = ix + incx


def idamax(n, dx, incx):
    ii = 1
    xmax = abs(dx[1])
    if incx == 1:
        for i in range(0, n):
            if abs(dx[i]) > xmax:
                xmax = abs(dx[i])
                ii = i
    else:
        ix = 0
        for i in range(0, n):
            ix = ix + incx
            if abs(dx[ix]) > xmax:
                xmax = abs(dx[ix])
                ii = i
    return ii


def dcopy(n, dx, incx, dy, incy):
    if incx == 1 and incy == 1:
        for i in range(0, n):
            dy[i] = dx[i]
    else:
        ix = 1
        iy = 1
        for i in range(0, n):
            dy[iy] = dx[ix]
            iy = iy + incy
            ix = ix + incx


def drot(n, dx, incx, dy, incy, c, s):
    if abs(c) > abs(s):
        c1 = c
        s1 = s
        s2 = -s1 / c1
        c2 = c1 - s2 * s1
    else:
        c1 = s
        s1 = c
        dswap(n, dx, incx, dy, incy)
        s2 = s1 / c1
        c2 = -s2 * s1 - c1
    dscal(n, c1, dx, incx)
    daxpy(n, s1, dy, incy, dx, incx)
    dscal(n, c2, dy, incy)
    daxpy(n, s2, dx, incx, dy, incy)
