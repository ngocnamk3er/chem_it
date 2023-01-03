import math
import numpy

def main():
    Planck = 6.6260693E-34
    NA = 6.0221415E23
    BOHR = 0.5291772108E-10
    J_cm = 1.98644561E-23
    GradS = [], GradT = [], FREQ = [], FREQ_MECP = [], Energy = []
    f = open("input.txt", "r")
    f1 = open("output.txt", "w+")
    f2 = open("energy.txt", "w+")
    f3 = open("Ncr(E).txt", "w+")
    f.readline()
    Natoms = f.readline()
    NumGrad = f.readline()
    H12 = f.readline()
    redmass = f.readline()
    EPS = f.readline()
    f.readline()
    E_min = f.readline()
    f.readline()
    EMECP = f.readline()
    f.readline()
    Ncalc = f.readline()
    i = 1
    while (i <= Ncalc):
        Energy[i] = f.readline()
        i += 1
    i = 1
    i = 1
    while (i <= NumGrad):
        GradS[i] = f.readline()
        i += 1
    i = 1
    f.readline()
    while (i <= NumGrad):
        GradT[i] = f.readline()
        i += 1
    NFREQ = 3*Natoms - 6
    NFREQ_MECP = 3*Natoms - 7
    f.readline()
    i = 1
    while (i <= NFREQ):
        FREQ[i] = f.readline()
        i += 1
    f.readline()
    i = 1
    while (i <= NFREQ_MECP):
        FREQ_MECP[i] = f.readline()
        i += 1
    H12 = H12*11.9626/NA
    redmass = redmass*1.66054E-27
    deltaF = 0
    i = 1
    while (i <= NumGrad):
        deltaF = deltaF + (GradS[i]-GradT[i])*(GradS[i]-GradT[i])
        i += 1
    coef = 62750.94717E-2*1000*4.184/NA/BOHR
    sqr = math.sqrt(deltaF)
    deltaF = sqr*coef
    f1.write("1")
    f1.write("2")
    f1.write("3")
    f1.write("4")

    i = 1
    while (i <= Ncalc):
        Energy_MECP = abs(Energy[i]-EMECP)
        Energy_Min = abs(Energy[i]-E_min)
        f2.write('\n')
        f2.write('Internal energy = ', Energy[i], ' kcal/mol\n')
        f2.write('Available energy above MECP = ', Energy_MECP, ' kcal/mol\n')
        f2.write('Available energy above minima = ', Energy_Min, ' kcal/mol\n')
        f3.write('\n')
        f2.write('Calculating the integrated density of states Ncr(E) ...\n')
        f2.write("\n")
        NMECP = qsimp(Energy_MECP, H12, redmass, deltaF,
                      FREQ_MECP, NFREQ_MECP, EPS)
        RHO_REACTANT = RHO_SADDLE_POINT(Energy_Min, FREQ, NFREQ)
        rate = NMECP/(RHO_REACTANT*(Planck/J_cm))
        P_sh = Psh(Energy_MECP, H12, redmass, deltaF)
        f2.write('NMECP(E) = ', NMECP, '\n')
        f2.write('RHO_R    = ', RHO_REACTANT, '\n')
        f2.write('k(E)     = ', rate, ' s-1', '\n')
        f2.write('\n')
        f1.write(Energy[i], rate, P_sh, NMECP, RHO_REACTANT, '\n')
        f2.write('Calculation finished! Open output file to get results.')
        i += 1

    f1.write("\n")
    f1.write("\n")
    f1.write('Calculation finished!')
    f.close()
    f1.close()
    f2.close()
    f3.close()


def qsimp(ENERGY, H12, RedMass, DeltaF, FREQ, NFREQ, EPS):
    a = 0
    b = 0
    ost = 0
    Cal_cm = 349.755
    st = -1.e30
    os = -1.e30
    JMAX = 200
    maxfreq = 1000
    i = 1
    while (i <= JMAX):
        if (i == 1):
            c = 0.5*(ENERGY-a)*Psh(a, H12, RedMass, DeltaF)
            d = RHO_SADDLE_POINT((ENERGY-a), FREQ, NFREQ)
            a1 = c*d
            b1 = Psh(ENERGY, H12, RedMass, DeltaF) * \
                RHO_SADDLE_POINT(b, FREQ, NFREQ)
            st = a1+b1
            st = st*Cal_cm
            sum = 0
        else:
            it = 2**(i-2)
            tnm = it
            del1 = (ENERGY-a)/tnm
            x = a+0.5*del1
            sum = 0
            j = 1
            while (j <= it):
                y0 = Psh(x, H12, RedMass, DeltaF)
                y1 = RHO_SADDLE_POINT((ENERGY-x), FREQ, NFREQ)
                sum = sum + y0*y1*Cal_cm
                x = x + del1
                j += 1
            st = 0.5*(st+(ENERGY-a)*sum/tnm)
        s1 = 4*st-ost
        s = s1/3
        qsimp = s
        print('Iter. ', i, '   NMECP(E)= ', s, '   Grad = ', abs(s-os))
        if (abs(s-os) < EPS):
            return
        if (s == 0 and os == 0 and i > 6):
            return
        os = s
        ost = st
        i += 1


def Psh(Eh, H12, redmass, deltaF):
    h = 6.6260693E-34
    NA = 6.0221415E23
    pi = 4*numpy.arctan(1)
    Psh
    if (Eh <= 0):
        Psh = 0
    else:
        deltaE = Eh*1000*4.184/NA
        A = -2*pi*H12*H12
        B = math.sqrt(redmass/(2*deltaE))
        C = h*deltaF
        P = math.exp(A*B/C)
        Psh = (1-P)*(1+P)
    return Psh


def RHO_SADDLE_POINT(ENERGY, FREQ, NFREQ):
    PI = numpy.arccos(-1)
    eng = ENERGY * 349.755
    RHO_SADDLE_POINT
    if (eng <= 1):
        RHO_SADDLE_POINT = 0
    else:
        x1 = 1.E-20
        x2 = 1.E+3
        xacc = 1.E-12
        b2 = ZRIDDR(x1, x2, xacc, eng, FREQ, NFREQ)
        sum = 0
        prod = 1
        i = 1
        while (i <= NFREQ):
            temp1 = FREQ[i]/(expsafe(b2*FREQ[i])-1)
            temp1 = temp1*temp1
            sum = sum+expsafe(b2*FREQ[i])*temp1
            prod = prod / (1.0 - expsafe(-b2*FREQ[i]))
        RHO_SADDLE_POINT = expsafe(b2*eng) * prod / math.sqrt(2.0 * PI * sum)
    return RHO_SADDLE_POINT


def ZRIDDR(x1, x2, xacc, eng, FREQ, NFREQ):
    MAXIT = 60
    UNUSED = -1.11E30
    maxfreq = 1000
    fl = funcrzo(x1, eng, FREQ, NFREQ)
    fh = funcrzo(x2, eng, FREQ, NFREQ)
    if ((fl > 0 and fh < 0) or (fl < 0 and fh > 0)):
        xl = x1
        xh = x2
        zriddr = UNUSED
        j = 1
        while (j <= MAXIT):
            xn = xl+xh
            xm = 0.5*xn
            fm = funcrzo(xm, eng, FREQ, NFREQ)
            s = math.sqrt(fm * fm - fl * fh)
            if (s == 0):
                return
            a = 1.0
            b = fl - fh
            xnew = xm+(xm-xl)*(numpy.sign(a, b)*fm/s)
            if (abs(xnew-zriddr) <= xacc):
                return
            zriddr = xnew
            fnew = funcrzo(zriddr, eng, FREQ, NFREQ)
            if (fnew == 0):
                return
            if (numpy.sign(fm, fnew) != fm):
                xl = xm
                fl = fm
                xh = zriddr
                fh = fnew
            elif (numpy.sign(fl, fnew) != fl):
                xh = zriddr
                fh = fnew
            elif (numpy.sign(fh, fnew) != fh):
                xl = zriddr
                fl = fnew
            else:
                print('never get here in zriddr')
            if (abs(xh-xl) <= xacc):
                return
            j += 1
            print('zriddr exceed maximum iterations')
    elif (fl == 0):
        zriddr = x1
    elif (fh == 0):
        zriddr = x2
    else:
        print('root must be bracketed in zriddr')


def funcrzo(x, eng, FREQ, NFREQ):
    sum = 0
    i = 1
    while (i <= NFREQ):
        y = x*FREQ(i)
        z = expsafe(y)-1.0
        sum = sum+FREQ(i)/divsafe(z)
        i += 1
    return eng-sum


def expsafe(arg):
    if (arg > 700):
        arg = 700
    return math.exp(arg)


def divsafe(arg):
    if (abs(arg) < 1.E-100):
        arg = 1.E-100
    return arg


if __name__ == "__main__":
    main()
