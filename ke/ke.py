import math
import numpy

def main():
    Planck = 6.6260693E-34
    NA = 6.0221415E23
    BOHR = 0.5291772108E-10
    J_cm = 1.98644561E-23
    GradS = [] 
    GradT = []
    FREQ = []
    FREQ_MECP = [] 
    Energy = []
    f = open("input.txt", "r")
    f1 = open("output.txt", "w+")
    f2 = open("energy.txt", "w+")
    f3 = open("Ncr(E).txt", "w+")
    f.readline()
    Natoms = int(f.readline())
    NumGrad = int(f.readline())
    H12 = float(f.readline())
    redmass = float(f.readline())
    EPS = float(f.readline())
    f.readline()
    E_min = float(f.readline())
    f.readline()
    EMECP = float(f.readline())
    f.readline()
    Ncalc = int(f.readline())
    f.readline()
    i = 1
    while (i <= Ncalc):
        Energy.append(float(f.readline()))
        i += 1
    i = 1
    f.readline()
    while (i <= NumGrad):
        GradS.append(float(f.readline()))
        i += 1
    i = 1
    f.readline()
    while (i <= NumGrad):
        GradT.append(float(f.readline()))
        i += 1
    NFREQ = 3*Natoms - 6
    NFREQ_MECP = 3*Natoms - 7
    f.readline()
    i = 1
    while (i <= NFREQ):
        FREQ.append(float(f.readline()))
        i += 1
    f.readline()
    i = 1
    while (i <= NFREQ_MECP):
        FREQ_MECP.append(float(f.readline()))
        i += 1
    H12 = H12*11.9626/NA
    redmass = redmass*1.66054E-27
    deltaF = 0
    i = 0
    while (i < NumGrad):
        deltaF = deltaF + (GradS[i]-GradT[i])*(GradS[i]-GradT[i])
        i += 1
    coef = 62750.94717E-2*1000*4.184/NA/BOHR
    sqr = math.sqrt(deltaF)
    deltaF = sqr*coef
    i = 0
    while (i < Ncalc):
        Energy_MECP = abs(Energy[i]-EMECP)
        Energy_Min = abs(Energy[i]-E_min)
        f2.write('\n')
        f2.write('Internal energy = ')
        f2.write(str(Energy[i]))
        f2.write(' kcal/mol')
        f2.write('\n')
        f2.write('Available energy above MECP = ')
        f2.write(str(Energy_MECP))
        f2.write(' kcal/mol')
        f2.write('\n\n')
        f2.write('Available energy above minima = ')
        f2.write(str(Energy_Min))
        f2.write(' kcal/mol')
        f3.write('\n\n')
        f3.write('\nCalculating the integrated density of states Ncr(E) ...\n')
        f2.write("\n")
        P_sh = Psh(Energy_MECP, H12, redmass, deltaF)
        NMECP = qsimp(Energy_MECP, H12, redmass, deltaF,FREQ_MECP, NFREQ_MECP, EPS, f3)
        RHO_REACTANT = RHO_SADDLE_POINT(Energy_Min, FREQ, NFREQ)
        rate = NMECP/(RHO_REACTANT*(Planck/J_cm))
        f2.write('NMECP(E) = ')
        f2.write(str(NMECP))
        f2.write('\n')
        f2.write('RHO_R    = ')
        f2.write(str(RHO_REACTANT))
        f2.write('\n')
        f2.write('k(E)     = ')
        f2.write(str(rate))
        f2.write(' s-1')
        
        f2.write('\n')
        f1.write(str(Energy[i]))
        f1.write('\n')
        f1.write(str(rate))
        f1.write('\n')
        f1.write(str(P_sh))
        f1.write('\n')
        f1.write(str(NMECP))
        f1.write('\n')
        f1.write(str(RHO_REACTANT))
        f1.write('\n\n')
        i += 1

    f1.write('Calculation finished!')
    f2.write('\nCalculation finished! Open output file to get results.\n')
    f.close()
    f1.close()
    f2.close()
    f3.close()


def qsimp(ENERGY, H12, RedMass, DeltaF, FREQ, NFREQ, EPS, f3):
    a = 0
    b = 0
    ost = 0
    Cal_cm = 349.755
    st = -1.e30
    os = -1.e30
    JMAX = 200
    maxfreq = 1000
    d = 0
    b1 = 0
    i = 1
    while (i <= JMAX):
        if (i == 1):
            c = 0.5*(ENERGY-a)*Psh(a, H12, RedMass, DeltaF)
            d = RHO_SADDLE_POINT(ENERGY-a, FREQ, NFREQ)
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
        f3.write('Iter. ')
        f3.write(str(i))
        f3.write('NMECP(E)= ')
        f3.write(str(s))
        f3.write('Grad = ')
        f3.write(str(abs(s-os)))
        f3.write('\n')
        
        if (abs(s-os) < EPS):
            break
        if (s == 0 and os == 0 and i > 6):
            break
        os = s
        ost = st
        i += 1
    return qsimp

def Psh(Eh, H12, redmass, deltaF):
    h = 6.6260693E-34
    NA = 6.0221415E23
    pi = 4*numpy.arctan(1)
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
    RHO_SADDLE_POINT = 0
    if (eng <= 1):
        RHO_SADDLE_POINT = 0
    else:
        x1 = 1.E-20
        x2 = 1.E+3
        xacc = 1.E-12
        b2 = ZRIDDR(x1, x2, xacc, eng, FREQ, NFREQ)
        sum = 0
        prod = 1
        i = 0
        while (i < NFREQ):
            sum = sum + expsafe(b2*FREQ[i])*(FREQ[i]/(expsafe(b2*FREQ[i])-1))
            prod = prod / (1 - expsafe(-b2*FREQ[i]))
            i += 1
        # RHO_SADDLE_POINT = (expsafe(b2*eng) * prod) / math.sqrt(2.0 * PI * sum)
        RHO_SADDLE_POINT = 1
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
                break
            a = 1.0
            b = fl - fh
            xnew = xm+(xm-xl)*(math.copysign(a, b)*fm/s)
            if (abs(xnew-zriddr) <= xacc):
                break
            zriddr = xnew
            fnew = funcrzo(zriddr, eng, FREQ, NFREQ)
            if (fnew == 0):
                break
            if (math.copysign(fm, fnew) != fm):
                xl = xm
                fl = fm
                xh = zriddr
                fh = fnew
            elif (math.copysign(fl, fnew) != fl):
                xh = zriddr
                fh = fnew
            elif (math.copysign(fh, fnew) != fh):
                xl = zriddr
                fl = fnew
            if (abs(xh-xl) <= xacc):
                break
            j += 1
    elif (fl == 0):
        zriddr = x1
    elif (fh == 0):
        zriddr = x2
    return zriddr

def funcrzo(x, eng, FREQ, NFREQ):
    sum = 0
    i = 0
    while (i < NFREQ):
        y = x*FREQ[i]
        z = expsafe(y)-1
        sum = sum+FREQ[i]/divsafe(z)
        i += 1
    return eng-sum


def expsafe(arg):
    if (arg > 700):
        arg = 700
    arg = numpy.exp(arg)
    return arg


def divsafe(arg):
    if (abs(arg) < 1.E-100):
        arg = 1.E-100
    return arg


if __name__ == "__main__":
    main()
