import math

def ke(file1, outputfile, energyfile, ncrfile):
    planck = 6.6260693e-34  # planck's constant (j * s)
    na = 6.0221415e23       # avogadro constant (mol-1)
    bohr = 0.5291772108e-10  # bohr radius (m)
    j_cm = 1.98644561e-23 	# j to cm-1 convertion factor

    freq, freq_mecp, energy = ([] for i in range(3))
    grads, gradt = ([] for i in range(2))
    inputarr = []
    lines = file1.readlines()
    for line in lines:
        inputarr.append(line.strip())
    natoms = float(inputarr[1])
    numgrad = float(inputarr[2])
    h12 = float(inputarr[3])
    redmass = float(inputarr[4])
    eps = float(inputarr[5])
    e_min = float(inputarr[7])
    emecp = float(inputarr[9])
    ncalc = float(inputarr[11])
    for i in range(0, int(ncalc)):
        energy.append(float(inputarr[13+i]))
    for i in range(0, int(numgrad)):
        grads.append(float(inputarr[20+i]))
    for i in range(0, int(numgrad)):
        gradt.append(float(inputarr[28+i]))
    nfreq = 3*natoms-6
    nfreq_mecp = 3*natoms-7
    for i in range(0, int(nfreq)):
        freq.append(float(inputarr[36+i]))
    for i in range(0, int(nfreq_mecp)):
        freq_mecp.append(float(inputarr[40+i]))
    h12 = h12*11.9626/na
    redmass = redmass*(1.66054e-27)
    deltaf = 0
    for i in range(0, int(numgrad)):
        deltaf = deltaf + (grads[i]-gradt[i])*(grads[i]-gradt[i])
    coef = 62750.94717e-2*1000*4.184/na/bohr
    sqr = math.sqrt(deltaf)
    deltaf = sqr*coef
    outputfile.write('the results of k(e) are listed below\n')
    outputfile.write('\n')
    outputfile.write('e(internal),kcal/mol'+'----'+'k(e),s-1' +
                     '----'+'p(sh)'+'----'+'nmecp'+'----'+'rho(reactant)\n')
    outputfile.write('\n')

    for i in range(0, int(ncalc)):
        energy_mecp = abs(energy[i]-emecp)
        energy_min = abs(energy[i]-e_min)
        energyfile.write("\n")
        energyfile.write('internal energy = %f kcal/mol\n' % energy[i])
        energyfile.write(
            'available energy above mecp = %f kcal/mol\n' % energy_mecp)
        energyfile.write(
            'available energy above minima = %f kcal/mol\n' % energy_min)
        ncrfile.write("\n")
        ncrfile.write(
            'calculating the integrated density of states ncr(e) ...\n')
        energyfile.write("\n")
        nmecp = qsimp(energy_mecp, h12, redmass, deltaf,
                      freq_mecp, nfreq_mecp, eps, ncrfile)

        rho_reactant = rho_saddle_point(energy_min, freq, nfreq)

        rate = nmecp/(rho_reactant*(planck/j_cm))
        p_sh = psh(energy_mecp, h12, redmass, deltaf)

        energyfile.write('nmecp(e) = %f\n' % nmecp)
        energyfile.write('rho_r    = %f\n' % rho_reactant)
        energyfile.write('k(e)     = %f s-1\n' % rate)
        energyfile.write("\n")
        outputfile.write(str(energy[i])+"         "+str(rate)+"         "+str(
            p_sh)+"         " + str(nmecp)+"         "+str(rho_reactant)+'\n')


def qsimp(energy, h12, redmass, deltaf, freq, nfreq, eps, ncrfile):
    # print ("qsimpxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    # print(energy, h12, redmass, deltaf, freq, nfreq, eps)
    jmax = 200
    cal_cm = 349.755
    a = 0.0
    b = 0.0
    st = -1.e30
    os = -1.e30
    ost = 0.0
    for i in range(1, jmax+1):
        if i == 1:
            c = 0.5*(energy-a)*psh(a, h12, redmass, deltaf)
            d = rho_saddle_point((energy-a), freq, nfreq)
            a1 = c*d
            b1 = psh(energy, h12, redmass, deltaf) * \
                rho_saddle_point(b, freq, nfreq)
            st = a1+b1
            st = st*cal_cm
        else:
            it = 2**(i-2)
            tnm = it
            del1 = (energy-a)/tnm
            # print("print(it): ", it)
            x = a+0.5*del1
            sum = 0.
            # print("for-----------------")
            # print(x,del1)
            for j in range(1, it+1):
                # print ("--------------------: ", j)
                y0 = psh(x, h12, redmass, deltaf)
                # print(psh(x, h12, redmass, deltaf),"psh")
                # print(y0,x, h12, redmass, deltaf, "xxx")
                # print ("--------------------")
                # print (y0)
                # print ("--------------------")
                y1 = rho_saddle_point((energy-x), freq, nfreq)
                sum = sum + y0*y1*cal_cm
                x = x + del1
            # print(x,y0,sum,y1)
            # print("for-----------------")
            st = 0.5*(st+(energy-a)*sum/tnm)
        s1 = 4*st-ost
        s = s1/3
        qsimp = s
        ncrfile.write('Iter. ' + str(i) + '   NMECP(E)= ' +
                      str(s)+'   Grad = ' + str(abs(s-os))+'\n')
        if abs(s-os) < eps:
            return qsimp
        if (s == 0 and os == 0 and i > 6):
            return qsimp
        os = s
        ost = st
    # print ("qsimpxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    return qsimp
    pass


def psh(eh, h12, redmass, deltaf):
    h = 6.6260693E-34
    na = 6.0221415E23
    pi = 3.14
    # print(eh, h12, redmass, deltaf)
    if eh <= 0:
        psh = 0.0
    else:
        deltae = eh*1000*4.184/na
        # print ("print (deltae)")
        # print (deltae)
        # print ("print (deltae)")

        a = -2*pi*h12*h12
        b = math.sqrt(redmass/(2*deltae))
        c = h*deltaf
        p = math.exp(a*b/c)
        # print ("-------------")
        # print (a,b,c)
        # print ("-------------")
        # print
        psh = (1-p)*(1+p)
    # print (psh)
    return psh
    return 10


def rho_saddle_point(energy, freq, nfreq):
    pi = 3.14
    eng = energy * 349.755
    if (eng <= 1):
        rho_saddle_point = 0.0
    else:
        x1 = 1.e-20
        x2 = 1.e+3
        xacc = 1.e-12
        b2 = zriddr(x1, x2, xacc, eng, freq, nfreq)
        sum = 0.
        prod = 1.
        for i in range(0, int(nfreq)):
            temp1 = freq[i]/(expsafe(b2*freq[i])-1.)
            temp1 = temp1 * temp1
            sum = sum+expsafe(b2*freq[i])*temp1
            prod = prod / (1.0 - expsafe(-b2*freq[i]))
        rho_saddle_point = expsafe(b2*eng) * prod / math.sqrt(2.0 * pi * sum)
        # rho_saddle_point = 10
    # print (rho_saddle_point)
    return rho_saddle_point
    return 10


def zriddr(x1, x2, xacc, eng, freq, nfreq):
    maxit = 60
    unused = -1.11e30
    fl = funcrzo(x1, eng, freq, nfreq)
    fh = funcrzo(x2, eng, freq, nfreq)
    zriddr = 0
    if ((fl > 0 and fh < 0) or (fl < 0 and fh > 0)):
        xl = x1
        xl = x1
        xh = x2
        zriddr = unused
        for j in range(1, maxit+1):
            xn = xl+xh
            xm = 0.5*xn
            fm = funcrzo(xm, eng, freq, nfreq)
            s = math.sqrt(fm * fm - fl * fh)
            if (s == 0.0):
                return zriddr

            a = 1.0
            b = fl - fh
            xnew = xm+(xm-xl)*(math.copysign(a, b)*fm/s)
            if (abs(xnew-zriddr) <= xacc):
                return zriddr
            zriddr = xnew
            fnew = funcrzo(zriddr, eng, freq, nfreq)

            if (fnew == 0.):
                return zriddr

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
            else:
                print('never get here in zriddr')
            if (abs(xh-xl) <= xacc):
                return zriddr

        print('zriddr exceed maximum iterations')
    elif (fl == 0.):
        zriddr = x1
    elif (fh == 0.):
        zriddr = x2
    else:
        print('root must be bracketed in zriddr')
    return zriddr


def funcrzo(x, eng, freq, nfreq):
    sum = 0
    for i in range(0, int(nfreq)):
        y = x*freq[i]
        z = expsafe(y) - 1
        sum = sum + freq[i]/divsafe(z)
    return eng - sum


def expsafe(arg):
    if (arg > 700.):
        arg = 700.
    expsafe = math.exp(arg)
    return expsafe


def divsafe(arg):
    if (abs(arg) < 1.e-100):
        arg = 1.e-100
    divsafe = arg
    return divsafe