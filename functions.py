import math

emin, gamma, v1, v22, delta, rt = (0 for i in range(0, 6))


def eckart(t, e0, e1, vimag, qtun):
    e0, emax, pi = (0 for i in range(0, 3))

    pi = 3.14
    H = 9.5377077e-14
    R = 1.9872e-3
    CL = 2.99792e10

    EPSIL = 1.e-6
    freq = CL*vimag
    V1 = e0-e1
    V22 = math.sqrt(e1)+math.sqrt(e0)

    gamma = 2*math.sqrt(e0*e1)/(H*freq)

    if gamma <= 0.5:
        print("    E1*E0 IS TOO SMALL. THE PROGRAM WILL STOP")
        exit()
    delta = math.sqrt(math.pow(gamma, 2)-0.25)
    rt = R*t
    emax = e0 + 14 * rt
    n = 20
    if e0 <= e1:
        emin = 0
    else:
        emin = e0 - e1
    st = (emax - emin)/n
    alfa = gamma * math.sqrt(emax)/v22
    beta = gamma * math.sqrt(emax-v1)/v22
    x = 2*math.pi(alfa+beta)
    y = 2*math.pi(alfa-beta)
    z = 2*math.pi*delta

    c = (math.exp(x) + math.exp(-x))/2
    d = (math.exp(y) + math.exp(-y))/2
    e = (math.exp(z) + math.exp(-z))/2

    edg = (c-d)/2
    edg = edg/(c+e)
    edg = edg*math.exp((e0 - emax))/rt

    sum1 = 0
    sum1 = total(n-1, st, e0, sum1)
    sum1 = sum1 + edg
    x0 = st/2

    sum2 = 0

    check  = 0
    while check == 0:
        sum2 = total(n, x0, st, e0, sum2)
        sum2 = sum2 + sum1
        n = n*2
        st = (emax - emin)/n
        x0 = st/2
        if (abs((sum2*st-sum1*2*st)/3) >= EPSIL):
            sum1 = sum2
            check = 1
    qtun = sum2 * st
    return qtun


def partif(t, q, nira, dira, air, nirb, dirb, bir, qx):
    return qx


def total(n, x0, st, e0, sum):
    return sum
