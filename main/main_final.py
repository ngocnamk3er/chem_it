import csv
import math
import ke
import cvtst
import inputfile
import os
emin, gamma, v1, v22, delta, rt = (0 for i in range(0, 6)) # /eint/
fva, fvb, fvx, air, dira = ([] for i in range(5)) #/FVS/
bir, dirb, xir, dirx, fvxir, virx =  ([] for i in range(6)) #  /HRA/
na, nb, nx, nirx = ( 0 for i in range(4))
h,r,an,rt,wa,wb,ea = ( 0 for i in range(7))
ai1,ai2,ai3,bi1,bi2,bi3,xi1,xi2,xi3 = ( 0 for i in range(9))

def partif(t, q, nira, dira, air, nirb, dirb, bir, qx):
    global fvx,fva,fvb,na,nb,nx,nirx,fvxir,dirx,virx
    global h,r,an,rt,wa,wb,ea
    global ai1,ai2,ai3,bi1,bi2,bi3,xi1,xi2,xi3,xir
    cir = 0
    cf3 = 0
    ct=3.11985E-04
    cr=2.48320E-02
    crp=6.93580E-03
    qta=ct*(wa*t)**1.5
    qtb=ct*(wb*t)**1.5
    if qtb<=0 :
        qtb = 1.0
    qtx = ct * ((wa + wb)*t)**1.5
    if xi3<=0:
        qtx = cr * math.sqrt(xi1*xi2)*t
    else:
        qrx = crp * math.sqrt(xi1*xi2*xi3)*t**1.5
    qva = 1.0
    qvb = 1.0
    qrb = 1.0
    qvx = 1.0

    for i in range(0,int(na)):
        qva = qva / (1.0 - math.exp(-h*fva[i]/rt))
    if(ai1<=0):
        qra = cr*math.sqrt(ai3*ai2)*t 
    else:
        qra = crp * math.sqrt(ai1*ai2*ai3)*t**1.5
    if nb != 0:
        for i in range(1, nb+1):
            qvb = qvb/(1.0 - math.exp(-h*fvb[i]/rt))
        if bi3<= 0:
            qrb = cr * math.sqrt(bi1 * bi2) * t
        else:  
            qrb =  crp * math.sqrt(bi1*bi2*bi3)*t**1.5
    for i in range(0,int(nx)):
        qvx = qvx/(1.0-math.exp(-h*fvx[i]/rt))
    q = qtx/(qta*qtb)*qrx/(qra*qrb)*qvx/(qva*qvb)
    cir = 0.279321
    if(nirx != 0):
        for i in range(0,nirx):
            qx = cir * math.sqrt(xir[i]*t)/dirx[i]
            if fvxir[i] <= 0:
                virx[i] = 0
                cf = 1.0
            else :
                virx[i] = 1.0215E-04*xir[i]*fvxir[i]**2/dirx
                [i]**2
                cf1 = (1.0-math.exp(-rt/(virx[i])))**1.2
                cf2 = math.exp(-1.2*rt/virx[i])
                cf3 = (qx*(1.0-math.exp(-dirx[i]/qx*math.sqrt(math.pi*virx[i]/rt))))
                cf = cf1 + (cf2/cf3)
            q = q*qx*cf
    
    if(nira != 0):
        for i in range(0,nira):
            qa = (cir*(math.sqrt(air[i]*t))/dira[i])
            q = q/qa
    
    if nirb != 0:
        for i in range(0,nirb):
            qb = (cir*(math.sqrt(bir[i]*t))/dirb[i])
            q = q / qb
    
    return [q,qx]

def eckart(t, e0, e1, vimag, qtun):
        global emin,gamma,v1,v22,delta,rt  
        emax = 0
        h = 9.5377077e-14
        r = 1.9872e-3
        cl = 2.99792e10

        EPSIL = 1.e-6
        freq =cl*vimag
        v1 = e0-e1
        v22 = math.sqrt(e1)+math.sqrt(e0)
        gamma = 2*math.sqrt(e0*e1)/(h*freq)
        if gamma < 0.5:
            exit()
        delta = math.sqrt(gamma**2-0.25)
        rt = r*t
        emax = e0 + 14 * rt
        n = 20
        if e0 <= e1:
            emin = 0
        else:
            emin = e0 - e1
        st = (emax - emin)/n
        alfa = gamma * math.sqrt(emax)/v22
        beta = gamma * math.sqrt(emax-v1)/v22
        x = 2*math.pi*(alfa+beta)
        y = 2*math.pi*(alfa-beta)
        z = 2*math.pi*delta
        # c = (np.exp(x) + np.exp(-x))/2
        # d = (np.exp(y) + np.exp(-y))/2
        # e = (np.exp(z) + np.exp(-z))/2
        c = 2
        d = 2
        e = 2
        edg = (c-d)/2
        edg = edg/(c+e)
        edg = edg*math.exp((e0 - emax)/rt)/rt

        sum1 = 0
        check = 0
        while check == 0:
            sum1 = total(n-1,st,  st, e0, sum1)
            sum1 = sum1 + edg
            x0 = st/2
            sum2 = 0
            sum2 = total(n, x0, st, e0, sum2)
            sum2 = sum2 + sum1
            n = n*2
            st = (emax - emin)/n
            x0 = st/2
            if (abs((sum2*st-sum1*2*st)/3) <= EPSIL):
                break
            else:
                sum1 = sum2
        qtun = sum2 * st
        return qtun

def total(n, x0, st, e0, sum):
    global emin,gamma,v1,v22,delta,rt
    sum = 0
    e = emin + x0
    for i in range(0,n):
        alfa=gamma*math.sqrt(e)/v22
        beta=gamma*math.sqrt(e-v1)/v22
        x1 = 2*math.pi*(alfa+beta)
        y1 = 2*math.pi*(alfa-beta)
        z1 = 2*math.pi*delta
        x1 = x1/10
        y1 = y1/10
        z1 = z1/10
        c1 = (math.exp(x1)+math.exp(-x1))/2 
        d1 = (math.exp(y1)+math.exp(-y1))/2
        f1 = (math.exp(z1)+math.exp(-z1))/2
        s = (c1 - d1)/(c1+f1)
        s = s*math.exp((e0 - e)/rt)/rt
        sum = sum + s
        e = e + st
    return sum
# def main():
if __name__ == "__main__":
        # Duong dan
    duongdan = os.getcwd()
    # Tao thu muc theo duong dan bang ham mkdir()
    if(os.path.exists(duongdan + '\output') == False):
      os.mkdir(duongdan + '\output')
      
    # run file ke
    file1 = open('input//ke_input.txt', 'r')
    outputfile = open('output//output_ke.csv', 'w+', encoding='UTF8', newline='')
    energyfile = open('output//energy.csv', 'w+', encoding='UTF8', newline='')
    ncrfile = open('output//ncr(e).csv', 'w+', encoding='UTF8', newline='')
    ke.ke(file1, outputfile, energyfile, ncrfile)
    file1.close()
    outputfile.close()
    energyfile.close()
    ncrfile.close()
    
    # run file cvtst
    f = open('input//cvtst_input.txt', 'r')
    f1 = open('output//bimoout.csv', 'w+', encoding='UTF8', newline='')
    f2 = open('output//bimopaste.csv', 'w+', encoding='UTF8', newline='')  
    cvtst.cvtst(f, f1, f2)
    f.close()
    f1.close()
    f2.close()
    
    inputArr = []
    with open("input//ktst_input.txt", 'r') as file5:
        lines = file5.readlines()
        for line in lines:
            inputArr.append(line.strip())
    
    output_ktst = open('output//output_ktst.csv', 'w+', encoding='UTF8', newline='')
    writer = csv.writer(output_ktst)
            
    t,eaa,vimag,qtun = (0 for i in range(4))
    nira = 0
    qx = 0
    q = 0
    
    # khai báo hằng số
    an=6.0225E+23
    h=2.8592E-03
    r=1.9872E-03

    Title = inputArr[0]
    na,nb,nx = list(map(float, inputArr[2].split(',')))

    dgn = (list(map(float, inputArr[4].split(','))))[0]

    vimag = (list(map(float, inputArr[6].split(','))))[0]

    fhr,dihr,rshr = list(map(float, inputArr[8].split(',')))


    if na != 0:
        wa,wb = list(map(float, inputArr[10].split(',')))

        ai1,ai2,ai3 = list(map(float, inputArr[12].split(',')))

        bi1,bi2,bi3 = list(map(float, inputArr[14].split(',')))

        xi1,xi2,xi3 = list(map(float, inputArr[16].split(',')))

        dl = (list(map(float, inputArr[18].split(','))))[0]

        ea,eaa = list(map(float, inputArr[20].split(',')))

        nirx,nira,nirb = list(map(float, inputArr[22].split(',')))

        if nirx != 0 :
            xir = (list(map(float, inputArr[23].split(','))))
            fvxir = (list(map(float, inputArr[23].split(','))))
        if nira != 0 :
            air = (list(map(float, inputArr[23].split(','))))
            dira = (list(map(float, inputArr[23].split(','))))
        if nirb != 0 :
            bir = (list(map(float, inputArr[23].split(','))))
            dirb = (list(map(float, inputArr[23].split(','))))
        fvx = (list(map(float, inputArr[24].split(','))))

        fva = (list(map(float, inputArr[26].split(','))))


        if nb != 0:
            FVB = (list(map(float, inputArr[28].split(','))))
        
        # # in ra màn hình console
        # print("----------"+Title+"---------------")
        # print("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f" % (na,nb,nx,dgn,vimag,fhr,dihr,rshr))
        # print("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f" % (wa,wb,ai1,ai2,ai3,bi1,bi2,bi3))
        # print("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f" % (dl,xi1,xi2,xi3,ea,eaa))
        # print("%10.4f%10.4f%10.4f" % (nirx,nira,nirb))
        vhr = 0.0
        writer.writerow(['THE RATE CONSTANTS ARE LISTED BELOW'])
        writer.writerow([])
        writer.writerow(['Temperature','k(T),cm^3/mole.s','k(T),cm^3/molecule.s','Eckart coefficients'])    
        k = list(map(int, inputArr[28].split(',')))[0]
        for i in range(1,k+1):

            t = list(map(float, inputArr[29+i].split(',')))[0]
            rt = r*t
            qhr = 1.0
            if vimag <= 0:
                qtun = 1.0
            qtun = eckart(t, ea,eaa,vimag,qtun)
            t1 = t
            q,qx = partif(t1,q,nira,dira,air,nirb,dirb,bir,qx)
            rkcal = dl * (2.083E10) * t * q * math.exp(-ea/rt) *qtun * qhr
            rkcal2 = rkcal / 6.022E+23
            writer.writerow([t,str(rkcal),str(rkcal2),qtun])