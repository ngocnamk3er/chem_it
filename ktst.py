import math
# định nghĩa các hàm
def eckart(t,e0,e1,vimag,qtun):
    return qtun

def partif(t,q,nira,dira,air,nirb,dirb,bir,qx):
    return qx

def total(n, x0, st, e0, sum):
    return sum

# đọc dữ liệu từ file input.txt và lưu vào mảng inputArr
inputArr = []
file1 = open('input.txt', 'r')
Lines = file1.readlines()
for line in Lines:
	inputArr.append(line.strip())
#khai báo các biến

fva, fvb, fvx, air, dira = ([] for i in range(5)) #/FVS/
bir, dirb, xir, dirx, fvxir, virx =  ([] for i in range(6)) #  /HRA/
na, nb, nx, nirx = ( 0 for i in range(4))
h,r,an,rt,wa,wb,ea = ( 0 for i in range(7))
ai1,ai2,ai3,bi1,bi2,bi3,xi1,xi2,xi3 = ( 0 for i in range(9))
T,eaa,vimag,qtun = (0 for i in range(4))
#khai báo các biến chưa được khai báo và không được gán trong fortran, giá trị mặc định là 0
#khai báo biến được khai báo ở tham số trong subroutine fortran
nira = 0
qx = 0
q = 0
#format

#-----

# khai báo hằng số
AN=6.0225E+23
H=2.8592E-03
R=1.9872E-03

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

    nirx,a,nirb = list(map(float, inputArr[22].split(',')))

    if nirx != 0 :
        xir = (list(map(float, inputArr[23].split(','))))
        fvxir = (list(map(float, inputArr[23].split(','))))
        # read (5,*) (XIR(I), FVXIR(I),I=1,NIRx) 
        # READ (5,*) (DIRx(I),I=1,NIRx) 
    if nira != 0 :
        air = (list(map(float, inputArr[23].split(','))))
        dira = (list(map(float, inputArr[23].split(','))))
        # read (5,*) (aIR(I),I=1,NIRa)
        # READ (5,*) (DIRa(I),I=1,NIRa) 
    if nirb != 0 :
        bir = (list(map(float, inputArr[23].split(','))))
        dirb = (list(map(float, inputArr[23].split(','))))
        # read (5,*) (bIR(I),I=1,NIRb)
        # READ (5,*) (DIRb(I),I=1,NIRb) 
    fvx = (list(map(float, inputArr[24].split(','))))

    fva = (list(map(float, inputArr[26].split(','))))


    if nb != 0:
        FVB = (list(map(float, inputArr[28].split(','))))
    
    # in ra màn hình console
    print("----------"+Title+"---------------")
    print("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f" % (na,nb,nx,dgn,vimag,fhr,dihr,rshr))
    print("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f" % (wa,wb,ai1,ai2,ai3,bi1,bi2,bi3))
    print("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f" % (dl,xi1,xi2,xi3,ea,eaa))
    print("%10.4f%10.4f%10.4f" % (nirx,nira,nirb))
    #-----------------------
    if nirx !=0 :
        print(xir)
        print(dirx)
    if nira != 0:
        print(air)
        print(dira)
    if nirb != 0:
        print(bir)
        print(dirb)
    print(fvx)
    print(fva)
    print(fvb)
    vhr = 0.0
    print("VHR = %.4f kcal/mol -- EA = %.4f kcal/mol" % (vhr,ea))

    # ghi tiêu đề ra file
    L = ["                         THE RATE CONSTANTS ARE LISTED BELOW\n", 
        "-----------------------------------------------------------------------------------------------------\n", 
        "Temperature          k(T),cm^3/mole.s          k(T),cm^3/molecule.s          Eckart coefficients\n",
        "-----------------------------------------------------------------------------------------------------"]

    # writing to file
    file1 = open('output.txt', 'w')
    file1.writelines(L)
    # file1.close()

    # đọc từ file input tiếp 
 
    k = list(map(int, inputArr[28].split(',')))[0]

    # vòng lặp sau đọc đoạn The range of temperatures và xử lí kết quả:

    for i in range(1,k+1):

        T = list(map(float, inputArr[29+i].split(',')))[0]
        print(T)
        if T <=0 :
            print("Go to # đọc dữ liệu từ file input.txt")
        rt = R*T
        qhr = 1.0
        if vimag <= 0:
            qtun = 1.0
        qtun = eckart(T, ea,eaa,vimag,qtun)
        t1 = T
        qx = partif(t1,q,nira,dira,air,nirb,dirb,bir,qx)
        rkcal = dl * (2.083E10) * T * q * math.exp(-ea/rt) * qhr
        rkcal2 = rkcal / 6.022E+23
        # print(str(T)+"---"+str(rkcal)+"---"+str(rkcal2)+"---"+str(qtun))
        
        space = "                    "
        # file1.writelines("\n"+str(T)+space+str(rkcal)+space+str(rkcal2)+space+str(qtun)+"\n")
        file1.writelines("\n%10.4f%20.4f%30.4f%30.4f\n" % (T, rkcal, rkcal2, qtun))

        






    

        
        









    

    
                






