L = 1
while (L <= LP):
    title = f.readline()
    readGeoESpe = f.readline()
    # tách readGeoESpe lấy 3 giá trị
    E[L] = E[L]*627.50595E0
    f1.write(GEO(L), E(L)/627.50595, SPECIES, '\n')
    BE = 0.0
    AE = 0.0
    WE = 0.0
    WEX = 0.0
    if (SPECIES == "diatomic"):
        # READ(10,*) WT,G0,G1,EE
        # READ(10,*) BE,AE,WE,WEX
        pass
    else:
        # READ(10,*) WT,G0,G1,EE
        pass
    # READ(10,*) A,B,C,SN,SF,F,NF,NV
    A = A * FAC
    B = B * FAC
    C = C * FAC
    F = F * FAC
    #IF(NV.NE.0) READ(10,*)(W(j),j=1,NV)
    #IF(TUN.EQ."yes") READ(10,*) RCV(L)
    GIBBS (NT,TEMP,L,SPECIES,WT,G0,G1,EE,BE,AE,WE,WEX,SN,SF,A,B,C,F,NF,W,NV,DGT,E[L],RES,ZPE,EZ)
    L += 1
