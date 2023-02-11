import ke
import cvtst
import os
if __name__ == "__main__":
    # Duong dan
    duongdan = os.getcwd()
    # Tao thu muc theo duong dan bang ham mkdir()
    if(os.path.exists(duongdan + '\output') == False):
      os.mkdir(duongdan + '\output')
      
    # run file ke
    file1 = open('..//ke_py//input.txt', 'r')
    outputfile = open('output//output.txt', 'w+')
    energyfile = open('output//energy.txt', 'w+')
    ncrfile = open('output//ncr(e).txt', 'w+')
    ke.ke(file1, outputfile, energyfile, ncrfile)
    file1.close()
    outputfile.close()
    energyfile.close()
    ncrfile.close()
    
    # run file cvtst
    f = open('..//CVTST-m_py//bimo.txt', 'r')
    f1 = open('output//bimoout.txt', 'w+')
    f2 = open('output//bimopaste.txt', 'w+')
    cvtst.cvtst(f, f1, f2)
    f.close()
    f1.close()
    f2.close()
