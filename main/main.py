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
    outputfile = open('output//output_ke.csv', 'w+', encoding='UTF8', newline='')
    energyfile = open('output//energy.csv', 'w+', encoding='UTF8', newline='')
    ncrfile = open('output//ncr(e).csv', 'w+', encoding='UTF8', newline='')
    ke.ke(file1, outputfile, energyfile, ncrfile)
    file1.close()
    outputfile.close()
    energyfile.close()
    ncrfile.close()
    
    # run file cvtst
    f = open('..//CVTST-m_py//bimo.txt', 'r')
    f1 = open('output//bimoout.csv', 'w+', encoding='UTF8', newline='')
    f2 = open('output//bimopaste.csv', 'w+', encoding='UTF8', newline='')  
    cvtst.cvtst(f, f1, f2)
    f.close()
    f1.close()
    f2.close()
