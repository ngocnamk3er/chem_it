import os
def inputfile():
  # Duong dan
  duongdan = os.getcwd()
  # Tao thu muc theo duong dan bang ham mkdir()
  if(os.path.exists(duongdan + '\input') == False):
    os.mkdir(duongdan + '\input')
  
  ke_input = open('input//ke_input.txt','w+')
  ktst_input = open('input//ktst_input.txt','w+')
  cvtst_input = open('input//cvtst_input.txt','w+')
  line_total = 0
  line_cvtst = 0
  line_ktst = 0
  with open("input.txt", 'r') as f:
      for readlline in f:
        line_total += 1
        if(readlline == 'cvtst_input\n'):
          line_cvtst = line_total
        if(readlline == 'ktst_input\n'):
          line_ktst = line_total
  line = 0
  with open("input.txt", 'r') as f:
      for readlline in f:
        line += 1
        if(line < line_cvtst):
          ke_input.write(readlline)
        elif(line > line_cvtst and line < line_ktst):
          cvtst_input.write(readlline)
        elif(line > line_ktst and line <= line_total):
          ktst_input.write(readlline)
  ke_input.close()
  cvtst_input.close()
  ktst_input.close()

inputfile()    

