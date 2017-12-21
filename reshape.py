import os
import numpy as np
import codecs
#输入相对路径
os.chdir("F:\\Training\\amsr2workspace\\")
array = np.zeros((),dtype=int)
array = np.loadtxt(open("quweima2.txt"))
#array.dtype = 'float'
array2 = array.astype(int)
#array2 = array.reshape(361,625)
#print(array2)
print(array2)
#f = codecs.open("data2_20160913_China_rectangle_retrieval_shape.txt", "w")
#for i in array2.shape[0]:
#    f.write(str(i)+'\r\n')
#f.close()
file_name = "quweima2reshape.txt"
np.savetxt(file_name, array2, newline='\n')#, fmt=['%s']*array2.shape[1],newline='\n')

