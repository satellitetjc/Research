import os
import numpy as np
import codecs
#此程序首先将aoi区域二维矩阵转化为一维矩阵，之后去除矩阵中的0值，生成一条有效的数据链，作为神经网络的训练数据，可有效将smc与tb的每个栅格一一对应，输出的是txt格式
os.chdir("F:\\Training\\amsr2workspace\\")
array = np.zeros((264,335),dtype=float)
array = np.loadtxt(open("gw1am2_20160828_01d_eqma_l3sgsmchf2210210_clip.csv"),delimiter=",",skiprows=0)
array2 = array.reshape(88440)
#print array2
array3 = array2[array2 > 0]

print(array3)
f = codecs.open("gw1am2_20160828_01d_eqma_l3sgsmchf2210210_clip_chu0.txt","w")
file_name = "gw1am2_20160828_01d_eqma_l3sgsmchf2210210_clip_chu0.txt"
np.savetxt(file_name,array3,newline='\n')
