import os
import numpy as np
import codecs
#此代码用于拼接数据，因为训练数据是多天，需将多天数据拼接一起，一并输入
os.chdir("F:\\Training\\amsr2workspace\\")
a = np.zeros((88440),dtype=float)
a = np.loadtxt(open("gw1am2_20160828_01d_eqma_l3sgsmchf2210210_clip_chu0.txt"))
b = np.zeros((88440),dtype=float)
b = np.loadtxt(open("smc_training_shuju.txt"))
c = np.concatenate((a,b),axis=0)
#拼接用concatenate函数，axis : 沿着某个轴拼接，默认为列方向
#合并用np.row_stack()函数，这是列合并，行合并是np.column_stack
f = codecs.open("concatenate.txt","w")
file_name = "concatenate.txt"
np.savetxt(file_name,c,newline='\n')
