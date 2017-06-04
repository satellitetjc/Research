import numpy as np

channels = 7  #通道数目
rows = 100   #行数
cols = 100   #列数

groundTruth = np.random.random((rows, cols))  # 随机真实值

obs = np.random.random((rows, cols, channels))  #随机观测值

n = 5  #设置扩散的范围，上下左右扩散这么大
squareDim = 2 * n + 1  #正方形的大小

trainIn = np.zeros((groundTruth.shape[0] - 2 * n, groundTruth.shape[0] - 2 * n, squareDim, squareDim, channels)) #训练的输入
trainOut = np.zeros((groundTruth.shape[0] - 2 * n, groundTruth.shape[0] - 2 * n)) #训练的输出
#groundTruth = np.loadtxt(open("E:\\TJC\\2016SMC_China_ascii\\gw1am2_20160100_01m_eqma_l3sgsmchf2210210_chip.csv","rb"),delimiter=",",skiprows=0)
#obs = np.loadtxt(open("E:\\TJC\\2016LW_China_ascii\\06h\\gw1am2_20160100_01m_eqma_l3sgt06ha2220220_h_chip.csv","rb"),delimiter=",",skiprows=0)  #随机观测值


for i in range(n, groundTruth.shape[0] - n):
    for j in range(n, groundTruth.shape[1] - n):
        trainIn[i - n, j - n, :] = obs[i - n: i + n + 1, j - n: j + n + 1, :]
        trainOut[i - n, j - n] = groundTruth[i, j]
        
trainIn = trainIn.reshape(trainIn.shape[0] * trainIn.shape[1], squareDim, squareDim, channels)
trainOut = trainOut.reshape(trainOut.shape[0] * trainOut.shape[1], -1)
