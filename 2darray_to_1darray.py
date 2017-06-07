import numpy as np

a = np.array([[1,1,1],[2,2,2]])
b = np.array([[6,5,4],[14,15,16]])
c = np.array([[111,111,111],[7,8,9]])

d = a.reshape(6)
e = b.reshape(6)
f = c.reshape(6)

#pixel = 6
channel =3
value = 1

trainin = np.zeros((d.shape[0],channel,value),dtype=float) #6个点，abc3个通道，每个通道1个值
for i in range(0,d.shape[0]):
    trainin[i] = [[d[i]],[e[i]],[f[i]]]
print(trainin)
