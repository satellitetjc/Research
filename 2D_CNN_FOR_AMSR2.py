from __future__ import print_function, division

#import os
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
#from keras.constraints import maxnorm
from keras.optimizers import SGD
from keras.layers.convolutional import Conv2D
from keras.layers.convolutional import MaxPooling2D
#from keras.utils import np_utils
from keras import backend as K
K.set_image_dim_ordering('th')

# input suqareDim dimensions
sd_rows, sd_cols = 11, 11

#load data
def load_data():
    #path = "E:\\TJC\\2016LW_NC_ASCII\\JULY\\"
    smc = np.loadtxt(open("E:\\TJC\\2016LW_NC_ASCII\\20160703_smc_shenggui_mean_area.csv","rb"),delimiter=",",skiprows=0) 
    tb = np.loadtxt(open("E:\\TJC\\2016LW_NC_ASCII\\JULY\\gw1am2_20160703_01d_eqma_l3sgt36ha2220220_v_clip.csv","rb"),delimiter=",",skiprows=0)   
    n = 5  #设置扩散的范围，上下左右扩散这么大
    squareDim = 2 * n + 1  #正方形的大小
    trainIn = np.zeros((smc.shape[0] - 2 * n, smc.shape[1] - 2 * n, squareDim, squareDim)) #训练的输入
    trainOut = np.zeros((smc.shape[0] - 2 * n, smc.shape[1] - 2 * n)) #训练的输出
    for i in range(n, tb.shape[0] - n):
        for j in range(n, tb.shape[1] - n):
            trainIn[i - n, j - n] = tb[i - n: i + n + 1, j - n: j + n + 1]
            trainOut[i - n, j - n] = smc[i, j]
        
    trainIn = trainIn.reshape(trainIn.shape[0] * trainIn.shape[1], squareDim, squareDim)
    trainOut = trainOut.reshape(trainOut.shape[0] * trainOut.shape[1], 1)
    
    train_lw = trainIn[:2940]
    train_smc = trainOut[:2940]
    #valid_lw = trainIn[2940:3310]
    #valid_smc = trainOut[2940:3310]
    test_lw = trainIn[3310:]
    test_smc = trainOut[3310:]
   
    return train_lw, train_smc, test_lw, test_smc
    
def smc_retrieval(lr=0.005,decay=1e-6,momentum=0.9):    
    # number of convolutional filters to use  
    nb_filters1, nb_filters2 = 5, 10  
    # size of pooling area for max pooling  
    nb_pool = 2  
    # convolution kernel size  
    nb_conv = 3  
    
    nb_output_smc = 1 
    model = Sequential()
    model.add(Conv2D(nb_filters1, nb_conv, nb_conv, border_mode='valid',
                     input_shape=(nb_output_smc, sd_rows, sd_cols)))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
    model.add(Conv2D(nb_filters2, nb_conv, nb_conv))
    model.add(MaxPooling2D(pool_size=(nb_pool,nb_pool)))
    model.add(Dropout(0.25))
    
    model.add(Flatten())
    model.add(Dense(1000))
    model.add(Activation('relu'))
    model.add(Dropout(0.5))
    model.add(Dense(nb_output_smc))
    model.add(Activation('linear'))
    
    sgd = SGD(lr=lr, decay=decay,momentum=momentum, nesterov=True)
    model.compile(loss='mse', optimizer = sgd, metrics=['mae'])
    
    return model
    
def evaluate_train():
    model = smc_retrieval()
    model.summary()
    
    train_lw, train_smc, test_lw, test_smc = load_data()
    model.fit(train_lw, train_smc, nb_epoch=1000, batch_size=200, validation_data=(test_lw, test_smc))
    
    pred = model.predict(test_lw)
    print('\n\nactual', 'predicted', sep='\t')
    for test_lw, actual, predicted in zip(test_lw, test_smc, pred.squeeze()):
        print(test_lw.squeeze(), actual.squeeze(), predicted, sep='\t') 
    
def main():
    np.set_printoptions(threshold=1) #输出起始点
    evaluate_train()
    
if __name__ == '__main__':
    main()  
