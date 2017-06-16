# 记得打开spyder时记得打开tensorflow环境
# 这个程序和lw、smc数据都放在E:\\TJC\\2016LW_NC_ASCII\\下，亮温数据同时在E:\\TJC\\2016LW_NC_ASCII\\JULY文件夹下
from __future__ import print_function, division

import os
import numpy as np
from keras.layers import Convolution1D, Dense, MaxPooling1D, Flatten
from keras.models import Sequential

def smc_retrieval(channel_size, filter_length, nb_input_lw, nb_output_smc, 
                  nb_filter):
    model = Sequential((
        # The first conv layer learns `nb_filter` filters (aka kernels), each of size ``(filter_length, nb_input_series)``.
        # Its output will have shape (None, window_size - filter_length + 1, nb_filter), i.e., for each position in
        # the input timeseries, the activation of each filter at that position.
        Convolution1D(nb_filter=nb_filter, filter_length=filter_length, 
                      activation='relu', input_shape=(channel_size, 
                                                      nb_input_lw)),
        MaxPooling1D(),     # Downsample the output of convolution by 2X.
        Convolution1D(nb_filter=nb_filter, filter_length=filter_length, 
                      activation='relu'),
        MaxPooling1D(),
        Flatten(),
        Dense(nb_output_smc, activation='linear'),     # For binary classification, change the activation to 'sigmoid'
    ))
    model.compile(loss='mse', optimizer='adam', metrics=['mae'])
    # To perform (binary) classification instead:
    # model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['binary_accuracy'])
    return model

def load_save():
    path = "E:\\TJC\\2016LW_NC_ASCII\\JULY\\"
    j = 0
    lw_csv = np.empty((101,171,14),dtype=float)
    for i in os.listdir(path):
        lw_csv[:,:,j] = np.loadtxt(open(i,"rb"),delimiter=",",skiprows=0)
        j = j + 1
        #print(lw_csv)
    relw_csv = lw_csv.reshape(17271,14)
    #print(relw_csv)
    X = relw_csv.reshape(17271,14,1)
    smc_csv = np.loadtxt(open("E:\\TJC\\2016LW_NC_ASCII\\smc_201607.csv","rb"),delimiter=",",skiprows=0) 
    Y = smc_csv.reshape(17271,1)
   # print(X)
   # print(Y)
    return X, Y       
#load_save()

def evaluate_train():
    filter_length = 4
    nb_filter = 4
    model = smc_retrieval(channel_size=14, filter_length=filter_length,
                          nb_input_lw=1, nb_output_smc=1, nb_filter=nb_filter)
    model.summary()
    
    X, Y = load_save()
    test_size = int(0.2 * 17271)
    X_train, X_test, Y_train, Y_test = X[:-test_size], X[-test_size:], Y[:-test_size], Y[-test_size:]
    model.fit(X_train, Y_train, nb_epoch=25, batch_size=2, validation_data=(X_test, Y_test))
    pred = model.predict(X_test)
    print('\n\nactual', 'predicted', sep='\t')
    for actual, predicted in zip(Y_test, pred.squeeze()):
        print(actual.squeeze(), predicted, sep='\t')
    #print('next', model.predict(q).squeeze(), sep='\t')
    
def main():
    np.set_printoptions(threshold=25)
    evaluate_train()
    
if __name__ == '__main__':
    main()

