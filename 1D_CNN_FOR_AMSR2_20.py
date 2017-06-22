#说明：通道维卷积的新版本，加了batchnormalzation层，加了dropout层以防过拟合，2016年7月28日的数据，加了全连接层1000个节点，打乱了数据，验证样本和测试样本不是同一个
#2017.6.21更：似乎不要batchnormalzation更好，8.28号的数据，中国中部地区（2537个像元），loss层mse，metrics层只是评估性能，不用来训练网络，epoch8000，batchsize300
#R方能达到0.65左右
from __future__ import print_function, division
import os
import numpy as np
from keras.layers import Convolution1D, Dense, MaxPooling1D, Flatten, BatchNormalization, Dropout, Activation
from keras.models import Sequential

def smc_retrieval(channel_size, filter_length, nb_input_lw, nb_output_smc, 
                  nb_filter):
        model = Sequential()
        # The first conv layer learns `nb_filter` filters (aka kernels), each of size ``(filter_length, nb_input_series)``.
        # Its output will have shape (None, window_size - filter_length + 1, nb_filter), i.e., for each position in
        # the input timeseries, the activation of each filter at that position.
        model.add(Convolution1D(nb_filter=nb_filter, filter_length=filter_length, 
                      activation='relu', input_shape=(channel_size, 
                                                      nb_input_lw)))
        model.add(MaxPooling1D())    # Downsample the output of convolution by 2X.
        model.add(Convolution1D(nb_filter=nb_filter, filter_length=filter_length, 
                      activation='relu'))
        model.add(MaxPooling1D())
        model.add(Dropout(0.05))
        model.add(Flatten())
        model.add(Dense(1000))
        model.add(Activation('relu'))
        #model.add(Dropout(0.5))
        #model.add(BatchNormalization())
        model.add(Dense(nb_output_smc, activation='linear'))     # For binary classification, change the activation to 'sigmoid'
        model.compile(loss='mse', optimizer='adam', metrics=['mae'])
    # To perform (binary) classification instead:
    # model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['binary_accuracy'])
        return model

def load_save():
    path = "E:\\TJC\\workspace\\shuju\\"
    j = 0
    all_csv = np.zeros((52,100,15),dtype=float)
    for i in os.listdir(path):
        all_csv[:,:,j] = np.loadtxt(open(i,"rb"),delimiter=",",skiprows=0)
        j = j + 1
    #print(lw_csv)
    #print(lw_csv)
    reall_csv = all_csv.reshape(5200,15)
    #print(relw_csv)
    A = reall_csv.reshape(5200,15,1)
    #smc_csv = np.loadtxt(open("E:\\TJC\\2016LW_NC_ASCII\\20160703_smc_region.csv","rb"),delimiter=",",skiprows=0) 
    #Y = smc_csv.reshape(12972,1)
    #X_test = X[:3]
    #X_luan = np.random.shuffle(X_test)
    #print(X_luan)
    #print(Y)
    #S = np.zeros((5200,15,1),dtype=float)
    S = np.random.permutation(A)
    X = S[:, :14, :]
    Y_L = S[:, 14:, :]
    Y = Y_L.reshape(5200,1)
    #print(Y_L)
    #print(Y)
    return X,Y
load_save()

def evaluate_train():
    filter_length = 4
    nb_filter = 30
    model = smc_retrieval(channel_size=14, filter_length=filter_length,
                          nb_input_lw=1, nb_output_smc=1, nb_filter=nb_filter)
    model.summary()
    
    X, Y = load_save()
    #test_size = int(0.1 * 5200)
    X_train, X_val, X_test, Y_train, Y_val, Y_test = X[:4160], X[4160:4680], X[4680:], Y[:4160], Y[4160:4680], Y[4680:]
    model.fit(X_train, Y_train, nb_epoch=4000, batch_size=500, validation_data=(X_val, Y_val))
    pred = model.predict(X_test)
    print('\n\nactual', 'predicted', sep='\t')
    for X_test, actual, predicted in zip(X_test, Y_test, pred.squeeze()):
        print(X_test.squeeze(), actual.squeeze(), predicted, sep='\t')
    #score = model.accuracy_score(actual.squeeze(), predicted)
    #print('Accuracy: %.2f%%' % (score * 100))
    #print('next', model.predict(q).squeeze(), sep='\t')
    
def main():
    np.set_printoptions(threshold=1) 
    evaluate_train()
    
if __name__ == '__main__':
    main()
