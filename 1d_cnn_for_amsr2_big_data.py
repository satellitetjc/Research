from __future__ import print_function, division
import os
import numpy as np
from keras.layers import Convolution1D, Dense, MaxPooling1D, Flatten, BatchNormalization, Dropout, Activation
from keras.models import Sequential
from keras.utils import plot_model, normalize

def smc_retrieval(channel_size, filter_length, nb_input_lw, nb_output_smc, 
                  nb_filter):
    model = Sequential((
        # The first conv layer learns `nb_filter` filters (aka kernels), each of size ``(filter_length, nb_input_series)``.
        # Its output will have shape (None, window_size - filter_length + 1, nb_filter), i.e., for each position in
        # the input timeseries, the activation of each filter at that position.
        
        Convolution1D(nb_filter=nb_filter, filter_length=filter_length, 
                      activation='relu', input_shape=(channel_size, 
                                                      nb_input_lw)),
        MaxPooling1D(), 
        Convolution1D(60, filter_length=filter_length, 
                      activation='relu'),
        MaxPooling1D(),
        Dropout(0.05),
        Flatten(),
        Dense(1000),
        #BatchNormalization(),
        Activation('relu'),      
        #BatchNormalization(), 
        #Dropout(0.005),        
        Dense(nb_output_smc, activation='relu'),     # For binary classification, change the activation to 'sigmoid'
    ))
    model.compile(loss='mse', optimizer='adam', metrics=['mae'])
    # To perform (binary) classification instead:
    # model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['binary_accuracy'])
    #plot_model(model, to_file="E:\\TJC\\workspace\\shuju\\model.png")
    return model

def load_save():
    os.chdir("E:\\TJC\\workspace\\shuju\\")
    j = 0
    all_csv = np.zeros((177,43,15),dtype=float)
    for i in os.listdir():
        all_csv[:,:,j] = np.loadtxt(open(i,"rb"),delimiter=",",skiprows=0)
        j = j + 1
    #print(lw_csv)
    #print(lw_csv)
    reall_csv = all_csv.reshape(7611,15)
    #print(relw_csv)
    A = reall_csv.reshape(7611,15,1)
    #smc_csv = np.loadtxt(open("E:\\TJC\\2016LW_NC_ASCII\\20160703_smc_region.csv","rb"),delimiter=",",skiprows=0) 
    #Y = smc_csv.reshape(12972,1)
    #X_test = X[:3]
    #X_luan = np.random.shuffle(X_test)
    #print(X_luan)
    #print(Y)
    #S = np.zeros((5200,15,1),dtype=float)
    S = np.random.permutation(A)
    X = S[:, :14, :]
    X_NOL = normalize(X, axis=1, order=2)
    Y_L = S[:, 14:, :]
    Y = Y_L.reshape(7611,1)
    #print(Y_L)
    #print(Y)
    return X_NOL,Y
load_save()

def test_save():
    os.chdir("E:\\TJC\\workspace\\test\\")
    k = 0
    all_csv2 = np.zeros((59,43,14),dtype=float)
    for l in os.listdir():
        all_csv2[:,:,k] = np.loadtxt(open(l,"rb"),delimiter=",",skiprows=0)
        k = k + 1
    reall_csv2 =all_csv2.reshape(2537,14)
    B = reall_csv2.reshape(2537,14,1)
    X_test = normalize(B,axis=1,order=2)
    print(X_test)
    return X_test
test_save()

def evaluate_train():
    filter_length = 4
    nb_filter = 30
    model = smc_retrieval(channel_size=14, filter_length=filter_length,
                          nb_input_lw=1, nb_output_smc=1, nb_filter=nb_filter)
    model.summary()
    
    X_NOL, Y = load_save()
    X_test = test_save()
    #test_size = int(0.1 * 5200)
    X_train, X_val, Y_train, Y_val = X_NOL[:6900], X_NOL[6900:], Y[:6900], Y[6900:]
    model.fit(X_train, Y_train, nb_epoch=20000, batch_size=3000, validation_data=(X_val, Y_val))
    pred = model.predict(X_test)
    f = open("log.txt","w")
    print('retrieval', sep='\t',file = f)
    for retrieval in pred.squeeze():
        print(retrieval, sep='\t',file = f)
        #squeeze()有压缩的意思，所以不会全部显示？
        #print(X_test.squeeze(), actual.squeeze(), predicted, sep='\t')
    #score = model.accuracy_score(actual.squeeze(), predicted)
    #print('Accuracy: %.2f%%' % (score * 100))
    #print('9.13', model.predict(X_test).squeeze(), sep='\t')
    
def main():
    #np.set_printoptions(threshold=1) 
    evaluate_train()
    #plot_model(smc_retrieval(), to_file="E:\\TJC\\workspace\\shuju\\model.png")
    
if __name__ == '__main__':
    main()
