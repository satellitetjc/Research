# Research
practice
from __future__ import print_function, division

import os
import sys
import numpy as np
from keras.layers import Convolution1D, Dense, MaxPooling1D, Flatten
from keras.models import Sequential

def smc_retrieval(channel_size=14, filter_length, nb_input_lw=1, nb_output_smc=1, 
                  nb_filter=4):
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
 
def make_lw_pixel(lw_csv,lw_list):
   for file in os.listdir("E:\\TJC\\2016LW_NC_ASCII\\JULY"):
      lw_csv[file] = np.loadtxt(open(file,"rb"),delimiter=",",skiprows=0)
      lw_list[file] = lw_csv[file].reshape(17271) #每个csv文件171*101大小，171*101=17271
   
   channel = 3
