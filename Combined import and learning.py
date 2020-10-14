# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 14:52:11 2019

@author: s1655357
"""

import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
'''import eli5'''
'''from eli5.sklearn import PermutationImportance'''
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
import scipy.sparse as sparse
import re
########################################################################################################
os.chdir('M:\Machine_Learning') 

#Import 8kb regions generated in R which contain either a random sample of the genome or a chromatin transition region
random_8kb_sequences = pd.read_csv('1000_8kb_random_samples_upper.tab',sep='\t',names=['Region','Sequence'])
sequence_extract=pd.Series(random_8kb_sequences['Sequence'])
sliced=sequence_extract.str.slice(start=3500,stop=4500)
random_sequences=pd.DataFrame(sliced)

#Slice out the 1000bp inuslator chunk at the core of each CT
ctr_sequences = pd.read_csv('insulator_library.tab',sep='\t',names=['Region','Sequence'])
sequence_extract=pd.Series(ctr_sequences['Sequence'])
sliced=sequence_extract.str.slice(start=3500,stop=4500)
insulator_sequences=pd.DataFrame(sliced)





sequence_extract=pd.Series(ctr_sequences['Sequence'])
sliced=sequence_extract.str.slice(start=3500,stop=4500)
insulator_sequences=pd.DataFrame(sliced)




combined_sequences=insulator_sequences.append(random_sequences)
#print(combined_sequences.head())

insulator_dataset= combined_sequences.copy()
## Running optimised model on whole dataset

##############################################################################
## Initialise
##############################################################################
# Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn import preprocessing
from sklearn.utils import shuffle
from sklearn.metrics import confusion_matrix
from sklearn.metrics import fbeta_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import cross_val_score
import random

import keras
np.random.seed(1337)
from keras.preprocessing import sequence
from keras.optimizers import RMSprop
from keras.models import Sequential
from keras.layers.core import Dense
from keras.layers.core import Dropout
from keras.layers.core import Activation
from keras.layers.core import Flatten
from keras.layers.convolutional import Conv1D
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import CSVLogger
import keras.backend as K



# One hot encoder
def one_hot_encode(df, col, seq_len):
    # Dictionary returning one-hot encoding of nucleotides. 
    nuc_d = {'a':[1,0,0,0],'c':[0,1,0,0],'g':[0,0,1,0],'t':[0,0,0,1], 'n':[0,0,0,0]}
    
    # Creat empty matrix.
    vectors=np.empty([len(df),seq_len,4])
    
    # Iterate through sequences and one-hot encode
    for i,seq in enumerate(df[col].str[:seq_len]): 
        seq = seq.lower()
        a = np.array([nuc_d[x] for x in seq])
        vectors[i] = a
    return vectors

one_hot_encoded=one_hot_encode(insulator_dataset,col='Sequence',seq_len=1000)
#print(one_hot_encoded[1][1])

##############################################################################
# Generating the output labels
ins=[1]*len(insulator_sequences)
data=pd.DataFrame(ins)
rando=[0]*len(random_sequences)
scores=data.append(rando)
y = np.ravel(scores)

########################

x=one_hot_encoded
x=shuffle(x)
from keras.utils import to_categorical
y_binary = to_categorical(y)


train_x,test_x,train_y,test_y=train_test_split(x,y, random_state=1)



################################################################################################################

def train_model(x, y, nb_epoch, conv_layers=2, inp_len=1000, border_mode='same', nbr_filters=128, filter_len=8, 
                conv_dropout2=0.4, conv_dropout3=0, dense_layers=1, nodes1=100, nodes2=40, nodes3=70, dense_dropout1=0, 
                dense_dropout2=0, dense_dropout3=0, loss='binary_crossentropy'):
     
    model = Sequential()
    
    # Convolutional layers
    if conv_layers >= 1:
        model.add(Conv1D(activation="relu", input_shape=(inp_len, 4), padding=border_mode, filters=nbr_filters, 
                         kernel_size=filter_len))
        
    if conv_layers >= 2:
        model.add(Conv1D(activation="relu", input_shape=(inp_len, 1), padding=border_mode, filters=nbr_filters, 
                         kernel_size=filter_len))
        model.add(Dropout(conv_dropout2))
        
    if conv_layers >= 3:
        model.add(Conv1D(activation="relu", input_shape=(inp_len, 1), padding=border_mode, filters=nbr_filters, 
                         kernel_size=filter_len))
        model.add(Dropout(conv_dropout3))
    
    # Flatten
    model.add(Flatten())

    # Dense layers
    if dense_layers >= 1:
        model.add(Dense(nodes1, activation='relu'))
        model.add(Dropout(dense_dropout1))
        
    if dense_layers >= 2:
        model.add(Dense(nodes2, activation='relu'))
        model.add(Dropout(dense_dropout2))
        
    if dense_layers >= 3:
        model.add(Dense(nodes3, activation='relu'))
        model.add(Dropout(dense_dropout3))
        
    # Output layer
    model.add(Dense(1, activation='sigmoid'))
    
    # Compile the model
    adam = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    model.compile(loss=loss, optimizer='adam', metrics=['accuracy'])
    
    model.fit(x, y, batch_size=10, epochs=nb_epoch, verbose=1)
    return model

# Run model
model = train_model(train_x,train_y, nb_epoch=10)
predictions = model.predict(test_x)
######################################################################


#model.save('insulator_training_parameters.h5')

##############################################################################################################################
predict_bin= np.where(predictions>=0.5, 1, 0)

# Analysis
cm_all = confusion_matrix(test_y, predict_bin)
print(cm)


test_acc = (cm_all[0,0]+cm_all[1,1])/np.sum(cm_all) # accuracy achieved on test data; number predicted correctly / total
pred_one_acc = cm_all[1,1]/(cm_all[0,1]+cm_all[1,1]) # proportion of ones predicted correctly on test data
pred_zero_acc = cm_all[0,0]/(cm_all[0,0]+cm_all[1,0]) # prcmoportion of zeros predicted correctly on test data
f0_5 = fbeta_score(predict_bin,test_y, beta=0.5)  # punishes false positives more (i.e. promoters that are predicted as being stress activated but aren't)


print("Test accuracy: " + str(round(test_acc,4)))
print("Ones accuracy: " + str(round(pred_one_acc,4)))  
print("Zeros accuracy: " + str(round(pred_zero_acc,4)))  
print("F0.5 score: " + str(round(f0_5,4)))

# Export results
results = {'Test_accuracy': [test_acc],
           'Ones_accuracy': [pred_one_acc], 'Zeroes_accuracy': [pred_zero_acc],
           'F0.5_score': [f0_5]}
results_df = pd.DataFrame(data=results)
results_df.to_csv('Chow_NN_7_predict_all_results.csv', index=False)
ml_test.to_csv('Chow_NN_7_predict_all_results_2.csv', index=False)
pd.DataFrame(cm_all).to_csv('Chow_NN_7_predict_all_results_cm.csv')


