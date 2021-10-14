from numpy import cos, sin, exp, pi, sqrt, arcsin
import numpy as np
from random import uniform, betavariate
import cmath
import matplotlib.pyplot as plt
import sys
from numass import nuM_range

if len( sys.argv ) > 1:
    number = 0
    number = int(sys.argv[1])

else:
    print("Run as $python getArrays.py <NH:0 or IH:1>, for example $python getArrays.py for NH. ")
    sys.exit()

########################################################################### 
###  best fit values,       best fit - 3sigma,   best fit + 3sigma      ###
########################################################################### 
## Values from Three-neutrino fit based on data available in July 2020: JHEP09(2020)178 (Table with SK athmospheric data) ##

T12 = [arcsin(sqrt(0.304)), arcsin(sqrt(0.269)), arcsin(sqrt(0.343))]
T13_N = [arcsin(sqrt(0.02219)), arcsin(sqrt(0.02032)), arcsin(sqrt(0.02410))]
T13_I = [arcsin(sqrt(0.02238)), arcsin(sqrt(0.02052)), arcsin(sqrt(0.02428))] 
## all in meV**2
dM2_12 = [7.42E1, 6.82E1, 8.04E1]
dM2_N = [2.517E3, 2.435E3, 2.598E3]
dM2_13_N = [i+0.5*j for i,j in zip(dM2_N, dM2_12)]
dM2_I = [2.498E3, 2.414E3, 2.581E3]
dM2_13_I = [i-0.5*j for i,j in zip(dM2_I, dM2_12)]

########### STERILE Neutrino Mixing ##########
## Outdated values (as of October 2021), need to update ##

T14 = [arcsin(sqrt(0.09))/2, arcsin(sqrt(0.09))/2, arcsin(sqrt(0.09))/2]
dM2_14 = [1.78E6,1.78E6,1.78E6]  

ml = np.linspace(0.01,1,100)
ml = np.append(ml, np.linspace(1.1,10,90))
ml = np.append(ml, np.linspace(11,100,90))
ml = np.append(ml, np.linspace(110,1000,90))
N =len(ml)

mbbN = np.zeros((N,4)) 
mbbI = np.zeros((N,4))
mbN = np.zeros((N,3))
mbI = np.zeros((N,3))
sumMN = np.zeros((N,3))
sumMI = np.zeros((N,3))

N_samples = 50000

if(int(sys.argv[1]))==0:    
    print("Normal ordering")
    for i in range (len(ml)):
        mbbN[i], mbN[i], sumMN[i] = \
                nuM_range(ml[i], 'NH', False, [T12, dM2_12], [T13_N, dM2_13_N],
                [T14, dM2_14], N_samples)
        
    with open('normalArray.npy', 'wb') as f:
        np.save(f, mbbN)
        np.save(f, mbN)
        np.save(f, sumMN)

    
else:
    print("Inverted ordering")
    for i in range (len(ml)):
        mbbI[i], mbI[i], sumMI[i] = \
                nuM_range(ml[i], 'IH', False, [T12, dM2_12], [T13_I, dM2_13_I],
                [T14, dM2_14], N_samples)
        
    with open('invertedArray.npy', 'wb') as f:
        np.save(f, mbbI)
        np.save(f, mbI)
        np.save(f, sumMI)