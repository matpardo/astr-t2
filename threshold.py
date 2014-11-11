#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import pylab as pl
import numpy as np
from scipy import random

# Realiza thresholding en data con un valor k.
# Entrega las posiciones en data en las cuales se encuentran los datos que sobrepasan el threshold
# Y el valor del threshold
def thresholding(data,k):

    mean_n = np.mean(data)
    std_n = np.std(data)
    threshold    = mean_n + k * std_n
    source_pos   = []
    for i in range(len(data)):
        if data[i] > threshold:
           source_pos.append(i)

    # Graficar
    val_sources = []
    FDR = 0
    for pos in source_pos:
        val_sources.append(data[pos])
        if not 65 <= pos <= 136 and pos != 180 and pos != 220:
           FDR += 1

    print_args = (k,FDR,len(source_pos),float(FDR)/len(source_pos))
    print "Threshold : k = %.1f & FDR = %d/%d = %f" % print_args

    pl.clf()
    pl.plot(data,'b')
    pl.plot(source_pos,val_sources,'ok')
    pl.plot(range(0,len(data)),[threshold]*N,'-g')
    ymin = np.min(data)
    ymax = np.max(data)
    pl.ylim(ymin-1.0,ymax+1.0)
    pl.show()

    return source_pos,threshold

# Para test

def noise_values(signal):

    data = signal[:]
    while True:
       noise = []
       sigma = np.std(data)
       for val in data:
           if val < 3 * sigma:
              noise.append(val)
       if len(noise) == len(data): break
       data = noise[:]
    
    return np.mean(data),np.std(data)

def gaussian(x,stdev=1.0,mean=0.0):

    return np.exp(-1.0*(x-mean)**2.0/(2*stdev**2.0))/(stdev*np.sqrt(2*np.pi))

# Señal: se tiene una gaussiana y dos deltas como fuentes

N = 250
signal = np.zeros(N)
for i in range(65,136):
    signal[i] = 250.0 * gaussian(i,10.0,100.0) # Fuente (intensa) entre pixeles 65 y 135
signal[180] = 10.0 # Fuente (intensa) en pixel 180
signal[220] = 1.5  # Fuente (débil) en pixel 220
n_signal = signal + random.standard_normal(N) # Ruido con media 0 y desviación estándar 1

# Ejemplos: se usa thresholding con k = {2.0,2.5,3.0}

k = [2.0,2.5,3.0]
for i in range(len(k)):
    pos_sources,threshold = thresholding(n_signal,k[i])