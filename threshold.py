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