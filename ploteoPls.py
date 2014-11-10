#!/usr/bin/env python

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal, misc

arch   = pyfits.open('stellar.fits')
hdr    = arch[0].header
img    = arch[0].data
flux20 = hdr['FLUX20']

ref_pix_hor   = float(hdr['CRPIX1']) # Column Pixel Coordinate of Ref. Pixel  
ref_pix_ver   = float(hdr['CRPIX2']) # Row Pixel Coordinate of Ref. Pixel
ref_pix_ra    = float(hdr['CRVAL1']) # RA at Reference Pixel
ref_pix_dec   = float(hdr['CRVAL2']) # DEC at Reference Pixel

delta_ra_col  = float(hdr['CD1_1'])  # RA  degrees per column pixel
delta_ra_row  = float(hdr['CD1_2'])  # RA  degrees per row pixel
delta_dec_col = float(hdr['CD2_1'])  # DEC degrees per column pixel
delta_dec_row = float(hdr['CD2_2'])  # DEC degrees per row pixel

# De pixeles (ver,hor) a (ra,dec)
def pix_2_ra_dec(ver,hor):

    diff_ver = ver - ref_pix_ver
    diff_hor = hor - ref_pix_hor

    diff_ra  = (diff_hor * delta_ra_col) + diff_ver * delta_ra_row
    diff_dec =  diff_hor * delta_dec_col  + (diff_ver * delta_dec_row)

    calc_ra  = ref_pix_ra  + diff_ra
    calc_dec = ref_pix_dec + diff_dec

    return (calc_ra,calc_dec)

# De (ra,dec) a pixeles (ver,hor)
def ra_dec_2_pix(ra,dec):

    diff_ra  = ra  - ref_pix_ra
    diff_dec = dec - ref_pix_dec

    diff_ver = (diff_dec/delta_dec_col - diff_ra/delta_ra_col)/(delta_dec_row/delta_dec_col - delta_ra_row/delta_ra_col)
    diff_hor = (diff_dec/delta_dec_row - diff_ra/delta_ra_row)/(delta_dec_col/delta_dec_row - delta_ra_col/delta_ra_row)

    calc_ver = ref_pix_ver + diff_ver
    calc_hor = ref_pix_hor + diff_hor

    return (calc_ver,calc_hor)

def plot_image(image,interpolation="nearest",log_scale=False,title=None):
    ''' buenos colores
	'hot' 'gray' 'Greys' 'bone' 'copper' 'gist_heat' 'pink' 'summer' 'afmhot'
    '''
    if log_scale:
       min_value = np.min(image)
       if min_value <= 0: image2 = image + np.abs(min_value) + 10.0
       else: image2 = image - min_value + 10.0
       plt.imshow(np.log(image2),interpolation=interpolation)
       plt.set_cmap('pink')
       plt.colorbar()
    else:
       plt.imshow(image,interpolation=interpolation)
       plt.set_cmap('pink') 
       plt.colorbar()

    if title != None:
       plt.title(str(title))

    plt.show()
    return

# Entrega flujo usando la formula m-m0=-2.5*log(F/F0)
def mToCounts(m,m0,F0):
    F = F0*math.exp((m0-m)/2.5)
    return F

#Funcion para calcular radio eliptico dados x e y
def getER(x,y,xc,yc,el,theta):
    return math.sqrt(((x-xc)*math.cos(theta)-(y-yc)*math.sin(theta))**2+( (((x-xc)*math.sin(theta)-(y-yc)*math.cos(theta))**2)/((1-el)**2)))

plot_image(img, log_scale=True)
print(hdr)
