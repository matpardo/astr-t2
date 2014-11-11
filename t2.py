#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import math


arch   = pyfits.open('blank.fits')
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

#Inserta una estrella con los parametros m,RA,DEC en la imagen hdu
def addStar(hdu,m,RA,DEC):
    ver,hor = ra_dec_2_pix(RA,DEC)
    max_ver,max_hor = hdu.shape
    if 0<=ver<max_ver and 0<hor<max_hor:
        cuentas = mToCounts(m,20,flux20)
        hdu[ver][hor]=cuentas
    return

#Lee un catalogo de estrellas y los agrega a la imagen fits
def addStellarCatalog (hdu, catalog):
    fp = open(catalog, 'r')
    for line in fp:
        obj,RA,DEC,r_mag,SED,indx,star = line.split()
        addStar(hdu,float(r_mag),float(RA),float(DEC))
    return

#Funcion para calcular radio eliptico dados x e y
def getER(x,y,xc,yc,el,theta):
    return math.sqrt(((x-xc)*math.cos(theta)-(y-yc)*math.sin(theta))**2+( (((x-xc)*math.sin(theta)-(y-yc)*math.cos(theta))**2)/((1-el)**2)))

#Agrega una galaxia a la imagen hdu
def addGalaxy (hdu, m, RA, DEC, n, Re, el, theta):
    #obtener coordenadas
    ver,hor = ra_dec_2_pix(RA,DEC)
    max_ver,max_hor = hdu.shape
    if 0<=ver<max_ver and 0<=hor<max_hor:
        #obtener cuentas
        cuentas = mToCounts(m,20,flux20)
        #calcular I0
        bn=2*n-0.324
        I0 = (cuentas*(math.pow(bn,2*n)))/(math.pow(Re,2)*2*math.pi*n*(1-el)*math.gamma(2*n))
        hdu[ver][hor]=I0
        #sacar delta max para delimitar cuadro a pintar
        #dmax=int(np.maximum(Re,Re*(1-el))*(2))
        #Para optimizar, se ocupan coordenadas polares, dejando de dibujar cuando las magnitudes sean cero.
        ang=0
        dang=2*math.pi/50
        while ang<2*math.pi:
            #print("ang={}".format(ang))
            d=1
            Ier=I0
            while Ier>0.1:
                x=hor+d*math.cos(ang)
                y=ver+d*math.sin(ang)
                if(d>1000):
                    print("ang={} d={} ({},{})".format(ang,d,x,y))
                
                d+=1
                if not (0<=y<max_ver and 0<=x<max_hor):
                    break
                Er=getER(x,y,hor,ver,el,theta)
                Ier=I0*math.exp(-bn*(math.pow((Er/Re),(1/n))))
                if hdu[y][x]<Ier:
                    hdu[y][x]=Ier
            ang+=dang
        '''
        for x in range(int(hor-dmax),int(hor+dmax)):
            for y in range(int(ver-dmax),int(ver+dmax)):
                #si esta fuera de la imagen, gg
                if not (0<=y<max_ver and 0<=x<max_hor):
                    continue
                Er=getER(x,y,hor,ver,el,theta)
                Ier=I0*math.exp(-bn*(math.pow((Er/Re),(1/n))))
                if hdu[y][x]<Ier:
                    hdu[y][x]=Ier
        '''
        
    return

#Lee un catalogo de galaxias y los agrega a la imagen fits
def addGalaxyCatalog (hdu, catalog):
    #contador solo para no cargar el catalogo completo en las pruebas
    fp = open(catalog, 'r')
    counter=0
    for line in fp:
        counter+=1
        if counter%1000==0:
            print("counter={}".format(counter))
        #print(line)
        obj,RA,DEC,r_mag,SED,redshift,galaxy,n,Re,el,theta=line.split()
        #if float(n)<5 and float(r_mag)>5:
        #    continue
        addGalaxy(hdu,float(r_mag),float(RA),float(DEC),float(n),float(Re),float(el),float(theta))
    return counter

#addStellarCatalog(img,"stellar.dat")
#plot_image(img, log_scale=True)  

cc=addGalaxyCatalog(img,"galaxy.dat")
print(cc)
plot_image(img,log_scale=True)

#addStellarCatalog(img,"stellar.dat")
#plot_image(img, log_scale=True)    

#Sube un nivel de intensidad constante a una imagen
def addBackground (hdu, background):
    return -1

#Realiza convolución de una imagen con un point spread function (PSF) Gaussiano
# con dev std sigma_PSF
def convolvePSF (hdu, sigma_PSF):
    return -1

#Agrega ruido Poissoniano y Gaussiano con desv std sigma_noise
def addNoise (hdu, sigma_noise):
    return -1

#Filtro para reducir ruido de la imagen.
def filterImage (hdu, omega_1, omega_2):
    #La idea es aplicar un filtro pasa-bajos usando la máscara 3x3 Genérica
    #Para esto debemos aplicar la fórmula a cada pixel y considerar casos bordes
    #obtener coordenadas
    ver,hor = ra_dec_2_pix(RA,DEC)
    max_ver,max_hor = hdu.shape

    #valores que se usarán con frecuencia
    dos_cos_omega1 = 2 * cos(omega_1)
    dos_cos_omega2 = 2 * cos(omega_2)

    suma_2cos = dos_cos_omega1 + dos_cos_omega2
    mult_2cos = dos_cos_omega1 * dos_cos_omega2

    #Para cada Fila
    for y in xrange (0, max_ver):
        #Para cada columna
        for x in xrange (0, max_hor):
            #usaremos variables para saber si hay que aplicar operación en pixel o no
            upper_left = 1
            upper_center = 1
            upper_right = 1
            left = 1
            right = 1
            down_left = 1
            down_center = 1
            down_left = 1
            #si está en la parte superior
            if y == 0:
                upper_left = 0
                upper_center = 0
                upper_right = 0
            #si está en la parte de abajo
            elif y == (max_ver - 1):
                down_left = 0
                down_center = 0
                down_right = 0
            #si está en la parte izq
            if x == 0:
                upper_left = 0
                left = 0
                down_left = 0
            #si está en la parte derecha
            elif x == (max_hor - 1):
                upper_right = 0
                right = 0
                down_right = 0
            #calcular
            if upper_left:
                hdu[x - 1][y - 1] = mult_2cos * hdu[x - 1][y - 1]
            if upper_center:
                hdu[x][y - 1] = suma_2cos * hdu[x][y - 1]
            if upper_right:
                hdu[x + 1][y - 1] = mult_2cos * hdu[x + 1][y - 1]
            if left:
                hdu[x - 1][y] = suma_2cos * hdu[x - 1][y]
            if right:
                hdu[x + 1][y] = suma_2cos * hdu[x + 1][y]
            if down_left:
                hdu[x - 1][y + 1] = mult_2cos * hdu[x - 1][y + 1]
            if down_center:
                hdu[x][y + 1] = suma_2cos * hdu[x][y + 1]
            if down_right:
                hdu[x + 1][y + 1] = mult_2cos * hdu[x + 1][y + 1]
    
    return hdu
