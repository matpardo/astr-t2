#!/usr/bin/env python
# -*- coding: utf-8 -*- 

def getData(archivo):
	f=open(archivo)
	cant = f.readline().split()
	dataTable = []
	for line in f:
		datos = line.split()
		data = []
		for dato in datos:
			data.append(float(dato))
		dataTable.append(data)
	return float(cant[0]),float(cant[1]),dataTable

#estrellas,galaxias,tabla = getData("datos.dat")
#print(str(estrellas)+"\n"+str(galaxias))
#print(tabla)
