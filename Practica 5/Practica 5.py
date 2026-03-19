from netCDF4 import Dataset #Módulo para leer el archivo de datos
import numpy as np
import matplotlib.pyplot as plt #Módulo para representar datos en gráficas
import scipy.stats as sps
import pymannkendall as mk
from scipy.stats import linregress

datillos = Dataset("Data_P5.nc", "r", format="NETCDF4") #Leemos los datos del archivo de de ERA5

temp, lon, lat = datillos.variables["t2m"], datillos.variables["longitude"][:], datillos.variables["latitude"][:]
data = temp[:]

ristra_datos = []

for i in data:
	ristra_datos.append(float(i[0][0]-273.15))

meses = list(range(1940,2026))

a, b, cosa1, cosa2 = sps.mstats.theilslopes(ristra_datos, meses)
prediccion_theisen = np.array(meses)*a +b
a_l, b_l, cosa, cosa, cosa = linregress(meses, ristra_datos)
prediccion_lin = np.array(meses)*a_l + b_l


plt.scatter(meses, ristra_datos, s=10, c="red")
plt.plot(meses,prediccion_theisen,c="blue")
plt.plot(meses,prediccion_lin,c="green")


resultao = mk.original_test(ristra_datos)
print(resultao)

print("pendiente Teil-Sen"+str(a),"intercept"+str(b))
print("pendiente Mann-Kendall"+str(resultao[-2]),"intercept"+str(resultao[-1]))



plt.show()