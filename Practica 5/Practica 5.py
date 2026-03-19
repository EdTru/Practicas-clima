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


a_l, b_l, cosa, cosa, cosa3 = linregress(meses, ristra_datos)
print(linregress(meses, ristra_datos))
prediccion_lin = np.array(meses)*a_l + b_l

plt.scatter(meses, ristra_datos, s=10, c="red",label="Datos ERAS5")
plt.plot(meses,prediccion_theisen,c="blue", label=rf"Ajuste Thei-Sen $m=${a.round(4)}")

plt.plot(meses,prediccion_lin,c="green", label=rf"Ajuste Lineal $m=${a_l.round(4)}")

plt.legend()

plt.xlabel("Años")
plt.ylabel("Promedio temperatura (ºC)")
plt.title("Promedio temperatura en el mes Julio en Madrid en los años 1940-2025")

resultao = mk.original_test(ristra_datos)
print(resultao)
print(a_l, b_l)
print("pendiente Teil-Sen"+str(a),"intercept"+str(b))
print("pendiente Mann-Kendall"+str(resultao[-2]),"intercept"+str(resultao[-1]))

plt.show()