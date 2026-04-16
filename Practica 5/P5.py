from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
import pymannkendall as mk
from scipy.stats import linregress
import os 

os.chdir(os.path.dirname(os.path.realpath(__file__)))

datillos = Dataset("Data_P5.nc", "r", format="NETCDF4")
temp, lon, lat = datillos.variables["t2m"], datillos.variables["longitude"][:], datillos.variables["latitude"][:]
data = temp[:]
datillos.close()

ristra_datos = []

for i in data:
	ristra_datos.append(float(i[0][0]-273.15))

meses = list(range(1940,2026))

pendiente_thei, ordenada_thei, cosa1, cosa2 = sps.mstats.theilslopes(ristra_datos, meses)
prediccion_theisen = np.array(meses)*pendiente_thei +ordenada_thei


pendiente_lin, ordenada_lin, cosa, cosa, cosa3 = linregress(meses, ristra_datos)
prediccion_lin = np.array(meses)*pendiente_lin + ordenada_lin

plt.scatter(meses, ristra_datos, s=10, c="red",label="Datos ERAS5")
plt.plot(meses,prediccion_theisen,c="blue", label=rf"Ajuste Thei-Sen $m=${pendiente_thei.round(4)}")

plt.plot(meses,prediccion_lin,c="green", label=rf"Ajuste Lineal $m=${pendiente_lin.round(4)}")
ticks = [1940,2026,10]
plt.set_ticks(ticks)

plt.legend()

plt.xlabel("Años")
plt.ylabel("Promedio temperatura (ºC)")
plt.title("Promedio temperatura en el mes Julio en Madrid en los años 1940-2025")

resultado_mannkendall = mk.original_test(ristra_datos)

print("pendiente Teil-Sen \t "+str(pendiente_thei),"\t intercept \t "+str(ordenada_thei))
print("pendiente regresion lin \t "+str(pendiente_lin),"\t intercept \t "+str(ordenada_lin))
print("pendiente Mann-Kendall \t "+str(resultado_mannkendall[-2]),"\t intercept \t "+str(resultado_mannkendall[-1]))

print(resultado_mannkendall)

plt.show()