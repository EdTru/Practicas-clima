from netCDF4 import Dataset #Módulo para leer el archivo de datos
import numpy as np
import matplotlib.pyplot as plt #Módulo para representar datos en gráficas

datillos = Dataset("Data_P5.nc", "r", format="NETCDF4") #Leemos los datos del archivo de de ERA5

temp, lon, lat = datillos.variables["t2m"], datillos.variables["longitude"][:], datillos.variables["latitude"][:]
data = temp[:]

print(data)
