from netCDF4 import Dataset #Módulo para leer el archivo de datos
import numpy as np #Módulo para trabajar con matrices
import matplotlib
import matplotlib.pyplot as plt #Módulo para representar datos en gráficas
import cartopy.crs as ccrs #Módulo para el formato de mapa global


datillos = Dataset("P1/P1Data.nc", "r", format="NETCDF4") #Leemos los datos del archivo de de ERA5


temp, lon, lat = datillos.variables["t2m"], datillos.variables["longitude"][:], datillos.variables["latitude"][:]
data = temp[:]
#Extraemos los datos correspondientes a cada medida en arrays.


celsius = data-np.ones_like(data)*273.15 #Pasamos las medidas de temperatura a celsius. 
#Lo hacemos restando una matriz llena de 273.15 a la matriz que almacena los datos de temperatura.

datillos.close() #Cerramos el archivo de datos.


ax = plt.axes(projection=ccrs.Robinson()) #Creamos una figura con la forma y los ejes correspondientes a un mapa global.
#También se puede en PlateCarree(mapa rectangular) o Orthographic(mapa circular).

plt.contourf(lon, lat, celsius[0] ,60 ,transform=ccrs.PlateCarree(), cmap="turbo", vmin=-60, vmax=50) #Coloreamos el mapa asociando los valores de temperatura en cada punto del globo con una escala de colores.
#Con cmap se especifíca la paleta de colores, se puede cambiar por coolwarm, Spectral_r, etc.
#Con vmin y vmax se modifica rango del espectro de colores
#Con contourf se interpola entre los puntos en la figura para obtener un coloreado suave y continuo, pero se puedo usar 
#plt.pcolormesh y cambiar 60 por shading="auto" y el programa no interpola representa los datos en bruto, cuadrado por cuadrado. 

plt.title("Temperatura a 2 m (°C) – ERA5 – 2010/07/15 16:00 UTC") #Título de la figura

ax.coastlines(resolution="110m") #Añadimos las líneas de costa
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False) #Adición de los meridianos y latitudes.


cbar = plt.colorbar(orientation="horizontal", pad=0.05, shrink=0.8) #Creamos la barra que 
#indica la escala de color con la temperatura, pad indica la separación al mapa y shrink el tamaño.
cbar.set_label("Temperatura a 2 m (°C)") #Texto que acompaña a la barra
cbar.ax.tick_params(labelsize=8) #Modifica el tamaño de los múmeros en la barra de color


plt.show() #Mostramos la figura
