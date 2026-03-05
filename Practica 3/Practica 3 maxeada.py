from netCDF4 import Dataset #Módulo para leer el archivo de datos
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt #Módulo para representar datos en gráficas
import cartopy.crs as ccrs #Módulo para el formato de mapa global


datillos = Dataset("Data_P3.nc", "r", format="NETCDF4") #Leemos los datos del archivo de de ERA5

#EXTRACCIÓN DE LAS DIFERENTES VARIABLES:

temp, lon, lat = datillos.variables["t2m"], datillos.variables["longitude"][:], datillos.variables["latitude"][:]
data = temp[:]

datos_media = data[1:31] #Guardamos las matrices que contienen los datos de referencia y los actuales por separado.
datos_2026 = data[31]

#CÁLCULO DE LA MEDIA DESDE 1961 HASTA 1990:

matriz_media = np.zeros_like(datos_2026)

for i in range(len(datos_media)):
    matriz_media = matriz_media + datos_media[i]

matriz_media = 1/(len(datos_media))*matriz_media

#CÁLCULO DE LA ANOMALÍA:

anomalia = datos_2026-matriz_media

datillos.close()

#CONSTRUCCIÓN DEL MAPA A COLORES:

ax = plt.axes(projection=ccrs.Robinson()) #Creamos una figura con la forma y los ejes correspondientes a un mapa global.
#También se puede en PlateCarree(mapa rectangular) o Orthographic(mapa circular).

levels = np.linspace(-8, 8, 100) #Creamos una discretización de puntos que servirá como base para los niveles de colores del mapa.
graf = plt.contourf(lon, lat, anomalia ,levels=levels ,transform=ccrs.PlateCarree(), cmap="coolwarm", extend='both') #Coloreamos el mapa asociando los valores de temperatura en cada punto del globo con una escala de colores.
#Con cmap se especifíca la paleta de colores, se puede cambiar por coolwarm, Spectral_r, etc.
#Con vmin y vmax se modifica rango del espectro de colores
#Con contourf se interpola entre los puntos en la figura para obtener un coloreado suave y continuo, pero se puedo usar
#plt.pcolormesh y cambiar 60 por shading="auto" y el programa no interpola representa los datos en bruto, cuadrado por cuadrado.

plt.title("Anomalía de Temperatura a 2 m (°C) – ERA5 – Enero 2026") #Título de la figura

ax.coastlines(resolution="110m") #Añadimos las líneas de costa
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False) #Adición de los meridianos y latitudes.

graf.set_clim(-8, 8) #Establece los valores mínimo y máximo para los colores límites en la escala(si no se añade esta línea se colocan los límites del levels por defecto).

ticks = list(np.arange(-8, 8, 2)) + [8] #Creamos una lista que indicará en que lugares de la barra se colocan las marcas qué indican el valor de la temperatura.
cbar = plt.colorbar(graf, orientation="horizontal", pad=0.05, shrink=0.8) #Creamos la barra que
#indica la escala de color con la temperatura, pad indica la separación al mapa y shrink el tamaño.
cbar.set_ticks(ticks)
cbar.set_label("Temperatura a 2 m (°C)") #Texto que acompaña a la barra
cbar.ax.tick_params(labelsize=8) #Modifica el tamaño de los múmeros en la barra de color

plt.show() #Mostramos la figura

