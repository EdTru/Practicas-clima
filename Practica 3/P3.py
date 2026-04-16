from netCDF4 import Dataset 
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs 
import os 

os.chdir(os.path.dirname(os.path.realpath(__file__)))

datillos = Dataset("Data_P3.nc", "r", format="NETCDF4")
temp, lon, lat = datillos.variables["t2m"], datillos.variables["longitude"][:], datillos.variables["latitude"][:]


data = temp[:]
datillos.close()


datos_media = data[1:31] 
datos_2026 = data[31]

matriz_media = np.zeros_like(datos_2026)

for i in range(len(datos_media)):
	matriz_media = matriz_media + datos_media[i]

matriz_media = 1/(len(datos_media))*matriz_media

anomalia = datos_2026-matriz_media



ax = plt.axes(projection=ccrs.Robinson())

levels = np.linspace(-8, 8, 50)
graf = plt.contourf(lon, lat, anomalia ,levels=levels ,transform=ccrs.PlateCarree(), cmap="Pastel2", extend='both') 

plt.title("Anomalía de Temperatura a 2 m (°C) – ERA5 – Enero 2026")

ax.coastlines(resolution="110m")
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

graf.set_clim(-8, 8)

ticks = list(np.arange(-8, 8, 1)) + [8]
cbar = plt.colorbar(graf, orientation="horizontal", pad=0.05, shrink=0.8) 

cbar.set_ticks(ticks)
cbar.set_label("Temperatura a 2 m (°C)") 
cbar.ax.tick_params(labelsize=8) 

plt.show()

