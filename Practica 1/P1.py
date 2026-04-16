from netCDF4 import Dataset 
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs 
import os 

os.chdir(os.path.dirname(os.path.realpath(__file__)))

datillos = Dataset("Data_P1.nc", "r", format="NETCDF4") 
temp, lon, lat = datillos.variables["t2m"], datillos.variables["longitude"][:], datillos.variables["latitude"][:]
data = temp[:]
datillos.close()

celsius = data-np.ones_like(data)*273.15 

ax = plt.axes(projection=ccrs.PlateCarree()) 

plt.contourf(lon, lat, celsius[0] ,60 ,transform=ccrs.PlateCarree(), cmap="turbo", vmin=-60, vmax=50) 


plt.title("Temperatura a 2 m (°C) – ERA5 – 2010/07/15 16:00 UTC") 

ax.coastlines(resolution="110m") 
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False) 


cbar = plt.colorbar(orientation="horizontal" ,pad=0.05, shrink=0.8)
cbar.set_label("Temperatura a 2 m (°C)") 
cbar.ax.tick_params(labelsize=8) 


plt.show() 
