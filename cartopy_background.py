import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from netCDF4 import Dataset
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Load your data
data = Dataset("../DATA/PPE_vars_accum.nc")
lon = data["lon"][:]
lat = data["lat"][:]
pres = data["p"][:]

extent = [-65, 45, 30, 85]

# Create figure and axis with PlateCarree projection
fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(extent, crs=ccrs.PlateCarree())
ax.coastlines(resolution='50m', color='black', linewidth=1)

# Set tick marks every 10 degrees for both longitude and latitude
ax.set_xticks(np.arange(extent[0], extent[1] + 1, 10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(extent[2], extent[3] + 1, 10), crs=ccrs.PlateCarree())

# Format tick labels as longitude and latitude
lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

# Plot your data with adjustable color limits
# (Here we use the min and max of the data; adjust vmin/vmax as needed)
vmin = pres.min()  # or set your own lower limit, e.g., 0
vmax = pres.max()  # or set your own upper limit, e.g., 100
sc = ax.scatter(lon, lat, s=0.1, c=pres, cmap='viridis', vmin=vmin, vmax=vmax)

bbox_ax = ax.get_position()
cbar_ax = fig.add_axes([0.94,bbox_ax.y0+0.01, 0.015, bbox_ax.y1-bbox_ax.y0 - 0.03])

# Add a colorbar with a label
cb = fig.colorbar(sc, cax=cbar_ax, orientation='vertical', pad=0.05)
cb.set_label('Pressure [hPa]')


plt.savefig("../PLOTS/lonlat.png",
        dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.close()
