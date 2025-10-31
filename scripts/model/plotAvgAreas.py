import xarray as xr
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
# Set font style to match latex document----------
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':16})
# ------------------------------------------------
import cartopy.crs as ccrs
import cartopy
from shapely.geometry.polygon import LinearRing
import matplotlib.patches as mpatches
import functions

wpath="/projects/NS9600K/astridbg/INP-atm-present/figures/model/"

#------------------------------
# Squares to plot
#------------------------------

# square = [lat, lon]

NYA = [78.63, 12.5] # Ny-Ålesund
ALE = [82.42, -62.5] # Alert
BAR = [71.05, -157.5] # Utqiavik/Barrow

coords = [NYA,ALE,BAR]

#------------------------------
# Plotting Arctic area
#------------------------------

fig = plt.figure(1, figsize=[5,5],dpi=300)

	
# Set the projection to use for plotting
proj = ccrs.Orthographic(central_longitude=0, central_latitude=90)

ax = plt.axes(projection=ccrs.Orthographic(0, 90))
ax.add_feature(cartopy.feature.OCEAN, zorder=0)
ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
functions.polarCentral_set_latlim([64,90], ax)
ax.coastlines()
ax.gridlines(color="darkblue",zorder=1)


# Add circle for the whole Arctic average

lat = 90
lon = 0
r = 23.5

def compute_radius(ortho, radius_degrees):
    phi1 = lat + radius_degrees if lat <= 0 else lat - radius_degrees
    _, y1 = ortho.transform_point(lon, phi1, ccrs.PlateCarree())
    return abs(y1)
r_ortho = compute_radius(proj,r)

ax.add_patch(mpatches.Circle(xy=[lon, lat], radius=r_ortho, facecolor=None, fill=False,edgecolor='red',transform=proj,lw=2,zorder=2))


# Add circle for the North Pole

lat = 90
lon = 0
r = 5

r_ortho = compute_radius(proj,r)

ax.add_patch(mpatches.Circle(xy=[lon, lat], radius=r_ortho, facecolor=None, fill=False,edgecolor='red',transform=proj,lw=2,zorder=2))

for coord in coords:
        lat = coord[0]; lat_buffer=0.95
        lon = coord[1]; lon_buffer=1.25
        LONS = [lon-lon_buffer, lon-lon_buffer, lon+lon_buffer, lon+lon_buffer]
        LATS = [lat-lat_buffer, lat+lat_buffer, lat+lat_buffer, lat-lat_buffer]
        ring = LinearRing(list(zip(LONS, LATS)))
        ax.add_geometries([ring], ccrs.PlateCarree(), facecolor='none', edgecolor='red', lw=2, zorder=2)
	
plt.savefig(wpath+"avgareas.pdf",bbox_inches="tight")
plt.savefig(wpath+"avgareas.png",bbox_inches="tight")
	
plt.clf()
