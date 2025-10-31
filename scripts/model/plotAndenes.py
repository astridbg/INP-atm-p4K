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

wpath="/projects/NS9600K/astridbg/master/figures/model/"

#------------------------------
# Andenes coordinates
#------------------------------

lat = 69.18
lon = 16.07

#------------------------------
# Plotting Arctic area
#------------------------------

fig = plt.figure(1, figsize=[5,5],dpi=300)

	
# Set the projection to use for plotting

ax = plt.axes(projection=ccrs.Orthographic(0, 90))
ax.add_feature(cartopy.feature.OCEAN, zorder=0)
ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
functions.polarCentral_set_latlim([65,90], ax)
ax.coastlines()
ax.gridlines()

ax.scatter(lon, lat, c='r',s=100, transform=ccrs.PlateCarree())

plt.savefig(wpath+"andenes.pdf",bbox_inches="tight")
plt.savefig(wpath+"andenes.png",bbox_inches="tight",transparent=True)
	
plt.clf()
