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
import functions

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/INP-atm-present/figures/model/diff_2D/"

# Default cases----------------
#case1 = "def_20210126"; case1nm = "CAM6"
#case1 = "meyers92_20220210"; case1nm = "M92"
case1 = "M92_20240522"; case1nm = "M92"
#case1 = "andenes21_20220222"; case1nm = "A21"
# Modified cases---------------
#case2 = "meyers92_20220210"; case2nm = "M92"
#case2 = "andenes21_20220222"; case2nm = "A21"
#case2 = "A21_20240522"; case2nm = "A21"
#case2 = "andenes21_20240322"; case2nm = "A21_nohet"
#case2 = "andenes21_20240322_biggsoff"; case2nm = "A21_nobigg"
#case2 = "M92_4K_20240607"; case2nm = "M92_4K"
case2 = "M92_20240610"; case2nm = "M92_NEW"
#------------------------------	
date1 = "2007-04-15_2010-03-15"
#date1 = "2007-04-15"

#date2 = "2007-04-15_2010-03-15"
#date2 = "2007-04-15_2008-03-15"
#date2 = "2007-04-15_2007-12-15"
date2 = "2007-04-15"

#------------------------------
# Two-dimensional fields
#------------------------------

variables = ["SWCF","LWCF","SWCFS","LWCFS","NETCFS","CLDTOT","CLDHGH","CLDMED","CLDLOW","TGCLDIWP","TGCLDLWP","TREFHT"]
variables = ["TREFHT"]

#------------------------------
# Shaping and plotting fields
#------------------------------
for var in variables:
    print(var)
    ds1 = xr.open_dataset(rpath+var+"_"+case1+"_"+date1+".nc")
    ds2 = xr.open_dataset(rpath+var+"_"+case2+"_"+date2+".nc")

    # Get difference between cases time averaged over the whole period
    diff = ds2[var].isel(time=0)-ds1[var].isel(time=0)

    fig = plt.figure(1, figsize=[9,10],dpi=300)

    fig.suptitle(ds1[var].long_name+" "+case2nm+r"$-$"+case1nm, fontsize=26)
	
    lev_extent = round(max(abs(np.min(diff.sel(lat=slice(66.5,90)).values)), abs(np.max(diff.sel(lat=slice(66.5,90)).values))),2)
    if lev_extent < 0.004:
        lev_extent = 0.004
    levels = np.linspace(-lev_extent,lev_extent,25)
	
    # Set the projection to use for plotting
    ax = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(0, 90))
    functions.polarCentral_set_latlim([65,90], ax)
    map = diff.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), 
						cmap='coolwarm',levels=levels,
						add_colorbar=False)

    ax.coastlines()
	
#   cb_ax = fig.add_axes([0.15, 0.05, 0.7, 0.04])
    cb_ax = fig.add_axes([0.15, 0.07, 0.7, 0.04])

    cbar = plt.colorbar(map, cax=cb_ax, spacing = 'uniform', extend='both', orientation='horizontal', fraction=0.046, pad=0.06)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_xlabel(ds1[var].units, fontsize=23)

    if lev_extent >= 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places        
    elif 0.4 <= lev_extent < 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # One decimal place
    elif 0.04 <= lev_extent < 0.4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # Two decimal places     
    elif 0.004 <= lev_extent < 0.04:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}')) # Three decimal places

    plt.savefig(wpath+"pdf/quicktest"+var+"_"+case1+"_"+case2+".pdf",bbox_inches="tight")
    plt.savefig(wpath+"png/quicktest"+var+"_"+case1+"_"+case2+".png",bbox_inches="tight")	
    
    plt.clf()
