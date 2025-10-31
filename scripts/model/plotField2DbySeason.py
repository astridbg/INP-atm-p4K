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
wpath="/projects/NS9600K/astridbg/INP-atm-present/figures/model/spatial/"

# Case ------------------------
#case = "def_20210126"; casenm = "CAM6"
case = "meyers92_20220210"; casenm = "M92"
#case = "andenes21_20220222"; casenm = "A21"
#case = "M92_20240522"; casenm="M92"
#------------------------------	
date = "2007-04-15_2010-03-15"
#------------------------------
# Add seasonal open ocean mask
#------------------------------
ocean_mask = False
if ocean_mask:
    ds_ocn = xr.open_dataset(rpath+"OCNFRAC"+"_"+case+"_"+date+".nc")
    open_sea = ds_ocn > 0.85
    open_sea = open_sea.groupby("time.season").mean("time")
    open_sea = open_sea >= 0.5


#------------------------------
# Two-dimensional fields
#------------------------------

variables = ["SWCF","LWCF","SWCFS","LWCFS","NETCFS","CLDTOT","CLDHGH","CLDMED","CLDLOW","TGCLDIWP","TGCLDLWP","TREFHT"]
variables = ["TGCLDIWP"]

#------------------------------
# Shaping and plotting fields
#------------------------------
for var in variables:
    print(var)
    ds = xr.open_dataset(rpath+var+"_"+case+"_"+date+".nc")

    # Get start and end date of period
    date_start = str(ds.time[0].values).split(" ")[0]  
    date_end = str(ds.time[-1].values).split(" ")[0]

    # Group cases by season and mean over the period by season
    ds_seas = ds.groupby("time.season").mean("time")

    lev_extent = round(np.max(ds_seas[var].sel(lat=slice(66.5,90)).values),10)
    print(lev_extent)
    if lev_extent < 0.004:
        lev_extent = 0.004
    levels = np.linspace(0,lev_extent,25)


    fig = plt.figure(1, figsize=[9,10],dpi=300)
    title = ds[var].long_name+"\n"+casenm+"\n"
    #fig.suptitle(title, fontsize=22)
	
    # Set the projection to use for plotting
    ax1 = plt.subplot(2, 2, 1, projection=ccrs.Orthographic(0, 90))
    ax2 = plt.subplot(2, 2, 2, projection=ccrs.Orthographic(0, 90))
    ax3 = plt.subplot(2, 2, 3, projection=ccrs.Orthographic(0, 90))
    ax4 = plt.subplot(2, 2, 4, projection=ccrs.Orthographic(0, 90))

    for ax,season in zip([ax1, ax2, ax3, ax4], ["DJF", "MAM","JJA","SON"]):
    	
        functions.polarCentral_set_latlim([66.5,90], ax)
        map = ds_seas[var].sel(season=season).plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), 
                                           cmap=plt.cm.get_cmap("Reds"),#cmap=plt.cm.get_cmap('Blues').reversed() 
                                           levels=levels,
                                           add_colorbar=False)

        if ocean_mask:
            ax.contourf(open_sea.lon, open_sea.lat, open_sea["OCNFRAC"].sel(season=season), transform=ccrs.PlateCarree(), colors='none',hatches=['..'],levels=[.5, 1.5])
        ax.set_title(season, fontsize=22)
        ax.coastlines()

	
    cb_ax = fig.add_axes([0.15, 0.07, 0.7, 0.04])

    cbar = plt.colorbar(map, cax=cb_ax, spacing = 'uniform', extend='both', orientation='horizontal', fraction=0.046, pad=0.06)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_xlabel(ds[var].units, fontsize=18)

    if lev_extent >= 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places        
    elif 0.4 <= lev_extent < 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # One decimal place
    elif 0.04 <= lev_extent < 0.4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # Two decimal places     
    elif 0.004 <= lev_extent < 0.04:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}')) # Three decimal places
    
    plt.savefig(wpath+"pdf/"+var+"_byseason_"+case+".pdf", bbox_inches='tight')
    plt.savefig(wpath+"png/"+var+"_byseason_"+case+".png", bbox_inches='tight')

    plt.clf()
