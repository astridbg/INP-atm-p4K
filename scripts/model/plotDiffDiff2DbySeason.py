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
wpath="/projects/NS9600K/astridbg/INP-atm-p4K/figures/model/diffdiff_byseason/"

### PRESENT CLIMATE
# Default cases----------------
#case1 = "meyers92_20220210"; case1nm = "M92"
#case1 = "andenes21_20220222"; case1nm = "A21"
case1 = "M92_20240612"; case1nm = "M92"
# Modified cases---------------
#case2 = "meyers92_20220210"; case2nm = "M92"
#case2 = "andenes21_20220222"; case2nm = "A21"
case2 = "A21_20240612"; case2nm = "A21"

### +4K CLIMATE
# Default cases----------------
case1_4K = "M92_SSTp4K_20240611"; case1nm_4K = "M92_4K"
#case1_4K = "M92_SSTp4K_NoIce_20240611"; case1nm_4K = "M92_4K_NoIce"
# Modified cases---------------
case2_4K = "A21_SSTp4K_20240612"; case2nm_4K = "A21_4K_"
#case2_4K = "A21_SSTp4K_NoIce_20240617"; case2nm_4K = "A21_4K_NoIce"

#------------------------------	
date1 = "2007-04-15_2010-03-15"
date2 = "2007-04-15_2010-03-15"
#date2 = "2007-04-15_2007-12-15"

#------------------------------
# Add seasonal open ocean mask
#------------------------------
ocean_mask = True
if ocean_mask:
    ds_ocn = xr.open_dataset(rpath+"OCNFRAC"+"_"+case2+"_"+date2+".nc")
    open_sea = ds_ocn > 0.85
    open_sea = open_sea.groupby("time.season").mean("time")
    open_sea = open_sea >= 0.5

#------------------------------
# Two-dimensional fields
#------------------------------

variables = ["SWCFS","LWCFS","TGCLDIWP","TGCLDLWP", "CLDTOT", "CLDLOW", "CLDMED", "CLDHGH", "TREFHT", "NETCFS", "CLDLWEMS"]
colorbar_extents = [8, 8, 5, 15, 0.1, 0.1, 0.1, 0.1, 1.5, 8, 1]

#------------------------------
# Shaping and plotting fields
#------------------------------
for var, extent in zip(variables, colorbar_extents):
    print(var)
    ds1 = xr.open_dataset(rpath+var+"_"+case1+"_"+date1+".nc")
    ds2 = xr.open_dataset(rpath+var+"_"+case2+"_"+date2+".nc")
    ds1_4K = xr.open_dataset(rpath+var+"_"+case1_4K+"_"+date1+".nc")
    ds2_4K = xr.open_dataset(rpath+var+"_"+case2_4K+"_"+date2+".nc")

    # Get start and end date of period
    date_start = str(ds1.time[0].values).split(" ")[0]  
    date_end = str(ds1.time[-1].values).split(" ")[0]

    # Group cases by season and mean over the period by season
    ds1_seas = ds1.groupby("time.season").mean("time")  
    ds2_seas = ds2.groupby("time.season").mean("time")
    ds1_4K_seas = ds1_4K.groupby("time.season").mean("time")  
    ds2_4K_seas = ds2_4K.groupby("time.season").mean("time")

    diff = (ds2_4K_seas[var]-ds2_seas[var])-(ds1_4K_seas[var]-ds1_seas[var])

    lev_extent = round(max(abs(np.min(diff.sel(lat=slice(66.5,90)).values)), 
                            abs(np.max(diff.sel(lat=slice(66.5,90)).values))),2)
    print("Maximum absolute change: ", lev_extent, " Colorbar extent: ", extent)
    #if lev_extent < 0.004:
    #    lev_extent = 0.004
    #lev_extent = 10
    levels = np.linspace(-extent,extent,25)


    fig = plt.figure(1, figsize=[9,10],dpi=300)
    title = ds1[var].long_name+"\n"+case2nm+r"$-$"+case1nm
    #fig.suptitle(title, fontsize=22)
	
    # Set the projection to use for plotting
    ax1 = plt.subplot(2, 2, 1, projection=ccrs.Orthographic(0, 90))
    ax2 = plt.subplot(2, 2, 2, projection=ccrs.Orthographic(0, 90))
    ax3 = plt.subplot(2, 2, 3, projection=ccrs.Orthographic(0, 90))
    ax4 = plt.subplot(2, 2, 4, projection=ccrs.Orthographic(0, 90))
    #plt.subplots_adjust(top=0.85)

    for ax,season,label in zip([ax1, ax2, ax3, ax4], ["DJF", "MAM","JJA","SON"], ["(a)", "(b)", "(c)", "(d)"]):
    	
        functions.polarCentral_set_latlim([66.5,90], ax)
        map = diff.sel(season=season).plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), 
                                           cmap='coolwarm', levels=levels,
                                           add_colorbar=False)
        if ocean_mask:
            ax.contourf(open_sea.lon, open_sea.lat, open_sea["OCNFRAC"].sel(season=season), transform=ccrs.PlateCarree(), colors='none',hatches=['..'],levels=[.5, 1.5])
        ax.set_title(label+" "+season, fontsize=22)
        ax.coastlines()

	
    cb_ax = fig.add_axes([0.15, 0.07, 0.7, 0.04])

    cbar = plt.colorbar(map, cax=cb_ax, spacing = 'uniform', extend='both', orientation='horizontal', fraction=0.046, pad=0.06)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_xlabel(r"$\Delta$"+ds1[var].units, fontsize=18)
    
    # Customize number of decimal points for each colorbar extent
    if lev_extent >= 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places        
    elif 0.4 <= lev_extent < 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # One decimal place
    elif 0.04 <= lev_extent < 0.4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # Two decimal places     
    elif 0.004 <= lev_extent < 0.04:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}')) # Three decimal places
    
    plt.savefig(wpath+"pdf/"+var+"_"+case1+"_"+case2+"_"+case1_4K+"_"+case2_4K+".pdf", bbox_inches='tight')
    plt.savefig(wpath+"png/"+var+"_"+case1+"_"+case2+"_"+case1_4K+"_"+case2_4K+".png", bbox_inches='tight')

    plt.clf()
