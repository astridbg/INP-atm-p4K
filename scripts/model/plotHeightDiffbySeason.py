import xarray as xr
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
# Set font style to match latex document----------
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':20})
# ------------------------------------------------
import cartopy.crs as ccrs
import functions

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/INP-atm-present/figures/model/diff_byseason/"

# Default cases----------------
#case1 = "def_20210126"; case1nm = "CAM6"
#case1 = "meyers92_20220210"; case1nm = "M92"
case1 = "M92_20240612"; case1nm = "M92"
#case1 = "andenes21_20220222"; case1nm = "A21"

# Modified cases---------------
#case2 = "meyers92_20220210"; case2nm = "CAM5"
#case2 = "andenes21_20220222"; case2nm = "A21"
#case2 = "A21_20240522"; case2nm = "A21"
case2 = "M92_SSTp4K_20240611"; case2nm = "M92_4K"
#case2 = "A21_SSTp4K_20240612"; case2nm = "A21_4K"

#------------------------------	
date1 = "2007-04-15_2010-03-15"
date2 = "2007-04-15_2010-03-15"


#------------------------------
# Three-dimensional fields
# for specific level
#------------------------------

#variables = ["AWNI", "FREQI","CLDICE"]
#variables = ["NIMEY","AWNI", "FREQI","CLDICE"]
variables = ["ICIMR","ICWMR"]
xlims = [[-0.05,0.5],[-0.0002,0.0025],[-0.01,0.05]]

for var,xlim in zip(variables,xlims):
    print(var)
    ds1 = xr.open_dataset(rpath+var+"_"+case1+"_"+date1+".nc")
    ds2 = xr.open_dataset(rpath+var+"_"+case2+"_"+date2+".nc")
        
    # Get start and end date of period
    date_start = str(ds1.time[0].values).split(" ")[0]
    date_end = str(ds1.time[-1].values).split(" ")[0]

    # Group cases by season and mean over the period by season
    ds1m = ds1.groupby("time.season").mean("time")
    ds2m = ds2.groupby("time.season").mean("time")
	
    # Make horizontal averages:
    # - for the Arctic
    ds1_arct = functions.computeWeightedMean(ds1m[var].sel(lat=slice(66.5,90)))
    ds2_arct = functions.computeWeightedMean(ds2m[var].sel(lat=slice(66.5,90)))
	
    fig,axs = plt.subplots(nrows=1, ncols=4, sharey=True, figsize=[10,5],dpi=300, constrained_layout=True)
    fig.suptitle(ds1[var].long_name+", Arctic average", fontsize=22)
    
    if var == "TH":       	
        levels = ds1.ilev.values
        print(ds1.ilev)
    else:
        levels = ds1.lev.values
        
    for ax, season in zip(axs.flatten(),["DJF", "MAM","JJA","SON"]):
        
        if season == "DJF":
            ax.plot(ds1_arct.sel(season=season), levels, label=case1nm,ls="--")
            ax.plot(ds2_arct.sel(season=season), levels, label=case2nm)
            ax.set_ylabel("hPa")
            ax.invert_yaxis()
        else:
            ax.plot(ds1_arct.sel(season=season), levels, ls="--")
            ax.plot(ds2_arct.sel(season=season), levels)
        
        ax.set_title(season)
        #ax.set_xlim(xlim)
        ax.set_xlabel(ds1[var].units)
        ax.grid(alpha=0.5)
    fig.legend(loc="lower center", bbox_to_anchor=(0.5, -0.13),ncol=2)
    
    plt.savefig(wpath+var+"_heightdiff_arctic_"+case1+"_"+case2+".pdf",bbox_inches="tight")
    plt.clf()	

