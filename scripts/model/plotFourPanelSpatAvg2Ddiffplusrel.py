import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
# Set font style to match latex document----------
#plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':18})
# ------------------------------------------------
import datetime
from functions import *

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/INP-atm-present/figures/model/spatavg_reldiff/"

# Default cases----------------
#case1 = "def_20210126"; case1nm = "CAM6"
case1 = "meyers92_20220210"; case1nm = "M92"
# Modified cases---------------
#case2 = "meyers92_20220210"; case2nm = "CAM5"
case2 = "andenes21_20220222"; case2nm = "A21"
#------------------------------	
date1 = "2007-04-15_2010-03-15"
date2 = "2007-04-15_2010-03-15"

#------------------------------
# Areas to analyse
#------------------------------

NYA = [78.9227, 11.9273] # Ny-Ålesund
ALE = [82.5, 297.6] # Alert
BAR = [71.17, 203.55] # Barrow
npole = [[0,360],[85,90]] # North Pole

#------------------------------
# Pairs of two-dimensional fields
#------------------------------

variables = ["SWCF","LWCF","SWCFS","LWCFS","CLDTOT","CLDHGH","CLDMED","CLDLOW","TGCLDIWP","TGCLDLWP","TREFHT", "NETCFS"]
variables = [["TGCLDIWP","TGCLDLWP"], ["SWCFS", "LWCFS"], ["NETCFS", "TREFHT"]]

#------------------------------
# Plotting monthly mean averages
#------------------------------
for pair in variables:
    
    ds1 = xr.open_mfdataset([rpath+pair[0]+"_"+case1+"_"+date1+".nc", rpath+pair[1]+"_"+case1+"_"+date1+".nc"])
    ds2 = xr.open_mfdataset([rpath+pair[0]+"_"+case2+"_"+date2+".nc", rpath+pair[1]+"_"+case2+"_"+date2+".nc"])
    
    # Get start and end date of period
    date_start = str(ds1.time[0].values).split(" ")[0]
    date_end = str(ds1.time[-1].values).split(" ")[0]

    # Get monthly mean
    ds1m = ds1.groupby("time.month").mean("time")
    ds2m = ds2.groupby("time.month").mean("time")
    
    # Make time array
    months = []
    for i in range(len(ds1m.month.values)):
        datetime_object = datetime.datetime.strptime(str(ds1m.month.values[i]), "%m")
        months.append(datetime_object.strftime("%b"))
    
    # Get spatial average over Arctic
    ds1_arct = computeWeightedMean(ds1m.sel(lat=slice(66.5,90)))
    ds2_arct = computeWeightedMean(ds2m.sel(lat=slice(66.5,90)))
    print("Arctic")
    print(pair[0])
    print((ds2_arct[pair[0]]-ds1_arct[pair[0]]).mean("month").values)
    print(pair[1])
    print((ds2_arct[pair[1]]-ds1_arct[pair[1]]).mean("month").values)

    # Get values over Ny-Ålesund
    ds1_NYA = ds1m.sel(lat=NYA[0], lon=NYA[1], method="nearest")
    ds2_NYA = ds2m.sel(lat=NYA[0], lon=NYA[1], method="nearest")

    # Get values over Alert
    ds1_ALE = ds1m.sel(lat=ALE[0], lon=ALE[1], method="nearest")
    ds2_ALE = ds2m.sel(lat=ALE[0], lon=ALE[1], method="nearest")

    # Get values over Barrow
    ds1_BAR = ds1m.sel(lat=BAR[0], lon=BAR[1], method="nearest")
    ds2_BAR = ds2m.sel(lat=BAR[0], lon=BAR[1], method="nearest")
        
    # Get spatial average over North Pole
    ds1_npol = computeWeightedMean(ds1m.sel(lon=slice(npole[0][0],npole[0][1]),
                                lat=slice(npole[1][0],npole[1][1])))
    ds2_npol = computeWeightedMean(ds2m.sel(lon=slice(npole[0][0],npole[0][1]),
                                lat=slice(npole[1][0],npole[1][1])))
    
    fig,axs = plt.subplots(ncols=2,nrows=2, gridspec_kw={'width_ratios': [3, 1]}, figsize=[13,7],dpi=300,constrained_layout=True)
    #fig.suptitle(ds1[var].long_name+", "+case2nm+r"$-$"+case1nm, fontsize=26)

    left_label = ["(a)", "(c)"]
    right_label = ["(b)", "(d)"]

    for i in range(2):
        axs[i, 0].plot(months, ds2_arct[pair[i]]-ds1_arct[pair[i]], color="cornflowerblue", label="Arctic")
        axs[i, 0].plot(months, ds2_NYA[pair[i]]-ds1_NYA[pair[i]], color="mediumseagreen", label="Ny-Ålesund")
        axs[i, 0].plot(months, ds2_ALE[pair[i]]-ds1_ALE[pair[i]], color="magenta", label="Alert")
        axs[i, 0].plot(months, ds2_BAR[pair[i]]-ds1_BAR[pair[i]], color="crimson", label="Utqiagvik")
        axs[i, 0].plot(months, ds2_npol[pair[i]]-ds1_npol[pair[i]], color="olive", label="North Pole")
        axs[i, 0].set_ylabel(r"$\Delta$"+ds1[pair[i]].units)
        axs[i, 0].tick_params(axis="x",labelsize=20)
        axs[i, 0].grid(alpha=0.5)
        axs[i, 0].annotate(left_label[i],fontsize=20,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 30), textcoords='offset points',
                ha='left', va='top')

        if pair[i] == "TREFHT":
            # Get average absolute change
            
            abs_arct = ds2_arct[pair[i]]-ds1_arct[pair[i]]
            abs_nya = ds2_NYA[pair[i]]-ds1_NYA[pair[i]]
            abs_ale = ds2_ALE[pair[i]]-ds1_ALE[pair[i]]
            abs_bar = ds2_BAR[pair[i]]-ds1_BAR[pair[i]]
            abs_npol = ds2_npol[pair[i]]-ds1_npol[pair[i]]
        
        
            change_all = pd.DataFrame({"Arctic":abs_arct,"Ny-Ålesund":abs_nya,"Alert":abs_ale,"Utqiagvik":abs_bar,"North Pole":abs_npol})
            bplot= axs[i, 1].boxplot(change_all,patch_artist=True,medianprops={"color":"black"})
            axs[i, 1].set_ylabel(r"$\Delta$"+ds1[pair[i]].units)
        
        elif pair[i] == "TGCLDLWP" or pair[i] == "CLDLWEM":

            # Get average relative change
            # Shrink y axis due to extreme values

            rel_arct = ((ds2_arct[pair[i]]-ds1_arct[pair[i]])/ds1_arct[pair[i]].where(ds1_arct[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_arct[pair[i]])
            rel_nya = ((ds2_NYA[pair[i]]-ds1_NYA[pair[i]])/ds1_NYA[pair[i]].where(ds1_NYA[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_NYA[pair[i]])
            rel_ale = ((ds2_ALE[pair[i]]-ds1_ALE[pair[i]])/ds1_ALE[pair[i]].where(ds1_ALE[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_ALE[pair[i]])
            rel_bar = ((ds2_BAR[pair[i]]-ds1_BAR[pair[i]])/ds1_BAR[pair[i]].where(ds1_BAR[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_BAR[pair[i]])
            rel_npol = ((ds2_npol[pair[i]]-ds1_npol[pair[i]])/ds1_npol[pair[i]].where(ds1_npol[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_npol[pair[i]])
        
        
            rel_all = pd.DataFrame({"Arctic":rel_arct,"Ny-Ålesund":rel_nya,"Alert":rel_ale,"Utqiagvik":rel_bar,"North Pole":rel_npol})
            bplot=axs[i, 1].boxplot(rel_all,patch_artist=True,medianprops={"color":"black"},showfliers=False)
            axs[i, 1].set_ylabel("% change")
            axs[i, 1].set_yscale("log")


        else:
            # Get average relative change
        
            rel_arct = ((ds2_arct[pair[i]]-ds1_arct[pair[i]])/ds1_arct[pair[i]].where(ds1_arct[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_arct[pair[i]])
            rel_nya = ((ds2_NYA[pair[i]]-ds1_NYA[pair[i]])/ds1_NYA[pair[i]].where(ds1_NYA[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_NYA[pair[i]])
            rel_ale = ((ds2_ALE[pair[i]]-ds1_ALE[pair[i]])/ds1_ALE[pair[i]].where(ds1_ALE[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_ALE[pair[i]])
            rel_bar = ((ds2_BAR[pair[i]]-ds1_BAR[pair[i]])/ds1_BAR[pair[i]].where(ds1_BAR[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_BAR[pair[i]])
            rel_npol = ((ds2_npol[pair[i]]-ds1_npol[pair[i]])/ds1_npol[pair[i]].where(ds1_npol[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_npol[pair[i]])
    
        
            rel_all = pd.DataFrame({"Arctic":rel_arct,"Ny-Ålesund":rel_nya,"Alert":rel_ale,"Utqiagvik":rel_bar,"North Pole":rel_npol})
            bplot=axs[i, 1].boxplot(rel_all,patch_artist=True,medianprops={"color":"black"},showfliers=False)
            axs[i, 1].set_ylabel("% change")
        
        colors=["cornflowerblue","mediumseagreen","magenta","crimson","olive"]
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
        axs[i, 1].grid(alpha=0.5)
        axs[i, 1].set_xticklabels([])
        axs[i, 1].annotate(right_label[i],fontsize=20,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 30), textcoords='offset points',
                ha='left', va='top')

        # Shrink current axis's height by 15% on the bottom
        box = axs[i, 0].get_position()
        axs[i, 0].set_position([box.x0, box.y0 + box.height * 0.15,
                    box.width, box.height * 0.85])
        box = axs[i, 1].get_position()
        axs[i, 1].set_position([box.x0, box.y0 + box.height * 0.15,
                    box.width, box.height * 0.85])

    # Put a legend below current axis
    axs[1, 0].legend(loc='upper center', bbox_to_anchor=(0.75, -0.12), ncol=5, columnspacing=0.5, handlelength=1,handletextpad=0.4)
	
    plt.grid(alpha=0.5)
    plt.savefig(wpath+"monthlymean/pdf/"+pair[0]+"_"+pair[1]+"_avg_"+case1+"_"+case2+".pdf",bbox_inches="tight")
    plt.savefig(wpath+"monthlymean/png/"+pair[0]+"_"+pair[1]+"_avg_"+case1+"_"+case2+".png",bbox_inches="tight")
    plt.clf()
