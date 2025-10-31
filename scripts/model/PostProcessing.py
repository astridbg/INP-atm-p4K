import xarray as xr
import pandas as pd
import numpy as np
import glob
import functions

rpath="/projects/NS9600K/astridbg/data/model/noresm_rawdata/cases/"
wpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"

#case = "meyers92_20220210"
#case = "andenes21_20220222"
#case = "M92_20240612"
#case = "M92_SSTp4K_20240611"
#case = "M92_SSTp4K_NoIce_20240611"
#case = "A21_20240612"
#case = "A21_SSTp4K_20240612"
case = "A21_SSTp4K_NoIce_20240617"

casefolder="NF2000climo_f19_tn14_"+case

all_files = glob.glob(rpath+casefolder+"/atm/hist/"+casefolder+".cam.h0.*")
all_files.sort()
print("Files found")

ds = xr.open_mfdataset(all_files)
print("Dataset created")

#-----------------------------
# Postprocessing of model data
#-----------------------------

# Fix timestamp of model data
ds = functions.fix_cam_time(ds)

# Remove spinup months of data set
ds = ds.isel(time=slice(3,len(ds.time)))

#-----------------------------

print("Postprocessing completed")

#--------------------------------------------------------------------
# Store relevant variables intermediately to save time when plotting,
# change to desired units and create combined variables 
#--------------------------------------------------------------------
date = "2007-04-15_2010-03-15"
#date = "2007-04-15_2008-03-15"
#date = "2007-04-15_2007-12-15"
#date = "2007-04-15"


# For cases meyers92 and andenes21
#variables = ["NIMEY","AWNI", "FREQI","CLDICE","SWCF","LWCF","LWCFS","SWCFS","NETCFS","AWNICC","CLDTOT","CLDHGH","CLDMED","CLDLOW","TGCLDIWP","TGCLDLWP","TREFHT"]

# For case def
variables = ["AWNI", "FREQI","NUMICE","LWCFS","SWCFS","NETCFS","CLDTOT","CLDHGH","CLDMED","CLDLOW","CLOUD","TGCLDIWP","TGCLDLWP","TREFHT", "CLDLWEMS", "OCNFRAC"]
#variables = ["OCNFRAC"]
#variables = ["FSNT","FSNTC","FSNTOA","FSNTOAC","FSUTOA", "FLNT", "FLNTC", "FLUT", "FLUTC"]

for var in variables:
    print("Started writing variable:")
	
    # Change to desired units
    if var == "NIMEY" or var == "AWNI":
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre and name
        ds[var].attrs["units"] = "1/L"

    if var == "T" or var == "TREFHT" or var=="TH":
        ds[var].values = ds[var].values - 273.15 # Change unit to degrees Celsius
        ds[var].attrs["units"] = r"$^{\circ}$C"

    if var == "CLDICE" or var == "CLDLIQ" or var == "Q" or var == "ICIMR" or var == "ICWMR": 
        ds[var].values = ds[var].values*1e+3 # Change unit to grams per kilograms
        ds[var].attrs["units"] = "g/kg"
	
    if var == "TGCLDIWP" or var == "TGCLDLWP":
        ds[var].values = ds[var].values*1e+3 # Change unit to grams per squared meter
        ds[var].attrs["units"] = "g/m$^2$"

    if var == "IWC":
        ds[var].values = ds[var].values*1e+3 # Change unit to grams per cubic meter
        ds[var].attrs["units"] = "g/m$^3$"

    # Change to more meaningful name

    if var == "TREFHT":
        ds[var].attrs["long_name"] = "Surface (2m) Temperature"
    if var == "NIMEY":
        ds[var].attrs["long_name"] = "Activated Ice Number Concentation due to Meyers' parameterisation"
    if var == "AWNI":
        ds[var].attrs["long_name"] = "Average cloud ice number concentration"
    if var == "ICIMR":
        ds[var].attrs["long_name"] = "In-cloud ice mixing ratio"
    if var == "ICWMR":
        ds[var].attrs["long_name"] = "In-cloud liquid water mixing ratio"

    # Make combined data variables
    if var == "LWCFS":
        ds = ds.assign(LWCFS=ds["FLNSC"]-ds["FLNS"])
        ds[var].attrs["units"] = "W/m$^2$"
        ds[var].attrs["long_name"] = "Longwave cloud radiative effect at surface"

    if var == "SWCFS":
        ds = ds.assign(SWCFS=ds["FSNS"]-ds["FSNSC"])
        ds[var].attrs["units"] = "W/m$^2$"
        ds[var].attrs["long_name"] = "Shortwave cloud radiative effect at surface"
	
    if var == "AWNINONIMEY":
        ds = ds.assign(AWNINONIMEY=ds["AWNI"]-ds["NIMEY"])
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Average cloud ice number concentration minus Meyer's contribution"
    
    if var == "AWNICC":
        ds = ds.assign(AWNICC=ds["AWNI"]/ds["FREQI"].where(ds["FREQI"]>0))
        ds[var] = ds[var].fillna(0)
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Average cloud ice number concentration in cold clouds"

    if var == "AWNICLD":
        ds = ds.assign(AWNICLD=ds["AWNI"]/ds["CLOUD"].where(ds["CLOUD"]>0))
        ds[var] = ds[var].fillna(0)
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Average cloud ice number concentration in clouds"
    
    if var == "NIMEYCC":
        ds = ds.assign(NIMEYCC=ds["NIMEY"]/ds["FREQI"].where(ds["FREQI"]>0))
        ds[var] = ds[var].fillna(0)
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Activated Ice Number Concentation due to Meyers' parameterisation in cold clouds"

    if var == "NIMEYCLD":
        ds = ds.assign(NIMEYCLD=ds["NIMEY"]/ds["CLOUD"].where(ds["CLOUD"]>0))
        ds[var] = ds[var].fillna(0)
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Activated Ice Number Concentation due to Meyers' parameterisation in clouds"
        
    if var == "NETCFS":
        ds = ds.assign(NETCFS=ds["FSNS"]-ds["FSNSC"]-ds["FLNS"]-(-ds["FLNSC"]))
        ds[var].attrs["units"] = "W/m$^2$"
        ds[var].attrs["long_name"] = "Net cloud radiative effect at surface"
    
    if var == "CLDLWEM":
        ds = ds.assign(CLDLWEM=1-np.exp(-0.13*ds["TGCLDLWP"]*1e+3))
        ds[var].attrs["units"] = "Emissivity"
        ds[var].attrs["long_name"] = "Cloud Longwave Emissivity"

    if var == "CLDLWEMS":
        ds = ds.assign(CLDLWEMS=1-np.exp(-0.158*ds["TGCLDLWP"]*1e+3))
        ds[var].attrs["units"] = "Emissivity"
        ds[var].attrs["long_name"] = "Cloud Longwave Emissivity"

    print(ds[var].attrs["long_name"])
    print("Units: ", ds[var].attrs["units"])
	
    ds[var].to_netcdf(wpath+var+"_"+case+"_"+date+".nc")
