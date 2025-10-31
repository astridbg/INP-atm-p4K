

def fix_cam_time(ds):
    # Author: Marte Sofie Buraas / Ada Gjermundsen    

    """ NorESM raw CAM h0 files has incorrect time variable output,
    thus it is necessary to use time boundaries to get the correct time
    If the time variable is not corrected, none of the functions involving time
    e.g. yearly_avg, seasonal_avg etc. will provide correct information

    Parameters
    ----------
    ds : xarray.DaraSet

    Returns
    -------
    ds_weighted : xarray.DaraSet with corrected time
    """
    from cftime import DatetimeNoLeap

    months = ds.time_bnds.isel(nbnd=0).dt.month.values
    years = ds.time_bnds.isel(nbnd=0).dt.year.values
    dates = [DatetimeNoLeap(year, month, 15) for year, month in zip(years, months)]
    ds = ds.assign_coords(time=dates)
    return ds

def computeWeightedMean(ds):

    # Author: Anne Fouilloux
    import numpy as np

    # Compute weights based on the xarray you pass
    weights = np.cos(np.deg2rad(ds.lat))
    weights.name = "weights"
    # Compute weighted mean
    air_weighted = ds.weighted(weights)
    weighted_mean = air_weighted.mean(("lon", "lat"))
    return weighted_mean


def polarCentral_set_latlim(lat_lims, ax):
    
    # Author: Anne Fouilloux
    import numpy as np
    import cartopy.crs as ccrs
    import matplotlib.path as mpath

    ax.set_extent([-180, 180, lat_lims[0], lat_lims[1]], ccrs.PlateCarree())
    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)
