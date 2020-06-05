# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import pysftp
import xarray as xr
import numpy as np
#import cartopy as cy
import matplotlib.pyplot as plt

from datetime import datetime
#from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

import functions as fct
# %load_ext autoreload
# %autoreload 2

# + active=""
# srv = pysftp.Connection(host="sftp://franzihe@login.uio.no/uio/kant/geo-metos-u1/franzihe/", username="franzihe", private_key="/home/franzihe/.ssh/id_rsa")
# -

# Set figure size for all our plots
plt.rcParams['figure.figsize'] = [10., 8.]
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams.update({'font.size': 22})

# datetime object containing current date and time
# needed if used latest 
now = datetime.now()
# dd/mm/YY H:M:S
date = now.strftime("%Y%m%d")
time = now.strftime('%H')
print("date and time =", date, time)	


ini_time = '00'

# +
savefig = 0   # 1=yes, 0=no
form = 'png'

# File direction where figures should be saved (to be changed)
#figdir = 'sftp://franzihe@login.uio.no/uio/kant/geo-metos-u1/franzihe/www_docs'
figdir = '/home/franzihe/Documents/Figures/Weathermaps'

# +
# Pick out area around 400km to Andenes
# Andøya Space Center Coordinates: 69.2950N, 16.0300E
# to be changed to the location of interest, define the lower left latitude and longitude,
# upper right corner latitude and longitude

#andenes_lat = 69.2950; andenes_lon = 16.03
#lower_lat = 65.69; lower_lon = 5.8
#upper_lat = 72.9;  upper_lon = 26.26
#lower_lat = 65.; lower_lon = 1.24
#upper_lat = 75;  upper_lon = 26.26
# -

# Nordmela Coordinates: 69.1358N, 15.6776E
andenes_lat = 69.135840; andenes_lon = 15.677645
lower_lat = 64.83; lower_lon = 0.98
upper_lat = 74.85;  upper_lon = 25.85

# +
#  Open the netCDF file containing the input data.
# MEPS forecasts (Norway without Svalbard) from recent initialisations: https://thredds.met.no/thredds/catalog/mepslatest/catalog.html
# AROME Arctic forecast (Arctic including Svalbard) can be found in the 
# archive: https://thredds.met.no/thredds/catalog/aromearcticarchive/catalog.html
# latests: https://thredds.met.no/thredds/catalog/aromearcticlatest/catalog.html

#thredds = 'https://thredds.met.no/thredds/dodsC/mepslatest/meps_det_2_5km_%sT%sZ.ncml' %(date, ini_time)   # deterministic forecast
year = '2019'
month = '09'
day = '27'
date = year+month+day

thredds = 'https://thredds.met.no/thredds/dodsC/aromearcticarchive/%s/%s/%s/arome_arctic_extracted_2_5km_%sT%sZ.nc' %(year, month, day, date, ini_time)
fnx = xr.open_dataset(thredds, decode_times  = True, use_cftime = True)
# -

lower_x, lower_y = fct.find_yx(fnx.latitude, fnx.longitude, lower_lat, lower_lon)
upper_x, upper_y = fct.find_yx(fnx.latitude, fnx.longitude, upper_lat, upper_lon)
andenes_x, andenes_y = fct.find_yx(fnx.latitude, fnx.longitude, andenes_lat, andenes_lon)



# + active=""
# forecast_in_hours = '03'

# +
for forecast_in_hours in range(0,39,3):
    
    if forecast_in_hours < 10:
        forecast_in_hours = '0%s' %forecast_in_hours
#    print(int(forecast_in_hours))
    fig_name = 'meps_2_5_km_%sT%sZ_%s.%s' %(date,ini_time,forecast_in_hours, form)
#################################################################################
    map_area = 'Andoya'
    
#################################################################################
    ### 700hP - Temp - RH
    fct.plt_700_humidity(fnx, fnx.geopotential_pl.sel(pressure = 700.,
                                             x = slice(lower_x, upper_x),
                                             y = slice(lower_y, upper_y)).isel(time=int(forecast_in_hours))/100, 
                     fnx.air_temperature_pl.sel(pressure = 700.,
                                                x = slice(lower_x, upper_x),
                                                y = slice(lower_y, upper_y)).isel(time = int(forecast_in_hours)) - 273.15, 
                     fnx.relative_humidity_pl.sel(pressure = 700.,
                                                  x = slice(lower_x, upper_x),
                                                  y = slice(lower_y, upper_y)).isel(time = int(forecast_in_hours))*100, 
                     andenes_x, andenes_y)

    if savefig == 1:
        fct.createFolder('%s/700hPa_RH_T/%s/' %(figdir,map_area))
        plt.savefig('%s/700hPa_RH_T/%s/%s' %(figdir, map_area, fig_name), format = form, bbox_inches='tight', transparent=True)
        print('plot saved: %s/700hPa_RH_T/%s/%s' %(figdir, map_area, fig_name))
    plt.close()

#################################################################################    

    # convert Geopotential to height
    # https://en.wikipedia.org/wiki/Geopotential
    a = 6.378*10**6     # average radius of the earth  [m]
    G = 6.673*10**(-11) # gravitational constant       [Nm2/kg2]
    ma = 5.975*10**24   #  mass of the earth           [kg]

    Z_1000 = (-a**2 * fnx.geopotential_pl.sel(pressure = 1000.,).isel(time = int(forecast_in_hours)))/\
             (a * fnx.geopotential_pl.sel(pressure = 1000.,).isel(time = int(forecast_in_hours)) - G * ma)

    Z_500 = (-a**2 * fnx.geopotential_pl.sel(pressure = 500.,).isel(time = int(forecast_in_hours)))/\
             (a * fnx.geopotential_pl.sel(pressure = 500.,).isel(time = int(forecast_in_hours)) - G * ma)

    Z_thickness = (Z_500 - Z_1000)/10

    ########################################    
    ### Jet - Thickness - MSLP    
    fct.plt_Jet_Thick_MSLP(fnx, np.sqrt(fnx.x_wind_pl.sel(pressure = 250.,x = slice(lower_x, upper_x),y = slice(lower_y, upper_y)).isel(time = int(forecast_in_hours))**2 + \
                                        fnx.y_wind_pl.sel(pressure = 250.,x = slice(lower_x, upper_x),y = slice(lower_y, upper_y)).isel(time = int(forecast_in_hours))**2),
                       fnx.air_pressure_at_sea_level.sel(height_above_msl = 0, x = slice(lower_x, upper_x),y = slice(lower_y, upper_y)).isel(time = int(forecast_in_hours))/100, 
                       Z_thickness.sel(x = slice(lower_x, upper_x),y = slice(lower_y, upper_y)),
                       fnx.lwe_thickness_of_atmosphere_mass_content_of_water_vapor.sel(surface = 0, x = slice(lower_x, upper_x),y = slice(lower_y, upper_y)).isel(time = int(forecast_in_hours)),
                       andenes_x, andenes_y)
    if savefig == 1:
        fct.createFolder('%s/MSLP_Thickness_Jet/%s/' %(figdir,map_area))
        plt.savefig('%s/MSLP_Thickness_Jet/%s/%s' %(figdir, map_area, fig_name), format = form, bbox_inches='tight', transparent=True)
        print('plot saved: %s/MSLP_Thickness_Jet/%s/%s' %(figdir, map_area, fig_name))
    plt.close()
#################################################################################    

    ### 850hP - temperature - wind
    XX, YY = np.meshgrid(fnx.x.sel(x = slice(lower_x, upper_x)), 
                         fnx.y.sel(y = slice(lower_y, upper_y)))
    fct.plt_temp_wind_850(fnx, fnx.air_temperature_pl.sel(pressure = 850.,x = slice(lower_x, upper_x),y = slice(lower_y, upper_y)).isel(time = int(forecast_in_hours)) - 273.15,
                          fnx.x_wind_pl.sel(pressure = 850.,x = slice(lower_x, upper_x),y = slice(lower_y, upper_y)).isel(time = int(forecast_in_hours)),
                          fnx.y_wind_pl.sel(pressure = 850.,x = slice(lower_x, upper_x),y = slice(lower_y, upper_y)).isel(time = int(forecast_in_hours)),
                          XX, YY,
                          andenes_x, andenes_y)

    if savefig == 1:
        fct.createFolder('%s/850hPa_U_T/%s/' %(figdir,map_area))
        plt.savefig('%s/850hPa_U_T/%s/%s' %(figdir, map_area, fig_name), format = form, bbox_inches='tight', transparent=True)
        print('plot saved: %s/850hPa_U_T/%s/%s' %(figdir, map_area, fig_name))
    plt.close()
#################################################################################    


    map_area = 'Norway'
#################################################################################    

    ### 700hP - Temp - RH
    fct.plt_700_humidity(fnx, fnx.geopotential_pl.sel(pressure = 700.,).isel(time=int(forecast_in_hours))/100, 
                     fnx.air_temperature_pl.sel(pressure = 700.,).isel(time = int(forecast_in_hours)) - 273.15, 
                     fnx.relative_humidity_pl.sel(pressure = 700.,).isel(time = int(forecast_in_hours))*100, 
                     andenes_x, andenes_y)

    if savefig == 1:
        fct.createFolder('%s/700hPa_RH_T/%s/' %(figdir,map_area))
        plt.savefig('%s/700hPa_RH_T/%s/%s' %(figdir, map_area, fig_name), format = form, bbox_inches='tight', transparent=True)
        print('plot saved: %s/700hPa_RH_T/%s/%s' %(figdir, map_area, fig_name))
    plt.close()
#################################################################################    

    fct.plt_Jet_Thick_MSLP(fnx, np.sqrt(fnx.x_wind_pl.sel(pressure = 250.,).isel(time = int(forecast_in_hours))**2 + fnx.y_wind_pl.sel(pressure = 250.,).isel(time = int(forecast_in_hours))**2), 
                       fnx.air_pressure_at_sea_level.sel(height_above_msl = 0).isel(time = int(forecast_in_hours))/100, 
                       Z_thickness, 
                       fnx.lwe_thickness_of_atmosphere_mass_content_of_water_vapor.sel(surface = 0).isel(time = int(forecast_in_hours)), 
                       andenes_x, andenes_y)

    if savefig == 1:
        fct.createFolder('%s/MSLP_Thickness_Jet/%s/' %(figdir,map_area))
        plt.savefig('%s/MSLP_Thickness_Jet/%s/%s' %(figdir, map_area, fig_name), format = form, bbox_inches='tight', transparent=True)
        print('plot saved: %s/MSLP_Thickness_Jet/%s/%s' %(figdir, map_area, fig_name))
    plt.close()
#################################################################################

    XX, YY = np.meshgrid(fnx.x, fnx.y)
    fct.plt_temp_wind_850(fnx, fnx.air_temperature_pl.sel(pressure = 850,).isel(time = int(forecast_in_hours)) - 273.15,
                          fnx.x_wind_pl.sel(pressure = 850,).isel(time = int(forecast_in_hours)),
                          fnx.y_wind_pl.sel(pressure = 850,).isel(time = int(forecast_in_hours)),
                          XX, YY,
                          andenes_x, andenes_y)

    if savefig == 1:
        fct.createFolder('%s/850hPa_U_T/%s/' %(figdir,map_area))
        plt.savefig('%s/850hPa_U_T/%s/%s' %(figdir, map_area, fig_name), format = form, bbox_inches='tight', transparent=True)
        print('plot saved: %s/850hPa_U_T/%s/%s' %(figdir, map_area, fig_name))
    plt.close()
# -

fnx.close()
