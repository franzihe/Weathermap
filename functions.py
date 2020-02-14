# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---
# https://gist.github.com/ajdawson/dd536f786741e987ae4e
from copy import copy
import cartopy.crs as ccrs
import numpy as np
import shapely.geometry as sgeom
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import os



### Define colorbar colors
champ = 255.
# tot. precipitable water (grey scale)
no1 = np.array([255,255,255])/champ
no2 = np.array([231,231,231])/champ
no3 = np.array([201,201,201])/champ
no4 = np.array([171,171,171])/champ
no5 = np.array([140,140,140])/champ
no6 = np.array([110,110,110])/champ
no7 = np.array([80,80,80])/champ

# 250 hPa wind speed (colored scale)
no11 = np.array([255,255,255])/champ
no12 = np.array([196,225,255])/champ
no13 = np.array([131,158,255])/champ
no14 = np.array([255,209,177])/champ
no15 = np.array([255,118,86])/champ
no16 = np.array([239,102,178])/champ
no17 = np.array([243,0,146])/champ

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
        
        

def find_yx(lat, lon, point_lat, point_lon):
    abs_lat = abs(lat - point_lat)
    abs_lon = abs(lon - point_lon)

    c = np.maximum(abs_lat, abs_lon)

    y, x = np.where(c == c.min())
    y = y[0]
    x = x[0]
    
    xx = lat[y, x].x
    yy = lon[y, x].y
    return(xx, yy)


def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.
    
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    te = lambda xy: xy[0]
    lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])
    

def lambert_yticks(ax, ticks):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for an axis of a Lambert Conformal projection."""
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:    
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels


def plt_Jet_Thick_MSLP(fnx, u_250, mslp, Z_thickness, pw, andenes_x, andenes_y):
    
    projection = ccrs.LambertConformal(central_longitude =fnx.projection_lambert.longitude_of_central_meridian,
                                       central_latitude  =fnx.projection_lambert.latitude_of_projection_origin,
                                       standard_parallels = fnx.projection_lambert.standard_parallel)
    f, ax = plt.subplots(subplot_kw={'projection' : projection}, )
    
    #ax.set_title('Ensemble mean %s')# %s' %(_dm.time))
    ax.coastlines(resolution = '50m')

    
###################################################
    
    
    levels_u = np.arange(40,120,10)
    levels_th1 = np.arange(402,546,6)
    levels_th2 = np.arange(546,650,6)
    levels_p = np.arange(800,1100,4)
#    levels_pw = np.arange(22,78,8)
#    levels_pw = np.arange(14,62,8)
    levels_pw = np.arange(0,62,8)
###################################################
        
    # Plot contour lines for 250-hPa wind and fill
    U_map = colors.ListedColormap([no11, no12, no13, no14, no15, no16, no17])
    norm = colors.BoundaryNorm(boundaries = levels_u, ncolors=U_map.N)
    _U_250 = u_250.plot.pcolormesh(ax = ax,
                                       transform = projection,
                                       levels = levels_u,
                                       cmap = U_map,
                                       norm = norm,
                                       add_colorbar = False,
                                       extend = 'both'
                                      )
    
    cb_U_250 = plt.colorbar(_U_250, ax=ax, orientation="vertical",extend='both', 
                         shrink = 0.5)
    cb_U_250.set_label(label='250$\,$hPa Wind (m$\,$s$^{-1}$)', #size='large', 
                    weight='bold')
    
###################################################

    # Plot MSL pressure every 4 hPa
    CS_p = mslp.plot.contour(ax = ax,
                             transform = projection,
                             levels = levels_p,
                             colors = 'k',
                             linewidths = 1.8)
    ax.clabel(CS_p, levels_p[::2], 
              inline=1, fmt='%1.0f', #fontsize=10
             )
###################################################
    
    # Plot the 1000-500 hPa thickness
    CS_th1 = Z_thickness.where(Z_thickness < 546).plot.contour(ax = ax,
                                      transform = projection,
                                      levels = levels_th1,
                                      colors = 'b',
                                      linewidths = 2.,
                                      linestyles = '--')
    ax.clabel(CS_th1, levels_th1, inline = 1, fmt = '%1.0f')
    
    try:
        CS_th2 = Z_thickness.plot.contour(ax = ax,
                                      transform = projection,
                                      levels = levels_th2,
                                      colors = 'r',
                                      linewidths = 2.,
                                      linestyles = '--')
        ax.clabel(CS_th2, levels_th2, inline = 1, fmt = '%1.0f')
    except (ValueError):
        pass
#    labels = ['line1']
 #   CS_th1.collections[0].set_label(labels[0])
  #  plt.legend(bbox_to_anchor=(1.1, 1.05))
###################################################
    
    # Plot contourf for precipitable water
    PW_map = colors.ListedColormap([no1, no2, no3, no4, no5, no6, no7])
    PW_norm = colors.BoundaryNorm(boundaries = levels_pw, ncolors=PW_map.N)
    _PW = pw.plot.pcolormesh(ax = ax,
                                       transform = projection,
                                       levels = levels_u,
                                       cmap = PW_map,
                                       norm = PW_norm,
                                       add_colorbar = False,
                                       extend = 'both'
                                      )
    cb_PW = plt.colorbar(_PW, ax=ax, orientation="vertical",extend='both', 
                         shrink = 0.5)
    cb_PW.set_label(label='Precipitable water (m)', #size='large', 
                    weight='bold')
    
    
    
###################################################

    ax.plot([andenes_x], [andenes_y], color = 'red', marker = "^", transform=projection, markersize = 22 )
###################################################
    map_design(f,ax)
    
    
def plt_700_humidity(fnx, geop, temp, RH,andenes_x, andenes_y):
    projection = ccrs.LambertConformal(central_longitude =fnx.projection_lambert.longitude_of_central_meridian,
                                       central_latitude  =fnx.projection_lambert.latitude_of_projection_origin,
                                       standard_parallels = fnx.projection_lambert.standard_parallel)
    f, ax = plt.subplots(subplot_kw={'projection' : projection}, )
    
    #ax.set_title('Ensemble mean %s')# %s' %(_dm.time))
    ax.coastlines(resolution = '50m')

    

###################################################


    # Geopotential
    levels_g = np.arange(200,300,4)
    levels_T = np.arange(temp.min(), temp.max(),2)
    CS_g = geop.plot.contour(ax = ax, 
                             transform = projection, 
                             levels = levels_g,
                             colors = 'k',
                             linewidths=2)
    ax.clabel(CS_g, levels_g[::2], 
              inline=1, fmt='%1.0f', #fontsize=10
             )
###################################################
    # Temperature
    CS_T = temp.plot.contour(ax = ax,
                             transform = projection,
                             levels = levels_T,
                        #vmin=-40, vmax=10,
                        cmap=plt.get_cmap('plasma'),
                                                              linewidths=2,
                       # extend='both',
                        add_colorbar = False)
    ax.clabel(CS_T, inline = 1, fmt ='%1.0f',)
    cb_T = plt.colorbar(CS_T,ax = ax, format ='%1.0f', orientation = 'vertical', #extend = 'both',
                        shrink = 0.5)
    cb_T.lines[0].set_linewidth(25)
    cb_T.set_label(label='Temperature ($^{o}C$)', #size='large', 
                   weight='bold')
    #cb_T.ax.tick_params(labelsize='large')
###################################################

    # relative Humidity
    CS_rh = RH.where(RH >= 70.).plot.pcolormesh(ax = ax,
                                                transform = projection,
                                                levels = np.arange(70,110,10),
                                                cmap = plt.get_cmap('YlGn'),
                                               # extend = 'both',
                                                add_colorbar = False)

    cb_rh = plt.colorbar(CS_rh, orientation="vertical",#extend='both', 
                         shrink = 0.5)

    cb_rh.set_label(label='Relative humidity (%)', #size='large', 
                    weight='bold')
    #cb_rh.ax.tick_params(labelsize='large')

###################################################

    ax.plot([andenes_x], [andenes_y], color = 'red', marker = "^", transform=projection, markersize = 22 )

###################################################
    map_design(f,ax)

    
def plt_temp_wind_850(fnx, temp, u_wind, v_wind, XX, YY, andenes_x, andenes_y):
    
    projection = ccrs.LambertConformal(central_longitude =fnx.projection_lambert.longitude_of_central_meridian,
                                       central_latitude  =fnx.projection_lambert.latitude_of_projection_origin,
                                       standard_parallels = fnx.projection_lambert.standard_parallel)
    f, ax = plt.subplots(subplot_kw={'projection' : projection}, )
    
    #ax.set_title('Ensemble mean %s')# %s' %(_dm.time))
    ax.coastlines(resolution = '50m')

    

###################################################


    # Temperature
    levels_T = np.arange(-30, 31)

###################################################
    # Temperature
    CB_T = temp.plot.pcolormesh(ax = ax,
                             transform = projection,
                             levels = levels_T[::2],
                             cmap=plt.get_cmap('BrBG_r'),
                                                              
                        extend='both',
                        add_colorbar = False
                               )
    cb_T = plt.colorbar(CB_T,orientation="vertical",#extend='both', 
                         shrink = 0.5)
    cb_T.set_label(label='Temperature ($^{o}C$)', #size='large', 
                    weight='bold')
#####    
    CS_T = temp.plot.contour(ax = ax,
                      transform = projection,
                      levels = levels_T[::5],
                      colors = 'k',
                      linewidths = 2,
                            linestyles = '-')
    ax.clabel(CS_T, levels_T[::5], 
              inline=1, fmt='%1.0f', )
    
#####
    CS_T2 = temp.plot.contour(ax = ax,
                             transform = projection,
                             levels = np.array([ -38, -36, -34, -32,  -28, -26, -24, -22,  -18, -16,
       -14, -12,  -8,  -6,  -4,  -2,   2,   4,   6,   8,
        12,  14,  16,  18,  22,  24,  26,  28,  32,  34,  36,
        38, ]),
                             colors = 'lightgray',
                              linewidths = 1.2,
                             linestyles = '--')
###################################################
    # Wind barbs
    
    ax.barbs(XX[::70,::70], YY[::70,::70], u_wind[::70,::70], v_wind[::70,::70], barbcolor=[0.31372549, 0.31372549, 0.31764706])

###################################################

    ax.plot([andenes_x], [andenes_y], color = 'red', marker = "^", transform=projection, markersize = 22 )

    map_design(f,ax)
###################################################

def map_design(f, ax,):
# *must* call draw in order to get the axis boundary used to add ticks:
    f.canvas.draw()

    # Define gridline locations and draw the lines using cartopy's built-in gridliner:
    xticks = [-110, -50, -40, -30, -20, -11, 0, 10, 20, 30, 40, 50]
    yticks = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
    ax.gridlines(xlocs=xticks, ylocs=yticks)
    #ax.add_feature(cy.feature.OCEAN)
    #ax.add_feature(cy.feature.LAND)

    # Label the end-points of the gridlines using the custom tick makers:
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)

    lambert_xticks(ax, xticks)
    lambert_yticks(ax, yticks)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.tight_layout();   