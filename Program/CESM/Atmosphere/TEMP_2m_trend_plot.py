#Program plots the 2-meter surface temperature trend

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh 		= netcdf.Dataset(directory+'Atmosphere/TEMP_2m_trend_year_1750-1850.nc', 'r')

lon			        = fh.variables['lon'][:]
lat			        = fh.variables['lat'][:] 			
temp_trend_1		= fh.variables['TEMP_trend'][:] 
temp_trend_1_sig	= fh.variables['TEMP_trend_sig'][:]

fh.close()

#-----------------------------------------------------------------------------------------
fh 			= netcdf.Dataset(directory+'Atmosphere/TEMP_2m_trend_year_4090-4190.nc', 'r')
			
temp_trend_2		= fh.variables['TEMP_trend'][:] 
temp_trend_2_sig	= fh.variables['TEMP_trend_sig'][:]

fh.close()

#-----------------------------------------------------------------------------------------

temp_ratio		= -temp_trend_2 / temp_trend_1

#-----------------------------------------------------------------------------------------
#Rescale the temperature plot
scale	= 5.0
cut_off	= 5

temp_trend_1[temp_trend_1 < -cut_off]			= (temp_trend_1[temp_trend_1 < -cut_off] - -cut_off) / scale - cut_off
temp_trend_1[temp_trend_1 > cut_off]			= (temp_trend_1[temp_trend_1 > cut_off] - cut_off) / scale + cut_off
temp_trend_2[temp_trend_2 < -cut_off]			= (temp_trend_2[temp_trend_2 < -cut_off] - -cut_off) / scale - cut_off
temp_trend_2[temp_trend_2 > cut_off]			= (temp_trend_2[temp_trend_2 > cut_off] - cut_off) / scale + cut_off

#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS      = ax.contourf(lon, lat, temp_trend_1, levels = np.arange(-8, 8.01, 0.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [-8, -6, -4, -2, 0, 2, 4, 6, 8], cax=ax_cb)
cbar.ax.set_yticklabels([-20, -10, -4, -2, 0, 2, 4, 10, 20])
cbar.set_label('2-meter temperature trend ($^{\circ}$C per century)')

ax.set_global()
ax.gridlines()
ax.coastlines()

for lat_i in range(0, len(lat), 3):
	for lon_i in range(0, len(lon), 3):
		#Determine significant difference

		if temp_trend_1_sig[lat_i, lon_i] <= 0.95:
			#Non-significant difference
			ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

ax.set_title('Yearly 2-meter temperature trend (1750 - 1850)')
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS      = ax.contourf(lon, lat, temp_trend_2, levels = np.arange(-8, 8.01, 0.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [-8, -6, -4, -2, 0, 2, 4, 6, 8], cax=ax_cb)
cbar.ax.set_yticklabels([-20, -10, -4, -2, 0, 2, 4, 10, 20])
cbar.set_label('2-meter temperature trend ($^{\circ}$C per century)')

ax.set_global()
ax.gridlines()
ax.coastlines()

for lat_i in range(0, len(lat), 3):
	for lon_i in range(0, len(lon), 3):
		#Determine sig2ificant difference

		if temp_trend_1_sig[lat_i, lon_i] > 0.95 and temp_trend_2_sig[lat_i, lon_i] > 0.95 and temp_ratio[lat_i, lon_i] >= 1.5:
			#Significant difference and 1.5 faster
			ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

		if temp_trend_1_sig[lat_i, lon_i] > 0.95 and temp_trend_2_sig[lat_i, lon_i] > 0.95 and temp_ratio[lat_i, lon_i] <= 0.5:
			#Significant difference and 0.5 slower
			ax.scatter(lon[lon_i], lat[lat_i], marker = 's', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

ax.set_title('b) 2-meter temperature trend (4090 - 4190)')

show()

