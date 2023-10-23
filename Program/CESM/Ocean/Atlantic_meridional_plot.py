#Program plots the temperature, salinity and potential density in the Atlantic basin

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

year_start	= 3501
year_end	= 3550

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Data/Atlantic_meridional_profiles/CESM_year_'+str(year_start).zfill(4)+'-'+str(year_end).zfill(4)+'.nc', 'r')

depth	= fh.variables['depth'][:] 		
lat	    = fh.variables['lat'][:] 		
temp	= np.mean(fh.variables['TEMP'][:], axis = 0)	
salt	= np.mean(fh.variables['SALT'][:], axis = 0) 		
dens	= np.mean(fh.variables['PD'][:], axis = 0) 		

fh.close()

#-----------------------------------------------------------------------------------------

depth_crop			        = 1000
factor_depth_crop		    = 4
depth[depth > depth_crop] 	= ((depth[depth > depth_crop] - depth_crop) / factor_depth_crop) + depth_crop

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-60, 70], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

CS	= contourf(lat, depth, temp, levels = np.arange(0, 30.01, 1), extend = 'both', cmap = 'Spectral_r')
cbar	= colorbar(CS, ticks = np.arange(0, 30.01, 5))
cbar.set_label('Temperature ($^{\circ}$C)')

ax.set_xlim(-30, 62)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')	

ax.set_xticks(np.arange(-20, 60.1, 20))
ax.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

if year_start == 1 and year_end == 50:
	ax.set_title('a) Temperature (1 - 50)')

if year_start == 1701 and year_end == 1750:
	ax.set_title('d) Temperature (1701 - 1750)')

if year_start == 3101 and year_end == 3150:
	ax.set_title('g) Temperature (3101 - 3150)')

if year_start == 3501 and year_end == 3550:
	ax.set_title('j) Temperature (3501 - 3550)')

if year_start == 3901 and year_end == 3950:
	ax.set_title('m) Temperature (3901 - 3950)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-60, 70], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

CS	= contourf(lat, depth, salt, levels = np.arange(32.5, 37.51, 0.25), extend = 'both', cmap = 'BrBG_r')
cbar	= colorbar(CS, ticks = np.arange(32.5, 37.51, 0.5))
cbar.set_label('Salinity (g kg$^{-1}$)')

ax.set_xlim(-30, 62)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')	

ax.set_xticks(np.arange(-20, 60.1, 20))
ax.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

if year_start == 1 and year_end == 50:
	ax.set_title('b) Salinity (1 - 50)')

if year_start == 1701 and year_end == 1750:
	ax.set_title('e) Salinity (1701 - 1750)')

if year_start == 3101 and year_end == 3150:
	ax.set_title('h) Salinity (3101 - 3150)')

if year_start == 3501 and year_end == 3550:
	ax.set_title('k) Salinity (3501 - 3550)')

if year_start == 3901 and year_end == 3950:
	ax.set_title('n) Salinity (3901 - 3950)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-60, 70], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

CS	= contourf(lat, depth, dens, levels = np.arange(1024, 1028.01, 0.2), extend = 'both', cmap = 'PuOr_r')
CS_2	= contour(lat, depth, dens, levels = [1027], linewidths = 2.0, colors = 'k')
cbar	= colorbar(CS, ticks = np.arange(1024, 1028.01, 1))
cbar.set_label('Potential density (kg m$^{-3}$)')

ax.set_xlim(-30, 62)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')	

ax.set_xticks(np.arange(-20, 60.1, 20))
ax.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

if year_start == 1 and year_end == 50:
	ax.set_title('c) Potential density (1 - 50)')

if year_start == 1701 and year_end == 1750:
	ax.set_title('f) Potential density (1701 - 1750)')

if year_start == 3101 and year_end == 3150:
	ax.set_title('i) Potential density (3101 - 3150)')

if year_start == 3501 and year_end == 3550:
	ax.set_title('l) Potential density (3501 - 3550)')

if year_start == 3901 and year_end == 3950:
	ax.set_title('o) Potential density (3901 - 3950)')

show()
