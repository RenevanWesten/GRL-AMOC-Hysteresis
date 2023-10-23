#Program plots the Atlantic meridional profiles between 250 - 500 m depths

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
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min	= 250
depth_max	= 500

#-----------------------------------------------------------------------------------------


fh = netcdf.Dataset(directory+'Ocean/Atlantic_'+str(depth_min)+'-'+str(depth_max)+'m_year_0001-2200.nc', 'r')

time_1	= fh.variables['time'][:]
lat	    = fh.variables['lat'][:]
lat_vol	= fh.variables['lat_volume'][:]	
temp_1	= fh.variables['TEMP'][:]	
salt_1	= fh.variables['SALT'][:] 		
dens_1	= fh.variables['PD'][:]	


fh.close()

fh = netcdf.Dataset(directory+'Ocean/Atlantic_'+str(depth_min)+'-'+str(depth_max)+'m_year_2200-4400.nc', 'r')

time_2	= fh.variables['time'][:]
temp_2	= fh.variables['TEMP'][:]	
salt_2	= fh.variables['SALT'][:] 		
dens_2	= fh.variables['PD'][:]	

fh.close()

#-----------------------------------------------------------------------------------------

lat_min_1	= -35
lat_max_1	= -30
lat_min_2	= 50
lat_max_2	= 55

#Get the indices
lat_min_index_1	= np.argmin(np.abs(lat - lat_min_1))
lat_max_index_1	= np.argmin(np.abs(lat - lat_max_1))+1
lat_min_index_2	= np.argmin(np.abs(lat - lat_min_2))
lat_max_index_2	= np.argmin(np.abs(lat - lat_max_2))+1

#Get the volumes to normalise the potential densities
volume_1	= lat_vol[lat_min_index_1:lat_max_index_1]
volume_2	= lat_vol[lat_min_index_2:lat_max_index_2]
volume_1	= volume_1 / np.sum(volume_1)
volume_2	= volume_2 / np.sum(volume_2)

dens_grad_1	= np.zeros(len(time_1))
dens_grad_2	= np.zeros(len(time_2))

for time_i in range(len(time_1)):
	#Take the meridional potential density difference
	dens_grad_1[time_i]	= np.sum(dens_1[time_i, lat_min_index_2:lat_max_index_2] * volume_2) - np.sum(dens_1[time_i, lat_min_index_1:lat_max_index_1] * volume_1)

for time_i in range(len(time_2)):
	#Take the meridional potential density difference
	dens_grad_2[time_i]	= np.sum(dens_2[time_i, lat_min_index_2:lat_max_index_2] * volume_2) - np.sum(dens_2[time_i, lat_min_index_1:lat_max_index_1] * volume_1)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(4400 - time_2, temp_2[:, 32], linestyle = '-', color = 'r', linewidth = 1.5, label = 'Temperature')
ax.plot(time_1, temp_1[:, 32], linestyle = '-', color = 'k', linewidth = 1.5, label = 'Temperature')

ax.set_xlabel('Model year')
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.set_xlim(1, 2200)
ax.set_ylim(10, 20)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 2, label = '1 - 2200')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 2, label = '2201 - 4400')

graphs	      	= graph_1 + graph_2
legend_labels = [l.get_label() for l in graphs]
legend        = ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)


ax.plot([4400-4168, 4400-4168], [10.8, 100], '--k', linewidth = 2.0)
ax.text(4400-4168, 10.1, 'AMOC$_{\mathrm{max}}$', verticalalignment='bottom', horizontalalignment='center', color = 'k', fontsize=12)

#-----------------------------------------------------------------------------------------

box1 = TextArea("1/", textprops=dict(color="k"))
box2 = TextArea("4400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(-0.023, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("500/", textprops=dict(color="k"))
box2 = TextArea("3900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.165, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("1000/", textprops=dict(color="k"))
box2 = TextArea("3400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.379, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("1500/", textprops=dict(color="k"))
box2 = TextArea("2900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.605, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("2000/", textprops=dict(color="k"))
box2 = TextArea("2400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.835, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

ax.set_title('b) Temperature (250 - 500 m) at 30$^{\circ}$N')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(4400 - time_2, salt_2[:, 32], linestyle = '-', color = 'r', linewidth = 1.5, label = 'Salinity')
ax.plot(time_1, salt_1[:, 32], linestyle = '-', color = 'k', linewidth = 1.5, label = 'Salinity')

ax.set_xlabel('Model year')
ax.set_ylabel('Salinity (g kg$^{-1}$)')
ax.set_xlim(1, 2200)
ax.set_ylim(32.5, 37.5)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 2, label = '1 - 2200')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 2, label = '2201 - 4400')

graphs	      	= graph_1 + graph_2
legend_labels = [l.get_label() for l in graphs]
legend        = ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)


ax.plot([4400-4168, 4400-4168], [32.85, 100], '--k', linewidth = 2.0)
ax.text(4400-4168, 32.55, 'AMOC$_{\mathrm{max}}$', verticalalignment='bottom', horizontalalignment='center', color = 'k', fontsize=12)

#-----------------------------------------------------------------------------------------

box1 = TextArea("1/", textprops=dict(color="k"))
box2 = TextArea("4400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(-0.023, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("500/", textprops=dict(color="k"))
box2 = TextArea("3900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.165, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("1000/", textprops=dict(color="k"))
box2 = TextArea("3400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.379, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("1500/", textprops=dict(color="k"))
box2 = TextArea("2900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.605, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("2000/", textprops=dict(color="k"))
box2 = TextArea("2400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.835, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

ax.set_title('d) Salinity (250 - 500 m) at 30$^{\circ}$N')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(4400 - time_2, dens_grad_2, linestyle = '-', color = 'r', linewidth = 1.5, label = 'Salinity')
ax.plot(time_1, dens_grad_1, linestyle = '-', color = 'k', linewidth = 1.5, label = 'Salinity')

ax.set_xlabel('Model year')
ax.set_ylabel('Potential density difference (kg m$^{-3}$)')
ax.set_xlim(1, 2200)
ax.set_ylim(0, 2)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 2, label = '1 - 2200')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 2, label = '2201 - 4400')

graphs	      	= graph_1 + graph_2
legend_labels = [l.get_label() for l in graphs]
legend        = ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)


ax.plot([4400-4168, 4400-4168], [0.14, 100], '--k', linewidth = 2.0)
ax.text(4400-4168, 0.02, 'AMOC$_{\mathrm{max}}$', verticalalignment='bottom', horizontalalignment='center', color = 'k', fontsize=12)

#-----------------------------------------------------------------------------------------

box1 = TextArea("1/", textprops=dict(color="k"))
box2 = TextArea("4400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(-0.023, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("500/", textprops=dict(color="k"))
box2 = TextArea("3900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.165, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("1000/", textprops=dict(color="k"))
box2 = TextArea("3400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.379, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("1500/", textprops=dict(color="k"))
box2 = TextArea("2900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.605, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("2000/", textprops=dict(color="k"))
box2 = TextArea("2400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.835, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

ax.set_title('Potential density (250 - 500 m) difference (north minus south)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CS	= contourf(lat, time_2, temp_2, levels = np.arange(0, 30.01, 1), extend = 'both', cmap = 'Spectral_r')
cbar	= colorbar(CS, ticks = np.arange(0, 30.01, 5))
cbar.set_label('Temperature ($^{\circ}$C)')

ax.set_xlim(-30, 62)
ax.set_xticks(np.arange(-20, 60.1, 20))
ax.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])
ax.set_ylabel('Model year')
ax.set_ylim(2200, 4400)
ax.set_yticks([2400, 2900, 3400, 3900, 4400])

ax.plot([-10, 100], [4168, 4168], '--k', linewidth = 2.0)
ax.text(-29, 4168, 'AMOC$_{\mathrm{max}}$', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=12)

ax.set_title('a) Temperature (250 - 500 m)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CS	= contourf(lat, time_2, salt_2, levels = np.arange(32.5, 37.51, 0.25), extend = 'both', cmap = 'BrBG_r')
cbar	= colorbar(CS, ticks = np.arange(32.5, 37.51, 0.5))
cbar.set_label('Salinity (g kg$^{-1}$)')

ax.set_xlim(-30, 62)
ax.set_xticks(np.arange(-20, 60.1, 20))
ax.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])
ax.set_ylabel('Model year')
ax.set_ylim(2200, 4400)
ax.set_yticks([2400, 2900, 3400, 3900, 4400])

ax.plot([-10, 100], [4168, 4168], '--k', linewidth = 2.0)
ax.text(-29, 4168, 'AMOC$_{\mathrm{max}}$', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=12)

ax.set_title('c) Salinity (250 - 500 m)')
#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

CS	= contourf(lat, time_2, dens_2, levels = np.arange(1024, 1028, 0.25), extend = 'both', cmap = 'PuOr_r')
cbar	= colorbar(CS, ticks = np.arange(1024, 1028.01, 1))
cbar.set_label('Potential density (kg m$^{-3}$)')

ax.set_xlim(-30, 62)
ax.set_xticks(np.arange(-20, 60.1, 20))
ax.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])
ax.set_ylabel('Model year')
ax.set_ylim(2200, 4400)
ax.set_yticks([2400, 2900, 3400, 3900, 4400])

ax.plot([-10, 100], [4168, 4168], '--k', linewidth = 2.0)
ax.text(-29, 4168, 'AMOC$_{\mathrm{max}}$', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=12)

ax.set_title('Potential density (250 - 500 m)')

show()



