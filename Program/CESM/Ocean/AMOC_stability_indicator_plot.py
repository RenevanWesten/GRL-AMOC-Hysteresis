#Program plots the AMOC stability indicator

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from matplotlib.collections import LineCollection
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinDataFOV(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_OV'][:]	   #Freshwater (Sv)

	fh.close()

	return time, FOV

def ReadinDataAMOC(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	AMOC	= fh.variables['Transport'][:]	#AMOC strength

	fh.close()

	return time, AMOC

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

year_average 	= 50
#-----------------------------------------------------------------------------------------

time, FOV_34S	= ReadinDataFOV(directory+'Ocean/FOV_index_section_34S.nc')
time, FOV_60N	= ReadinDataFOV(directory+'Ocean/FOV_index_section_60N.nc')
time, AMOC	    = ReadinDataAMOC(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')

#Determine the freshwater convergence, only the AMOC component
FOV			        = FOV_34S - FOV_60N

time_average    	= int(len(time) / year_average)
time_average		= np.zeros(time_average)

FOV_average     	= np.zeros(len(time_average))
AMOC_average     	= np.zeros(len(time_average))


for time_i in range(len(time_average)):
	#Loop over each time slice
	time_average[time_i]	= np.mean(time[time_i*year_average:(time_i+1)*year_average])
	FOV_average[time_i] 	= np.mean(FOV[time_i*year_average:(time_i+1)*year_average])
	AMOC_average[time_i] 	= np.mean(AMOC[time_i*year_average:(time_i+1)*year_average])

#Determine the stability indicator
stability	= np.diff(FOV_average) / np.diff(AMOC_average)
time_average 	= time_average[:-1] + year_average / 2.0
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_2	= ax.plot(4400 - time_average, stability, linestyle = '-', color = 'r', linewidth = 2.0, label = '2201 - 4400')
graph_1	= ax.plot(time_average, stability, linestyle = '-', color = 'k', linewidth = 2.0, label = '1 - 2200')

ax.set_xlabel('Model year')
ax.set_ylabel('Stability indicator')
ax.set_xlim(1, 2200)
ax.set_ylim(-0.5, 0.5)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

graphs	      = graph_1 + graph_2
legend_labels = [l.get_label() for l in graphs]
legend        = ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

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

ax.set_title('d) AMOC stability indicator')

show()

#-----------------------------------------------------------------------------------------

