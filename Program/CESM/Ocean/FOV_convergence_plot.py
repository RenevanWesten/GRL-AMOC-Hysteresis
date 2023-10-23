#Program plots the freshwater transport at 34S and 60N

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinDataFOV(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_OV'][:]	#Fresh water

	fh.close()

	return time, FOV

def ReadinDataGyre(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_gyre'][:]	#Fresh water

	fh.close()

	return time, FOV

def ReadinDataFW(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['FW_transport'][:]	#Fresh water

	fh.close()

	return time, FOV


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time, FW_34S		= ReadinDataFW(directory+'Ocean/Freshwater_transport_section_34S.nc')
time, FOV_34S		= ReadinDataFOV(directory+'Ocean/FOV_index_section_34S.nc')
time, FW_gyre_34S 	= ReadinDataGyre(directory+'Ocean/FW_gyre_section_34S.nc')

time, FW_60N		= ReadinDataFW(directory+'Ocean/Freshwater_transport_section_60N.nc')
time, FOV_60N		= ReadinDataFOV(directory+'Ocean/FOV_index_section_60N.nc')
time, FW_gyre_60N 	= ReadinDataGyre(directory+'Ocean/FW_gyre_section_60N.nc')

#-----------------------------------------------------------------------------------------

year_average 	= 25

time_average    	    = int(len(time) / year_average)
time_average		    = np.zeros(time_average)
FW_34S_average     	    = np.zeros(len(time_average))
FOV_34S_average     	= np.zeros(len(time_average))
FW_gyre_34S_average     = np.zeros(len(time_average))

FW_60N_average     	    = np.zeros(len(time_average))
FOV_60N_average     	= np.zeros(len(time_average))
FW_gyre_60N_average     = np.zeros(len(time_average))



for time_i in range(len(time_average)):
	#Loop over each time slice
    time_average[time_i]		= np.mean(time[time_i*year_average:(time_i+1)*year_average])
    FW_34S_average[time_i] 		= np.mean(FW_34S[time_i*year_average:(time_i+1)*year_average])
    FOV_34S_average[time_i] 	= np.mean(FOV_34S[time_i*year_average:(time_i+1)*year_average])
    FW_gyre_34S_average[time_i] 	= np.mean(FW_gyre_34S[time_i*year_average:(time_i+1)*year_average])
    
    FW_60N_average[time_i] 		= np.mean(FW_60N[time_i*year_average:(time_i+1)*year_average])
    FOV_60N_average[time_i] 	= np.mean(FOV_60N[time_i*year_average:(time_i+1)*year_average])
    FW_gyre_60N_average[time_i] = np.mean(FW_gyre_60N[time_i*year_average:(time_i+1)*year_average])

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(4400 - time_average, FW_34S_average, linestyle = '-', color = 'r', linewidth = 2.0, label = r'$F_{\nabla\mathrm{S}}$')
ax.plot(4400 - time_average, FOV_34S_average, linestyle = '--', color = 'r', linewidth = 2.0, label = '$F_{\mathrm{ovS}}$')
ax.plot(4400 - time_average, FW_gyre_34S_average, linestyle = ':', color = 'r', linewidth = 2.0, label = '$F_{\mathrm{azS}}$')

graph_1	= ax.plot(time_average, FW_34S_average, linestyle = '-', color = 'k', linewidth = 2.0, label = r'$F_{\nabla\mathrm{S}}$')
graph_2	= ax.plot(time_average, FOV_34S_average, linestyle = '--', color = 'k', linewidth = 2.0, label = '$F_{\mathrm{ovS}}$')
graph_3	= ax.plot(time_average, FW_gyre_34S_average, linestyle = ':', color = 'k', linewidth = 2.0, label = '$F_{\mathrm{azS}}$')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-0.65, 0.65)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

graph_4		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 2, label = '1 - 2200')
graph_5		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 2, label = '2201 - 4400')

graphs	      	= graph_4 + graph_5
legend_labels = [l.get_label() for l in graphs]
legend_2      = ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.add_artist(legend_1)

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

ax.set_title('a) Freshwater transport at 34$^{\circ}$S')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(4400 - time_average, FW_60N_average, linestyle = '-', color = 'r', linewidth = 2.0, label = r'$F_{\nabla\mathrm{N}}$')
ax.plot(4400 - time_average, FOV_60N_average, linestyle = '--', color = 'r', linewidth = 2.0, label = '$F_{\mathrm{ovN}}$')
ax.plot(4400 - time_average, FW_gyre_60N_average, linestyle = ':', color = 'r', linewidth = 2.0, label = '$F_{\mathrm{azN}}$')

graph_1	= ax.plot(time_average, FW_60N_average, linestyle = '-', color = 'k', linewidth = 2.0, label = r'$F_{\nabla\mathrm{N}}$')
graph_2	= ax.plot(time_average, FOV_60N_average, linestyle = '--', color = 'k', linewidth = 2.0, label = '$F_{\mathrm{ovN}}$')
graph_3	= ax.plot(time_average, FW_gyre_60N_average, linestyle = ':', color = 'k', linewidth = 2.0, label = '$F_{\mathrm{azN}}$')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-0.65, 0.65)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

graph_4		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 2, label = '1 - 2200')
graph_5		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 2, label = '2201 - 4400')

graphs	      	= graph_4 + graph_5
legend_labels = [l.get_label() for l in graphs]
legend_2      = ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.add_artist(legend_1)

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

ax.set_title('b) Freshwater transport at 60$^{\circ}$N')

show()