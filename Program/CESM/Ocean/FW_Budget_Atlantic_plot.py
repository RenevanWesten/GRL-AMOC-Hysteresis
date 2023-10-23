#Program plots the Atlantic freshwater content and budget

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Ocean/Freshwater_content_Atlantic.nc', 'r')

time		= fh.variables['time'][:]     
FW_0_1000	= fh.variables['FW_0-1000m'][:]    	
FW_1000_6000= fh.variables['FW_1000-6000m'][:]    
FW_0_6000	= fh.variables['FW_0-6000m'][:]    	
FW_surf		= fh.variables['FW_surface'][:]    				

fh.close()

#Convert the freshwater content to Sv
FW_Atlantic	    = ma.masked_all(len(time))
FW_Atlantic[:-1]= np.diff(FW_0_6000) / (86400.0 * 365.0) / 10**6.0

#-----------------------------------------------------------------------------------------

fh 			= netcdf.Dataset(directory+'Ocean/Freshwater_transport_section_34S.nc', 'r')

transport_salt_34S	= fh.variables['FW_transport'][:]		

fh.close()


fh 			= netcdf.Dataset(directory+'Ocean/Freshwater_transport_section_60N.nc', 'r')

transport_salt_60N	= fh.variables['FW_transport'][:]		

fh.close()

fh 			= netcdf.Dataset(directory+'Ocean/Freshwater_transport_Mediterranean.nc', 'r')

transport_salt_med	= fh.variables['FW_transport'][:]		

fh.close()


#-----------------------------------------------------------------------------------------

#Determine the freshwater transport into the Atlantic Ocean
#A positive transport means net freshwater through the section, thus 60N and the Mediterranean are subtracted
FW_transport	= transport_salt_34S - transport_salt_60N - transport_salt_med

FW_transport_diff	    = np.mean(FW_transport[1700:1750])-np.mean(FW_transport[0:50])
FW_transport_34S_diff	= np.mean(transport_salt_34S[1700:1750])-np.mean(transport_salt_34S[0:50])
FW_transport_60N_diff	= np.mean(transport_salt_60N[1700:1750])-np.mean(transport_salt_60N[0:50])
FW_transport_surf_diff	= np.mean(FW_surf[1700:1750])-np.mean(FW_surf[0:50])

#-----------------------------------------------------------------------------------------

year_average 		= 25

time_average    	= int(len(time) / year_average)
time_average		= ma.masked_all(time_average)
FW_average     		= ma.masked_all(len(time_average))
FW_surf_average		= ma.masked_all(len(time_average))
FW_transport_average	= ma.masked_all(len(time_average))

for time_i in range(len(time_average)):
	#Loop over each time slice
	time_average[time_i]		    = np.mean(time[time_i*year_average:(time_i+1)*year_average])
	FW_average[time_i] 		        = np.mean(FW_Atlantic[time_i*year_average:(time_i+1)*year_average])
	FW_surf_average[time_i] 	    = np.mean(FW_surf[time_i*year_average:(time_i+1)*year_average])
	FW_transport_average[time_i] 	= np.mean(FW_transport[time_i*year_average:(time_i+1)*year_average])

#-----------------------------------------------------------------------------------------

#Now determine the residual term
FW_mix		= FW_Atlantic - FW_transport - FW_surf

#Determine the freshwater content w.r.t. the first 50 years
FW_0_1000	= FW_0_1000 - np.mean(FW_0_1000[0:50])
FW_1000_6000= FW_1000_6000 - np.mean(FW_1000_6000[0:50])
FW_0_6000	= FW_0_6000 - np.mean(FW_0_6000[0:50])

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()


ax.plot(4400 - time_average, FW_average, linestyle = '-', color = 'r', linewidth = 2)
ax.plot(4400 - time_average, FW_surf_average, linestyle = '--', color = 'r', linewidth = 2)
ax.plot(4400 - time_average, FW_transport_average, linestyle = ':', color = 'r', linewidth = 2)

graph_1	= ax.plot(time_average, FW_average, linestyle = '-', color = 'k', linewidth = 2, label = r'$\frac{\mathrm{d}\overline{W}}{\mathrm{d}t}$')
graph_2	= ax.plot(time_average, FW_surf_average, linestyle = '--', color = 'k', linewidth = 2, label = '$F_{\mathrm{surf}}$')
graph_3	= ax.plot(time_average, FW_transport_average, linestyle = ':', color = 'k', linewidth = 2, label = r'$F_{\nabla}^b$')

ax.set_xlabel('Model year')
ax.set_ylabel(r'Freshwater transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-1, 1)
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

ax.set_title('d) Atlantic freshwater budget')

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
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()


ax.plot(4400 - time, FW_0_6000 / 10**14.0, linestyle = '-', color = 'r', linewidth = 2)
ax.plot(4400 - time, FW_0_1000 / 10**14.0, linestyle = '--', color = 'r', linewidth = 2)
ax.plot(4400 - time, FW_1000_6000 / 10**14.0, linestyle = ':', color = 'r', linewidth = 2)

graph_1	= ax.plot(time, FW_0_6000 / 10**14.0, linestyle = '-', color = 'k', linewidth = 2, label = '$\overline{W}$')
graph_2	= ax.plot(time, FW_0_1000 / 10**14.0, linestyle = '--', color = 'k', linewidth = 2, label = r'$\overline{W}_{1000\uparrow}$')
graph_3	= ax.plot(time, FW_1000_6000 / 10**14.0, linestyle = ':', color = 'k', linewidth = 2, label = r'$\overline{W}_{1000\downarrow}$')

ax.set_xlabel('Model year')
ax.set_ylabel(r'Freshwater content difference ($\times 10^{14}$ m$^3$)')
ax.set_xlim(1, 2200)
ax.set_ylim(-25, 25)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')


graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

graph_4		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 2, label = '1 - 2200')
graph_5		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 2, label = '2201 - 4400')

graphs	      	= graph_4 + graph_5
legend_labels = [l.get_label() for l in graphs]
legend_2      = ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.add_artist(legend_1)

ax.set_title('c) Atlantic freshwater content')

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
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Data/Atlantic_meridional_profiles/CESM_year_4001-4050.nc', 'r')

depth	= fh.variables['depth'][:] 		
lat	    = fh.variables['lat'][:] 		
salt	= np.mean(fh.variables['SALT'][:], axis = 0) 		

fh.close()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

depth_crop			        = 1000
factor_depth_crop		    = 4
depth[depth > depth_crop] 	= ((depth[depth > depth_crop] - depth_crop) / factor_depth_crop) + depth_crop

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.49, 0.163, 0.46, 0.25])

ax2.fill_between([-60, 70], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

CS	= ax2.contourf(lat, depth, salt, levels = np.arange(32.5, 37.51, 0.25), extend = 'both', cmap = 'BrBG_r')
cbar	= colorbar(CS, ticks = np.arange(33, 37.01, 1))
cbar.set_label('Salinity (g kg$^{-1}$)')

ax2.set_xlim(-30, 62)
ax2.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax2.set_ylabel('Depth (m)')	

ax2.set_xticks(np.arange(-20, 60.1, 20))
ax2.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

ax2.set_yticks([0, 500, 1000, 1500, 2000])
ax2.set_yticklabels([0, 500, 1000, 3000, 5000])

ax2.set_title('Salinity (4001 - 4050)', fontsize = 10)

show()

