#Time series of the Irminger basin

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

month_ice	= 3


#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Ocean/Irminger_basin_properties_month_'+str(month_ice)+'.nc', 'r')

time		= fh.variables['time'][:]
fraction	= fh.variables['Fraction'][:] 		
mixed		= fh.variables['MXL'][:] 			
heat_flux	= fh.variables['SHF'][:] 			

fh.close()


fh = netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_0-1000m.nc', 'r')

time		= fh.variables['time'][:]		
AMOC		= fh.variables['Transport'][:]	#AMOC strength

fh.close()


year_average 	= 5

time_average    	= int(len(time) / year_average)
time_average		= np.zeros(time_average)
fraction_average    = np.zeros(len(time_average))
mixed_average     	= np.zeros(len(time_average))
heat_flux_average   = np.zeros(len(time_average))

for time_i in range(len(time_average)):
	#Loop over each time slice
	time_average[time_i]		= np.mean(time[time_i*year_average:(time_i+1)*year_average])
	fraction_average[time_i] 	= np.mean(fraction[time_i*year_average:(time_i+1)*year_average])
	mixed_average[time_i] 		= np.mean(mixed[time_i*year_average:(time_i+1)*year_average])
	heat_flux_average[time_i] 	= np.mean(heat_flux[time_i*year_average:(time_i+1)*year_average])

#-----------------------------------------------------------------------------------------

fig, ax		= subplots()

graph_ice = ax.plot(time_average, fraction_average, '-b', linewidth = 1.0, label = 'Sea-ice fraction')

ax.set_xlim(1, 2200)
ax.set_ylim(0, 100)
ax.set_xlabel('Model year')
ax.set_ylabel('Sea-ice fraction ($\%$)')
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])

ax2 = ax.twinx()

graph_mixed = ax2.plot(time_average, mixed_average, '-k', linewidth = 1.0, label = 'Mixed layer depth')
graph_heat  = ax2.plot(time_average, -heat_flux_average, '-r', linewidth = 1.0, label = 'Surface heat flux')
ax2.set_ylabel('Mixed layer depth (m) and surface heat flux (W m$^{-2}$)')

ax2.set_ylim(800, 0)

graph_ice 	= ax.plot([-100,-100], [-100,-100], '-', color = 'b', linewidth = 2.0, label = 'Sea-ice fraction')
graph_mixed 	= ax.plot([-100,-100], [-100,-100], '-', color = 'k', linewidth = 2.0, label = 'Mixed layer depth')
graph_heat 	= ax.plot([-100,-100], [-100,-100], '-', color = 'r', linewidth = 2.0, label = 'Surface heat flux')
graphs	      	= graph_ice + graph_mixed + graph_heat

legend_labels 	= [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc = 'lower right', ncol=1, framealpha = 1.0)

ax.set_title('c) Irminger basin sea-ice fraction, MXL and SHF (March)')
#-----------------------------------------------------------------------------------------
fig, ax		= subplots()

graph_ice = ax.plot(time_average, fraction_average, '-b', linewidth = 1.0, label = 'Sea-ice fraction')

ax.set_xlim(2200, 4400)
ax.set_ylim(0, 100)
ax.set_xlabel('Model year')
ax.set_ylabel('Sea-ice fraction ($\%$)')
ax.grid()

ax.set_xticks([2400, 2900, 3400, 3900, 4400])

ax2 = ax.twinx()

graph_mixed = ax2.plot(time_average, mixed_average, '-k', linewidth = 1.0, label = 'Mixed layer depth')
graph_heat  = ax2.plot(time_average, -heat_flux_average, '-r', linewidth = 1.0, label = 'Surface heat flux')
ax2.set_ylabel('Mixed layer depth (m) and surface heat flux (W m$^{-2}$)')

ax2.set_ylim(800, 0)

#-----------------------------------------------------------------------------------------
ax3 	= fig.add_axes([0.28, 0.37, 0.40, 0.35])

ax3.plot([4091, 4091], [-5, 105], '--', linewidth = 2.0, color = 'dodgerblue')
ax3.plot(time, fraction, '-b', linewidth = 2.0, label = 'Sea-ice fraction')

ax3.set_xlabel('Model year')
ax3.set_ylabel('Sea-ice fraction ($\%$)')
ax3.set_ylim(0, 100)
ax3.set_xlim(4075, 4125)
ax3.grid()

ax3_2 = ax3.twinx()

graph_AMOC = ax3_2.plot(time, AMOC, '-', color = 'gray', linewidth = 2.0, label = 'AMOC')
ax3_2.set_ylabel('Volume transport (Sv)')

ax3_2.set_ylim(-2, 15)

ax3.set_title('Sea-ice fraction and AMOC strength')


graph_ice 	= ax.plot([-100,-100], [-100,-100], '-', color = 'b', linewidth = 2.0, label = 'Sea-ice fraction')
graph_mixed 	= ax.plot([-100,-100], [-100,-100], '-', color = 'k', linewidth = 2.0, label = 'Mixed layer depth')
graph_heat 	= ax.plot([-100,-100], [-100,-100], '-', color = 'r', linewidth = 2.0, label = 'Surface heat flux')
graph_AMOC 	= ax.plot([-100,-100], [-100,-100], '-', color = 'gray', linewidth = 2.0, label = 'AMOC')
graphs	      	= graph_ice + graph_mixed + graph_heat + graph_AMOC

legend_labels 	= [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc = 'lower left', ncol=1, framealpha = 1.0)


ax.set_title('d) Irminger basin sea-ice fraction, MXL and SHF (March)')

show()
