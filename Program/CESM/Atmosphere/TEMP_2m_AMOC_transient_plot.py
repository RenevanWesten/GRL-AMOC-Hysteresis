#Program plots the 2-meter surface temperature and the AMOC response

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
from scipy import stats
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker


#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

	fh.close()

	return time, transport

def SignificantTrend(time, data):
	"""Finds whether trend is significant
	Returns the trend and if it significant (= 1)"""

	#Set time similar to Santer et al. (2000), time array from 1 till N
	#Statistical significance of trends and trend differences in layer-average atmospheric salterature time series
	time		= np.arange(1, len(time) + 1)

	#Determine the detrended time series
	trend, base 	= polyfit(time, data, 1)
	data_res	= data - ((trend * time) + base)

	#Effective sample size, based on the lag-1 correlation
	corr_1		= np.corrcoef(data_res[:-1], data_res[1:])[0, 1]
	N_eff		= int(len(time) * (1.0 - corr_1) / (1.0 + corr_1))

	#Determine the variance of the anomalies
	data_var	= np.sum(data_res**2.0) / (N_eff - 2.0)

	#Determine the standard error
	standard_error	=  np.sqrt(data_var) / np.sqrt(np.sum((time - np.mean(time))**2.0))

	#Determine the Student-T value
	t_value		= trend / standard_error

	#Get the significance levels and the corresponding critical values (two-sided)
	sig_levels 	= np.arange(50, 100, 0.5) / 100.0
	t_crit 		= stats.t.ppf((1.0 + sig_levels) / 2.0, N_eff - 2)

	#Get the indices where the significance is exceeding the critical values
	sig_index	= np.where(fabs(t_value) > t_crit)[0]
	significant	= 0.0

	if len(sig_index) > 0:
		#If there are significance values, take the highest significant level
		significant = sig_levels[sig_index[-1]]

	return trend, np.sqrt(standard_error), significant


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

fh 		= netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc', 'r')

time	    = fh.variables['time'][:]		
transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

fh.close()

fh 		    = netcdf.Dataset(directory+'Atmosphere/TEMP_2m_global.nc', 'r')

time		= fh.variables['time'][:] 		
temp_global	= fh.variables['TEMP_global'][:] #2-meter temperature (deg C)

fh.close()

#-----------------------------------------------------------------------------------------

trend_1, error_1, sig_1	= SignificantTrend(time[1749:1850], transport[1749:1850])
trend_2, error_2, sig_2	= SignificantTrend(time[4089:4190], transport[4089:4190])
trend_3, error_3, sig_3	= SignificantTrend(time[1749:1850], temp_global[1749:1850])
trend_4, error_4, sig_4	= SignificantTrend(time[4089:4190], temp_global[4089:4190])

print('AMOC ratio:', -trend_2 / trend_1)
print('TEMP ratio:', -trend_4 / trend_3)

#-----------------------------------------------------------------------------------------

trend_length	= 20

year_start	= np.arange(4090, 4191)
trend_period_1	= np.zeros(len(year_start))
trend_period_2	= np.zeros(len(year_start))

for year_i in range(len(year_start)):
	#Take 20 year trends
	trend, base 		= polyfit(time[1749+year_i:1749+trend_length+year_i], transport[1749+year_i:1749+trend_length+year_i], 1)
	trend_period_1[year_i]	= trend

	trend, base 		= polyfit(time[4089+year_i:4089+trend_length+year_i], transport[4089+year_i:4089+trend_length+year_i], 1)
	trend_period_2[year_i]	= trend

print()
print('Sliding window for AMOC trend')
print('Maximum trend (1750 - 1850):', np.max(-trend_period_1))
print('Maximum trend (4090 - 4190):', np.max(trend_period_2))

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

ax.fill_between([1750, 1850], -30, 30, alpha=0.25, edgecolor='orange', facecolor='orange')

graph_1	= ax.plot(time, temp_global - temp_global[1599], linestyle = '-', color = 'gray', linewidth = 1.5, label = 'Temperature (1600 - 2000)')

ax.set_xlabel('Model year')
ax.set_ylabel('Temperature difference ($^{\circ}$C)')
ax.set_xlim(1600, 2000)
ax.set_ylim(-2.5, 2.5)
ax.grid()

ax2 	 = ax.twiny()
graph_2  = ax2.plot(time, temp_global - temp_global[3939], linestyle = '-', color = 'firebrick', linewidth = 1.5, label = 'Temperature (3940 - 4340)')

ax2.set_xticks([3940, 3990, 4040, 4090, 4140, 4190, 4240, 4290, 4340])
ax2.set_xlim(3940, 4340)

for tl in ax2.get_xticklabels():
    tl.set_color('r')

ax3 	= ax.twinx()
graph_3	= ax3.plot(time, transport, linestyle = '-', color = 'k', linewidth = 1.5, label = 'AMOC (1600 - 2000)')
graph_4	= ax3.plot(time-2340, transport, linestyle = '-', color = 'r', linewidth = 1.5, label = 'AMOC (3940 - 4340)')

ax3.set_ylim(-2, 35)
ax3.set_ylabel('Volume transport (Sv)')

graphs	      = graph_1 + graph_2 + graph_3 + graph_4

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

fig.suptitle('a) Global mean surface temperature and AMOC strength')

show()
#-----------------------------------------------------------------------------------------
