#Program plots one particular depth level of AMOC

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy.interpolate import interp1d

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinDataAMOC(filename):

    fh = netcdf.Dataset(filename, 'r')

    year	= fh.variables['time'][:] 		
    depth	= fh.variables['depth'][:] 		
    lat	    = fh.variables['lat'][:] 		
    AMOC	= fh.variables['AMOC'][:]	

    fh.close()

    return year, lat, depth, AMOC

def ReadinDataMHT(filename):

    fh = netcdf.Dataset(filename, 'r')

    year	= fh.variables['time'][:] 		
    lat	    = fh.variables['lat'][:] 		
    MHT  	= fh.variables['MHT_Atlantic'][:]	

    fh.close()

    return year, lat, MHT
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_level = 500

files = glob.glob(directory+'Data/AMOC_structure/CESM_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

year, lat, depth, AMOC = ReadinDataAMOC(files[0])

#-----------------------------------------------------------------------------------------

time_all        = np.zeros(len(files) * len(year))
AMOC_all		= ma.masked_all((len(time_all), len(lat)))


for file_i in range(len(files)):
    #Now determine for AMOC strength
    print(file_i)
	    
    year, lat, depth, AMOC = ReadinDataAMOC(files[file_i])

    for year_i in range(len(year)):
        #Determine the meridional transport
        time_all[file_i*len(year)+year_i] = year[year_i]

        for lat_i in range(len(lat)):
            #Interpolate to given grid	
            AMOC_all[file_i*len(year)+year_i, lat_i]	= interp1d(depth, AMOC[year_i, :, lat_i])(depth_level)

#-----------------------------------------------------------------------------------------
#Read in the meridional heat transport files
files = glob.glob(directory+'Data/AMOC_MHT/CESM_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

year, lat_MHT, MHT = ReadinDataMHT(files[0])

#-----------------------------------------------------------------------------------------

MHT_all		= ma.masked_all((len(time_all), len(lat_MHT)))

for file_i in range(len(files)):
    #Now determine for AMOC strength
    print(file_i)
	    
    year, lat_MHT, MHT = ReadinDataMHT(files[file_i])

    for year_i in range(len(year)):
        #Save the meridional heat transport
        MHT_all[file_i*len(year)+year_i]    = MHT[year_i]


year_average		= 25
time_average    	= int(len(time_all) / year_average)
time_average		= np.zeros(time_average)
MHT_average     	= np.zeros((len(time_average), len(lat_MHT)))

for time_i in range(len(time_average)):
    #Loop over each time slice
    time_average[time_i]		= np.mean(time_all[time_i*year_average:(time_i+1)*year_average])
    MHT_average[time_i] 		= np.mean(MHT_all[time_i*year_average:(time_i+1)*year_average], axis = 0)

#Read in one file from a particular year
fh = netcdf.Dataset(directory+'/Data/AMOC_MHT/CESM_year_4101-4150.nc', 'r')

MHT_atlantic_4101_4150	= np.mean(fh.variables['MHT_Atlantic'][:], axis = 0)	#Meridional heat transport (PW)

fh.close()

#-----------------------------------------------------------------------------------------    

#Rescale the AMOC plot
scale	= 3
cut_off	= 3

AMOC_all[AMOC_all < -cut_off]	= (AMOC_all[AMOC_all < -cut_off] - -cut_off) / scale - cut_off
AMOC_all[AMOC_all > cut_off]	= (AMOC_all[AMOC_all > cut_off] - cut_off) / scale + cut_off
#-----------------------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1], 'hspace': 0.05})

CS	= ax1.contourf(lat, time_all, AMOC_all, levels = np.arange(-9, 9.01, 0.5), extend = 'both', cmap = 'RdBu_r')

CS_2	= ax1.contour(lat_MHT, time_average, MHT_average, levels = [0.25], colors = 'royalblue', linewidths = 2)
CS_2	= ax1.contour(lat_MHT, time_average, MHT_average, levels = [1.0], colors = 'cyan', linewidths = 2)

fig.subplots_adjust(right=0.8)
cbar 	= fig.add_axes([0.82, 0.278, 0.030, 0.60])
fig.colorbar(CS, cax = cbar, ticks = np.arange(-9, 9.01, 1))
cbar.set_yticklabels([-21, -18, -15, -12, -9, -6, -3, -2, -1, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21])
cbar.set_ylabel('Atlantic Meridional Overturning Circulation (Sv)')

ax1.plot([-100, 100], [3357, 3357], '--k', linewidth = 2.0)
ax1.text(68, 3357 + 25, '$F_{\mathrm{ovS}}^{+} > 0$', verticalalignment='bottom', horizontalalignment='right', color = 'k', fontsize=12)
ax1.text(68, 3357 - 25, '$F_{\mathrm{ovS}}^{+} < 0$', verticalalignment='top', horizontalalignment='right', color = 'k', fontsize=12)

ax1.set_xlim(-30, 70)
ax1.set_ylabel('Model year')
ax1.set_ylim(2200, 4400)
ax1.set_yticks([2400, 2900, 3400, 3900, 4400])
ax1.set_xticklabels([])

CS_1		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'cyan', linewidth = 2, label = '1 PW')
CS_2		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'royalblue', linewidth = 2, label = '0.25 PW')

graphs		= CS_1 + CS_2
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax1.legend(graphs, legend_labels, loc = 'lower right', ncol=1, numpoints = 1, framealpha = 1.0)


ax2.plot(lat_MHT,MHT_atlantic_4101_4150, '-k', linewidth = 2, label = '4101 - 4150')

ax2.set_ylabel('MHT (PW)')
ax2.set_xlim(-30, 70)
ax2.set_ylim(-1.2, 1.2)
ax2.grid()

ax2.set_xticks(np.arange(-20, 60.1, 20))
ax2.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

legend_2      = ax2.legend(loc='lower center', ncol=1, framealpha = 1.0, numpoints = 1)

ax1.set_title('f) Atlantic overturning circulation (500 m) and MHT')

show()
