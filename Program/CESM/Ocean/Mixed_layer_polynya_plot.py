#Program plots the mixed layer depth during polynya formation

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

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'


def ConverterField(index_break, field):
	"""Shifts field, where it starts at 0E and ends at 360E"""

	new_field	= ma.masked_all(shape(field))
	length_section	= len(field[0]) - index_break

	#Shift the first part
	new_field[:, :length_section] = field[:, index_break:]

	#Shift the last part
	new_field[:, length_section:] = field[:, :index_break] 

	return new_field

def PeriodicBoundaries3D(lon, lat, field, lon_grids = 1):
        """Add periodic zonal boundaries for 2D field"""

        #Empty field with additional zonal boundaries
        lon_2                   = np.zeros((len(lat), len(lon[0]) + lon_grids * 2))
        lat_2                   = np.zeros((len(lat), len(lon_2[0])))
        field_2                 = ma.masked_all((len(field), len(lat), len(lon_2[0])))

        #Get the left boundary, which is the right boundary of the original field
        lon_2[:, :lon_grids]       	= lon[:, -lon_grids:]
        lat_2[:, :lon_grids]       	= lat[:, -lon_grids:]
        field_2[:, :, :lon_grids]	= field[:, :, -lon_grids:]

        #Same for the right boundary
        lon_2[:, -lon_grids:]        	= lon[:, :lon_grids]
        lat_2[:, -lon_grids:]        	= lat[:, :lon_grids]
        field_2[:, :, -lon_grids:]	= field[:, :, :lon_grids]

        #And the complete field
        lon_2[:, lon_grids:-lon_grids]          = lon
        lat_2[:, lon_grids:-lon_grids]          = lat
        field_2[:, :, lon_grids:-lon_grids]     = field

        return lon_2, lat_2, field_2

def LowCESMPlot(lon, lat, field):
	"""Returns 4 array's to plot on a global projection"""

	#Left of pole
	lon[lon > 180]	= lon[lon > 180] - 360.0

	lon_1		= lon[:, :160]
	lat_1		= lat[:, :160]
	field_1		= field[:, :160]

	#Right of pole
	lon_2		= lon[:, 159:]
	lat_2		= lat[:, 159:]
	field_2		= field[:, 159:]

	lat_3		= ma.masked_where(lon_2 > 0.0, lat_2)
	field_3		= ma.masked_where(lon_2 > 0.0, field_2)
	lon_3		= ma.masked_where(lon_2 > 0.0, lon_2)

	lon_2[lon_2 < -160] 	= lon_2[lon_2 < -160] + 360
	lat_2			= ma.masked_where(lon_2 < 0.0, lat_2)
	field_2			= ma.masked_where(lon_2 < 0.0, field_2)
	lon_2			= ma.masked_where(lon_2 < 0.0, lon_2)

	#To match at 40W
	index_1		= (fabs(lon[40] - 0.0)).argmin()

	lon_4		= ConverterField(index_1, lon)
	lat_4		= ConverterField(index_1, lat)
	field_4		= ConverterField(index_1, field)

	lon_4		= lon_4[:, 280:300]
	lat_4		= lat_4[:, 280:300]
	field_4		= field_4[:, 280:300]

	return lon_1, lat_1, field_1, lon_2, lat_2, field_2, lon_3, lat_3, field_3, lon_4, lat_4, field_4

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

month_ice	= 3

#-----------------------------------------------------------------------------------------

fh 	= netcdf.Dataset(directory+'Data/Polynya/Mixed_layer_depth_month_'+str(month_ice)+'_year_4051-4100.nc', 'r')

lon		= fh.variables['lon'][:]
lat		= fh.variables['lat'][:] 			
mixed	= fh.variables['MXL'][:]

fh.close()

#-----------------------------------------------------------------------------------------

lon_1, lat_1, mixed_1_plot, lon_2, lat_2, mixed_2_plot, lon_3, lat_3, mixed_3_plot, lon_4, lat_4, mixed_4_plot	= LowCESMPlot(lon, lat, mixed[40])

fig, ax         = plt.subplots(subplot_kw={'projection': ccrs.NearsidePerspective(-30, 60, 3500000)})

CS      = ax.contourf(lon_1, lat_1, mixed_1_plot, levels = np.arange(0, 1600.01, 100), cmap = 'Reds', extend = 'max', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_2, lat_2, mixed_2_plot, levels = np.arange(0, 1600.01, 100), cmap = 'Reds', extend = 'max', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_3, lat_3, mixed_3_plot, levels = np.arange(0, 1600.01, 100), cmap = 'Reds', extend = 'max', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_4, lat_4, mixed_4_plot, levels = np.arange(0, 1600.01, 100), cmap = 'Reds', extend = 'max', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(0, 1601, 200), cax=ax_cb)
cbar.set_label('Mixed layer depth (m)', fontsize = 12)

#Irminger Basin
lon_irm	= [-44, -29.4, -20, -44, -44]
lat_irm = [60, 60, 66, 66, 60]

lon_irm_1 	= np.arange(-44, -29.4, 0.1)
lat_irm_1 	= np.zeros(len(lon_irm_1))+60
lon_irm_2 	= np.arange(-29.4, -19.99, 0.1)
lat_irm_2	= np.linspace(60, 66, len(lon_irm_2), endpoint = True)
lon_irm_3 	= np.arange(-44, -20, 0.1)[::-1]
lat_irm_3	= np.zeros(len(lon_irm_3))+66
lat_irm_4	= np.arange(60.1, 66, 0.1)[::-1]
lon_irm_4	= np.zeros(len(lat_irm_4))-44

lon_irm		= np.append(lon_irm_1, lon_irm_2)
lon_irm		= np.append(lon_irm, lon_irm_3)
lon_irm		= np.append(lon_irm, lon_irm_4)
lat_irm		= np.append(lat_irm_1, lat_irm_2)
lat_irm		= np.append(lat_irm, lat_irm_3)
lat_irm		= np.append(lat_irm, lat_irm_4)

ax.plot(lon_irm, lat_irm, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 15)

ax.gridlines(zorder=10)
ax.add_feature(cfeature.LAND, zorder=10)
ax.coastlines()
ax.set_global()

ax.set_title('b) Mixed layer depth, March (4091)')
show()



