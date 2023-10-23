#Program determines the MOV index

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

files	= [directory+'Data/Polynya/Arctic_ice_fraction_month_10_year_4090.nc',
directory+'Data/Polynya/Arctic_ice_fraction_month_11_year_4090.nc',
directory+'Data/Polynya/Arctic_ice_fraction_month_12_year_4090.nc',
directory+'Data/Polynya/Arctic_ice_fraction_month_1_year_4091.nc',
directory+'Data/Polynya/Arctic_ice_fraction_month_2_year_4091.nc',
directory+'Data/Polynya/Arctic_ice_fraction_month_3_year_4091.nc',
directory+'Data/Polynya/Arctic_ice_fraction_month_4_year_4091.nc',
directory+'Data/Polynya/Arctic_ice_fraction_month_5_year_4091.nc',
directory+'Data/Polynya/Arctic_ice_fraction_month_6_year_4091.nc']


for file_i in range(len(files)):

	fh 		= netcdf.Dataset(files[file_i], 'r')

	lon		= fh.variables['lon'][:]
	lat		= fh.variables['lat'][:] 			
	fraction_1	= fh.variables['Fraction'][0]

	fh.close()

	#-----------------------------------------------------------------------------------------

	lon_1, lat_1, fraction_1_plot, lon_2, lat_2, fraction_2_plot, lon_3, lat_3, fraction_3_plot, lon_4, lat_4, fraction_4_plot	= LowCESMPlot(lon, lat, fraction_1)

	fig, ax         = plt.subplots(subplot_kw={'projection': ccrs.NearsidePerspective(-30, 60, 3500000)})

	CS      = ax.contourf(lon_1, lat_1, fraction_1_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
	CS      = ax.contourf(lon_2, lat_2, fraction_2_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
	CS      = ax.contourf(lon_3, lat_3, fraction_3_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
	CS      = ax.contourf(lon_4, lat_4, fraction_4_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())

	divider = make_axes_locatable(ax)
	ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
	fig.add_axes(ax_cb)

	cbar    = colorbar(CS, ticks = [15, 40, 60, 80, 100], cax=ax_cb)
	cbar.set_label('Sea-ice fraction ($\%$)', fontsize = 12)

	ax.gridlines(zorder=10)
	ax.add_feature(cfeature.LAND, zorder=10)
	ax.coastlines()
	ax.set_global()

	if file_i == 0: ax.set_title('a) Arctic sea-ice extent, October (4090)')
	if file_i == 1: ax.set_title('b) Arctic sea-ice extent, November (4090)')
	if file_i == 2: ax.set_title('c) Arctic sea-ice extent, December (4090)')
	if file_i == 3: ax.set_title('d) Arctic sea-ice extent, January (4091)')
	if file_i == 4: ax.set_title('e) Arctic sea-ice extent, February (4091)')
	if file_i == 5: ax.set_title('f) Arctic sea-ice extent, March (4091)')
	if file_i == 6: ax.set_title('g) Arctic sea-ice extent, April (4091)')
	if file_i == 7: ax.set_title('h) Arctic sea-ice extent, May (4091)')
	if file_i == 8: ax.set_title('i) Arctic sea-ice extent, June (4091)')


#-----------------------------------------------------------------------------------------

month_ice	= 3
#-----------------------------------------------------------------------------------------

fh 	= netcdf.Dataset(directory+'Data/Polynya/Arctic_ice_fraction_month_'+str(month_ice)+'_year_4091.nc', 'r')

lon		    = fh.variables['lon'][:]
lat		    = fh.variables['lat'][:] 
area		= fh.variables['AREA'][:] 						
fraction_1	= fh.variables['Fraction'][0]

fh.close()

fh 		= netcdf.Dataset(directory+'Data/Polynya/Arctic_ice_fraction_month_'+str(month_ice)+'_year_3901-3950.nc', 'r')
fraction_2	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

fh 		= netcdf.Dataset(directory+'Data/Polynya/Arctic_ice_fraction_month_'+str(month_ice)+'_year_4201-4250.nc', 'r')
fraction_3	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

#-----------------------------------------------------------------------------------------

#Get the region around the polynya
index_1			    = (fabs(lon[40] - 0.0)).argmin()
lon_polynya		    = ConverterField(index_1, lon)
lat_polynya		    = ConverterField(index_1, lat)
fraction_polynya	= ConverterField(index_1, fraction_1)
area_polynya		= ConverterField(index_1, area)

lon_polynya		    = lon_polynya[70:100, 275:310]
lat_polynya		    = lat_polynya[70:100, 275:310]
area_polynya		= area_polynya[70:100, 275:310]
fraction_polynya	= fraction_polynya[70:100, 275:310]

area_polynya	    = ma.masked_where(fraction_polynya >= 15.0, area_polynya)

print('Area polynya:', np.sum(area_polynya) /  10**12, 'million km^2')
#contourf(np.arange(len(lon_polynya[0])), np.arange(len(lat_polynya)), fraction_polynya)

fig, ax = subplots()
contourf(lon_polynya, lat_polynya, area_polynya, cmap = 'Spectral_r')

#-----------------------------------------------------------------------------------------

lon_1, lat_1, fraction_1_plot, lon_2, lat_2, fraction_2_plot, lon_3, lat_3, fraction_3_plot, lon_4, lat_4, fraction_4_plot	= LowCESMPlot(lon, lat, fraction_1)

fig, ax         = plt.subplots(subplot_kw={'projection': ccrs.NearsidePerspective(-30, 60, 3500000)})

CS      = ax.contourf(lon_1, lat_1, fraction_1_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_2, lat_2, fraction_2_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_3, lat_3, fraction_3_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_4, lat_4, fraction_4_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())

lon_1, lat_1, fraction_1_plot, lon_2, lat_2, fraction_2_plot, lon_3, lat_3, fraction_3_plot, lon_4, lat_4, fraction_4_plot	= LowCESMPlot(lon, lat, fraction_2)

CS_1	= ax.contour(lon_1, lat_1, fraction_1_plot, levels = [15], colors = 'darkblue', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_2, lat_2, fraction_2_plot, levels = [15], colors = 'darkblue', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_3, lat_3, fraction_3_plot, levels = [15], colors = 'darkblue', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_4, lat_4, fraction_4_plot, levels = [15], colors = 'darkblue', linewidths = 2, transform=ccrs.PlateCarree())

lon_1, lat_1, fraction_1_plot, lon_2, lat_2, fraction_2_plot, lon_3, lat_3, fraction_3_plot, lon_4, lat_4, fraction_4_plot	= LowCESMPlot(lon, lat, fraction_3)

CS_1	= ax.contour(lon_1, lat_1, fraction_1_plot, levels = [15], colors = 'darkturquoise', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_2, lat_2, fraction_2_plot, levels = [15], colors = 'darkturquoise', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_3, lat_3, fraction_3_plot, levels = [15], colors = 'darkturquoise', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_4, lat_4, fraction_4_plot, levels = [15], colors = 'darkturquoise', linewidths = 2, transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [15, 40, 60, 80, 100], cax=ax_cb)
cbar.set_label('Sea-ice fraction ($\%$)', fontsize = 12)

ax.gridlines(zorder=10)
ax.add_feature(cfeature.LAND, zorder=10)
ax.coastlines()
ax.set_global()

graph_1	= ax.plot([50, 50], [90, 90], '-', color = 'darkblue', linewidth = 2.0, label = '3901 - 3950', transform=ccrs.PlateCarree(), zorder =0)
graph_2	= ax.plot([50, 50], [90, 90], '-', color = 'darkturquoise', linewidth = 2.0, label = '4201 - 4250', transform=ccrs.PlateCarree(), zorder =0)

graphs	      = graph_1 + graph_2

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(12)


ax.set_title('a) Arctic sea-ice extent, March (4091)')
show()






