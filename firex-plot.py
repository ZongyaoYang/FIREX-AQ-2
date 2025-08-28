#!/usr/bin/env python3

# firex-plot.py
# Python script to create a maps of NetCDF data with global coastlines.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: January 2019
# Copyright 2019 by the University Corporation for Atmospheric Research


# set up environment variables needed by the imports
import os
os.environ['MPLCONFIGDIR']="/tmp"
os.environ['NCARG']="/usr/local/ncarg"
os.environ['NCARG_ROOT']="/usr/local/ncarg"

import warnings
import matplotlib

import random
import numpy
import math
import sys
import datetime
import calendar
import utilsLite

import netCDF4
import glob
import re
import traceback

# set up for PNG image output
matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg 
import matplotlib.pyplot

import cartopy
import mapDataFirex3
import wrf
import utilsCartopy
import warnings



BASE_DIRECTORY = "/data14a/WRF-Chem/"	# for model output on modeling1
					# change for different machines
BASE_DIRECTORY_AQ_WATCH = "/glade/campaign/acom/acom-da/shawnh/AQ_WATCH/"	# output on derecho

# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Convert string to number, or 0.0 if not numeric.
# numString = string that probably can be converted
def safeFloat(numString):
   result = -1.0
   try:
      result = float(numString)
   except ValueError:
      result = 0.0

   return result



# Print out html header.
def htmlHeader():
   print("Content-type: text/html\n")
   print("<html>")
   print("<head>")
   print("<title>Plot</title>")
   print("</head>")
   print("<body>")



# Print out html footer.
def htmlFooter():
   print("</body>")
   print("</html>")



# Print a prominent header for error messages.
def errorHeader():
   print("<h2>Error:</h2><br>")



# Print out error message in html.
def errorHtml(errorMessage):
   htmlHeader()
   errorHeader()

   htmlMessage = errorMessage.replace("\n", "<br>")
   print(htmlMessage + "\n")
   htmlFooter()



def isSurfaceVariable(varName):
   if (varName == "PM25_SRF"
      or varName == "PS"
      or varName == "AODVISdn"
      or varName == "PBLH"
      or varName == "SWDOWN"
      or varName == "NOx"	# derived from surface variables
      or varName == "AOD_550"
      or varName == "PM10_SFC"
      or varName == "PM2_5_DRY_SFC"
      or varName == "vent_rate"	# derived from surface variables
      ):
      return(True)

   return(False)



# Find the nearest value in an array.
# Return nearest value and index where it was found
def findNearest(array, value):
    index = (numpy.abs(array - value)).argmin()
    return [array[index], index]



# Read NetCDF data file and retrieve some 2-D variable.
# filename = directory and base filename from which to extract
# varName = name of at least 2-D variable in that file
# plotHeight = requested distance above surface (km)
#	Use negative value for surface.
# hour = index into the time dimension
# varUnits = units requested for values
# return tuple of [varValues, lats, lons, units]
def read2Dvar(filename, varName, plotHeight=-1.0, hour=0, varUnits=None):

   progress("read2Dvar {} from {} ".format(varName, filename))

   # open the NetCDF file for reading
   wrfData = netCDF4.Dataset(filename)

   # extract the height for interpolation
   heights = None
   if (plotHeight >= 0.0):
      heights = wrf.getvar(wrfData, "z", units="km").data
      #progress("{} shape0 = {}".format("z", heights.shape))

   # extract the requested variable
   # The indexes are: time, level, lat, lon
   if (varUnits is not None):
      varValues = wrf.getvar(wrfData, varName, units=varUnits, squeeze=False, meta=True)
   else:
      varValues = wrf.getvar(wrfData, varName, squeeze=False, meta=True)

   actualUnits = varValues.attrs["units"]
   progress("actualUnits0 = {}".format(actualUnits))

   varValues = varValues.data
   #progress("{} shape1 = {}".format(varName, varValues.shape))

   # select one time if several are supplied
   if (len(varValues.shape) == 4):
      varValues = varValues[hour, :, :, :]
   elif (len(varValues.shape) == 3):
      varValues = varValues[hour, :, :]
   #progress("{} shape2 = {}".format(varName, varValues.shape))

   if (len(varValues.shape) == 2):
      # 2D variable; ignore plot height
      {}
   else:
      # 3D variable; select some level
      if (plotHeight < 0.0):
         # surface request
         varValues = (varValues[0,:,:])
      else:
         # request is for some specified altitude
         #progress("interpolating to level {}".format(plotHeight))
         varValues = wrf.interplevel(varValues, heights, plotHeight).data

   #progress("{} shape3 = {}".format(varName, varValues.shape))

   # collapse remaining level into 2-D array
   if (len(varValues.shape) == 3):
      varValues = varValues[0, :, :]

   #progress("{} shape4 = {}".format(varName, varValues.shape))

   # retrieve latitudes and longitudes
   lats = wrf.getvar(wrfData, "XLAT").data
   lons = wrf.getvar(wrfData, "XLONG").data

   wrfData.close()

   return [varValues, lats, lons, actualUnits]



# Read NetCDF data file and retrieve wind values,
# which are rotated from grid to earth coordinates.
# filename = directory and base filename from which to extract
# plotHeight = requested distance above surface (km)
#	Use negative value for surface.
# hour = index into the time dimension
# varUnits = units requested for values
# return tuple of [uWind, vWind, lats, lons]
def read2Dwind(filename, plotHeight=-1.0, hour=0, varUnits=None):

   progress("read2Dwind from {} ".format(filename))

   # open the NetCDF file for reading
   wrfData = netCDF4.Dataset(filename)

   # extract the height for interpolation
   heights = None
   windVarName = "uvmet10"
   wantSurface = True
   if (plotHeight >= 0.0):
      heights = wrf.getvar(wrfData, "z", units="km").data
      windVarName = "uvmet"
      wantSurface = False

   # extract the requested variable
   # The indexes are: uv, time, level, lat, lon
   if (varUnits is not None):
      varValues = wrf.getvar(wrfData, windVarName, units=varUnits, squeeze=False).data
   else:
      varValues = wrf.getvar(wrfData, windVvarName, squeeze=False).data
   #progress("{} shape1 = {}".format(windVarName, varValues.shape))

   # select one time if several are supplied
   if (wantSurface):
      if (len(varValues.shape) == 4):
         varValues = varValues[:, hour, :, :]
   else:
      if (len(varValues.shape) == 5):
         varValues = varValues[:, hour, :, :, :]
   #progress("{} shape2 = {}".format(windVarName, varValues.shape))

   if (not wantSurface):
      # match size of the heights array
      #if (varValues.shape[1] > heights.shape[1]):
      #   varValues = varValues[:, :heights.shape[1], :]
      #if (varValues.shape[2] > heights.shape[2]):
      #   varValues = varValues[:, :, :heights.shape[2]]

      # request is for some specified altitude
      varValues = wrf.interplevel(varValues, heights, plotHeight).data

   #progress("{} shape3 = {}".format(windVarName, varValues.shape))

   # collapse remaining level into 3-D array
   if (len(varValues.shape) == 4):
      varValues = varValues[:, 0, :, :]

   #progress("{} shape4 = {}".format(windVarName, varValues.shape))

   # extract wind components
   uWind = varValues[0, :, :]
   vWind = varValues[1, :, :]
   #progress("{} shape5 = {}".format(windVarName, uWind.shape))

   # retrieve latitudes and longitudes
   lats = wrf.getvar(wrfData, "XLAT").data
   lons = wrf.getvar(wrfData, "XLONG").data

   wrfData.close()

   return [uWind, vWind, lats, lons]



# Return the date-dependent filename prefix of data archive.
# yearMonthDay = formatted with hyphens as 2018-09-22
# source = "firex-aq" or "aq-watch"
# domain = d01 for CONUS, d02 for Colorado region
def getNamePrefix(yearMonthDay, source="firex-aq", domain="d01"):
   # WRF-Chem
   if (source == "aq-watch" and domain== "d03"):
      # As of February 22, 2022, aq-watch domain d03 is tagged as d02.
      return("wrfmet_hourly_{}_".format("d02"))

   return("wrfout_hourly_{}_".format(domain))



# Convert units after scaling.
# originalUnits = extracted from NetCDF file
# scaledUnits = returned from scaleData()
# return scaled version of original
def convertUnits(originalUnits, scaledUnits):
   if (originalUnits == "kg/m3"):
      if (scaledUnits == "ppb"):
         return(u'\u03BC' + "g/m3")	# Greek symbol for micron

   if (originalUnits == "K"
      or originalUnits == "Pa"):
      return(originalUnits)

   if (scaledUnits is None):
      return(originalUnits)

   return(scaledUnits)



# Look up and return contour levels to use.
# varName = determines the levels used
def getContourLevelsIncremental(varName):
   # set up the default
   contourLevels = [-0.001, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15]
   if (mapDataFirex.useLogScale(varName)):
      contourLevels = [0.0009, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15]

   contourLevels = [x / 2 for x in contourLevels]
   numLevels = len(contourLevels)

   if (varName == "NOx"):
      # the typical range is 1e-5 to 1e-1
      contourLevels = [1.9e-5,
         2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3,
         1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0, 1.1]
      contourLevels = [0.9e-5,
         1e-5, 2e-5, 5e-5, 1e-4, 5e-4, 1e-3, 2e-3, 4e-3, 6e-3,
         1e-2, 2e-2, 3e-2, 7e-2, 1e-1, 2e-1, 2.1e-1] 

   if (varName == "tr17_5" or varName == "tr17_6"):
      contourLevels = [x * 2 for x in contourLevels]

   if (varName ==  "tr17_7" or varName == "tr17_8"):
      contourLevels  = contourLevels	# placeholder code from NCL

   if (varName == "tr17_1_age" or varName == "tr17_5_age"):
      contourLevels = [0.49, 0.5, 1, 2, 3, 5, 6, 8, 10, 12, 15, 20, 24, 36, 48, 72, 73]

   if (varName == "PBLH"):
      contourLevels  = [49.0, 50., 100.,200., 350, 500. ,750., 1000., 1500., 2000., 2500., 3000., 4000., 5000.,6000.,7000.,7001]

   if (varName == "SWDOWN"):
      contourLevels  = [19.0, 20, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1001]

   if (varName == "REFL_10cm"):
      contourLevels  = [-41,-40., -30., -20., -10., -5., 0., 5., 10., 20., 30., 40., 50., 60., 70., 80.,81]

   if (varName == "PM2_5_DRY_SFC" or varName == "PM10_SFC"):
      contourLevels  = [0.9, 1, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 101]

   if (varName == "AOD_550"):
      contourLevels  = [0.009, 0.01, .05, .10, .15, .20, .25, .3, .35, .4, .5, .6, .7, .8, .9, 1.0, 1.1]

   if (varName == "co"):
      contourLevels  = [0.0049, 0.005, .01, .02, .05, .1, .12, .14, .17, .2, .25, .3, .35, .4, .45, 0.5, 0.51]

   if (varName == "co_bdry"):
      contourLevels  = [0.009, 0.01, .02, .03, .04, .05, .06, .07, .08, .09, .10, .12, .14, .16, .18, 0.2, 0.21]

   if (varName == "co_chem"):
      # logarithmic
      contourLevels  = [1.9e-8, 2e-8, 4e-8, 1e-7, 2e-7, 4e-7, 1e-6, 2e-6, 4e-6,
         1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 1.1e-3]

   return(contourLevels)



# Return tick mark locations for logarithmic color scale.
# logValues = log10() has already been applied to the data
# minValue, maxValue = user-requested log bounds for color scale, or None
def logTicks(logValues, minValue, maxValue):
   # determine numerical range for the color scale
   minimum = numpy.nanmin(logValues)
   if (minValue != None):
      minimum = minValue
   maximum = numpy.nanmax(logValues)
   if (maxValue != None):
      maximum = maxValue
   span = math.ceil(maximum) - math.floor(minimum)

   if (span > 20):
      # triple powers of 10
      return numpy.arange(math.floor(minimum), math.ceil(maximum) + 2, 3)
   if (span > 5):
      # powers of 10
      return numpy.arange(math.floor(minimum), math.ceil(maximum))

   # space ticks at multiples of 1, 2, and 5
   ticks = numpy.zeros(int(span * 3))
   exponent = math.floor(minimum)
   for ti in range(0, len(ticks), 3):
      ticks[ti] = exponent
      ticks[ti+1] = exponent + 0.30103
      ticks[ti+2] = exponent + 0.69897
      exponent += 1

   return ticks



# Format a linear tick label.
# tickValue = to show on axis
# pos = tick position; index from 0 through greatest tick value
def formatLinear(tickValue, pos=None):
   return("{:.3g}".format(tickValue))

# Format an exponential tick label.
# expValue = 10 to this power
# pos = tick position; index from 0 through greatest tick value
def formatExponent(expValue, pos=None):
   return("{:.3g}".format(math.pow(10, expValue)))



degreeSymbol = chr(176)
microSymbol = u'\u03BC'

# Create a filename for saving the image.
def makeSaveName(plotTitle):
   temp = plotTitle.replace(' at', '')
   temp = temp.replace(' time', '')
   temp = temp.replace('Time: ', '')
   temp = temp.replace(' hPa', '')

   temp = temp.replace(' (ppm)', '')
   temp = temp.replace(' (ppb)', '')
   temp = temp.replace(' (ppt)', '')

   temp = temp.replace(degreeSymbol, '')
   temp = temp.replace(microSymbol, 'u')
   temp = temp.replace(',', '_')

   return temp.replace(' ', '-')



# Plot data somewhere on a coastline map of the world.
# var2D = some lat-lon variable to plot
# lats = latitude values of var2D
# lons = longitude values of var2D
# latBounds = lower and upper latitude (-90:90)
# lonBounds = west and east longitude (-180:180)
# altitude = height in km above surface, or -1.0 for surface
# plotTitle = appears above the plot
# varName = friendly name of species to plot
# subDir = sub-directory within images/, probably year-month-day
# saveToFilename = name of image file to create in data directory
# uWind, vWind = west-east and south-north components of wind
# initDateStr = forecast initialized in this format: 2019-07-19
# source = "firex-aq" or "aq-watch"
# varUnits = units for this variable
def plotOnWorld(var2D, lats, lons, latBounds, lonBounds, altitude,
   plotTitle, varName, subDir, saveToFilename=None, uWind=None, vWind=None,
   initDateStr=None, source="firex-aq", varUnits="m"):

   #progress("lons = {}".format(lons))
   #progress("lats = {}".format(lats))

   progress("plotOnWorld save to {}".format(saveToFilename))
   numLats = lats.shape[0]
   numLons = lats.shape[1]

   # get the actual lat-lon bounds of the plot
   lat1 = latBounds[0]
   lat2 = latBounds[1]
   lon1 = lonBounds[0]
   lon2 = lonBounds[1]
   #progress("Region lat bounds = {} to {} lon bounds = {} to {}".format(lat1, lat2, lon1, lon2))

   if (lon1 > lon2):
      # map crosses International Date Line; shift left
      # Need the cartopy equivalent here. Carl Drews - December 15, 2020.
      var2D, lons = TKBM.shiftgrid(lon1, var2D, lons, start=False)
      lon1 -= 360

   # create Cartopy instance for projection.
   matplotlib.pyplot.figure(dpi=150)
   map = cartopy.crs.LambertConformal(
      central_longitude=-97.0, central_latitude=39.0,
      standard_parallels=(30.0, 60.0))
   ax = matplotlib.pyplot.axes(projection=map)
   ax.set_extent([lon1, lon2, lat1, lat2])

   #progress("plotOnWorld(): range of var2D = {} to {}"
   #   .format(numpy.nanmin(var2D), numpy.nanmax(var2D)))

   # Make filled contour plot.
   # begin with the pre-set contour levels
   contourLevels = mapDataFirex3.getContourLevels(varName, units=varUnits)
   numContours = len(contourLevels)
   progress("contourLevels1 = {}".format(contourLevels))

   # determine if data exceeds color scale
   extendAdded = "neither"
   if (contourLevels[numContours-1] < numpy.nanmax(var2D)):
      extendAdded = "max"
   if (numpy.nanmin(var2D) < contourLevels[0]):
      if (extendAdded == "neither"):
         extendAdded = "min"
      else:
         extendAdded = "both"

   with warnings.catch_warnings():
      # warnings mess up our script output to the browser
      warnings.filterwarnings("ignore", category=DeprecationWarning)

      progress("plotOnWorld(): range of var2D = {} to {}"
         .format(numpy.nanmin(var2D), numpy.nanmax(var2D)))
      contourMap = ax.contourf(lons, lats, var2D, contourLevels,
         cmap=mapDataFirex3.getColorMap(varName, numContours-1),
         transform=cartopy.crs.PlateCarree(),
         alpha=0.8, extend=extendAdded)

   continentalScale = (lat2 - lat1 >= 15.0)
   stateScale = (lat2 - lat1 >= 4.4)
   colorbarShrink = 0.8
   if ((not continentalScale) and (not stateScale)):
      colorbarShrink = 0.6

   # set up line width
   width = 1.0
   stateWidth = width/2
   countyWidth = stateWidth/2

   if ((not continentalScale)):
      stateWidth *= 2
      countyWidth *= 2
      if ((not stateScale)):
         countyWidth *= 2

   ax.coastlines(resolution='50m', linewidth=width)
   ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=width)
   ax.add_feature(cartopy.feature.STATES.with_scale("50m"), linewidth=stateWidth)

   # Mexico has states (provinces) too!
   ax.add_geometries(utilsCartopy.mexicanStates(), cartopy.crs.PlateCarree(),
      edgecolor="black", facecolor='none', linewidth=stateWidth)

   if (not continentalScale):
      ax.add_geometries(utilsCartopy.usCounties(), cartopy.crs.PlateCarree(),
         edgecolor="black", facecolor='none', linewidth=countyWidth)
      progress("stateWidth = {}   countyWidth = {}".format(stateWidth, countyWidth))

   xTickLoc = utilsCartopy.getXTickLocator(lonBounds)
   yTickLoc = utilsCartopy.getXTickLocator(latBounds)

   # draw geographic lines on map
   geoLines = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False,
      xlocs=xTickLoc, ylocs=yTickLoc,
      xformatter=matplotlib.ticker.FuncFormatter(utilsCartopy.formatLongitudeConus),
      yformatter=matplotlib.ticker.FuncFormatter(utilsCartopy.formatLatitudeConus),
      linewidth=width/2, linestyle=':', color="black")
   geoLines.right_labels = False
   geoLines.top_labels = False
   geoLines.rotate_labels = False

   # adjust aspect ratio for this latitude
   ax.set_aspect(1.0)

   progress("wind barbs...", endString="\n")
   # draw wind barbs
   if (uWind is not None):
      if (continentalScale):
         stride = 10
      elif (stateScale):
         stride = 2
      else:
         stride = 1

      # adjust barb spacing for high-resolution grids (< 12 km)
      halfLatIndex = int(lats.shape[0] / 2)
      halfLonIndex = int(lats.shape[1] / 2)
      gridResolution = (lats[halfLatIndex + 1, halfLonIndex]
         - lats[halfLatIndex, halfLonIndex]) * 111.1	# km
      progress("gridResolution = {} km".format(gridResolution))
      stride *= int(round(12.3 / gridResolution))
      progress("stride = {}".format(stride))

      barbUWind = uWind[0::stride, 0::stride]
      barbVWind = vWind[0::stride, 0::stride]
      numRows = min(barbUWind.shape[0], barbVWind.shape[0])
      numColumns = min(barbUWind.shape[1], barbVWind.shape[1])
      #progress("Rows: {} Columns: {} for wind barbs.".format(numRows, numColumns))

      barbX = lons[0::stride, 0::stride]
      barbY = lats[0::stride, 0::stride]
      #progress("Shapes: {} {} {} {}".
      #   format(str(barbX.shape), str(barbY.shape), str(uWind.shape), str(vWind.shape)))

      #progress("barbUWind = " + str(barbUWind))
      #progress("barbVWind = " + str(barbVWind))

      # At 3 km height barbs() produces this runtime warning:
      # /usr/local/lib64/python3.6/site-packages/cartopy/mpl/geoaxes.py:1913: RuntimeWarning: invalid value encountered in less / greater
      # There are NaNs in the wind fields because of high mountains in Colorado and Wyoming.
      # Carl Drews - November 16, 2020
      with warnings.catch_warnings():
         warnings.filterwarnings("ignore", category=RuntimeWarning)
         ax.barbs(barbX, barbY,
            barbUWind[:numRows, :numColumns], barbVWind[:numRows, :numColumns],
            transform=cartopy.crs.PlateCarree(),
            length=4, linewidth=0.2)

   # place title and color scale
   if (initDateStr is None):
      matplotlib.pyplot.title(plotTitle)
   else:
      matplotlib.pyplot.suptitle(plotTitle, fontsize=14, y=0.98)
      matplotlib.pyplot.title(
         "Forecast initialized at: {} 00:00 UTC".format(initDateStr),
         fontsize=10)

   formatter = formatLinear
   if mapDataFirex3.useLogScale(varName):
      formatter = formatExponent

   barTicks = contourLevels[0:numContours:2]
   if (varName == "vent_rate"):
      barTicks = contourLevels[0:numContours]

   cbar = matplotlib.pyplot.colorbar(contourMap,
      shrink=colorbarShrink, ticks=barTicks,
      format=matplotlib.ticker.FuncFormatter(formatter),
      orientation='horizontal', pad=0.08, drawedges=True)
   cbar.ax.tick_params(labelsize=8)

   # stamp the image with the chemical species
   #matplotlib.pyplot.figtext(0.149, 0.568, varName, ha="left", fontsize=11)

   # make an image file on disk
   imageDir = "/webt/firex-aq/"
   imageDir = "/home/drews/FIREX-AQ/webt/firex-aq/"
   #imageDir = "/home/drews/FIREX-AQ/images/"		# bogus - for testing
   if (source == "aq-watch"):
      imageDir = "/glade/derecho/scratch/drews/FIREX-AQ/webt/aq-watch/"		# derecho
   imageDir += subDir

   # ensure that year-month-day sub-directory exists
   if (not os.path.isdir(imageDir)):
      progress("Create year-month-day directory {}".format(imageDir))
      os.mkdir(imageDir)

   fullSavePath = imageDir + saveToFilename
   progress("Saving image to file {} ...".format(fullSavePath))

   try:
      # trim whitespace but leave a narrow margin
      matplotlib.pyplot.savefig(fullSavePath,
         bbox_inches="tight", pad_inches=0.04)

   except RuntimeError as oops:
      progress("{}".format(traceback.format_exc()))

   matplotlib.pyplot.clf()
   matplotlib.pyplot.close()

   return(0)	# no error code



# Generate and return suitable tick marks for vertical log scale
# minHeight, maxHeight = vertical range of plot in hPa (1000.0 : 1.0)
def getYTicks(minHeight, maxHeight):
   allTicks = [
         1.0e-5, 2.0e-5, 5.0e-5,
         1.0e-4, 2.0e-4, 5.0e-4,
         1.0e-3, 2.0e-3, 5.0e-3,
         1.0e-2, 2.0e-2, 5.0e-2,
         1.0e-1, 2.0e-1, 5.0e-1,
         1.0, 2.0, 5.0,
         10.0, 20.0, 50.0,
         100.0, 200.0, 500.0, 1000.0]

   # carve off tick marks at both ends if not within range
   while (allTicks[0] < maxHeight):
      del allTicks[0]
   while (allTicks[len(allTicks)-1] > minHeight):
      del allTicks[len(allTicks)-1]

   return allTicks



# Return a more readable name for a chemical species.
def friendlyName(netcdfName):
   betterName = netcdfName
   suffix1 = "_VMR_inst"
   suffix2 = "_VMR_avrg"

   # remove the suffix, if any
   if (netcdfName.endswith(suffix1)):
      betterName = netcdfName[:-len(suffix1)]
   elif (netcdfName.endswith(suffix2)):
      betterName = netcdfName[:-len(suffix2)]

   return betterName



# Return lat-lon bounds for the requested region.
# return tuple of [gotRegion, latBounds[], lonBounds[]]
def getRegion(regionName):
   latBounds = None
   lonBounds = None

   if (regionName == "Global"):
      latBounds = [-90, 90]	# world
      lonBounds = [-180, 180]

   if (regionName == "ConUS".lower()):
      latBounds = [23.50, 50.00]	# continental United States
      lonBounds = [-119.90, -74.00]

   if (regionName == "Colorado".lower()):
      latBounds = [36.80, 41.20]	# state of Colorado
      lonBounds = [-110.20, -100.92]

   if (regionName == "Front-Range".lower()):
      latBounds = [39.33, 40.82]	# Colorado Front Range
      lonBounds = [-106.03, -103.66]

   if (regionName == "TropicalAmerica"):
      latBounds = [0, 32]	# central America
      lonBounds = [-120, -70]

   if (regionName == "SouthAmerica"):
      latBounds = [-60, 20]	# South America
      lonBounds = [-90, -30]

   if (regionName == "Asia"):
      latBounds = [0, 70]	# Asia
      lonBounds = [30, 150]

   if (regionName == "Australia"):
      latBounds = [-45, 5]	# Australia
      lonBounds = [105, 180]

   if (regionName == "Europe"):
      latBounds = [30, 70]	# Europe
      lonBounds = [-15, 60]

   if (regionName == "Africa"):
      latBounds = [-40, 40]	# Africa
      lonBounds = [-20, 55]

   if (regionName == "Pacific"):
      latBounds = [-60, 60]	# Pacific Ocean
      lonBounds = [120, -70]

   if (regionName == "NorthAmerica"):
      latBounds = [5, 75]	# North America
      lonBounds = [-170, -50]

   return [latBounds != None, latBounds, lonBounds]



# Is this variable derived by calculation from other variables?
def isDerived(myVarName):
   if (myVarName == "NOx"
      or myVarName == "tr17_1_age" or myVarName == "tr17_5_age"
      or myVarName == "vent_rate" or myVarName == "RH"
      or myVarName == "PBLH"):          # converted to feet
      return(True)

   return(False)



# Plot all WRF output found in given directory for today.
# today = UTC date-time at which script began running
# wrfDir = directory containing WRF-Chem output
# regions = list of regions to plot
# chemical = list of species to plot
# heights = list of vertical levels to plot
# source = "firex-aq" or "aq-watch"
# domain = "d01" for CONUS, "d02" for AQ_WATCH Colorado and Front Range
def plotToday(today, wrfDir, regions, chemicals, heights,
   source="firex-aq", domain="d01"):
   retValue = 0		# no error

   # look in the data sub-directory for today
   year = today.year
   month = today.month
   day = today.day

   wrfSubDir = "wrf/"
   namePrefix = "wrfout"
   if (source == "aq-watch" and domain == "d03"):
      wrfSubDir = "wrf_met_frange/"
      namePrefix = "wrfmet"

   yearMonthDay = "{:04d}{:02d}{:02d}".format(year, month, day)
   todayDir = wrfDir + yearMonthDay + "/" + wrfSubDir
   progress("todayDir = {}".format(todayDir))
   todayStr = "{:04d}-{:02d}-{:02d}".format(year, month, day)

   # pattern-match the file names
   filePattern = getNamePrefix(yearMonthDay, source, domain) + "*:00:00"
   progress("filePattern = {}".format(filePattern))
   fileNames = glob.glob(todayDir + filePattern)
   fileNames.sort()

   # loop through the hourly files we found
   plotResult = 0
   for fileName in fileNames:
      # the full filenames look like this:
      # /WRF-Chem/20190418/wrf/wrfout_hourly_d01_2019-04-19_16:00:00
      progress("wrfOutput = {}".format(fileName))

      # extract the date and time
      justName = os.path.basename(fileName)
      fileStamp = (datetime.datetime.strptime(justName,
         getNamePrefix(yearMonthDay, source, domain) + "%Y-%m-%d_%H:00:00"))
      forecastYear = fileStamp.year
      forecastMonth = fileStamp.month
      forecastDay = fileStamp.day
      forecastHour = fileStamp.hour

      # loop through the requested heights
      for height in heights:
         # read the wind fields for this time and height
         yearMonthDayHour = ("{:04}-{:02}-{:02}_{:02}".
            format(forecastYear, forecastMonth, forecastDay, forecastHour))
         windFilename = getNamePrefix(yearMonthDayHour, source, domain) + yearMonthDayHour + ":00:00"

         todayWind = todayDir + windFilename
         if (not os.path.exists(todayWind)):
            progress("WRF file {} does not exist.".format(todayWind))
            continue

         try:
            # these are hourly files containing only one time step
            windData = read2Dwind(todayWind, height, 0, "kt")

         except FileNotFoundError as oops:
            progress(oops)
            continue

         for species in chemicals:
            speciesDir = todayDir
            if (isDerived(species)):
               # derived files have to go in a writable directory on derecho
               speciesDir = speciesDir.replace("shawnh", "drews")

            if (height >= 0.0):
               if (isSurfaceVariable(species)):
                  # skip the higher heights for surface-only variables
                  continue

            nowStr = ("{:04}-{:02}-{:02}"
               .format(forecastYear, forecastMonth, forecastDay))
            if (todayStr == nowStr and forecastHour == 0):
               if (species == "PBLH" or species == "vent_rate"):
                  progress("Skipping initialization hour for {}.".format(species))
                  continue

            plotResult = createPlot(speciesDir,
               forecastYear, forecastMonth, forecastDay, forecastHour,
               regions, species, height, windData, todayStr,
               source, domain)

   return(plotResult)



# Make regional plots of WRF-Chem data.
# dataDir = where to look for WRF-Chem model output
# year, month, day, hour = time of model output
# regions = list of named regions on globe or continent
# species = name of the variable to plot
# height = km above surface, or negative for at surface
# windData = already loaded by read2Dwind() for this time and height
# initDateStr = date forecast initialized as a string: 2019-07-19
# source = "firex-aq" or "aq-watch"
# domain = "d01" or "d02"
def createPlot(dataDir, year, month, day, hour,
   regions, species, height, windData,
   initDateStr=None, source="firex-aq", domain="d01"):

   dataDirDerived = dataDir

   ### Read data values from NetCDF data file. ###
   yearMonthDayHour = "{:04}-{:02}-{:02}_{:02}".format(year, month, day, hour)
   filename = getNamePrefix(yearMonthDayHour, source, domain) + yearMonthDayHour + ":00:00"
   if (isDerived(species)):
      filename += "-derived.nc"

   if (source == "aq-watch" and domain == "d03"):
      # different naming for derived files
      dataDirDerived = dataDirDerived.replace("wrf_met_frange/", "wrf/")
      filename = filename.replace("wrfmet_hourly", "wrfout_hourly")
      filename = filename.replace("d02", "d03")

   if (not os.path.exists(dataDirDerived + filename)):
      progress("WRF file {} does not exist.".format(dataDirDerived + filename))
      return(1)		# error

   try:
      # these are hourly files containing only one time step
      waccmData = read2Dvar(dataDirDerived + filename, species, height, 0)

   except (FileNotFoundError, IOError) as oops:
      htmlHeader()
      errorHeader()
      print("Could not access file {}<br>".format(filename))
      print("{}<br>".format(oops))
      print("Please check the date range.")
      htmlFooter()
      sys.exit(1)

   # collect parts of the tuple that was returned
   values = waccmData[0]
   lats = waccmData[1]
   lons = waccmData[2]
   units = waccmData[3]
   progress("units0 = {}".format(units))

   # fill in missing units
   if ((species == "c2h6" or species == "isopr")
      and units == ""):
      units = "ppmv"
      progress("Setting {} units to {}".format(species, units))

   uWind = windData[0]
   vWind = windData[1]

   #progress("lats = " + str(lats))
   #progress("lons = " + str(lons))

   # create hour portions of the plot title
   hoursTitle = "{:02}:{:02}".format(hour, 0)

   heightTitle = "surface"
   if (height >= 0.0):
      heightTitle = "{0:.1f} km".format(height)

   varName = friendlyName(species)

   # scale the data for plotting
   progress("createPlot(): range of values = {} to {}"
      .format(numpy.nanmin(values), numpy.nanmax(values)))
   scaleUnits = mapDataFirex3.scaleData(varName, values, False)
   progress("scaleUnits1 = {}".format(scaleUnits))

   ### Create the Map plot type. ###
   scaleUnits = convertUnits(units, scaleUnits)
   progress("scaleUnits2 = {}".format(scaleUnits))

   if (scaleUnits == "Pa"):
      # convert Pascals to hectoPascals
      scaleUnits = "hPa"
      values /= 100

   if (scaleUnits == "ppmv"):
      # convert parts per million to per billion
      scaleUnits = "ppbv"
      if (mapDataFirex3.useLogScale(varName)):
         values += 3
      else:
         values *= 1000

   if (varName == "PBLH" and units == "m"):
      # show boundary layer height in km
      scaleUnits = "km"
      values /= 1000

   if (varName == "PBLH" and units == "ft"):
      # show boundary layer height in k-ft
      scaleUnits = "kft"
      values /= 1000

   if (varName == "vent_rate" and units == "kt ft"):
      # convert to a manageable scale
      scaleUnits = "1e3 kt ft"
      values /= 1000

   # create plot title
   yearMonthDay = "{:04}-{:02}-{:02}".format(year, month, day)
   plotTitle = "{} at {}".format(varName, heightTitle)
   if (source == "aq-watch"):
      plotTitle = "{} {} at {}".format(varName, domain, heightTitle)
   plotTitle += " {} {}".format(yearMonthDay, hoursTitle)
   if (scaleUnits != None):
      plotTitle += " (" + scaleUnits + ")"

   # set default map boundaries
   latBounds = [18.0, 55.0]
   lonBounds = [-125.0, -50.0]

   for region in regions:
      if (region.lower() == "front-range" and height > 0.0):
         # skip Front Range of Colorado for above surface
         progress("Skipping region {} for height {}.".format(region, height))
         response = 0	# no error code
         continue

      # get the preset region, if any
      temp = getRegion(region)
      if (temp[0]):
         latBounds = temp[1]
         lonBounds = temp[2]

      yearMonthDayDir = ("{:04d}{:02d}{:02d}/"
         .format(year, month, day))
      plotFilename = ("{:04d}-{:02d}-{:02d}_{:02d}00_{}_{}_{}"
         .format(year, month, day, hour, species,
         heightString(height), region))
      if (source == "aq-watch"):
         # include domain for AQ_WATCH
         plotFilename += "_" + domain
      plotFilename += ".png"
      response = plotOnWorld(values, lats, lons, latBounds, lonBounds, height,
         plotTitle, varName, yearMonthDayDir, plotFilename, uWind, vWind,
         initDateStr, source, varUnits=scaleUnits)
      progress("")

   return(response)



# Convert height as float kilometers to string.
def heightString(floatHeight):
   if (floatHeight >= 0.0):
      return("{:.0f}km".format(floatHeight))

   else:
      return("Sfc")



# Main routine begins here.
def main():
   ### Set up the default parameters. ###
   regions = ["conus", "colorado", "front-range"]
   source = "firex-aq"
   domains = ["d01"]
   variables = ["NOx",
      "PBLH", "SWDOWN",
      # Inert tracers are not calculated beginning on 22 April 2020.
      #"tr17_1_age", "tr17_2", "tr17_4", "tr17_5_age", "tr17_6", "tr17_8",
      "o3", "PM2_5_DRY_SFC", "PM10_SFC", "AOD_550",
      "co", "co_anth", "co_fire", "co_bdry", "co_bdry_fire", "co_chem",
      "co_asia", "isopr",
      "vent_rate", "RH"
   ]
   heights = [-1.0, 3.0, 5.0, 8.0]		# km above surface, or -1 for surface

   today = False
   year = 2020
   month = 8
   days = None
   hours = None

   # The forecast is generated on a certain date for several days ahead.
   # specify the forecast date like this: 20190428
   # If no forecast specified, assume the same date as requested.
   forecast = None

   # get command-line arguments, if any
   for argPair in sys.argv:
      progress("argPair = " + argPair)
      # the arguments are: arg=value
      pairValue = argPair.split('=')
      if (len(pairValue) < 2):
         continue

      if (pairValue[0] == "forecast"):
         forecast = re.sub("[^0-9]", "", pairValue[1])

      if (pairValue[0] == "year"):
         if (utilsLite.sanitize(pairValue[1]).lower() == "today"):
            today = True
         else:
            year = int(pairValue[1])

      if (pairValue[0] == "month"):
         month = int(pairValue[1])

      # Days and hours can be specified singly;
      # otherwise user wants all days and hours.
      if (pairValue[0] == "day"):
         days = [int(pairValue[1])]
      if (pairValue[0] == "hour"):
         hours = [int(pairValue[1])]

      if (pairValue[0] == "region"):
         regions = pairValue[1].split(',')
         for ri in range(len(regions)):
            regions[ri] = utilsLite.sanitize(regions[ri]).lower().strip()
      if (pairValue[0] == "source"):
         source = utilsLite.sanitize(pairValue[1])
      if (pairValue[0] == "domain"):
         domains = pairValue[1].split(',')
         for di in range(len(domains)):
            domains[di] = utilsLite.sanitize(domains[di]).lower().strip()
      if (pairValue[0] == "species"):
         variables = pairValue[1].split(',')
         for si in range(len(variables)):
            variables[si] = utilsLite.sanitize(variables[si]).strip()
      if (pairValue[0] == "height"):
         heights = [safeFloat(pairValue[1])]

   baseDirectory = BASE_DIRECTORY
   if (source == "aq-watch"):
      baseDirectory = BASE_DIRECTORY_AQ_WATCH

   if (today):
      # Record date-time now in case domain plotting extends into UTC tomorrow.
      dayDelta = 0		# use this to offset today (yesterday, tomorrow, last week...)
      nowScript = datetime.datetime.utcnow()
      nowScript += datetime.timedelta(days=dayDelta)

      # process all files in today's directory
      for domain in domains:
         retValue = plotToday(nowScript, baseDirectory, regions, variables,
            heights, source, domain)
      return(retValue)

   if (days is None):
      # expand into all days that month
      lastDay = calendar.monthrange(year, month)[1]
      days = [di for di in range(1, lastDay+1)]

   # expand into all hours of the day
   if (hours is None):
      # Consider running this loop backwards to plot recent data first. Carl Drews - April 12, 2019
      hours = [hi for hi in range(0, 24)]

   progress("year = {} month = {}".format(year, month))
   progress("days = {}".format(days))
   progress("hours = {}".format(hours))
   progress("heights = {}".format(heights))
   progress("species = {}".format(variables))
   progress("regions = {}".format(regions))
   progress("source = {}".format(source))
   progress("domains = {}".format(domains))

   plotResult = 0	# no error code
   # loop through the lists of days and hours
   for day in days:
      for hour in hours:
         forecastStr = None

         ### Read data values from NetCDF data file. ###
         if (forecast is not None):
            dataDir = "{}{}/wrf/".format(baseDirectory, forecast)
            forecastStr = "{}-{}-{}".format(
               forecast[0:4], forecast[4:6], forecast[6:8])
         else:
            dataDir = "{}{:04}{:02}{:02}/wrf/".format(baseDirectory, year, month, day)
            forecastStr = "{:04}-{:02}-{:02}".format(year, month, day)
         progress("dataDir = {}".format(dataDir))

         for height in heights:
            for domain in domains:
               progress("\ndomain = {}".format(domain))

               # read the wind fields for this time and height
               yearMonthDayHour = "{:04}-{:02}-{:02}_{:02}".format(year, month, day, hour)
               windFilename = getNamePrefix(yearMonthDayHour, source, domain) + yearMonthDayHour + ":00:00"

               # As of February 22, 2022, the WRF-Chem d03 output files are tagged "d02".
               if (domain == "d03"):
                  windFilename = windFilename.replace("d03", "d02")

               dataWindFile = dataDir + windFilename
               if (source == "aq-watch" and domain == "d03"):
                  # different naming for this model output
                  dataWindFile = dataWindFile.replace("wrf/", "wrf_met_frange/")

               if (not os.path.exists(dataWindFile)):
                  progress("WRF file {} does not exist.".format(dataWindFile))
                  continue
   
               try:
                  # these are hourly files containing only one time step
                  windData = read2Dwind(dataWindFile, height, 0, "kt")
   
               except FileNotFoundError as oops:
                  progress(oops)
                  continue
   
               for species in variables:
                  speciesDir = dataDir
                  if (isDerived(species)):
                     # derived files have to go in a writable directory on derecho
                     speciesDir = speciesDir.replace("cheyenne", "derecho")
                     speciesDir = speciesDir.replace("shawnh", "drews")

                  if (height >= 0.0):
                     if (isSurfaceVariable(species)):
                        # skip the higher heights for surface-only variables
                        continue
   
                  yearMonthDay = "{:04}-{:02}-{:02}".format(year, month, day)
                  progress("forecastStr = {}   yearMonthDay = {}".format(forecastStr, yearMonthDay))
                  if (forecastStr == yearMonthDay and hour == 0):     # disable for backfill
                     if (species == "PBLH" or species == "vent_rate"):
                        progress("Skipping initialization hour for {}.".format(species))
                        continue
   
                  plotResult = createPlot(speciesDir, year, month, day, hour,
                     regions, species, height, windData, forecastStr,
                     source, domain)

   return(plotResult)



# call the main
progress("{}".format(__file__))
progress("Start time: {}".format(datetime.datetime.now()))
retValue = main()
progress("End time: {}".format(datetime.datetime.now()))

sys.exit(retValue)

