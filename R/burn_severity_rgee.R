# this script adapts Java Script GEE code 
# originally authored by Kyle C. Rodman, Ecological Restoration Institute, NAU 
# translational author: Zachary Hoylman, zachary.hoylman@umontana.edu 

# set up conda enviorment
# conda create --name gee      # Create a conda environment
# conda activate gee# Activate the environment
# conda install python=3.7
# conda install -c conda-forge earthengine-api # Install the Earth Engine Python API
# earthengine authenticate          # Authenticate your access with the command line tool
# conda install pandas
# conda install numpy

library(reticulate)
library(cptcity)
library(raster)
library(stars)
library(sf)
library(tidyverse)

#set up enviorment
use_condaenv("gee", conda = "auto",required = TRUE)
ee = import("ee")
#crashes rstudio if you load this before initiating the conda enviorment??
library(rgee)
ee_Initialize(drive = TRUE)

# Pulling in Landsat data, creating composite, and identifying survivors
# Harmonization function and image composite functions from LandTrendr source code 
# (Yang, Braaten, and Kennedy; users/emaprlab/public)

harmonizationRoy = function(oli) {
  slopes = ee$Image$constant(c(0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949))           # create an image of slopes per band for L8 TO L7 regression line - David Roy
  itcp = ee$Image$constant(c(-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029))        # create an image of y-intercepts per band for L8 TO L7 regression line - David Roy
  y = oli$select(c('B2','B3','B4','B5','B6','B7'),c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))$ # select OLI bands 2-7 and rename them to match L7 band names
    subtract(itcp$multiply(10000))$divide(slopes)$                                        # ...multiply the y-intercept bands by 10000 to match the scale of the L7 bands then apply the line equation - subtract the intercept and divide by the slope
    set('system:time_start', oli$get('system:time_start'))                                # ...set the output system:time_start metadata to the input image time_start otherwise it is null
  return(y$toShort())                                                                     # return the image as short to match the type of the other data
}

# Function to extract Landsat data within specified date range

getSRcollection = function(year, startDay, endDay, sensor, aoi) {
  # get a landsat collection for given year, day range, and sensor
  srCollection = ee.ImageCollection(paste0('LANDSAT/', sensor, '/C01/T1_SR'))$ # get surface reflectance images
    filterBounds(aoi)$                                  # ...filter them by intersection with AOI
    filterDate(paste0(year,'-',startDay), paste0(year,'-',endDay))    # ...filter them by year and day range
  
  # apply the harmonization function to LC08 (if LC08), subset bands, unmask, and resample           
  srCollection = srCollection$map(
    function(img) {
      ee_utils_pyfunc(
        dat = ee$Image(
          ee$Algorithms$If(
            sensor == 'LC08',                                               # condition - if image is OLI
            harmonizationRoy(img.unmask()),                                 # true - then apply the L8 TO L7 alignment function after unmasking pixels that were previosuly masked (why/when are pixels masked)
            img$select(c('B1', 'B2', 'B3', 'B4', 'B5', 'B7')$               # false - else select out the reflectance bands from the non-OLI image
            unmask()$                                                       # ...unmask any previously masked pixels 
            set('system:time_start', img$get('system:time_start'))          # ...set the output system:time_start metadata to the input image time_start otherwise it is null
            )
          )
        ))
    # make a cloud, shadow, snow, and water mask from fmask band
    qa = img$select('pixel_qa')                                # select the fmask band
    mask = qa$bitwiseAnd(8)$eq(0)$and(                         # exclude shadow
        qa$bitwiseAnd(16)$eq(0))$and(                          # exclude snow
        qa$bitwiseAnd(32)$eq(0)$and(                           # exclude clouds
        qa$bitwiseAnd(4)$eq(0)))                               # exclude water
    
    # apply the mask to the image and return it
    return(dat$mask(mask)) #apply the mask - 0's in mask will be excluded from computation and set to opacity=0 in display
  })
  return(srCollection); # return the prepared collection
}

# Function to merge LS collections to create mosaic by year
getCombinedSRcollection = function(year, startDay, endDay, aoi) {
  var lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi)       # get TM collection for a given year, date range, and area
  var lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi)
  if((year < 2003)||(year == 2012)){                                      # Getting LS7 only if pre-SLC failure or if 2012 is needed (gap b/w 5 and 8)
    var le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi)       # get ETM+ collection for a given year, date range, and area
    var mergedCollection = ee.ImageCollection(lt5.merge(le7).merge(lc8)) # merge the individual sensor collections into one imageCollection object
  }else{
    var mergedCollection = ee.ImageCollection(lt5.merge(lc8))
  }
  return mergedCollection                                              # return the Imagecollection
}