

#load libraries
library(tidyverse)
library(MODIStsp)
library(rgeoboundaries)
library(sf)
library(raster)
library(rasterVis)
library(here)
library(viridis)


#get land cover layers

MODIStsp_get_prodlayers("MCD12Q1")

map_boundary <- geoboundaries("France")

plot(map_boundary)

# Defining filepath to save downloaded spatial file
spatial_filepath <- "LandCoverData/france.shp"
# Saving downloaded spatial file on to our computer
st_write(map_boundary, paste0(spatial_filepath))

MODIStsp(gui             = FALSE,
         out_folder      = "LandCoverData",
         out_folder_mod  = "LandCoverData",
         selprod         = "LandCover_Type_Yearly_500m (MCD12Q1)",
         bandsel         = "LC1", 
         user            = "briancalhoon" ,
         password        = "00h0OqKWDw$67R",
         start_date      = "2020.01.01", 
         end_date        = "2020.12.31", 
         verbose         = FALSE,
         spatmeth        = "file",
         spafile         = spatial_filepath,
         out_format      = "GTiff")

IGBP_raster <- raster(here::here("LandCoverData/france/LandCover_Type_Yearly_500m_v6/LC1/MCD12Q1_LC1_2020_001.tif"))
