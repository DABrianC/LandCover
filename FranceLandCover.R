

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

map_boundary <- geoboundaries("Haiti")

# Defining filepath to save downloaded spatial file
spatial_filepath <- "LandCoverData/haiti.shp"
# Saving downloaded spatial file on  our computer
st_write(map_boundary, paste0(spatial_filepath))

MODIStsp(gui = FALSE
         , out_folder = "LandCoverData"
         , out_folder_mod = "LandCoverData"
         , selprod = "LandCover_Type_Yearly_500m (MCD12Q1)"
         , bandsel = "LC1"
         , use = "briancalhoon"
         , password = "00h0OqKWDw$67R"
         , start_date = "2020.01.01"
         , end_date = "2020.12.31"
         , verbose = FALSE
         , spatmeth = "file"
         , spafile = spatial_filepath
         , out_format = "GTIFF")


# Downloading the boundary of Zimbabwe
map_boundary <- geoboundaries("Haiti")

# Reading in the downloaded landcover raster data
IGBP_raster <- raster(here::here("LandCoverData/haiti/LandCover_Type_Yearly_500m_v6/LC1/MCD12Q1_LC1_2020_001.tif"))

# Transforming data
IGBP_raster <- projectRaster(IGBP_raster, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Cropping data
IGBP_raster <- raster::mask(IGBP_raster, as_Spatial(map_boundary))

# Converting the raster object into a dataframe and converting the IGBP classification into a factor
IGBP_df <- as.data.frame(IGBP_raster, xy = TRUE, na.rm = TRUE) %>%
  mutate(MCD12Q1_LC1_2020_001 = as.factor(round(MCD12Q1_LC1_2020_001)))
rownames(IGBP_df) <- c()
# Renaming IGBP classification levels
levels(IGBP_df$MCD12Q1_LC1_2020_001) <- c( "Evergreen needleleaf forests",
                                           "Evergreen broadleaf forests",
                                           "Deciduous needleleaf forests",
                                           "Deciduous broadleaf forests",
                                           "Mixed forests",
                                           "Closed shrublands",
                                           "Open shrublands",
                                           "Woody savannas",
                                           "Savannas",
                                           "Grasslands",
                                           "Permanent wetlands",
                                           "Croplands",
                                           "Urban and built-up lands",
                                           "Cropland/natural vegetation mosaics",
                                           "Snow and ice",
                                           "Barren",
                                           "Water bodies")
# Visualising using ggplot2
ggplot() + 
  geom_raster(data = IGBP_df,
              aes(x = x, y = y, fill = MCD12Q1_LC1_2020_001)) +
  geom_sf(data = map_boundary, inherit.aes = FALSE, fill = NA) +
  scale_fill_viridis(name = "Land Cover Type", discrete = TRUE) +
  labs(title = "Land Cover Classification in Haiti",
       subtitle = "January 1, 2020 - December 31, 2020",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()