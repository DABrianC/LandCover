

#load libraries
library(tidyverse)
library(rgeoboundaries)
library(sf)
library(raster)
library(rasterVis)
library(here)
library(viridis)
library(MODIStsp)
library(extrafont)
library(extrafontdb)

#get land cover layers

MODIStsp_get_prodlayers("MCD12Q1")

map_boundary<- geoboundaries("HTI")



# Defining filepath to save downloaded spatial file
spatial_filepath <- "LandCoverData/haiti.shp"
# Saving downloaded spatial file on  our computer
st_write(map_boundary, paste0(spatial_filepath), append = FALSE)

MODIStsp(gui = FALSE
         , out_folder = "LandCoverData"
         , out_folder_mod = "LandCoverData"
         , selprod = "LandCover_Type_Yearly_500m (MCD12Q1)"
         , bandsel = "LC1"
         , user = "" #insert your own username
         , password = "" #insert your own password
         , start_date = "2019.01.01"
         , end_date = "2019.12.31"
         , verbose = FALSE
         , spatmeth = "file"
         , spafile = spatial_filepath
         , out_format = "GTiff")


# Reading in the downloaded landcover raster data
IGBP_raster <- raster(here::here("LandCoverData/haiti/LandCover_Type_Yearly_500m_v6/LC1/MCD12Q1_LC1_2019_001.tif"))

# Transforming data
IGBP_raster <- projectRaster(IGBP_raster, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Cropping data
IGBP_raster <- raster::mask(IGBP_raster, as_Spatial(map_boundary))

# Converting the raster object into a dataframe and converting the IGBP classification into a factor
IGBP_df <- as.data.frame(IGBP_raster, xy = TRUE, na.rm = TRUE) %>%
  mutate(MCD12Q1_LC1_2019_001 = as.factor(round(MCD12Q1_LC1_2019_001)))
rownames(IGBP_df) <- c()
# Renaming IGBP classification levels
levels(IGBP_df$MCD12Q1_LC1_2019_001) <- c( "Evergreen needleleaf forests",
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
p <- ggplot() + 
    geom_raster(data = IGBP_df,
              aes(x = x, y = y, fill = MCD12Q1_LC1_2019_001)) +
  geom_sf(data = map_boundary, inherit.aes = FALSE, fill = NA) +
  scale_fill_viridis(name = "Land Cover Type", discrete = TRUE) +
  labs(title = "Land Cover\nClassification\nin Haiti",
      subtitle = "January 1, 2019 - December 31, 2019",
      x = "Longitude",
      y = "Latitude"
      , caption = "Author: Brian Calhoon | Data: MODIS | Packages: rgeoboundaries, MODIStsp") +
  theme_void() +
  theme(panel.background = element_rect(fill = "#E7CEC8",
                                        colour = "#E7CEC8",
                                        size = 0.5, linetype = "solid")
        , plot.title.position = "plot"
        , plot.title = element_text(size = 44
                                    , family = "Corbel"
                                    , color = "#04366D"
                                    , vjust = - 20)
        , plot.subtitle = element_text(size = 28
                                       , family = "Corbel"
                                       , color = "#0960BE"
                                       , vjust = -35)
        , legend.title = ggplot2::element_text(size = 20
                                               , family = "Corbel"
                                               , color = "#000000")
        , legend.key.size = unit(1, "cm")
        , legend.text = ggplot2::element_text(size = 16
                                             , family = "Corbel"
                                             , color = "#000000")
        , plot.caption = ggplot2::element_text(size = 16, family = "Corbel", color = "#04366D")) +
  guides(color=guide_legend(ncol=2))
  png("Haiti.png"
      , width = 1080
      , height = 920
      , unit = "px")
  print(p)
  dev.off()

