

library(RNetCDF)
library(readr)
library(dplyr)
library(glue)
library(raster)

proWG<-CRS("+proj=longlat +datum=WGS84")
proUTM <- CRS("+proj=utm +zone=31 +ellps=GRS80 +units=m +no_defs")


# dataset = raster()
# load('product/rf_model_Sabellaria.Rdata')

#== make Specieslist ================================ 

source('make_taxon_list.R')
taxon <- taxon %>% 
  mutate(
    filename = case_when(
      grepl("sabellaria", Species, ignore.case = TRUE) ~ 'RF_pred_Sabellaria.tif',
      grepl("lanice", Species, ignore.case = TRUE) ~ 'RF_pred_Lanice.tif',
      grepl("modiolus", Species, ignore.case = TRUE) ~ 'RF_pred_Modiolus.tif',
      grepl("ostrea", Species, ignore.case = TRUE) ~ 'RF_pred_Ostrea.tif'
    )
  )


# or all in one go

dataset.list <- brick(lapply(taxon$filename, function(x) raster(file.path('product', x))))
names(dataset.list) <- taxon$Species

dataset.list.wg <- raster::projectRaster(from = dataset.list, crs = proWG)
plot(dataset.list.wg)
plot(dataset.list)

# Extract the unique and sorted values of the 4 dimensions
lon <- sort(unique(raster::coordinates(dataset.list.wg)[,1]))
lat <- rev(sort(unique(raster::coordinates(dataset.list.wg)[,2])))

# no time dimension

AphiaID <- taxon$AphiaID_species

# # Use expand.grid() to create a data frame with all the possible 
# # combinations of the 4 dimensions
# longer <- expand.grid(
#   lon = lon, 
#   lat = lat, 
#   # time = time, 
#   aphiaid = AphiaID, 
#   stringsAsFactors = FALSE)

# Define unique identifier again and merge the variable occurrenceStatus 
# with presences and absences

dataset2 <- lapply(
  1:4, 
  function(ii) 
    as.data.frame(dataset.list.wg[[ii]], xy=TRUE) %>%
    mutate(AphiaID = taxon$AphiaID_species[ii]) %>%
    rename(lon = x,
           lat = y,
           value = 3)
) %>%
  bind_rows()

# Save for later
# write_csv(dataset2, "./data/derived/all_products.csv")


#== create 3d array =============================================

# The product of the lengths of all dimensions
length(unique(dataset2$lon)) * length(unique(dataset2$lat)) * length(unique(dataset2$AphiaID))
#> 

# Is the same as the length of the variable of interest, including all 
# possible combinations of the dimensions even if this coerce NA's

# Create array
array <- array(
  data = dataset2$value,
  dim = c(
    length(unique(dataset2$lon)), 
    length(unique(dataset2$lat)), 
    length(unique(dataset2$AphiaID)))
)

#== create netcdf file =========================================

nc <- create.nc("product/foo.nc") 

#== add lon coordinates =================

# Define lon dimension
dim.def.nc(nc, dimname = "lon", dimlength = length(lon)) 

# Define lon variable
var.def.nc(nc, varname = "lon", vartype = "NC_DOUBLE", dimensions = "lon")

# Add attributes
att.put.nc(nc, variable = "lon", name = "units", type = "NC_CHAR", value = "degrees_east")
att.put.nc(nc, variable = "lon", name = "standard_name", type = "NC_CHAR", value = "longitude")
att.put.nc(nc, variable = "lon", name = "long_name", type = "NC_CHAR", value = "Longitude")

# Put data
var.put.nc(nc, variable = "lon", data = lon) 

# Check
var.get.nc(nc, variable = "lon")

# Define lat dimension
dim.def.nc(nc, dimname = "lat", dimlength = length(lat)) 

# Define lat variable
var.def.nc(nc, varname = "lat", vartype = "NC_DOUBLE", dimensions = "lat")

# Add attributes
att.put.nc(nc, variable = "lat", name = "units", type = "NC_CHAR", value = "degrees_north")
att.put.nc(nc, variable = "lat", name = "standard_name", type = "NC_CHAR", value = "latitude")
att.put.nc(nc, variable = "lat", name = "long_name", type = "NC_CHAR", value = "Latitude")

# Put data
var.put.nc(nc, variable = "lat", data = lat) 

# Check
var.get.nc(nc, variable = "lat")


# # Define time dimension
# dim.def.nc(nc, dimname = "time", dimlength = length(time)) 
# 
# # Define time variable
# var.def.nc(nc, varname = "time", vartype = "NC_DOUBLE", dimensions = "time")
# 
# # Add attributes
# att.put.nc(nc, variable = "time", name = "standard_name", type = "NC_CHAR", value = "time")
# att.put.nc(nc, variable = "time", name = "long_name", type = "NC_CHAR", value = "Time")
# att.put.nc(nc, variable = "time", name = "units", type = "NC_CHAR", value = "days since 1970-01-01 00:00:00")
# att.put.nc(nc, variable = "time", name = "calendar", type = "NC_CHAR", value = "gregorian")
# 
# # Put data
# var.put.nc(nc, variable = "time", data = time)
# 
# # Check
# var.get.nc(nc, variable = "time")

# Define the aphia and string80 dimensions
dim.def.nc(nc, dimname = "aphiaid", dimlength = nrow(taxon))
dim.def.nc(nc, dimname = "string80", dimlength = 80)

# Add aphiaid variable and attributes 
var.def.nc(nc, varname = "aphiaid", vartype = "NC_INT", dimensions = "aphiaid")
att.put.nc(nc, variable = "aphiaid", name = "long_name", type = "NC_CHAR", value = "Life Science Identifier - World Register of Marine Species")

# Put aphiaid data
var.put.nc(nc, variable = "aphiaid", data = taxon$AphiaID_species)

# Check
var.get.nc(nc, variable = "aphiaid")

# Add taxon_name variable and attributes
var.def.nc(nc, varname = "taxon_name", vartype = "NC_CHAR", dimension = c("string80", "aphiaid"))
att.put.nc(nc, variable = "taxon_name", name = "standard_name", type = "NC_CHAR", value = "biological_taxon_name")
att.put.nc(nc, variable = "taxon_name", name = "long_name", type = "NC_CHAR", value = "Scientific name of the taxa")

# Put taxon_name data
var.put.nc(nc, variable = "taxon_name", data = taxon$Species)

# Check
var.get.nc(nc, variable = "taxon_name")


# Add taxon_lsid variable and attributes
var.def.nc(nc, varname = "taxon_lsid", vartype = "NC_CHAR", dimension = c("string80", "aphiaid"))
att.put.nc(nc, variable = "taxon_lsid", name = "standard_name", type = "NC_CHAR", value = "biological_taxon_lsid")
att.put.nc(nc, variable = "taxon_lsid", name = "long_name", type = "NC_CHAR", value = "Life Science Identifier - World Register of Marine Species")

# Put taxon_name data
var.put.nc(nc, variable = "taxon_lsid", data = taxon$Species)

# Check
var.get.nc(nc, variable = "taxon_lsid")


# Define non-dimensional crs variable 
var.def.nc(nc, varname = "crs", vartype = "NC_CHAR", dimensions = NA)

# Add attributes
att.put.nc(nc, variable = "crs", name = "long_name", type = "NC_CHAR", value = "Coordinate Reference System")
att.put.nc(nc, variable = "crs", name = "geographic_crs_name", type = "NC_CHAR", value = "WGS 84")
att.put.nc(nc, variable = "crs", name = "grid_mapping_name", type = "NC_CHAR", value = "latitude_longitude")
att.put.nc(nc, variable = "crs", name = "reference_ellipsoid_name", type = "NC_CHAR", value = "WGS 84")
att.put.nc(nc, variable = "crs", name = "horizontal_datum_name", type = "NC_CHAR", value = "WGS 84")
att.put.nc(nc, variable = "crs", name = "prime_meridian_name", type = "NC_CHAR", value = "Greenwich")
att.put.nc(nc, variable = "crs", name = "longitude_of_prime_meridian", type = "NC_DOUBLE", value = 0.)
att.put.nc(nc, variable = "crs", name = "semi_major_axis", type = "NC_DOUBLE", value = 6378137.)
att.put.nc(nc, variable = "crs", name = "semi_minor_axis", type = "NC_DOUBLE", value = 6356752.314245179)
att.put.nc(nc, variable = "crs", name = "inverse_flattening", type = "NC_DOUBLE", value = 298.257223563)
att.put.nc(nc, variable = "crs", name = "spatial_ref", type = "NC_CHAR", value = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]')
att.put.nc(nc, variable = "crs", name = "GeoTransform", type = "NC_CHAR", value = '-180 0.08333333333333333 0 90 0 -0.08333333333333333 ')


# Create the presence_absence variable defined by the four dimensions
var.def.nc(nc, varname = "habitat_suitability", vartype = "NC_DOUBLE", dimensions = c("lon", "lat", "aphiaid"))

# Add attributes
att.put.nc(nc, variable = "habitat_suitability", name = "_FillValue", type = "NC_DOUBLE", value = -99999)
att.put.nc(nc, variable = "habitat_suitability", name = "long_name", type = "NC_CHAR", value = "Modelled habitat suitability for a species.")

var.put.nc(nc, variable = "habitat_suitability", data = array) 

sync.nc(nc)
print.nc(nc)
close.nc(nc)

