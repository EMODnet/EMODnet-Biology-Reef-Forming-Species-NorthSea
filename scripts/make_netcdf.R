

library(RNetCDF)
library(readr)
library(dplyr)
library(glue)
library(raster)

# projection strings
proWG  <- CRS("+proj=longlat +datum=WGS84")
proUTM <- CRS("+proj=utm +zone=31 +ellps=GRS80 +units=m +no_defs")

#== make Specieslist ================================ 

source('scripts/make_taxon_list.R')

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

# Extract the unique and sorted values of the 4 dimensions for later use 
lon <- sort(unique(raster::coordinates(dataset.list.wg)[,1]))
lat <- sort(unique(raster::coordinates(dataset.list.wg)[,2]))

# no time dimension

AphiaID <- taxon$AphiaID_species

# convert rasters to dataframes
dataset2 <- lapply(
  1:4, 
  function(ii) 
    as.data.frame(dataset.list.wg[[ii]], xy=TRUE) %>%
    mutate(AphiaID = taxon$AphiaID_species[ii]) %>%
    rename(lon = x,
           lat = y,
           value = 3) %>%
    arrange(lat, lon)
) %>%
  bind_rows()


#== create 3d array =============================================

# The product of the lengths of all dimensions
length(unique(dataset2$lon)) * length(unique(dataset2$lat)) * length(unique(dataset2$AphiaID)) == nrow(dataset2)

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

dim(array)


#== create netcdf file =========================================================

nc <- create.nc("product/habitat_suitability_reef_forming_species_north_sea.nc") 

#== add lon and Lat coordinates ================================================

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

#=== Define time dimension =====================================================
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
att.put.nc(nc, variable = "aphiaid", name = "units", type = "NC_CHAR", value = "level")

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
var.put.nc(nc, variable = "taxon_lsid", data = paste0("urn:lsid:marinespecies.org:taxname:", taxon$Species))

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
att.put.nc(nc, variable = "habitat_suitability", name = "long_name", type = "NC_CHAR", value = "Habitat suitability of biological entity specified elsewhere per unit area of the bed")

var.put.nc(nc, variable = "habitat_suitability", data = array) 

#== add global attributes ======================================================

attributes <- list(
  title = "Habitat suitability for reef-forming species in the North Sea.",
  summary = "Due to fishing and other human activities, reef forming species have almost completely disappeared over roughly the past century. They are important structures that accomodate juvenile fish and other small organisms. For protection of areas where such reefs could possibly be reintroduced, it is important to define areas that are suitable habitats. This product aims to classify areas in the North Sea based on current occurence in combination with environmental variables that are particularly suitable for these organisms.",                       
  Conventions = "CF-1.8",
  # id = "",
  naming_authority = "emodnet-biology.eu",
  history = "https://github.com/EMODnet/EMODnet-Biology-Reef-Forming-Species-NorthSea",
  source = "https://github.com/EMODnet/EMODnet-Biology-Reef-Forming-Species-NorthSea",
  # processing_level = "",
  # comment = "", 
  # acknowledgment = "",
  license = "CC-BY",
  standard_name_vocabulary = "CF Standard Name Table v1.8",
  date_created = as.character(Sys.Date()),
  creator_name = "Willem Stolte",
  creator_email = "willem.stolte@deltares.nl",
  creator_url = "www.deltares.nl",
  institution = "Deltares",
  project = "EMODnet-Biology",
  publisher_name = "EMODnet-Biology",                 
  publisher_email = "bio@emodnet.eu",                
  publisher_url = "www.emodnet-biology.eu",                  
  # geospatial_bounds = "",              
  # geospatial_bounds_crs = "",          
  # geospatial_bounds_vertical_crs = "", 
  geospatial_lat_min = min(lat),
  geospatial_lat_max = max(lat),
  geospatial_lon_min = min(lon),
  geospatial_lon_max = max(lon),
  # geospatial_vertical_min = "",        
  # geospatial_vertical_max = "",        
  # geospatial_vertical_positive = "",  
  # time_coverage_start = "1911",            
  # time_coverage_end = "2016",              
  # time_coverage_duration = "",         
  # time_coverage_resolution = "",       
  # uuid = "",                           
  # sea_name = "",                       
  # creator_type = "",                   
  creator_institution = "Deltares",            
  # publisher_type = "",                 
  publisher_institution = "Flanders Marine Institute (VLIZ)",        
  # program = "",                        
  # contributor_name = "",               
  # contributor_role  = "",              
  geospatial_lat_units = "degrees_north",           
  geospatial_lon_units = "degrees_east",           
  # geospatial_vertical_units   = "",    
  # date_modified = "",               
  # date_issued = "",                    
  # date_metadata_modified   = "",       
  # product_version = "",            
  # keywords_vocabulary = "",          
  # platform  = "",              
  # platform_vocabulary = "",          
  # instrument = "",          
  # instrument_vocabulary  = "",        
  # featureType = "Point",                  
  # metadata_link = "",                  
  # references = "",
  comment = "Uses attributes recommended by http://cfconventions.org",
  license = "CC-BY", 
  publisher_name = "EMODnet Biology Data Management Team",
  citation = "P.M.J.Herman, and F.F. van Rees. 2022. “Mapping Reef Forming North Sea Species.” Delft. https://pub.kennisbank.deltares.nl/Details/fullCatalogue/1000020826.",
  acknowledgement = "European Marine Observation Data Network (EMODnet) Biology project (EMFF/2019/1.3.1.9/Lot 6/SI2.837974), funded by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund"
)

# Define function that detects if the data type should be character of 
# integer and add to global attributes
add_global_attributes <- function(nc, attributes){
  
  stopifnot(is.list(attributes))
  
  for(i in 1:length(attributes)){
    if(is.character(attributes[[i]])){
      type <- "NC_CHAR"
    }else if(is.numeric(attributes[[i]])){
      type <- "NC_DOUBLE"
    }
    att.put.nc(nc, variable = "NC_GLOBAL", name = names(attributes[i]), type = type, value = attributes[[i]])
  }
  sync.nc(nc)
}

# Add attributes
add_global_attributes(nc, attributes)

sync.nc(nc)
print.nc(nc)
close.nc(nc)

