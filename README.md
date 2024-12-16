# Habitat suitability of reef-froming species in the North Sea

## Introduction

Due to fishing and other human activities, reef forming species have almost completely disappeared over roughly the past century. They are important structures that accomodate juvenile fish and other small organisms. For protection of areas where such reefs could possibly be reintroduced, it is important to define areas that are suitable habitats. This product aims to classify areas in the North Sea based on current occurence in combination with environmental variables that are particularly suitable for these organisms.

## Directory structure

```         
EMODnet-Biology-Reef-Forming-Species-NorthSea/
├── analysis
├── data/
│   ├── derived_data/
│   ├── environment/
│   └── raw_data/
├── docs/
├── product/
└── scripts/
```

-   **analysis** - Markdown or Jupyter notebooks
-   **data** - Raw and derived data
-   **docs** - Rendered reports
-   **product** - Output product files
-   **scripts** - Reusable code

## Data series

Environmental information is needed as a basis for species distribution models. For this project, we rely heavily on a recent compilation of North Sea wide environmental information by van der Reijden et al. (2018). These authors have compiled their datasets on bathymetry, grain size distribution, temperature and salinity from diverse literature sources. They have made their data available in the form of geo-tiff files, that we have downloaded for use in the present project. In the files, there is also information on bottom shear stress, but this is based on a rather coarse model. We have replaced it with results of the Deltares DCSM-FM model for the greater North Sea. The datasets used are listed in Table II. Sources of the data are van der Reijden et al. (2018) for fisheries and calculations of 'Bathymetric Position Index' values based on bathymetry, Stephens (2015) for grain size data, Copernicus marine services (www.marine.copernicus.eu) for salinity and temperature, EMODnet bathymetry ([http://portal.emodnet-bathymetry.eu/)](http://portal.emodnet-bathymetry.eu/)) for basic bathymetry, Deltares for bottom shear stress calculated with DCSM-FM.

The 'BPI' (Bathymetric position index) calculates for each point, the difference of the depth of the point with the average depth of the surrounding area, where the surrounding area is a circle with a fixed radius. BPI5 uses 5 km as a radius for the surroundings, and similar for the other BPI variables. van der Reijden et al. (2018) also define a weighted average BPI, but we did not use that in our analysis.

Temperature difference is a measure for the change in temperature between 2008 and 2013. This is not distributed homogeneously over the North Sea. Atlantic water has warmed very little, whereas the North Sea has been warming considerably over the past decades. Consequently, the largest temperature differences are seen in the eastern and north-eastern parts of the North Sea.

No temporal (e.g. seasonal) variance of salinity and temperature has been used in the present study. It is known that variation of these variables is often very important in estuarine conditions. However, in the North Sea the ranges are much more limited. It is unlikely that any of these parameters would fall outside of the tolerance of the species, with the probable exception of temperature for the boreal species Modiolus modiolus. However, also mean temperature appeared to be a very useful variable in predicting the range of this species, and obviously there is a tight correlation between mean temperature and yearly temperature range in the North Sea.

Datasets used in analyses

| Env.Variable    | Explanation                           | Souce         |
|-----------------|---------------------------------------|---------------|
| Depth           | Depth at 178 m resolution             | EMODnet       |
| BPI5            | Bathymetric Position Index 5 km       | vdReijden2018 |
| BPI10           | Bathymetric Position Index 10 km      | vdReijden2018 |
| BPI75           | Bathymetric Position Index 75 km      | vdReijden2018 |
| Bott.shr.stress | Bottom shear stress from currents     | DCSM-FM       |
| Salinity        | Mean Salinity                         | Copernicus    |
| Temperature     | Mean Temperature                      | Copernicus    |
| Temp.diff       | Temperature differernce over the year | Copernicus    |
| Gravel          | Fraction gravel in the sediment       | Stephens2015  |
| Mud             | Fraction mud in the sediment          | Stephens2015  |
| Sand            | Fraction sand in the sediment         | Stephens2015  |
| Beam_plaice     | Intensity beam trawling for plaice    | vdReijden2018 |
| Beam_sole       | Intensity beam trawling for sole      | vdReijden2018 |
| Otter_mix       | Intensity trawling of mixed species   | vdReijden2018 |

## Extracting species presence/absence information from the available data sets

The EMODnet Biology product on presence/absence of species in samples in the Greater North Sea is delivered as a binary R file. Alternatively, it is also available as a .csv file, but this takes longer to read in.

In code chunk#3 a function is defined that retrieves the data for a particular species and writes the results as a shape file for use in GIS. Species are identified using their AphiaID, which is their unique identification in WoRMS, the World Register of Marine Species ([https://marinespecies.org)](https://marinespecies.org))

### Retrieving data from Wageningen Marine Research fisheries database

Wageningen Marine Research has made available all data in their 'Frisbee' database on the concerned species (called 'DATRAS' data). The data set is composed of all hauls with a diversity of instruments, including beam trawls, otter trawls, plankton nets and others. The species concerned were never retrieved from some of these instruments, probably because some instruments (e.g. plankton nets) are not able to catch them. In order to avoid excess zeroes, suggesting absence of the species whereas presence could not have been established, we restricted the database to those instruments that had at least once caught one of the concerned species. These are beam trawls, otter trawls and an instrument called

'GOV'. Closer examination showed that Lanice conchilega, one of the most frequently found species of macrobenthos in the North Sea, was only found 12 times in total in this database. We concluded that inclusion of the database for this species would lead to too many false zeroes, and restricted use of the database to Modiolus and Sabellaria only. The oyster was not reported from this database. However, in the retrieval code illustrated here, all four species are looked after in the DATRAS data base and illustrating shapefiles for all four are produced.

### Retrieving historical data on oyster distribution

Historical data on the distribution of oysters in the North Sea, and more particularly in the Dutch waters, during the nineteenth century were derived from Bennema et al. (2020), and courteously made available to us by Floris Bennema. These authors discuss two different sources of data in their paper. One source are historical expeditions in the North Sea, the data of which have been digitized. We received these data in two files: one file describing finds by the Huxley_Wodan expeditions, and one by the Poseidon expeditions. These data have been read in and converted to spatial files. The other source were old maps, that have been critically evaluated by the authors and compiled into an overall map indicating the area of high oyster occurrence in the region around the Oyster Grounds. We digitized this map into a polygon using QGIS and used it as a basis to generate pseudo-absences and pseudo- presences. Random points were generated in the North Sea, and points within the map polygon were attributed a probability of 0.7 to contain oysters, whereas points outside of the polygon had absence. In order to complement this data base with information on the Flemish Banks, that could also be of importance to the Dutch waters off Zeeland, we used the report by Houziaux et al. (2008) on the findings of the extensive set of dredge surveys by Gilson in the beginning of the twentieth century. We digitized all sample points of Gilson from the figures in the report of Houziaux, indicating presence of oysters where this was recorded. The points were saved as a shapefile and read in to extract the points with absence and presence information. The data provided by Wageningen Marine Research also contain findings of flat oysters in the Voordelta, the Rotterdam harbour and the Wadden Sea. All three of these populations fall outside of the environmental rasters available in the present project. Two of them seem to depend on artificial hard substrate, although it remains to be seen if that is only a transition phase or not. It is also likely that in these estuarine or near-estuarine conditions, other environmental factors (e.g. salinity) will have an influence on habitat suitability than in the open North Sea. For these reasons, information from these populations was not used in the present analysis, which was restricted to historical data of oyster occurrence on natural substrates.

## Data product

The random forest regression analysis generated maps with suitable habitats for the 4 species regarded in this analysis.

## More information:

The Rmarkdown file "Description_analysis.Rmd" describes the complete analysis, with code.

### References

### Code and methodology

#### Regression analysis

Species distribution models have been prepared using random forest regression. The fisheries intensity was not relevant as a predictor for the oyster, as the oyster data are historical nineteenth-century reconstructions. It turned out that for the other three species, the predictive power of the three fisheries intensities was very low. The variables have been removed from the analysis. Furthermore, sand fraction has also been removed from the analysis, because it is fully collinear with mud and gravel fractions: the three together always sum to 1. From the BPI variables, we only retained BPI at 5, 10 and 75 km, as the other classes (30 and 50 km) were usually redundant with these three. The remaining variables all had at least some importance in almost all regression models. If a single factor was occasionally not statistically relevant, it was still maintained in the analysis in order to keep consistency between the different species. n the random forest model, there were clear signs of overfitting in the *Ostrea* model, when only the expedition data were used. Overfitting was manifested because the prediction model only predicted occurrence in a very narrow band around the positive observations, not in between them. This defect was much less apparent after we added the pseudo-data based on the historical maps. For the other random forest models, no clear signs of overfitting were apparent, although it might sometimes be the case in the *Modiolus* map.

### Citation and download link

This product should be cited as:

P.M.J.Herman, and F.F. van Rees. 2022. “Mapping Reef Forming North Sea Species.” Delft. <https://pub.kennisbank.deltares.nl/Details/fullCatalogue/1000020826>.

Available to download in:

{{link_download}}

### Authors

Peter Herman & Floris van Rees

contact: willem.stolte\@deltares.nl
