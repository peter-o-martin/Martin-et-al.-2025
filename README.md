### Spatial-Temporal Patterns of Per- and Polyfluoroalkyl Substances (PFAS) in the Biota of the Laurentian Great Lakes: A Meta-Analysis
## Leyerle Martin et al., 2025 

This GitHub repository contains the R code, finalized dataset (i.e., the dataset produced after performing imputation of missing values and adding in additional variables), and the generalized additive models produced through the methodological workflow described in detail in Leyerle Martin et al. (2025). The different parts of the repository are as follows:

1. **Great Lakes Shapefiles** folder: this folder contains the finalized shapefile, GL_Watershed_shapefile.shp, that was used to assign data points to one of the five watersheds of the Great Lakes, as well as to remove any points that fell outside the Great Lakes watersheds. This finalized vector layer, created with QGIS software version 3.36.0 – Maidenhead (QGIS Development Team, 2024), cobmines information from shapefiles characterizing the Great Lakes subbasins (US Geological Survey, 2010) and the Saint Lawrence River (Flanders Marine Institute, 2017; Natural Resources Canada & US Geological Survey, 2010). We also buffered the subbasins shapefile in QGIS by approximately 48 km (0.311°), a distance that corresponds to the estimated summer range limits of the bald eagle (Haliaeetus leucocephalus) (Mandernack et al., 2012), the most mobile species included in the meta-analysis.





## References
Flanders Marine Institute. (2017). *Gulf of Saint Lawrence* [Shapefile]. Marine Regions (the VLIMAR Gazetteer and the VLIZ Maritime Boundaries Geodatabase). http://marineregions.org/mrgid/4290
Mandernack, B. A., Solensky, M., Martell, M., & Schmitz, R. T. (2012). Satellite Tracking of Bald Eagles in the Upper Midwest. *Journal of Raptor Research, 46*(3), 258–273. https://doi.org/10.3356/JRR-10-77.1
Natural Resources Canada & US Geological Survey. (2010). *North American Atlas – Basin Watersheds* [Shapefile]. Government of Canada. http://www.cec.org/north-american-environmental-atlas/watersheds/
QGIS Development Team. (2024). *QGIS Geographic Information System* (Version 3.36.0 – Maidenhead) [Computer software]. QGIS Association. https://www.qgis.org
US Geological Survey. (2010). *Great Lakes and Watersheds Shapefiles* [Shapefile]. USGS ScienceBase-Catalog. https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd
