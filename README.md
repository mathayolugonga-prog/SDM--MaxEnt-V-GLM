
# Species Distribution Modelling of Eurasian Badger (Meles meles)

#Project Overview
This project applies species distribution modelling (SDM) techniques to analyse and predict the spatial distribution of the Eurasian badger (*Meles meles*) across Scotland. The study integrates environmental predictors including land cover and elevation to understand habitat suitability patterns.

Two modelling approaches are used:
- Generalised Linear Model (GLM)
- Maximum Entropy Model (MaxEnt)

The project demonstrates practical applications of spatial ecology, GIS analysis, and statistical modelling in R.

---

#Data Used
- Species occurrence data (CSV format)
- Land Cover Map (LCM raster)
- Digital Elevation Model (DEM raster)
- Study area boundary (Shapefile)

- Environmental predictors included:
- Proportion of broadleaf woodland (1800m radius)
- Proportion of urban land cover (2300m radius)
- Elevation (metres)

# Methods

#1. Data Preparation
- Cleaned occurrence records based on coordinate uncertainty
- Converted tabular data into spatial objects (sf)
- Cropped and masked raster datasets to study area
- Standardised coordinate reference systems (CRS)

#2. Feature Engineering
- Reclassified land cover into:
  - Broadleaf woodland (binary)
  - Urban land cover (binary)
- Applied focal analysis to compute landscape composition at different spatial scales
- Generated environmental raster layers for modelling

- #3. Background Data Generation
- Created 2000 random background (pseudo-absence) points
- Extracted environmental variables for presence and background locations

#4. Model Development

#GLM (Logistic Regression)
- Included polynomial terms to capture non-linear species responses
- Used spatial cross-validation (mlr package)
- Evaluated performance using AUC

#MaxEnt Model
- Implemented using `maxnet`
- Used spatial block cross-validation
- Model evaluated using AUC values across folds

- 
#Outputs
The project produces:
- Habitat suitability maps (GLM and MaxEnt)
- Response curves for environmental variables
- Spatial cross-validation visualisations
- Comparative modelling results
---

#How to Run

#Requirements
Install required R packages:

```r
install.packages(c("terra", "sf", "ggplot2", "cowplot", "mlr", "dismo", "maxnet", "glmnet", "precrec"))
