#SPECIES DISTRIBUTION MODELLING: ASSESSMENT ONE

# Setting working directory
setwd("D:/MSc. GIS/GEOG71922- Spatial Ecology/Ass.1/Assessment/Sciurus")

# Installing required packages
install.packages("terra")
install.packages("maxnet")
install.packages("cowplot")
install.packages("ggplot2")
install.packages("precrec")
install.packages("glmnet")
install.packages("mlr")

# Load core spatial libraries
library(terra)    # for spatial raster data
library(sf)       # for spatial vector data (data frames with geometry)
library(ggplot2)  # for response plots
library(cowplot)  # for arranging multiple plots



# Loading of occurrence data in CSV formart
melesmeles = read.csv("Melesmeles.csv")

# Data preparation and cleaning; Subset to records with coordinate uncertainty under/equal 1000 metres
melesmeles = melesmeles[melesmeles$Coordinate.uncertainty_m <= 1000, ]

# Extracting coordinate pairs from occurance data (longitude and latitude columns)
melesmeles.latlong = data.frame(x = melesmeles$Longitude, y = melesmeles$Latitude)

# Creating spatial presentation in points and assign data to crs (WGS84 (lat/lon))
melesmeles.sp = st_as_sf(melesmeles.latlong, coords = c("x", "y"), crs = "epsg:4326")

# Load a study area (Scotland border in .shp) and environment data (Land cover in raster
scot = st_read('scotSamp.shp')

# Load the UK Land Cover Map raster
LCM = rast("LCMUK.tif")

# Crop LCM to a buffered study area (buffer avoids edge data loss in next step)
LCM = crop(LCM, st_buffer(scot, dist = 1000))



#Initial data processing

#Raster layer aggregation to enhance PC processing efficiency
LCM = aggregate(LCM$LCMUK_1, fact = 4, fun = "modal") # Increase pixel size 25m>>100m and apply 'modal' to keep land cover categories

melesmeles.sp = st_transform(melesmeles.sp, crs(LCM)) # Reproject points to match the CRS of the LCM raster

melesmelesFin = melesmeles.sp[scot, ]                 # Crop points to the study area boundary
cat("Records after cropping to study area:", nrow(melesmelesFin), "\n")

LCM = crop(LCM, scot, mask = TRUE)                    # Mask the LCM to the study area boundary



# Preparing environment predictors from the land cover.

# Treating LCM as categorical (factor) for reclassification
LCM = as.factor(LCM$LCMUK_1)
levels(LCM)

# Creating reclassification vector: 1 for broadleaf class, 0 for all others
reclass = c(0, 1, rep(0, 20))


# Building reclassification matrix from LCM levels and new values
RCmatrix = cbind(levels(LCM)[[1]], reclass)
RCmatrix = RCmatrix[, 2:3]
RCmatrix = apply(RCmatrix, 2, FUN = as.numeric)

# Applying reclassification to produce binary broadleaf woodland raster
broadleaf = classify(LCM, RCmatrix)



# Focal Analysis- at 1800m scale 

# Number of pixels needed to cover 1800m radius
nPix = round(1800 / res(LCM)[1])
# Converting radius to diameter and ensure matrix dimensions are odd
nPix = (nPix * 2) + 1

# Building a base weights matrix
weightsMatrix = matrix(1:nPix^2, nrow = nPix, ncol = nPix)

# Setting up the focal (central) cell
x = ceiling(ncol(weightsMatrix) / 2)
y = ceiling(nrow(weightsMatrix) / 2)
focalCell = weightsMatrix[x, y]
indFocal = which(weightsMatrix == focalCell, arr.ind = TRUE)

# Calculating Euclidean distance from each cell to the focal cell
distances = list()
for(i in 1:nPix^2){
  ind.i = which(weightsMatrix == i, arr.ind = T)
  diffX = abs(ind.i[1,1] - indFocal[1,1]) * res(LCM)[1]
  diffY = abs(ind.i[1,2] - indFocal[1,2]) * res(LCM)[1]
  dist.i = sqrt(diffX^2 + diffY^2)
  distances[[i]] = dist.i
}

# Filling up matrix with computed distances
weightsMatrix[] = unlist(distances)

# Excluding cells beyond 1800m by assigning them to 'NA'
weightsMatrix[weightsMatrix > 1800] = NA

# Visualise the circular neighbourhood weights matrix
plot(rast(weightsMatrix), main = "Broadleaf Neighbourhood Weights Matrix (1800m radius)")

# Standardising matrix weights so they sum to 1 (results will be proportions)
weightsMatrixNorm = weightsMatrix
weightsMatrixNorm[!is.na(weightsMatrixNorm)] = 1 / length(weightsMatrixNorm[!is.na(weightsMatrixNorm)])

# Visualise the normalised weights matrix
plot(rast(weightsMatrixNorm), main = "Normalised Broadleaf Weights Matrix")

# Applying focal analysis to compute proportion of broadleaf per raster cell
lcm_wood_1800 = focal(broadleaf, w = weightsMatrixNorm, fun = "sum")

# Plot broadleaf proportion surface with presence points
plot(lcm_wood_1800, main = "Proportion of Broadleaf Woodland within 1800m")
plot(melesmelesFin$geometry, add = T, pch = 16, cex = 0.5, col = "red")



# Focal at 2300m scale for Urban cover

# Reclassifying urban and suburban class with 1 and 0 for all others
reclassUrban = c(rep(0, 19), 1, 1)

# Building reclassification matrix
RCmatrixUrban = cbind(levels(LCM)[[1]], reclass)
RCmatrixUrban = RCmatrixUrban[, 2:3]
RCmatrixUrban = apply(RCmatrixUrban, 2, FUN = as.numeric)
urban = classify(LCM, RCmatrixUrban)

# Building weights matrix for 2300m characteristic scale
nPixUrban = round(2300 / res(LCM)[1])
nPixUrban = (nPixUrban * 2) + 1
weightsMatrixUrban = matrix(1:nPixUrban^2, nrow = nPixUrban, ncol = nPixUrban)

#Deciding on the focal cell
x = ceiling(ncol(weightsMatrixUrban) / 2)
y = ceiling(nrow(weightsMatrixUrban) / 2)
focalCell = weightsMatrixUrban[x, y]
indFocal = which(weightsMatrixUrban == focalCell, arr.ind = TRUE)

# Calcultating eucledian distances
distancesUrban = list()
for(i in 1:nPixUrban^2){
  ind.i = which(weightsMatrixUrban == i, arr.ind = T)
  diffX = abs(ind.i[1,1] - indFocal[1,1]) * res(LCM)[1]
  diffY = abs(ind.i[1,2] - indFocal[1,2]) * res(LCM)[1]
  dist.i = sqrt(diffX^2 + diffY^2)
  distancesUrban[[i]] = dist.i
}

# Filling weights along the distances
weightsMatrixUrban[] = unlist(distancesUrban)
weightsMatrixUrban[weightsMatrixUrban > 2300] = NA
weightsMatrixUrban[!is.na(weightsMatrixUrban)] = 1 / length(weightsMatrixUrban[!is.na(weightsMatrixUrban)])

# Applying focal analysis to compute proportion of urban cover per cell
lcm_urban_2300 = focal(urban, w = weightsMatrixUrban, fun = "sum")

# Plotting urban proportion surface with presence points
plot(lcm_urban_2300, main = "Proportion of Urban Cover within 2300m")
plot(melesmelesFin$geometry, add = T, pch = 16, cex = 0.5, col = "red")



# Adding elevation data and processing it

# Reading in elevation layer
demScot = rast('demScotland.tif')

# Applying resample to a DEM layer to match resolution and extent of a land cover layer
demScot = terra::resample(demScot, lcm_wood_1800)

# Stacking all environmental data to prepare predictors
allEnv = c(lcm_wood_1800, lcm_urban_2300, demScot)
names(allEnv) = c("broadleaf", "urban", "elevation")



# Generating pseudo absence/ background points

set.seed(11)

# Randomly sample 2000 background points across the study area
back = spatSample(allEnv, size = 2000, as.points = TRUE, method = "random", na.rm = TRUE)
back = back[!is.na(back$broadleaf), ]
back = st_as_sf(back, crs = "EPSG:27700")

# Plotting presence points and background points
plot(allEnv$broadleaf, main = "Presence Points (red) and Background Points (blue)")
plot(melesmelesFin$geometry, add = T, pch = 16, cex = 0.7, col = "black")
plot(back$geometry, add = T, pch = 16, cex = 0.7, col = "red")



# Extraction of predictors values from presence/absence occurances

# Extraction of the environmental values at presence locations
eP = terra::extract(allEnv, melesmelesFin)

# Combining extracted values with presence geometry and assign presence label
Pres.cov = st_as_sf(cbind(eP, melesmelesFin))
Pres.cov$Pres = 1

Pres.cov = Pres.cov[, -1]   # remove the auto-generated ID column

# Retaining coordinates for spatial cross-validation for model evaluation
coordsPres = st_coordinates(Pres.cov)

# Assigning absence to background points
Back.cov = st_as_sf(data.frame(back, Pres = 0))
coordsBack = st_coordinates(back)

# Combining presence and background coordinates
coords = data.frame(rbind(coordsPres, coordsBack))
colnames(coords) = c("x", "y")

# Combining all presence and background data and storing them into one data frame
all.cov = rbind(Pres.cov, Back.cov)
all.cov = cbind(all.cov, coords)
all.cov = na.omit(all.cov)

# Droping off geometry column to work with tabular data
all.cov = st_drop_geometry(all.cov)


# Displaying distributions of covariate values for presence/background
par(mfrow = c(1, 3))
boxplot(broadleaf ~ Pres, data = all.cov,
        names = c("Background", "Presence"),
        main = "Broadleaf: Presence vs Background",
        ylab = "Proportion Broadleaf (1800m)",
        col = c("grey", "darkgreen"))

boxplot(urban ~ Pres, data = all.cov,
        names = c("Background", "Presence"),
        main = "Urban: Presence vs Background",
        ylab = "Proportion Urban (2300m)",
        col = c("grey", "darkred"))

boxplot(elevation ~ Pres, data = all.cov,
        names = c("Background", "Presence"),
        main = "Elevation: Presence vs Background",
        ylab = "Elevation (metres)",
        col = c("grey", "steelblue"))
par(mfrow = c(1, 1))



# Set up of MLR for Spatial Cross- Validation 

library(mlr)   # Loading MLR library

# Converting response variable to factor for creating task
task = all.cov
task$Pres = as.factor(task$Pres)

# Creating the MLR classification task, specifying target column and coordinates
task = makeClassifTask(data = task[, c(1:4)], target = "Pres",
                       positive = "1", coordinates = task[, 5:6])

# Defining the binomial (logistic regression) learner
lrnBinomial = makeLearner("classif.binomial",
                          predict.type = "prob",
                          fix.factors.prediction = TRUE)

# Creating 5 folds for random points separation to evaluate the model (GLM)
perf_level_spCV = makeResampleDesc(method = "SpRepCV", folds = 5, reps = 5)

# Validating the model in binomial regression
sp_cvBinomial = resample(learner = lrnBinomial, task = task,
                         resampling = perf_level_spCV,
                         measures = mlr::auc,
                         show.info = FALSE)

# Visualising the spatial fold partitioning on a map
plotsSP = createSpatialResamplingPlots(task, resample = sp_cvBinomial,
                                       crs = crs(allEnv), datum = crs(allEnv),
                                       color.test = "red", point.size = 1)
cowplot::plot_grid(plotlist = plotsSP[["Plots"]], ncol = 3, nrow = 2,
                   labels = plotsSP[["Labels"]])




# Fitting the GLM model

# Using poly() to consider non-linear species responses
glm.melesmeles = glm(Pres ~ poly(broadleaf, 3) + urban + elevation,
                     binomial(link = 'logit'),
                     data = all.cov)

# Printing model summary
summary(glm.melesmeles)

# Predicting across the environmental raster stack
prGLM = predict(allEnv, glm.melesmeles, type = "response")

# Plotting final GLM predicted distribution map
plot(prGLM, main = "GLM Predicted Distribution - Meles meles")
plot(melesmelesFin$geometry, add = T, pch = 16, cex = 0.7, col = "red")


# Displaying Broadleaf responses
glmNew = data.frame(broadleaf = seq(0, max(all.cov$broadleaf), length = 1000),
                    elevation = mean(all.cov$elev),
                    urban = mean(all.cov$urban))

# Predicting meles meles responce
preds = predict(glm.melesmeles, newdata = glmNew, type = "response", se.fit = TRUE)
glmNew$fit = preds$fit
glmNew$se = preds$se.fit

# Plotting responce accross broadleaf woodland
glm_p1 = ggplot(glmNew, aes(x = broadleaf, y = fit)) +
  geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "darkgreen", alpha = 0.3) +
  geom_line(colour = "darkgreen", linewidth = 1) +
  labs(title = "GLM-Broadleaf",
       x = "Proportion of Broadleaf Woodland (1800m radius)",
       y = "Predicted Occurrence Probability") +
  theme_bw()

print(glm_p1)

# Displaying Urban response
glmNewUrban = data.frame(urban = seq(0, max(all.cov$urban), length = 1000),
                         elevation = mean(all.cov$elev),
                         broadleaf = mean(all.cov$broadleaf))

predUrban = predict(glm.melesmeles, newdata = glmNewUrban, type = "response", se.fit = TRUE)
glmNewUrban$fit = predUrban$fit
glmNewUrban$se = predUrban$se.fit

# Plotting responce accross urban land cover
glm_p2 = ggplot(glmNewUrban, aes(x = urban, y = fit)) +
  geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "darkred", alpha = 0.3) +
  geom_line(colour = "darkred", linewidth = 1) +
  labs(title = "GLM-Urban",
       x = "Proportion of Urban Cover (2300m radius)",
       y = "Predicted Occurrence Probability") +
  theme_bw()

print(glm_p2)

# Displaying Elevation response
glmNewElev = data.frame(elevation = seq(0, max(all.cov$elev), length = 1000),
                        urban = mean(all.cov$urban),
                        broadleaf = mean(all.cov$broadleaf))

predElev = predict(glm.melesmeles, newdata = glmNewElev, type = "response", se.fit = TRUE)
glmNewElev$fit = predElev$fit
glmNewElev$se = predElev$se.fit

# Plotting responce against elevation
glm_p3 = ggplot(glmNewElev, aes(x = elevation, y = fit)) +
  geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "steelblue", alpha = 0.3) +
  geom_line(colour = "steelblue", linewidth = 1) +
  labs(title = "GLM-Elevation",
       x = "Elevation (metres)",
       y = "Predicted Occurrence Probability") +
  theme_bw()

print(glm_p3)



# Evaluating MaxEnt Model

library(dismo)   # loading dismo to perform kfold() and evaluate() functions
library(maxnet)  # loading maxnet for fitting models
library(glmnet)  # loading glmnet for robust regression model
library(precrec) # loading precrec for AUC calculation

# Separating presence and background data
Pres.cov = all.cov[all.cov$Pres == 1, ]
Back.cov = all.cov[all.cov$Pres == 0, ]


# Partitioning of data based on geographic location using a spatial grid for Spatial CV.
# Deving study area into six grids for Training and testing data in a spatially separated way.

area_grid = st_make_grid(melesmelesFin, c(50000, 50000),
                         what = "polygons", square = T)

# Converting grid to sf object and assign an ID to each square
area_grid_sf = st_as_sf(area_grid)
area_grid_sf$grid_id = 1:length(lengths(area_grid))


# Setting folds from the grids
folds = area_grid_sf$grid_id

# Converting the full dataset back to spatial object for intersection operations
dataPoints = st_as_sf(all.cov, coords = c("x", "y"))
st_crs(dataPoints) = crs(area_grid_sf)

maxEvalList = list()

for (i in folds) {
  # Selecting all grid squares except i as the training region
  gridTrain = subset(area_grid_sf, area_grid_sf$grid_id != i)
  
  # Spatially subsetting data points falling within the training region
  train = data.frame(st_drop_geometry(st_intersection(gridTrain, dataPoints)))
  
  # Selecting grid square i as the test region
  gridTest = subset(area_grid_sf, area_grid_sf$grid_id == i)
  
  # Spatially subsetting data points falling within the test region
  test = data.frame(st_drop_geometry(st_intersection(gridTest, dataPoints)))
  
  # Fitting MaxEnt model on spatially defined training data
  # lqph = linear, quadratic, product and hinge features for smooth responses
  maxnetMod = maxnet(train$Pres, train[1:3], classes = "lqph")
  
  # Predicting to spatially held-out test data
  pred = predict(maxnetMod, test, type = "cloglog")
  
  # Calculating AUC using precrec
  modauc = precrec::auc(precrec::evalmod(scores = pred,
                                         labels = test$Pres))
  # Storing AUC for this fold
  maxEvalList[[i]] = modauc$aucs[1]
  
}






# fitting MaxEnt model accorss the for best prediction
# Combining presence and background into one training dataset
dataTrain = rbind(Pres.cov, Back.cov)

# Fitting final MaxEnt model
finalMaxnet = maxnet(dataTrain$Pres, dataTrain[, 1:3], classes = "lqph")


# Prediction of species responce to identify suitable habitats

# Converting raster stack to data frame keeping track of which cells are NA
envDF = as.data.frame(allEnv, na.rm = FALSE)

# Identifying rows that have complete data (no NAs across all covariates)
validRows = complete.cases(envDF)

# Predicting for valid (non-NA) rows to avoid misalignment
predMaxnet = rep(NA, nrow(envDF))
predMaxnet[validRows] = predict(finalMaxnet,
                                envDF[validRows, ],
                                type = "cloglog")

# Applying predictions back into raster template
prMaxnet = allEnv[[1]]
values(prMaxnet) = predMaxnet

# Plotting the MaxEnt predicted distribution map
plot(prMaxnet, main = "MaxEnt Predicted Distribution - Meles meles")
plot(melesmelesFin$geometry, add = T, pch = 16, cex = 0.7, col = "red")


# Displaying MaxEnt model responces.

# Broadleaf response
maxNew = data.frame(
  broadleaf = seq(0, max(all.cov$broadleaf), length = 1000),
  elevation = mean(all.cov$elev),
  urban = mean(all.cov$urban)
)
maxNew$fit = predict(finalMaxnet, maxNew, type = "cloglog")

max_p1 = ggplot(maxNew, aes(x = broadleaf, y = fit)) +
  geom_line(colour = "darkgreen", linewidth = 1) +
  labs(title = "MaxEnt-Woodland",
       x = "Proportion of Broadleaf Woodland (1800m radius)",
       y = "Predicted Occurrence Probability (cloglog)") +
  theme_bw()

print(max_p1)

# Urban response
maxNewUrban = data.frame(
  urban = seq(0, max(all.cov$urban), length = 1000),
  elevation = mean(all.cov$elev),
  broadleaf = mean(all.cov$broadleaf)
)
maxNewUrban$fit = predict(finalMaxnet, maxNewUrban, type = "cloglog")

max_p2 = ggplot(maxNewUrban, aes(x = urban, y = fit)) +
  geom_line(colour = "darkred", linewidth = 1) +
  labs(title = "MaxEnt Response to Urban Land Cover",
       x = "Proportion of Urban Cover (2300m radius)",
       y = "Predicted Occurrence Probability (cloglog)") +
  theme_bw()

print(max_p2)

# Elevation response
maxNewElev = data.frame(
  elevation = seq(0, max(all.cov$elev), length = 1000),
  urban = mean(all.cov$urban),
  broadleaf = mean(all.cov$broadleaf)
)
maxNewElev$fit = predict(finalMaxnet, maxNewElev, type = "cloglog")

max_p3 = ggplot(maxNewElev, aes(x = elevation, y = fit)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  labs(title = "MaxEnt Response to Elevation",
       x = "Elevation (metres)",
       y = "Predicted Occurrence Probability (cloglog)") +
  theme_bw()

print(max_p3)


# Plot side by side predictions from  both GLM and MaxEnt models for comparison

par(mfrow = c(1, 2))

plot(prGLM, main = "GLM SDM - (A)")
plot(melesmelesFin$geometry, add = T, pch = 16, cex = 0.5, col = "red")

plot(prMaxnet, main = "MaxEnt SDM - (B)")
plot(melesmelesFin$geometry, add = T, pch = 16, cex = 0.5, col = "red")

par(mfrow = c(1, 1))

