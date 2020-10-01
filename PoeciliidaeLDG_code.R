### This code replicates all the analyses done in the paper entitled: 
# "Evolutionary and environmental drivers of viviparous freshwater fishes richness across the Americas"
# by Garcia-Andrade et al. Submitted to Global Ecology and Biogeography


#--------------------------------------------------------------#
####                1. EXTENT OF OCCURRENCE                #####
#--------------------------------------------------------------#

# Load packages 
library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(alphahull)

# Set your working directory
setwd("~/yourpath")

# Load function that convert the alphahull in a spatial object
# Download this function by Andrew Bevan from (https://stat.ethz.ch/pipermail/r-sig-geo/2012-March/014409.html)
source(file = "ah2sp.R")

### Read and prepare data ###
# Read the shapefile of hydrobasins level 8 layer (Lehner & Grill, 2013)
hidro <- readOGR(".", layer = "America_lev8")
# Add a unique code to each sub-basin to avoid duplicated names
hidro@data$unico <- c(1:60838)
# Read the clean species occurrence database
data_ocurr <- read.csv("poeciliidae_occ_june2019.csv", header = TRUE)
# Create spatial points with the occurrences
coord <- data.frame(coord_x = data_ocurr$long, coord_y= data_ocurr$lat)
occurrences <- SpatialPoints(coord)
# Add coordinate reference system (crs) to the spatial points object
wgs_84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(occurrences) <- wgs_84
# Add attributes to spatial points
df <- data.frame(genus = data_ocurr$Genus, species= data_ocurr$sp)
occurrences_df <- SpatialPointsDataFrame (occurrences, df)
# Re-project both occurrences and basin shapefile with the wgs 84 projection
occu_reproj <- spTransform(occurrences_df, wgs_84)
basin_reproj <- spTransform(hidro, wgs_84)

### Create the EOOs ###
# Note: all species must have at least three occurrences to allow creation of the alpha convex hull polygon
# Also, errors could happen if the occurrences are too geographically close

spp_a <- data.frame(table(df$species))
species <- subset(spp_a, spp_a$Freq >= 3)

# Create the polygons (EOO)
for(i in 1:dim(species)[1]){
# Select species
	sp <- species[i,1]
# Select species occurrences
	occ_sp_data <- data_ocurr[data_ocurr$sp==sp,]
	occ_sp_coord <- cbind(occ_sp_data$long, occ_sp_data$lat)
# Create the alpha convexhull with an alpha of six (as suggested by Pelayo-Villamil et al. 2015)
	ahull_sp <- ahull(occ_sp_coord, alpha = 6)
# Tranform the ahull object to spatial object
	ahull_polygon <- ah2sp(ahull_sp)
	slot(ahull_polygon, "polygons") <- lapply(slot(ahull_polygon, "polygons"),
	checkPolygonsHoles)
# Join all polygons
	ahull_polygon<- unionSpatialPolygons (ahull_polygon,
	IDs = as.character(ahull_polygon$HID))
# Project in WGS 84
	proj4string(ahull_polygon) <- wgs_84
	ahull_reproj <- spTransform(ahull_polygon, wgs_84)
	ahull_reproj <- gBuffer(ahull_polygon, byid = TRUE, width = 0)
# Overlap the ahull convex hull on sub-basins
	hydbas_map_sp <- gIntersects(ahull_reproj, basin_reproj, byid = TRUE)
	hydbas_map_sp_1 <- basin_reproj[as.vector(hydbas_map_sp), ]
	hydbas_map_sp_final <- gUnionCascaded(hydbas_map_sp_1)
# Add attributes to spatial object
	datfram <- data.frame(sciname = sp)
	hydbas_map_sp_final_df <- SpatialPolygonsDataFrame(hydbas_map_sp_final,
	datfram)
# Save EOO as a shapefile
	writeOGR(hydbas_map_sp_final_df, dsn = ".", layer = sp,
	driver="ESRI Shapefile")
}

### Calculate the sub-basin distribution ### 
# This procedure was used for species with less than three available occurrences

species_sub <- subset(spp_a, spp_a$Freq <3)

for (i in 1:dim(species_sub)[1]){
# Select species
	sp <- as.vector(species_sub[i,1])
# Select species occurrences
	occ_sp <- occu_reproj[occu_reproj$species == sp,]
# Select the basins where the occurrences are
	subb_spp <- gIntersects(occ_sp, basin_reproj, byid = TRUE)
	subb_spp2 <- ifelse(subb_spp == TRUE, 1L, 0L)
	subb_spp2_rows <- data.frame(presence= rowSums(subb_spp2))
	spp_prev1 <- cbind(basin_reproj@data$unico, subb_spp2_rows)
	spp_prev2 <- spp_prev1[,1:2]; colnames(spp_prev2) <- c("unico", "presence")
	spp_dist <- merge(basin_reproj, spp_prev2, by.y = "unico")
	spp_dist_final <-  spp_dist[which(spp_dist$presence >= 1),]
	spp_dist_final1 <- gUnionCascaded(spp_dist_final)
# Add attributes to spatial object
	datfram <- data.frame(sciname = sp)
	spp_dist_final_df <- SpatialPolygonsDataFrame(spp_dist_final1, datfram)
# Save sub-basin distribution
	writeOGR(spp_dist_final_df, dsn =".", layer = sp, driver="ESRI Shapefile")
}

#--------------------------------------------------------------#
####              2. PRESENCE-ABSENCE MATRIX               #####
#--------------------------------------------------------------#

# Load packages
library(raster)
library(rgdal)
library(maptools)
library(letsR)

# Set your working directory
setwd("~/yourpath") 

### Read and prepare data ###
# Read and join all the species' range
# Note: All shapefiles must be in the same folder
spp_shps <- list.files(path = ".", pattern = "*.shp$", full.names = TRUE)
all_shps <- lapply(spp_shps, readOGR)
poec_shps <- do.call(rbind, all_shps)
# Reproject to equal area Mollweide projection (km)
poe_shps_reproj <- rgdal::spTransform(poec_shps, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs")

### Create a PAM of grid cells resolution of 50X50 km ###
poe_pam <- letsR::lets.presab(poe_shps_reproj, xmn = -11625.17, xmx = -2406.773, ymn = -5321.016, ymx = 5864.784, resol = 50, 
remove.cells = TRUE, remove.sp = TRUE, show.matrix = FALSE, crs= "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs", count = TRUE)

### Save the PAM ###
save(poe_pam, file= "PAM.AHULL.MEDG_256spp.RData")

#--------------------------------------------------------------#
####                   3. EVOLUTIONARY TIME                #####
#--------------------------------------------------------------#

# Load packages
library(picante)
library(phytools)
library(raster)

# Set your working directory
setwd("~/yourpath")

### Read and prepare data (phylogeny and PAM) ###
# Read the expanded multiphylo from SUNPLIN (Martins et al. 2013)
trees <- read.nexus("node1000trees.nex", force.multi = TRUE)
# Load the presence absence matrix
load("PAM.AHULL.MEDG_256spp.RData")
# Separate the matrix from the PAM object
poe_pam <- as.data.frame(PAM.AHULL.MEDG$Presence_and_Absence_Matrix)
# Replace the space by undescore to match with the phylogeny tip names
colnames(poe_pam ) <- gsub(" ", "_", colnames(poe.pam))

### Calculate MPD across the multiphylo object ###
# Create an empty object for save mpd values
mpd_trees <- NULL
# Loop to calculate MPD across the multiphylo object
	for (i in 1:length(trees)){
  # Match PAM species & tree species
  	poe_data_spp <- match.phylo.comm(trees[[i]], poe_pam )
  # Clean community matrix & tree
  	poe_pam_spp <- poe_data_spp$comm
  	poe_tree_spp <- poe_data_spp$phy
  # Calculate cophenetic distances
  	poe_dist_spp <- cophenetic.phylo(poe.tree.spp)
  # Calculate MPD
  	poe_mpd_spp <- as.vector(mpd(poe_pam_spp, poe_dist_spp, abundance.weighted=FALSE))
  # Save mpd for each tree[i]
  	mpd_trees <- cbind(mpd_trees, poe_mpd_spp)
	}

# Average MPD for the 1000 expanded trees
mpd_mean <- as.data.frame(rowMeans(mpd_trees))

# Add coordinates to mean MPD values. Note that the value is divided by two to get the real age
mpd_ahull <- cbind(poe_pam [,1:2], (mpd.mean$`rowMeans(mpd_trees)`/2)); colnames(mpd_ahull) <- c("longitude", "latitude", "mpd")

### Calculate the median edge length for grid cells with species richness of one across the multiphylo object ###
# Code from Liam Revell (http://blog.phytool.org/2013/10/finding-edge-lengths-of-all-terminal.html)
edge_length <- NULL
for (i in 1:length(trees)){
	tree <- trees[[i]]
	n <- length(tree$tip.label)
	ee <- as.data.frame(setNames(tree$edge.length[sapply(1:n,function(x,y) which(y==x),y=tree$edge[,2])],tree$tip.label))
	edge.length <- cbind(ee[,1], edge.length)
}
# Compute median edge length for each species
median_elength <- apply(edge_length, 1, median, na.rm = TRUE)
# Join species name with its median edge length 
tips <- row.names(ee)
spp_elength<- cbind(tips, median_elength)

### Add the median edge length to grid cells with one species ###
# Filter grid cells with one species and get the species name
poe_pam_spp$sr_phylo <- rowSums(poe_pam_spp) # Species richness in grid cells
spp_grid <- vector()
for (i in 1:nrow(poe_pam_spp)){
	if (poe_pam_spp[i, 257] == 1) {
	spp_grid[i] <- colnames(poe_pam_spp[which(poe_pam_spp[i, 1:256] == 1)])
	}
	else{
	spp_grid[i] <- NA
	}
}
# Add a column with the species name of those species in grid cells with richness of one
poe_pam_spp$spp_grid <- spp_grid

### Assign the median edge value in the MPD column ###
for (i in 1:nrow(poe_pam_spp)){
# If SR values is equal to 1 replace the row value for the especies edge length
	if (poe_pam_spp[i, 257] == 1) {
	sp <- poe_pam_spp[i, 258]
	sp_data <- subset(spp_elength, tips %in% sp)
	mpd_ahull[i,3] <- as.numeric(sp_data[,2])
	}
# If the value is different of 1 pass to the next row
	else{
	next;
	}
}

### Prepare MPD results ###
# MPD values were multiple by 100 because the edge length in the trees was in decimal format
mpd_ahull$mpd <- ((mpd_ahull$mpd)*100)
# Delete NAs
mpd_ahull <- na.omit(mpd_ahull)
# Razterize MPD data
mpd_ras <- rasterize(mpd_ahull[,c(1,2)], PAM.AHULL.MEDG$Richness_Raster, mpd_ahull[,3])
# Save MPD raster as tif
writeRaster(mpd_ras, filename = "mpd_ahull05_1000sun_MLprj_256spp.tif", format = "GTiff")

#--------------------------------------------------------------#
####                  4. SPECIATION RATES                  #####
#--------------------------------------------------------------#

# Load packages
library(picante)
library(psych)

# Set your working directory
setwd("~/yourpath")

### Load and prepare data (phylogeny and PAM) ###
# Load expanded trees from SUNPLIN
trees <- read.nexus("node1000trees.nex", force.multi = TRUE)
# Load the presence absence object
load("PAM.AHULL.MEDG_256spp.RData")
# Separate the presence absence
poe_pam <- as.data.frame(PAM.AHULL.MEDG$Presence_and_Absence_Matrix)
# Replace the space by undescore to match with the tip names of the tree
colnames(poe_pam) <- gsub(" ", "_", colnames(poe_pam))
# Match species from the pam with the species in the tree
poe_data_spp <- match.phylo.comm(trees[[1]], poe_pam)
poe_pam_spp <- poe_data_spp$comm

### Function that calculateS the mDR for a multiphylo object, mDR sensu Jetz et al. (2012) ###
# trees = must be a multiphylo object
# pam.spp = a presence absence matrix that only contains the species in the phylo
mDR <- function(trees, pam_spp){
# Empty object to save the mDR results
	tip_dr_trees <- NULL
# Loop to calculate the mean diversification rate (mDR)
	for (i in 1:length(trees)){
		print(i)
# Calculate the mean equal splits based in Redding & Mooers 2005
		tree <- trees[[i]]
		tree$edge.length <- tree$edge.length*100
# Note the multiplication by 100 because the phylogeny of Reznick et al. (2017) is in decimal units
		equal_splits <- evol.distinct(tree, 
                                  type = "equal.splits",
                                  scale = FALSE,
                                  use.branch.lengths = TRUE)
# Calculate the tipDR for each species across the 1000 trees
		tip_dr <- 1/equal_splits$w # inverse of equal splits
		tip_dr_trees <- cbind(tip_dr_trees, tip_dr)
					}
# Calculate the harmonic mean of tip DR estimate for each species across the 1000 trees
		dr_spp <- data.frame(species = equal_splits$Species)
		dr_spp$tip_dr<- apply(tip_dr_trees, 1, harmonic.mean) 

### Calculate the mDR in each grid cell
		mean_dr <- NULL 
		for (j in 1:nrow(poe_pam_spp)){
# Get the list of species in the grid cell
			spp_pixel <- colnames(poe_pam_spp[which(poe_pam_spp[j, ] == 1)])
# Get the mean tip DR for all the species in the grid cell
			sp_data <- subset(dr_spp, species %in% spp_pixel)
# Calculate the harmonic mean
			mean_dr[j] <- psych::harmonic.mean(sp_data$tip_dr)
										}
		return(mean.dr)
			}

### Calculate mDR 
mDR_ahull <- mDR(trees = trees, pam_spp = poe_pam_spp)
# Save mDR
mean_dr_final <- data.frame(longitude = poe_pam[,1], latitude = poe_pam[,2], mDR = mDR_ahull)
write.csv(mean_dr_final, file = "mDR_ahull05_1000trees.csv", row.names = F)

#--------------------------------------------------------------#
####             5. ENVIRONMENTAL VARIABLES                #####
#--------------------------------------------------------------#

#Load packages
library(raster)
library(rgdal)
library(plotKML)

# Set your working directory
setwd("~/yourpath")

### Prepare the current temperature and precipitation seasonality ###
# Load bioclimatic variables at 10-minutes from WorldClim v1.4, including the 
bios <- stack(list.files(path = "~/bio_10m_bil/", pattern= '.bil$', full.names = T))
# Reprojection to equal-area grid based on the Mollweide projection
ml_km <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
bios_ml <- plotKML::reproject(bios, ml_km)
# Mask raster: richness raster
load("PAM.AHULL.MEDG_256spp.RData") # Prescence-absence matrix
sr_ras <- PAM.AHULL.MEDG$Richness_Raster
ext_mask <- extent(sr_ras)
# Crop climatic variables
bios_crop <- crop(bios_ml, ext_mask, snap="near")
# Resample with "bilinear interpolation" to the resolution of the richness raster
bios_crop_resam <- resample(bios_crop, sr_ras, method = "bilinear")

# Select temperature and precipitation seasonality
TempSeas <- bios_crop_resam$bio4
PrepSeas <- bios_crop_resam$bio15

### Create the annual temperature and precipitation anomaly ###
# Select current annual temperature and precipitation
bio01_current <- bios_crop_resam$bio1
bio12_current <- bios_crop_resam$bio12

# Read paleoclimatic variables for the LGM by the CCSM4 model from WorldClim v1.4 at 10-minute resolution
ccsm4 <- stack(list.files(path = "~/cclgmbi_10m/", pattern= '.tif$', full.names = T))
# Reprojection to equal-area grid based on the Mollweide projection 
ccsm4_ml <- plotKML::reproject(ccsm4, ml_km)
# Crop climatic variables
ccsm4_crop <- crop(ccsm4_ml, ext_mask, snap="near")
# Resample with "bilinear interpolation" to the resolution of the richness raster
ccsm4_crop_resam <- resample(ccsm4_crop, sr_ras, method = "bilinear")
# Select annual temperature and precipitation for LGM by CCSM4 model
bio01_ccsm4 <- ccsm4_crop_resam$cclgmbi1
bio12_ccsm4 <- ccsm4_crop_resam$cclgmbi12

# Load paleoclimatic variables for the LGM by MIROC model from WorlClim v1.4 at 10-minute resolution
miroc <- stack(list.files(path = "~/mrlgmbi_10m/", pattern= '.tif$', full.names = T))
# Reprojection to equal-area grid based on the Mollweide projection 
miroc_ml <- plotKML::reproject(miroc, ml_km)
# Crop climatic variables
miroc_crop <- crop(miroc_ml, ext_mask, snap="near")
# Resample with "bilinear interpolation" to the resolution of the richness raster
miroc_crop_resam <- resample(miroc_crop, sr_ras, method = "bilinear")

# Select annual temperature and precipitation for LGM by MIROC model
bio01_miroc <- miroc_crop_resam$mrlgmbi1 # BIO 01
bio12_miroc <- miroc_crop_resam$mrlgmbi12 #BIO 12

# Create annual temperature anomaly
bios01_LGM <- stack(bio01_ccsm4, bio01_miroc)
bios01_LGM_avg <- sum(bios01_LGM)/2
TempAnom <- bio01_current-bios01_LGM_avg; names(TempAnom) <- "TempAnom"

# Create annual precipitation anomaly
bios12_LGM <- stack(bio12_ccsm4, bio12_miroc)
bios12_LGM_avg <- sum(bios12_LGM)/2
PrepAnom <- bio12_current-bios12_LGM_avg; names(PrepAnom) <- "PrepAnom"

### Drainage basin area ###
# Load and prepare the 30-sec HydroSHEDS layer (Lehner, Verdin, & Jarvis, 2006)
basins <- readOGR(".", layer = "hydrobasin")
# Reproject to equal-area grid based on the Mollweide projection 
basins_ml <- spTransform(basins, ml_km)
# Extract area values
basins_ras <- rasterize(basins_ml, sr_ras, field="area")

### Topography heterogeneity ###
# Load the topographical heterogeneity index (TH8; www.ipez.es/ModestR/) 
TH8 <- raster("Topographic heterogeneity8.ASC")
# Reprojection to equal-area grid based on the Mollweide projection 
crs(TH8) <- wgs84
TH8_ml <- plotKML::reproject(TH8, ml_km)
# Crop th8 variable
TH8_crop <- crop(TH8_ml, ext_mask, snap= "near")
# Resample with "bilinear interpolation" to the resolution of the richness raster
TH8_crop_resam <- resample(TH8_crop, sr_ras, method="bilinear")

### Net Primary Productivity (MOD17A3 v55 layer from NASA) for 2000 to 2014 ###
# Load NPP raster
npp <- stack(list.files(path = "~/NPP/", pattern= '.tif$', full.names = TRUE))
# Reprojection to equal-area grid based on the Mollweide projection 
crs(npp) <- wgs84
npp_ml <- plotKML::reproject(npp, ml_km)
# Crop NPP layers
npp_crop <- crop(npp_ml, ext_mask, snap="near")
npp_crop_resam <- resample(npp_crop, sr_ras, method="bilinear")
# Calculate NPP mean for the year 2000 to 2004
npp_avg <- mean(npp_crop_resam); names(npp_avg) <- "npp_avg"

### Create a data frama with all predictor variables 
# Load the mpd raster 
EvolTime <- raster("mpd_ahull05_1000sun_MLprj_256spp.tif")

# Read mDR values 
SpecRate <- read.csv("mDR_ahull05_1000trees.csv", header=TRUE)
colnames(SpecRate) <- c("x","y", "SpecRate")
SpecRate <- na.omit(SpecRate)
mDR <- rasterize(SpecRate[,c(1,2)], sr_ras, SpecRate[,3])

# Stack all variables: response variable & predictors variables (environmental and evolutionary)
vars_pred <- stack(sr_ras,TempSeas,PrepSeas,TempAnom,PrepAnom,TH8_crop_resam,basins_ras,npp_avg, EvolTime, mDR); names(vars_pred) <- c("SR","TempSeas", "PrepSeas", "TempAnom", "PrepAnom", "TH8", "BasinArea", "NPP", "EvolTime", "SpecRate")

# Create a data frame with all variables
df_vars <- as.data.frame(rasterToPoints(vars_pred))
# Clean data
df_vars_comp <- subset(df_vars, SR >0)
# Save data frame
write.csv(df_vars_comp, file = "vars_mollprj_50km.csv", row.names=FALSE)

# Save raster variables 
writeRaster(vars_pred, filename=names(vars_pred), bylayer=TRUE, format="GTiff")


#----------------------------------------------------------------------------------------------#
#### 6. Piecewise Structural Equation modelling using Spatial Autoregressive Error models  #####
#----------------------------------------------------------------------------------------------#
# This code is based on Skeels et al. 2020 Global Ecology and Biogeography ()

# Load packages
library(raster)
library(rgdal)
library(ncf)
library(spdep)
library(spatialreg) # version 1.1-3
library(dplyr)
library(piecewiseSEM) # version 2.0.2

# Set your working directory
dir_work <- ("~/yourpath") 
setwd(dir_work)

### Preparing data ###
# Read data (data frame that contain response and predictor variables)
vars <- read.csv("vars_mollprj_50km_256spp_final.csv", header = TRUE);
colnames(vars) <- c("x","y","BasinArea","EvolTime","NPP", "PrepAnom","PrepSeas",
					"SpecRate","SR", "TempAnom", "TempSeas", "TH8")
# Log transformation of the response variable (species richness)
vars$SR <- log(vars$SR+1)
# Scale all predictors variables to mean zero and standard deviation of 1 
# to allow coeficient comparison
vars$BasinArea <- scale(vars$BasinArea, center = T)
vars$PrepSeas <- scale(vars$PrepSeas, center = T)
vars$TempSeas <- scale(vars$TempSeas, center = T)
vars$NPP <- scale(vars$NPP, center = T)
vars$TempAnom <- scale(vars$TempAnom, center = T)
vars$PrepAnom <- scale(vars$PrepAnom, center = T)
vars$TH8 <- scale(vars$TH8, center = T)
vars$EvolTime <- scale(vars$EvolTime, center = T)
vars$SpecRate <- scale(vars$SpecRate, center = T)

### Model selection for SAR neighbourhood distances and weighting schemes ###
# This code fits a series of spatial autoregressive models,
# using different neighbour distances and weights scheme to find the best fit for each path 

# Matrix weighting schemes to prove
scheme <- c("W", "C", "S")

# Make a matrix of spatial coordinates (X and Y coordinates)
sp <- SpatialPoints(data.frame(x=vars$x, y=vars$y))
crs(sp) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
coords <- coordinates(sp)

# find min and maximum nearest neighbour
k1 = knn2nb(knearneigh(coords, k=1)) 
dist <- unlist(nbdists(k1, coords))
max1 <- max(dist)
min1 <- min(dist)

# create a series of neighbour matrices based on different distances (distances are in km)
d1 <- dnearneigh(sp, longlat = T, d1=0, d2=min1)
d2 <- dnearneigh(sp, longlat = T, d1=0, d2=max1) 

# Create a data frame with path number and equations that will be tested in the pSEM
path_equation <- data.frame(path_number = paste0("path",1:6))
path_equation$equation <-c(
"SR ~ BasinArea + TempSeas + PrepSeas + NPP + TempAnom + PrepAnom + TH8 + EvolTime + SpecRate",
"SpecRate ~ BasinArea + NPP + TempAnom + PrepAnom + TH8 + EvolTime",
"EvolTime ~ BasinArea + TempAnom + PrepAnom + NPP + TH8",
"NPP ~ TempSeas + PrepSeas", "TempSeas ~ TH8 + TempAnom",
"TempSeas ~ TH8 + TempAnom",
"PrepSeas ~ TH8 + PrepAnom")

### Create models for each path and weighting scheme ###
# We tested three schemes ("W", "C", "S"), two distance (d1, d2) and two models (OLS and SAR)
# To find the best combination that minimize the autocorrelation and increase the fit for each path equation

for (i in 1:length(scheme)){
# Create a directory by scheme
	setwd(dir_work)
	dir_scheme <- scheme[i]
	if(!dir.exists(dir_scheme)) dir.create(dir_scheme)
	setwd(dir_scheme)
# Create weighting scheme distance matrix
	spatial_weights_d1 <- nb2listw(d1, zero.policy=TRUE, style=scheme[i])
	spatial_weights_d2 <- nb2listw(d2, zero.policy=TRUE, style=scheme[i])
	for (p in 1:nrow(path_equation)) {
		schemes <- rep(c(paste(scheme[i])), times=3)
		paths <- rep(c(paste(path_equation[p,1])), times=3)
		models <- c("lm_mod", "error_d1", "error_d2")
		results_path <- data.frame(scheme=schemes, path=paths, model=models)
# OLS model 
		lm_mod <- lm(paste(path_equation[p,2]),
                 data = vars)
		lm_mod_s <- summary(lm_mod)
		R2_lm <- lm_mod_s[["adj.r.squared"]]
# SAR error models 
# Minimum distance
		error_d1 <- spatialreg::errorsarlm(paste(path_equation[p,2]),
		data = vars,listw = spatial_weights_d1,tol=1e-12,zero.policy=T)
		error_d1_s <-summary(error_d1, Nagelkerke=TRUE)
		R2_error_d1 <- error_d1_s$NK
# Maximum distance
		error_d2 <- spatialreg::errorsarlm(paste(path_equation[p,2]),
		data = vars,listw = spatial_weights_d2, tol=1e-12,zero.policy=T)
		error_d2_s <-summary(error_d2, Nagelkerke=TRUE)
		R2_error_d2 <- error_d2_s$NK
# Save R2 (pseudo R2 for errorsar) and AIC
		results_path$R2_models <- c(R2_lm, R2_error_d1, R2_error_d2)
		results_path$AIC <- AIC(lm_mod, error_d1, error_d2)
		write.csv(results_path, file=paste(path_equation[p,1], ".csv"),
		row.names=F)
# Make correlograms of residual autocorrelation
		cor.ols1.res<-correlog(vars$x, vars$y, z=residuals(lm_mod), na.rm=TRUE, increment=1, resamp=1)
		cor.sar1.res<-correlog(vars$x, vars$y, z=residuals(error_d1),na.rm=TRUE, increment=1, resamp=1)
		cor.sar2.res<-correlog(vars$x, vars$y, z=residuals(error_d2), na.rm=TRUE, increment=1, resamp=1)
# Save correlograms in a jpeg file
		jpeg(filename = paste(path_equation[p,1],".jpg"),
		width = 215, height = 279, units = "mm",res = 600)
		par(mfrow = c(3, 1))
		plot(cor.ols1.res, xlab = "Distance (Km)", ylab = "Moran's I", ylim=c(-1,1), type = "l",
		lwd= 2, main = paste(models[1]), cex.main=2, cex.lab=1.8, cex.axis=1.5)
		abline(h=0, lty=5)
		plot(cor.sar1.res, xlab = "Distance (Km)", ylab = "Moran's I", ylim=c(-1,1), type = "l",
		lwd= 2, main =paste(models[2]), cex.main=2, cex.lab=1.8, cex.axis=1.5)
		abline(h=0, lty=5)
		plot(cor.sar2.res, xlab = "Distance (Km)", ylab = "Moran's I", ylim=c(-1,1), type = "l",
		lwd= 2, main = paste(models[3]), cex.main=2, cex.lab=1.8, cex.axis=1.5)
		abline(h=0, lty=5)
		dev.off()
# Save path results
		save.image(file=paste(path_equation[p,1], ".RData"))
				}
}

### Summarize the results in a data frame easily accessable ###
results <- data.frame()
for (i in 1:length(scheme)){
	setwd(dir_work)
# Extract results from each scheme directory
	dir_scheme <- scheme[i]
	setwd(dir_scheme)
	scheme_results_list <- list.files(path = ".", pattern= '.csv$')
# Join all results
	scheme_results <- lapply(scheme_results_list, read.csv)
	scheme_results2 <- do.call("rbind", scheme_results)
	results <- rbind(results, scheme_results2)
				}
# Back to your working directory
setwd(dir_work)

# Save the model selection results in a csv 
write.csv(results, file="results_model_selection.csv", row.names=F)

# Check results and select the better fitted scheme and model
# by the lower AIC and the higher pseudo-R for each path

### Piecewise Structural Equation modelling evaluation ###
# After the model selection we choose the "W" spatial weighting matrix with d2 distance (maximum distance) for all paths 

# Create spatial weighting matrix
spatial_weights_d2 <- nb2listw(d2, zero.policy=TRUE, style="W")

# Testing the theoretical pSEM model
sem_sar_model_teorico <- piecewiseSEM::psem(
# 1 # Equation 1: species richness as response
spatialreg::errorsarlm(SR ~ BasinArea + TempSeas + PrepSeas + NPP + TempAnom + PrepAnom + TH8 + EvolTime + SpecRate,
                           data = vars,
                           listw = spatial_weights_d2,
                           tol=1e-12,zero.policy=TRUE),
# 2 # Equation 2: speciation rate as response
spatialreg::errorsarlm(SpecRate ~ BasinArea + NPP + TempAnom + PrepAnom + TH8 + EvolTime,
                           data = vars,
                           listw = spatial_weights_d2,
                           tol=1e-12,zero.policy=TRUE),
# 3 # Equation 3: evolutionary time as response  
spatialreg::errorsarlm(EvolTime ~ BasinArea + TempAnom + PrepAnom + TH8 + NPP,
                           data = vars,
                           listw = spatial_weights_d2,
                           tol=1e-12,zero.policy=TRUE),
# 4 # Equation 4: net primary productivity as response  
spatialreg::errorsarlm(NPP ~ TempSeas + PrepSeas + TH8,
                           data = vars,
                           listw = spatial_weights_d2,
                           tol=1e-12,zero.policy=TRUE),
# 5 # Equation 5: temperature seasonality as response  
spatialreg::errorsarlm(TempSeas ~ TH8 + TempAnom,
                           data = vars,
                           listw = spatial_weights_d2,
                           tol=1e-12,zero.policy=TRUE),
# 6 # Equation 6: precipitacion seasonality as response  
spatialreg::errorsarlm(PrepSeas ~ TH8 + PrepAnom,
                           data = vars,
                           listw = spatial_weights_d2,
                           tol=1e-12,zero.policy=TRUE), 
data=vars)

# Run summary model 
summary_sar_model_teorico <- summary(sem_sar_model_teorico)

# Save RData with the pSEM model and summary model
save(sem_sar_model_teorico, summary_sar_model_teorico, file="psem_sar_teorico_mollprj_50km_ahull.RData")

# Extract path coefficients
coefs_sar_model_teorico <- coefs(sem_sar_model_teorico)

# Save coefficients
write.csv(coefs_sar_model_teorico, file = "coefs_psem_sar_modelo_teorico_mollprj_50km.csv")

### Running final pSEM 
# Delete non-significant paths and add missing ones (check the summary model)

sem_sar_model_final <- piecewiseSEM::psem(
# 1 # Equation 1: species richness as response of all predictor variables
spatialreg::errorsarlm(SR ~ BasinArea + TempSeas + NPP + TempAnom + PrepAnom + EvolTime + SpecRate,
                         data = vars,
                         listw = spatial_weights_d2,
                         tol=1e-12,zero.policy=TRUE),
# 2 # Equation 2: speciation rate as response of evolutionary time
spatialreg::errorsarlm(SpecRate ~ EvolTime,
                         data = vars,
                         listw = spatial_weights_d2,
                         tol=1e-12,zero.policy=TRUE),
# 3 # Equation 3: evolutionary time as response of basin area + productivity + temperature seasonality + precipitation seasonality
spatialreg::errorsarlm(EvolTime ~ BasinArea + NPP + TempSeas + PrepSeas,
                         data = vars,
                         listw = spatial_weights_d2,
                         tol=1e-12,zero.policy=TRUE),
# 4 # Equation 4: net primary productivity as response of precipitacion seasonality + precipitacion anomaly
spatialreg::errorsarlm(NPP ~ TempSeas + PrepSeas + PrepAnom,
                         data = vars,
                         listw = spatial_weights_d2,
                         tol=1e-12,zero.policy=TRUE),
# 5 # Equation 5: temperature seasonality as response of temperature anomaly 
spatialreg::errorsarlm(TempSeas ~ TempAnom,
                         data = vars,
                         listw = spatial_weights_d2,
                         tol=1e-12,zero.policy=TRUE),
# 6 # Equation 6: precipitacion seasonality as response of precipitacion anomaly
spatialreg::errorsarlm(PrepSeas ~ TH8 + PrepAnom,
                         data = vars,
                         listw = spatial_weights_d2,
                         tol=1e-12,zero.policy=TRUE), 
# 7 # Correlate errors
PrepSeas %~~% TempSeas,
data=vars)

# Run summary model 
summary_sar_model_final <- summary(sem_sar_model_final)

# Save pSEM model and summary model
save(sem_sar_model_final, summary_sar_model_final, file = "psem_sar_final_mollprj_50km_ahull.RData")

# Extract path coefficients
coefs_sar_model_final<- coefs(sem_sar_model_final)

# Save coefficients
write.csv(coefs_sar_model_final, file = "coefs_psem_sar_modelo_final_mollprj_50km.csv")