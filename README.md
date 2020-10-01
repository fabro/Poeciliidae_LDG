# Poeciliidae Latitudinal Diversity Gradient
R code for the analyses from "Evolutionary and environmental drivers of viviparous freshwater fishes richness across the Americas" by García-Andrade et al. Submitted to Global Ecology and Biogeography

R code to replicate all data treatment and analyses from the paper

## Code is composed of six sections:
(1) Extents of occurrence: This section generates the species' extents of occurrence (EOO, range maps) using alpha convex hull polygons overlapped on sub-basins, and/or the sub-basin distribution of each species.

(2) Presence-absence matrix: This section creates the presence-absence matrix (PAM) that will be used to calculate all modeled variables (response and predictors),  using the species' extents of occurrence. 

(3) Evolutionary time: This section calculates the mean pairwise distance (MPD) across a multiphylo object

(4) Speciation rates: This section calculates the mean diversification rate (mDR)

(5) Environmental variables: This section prepares all environmental variables to the resolution of 50x50 km, Mollweide projection and the extent of study area

(6) Piecewise structural equation modeling:  This code section fits a series of spatial autoregressive models, using different  neighbour distances and weights scheme to find the best fit for each pathway. Also, it runs and fits a piecewise structural equation model (pSEM) for the data.


## References

## Contact
aberenicega@gmail.com, fabricio.villalobos@gmail.com
