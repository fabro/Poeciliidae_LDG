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
Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. Ø. (2012). The global diversity of birds in space and time. Nature, 491, 444–448.

Lehner, B., & Grill, G. (2013). Global river hydrography and network routing: Baseline data and new approaches to study the world’s large river systems. Hydrological Processes, 27, 2171–2186.

Lehner, B., Verdin, K., & Jarvis, A. (2006). HydroSHEDS technical documentation. Washington DC.

Pelayo-Villamil, P., Guisande, C., Vari, R. P., Manjarrés-Hernández, A., García-Roselló, E., González-Dacosta, J., … Lobo, J. M. (2015). Global diversity patterns of freshwater fishes - potential victims of their own success. Diversity and Distributions, 21, 345–356.

Reznick, D. N., Furness, A. I., Meredith, R. W., & Springer, M. S. (2017). The origin and biogeographic diversification of fishes in the family Poeciliidae. PLoS ONE, 12, e0172546.

Skeels, A., Esquerré, D., & Cardillo, M. (2020). Alternative pathways to diversity across ecologically distinct lizard radiations. Global Ecology and Biogeography, 29, 454–469.

## Contact
aberenicega@gmail.com, fabricio.villalobos@gmail.com
