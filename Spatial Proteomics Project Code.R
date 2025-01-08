if (!require("MSnbase", quietly = TRUE))
  install.packages("MSnbase")
BiocManager::install("MSnbase")

if (!require("pRoloc", quietly = TRUE))
  install.packages("pRoloc")
BiocManager::install("pRoloc")

if (!require("pRolocdata", quietly = TRUE))
  install.packages("pRolocdata")
BiocManager::install("pRolocdata")

if (!require("pRolocGUI", quietly = TRUE))
  install.packages("pRolocGUI")
BiocManager::install("pRolocGUI")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rpx")

library("MSnbase")
library("pRolocdata")
library("pRoloc")
library("pRolocGUI")
library("rpx")

## get and view data
data(hyperLOPIT2015)
dim(hyperLOPIT2015)
table(fData(hyperLOPIT2015)$markers)

## view Subcellular Locations
hlm <- pRoloc::markerMSnSet(hyperLOPIT2015[, hyperLOPIT2015$Replicate == 1])
dim(hlm)
getMarkerClasses(hlm)

# Single: Spatial Map and Profile Plot
pRolocVis(object = hyperLOPIT2015, fcol = "markers")

## Compare: Spatial Maps and Profile Plots
data(hyperLOPIT2015ms3r1)
data(hyperLOPIT2015ms3r2)
spatial_proteomics_2 <- MSnSetList(list(hyperLOPIT2015ms3r1, hyperLOPIT2015ms3r2))
pRolocVis(spatial_proteomics_2, app = "compare", fcol = "markers")

data("hyperLOPITU2OS2018")
data("lopitdcU2OS2018")
xx <- MSnSetList(list(hyperLOPITU2OS2018, lopitdcU2OS2018))
if (interactive()) {
  pRolocVis(xx, app = "compare", fcol = c("markers", "final.assignment"))
}


## PCA Scree
plot2D(hyperLOPIT2015, method = "scree")


## Heatmap Locations vs Locations
hlq <- QSep(hyperLOPIT2015)
levelPlot(hlq)

fmat <- mrkConsProfiles(hyperLOPIT2015)
plotConsProfiles


## Heatmap Proteins vs Locations
mm <- getMarkerClasses(hyperLOPIT2015)
m_order <- levels(factor(mm))[order.dendrogram(hc)]


## Subcellular Locations Dendogram
hc <- mrkHClust(hyperLOPIT2015)

par(cex=0.55, mar=c(20, 10, 4, 2))
plot(hc, xlab="", ylab="", main="", sub="", axes=FALSE)
par(cex=1)
title(main="Subcellular Location Dendogram")
axis(2)


## Protein Abundance Heatmap
heatmap(MSnbase::exprs(hyperLOPIT2015))