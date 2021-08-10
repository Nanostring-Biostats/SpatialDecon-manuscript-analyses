###### analysis of grid over NSCLC
# study design: grid of ROIs over a single tumor, each split into Tumor/TME


rm(list = ls())
library(pheatmap)
library(scales)
library(viridis)
library(readxl)
library(dspNgs)
library(ComplexHeatmap)
library(SpatialDecon)
library(ggthemes)
library(logNormReg)
library(RColorBrewer)
library(dendextend)
source("spaceplot utils.R")


#### data loading ----------------------------------------
name = "ICP20th"
# load segment/AOI annotation:
annot = read.table(paste0("data/", name, "_SegmentProperties.txt"), 
                   header = T, sep = "\t", quote = "", stringsAsFactors = F, comment.char = "")
rownames(annot) = make.names(annot$Sample_ID)
head(annot)

# load xy coords and merge into annot:
xy = read.csv("data/ICP20th_Lung11_grid_SampleSheet_20200207 - xy and nuclei.csv")
rownames(xy) = make.names(xy$Sample_ID)
xy$roiid = c()
annot$x = xy[rownames(annot), "X.coordinate"]
annot$y = xy[rownames(annot), "Y.coordinate"]
annot$nuclei =  xy[rownames(annot), "Cell.number"]

# load gene/target annotations:
geneannot = read.table(paste0("data/", name, "_TargetProperties.txt"),
                       header = T, sep = "\t", quote = "", stringsAsFactors = F, comment.char = "")
rownames(geneannot) = geneannot$TargetName

# extract geneannot pathways:
allgenesets = unique(unlist(strsplit(geneannot$TargetGroup, ";")))
genesets = matrix(NA, nrow(geneannot), length(allgenesets), dimnames = list(rownames(geneannot), allgenesets))

# load raw counts:
raw = as.matrix(read.table(paste0("data/", name, "_TargetCountMatrix.txt"),
                           header = T, row.names = 1))
raw[1:5,1:5]

# load normalized data:
negnorm = as.matrix(read.table(paste0("data/", name, "_NegNorm_TargetCountMatrix.txt"),
                            header = T, row.names = 1))
negnorm[1:5,1:5]

# get background:
bg = derive_GeoMx_background(norm = negnorm,  
                             probepool = geneannot[rownames(negnorm), "Pooling"], 
                             negnames = rownames(geneannot)[grepl("Neg", rownames(geneannot))])
bg = bg[rownames(negnorm), match(colnames(negnorm), make.names(colnames(bg)))]
colnames(bg) = colnames(negnorm)

# get SNR (signal-to-noise ratio):
snr = negnorm / bg


## try new normalization: bg sub then Q3 scale:
# bg-subtracted snr:
snr.sub = pmax(snr - 1, 0)
hist(apply(snr.sub, 2, quantile, 0.75), breaks = 20)
# scaled bg subtracted snr:
snr.sub.scale = sweep(snr.sub, 2, apply(snr.sub, 2, quantile, 0.85, na.rm = T), "/")


#### Custom annotation handling and color definitions ------------------------

aoicols = c("chartreuse2", "darkcyan")  
names(aoicols) = c("Tumor", "TME")

annot$aoicol = aoicols[annot$AOI.name]

#### data alignment and subsetting ------------------------------

# evaluate alignment between sample notes and raw
print("AOIs present in expression data but missing from segment notes: ")
print(setdiff(colnames(raw), rownames(annot)))
print("AOIs present in segment notes but missing from expression data: ")
print(setdiff(rownames(annot), colnames(raw)))

# use only aois present in both
shared.aois = intersect(rownames(annot), colnames(raw))

# define a subset of AOIs for this study: 
# (in the below, we're just removing a glass AOI)
include.aois = colnames(raw)

# and keep only that subset of the data:
annot = annot[include.aois, ]
raw = raw[, include.aois]
#norm = norm[, include.aois]
negnorm = negnorm[, include.aois]
snr = snr[, include.aois]
snr.sub = snr.sub[, include.aois]
snr.sub.scale = snr.sub.scale[, include.aois]


# sanity check
stopifnot(identical(colnames(raw), rownames(annot)))


#### QC genes and AOIs: ------------------------

negs = rownames(raw)[grepl("Neg", rownames(raw))]

# flag genes that rise above 1.5 * background at least 2 times:
genes.ever.above.lod = which(rowSums(snr > 1.5) >= 2)

# exclude AOIs with outlier low signal:
signal.metric = apply(snr, 2, quantile, 0.85)
hist(signal.metric, breaks = 30, col = "grey50", xlim = c(0, max(signal.metric)))
signal.thresh = 3.25
abline(v = signal.thresh, col = 2, lwd = 2)
aois.with.good.signal = rownames(annot)[which(signal.metric > signal.thresh)]
raw = raw[, aois.with.good.signal]
annot = annot[aois.with.good.signal, ]
negnorm = negnorm[, aois.with.good.signal]
snr = snr[, aois.with.good.signal]
snr.sub = snr.sub[, aois.with.good.signal]
snr.sub.scale = snr.sub.scale[, aois.with.good.signal]
bg = bg[, aois.with.good.signal]

# sanity check
stopifnot(identical(colnames(raw), rownames(annot)))
stopifnot(identical(colnames(snr.sub.scale), colnames(bg)))


#### choose a normalization method ---------------------------------------
# signal strength:
pairs(annot[, c("NormFactorHK", "NormFactorQ3", "NormFactorNeg_01")], pch = 16, col = annot$aoicol)
# all 3 methods look good for both tumor and TME, but it'll be risky to compare tumor vs. TME directly.

# decision: go with neg normalization, in the theory that most genes are driven by background:
norm = snr
# expected background for norm matrix:
bg.norm = replace(norm, TRUE, 1)


#### split AOI-level data into ROI-level data: --------------------------------------

# (below code assumes there are <= 1 AOI of each type per ROI)
rannot = list(ROI = unique(annot$ROI))
aoitypes = unique(annot$AOI.name)
roinorm = roinorm.bgsub = list()
for (atype in aoitypes) {
  roinorm[[atype]] = roinorm.bgsub[[atype]] = c()
}

## build a data matrix for each AOI type:
for (roi in rannot$ROI) {
  ids = match(paste0(roi, aoitypes), paste0(annot$ROI, annot$AOI.name))
  names(ids) = aoitypes
  for (atype in aoitypes) {
    roinorm[[atype]] = cbind(roinorm[[atype]], norm[, ids[atype]])
    roinorm.bgsub[[atype]] = cbind(roinorm.bgsub[[atype]], pmax(norm - bg.norm, 0)[, ids[atype]])
  }
}
for (atype in names(roinorm)) {
  colnames(roinorm[[atype]]) = rannot$ROI
}

## record indices of norm in the roi-level datasets:
roiindices = list()
for (atype in names(roinorm)) {
  roiindices[[atype]] = rownames(annot)[match(paste0(rannot$ROI, atype), paste0(annot$ROI, annot$AOI.name))]
}

  
#### flag genes as having good enough data to analyze ---------------------------
# (though all genes get included in decon)
use.genes = (apply(norm, 1, sd, na.rm = T) > 0.2) & (apply(norm, 1, max) > 2.5)
table(use.genes)
use.gene.names = names(which(use.genes))

# genes to use specific to each AOI type (tumor and TME)
aoitype.use.genes = list()
for (atype in aoitypes) {
  tempind = which(annot$AOI.name == atype)
  aoitype.use.genes[[atype]] = names(which((apply(norm[, tempind, drop = F], 1, sd, na.rm = T) > 0.2) & (apply(norm[, tempind, drop = F], 1, max) > 2.5)))
}




#### Run decon --------------------------------------------

# load in fine-grained profile matrices:
load("data/zilionis.minor.RData")

# full version:

usegenes = intersect(intersect(rownames(norm), rownames(zilionis.minor)),                     rownames(safeTME))
X = cbind(zilionis.minor[usegenes, ], safeTME[usegenes, c("fibroblasts", "endothelial.cells")])

res = spatialdecon(norm = negnorm, 
                   raw = raw,
                   bg = bg,
                   X = X,
                   cell_counts = annot$nuclei,
                   is_pure_tumor = annot$AOI.name == 'Tumor',
                   maxit = 5000)
res$cell.counts = res$cell.counts$cell.counts
save(res, file = "decon results from zilionis.minor plus gentles.RData")

#### plot cells in space ------------------------------------------

# infer polygon:
bound = getBoundary(annot$x, annot$y, marg = 0.1)


# plots one cell at a time:
svg(paste0("cells in space zilionis.minor.svg"), height = 15, width = 15)
par(mfrow = c(7, 7))
for (cell in colnames(X)) {
  par(mar = c(0,0,0,0))
  plot(x = annot[roiindices$TME, "x"],
       y = annot[roiindices$TME, "y"],
       cex = res$cell.counts[cell, roiindices$TME] * .05 + 0.1,
       col = alpha("dodgerblue4", 0.7), pch = 16,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
       xlim = c(0, 9000), ylim = c(0, 9000), asp = 1)
  text(4500, 8400, cell, cex = 2)
  polygon(bound$x, bound$y, col = rgb(0,0,0,0.07), border = NA)
}
dev.off()

# stats:
round(rowSums(res$beta[paste0("hMac", 1:9), ]) / sum(rowSums(res$beta[paste0("hMac", 1:9), ])), 2)
round(rowSums(res$beta[paste0("hDC", 1:3), ]) / sum(rowSums(res$beta[paste0("hDC", 1:3), ])), 2)
round(rowSums(res$beta[paste0("hT", 1:7), ]) / sum(rowSums(res$beta[paste0("hT", 1:7), ])), 2)
