###### analysis of grid over NSCLC
# study design: grid of ROIs over a single tumor, each split into Tumor/TME


rm(list = ls())
library(pheatmap)
library(scales)
library(viridis)
library(readxl)
library(ComplexHeatmap)
library(SpatialDecon)
library(ggthemes)
library(logNormReg)
library(RColorBrewer)
library(dendextend)
library(ggplot2)
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

#norm = sweep(raw, 2, annot$NormFactorHK, "/")
negnorm[1:5,1:5]

# get background:
#geneannot$Pooling2 = c("")
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
# all 3 methods look good for both tumor and TME, though risky to compare them directly.

# decision: go with neg normalization, in the theory that most genes are driven by background:
norm = snr
# expected background for norm matrix:
bg.norm = replace(norm, TRUE, 1)


#### split AOI-level data into ROI-level data: --------------------------------------

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

# full version:
res = spatialdecon(norm = negnorm, 
                   raw = raw,
                   bg = bg,
                   cell_counts = annot$nuclei,
                   is_pure_tumor = annot$AOI.name == 'Tumor',
                   maxit = 500)

res.collapsed = collapseCellTypes(fit = res, 
                                  matching = safeTME.matches)
res$cell.counts = res$cell.counts$cell.counts


# and without using tumor-specific information, to show the impact of tumor modelling
res.notumor = spatialdecon(norm = negnorm, 
                           raw = raw,
                           bg = bg,
                           cell_counts = annot$nuclei,
                           maxit = 500)

immcells = rownames(res$beta)[!grepl("tumor", rownames(res$beta))]


#### results in text: how many immune cells in tumor vs. TME?

mean(colSums(res$cell.counts[, annot$AOI.name == "Tumor"]))
mean(colSums(res$cell.counts[, annot$AOI.name == "TME"]))

sum(res.collapsed$cell.counts$cell.counts[c("T.CD8.naive", "T.CD8.memory"), ])
sum(res.collapsed$cell.counts$cell.counts[c("T.CD4.naive", "T.CD4.memory", "Treg"), ])
sum(res.collapsed$cell.counts$cell.counts[c("macrophages"), ])



#### supp figure: compare results with/without tumor:


# convert to data frame for ggplot:
betay = betan = tempcell = tempaoi = c()
for (cell in immcells) {
  betay = c(betay, res$beta[cell, ])
  betan = c(betan, res.notumor$beta[cell, ])
  tempcell = c(tempcell, rep(cell, ncol(res$beta)))
  tempaoi = c(tempaoi, annot$AOI.name)
}
plotdf = data.frame(y = betay, n = betan, cell = tempcell, aoi = tempaoi)
g = ggplot(plotdf, aes(x = n, y = y, col = aoi)) + 
  geom_point(aes(alpha = I(0.7))) +
  facet_wrap(~cell, nrow = 5) +
  theme_few() + 
  scale_color_manual(values = aoicols) +
  labs(x = "Estimate in model without tumor", y = "Estimate from model accounting for tumor") + 
  theme(axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        strip.text.x = element_text(size = 10), 
        legend.title = element_blank())

svg("results/supp figure - decon results with vs. without tumor.svg", height = 10, width = 7.2)
print(g)
dev.off()



#### gini coefficients for each cell type --------------------------------

gini = function(x) {
  x1 = matrix(rep(x, each = length(x)), nrow = length(x))
  x2 = t(x1)
  sum(abs(x1 - x2)) / (2 * mean(x) * length(x)^2)
}

# plot gini coefs vs. means:
plot(apply(res$cell.counts[, annot$AOI.name == "TME"], 1, mean),
     apply(pmax(res$cell.counts[, annot$AOI.name == "TME"], 1), 1, gini))
text(apply(res$cell.counts[, annot$AOI.name == "TME"], 1, mean),
     apply(pmax(res$cell.counts[, annot$AOI.name == "TME"], 1), 1, gini),
     rownames(res$cell.counts))

# barplot of gini coefs:
ginis = apply(pmax(res$cell.counts[, annot$AOI.name == "TME"], 1), 1, gini)

print(round(ginis, 2))


#### plot cells in space ------------------------------------------


# infer polygon:
bound = getBoundary(annot$x, annot$y, marg = 0.1)


# plots one cell at a time:
for (cell in immcells) {
  svg(paste0("results/cells in space granular - ", cell, ".svg"), width = 2, height = 2)
  par(mar = c(0,0,0,0))
  plot(x = annot[roiindices$TME, "x"],
       y = annot[roiindices$TME, "y"],
       cex = res$cell.counts[cell, roiindices$TME] / max(res$cell.counts[cell, ]) * 6 + 0.1,
       col = alpha(cellcols[cell], 0.7), pch = 16,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
       xlim = c(0, 9000), ylim = c(0, 9000))
  text(4500, 8400, cell, cex = 0.7)
  polygon(bound$x, bound$y, col = rgb(0,0,0,0.07), border = NA)
  dev.off()
}


### plot cell scores ~ space
svg("results//florets - cells in space - legend.svg", width = 4.75, height = 4)
par(mar = c(0,0,0,0))
# legend:
plot(0, 0, col = 0, xlim = c(-2, 2), ylim = c(-1.75, 1.75), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     bty = "n")
b = rep(1, length(immcells))
names(b) = immcells
angles = seq(0, 2 * pi, length.out = length(b) + 1)
for (j in 1:length(b)) {
  tempangles = seq(angles[j] + 0.05, angles[j+1] - 0.05, length.out = 20)
  xt = b[j] * cos(tempangles)
  yt = b[j] *  sin(tempangles)
  polygon(c(0, xt), c(0, yt), col = cellcols[names(b)[j]], border = "white")
  if (names(b)[j] == "monocytes.NC.I") {
    yt = yt - 0.08
  }
  if (names(b)[j] == "plasma") {
    yt = yt + 0.08
  }
  text(median(xt) * 1.4, median(yt) * 1.2, names(b)[j], cex = 1.1)
}
dev.off()

svg("results//florets - cells in space.svg")
par(mar = c(0,0,0,0))
florets(x = annot[roiindices$TME, "x"],
        y = annot[roiindices$TME, "y"],
        b = res$cell.counts[immcells, roiindices$TME] * 450,
        legendwindow = F,
        xaxt = "n", yaxt = "n", xlab = "", ylab = "",
        #border = "grey70",
        border = NA,
        xlim = c(0, 9000), ylim = c(0, 9000), bty = "n",
        cex = 1)
polygon(bound$x, bound$y, col = rgb(0,0,0,0.07), border = NA)
dev.off()


#### heatmap of cell abundances ------------------------------------


# cell score heatmap in TME only:
use = annot$AOI.name == "TME"
p1 = pheatmap(pmin(res$cell.counts[, use], 100),
              col = viridis_pal(option = "B")(100),
              show_colnames = F)
o = p1$tree_col$order

# plot clusters in space:
set.seed(0)
clusts = cutree(p1$tree_col, k = 7)
# use a heatmap to name the clusters:
pheatmap(pmin(res$cell.counts[p1$tree_row$order, use][, o], 100),
         col = viridis_pal(option = "B")(100),
         show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         annotation_col = data.frame(cl = as.factor(clusts)))
clusts[clusts == 3] = "O"
clusts[clusts == 2] = "M"
clusts[clusts == 5] = "T8"
clusts[clusts == 1] = "T4"
clusts[clusts == 6] = "LT"
clusts[clusts == 7] = "LB"
clusts[clusts == 4] = "LM"

# declare color of TME clusters
clustcols =c("grey30", "forestgreen", "orange", "red", "darkorchid3", "darkblue", "blue")
names(clustcols) = c("O", "M", "T8", 'T4', "LT", "LB", "LM")


dev.off()
svg("results/decon heatmap - TME only.svg", height = 3, width = 8)
pheatmap(pmin(res$cell.counts[p1$tree_row$order, use][, o], 100),
         col = viridis_pal(option = "B")(100),
         show_colnames = F,
         cluster_rows = F, cluster_cols = F)
dev.off()


svg("results/decon barplot - props - TME only.svg", height = 3, width = 8.1)
layout(mat = c(1,2,3), height = c(1.5, 2.5, 1.5))
par(mar = c(.2,0,0,3))
plot(color_branches(as.dendrogram(p1$tree_col), 
                    k = 7,
                    col = clustcols), 
     xlab = "", 
     ylab = "", 
     #edgePar = list(lwd = 2),
     leaflab = "none", 
     yaxt = "n", 
     ylim = c(-20, 250))
x0 = 0.5
for (cid in c("O", "M", "T8", "T4", "LT", "LB", "LM")) {
  rect(x0, -20, x0 + sum(clusts == cid), 0, col = clustcols[cid], border = NA)
  x0 = x0 + sum(clusts == cid)
}

# proportions:
barplot(sweep(res$cell.counts, 2, colSums(res$cell.counts), "/")[, use][, o], col = cellcols[immcells],
        border = NA,
        names.arg = rep("", length(o)),
        yaxt = "n",
        ylab = "Proportion")

# abundances:
barplot(res$cell.counts[, use][, o], col = cellcols[immcells],
        border = NA,
        names.arg = rep("", length(o)),
        yaxt = "n", 
        ylab = "Total cells")
axis(4, at = c(0, 400, 800, 1200), las = 2)
dev.off()


svg(paste0("results/clusters in space.svg"), width = 5, height = 4)
layout(mat = matrix(c(1,2), nrow = 1), widths = c(8, 2))
par(mar = c(0,0,0,0))
plot(annot[names(clusts), "x"], annot[names(clusts), "y"],
     pch = 16, 
     col = alpha(clustcols[clusts], 0.7),
     cex = colSums(res$cell.counts[, names(clusts)]) * 0.005 + 0.1,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
     xlim = c(0, 9000), ylim = c(0, 9000))
polygon(bound$x, bound$y, col = rgb(0,0,0,0.07), border = NA)
frame()
legend("right", pch = c(NA, rep(16, length(clustcols))), 
       col = c(NA, alpha(clustcols, 0.7)), 
       legend = c("Cluster", names(clustcols)), 
       cex = 1, bty = "n")
dev.off()



