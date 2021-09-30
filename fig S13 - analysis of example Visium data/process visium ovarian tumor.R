rm(list = ls())
library(SpatialDecon)
library(scales)
source("spaceplot utils.R")

#### load and process visium data -------------------------------------

load("filtered_data.Rdata")
# extract counts matrix:
raw = as.matrix(filtered_data@assays$Spatial@counts)
sharedgenes = intersect(rownames(raw), rownames(safeTME))

# extract positions:
x = filtered_data@images$slice1@coordinates[colnames(raw), "row"]
y = filtered_data@images$slice1@coordinates[colnames(raw), "col"]
keep = !is.na(x) & (colSums(raw[sharedgenes, ]) > 0)
raw = raw[, keep]
x = x[keep]
y = y[keep]


#### run decon: -------------------------------------

# find all tumor regions (low counts of safeTME vs. other genes):
plot(colSums(raw), colSums(raw[sharedgenes, ]), log = "xy")
hist(colSums(raw[sharedgenes, ]) / colSums(raw), breaks = 60)
propsafeTMEcounts = colSums(raw[sharedgenes, ]) / colSums(raw)
alltumor = propsafeTMEcounts <= propsafeTMEcounts[order(propsafeTMEcounts)[100]]

# run decon:
rawmat = as.matrix(raw[sharedgenes, ])
use = TRUE
if (FALSE) {
  res = spatialdecon(norm = rawmat[,use], raw = rawmat[,use], lower_thresh = 0.5,
                   bg = 0.01, X = safeTME[sharedgenes, ], is_pure_tumor = alltumor[use])
  save(res, file = "visium example data tumor results.RData")
}
load("visium example data tumor results.RData")
#### plot cells in space ------------------------------------------

# infer polygon:
bound = getBoundary(x, y, marg = 0.1)

# plots one cell at a time:
svg(paste0("Visium decon results.svg"), height = 6.5, width = 8)
par(mfrow = c(5,4))
for (cell in c("T.CD4.naive","T.CD4.memory","T.CD8.naive",
               "T.CD8.memory","Treg","NK","B.naive","B.memory","plasma","pDCs",
               "macrophages","mDCs","monocytes.C","monocytes.NC.I","neutrophils","mast",
               "endothelial.cells" ,"fibroblasts")) {
  par(mar = c(0,0,0,0))
  plot(y, x,
       cex = pmax((res$beta[cell, ] - 0.1) * 5, 0),
       col = alpha("dodgerblue4", 0.7), pch = 16,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
       asp = 1)
  text(mean(range(y)), max(x)+6, cell, cex = 1.55)
  polygon(bound$y, bound$x, col = rgb(0,0,0,0.07), border = NA)
}
dev.off()


florets(x = x, y = y, b = res$beta, col = cellcols[rownames(res$beta)])
polygon(bound$x, bound$y, col = rgb(0,0,0,0.07), border = NA)


