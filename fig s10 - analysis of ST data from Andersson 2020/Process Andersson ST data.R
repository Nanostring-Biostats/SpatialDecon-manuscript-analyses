# script to run spatialdecon on the ST data for Andersson et al 2020

rm(list = ls())
library(SpatialDecon)
library(scales)
source("spaceplot utils.R")

#### load data from sample G1, which was highlighted in paper's figures: ---------------------------
temp = read.table(gzfile("G1.tsv.gz"))
temp[1:5,1:5]

raw = t(as.matrix(temp))
x = as.numeric(substr(rownames(temp), 1, unlist(gregexpr("x", rownames(temp))) - 1))
y = -as.numeric(substr(rownames(temp), unlist(gregexpr("x", rownames(temp))) + 1, nchar(rownames(temp))))

#### run decon: ---------------------------------------------
# find all tumor regions (low counts of safeTME vs. other genes):
sharedgenes = intersect(rownames(raw), rownames(safeTME))
plot(colSums(raw), colSums(raw[sharedgenes, ]), log = "xy")
hist(colSums(raw[sharedgenes, ]) / colSums(raw), breaks = 20)
alltumor = colSums(raw[sharedgenes, ]) / colSums(raw) < 0.03  # for alma data

table(alltumor)


# run decon:
res = spatialdecon(norm = raw, raw = raw, lower_thresh = 1,
                   bg = 0.01, X = safeTME, is_pure_tumor = alltumor)
save(res, file = "spatialdecon results in Andersson data.Rdata")



#### plot cells in space ------------------------------------------


# infer polygon:
bound = getBoundary(x, y, marg = 0.1)


# plots one cell at a time:
svg(paste0("Andersson decon results.svg"), height = 10, width = 10)
par(mfrow = c(5,4))
for (cell in c("T.CD4.naive","T.CD4.memory","T.CD8.naive",
               "T.CD8.memory","Treg","NK","B.naive","B.memory","plasma","pDCs",
               "macrophages","mDCs","monocytes.C","monocytes.NC.I","neutrophils","mast",
               "endothelial.cells" ,"fibroblasts")) {
  par(mar = c(0,0,0,0))
  plot(x = x,
       y = y,
       cex = res$beta[cell, ] * 25,
       col = alpha("dodgerblue4", 0.7), pch = 16,
       #col = alpha(cellcols[cell], 0.7), pch = 16,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
       asp = 1)
  text(mean(range(x)), max(y), cell, cex = 2)
  polygon(bound$x, bound$y, col = rgb(0,0,0,0.07), border = NA)
}
dev.off()


florets(x = x, y = y, b = res$beta, col = cellcols[rownames(res$beta)])
polygon(bound$x, bound$y, col = rgb(0,0,0,0.07), border = NA)


#### examine impact of background variable: -----------------------------

res.01 = spatialdecon(norm = raw, raw = raw, lower_thresh = 1,
                      bg = .01, X = safeTME, is_pure_tumor = alltumor)

res1 = spatialdecon(norm = raw, raw = raw, lower_thresh = 1,
                    bg = 1, X = safeTME, is_pure_tumor = alltumor)


cor(as.vector(res.01$beta[colnames(safeTME), ]), as.vector(res1$beta[colnames(safeTME), ]))
