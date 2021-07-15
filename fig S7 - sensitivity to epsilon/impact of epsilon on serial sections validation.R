#### validation in serial sections:

rm(list = ls())
library(scales)
library(SpatialDecon)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)

#### load data --------------------------------

rna = as.matrix(read.csv("data/raw RNA counts.csv", row.names = 1, header = T, stringsAsFactors = F))
prot = as.matrix(read.csv("data/raw prot counts.csv", row.names = 1, header = T, stringsAsFactors = F))
annot = read.csv("data/AOI annotations.csv", row.names = 1, header = T, stringsAsFactors = F)


# color by AOI type (tumor vs TME):
annot$AOI_type = c(NA, "Tumor", "TME", "Mixed")[1 + grepl("Tumor", annot$aoi.id) + 
                                                2 * grepl("TME", annot$aoi.id) + 
                                                3 * grepl("Geometric", annot$aoi.id)]
aoicols = c("#990033", "#003399", "grey30")
names(aoicols) = c("TME", "Tumor", "Mixed")
annot$aoicol = aoicols[annot$AOI_type]

# color tissues:
tissuecols = brewer.pal(8, "Set3")[1:(length(unique(annot$tissue))+1)][-2]
names(tissuecols) = unique(annot$tissue)
annot$tissuecol = tissuecols[annot$tissue]

#### QC rna and prot data and remove AOIs with failed signal ----------------------------

# qc rna:
negprobes = rownames(rna)[grepl("Neg", rownames(rna))]
print(negprobes)
annot$neg = colMeans(rna[negprobes, , drop = F])
hist(annot$neg, breaks = 20)
fail.rna = rownames(annot)[annot$neg < 10]

# qc protein:
iggs = colnames(prot)[grepl("IgG", colnames(prot))]
annot$iggfactor = exp(rowMeans(log(prot[, iggs])))
pairs(prot[, iggs], pch =16, col = annot$tissuecol)
pairs(prot[, iggs], log = "xy", col = annot$tissuecol, pch = 16)
plot(annot$iggfactor, col = annot$tissuecol, pch = 16)

boxplot(annot$iggfactor ~ annot$tissue, outline = F, ylim = c(0, max(annot$iggfactor)))
points(jitter(as.numeric(as.factor(annot$tissue))), annot$iggfactor, col = annot$aoicol, pch = 16, cex = 1.3)
abline(h = 10)
fail.prot = rownames(annot)[annot$iggfactor < 15]

keep = setdiff(rownames(annot), c(fail.rna, fail.prot))

annot = annot[keep, ]
rna = rna[, keep]
prot = prot[keep, ]

#### normalize protein ----------------------------

pnorm = sweep(prot, 1, annot$iggfactor, "/")


#### prep RNA data for decon run -----------------------------

# load ICP gene list:
icpgenes = read.csv("data/ICP gene list.csv", stringsAsFactors = F, header = F)[, 1]

# get expected background matrix
bg = rna * 0
bg[is.element(rownames(bg), icpgenes), ] = sweep(bg[is.element(rownames(bg), icpgenes), ], 2, rna["NegProbe-CTP01", ], "+")
bg[!is.element(rownames(bg), icpgenes), ] = sweep(bg[!is.element(rownames(bg), icpgenes), ], 2, rna["NegProbe-Kilo", ], "+")

# normalized data: signal to noise:
norm = rna / bg

#### modified algorithm2 where epsilon is exposed as an argument: ------------------------------
alg2e <- function(Y, X, bg = 0, weights = NULL, resid_thresh = 3, lower_thresh = 0.5, 
          align_genes = TRUE, maxit = 1000, epsilon = NULL) 
{
  if (align_genes) {
    sharedgenes <- intersect(rownames(X), rownames(Y))
    Y <- Y[sharedgenes, ]
    X <- X[sharedgenes, ]
    if (is.matrix(bg)) {
      bg <- bg[sharedgenes, ]
    }
    if (is.matrix(weights)) {
      weights <- weights[sharedgenes, ]
    }
  }
  tidied <- SpatialDecon:::tidy_X_and_Y(X, Y)
  X <- tidied$X
  Y <- tidied$Y
  if ((length(bg) > 0) & (is.vector(bg))) {
    bg <- matrix(bg, nrow = length(bg))
  }
  out0 <- SpatialDecon:::deconLNR(Y = Y, X = X, bg = bg, weights = weights, 
                   epsilon = epsilon, maxit = maxit)
  out0$yhat <- X %*% out0$beta + bg
  out0$resids <- log2(pmax(Y, lower_thresh)) - log2(pmax(out0$yhat, 
                                                         lower_thresh))
  outliers <- SpatialDecon:::flagOutliers(Y = Y, yhat = out0$yhat, wts = weights, 
                           resids = out0$resids, resid_thresh = resid_thresh)
  Y.nooutliers <- replace(Y, outliers, NA)
  out <- SpatialDecon:::deconLNR(Y = Y.nooutliers, X = X, bg = bg, weights = weights, 
                  epsilon = epsilon)
  out$yhat <- X %*% out$beta + bg
  out$resids <- log2(pmax(Y.nooutliers, 0.5)) - log2(pmax(out$yhat, 
                                                          0.5))
  tempbeta <- out$beta
  tempse <- tempp <- tempt <- tempbeta * NA
  for (i in seq_len(ncol(tempse))) {
    tempse[, i] <- suppressWarnings(sqrt(diag(out$sigmas[, 
                                                         , i])))
  }
  tempt <- (tempbeta/tempse)
  tempp <- 2 * (1 - stats::pnorm(tempt))
  out$p <- tempp
  out$t <- tempt
  out$se <- tempse
  return(out)
}



#### run decon over a range of epsilons --------------------------------------

usegenes = intersect(rownames(norm), rownames(safeTME))
betas = list()

# default epsilon: 
eps.default = min(norm[usegenes, ][norm[usegenes, ] > 0])
ratios = c(0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 5, 10, 20)
epsilons = eps.default * ratios

for (i in 1:length(epsilons)) {
  print(i)
  res = alg2e(Y = norm[usegenes, ], 
              X = safeTME[usegenes, ], 
              bg = replace(norm[usegenes, ], TRUE, 1), 
              epsilon = epsilons[i])
  betas[[i]] = res$beta
}
save(betas, file = "alg2e results.RData")


#### evaluate accuracy: -------------------------------

# metric 1: mean cor with protein (show a heatmap for each eps)


#### compare spatialdecon vs. protein in detail: -------------------------------------------

# define cell-protein matches:
cpmatch = list()
cpmatch[["CD3"]] = c("T.CD4.naive", "T.CD4.memory", "T.CD8.naive", "T.CD8.memory", "Treg")
cpmatch[["CD8"]] = c("T.CD8.naive", "T.CD8.memory")
cpmatch[["CD68"]] = "macrophages"
cpmatch[["CD66b"]] = "neutrophils"
cpmatch[["CD20"]] = c("B.naive", "B.memory")
cpmatch[["SMA"]] = "fibroblasts"




#### compare results of all methods:

cors = spearmans = list()
for (i in 1:length(betas)) {
  cors[[i]] = spearmans[[i]] = matrix(NA, nrow = length(unique(annot$tissue)), ncol = length(cpmatch),
                        dimnames = list(unique(annot$tissue), names(cpmatch)))
}

for (tiss in unique(annot$tissue)) {
  for (pname in names(cpmatch)) {
    for (i in 1:length(betas)) {
      tempcells = cpmatch[[pname]]
      tempbeta = colSums(betas[[i]][tempcells, rownames(pnorm), drop = F])
      tempprot = pnorm[, pname]
      cors[[i]][tiss, pname] = cor(tempbeta[annot$tissue == tiss], tempprot[annot$tissue == tiss], method = "pearson", use = "complete")
      spearmans[[i]][tiss, pname] = cor(tempbeta[annot$tissue == tiss], tempprot[annot$tissue == tiss], method = "spearman", use = "complete")
    }
  }
}


svg("epsilon vs. accuracy.svg", height = 8, width = 8)
o = order(epsilons)
par(mfrow = c(2,3))
for (tiss in rownames(cors[[1]])) {
  cols = brewer.pal(7, "Set1")[-6]
  names(cols) = colnames(cors[[1]])
  plot(c(0,0), xlim = range(ratios), ylim = range(unlist(cors)), col = 0, log = "x",
       xlab = "Ratio to default epsilon", ylab = "Correlation with protein",
       cex.lab = 1.5, main = paste0("Sample ", tiss))
  for (prot in colnames(cors[[1]])) {
    tempcors = c()
    for (i in 1:length(cors)) {
      tempcors[i] = cors[[i]][tiss, prot]
    }
    lines(ratios[o], tempcors[o], col = cols[prot])
    points(ratios[o], tempcors[o], col = cols[prot], pch = 16)
  }
  abline(v = 1, lty = 2)
}
frame()
legend("center", lty = 2, pch = c(rep(16, length(cols)), NA), 
                                  col = c(cols, 1), 
                                  legend = c(names(cols), "default epsilon"), cex = 1.5)
dev.off()

