#### validation in serial sections:

rm(list = ls())
library(scales)
library(SpatialDecon)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)

#### load and process serial sections data --------------------------------

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

usegenes = intersect(rownames(norm), rownames(safeTME))

#### for each sample, run spatialdecon under varying G --------------------------


# define cell-protein matches:
cpmatch = list()
cpmatch[["CD3"]] = c("T.CD4.naive", "T.CD4.memory", "T.CD8.naive", "T.CD8.memory", "Treg")
cpmatch[["CD8"]] = c("T.CD8.naive", "T.CD8.memory")
cpmatch[["CD68"]] = "macrophages"
cpmatch[["CD66b"]] = "neutrophils"
cpmatch[["CD20"]] = c("B.naive", "B.memory")
cpmatch[["SMA"]] = "fibroblasts"

# set up decon options:
tissues = unique(annot$tissue)
Gs = c(2,3,5,8,10,15,20,25,30)
cors = array(NA, dim = c(length(unique(tissues)), length(cpmatch), length(Gs)), 
             dimnames = list(tissues, names(cpmatch), as.character(Gs)))

for (tiss in tissues) {
  use = annot$tissue == tiss
  print(tiss)
  for (G in Gs) {
    print(G)
    # run decon:
    res = spatialdecon(norm = norm[usegenes, use],
                       raw = rna[usegenes, use],
                       bg = replace(norm[usegenes, use], TRUE, 1), 
                       is_pure_tumor = (annot$AOI_type == "Tumor")[use],
                       n_tumor_clusters = G) 
    # get cor:
    for (pname in names(cpmatch)) {
      tempcells = cpmatch[[pname]]
      tempbeta = colSums(res$beta[tempcells, rownames(pnorm)[use], drop = F])
      tempprot = pnorm[use, pname]
      cors[tiss, pname, as.character(G)] = cor(tempprot, tempbeta)  
    }
  }
}


#### display results -------------------------

meancors = apply(cors, c(1,3), mean)
pheatmap(meancors, cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white", "dodgerblue"))(100),
         display_numbers = T)

# convert cors to list:
lcors = list()
for (G in Gs) {
  lcors[[as.character(G)]] = cors[,,as.character(G)]
}

  
svg("G vs. accuracy.svg", height = 8, width = 8)
o = order(Gs)
par(mfrow = c(2,3))
for (tiss in rownames(lcors[[1]])) {
  cols = brewer.pal(7, "Set1")[-6]
  names(cols) = colnames(lcors[[1]])
  plot(c(0,0), xlim = range(Gs), ylim = range(unlist(lcors)), col = 0, log = "x",
       xlab = "G (# of tumor clusters)", ylab = "SpatialDecon / protein correlation",
       cex.lab = 1.5, main = paste0("Sample ", tiss))
  for (prot in colnames(lcors[[1]])) {
    tempcors = c()
    for (i in 1:length(lcors)) {
      tempcors[i] = lcors[[i]][tiss, prot]
    }
    lines(Gs[o], tempcors[o], col = cols[prot])
    points(Gs[o], tempcors[o], col = cols[prot], pch = 16)
  }
  abline(v = 1, lty = 2)
}
frame()
legend("center", lty = 2, pch = c(rep(16, length(cols))), 
                                  col = c(cols, 1), 
                                  legend = c(names(cols)), cex = 1.5)
dev.off()




#### multi-tumor analysis --------------------------------

cors2 = list()
for (G in Gs) {
  cors2[[as.character(G)]] = matrix(NA, length(unique(annot$tissue)), length(cpmatch),
                                    dimnames = list(unique(annot$tissue), names(cpmatch)))
  print(G)
  # run decon:
  res = spatialdecon(norm = norm[usegenes, ],
                     raw = rna[usegenes, ],
                     bg = replace(norm[usegenes, ], TRUE, 1), 
                     is_pure_tumor = (annot$AOI_type == "Tumor"),
                     n_tumor_clusters = G) 
  # get cor:
  for (tiss in unique(annot$tissue)) {
    use = annot$tissue == tiss
    for (pname in names(cpmatch)) {
      tempcells = cpmatch[[pname]]
      tempbeta = colSums(res$beta[tempcells, rownames(pnorm)[use], drop = F])
      tempprot = pnorm[use, pname]
      cors2[[as.character(G)]][tiss, pname] = cor(tempprot, tempbeta)  
    }
  }
}

## display results:
svg("G vs. accuracy - multi-sample.svg", height = 8, width = 8)
o = order(Gs)
par(mfrow = c(2,3))
for (tiss in rownames(cors2[[1]])) {
  cols = brewer.pal(7, "Set1")[-6]
  names(cols) = colnames(cors2[[1]])
  plot(c(0,0), xlim = range(Gs), ylim = range(unlist(cors2)), col = 0, log = "x",
       xlab = "G (# of tumor clusters)", ylab = "SpatialDecon / protein correlation",
       cex.lab = 1.5, main = paste0("Sample ", tiss))
  for (prot in colnames(cors2[[1]])) {
    tempcors = c()
    for (i in 1:length(cors2)) {
      tempcors[i] = cors2[[i]][tiss, prot]
    }
    lines(Gs[o], tempcors[o], col = cols[prot])
    points(Gs[o], tempcors[o], col = cols[prot], pch = 16)
  }
  abline(v = 1, lty = 2)
}
frame()
legend("center", lty = 2, pch = c(rep(16, length(cols))), 
       col = c(cols, 1), 
       legend = c(names(cols)), cex = 1.5)
dev.off()