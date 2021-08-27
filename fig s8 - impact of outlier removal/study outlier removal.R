

rm(list = ls())
library(scales)
library(SpatialDecon)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)

#### load serial sections data --------------------------------

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


#### perturb the data: ---------------------------------------------

# save original normalized data then perturb norm:
norm0 = norm


usegenes = intersect(rownames(norm), rownames(safeTME))

# add noise to a subset of genes:
set.seed(0)
to.perturb = sample(usegenes, 50, replace = F)
perturbation = matrix((rnorm(length(norm[to.perturb, ]), mean = 0, sd = 3)), length(to.perturb))
norm[to.perturb, ] = norm[to.perturb, ] * 2 ^ perturbation

#### run decon on perturbed dataset with and without outlier removal: ------------------------

# with outlier removal:
res = spatialdecon(norm = norm[usegenes, ],
                   raw = rna[usegenes, ],
                   bg = replace(norm[usegenes, ], TRUE, 1), 
                   is_pure_tumor = annot$AOI_type == "Tumor",
                   n_tumor_clusters = 10,
                   resid_thresh = 3)   # default value for outlier removal
# without:
res2 = spatialdecon(norm = norm[usegenes, ],
                   raw = rna[usegenes, ],
                   bg = replace(norm[usegenes, ], TRUE, 1), 
                   is_pure_tumor = annot$AOI_type == "Tumor",
                   n_tumor_clusters = 10,
                   resid_thresh = Inf)  # no outlier removal


#### compare genes' freq. removed: -------------------------------

# simplest:
flagged = 1 * is.na(res$resids)
nflagged = rowMeans(flagged[usegenes, ])
svg("rate of outlier removal.svg", width = 4.5)
boxplot(nflagged ~ is.element(usegenes, to.perturb), ylab = "Proportion of AOIs where gene is flagged",
        xlab = "", names = c("Unperturbed\n(494 genes)", "Added noise\n(50 genes)"), cex.lab = 1.5, outline = F)
points(1 + jitter(as.numeric(is.element(usegenes, to.perturb))), nflagged, pch = 16, col = alpha("darkblue", 0.7))
dev.off()

#### compare performance vs. protein: ------------------------

betas = list(outliers.removed = res$beta,
             outliers.retained = res2$beta)

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
for (name in names(betas)) {
  cors[[name]] = spearmans[[name]] = matrix(NA, nrow = length(unique(annot$tissue)), ncol = length(cpmatch),
                                      dimnames = list(unique(annot$tissue), names(cpmatch)))
}

for (tiss in unique(annot$tissue)) {
  for (pname in names(cpmatch)) {
    for (name in names(betas)) {
      tempcells = cpmatch[[pname]]
      tempbeta = colSums(betas[[name]][tempcells, rownames(pnorm), drop = F])
      tempprot = pnorm[, pname]
      cors[[name]][tiss, pname] = cor(tempbeta[annot$tissue == tiss], tempprot[annot$tissue == tiss], method = "pearson", use = "complete")
      spearmans[[name]][tiss, pname] = cor(tempbeta[annot$tissue == tiss], tempprot[annot$tissue == tiss], method = "spearman", use = "complete")
      cors[[name]] = replace(cors[[name]], is.na(cors[[name]]), 0)
    }
  }
}

pheatmap(cors[[1]], display_numbers = T, cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white", "dodgerblue"))(100), breaks = seq(0, 1, length.out = 101))
pheatmap(cors[[2]], display_numbers = T, cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white", "dodgerblue"))(100), breaks = seq(0, 1, length.out = 101))






## plot results: -------------------

# assemble data frame for ggplot:
method = tissue = protein = stat = value = c()
for (obj in c("cors", "spearmans")) {
  temp = get(obj)
  for (name in names(temp)) {
    for (tiss in rownames(temp[[name]])) {
      method = c(method, rep(name, ncol(temp[[name]])))
      tissue = c(tissue, rep(tiss, ncol(temp[[name]])))
      protein = c(protein, colnames(temp[[name]]))
      stat = c(stat, rep(obj, ncol(temp[[name]])))
      value = c(value, temp[[name]][tiss, ])
    }
  }
}

plotdf = data.frame(method = method, tissue = tissue, protein = protein, stat = stat, value = value)

# custom names
plotdf$protein[plotdf$protein == "CD20"] = "CD20 protein\nvs. B-cells"
plotdf$protein[plotdf$protein == "CD3"] = "CD3 protein\nvs. total T-cells"
plotdf$protein[plotdf$protein == "CD8"] = "CD8 protein\nvs. CD8 T-cells"
plotdf$protein[plotdf$protein == "CD68"] = "CD68 protein\nvs. macrophages"
plotdf$protein[plotdf$protein == "CD66b"] = "CD66b protein\nvs. neutrophils"
plotdf$protein[plotdf$protein == "SMA"] = "SMA protein\nvs. fibroblasts"
plotdf$tumor = paste0("Tumor ", 1:5)[as.numeric(as.factor(plotdf$tissue))]
plotdf$protein = factor(plotdf$protein, 
                        levels = c("CD3 protein\nvs. total T-cells",
                                   "CD8 protein\nvs. CD8 T-cells",
                                   "CD20 protein\nvs. B-cells",
                                   "CD68 protein\nvs. macrophages",
                                   "CD66b protein\nvs. neutrophils",
                                   "SMA protein\nvs. fibroblasts"))

plotdf$method[plotdf$method == "outliers.removed"] = "Outliers removed"
plotdf$method[plotdf$method == "outliers.retained"] = "Outliers retained"


plotdf$method = factor(plotdf$method, levels = c("Outliers removed",
                                                 "Outliers retained"))

colmap = c(
  "Outliers removed" = "cornflowerblue",
  "Outliers retained" = "orange")


use = (plotdf$stat == "cors") 
svg("serial sections - with vs. without outliers.svg", width = 9, height = 9)
g = ggplot(data = plotdf[use, ], aes(x = method, y = value, fill = method, col = method, )) +  
  geom_bar(alpha = 0.5, size = 1, stat = "identity") +
  theme_few() + 
  facet_grid(tumor ~ protein, scales = "free") +
  scale_y_continuous(limits=c(-0.45, 1)) +
  scale_x_discrete(labels = NULL) +
  scale_fill_manual(values = colmap) + 
  scale_color_manual(values = colmap, ) + 
  labs(x = "Method", y = "Correlation") + 
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17),
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), 
        legend.title = element_blank()) +
  theme(legend.position="bottom") + theme(legend.text = element_text(size=16))
print(g)
dev.off()

mean(cors[[1]])
mean(cors[[2]])



#### now compare performance under various outlier thresholds -------------------------------


# try other outlier thresholds:
reslist = list(t3 = res)
reslist[["t2"]] = spatialdecon(norm = norm[usegenes, ],
                               raw = rna[usegenes, ],
                               bg = replace(norm[usegenes, ], TRUE, 1), 
                               is_pure_tumor = annot$AOI_type == "Tumor",
                               n_tumor_clusters = 10,
                               resid_thresh = 2) 
reslist[["t4"]] = spatialdecon(norm = norm[usegenes, ],
                               raw = rna[usegenes, ],
                               bg = replace(norm[usegenes, ], TRUE, 1), 
                               is_pure_tumor = annot$AOI_type == "Tumor",
                               n_tumor_clusters = 10,
                               resid_thresh = 4) 
reslist[["t5"]] = spatialdecon(norm = norm[usegenes, ],
                               raw = rna[usegenes, ],
                               bg = replace(norm[usegenes, ], TRUE, 1), 
                               is_pure_tumor = annot$AOI_type == "Tumor",
                               n_tumor_clusters = 10,
                               resid_thresh = 5) 
reslist[["t6"]] = spatialdecon(norm = norm[usegenes, ],
                               raw = rna[usegenes, ],
                               bg = replace(norm[usegenes, ], TRUE, 1), 
                               is_pure_tumor = annot$AOI_type == "Tumor",
                               n_tumor_clusters = 10,
                               resid_thresh = 5) 
reslist[["tinf"]] = res2

### compare performance:

cors = spearmans = list()
for (name in names(reslist)) {
  cors[[name]] = spearmans[[name]] = matrix(NA, nrow = length(unique(annot$tissue)), ncol = length(cpmatch),
                                            dimnames = list(unique(annot$tissue), names(cpmatch)))
}

for (tiss in unique(annot$tissue)) {
  for (pname in names(cpmatch)) {
    for (name in names(reslist)) {
      tempcells = cpmatch[[pname]]
      tempbeta = colSums(reslist[[name]]$beta[tempcells, rownames(pnorm), drop = F])
      tempprot = pnorm[, pname]
      cors[[name]][tiss, pname] = cor(tempbeta[annot$tissue == tiss], tempprot[annot$tissue == tiss], method = "pearson", use = "complete")
      spearmans[[name]][tiss, pname] = cor(tempbeta[annot$tissue == tiss], tempprot[annot$tissue == tiss], method = "spearman", use = "complete")
      cors[[name]] = replace(cors[[name]], is.na(cors[[name]]), 0)
    }
  }
}

round(sapply(cors, mean), 3)

pheatmap(cors[[1]], display_numbers = T, cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white", "dodgerblue"))(100), breaks = seq(0, 1, length.out = 101))
pheatmap(cors[[2]], display_numbers = T, cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white", "dodgerblue"))(100), breaks = seq(0, 1, length.out = 101))




# assemble data frame for ggplot:
method = tissue = protein = stat = value = c()
for (obj in c("cors", "spearmans")) {
  temp = get(obj)
  for (name in names(temp)) {
    for (tiss in rownames(temp[[name]])) {
      method = c(method, rep(name, ncol(temp[[name]])))
      tissue = c(tissue, rep(tiss, ncol(temp[[name]])))
      protein = c(protein, colnames(temp[[name]]))
      stat = c(stat, rep(obj, ncol(temp[[name]])))
      value = c(value, temp[[name]][tiss, ])
    }
  }
}

plotdf = data.frame(method = method, tissue = tissue, protein = protein, stat = stat, value = value)

# custom names
plotdf$protein[plotdf$protein == "CD20"] = "CD20 protein\nvs. B-cells"
plotdf$protein[plotdf$protein == "CD3"] = "CD3 protein\nvs. total T-cells"
plotdf$protein[plotdf$protein == "CD8"] = "CD8 protein\nvs. CD8 T-cells"
plotdf$protein[plotdf$protein == "CD68"] = "CD68 protein\nvs. macrophages"
plotdf$protein[plotdf$protein == "CD66b"] = "CD66b protein\nvs. neutrophils"
plotdf$protein[plotdf$protein == "SMA"] = "SMA protein\nvs. fibroblasts"
plotdf$tumor = paste0("Tumor ", 1:5)[as.numeric(as.factor(plotdf$tissue))]
plotdf$protein = factor(plotdf$protein, 
                        levels = c("CD3 protein\nvs. total T-cells",
                                   "CD8 protein\nvs. CD8 T-cells",
                                   "CD20 protein\nvs. B-cells",
                                   "CD68 protein\nvs. macrophages",
                                   "CD66b protein\nvs. neutrophils",
                                   "SMA protein\nvs. fibroblasts"))

plotdf$method[plotdf$method == "t2"] = "Threshold = 2"
plotdf$method[plotdf$method == "t3"] = "Threshold = 3"
plotdf$method[plotdf$method == "t4"] = "Threshold = 4"
plotdf$method[plotdf$method == "t5"] = "Threshold = 5"
plotdf$method[plotdf$method == "t6"] = "Threshold = 6"
plotdf$method[plotdf$method == "tinf"] = "Outliers retained"


plotdf$method = factor(plotdf$method, levels = c(paste0("Threshold = ", 2:6), "Outliers retained"))

colmap = c(
  "Threshold = 2" = "darkviolet",
  "Threshold = 3" = "red",
  "Threshold = 4" = "darkgoldenrod3",
  "Threshold = 5" = "cornflowerblue",
  "Threshold = 6" = "dodgerblue4",
  "Outliers retained" = "grey30")


use = (plotdf$stat == "cors") 
svg("serial sections - accuracy vs outlier threshold.svg", width = 9, height = 9)
g = ggplot(data = plotdf[use, ], aes(x = method, y = value, fill = method, col = method, )) +  
  geom_bar(alpha = 0.5, size = 1, stat = "identity") +
  theme_few() + 
  facet_grid(tumor ~ protein, scales = "free") +
  scale_y_continuous(limits=c(-0.45, 1)) +
  scale_x_discrete(labels = NULL) +
  scale_fill_manual(values = colmap) + 
  scale_color_manual(values = colmap, ) + 
  labs(x = "Method", y = "Correlation") + 
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17),
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), 
        legend.title = element_blank()) +
  theme(legend.position="bottom") + theme(legend.text = element_text(size=16))
print(g)
dev.off()



t.test(as.vector(cors[["t3"]]) - as.vector(cors[["t2"]]))
t.test(as.vector(cors[["t3"]]) - as.vector(cors[["t4"]]))
t.test(as.vector(cors[["t3"]]) - as.vector(cors[["t5"]]))
t.test(as.vector(cors[["t3"]]) - as.vector(cors[["t6"]]))
t.test(as.vector(cors[["t3"]]) - as.vector(cors[["tinf"]]))
