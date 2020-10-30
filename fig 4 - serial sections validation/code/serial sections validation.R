#### validation in serial sections:

rm(list = ls())
library(scales)
library(SpatialDecon)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(here)

#### load data --------------------------------
folderpath = "fig 4 - serial sections validation"
rna = as.matrix(read.csv(here(folderpath, "data", "raw RNA counts.csv"), row.names = 1, header = T, stringsAsFactors = F))
prot = as.matrix(read.csv(here(folderpath, "data", "raw prot counts.csv"), row.names = 1, header = T, stringsAsFactors = F))
annot = read.csv(here(folderpath, "data", "AOI annotations.csv"), row.names = 1, header = T, stringsAsFactors = F)



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
icpgenes = read.csv(here(folderpath, "data", "ICP gene list.csv"), stringsAsFactors = F, header = F)[, 1]


# get expected background matrix
bg = rna * 0
bg[is.element(rownames(bg), icpgenes), ] = sweep(bg[is.element(rownames(bg), icpgenes), ], 2, rna["NegProbe-CTP01", ], "+")
bg[!is.element(rownames(bg), icpgenes), ] = sweep(bg[!is.element(rownames(bg), icpgenes), ], 2, rna["NegProbe-Kilo", ], "+")

# normalized data: signal to noise:
norm = rna / bg


#### save data for benchmarking dataset --------------------------------

write.csv(annot[, c("tissue", "AOI_type")], file = here(folderpath, "benchmarking_data", "segment_annotations.csv"))
write.csv(norm, file = here(folderpath, "benchmarking_data", "normalized_rna.csv"))
write.csv(rna, file = here(folderpath, "benchmarking_data", "raw_rna.csv"))
write.csv(pnorm, file = here(folderpath, "benchmarking_data", "normalized_protein.csv"))
write.csv(prot, file = here(folderpath, "benchmarking_data", "raw_protein.csv"))



#### run decon --------------------------------------
# run initial decon without tumor:
res0 = spatialdecon(norm = norm,
                    raw = rna,
                    bg = replace(norm, TRUE, 1), 
                    cellmerges = safeTME.matches,
                    is_pure_tumor = NULL) 

# while modelling tumor columns:
res = spatialdecon(norm = norm,
                   raw = rna,
                   bg = replace(norm, TRUE, 1), 
                   is_pure_tumor = annot$AOI_type == "Tumor", 
                   cellmerges = safeTME.matches,
                   n_tumor_clusters = 10) 

rlist = list(res0 = res0, res = res)

write.csv(res$beta, file = here(folderpath, "results", "cell abundance estimates.csv"))


#### compare decon vs. prot -------------------------------------------

# scatterplots stratified by tissue, and scatterplots colored by tissue

# define cell-prot matches:
cpmatch = list()
cpmatch[["CD3"]] = c("CD4.T.cells", "CD8.T.cells", "Treg")
cpmatch[["CD8"]] = "CD8.T.cells"
cpmatch[["CD68"]] = "macrophages"
cpmatch[["CD66b"]] = "neutrophils"
cpmatch[["CD20"]] = "B"
cpmatch[["SMA"]] = "fibroblasts"



#### format as data frame for ggplot:
tempbeta = res$beta
betahat = protexpr = aoiid = protid = c()
for (name in names(cpmatch)) {
  # get matching cell estimates:
  tempcells = cpmatch[[name]]
  # fill in data frame:
  betahat = c(betahat, colSums(tempbeta[tempcells, rownames(pnorm), drop = F]))
  protexpr = c(protexpr, pnorm[rownames(pnorm), name])
  aoiid = c(aoiid, rownames(pnorm))
  protid = c(protid, rep(name, nrow(pnorm)))
}
plotdf = data.frame(betahat = betahat, protexpr = protexpr, aoiid = aoiid, protid = protid, stringsAsFactors = F)
plotdf$tissue = annot[plotdf$aoiid, "tissue"]
plotdf$AOI_type = annot[plotdf$aoiid, "AOI_type"]

# custom names
plotdf$protid[plotdf$protid == "CD20"] = "CD20 protein\nvs. B-cells"
plotdf$protid[plotdf$protid == "CD3"] = "CD3 protein\nvs. total T-cells"
plotdf$protid[plotdf$protid == "CD8"] = "CD8 protein\nvs. CD8 T-cells"
plotdf$protid[plotdf$protid == "CD68"] = "CD68 protein\nvs. macrophages"
plotdf$protid[plotdf$protid == "CD66b"] = "CD66b protein\nvs. neutrophils"
plotdf$protid[plotdf$protid == "SMA"] = "SMA protein\nvs. fibroblasts"
plotdf$tumor = paste0("Tumor ", 1:5)[as.numeric(as.factor(plotdf$tissue))]
plotdf$protid = factor(plotdf$protid, 
                       levels = c("CD3 protein\nvs. total T-cells",
                                  "CD8 protein\nvs. CD8 T-cells",
                                  "CD20 protein\nvs. B-cells",
                                  "CD68 protein\nvs. macrophages",
                                  "CD66b protein\nvs. neutrophils",
                                  "SMA protein\nvs. fibroblasts"))

#### plot prot vs. cell estimate, faceted by patient/ cell:
svg(here(folderpath, "results", "serial sections - prot vs cell scores.svg"), width = 10, height = 9)
g = ggplot(data = plotdf, aes(x = protexpr, y = betahat, col = AOI_type)) + 
  geom_point(alpha = 0.5, size = 2) +
  facet_grid(tumor ~ protid, scales = "free") +
  scale_y_continuous(limits=c(0, 10)) +
  theme_few() + 
  scale_color_manual(values = aoicols) +
  labs(x = "Protein expression", y = "Cell abundance estimated from RNA") + 
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17),
        strip.text.x = element_text(size = 12), legend.title = element_blank())
print(g)
dev.off()



#### compare results with, without tumor adjustment:

cors = list()
cors[["res"]] = cors[["res0"]] = matrix(NA, nrow = length(unique(annot$tissue)), ncol = length(cpmatch),
                                        dimnames = list(unique(annot$tissue), names(cpmatch)))

for (tiss in unique(annot$tissue)) {
  for (pname in names(cpmatch)) {
    for (name in names(rlist)) {
      tempcells = cpmatch[[pname]]
      tempbeta = colSums(rlist[[name]]$beta[tempcells, rownames(pnorm), drop = F])
      tempprot = pnorm[, pname]
      cors[[name]][tiss, pname] = cor(tempbeta[annot$tissue == tiss], tempprot[annot$tissue == tiss])
    }
  }
}

write.csv(round(sapply(cors, colMeans), 3), file = here(folderpath, "results", "supp table - prot vs. cell cor with and without tumor modelling.csv"))

write.csv(cors$res, file = here(folderpath, "results", "correlations per tissue and per cell.csv"))

