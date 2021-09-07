#### validation in serial sections:

rm(list = ls())
library(scales)
library(SpatialDecon)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(Giotto)

source("DWLS_functions.R")
source("deconvolution functions - NNLS and vSVR and DWLS.R")

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

#### save data for benchmarking dataset --------------------------------

write.csv(annot[, c("tissue", "AOI_type")], file = "benchmarking_data/segment_annotations.csv")
write.csv(norm, file = "benchmarking_data/normalized_rna.csv")
write.csv(rna, file = "benchmarking_data/raw_rna.csv")
write.csv(pnorm, file = "benchmarking_data/normalized_protein.csv")
write.csv(prot, file = "benchmarking_data/raw_protein.csv")

#### run decon --------------------------------------

usegenes = intersect(rownames(norm), rownames(safeTME))

betas = list()

# run initial decon without tumor:
res0 = spatialdecon(norm = norm[usegenes, ],
                    raw = rna[usegenes, ],
                    bg = replace(norm[usegenes, ], TRUE, 1), 
                    is_pure_tumor = NULL) 

betas$spatialdecon.ignoretumor = res0$beta

# while modeling tumor columns:
res = spatialdecon(norm = norm[usegenes, ],
                      raw = rna[usegenes, ],
                      bg = replace(norm[usegenes, ], TRUE, 1), 
                      is_pure_tumor = annot$AOI_type == "Tumor",
                      n_tumor_clusters = 10) 
betas$spatialdecon.modeltumor = res$beta

# run NNLS:
nnls = deconUsingConstrainedLM(Y = norm[usegenes, ], 
                               X = SpatialDecon::safeTME[usegenes, ], 
                               bg = replace(norm[usegenes, ], TRUE, 1))
betas$NNLS = nnls$beta

# run vSRV:
vsvr = deconSVR(Y = norm[usegenes, ], 
                X = SpatialDecon::safeTME[usegenes, ], 
                bg = replace(norm[usegenes, ], TRUE, 1))
colnames(vsvr$beta) = colnames(norm)
betas$VSVR = vsvr$beta

# run DWLS:
dw = deconUsingDWLS(Y = norm[usegenes, ], 
                    X = SpatialDecon::safeTME[usegenes, ], 
                    bg = replace(norm[usegenes, ], TRUE, 1))
betas$DWLS = dw$beta

save(betas, file = "all decon results.RData")

# load stereoscope results and align:
stereo.safeTME = t(read.table("stereoscope_results.tsv", sep = "\t", header = T, row.names = 1))
betas[["stereo"]] = stereo.safeTME[rownames(betas[[1]]), colnames(betas[[1]])]


# load stereoscope results from lung initialization:
stereo.lung = t(read.table("stereoscope_lung_params.tsv", sep = "\t", header = T, row.names = 1))
betas[["stereo.lung"]] = betas[["stereo"]]*0
betas[["stereo.lung"]]["fibroblasts", ] = stereo.lung["Fibroblasts", colnames(betas[[1]])]
betas[["stereo.lung"]]["B.naive", ] = stereo.lung["tB.cells", colnames(betas[[1]])]
betas[["stereo.lung"]]["T.CD4.naive", ] = stereo.lung["tT.cells", colnames(betas[[1]])]
betas[["stereo.lung"]]["neutrophils", ] = stereo.lung["tNeutrophils", colnames(betas[[1]])]

## run spatialDWLS:
if (FALSE) {
  # format as giotto object:
  sign_matrix = SpatialDecon::safeTME[usegenes, ]
  save(rna, norm, annot, sign_matrix, usegenes, file = "for debugging dwls.RData")
  load("for debugging dwls.RData")
  exprs_list = list(raw = rna[usegenes, ], normalized = norm[usegenes, ])
  # giotto master version:
  gobj = Giotto::createGiottoObject(
    raw_exprs = exprs_list$raw,
    #expression_feat = c('rna'),
    spatial_locs = data.frame(sdimx = rnorm(nrow(annot)), 
                              sdimy = rnorm(nrow(annot)),
                              cell_ID = rownames(annot)),
    cell_metadata = annot)
  gobj@cell_metadata$rna$cluster = substr(annot$aoi.id, nchar(annot$aoi.id), nchar(annot$aoi.id))
  
  gobj <- normalizeGiotto(gobject = gobj)
  gobj <- runPCA(gobject = gobj, feats_to_use = usegenes, scale_unit = F)
  signPCA(gobj, genes_to_use = usegenes, scale_unit = F)
  gobj <- createNearestNetwork(gobject = gobj, dimensions_to_use = 1:10, k = 10)
  gobj <- doLeidenCluster(gobject = gobj, resolution = 0.4, n_iterations = 1000)
  # run spatialdwls:
  sdw = runDWLSDeconv(gobject = gobj,  
                      expression_values = c("normalized"),
                      logbase = 2,
                      cluster_column = "leiden_clus",
                      sign_matrix = sign_matrix,
                      n_cell = 100,
                      cutoff = 2,
                      name = NULL,
                      return_gobject = TRUE)
  tempmat = (sdw@spatial_enrichment$DWLS)[,-1]
  tempmat = apply(tempmat, 2, as.numeric)
  rownames(tempmat) = sdw@spatial_enrichment$DWLS$cell_ID
  save(tempmat, file = "spatialDWLS results.RData")
}
load("spatialDWLS results.RData")
betas[["spatialDWLS"]] = t(as.matrix(tempmat))


#### compare spatialdecon vs. protein in detail: -------------------------------------------

# define cell-protein matches:
cpmatch = list()
cpmatch[["CD3"]] = c("T.CD4.naive", "T.CD4.memory", "T.CD8.naive", "T.CD8.memory", "Treg")
cpmatch[["CD8"]] = c("T.CD8.naive", "T.CD8.memory")
cpmatch[["CD68"]] = "macrophages"
cpmatch[["CD66b"]] = "neutrophils"
cpmatch[["CD20"]] = c("B.naive", "B.memory")
cpmatch[["SMA"]] = "fibroblasts"


#### format as data frame for ggplot:
tempbeta = betas$spatialdecon.modeltumor
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
svg("serial sections - prot vs cell scores.svg", width = 9, height = 9)
g = ggplot(data = plotdf, aes(x = protexpr, y = betahat, col = AOI_type)) + 
  geom_point(alpha = 0.5, size = 2) +
  facet_grid(tumor ~ protid, scales = "free") +
  scale_y_continuous(limits=c(0, 10)) +
  theme_few() + 
  scale_color_manual(values = aoicols) +
  labs(x = "Protein expression", y = "Cell abundance estimated from RNA") + 
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17),
        strip.text.x = element_text(size = 12), legend.title = element_blank()) +
  theme(legend.position="bottom") + theme(legend.text = element_text(size=16))
print(g)
dev.off()



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
      # set NA's (from decon results with 0 variability) to zero:
      if (name == "spatialDWLS") {
        cors[[name]][is.na(cors[[name]])] = 0
        spearmans[[name]][is.na(spearmans[[name]])] = 0
      }
    }
  }
}


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

plotdf$method[plotdf$method == "spatialdecon.modeltumor"] = "SpatialDecon modelling tumor"
plotdf$method[plotdf$method == "spatialdecon.ignoretumor"] = "SpatialDecon ignoring tumor"
plotdf$method[plotdf$method == "stereo"] = "Stereoscope + safeTME"
plotdf$method[plotdf$method == "stereo.lung"] = "Stereoscope + scRNAseq lung profiles"
plotdf$method[plotdf$method == "VSVR"] = "v-SVR"
plotdf$method[plotdf$method == "spatialDWLS"] = "SpatialDWLS"

plotdf$method = factor(plotdf$method, levels = c("SpatialDecon modelling tumor",
                                                 "SpatialDecon ignoring tumor",
                                                 "NNLS", "v-SVR", "DWLS", "SpatialDWLS",
                                                 "Stereoscope + safeTME",
                                                 "Stereoscope + scRNAseq lung profiles"))

colmap = c(
    "SpatialDecon modelling tumor" = "red",
    "SpatialDecon ignoring tumor" = "orange",
    "NNLS" = "grey30", 
    "v-SVR" = "chartreuse3", 
    "DWLS" = "darkviolet",
    "SpatialDWLS" = "#330099",
    "Stereoscope + safeTME" = "dodgerblue4",
    "Stereoscope + scRNAseq lung profiles" = "cornflowerblue")
  

use = (plotdf$stat == "cors") #& (is.element(plotdf$method, c("DWLS", "spatialdecon.ignoretumor", "spatialdecon.modeltumor")))
svg("serial sections - method comparison.svg", width = 9, height = 9)
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
        strip.text.x = element_text(size = 12), legend.title = element_blank()) +
  theme(legend.position="bottom") + theme(legend.text = element_text(size=10))
print(g)
dev.off()


# save betas:
write.csv(plotdf, file = "correlations per tissue and per cell.csv")

# statistics:
ps = deltas = means = cis = dfs = c()
for (name in names(cors)) {
  mod = t.test(as.vector(cors$spatialdecon.modeltumor) - as.vector(cors[[name]]))
  ps[name] = mod$p.value
  deltas[name] = mod$estimate
  means[name] = mean(cors[[name]])
  cis[name] = sd(cors[[name]]) / sqrt(length(cors[[name]]))
}


# look at p-values vs. spatialdecon:
print(ps)

t.test(as.vector(cors$spatialdecon.modeltumor) - as.vector(cors$VSVR))
t.test(as.vector(cors$spatialdecon.modeltumor) - as.vector(cors$NNLS))
t.test(as.vector(cors$spatialdecon.modeltumor) - as.vector(cors$spatialdecon.ignoretumor))

# compare two stereoscope runs:
t.test(as.vector(cors$stereo) - as.vector(cors$stereo.lung))


# plot summary statistics:
plotdf = data.frame(mean = means, ci = cis, method = names(means))
# remove stereoscope + lung scRNAseq, which isn't comparable to the rest since it's missing cell types:
plotdf = plotdf[setdiff(rownames(plotdf), "stereo.lung"), ]
# format names:
plotdf$method[plotdf$method == "spatialdecon.modeltumor"] = "SpatialDecon modelling tumor"
plotdf$method[plotdf$method == "spatialdecon.ignoretumor"] = "SpatialDecon ignoring tumor"
plotdf$method[plotdf$method == "stereo"] = "Stereoscope + safeTME"
plotdf$method[plotdf$method == "VSVR"] = "v-SVR"
plotdf$method[plotdf$method == "spatialDWLS"] = "SpatialDWLS"

plotdf$method = factor(plotdf$method, levels = c("SpatialDecon modelling tumor",
                                                 "SpatialDecon ignoring tumor",
                                                 "NNLS", "v-SVR", "DWLS", "SpatialDWLS",
                                                 "Stereoscope + safeTME")) 


svg("serial sections - mean correlation.svg", height = 9, width = 2.5)
g = ggplot(data = plotdf, aes(x = method, y = mean, fill = method, col = method)) +  
  geom_bar(alpha = 0.5, size = 1, stat = "identity", show.legend = FALSE) +
  theme_few() + 
  scale_y_continuous(limits=c(0, 0.85)) +
  #scale_x_discrete(labels = NULL) +
  scale_fill_manual(values = colmap) + 
  scale_color_manual(values = colmap, ) + 
  labs(x = "", y = "Mean correlation") + 
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17),
        strip.text.x = element_text(size = 16), legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle=-90, size = 15)) +
  geom_errorbar(aes(ymin=mean-2*ci, ymax=mean + 2*ci), width=.2,
                position=position_dodge(.9), col = I("black"))
print(g)
dev.off()

