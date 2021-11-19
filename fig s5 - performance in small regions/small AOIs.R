#### compare spatialdecon in small vs. large AOIs:

library(SpatialDecon)
library(pheatmap)
library(readxl)

### load data ----------------------

load("small AOI data.RData")
annot$area = as.numeric(annot$area)

# remove NTCs:
use = rownames(annot)[annot$segment != "NTC"]

# remove AOIS with insufficient signal:
annot$q9 = apply(raw, 2, quantile, 0.9)
use = intersect(use, rownames(annot)[(annot$q9 > 0)])

annot = annot[use, ]
raw = raw[, use]

# load xy positions:
xannot = as.data.frame(read_xlsx("data/KN136_with_small-annotations.xlsx"))
rownames(xannot) = xannot$Sample_ID
annot = cbind(annot, xannot[rownames(annot), c("nuclei","stage_x","stage_y")])
annot$x = as.numeric(annot$stage_x)
annot$y = as.numeric(annot$stage_y)


#### run spatial decon ------------------------------

# normalize:
norm = sweep(raw, 2, annot$q9, "/") * mean(annot$q9)
bg = sweep(norm*0, 2, norm["NegProbe-WTX", ], "+")


res = spatialdecon(norm = norm, 
                   bg = bg, 
                   X = safeTME, 
                   is_pure_tumor = annot$segment == "Tumor")

#### analyze decon results -----------------------------

# shorthand for tumor and immuen cells and AOIs:
immcells = rownames(res$beta)[!grepl("tumor", rownames(res$beta))]
tumorcells = rownames(res$beta)[grepl("tumor", rownames(res$beta))]

aois.tumor = annot$segment == "Tumor"
aois.tme = annot$segment == "Immune"


# explore CRC 3, the only slide with truly small AOIs:
slide = "CRC 3"
use = (annot$slide_name == slide) & (aois.tme)
o = order(annot$area[use])


svg("supp fig 3 - range of AOI sizes.svg", width = 10)
graphics::layout(mat = matrix(c(1:2, 3, 3), 2), 
       widths = c(6,2,6,2),
       heights = rep(c(1,3), 2))
par(mar = c(1,5,1.3,0))
barplot(2*sqrt(annot$area[use][o] / pi), xaxt = "n", 
        col = "dodgerblue3",
        ylab = "Effective diameter (um)", cex.lab = 1.2)
abline(h = 55, col = "darkred")
text(1,75,"55 um")
par(mar = c(4,5,0,0))
barplot(res$prop_of_nontumor[immcells, use][, o], xaxt="n",
        col = cellcols[immcells],
        xlab = "Regions",
        ylab = "Proportion of immune cells", cex.lab = 1.2,
        names.arg = rep("", sum(use)))
par(mar = c(4,0,0,0))
frame()
legend("center", fill = rev(cellcols[immcells]), legend = rev(immcells), cex = 1.25)
dev.off()
