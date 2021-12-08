library(svglite)
library(tidyr)
library(dplyr)
library(ggplot2)
library(here)
library(RColorBrewer)
library(readxl)
require(ggplotify)
library(ggridges)
library(ComplexHeatmap)
library(patchwork)

# retrieve commonly shared genes across all TCGA cancer types 
filesList <- list.files(here("fig 3 - gene rates in tumors", "output"), ".rdata")
uniqueGene <- c()

for(i in seq_along(filesList)) { 
  fileName <- strsplit(filesList[i], "_res.rdata")[[1]]
  load(here("fig 3 - gene rates in tumors", "output", filesList[i]))
  if(length(uniqueGene)==0){
    uniqueGene <- c(uniqueGene, as.character(res$gene))
  }
  uniqueGene <- intersect(uniqueGene, as.character(res$gene))
}

# combine expression data from TCGA cancer types 
## assign each saved data into res
for(i in seq_along(filesList)) { 
  fileName <- strsplit(filesList[i], "_res.rdata")[[1]]
  load(here("fig 3 - gene rates in tumors", "output", filesList[i]))
  rownames(res) <- res$gene
  newres <- res[uniqueGene,c(2,3)]
  colnames(newres) <- c("gene", fileName)
  newres$gene <- as.character(newres$gene)
  assign(paste("res", i, sep = ""), newres)
}

## combine results from each TCGA into one data.frame called res_combined
res_combined <- res1
for (i in setdiff(seq_along(filesList), 1)){
  res_combined <- inner_join(res_combined, get(paste("res", i, sep = "")), by = "gene") 
}

rownames(res_combined) <- res_combined[,1]
res_combined <- res_combined[,-1]

# Data visualization
## generate the heatmap to obtain the clustering 
set.seed(2020)
p <- pheatmap::pheatmap(t(res_combined),
                        fontsize = 10, show_colnames = FALSE,
                        treeheight_row = 0, treeheight_col = 0)
## get the data ready for heatmap
res_combined$gene <- rownames(res_combined)
res_long <- pivot_longer(res_combined, -gene, values_to = "value",
                         names_to = "cancer")
res_long$cancer <- factor(res_long$cancer,
                          levels = rev(colnames(res_combined)[p$tree_row$order]))
res_long$gene <- factor(res_long$gene,
                        levels = rownames(res_combined)[p$tree_col$order])

## prepare the color palette for heatmap 
col_fun = circlize::colorRamp2(seq(0, 100), 
                               colorRampPalette(c("#F2F0F7", "#4A1486"))(101))
mat <- as.matrix(t(res_combined[, -which(colnames(res_combined)=="gene")]))
mat <- mat[colnames(res_combined)[p$tree_row$order], 
           rownames(res_combined)[p$tree_col$order]]
## list genes to annotate 
geneannot <- c("CD3D", "CD3E", "CD8A", "CD8B",
               "CD163", "CD68", "CD34", "PTPRC", "MS4A1",
               "CD19", "CD4", "FOXP3", "EPCAM", 'KRT18')
locs <-  which(colnames(mat) %in% geneannot)
labels <- colnames(mat)[locs]
ca <- columnAnnotation(mark = anno_mark(at = locs, labels = labels,
                                        labels_gp = gpar(fontsize = 10)))
## create the heatmap
pr <- Heatmap(mat, name = " ",
              show_column_names=FALSE, 
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              col = col_fun,
              top_annotation = ca,
              row_names_gp = gpar(fontsize = 10) 
)
pr

## create the ridgeline plot version of the data in the heatmap
pl <- ggplot(res_long, aes(x = value, y = cancer, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.001) +
  scale_fill_gradient(low = "#F2F0F7", high = "#4A1486", 
                      breaks = seq(0, 100, 20)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        axis.line = element_line(colour = 'black'))

## generate the scatterplot for eight panels   
### filter out cancer types 
if(FALSE){ ## not run due to unavailability of the TCGA raw data
  filesList <- list.files(here("fig 3 - gene rates in tumors", "output"), ".RData")
  mergedatList <- NULL
  for(file in filesList){
    load(here("fig 3 - gene rates in tumors", "output", file))
    load(here("fig 3 - gene rates in tumors", "data", gsub("_res.rdata", ".rdata", file)))
    cancer <- gsub("_res.rdata", "", file)
    
    ## exclude several cancer types for plotting 
    if(file %in% c("LAML_res.rdata", "DLBC_res.rdata", "THYM_res.rdata")){
      next
    }
    
    dat <- get(paste0(cancer,".dat.subset"))$e
    dat <- dat[, -which(colSums(dat)==0)]
    
    if(length(which(duplicated(colnames(dat))))>0){
      dat <- dat[, -which(duplicated(colnames(dat)))]
    }
    
    ## collect the average expression level across cancer types 
    meandat <- data.frame(gene = colnames(dat), 
                          mean = apply(log2(dat+1), 2, mean),
                          count = nrow(dat), 
                          stringsAsFactors = FALSE)
    
    ## combine the average expression level with the intercept terms (beta0) from the log normal models 
    res$gene <- as.character(res$gene)
    colnames(res)[3] <- "beta0"
    
    mergedat <- inner_join(meandat, res[, c("gene", "beta0")], by = "gene")
    rownames(mergedat) <- mergedat$gene
    mergedat$type <- cancer
    mergedatList[[cancer]] <- mergedat[uniqueGene,]
  }
}

# load the gene lists from eight panels and the color scheme
load(here("fig 3 - gene rates in tumors", "data", "gene", "gene.rdata"))
load(here("fig 3 - gene rates in tumors", "data", "gene", "cellcols.RData"))
load(here("fig 3 - gene rates in tumors", "data", "mergedatList.rdata"))

## beta0 across TCGA cancer type
ntotal <- unique(rowSums(sapply(mergedatList, 
                                function(list) list$count)))
mat2 <- mergedatList[[1]]
mat2$gene <- mergedatList[[1]]$gene
mat2$mean <- rowSums(sapply(mergedatList, 
                            function(list) list$mean*list$count))/ntotal
mat2$beta0 <- rowMeans(sapply(mergedatList, function(list) list$beta0))
genecolor <- alpha("#E74C3C", 0.5) #alpha("orange", 0.5)


## prepare plotting specifications for scatter plots of eight panels 
plot_list <- rep(list(mat2), 8)
for(k in 1:8){
  plot_list[[k]]$title <- c("CIBERSORT", "EPIC", "MCP-counter",
                            "quanTIseq", "Timer", "xCell", 
                            "Danaher (2017)", "SafeTME")[k]
  
  
  #transgrey <- adjustcolor("#e6e6e6", alpha.f = 0.25) #adjustcolor("grey", alpha.f = 0.1)
  transgrey <- rgb(0,0,0,0)
  plot_list[[k]]$color <- transgrey
  #plot_list[[k]]$color <- rgb(0,0,0,0)
  
  if (k == 7){
    geneshighlight <- read_excel(here("fig 3 - gene rates in tumors", 
                                      "data", "gene", "Nanostring",
                                      "40425_2017_215_MOESM1_ESM.xlsx"), 
                                 sheet = "S4. Selected markers", 
                                 col_names = c("gene", "cell"), skip = 1)
    
    plot_list[[k]]$color[which(plot_list[[k]]$gene %in% geneshighlight$gene)] <- genecolor  
  } else if (k == 8 ) {
    load(here("fig 3 - gene rates in tumors", "data", "gene", "safeTME.RData"))
    geneshighlight <- rownames(safeTME)
    
    plot_list[[k]]$color[which(plot_list[[k]]$gene %in% geneshighlight)] <- genecolor 
  } else {
    plot_list[[k]]$color[which(plot_list[[k]]$gene %in% geneList[[k]]$gene)] <- genecolor
  }
  
  plot_list[[k]] <- plot_list[[k]][order(plot_list[[k]]$color, decreasing = FALSE), ]
}

mat2 <- do.call(rbind, plot_list)
mat2$title <- factor(mat2$title, 
                     levels = c("CIBERSORT", "EPIC", "MCP-counter",
                                "quanTIseq", "Timer", "xCell", 
                                "Danaher (2017)", "SafeTME"))
# remove genes outside of gene lists:
mat2 <- mat2[mat2$color == genecolor, ]

## combine all eight scatter plots into one plot
pb <- ggplot(mat2, aes(x = mean, y = beta0)) +
  geom_point(aes(color = color), alpha = 0.75,
             size = ifelse(mat2$color==transgrey, 0.25, 1))+
  scale_color_identity(guide = "legend")+
  guides(size = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(color = "Included in\ngene lists") +
  facet_wrap(~title, nrow = 2) +
  xlab(expression("Mean "*log[2]*"(expression)")) +
  ylab("Percent of transcripts attributed to cancer cells") +
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 13), 
        plot.margin = unit(c(-10,0,0,0), "pt"))


# Combine the heatmap, ridgeline plot, and the scatterplots into one plot
gb1 <- pl + theme(plot.margin = unit(c(1.2,0.5,0,0), "cm"))
gb2 <- grid.grabExpr(draw(pr, padding  = unit(c(1.1, 0, 0, 1), "cm"))) #+ theme(plot.margin = unit(c(1.2,0.5,0,1), "cm"))

p <- (gb1 + gb2 + plot_layout(widths = c(1, 2)))/pb +  
  plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = c('a'))

#ggsave(file = here("fig 3 - gene rates in tumors", "output", paste0("combined.png")),
#       dpi = 300, type = "cairo",
#       width = 10, height = 11)


ggsave(file = here("fig 3 - gene rates in tumors", "output", paste0("combined.pdf")),
       width = 10, height = 11)


#ggsave(file = here("fig 3 - gene rates in tumors", "output", paste0("combined.svg")),
#       width = 10, height = 11)


