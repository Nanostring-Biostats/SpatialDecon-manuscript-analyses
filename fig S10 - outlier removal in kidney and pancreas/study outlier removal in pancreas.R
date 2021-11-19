

rm(list = ls())
library(scales)
library(SpatialDecon)
library(readxl)

#### load data --------------------------------

annot = as.data.frame(read_xlsx("data/pancreas Q3 Norm.xlsx", sheet = "SegmentProperties"))
raw0 = read_xlsx("data/pancreas Q3 Norm.xlsx", sheet = "TargetCountMatrix")
raw = as.matrix(raw0[, -1])
rownames(raw) = raw0$TargetName
rm(raw0)

# normalize:
norm = sweep(raw, 2, apply(raw, 2, quantile, 0.9), "/") *  mean(apply(raw, 2, quantile, 0.9))

#### get pancreas profile matrix: -------------------------------

refmat = download_profile_matrix(species = "Human", age_group = "Adult", matrixname = "Pancreas_HCA")
sharedgenes = intersect(rownames(refmat), rownames(norm))

#### perturb the data: ---------------------------------------------

# save original normalized data then perturb norm:
norm0 = norm

# add noise to a subset of genes:
set.seed(0)
to.perturb = sample(sharedgenes, 1000, replace = F)
perturbation = matrix((rnorm(length(norm[to.perturb, ]), mean = 0, sd = 3)), length(to.perturb))
norm[to.perturb, ] = norm[to.perturb, ] * 2 ^ perturbation


#### run decon with no outlier removal on normalized data --------------------------------
res_truth_no = spatialdecon(norm = norm0[sharedgenes, ], 
                         raw = raw[sharedgenes, ],
                         bg = sweep(norm0*0, 2, norm["NegProbe-WTX", ], "+"), 
                         X = refmat[sharedgenes, ],
                         resid_thresh = Inf)  # no outlier removal
res_truth_yes = spatialdecon(norm = norm0[sharedgenes, ], 
                            raw = raw[sharedgenes, ],
                            bg = sweep(norm0*0, 2, norm["NegProbe-WTX", ], "+"), 
                            X = refmat[sharedgenes, ],
                            resid_thresh = 3)  # default outlier removal


#### run decon on perturbed dataset with and without outlier removal: ------------------------

# with outlier removal:
res_yes = spatialdecon(norm = norm[sharedgenes, ], 
                       raw = raw[sharedgenes, ],
                       bg = sweep(norm*0, 2, norm["NegProbe-WTX", ], "+"), 
                       X = refmat[sharedgenes, ],
                       resid_thresh = 3)   # default value for outlier removal
# without:
res_no = spatialdecon(norm = norm[sharedgenes, ], 
                      raw = raw[sharedgenes, ],
                      bg = sweep(norm*0, 2, norm["NegProbe-WTX", ], "+"), 
                      X = refmat[sharedgenes, ],
                      resid_thresh = Inf)  # no outlier removal

save(res_truth_yes, res_truth_no, res_no, res_yes, file = "results/decon results - pancreas.RData")

#### compare genes' freq. removed: -------------------------------

# simplest:
flagged = 1 * is.na(res_yes$resids)
nflagged = rowMeans(flagged[sharedgenes, ])
png("results/rate of outlier removal - pancreas.png", width = 4.5, height = 6, units = "in", res= 300)
boxplot(nflagged ~ is.element(sharedgenes, to.perturb), ylab = "Proportion of regions where gene is flagged",
        xlab = "", names = c("Unperturbed\n(14975 genes)", "Added noise\n(1000 genes)"), cex.lab = 1.5, outline = F)
points(1 + jitter(as.numeric(is.element(sharedgenes, to.perturb))), nflagged, pch = 16, col = alpha("darkblue", 0.7))
dev.off()


#### compare to "truth" decon results -------------------------------

cors = data.frame("no" = diag(cor(t(res_truth_no$beta), t(res_no$beta))),
                  "yes" = diag(cor(t(res_truth_yes$beta), t(res_yes$beta))))
mses = data.frame("no" = rowMeans((res_truth_no$beta - res_no$beta)^2),
                  "yes" = rowMeans((res_truth_yes$beta - res_yes$beta)^2))


svg("results/concordance with truth - pancreas.svg", width = 6, height = 4)
par(mfrow = c(1,2))
par(mar = c(5,6,2,1))
plot(rowMeans(res_truth_no$beta)[rownames(cors)], cors$no, pch = 16, 
     ylab = "Correlation between results from\nperturbed and unperturbed data", 
     xlab = "Mean cell type abundance")
points(rowMeans(res_truth_yes$beta)[rownames(cors)], cors$yes, col = 2, pch = 16)
lines(lowess(cors$no ~ rowMeans(res_truth_no$beta)[rownames(cors)]))
lines(lowess(cors$yes ~ rowMeans(res_truth_yes$beta)[rownames(cors)]), col = 2)
legend("bottomright", pch = 16, col = c(1,2), legend = c("Outliers retained", "Outliers removed"), cex = 0.8)

plot(rowMeans(res_truth_no$beta), mses$no, pch = 16, 
     ylab = "Mean Squared Error between results from\nperturbed and unperturbed data", 
     xlab = "Mean cell type abundance")
points(rowMeans(res_truth_yes$beta), mses$yes, col = 2, pch = 16)
lines(lowess(mses$no ~ rowMeans(res_truth_no$beta)))
lines(lowess(mses$yes ~ rowMeans(res_truth_yes$beta)), col = 2)
dev.off()


