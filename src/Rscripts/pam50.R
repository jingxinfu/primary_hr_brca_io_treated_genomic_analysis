library(genefu)
library(reshape2)
library(ggdendro)
library(ggplot2)

# Steps based on the discussion with Jorge:
## Subytpe codes: https://unclineberger.org/peroulab/algorithms/
## 1. pre-normalize the tpm using the reference dataset following the code: https://unclineberger.org/peroulab/wp-content/uploads/sites/1008/2022/02/JCO-2020-Subgroup-specific-gene-centering-method-AFM.zip
## 2. Forge Genefu to not normalize the data
## 3. Taking the correlation scores from genefu other than the categorical subtype variable.
## 4. purity <0.1 should be normal-like subtypes.
#### - Because it will assign sample a subtype even if the Rho < 0.1, which should be categorized into NC (non-classified)
#### - Be careful of the scores, it will also assign sample to the subtype with largest probability when prob of lumA equals to 0.25 and prob of the lumB equals to 0.24.

# https://github.com/RMolania/TCGA_PanCancer_UnwantedVariation/blob/master/Scripts/TCGA_BRCA_RNAseq.Rmd

################################################### Load pre-normalized expression data ##################################################
mrna <- t(read.delim("data/processed/bulkRNA_PAM50_Normalized_Expression.csv", sep = ",", row.names = 1))

################################################### Make the gene name consistent with the one used in model ##################################################
data(pam50)
# colnames(mrna)["NUF2" == colnames(mrna)] <- "CDCA1"
# colnames(mrna)["NDC80" == colnames(mrna)] <- "KNTC2"
# colnames(mrna)["ORC6" == colnames(mrna)] <- "ORC6L"

################################################### Run PAM50 subtyping##################################################
# PAM50Preds <- molecular.subtyping(
#     sbt.model = "pam50.robust",
#     data = as.matrix(mrna[, pam50.robust$centroids.map$probe]),
#     annot = pam50.robust$centroids.map,
#     do.mapping = TRUE
# )


# sbt.model=pam50 means pam50 without data normalization. “pam50.robust” and “pam50.scale” will re-normalize data.
PAM50Preds <- intrinsic.cluster.predict(sbt.model = pam50, data = mrna, annot = pam50$centroids.map, do.mapping = TRUE) # This is using spearman

################################################### Run Clustering ##################################################
pam50.dendro <- as.dendrogram(hclust(d = dist(x = t(mrna))))
# Create dendrogram plot
dendro.plot <- ggdendrogram(
    data = pam50.dendro, rotate = TRUE
) + theme(axis.text.y = element_text(size = 6))
# melt data to format for ggplot
mrna.long <- melt(t(mrna))
# extract the column order in the dendrogram
pam50.order <- order.dendrogram(pam50.dendro)

# order the levels according to their position in the cluster
mrna.long$Var1 <- factor(
    x = mrna.long$Var1,
    levels = mrna.long$Var1[pam50.order],
    ordered = TRUE
)
# create heatmap plot
heatmap.plot <- ggplot(
    data = mrna.long, aes(x = Var2, y = Var1)
) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2() +
    theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top"
    )

# plot dendrogram and heatmap together
library("grid")
pdf("report/figure/bulkRNA_PAM50.pdf", width = 6, height = 6)
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 0.98))
print(dendro.plot, vp = viewport(x = 0.87, y = 0.416, width = 0.2, height = 0.98))
dev.off()
################################################### Output Subtype ##################################################
subtypes <- data.frame(Sample = gsub(".", "-", names(PAM50Preds$subtype), fixed = T), PAM50 = PAM50Preds$subtype)
subtypes$PAM50[apply(PAM50Preds$cor, 1, max) < 0.1] <- "NC"
write.table(
    subtypes,
    "data/processed/bulkRNA_PAM50.csv",
    sep = ",", quote = F, row.names = F
)