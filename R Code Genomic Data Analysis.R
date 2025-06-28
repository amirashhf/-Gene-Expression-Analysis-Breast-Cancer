## R Code Genomic Data Analysis 

## 1. Install & Load Library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase", "GEOquery", "limma", "genefilter", "annotate", "hgu133a.db", "GO.db"))

library(Biobase)
library(GEOquery)
library(limma)
library(genefilter)
library(annotate)
library(hgu133a.db)
library(GO.db)

## 2. Download Data dari GEO
dtgeo <- getGEO('GDS1329', destdir = ".")
dtgeo

## 3. Convert ke ExpressionSet
eset <- GDS2eSet(dtgeo, do.log2 = TRUE)
eset

## 4. Ambil Data Fenotipe
phdtgeo <- pData(eset)
head(phdtgeo)

## 5. Ambil Data Ekspresi Gen
expdtgeo <- exprs(eset)
dim(expdtgeo)
head(expdtgeo)

## 6. Load Platform Annotation
Meta(dtgeo)$platform
annotation(eset) <- "hgu133a"

## 7. Filter Gene (Menghapus Gen Varians Rendah dll)
esetFilt <- nsFilter(eset)$eset
expdtgeoFilt <- exprs(esetFilt)
dim(expdtgeoFilt)

## 8. Visualisasi Distribusi Data
par(mfrow = c(1, 2))
hist(expdtgeo, main = "Original Data")
hist(expdtgeoFilt, main = "Filtered Data")

## 9. Buat Variabel Grup
# Adjust sesuai data: misal 'apocrine', 'basal', 'luminal'
vargrp <- phdtgeo$disease.state
table(vargrp)

vargrp_clean <- gsub(" tumor", "", vargrp)
table(vargrp_clean)

# Assign grup kode
group <- ifelse(vargrp_clean == "apocrine", 0,
                ifelse(vargrp_clean == "basal", 1,
                       ifelse(vargrp_clean == "luminal", 2, NA)))
group <- factor(group)


## 9. Install ggplot2 (jika belum diinstal)
library(ggplot2)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)

# Pakai vargrp_clean sebagai faktor grup baru
group <- as.factor(vargrp_clean)

# Bikin desain matriks baru
design <- model.matrix(~ 0 + group)  # tanpa intercept
colnames(design) <- levels(group)

# Fit model
fit <- lmFit(expdtgeoFilt, design)

# Misal mau bandingkan 'basal' vs 'luminal'
contrast.matrix <- makeContrasts(BasalVsLuminal = basal - luminal, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Top table hasil basal vs luminal
topTableRes <- topTable(fit2, coef = "BasalVsLuminal", number = Inf)
head(topTableRes)

# Volcano plot
EnhancedVolcano(topTableRes,
                lab = rownames(topTableRes),
                x = 'logFC',
                y = 'P.Value',
                title = 'Basal vs Luminal',
                pCutoff = 0.05,
                FCcutoff = 1)


## 9. Buat Variabel Grup
# Adjust sesuai data: misal 'apocrine', 'basal', 'luminal'
vargrp <- phdtgeo$disease.state
table(vargrp)

vargrp_clean <- gsub(" tumor", "", vargrp)
table(vargrp_clean)

# Assign grup kode
group <- ifelse(vargrp_clean == "apocrine", 0,
                ifelse(vargrp_clean == "basal", 1,
                       ifelse(vargrp_clean == "luminal", 2, NA)))
group <- factor(group)


## 10. DE Analysis pakai Limma
design <- model.matrix(~group)
fit <- eBayes(lmFit(expdtgeoFilt, design))
summary(decideTests(fit))

## 11. Tampilkan Top Gene
topResult <- topTable(fit, coef = 2, number = 50)
topResult

## 12. Heatmap dari Top Genes
selected <- rownames(expdtgeoFilt) %in% rownames(topResult)
expdtgeosel <- expdtgeoFilt[selected, ]

heatmap(expdtgeosel)

## 13. Boxplot 4 Gene Teratas
par(mfrow = c(2, 2))
for (i in 1:4) {
  plot(vargrp, expdtgeosel[i, ], main = rownames(expdtgeosel)[i])
}

## 14. Gene Annotation
ids <- rownames(topResult)
GeneSelected <- AnnotationDbi::select(hgu133a.db, keys = ids, columns = c("SYMBOL", "ENTREZID", "GENENAME", "GO"), keytype = "PROBEID")
head(GeneSelected)

## 15. Gene Ontology
GOselected <- AnnotationDbi::select(GO.db, keys = GeneSelected$GO, columns = c("TERM", "GOID"), keytype = "GOID")
head(GOselected)

## 16. Combine Annotation
finalres <- cbind(GeneSelected, GOselected)
head(finalres)

## (Optional) Simpan hasil
write.csv(finalres, file = "C:/Users/LENOVO/DOwnloads/GeneAnnotation_GDS1329.csv", row.names = FALSE)


## 17. Pilih Top Variance Genes untuk Klastering
genes.var <- apply(expdtgeoFilt, 1, var)
genes.var.select <- order(-genes.var)[1:100]

data.s <- expdtgeoFilt[genes.var.select, ]
dim(data.s)

## 18. Cari Jumlah Klaster Optimal (Gap Statistic)
set.seed(77)
library(cluster)

gap_stat <- clusGap(t(data.s), FUN = kmeans, nstart = 25, K.max = 7, B = 50)
library(factoextra)
fviz_gap_stat(gap_stat)

## 19. Hierarchical Clustering
# Hitung jarak Euclidean
d <- dist(t(data.s), method = "euclidean")

# Hierarchical Clustering dengan Complete Linkage
hc <- hclust(d, method = "complete")

# Hierarchical Clustering Dendrogram tanpa label
plot(hc, labels = FALSE, main = "Hierarchical Clustering Dendrogram")

# Warnai branches
dend <- as.dendrogram(hc)
dend_colored <- color_branches(dend, k = 3)

# Plot tanpa label
plot(dend_colored, main = "Colored Dendrogram (k=3)", axes = TRUE, leaflab = "none")


# Cut Tree
k <- 3
groups_hc <- cutree(hc, k = k)

# Tabel Hasil Cluster
table(groups_hc, vargrp_clean)

## 20. K-Means Clustering
set.seed(77)
k_means <- kmeans(t(expdtgeoFilt), centers = k)

# Tabel Hasil Cluster
table(k_means$cluster, vargrp_clean)

# Visualisasi Clustering
fviz_cluster(list(data = t(expdtgeoFilt), cluster = k_means$cluster),
             geom = "point", stand = FALSE, main = "K-Means Clustering (PCA)")

## 21. PAM Clustering
pam_result <- pam(t(expdtgeoFilt), k = k)

# Tabel Hasil Cluster
table(pam_result$clustering, vargrp_clean)

# Visualisasi Clustering
fviz_cluster(pam_result, geom = "point", stand = FALSE,
             main = "PAM Clustering (PCA)")

## 22. Biclustering BCBimax
library(biclust)

bcbimax_biclust <- biclust(data.s, method = BCBimax())

# Tampilkan Jumlah Bicluster
bcbimax_biclust

# Heatmap dari Bicluster
drawHeatmap2(data.s, bcbimax_biclust, number = 1)

## 23. Biclustering Plaid Model
set.seed(57)

plaid_biclust <- biclust(data.s, method = BCPlaid(),
                         fit.model = ~m + a + b, background = TRUE, cluster = "b",
                         iter.startup = 10, iter.layer = 50)

# Ringkasan Plaid Model
summary(plaid_biclust)

# Ekstraksi Data Bicluster
bicluster_rows <- which(plaid_biclust@RowxNumber[,1])
bicluster_cols <- which(plaid_biclust@NumberxCol[1,])

bicluster_data <- expdtgeoFilt[bicluster_rows, bicluster_cols]

# Visualisasi Heatmap Bicluster
library(pheatmap)
pheatmap(bicluster_data, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE)

## 24. Export Data untuk Klasifikasi
write.csv(expdtgeoFilt, file = "gene_exp_data.csv", row.names = TRUE)
write.csv(finalres, file = "final_annotation.csv", row.names = FALSE)
getwd()

write.csv(expdtgeoFilt, file = "C:/Users/LENOVO/DOwnloads/gene_exp_data.csv", row.names = TRUE)
write.csv(finalres, file = "C:/Users/LENOVO/DOwnloads/final_annotation.csv", row.names = FALSE)


# 1. Buat faktor grup
group <- as.factor(vargrp_clean)
levels(group)  # harus apocrine, basal, luminal

# 2. Loop untuk setiap kombinasi grup
pairs <- combn(levels(group), 2, simplify = FALSE)

for (p in pairs) {
  
  # Ambil dua grup
  idx <- group %in% p
  group2 <- droplevels(group[idx])
  data2 <- expdtgeoFilt[, idx]
  
  # Hitung p-value (t-test per gen)
  pvalues <- apply(data2, 1, function(x) {
    t.test(x ~ group2)$p.value
  })
  
  # Hitung mean difference
  mean_diff <- apply(data2, 1, function(x) {
    tapply(x, group2, mean)[1] - tapply(x, group2, mean)[2]
  })
  
  # Plot volcano
  plot(mean_diff, -log10(pvalues),
       pch = 20, cex = 0.5,
       main = paste("Volcano Plot:", p[1], "vs", p[2]),
       xlab = "Mean Difference",
       ylab = "-log10(P-value)",
       col = ifelse(pvalues < 0.05, "red", "black"))
  
  abline(h = -log10(0.05), col = "blue", lty = 2)
  
  legend("topright",
         legend = c("Significant", "Not Significant"),
         col = c("red", "black"),
         pch = 20,
         cex = 0.8)
}
