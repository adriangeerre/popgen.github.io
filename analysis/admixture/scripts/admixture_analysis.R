########################
## Admixture analysis ##
########################

# Load libraries
library(tidyverse)
library(SNPRelate)

# Set working directory
setwd("/datos/analisis/admixture")

# Transform VCF to GDS (Run one time!!)
vcf_file <- "chr22.phase3.vcf.gz"
snpgdsVCF2GDS(vcf_file, "chr22.phase3.gds", method="biallelic.only")

# Load data
metadata <- read.csv("igsr_samples_phase3.csv", sep = '\t') # Metadata
variants <- snpgdsOpen("chr22.phase3.gds") # Open gds file (Variants)

# PCA
pca <- snpgdsPCA(variants)
summary(pca)

# Create dataset matching metadata and PCA
eigenvectors <- as.data.frame(pca$eigenvect)
eigenvectors$region <- metadata[match(pca$sample.id, metadata$Sample.name),]$Superpopulation.code
eigenvectors$population <- metadata[match(pca$sample.id, metadata$Sample.name),]$Population.name

# PCA variance
pca_var <- pca$varprop[!is.nan(pca$varprop)] * 100 # It contains NAN because the variance is below the representable limit
pcs <- seq(1, length(pca_var))
plot(pcs, pca_var, type = "b", col = "red", pch = c(16), xlab="PCs", ylab="Variance (%)") # The first two PCs contains the most variance (2.32%)

# Plot of the first two PCs
ggplot(data = eigenvectors, aes(x=V1, y=V2, col=region)) +
  geom_point(size=3, alpha=0.7) +
  labs(x="PC2", y="PC1") +
  theme_bw()

'With the two first PCs we have 5 regions: Africa, America, East Asia, Europe and South Asia. The clusters for East Asia, Europe and South Asia are clearly divided. South Asia and Europe are close but only a few samples intercept. The African cluster is separated from the others with a large distance but a few individuals draw a line approaching Europe. Finally, America has a spread cluster with a triangle shape were the individuals are mainly group around South Asia and Europe. In addition, one of the peaks of the triangle approaches the african cluster.'

# Correlation between PCs
cor(eigenvectors %>% select(V1,V2), method = "pearson")

'The correlation is 5.0e-16. Is a small correlation between the PC1 and PC2'

# GDS to BED to use in ADMIXTURE
snpgdsGDS2BED(variants, "chr22.phase3.gds", sample.id = NULL)

# Q-matrix (Basic plot)
q_matrix <- read.table("chr22.phase3.gds.5.Q")
ord <- q_matrix[order(q_matrix$V1, q_matrix$V2, q_matrix$V3, q_matrix$V4, q_matrix$V5),] # Order min to max value
admixture_plot <- barplot(t(as.matrix(ord)), space=c(0.2), col=rainbow(5), xlab="Individual #", ylab="Ancestry", border=NA, las=2)

# Q-matrix (Tidyverse)
tidy_q <- q_matrix %>% as.data.frame() %>% mutate(sample = 1:nrow(q_matrix)) %>% gather(key = region, values, 1:5)
tidy_q$sample <- as.character(tidy_q$sample)
tidy_q %>% arrange(region) %>% ggplot() + geom_col(aes(x=sample,y=values, fill=region)) + theme_bw()

#######################
## Advanced analysis ##
#######################

# LD prunning

'The SNPs could be related to each other because of LD. This situation could mask the real clustering because not all SNPs have the same relevance in the evolutionary history of the populations. By prunning the SNPs by LD we remove data that clear our visualization.'

spns <- snpgdsLDpruning(variants, ld.threshold = 0.5, num.thread = 4)
snpset.id <- unlist(snpset)