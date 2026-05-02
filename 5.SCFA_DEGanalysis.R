## RNAseq Analysis on RNAseq from microglia isolated 24hr after LPS. 4 weeks SCFA diet
#RNA Seq using edgeR and limma voom
# AVC April 2026
## following: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

library(readr)
library(dplyr)
library(ggplot2)
library(limma)
library(edgeR)
library(tidyverse)


##LOADING IN THE DATA##

path = "/Users/aciernia/Sync/CierniaLabMembers/AnnieCiernia/Experiments/SCFA_Invivo2025/RNAseq_April2026"
setwd(path)


#Counts 1

counts <- read.csv("R696_raw_counts.csv")
head(counts)

tail(counts)

rownames(counts) <- counts$gene_id

#samples
Samples <- colnames(counts)[2:ncol(counts)]
Samples<-as.data.frame(Samples)
head(Samples)


#Counts 2

counts2 <- read.csv("R799_raw_counts_updated.csv")
head(counts2)

tail(counts2)

#remove last 5 rows
counts2 <- head(counts2, -5)

rownames(counts2) <- counts2$gene_id

Samples2 <- colnames(counts2)[2:ncol(counts2)]
Samples2 <- as.data.frame(Samples2)
head(Samples2)


#combine matrixs
counts_master <- merge(counts, counts2, by= c("gene_id"))

rownames(counts_master) <- counts_master$gene_id

counts_master <- counts_master %>% dplyr::select(-gene_id, -gene_name.x,-gene_name.y)
names(counts_master)


#matrix gene names and ensembl IDs (column bind) 
genenames = as.data.frame(rownames(counts_master))
colnames(genenames) <- c("EnsemblID")

############################################################################################
## add gene length
############################################################################################

# BiocManager::install(c("GenomicFeatures","TxDb.Mmusculus.UCSC.mm10.knownGene","org.Mm.eg.db"))
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(GenomicRanges)
library(dplyr)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

ex_by_gene <- exonsBy(txdb, by = "gene")                 # names: Entrez IDs
gene_len_bp <- sapply(reduce(ex_by_gene), function(gr) sum(width(gr)))

len_df <- tibble(ENTREZID = names(gene_len_bp),
                 length_bp = as.integer(gene_len_bp))

# Map Entrez -> Ensembl + Symbol
map_df <- AnnotationDbi::select(org.Mm.eg.db,
                                keys = len_df$ENTREZID,
                                keytype = "ENTREZID",
                                columns = c("ENSEMBL","SYMBOL")) %>%
  distinct(ENTREZID, .keep_all = TRUE)

lengths_df <- len_df %>%
  left_join(map_df, by = "ENTREZID") %>%
  dplyr::select(ENSEMBL, SYMBOL, ENTREZID, length_bp)

lengths_df$length_kb = lengths_df$length_bp/1000

symbol_length <- lengths_df %>% dplyr::select(ENSEMBL,SYMBOL,length_kb ) %>% distinct()

#add to genenames
genenames <- merge(genenames, symbol_length, by.x="EnsemblID",by.y="ENSEMBL")

#reorder to match orginal
idx <- match(rownames(counts_master), genenames$EnsemblID)

genenames_aligned <- genenames[idx, , drop = FALSE]

############################################################################################
## metadata
############################################################################################

#read in metadata 1
metadata <- read.csv("SCFA_rna_seq_metadata1.csv")
head(metadata)

metadata$Sample.ID <- paste(metadata$Cohort_1, metadata$Sample.number, sep=".")

#remove SCFA2 11 > appears to be wrong LPS dose
metadata <- metadata %>% filter(Sample.ID != "SCFA2.11")

#match samples 1
metadata$Sample.ID

#filter counts 1
counts1 <- counts_master[, which(colnames(counts_master) %in% metadata$Sample.ID)]
colnames(counts1)


# Step 1: keep only samples present in counts1
metadata <- metadata[metadata$Sample.ID %in% colnames(counts1), ]

# Step 2: reorder metadata to match counts1 column order
metadata <- metadata[match(colnames(counts1), metadata$Sample.ID), ]

colnames(counts1)
metadata$Sample.ID

#read in metadata 2
metadata2 <- read.csv("2025-12-19_SCFA7.csv")
head(metadata2)

metadata2$Sample.ID <- metadata2$Number

#clean up column names
colnames(counts_master) <- sub("\\.MT.*", "", colnames(counts_master))
colnames(counts_master)


#filter counts 2
counts2 <- counts_master[, which(colnames(counts_master) %in% metadata2$Sample.ID)]
colnames(counts2)


#make metadata match
# Step 1: keep only samples present in counts1
metadata2 <- metadata2[metadata2$Sample.ID %in% colnames(counts2), ]

# Step 2: reorder metadata to match counts1 column order
metadata2 <- metadata2[match(colnames(counts2), metadata2$Sample.ID), ]

colnames(counts2)
metadata2$Sample.ID






#combine metadata
names(metadata)
names(metadata2)

metadata$Diet <- metadata$Condition
metadata$Sac.Date <- metadata$Sac.date

m1 <- metadata %>% dplyr::select(Sample.ID,Cohort,Cohort_1, Sex,Diet,Injection)
m2 <- metadata2 %>% dplyr::select(Sample.ID, Cohort,Cohort_1, Sex,Diet,Injection)


meta_master <- rbind(m1,m2)
meta_master$Sex <- gsub("female", "F", meta_master$Sex)
meta_master$Sex <- gsub("male", "M", meta_master$Sex)

meta_master$Diet <- factor(meta_master$Diet, levels = c("Saline", "Butyrate","Propionate","Acetate"))

meta_master$Injection <- factor(meta_master$Injection, levels = c("PBS","LPS"))

unique(meta_master$Cohort_1)

#order meta_master

meta_master <- meta_master %>% arrange(Sex, Diet, Injection)



#combne counts
counts_final <- cbind(counts1,counts2)
  
  
#match to metadata
counts_final <- counts_final[, meta_master$Sample.ID]

meta_master <- meta_master[match(colnames(counts_final), meta_master$Sample.ID), ]

colnames(counts_final)
meta_master$Sample.ID



write.csv(counts_final,"FinalRawCountsMatrix.csv")

write.csv(meta_master,"Final_Metadata.csv")

#new treatment variable that is LPS+sex

meta_master$condition <- paste(meta_master$Sex, meta_master$Diet, meta_master$Injection, sep="_")

meta_master$condition = factor(meta_master$condition)

table(meta_master$condition)

#Make DEG list object
d0 <- DGEList(counts_final, group = meta_master$condition)
d0$genes = genenames_aligned
dim(d0)

save(d0, counts_final,meta_master, file="RawInputs.Rdata")

write.csv(d0$samples,file="sample_counts.csv")


metadata <- meta_master

############################################################################################
##getting rid of low expressed genes (cutoff for the number of samples that we want above 1,getting rid of things that are less than one group)
############################################################################################

keep <- filterByExpr(
  d0,
  group = meta_master$condition,
  min.prop = 0.5
)

d <- d0[keep, , keep.lib.sizes = FALSE]

dim(d)
#21871    83

cat("Genes before filtering:", nrow(d0$counts), "\n")
cat("Genes after filtering:", sum(keep), "\n")

############################################################################################
#Making a filtering plot
############################################################################################
library(RColorBrewer)
nsamples <- ncol(d0)

colourCount = nsamples
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
fill=getPalette(colourCount)

#plot:
pdf('FilteringCPM_plots.pdf', h=4, w=6)
par(mfrow=c(1,2))

#prefilter:
lcpm <- cpm(d0, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}

#filtered data
#og-CPM of zero threshold (equivalent to a CPM value of 1) used in the filtering ste
lcpm <- cpm(d, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")
dev.off()


############################################################################################
##Plot library size##
############################################################################################
#plot library sizes
pdf('LibrarySizes.pdf',w=30,h=8)
barplot(d$samples$lib.size,names=colnames(d),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()

############################################################################################
####Design model matrix####
############################################################################################


# 2. build your design matrix
design <- model.matrix(~ 0 + meta_master$condition + meta_master$Cohort_1  )  # or whatever your formula is

rownames(design) <- meta_master$Sample.ID

##Run Calculation Normalization
DGE <- calcNormFactors(d, method = "TMM")

# Verify normalization happened
print(head(DGE$samples$norm.factors))

#Compute normalized logCPM   
lcpm_norm <- cpm(DGE, log = TRUE, normalized.lib.sizes = TRUE)

#voom
pdf('Voom.pdf',w=6,h=4)
v <- voom(DGE, design, plot = TRUE)
dev.off()


############################################################################################
#MDS Plots no corrections
############################################################################################

library(RColorBrewer)
colors_16 <- colorRampPalette(brewer.pal(8, "Set3"))(16)


pdf(file = "MDSplots_NoCorrection.pdf", width = 16, height = 16)  # Wider PDF

par(mfrow = c(4, 2),           # 4 rows, 2 columns of plots
    mar = c(5, 5, 4, 8),       # Extra space on the right
    oma = c(1, 1, 1, 4),       # Outer margin space
    xpd = TRUE)               # Allow plotting outside the box

### Panel 1: LPS x Sex
col.cell <- colors_16[metadata$condition]
plotMDS(DGE, col = col.cell, dim.plot = c(1, 2))
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$condition), cex = 0.8, bty = "n")
title("Diet x LPS x Sex")

### Panel 2: LPS Treatment

metadata$col.cell <- c("blue","green")[metadata$Injection]
plotMDS(DGE, col = metadata$col.cell, dim.plot = c(1, 2))
legend("topleft", fill = unique(metadata$col.cell),
       legend = unique(metadata$Injection), cex = 0.8, bty = "n")
title("LPS Treatment")

### Panel 3: Sex (dim 1 vs 2)
metadata$Sex <- as.factor(metadata$Sex)
col.cell <- c("black","red")[metadata$Sex]
plotMDS(DGE, col = col.cell, dim.plot = c(1, 2), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Sex), cex = 0.8, bty = "n")
title("Sex (Dim 1 vs 2)")

### Panel 4: Sex (dim 3 vs 4)
plotMDS(DGE, col = col.cell, dim.plot = c(3, 4), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Sex), cex = 0.8, bty = "n")
title("Sex (Dim 3 vs 4)")

### Panel 5: Pups per litter
col.cell <- c("black","red","purple","darkgreen","orange")[metadata$Diet]
plotMDS(DGE, col = col.cell, dim.plot = c(1, 2), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Diet), cex = 0.8, bty = "n")
title("Diet (Dim 2 vs 3")

### Panel 6: Pups per litter
col.cell <- c("black","red","purple","darkgreen","orange")[metadata$Diet]
plotMDS(DGE, col = col.cell, dim.plot = c(3, 4), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Diet), cex = 0.8, bty = "n")
title("Diet (Dim 3 vs 4")

### Panel 7: Batch
metadata$Cohort_1 <- as.factor(metadata$Cohort_1)
col.cell <- colors_16[metadata$Cohort_1]
plotMDS(DGE, col = col.cell, dim.plot = c(1, 2), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Cohort_1), cex = 0.8, bty = "n")
title("Collection Date (Dim 1 vs 2)")


### Panel 7: Batch

col.cell <- colors_16[metadata$Cohort_1]
plotMDS(DGE, col = col.cell, dim.plot = c(3, 4), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Cohort_1), cex = 0.8, bty = "n")
title("Collection Date (Dim 3 vs 4)")


dev.off()


############################################################################################
#MDS Plots with batch corrections
############################################################################################

DGE1<-removeBatchEffect(DGE$counts, group = metadata$condition, batch = metadata$Cohort_1)



library(RColorBrewer)
colors_16 <- colorRampPalette(brewer.pal(8, "Set3"))(16)


pdf(file = "MDSplots_Correction.pdf", width = 16, height = 16)  # Wider PDF

par(mfrow = c(4, 2),           # 4 rows, 2 columns of plots
    mar = c(5, 5, 4, 8),       # Extra space on the right
    oma = c(1, 1, 1, 4),       # Outer margin space
    xpd = TRUE)               # Allow plotting outside the box

### Panel 1: LPS x Sex
col.cell <- colors_16[metadata$condition]
plotMDS(DGE1, col = col.cell, dim.plot = c(1, 2))
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$condition), cex = 0.8, bty = "n")
title("Diet x LPS x Sex")

### Panel 2: LPS Treatment

metadata$col.cell <- c("blue","green")[metadata$Injection]
plotMDS(DGE1, col = metadata$col.cell, dim.plot = c(1, 2))
legend("topleft", fill = unique(metadata$col.cell),
       legend = unique(metadata$Injection), cex = 0.8, bty = "n")
title("LPS Treatment")

### Panel 3: Sex (dim 1 vs 2)
metadata$Sex <- as.factor(metadata$Sex)
col.cell <- c("black","red")[metadata$Sex]
plotMDS(DGE1, col = col.cell, dim.plot = c(1, 2), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Sex), cex = 0.8, bty = "n")
title("Sex (Dim 1 vs 2)")

### Panel 4: Sex (dim 3 vs 4)
plotMDS(DGE1, col = col.cell, dim.plot = c(3, 4), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Sex), cex = 0.8, bty = "n")
title("Sex (Dim 3 vs 4)")

### Panel 5: Pups per litter
col.cell <- c("black","red","purple","darkgreen","orange")[metadata$Diet]
plotMDS(DGE1, col = col.cell, dim.plot = c(1, 2), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Diet), cex = 0.8, bty = "n")
title("Diet (Dim 2 vs 3")

### Panel 6: Pups per litter
col.cell <- c("black","red","purple","darkgreen","orange")[metadata$Diet]
plotMDS(DGE1, col = col.cell, dim.plot = c(3, 4), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Diet), cex = 0.8, bty = "n")
title("Diet (Dim 3 vs 4")

### Panel 7: Batch
metadata$Cohort_1 <- as.factor(metadata$Cohort_1)
col.cell <- colors_16[metadata$Cohort_1]
plotMDS(DGE1, col = col.cell, dim.plot = c(1, 2), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Cohort_1), cex = 0.8, bty = "n")
title("Collection Date (Dim 1 vs 2)")


### Panel 7: Batch

col.cell <- colors_16[metadata$Cohort_1]
plotMDS(DGE1, col = col.cell, dim.plot = c(3, 4), pch = 16)
legend("topleft", fill = unique(col.cell),
       legend = levels(metadata$Cohort_1), cex = 0.8, bty = "n")
title("Collection Date (Dim 3 vs 4)")


dev.off()






# 
# ##########################################
# ### Model with Cohort as a covariate
# ##########################################
# 
# 2. build your design matrix
design <- model.matrix(~ 0 + meta_master$condition + meta_master$Cohort_1  ) 

rownames(design) <- meta_master$Sample.ID

##Run Calculation Normalization
DGE <- calcNormFactors(d, method = "TMM")

# Verify normalization happened
print(head(DGE$samples$norm.factors))


#voom
pdf('Voom.pdf',w=6,h=4)
v <- voom(DGE, design, plot = TRUE)
dev.off()


##########################################
### Fit the linear model using block + correlation
##########################################

fit <- lmFit(v, design) 
        

colnames(design) <- gsub("meta_master\\$condition","",colnames(design))
colnames(design) <- gsub("meta_master\\$","",colnames(design))
############################################################################################
### make contrast matrix
############################################################################################
contr.matrix <- makeContrasts(
  
  # F: Diet x LPS
  F_Acetate_LPSvsPBS    = F_Acetate_LPS    - F_Acetate_PBS,
  F_Butyrate_LPSvsPBS   = F_Butyrate_LPS   - F_Butyrate_PBS,
  F_Propionate_LPSvsPBS = F_Propionate_LPS - F_Propionate_PBS,
  F_Saline_LPSvsPBS     = F_Saline_LPS     - F_Saline_PBS,
  
  # M: Diet x LPS 
  M_Acetate_LPSvsPBS      = M_Acetate_LPS      - M_Acetate_PBS,
  M_Butyrate_LPSvsPBS     = M_Butyrate_LPS     - M_Butyrate_PBS,
  M_Propionate_LPSvsPBS   = M_Propionate_LPS   - M_Propionate_PBS,
  M_Saline_LPSvsPBS       = M_Saline_LPS       - M_Saline_PBS,
  
  # Sex-collapsed: Diet x LPS (average of M + F)
  BothSex_Acetate_LPSvsPBS    = 0.5*(F_Acetate_LPS    - F_Acetate_PBS)    +
    0.5*(M_Acetate_LPS      - M_Acetate_PBS),
  BothSex_Butyrate_LPSvsPBS   = 0.5*(F_Butyrate_LPS   - F_Butyrate_PBS)   +
    0.5*(M_Butyrate_LPS     - M_Butyrate_PBS),
  BothSex_Propionate_LPSvsPBS = 0.5*(F_Propionate_LPS - F_Propionate_PBS) +
    0.5*(M_Propionate_LPS   - M_Propionate_PBS),
  BothSex_Saline_LPSvsPBS     = 0.5*(F_Saline_LPS     - F_Saline_PBS)     +
    0.5*(M_Saline_LPS       - M_Saline_PBS),
  
  # F: Diet within LPS 
  F_LPS_AcetatevsSaline    = F_Acetate_LPS    - F_Saline_LPS, 
  F_LPS_ButyratevsSaline   = F_Butyrate_LPS   - F_Saline_LPS, 
  F_LPS_PropionatevsSaline = F_Propionate_LPS - F_Saline_LPS, 
  
  # F: Diet within PBS
  F_PBS_AcetatevsSaline    = F_Acetate_PBS    - F_Saline_PBS, 
  F_PBS_ButyratevsSaline   = F_Butyrate_PBS   - F_Saline_PBS, 
  F_PBS_PropionatevsSaline = F_Propionate_PBS - F_Saline_PBS, 
  
  # M: Diet within LPS
  M_LPS_AcetatevsSaline    = M_Acetate_LPS    - M_Saline_LPS, 
  M_LPS_ButyratevsSaline   = M_Butyrate_LPS   - M_Saline_LPS, 
  M_LPS_PropionatevsSaline = M_Propionate_LPS - M_Saline_LPS, 
  
  # M: Diet within PBS
  M_PBS_AcetatevsSaline    = M_Acetate_PBS    - M_Saline_PBS, 
  M_PBS_ButyratevsSaline   = M_Butyrate_PBS   - M_Saline_PBS, 
  M_PBS_PropionatevsSaline = M_Propionate_PBS - M_Saline_PBS, 
  
  # Sex-collapsed: Diet within LPS
  BothSex_LPS_AcetatevsSaline    = 0.5*(F_Acetate_LPS    - F_Saline_LPS)    +
    0.5*(M_Acetate_LPS      - M_Saline_LPS),
  BothSex_LPS_ButyratevsSaline   = 0.5*(F_Butyrate_LPS   - F_Saline_LPS)   +
    0.5*(M_Butyrate_LPS     - M_Saline_LPS),
  BothSex_LPS_PropionatevsSaline = 0.5*(F_Propionate_LPS - F_Saline_LPS) +
    0.5*(M_Propionate_LPS   - M_Saline_LPS),
  
  # Sex-collapsed: Diet within PBS
  BothSex_PBS_AcetatevsSaline    = 0.5*(F_Acetate_PBS    - F_Saline_PBS)    +
    0.5*(M_Acetate_PBS      - M_Saline_PBS),
  BothSex_PBS_ButyratevsSaline   = 0.5*(F_Butyrate_PBS   - F_Saline_PBS)   +
    0.5*(M_Butyrate_PBS     - M_Saline_PBS),
  BothSex_PBS_PropionatevsSaline = 0.5*(F_Propionate_PBS - F_Saline_PBS) +
    0.5*(M_Propionate_PBS   - M_Saline_PBS),
  
  levels = colnames(design)
)

contr.matrix

#apply contrast matrix

fit2 <- contrasts.fit(fit, contr.matrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

pdf('PlotSA_VoomTrend.pdf',w=6,h=4)
plotSA(fit2, main="Final model: Mean variance trend")
dev.off()

#run tests
dt <- decideTests(fit2)
summary(dt)
write.csv(summary(dt), file = "DEGcounts.csv")

############################################################################################
###Normalize Log counts###
############################################################################################


# Order samples by condition
ord <- order(metadata$condition)

# Reorder condition vector
cond_ord <- metadata$condition[ord]

# Define colors per condition
cond_levels <- levels(cond_ord)
cond_cols <- setNames(rainbow(length(cond_levels)), cond_levels)
sample_cols <- cond_cols[cond_ord]



pdf("NormalizationPlot.pdf", width = 20, height = 24)

# Two stacked panels
layout(matrix(c(1, 2), nrow = 2), heights = c(1, 1.2))

##########################################################
### Panel C: Raw (unnormalized) log2 CPM
##########################################################

par(mar = c(10, 6, 6, 2) + 0.1)

logcounts_raw <- cpm(d,
                     log = TRUE,
                     normalized.lib.sizes = FALSE)

bp_raw <- boxplot(logcounts_raw[, ord],
                  xlab = "",
                  ylab = "Log2 Counts Per Million (raw)",
                  las = 2,
                  cex.names = 0.75,
                  col = sample_cols,
                  border = sample_cols,
                  medlwd = 2,        # thicker sample median bars
                  outline = FALSE)

# Overall median (global)
abline(h = median(logcounts_raw), col = "blue", lty = 2, lwd = 2)

# Explicit sample-median markers
points(seq_along(bp_raw$stats[3, ]),
       bp_raw$stats[3, ],
       pch = 16,
       cex = 0.6)

mtext("C. Raw Log2CPM (unnormalized)",
      side = 3, line = 3, adj = 0, cex = 2)

##########################################################
### Panel D: Post-voom log2 CPM
##########################################################

par(mar = c(10, 6, 6, 2) + 0.1)

bp_voom <- boxplot(v$E[, ord],
                   xlab = "",
                   ylab = "Log2 Counts Per Million (post-voom)",
                   las = 2,
                   cex.names = 0.75,
                   col = sample_cols,
                   border = sample_cols,
                   medlwd = 2,        # thicker sample median bars
                   outline = FALSE)

# Overall median (global)
abline(h = median(v$E), col = "blue", lty = 2, lwd = 2)

# Explicit sample-median markers
points(seq_along(bp_voom$stats[3, ]),
       bp_voom$stats[3, ],
       pch = 16,
       cex = 0.6)

mtext("D. Post-voom Log2CPM (final expression values)",
      side = 3, line = 3, adj = 0, cex = 2)

##########################################################
### Legend
##########################################################

legend("topright",
       legend = cond_levels,
       fill = cond_cols,
       border = cond_cols,
       cex = 1.2,
       bty = "n")

dev.off()





############################################################################################
###GLIMMA interactive plot building###
############################################################################################
#BiocManager::install("Glimma")
#http://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/Glimma.pdf
library(Glimma)
#writes out html file: ** create MDS glimma from DGE1***
glMDSPlot(DGE1, 
          groups = metadata$condition,
          labels = metadata$Sample.ID,
          launch = TRUE)



#batch corrected log2CPM
#nc_norm <-removeBatchEffect(v$E, group = metadata$treatment, batch = metadata$date)

#glMDSPlot(nc_norm, groups=group)

for (COEF in 1:ncol(fit2)) {
  glMDPlot(fit2, counts=v$E,transform=FALSE,anno=v$genes,
           coef=COEF, status=dt, main=colnames(fit2)[COEF],
           groups=metadata$condition, folder="glimma_results", 
           launch=FALSE, html = paste("MD-Plot",colnames(contr.matrix)[COEF]))}


############################################################################################
#####Average CPM for each condition#####
############################################################################################

###get log2CPM counts from voom and put in dataframe:
library(plotrix)
#average log2 CPM and sem
countdf <- as.data.frame(v$E)
countdf$GeneID <- rownames(v$E)

#add gene names
DF <- merge(countdf,genenames_aligned, by.x ="GeneID",by.y="EnsemblID")
#write as csv
write.csv(DF,file="log2CPMvalues.csv")

#summarize 
countdf2 <- DF %>% group_by(GeneID,SYMBOL) %>% tidyr::gather(sample,log2CPM, 2:84) 
countdf2 <- as.data.frame(countdf2)

countdf3 <-merge(countdf2, metadata,by.x="sample",by.y="Sample.ID",all=T)


#save
write.csv(GeneSummary, file = "AverageLog2CPM.csv")

# Wide matrix: rows = genes, columns = samples, values = log2CPM
log2cpm_mat <- countdf3 %>%
  dplyr::select(GeneID, sample, log2CPM) %>%   # drop length_kb here
  distinct() %>%                         # in case of duplicates
  pivot_wider(
    names_from  = sample,
    values_from = log2CPM
  )

write.csv(log2cpm_mat, "Log2CPM_Data.csv")



#rpkm
rpkm_mat   <- edgeR::rpkm(DGE, gene.length = v$genes$length_kb, log = FALSE, normalized.lib.sizes = TRUE)
rpkm_mat   <- as.data.frame(rpkm_mat)
rpkm_mat$Ensembl <- rownames(rpkm_mat)

rpkm_wide <- rpkm_mat %>% group_by(Ensembl) %>% gather(sample, RPKM,1:(ncol(rpkm_mat)-1))

#add genesymboles
countdf3$Ensembl <- countdf3$GeneID
DF <- merge(rpkm_wide, countdf3, by=c("sample","Ensembl"))

write.csv(DF,"RPKM_Data.csv")



############################################################################################
##### GET DEGS #####
############################################################################################

library(calibrate)
library(dplyr)

####Make contrasts####
comparisons=(coef(fit2))
comparisons=colnames(comparisons)
comp_out <- as.data.frame(rownames(v$E))
names(comp_out) <- c("GeneID")
nrowkeep <- nrow(comp_out)

SumTableOut <- NULL

for(i in 1:length(comparisons)){
  #comparison name
  comp=comparisons[i]
  print(comp)
  #make comparisons 
  
  topT=topTreat(fit2,coef=i,number=nrowkeep,adjust.method="BH")
  print(nrow(topT[(topT$adj.P.Val<0.05),]))

  #LogFC values:https://support.bioconductor.org/p/82478/
  topT$direction <- c("none")
  topT$direction[which(topT$logFC > 0)] = c("Increase")
  topT$direction[which(topT$logFC < 0)] = c("Decrease")
  
  topT$significance <- c("nonDE")
  topT$significance[which(topT$adj.P.Val <0.05)] <- c("DE")

  #summary counts table based on Ensemble Gene ID counts:
  SumTable <- table(topT$significance,topT$direction)
  SumTable <- as.data.frame(SumTable)
  SumTable$comparison <- paste(comp)
  SumTableOut <- rbind(SumTable,SumTableOut)


  
  #gene gene names and expression levels
  topT2 <- topT
  topT2$comparison <- paste(comp)
  write.csv(topT2,file = paste(comp,"_DEgenes.csv"))
  
  #get master file:
  colnames(topT)[3:ncol(topT)] <- paste(colnames(topT)[3:ncol(topT)],comp)
  comp_out <- merge(comp_out,topT, by.x = "GeneID" , by.y= "EnsemblID")
  
  #data for plot with gene names:
  genenames <- topT2 %>% dplyr::select(adj.P.Val,logFC,SYMBOL) %>% distinct()
  
  #names for plots
  plotname <- gsub("\\."," ",comp)
  plotname <- gsub("vs"," vs ",plotname)
  
  #volcano plot
  pdf(file = paste(comp,"_Volcano.pdf", sep=""), wi = 9, he = 6, useDingbats=F)
  
  with(genenames, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main=paste(plotname,"\nVolcano plot", sep=" "), ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2 
  with(subset(genenames, logFC < -2 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(genenames, logFC > 2 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-2,2), col = "black", lty = 2, lwd = 1)
  
  #Label points with the textxy function from the calibrate plot
  library(calibrate)
  with(subset(genenames, adj.P.Val<0.05 & abs(logFC)>2), textxy(logFC, -log10(adj.P.Val), labs=SYMBOL, cex=.5))
  
  dev.off()
  
}

write.csv(SumTableOut,"SummaryTableDEgenes.csv")
#master outfile to get log2CPM values


#clean comp_out
# remove all genesymbol.x and genesymbol.y columns
comp_out_clean <- comp_out%>%
  dplyr::select(
    -dplyr::matches("^SYMBOL\\.x$"),
    -dplyr::matches("^SYMBOL\\.y$"),
    -dplyr::matches("^length_kb"),
    -dplyr::matches("^AveExp")
  )



mout <- merge(comp_out_clean, log2cpm_mat, by="GeneID")

write.csv(mout,"AllDEG_AllConditions_log2CPM.csv")

mout2 <- merge(comp_out_clean, rpkm_mat, by.x="GeneID", by.y="Ensembl")

write.csv(mout2,"AllDEG_AllConditions_RPKM.csv")



###########################################################
# Compare BothSex to Female and Male for same contrasts
# (per diet, LPS > PBS and LPS < PBS)
###########################################################
library(UpSetR)
## UP-REGULATED (== 1)
dt2 <- as.data.frame(dt)
# Acetate
LPS_up_Acetate_SexCompare <- list(
  'Female_Acetate_LPS>PBS'   = rownames(dt2[dt2$F_Acetate_LPSvsPBS    == 1,]),
  'male_Acetate_LPS>PBS'     = rownames(dt2[dt2$M_Acetate_LPSvsPBS      == 1,]),
  'BothSex_Acetate_LPS>PBS'  = rownames(dt2[dt2$BothSex_Acetate_LPSvsPBS   == 1,])
)

pdf('Upset_SexCompare_Acetate_LPSDEGsup.pdf', w = 10, h = 6)
UpSetR::upset(
  fromList(LPS_up_Acetate_SexCompare),
  nintersects = NA,
  nsets       = length(LPS_up_Acetate_SexCompare),
  text.scale  = 2,
  keep.order  = TRUE
)
dev.off()

# Propionate
LPS_up_Propionate_SexCompare <- list(
  'Female_Propionate_LPS>PBS'   = rownames(dt2[dt2$F_Propionate_LPSvsPBS    == 1,]),
  'male_Propionate_LPS>PBS'     = rownames(dt2[dt2$M_Propionate_LPSvsPBS      == 1,]),
  'BothSex_Propionate_LPS>PBS'  = rownames(dt2[dt2$BothSex_Propionate_LPSvsPBS   == 1,])
)

pdf('Upset_SexCompare_Propionate_LPSDEGsup.pdf', w = 10, h = 6)
UpSetR::upset(
  fromList(LPS_up_Propionate_SexCompare),
  nintersects = NA,
  nsets       = length(LPS_up_Propionate_SexCompare),
  text.scale  = 2,
  keep.order  = TRUE
)
dev.off()

# Butyrate
LPS_up_Butyrate_SexCompare <- list(
  'Female_Butyrate_LPS>PBS'   = rownames(dt2[dt2$F_Butyrate_LPSvsPBS    == 1,]),
  'male_Butyrate_LPS>PBS'     = rownames(dt2[dt2$M_Butyrate_LPSvsPBS      == 1,]),
  'BothSex_Butyrate_LPS>PBS'  = rownames(dt2[dt2$BothSex_Butyrate_LPSvsPBS   == 1,])
)

pdf('Upset_SexCompare_Butyrate_LPSDEGsup.pdf', w = 10, h = 6)
UpSetR::upset(
  fromList(LPS_up_Butyrate_SexCompare),
  nintersects = NA,
  nsets       = length(LPS_up_Butyrate_SexCompare),
  text.scale  = 2,
  keep.order  = TRUE
)
dev.off()

# Saline
LPS_up_Saline_SexCompare <- list(
  'Female_Saline_LPS>PBS'   = rownames(dt2[dt2$F_Saline_LPSvsPBS    == 1,]),
  'male_Saline_LPS>PBS'     = rownames(dt2[dt2$M_Saline_LPSvsPBS      == 1,]),
  'BothSex_Saline_LPS>PBS'  = rownames(dt2[dt2$BothSex_Saline_LPSvsPBS   == 1,])
)

pdf('Upset_SexCompare_Saline_LPSDEGsup.pdf', w = 10, h = 6)
UpSetR::upset(
  fromList(LPS_up_Saline_SexCompare),
  nintersects = NA,
  nsets       = length(LPS_up_Saline_SexCompare),
  text.scale  = 2,
  keep.order  = TRUE
)
dev.off()


## DOWN-REGULATED (== -1)

# Acetate
LPS_down_Acetate_SexCompare <- list(
  'Female_Acetate_LPS<PBS'   = rownames(dt2[dt2$F_Acetate_LPSvsPBS    == -1,]),
  'male_Acetate_LPS<PBS'     = rownames(dt2[dt2$M_Acetate_LPSvsPBS      == -1,]),
  'BothSex_Acetate_LPS<PBS'  = rownames(dt2[dt2$BothSex_Acetate_LPSvsPBS   == -1,])
)

pdf('Upset_SexCompare_Acetate_LPSDEGsdown.pdf', w = 10, h = 6)
UpSetR::upset(
  fromList(LPS_down_Acetate_SexCompare),
  nintersects = NA,
  nsets       = length(LPS_down_Acetate_SexCompare),
  text.scale  = 2,
  keep.order  = TRUE
)
dev.off()

# Propionate
LPS_down_Propionate_SexCompare <- list(
  'Female_Propionate_LPS<PBS'   = rownames(dt2[dt2$F_Propionate_LPSvsPBS    == -1,]),
  'male_Propionate_LPS<PBS'     = rownames(dt2[dt2$M_Propionate_LPSvsPBS      == -1,]),
  'BothSex_Propionate_LPS<PBS'  = rownames(dt2[dt2$BothSex_Propionate_LPSvsPBS   == -1,])
)

pdf('Upset_SexCompare_Propionate_LPSDEGsdown.pdf', w = 10, h = 6)
UpSetR::upset(
  fromList(LPS_down_Propionate_SexCompare),
  nintersects = NA,
  nsets       = length(LPS_down_Propionate_SexCompare),
  text.scale  = 2,
  keep.order  = TRUE
)
dev.off()

# Butyrate
LPS_down_Butyrate_SexCompare <- list(
  'Female_Butyrate_LPS<PBS'   = rownames(dt2[dt2$F_Butyrate_LPSvsPBS    == -1,]),
  'male_Butyrate_LPS<PBS'     = rownames(dt2[dt2$M_Butyrate_LPSvsPBS      == -1,]),
  'BothSex_Butyrate_LPS<PBS'  = rownames(dt2[dt2$BothSex_Butyrate_LPSvsPBS   == -1,])
)

pdf('Upset_SexCompare_Butyrate_LPSDEGsdown.pdf', w = 10, h = 6)
UpSetR::upset(
  fromList(LPS_down_Butyrate_SexCompare),
  nintersects = NA,
  nsets       = length(LPS_down_Butyrate_SexCompare),
  text.scale  = 2,
  keep.order  = TRUE
)
dev.off()

# Saline
LPS_down_Saline_SexCompare <- list(
  'Female_Saline_LPS<PBS'   = rownames(dt2[dt2$F_Saline_LPSvsPBS    == -1,]),
  'male_Saline_LPS<PBS'     = rownames(dt2[dt2$M_Saline_LPSvsPBS      == -1,]),
  'BothSex_Saline_LPS<PBS'  = rownames(dt2[dt2$BothSex_Saline_LPSvsPBS   == -1,])
)

pdf('Upset_SexCompare_Saline_LPSDEGsdown.pdf', w = 10, h = 6)
UpSetR::upset(
  fromList(LPS_down_Saline_SexCompare),
  nintersects = NA,
  nsets       = length(LPS_down_Saline_SexCompare),
  text.scale  = 2,
  keep.order  = TRUE
)
dev.off()





############################################################################################
#upset plots to compare gene lists for males vs females
############################################################################################

library("VennDiagram")


#make list for Females 
LPS_up_Female <- list( 
  'Female_Acetate_LPS>PBS' = rownames(dt2[dt2$F_Acetate_LPSvsPBS == 1,]),
  'Female_Proprionate_LPS>PBS' = rownames(dt2[dt2$F_Propionate_LPSvsPBS == 1,]),
  'Female_Butyrate_LPS>PBS' = rownames(dt2[dt2$F_Butyrate_LPSvsPBS == 1,]),
  'Female_Saline_LPS>PBS' = rownames(dt2[dt2$F_Saline_LPSvsPBS == 1,]))

LPS_down_Female <- list( 
  'Female_Acetate_LPS<PBS' = rownames(dt2[dt2$F_Acetate_LPSvsPBS == -1,]),
   'Female_Propionate_LPS<PBS' = rownames(dt2[dt2$F_Propionate_LPSvsPBS == -1,]),
 'Female_Butyrate_LPS<PBS' = rownames(dt2[dt2$F_Butyrate_LPSvsPBS == -1,]),
'Female_Saline_LPS<PBS' = rownames(dt2[dt2$F_Saline_LPSvsPBS == -1,]))



#plot
library(UpSetR)


#barcolor <- c("orange","orange","purple","purple","blue","blue")

pdf('Upset_Female_LPSDEGsup.pdf',w=12,h=8)
UpSetR::upset(fromList(LPS_up_Female), nintersects =NA,nsets = 8, 
              text.scale = 2, 
              #sets.bar.color = barcolor,
              keep.order = T)
#plot(ven, type = "upset", nintersects=NA, keep.order = T, order.by = "degree")
dev.off()

pdf('Upset_Female_LPSDEGsdown.pdf',w=12,h=8)
UpSetR::upset(fromList(LPS_down_Female), nintersects =NA,nsets = 8, 
              text.scale = 2, 
              #sets.bar.color = barcolor,
              keep.order = T)
#plot(ven, type = "upset", nintersects=NA, keep.order = T, order.by = "degree")
dev.off()



#make list for males 
LPS_up_male <- list( 
  'male_Acetate_LPS>PBS' = rownames(dt2[dt2$M_Acetate_LPSvsPBS == 1,]),
  'male_Proprionate_LPS>PBS' = rownames(dt2[dt2$M_Propionate_LPSvsPBS == 1,]),
  'male_Butyrate_LPS>PBS' = rownames(dt2[dt2$M_Butyrate_LPSvsPBS == 1,]),
  'male_Saline_LPS>PBS' = rownames(dt2[dt2$M_Saline_LPSvsPBS == 1,]))

LPS_down_male <- list( 
  'male_Acetate_LPS<PBS' = rownames(dt2[dt2$M_Acetate_LPSvsPBS == -1,]),
  'male_Propionate_LPS<PBS' = rownames(dt2[dt2$M_Propionate_LPSvsPBS == -1,]),
  'male_Butyrate_LPS<PBS' = rownames(dt2[dt2$M_Butyrate_LPSvsPBS == -1,]),
  'male_Saline_LPS<PBS' = rownames(dt2[dt2$M_Saline_LPSvsPBS == -1,]))


pdf('Upset_male_LPSDEGsup.pdf',w=12,h=8)
UpSetR::upset(fromList(LPS_up_male), nintersects =NA,nsets = 8, 
              text.scale = 2, 
              #sets.bar.color = barcolor,
              keep.order = T)
#plot(ven, type = "upset", nintersects=NA, keep.order = T, order.by = "degree")
dev.off()

pdf('Upset_male_LPSDEGsdown.pdf',w=12,h=8)
UpSetR::upset(fromList(LPS_down_male), nintersects =NA,nsets = 8, 
              text.scale = 2, 
              #sets.bar.color = barcolor,
              keep.order = T)
#plot(ven, type = "upset", nintersects=NA, keep.order = T, order.by = "degree")
dev.off()



  ################################################################################################
  ##eular plots
  ################################################################################################

  # ─────────────────────────────────────────────────────────────
  # 0. Load required libraries (install once if needed)
  # ─────────────────────────────────────────────────────────────
  # install.packages(c("eulerr", "gridExtra", "grid"))
  library(eulerr)
  library(gridExtra)
  library(grid)
  
  # ─────────────────────────────────────────────────────────────
  # 1. Define set colors
  # ─────────────────────────────────────────────────────────────
  sex_colors <- c(
    Male = "#4B9CD3",     # UBC Blue
    Female = "#F59FB1"    # Soft pink
  )
  
  # ─────────────────────────────────────────────────────────────
  # 2. Build gene lists: up- and downregulated
  # ─────────────────────────────────────────────────────────────
  LPS_up_combined <- list(
    'Male_Acetate_LPS>PBS'     = rownames(dt2[dt2$M_Acetate_LPSvsPBS     ==  1, ]),
    'Female_Acetate_LPS>PBS'   = rownames(dt2[dt2$F_Acetate_LPSvsPBS   ==  1, ]),
    'Male_Proprionate_LPS>PBS' = rownames(dt2[dt2$M_Propionate_LPSvsPBS  ==  1, ]),
    'Female_Proprionate_LPS>PBS' = rownames(dt2[dt2$F_Propionate_LPSvsPBS == 1, ]),
    'Male_Butyrate_LPS>PBS'    = rownames(dt2[dt2$M_Butyrate_LPSvsPBS    ==  1, ]),
    'Female_Butyrate_LPS>PBS'  = rownames(dt2[dt2$F_Butyrate_LPSvsPBS  ==  1, ]),
    'Male_Saline_LPS>PBS'      = rownames(dt2[dt2$M_Saline_LPSvsPBS      ==  1, ]),
    'Female_Saline_LPS>PBS'    = rownames(dt2[dt2$F_Saline_LPSvsPBS    ==  1, ])
  )
  
  LPS_down_combined <- list(
    'Male_Acetate_LPS<PBS'     = rownames(dt2[dt2$M_Acetate_LPSvsPBS     == -1, ]),
    'Female_Acetate_LPS<PBS'   = rownames(dt2[dt2$F_Acetate_LPSvsPBS   == -1, ]),
    'Male_Proprionate_LPS<PBS' = rownames(dt2[dt2$M_Propionate_LPSvsPBS  == -1, ]),
    'Female_Proprionate_LPS<PBS' = rownames(dt2[dt2$F_Propionate_LPSvsPBS == -1, ]),
    'Male_Butyrate_LPS<PBS'    = rownames(dt2[dt2$M_Butyrate_LPSvsPBS    == -1, ]),
    'Female_Butyrate_LPS<PBS'  = rownames(dt2[dt2$F_Butyrate_LPSvsPBS  == -1, ]),
    'Male_Saline_LPS<PBS'      = rownames(dt2[dt2$M_Saline_LPSvsPBS      == -1, ]),
    'Female_Saline_LPS<PBS'    = rownames(dt2[dt2$F_Saline_LPSvsPBS    == -1, ])
  )
  
  # ─────────────────────────────────────────────────────────────
  # 3. Define function to generate labeled Euler plots
  # ─────────────────────────────────────────────────────────────
  make_euler_grob <- function(male_vec, female_vec, title){
    plot(
      euler(list(Male = male_vec, Female = female_vec)),
      fills = list(fill = sex_colors, alpha = 0.6),
      quantities = TRUE,
      main = title
    )
  }
  
  # ─────────────────────────────────────────────────────────────
  # 4. Generate plots for LPS-upregulated genes (A–D)
  # ─────────────────────────────────────────────────────────────
  g_up_A <- make_euler_grob(LPS_up_combined[["Male_Acetate_LPS>PBS"]],
                            LPS_up_combined[["Female_Acetate_LPS>PBS"]],
                            "Acetate  (LPS>VEH)")
  g_up_B <- make_euler_grob(LPS_up_combined[["Male_Proprionate_LPS>PBS"]],
                            LPS_up_combined[["Female_Proprionate_LPS>PBS"]],
                            "Propionate  (LPS>VEH)")
  g_up_C <- make_euler_grob(LPS_up_combined[["Male_Butyrate_LPS>PBS"]],
                            LPS_up_combined[["Female_Butyrate_LPS>PBS"]],
                            "Butyrate  (LPS>VEH)")
  g_up_D <- make_euler_grob(LPS_up_combined[["Male_Saline_LPS>PBS"]],
                            LPS_up_combined[["Female_Saline_LPS>PBS"]],
                            "Saline  (LPS>VEH)")
  
  # ─────────────────────────────────────────────────────────────
  # 5. Generate plots for LPS-downregulated genes (E–H)
  # ─────────────────────────────────────────────────────────────
  g_dn_E <- make_euler_grob(LPS_down_combined[["Male_Acetate_LPS<PBS"]],
                            LPS_down_combined[["Female_Acetate_LPS<PBS"]],
                            "Acetate  (LPS<VEH)")
  g_dn_F <- make_euler_grob(LPS_down_combined[["Male_Proprionate_LPS<PBS"]],
                            LPS_down_combined[["Female_Proprionate_LPS<PBS"]],
                            "Propionate  (LPS<VEH)")
  g_dn_G <- make_euler_grob(LPS_down_combined[["Male_Butyrate_LPS<PBS"]],
                            LPS_down_combined[["Female_Butyrate_LPS<PBS"]],
                            "Butyrate  (LPS<VEH)")
  g_dn_H <- make_euler_grob(LPS_down_combined[["Male_Saline_LPS<PBS"]],
                            LPS_down_combined[["Female_Saline_LPS<PBS"]],
                            "Saline  (LPS<VEH)")
  
  # ─────────────────────────────────────────────────────────────
  # 6. Function to add panel labels
  # ─────────────────────────────────────────────────────────────
  label_panel <- function(p, lab){
    arrangeGrob(
      p,
      top = textGrob(lab, x = unit(0, "npc"), y = unit(1, "npc"),
                     hjust = -0.2, vjust = 1.2,
                     gp = gpar(fontface = "bold", cex = 1.4))
    )
  }
  
  # Label plots
  gA <- label_panel(g_up_A, "A");  gB <- label_panel(g_up_B, "B")
  gC <- label_panel(g_up_C, "C");  gD <- label_panel(g_up_D, "D")
  gE <- label_panel(g_dn_E, "E");  gF <- label_panel(g_dn_F, "F")
  gG <- label_panel(g_dn_G, "G");  gH <- label_panel(g_dn_H, "H")
  
  # ─────────────────────────────────────────────────────────────
  # 7. Save all 8 plots to one PDF (2 rows × 4 columns)
  # ─────────────────────────────────────────────────────────────
  pdf("LPS_up_down_Euler_plots.pdf", width = 14, height = 9)
  
  grid.arrange(
    arrangeGrob(gA, gB, gC, gD, ncol = 4),
    arrangeGrob(gE, gF, gG, gH, ncol = 4),
    heights = c(1, 1)
  )
  
  dev.off()
  
  # ─────────────────────────────────────────────────────────────
  # 8. Confirmation message
  # ─────────────────────────────────────────────────────────────
  cat("✅ PDF saved as 'LPS_up_down_Euler_plots.pdf' in the working directory\n")
  
  


  

#stop here
save.image("Image.RData")

#load("Image.RData")



