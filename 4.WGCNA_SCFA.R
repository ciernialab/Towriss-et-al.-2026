#####################################################################################
# Gene expression from isolated microglia after 4 wks SCFA treatment then PBS/LPS 24hr sac
####################################################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)

#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gplots)


#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

######################################################################

############################################################################################################################
#Initial data processing, sample cleaning and differential expression analysis 
############################################################################################################################

############################################################################################################################
##LOADING IN THE DATA##
path = "/Users/aciernia/Sync/CierniaLabMembers/AnnieCiernia/Experiments/SCFA_Invivo2025/RNAseq_April2026/WGCNA"
setwd(path)


#raw counts
counts <- read.csv("FinalRawCountsMatrix.csv")
rownames(counts) <- counts$X
counts <- counts %>% dplyr::select(-X)
counts <- as.matrix(counts)


#read in metadata
metadata <- read.csv("Final_Metadata.csv")

#remove SCFA2 11 > appears to be wrong LPS dose

metadata <- metadata %>% filter(Sample.ID != "SCFA2.11")

#match samples
metadata$Sample.ID
colnames(counts)

#match
counts2 <- counts[, which(colnames(counts) %in% metadata$Sample.ID)]
counts3 <- counts2[,metadata$Sample.ID]
colnames(counts3)

counts <- as.matrix(counts3)

#new treatment variable that is LPS+sex

metadata$Diet <- factor(metadata$Diet, levels = c("Saline"  ,   "Butyrate"  , "Propionate", "Acetate"))

metadata$LPS.Treatment <- factor(metadata$Injection, levels = c("PBS","LPS"))


metadata$condition <- paste(metadata$Sex, metadata$Diet, metadata$LPS.Treatment, sep="_")

metadata$condition = factor(metadata$condition)



#Make DEG list object
d0 <- DGEList(counts, group = metadata$condition)
#d0$genes = genenames
dim(d0)
############################################################################################
##getting rid of low expressed genes (cutoff for the number of samples that we want above 1,getting rid of things that are less than one group)
############################################################################################

#filter by group
#filter by group
keep <- filterByExpr(
  d0,
  group = metadata$condition,
  min.prop = 0.5
)

d <- d0[keep, , keep.lib.sizes = FALSE]

dim(d)

#21871    83



# Log-transformed CPM
cpm_data <- cpm(d, log = TRUE, prior.count = 1)

#remove low variance genes (keeps top 75% most variable)
var_genes <- apply(cpm_data, 1, var)
expr_filtered <- cpm_data[var_genes > quantile(var_genes, 0.25), ]


#remove batch effects

library(limma)

# expr: genes × samples (log2 normalized)
# batch: factor
# design: optional biological variables to preserve

expr_corrected <- removeBatchEffect(
  expr_filtered,
  batch = metadata$Cohort_1,
  design = model.matrix(~ Diet * Injection * Sex, data = metadata)
)

dim(expr_corrected)
#16403    83


############################################################################################################################
############################################################################################################################
############################################################################################################################
#WGCNA


# Allow WGCNA functions to use multiple threads
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
library(WGCNA)
############################################################################################################################
############################################################################################################################
############################################################################################################################


#make gene matrix with genes as columns and rows as samples

datExpr <- as.data.frame(t(expr_corrected)) # Transpose for WGCNA input (genes as columns)


################## process data ########################
#detect genes with missing values or samples with missing values
# Assuming `datExpr` is your log-CPM data with genes as columns and samples as rows
# Check for missing data
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- expr_corrected[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK # TRUE


write.csv(datExpr,file="Filtered_inputWGCNAdata.csv")
# 

filteredData <- datExpr
#=====================================================================================
#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(filteredData), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)


pdf(file = "sampleClustering_preCut.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
#abline(h = 160, col = "red");

dev.off()

# # Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)
# table(clust)
# # clust 1 contains the samples we want to keep. > all in this case
# keepSamples = (clust==1)
# datExpr = DF4[keepSamples, ]
nGenes = ncol(filteredData)
nSamples = nrow(filteredData)

#=====================================================================================
#phenotype data:
info <- metadata

datTraits <- info


#match sample names and trait data:
Samples = rownames(filteredData)
traitRows = match(Samples, info$Sample.ID)
datTraits = info[traitRows, ]

datTraits<- as.data.frame(datTraits)
datTraits$Cohort <- as.factor(datTraits$Cohort)
datTraits$Sample.ID

#change genotype to numeric #M =2 = red,  F= 1 = white
datTraits$sex2 <- as.numeric(as.factor(datTraits$Sex)) 
datTraits$LPS.Treatment2 <- as.numeric(datTraits$LPS.Treatment)
datTraits$Diet2 <- as.numeric(datTraits$Diet)
datTraits$Cohort2 <- as.numeric(datTraits$Cohort)



#numeric traits only
numerictraits <- datTraits %>% dplyr::select(sex2, LPS.Treatment2,Diet2,Cohort2)
rownames(numerictraits) <- datTraits$Sample.ID

#We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable datTraits. Before we continue with network construction and module detection, we visualize how the clinical traits relate to the sample dendrogram.
# Re-cluster samples
sampleTree2 = hclust(dist(filteredData), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors1 = numbers2colors(numerictraits$Diet2,signed = FALSE);
  traitColors2 = numbers2colors(numerictraits$LPS.Treatment2,signed = FALSE);
  traitColors3 = numbers2colors(numerictraits$sex2,signed = FALSE);
  traitColors4 = numbers2colors(numerictraits$Cohort2,signed = FALSE);
  traitColors <-cbind(traitColors1,traitColors2,traitColors3,traitColors4)
  colnames(traitColors) <- c("Diet","LPS", "sex","cohort")
# Plot the sample dendrogram and the colors underneath.

pdf(file = "sampleClustering_withtraits.pdf", width = 12, height = 9)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = colnames(traitColors),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(numerictraits,datTraits,filteredData, file = "WTWGCNAnetworkConstruction-inputdata.RData")


#=====================================================================================
#choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency 
# Choose a set of soft-thresholding powers
powers = c(c(1:30))
# Call the network topology analysis function
sft = pickSoftThreshold(filteredData, powerVector = powers, networkType = "signed",corFnc="bicor" ,verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="blue");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="black")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")

write.csv(sft, "SoftThresholding_powercalc.csv")
#lowest power for which the scale-free topology fit index of > 0.8 ->10





#=====================================================================================
#Run this section on an external server, save the network data and load into R studio
#=====================================================================================
#Co-expression similarity and adjacency
library(WGCNA)
allowWGCNAThreads()

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#load data


Allnet = blockwiseModules(filteredData, maxBlockSize = 50000,
                          power = 8, TOMType = "signed", minModuleSize = 25,
                          networkType = "signed",
                          corType = "bicor", #biweight midcorrelation
                          maxPOutliers = 0.05, #forces bicor to never regard more than the specified proportion of samples as outliers.
                          reassignThreshold = 0,
                          numericLabels = TRUE,
                          saveTOMs = FALSE,
                          nThreads = 12,
                          #saveTOMFileBase = "TOM-BilboMGExpression",
                          verbose = 3)

save(Allnet, file = "NetworkConstruction-auto_WTWGCNA.RData")
#=====================================================================================
table(Allnet$colors)

# Convert labels to colors for plotting
AllnetModuleColors = labels2colors(Allnet$colors)


AllnetTree = Allnet$dendrograms[[1]]

#plot the gene dendrogram and the corresponding module colors
sizeGrWindow(8,6);

pdf(file = "AllnetDendrogram_pretrim.pdf", wi = 8, he = 6)
plotDendroAndColors(AllnetTree, AllnetModuleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "All Samples Cluster Dendrogram")
dev.off()


#=====================================================================================

# Calculate eigengenes
MEList = moduleEigengenes(filteredData, colors = AllnetModuleColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

pdf(file = "Pretrim_allnet_dendrogram.pdf", wi = 9, he = 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=0.1, col = "red")

dev.off()



#=====================================================================================
#merge down to fewer modules based on plot height
MEDissThres = 0.1

# Call an automatic merging function
merge = mergeCloseModules(filteredData, AllnetModuleColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#replot
sizeGrWindow(12, 9)
pdf(file = "Allnet_postmerg.4MEClustering.pdf", wi = 9, he = 6)

plotDendroAndColors(AllnetTree, mColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "All Samples Cluster Dendrogram")
dev.off()

#=====================================================================================
# Rename to moduleColors
moduleColors = mColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

#=====================================================================================
#repeat clustering of new MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
pdf(file = "allnet_dendrogram_postmerge.pdf", wi = 9, he = 6)
plot(METree, main = "Clustering of module eigengenes after merging modules",
     xlab = "", sub = "")

dev.off()

#=====================================================================================
#correlate with diet, lps and sex
#=====================================================================================
#correlate eigengenes with external traits and look for the most signicant associations
# Define numbers of genes and samples
nGenes = ncol(filteredData);
nSamples = nrow(filteredData);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(filteredData, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Check alignment of samples between MEs and traits
rownames(datTraits) <- datTraits$Sample.ID

datTraits$Interaction <- paste(datTraits$Diet, datTraits$Injection, sep="_")
datTraits$Interaction2 <- as.numeric(as.factor(datTraits$Interaction))

datTraitscort <- datTraits %>% dplyr::select(Diet2,LPS.Treatment2, Interaction2, sex2)
all(rownames(MEs) == rownames(datTraitscort)) # Should return TRUE


#correlate
moduleTraitCor = cor(MEs,datTraitscort, use = "p");
moduleTraitCor <- moduleTraitCor[-nrow(moduleTraitCor), ]
moduleTraitCor <- moduleTraitCor[-nrow(moduleTraitCor), ]

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#FDR correct the pvalues
#sapply(pval,p.adjust,method="fdr") #per column

#FDR correction on entire matrix
FDR <- matrix(p.adjust(as.vector(moduleTraitPvalue), method='fdr'),ncol=ncol(moduleTraitPvalue))
colnames(FDR) <- colnames(moduleTraitPvalue)
rownames(FDR) <- rownames(moduleTraitPvalue)


corout <- cbind(moduleTraitCor,FDR)
write.csv(corout, "PearsonsCorrelations_ME_traits.csv",row.names = T)

#We color code each association by the correlation value:

names(datTraitscort) <- c("Diet","LPS","Diet x LPS","Sex")

sizeGrWindow(10,10)


pdf(file = "module-traits_relationshipsFDRcorrected.pdf", wi = 16, he = 14)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\ (", #pearson correlation coefficient, space, (FDR corrected pvalue)
                    signif(FDR, 4), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(12, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraitscort),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


#=====================================================================================
#graphs of ME
#=====================================================================================


datTraits <- as.data.frame(datTraits)
datTraits$subject <- datTraits$Sample.ID


MEs$subject <- rownames(MEs)
datComb <- merge(datTraits,MEs, by="subject")
datComb2 <- datComb %>% gather(module, ME, 17:ncol(datComb))


mod <- unique(datComb2$module)
#remove grey
mod <- head(mod, -1)

anova_out <- NULL
posthoc1_out <- NULL
posthoc2_out <- NULL

for (i in mod) {
  print(i)
  #select data
  dattmp <- datComb2 %>% dplyr::filter(module == i)
  
  
  # Set color palette
  cbPalette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")
  
  # Plot
  p <- ggplot(dattmp, aes(x = LPS.Treatment, y = ME, fill = Diet)) + 
    facet_wrap(~Sex) +
    stat_summary(
      geom = "boxplot", 
      fun.data = function(x) setNames(
        quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), 
        c("ymin", "lower", "middle", "upper", "ymax")
      ), 
      position = position_dodge(width = 0.9),
      aes(group = interaction(Diet, LPS.Treatment))
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
      aes(group = interaction(Diet, LPS.Treatment), color = Sex),
      size = 2, alpha = 0.7
    ) +
    scale_fill_manual(values = cbPalette) +
    scale_color_manual(values = c("female" = "purple", "male" = "blue")) +
    theme_cowplot(font_size = 15) +
   # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(name = "LPS Treatment") +
    scale_y_continuous(name = "Eigengene Expression") +
    ggtitle(paste(i)) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black")
    )
  
  # Display plot
  p
  
  
  pdf(file=paste("boxplot_",i,".pdf",sep=""), height = 5, width =6)   
  print(p) 
  dev.off()
  
}




#=====================================================================================

#We have found modules with high association with our trait of interest, and have identified their central players by the Module Membership measure.
#We now merge this statistical information with gene annotation and write out a file that summarizes the most important results and can be inspected in standard spreadsheet software 
# Create the starting data frame
Ensemble = colnames(filteredData)
Ensemble <- Ensemble[!is.na(Ensemble) & !is.infinite(Ensemble)]


library(org.Mm.eg.db)

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = Ensemble,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Map Ensembl IDs to full names
gene_descriptions <- mapIds(
  org.Mm.eg.db,
  keys = Ensemble,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Combine into a data frame
gene_info <- data.frame(
  ensembl_gene_id = names(gene_symbols),
  external_gene_name = gene_symbols,
  description = gene_descriptions
)



#combine

tmp = data.frame(ensembl_gene_id = colnames(filteredData),
                 moduleColor = moduleColors,
                 t(filteredData))

#merge
geneInfo0 <- merge(gene_info,tmp, by= "ensembl_gene_id")

write.csv(geneInfo0, file = "geneInfo_ME.csv")





### -------------------------
### MGENrichment plots
### -------------------------


library(tidyverse)
library(cowplot)


# -------------------------------
# 1. Read enrichment results
# -------------------------------
blue_enrich <- read.csv(
  "MG_Enrichment_Results_bluemodule.csv",
  stringsAsFactors = FALSE
) %>%
  mutate(module = "Blue module")

black_enrich <- read.csv(
  "MG_Enrichment_Results_BlackModuleGenes.csv",
  stringsAsFactors = FALSE
) %>%
  mutate(module = "Black module")

# -------------------------------
# 2. Combine and select top pathways
#    (lowest FDR within each module)
# -------------------------------
top_n_pathways <- 10

#select top pathways by OR and fix FDR = 0
plot_df <- bind_rows(blue_enrich, black_enrich) %>%
  filter(!is.na(FDR), !is.na(OR)) %>%
  group_by(module) %>%
  arrange(desc(OR)) %>%          # ← sort by OR (largest first)
  slice_head(n = top_n_pathways) %>%
  ungroup() %>%
  mutate(
    FDR_plot = ifelse(FDR == 0, .Machine$double.xmin, FDR),
    neglog10_FDR = -log10(FDR_plot),
    listname = factor(listname, levels = rev(unique(listname)))
  )




# -------------------------------
# 3. Dot plot (size = OR)
# -------------------------------
p_modules <- ggplot(
  plot_df,
  aes(
    x = module,
    y = listname,
    size = OR,
    fill = neglog10_FDR
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    alpha = 0.85
  ) +
  
  scale_size_continuous(
    name = "Odds ratio (OR)"
  ) +
  
  scale_fill_viridis_c(
    name = expression(-log[10]~FDR),
    option = "C"
  ) +
  
  theme_cowplot() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey90")
  ) +
  
  labs(
    title = "Top enriched pathways in Blue and Black modules",
    x = "Module",
    y = "Pathway"
  )

# -------------------------------
# 4. Save
# -------------------------------
ggsave(
  "Module_Enrichment_Blue_vs_Black_dotplot_ORsize.pdf",
  p_modules,
  width = 8,
  height = 5
)

############################################################
# GO TERM ENRICHMENT FOR WGCNA MODULES
# - Biological Process (BP)
# - Molecular Function (MF)
# - Cellular Component (CC)
#
# INPUT:
#   geneInfo_ME.csv
#     - external_gene_name : mouse gene symbols
#     - moduleColor        : module assignment
#

############################################################


#############################
# 1. LOAD REQUIRED PACKAGES
#############################

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})


#############################
# 2. READ INPUT DATA
#############################

# Read the gene-level table
gene_info <- read.csv(
  "geneInfo_ME.csv",
  stringsAsFactors = FALSE
)


# Basic sanity filtering:
#  - remove missing gene symbols
#  - remove empty symbols
gene_info <- gene_info %>%
  filter(
    !is.na(external_gene_name),
    external_gene_name != ""
  )


#############################
# 3. DEFINE GENE UNIVERSE
#############################

# Background universe = all genes used in the analysis
gene_universe <- unique(gene_info$external_gene_name)


#############################
# 4. FUNCTION: RUN GO ENRICHMENT FOR ONE MODULE
#############################

run_go_for_module <- function(gene_symbols, universe, module_name) {
  
  # GO domains to test
  go_ontologies <- c("BP", "MF", "CC")
  
  # Run enrichGO separately for each ontology
  go_results <- lapply(go_ontologies, function(ont) {
    
    ego <- enrichGO(
      gene          = gene_symbols,
      universe      = universe,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,   # keep all terms; filter later
      qvalueCutoff  = 1,
      readable      = TRUE
    )
    
    # If no enrichment results, return NULL
    if (is.null(ego) || nrow(ego@result) == 0) {
      return(NULL)
    }
    
    # Convert result slot to tidy data frame
    ego@result %>%
      mutate(
        module   = module_name,
        ontology = ont
      )
  })
  
  # Combine BP, MF, CC into one data frame
  bind_rows(go_results)
}


#############################
# 5. RUN GO ENRICHMENT FOR ALL MODULES
#############################

# Split genes by moduleColor
# Run GO enrichment per module
# Combine all modules into one table
all_go_results <- gene_info %>%
  split(.$moduleColor) %>%
  purrr::imap_dfr(function(df, module_name) {
    
    module_genes <- unique(df$external_gene_name)
    
    message(
      "Running GO enrichment for module: ",
      module_name,
      " (n = ", length(module_genes), " genes)"
    )
    
    run_go_for_module(
      gene_symbols = module_genes,
      universe     = gene_universe,
      module_name  = module_name
    )
  })


#############################
# 6. CLEAN AND FORMAT OUTPUT
#############################

# Select and rename columns for clarity
all_go_results <- all_go_results %>%
  select(
    module,
    ontology,
    ID,
    Description,
    GeneRatio,
    BgRatio,
    Count,
    pvalue,
    p.adjust,
    qvalue,
    geneID
  ) %>%
  rename(
    FDR = p.adjust
  ) %>%
  arrange(module, ontology, FDR)


#############################
# 7. WRITE OUTPUT FILE
#############################

write.csv(
  all_go_results,
  file = "GO_Enrichment_AllModules_BP_MF_CC.csv",
  row.names = FALSE
)

#############################
#plot
#############################

library(tidyverse)
library(cowplot)
library(scales)


# Helper function to convert "a/b" to numeric values
parse_ratio <- function(x) {
  sapply(strsplit(x, "/"), function(y) as.numeric(y[1]) / as.numeric(y[2]))
}

go_df <- all_go_results %>%
  mutate(
    gene_ratio = parse_ratio(GeneRatio),
    bg_ratio   = parse_ratio(BgRatio),
    
    # Approximate Odds Ratio
    OR = (gene_ratio / (1 - gene_ratio)) /
      (bg_ratio   / (1 - bg_ratio)),
    
    # Handle zero FDR for plotting
    FDR_plot = ifelse(p.adjust == 0, .Machine$double.xmin, p.adjust),
    neglog10_FDR = -log10(FDR_plot)
  ) %>% filter(p.adjust > 0.05)

#filter
top_n_terms <- 10


plot_go_df <- go_df %>%
  filter(
    module %in% c("blue", "black"),
    ontology %in% c("BP", "MF", "CC"),
    is.finite(OR),
    !is.na(OR)
  ) %>%
  group_by(module) %>%
  arrange(desc(OR), p.adjust) %>%
  slice_head(n = top_n_terms) %>%
  ungroup() %>%
  mutate(
    Description = factor(Description, levels = rev(unique(Description))),
    module = factor(module, levels = c("blue", "black"))
  )


p_go_modules <- ggplot(
  plot_go_df,
  aes(
    x = module,
    y = Description,
    size = OR,
    fill = neglog10_FDR
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    alpha = 0.85
  ) +
  
  scale_size_continuous(
    name = "Odds ratio (OR)",
    trans = "log10"
  ) +
  
  scale_fill_viridis_c(
    name = expression(-log[10]~FDR),
    option = "C"
  ) +
  
  theme_cowplot() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 11),
    panel.grid.major.x = element_line(color = "grey90")
  ) +
  
  labs(
    title = "Top GO terms enriched in Blue and Black modules",
    x = "Module",
    y = "GO term"
  )


ggsave(
  "GO_Enrichment_Blue_vs_Black_dotplot_ORsize.pdf",
  p_go_modules,
  width = 9,
  height = 7
)



 save.image("Workspace.RData")
 
 #load("Workspace.RData")
