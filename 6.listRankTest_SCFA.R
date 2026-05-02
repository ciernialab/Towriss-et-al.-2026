
#analysis to compare ranked lists of genes
#AVC April 2025

#https://htmlpreview.github.io/?https://github.com/RRHO2/RRHO2/blob/master/vignettes/RRHO2.html

library(devtools)
#install_github("RRHO2/RRHO2", build_opts = c("--no-resave-data", "--no-manual"))
library(RRHO2)


# 1.5 Format of the input data
# The input gene list should have 2 columns, 1st column is the gene symbol, 2nd column is the input score.
# The score should be calculated as -log10(pvalue) * sign(effectSize)
# NA is not allowed
# Gene symbols of the two gene lists should be identical, but don’t have to be in the same order

#read in gene lists
path = "/Users/aciernia/Sync/CierniaLabMembers/AnnieCiernia/Experiments/SCFA_Invivo2025/RNAseq_April2026"
setwd(path)
lists <- read.csv("AllDEG_AllConditions_log2CPM.csv")

#males
#Compare butyrate to saline
#make lists
lists$Male_Sal_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.M_Saline_LPSvsPBS, base = 10)

lists$Male_Sal_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.M_Saline_LPSvsPBS == "Decrease" ~ lists$Male_Sal_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Male_Sal_LPSvsPBS_neglogPvalue *1)

lists$Male_Butyrate_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.M_Butyrate_LPSvsPBS, base = 10)

lists$Male_Butyrate_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.M_Butyrate_LPSvsPBS == "Decrease" ~ lists$Male_Butyrate_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Male_Butyrate_LPSvsPBS_neglogPvalue *1)


list1 <- cbind(lists$GeneID, lists$Male_Sal_LPSvsPBS_neglogPvalue2)

list2 <- cbind(lists$GeneID, lists$Male_Butyrate_LPSvsPBS_neglogPvalue2)


#Create the RRHO2 object
RRHO_obj <-  RRHO2_initialize(list1, list2, labels = c("Male Saline", "Male Butyrate"), log10.ind=TRUE, multipleTesting = "none")

#Visualize the heatmap
pdf("RRHO2_heatmap_Male_ButyratevsSal.pdf", width=6, height=6)
RRHO2_heatmap(RRHO_obj)

dev.off()


#Females
#Compare butyrate to saline
#make lists
lists$Female_Sal_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.F_Saline_LPSvsPBS, base = 10)

lists$Female_Sal_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.F_Saline_LPSvsPBS == "Decrease" ~ lists$Female_Sal_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Female_Sal_LPSvsPBS_neglogPvalue *1)

lists$Female_Butyrate_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.F_Butyrate_LPSvsPBS, base = 10)

lists$Female_Butyrate_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.F_Butyrate_LPSvsPBS == "Decrease" ~ lists$Female_Butyrate_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Female_Butyrate_LPSvsPBS_neglogPvalue *1)


list1 <- cbind(lists$GeneID, lists$Female_Sal_LPSvsPBS_neglogPvalue2)

list2 <- cbind(lists$GeneID, lists$Female_Butyrate_LPSvsPBS_neglogPvalue2)


#Create the RRHO2 object
RRHO_obj <-  RRHO2_initialize(list1, list2, labels = c("Female Saline", "Female Butyrate"), log10.ind=TRUE, multipleTesting = "none")

#Visualize the heatmap
pdf("RRHO2_heatmap_Female_ButyratevsSal.pdf", width=6, height=6)
RRHO2_heatmap(RRHO_obj)

dev.off()


#males
#Compare Acetate to saline
#make lists
lists$Male_Sal_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.M_Saline_LPSvsPBS, base = 10)

lists$Male_Sal_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.M_Saline_LPSvsPBS == "Decrease" ~ lists$Male_Sal_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Male_Sal_LPSvsPBS_neglogPvalue *1)

lists$Male_Acetate_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.M_Acetate_LPSvsPBS, base = 10)

lists$Male_Acetate_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.M_Acetate_LPSvsPBS == "Decrease" ~ lists$Male_Acetate_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Male_Acetate_LPSvsPBS_neglogPvalue *1)


list1 <- cbind(lists$GeneID, lists$Male_Sal_LPSvsPBS_neglogPvalue2)

list2 <- cbind(lists$GeneID, lists$Male_Acetate_LPSvsPBS_neglogPvalue2)


#Create the RRHO2 object
RRHO_obj <-  RRHO2_initialize(list1, list2, labels = c("Male Saline", "Male Acetate"), log10.ind=TRUE, multipleTesting = "none")

#Visualize the heatmap
pdf("RRHO2_heatmap_Male_AcetatevsSal.pdf", width=6, height=6)
RRHO2_heatmap(RRHO_obj)

dev.off()


#Females
#Compare Acetate to saline
#make lists
lists$Female_Sal_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.F_Saline_LPSvsPBS, base = 10)

lists$Female_Sal_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.F_Saline_LPSvsPBS == "Decrease" ~ lists$Female_Sal_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Female_Sal_LPSvsPBS_neglogPvalue *1)

lists$Female_Acetate_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.F_Acetate_LPSvsPBS, base = 10)

lists$Female_Acetate_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.F_Acetate_LPSvsPBS == "Decrease" ~ lists$Female_Acetate_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Female_Acetate_LPSvsPBS_neglogPvalue *1)


list1 <- cbind(lists$GeneID, lists$Female_Sal_LPSvsPBS_neglogPvalue2)

list2 <- cbind(lists$GeneID, lists$Female_Acetate_LPSvsPBS_neglogPvalue2)


#Create the RRHO2 object
RRHO_obj <-  RRHO2_initialize(list1, list2, labels = c("Female Saline", "Female Acetate"), log10.ind=TRUE, multipleTesting = "none")

#Visualize the heatmap
pdf("RRHO2_heatmap_Female_AcetatevsSal.pdf", width=6, height=6)
RRHO2_heatmap(RRHO_obj)

dev.off()


#males
#Compare Propionate to saline
#make lists
lists$Male_Sal_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.M_Saline_LPSvsPBS, base = 10)

lists$Male_Sal_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.M_Saline_LPSvsPBS == "Decrease" ~ lists$Male_Sal_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Male_Sal_LPSvsPBS_neglogPvalue *1)

lists$Male_Propionate_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.M_Propionate_LPSvsPBS, base = 10)

lists$Male_Propionate_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.M_Propionate_LPSvsPBS == "Decrease" ~ lists$Male_Propionate_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Male_Propionate_LPSvsPBS_neglogPvalue *1)


list1 <- cbind(lists$GeneID, lists$Male_Sal_LPSvsPBS_neglogPvalue2)

list2 <- cbind(lists$GeneID, lists$Male_Propionate_LPSvsPBS_neglogPvalue2)


#Create the RRHO2 object
RRHO_obj <-  RRHO2_initialize(list1, list2, labels = c("Male Saline", "Male Propionate"), log10.ind=TRUE, multipleTesting = "none")

#Visualize the heatmap
pdf("RRHO2_heatmap_Male_PropionatevsSal.pdf", width=6, height=6)
RRHO2_heatmap(RRHO_obj)

dev.off()


#Females
#Compare Propionate to saline
#make lists
lists$Female_Sal_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.F_Saline_LPSvsPBS, base = 10)

lists$Female_Sal_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.F_Saline_LPSvsPBS == "Decrease" ~ lists$Female_Sal_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Female_Sal_LPSvsPBS_neglogPvalue *1)

lists$Female_Propionate_LPSvsPBS_neglogPvalue <- -log(lists$adj.P.Val.F_Propionate_LPSvsPBS, base = 10)

lists$Female_Propionate_LPSvsPBS_neglogPvalue2 <- case_when(lists$direction.F_Propionate_LPSvsPBS == "Decrease" ~ lists$Female_Propionate_LPSvsPBS_neglogPvalue *-1, TRUE ~ lists$Female_Propionate_LPSvsPBS_neglogPvalue *1)


list1 <- cbind(lists$GeneID, lists$Female_Sal_LPSvsPBS_neglogPvalue2)

list2 <- cbind(lists$GeneID, lists$Female_Propionate_LPSvsPBS_neglogPvalue2)


#Create the RRHO2 object
RRHO_obj <-  RRHO2_initialize(list1, list2, labels = c("Female Saline", "Female Propionate"), log10.ind=TRUE, multipleTesting = "none")

#Visualize the heatmap
pdf("RRHO2_heatmap_Female_PropionatevsSal.pdf", width=6, height=6)
RRHO2_heatmap(RRHO_obj)

dev.off()