#######################################################################
# FINAL, WORKING GSEA PIPELINE
#######################################################################

library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

set.seed(123)

#######################################################################
# 1.load in data and MAP ENSEMBL → SYMBOL (WITHOUT SHRINKING GENE UNIVERSE)
#######################################################################

mout <- read.csv("AllDEG_AllConditions_log2CPM.csv")

ensembl_to_symbol <- bitr(
  mout$GeneID,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Mm.eg.db
)

mout_sym <- mout %>%
  left_join(ensembl_to_symbol, by = c("GeneID" = "ENSEMBL"))


#######################################################################
# 1. load MGEnrichment dataset > symbols
#######################################################################

load("/Users/aciernia/Sync/CierniaLabMembers/AnnieCiernia/MGErichement_newlists_2025/ScienceDirect_files_09Dec2025_00-16-47.787/Microglia.Mouse.GeneListDatabase.Dec2025.Rdata")


head(mouse.master4)

unique(mouse.master4$groups)
mouse.master4 <- mouse.master4 %>% filter(groups != "ASD genetics") %>% filter(groups !=  "ASD regulators")

#make gene sets
gene_sets <- mouse.master4 %>%
  dplyr::select(listname, mgi_symbol) %>%
  dplyr::filter(!is.na(mgi_symbol)) %>%
  dplyr::distinct() %>%                     # remove duplicates
  dplyr::group_by(listname) %>%
  dplyr::summarise(
    genes = list(unique(mgi_symbol)),
    .groups = "drop"
  ) %>%
  tibble::deframe()


#Clean gene sets to remove genes not in our RNAseq
all_genes <- unique(mout_sym$SYMBOL)

gene_sets <- lapply(
  gene_sets,
  function(gs) intersect(gs, all_genes)
)

#remove small genelists
gene_sets <- gene_sets[lengths(gene_sets) >= 10]

#check
length(gene_sets)
gene_sets[[1]]
names(gene_sets)[1:5]



#######################################################################
# 3. PREPARE RANKED GENE LISTS
#######################################################################

geneLists <- mout_sym %>% dplyr::select(contains("logFC"))
geneLists$SYMBOL <- mout_sym$SYMBOL

# Restrict to biologically meaningful contrasts
conditionnames <- grep(
  "(F|M)_(Saline|Acetate|Butyrate|Propionate)_LPSvsPBS",
  colnames(geneLists),
  value = TRUE
)

conditionnames2 <- grep(
  "(F|M)_(LPS|PBS)_(Saline|Acetate|Butyrate|Propionate)",
  colnames(geneLists),
  value = TRUE
)


conditionsnames_master <- append(conditionnames, conditionnames2)

####################################################
#run enrichments for each list x DB
####################################################

library(fgsea)

test1 <- geneLists[,i]
names(test1) <- geneLists$SYMBOL

# Sort
test1 <- sort(test1, decreasing = TRUE)

#remove duplicates
test1 <- test1[!duplicated(names(test1))]

# Sort
test1 <- sort(test1, decreasing = TRUE)

fgsea_res <- fgsea(
  pathways = gene_sets,
  stats = test1,
  minSize = 10,
  maxSize = 500
)

fgsea_res %>%
  dplyr::arrange(padj) %>%
  head()


#loop through each LPS condition vs PBS log2FC and run GSEA enrichment for different MsigDBs
library(dplyr)
library(fgsea)

# initialize as empty data frame
MGEnrichmentGSEA <- data.frame()

for (i in conditionsnames_master) {
  
  # extract ranking vector
  test1 <- geneLists[[i]]
  names(test1) <- geneLists$SYMBOL
  
  # remove NA values
  test1 <- test1[!is.na(test1)]
  
  # remove duplicated gene names
  test1 <- test1[!duplicated(names(test1))]
  
  message("Running GSEA for: ", i)
  
  # skip if too few genes
  if (length(test1) < 20) next
  
  # sort decreasing for GSEA
  test1 <- sort(test1, decreasing = TRUE)
  
  # run fgsea
  MGEnrichmentGSEA1 <- fgsea(
    pathways = gene_sets,
    stats    = test1,
    minSize  = 10,
    maxSize  = 500
  )
  
  # convert to data frame
  MGEnrichmentGSEA1 <- as.data.frame(MGEnrichmentGSEA1)
  
  # skip if no results
  if (nrow(MGEnrichmentGSEA1) == 0) next
  
  # add comparison label (FIXED)
  MGEnrichmentGSEA1$comparision <- i
  
  # append
  MGEnrichmentGSEA <- bind_rows(
    MGEnrichmentGSEA,
    MGEnrichmentGSEA1
  )
}

# Did all comparisons run?
length(unique(MGEnrichmentGSEA$comparision))
length(conditionsnames_master)

# Any missing comparison labels?
table(is.na(MGEnrichmentGSEA$comparision))

# Quick peek
head(MGEnrichmentGSEA)


#add back descriptions
mouse.master4.sum <- mouse.master4 %>% dplyr::select(-ensembl_gene_id,-mgi_symbol,-entrezgene_id) %>% distinct()
dim(mouse.master4.sum)

MG_GSEA <- merge(MGEnrichmentGSEA, mouse.master4.sum, by.x = "pathway",by.y ="listname")


MG_GSEA <- MG_GSEA %>%
  filter(!(is.na(padj) & is.na(NES)))

MG_GSEA_df <- MG_GSEA
library(purrr)

MG_GSEA_df$leadingEdge <- map_chr(MG_GSEA_df$leadingEdge, ~ if (length(.x) == 0) NA_character_ else paste(.x, collapse = ","))


#write out full files
write.csv(MG_GSEA_df,"MGEnrichment_GSEAenrichments.csv")

write.csv(gsea_clean,"MGEnrichment_GSAenrichments_significant.csv")

#clean and filter
library(dplyr)

gsea_clean <- MG_GSEA_df %>%
  filter(
    !is.na(NES),
    !is.na(padj),
    padj < 1        # drop totally uninformative results
  ) %>%
  mutate(
    neglog10_padj = -log10(padj),
    direction = ifelse(NES > 0, "Positive NES", "Negative NES")
  )

#Pick top 5 enrichments per comparison

top_gsea <- gsea_clean %>%
  group_by(comparision) %>%
  arrange(padj) %>%
  slice_head(n = 5) %>%
  ungroup()

#make pathwaynames shorter for graphing
top_gsea <- top_gsea %>%
mutate(
  pathway_short = factor(
    pathway,
    levels = unique(pathway)
  )
)

#Order pathways by average NES
pathway_levels <- with(
  top_gsea,
  tapply(NES, pathway_short, mean, na.rm = TRUE)
) |> sort() |> names()

top_gsea$pathway_short <- factor(top_gsea$pathway_short,
                                 levels = pathway_levels)


#filter
top_gsea <- top_gsea %>%
  filter(padj < 0.05)


#Dot plot
library(ggplot2)
library(cowplot)

p <- ggplot(
  top_gsea,
  aes(
    x = comparision,
    y = pathway_short,
    size = abs(NES),
    color = neglog10_padj
  )
) +
  geom_point(alpha = 0.85) +
  
  scale_size_continuous(
    name = "|NES|",
    range = c(2, 10)
  ) +
  
  scale_color_viridis_c(
    name = expression(-log[10](adj~p)),
    option = "C"
  ) +
  
  theme_cowplot(font_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Comparison",
    y = "Gene set",
    title = "GSEA enrichment across comparison lists"
  )

p + facet_wrap(~ direction, scales = "free_y")



ggsave(
  "GSEA_DotPlot_TopPathways_PerComparison.pdf",
  plot = p,
  width = 14,
  height = 10
)

#######################################################################

# You want separate panels (or separate plots) for:
#   
#   all LPSvsPBS
# all _LPS_ (diet effects within LPS)
# all _PBS_ (diet effects within PBS)


library(dplyr)
library(stringr)

top_gsea <- top_gsea %>%
  mutate(
    comparison_type = case_when(
      str_detect(comparision, "LPSvsPBS")              ~ "LPS vs PBS",
      str_detect(comparision, "_LPS_")                 ~ "Within LPS (Diet)",
      str_detect(comparision, "_PBS_")                 ~ "Within PBS (Diet)",
      TRUE                                             ~ "Other"
    )
  )

table(top_gsea$comparison_type)

top_gsea$comparision <- factor(
  top_gsea$comparision,
  levels = unique(top_gsea$comparision)
)


#plot
library(ggplot2)
library(cowplot)

p <- ggplot(
  top_gsea,
  aes(
    x = comparision,
    y = pathway_short,
    size = abs(NES),
    color = neglog10_padj
  )
) +
  geom_point(alpha = 0.85) +
  
  facet_wrap(
    ~ comparison_type,
    scales = "free_x",
    nrow = 1
  ) +
  
  scale_size_continuous(
    name = "|NES|",
    range = c(2, 9)
  ) +
  
  scale_color_viridis_c(
    name = expression(-log[10]~adj~p),
    option = "C"
  ) +
  
  theme_cowplot(font_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Comparison",
    y = "Gene set",
    title = "GSEA enrichment separated by comparison type"
  )

p

ggsave(
  "GSEA_Dotplot_ByComparisonType.pdf",
  plot = p,
  width = 16,
  height = 8
)

#separate plots

#LPS vs PBS comparisons
p_LPSvsPBS <- top_gsea %>%
  filter(comparison_type == "LPS vs PBS") %>%
  ggplot(aes(
    x = comparision,
    y = pathway_short,
    size = abs(NES),
    color = neglog10_padj
  )) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(2, 9)) +
  scale_color_viridis_c(option = "C") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "GSEA: LPS vs PBS comparisons",
    x = "Comparison",
    y = "Gene set"
  )

ggsave(
  "GSEA_Dotplot_ByComparisonType_onlyLPSvsPBS.pdf",
  plot = p_LPSvsPBS,
  width = 12,
  height = 8
)

#Diet effects within LPS
p_LPS <- top_gsea %>%
  filter(comparison_type == "Within LPS (Diet)") %>%
  ggplot(aes(
    x = comparision,
    y = pathway_short,
    size = abs(NES),
    color = neglog10_padj
  )) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(2, 9)) +
  scale_color_viridis_c(option = "C") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "GSEA: Diet effects within LPS",
    x = "Comparison",
    y = "Gene set"
  )

ggsave(
  "GSEA_Dotplot_ByComparisonType_onlyLPS.pdf",
  plot = p_LPS,
  width = 10,
  height = 8
)

#Diet effects within PBS
p_PBS <- top_gsea %>%
  filter(comparison_type == "Within PBS (Diet)") %>%
  ggplot(aes(
    x = comparision,
    y = pathway_short,
    size = abs(NES),
    color = neglog10_padj
  )) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(2, 9)) +
  scale_color_viridis_c(option = "C") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "GSEA: Diet effects within PBS",
    x = "Comparison",
    y = "Gene set"
  )

ggsave(
  "GSEA_Dotplot_ByComparisonType_onlyPBS.pdf",
  plot = p_PBS,
  width = 10,
  height = 8
)







#######################################################################
#top 5 positive and top 5 negative NES per comparison
#######################################################################

top_gsea <- gsea_clean %>%
  filter(padj < 0.05) %>%                     # significance first
  group_by(comparision) %>%
  
  # top 5 positive NES
  slice_max(NES, n = 5, with_ties = FALSE) %>%
  
  bind_rows(
    gsea_clean %>%
      filter(padj < 0.05) %>%
      group_by(comparision) %>%
      slice_min(NES, n = 5, with_ties = FALSE)
  ) %>%
  
  ungroup() %>%
  
  mutate(
    direction = case_when(
      NES > 0 ~ "Positive NES (Up in numerator)",
      NES < 0 ~ "Negative NES (Down in numerator)"
    )
  )



top_gsea <- top_gsea %>%
  mutate(pathway_short = pathway)

pathway_levels <- with(
  top_gsea,
  tapply(NES, pathway_short, mean, na.rm = TRUE)
) |>
  sort() |>
  names()

top_gsea$pathway_short <- factor(
  top_gsea$pathway_short,
  levels = pathway_levels
)


top_gsea <- top_gsea %>%
  mutate(
    comparison_type = case_when(
      str_detect(comparision, "LPSvsPBS") ~ "LPS vs PBS",
      str_detect(comparision, "_LPS_")    ~ "Within LPS (Diet)",
      str_detect(comparision, "_PBS_")    ~ "Within PBS (Diet)",
      TRUE                                ~ "Other"
    )
  )

library(stringr)
library(dplyr)

#set order
top_gsea <- top_gsea %>%
  mutate(
    sex = case_when(
      str_detect(comparision, "^logFC\\.F_") ~ "F",
      str_detect(comparision, "^logFC\\.M_") ~ "M",
      TRUE ~ NA_character_
    ),
    
    treatment = case_when(
      str_detect(comparision, "Saline")      ~ "Saline",
      str_detect(comparision, "Butyrate")    ~ "Butyrate",
      str_detect(comparision, "Propionate") ~ "Propionate",
      str_detect(comparision, "Acetate")     ~ "Acetate",
      TRUE ~ NA_character_
    )
  )


desired_levels <- top_gsea %>%
  distinct(sex, treatment, comparision) %>%
  mutate(
    sex = factor(sex, levels = c("F", "M")),
    treatment = factor(
      treatment,
      levels = c("Saline", "Butyrate", "Propionate", "Acetate")
    )
  ) %>%
  arrange(sex, treatment) %>%
  pull(comparision)



top_gsea$comparision <- factor(
  top_gsea$comparision,
  levels = desired_levels
)




top_gsea <- top_gsea %>%
  mutate(
    direction = case_when(
      NES > 0 ~ "Positive NES",
      NES < 0 ~ "Negative NES"
    )
  )


#set legend for circles and triangles

top_gsea <- top_gsea %>%
  mutate(
    direction = factor(
      ifelse(NES > 0, "Higher in LPS", "Higher in PBS"),
      levels = c("Higher in LPS", "Higher in PBS")
    )
  )



#separate plots

#LPS vs PBS comparisons
p_LPSvsPBS <- top_gsea %>%
  filter(comparison_type == "LPS vs PBS") %>%
  ggplot(aes(
    x = comparision,
    y = pathway_short,
    size = abs(NES),
    fill = neglog10_padj,
    shape = direction
  )) +
  geom_point(alpha = 0.85, color = "black") +
  
  scale_shape_manual(
    name = "Pathway enrichment",
    values = c(
      "Higher in LPS"   = 21,  # filled circle
      "Higher in PBS" = 25   # filled triangle
    )
  ) +
  
  scale_size_continuous(
    name = "|NES|"
  ) +
  
  scale_fill_viridis_c(
    name = expression(-log[10]~adj~p),
    option = "C"
  ) +
  
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "GSEA: LPS vs PBS comparisons",
    x = "Comparison",
    y = "Gene set"
  )

ggsave(
  "GSEA_LPSvsPBS_final_withInterpretation.pdf",
  p_LPSvsPBS,
  width = 12,
  height = 8
)

#second version of plot
# LPS vs PBS comparisons
p_LPSvsPBS <- top_gsea %>%
  filter(comparison_type == "LPS vs PBS") %>%
  ggplot(aes(
    x = comparision,
    y = pathway_short,
    size = neglog10_padj,
    fill = NES
  )) +
  
  geom_point(
    shape = 21,            # filled circle
    color = "black",       # outline
    alpha = 0.85
  ) +
  
  scale_size_continuous(
    name = expression(-log[10]~adj~p)
  ) +
  
  scale_fill_gradient2(
    name = "NES",
    low = "#2166AC",       # blue
    mid = "white",
    high = "#B2182B",      # red
    midpoint = 0
  ) +
  
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  
  labs(
    title = "GSEA: LPS vs PBS comparisons",
    x = "Comparison",
    y = "Gene set"
  )

ggsave(
  "GSEA_LPSvsPBS_final_withInterpretation_dots.pdf",
  p_LPSvsPBS,
  width = 12,
  height = 8
)


#Diet effects within LPS

order_diet_comparisons <- function(df) {
  df %>%
    mutate(
      sex = case_when(
        str_detect(comparision, "^logFC\\.F_") ~ "F",
        str_detect(comparision, "^logFC\\.M_") ~ "M"
      ),
      diet = case_when(
        str_detect(comparision, "ButyratevsSaline")    ~ "Butyrate",
        str_detect(comparision, "PropionatevsSaline") ~ "Propionate",
        str_detect(comparision, "AcetatevsSaline")    ~ "Acetate"
      )
    ) %>%
    mutate(
      sex = factor(sex, levels = c("F", "M")),
      diet = factor(
        diet,
        levels = c("Butyrate", "Propionate", "Acetate")
      )
    ) %>%
    arrange(sex, diet) %>%
    pull(comparision) %>%
    unique()
}



top_gsea <- top_gsea %>%
  mutate(
    direction = factor(
      ifelse(NES > 0, "Higher in LPS+SCFA", "Higher in LPS+saline"),
      levels = c("Higher in LPS+SCFA", "Higher in LPS+saline")
    )
  )



p_LPS <- top_gsea %>%
  filter(comparison_type == "Within LPS (Diet)") %>%
  mutate(comparision = factor(comparision, levels = lps_levels)) %>%
  ggplot(aes(
    x = comparision,
    y = pathway_short,
    size = abs(NES),
    fill = neglog10_padj,
    shape = direction
  )) +
  geom_point(alpha = 0.85, color = "black") +
  
  scale_shape_manual(
    name = "Pathway enrichment",
    values = c(
      "Higher in LPS+SCFA"   = 21,  # filled circle
      "Higher in LPS+saline" = 25   # filled triangle
    )
  )+
  
  scale_size_continuous(
    name = "|NES|",
    range = c(2, 9)
  ) +
  
  scale_fill_viridis_c(
    name = expression(-log[10]~adj~p),
    option = "C"
  ) +
  
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "GSEA: Diet effects within LPS",
    x = "Comparison",
    y = "Gene set",
    shape = "NES direction"
  )

ggsave(
  "GSEA_Dotplot_ByComparisonType_onlyLPS_final.pdf",
  p_LPS,
  width = 10,
  height = 8
)

#second plot option
p_LPS <- top_gsea %>%
  filter(comparison_type == "Within LPS (Diet)") %>%
  mutate(comparision = factor(comparision, levels = lps_levels)) %>%
  ggplot(aes(
    x = comparision,
    y = pathway_short,
    size = neglog10_padj,
    fill = NES
  )) +
  
  geom_point(
    shape = 21,        # filled circle only
    color = "black",   # outline
    alpha = 0.85
  ) +
  
  scale_size_continuous(
    name = expression(-log[10]~adj~p),
    range = c(2, 9)
  ) +
  
  scale_fill_gradient2(
    name = "NES",
    low = "#2166AC",   # blue (negative NES)
    mid = "white",
    high = "#B2182B",  # red (positive NES)
    midpoint = 0
  ) +
  
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  
  labs(
    title = "GSEA: Diet effects within LPS",
    x = "Comparison",
    y = "Gene set"
  )

ggsave(
  "GSEA_Dotplot_ByComparisonType_onlyLPS_dots.pdf",
  p_LPS,
  width = 10,
  height = 12
)






#Diet effects within PBS


pbs_levels <- top_gsea %>%
  filter(comparison_type == "Within PBS (Diet)") %>%
  order_diet_comparisons()

top_gsea <- top_gsea %>%
  mutate(
    direction = factor(
      ifelse(NES > 0, "Higher in PBS+SCFA", "Higher in PBS+saline"),
      levels = c("Higher in PBS+SCFA", "Higher in PBS+saline")
    )
  )

p_PBS <- top_gsea %>%
  filter(comparison_type == "Within PBS (Diet)") %>%
  mutate(comparision = factor(comparision, levels = pbs_levels)) %>%
  ggplot(aes(
    x = comparision,
    y = pathway_short,
    size = abs(NES),
    fill = neglog10_padj,
    shape = direction
  )) +
  geom_point(alpha = 0.85, color = "black") +
  
  scale_shape_manual(
    name = "Pathway enrichment",
    values = c(
      "Higher in PBS+SCFA"   = 21,  # filled circle
      "Higher in PBS+saline" = 25   # filled triangle
    )
  ) +
  
  scale_size_continuous(
    name = "|NES|",
    range = c(2, 9)
  ) +
  
  scale_fill_viridis_c(
    name = expression(-log[10]~adj~p),
    option = "C"
  ) +
  
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "GSEA: Diet effects within PBS",
    x = "Comparison",
    y = "Gene set",
    shape = "NES direction"
  )

ggsave(
  "GSEA_Dotplot_ByComparisonType_onlyPBS_final.pdf",
  p_PBS,
  width = 10,
  height = 8
)


#second plot option
p_PBS <- top_gsea %>%
  filter(comparison_type == "Within PBS (Diet)") %>%
  mutate(comparision = factor(comparision, levels = pbs_levels)) %>%
  ggplot(aes(
    x = comparision,
    y = pathway_short,
    size = neglog10_padj,
    fill = NES
  )) +
  
  geom_point(
    shape = 21,        # filled circle only
    color = "black",   # outline
    alpha = 0.85
  ) +
  
  scale_size_continuous(
    name = expression(-log[10]~adj~p),
    range = c(2, 9)
  ) +
  
  scale_fill_gradient2(
    name = "NES",
    low = "#2166AC",   # blue (negative NES)
    mid = "white",
    high = "#B2182B",  # red (positive NES)
    midpoint = 0
  ) +
  
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  
  labs(
    title = "GSEA: Diet effects within PBS",
    x = "Comparison",
    y = "Gene set"
  )

ggsave(
  "GSEA_Dotplot_ByComparisonType_onlyPBS_dots.pdf",
  p_PBS,
  width = 10,
  height = 10
)







