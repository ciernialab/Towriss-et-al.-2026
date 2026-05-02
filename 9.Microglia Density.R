title: "microglial density"
author: "Shreya Gandhi"
date: "`r Sys.Date()`"

## Analysis of microglia density
#Iba1 staining in both male and female mice 

# Load libraries
library(tidyverse)
library(stringr)

# Set working directory
setwd("/Users/aciernia/Sync/CierniaLabMembers/AnnieCiernia/Experiments/SCFA_Invivo2025/microglia_density")

#run for each timmepoint:


#Path to single-cell TIFF images
tiff_dir <- "/Volumes/NINC/CierniaLab/Clianta/SCFA4/SCFA4\ Microglia\ Morphology/Density/SingleCells\ 2"

#List TIFF files
tiff_files <- list.files(tiff_dir, pattern = "\\.tif$", full.names = FALSE)

# Create dataframe from file names
# Assumes filenames like: MouseID_Region_Subregion_Sex_Diet_Stain_123.tif
metadata_df <- data.frame(File = tiff_files) %>%
  mutate(File_no_ext = str_remove(File, "\\.tif$")) %>%
  separate(
    File_no_ext,
    into = c("MouseID", "Sex", "Region","Cohort","Antibody"),
    sep = "_",
    remove = FALSE,
    fill = "right"
  ) %>%
  mutate(
    MouseID = str_to_title(MouseID),
  ) 

#remove .tif from region
unique(metadata_df$Region )
unique(metadata_df$Sex)

#get subfields for HC
metadata_df <- metadata_df %>%
  mutate(HPC_region = str_match(File, "_(DG|CA1|CA2|CA3)\\.tif")[,2])

unique(metadata_df$HPC_region)

#replace DH with subregion
metadata_df <- metadata_df %>%
  mutate(Region = if_else(Region == "DH", HPC_region, Region))

unique(metadata_df$Region )

metadata_df$MouseID <- paste(metadata_df$MouseID,metadata_df$Sex,sep="_")
unique(metadata_df$MouseID)

# Count number of TIFFs (cells) per mouse-subregion group
microglia_counts <- metadata_df %>%
  group_by(MouseID, Sex, Region) %>%
  summarise(num = n(), .groups = "drop")

# Load Areas.csv
areas <- read.csv("/Volumes/NINC/CierniaLab/Clianta/SCFA4/SCFA4\ Microglia\ Morphology/Density/Areas.csv")

# Extract metadata from Label column in Areas.csv
areas <- areas %>%
  mutate(Label_no_ext = str_remove(Label, "\\.tif$")) %>%
  separate(Label_no_ext, into = c("MouseID", "Sex",  "Region","antibody","cohort","HC_region"), sep = "_", remove = TRUE) %>%
  mutate(
    MouseID = paste0(toupper(substr(MouseID, 1, 1)), substr(MouseID, 2, nchar(MouseID))),
  )

#clean up subregions

areas$HC_region <- gsub("\\..*", "", areas$HC_region)   # remove everything from the first '.' onward

#replace DH with subregion
areas <- areas %>%
  mutate(Region = if_else(Region == "DH", HC_region, Region))


unique(areas$Region)

areas$MouseID <- paste(areas$MouseID,areas$Sex,sep="_")
unique(areas$MouseID)


# Merge area info with microglia counts
microglia_counts$num <- as.numeric(microglia_counts$num)
marea <- left_join(areas, microglia_counts, by = c("MouseID", "Region"))

marea$Density <- marea$num/marea$Area

#######################################################################################################################
#read in metadata
#######################################################################################################################
md <- read.csv("SCFA4_metadata_final.csv")


#make mouse id 
md <- md %>%
  mutate(
    Cohort_suffix = case_when(
      grepl("Females 1", Cohort) ~ "F1",
      grepl("Females 2", Cohort) ~ "F2",
      grepl("Males 1", Cohort)   ~ "M1",
      grepl("Males 2", Cohort)   ~ "M2",
      TRUE ~ NA_character_
    ),
    MouseID = paste0("Mouse", ID, "_", Cohort_suffix)
  )

unique(md$MouseID)

#merge in data
DF_ALL <- merge(md,marea,by="MouseID")

write.csv(DF_ALL,"Density_Finaldata.csv")
#######################################################################################################################
#plot
#######################################################################################################################
library(lemon)
library(cowplot)

DF_ALL <- DF_ALL %>% filter(Region != "CC")


# --- data prep ---
dat <- DF_ALL %>%
  mutate(
    Region = factor(Region, levels = c("CA1","CA3","DG","PFC","CP")),
    Diet = factor(Condition, levels = c("Saline","Butyrate"  ,  "Propionate", "Acetate")),
    Sex = as.factor(Sex),
    Injection = factor(Injection, levels= c("PBS","LPS")),
    MouseID = as.factor(MouseID)
  ) 

dat$Sex <- toupper(dat$Sex)

pdf("Boplot_Density.pdf", height = 12, width =12)    # create PNG for the heat map       

cbPalette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")

p <- dat %>%
  ggplot(aes(x = Injection, y = Density, fill = Diet)) +
  facet_rep_wrap(~ Region * Sex,
                 repeat.tick.labels = "bottom",
                 ncol = 4
                # scales = "free_y"
                ) +

  stat_summary(
    geom = "boxplot",
    fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)),
                                    c("ymin", "lower", "middle", "upper", "ymax")),
    position = position_dodge(width = 0.9)
  ) +
  
  geom_point(aes(group = Diet),
             position = position_dodge(width = 0.9),
             shape = 21, size = 2, alpha = 0.7, stroke = 0.3) +
  
  scale_fill_manual(values = cbPalette) +
  scale_x_discrete(name = "Diet") +
  scale_y_continuous(name = "Microglial Density") +
  ggtitle("Microglial Density") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

dev.off();


#######################################################################################################################
#stats: ANOVA+posthocs
#######################################################################################################################

library(nlme)
library(emmeans)
library(car)        # Anova(type = "III")
library(dplyr)
library(readr)
library(stringr)
library(tibble)

options(contrasts = c("contr.sum", "contr.poly"))  # for Type III


regions <- levels(dat$Region)


ANOVAout <- NULL
postout <- NULL

for (rg in regions) {
  dd <- filter(dat, Region == rg)
  
  fit <- lme(
    Density ~ Diet * Injection *Sex,
    random = ~ 1 | MouseID,
    data = dd,
    na.action = na.omit,
    method = "REML"
  )
  
  # Type III ANOVA
  a3 <- car::Anova(fit, type = "III")
  a3_df <- as.data.frame(a3) |>
    rownames_to_column("Effect") |>
    mutate(Region = rg, .before = 1)

  # Marginal means & posthocs
  emm_grid <- emmeans(fit, ~ Diet | Sex| Injection)

  pairs_full <- pairs(emm_grid, adjust = "tukey") |>
    as.data.frame() |>
    mutate(Region = rg, .before = 1)

#save to output
  ANOVAout <- rbind(ANOVAout, a3_df)
  postout <- rbind(postout, pairs_full)
  
  
}

#add significance column
ANOVAout <- ANOVAout  %>%
  mutate(Significance = ifelse(`Pr(>Chisq)` < 0.05, "significant", "ns"))

postout <- postout  %>%
  mutate(Significance = ifelse(`p.value` < 0.05, "significant", "ns"))

write.csv(ANOVAout,"ANOVA_density.csv")
write.csv(postout,"TukeyPosthocs_density.csv")

save.image("image.RData")
