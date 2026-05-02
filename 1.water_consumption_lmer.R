#water consumption analysis

library(tidyverse)
library(cowplot)
library(lubridate) 
library(readxl)
setwd("/Users/aciernia/Sync/CierniaLabMembers/AnnieCiernia/Experiments/SCFA_Invivo2025/water")

#change in water in ml every 3-4 days

all_data <- read.csv("scfa1-4-master water.csv")

#define factors
all_data$Sex <- as.factor(all_data$Sex)

all_data$Diet <- all_data$Condition
all_data$Diet <- gsub("saline","Saline",all_data$Diet)
all_data$Diet <- gsub("butyrate","Butyrate",all_data$Diet)
all_data$Diet <- gsub("propionate","Propionate",all_data$Diet)
all_data$Diet <- gsub("acetate","Acetate",all_data$Diet)


all_data$Diet <- factor(all_data$Diet, levels = c("Saline"  ,   "Butyrate"  , "Propionate", "Acetate"))

all_data$Cage <- factor(all_data$Cage)

all_data <- all_data %>% gather(Day_Group, water_weight_diff, 6:13)

all_data$Day_Group <- factor(all_data$Day_Group)


table(all_data$Diet)
table(all_data$Sex)
table(all_data$Day_Group)

#stats

library(lsmeans)
library(nlme)
library(lmerTest)


#filter out outliers by sex and diet
data_no_outliers <- all_data %>%
  group_by(Diet,Sex) %>%
  filter(
    water_weight_diff >= quantile(water_weight_diff, 0.25, na.rm = TRUE) - 1.5 * IQR(water_weight_diff, na.rm = TRUE) &
      water_weight_diff <= quantile(water_weight_diff, 0.75, na.rm = TRUE) + 1.5 * IQR(water_weight_diff, na.rm = TRUE)
  ) %>%
  ungroup()




# Convert Day_Group to numeric for emtrends

all_data_no_outliers$Day_Group <- as.character(all_data_no_outliers$Day_Group)

all_data_no_outliers$Day_Group <- gsub( "W","",all_data_no_outliers$Day_Group)
all_data_no_outliers$Day_Group <- as.numeric(all_data_no_outliers$Day_Group )

# Refit the model using numeric Day_Group
model_trend <- lmer(
  water_weight_diff ~ Diet * Sex * Day_Group + Cohort + (1 | Cage),
  data = model_data,
  na.action = na.omit
)

anova <- as.data.frame(anova(model_trend))

write.csv(anova, "water_lmer_anova.csv", row.names = T)



# Extract estimated slopes across Day_Group per Diet, Sex
library(emmeans)
em_trends <- emtrends(model_trend, ~ Diet|Sex, var = "Day_Group")

# Pairwise comparisons of slopes
trend_comparisons <- pairs(em_trends, adjust = "tukey")


# Convert to data frame and save to CSV
tukey_df <- as.data.frame(trend_comparisons)


# 4. Write to CSV
write.csv(tukey_df, "water_Diet_DayGroup_Tukey_Comparisons.csv", row.names = FALSE)



library(ggplot2)
library(lemon)
library(cowplot)

# Set colors
cbPalette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "purple")

# Output to PDF
pdf("Boxplot_water_Diet_Sex_Day_Group.pdf", height = 5, width = 10)

p <- ggplot(all_data_no_outliers, aes(x = as.factor(Day_Group), y = water_weight_diff)) + 
  facet_wrap(~Sex) +
  stat_summary(
    geom = "boxplot", 
    fun.data = function(x) setNames(
      quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), 
      c("ymin", "lower", "middle", "upper", "ymax")
    ), 
    position = position_dodge(width = 0.9),
    aes(fill = Diet)
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
    aes(fill = Diet, group = Diet),  # <-- corrected
    shape = 21,
    color = "black",
    size = 2,
    alpha = 0.5
  ) +
  scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 15) +
  scale_x_discrete(name = "Measurement Day") +
   scale_y_continuous(name = "Water Consumption") +
  theme(strip.background = element_rect(fill = "white", color = "black"))

print(p)

dev.off()



save.image("Data.RData")

all_data_no_outliers_count <- all_data_no_outliers %>% select(Sex,Cage,Diet) %>% distinct()
 
n <- table( all_data_no_outliers_count$Sex, all_data_no_outliers_count$Diet)

write.csv(n, "npercondition.csv")

