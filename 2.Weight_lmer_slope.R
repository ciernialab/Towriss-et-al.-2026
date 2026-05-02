#analysis of changes in body weight over time
#August 2025
#AVC


library(tidyverse)
library(cowplot)
library(lubridate) 
library(readxl)
library(lsmeans)
library(nlme)
library(lmerTest)
library(emmeans)

setwd("/Users/aciernia/Sync/CierniaLabMembers/AnnieCiernia/Experiments/SCFA_Invivo2025/weight")



# Set file path
file_path <- "scfa1-4_weights.xlsx"

# Get all sheet names
sheet_names <- excel_sheets(file_path)


# Initialize an empty list to hold each sheet's processed data
list_of_data <- list()

# For loop over each sheet
for (sheet in sheet_names) {
  
  # 1. Read one sheet
  df <- read_excel(file_path, sheet = sheet)
  
  
  # 1. Pivot longer: move all date columns into 'Date' and 'Weight' columns
  df_long <- df %>%
    pivot_longer(
      cols = matches("^\\d"),  # columns that start with a number (your dates)
      names_to = "Date",
      values_to = "Weight"
    ) %>%
    filter(!is.na(Weight))  # remove rows where Weight is NA
  
  # 2. Convert Date from text to actual date format
  
  # ✅ Fix: convert Date to Date class before anything else
  df_long <- df_long %>%
    mutate(
      Date = mdy(str_replace_all(Date, "\\.", "/")),  # Convert "12.9.2024" -> "12/9/2024" first
      Cohort = sheet
    ) 
  
  df_long <- df_long %>%
    group_by(Animal) %>%
    arrange(Date, .by_group = TRUE) %>%   # ✅ Make sure dates are ordered
    mutate(
      First_Date = min(Date, na.rm = TRUE),
      Starting_Weight = Weight[Date == First_Date][1],
      Days_Since_First = as.numeric(Date - First_Date),
      Weight_Percent = (Weight / Starting_Weight) * 100,
      Measurement_Number = row_number()  # ✅ new: count 1,2,3,... per animal
    ) %>%
    ungroup()
  
  
  # 5. Save this processed sheet into the list
  list_of_data[[sheet]] <- df_long
}

# 6. Combine all sheets into one big data frame
all_data <- bind_rows(list_of_data)

# LPS <- read_excel(file_path, sheet = "SCFA4 weight")
# LPS <- LPS %>% select(Animal,Injection, Cage)
# LPS <- distinct(LPS)
# 
# write.csv(LPS,"LPSmetadata.csv")

# Done!


#define factors


all_data$Sex <- as.factor(all_data$Sex)

all_data$Diet <- all_data$Condition
all_data$Diet <- gsub("saline","Saline",all_data$Diet)
all_data$Diet <- gsub("butyrate","Butyrate",all_data$Diet)
all_data$Diet <- gsub("propionate","Propionate",all_data$Diet)
all_data$Diet <- gsub("acetate","Acetate",all_data$Diet)


all_data$Diet <- factor(all_data$Diet, levels = c("Saline"  ,   "Butyrate"  , "Propionate", "Acetate"))




all_data <- all_data %>%
  mutate(
    Day_Group = cut(
      Days_Since_First,
      breaks = seq(0, max(Days_Since_First, na.rm = TRUE) + 4, by = 4),
      include.lowest = TRUE,
      right = FALSE,
      labels = FALSE
    )
  )
# 2. Then, for multiple points per animal per bin, take their mean
library(dplyr)

all_data <- all_data %>%
  filter(!is.na(Day_Group)) %>%  # Only keep non-NA binned points
  dplyr::group_by(Animal, Sex, Diet, Injection,Cohort, Day_Group) %>%
  dplyr::summarize(
    Weight_Percent = mean(Weight_Percent, na.rm = TRUE),
    .groups = "drop"
  )


#fix variables

#remove day 12 for not enough data
all_data <- all_data %>% filter(Day_Group != 9)
all_data$Day_Group <- factor(all_data$Day_Group)



#check variables
table(all_data$Diet)
table(all_data$Sex)
table(all_data$Day_Group)
all_data$Day_Group <- as.numeric(as.character(all_data$Day_Group))



#female
Fdata <- all_data %>% filter(Sex == "female")

#filter out outliers by sex
Fdata_no_outliers <- Fdata %>%
  group_by(Diet,Injection) %>%
  filter(
    Weight_Percent >= quantile(Weight_Percent, 0.25, na.rm = TRUE) - 1.5 * IQR(Weight_Percent, na.rm = TRUE) &
      Weight_Percent <= quantile(Weight_Percent, 0.75, na.rm = TRUE) + 1.5 * IQR(Weight_Percent, na.rm = TRUE)
  ) %>%
  ungroup()

#Male
Mdata <- all_data %>% filter(Sex == "male")

#filter out outliers by sex
Mdata_no_outliers <- Mdata %>%
  group_by(Diet,Injection) %>%
  filter(
    Weight_Percent >= quantile(Weight_Percent, 0.25, na.rm = TRUE) - 1.5 * IQR(Weight_Percent, na.rm = TRUE) &
      Weight_Percent <= quantile(Weight_Percent, 0.75, na.rm = TRUE) + 1.5 * IQR(Weight_Percent, na.rm = TRUE)
  ) %>%
  ungroup()



Alldata <- rbind(Fdata_no_outliers, Mdata_no_outliers)
Alldata$Sex

# Model with only random intercept
model_simple <- lmer(
  Weight_Percent ~ Diet * Sex * Injection * Day_Group + Cohort + (1 + Day_Group | Animal),
  data = Alldata ,
  na.action = na.omit
)

summary(model_simple)

anova <- as.data.frame(anova(model_simple))

write.csv(anova, "lmer_anova.csv", row.names = T)

#To compare how weight changes over time across diets (i.e., is the slope of Day_Group different between diets?), use emtrends():

# 1. Estimate Day_Group slope for each Diet
diet_slopes <- emtrends(model_simple, ~ Diet|Sex*Injection, var = "Day_Group")

# 2. Get Tukey-adjusted pairwise comparisons
tukey_slopes <- pairs(diet_slopes, adjust = "tukey")

# 3. Convert to data frame and save to CSV
tukey_df <- as.data.frame(tukey_slopes)

# 4. Write to CSV
write.csv(tukey_df, "Diet_DayGroup_Slope_Tukey_Comparisons.csv", row.names = FALSE)




# 6. Plot % weight over time!

library(ggplot2)
library(lemon)
library(cowplot)

# Set colors
cbPalette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "purple")

Alldata$Injection <- factor(Alldata$Injection, levels = c("PBS","LPS"))

# Output to PDF
pdf("Boxplot_WeightPercent_Diet_Sex_Day_Group.pdf", height = 10, width = 10)

p <- ggplot(Alldata, aes(x = as.factor(Day_Group), y = Weight_Percent)) + 
  facet_wrap(~Sex*Injection) +
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
   scale_y_continuous(name = "Weight (% of Starting Weight)") +
  theme(strip.background = element_rect(fill = "white", color = "black"))

print(p)

dev.off()


all_data_no_outliers_count <- Alldata %>% select(Sex,Animal,Diet,Injection) %>% distinct()

n <- table( all_data_no_outliers_count$Sex, all_data_no_outliers_count$Diet, all_data_no_outliers_count$Injection)

write.csv(n, "npercondition.csv")


save.image("Data.RData")


