setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/paper/optimise-ar")
source("functions_synthesise_data.R")
library(tidyverse)
library(lubridate)
library(ggsurvfit)
library(survival)
library(ggsurvfit)
library(gridExtra)
library(grid)
#library(tidycmprsk)
library(survival)

#import dummy data
ordinal_data <- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/paper/optimise-ar/dummy_ordinal.csv")[,-1]


#### FIGURE GENERATION
ordinal_data<-ordinal_data%>%group_by(id, toxicity)%>%mutate(base_adj = grade-grade[timepoint==0])

ordinal_data_available<-ordinal_data%>%filter(is.na(base_adj)==FALSE)

summary_df <- ordinal_data_available %>%
  group_by(id) %>%
  summarise(
    time = if (any(base_adj == 3)) min(timepoint[base_adj == 3]) else max(timepoint),
    event = as.integer(any(base_adj == 3))  # 1 if grade 3 occurred, 0 if censored
  )

summary_df <- cbind(summary_df, ordinal_data%>%group_by(id)%>%dplyr::summarise(max(dose)))
colnames(summary_df)[5]<- "dose"

summary_df$dose<- as.factor(summary_df$dose)
levels(summary_df$dose) <- c("Dose 1", "Dose 2", "Dose 3")

fit<- survfit2(Surv(time,event) ~ dose, data=summary_df)

median_df <- survminer::surv_median(fit)%>%
  dplyr::select(strata, median) %>%
  rename(Dose = strata,
         `Median time \n to event` = median)
median_df[,1]<-c("Dose 1", "Dose 2", "Dose 3")

# tbl_grob <- tableGrob(median_df, rows = NULL, theme = ttheme_default(
#   core = list(fg_params = list(cex = 0.7)),
#   colhead = list(fg_params = list(cex = 0.8, fontface = "bold"))
# ))

tbl_grob <- tableGrob(
  median_df,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.7, col = "black"),     # text
      bg_params = list(fill = NA, col = "black")      # no fill, black borders
    ),
    colhead = list(
      fg_params = list(cex = 0.8, fontface = "bold", col = "black"),
      bg_params = list(fill = NA, col = "black")      # no fill, black borders
    )
  )
)

setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/paper/optimise-ar/for_paper")
pdf("Figure4.pdf", width=10, height=6)
  ggsurvfit(fit, linewidth = 1.5) + 
  add_censor_mark() +
  add_confidence_interval() +
  scale_ggsurvfit() + add_risktable()+
  scale_color_manual(
    values = c(
      "Dose 1" = "#66C2A5",   # Dose level 1
      "Dose 2" = "#8DA0CB",   # Dose level 2
      "Dose 3" = "#FC8D62"    # Dose level 3
    ))+
  scale_fill_manual(
    values = c(
      "Dose 1" = "#66C2A5",   # Dose level 1
      "Dose 2" = "#8DA0CB",   # Dose level 2
      "Dose 3" = "#FC8D62"    # Dose level 3
    )) +
    scale_y_continuous(
      name = expression(
        atop(
          "Proportion without first baseline-adjusted",
          "Grade" >= 3 ~ " PRO-CTCAE score"
        )
      )
    )+
    annotation_custom(
      tbl_grob,
      xmin = -0.5, xmax = 1.75,  # x-position range for table
      ymin = 0.05, ymax = 0.35  # y-position range for table
    )
dev.off()

