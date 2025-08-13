library(RColorBrewer)
library(tidyverse)
library(MASS)
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")
#### FIGURE GENERATION
continuous_data <- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/dummy_continuous.csv")[,-1]

ice_n_pat<-continuous_data %>%
  distinct(id, dose) %>%
  anti_join(
    continuous_data %>% filter(timepoint == 7) %>% distinct(id),
    by = "id"
  ) %>%
  count(dose, name = "n.pat")

continuous_data <- continuous_data %>%
  mutate(score = ifelse(is.na(score), "missing", score))


df <- continuous_data %>%
  filter(timepoint %in% c(0, 7)) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = score,
    names_prefix = "tp_"
  ) %>%
  filter(!is.na(tp_7)) %>%
  mutate(
    tp_0 = as.character(tp_0),
    tp_7 = as.character(tp_7),
    score = ifelse(tp_0 == "missing" | tp_7 == "missing",
                   "missing",
                   as.character(as.numeric(tp_7) - as.numeric(tp_0)))
  ) %>%
  dplyr::select(id, score, dose, functioning)

df<-df %>%group_by(dose, functioning) %>%
  mutate(non_na_count = sum(score!="missing")) %>%
  mutate(na_count = sum(score=="missing")) %>%
  ungroup()

unrep <- df %>%
  group_by(dose, functioning) %>%
  summarise(
    n_missing = sum(score == "missing"),
    .groups = "drop"
  )

df<- df%>%subset(score !="missing")
df$score<- as.numeric(df$score)

sum.df<-df%>%group_by(dose, functioning)%>%summarise(mean=mean(score, na.rm=TRUE),non_na_count = unique(non_na_count), se = sd(score, na.rm=TRUE)/sqrt(non_na_count))
sum.df<-sum.df%>%mutate(l.ci= mean-(1.96*se))
sum.df<-sum.df%>%mutate(u.ci= mean+(1.96*se))

