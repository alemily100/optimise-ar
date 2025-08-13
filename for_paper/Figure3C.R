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

#unrep<-df%>%group_by(dose, functioning)%>%distinct(na_count)%>%ungroup()

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

label_names <- c("1" = "Global health \n status", "2" = "Physical \n functioning", "3" = "Role \n functioning",
                 "4" = "Emotional \n functioning", "5" = "Cognitive \n functioning", "6" = "Social \n functioning")

thresholds <- data.frame(
  ymin = c(10, -25),
  ymax = c(25, -10),
  Threshold = factor(c("Improvement", "Worsening"))
)

sum.df<-sum.df %>%
  left_join(unrep, by = c("dose", "functioning"))%>%
  mutate(custom_x = paste0("Dose ", dose, "\nUR=", n_missing))

arrowdf <- tibble(functioning = "1")
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/for_paper")
pdf("Figure3C.pdf", width=13, height=6)
data.frame(sum.df)%>% subset(functioning %in% 1:3) %>% ggplot(aes(x = custom_x, y = mean, colour=as.factor(dose))) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = l.ci, ymax = u.ci), width = 0.2, size=1.2) +
  theme_minimal(base_size=14) + ylim(c(-25,25))+
  geom_rect(data = thresholds,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = Threshold),
            alpha = 0.3, inherit.aes = FALSE)+
  scale_fill_manual(
    name = "Threshold",
    values = c("Improvement" = "palegreen1", "Worsening" = "red")
  )+ 
    labs(
    y = "Mean change from baseline to timepoint 7(with 95% CI) \n for subset of EORTC QLQ-C30 functional domains",
    x = paste0("Dose (N=", paste(10-ice_n_pat$n, collapse = "; "), ")"))+
  facet_wrap(~ as.factor(functioning), nrow=1,labeller = as_labeller(label_names), scales = "free_x")+
  scale_x_discrete(expand = expansion(add = c(1,1)))+
  scale_color_manual(
    values = c(
      "1" = "#66C2A5",   # Dose level 1
      "2" = "#8DA0CB",   # Dose level 2
      "3" = "#FC8D62"    # Dose level 3
    ),
    guide = "none")+ 
  theme(fill="Threshold")+
  geom_segment(data = arrowdf, 
               aes(x = 0.5, xend = 0.5, y = 0, yend = 25), 
               colour = "darkgrey", size = 0.8, arrow = arrow(type="closed",
                                                                      length = unit(0.1, "inches")))+
  geom_segment(data = arrowdf, 
               aes(x = 0.5, xend = 0.5, y = 0, yend = -25), 
               colour = "darkgrey", size = 0.8, arrow = arrow(type="closed",
                                                              length = unit(0.1, "inches")))+
  geom_text(data = arrowdf, aes(x = 1.2, y = 20, label = "Improvement"),
            inherit.aes = FALSE)+
  geom_text(data = arrowdf, aes(x = 1.2, y = -20, label = "Worsening"),
            inherit.aes = FALSE)
dev.off()




  
