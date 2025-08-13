setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")

continuous_data <- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/dummy_continuous.csv")[,-1]

#### FIGURE GENERATION
continuous_data<- continuous_data%>%filter(functioning==1)%>%dplyr::select(-functioning)

diff<-continuous_data %>%
  group_by(id) %>%  
  mutate(baseline = score[timepoint == 0], 
         adjusted_value = score - baseline) %>%  
  ungroup()


summary<-diff%>%group_by(dose, timepoint)%>%summarise(mean=mean(adjusted_value, na.rm=TRUE), sd = sd(adjusted_value, na.rm = TRUE) / sqrt(n()))
summary<-summary%>%mutate(l.ci=mean-(1.96*sd))%>%mutate(u.ci=mean+(1.96*sd))
pd <- position_dodge(width = 0.2)

#### FIGURE GENERATION
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/for_paper")
pdf("Figure2C.pdf", width=10, height=6)
summary %>% ggplot(aes(x = as.factor(timepoint), y = mean, group=as.factor(dose), shape=as.factor(dose),colour = as.factor(dose))) +
          scale_colour_manual(
             values = c("1" = "#66C2A5",
                        "2" = "#8DA0CB",
                        "3" = "#FC8D62",
                        label = c("Dose level 1", "Dose level 2", "Dose level 3")
           ), name=element_blank()) +
          geom_ribbon(aes(ymin = 10, ymax = 30, fill = "palegreen1"), alpha = 0.3, colour = NA) + 
          geom_ribbon(aes(ymin = -30, ymax = -10, fill = "red"), alpha = 0.1, colour = NA) +
          geom_point(position = pd, size=3) + geom_line(size=1.5)+
  geom_errorbar(position = pd, aes(ymin = l.ci, ymax = u.ci), width = 0.2, size=1.2)+
  theme_minimal(base_size=14) +
  guides(shape = guide_legend(title = "Dose"),
         colour = guide_legend(title = "Dose"))+
          theme_minimal(base_size=14)+
          xlab("Timepoint") + ylab("Mean change from baseline across timepoints (with 95% CI) \n for EORTC QLQ-C30 Global health status score") +
  scale_fill_manual(name = "Threshold", values=c("palegreen1", "red"), label=c("Improvement", "Worsening"))+
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0, yend = 30), 
               colour = "darkgrey", size = 0.8, arrow = arrow(type="closed",
                                                              length = unit(0.1, "inches")))+
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0, yend = -30), 
               colour = "darkgrey", size = 0.8, arrow = arrow(type="closed",
                                                              length = unit(0.1, "inches")))+
  geom_text(aes(x = 1.6, y = 17, label = "Improvement"),
            inherit.aes = FALSE)+
  geom_text(aes(x = 1.6, y = -17, label = "Worsening"),
            inherit.aes = FALSE)
dev.off()

#find intercurrent events 
all_ids_doses <- continuous_data %>%
  distinct(id, dose)

timepoints <- unique(continuous_data$timepoint)

full_grid <- all_ids_doses %>%
  crossing(timepoint = timepoints)

missing <- full_grid %>%
  anti_join(
    continuous_data %>% distinct(id, timepoint),
    by = c("id", "timepoint")
  )

ice_n_pat <- missing %>%
  count(dose, timepoint, name = "n_missing_patients")

ice_n_pat


#find unreporting 
continuous_data <- continuous_data %>%
  mutate(score = ifelse(is.na(score), "missing", score))

unrep<- continuous_data %>%
  group_by(timepoint, dose)%>%
  summarise(n_missing_score = sum(score=="missing"))

