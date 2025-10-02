setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")
source("functions_synthesise_data.R")
library(lme4)
library(broom.mixed)
### CONTINUOUS DATA GENERATION 
#systematic review where study was identified: https://www.sciencedirect.com/science/article/pii/S0959804921011485#bib17
#CheckMate066: https://pubmed.ncbi.nlm.nih.gov/27405322/
CM <- read.csv("CheckMate066.csv")
#we utilise data for Dacarbazine only
CM<-data.frame(CM) %>% filter(Drug == "D")
n.pat<-rep(30, times=3)
n.dose<-3
pat.sd<-5
baseline_mean<- 66.2
baseline_sd<- 25.1
mean_change<- CM[,3]
se_change<- CM[,4]
n.timepoints<- length(mean_change)
dose_mean_vector<-c(0,-5,-10)
n.timepoints<-8
set.seed(10014)
time_last_report<-pmin(round(rexp(sum(n.pat), log(2)/36),0), n.timepoints)
time<- data.frame(cbind(1:sum(n.pat), time_last_report))
colnames(time)<- c("id", "time_last_report")

set.seed(1250)
continuous_data<-final_cont(n.dose,n.pat, pat.sd, baseline_mean, baseline_sd,  mean_change, se_change, n.timepoints, dose_mean_vector)

continuous_data<-data.frame(continuous_data)%>% left_join(time, by = "id")
## add non-treatment related death to dataset 
continuous_data<-data.frame(continuous_data)%>%group_by(id)%>%mutate(ICE = ifelse(timepoint<=time_last_report, FALSE, TRUE))

missing<-rbinom(nrow(continuous_data), 1, 0.06)
continuous_data$score[missing==1]<- NA
sample<- continuous_data %>%mutate(missing=is.na(score))
## change in score 

continuous_data<-continuous_data%>%group_by(id)%>%subset(ICE==FALSE)

diff<-continuous_data %>%
  group_by(id) %>%  
  mutate(baseline = score[timepoint == 0], 
         adjusted_value = score - baseline) %>%  
  ungroup()


summary<-diff%>%group_by(dose, timepoint)%>%summarise(mean=mean(adjusted_value, na.rm=TRUE), sd = sd(adjusted_value, na.rm = TRUE) / sqrt(n()))
summary<-summary%>%mutate(l.ci=mean-(1.96*sd))%>%mutate(u.ci=mean+(1.96*sd))

pd <- position_dodge(width = 0.2)

fm1 <- lmer(adjusted_value ~  timepoint + factor(dose) + (1 | id) + baseline, diff, na.action = na.omit)
tidy_model <- broom.mixed::tidy(fm1, effects = "fixed", conf.level=0.95, conf.int = TRUE)

summary<-summary%>%mutate(pred = case_when(
  dose == 1 ~ tidy_model$estimate[1]+tidy_model$estimate[2]*timepoint,
  dose == 2 ~ tidy_model$estimate[1]+tidy_model$estimate[2]*timepoint+tidy_model$estimate[3],
  dose == 3 ~ tidy_model$estimate[1]+tidy_model$estimate[2]*timepoint+tidy_model$estimate[4]
))
  
#### FIGURE GENERATION
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/paper/optimise-ar/for_paper")

pdf("case_study.pdf", width=10, height=5)
summary %>% ggplot(aes(x = as.factor(timepoint), y = mean, group=as.factor(dose), shape=as.factor(dose),colour = as.factor(dose))) +
          scale_colour_manual(
             values = c("1" = "#66C2A5",
                        "2" = "#8DA0CB",
                        "3" = "#FC8D62",
                        label = c("Dose level 1", "Dose level 2", "Dose level 3")
           ), name=element_blank()) +
          geom_ribbon(aes(ymin = 10, ymax = 20, fill = "#009E73"), alpha = 0.1, colour = NA) + 
          geom_ribbon(aes(ymin = 5, ymax = 10, fill = "#56B4E9"), alpha = 0.1, colour = NA) + 
          geom_ribbon(aes(ymin = -5, ymax = -10, fill = "#E69F00"), alpha = 0.05, colour = NA) +
          geom_ribbon(aes(ymin = -10, ymax = -20, fill = "#D55E00"), alpha = 0.1, colour = NA) + 
  geom_line(size=1.5, alpha=0.4)+
  theme_minimal(base_size=14) +
  guides(shape = guide_legend(title = "Dose"),
         colour = guide_legend(title = "Dose"))+
          theme_minimal(base_size=14)+
          xlab("Weeks from baseline") + ylab("Mean change from baseline across timepoints \n for EORTC QLQ-C30 Global health status score") +
  scale_fill_manual(name = "Threshold", values=c("#009E73","#56B4E9","#E69F00", "#D55E00"), label=c("Improvement of at least 10","Improvement of at least 5","Worsening of at least 5", "Worsening of at least 10"))+
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0, yend = 20), 
               colour = "darkgrey", size = 0.8, arrow = arrow(type="closed",
                                                              length = unit(0.1, "inches")))+
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0, yend = -20), 
               colour = "darkgrey", size = 0.8, arrow = arrow(type="closed",
                                                              length = unit(0.1, "inches")))+
  geom_text(aes(x = 1.6, y = 17, label = "Improvement"),
            inherit.aes = FALSE)+
  geom_text(aes(x = 1.6, y = -17, label = "Worsening"),
            inherit.aes = FALSE)+
  geom_line(aes(x=as.factor(timepoint), y=pred), linetype=2, size=1.5)
dev.off()
  
no.pat<-sample%>%group_by(timepoint, dose)%>%count(ICE, missing)%>%filter(dose==3)
print(no.pat, n=23)
