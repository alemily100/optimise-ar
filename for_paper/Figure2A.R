setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")
library(ggh4x)
library(tidyverse)
#import dummy data
ordinal_data <- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/dummy_ordinal.csv")[,-1]

#create plot
ordinal_data$grade[is.na(ordinal_data$grade)]<- "Unreported"

n_pat <- ordinal_data %>%
  filter(timepoint == 0, toxicity == 1) %>%
  group_by(dose) %>%
  summarise(n.pat = n())

prop<- ordinal_data%>%group_by(toxicity, timepoint, dose)%>%count(grade) %>%
  left_join(n_pat, by = "dose") %>%
  mutate(proportion = n / n.pat)%>%dplyr::select(-n)
colnames(prop)<- c("toxicity","timepoint", "dose","grade","n.pat", "proportion")

prop<- prop %>% dplyr::select(-n.pat)

missing_layer <- prop %>%
  group_by(toxicity, timepoint, dose)%>%
  summarise(proportion = 1 - sum(proportion), .groups = "drop") %>%
  filter(proportion > 0) 

missing_layer<- cbind(missing_layer, "Discontinued")
colnames(missing_layer)<-c("toxicity","timepoint","dose","proportion","grade")

missing_layer <- missing_layer %>%
  dplyr::select(names(prop)) 

prop$grade<- as.character(prop$grade)
prop <- bind_rows(prop, missing_layer)

prop$toxicity<- as.factor(prop$toxicity)
prop$grade<- as.factor(prop$grade)
prop$dose<- as.factor(prop$dose)
prop$timepoint<- as.factor(prop$timepoint)

prop$dose <- factor(prop$dose, levels = c("1", "2", "3"))
prop$grade<- factor(prop$grade, levels = c("Discontinued", "Unreported", 5:1))
facet_labels<-c("1"="Dose 1 (N=10)", "2"="Dose 2 (N=10)", "3"="Dose 3 (N=10)")

#Final figure 
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/paper/optimise-ar/for_paper")
pdf("Figure2A.pdf")
prop%>%filter(toxicity==2)%>%
  ggplot(aes(x=as.factor(timepoint)))+
  geom_bar(aes(x=as.factor(timepoint),y = proportion, fill = as.factor(grade)), stat = "identity", position = "stack")+
  scale_fill_manual(values=rev(c("#009E73", "#56B4E9","#F0E442","#E69F00", "#D55E00", "grey","grey30")), labels =rev(c("None", "Mild", "Moderate", "Severe", "Very Severe", "Unreported", "Discontinued")),
                    guide = guide_legend(reverse = FALSE))+
  theme_bw(base_size=14) + 
  facet_wrap2(~as.factor(dose), nrow=3, strip.position = "left", axes="all", labeller = as_labeller(facet_labels))+
  labs(x = "Timepoint",                       
       y = "Proportion",                       
       fill = "Severity")
dev.off()