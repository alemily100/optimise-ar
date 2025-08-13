setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")

continuous_data <- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/dummy_continuous.csv")[,-1]

#### FIGURE GENERATION
continuous_data<- continuous_data%>%filter(functioning==1)%>%dplyr::select(-functioning)
continuous_data<-data.frame(continuous_data)  
continuous_data <- continuous_data %>% subset(id %in% c(1,2,3,11,12,13,21,22,23)) 
continuous_data <- continuous_data %>% mutate(id = (1:9)[match(id, unique(continuous_data$id))])
continuous_data$label <- paste0("ID ", continuous_data$id, " (Dose ", continuous_data$dose, ")")

setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/for_paper")
pdf("Figure2B.pdf", width=12, height=6)
data.frame(continuous_data) %>% ggplot(aes(x=as.factor(timepoint), y=score, group=id,colour=label,
                                           linetype = label))+
  geom_line(size = 1)+
  scale_colour_manual(values = rep(c("#66C2A5","#8DA0CB","#FC8D62"), each=3),  guide = "none" )+
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotdash"), times=3),  guide = "none" )+
  theme_minimal(base_size = 16)+ylim(c(0,100))+
  geom_point()+
  labs(
    y="EORTC QLQ-C30 global health status score",
    x = "Timepoint",         
    color = "Patient (Dose)"       
  )+
  guides(colour = guide_legend(override.aes = list(linetype = rep(c("solid", "dashed", "dotted"), times = 3))))+
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0, yend = 100), 
               colour = "darkgrey", size = 0.8, arrow = arrow(type="closed",
                                                              length = unit(0.1, "inches")))+ylim(c(0,100))+
  scale_x_discrete(expand = expansion(add = c(1, 0)))+
  geom_text(aes(x = 0.4, y = 90, label = "Better \n tolerability"), col="black")
dev.off()

