setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")

continuous_data <- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/dummy_continuous.csv")[,-1]

#### FIGURE GENERATION
continuous_data<- continuous_data%>%filter(functioning==1)%>%dplyr::select(-functioning)
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/for_paper")
pdf("Figure2D.pdf", width=10, height=6)
ggplot(continuous_data, aes(x = as.factor(timepoint),y = score, fill = as.factor(dose))) + 
          geom_boxplot() + 
          scale_fill_manual(
                      values = c("1" = "#66C2A5",   # Dose level 1
                                 "2" = "#8DA0CB",   # Dose level 2
                                 "3" = "#FC8D62"    # Dose level 3
                                 ),
                      name=element_blank(),label = c("Dose 1", "Dose 2", "Dose 3"))+
          theme_minimal(base_size=14)+
          xlab("Timepoint") + ylab("EORTC QLQ-C30 global health status score")+
  geom_segment(aes(x = 0.2, xend = 0.2, y = 0, yend = 100), 
               colour = "darkgrey", size = 0.8, arrow = arrow(type="closed",
                                                              length = unit(0.1, "inches")))+ylim(c(0,100))+
  scale_x_discrete(expand = expansion(add = c(2, 0)))+
  geom_text(aes(x = -0.45, y = 90, label = "Better \n tolerability"))
dev.off()
                 