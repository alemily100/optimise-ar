setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")


#import dummy data
ordinal_data <- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/dummy_ordinal.csv")[,-1]

#create plot
ordinal_data<-ordinal_data%>%group_by(id)%>%
  mutate(has_tp = any(timepoint == 7)) %>%
  filter(timepoint %in% c(0, 7))%>%
  pivot_wider(names_from = timepoint, values_from = grade, names_prefix = "tp_") %>%
  mutate(grade_change = tp_7 - tp_0)%>%
  mutate(base_adj = case_when(
    has_tp == FALSE ~ "ICE" ,
    TRUE ~  ifelse(grade_change>0, as.character(tp_7), as.character(1))))

ordinal_data$base_adj[is.na(ordinal_data$base_adj)==TRUE]<- "Unreported"

# Data Manipulation 
#prop<- ordinal_data%>%group_by(toxicity,dose)%>%count(base_adj) %>%
#  mutate(proportion = n / n.pat)%>%dplyr::select(-n)

prop<- ordinal_data%>%group_by(toxicity,dose)%>%
  count(base_adj)%>% 
  ungroup()%>%
  mutate(original = TRUE)%>%
  complete(toxicity, dose, base_adj, fill = list(n = 0))%>%
  group_by(toxicity, dose) %>%
  mutate(n.pat = sum(n[base_adj != "ICE"], na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(proportion = n / n.pat)%>%dplyr::select(-n)

prop$toxicity<- as.factor(prop$toxicity)
prop$base_adj<- as.factor(prop$base_adj)
prop$dose<- as.factor(prop$dose)
levels(prop$base_adj) <- c(levels(prop$base_adj), "Unreported")


#Final figure 
prop$dose <- factor(prop$dose, levels = rev(levels(prop$dose)))
prop<-prop%>%subset(base_adj != "ICE")
prop$base_adj <- factor(prop$base_adj, levels = c("Unreported",5,4,3,2,1))
facet_labels<- c("1"="Symptomatic Adverse Event 1","2"="Symptomatic Adverse Event 2", "3"="Symptomatic Adverse Event 3", "4"="Symptomatic Adverse Event 1", "5"="Symptomatic Adverse Event 2", "6"="Symptomatic Adverse Event 3", "7"="Symptomatic Adverse Event 7", "8"="Symptomatic Adverse Event 8", "9"="Symptomatic Adverse Event 9", "10"="Symptomatic Adverse Event 10")

ice_n_pat <- prop %>%
  distinct(dose, n.pat)

setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/paper/optimise-ar/for_paper")
pdf("Figure3B.pdf", height=4, width=7)
prop%>%subset(toxicity %in% 4:6)%>% ggplot(aes(fill=base_adj,x=dose, y=proportion*100)) + 
  geom_bar(position="stack", stat="identity") + coord_flip()+
  facet_grid(toxicity ~ .,scales = "free_y", labeller = as_labeller(facet_labels))+
  scale_fill_manual(values=rev(c("darkgreen", "lightgreen","yellow","orange", "red", "darkgrey")), labels =rev(c("None", "Mild", "Moderate", "Severe", "Very Severe", "Unreported")),
                    guide = guide_legend(reverse = TRUE))+
  labs(x = paste0("Dose (N=", paste(ice_n_pat$n.pat, collapse = "; "), ")"),                       
       y = "Percentage (%)",                       
       fill = "Severity")+
  theme(
    strip.text.y = element_text(angle = 0))
dev.off()



