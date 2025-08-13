library(dplyr)
library(tidyr)
library(ggplot2)
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")


#import dummy data
ordinal_data <- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/dummy_ordinal.csv")[,-1]

#create plot
ordinal_data<-ordinal_data%>%group_by(id)%>%
  mutate(has_tp = any(timepoint == 7)) %>%
  filter(timepoint %in% c(0, 7))%>%
  pivot_wider(names_from = timepoint, values_from = grade, names_prefix = "tp_") %>%
  mutate(grade_change = case_when(
    has_tp == FALSE ~ "ICE" ,
    TRUE ~  as.character(tp_0-tp_7)))


ordinal_data$grade_change[is.na(ordinal_data$grade_change)==TRUE]<- "Missing"
prop<- ordinal_data%>%group_by(toxicity,dose)%>%
  count(grade_change)%>% 
  ungroup()%>%
  mutate(original = TRUE)%>%
  complete(toxicity, dose, grade_change, fill = list(n = 0))%>%
  group_by(toxicity, dose) %>%
  mutate(n.pat = sum(n[grade_change != "ICE"], na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(proportion = n / n.pat)%>%dplyr::select(-n)

prop$toxicity<- as.factor(prop$toxicity)
prop$dose<- as.factor(prop$dose)



#Final figure 
prop$dose <- factor(prop$dose, levels = rev(levels(prop$dose)))
#prop$base_adj <- factor(prop$base_adj, levels = c("Unreported",5,4,3,2,1))
prop$grade_change <- as.character(prop$grade_change)
#prop <- prop %>%
#  replace_na(list(grade_change = "Unreported"))

missing<- prop%>%subset(grade_change=="Missing")%>%mutate(label = paste0("Unreported: ",round(proportion,2)*100,"%"))

missing<-cbind(missing, rep(0.7, times=nrow(missing)))
colnames(missing)[8]<- "loc"
prop <- prop %>%
  mutate(proportion = ifelse(grade_change>=0, proportion, -proportion))


prop<-prop%>%subset(grade_change!="Missing")%>%subset(grade_change!=0)
prop <- prop %>% filter(original)
val<-prop %>% 
  group_by(dose, toxicity) %>%
  arrange(dose, toxicity, grade_change) %>%
  mutate(percent = paste0(round(abs(proportion) * 100, 1), "%")) %>%
  
  # Create separate cumsum columns for pos/neg
  mutate(
    cumsum_pos = if_else(proportion > 0, cumsum(if_else(proportion > 0, proportion, 0)), NA_real_),
    cumsum_neg = if_else(proportion < 0, cumsum(if_else(proportion < 0, proportion, 0)), NA_real_)
  ) %>%
  mutate(
    pos = case_when(
      proportion > 0 ~ cumsum_pos - 0.5 * proportion,
      proportion < 0 ~ cumsum_neg - 0.5 * proportion,
      TRUE ~ 0  # fallback for proportion == 0
    )
  )


facet_labels<- c("1"="Symptomatic Adverse Event 1","2"="Symptomatic Adverse Event 2", "3"="Symptomatic Adverse Event 3", "4"="Symptomatic Adverse Event 1", "5"="Symptomatic Adverse Event 2", "6"="Symptomatic Adverse Event 3", "7"="Symptomatic Adverse Event 7", "8"="Symptomatic Adverse Event 8", "9"="Symptomatic Adverse Event 9", "10"="Symptomatic Adverse Event 10")
val<- val%>%subset(grade_change!="ICE")
val$grade_change <- factor(val$grade_change)
levels(val$dose) <-  c("Dose 3", "Dose 2", "Dose 1")

#arrowdf <- tibble(toxicity = "10")

arrowdf <- tibble(
  toxicity = factor("6", levels = levels(val$toxicity)),  # important!
  x = 0.5,
  xend = 0.5,
  y = 0,
  yend = 1
)

arrowdf.rev <- tibble(
  toxicity = factor("6", levels = levels(val$toxicity)),  # important!
  x = 0.5,
  xend = 0.5,
  y = 0,
  yend = -1
)

levels(missing$dose) <-  c("Dose 3", "Dose 2", "Dose 1")

ice_n_pat <- prop %>%
  distinct(dose, n.pat)


labels<- cbind(c(-5, -4, -3, -2, -1, 1, 2, 3, 4, 5),
                   c("red", "deeppink", "orange","gold","lightyellow3", "lightblue", "darkseagreen3", "darkolivegreen3", "chartreuse3", "chartreuse4"),
                   rev(c("Improve by 5", "Improve by 4", "Improve by 3","Improve by 2", "Improve by 1", "Worsen by 1", 
                  "Worsen by 2", "Worsen by 3", "Worsen by 4", "Worsen by 5")))

colnames(labels)<- c("score", "colour", "leg")
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/for_paper")

pdf("Figure3A.pdf", height=8, width=15)
ggplot(val%>%subset(toxicity %in% 4:6),aes(x = dose, y = proportion, fill = grade_change, group=dose)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_y_continuous(labels = abs) +
  theme_minimal(base_size=15)+theme(
    strip.text.y = element_text(angle = 0))+
  facet_grid(toxicity ~ .,scales = "free_y", labeller=as_labeller(facet_labels)) + ylim(c(-1,1))+
  scale_fill_manual(values=labels[match(as.numeric(levels(val$grade_change)), labels[,1]),2],
                    labels=labels[match(as.numeric(levels(val$grade_change)), labels[,1]),3])+
  geom_text(aes(y = pos, label = percent), size = 3, color = "black")+
  labs(x = paste0("Dose (N=", paste(ice_n_pat$n.pat, collapse = "; "), ")"),                       
       y = "Percentage (%)",                       
       fill = "Change in Severity ")+
  scale_y_continuous(breaks = c(0),
                     limits=c(-1,1))+
  geom_hline(yintercept = 0, color = "grey", size = 1)+
  geom_segment(data = arrowdf, 
               aes(x = x, xend = xend, y = y, yend = yend), 
               colour = "darkgrey", 
               size = 0.8, 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               inherit.aes = FALSE)+
  scale_x_discrete(expand = expansion(mult = c(1.5, 0)))+
  geom_segment(data = arrowdf.rev, 
               aes(x = x, xend = xend, y = y, yend = yend), 
               colour = "darkgrey", 
               size = 0.8, 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               inherit.aes = FALSE)+
  geom_text(data = arrowdf, aes(x = -0.25, y = 0.5, label = "Improvement"),
            inherit.aes = FALSE)+
  geom_text(data = arrowdf.rev, aes(x = -0.25, y = -0.5, label = "Worsening"),
            inherit.aes = FALSE)+
  geom_text(
    data = missing%>%filter(toxicity %in% 4:6),
    aes(x = dose, y = loc, label = label, group = dose),  # 'toxicity' must also be in 'missing'
    inherit.aes = FALSE,
    vjust = -0.5
  )
dev.off()


