setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")
source("functions_synthesise_data.R")

### CONTINUOUS DATA GENERATION 

n.toxicity<-15
mod_a<- seq(from=0.2, to=1.2, length.out=3)
mod_b<- seq(from=4, to=1.5, length.out=8)
n.pat<- rep(10, times=3)
n.grade<- 5
n.dose<-length(mod_a)
n.timepoints<-length(mod_b)

set.seed(10029)
time_last_report<-pmin(round(rexp(sum(n.pat), log(2)/36),0), n.timepoints)
ordinal_data<-final_ord(n.dose, n.grade,mod_a, mod_b, n.pat, n.timepoints, n.toxicity)


## add non-treatment related death to dataset 
time_df <- data.frame(id = 1:sum(n.pat), time_last_report = time_last_report)
ordinal_data <- data.frame(ordinal_data) %>%
  left_join(time_df, by = "id") %>%
  filter(timepoint <= time_last_report)
missing<-rbinom(nrow(ordinal_data), 1, 0.05)
ordinal_data$grade[missing==1]<- NA
ordinal_data <- ordinal_data%>% dplyr::select(-time_last_report)
write.csv(ordinal_data, "dummy_ordinal.csv")
