setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations")
source("functions_synthesise_data.R")

### CONTINUOUS DATA GENERATION 
#systematic review where study was identified: https://www.sciencedirect.com/science/article/pii/S0959804921011485#bib17
#CheckMate066: https://pubmed.ncbi.nlm.nih.gov/27405322/
CM <- read.csv("CheckMate066.csv")
#we utilise data for Dacarbazine only
CM<-data.frame(CM) %>% filter(Drug == "D")
n.pat<-rep(10, times=3)
n.dose<-3
pat.sd<-5
baseline_mean<- 66.2
baseline_sd<- 25.1
mean_change<- CM[,3]
se_change<- CM[,4]
n.timepoints<- length(mean_change)
dose_mean_vector<-c(0,-5,-10)
functioning<-1:6
mean<-c(0.9, -2.7, 3.6, 5.3, 1.0, 0.3)
se<-c(3.4, 2.6, 3.9, 2.7, 3.2, 3.7)

set.seed(10014)
time_last_report<-pmin(round(rexp(sum(n.pat), log(2)/36),0), n.timepoints+1)

set.seed(1234)
continuous_data<-final_cont(n.dose,n.pat, pat.sd, baseline_mean, baseline_sd,  mean_change, se_change, n.timepoints, dose_mean_vector)
continuous_data<- cbind(continuous_data, rep(1, times=nrow(continuous_data)))
colnames(continuous_data)[length(colnames(continuous_data))]<- "functioning"
for(i in 2:length(functioning)){
  new<-final_cont(n.dose,n.pat, pat.sd, baseline_mean, baseline_sd,  mean_change/mean_change[length(mean_change)]*mean[i], se_change/se_change[length(se_change)]*se[i], n.timepoints, dose_mean_vector)
  new<- cbind(new, rep(i, times=nrow(new)))
  colnames(new)[length(colnames(new))]<- "functioning"
  continuous_data<-
    rbind(continuous_data,new)
}

continuous_data<-data.frame(continuous_data)
id_list <- unique(continuous_data$id)
functioning_list <- unique(continuous_data$functioning)

time_last_df <- tibble(
  id = id_list,
  time_last_report = time_last_report
)

continuous_data <- continuous_data %>%
  left_join(time_last_df, by = c("id"))

continuous_data<- continuous_data %>%
  filter(timepoint <= time_last_report)%>%dplyr::select(-time_last_report)

missing<-rbinom(nrow(continuous_data), 1, 0.06)
continuous_data$score[missing==1]<- NA
write.csv(continuous_data, "dummy_continuous.csv")
