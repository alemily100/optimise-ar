library(RColorBrewer)
library(tidyverse)
library(dfcrm)

### FUNCTIONS
like_c<- function(gamma, data_frame,u){
  l<-0
  outcome<-data_frame[,3]
  dose<- data_frame[,2]
  for(i in 1:nrow(data_frame)){
    l<-l+(as.numeric(outcome[i]==0))*log(1-u[dose[i]]^gamma)+(as.numeric(outcome[i]==1)*gamma)*log(u[dose[i]])
  }
  l
}

like_p<- function(eta, data_frame,v){
  l<-0
  outcome<-data_frame[,4]
  dose<- data_frame[,2]
  for(i in 1:nrow(data_frame)){
    l<-l+(as.numeric(outcome[i]==0))*log(1-v[dose[i]]^eta)+(as.numeric(outcome[i]==1)*eta)*log(v[dose[i]])
  }
  l
}

time_to_dlt<- function(dose,clin_tox, pat_tox, phi){
  lambda_c<- -log(1-clin_tox[dose])
  lambda_p<- -log(1-pat_tox[dose])
  u1<-runif(1)
  time_c<- (-log(u1))/lambda_c
  u2<-runif(1)
  a<- u2/(u1^(-(phi+1)/phi))
  b<- u1^(-1/phi)-1
  time_p<- (phi/lambda_p)*log(a^(1/(-phi-1))-b)
  times<- c(time_c, time_p)
  return(times)
}

truncated_matrix<- function(no_patients, matrix_time_to_toxicity, previous_matrix, dose){
  new_mat<- matrix(nrow=no_patients, ncol=4)
  new_mat[,1]<-(max(previous_matrix[,1])+1):(max(previous_matrix[,1])+no_patients)
  new_mat[,2]<-rep(dose, times=no_patients)
  new_mat[,3]<-as.numeric(matrix_time_to_toxicity[,1]<=1)
  new_mat[,4]<-as.numeric(matrix_time_to_toxicity[,2]<=1)
  return(rbind(previous_matrix, new_mat))
}

stopping_rule<-function(matrix, a,b, target, min_p){
  m_1<- matrix[matrix[,2]==1,]
  sum_dlt<- sum((m_1[,3] | m_1[,4]))
  a_tilde<- a + sum_dlt
  b_tilde<- b + nrow(m_1) - sum_dlt
  ifelse(pbeta(target, a_tilde, b_tilde, lower.tail = FALSE)>min_p,
         "stop", "continue")
}

trial_sim_original<- function(u_skeleton, v_skeleton, sample,no_enrolled, phi, true_tox_c, true_tox_p, number_dosages, target_clin, target_pat,stopping_rule_prob, a, b, target){
  current_dose<-1
  val_mat<- NULL
  val<-t(sapply(rep(current_dose, times=3), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
  M<- matrix(nrow=3, ncol=4)
  M[,1]<-1:3
  M[,2]<-rep(current_dose, times=3)
  M[,3]<-val[,1]<=1
  M[,4]<-val[,2]<=1
  val_mat<-rbind(val_mat, val)
  while(sum(M[,3], M[,4])==0 && nrow(M)<sample){
    #use rule based design when there is no heterogeneity in either endpoint 
    if(stopping_rule(M, a,b, target, stopping_rule_prob)=="stop"){
      current_dose<- NA
      break
    }
    ifelse(current_dose<number_dosages,current_dose<-current_dose+1, current_dose<-number_dosages)
    val<-t(sapply(rep(current_dose, times=3), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
    val_mat<-rbind(val_mat, val)
    M<-truncated_matrix(no_enrolled, val, M, current_dose)
  }
  while(xor(sum(M[,3])>0,sum(M[,4])>0) & nrow(M)<sample){
    #use rule based design and likelihood CRM when there is heterogeneity in one endpoint only  
    if(stopping_rule(M, 0.1,0.9, target, stopping_rule_prob)=="stop"){
      current_dose<- NA
      break
    }
    if(sum(M[,3])>0){
      est_c<-optimise(like_c,c(0,10), data_frame=M,u=u_skeleton, maximum = TRUE)
      p_c<-u_skeleton^est_c$maximum 
      current_dose<-min(which.min(abs(p_c-target_clin)), current_dose+1)
      val<-t(sapply(rep(current_dose, times=no_enrolled), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
      val_mat<-rbind(val_mat, val)
      M<-truncated_matrix(no_enrolled, val, M, current_dose)
    }
    else{
      est_p<-optimise(like_p, c(0,10), data_frame=M,v=v_skeleton,maximum = TRUE)
      p_p <- v_skeleton^est_p$maximum
      current_dose<-min(which.min(abs(p_p-target_pat)), current_dose+1)
      val<-t(sapply(rep(current_dose, times=no_enrolled), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
      val_mat<-rbind(val_mat, val)
      M<-truncated_matrix(no_enrolled, val, M, current_dose)
    }
  }
  while(sum(M[,3])>0&sum(M[,4])>0&nrow(M)<sample){
    #use likelihood CRM when there is heterogeneity in both endpoints   
    if(stopping_rule(M, a,b, target, stopping_rule_prob)=="stop"){
      current_dose<- NA
      break
    }
    est_c<-optimise(like_c,c(0,10), u=u_skeleton, data_frame=M, maximum = TRUE)
    est_p<-optimise(like_p,c(0,10), v=v_skeleton, data_frame=M, maximum = TRUE)
    #dose-finding decision 
    p_c<-u_skeleton^est_c$maximum 
    p_p<-v_skeleton^est_p$maximum
    rec_dose<-min(which.min(abs(p_c-target_clin)), which.min(abs(p_p-target_pat)))
    ifelse(rec_dose>max(M[,2]), current_dose<-max(M[,2])+1, current_dose<-rec_dose)
    val<-t(sapply(rep(current_dose, times=no_enrolled), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
    val_mat<-rbind(val_mat, val)
    M<-truncated_matrix(no_enrolled, val, M, current_dose)
  }
  est_c<-optimise(like_c,c(0,10), u=u_skeleton, data_frame=M, maximum = TRUE)
  est_p<-optimise(like_p,c(0,10), v=v_skeleton, data_frame=M, maximum = TRUE)
  #dose-finding decision 
  p_c<-u_skeleton^est_c$maximum 
  p_p<-v_skeleton^est_p$maximum
  rec_dose<-min(which.min(abs(p_c-target_clin)), which.min(abs(p_p-target_pat)))
  if(is.na(current_dose)==FALSE){
    ifelse(rec_dose>max(M[,2]), current_dose<-max(M[,2])+1, current_dose<-rec_dose)
  }
  unlink(file.path("C:/Users/ealger/AppData/Local/Temp", "Rtmp*"), recursive = T)
  optimal_dose<- min(which.min(abs(true_tox_c-target_clin)), which.min(abs(true_tox_p-target_pat)))
  mtd_assigned<-sum(M[,2]==optimal_dose)/sample
  overdosed<- sum(M[,2]>optimal_dose)/sample
  return(list(M, val_mat))
}

### DATA GENERATION
u<- c(0.06, 0.14, 0.25, 0.38, 0.50)
v<- c(0.10, 0.21, 0.35, 0.49, 0.61)
target_c<- 0.25
target_p<- 0.35
sample<- 15
no_enrolled<-3
phi<-0.9
true_tox_c<-c(0.05,0.25,0.40)
true_tox_p<-c(0.10,0.15,0.35)
no.dosages<- length(true_tox_c)
target_clin<- 0.25
target_pat<- 0.35
stop_prob<-0.8
target_stop<- 0.5
all_sc<-NULL
a<- 0.35
b<- 0.65

set.seed(130)
study<-trial_sim_original(u, v, sample,no_enrolled, phi, true_tox_c, true_tox_p, no.dosages, target_clin, target_pat,stop_prob, a, b, target_stop)

### FIGURE GENERATION
time.dlt<- data.frame(study[[2]])
times<- data.frame(study[[1]])
times$start<-c(sapply(1:5, function(k) ((k-1)*5):(5*k-3)))
times$end<-c(sapply(1:5, function(k) ((5*k-2):(5*k))))
colnames(times)<-c("id", "dose", "cdlt", "pdlt", "start", "end")
cdlt_elements<-(times%>%filter(cdlt!=0))[,1]
pdlt_elements<-(times%>%filter(pdlt!=0))[,1]
times[cdlt_elements,3]<-(time.dlt[cdlt_elements,1]*3)+times[cdlt_elements,5]
times[pdlt_elements,4]<-(time.dlt[pdlt_elements,2]*3)+times[pdlt_elements,5]
times[(times%>%filter(cdlt==0))[,1],3]<-NA
times[(times%>%filter(pdlt==0))[,1],4]<-NA
times$id <- factor(times$id, levels = unique(times$id))

setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/for_paper")
pdf("case_study2_Figure8.pdf", width=12, height = 5)
ggplot(times, aes(y = id)) +
  geom_segment(aes(x = start, xend = end, yend = id, color = factor(dose)), size = 3) +
  geom_point(aes(x = cdlt, fill = "Clinician-DLT observation"), shape = 21, color = "black", size = 4, stroke = 1) +
  geom_point(aes(x = pdlt, fill = "Patient-DLT observation"), shape = 21, color = "black", size = 4, stroke = 1) +
  labs(x = "Trial timeline (weeks)", y = "Patient", color = "Dose", fill = "DLT Type") +
  scale_color_manual(
    values = c(
      "1" = "#66C2A5",   # Dose level 1
      "2" = "#8DA0CB",   # Dose level 2
      "3" = "#FC8D62"    # Dose level 3
    ),
    labels = c("1" = "Dose level 1", "2" = "Dose level 2", "3" = "Dose level 3")
  ) +
  scale_fill_manual(
    values = c("Clinician-DLT observation" = "yellow", "Patient-DLT observation" = "red")
  ) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
  scale_y_discrete(labels = paste("Patient", seq_along(levels(times$id)))) +
  theme_minimal(base_size = 14)+
  scale_x_continuous(breaks = 1:25)+
  theme(
    panel.grid.minor = element_blank())
dev.off()

### Credible intervals 
ci.data<-data.frame(study[[1]])
colnames(ci.data)<- c("id", "dose", "cdlt", "pdlt")
ci.data%>%filter(dose==1)

c_dlt<-crm(c(00.06, 0.14, 0.25), 0.25, study[[1]][,3], study[[1]][,2],
    conf.level = 0.9, method = "mle",
    model = "empiric", model.detail = TRUE,
    patient.detail = TRUE, var.est = TRUE)

p_dlt<-crm(c(0.10, 0.21, 0.35), 0.35, study[[1]][,4], study[[1]][,2],
    conf.level = 0.9, method = "mle",
    model = "empiric", model.detail = TRUE,
    patient.detail = TRUE, var.est = TRUE)
min(c_dlt$mtd, p_dlt$mtd)
