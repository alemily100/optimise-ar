setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/paper")
library(RColorBrewer)
library(tidyverse)
library(dfcrm)
library(cowplot)
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

trial<-trial_sim_original(u, v, sample,no_enrolled, phi, true_tox_c, true_tox_p, no.dosages, target_clin, target_pat,stop_prob, a, b, target_stop)
study<-data.frame(trial[[1]])
study$x<- seq(from =0.5, length.out=15, by=0.5)
colnames(study)<- c("pat", "dose", "cdlt", "pdlt", "loc")

study_dlt<- study%>%filter(cdlt==1 |pdlt==1)

study$status <- with(study, 
                     ifelse(cdlt == 1 & pdlt == 1, "Both",
                            ifelse(cdlt == 1, "CDLT only",
                                   ifelse(pdlt == 1, "PDLT only", "None"))))

layout(matrix(c(1,2,1,3), nrow=2))
par(mar=c(6,6,1,2)+0.1)
p<-ggplot(study, aes(loc, dose,  label = pat), show.legend = TRUE) +
  annotate('rect', xmin=0, xmax=8, ymin=1.5, ymax=2.5, alpha=.2, fill='lightgreen')+
  geom_point(size = 10, shape = 21, fill = "white", colour="black", stroke=2) +
  geom_point(data=study, aes(loc, dose, fill=status), size = 10, shape = 21, stroke=2)+
  geom_text(data=study, aes(label = pat), vjust = 0.4)+
  xlab("Cohort")+ ylab("Dose")+scale_x_continuous(breaks=seq(from=1, by=1.5, length.out=5), labels=1:5)+theme_bw()+ labs(caption="(a) Patient dose assignment in PRO-CRM trial case study with observed Clinician-DLT and Patient-DLT outcomes.")+ theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),
                                                                                                                                                                                                                                       axis.title.y=element_text(size=15),panel.grid.major = element_blank(),
                                                                                                                                                                                                                                       panel.grid.minor = element_blank(),legend.text = element_text(size = 15),
                                                                                                                                                                                                                                       panel.border = element_blank(),legend.title = element_text(size = 15),
                                                                                                                                                                                                                                       panel.background = element_blank(), plot.caption=element_text(size=15,hjust=0))+
  scale_fill_manual(
    name = "DLT observation",
    values = c(
      "CDLT only" = "yellow",
      "PDLT only" = "red",
      "Both" = "orange",
      "None" = "white"
    ),breaks = c("CDLT only", "PDLT only", "Both"),
    labels = c(
      "CDLT only" = "C-DLT only",
      "PDLT only" = "P-DLT only",
      "Both" = "Both C-DLT & P-DLT",
      "None" = NULL   
    ))+
  scale_y_continuous(breaks=seq(from=1, by=1, length.out=3), labels=1:3)
p

est_c<-optimise(like_c,c(0,10), u=c(0.06, 0.14, 0.25), data_frame=study, maximum = TRUE)
est_p<-optimise(like_p,c(0,10), v=c(0.10, 0.21, 0.35), data_frame=study, maximum = TRUE)
#dose-finding decision 
p_c<-c(0.06, 0.14, 0.25)^est_c$maximum 
p_p<-c(0.10, 0.21, 0.35)^est_p$maximum
#pl<-data.frame(x,y)
#pl15<- data.frame(x, sapply(x, function (k) uniroot(f, a=15, pat_prob=k, pat_target=0.35, clin_target=0.25,interval= c(1.e-14, 1e04),
#                                                    extendInt="yes")$root))
point<- data.frame(p_p, p_c, c("Dose 1", "Dose 2", "Dose 3"))
colnames(point)[3]<- "label"
#colnames(pl)<-c("x_val", "y_val")
#colnames(pl15)<-c("x_val", "y_val")
l<-ggplot()+ geom_point(data=point, aes(x=p_p, y=p_c), size=4, pch=4, stroke=1.5)+geom_text(data=point,aes(x = p_p, y = p_c,label=label), vjust=-2)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15), plot.caption=element_text(size=15, hjust=0),
                   axis.title.y=element_text(size=15))+labs(x="Estimated P-DLT rate", y="Estimated C-DLT rate")+labs(caption="(b) Estimated P-DLT and C-DLT rate at final analysis.")+
  geom_vline(xintercept=0.35, col="red", lwd=1.5)+geom_hline(yintercept=0.25, col="blue", lwd=1.5, lty=2) + xlim(0,0.5) + ylim(0,0.5)


setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/OPTIMISE-AR (PRO Guidance paper)/Aim 2/emily-optimisear-generate-recommendations/for_paper")
pdf("case_study2_Figure6.pdf", width=13, height=8)
plot_grid(
  p, l,
  ncol = 1,               # stack vertically
  rel_heights = c(0.6, 0.4)   # adjust as needed
)
dev.off()