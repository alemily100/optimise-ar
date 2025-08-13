library(truncnorm)
library(tidyverse)
library(MASS)

#synthesise continuous data
baseline<- function (total.pat, baseline_mean, baseline_sd, ran_eff){
  return(sapply(1:total.pat, function (k) rtruncnorm(1,a=0,b=100, baseline_mean, baseline_sd)+ran_eff[k]))
}

at_time<- function(dose, n.pat, mean_change, se_change, baseline_val, ran_eff, time_point){
  pat.id<- (((dose-1)*n.pat)+1):(n.pat*dose)
  at_timepoint<-baseline_val[pat.id]+sapply(pat.id, function (k) rnorm(1, mean_change[time_point], se_change[time_point])+ran_eff[k])
  return(at_timepoint)
}

over_time<- function(dose, n.pat, mean_change, se_change, baseline_val, ran_eff, n.timepoints){
  pat.id<- (((dose-1)*n.pat)+1):(n.pat*dose)
  function_matrix<-sapply(1:n.timepoints, function(j) at_time(dose, n.pat, mean_change, se_change, 
                                                              baseline_val, ran_eff,j))
  new_mat<-cbind(rep(pat.id, times=n.timepoints), c(function_matrix), rep(1:n.timepoints, each=n.pat), rep(dose, times=n.pat*n.timepoints))
  colnames(new_mat)<- c("id", "score", "timepoint", "dose")
  return(new_mat)
}
  

final_cont<- function(n.dose, n.pat, pat_sd, baseline_mean, baseline_sd, 
                      mean_change, se_change, n.timepoints, dose_mean_vector){
  total.pat<- sum(n.pat)
  random_eff<- rnorm(total.pat, 0, pat_sd)
  baseline_vals<- baseline(total.pat, baseline_mean,baseline_sd, random_eff)
  mat<- cbind(1:total.pat, baseline_vals, rep(0, times=total.pat), rep(1:n.dose,times=n.pat))
  colnames(mat)<-c("id", "score", "timepoint", "dose")
  for(i in 1:n.dose){
    new_mat<- over_time(i, n.pat[i], mean_change+dose_mean_vector[i], se_change, baseline_vals, random_eff, 7)
    colnames(new_mat)<-c("id", "score", "timepoint", "dose")
    mat<- rbind(mat, new_mat)
  }
  return(mat)
}

#synthesise ordinal data 
generate_matrix_beta<- function(ndose, ngrades, vec_a, a, b){
  u<-sort(runif(ngrades-1))
  dose1<- c(pbeta(u[1],a,b), pbeta(u[2],a,b)-pbeta(u[1],a,b), pbeta(u[3],a,b)-pbeta(u[2],a,b),pbeta(u[4],a,b)-pbeta(u[3],a,b),pbeta(u[4],a,b, lower.tail = FALSE))
  M<- rbind(dose1, t(sapply(1:(ngrades-1), function(k) dose_effect(a,b, dose1, vec_a[k]))))
  return(M)
}

dose_effect<- function(a,b, vec_initial_probs, modified_a){
  ngrades<- length(vec_initial_probs)
  C<- sapply(1:(ngrades-1), function(j) qbeta(sum(vec_initial_probs[1:j]), a, b))
  vec<- c()
  vec[1]<-pbeta(C[1], modified_a, b)
  for(i in 2:(ngrades-1)){
    vec[i]<- pbeta(C[i], modified_a, b)-pbeta(C[i-1],modified_a, b)
  }
  vec[ngrades]<- 1-pbeta(C[ngrades-1], modified_a,b) 
  return(vec)
}

generate_ls<- function(ndose, ngrades, n.toxicity,mod_a, mod_b){
  #M<- matrix(nrow=n.toxicity, ncol=5)
  ls<- vector(mode='list', length=n.toxicity)
  for(i in 1:n.toxicity){
    ls[[i]]<- generate_matrix_beta(ndose, ngrades, mod_a[-1],mod_a[1],4)
  }
  return(ls)
}

ae.score.cycle<- function(toxicity_matrix, a, vec_a,b, vec_modified_cycle_b){
  cyc<-lapply(c(b,vec_modified_cycle_b), function (k) time_effect(toxicity_matrix,b, vec_a,a,k))
  return(cyc)
}

time_effect<- function(matrix_tox, b, vec_modified_a, a, modified_b){
  vec<-c()
  ngrades=ncol(matrix_tox)
  dose<-0
  M<- matrix(nrow=nrow(matrix_tox), ncol=ncol(matrix_tox))
  while(dose< nrow(matrix_tox)){
    dose<- dose+1
    modified_a<-c(a, vec_modified_a)[dose] 
    C<- sapply(1:(ngrades-1), function(j) qbeta(sum(matrix_tox[dose,1:j]), modified_a, b))
    vec[1]<-pbeta(C[1], modified_a, modified_b)
    for(i in 2:(ngrades-1)){
      vec[i]<- pbeta(C[i], modified_a, modified_b)-pbeta(C[i-1], modified_a, modified_b)
    }
    vec[ngrades]<- 1-pbeta(C[ngrades-1], modified_a, modified_b) 
    M[dose,]<- vec
  }
  return(M)
}

ae_score_patient_wish<- function(cycle, dose, n.toxicity, toxicity_dose_grade_cycle_mat){
  sigma<- diag(15)
  z<- mvrnorm(mu=rep(0, times=n.toxicity), Sigma=sigma)
  u<-sapply(1:n.toxicity, function (k) pnorm(z[k], mean =0, sd=sqrt(diag(sigma)[k])))
  scores<- sapply(1:n.toxicity, function (j) min(which(sapply(1:5, function(k) sum(toxicity_dose_grade_cycle_mat[,j][[cycle]][dose,1:k]))>u[j])))
  return(scores)
}

at_time_ord<- function(dose, n.pat, timepoint, n.toxicity, prob_t){
  pat.id<- (((dose-1)*n.pat)+1):(n.pat*dose)
  M_score<-sapply(pat.id, function (k) ae_score_patient_wish(timepoint, dose, n.toxicity,prob_t))
  mat<-cbind(rep(pat.id,each=n.toxicity), rep(1:n.toxicity, times=n.pat), c(M_score), rep(timepoint-1, times=n.toxicity*n.pat), 
             rep(dose, times=n.toxicity*n.pat))
  colnames(mat)<-c("id","toxicity","grade", "timepoint", "dose")
  return(mat)
}

over_time_ord<- function(dose, n.pat, n.timepoints, n.toxicity, prob_t){
mat<-at_time_ord(dose, n.pat, 1, n.toxicity, prob_t)
for(k in 2:n.timepoints){
  mat<- rbind(mat, at_time_ord(dose, n.pat, k, n.toxicity, prob_t))
}
return(mat)
}

final_ord<- function(ndose, ngrades, mod_a, mod_b, n.pat, n.timepoints, n.toxicity){
matrix_tox<-generate_ls(ndose, ngrades, n.toxicity, mod_a, mod_b)
prob_t<-sapply(1:n.toxicity, function (k) ae.score.cycle(matrix_tox[[k]], mod_a[1], mod_a[-1],
                                                         mod_b[1], mod_b[-1]))  
mat<- over_time_ord(1, n.pat[1], n.timepoints, n.toxicity, prob_t)
for(j in 2:n.dose){
  mat<- rbind(mat, over_time_ord(j, n.pat[j], n.timepoints, n.toxicity, prob_t))
}
return(mat)
}





