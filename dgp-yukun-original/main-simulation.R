library(trust)
library(Matrix)
library(Rlab)
library(MASS)
##################################
######## Set parameters ##########
##################################
alpha=0.8
n=1000                         # number of obs
NUM_ITERATIONS = 2000;         # number of iteration
dgp=4;                         # dgp list
h=0.05;                        # bandwidth
sieve_regularization = 0;
K_con=5;                       #seive order K and k
k <- 5;
DELTA = 0.00000000001;         # precise of derivatiion
dimX = 6;
true_list = c(1.7387,1.47289,2.8483,3.2091); #alpha=0.8
true<-true_list[dgp];
##################################
######## generate data ###########
##################################

RESULTS <- NULL;
AllRESULTS <- NULL;
source("~/Desktop/Weak-DID/dgp-yukun-original/DGP_yang_ding.R")




##################################################
###define functions ###############################
##################################################


alpha2 <- function(d,Y1,Y0,gamma1,gamma2,x,k,mm,j,h){
  b=0.01
  gamma_ps = matrix(gamma1,,1)
  gamma_reg = matrix(gamma2,,1)
  reg = x%*%gamma_reg
  ps  = exp(x%*%gamma_ps)/(1+exp(x%*%gamma_ps))
  
  integral =0
  plist = 0:100/100
  for (p in plist) {
    if (p<=1-h){
      integral = integral+p/(1-p)*mean((1-d)*(Y1-Y0-reg)*dnorm((ps-p)/b)/b)
    }
    else{
        for (kappa in 1:k) {
          integral=integral+(1-p)^(kappa-1)*mean(dnorm((ps-p)/b)/b)/factorial(kappa)*mm[kappa,j]
        }
      }
    }
  
  return( integral * (plist[2]-plist[1]) )
}


##################################################
###define P_K(0)##################################
##################################################

a = 0;
dpda1 = matrix(c(0, 3^0.5*2, 5^0.5*(2*6*a^1-6), 7^0.5*(3*20*a^2-2*30*a^1+12), 9^0.5*(4*70*a^3-3*140*a^2+2*90*a^1-20), 11^0.5*(5*252*a^4-4*630*a^3+3*560*a^2-2*210*a^1+30)), 6,1);
dpda2 = matrix(c(0,       0,       5^0.5*(2*6),      7^0.5*(2*3*20*a^1-2*30),    9^0.5*(3*4*70*a^2-2*3*140*a^1+2*90),  11^0.5*(4*5*252*a^3-3*4*630*a^2+2*3*560*a^1-2*210)), 6,1);
dpda3 = matrix(c(0,       0,                 0,               7^0.5*(2*3*20),           9^0.5*(2*3*4*70*a^1-2*3*140),        11^0.5*(3*4*5*252*a^2-2*3*4*630*a^1+2*3*560)), 6,1);
dpda4 = matrix(c(0,       0,                 0,                            0,                       9^0.5*(2*3*4*70),                  11^0.5*(2*3*4*5*252*a^1-2*3*4*630)), 6,1);
dpda5 = matrix(c(0,       0,                 0,                            0,                                      0,                                11^0.5*(2*3*4*5*252)), 6,1);


for(iter in 1:NUM_ITERATIONS){
  
  
  ##*****remember to add  dpg=1:4
  #data.did<-dgps(dgp,n)
  data.did<-dgps(n,dgp,alpha)
  y0 <- data.did$y0
  y1 <- data.did$y1
  d <-  data.did$d
  x <-  data.did$X
  deltaY <- y1 - y0;
  int.cov <-cbind(1,x)

  ##########################################################################################
  ###reg coefficients for selection model and outcome regression model######################
  ##########################################################################################
  
  i.weights <- as.vector(rep(1, n))
  reg_result = lm( I(deltaY) ~ int.cov[,2:7], subset=(d==0) )
  ps_result  = glm( d ~ int.cov[,2:7], family=binomial(link='logit') )
  gamma_reg_hat = matrix(coef(reg_result),,1)
  gamma_ps_hat  = matrix(coef(ps_result),,1)
  gamma_hat = rbind(gamma_ps_hat,gamma_reg_hat)
  reg_hat = int.cov%*%gamma_reg_hat
  ps.fit  = exp(int.cov%*%gamma_ps_hat)/(1+exp(int.cov%*%gamma_ps_hat))
  

  ########################
  ###define A and B#######
  ########################
  A <- 1 - ps.fit;
  A2 <- matrix(A,n,1);
  B1 <- matrix(d*(deltaY-reg_hat),n,1);
  B2 <- ps.fit * (1-d) * (deltaY-reg_hat);
  B2 <- matrix (B2,n,1);
  B3 <- matrix (d,n,1);
  
  
  
  ##############################################################
  #########influence function for gamma.   #####################
  ##############################################################

  
  phi1_den=0
  phi1_num=NULL
  phi2_den=0
  phi2_num=NULL
  for( idx in 1:n ){
    phi1_den = phi1_den + matrix(int.cov[idx,],,1)%*%(ps.fit[idx]*(1-ps.fit[idx]))%*%matrix(int.cov[idx,],1,)
    phi1_num = cbind(phi1_num, matrix(int.cov[idx,],,1)%*%(d[idx]-ps.fit[idx]))
    phi2_den = phi2_den + matrix(int.cov[idx,],,1)%*%(1-d[idx])%*%matrix(int.cov[idx,],1,)
    phi2_num = cbind(phi2_num, matrix(int.cov[idx,],,1)%*%((1-d[idx])*(deltaY[idx])))
  }
  phi1_hat = ginv(phi1_den) %*% phi1_num
  phi2_hat = ginv(phi2_den) %*% phi2_num
  phi_hat = rbind(phi1_hat,phi2_hat)
  

  ###first def gamma_diff############
  
  gamma_ps_hat_diff<-matrix(0,7,7);
  
  gamma_ps_hat_diff[1,]=c(gamma_ps_hat[1]+DELTA,gamma_ps_hat[2:7]);
  gamma_ps_hat_diff[2,]=c(gamma_ps_hat[1],gamma_ps_hat[2]+DELTA,gamma_ps_hat[3:7]);
  gamma_ps_hat_diff[3,]=c(gamma_ps_hat[1:2],gamma_ps_hat[3]+DELTA,gamma_ps_hat[4:7]);
  gamma_ps_hat_diff[4,]=c(gamma_ps_hat[1:3],gamma_ps_hat[4]+DELTA,gamma_ps_hat[5:7]);
  gamma_ps_hat_diff[5,]=c(gamma_ps_hat[1:4],gamma_ps_hat[5]+DELTA,gamma_ps_hat[6:7]);
  gamma_ps_hat_diff[6,]=c(gamma_ps_hat[1:5],gamma_ps_hat[6]+DELTA,gamma_ps_hat[7]);
  gamma_ps_hat_diff[7,]=c(gamma_ps_hat[1:6],gamma_ps_hat[7]+DELTA)
  
  
  gamma_reg_hat_diff<-matrix(0,7,7);
  gamma_reg_hat_diff[1,]=c(gamma_reg_hat[1]+DELTA,gamma_reg_hat[2:7]);
  gamma_reg_hat_diff[2,]=c(gamma_reg_hat[1],gamma_reg_hat[2]+DELTA,gamma_reg_hat[3:7]);
  gamma_reg_hat_diff[3,]=c(gamma_reg_hat[1:2],gamma_reg_hat[3]+DELTA,gamma_reg_hat[4:7]);
  gamma_reg_hat_diff[4,]=c(gamma_reg_hat[1:3],gamma_reg_hat[4]+DELTA,gamma_reg_hat[5:7]);
  gamma_reg_hat_diff[5,]=c(gamma_reg_hat[1:4],gamma_reg_hat[5]+DELTA,gamma_reg_hat[6:7]);
  gamma_reg_hat_diff[6,]=c(gamma_reg_hat[1:5],gamma_reg_hat[6]+DELTA,gamma_reg_hat[7]);
  gamma_reg_hat_diff[7,]=c(gamma_reg_hat[1:6],gamma_reg_hat[7]+DELTA);
  
alpha20_diff<- matrix(0,length(gamma_ps_hat)+length(gamma_reg_hat),1);
  
for (i in 1:length(gamma_ps_hat)){
  pscore_aux<-int.cov%*%gamma_ps_hat_diff[i,];
  pscore_aux[pscore_aux>16] <- 16;
  alpha20_diff[i] <- (mean(B2/(1-plogis(pscore_aux)))-mean(B2/A2))/(gamma_ps_hat_diff[i,i]-gamma_ps_hat[i]);
}
  
for (i in (1+length(gamma_ps_hat)):(length(gamma_ps_hat)+length(gamma_reg_hat))){
  alpha20_diff[i] <- (mean(ps.fit*(1-d)*(deltaY-int.cov%*%gamma_reg_hat_diff[i-7,])/A2)-mean(B2/A2))/(gamma_reg_hat_diff[i-7,i-7]-gamma_reg_hat[i-7]);
}


  

  ########################################################
  ##########calculate sample mean and sd##################
  ########################################################
  sample_mean <- (mean(B1)-mean(B2/A2))/mean(B3);
  influ_sample <- matrix(c(B1 - matrix(colMeans(d*int.cov)%*%phi2_hat,n,1), B2/A2+t(phi_hat)%*%alpha20_diff, B3),n,3);
  jacobian_sample <- matrix(c(1/mean(B3),-1/mean(B3),-(mean(B1)-mean(B2/A2))/(mean(B3)^2)),1,3);
  phi_hat_sample <- influ_sample%*%t(jacobian_sample);
  SD_sample_mean <- (var(phi_hat_sample)/n)^0.5;
  cover_sample_mean <- abs(sample_mean - true) < qnorm(0.975) * SD_sample_mean;
  


  #######################################
  ##########seive estimation#############
  #######################################
  
 
  P0 = matrix(1,n,1);
  P1 = matrix(3^0.5*(2*A-1),n,1);
  P2 = matrix(5^0.5*(6*A^2-6*A+1),n,1);
  P3 = matrix(7^0.5*(20*A^3-30*A^2+12*A-1),n,1);
  P4 = matrix(9^0.5*(70*A^4-140*A^3+90*A^2-20*A+1),n,1);
  P5 = matrix(11^0.5*(252*A^5-630*A^4+560*A^3-210*A^2+30*A-1),n,1);
  P = cbind(P0,P1,P2,P3,P4,P5)[,1:(K_con+1)];
  beta2=solve(t(P)%*%P) %*% t(P)%*%B2;
  ################################################################
  ###calculate the (k)-th derivative of m evaluate at zero####
  ################################################################
  
  m = array(0,k);
  
  m[1] = t(dpda1[1:(K_con+1),]) %*% beta2;
  m[2] = t(dpda2[1:(K_con+1),]) %*% beta2;
  m[3] = t(dpda3[1:(K_con+1),]) %*% beta2;
  m[4] = t(dpda4[1:(K_con+1),]) %*% beta2;
  m[5] = t(dpda5[1:(K_con+1),]) %*% beta2;
  
  
  #####################################################
  ###bias estimator using m function defined before####
  #####################################################
  
  bias_estimator_1 = 0
  for(kappa in 1:k){
    bias_estimator_1 = bias_estimator_1 + mean(A2^(kappa-1)*(A2<h)) / factorial(kappa) * m[kappa];
  }
  
 
  
  ##########################################################
  ######mm report different m value for different gamma####
  ##########################################################
  
  AA<-array(0,n);
  beta22<-matrix(0,K_con+1,length(gamma_ps_hat)+length(gamma_reg_hat));
  mm<-matrix(0,K_con+1,length(gamma_ps_hat)+length(gamma_reg_hat));

  for (i in 1:length(gamma_ps_hat))
  {
    
    AA<- 1-plogis(int.cov%*%gamma_ps_hat_diff[i,]);
    P0 = matrix(1,n,1);
    P1 = matrix(3^0.5*(2*AA-1),n,1);
    P2 = matrix(5^0.5*(6*AA^2-6*AA+1),n,1);
    P3 = matrix(7^0.5*(20*AA^3-30*AA^2+12*AA-1),n,1);
    P4 = matrix(9^0.5*(70*AA^4-140*AA^3+90*AA^2-20*AA+1),n,1);
    P5 = matrix(11^0.5*(252*AA^5-630*AA^4+560*AA^3-210*AA^2+30*AA-1),n,1);
    PP = cbind(P0,P1,P2,P3,P4,P5)[,1:(K_con+1)];
    #beta1 = solve(t(P)%*%P + diag(array(sieve_regularization+1,K+1))) %*% t(P)%*%B1;
    beta22[,i] = solve(t(PP)%*%PP + diag(array(sieve_regularization,K_con+1))) %*% t(PP)%*%B2;
    
  }
  

  for (i in (length(gamma_ps_hat)+1):(length(gamma_ps_hat)+length(gamma_reg_hat))){
    BB <- ps.fit*(1-d)*(deltaY-int.cov%*%gamma_reg_hat_diff[i-7,]);
    beta22[,i] = solve(t(P)%*%P) %*% t(P)%*%BB;
  }
  
  for (i in 1:14){
    mm[1,i] = t(dpda1[1:(K_con+1),]) %*% beta22[,i];
    mm[2,i] = t(dpda2[1:(K_con+1),]) %*% beta22[,i];
    mm[3,i] = t(dpda3[1:(K_con+1),]) %*% beta22[,i];
    mm[4,i] = t(dpda4[1:(K_con+1),]) %*% beta22[,i];
    mm[5,i] = t(dpda5[1:(K_con+1),]) %*% beta22[,i];
  }
  
  ######################################################
  ### calculate alpha2_diff#############################
  #######################################################


  
  alpha2_diff<-matrix(0,14,1);
  for (i in 1:7) {
    alpha2_diff[i] = (alpha2(d,y1,y0,gamma_ps_hat_diff[i,],gamma_reg_hat,int.cov,k,mm,i,h)-alpha2(d,y1,y0,gamma_ps_hat,gamma_reg_hat,int.cov,k,mm,i,h))/DELTA
  }
  
  for (i in 8:14) {
    alpha2_diff[i] = (alpha2(d,y1,y0,gamma_ps_hat,gamma_reg_hat_diff[i-7,],int.cov,k,mm,i,h)-alpha2(d,y1,y0,gamma_ps_hat,gamma_reg_hat,int.cov,k,mm,i,h))/DELTA
  }
  #################################
  ###influence function omega2#####
  #################################
  
  
  ch_psi = 0;
  psi1 = P%*%solve(t(P)%*%P/n)%*%dpda1[1:(K_con+1),] * (B2-P%*%beta2);
  if(k>=1){
    ch_psi = ch_psi + matrix(mean(A2^0*(A2<h)/factorial(1))*psi1+A2^0*(A2<h)*m[1]/factorial(1),n,1)
  }
  psi2 = P%*%solve(t(P)%*%P/n)%*%dpda2[1:(K_con+1),] * (B2-P%*%beta2);
  if(k>=2){
    ch_psi = ch_psi + matrix(mean(A2^1*(A2<h)/factorial(2))*psi2+A2^1*(A2<h)*m[2]/factorial(2),n,1)
  }
  psi3 = P%*%solve(t(P)%*%P/n)%*%dpda3[1:(K_con+1),] * (B2-P%*%beta2);
  if(k>=3){
    ch_psi = ch_psi + matrix(mean(A2^2*(A2<h)/factorial(3))*psi3+A2^2*(A2<h)*m[3]/factorial(3),n,1)
  }
  psi4 = P%*%solve(t(P)%*%P/n)%*%dpda4[1:(K_con+1),] * (B2-P%*%beta2);
  if(k>=4){
    ch_psi = ch_psi + matrix(mean(A2^3*(A2<h)/factorial(4))*psi4+A2^3*(A2<h)*m[4]/factorial(4),n,1)
  }
  psi5 = P%*%solve(t(P)%*%P/n)%*%dpda5[1:(K_con+1),] * (B2-P%*%beta2);
  if(k>=5){
    ch_psi = ch_psi + matrix(mean(A2^4*(A2<h)/factorial(5))*psi5+A2^4*(A2<h)*m[5]/factorial(5),n,1)
  }
  ch_psi <- ch_psi+t(phi_hat)%*%alpha2_diff;
  #ch_psi <- ch_psi+t(rbind(phi1,phi2))%*%(alpha2_diff);
  
  
  influ_corrected <-matrix(c(B1-matrix(colMeans(d*int.cov)%*%phi2_hat,n,1) ,B2/A2*(A2>=h) + ch_psi,B3),n,3);
  jacobian_corrected <- matrix(c(1/mean(B3),-1/mean(B3),-(mean(B1)-mean(B2/A2*(A2>=h))-bias_estimator_1)/(mean(B3)^2)),1,3);
  phi_corrected_hat <-  influ_corrected %*%t(jacobian_corrected);
  SD_corrected_mean <-(var(phi_corrected_hat)/n)^0.5;
  sample_mean_corrected <- (mean(B1)-mean(B2/A2*(A2>=h))-bias_estimator_1)/mean(B3);
  cover_corrected_mean <- abs(sample_mean_corrected- true) < qnorm(0.975)*SD_corrected_mean;
  #self_normalized <- (sample_mean_corrected - true) / SD_corrected_mean;
  
  # #############################################
  # ################trimmed without BC###########
  # #############################################
  # # alpha2_diff_trimmed<- matrix(0,14,1);
  # # for (i in 1:7) {
  # #   alpha2_diff_trimmed[i] = (alpha2(d,y1,y0,gamma_ps_hat_diff[i,],gamma_reg_hat,int.cov,0,mm,i,h)-alpha2(d,y1,y0,gamma_ps_hat,gamma_reg_hat,int.cov,0,mm,i,h))/DELTA
  # # }
  # # 
  # # for (i in 8:14) {
  # #   alpha2_diff_trimmed[i] = (alpha2(d,y1,y0,gamma_ps_hat,gamma_reg_hat_diff[i-7,],int.cov,0,mm,i,h)-alpha2(d,y1,y0,gamma_ps_hat,gamma_reg_hat,int.cov,0,mm,i,h))/DELTA
  # # }
  # 
  # alpha2_diff_trimmed<- matrix(0,14,1);
  # for (i in 1:7) {
  #   alpha2_diff_trimmed[i] = (alpha_trimmed_without_BC(d,y1,y0,gamma_ps_hat_diff[i,],gamma_reg_hat,int.cov,h)-alpha_trimmed_without_BC(d,y1,y0,gamma_ps_hat,gamma_reg_hat,int.cov,h)/DELTA
  # }
  # 
  # for (i in 8:14) {
  #   alpha2_diff_trimmed[i] = (alpha2(d,y1,y0,gamma_ps_hat,gamma_reg_hat_diff[i-7,],int.cov,0,mm,i,h)-alpha2(d,y1,y0,gamma_ps_hat,gamma_reg_hat,int.cov,0,mm,i,h))/DELTA
  # }
  # 
  # influ_trimmed <-matrix(c(B1-matrix(colMeans(d*int.cov)%*%phi2_hat,n,1) ,B2/A2*(A2>=h) + t(phi_hat)%*%alpha2_diff_trimmed,B3),n,3);
  # jacobian_trimmed <- matrix(c(1/mean(B3),-1/mean(B3),-(mean(B1)-mean(B2/A2*(A2>=h)))/(mean(B3)^2)),1,3);
  # phi_trimmed <-  influ_trimmed %*%t(jacobian_trimmed);
  # SD_trimmed <-(var(phi_trimmed)/n)^0.5;
  # sample_mean_trimmed <- (mean(B1)-mean(B2/A2*(A2>=h)))/mean(B3);
  # cover_trimmed <- abs(sample_mean_trimmed- true) < qnorm(0.975)*SD_trimmed;
  # 
  # 
  
  ###############################
  ######display result###########
  ###############################
  results = c(true,sample_mean, cover_sample_mean, 2*qnorm(0.975)*SD_sample_mean, sample_mean_corrected, cover_corrected_mean,2*qnorm(0.975)*SD_corrected_mean);
  RESULTS = rbind(RESULTS, results);
  print(c(iter,colMeans(RESULTS)));
  
}



 # TRUE, MEAN, BIAS, sd,RMSE, COVER,95%length -- FOR SAMPLE MEAN
c(mean(RESULTS[,1]),mean(RESULTS[,2]),mean(RESULTS[,2])-mean(RESULTS[,1]),var(RESULTS[,2])^0.5,(mean(RESULTS[,2]-RESULTS[,1])^2+var(RESULTS[,2]))^0.5, mean(RESULTS[,3]),mean(RESULTS[,4]));


 #TRUE, MEAN, BIAS,sd,RMSE, COVER,95%length -- FOR BIAS-CORRECTED TRIMMED MEAN
c(mean(RESULTS[,1]),mean(RESULTS[,5]),mean(RESULTS[,5])-mean(RESULTS[,1]),var(RESULTS[,5])^0.5,(mean(RESULTS[,5]-RESULTS[,1])^2+var(RESULTS[,5]))^0.5,mean(RESULTS[,6]),mean(RESULTS[,7]));

  #TRUE, MEAN, BIAS,sd,RMSE, COVER,95%length -- FOR TRIMMED MEAN WITHOUT BIAS CORRECTION
# c(mean(RESULTS[,1]),mean(RESULTS[,8]),mean(RESULTS[,8])-mean(RESULTS[,1]),var(RESULTS[,8])^0.5,(mean(RESULTS[,8]-RESULTS[,1])^2+var(RESULTS[,8]))^0.5,mean(RESULTS[,9]),mean(RESULTS[,10]));



