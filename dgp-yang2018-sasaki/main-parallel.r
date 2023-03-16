library(foreach)
library(doParallel)

start_time <- proc.time()

#############
# Constants #
#############
main_dir = ""
dgp_path = "~/Desktop/simulation-20230311/simulation/dgp/dgp_y0_non_rdm.R"   # enter your own path here     
dimX = 6              # Dimension of X
dimZ = 6              # Dimension of Z
DELTA = 0.00000000001 # Precision of numerical derivative

##############
# Parameters #
##############
n <- 500                                      # Sample size
MC_ITERATIONS <- 500                         # Monte Carlo
alpha <- 0.8                                    #parameter in selection model
klist <- c(3)                             # list of k, K
hlist <- seq(from = 0, to = 0.1, by = 0.01)  # list of h, bandwidth
b = 0.01  # bandwidth for influence function estimation

dgp_list  <- list(1, 2)                # list of dgp
#true_list <- list(1.732, 1.50, 2.874, 3.215)        # list of dgp true values alpha=0.7
true_list <- list(1.73, 2.064, 2.905, 3.243)        # list of dgp true values alpha=0.8

selected_dgps <- c(1,2)                     # selected dgps

#######################
# LEGENDRE POLYNOMIAL #
#######################
shifted_legendre <- function(k, x){
  lp = 0
  for( kdx in 0:k ){
    lp = lp + choose(k,kdx)*choose(k+kdx,kdx)*(((2*x-1)-1)/2)^(kdx)
  }
  return( lp )
}

############################################################
# DERIVATIVE OF LEGENDRE POLYNOMIAL BY RECURSIVE ALGORITHM #
############################################################
deriv_shifted_legendre <- function(k,x,order){
  if(order == 0){
    return( shifted_legendre(k,x) )
  }else{
    ( deriv_shifted_legendre(k,x+DELTA,order-1) - deriv_shifted_legendre(k,x,order-1) ) / (DELTA)
  }
}

################################################
# ALPHA2 FUNCTION FOR JACOBIAN FOR FIRST STAGE #
################################################
alpha2 <- function(gamma,Z,D,Y0,Y1){
  gamma_ps = matrix(gamma[1:(dimZ+1)],,1)
  gamma_reg = matrix(gamma[(dimZ+2):(2*dimZ+2)],,1)
  reg = Z%*%gamma_reg
  ps  = exp(Z%*%gamma_ps)/(1+exp(Z%*%gamma_ps))
  B2_ = ps*(1-D)*((Y1-Y0) - reg)
  A2_ = 1 - ps
  pK_ = shifted_legendre(0,A2_)
  for( kdx in 1:K ){
    pK_ = cbind(pK_, shifted_legendre(kdx,A2_))
  }
  sieve_ = ginv(t(pK_)%*%pK_)%*%t(pK_)%*%B2_
  
  plist = 0:100/100
  integral = 0
  for( p in plist ){
    if(p<1-h){
      P_ = exp(Z%*%gamma_ps)/(1+exp(Z%*%gamma_ps))
      NU_ = X%*%gamma_reg
      integral = integral + mean((1-D)*((Y1-Y0)-NU_)*dnorm((P_-p)/b)/b)*p/(1-p)
    }else{
      if( k > 0 ){
        for( kappa in 1:k ){
          pKkappa_ = deriv_shifted_legendre(0,0,kappa)
          for( kdx in 1:K ){
            pKkappa_ = cbind(pKkappa_, deriv_shifted_legendre(kdx,0,kappa))
          }
          m_hat_kappa_ = pKkappa_%*%sieve_  
          integral = integral + (1-p)^(kappa-1)*mean(dnorm((P_-p)/b)/b)/factorial(kappa)*m_hat_kappa_
        }
      }
    }
  }
  return( integral * (plist[2]-plist[1]) )
}

########################################
# Start to Generate Simulation Results #
########################################

for(dgp_index in selected_dgps){
  
  dgp  <- dgp_list[[dgp_index]]
  true <- true_list[[dgp_index]]
  
  for(k in klist) {
    K = k                   # Sieve order
 
    # Prepare output files. Create path if not exists.
    # filename is like main_dir/results/dgp1/simu_k=3_MC=2000.csv
    sub_dir <- paste0("results/", "dgp", dgp)
    output_dir = paste0(main_dir, sub_dir)
    filename <- paste0(output_dir, "/simu_n=", n, "_k=", k, "_MC=", MC_ITERATIONS, "_alpha=", alpha, ".csv")
    ifelse(!dir.exists(file.path(output_dir)), dir.create(file.path(output_dir)), FALSE)
    file_conn <- file(filename, open = "w")
    header <- c("h", "BIAS", "SDEV", "SDE", "RMSE_old", "RMSE_new", "CF95", "CF95length") # write output file headers.
    cat(paste(header, collapse = ","), file = file_conn, "\n")
    
    for(h in hlist) { # h is bandwidth
      
      #####################
      # BEGIN MONTE CARLO #
      #####################
      number_of_cores = 8 # or detectCores().   
      # Set up parallel processing environment
      cl <- makeCluster(number_of_cores)
      registerDoParallel(cl)
      
      clusterEvalQ(cl, {
        library(MASS)
        library(Rlab)
      })
      
      RESULTS = NULL
      RESULTS <- foreach(iter = 1 : MC_ITERATIONS, .combine = rbind) %dopar% {
        ###################
        # DATA GENERATION #
        ###################
        source(dgp_path)
        data.did<-dgps(n, dgp, alpha)
        Y0 <- data.did$y0
        Y1 <- data.did$y1
        D  <- data.did$d
        X  <- data.did$X
        X  <- cbind(1, X)
        
        ####################################
        # ESTIMATION OBSERVING (Y0,Y1,D,Z) #
        ####################################
        reg_result = lm( I(Y1-Y0) ~ X[,2:7], subset=(D==0) )
        ps_result  = glm( D ~ X[,2:7], family=binomial(link='logit') )
        gamma_reg_hat = matrix(coef(reg_result),,1)
        gamma_ps_hat  = matrix(coef(ps_result),,1)
        reg_hat = X%*%gamma_reg_hat
        ps_hat  = exp(X%*%gamma_ps_hat)/(1+exp(X%*%gamma_ps_hat))
        
        B1 = D*((Y1-Y0) - reg_hat)
        B2 = ps_hat*(1-D)*((Y1-Y0) - reg_hat)
        A2 = 1 - ps_hat
        B3 = D
        
        pK = shifted_legendre(0,A2)
        for( kdx in 1:K ){
          pK = cbind(pK, shifted_legendre(kdx,A2))
        }
        sieve = ginv(t(pK)%*%pK)%*%t(pK)%*%B2
        
        BC = 0
        if( k > 0 ){
          for( kappa in 1:k ){
            pKkappa = deriv_shifted_legendre(0,0,kappa)
            for( kdx in 1:K ){
              pKkappa = cbind(pKkappa, deriv_shifted_legendre(kdx,0,kappa))
            }
            m_hat_kappa = pKkappa%*%sieve 
            BC = BC + mean(A2^(kappa-1) * (A2 < h))/factorial(kappa) * m_hat_kappa
          }
        }
        alpha2_hat = mean(B2/A2*(A2>=h)) + BC
        
        EST = (mean(B1) - alpha2_hat) / mean(B3)
        
        ##################################
        # FIRST STAGE INFLUENCE FUNCTION #
        ##################################
        phi1_den=0
        phi1_num=NULL
        phi2_den=0
        phi2_num=NULL
        for( idx in 1:n ){
          phi1_den = phi1_den + matrix(X[idx,],,1)%*%(ps_hat[idx]*(1-ps_hat[idx]))%*%matrix(X[idx,],1,)
          phi1_num = cbind(phi1_num, matrix(X[idx,],,1)%*%(D[idx]-ps_hat[idx]))
          phi2_den = phi2_den + matrix(X[idx,],,1)%*%(1-D[idx])%*%matrix(X[idx,],1,)
          phi2_num = cbind(phi2_num, matrix(X[idx,],,1)%*%((1-D[idx])*(Y1[idx]-Y0[idx])))
        }
        phi1_hat = ginv(phi1_den) %*% phi1_num
        phi2_hat = ginv(phi2_den) %*% phi2_num
        phi_hat = rbind(phi1_hat,phi2_hat)
        
        ############################
        # JACOBIAN FOR FIRST STAGE #
        ############################
        gamma_hat = c(gamma_ps_hat,gamma_reg_hat)
        deriv_alpha2_hat = NULL
        for( idx in 1:length(gamma_hat) ){
          gamma_hat_d = gamma_hat
          gamma_hat_d[idx] = gamma_hat_d[idx] + DELTA
          deriv_alpha2_hat = c(deriv_alpha2_hat, alpha2(gamma_hat_d,X,D,Y0,Y1))
        }
        deriv_alpha2_hat = ( deriv_alpha2_hat - rep(alpha2(gamma_hat,X,D,Y0,Y1),length(gamma_hat)) ) / DELTA
        
        ##############
        # OMEGA2 HAT #
        ##############
        omega2_hat1 = B2/A2*(A2>=h)
        omega2_hat2 = 0
        omega2_hat3 = 0
        if( k > 0 ){
          for( kappa in 1:k ){
            pKkappa = deriv_shifted_legendre(0,0,kappa)
            for( kdx in 1:K ){
              pKkappa = cbind(pKkappa, deriv_shifted_legendre(kdx,0,kappa))
            }
            m_hat_kappa = pKkappa%*%sieve 
            omega2_hat2 = omega2_hat2 + (A2^(kappa-1) * (A2 < h))/factorial(kappa) * c(m_hat_kappa)
            varphi2_hat_kappa = t( pKkappa%*%ginv(t(pK)%*%pK/n)%*%t(pK * matrix(B2-pK%*%sieve,n,K+1)) )
            omega2_hat3 = omega2_hat3 + mean(A2^(kappa-1)*(A2 < h))/factorial(kappa) * varphi2_hat_kappa
          }
        }
        omega2_hat4 = t( deriv_alpha2_hat%*%phi_hat )
        omega2_hat = omega2_hat1 + omega2_hat2 + omega2_hat3 + omega2_hat4
        
        ##############
        # VARPHI HAT #
        ##############
        varphi_hat = ( D*(Y1-Y0-reg_hat) - t(colMeans(matrix(D,n,dimZ+1)*X)%*%phi2_hat) ) / mean(D) - omega2_hat / mean(D) - c(mean( D*(Y1-Y0-reg_hat) ) - alpha2_hat) / mean(D)^2 * D
        
        ##################
        # STANDARD ERROR #
        ##################
        SE = ( var(varphi_hat)/n )^0.5
        
        ############
        # COVERAGE #
        ############
        COV90 = abs(EST-true) < SE * qnorm(0.95)
        COV95 = abs(EST-true) < SE * qnorm(0.975)
        
        if(!is.nan(EST)){
          result = cbind(EST,COV90,COV95,SE)
        }

        return(result)
      }
      
      # Stop the parallel processing environment
      stopCluster(cl)
      
      #################
      # WRITE RESULTS #
      #################
      BIAS = mean(RESULTS[,1]-true)
      SDEV  = var(RESULTS[,1])^0.5
      SDE = mean(RESULTS[,4])
      RMSE_old =(BIAS^2+var(RESULTS[,1]))^0.5
      RMSE_new = ( BIAS^2 + SDE^2 )^0.5
      CF95 = mean(RESULTS[,3])
      CF95length = 2 *qnorm(0.975)*mean(RESULTS[,4])
      result_row = c(h, BIAS, SDEV, SDE, RMSE_old, RMSE_new, CF95, CF95length) # Don't forget to modify headers too
      cat(paste(result_row, collapse = ","), "\n", file = file_conn)
      print(result_row)
    }
    close(file_conn) # close output file for given dgp, k, h, and MC
  }
}

# Record end time
end_time <- proc.time()

# Calculate elapsed time
elapsed_time <- end_time - start_time

# Print elapsed time
print(elapsed_time)