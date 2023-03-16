
######## Yang and Ding (2018) biametrika######
dgps<-function(n,dgp,alpha){
# Gen covariates
  if (dgp==1){
nu = matrix(c(0,0,0),3,1);
mul_var =  matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1),3,3);
x<-MASS::mvrnorm(n,nu,mul_var);
x1 <- x[,1];
x2 <- x[,2];
x3 <- x[,3];
x4 <- runif(n,-3,3);
x5<-rchisq(n,1);
x6<- rbern(n,0.5);

pp<- plogis(alpha*(x1+x2+x3+x4+x5+x6));
#d  <- as.numeric(runif(n) <= pp)
d<- rbern(n,pp);
a<-1;
Y11 <- a*(x1+x2+x3-x4+x5+x6)+stats::rnorm(n, mean = 0, sd = 1);
Y10 <- stats::rnorm(n, mean = 0, sd = 1);
Y1 <-d*Y11+(1-d)*Y10;
Y0 <- a*(x1+x2+x3-x4+x5+x6)+stats::rnorm(n, mean = 0, sd = 1);

#Y0 <- stats::rnorm(n, mean = 0, sd = 1)
att_unf <-(mean(d*Y11)-mean(d*Y10))/mean(d);
return(list(y1 = Y1, y0 = Y0,y11=Y11,y10=Y10, d = d, pp=pp,  att= att_unf, X = cbind(x1,x2,x3,x4,x5,x6)))
  }
  
  if (dgp==2){
    nu = matrix(c(0,0,0),3,1);
    mul_var =  matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1),3,3);
    x<-MASS::mvrnorm(n,nu,mul_var);
    x1 <- x[,1];
    x2 <- x[,2];
    x3 <- x[,3];
    x4 <- runif(n,-3,3);
    x5<-rchisq(n,1);
    x6<- rbern(n,0.5);
    
    pp<- plogis(alpha*(x1+x2^2+x3^2+x4+x5+x6));
    #d<- rbern(n,pp);
    d  <- as.numeric(runif(n) <= pp)
    a<-1;
    Y11 <- a*(x1+x2+x3+x4+x5+x6)+stats::rnorm(n, mean = 0, sd = 1);
    Y10 <- stats::rnorm(n, mean = 0, sd = 1);
    Y1 <-d*Y11+(1-d)*Y10;
    Y0 <- stats::rnorm(n, mean = 0, sd = 1)
    
    #Y0 <- stats::rnorm(n, mean = 0, sd = 1)
    att_unf <-(mean(d*Y11)-mean(d*Y10))/mean(d);
    return(list(y1 = Y1, y0 = Y0, y11=Y11,y10=Y10,d = d, pp=pp, X = cbind(x1,x2,x3,x4,x5,x6)))
    
    
  }

  
  if (dgp==3){
    nu = matrix(c(0,0,0),3,1);
    mul_var =  matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1),3,3);
    x<-MASS::mvrnorm(n,nu,mul_var);
    x1 <- x[,1];
    x2 <- x[,2];
    x3 <- x[,3];
    x4 <- runif(n,-3,3);
    x5<-rchisq(n,1);
    x6<- rbern(n,0.5);
    
    pp<- plogis(alpha*(x1+x2+x3+x4+x5+x6));
    #d<- rbern(n,pp);
    d  <- as.numeric(runif(n) <= pp)
    a<-1;
    Y11 <- a*(x1+x2+x3)^2+stats::rnorm(n, mean = 0, sd = 1);
    Y10 <- stats::rnorm(n, mean = 0, sd = 1);
    Y1 <-d*Y11+(1-d)*Y10;
    Y0 <- stats::rnorm(n, mean = 0, sd = 1)
    #Y0 <- a*(x1+x2+x3)^2+stats::rnorm(n, mean = 0, sd = 1)
    
    #Y0 <- stats::rnorm(n, mean = 0, sd = 1)
    att_unf <-(mean(d*Y11)-mean(d*Y10))/mean(d);
    return(list(y1 = Y1, y0 = Y0,y11=Y11,y10=Y10, d = d, pp=pp, X = cbind(x1,x2,x3,x4,x5,x6)))
    
    
  }
  if (dgp==4){
    nu = matrix(c(0,0,0),3,1);
    mul_var =  matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1),3,3);
    x<-MASS::mvrnorm(n,nu,mul_var);
    x1 <- x[,1];
    x2 <- x[,2];
    x3 <- x[,3];
    x4 <- runif(n,-3,3);
    x5 <- rchisq(n,1);
    x6 <- rbern(n,0.5);
    
    pp <- plogis(alpha*(x1+x2^2+x3^2+x4+x5+x6));
    d <- rbern(n,pp);
    #d  <- as.numeric(runif(n) <= pp)
    a <- 1;
    Y11 <- a*(x1+x2+x3)^2+stats::rnorm(n, mean = 0, sd = 1);
    #Y11 <- a*(x1+x2-x3)^2+stats::rnorm(n, mean = 0, sd = 1);
    Y10 <- stats::rnorm(n, mean = 0, sd = 1);
    Y1 <- d*Y11+(1-d)*Y10;
    
    Y0 <- stats::rnorm(n, mean = 0, sd = 1)
    att_unf <-(mean(d*Y11) - mean(d*Y10))/mean(d);
    return(list(y1 = Y1, y0 = Y0,y11=Y11,y10=Y10, d = d, pp=pp, X = cbind(x1,x2,x3,x4,x5,x6)))
  }
}