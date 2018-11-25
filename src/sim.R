# This code is converted from 'code13May.m'
# 
#
# Date          Changed
#----------------------------------------------------------------------------
# 2016-12-18     Original Code
# 2017-01-22     Change to mimic the MATLAB code 

# -----------------------------------------------------------
# 
#
# PURPOSE:
#     Conduct simulation studies in Lee and Shin
#
# INPUT:
#     ap.hat: (R x d) matrix; the estimates in each replication
#             where R is the number of replication
#                   d is the size of estimator
#     ap.0: (d x 1) vector; the true value of ap
#
# OUPUT:
#     mb: (d x 1) vector; median of bias distribution
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#   2016-12-18    Y. Shin             Original Code
#   2017-02-13    Y. shin             Change the structure  
#   2017-03-03    Y. Shin             Make it suitable for Sharcnet. No KO
#   2018-01-17
#   2018-06-07       
#   2018-09-20

rm(list=ls())
start.time = proc.time()[3]
library(MASS)
library(sandwich)
library(doMC)

seed='001'
set.seed(as.numeric(seed))

source('lib_csa.R')

# Set parameter values
N=100               # Sample Size
K=20                # The number of instruments
bt.0=0.1            # The true parameter value for bt
rho.eu=0.1          # Cov(e,u), the degree of endogeneity
R.sq.f=0.01         # Theoretical first stage R-squared
rho.z=0.5           # Correlation between instruments
rep=1000            # The number of replications
high.percentile=90  # High value for inter-percentile deviation
low.percentile=10   # Low value for inter-percentile deviation
pi.con=rep(sqrt(R.sq.f/(K*(1-R.sq.f))),K) # pi.con: constant pi
pi.dec=gen_pi_dec(K,rho.z,R.sq.f)         # pi.dec: decreasing pi
pi.irr=gen_pi_irr(K,rho.z,R.sq.f)
pi.0=pi.con              # The true parameter values for pi.0
R=1000
id = c(1:N)              # no-clustering in this simulation
I.N = diag(N)            # N-dim identity matrix

# parallel computing
use.par = TRUE   # TURE or FALSE

# Declare variables for simulation results
bt.OLS=matrix(NA,rep,2)     # OLS 
se.OLS=matrix(NA,rep,2)     # standard error
colnames(bt.OLS)=c('x.endogenous', 'const')

bt.TSLS=matrix(NA,rep,2)    # 2SLS - all IV's
se.TSLS=matrix(NA,rep,2)
#colnames(bt.TSLS)=c('const','x.endogenous')

bt.DN.M=matrix(NA,rep,2)    # 2SLS - Donald and Newey, 'Mallows' prelim est
se.DN.M=matrix(NA,rep,2)

bt.KO.M=matrix(NA,rep,2)    # 2SLS - Kuersteiner and Okui, 'Mallows' prelim est
se.KO.M=matrix(NA,rep,2)

bt.CSA.M=matrix(NA,rep,2)   # 2SLS - Complete Subset Average, 'Mallows' prelim est
se.CSA.M=matrix(NA,rep,2)

# Declare variables for checking mid-results
k.DN.M    = rep(NA, rep)      # Estimated number of IVs in Donald and Newey method, 'Mallows' prelim est
k.DN.inf  = rep(NA, rep)    # Estimated number of IVs in Donald and Newey method, infeasible true model

w.KO.M   = matrix(NA, nrow=rep, ncol=K)
w.KO.inf = matrix(NA, nrow=rep, ncol=K)

k.CSA.M   = rep(NA, rep)      
k.CSA.inf = rep(NA, rep)    

# Repeat the simulation for 'rep' times 
for (re in (1:rep)){
  
  cat('------------- \n')
  cat(re,'iterations \n')
  cat('------------- \n')

  
  # Data generation
  dat = dgp_LS( N=N, K=K, rho.eu=rho.eu, rho.z=rho.z, pi.0=pi.0, bt.0=bt.0 )
  y = dat$y
  z = dat$z
  x.end = as.matrix(dat$x.end)
  colnames(x.end)='x.end'
  x.exo = rep(1,N)
  x = cbind(x.end, x.exo)
  f = dat$f
  e = dat$e
  u = dat$u
  Sigma.z = dat$Sigma.z
  

  #--------------------------------------
  # Calculate aveP matrix 
  #--------------------------------------
  # Use multi cores / single core (no parallel)
  if (use.par==T){
    n.core = detectCores()-1  
    aveP=tsls_CSA_get_P(y, x.end, x.exo, z.excl=z, R=R, ld=NULL, sub.K=c(1:K.excl), use.par = T, n.core=n.core)
  } else {
    aveP=tsls_CSA_get_P(y, x.end, x.exo, z.excl=z, R=R, ld=NULL, sub.K=c(1:K.excl))
  }
  
  #--------------------------------------
  # OLS
  #--------------------------------------
  m.OLS = lm(y~x.end)
  bt.OLS[re,]=c(coef(m.OLS)['x.end'], coef(m.OLS)[1])
  se.OLS[re,]=cluster.se(id,bt.OLS[re,],x,y,diag(N))$se

  #--------------------------------------
  # TSLS: using all IV
  #--------------------------------------
  bt.TSLS[re,]=tsls_est( y=y, x.end=x.end, x.exo=x.exo, z=z )
  se.TSLS[re,]=cluster.se(id,bt.TSLS[re,],x,y,P=aveP[[K]])$se
  
  
  
  #--------------------------------------
  # 2SLS - Donald and Newey
  #--------------------------------------
  
  # Pick a pre.iv that has the highest correlation with x.end
  z.sort=z
  cor.order = rank(-abs(cor(x.end,z)))
  z.sort[,cor.order] = z
  z.pre.one=z.sort[,1]
  
  # Two-step: Mallows
  m.DN.M=tsls_DN(y=y, x.end=x.end, x.exo=x.exo, z.excl=z, ld=NULL, method.pre='two.step', z.pre=z.pre.one)
  bt.DN.M[re,]=m.DN.M$bt.hat
  se.DN.M[re,]=cluster.se(id,bt.DN.M[re,],x,y,P=m.DN.M$P.DN)$se
  k.DN.M[re]=m.DN.M$opt.k
  

  #--------------------------------------
  # 2SLS - Kruesteiner and Okui
  #--------------------------------------

  # Two-step: Mallows
  m.KO.M = tsls_KO( y=y, x.end=x.end, x.exo=x.exo, z.excl=z, method.pre='two.step', z.pre=z.pre.one)
  bt.KO.M[re,] = m.KO.M$bt.hat
  se.KO.M[re,] = cluster.se(id,bt.KO.M[re,],x,y,P=m.KO.M$P.KO)$se
  w.KO.M[re,] = m.KO.M$opt.w 
  

  #--------------------------------------
  # 2SLS - Complete Subset Average
  #--------------------------------------

  # Two-step
  m.CSA.M = tsls_CSA_aveP_in( y=y, x.end=x.end, x.exo=x.exo, z.excl=z, ld=NULL, aveP=aveP, method.pre='two.step', z.pre=z.pre.one)
  bt.CSA.M[re,] = m.CSA.M$bt.hat
  k.CSA.M[re] = m.CSA.M$opt.k
  se.CSA.M[re,] = cluster.se(id,bt.CSA.M[re,],x,y,P=aveP[[k.CSA.M[re]]])$se
  
}


# Print computation time & save the result    
runt <- proc.time()[3]-start.time
runt_h <- floor(runt/3600)
runt_m <- floor( (runt/3600 - runt_h) * 60 )
runt_s <- floor( (runt/60 - (runt_h*60) - runt_m) * 60 )
cat('\n',"runtime = ",runt_h,"h",runt_m,"m",runt_s,"s",'\n','\n')
cat('-----------------------------------------------------------','\n')
cat('runtime = ',runt,'\n')

library('xtable')

name.set = c('OLS','TSLS','DN.M','KO.M','CSA.M')

  bt.hat = cbind(bt.OLS[,1], bt.TSLS[,1], bt.DN.M[,1],bt.KO.M[,1],bt.CSA.M[,1])
  se = cbind(se.OLS[,1], se.TSLS[,1],se.DN.M[,1],se.KO.M[,1],se.CSA.M[,1])
  k.hat = cbind(NA, NA, k.DN.M, NA,  k.CSA.M )
  colnames(bt.hat)=name.set
  #median.bt = median(bt.hat)  # This code should be fixed 
  bias = apply(bt.hat,2,bias,ap.0=bt.0) 
  mb = apply(bt.hat,2,mb,ap.0=bt.0)
  mse = apply(bt.hat,2,mse,ap.0=bt.0)
  mad = apply(bt.hat,2,mad)
  range = apply(bt.hat,2,range,high=90, low=10)
  cover = sapply(c(1:5), function(i) {coverage(bt.hat[,i], se[,i], ap.0=bt.0)} )
  mean.k = apply(k.hat,2,mean)
  med.k = apply(k.hat,2, median)
  
  result=cbind(mse, bias, mad, mb, range, cover, mean.k, med.k)
  rownames(result)=name.set
  colnames(result)=c('MSE', 'Bias', 'MAD', 'Median Bias', 'Range', 'Coverage', 'Mean(K.hat)', 'Med(K.hat)')
  sink(file='../result/sim-table-1-panel-1.out')
  print(round(result,3))
  sink()

