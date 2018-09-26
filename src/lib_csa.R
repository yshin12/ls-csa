# Updates
#
# 2018-03-21: Implement the sampling procedure for the complete subsets (tsls_LS)
# 2017-12-14: Tidy-up 'get_weight'; changes in TSLS-LS Objective function form updated
#


gen_pi_con = function(K,rho.z,R.sq.f) {
  pi.con = rep(sqrt( R.sq.f / ( ( K + K * ( K - 1 ) * rho.z ) * ( 1 - R.sq.f ) ) ), K)
  return(pi.con)
}

gen_pi_dec = function(K,rho.z,R.sq.f) {
  
  sum1 = sum2 = 0
  for (ind.k in (1:K)) {
    sum1 = sum1 + ( 1 - ind.k / (K+1) )^8 
    Set.for.j = c(1:K)
    Set.for.j = Set.for.j[-ind.k]
    for (ind.j in Set.for.j) {
      sum2 = sum2 + ( 1 - ind.k / ( K + 1 ) )^4 * ( 1 - ind.j / (K+1) )^4 * rho.z	
    }
  }
  
  c.K.dec = sqrt( ( R.sq.f / ( 1-R.sq.f ) ) * ( 1 / (sum1 + sum2) ) )
  
  pi.dec = c.K.dec * ( 1 - c(1:K) / ( K + 1 ) )^4
  
  return(pi.dec)
  
}


gen_pi_irr = function(K,rho.z,R.sq.f) {
  
  half.K = K / 2
  sum1 = sum2 = 0
  for (ind.k in ((half.K+1):K)) {
    sum1 = sum1 + ( 1 - (ind.k-half.K) / (half.K+1) )^8 
    Set.for.j = c((half.K+1):K)
    Set.for.j = setdiff(Set.for.j,ind.k)
    for (ind.j in Set.for.j) {
      sum2 = sum2 + ( 1 - (ind.k - half.K) / ( half.K + 1 ) )^4 * ( 1 - (ind.j - half.K) / (half.K+1) )^4 * rho.z	
    }
  }
  
  c.K.irr = sqrt( ( R.sq.f / ( 1-R.sq.f ) ) * ( 1 / (sum1 + sum2) ) )
  
  K.vec = c(1:K)
  
  pi.irr = c.K.irr * ( 1 - (K.vec - half.K ) / ( half.K + 1 ) )^4 * ( K.vec > half.K ) 
  
  return(pi.irr)
  
}


# dgp_LS
# PURPOSE:
#     Generate data for simulation in Lee and Shin 
#
# USAGE: 
#     dgp_LS(N,K,rho.eu,rho.z,pi.0,bt.0)
#
# MODEL:
# 
#     Structural Eq:
#           y = x * beta + e
#     1st Stage Eq:
#           x = z' * pi  + u
#
# INPUT:
#     N: a scalar; the sample sizee
#     K: a scalar; the number of instruments
#     rho: a scalar; Cov(e, u)
#                            | 1        rho.eu |
#     Cov(e,u) = Sigma.eu =  | rho.eu   1      |
#      2 x 2
#
#     rho.z: a scalar; parameter for the correlation of z
#                          |   1       rho.z  ...   rho.z  |
#     Var(Z) = Sigma.Z  =  |   rho.z     1    ...   rho.z  |
#    ( K x K )             |   ...     ...    ...   ...    |
#                          |   rho.z   rho.z  rho.z  1     |
#     pi.0: K x 1 vector; pi 
#     bt.0: a scalar; beta
#     high: a scalar in [0,100] ; higher percentile
#     low: a scalar in [0,100] ; lower percentile
#
# OUPUT:
#     y: N x 1; a dependenvt variable
#     x.end: N x 1; an endogenous variable
#     z: N x K; instruments
#     
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#    2/11/17      Y. shin             Original Code       
# -----------------------------------------------------------
dgp_LS = function(N, K, rho.eu, rho.z, pi.0, bt.0){
  # Construct Cov(e.u) using rho.eu
  Sigma.eu = matrix(c(1,rho.eu,rho.eu,1),2,2)
  
  # Generate error terms e and u
  chol.eu = chol(Sigma.eu)
  eu = matrix(rnorm(N*2),N,2) %*% chol.eu
  e = eu[,1]   
  u = eu[,2]
  
  # Construct Var(Z)
  Sigma.z = matrix(rho.z,K,K)
  diag(Sigma.z) = 1
  
  # Generate instruments Z, N x K
  chol.z = chol(Sigma.z)
  z = matrix(rnorm(N*K),N,K) %*% chol.z
  
  # Construct a vector of endogenous variables x
  x = z %*% pi.0 + u
  
  # Construct a vector of dependent variable y
  y = x * bt.0 + e
  
  # Return a generated sample
  return( list( y=y, x.end=x, z=z ,e=e, u=u, Sigma.z=Sigma.z) )
}


# -----------------------------------------------------------
# USAGE: 
#     proj_m(z)
#
# PURPOSE:
#     Calculate the projection matrix of z
#
# INPUT:
#     z: N x K, a matrix 
#
# OUPUT:
#     proj.z: N x n; the projection matrix of z
#            proj.z = z(z'z)^{-1}z'
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#    2/11/17      Y. shin             Original Code       
# -----------------------------------------------------------
proj_m = function(z){
  #library('gpuR')
  #z = gpuMatrix(z)
  #inv.z = (z %*% as.matrix( solve( t(z) %*% z ) ) %*% t(z))
  #return(as.matrix(inv.z))
  (z %*% ( solve( t(z) %*% z ) ) %*% t(z))
}

tr = function(X) {
  sum(diag(X))
}

matmul = function(X) {
  X %*% X
}

# -----------------------------------------------------------
# USAGE: 
#     ls(y,x)
#
# PURPOSE:
#     Calculate the least squares estimator bt.hat for
#         y = x' bt + eps
#
# INPUT:
#     y: N x 1; a dependent variable
#     x: N x K; regressors
#
# OUPUT:
#     bt.hat: K x 1, the estimate for bt
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#    2/11/17      Y. shin             Original Code       
# -----------------------------------------------------------
ols = function(y,x){
  bt.hat = solve( t(x) %*% x ) %*% t(x) %*% y
  return(bt.hat)
}

# -----------------------------------------------------------
# USAGE: 
#     tsls_est(y,x.end,x.exo,z)
#
# PURPOSE:
#     Calculate the two-stage least squares estimator bt.hat for
#         y = x' bt + eps
#         where x = ( x.end, x.exo )'
#
# INPUT:
#     y: N x 1, a dependent variable
#     x.end: N x d1, endogenous regressors
#     x.exo: N x d2, exogenous regressors
#     z: N x K; instrudments excluding x.exo
#
# OUPUT:
#     bt.hat: (d1+d2) x 1, the estimate for bt
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#    2/11/17      Y. shin             Original Code       
# -----------------------------------------------------------
tsls_est = function(y, x.end, x.exo, z) {

  # Combine z and x.exo if x.exo exists
  if ( !is.null(x.exo) ) {
    z.big = cbind(z, x.exo)
  } else {
    z.big = z
  }
  x = cbind(x.end, x.exo)
  
  # Construct a projection matrix with z.big
  p=proj_m(z.big)
  
  # Calculate 2SLS estimator for b.hat
  bt.hat=solve(t(x)%*%p%*%x) %*% t(x)%*%p%*%y
  return(bt.hat)
}


# -----------------------------------------------------------
# USAGE: 
#     mb(ap.hat, ap.0)
#
# PURPOSE:
#     Calculate the mean squared error of the estimator
#
# INPUT:
#     ap.hat: (R x d) matrix; the estimates in each replication
#             where R is the number of replication
#                   d is the size of estimator
#     ap.0: (d x 1) vector; the true value of ap
#
# OUPUT:
#     mse: (d x 1) vector; mse of each parameter values
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#    2/11/17      Y. shin             Original Code       
# -----------------------------------------------------------
mse = function ( ap.hat, ap.0 ) {
  
  # Change it into a matrix just in case it is a numeric vector
  ap.hat = as.matrix(ap.hat)
  
  # Define constants
  R = nrow(ap.hat)   # The number of replications
  d = ncol(ap.hat)   # The number of parameters
  
  # Construct a (Rxd) matrix of true parameter values 
  ones = rep(1,R)               # (Rx1) vector of 1's
  mat.ap.0 = t(ap.0) %x% ones   # (Rxd) matrix of true values 
  
  # Calculate the squared error 
  sq.err = ( ap.hat - ap.0 ) ^ 2
  
  # Calculat the mean squared error 
  mse = apply (sq.err, 2, mean)
  return(mse)
}


# -----------------------------------------------------------
# USAGE: 
#     mad(ap.hat, ap.0)
#
# PURPOSE:
#     Calculate median abolute deviation between the estimator
#     and its true value
#
# INPUT:
#     ap.hat: (R x d) matrix; the estimates in each replication
#             where R is the number of replication
#                   d is the size of estimator
#     ap.0: (d x 1) vector; the true value of ap
#
# OUPUT:
#     mad: (d x 1) vector; median of asolute deviation distribution
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#    2/11/17      Y. shin             Original Code       
# -----------------------------------------------------------
mad = function ( ap.hat, ap.0 ) {
  
  # Change it into a matrix just in case it is a numeric vector
  ap.hat = as.matrix(ap.hat)
  
  # Error if the dimension does not match
  try( if ( ncol(ap.hat) != length(ap.0) ) stop('ERROR : Dimesion does not match!'))
  
  # Define constants
  R = nrow(ap.hat)   # The number of replications
  d = ncol(ap.hat)   # The number of parameters
  
  # Construct a (Rxd) matrix of true parameter values 
  ones = rep(1,R)               # (Rx1) vector of 1's
  mat.ap.0 = t(ap.0) %x% ones   # (Rxd) matrix of true values 
  
  # Calculate the absolute values of the bias
  abs.bias=abs(ap.hat-mat.ap.0)
  
  # Calculate te meidan of the abolute biases
  median.absolute.dev=apply(abs.bias,2,median)
  return(as.numeric(median.absolute.dev))
}


# -----------------------------------------------------------
# USAGE: 
#     range(ap.hat, high, low)
#
# PURPOSE:
#     Calculate inter-percentile range of the estimator
#
# INPUT:
#     ap.hat: (R x d) matrix; the estimates in each replication
#             where R is the number of replication
#                   d is the size of estimator
#     high: a scalar in [0,100] ; higher percentile
#     low: a scalar in [0,100] ; lower percentile
#
# OUPUT:
#     inter.per: (d x 1) vector; (high-low) percentile range
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#    2/11/17      Y. shin             Original Code       
# -----------------------------------------------------------

range = function (ap.hat, high, low) { 
  
  # Change it into a matrix just in case it is a numeric vector
  ap.hat = as.matrix(ap.hat)
  
  
  # Calculate the percentiles at high and low
  high.per = apply(ap.hat, 2, quantile, probs=high/100)
  low.per = apply(ap.hat, 2, quantile, probs=low/100)
  
  # Calculate and return the inter percentile
  inter.per = high.per - low.per
  return(inter.per)
}


# -----------------------------------------------------------
# USAGE: 
#     mb(ap.hat, ap.0)
#
# PURPOSE:
#     Calculate the median of the bias distribution
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
#    2/11/17      Y. shin             Original Code       
# -----------------------------------------------------------
mb = function ( ap.hat, ap.0 ) {
  
  # Change it into a matrix just in case it is a numeric vector
  ap.hat = as.matrix(ap.hat)
  
  # Define constants
  R = nrow(ap.hat)   # The number of replications
  d = ncol(ap.hat)   # The number of parameters
  
  # Construct a (Rxd) matrix of true parameter values 
  ones = rep(1,R)               # (Rx1) vector of 1's
  mat.ap.0 = t(ap.0) %x% ones   # (Rxd) matrix of true values 
  
  # Calculat the bias
  bias=ap.hat-ap.0
  
  # Calculate the median bias and return the value
  mb=apply(bias, 2, median) 
  return(mb)
}

bias = function ( ap.hat, ap.0 ) {
  
  # Change it into a matrix just in case it is a numeric vector
  ap.hat = as.matrix(ap.hat)
  
  # Define constants
  R = nrow(ap.hat)   # The number of replications
  d = ncol(ap.hat)   # The number of parameters
  
  # Construct a (Rxd) matrix of true parameter values 
  ones = rep(1,R)               # (Rx1) vector of 1's
  mat.ap.0 = t(ap.0) %x% ones   # (Rxd) matrix of true values 
  
  # Calculat the bias
  bias=ap.hat-ap.0
  
  # Calculate the median bias and return the value
  bias=apply(bias, 2, mean) 
  return(bias)
}



gen_p_triangle = function(K) {
  p.triangle = matrix(0,K,K)
  for (i in (1:K)) {
    for (j in (1:i)){
      p.triangle[i,j]=choose(i,j)
    }
  }
  return(p.triangle)
}


# -----------------------------------------------------------
# USAGE: 
#     tsls_CSA(y, x.end, x.exo, z, R)
#
# PURPOSE:
#     Calculate the two-stage least squares estimator using subsets for
#         y = x' bt + eps
#         where x = ( x.end, x.exo )'
#
# INPUT:
#     y: N x 1, a dependent variable
#     x.end: N x d1, endogenous regressors
#     x.exo: N x d2, exogenous regressors
#     z: N x K; instrudments excluding e.exo
#
# OUPUT:
#     bt.hat: (d1+d2) x 1, the estimate for bt
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#    2/11/17      Y. shin             Original Code       
#   [Should be revised later to reflect the maximum number of subsets. ]
#
# -----------------------------------------------------------

tsls_CSA = function(y, x.end, x.exo, z.excl, method.pre, z.pre, R){
  
  # Declare constants
  z = as.matrix(cbind(x.exo,z.excl))
  x = cbind(x.end, x.exo)
  iv.pre = cbind(x.exo, z.pre)
  N = nrow(z)
  K = ncol(z)
  K.excl = ncol(as.matrix(z.excl))
  d.2 =ncol(as.matrix(x.exo))  # dim of included exogenous regressros
  if (is.null(d.2)) d.2=0
  d.x = ncol(x)
  
  # If lambda (ld) is not defined and use the equal weight
  if (is.null(ld)){
    ld = rep(1/d.x, d.x)
  }

  # Construct projection matrices of each model and its weighted averages.
  P.all = array( NA, c(N, N, K.excl) )
  for ( i in ( 1 : K.excl ) ) {
    P.m = proj_m(z[,(1:(d.2+i))])    
    P.all[ , , i] = P.m
  }
  
  
  # Generate P(N,N,K) array, all projection matrices with averaging
  # Calculate u.hat.k
  aveP = array(0,c(N,N,K.excl))
  n_subset = gen_p_triangle(K.excl)[K.excl,]

  if ( R == Inf ){
    #
    #
    #   Correct the fortran code later 
    #   we need to include z.excl and x.exo separately 
    #
    #
    n_subset = gen_p_triangle(K.excl)[K.excl,]
    dyn.load("gen_aveP.so")
    gen_aveP = (.Fortran("gen_aveP", z=z, N=N, K=K.excl, aveP=aveP, n_subset=as.integer(n_subset)))
    aveP = gen_aveP$aveP
  } else {
    for (i.aveP in c(1,K.excl-1)){
      cat('No. of IVs =', i.aveP, '\n')
      cat('No. of Subsets =', n_subset[i.aveP], '\n')
      if (n_subset[i.aveP] <= R){
        sum.sub = matrix(0, N, N)
        subset_ind = combn(K.excl, i.aveP)
	      for ( i.subset in c(1:n_subset[i.aveP]) ){
          z.sub = z.excl[,subset_ind[,i.subset]]
          proj.sub = proj_m(cbind(x.exo, z.sub))
          sum.sub = sum.sub + proj.sub
        }
        aveP[ , , i.aveP] = sum.sub / n_subset[i.aveP]
      } else {  # if the number of subset if bigger than R, then we conduct only R radom selections
        cat('Random Selection', '\n')
        sum.sub = matrix(0, N, N)
        for (i.sub in (1:R)){
          sub.rnd = sort(sample(1:K.excl, i.aveP, replace=F))
          z.sub = z.excl[,sub.rnd]
          proj.sub = proj_m(cbind(x.exo, z.sub))
          sum.sub = sum.sub + proj.sub  
        }
        aveP[ , , i.aveP] = sum.sub / R
      }
    }
  }


    # Preliminary Regression to get sig2's
    prelim = pre.est(method.pre='one.step', y, x.end, x.exo, z.pre=Z.pre.all, z, P.all, ld, d.2)
    prelim = pre.est(method.pre='lasso', y, x.end, x.exo, z.pre=Z.pre.all, z, P.all, ld, d.2)
    sig2.eps = prelim$sig2.eps
    sig2.ld = prelim$sig2.ld
    sig.ld.eps = prelim$sig.ld.eps
    sig.u.eps = prelim$sig.u.eps
    H.inv = prelim$H.inv
    H.ld = H.inv %*% ld
    u.pre = prelim$u.pre
    f.pre = prelim$f.pre
    P.f.pre = proj_m(f.pre)
    
    # Calculate S.hat.ld.k for k=1,..., K
    S.hat.ld.k = array(0,K.excl-1)
    bias.k = array(0,K.excl-1)
    var.k = array(0,K.excl-1)
    
    
    for ( j in c(1:(K.excl-1)) ) {
    #for ( j in c(1:2) ) {
      P.K = aveP[ , , j]
      I.N = diag(N)
      # In the paper
      S.term.1 =  (t(x) %*% (I.N - P.K)%*% (I.N - P.K)  %*% x) / N - (t(u.pre) %*% (I.N - P.K) %*% (I.N - P.K) %*% u.pre ) / N
      S.term.2 =  (t(x) %*% (I.N - P.K)  %*% x) / N - (t(u.pre) %*% (I.N - P.K) %*% u.pre) / N 
      S.term.1.c =  ( (t(x) %*% (I.N - P.K)%*% (I.N - P.K)  %*% (x)) - (t(u.pre) %*% (I.N - P.K)%*% (I.N - P.K)  %*% (u.pre)) -2*t(f.pre)%*%(I.N - P.K)%*% (I.N - P.K) %*% u.pre) / N 
      S.term.2.c =  ( t(x) %*% (I.N - P.K)  %*% (x) - t(u.pre) %*% (I.N - P.K)  %*% (u.pre) -2*t(f.pre) %*% (I.N - P.K)  %*% (u.pre)    ) / N 
      # When we use all IVs. P.K becomes a projection matrix
      #S.term.1.c =  ( (t(x) %*% (I.N - P.K)%*% (I.N - P.K)  %*% (x)) ) / N
      #S.term.2.c =  ( t(x) %*% (I.N - P.K)  %*% (x) + t(u.pre) %*% (I.N - P.K) )/ N
      mid.term =  t(f.pre) %*% (I.N - P.K)%*%(I.N - P.f.pre)%*%(I.N - P.K) %*% f.pre / N 
      mid.term.2 = S.term.1 - S.term.2 %*% H.inv %*% S.term.2
      mid.term.2.c = S.term.1.c - S.term.2.c %*% H.inv %*% S.term.2.c 
      S.hat.ld.k[j] = sig.ld.eps^2 * j^2/N + sig2.eps * ( t(H.ld) %*%  mid.term %*% H.ld )
      #S.hat.ld.k[j] = sig.ld.eps^2 * j^2/N + sig2.eps * ( t(H.ld) %*%  mid.term.2 %*% H.ld )
      bias.k[j] = sig.ld.eps^2 * j^2/N
      var.k[j] = sig2.eps * ( t(H.ld) %*%  mid.term.2 %*% H.ld )
      cat('No of k = ', j, '\n')
      cat('bias term = ', sig.ld.eps^2 * j^2/N, '\n')
      cat('variance term (quad)      = ', sig2.eps * ( t(H.ld) %*%  mid.term %*% H.ld ), '\n')
      cat('variance term (subt)      = ', sig2.eps * ( t(H.ld) %*%  mid.term.2 %*% H.ld ), '\n')
      cat('variance term (subt.c)    = ', sig2.eps * ( t(H.ld) %*%  mid.term.2.c %*% H.ld ), '\n')
      cat('Diff b/w two var est=', abs(sig2.eps * ( t(H.ld) %*%  mid.term.2 %*% H.ld ) - sig2.eps * ( t(H.ld) %*%  mid.term %*% H.ld )), '\n')
    }
    plot(S.hat.ld.k, type='l', col='green', ylim=c(0,1))
    lines(bias.k, col='blue')
    lines(var.k, col='red')
    opt.k = which.min(S.hat.ld.k)
    cat('opt.k =', opt.k, '\n')
  

  bt.hat.LS =  solve( t(x) %*% aveP[ , , opt.k] %*% x ) %*% ( t(x) %*% aveP[ , , opt.k] %*%y )
  
  cat('bt.hat.LS = ',bt.hat.LS,'\n')
  # cat('bt.hat.inf.LS = ',bt.hat.inf.LS,'\n')
  
  return(list(bt.hat=bt.hat.LS, opt.k=opt.k, obj.v=S.hat.ld.k,aveP=aveP))
}

tsls_CSA_get_P = function(y, x.end, x.exo, z.excl, R=Inf, ld=NULL, sub.K){
  
  z = cbind(x.exo, z.excl)
  x = as.matrix(cbind(x.end, x.exo))
  d.x =ncol(x)
  d.2 =ncol(as.matrix(x.exo))
  N = nrow(z)
  K = ncol(z)
  
  # If lambda (ld) is not defined and use the equal weight
  if (is.null(ld)){
    ld = rep(1/d.x, d.x)
  }
  
  # Generate P(N,N,K) array, all projection matrices with averaging
  # Calculate u.hat.k
  K.excl = K - d.2
  aveP = array(0,c(N,N,K.excl))
  n_subset = gen_p_triangle(K.excl)[K.excl,]
  
  if ( R == Inf ){
    #
    #
    #   Correct the fortran code later 
    #   we need to include z.excl and x.exo separately 
    #
    #
    n_subset = gen_p_triangle(K.excl)[K.excl,]
    dyn.load("gen_aveP.so")
    gen_aveP = (.Fortran("gen_aveP", z=z, N=N, K=K.excl, aveP=aveP, n_subset=as.integer(n_subset)))
    aveP = gen_aveP$aveP
  } else {
    for (i.aveP in sub.K){
      cat('No. of IVs =', i.aveP, '\n')
      cat('No. of Subsets =', n_subset[i.aveP], '\n')
      if (n_subset[i.aveP] <= R){
        sum.sub = matrix(0, N, N)
        subset_ind = combn(K.excl, i.aveP)
        for ( i.subset in c(1:n_subset[i.aveP]) ){
          z.sub = z.excl[,subset_ind[,i.subset]]
          proj.sub = proj_m(cbind(x.exo, z.sub))
          sum.sub = sum.sub + proj.sub
        }
        aveP[ , , i.aveP] = sum.sub / n_subset[i.aveP]
      } else {  # if the number of subset if bigger than R, then we conduct only R radom selections
        cat('Random Selection', '\n')
        sum.sub = matrix(0, N, N)
        for (i.sub in (1:R)){
          sub.rnd = sort(sample(1:K.excl, i.aveP, replace=F))
          z.sub = z.excl[,sub.rnd]
          proj.sub = proj_m(cbind(x.exo, z.sub))
          sum.sub = sum.sub + proj.sub  
        }
        aveP[ , , i.aveP] = sum.sub / R
      }
    }
  }
  return(list(sub.K=sub.K, aveP=aveP))
}


tsls_CSA_aveP_in = function(y, x.end, x.exo, z.excl, ld=NULL, aveP, method.pre, z.pre){
  
  
  # Declare constants
  z = as.matrix(cbind(x.exo,z.excl))
  x = cbind(x.end, x.exo)
  iv.pre = cbind(x.exo, z.pre)
  N = nrow(z)
  K = ncol(z)
  K.excl = ncol(as.matrix(z.excl))
  d.2 =ncol(as.matrix(x.exo))  # dim of included exogenous regressros
  if (is.null(d.2)) d.2=0
  d.x = ncol(x)
  
  # If lambda (ld) is not defined and use the equal weight
  if (is.null(ld)){
    ld = rep(1/d.x, d.x)
  }
  
  # Construct projection matrices of each model and its weighted averages.
  P.all = array( NA, c(N, N, K.excl) )
  for ( i in ( 1 : K.excl ) ) {
    P.m = proj_m(z[,(1:(d.2+i))])    
    P.all[ , , i] = P.m
  }
  
  cat('ld = ', ld, '\n')
  cat('d.x = ', d.x, '\n')
  cat('ld = ', length(ld), '\n')


  # Preliminary Regression to get sig2's
  prelim = pre.est(method.pre=method.pre, y, x.end, x.exo, z.pre=z.pre, z, P.all, ld, d.2)

  sig2.eps = prelim$sig2.eps
  sig2.ld = prelim$sig2.ld
  sig.ld.eps = prelim$sig.ld.eps
  sig.u.eps = prelim$sig.u.eps
  H.inv = prelim$H.inv
  H.ld = H.inv %*% ld
  u.pre = prelim$u.pre
  f.pre = prelim$f.pre
  P.f.pre = proj_m(f.pre)
  
  # Calculate S.hat.ld.k for k=1,..., K
  S.hat.ld.k = array(0,K.excl-1)
  bias.k = array(0,K.excl-1)
  var.k = array(0,K.excl-1)
  
  
  for ( j in c(1:(K.excl-1)) ) {
    #for ( j in c(1:2) ) {
    P.K = aveP[ , , j]
    I.N = diag(N)
    # In the paper
    S.term.1 =  (t(x) %*% (I.N - P.K)%*% (I.N - P.K)  %*% x) / N - (t(u.pre) %*% (I.N - P.K) %*% (I.N - P.K) %*% u.pre ) / N
    S.term.2 =  (t(x) %*% (I.N - P.K)  %*% x) / N - (t(u.pre) %*% (I.N - P.K) %*% u.pre) / N 
    S.term.1.c =  ( (t(x) %*% (I.N - P.K)%*% (I.N - P.K)  %*% (x)) - (t(u.pre) %*% (I.N - P.K)%*% (I.N - P.K)  %*% (u.pre)) -2*t(f.pre)%*%(I.N - P.K)%*% (I.N - P.K) %*% u.pre) / N 
    S.term.2.c =  ( t(x) %*% (I.N - P.K)  %*% (x) - t(u.pre) %*% (I.N - P.K)  %*% (u.pre) -2*t(f.pre) %*% (I.N - P.K)  %*% (u.pre)    ) / N 
    # When we use all IVs. P.K becomes a projection matrix
    #S.term.1.c =  ( (t(x) %*% (I.N - P.K)%*% (I.N - P.K)  %*% (x)) ) / N
    #S.term.2.c =  ( t(x) %*% (I.N - P.K)  %*% (x) + t(u.pre) %*% (I.N - P.K) )/ N
    mid.term =  t(f.pre) %*% (I.N - P.K)%*%(I.N - P.f.pre)%*%(I.N - P.K) %*% f.pre / N 
    mid.term.2 = S.term.1 - S.term.2 %*% H.inv %*% S.term.2
    mid.term.2.c = S.term.1.c - S.term.2.c %*% H.inv %*% S.term.2.c 
    S.hat.ld.k[j] = sig.ld.eps^2 * j^2/N + sig2.eps * ( t(H.ld) %*%  mid.term %*% H.ld )
    #S.hat.ld.k[j] = sig.ld.eps^2 * j^2/N + sig2.eps * ( t(H.ld) %*%  mid.term.2 %*% H.ld )
    bias.k[j] = sig.ld.eps^2 * j^2/N
    var.k[j] = sig2.eps * ( t(H.ld) %*%  mid.term.2 %*% H.ld )
    cat('No of k = ', j, '\n')
    cat('bias term = ', sig.ld.eps^2 * j^2/N, '\n')
    cat('variance term (quad)      = ', sig2.eps * ( t(H.ld) %*%  mid.term %*% H.ld ), '\n')
    cat('variance term (subt)      = ', sig2.eps * ( t(H.ld) %*%  mid.term.2 %*% H.ld ), '\n')
    cat('variance term (subt.c)    = ', sig2.eps * ( t(H.ld) %*%  mid.term.2.c %*% H.ld ), '\n')
    cat('Diff b/w two var est=', abs(sig2.eps * ( t(H.ld) %*%  mid.term.2 %*% H.ld ) - sig2.eps * ( t(H.ld) %*%  mid.term %*% H.ld )), '\n')
  }
  plot(S.hat.ld.k, type='l', col='green', ylim=c(0,1))
  lines(bias.k, col='blue')
  lines(var.k, col='red')
  opt.k = which.min(S.hat.ld.k)
  cat('opt.k =', opt.k, '\n')
  
  
  bt.hat.LS =  solve( t(x) %*% aveP[ , , opt.k] %*% x ) %*% ( t(x) %*% aveP[ , , opt.k] %*%y )
  
  cat('bt.hat.LS = ',bt.hat.LS,'\n')
  # cat('bt.hat.inf.LS = ',bt.hat.inf.LS,'\n')
  
  return(list(bt.hat=bt.hat.LS, opt.k=opt.k, obj.v=S.hat.ld.k,aveP=aveP))
}



# -----------------------------------------------------------
# USAGE: 
#     R_hat_cv( K, x, z, H.tilde ) 
#
# PURPOSE:
#     Calculate the objective function value of CV 
#
# INPUT:
#     K: the number of IVs where the ftn will be evaluated
#     x: N x 1, endogenous variable
#     z: N x K; instrudments excluding e.exo
#
# OUPUT:
#     the evaluated funcation value
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#    2/12/17      Y. shin             Original Code       
#    [Consider a multidensional delta and add 'ld' laster]
#
# -----------------------------------------------------------
R_hat_cv = function( K, x, z, H.tilde) {

  # Declare some constants
  z = as.matrix(z)
  N = nrow(z)
  
  z.K = z[,(1:K)]
  P.K = proj_m(z.K)
  u.hat.K = ( diag(N) - P.K ) %*% x
  u.hat.ld.K = u.hat.K %*% solve(H.tilde)
  
  R.hat.cv.value=mean( u.hat.ld.K^2/((1-diag(P.K))^2) )
  
  return(R.hat.cv.value)
}



# -----------------------------------------------------------
# USAGE: 
#     R_hat_m( K, x, z, H.tilde, sig.hat.ld ) 
#
# PURPOSE:
#     Calculate the objective function value of the Mallows criterion 
#
# INPUT:
#     K: the number of IVs where the ftn will be evaluated
#     x: N x 1, endogenous variable
#     z: N x K; instrudments excluding e.exo
#
# OUPUT:
#     the evaluated funcation value
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#   2/12/2017     Y. Shin             Original Code       
#    [Consider a multidensional delta and add 'ld' laster]
#
# -----------------------------------------------------------
R_hat_mallows = function( K, x, z, H.tilde, sig.hat.ld, ld ) {

  # Declare constants
  z = as.matrix(z)
  N = nrow(z)
  
  z.K = z[,(1:K)]
  P.K = proj_m(z.K)
  u.hat.K = (diag(N)-P.K) %*% x
  u.hat.ld.K = u.hat.K %*% solve(H.tilde) %*% ld
  
  R.hat.m.value = ( ( t(u.hat.ld.K) %*% u.hat.ld.K ) / N ) + ( ( 2 * sig.hat.ld * K )/ N )
  return(R.hat.m.value)
}


# -----------------------------------------------------------
# USAGE: 
#     first_stage_GOF( y, x, z, criteria )
#     # first_stage_GOF( y, x.end, x.exo, z, criteria )
#
# PURPOSE:
#     Calculate the preliminary number of instruments K.opt 
#
# INPUT:
#     y: N x 1, a dependent variable
#     x: N x 1, endogenous regressors
#     # x.end: N x d1, endogenous regressors
#     # x.exo: N x d2, exogenous regressors
#     z: N x K; instrudments excluding e.exo
#     H.tilde: d x d, matrix, for fixed K.
#     criteria: characters, 'CV' (cross-validation) or 'Mallows'
#
# OUPUT:
#     opt.K: a scalar, the preliminary number of instruments
#     obj.value: a scalar, the objective function values. The 'opt.K'
#           is the minimizer of this function valuef
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#   2017-02-12    Y. shin             Original Code       
#   2017-10-17                        Knitro Part error corrected (get_weight)
#				
#    [Consider a multidensional delta and add 'ld' laster]
#
#
# -----------------------------------------------------------
first_stage_GOF <- function( y, x, z, H.tilde, sig.hat.ld, criteria, ld, d.2 ) {

  z = as.matrix(z)  # z matrix should be [x.exo, z.excl]
  N = nrow(z)
  K = ncol(z)
  
  obj.value = rep(NA,K)
  
  # When the criteria is Cross-validations
  if ( criteria == 'CV' ){
    for ( i in (1:K) ){
      obj.value[i] = R_hat_cv( K=i, x=x, z=z, H.tilde=H.tilde )
    }
  } else if ( criteria == 'Mallows' ) {
    for (i in (1:K)){
      obj.value[i] = R_hat_mallows( K=i, x=x, z=z, H.tilde=H.tilde, sig.hat.ld=sig.hat.ld, ld=ld )
    }
  }
  opt.K = which.min(obj.value[-c(1:d.2)]) + d.2   # Find the minimum value outside (1:d.2)
  
  return( list( opt.K=opt.K, obj.value=obj.value ) )
}


# -----------------------------------------------------------
# This functions returns estimates for residuals required for
# estimating a MSE estimate.
# -----------------------------------------------------------
get_sigma_DN <- function ( y, x, z, opt.K, H.tilde, ld ) {
  N = length(y)
  Z.1st = z[,(1:opt.K)]
  P.1st = proj_m(Z.1st)
  bt.1st = tsls_est(y=y, x.end=x, x.exo=NULL, z=Z.1st)
  eps.tilde = y - x %*% bt.1st
  u.tilde = ( diag( N ) - P.1st ) %*% x
  u.ld = u.tilde %*% solve(H.tilde) %*% ld
  sig.eps2 = ( t(eps.tilde) %*% eps.tilde ) / N
  sig.ld2 = ( t(u.ld) %*% u.ld ) / N
  sig.ld.eps2 = ( t(u.ld) %*% eps.tilde ) / N
  
  return( list( eps2=sig.eps2, ld2=sig.ld2, ld.eps2=sig.ld.eps2, bt.1st=bt.1st ) )
}


# -----------------------------------------------------------
# USAGE: 
#     tsls_est(y, x.end, x.exo, z, ... )
#
# PURPOSE:
#     Calculate the two-stage least squares estimator bt.hat for
#         y = x' bt + eps
#         where x = ( x.end, x.exo )'
#
# INPUT:
#     y: N x 1, a dependent variable
#     x.end: N x d1, endogenous regressors
#     x.exo: N x d2, exogenous regressors
#     z: N x K; instrudments excluding e.exo
#     ... 
#     pi.0: pi.0 for calculating theoretical S(K)
#     sigma.ue: Cov(eps,u) for calculating theoretical S(K)
#     sigma.eps.2: Var(eps) for calculating theoretical S(K) 
#
# OUPUT:
#     est: d x 1, the estimate for bt
#     k.opt: a scalar, the optimal K.star.hat, estimated value, arg.max S.hat(K)
#     k.opt.inf: a scalar, the optimal K.star, infeasible theretical value, arg.max S(K)
# 
# REVISION LOG:
#     Date        Programmer          Description
#     ====        ==========          ===========
#   2/13/2017     Y. Shin             Original Code       
#   2/16/2017     Y. Shin             For R.hat.ld.K, use updated sig.ld instead of sig.hat.ld.pre (for Mallows)
#
# -----------------------------------------------------------
tsls_DN <- function(y, x.end, x.exo, z.excl, ld=NULL, method.pre, z.pre) {
  # Required functions
  
  z = as.matrix(cbind(x.exo,z.excl))
  x = cbind(x.end,x.exo)
  iv.pre = cbind(x.exo, z.pre)
  N = nrow(z)
  K = ncol(z)
  K.excl = ncol(as.matrix(z.excl))
  d.2 =ncol(as.matrix(x.exo))  # dim of included exogenous regressros
  if (is.null(d.2)) d.2=0
  d.x = ncol(x)
  
  # If lambda (ld) is not defined and use the equal weight
  if (is.null(ld)){
    ld = rep(1/d.x, d.x)
  }

  # Collect all projection matrices
  P.all = array( NA, c(N, N, K.excl) )
  for ( i in ( 1 : K.excl ) ) {
    P.m = proj_m(z[,(1:(d.2+i))])    
    P.all[ , , i] = P.m
  }
  
  
  prelim = pre.est(method.pre, y, x.end, x.exo, z.pre, z, P.all, ld, d.2)
  sig2.eps = prelim$sig2.eps
  sig2.ld = prelim$sig2.ld
  sig.ld.eps = prelim$sig.ld.eps
  H.inv = prelim$H.inv

  # Evaluate the obj. fun for different k
  S.DN.k = rep(NA, K.excl)
  for (i.DN in (1:K.excl)){
    P.k = P.all[,,i.DN]
    u.pre.k = ( diag(N) - P.k ) %*% x
    u.ld.k = u.pre.k %*% H.inv %*% ld
    S.DN.k[i.DN] = sig.ld.eps^2 * (i.DN^2 / N) + sig2.eps^2 * ( (t(u.ld.k) %*% u.ld.k) / N + sig2.ld *(i.DN/N) )
  } 

  

  k.opt = which.min(S.DN.k)   
  z.opt = z[,(1:(d.2+k.opt))]
  P.DN = proj_m(z.opt)
  bt.DN = tsls_est( y=y, x.end=x.end, x.exo=x.exo, z=z.opt[,-c(1:d.2)])
  cat('k.opt = ', k.opt,'\n')
  cat('bt.hat.DN = \n', bt.DN,'\n')
  return(list(opt.k=k.opt, bt.hat=bt.DN, P.DN=P.DN))
}  


#------------------------------------------------------------------------------
# Functions required for KO

get_weights = function(y, x.end, x.exo, z.excl, H.tilde, sig.hat.ld, P.all, Gm, initial.z='Mallows', ld){
  library("quadprog")
  
  # Declare constants 
  z = cbind(x.exo, z.excl)
  K = ncol(z)
  N = nrow(z)
  d.2 =ncol(as.matrix(x.exo))
  
  x = cbind(x.end, x.exo)
  opt.w = rep(NA, K)
  
  # To get the preliminary estimate bt.tilde
  #r1 = first_stage_GOF(y=y, x=x, z=z, H.tilde=H.tilde, sig.hat.ld=sig.hat.ld, criteria='Mallows')
  #r1 = first_stage_GOF( y=y, x=x, z=z, H.tilde=H.tilde, sig.hat.ld=sig.hat.ld, criteria='Mallows', ld=ld, d.2=d.2 )
  
  #K.tilde = r1$opt.K
  #if (GF=='all') {K.tilde=K}
  #cat('KO estimator: K.tilde =', K.tilde, '\n')
  
  # Get preliminary K.opt using the first stage Goodness of fit criteria
  if (initial.z == 'all') {
    K.tilde=K
  } else {
    r1 = first_stage_GOF( y=y, x=x, z=z, H.tilde=H.tilde, sig.hat.ld=sig.hat.ld, criteria=initial.z, ld=ld, d.2=d.2 )
    K.tilde = r1$opt.K
  }
  cat('KO estimator: K.tilde = ',K.tilde,'\n')
  
  
  # Update H.tilde using K.tilde
  H.tilde.update = ( t(x) %*% P.all[ , , K.tilde] %*% x ) / N
  
  # Calculate sig.hat.eps, sig.hat.ld, sig.hat.ld.eps, H.hat using 'get_sigma_DN'
  #r2 = get_sigma_DN(y=y ,x=x ,z=z, opt.K=K.tilde, H.tilde=H.tilde.update)
  r2 = get_sigma_DN(y=y, x=x, z=z, opt.K=K.tilde, H.tilde=H.tilde.update, ld)
  
  sig.hat.eps = as.numeric(r2$eps2)
  sig.hat.ld = as.numeric(r2$ld2)
  sig.hat.ld.eps = as.numeric(r2$ld.eps2)
  
  
  # Calculate u.hat.ld and U.hat
  u.stack = matrix(0, N, K)
  for (i in (1:K)){
    u.stack[ , i]=( P.all[ , , K] - P.all[ , , i] ) %*% x %*% solve(H.tilde.update) %*% ld
  }
  U.hat = t(u.stack) %*% u.stack              # (K x K) matrix for S(K)
  
  
  
  # Generate a (K x 1) vector K.vec s.t. K.vec[i]=i
  K.vec = c(1:K)
  
  # Generate a (K x 1) vector of 1's
  one.K = rep(1,K)
  
  # Generate a (K x 1) vector of 0's
  zero.K = rep(0,K)
  
  
  #  This expression comes from Equation (2.6) in KO. 
  #  c = as.numeric(2 * sig.hat.eps * sig.hat.ld * t(K.vec)  ) /N
  #  H = 2 * ( sig.hat.ld.eps^2 * (K.vec %*% t(K.vec)  ) + (sig.hat.eps * U.hat  ) - (sig.hat.eps * sig.hat.ld * Gm )  ) / N
  
  # Equation (2.5) is rewritten: H is a matrix for the quadratic term; c is a matrix for the linear term
  c = as.numeric(-(B.ld.N) * t(K.vec) ) /N
  H = ( sig.hat.ld.eps^2 * (K.vec %*% t(K.vec)  ) +  (sig.hat.ld.eps^2 * Gm )  + (sig.hat.eps * U.hat ) ) / N
  
  Dmat = 2*H
  sc = norm(Dmat, "2")
  dvec = -c
  # Rescale the objective function to resolve the overflow issue
  Dmat = Dmat / sc
  dvec = dvec / sc
  Amat = rbind(one.K, diag(one.K), -diag(one.K))
  Amat = t(Amat)   
  bvec = c(1, zero.K, -one.K)
  
  r.qp=solve.QP(Dmat, dvec, Amat, bvec, meq=1)
  
  obj.val=r.qp$value
  opt.w = round(r.qp$solution,4)
  
  
  return(list(opt.w=opt.w, obj.val=obj.val))
}





#obj.fnt=function(x, c, H) {c %*% x + 0.5* t(x) %*% H %*% x }


get_weights_knitro = function(y, x.end, x.exo, z, H.tilde, sig.hat.ld, P.all, Gm){
  library("KnitroR")
  
  # Declare constants 
  z = as.matrix(z)
  K = ncol(z)
  N = nrow(z)
  
  x = cbind(x.end, x.exo)
  opt.w = rep(NA, K)
  
  # To get the preliminary estimate bt.tilde
  r1 = first_stage_GOF(y=y, x=x, z=z, H.tilde=H.tilde, sig.hat.ld=sig.hat.ld, criteria='Mallows')
  r1 = first_stage_GOF(y=y, x=x, z=z, H.tilde=1, sig.hat.ld=sig.hat.ld*H.tilde^2, criteria='Mallows')
  K.tilde = r1$opt.K
  # R.hat.all = r1$obj.value    # (K x 1) vector for R.hat (K), required
  
  # Update H.tilde using K.tilde
  H.tilde.update = ( t(x) %*% P.all[ , , K.tilde] %*% x ) / N
  
  # Calculate sig.hat.eps, sig.hat.ld, sig.hat.ld.eps, H.hat using 'get_sigma_DN'
  r2 = get_sigma_DN(y=y ,x=x ,z=z, opt.K=K.tilde, H.tilde=H.tilde.update)
  sig.hat.eps = as.numeric(r2$eps2)
  sig.hat.ld = as.numeric(r2$ld2)
  sig.hat.ld.eps = as.numeric(r2$ld.eps2)
  
  
  # Calculate u.hat.ld and U.hat
  u.stack = matrix(0, N, K)
  for (i in (1:K)){
    u.stack[ , i]=( P.all[ , , K] - P.all[ , , i] ) %*% x %*% solve(H.tilde.update)
  }
  U.hat = t(u.stack) %*% u.stack              # (K x K) matrix for S(K)
  
  
  
  # Generate a (K x 1) vector K.vec s.t. K.vec[i]=i
  K.vec = c(1:K)
  

  
  #--------------------------------
  # QP problem using Knitro
  #---------------------------------
  
  # These parameters and functions are dropping (H^{-1})^2 as in the Ox codes of KO.
  # They should be equivalent to the one that uses the original formular.
  # If you want to use these formular then change the objective function in the knitro call from obj.fnt to obj.fnt.KO
  #--------------------------------
  # c.KO = as.numeric(2 * sig.hat.eps * sig.hat.u * t(K.vec)  ) 
  # H.KO = 2 * ( sig.hat.u.eps^2 * (K.vec %*% t(K.vec)  ) + (sig.hat.eps * U.hat.KO  ) - (sig.hat.eps * sig.hat.u * Gm ) ) 
  # obj.fnt.KO=function(x) {t(c.KO) %*% x + 0.5* t(x) %*% H.KO %*% x }
  #---------------------------------------------
  
  c.fnt = function(x) {
    K.ones = rep(1,K)
    return(K.ones %*% x)
  }
#  This expression comes from Equation (2.6) in KO. 
#  c = as.numeric(2 * sig.hat.eps * sig.hat.ld * t(K.vec)  ) /N
#  H = 2 * ( sig.hat.ld.eps^2 * (K.vec %*% t(K.vec)  ) + (sig.hat.eps * U.hat  ) - (sig.hat.eps * sig.hat.ld * Gm )  ) / N

# Equation (2.5) is rewritten: H is a matrix for the quadratic term; c is a matrix for the linear term
  c = as.numeric(   -2*t(K.vec)*(sig.hat.eps*sig.hat.ld + 4*sig.hat.ld.eps^2)  + 2 * sig.hat.eps * sig.hat.ld * t(K.vec)  ) /N
  H = 2 * ( sig.hat.ld.eps^2 * (K.vec %*% t(K.vec)  ) + (sig.hat.eps * U.hat  ) + (sig.hat.ld.eps^2 * Gm )  ) / N


  obj.fnt=function(x) {t(c) %*% x + 0.5* t(x) %*% H %*% x }


  eval_grad_f <- function(x) {
    grad_f = H %*% x + c
    return( grad_f )
  }

# Jacobian callback
  eval_jac_g <- function( x ) {
    return( rep(1,K) )
  }

  cL = c(1.0)
  cU = rep(1.0)

  xL = rep(0,K)
  xU = rep(1,K)

  #
  # Knitro algorithms: 1 = Interior/Direct
  #                    2 = Interior/CG
  #                    3 = Active Set
  #                    4 = SQR
  # ms_enable: multiple starting points, T or F
  # URL: https://www.artelys.com/tools/knitro_doc/2_userGuide/algorithms.html
  #

  #r7=knitro(nvar=K, ncon=K+1, nnzJ=-1, x0=w0, objective = obj.fnt, constraints=c.fnt ,cL=cL, cU=cU, options=knitro.opt)
  r7=knitro(nvar=K, ncon=1, gradient=eval_grad_f,  jacobian=eval_jac_g, objective = obj.fnt, constraints=c.fnt ,cL=cL,  cU=cU, xL=xL, xU=xU)
  sol7=r7$x
  obj.val=obj.fnt(sol7)
  opt.w = round(sol7,4)

  
  return(list(opt.w=opt.w, obj.val=obj.val))
}





tsls_KO <- function(y, x.end, x.exo, z.excl, ld=NULL, method.pre, z.pre) {
  library('quadprog')

  # Declare constants
  z = as.matrix(cbind(x.exo,z.excl))
  x = cbind(x.end, x.exo)
  iv.pre = cbind(x.exo, z.pre)
  N = nrow(z)
  K = ncol(z)
  K.excl = ncol(as.matrix(z.excl))
  d.2 =ncol(as.matrix(x.exo))  # dim of included exogenous regressros
  if (is.null(d.2)) d.2=0
  d.x = ncol(x)
  
  # Generate a (K x 1) vector K.vec s.t. K.vec[i]=i
  K.vec = c(1:K.excl)
  
  # Generate a (K x 1) vector of 1's
  one.K = rep(1,K.excl)
  
  # Generate a (K x 1) vector of 0's
  zero.K = rep(0,K.excl)
  
  # Generate Gamma matrix such that GM[i,j] = min{i,j}
  Gm = matrix(0, K.excl, K.excl)
  for ( i in (1:K.excl) ){
    for ( j in (1:K.excl) ) {
      Gm[i,j] = min(c(i,j)) 
    }
  }
  
  # Construct projection matrices of each model and its weighted averages.
  P.all = array( NA, c(N, N, K.excl) )
  for ( i in ( 1 : K.excl ) ) {
    P.m = proj_m(z[,(1:(d.2+i))])    
    P.all[ , , i] = P.m
  }
  
  # If lambda (ld) is not defined and use the equal weight
  if (is.null(ld)){
    ld = rep(1/d.x, d.x)
  }
  
  # Preliminary Regression to get sig2's
  prelim = pre.est(method.pre, y, x.end, x.exo, z.pre, z, P.all, ld, d.2)
  sig2.eps = prelim$sig2.eps
  sig2.ld = prelim$sig2.ld
  sig.ld.eps = prelim$sig.ld.eps
  sig.u.eps = prelim$sig.u.eps
  H.inv = prelim$H.inv
  u.pre = prelim$u.pre
  f.pre = prelim$f.pre
  
  #-----------------------------------------------------------------------
  # Evaluate the obj. fun 
  #-----------------------------------------------------------------------
  # Calculate U.hat
  u.stack = matrix(0, N, K.excl)
  for (i in (1:K.excl)){
    u.stack[ , i]=( P.all[ , , K.excl] - P.all[ , , i] ) %*% x %*% H.inv %*% ld
  }
  U.hat = t(u.stack) %*% u.stack  
  
  # Calculate Sig.u
  Sig.u = (t(u.pre) %*% u.pre) / N
  
  # Calculate B.ld.N
  B.N.1 = sig2.eps * Sig.u
  B.N.2 = d.x * (sig.u.eps %*% t(sig.u.eps))
  B.N.3 = matrix(0, d.x, d.x)
  for (i.B in (1:N)){
    f.i = f.pre[i.B,]
    term.1 = f.i %*% t(sig.u.eps) %*% H.inv %*% sig.u.eps %*%t(f.i)
    term.2 = f.i %*% t(sig.u.eps) %*% H.inv %*% f.i %*% t(sig.u.eps)
    term.3 = sig.u.eps %*% t(f.i) %*% H.inv %*% sig.u.eps %*% t(f.i)
    B.N.3 = B.N.3 + term.1 + term.2 + term.3 
  }
  B.N.3 = B.N.3 / N
  B.N = 2*(B.N.1 + B.N.2 + B.N.3)
  B.ld.N = drop( t(ld) %*% H.inv %*% B.N %*% H.inv %*% ld )
  
  # Equation (2.5) is rewritten: H is a matrix for the quadratic term; c is a matrix for the linear term
  c = as.numeric(-(B.ld.N) * t(K.vec) ) /N
  H = ( sig.ld.eps^2 * (K.vec %*% t(K.vec)  ) +  (sig.ld.eps^2 * Gm )  + (sig2.eps * U.hat ) ) / N
  
  Dmat = 2*H
  sc = norm(Dmat, "2")
  dvec = -c
  # Rescale the objective function to resolve the overflow issue
  Dmat = Dmat / sc
  dvec = dvec / sc
  Amat = rbind(one.K, diag(one.K), -diag(one.K))
  Amat = t(Amat)   
  bvec = c(1, zero.K, -one.K)
  
  r.qp=solve.QP(Dmat, dvec, Amat, bvec, meq=1)
  
  obj.val=r.qp$value
  opt.w = round(r.qp$solution,4)

  cat('Optimal Weights = ', opt.w,'\n')

  # Calculate the averaged projection matrix with the optimal weight. 
  P.W = matrix(0, N, N)
  for ( i in ( 1 : K.excl ) ) {
    P.W = P.W + opt.w[i] * P.all[ , ,i]
  }
  
  # Calculate the 2SLS estimator with P.W and return the result
  bt.hat=solve(t(x)%*%P.W%*%x)%*%(t(x)%*%P.W%*%y)
  cat('bt.hat (KO) = ', bt.hat, '\n')

  return(list(bt.hat=bt.hat,opt.w=opt.w, obj.val=obj.val, P.KO=P.W))  
}


cluster.se = function(id, bt, x, y, P){
  eps.hat = y - x %*% bt
  id.index = list()
  P.g = list()
  eps.g.hat = list()
  
  G = length(table(id))
  d = ncol(x)
  mid.sum.G = matrix(0,d,d)
  
  for (i in (1:G)){
    id.index[[i]] = which(id==i)
    P.g[[i]] = P[id.index[[i]],]
    eps.g.hat[[i]] = eps.hat[id.index[[i]]]
    mid.sum.G = mid.sum.G + t(x) %*% t(P.g[[i]]) %*% eps.g.hat[[i]] %*% t(eps.g.hat[[i]]) %*% P.g[[i]] %*% x
  }
  
  outer.part = solve(t(x) %*% P %*% x)
  Sig.hat = N*outer.part %*% (mid.sum.G) %*% outer.part
  vcov = Sig.hat/N
  se = sqrt(diag(vcov))
  cat('hetero-cluster robust se =', se,'\n')
  return(list(se=se, vcov=vcov))
}


pre.est=function(method.pre, y, x.end, x.exo, z.pre, z, P.all, ld, d.2){
  cat('---------------------- \n')
  cat('Preliminary Estimation \n')
  cat('---------------------- \n')
  cat('Method = ', method.pre, '\n')
  
  x = cbind(x.end, x.exo)
  d.x = ncol(as.matrix(x))
  
  N = nrow(x)
  iv.pre = cbind(x.exo, z.pre)
  
  if (method.pre =='one.step'){
    P.pre = proj_m(iv.pre)
    f.pre = P.pre %*% x
    H.pre = ( t(x) %*% P.pre %*% x ) / N
    H.inv = solve(H.pre)
    u.pre = ( diag(N) - P.pre ) %*% x
    u.ld.pre = u.pre %*% H.inv %*% ld
    
    bt.pre = tsls_est(y,x.end,x.exo,z=z.pre)
    eps.pre = y - x %*% bt.pre
    sig2.eps = drop( ( t(eps.pre) %*% eps.pre ) / N )
    sig.ld.eps = drop( ( t(u.ld.pre) %*% eps.pre ) / N )
    sig2.ld = (t(u.ld.pre) %*% u.ld.pre) / N
    sig.u.eps = (t(u.pre) %*% eps.pre) / N
    
    cat('sig2.eps   = ', sig2.eps, '\n')
    cat('sig2.ld    = ', sig2.ld, '\n')
    cat('sig.ld.eps = ', sig.ld.eps, '\n')

  } else if (method.pre =='two.step'){
    P.pre = proj_m(iv.pre)
    H.pre = ( t(x) %*% P.pre %*% x ) / N
    H.inv = solve(H.pre)
    u.pre = ( diag(N) - P.pre ) %*% x
    u.ld.pre = u.pre %*% H.inv %*% ld
    sig2.ld = ( t(u.ld.pre) %*% u.ld.pre ) / N
    
    Mallows.k = rep(NA, K.excl)
    for (i.pre in (1:K.excl)){
      P.k = P.all[,,i.pre]
      u.pre.k = ( diag(N) - P.k ) %*% x
      u.ld.k = u.pre.k %*% H.inv %*% ld
      Mallows.k[i.pre] = (t(u.ld.k) %*% u.ld.k) / N + sig2.ld *(2*i.pre/N)
    }
    k.m = which.min(Mallows.k) # Choice of k after applying the Mallows criterion
    cat('k.m          =', k.m, '\n')
    
    # Re-estimate, eps, u, H with k.m
    z.m = z[,(1:(d.2+k.m))]
    P.m = P.all[,,k.m]
    f.pre = P.m %*% x
    H.m = ( t(x) %*% P.m %*% x ) / N
    H.inv = solve(H.m)
    u.m = u.pre = (diag(N) - P.m ) %*% x
    u.ld.m = u.m %*% H.inv %*% ld
    
    bt.m = tsls_est(y,x.end,x.exo,z=z.m[,-c(1:d.2)])
    eps.m = y - x %*% bt.m
    sig2.eps = drop( ( t(eps.m) %*% eps.m ) / N )
    sig.ld.eps = drop( ( t(u.ld.m) %*% eps.m ) / N )
    sig2.ld = t(u.ld.m) %*% u.ld.m / N
    sig.u.eps = (t(u.m) %*% eps.m) / N
    
    cat('sig2.eps   = ', sig2.eps, '\n')
    cat('sig.ld    = ', sig2.ld, '\n')
    cat('sig.ld.eps = ', sig.ld.eps, '\n')

  } else  if (method.pre=='lasso'){
    library('glmnet')
    if (ncol(as.matrix(x.end)==1)){
      cv.m = cv.glmnet(x=z, y=x.end)
      ld.1se = cv.m$lambda.1se
      m.first.stage = glmnet(x=z, y=x.end, family='gaussian', intercept=T, penalty.factor=c(rep(0,d.x),rep(1,K.excl)), lambda = ld.1se)  
      pi.1st = coef(m.first.stage)[-2]
      selected = which(pi.1st!=0)
      cat('1st Stage Lasso Selection =', selected, '\n')
      cat('No of selected z.excl     =', length(selected)-24, '\n')
      
      # Calculate eps, u, H
      z.sel = z[,selected]
      P.pre = proj_m(z.sel)
      f.pre = P.pre %*% x
      H.pre = ( t(x) %*% P.pre %*% x ) / N
      H.inv = solve(H.pre)
      u.pre = ( diag(N) - P.pre ) %*% x
      u.ld.pre = u.pre %*% H.inv %*% ld
      
      bt.pre = tsls_est(y,x.end,x.exo,z=z.sel[,-c(1:d.2)])
      eps.pre = y - x %*% bt.pre
      sig2.eps = drop( ( t(eps.pre) %*% eps.pre ) / N )
      sig.ld.eps = drop( ( t(u.ld.pre) %*% eps.pre ) / N )
      sig2.ld = (t(u.ld.pre) %*% u.ld.pre) / N
      sig.u.eps = (t(u.pre) %*% eps.pre) / N
      
      # Evaluate the obj. fun for different k
      cat('sig2.eps   = ', sig2.eps, '\n')
      cat('sig2.ld    = ', sig2.ld, '\n')
      cat('sig.ld.eps = ', sig.ld.eps, '\n')
    } else {
      cat('To be added \n')
      break
    }
    
  }
  
  return(list(sig2.eps=sig2.eps, sig2.ld=sig2.ld, sig.ld.eps=sig.ld.eps, sig.u.eps=sig.u.eps, H.inv=H.inv, u.pre=u.pre, f.pre=f.pre))
}

# Count the number of products that have inelastic demand
inelastic = function(ap, price, share){
  sum(ap*(price)*(1-share) > -1)
}
