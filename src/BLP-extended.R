rm(list=ls())

if (!require(hdm)) install.packages("hdm")
if (!require(sandwich)) install.packages("sandwich")
if (!require(lmtest)) install.packages("lmtest")
if (!require(R.matlab)) install.packages("R.matlab")
if (!require(doMC)) install.packages("doMC")

library(hdm)
library(sandwich)
library(lmtest)
library(R.matlab)
library(doMC)
source('lib_csa.R')

start.time=proc.time()[3]
# Set seed
set.seed(152732)

# Folder and File names
result.folder = "../result/"
data.folder = "../data/"
out.file.name = paste(result.folder, 'BLP-extended.out', sep='')

# Import the BLP data set
data(BLP)
dat = BLP$BLP
attach(dat)

# Design Selection
design = 'big' # small: dim(x)=5, dim(z)=10, big: dim(x)=24, dim(z)=48

# parallel computing
use.par = TRUE   # TURE or FALSE

# Number of random draws for the large complete subset
R = 1000

# Delcare Variables
# X: included exogenous variable
# Z: excluded exogenous variable
# price: the endogenous vriable

if (design =='small') {
  # X: Exogenous regressors
  X = as.matrix(cbind(1, dat[,c("hpwt", "air", "mpd","space")]))
  # Z: IVs
  Z = BLP$Z
} else if (design =='big') {
  tu <- trend/19
  mpdu <- mpd/7
  spaceu <- space/2
  augX = cbind(1, hpwt, air, mpdu, spaceu, tu, hpwt^2, hpwt^3, mpdu^2, mpdu^3, spaceu^2, spaceu^3, tu^2, tu^3, hpwt*air, mpdu*air, spaceu*air, tu*air,
               hpwt*mpdu, hpwt*spaceu, hpwt*tu, mpdu*spaceu, mpdu*tu, spaceu*tu)
  X = augX
  colnames(X) = c('const', 'hpwt', 'air', 'mpdu', 'spaceu', 'tu', 'hpwt^2', 'hpwt^3', 'mpdu^2', 'mpdu^3',
                  'spaceu^2', 'spaceu^3', 'tu^2', 'tu^3', 'hpwt*air', 'mpdu*air', 'spaceu*air', 'tu*air',
                  'hpwt*mpdu', 'hpwt*spaceu', 'hpwt*tu', 'mpdu*spaceu', 'mpdu*tu', 'spaceu*tu')
  Z = BLP$augZ
}
N = nrow(X)
K.excl = ncol(Z)

x.reg = cbind(price,X)

Z.sort = Z
cor.order = rank(-abs(cor(price,Z)))
Z.sort[,cor.order] = Z

Z.pre.one=Z.sort[,1]
Z.pre.all=Z


# Use multi cores / single core (no parallel)
if (use.par==T){
  n.core = detectCores()-1  
  aveP=tsls_CSA_get_P(y, x.end=price, x.exo=X, z.excl=Z.sort, R=R, ld=NULL, sub.K=c(1:K.excl), use.par = T, n.core=n.core)
} else {
  aveP=tsls_CSA_get_P(y, x.end=price, x.exo=X, z.excl=Z.sort, R=R, ld=NULL, sub.K=c(1:K.excl))
}

  
# OLS Estimation
sink(out.file.name)
cat('--------------------------------- \n')
cat('OLS \n')
cat('--------------------------------- \n')
m.ols = lm(y~ price + X-1, data=dat)
print(m.ols)
robust.inf.OLS = cluster.se(id=firm.id, bt=coef(m.ols), x=x.reg, y=y, P=diag(N))
se.OLS = robust.inf.OLS$se
cat('se ap = ', se.OLS[1],'\n')
sink()

# 2SLS Estimation
sink(out.file.name, append=T)
cat('--------------------------------- \n')
cat('2SLS \n')
cat('--------------------------------- \n')
m.tsls = tsls(x=X, d=price, y=y, z=Z, intercept =FALSE, homoscedastic = FALSE)
print(summary(m.tsls))
P.2SLS = aveP[[K.excl]] 
robust.inf.2SLS = cluster.se(id=firm.id, bt=coef(m.tsls), x=x.reg, y=y, P=P.2SLS)
se.2SLS = robust.inf.2SLS$se
cat('se ap = ', se.2SLS[1],'\n')
sink()

# DN
sink(out.file.name, append=T)
cat('--------------------------------- \n')
cat('Donald and Newey Estimator\n')
cat('--------------------------------- \n')
result_DN_two_step = tsls_DN( y=y, x.end=price, x.exo=X, z.excl=Z.sort, method.pre='two.step', z.pre= Z.pre.one)
bt.DN.two.step = result_DN_two_step$bt.hat
K.DN.two.step = result_DN_two_step$opt.k
robust.inf.DN.two.step = cluster.se(id=firm.id, bt=bt.DN.two.step, x=x.reg, y=y, P=result_DN_two_step$P.DN)
se.DN.two.step = robust.inf.DN.two.step$se
cat('se ap = ', se.DN.two.step[1],'\n')
sink()

# KO
sink(out.file.name, append=T)
cat('--------------------------------- \n')
cat('Kuersteiner and Okui Estimator \n')
cat('--------------------------------- \n')
# KO: two-step
result_KO_two_step = tsls_KO( y=y, x.end=price, x.exo=X, z.excl=Z.sort, method.pre='two.step', z.pre=Z.pre.one)
bt.KO.two.step = result_KO_two_step$bt.hat
opt.w.two.step = result_KO_two_step$opt.w 
P.KO.two.step = result_KO_two_step$P.KO
robust.inf.KO.two.step = cluster.se(id=firm.id, bt=bt.KO.two.step, x=x.reg, y=y, P=P.KO.two.step)
se.KO.two.step = robust.inf.KO.two.step$se
cat('se ap = ', se.KO.two.step[1],'\n')
sink()

# CSA: Complete Subset Average Estimation
sink(out.file.name, append=T)
cat('--------------------------------- \n')
cat('CSA \n')
cat('--------------------------------- \n')
# CSA: two-step
result_CSA_two_step = tsls_CSA_aveP_in( y=y, x.end=price, x.exo=X, z.excl=Z.sort, ld=NULL, aveP=aveP, method.pre='two.step', z.pre=Z.pre.one)
bt.CSA.two.step = result_CSA_two_step$bt.hat
K.CSA.two.step = result_CSA_two_step$opt.k
P.CSA.two.step = aveP[[K.CSA.two.step]]
robust.inf.CSA.two.step = cluster.se(id=firm.id, bt=bt.CSA.two.step, x=x.reg, y=y, P=P.CSA.two.step)
se.CSA.two.step = robust.inf.CSA.two.step$se
cat('se ap = ', se.CSA.two.step[1],'\n')
sink()

# Calculate the number of inelastic products
BLP_data=readMat(paste(data.folder,'BLP_data.mat',sep=""))
o_price = BLP_data$price
o_share = BLP_data$share

names.method=c('OLS', '2SLS', 'DN.two.step', 'KO.two.step', 'CSA.two.step')

# Construct the Table
ap.hat = c(coef(m.ols)['price'], coef(m.tsls)[1,1], bt.DN.two.step[1,1], bt.KO.two.step[1,1], bt.CSA.two.step[1,1])
No.Inelastic = sapply(ap.hat, inelastic, price=o_price, share=o_share)
K.opt = c(NA, NA, K.DN.two.step, NA, K.CSA.two.step)
robust.se = round(c(se.OLS[1], se.2SLS[1], se.DN.two.step[1], se.KO.two.step[1], se.CSA.two.step[1]),3 )


table.main = cbind(ap.hat, robust.se, K.opt, No.Inelastic)
row.names(table.main) = names.method

sink(out.file.name, append=T)
cat('\n \n')
cat('================== \n')
cat('BLP-extended-Table \n')
cat('================== \n')
print(table.main)
sink()

#save(table.main, file='BLP-original-table.RData')
#save.image(file=paste(result.folder,'BLP-original.RData',sep=''))
end.time=proc.time()[3]
exec.time = end.time-start.time
cat('Elapsed Time: ',exec.time,'\n')
