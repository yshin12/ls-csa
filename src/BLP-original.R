rm(list=ls())

if (!require(hdm)) install.packages("hdm")
if (!require(sandwich)) install.packages("sandwich")
if (!require(lmtest)) install.packages("lmtest")
if (!require(R.matlab)) install.packages("R.matlab")
if (!require(quadprog)) install.packages("quadprog")
if (!require(glmnet)) install.packages("glmnet")

library(hdm)
library(sandwich)
library(lmtest)
library(R.matlab)
source('lib_csa.R')

start.time=proc.time()[3]
# Set seed
set.seed(152732)

# Import the BLP data set
data(BLP)
dat = BLP$BLP
attach(dat)

# Design Selection
design = 'small' # small: dim(x)=5, dim(z)=10, big: dim(x)=24, dim(z)=48

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

# Options for method.pre
# one.step: Pick preliminary Z -> obj. fnt (Z will be picked as (i) min number of Z with the highest correlation or (ii) all IVs)
# two.step: preliminary Z -> Mallows -> obj.fnt (Z will be picked as (i) min number of Z with the highest correlation or (ii) all IVs)
# lasso: lasso selection for Z (no selection on the original eq) -> obj. fnt


# OLS Estimation

cat('--------------------------------- \n')
cat('OLS \n')
cat('--------------------------------- \n')
m.ols = lm(y~ price + X-1, data=dat)
print(m.ols)
robust.inf.OLS = cluster.se(id=firm.id, bt=coef(m.ols), x=x.reg, y=y, P=diag(N))
se.OLS = robust.inf.OLS$se
cat('se ap = ', se.OLS[1],'\n')


# 2SLS Estimation
cat('--------------------------------- \n')
cat('2SLS \n')
cat('--------------------------------- \n')
m.tsls = tsls(x=X, d=price, y=y, z=Z, intercept =FALSE, homoscedastic = FALSE)
print(summary(m.tsls))
P.2SLS = aveP[,,K.excl]
robust.inf.2SLS = cluster.se(id=firm.id, bt=coef(m.tsls), x=x.reg, y=y, P=P.2SLS)
se.2SLS = robust.inf.2SLS$se
cat('se ap = ', se.2SLS[1],'\n')


# DN
cat('--------------------------------- \n')
cat('Donald and Newey Estimator\n')
cat('--------------------------------- \n')
# DN: one-step 
result_DN_one_step = tsls_DN( y=y, x.end=price, x.exo=X, z.excl=Z.sort, method.pre='one.step', z.pre= Z.pre.one)
bt.DN.one.step = result_DN_one_step$bt.hat
K.DN.one.step = result_DN_one_step$opt.k
robust.inf.DN.one.step = cluster.se(id=firm.id, bt=bt.DN.one.step, x=x.reg, y=y, P=result_DN_one_step$P.DN)
se.DN.one.step = robust.inf.DN.one.step$se
cat('se ap = ', se.DN.one.step[1],'\n')

# DN: lasso
result_DN_lasso = tsls_DN( y=y, x.end=price, x.exo=X, z.excl=Z.sort, method.pre='lasso', z.pre= Z.pre.one)
bt.DN.lasso = result_DN_lasso$bt.hat
K.DN.lasso = result_DN_lasso$opt.k
robust.inf.DN.lasso = cluster.se(id=firm.id, bt=bt.DN.lasso, x=x.reg, y=y, P=result_DN_lasso$P.DN)
se.DN.lasso = robust.inf.DN.lasso$se
cat('se ap = ', se.DN.lasso[1],'\n')


# KO
cat('--------------------------------- \n')
cat('Kuersteiner and Okui Estimator \n')
cat('--------------------------------- \n')
# KO: one-step
result_KO_one_step = tsls_KO( y=y, x.end=price, x.exo=X, z.excl=Z.sort, method.pre='one.step', z.pre=Z.pre.one)
bt.KO.one.step = result_KO_one_step$bt.hat
opt.w.one.step = result_KO_one_step$opt.w 
P.KO.one.step = result_KO_one_step$P.KO
robust.inf.KO.one.step = cluster.se(id=firm.id, bt=bt.KO.one.step, x=x.reg, y=y, P=P.KO.one.step)
se.KO.one.step = robust.inf.KO.one.step$se
cat('se ap = ', se.KO.one.step[1],'\n')

# KO: lasso
result_KO_lasso = tsls_KO( y=y, x.end=price, x.exo=X, z.excl=Z.sort, method.pre='lasso', z.pre=NULL)
bt.KO.lasso = result_KO_lasso$bt.hat
opt.w.lasso = result_KO_lasso$opt.w
P.KO.lasso = result_KO_lasso$P.KO
robust.inf.KO.lasso = cluster.se(id=firm.id, bt=bt.KO.lasso, x=x.reg, y=y, P=P.KO.lasso)
se.KO.lasso = robust.inf.KO.lasso$se
cat('se ap = ', se.KO.lasso[1],'\n')


# CSA: Complete Subset Average Estimation
cat('--------------------------------- \n')
cat('CSA \n')
cat('--------------------------------- \n')
# CSA: one-step
result_CSA_one_step = tsls_CSA_aveP_in( y=y, x.end=price, x.exo=X, z.excl=Z.sort, ld=NULL, aveP=aveP, method.pre='one.step', z.pre=Z.pre.one)
bt.CSA.one.step = result_CSA_one_step$bt.hat
K.CSA.one.step = result_CSA_one_step$opt.k
P.CSA.one.step = aveP[,,K.CSA.one.step]
robust.inf.CSA.one.step = cluster.se(id=firm.id, bt=bt.CSA.one.step, x=x.reg, y=y, P=P.CSA.one.step)
se.CSA.one.step = robust.inf.CSA.one.step$se
cat('se ap = ', se.CSA.one.step[1],'\n')

# CSA: lasso
result_CSA_lasso = tsls_CSA_aveP_in( y=y, x.end=price, x.exo=X, z.excl=Z.sort, ld=NULL, aveP=aveP, method.pre='lasso', z.pre=NULL)
bt.CSA.lasso = result_CSA_lasso$bt.hat
K.CSA.lasso = result_CSA_lasso$opt.k
P.CSA.lasso = aveP[,,K.CSA.lasso]
robust.inf.CSA.lasso = cluster.se(id=firm.id, bt=bt.CSA.lasso, x=x.reg, y=y, P=P.CSA.lasso)
se.CSA.lasso = robust.inf.CSA.lasso$se
cat('se ap = ', se.CSA.lasso[1],'\n')


# Calculate the number of inelastic products
BLP_data=readMat('BLP_data.mat')
o_price = BLP_data$price
o_share = BLP_data$share

names.method=c('OLS', '2SLS', 'DN.one.step', 'DN.lasso', 'KO.one.step', 'KO.lasso', 'CSA.one.step', 'CSA.lasso')

# Construct the Table
ap.hat = c(coef(m.ols)['price'], coef(m.tsls)[1,1], bt.DN.one.step[1,1], bt.DN.lasso[1,1], bt.KO.one.step[1,1], bt.KO.lasso[1,1], 
           bt.CSA.one.step[1,1], bt.CSA.lasso[1,1])
No.Inelastic = sapply(ap.hat, inelastic, price=o_price, share=o_share)
K.opt = c(NA, NA, K.DN.one.step, K.DN.lasso, NA, NA
              , K.CSA.one.step, K.CSA.lasso)
robust.se = round(c(se.OLS[1], se.2SLS[1], se.DN.one.step[1], se.DN.lasso[1], se.KO.one.step[1], se.KO.lasso[1], se.CSA.one.step[1], se.CSA.lasso[1]),3 )


table.main = cbind(ap.hat, robust.se, K.opt, No.Inelastic)
row.names(table.main) = names.method

print(table.main)

end.time=proc.time()[3]
exec.time = end.time-start.time
cat('Elapsed Time: ',exec.time,'\n')