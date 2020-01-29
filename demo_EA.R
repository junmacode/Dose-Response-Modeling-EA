rm(list=ls())
#------ Begin Demo ------#
source("ea_curve_fitting.r")

## Data from hill-slope/log-logistic model##
load("drugExp.Rdata")


xvals<-c(0.04193, 0.01939,0.011437, 0.007041, 0.004522, 0.001693)   ### carboA dose
dat<-drugExp$CarboA[,3:8]
dat<-dat[60,]               #carboA response


# EA parameters are initial population size, stable/equilibrium populaton size, number of tournaments,
# tournament size, and number of generations

ea.params = c(1000,200,20,10,100)



# 5 parameter log logistic model

fit1 = eadrm( dat, xvals, model="h5",ea.params )
params<-as.numeric(fit1$params)
x1<-seq(0, 1, by = 0.00001)
x2<-params[2]+(params[1]-params[2])/((1+exp(params[4]*(log(x1)-log(params[3]))))^params[5])
plot(xvals,dat,xlab="Dose",ylab="Response")
lines(x1,x2)



# 4 parameter log logistic model

fit2 = eadrm( dat, xvals,  model="h4",ea.params )
m1_value<-as.numeric(fit2$params)
t1<-seq(0, 1, by = 0.00001)
t2<-m1_value[1] - (m1_value[1]-m1_value[2])/((1+(t1/m1_value[3])^(m1_value[4])))
plot(xvals,dat,xlab="Dose",ylab="Response")
lines(t1,t2)



# 3 parameter log logistic model

fit3 = eadrm( dat, xvals,  model="h3",ea.params )
params<-as.numeric(fit3$params)
x1<-seq(0, 1, by = 0.00001)
x2<-(params[1])/(1+exp(params[3]*(log(x1)-log(params[2]))))
plot(xvals,dat,xlab="Dose",ylab="Response")
lines(x1,x2)


# exponential model
fit4 = eadrm( dat, xvals,  model="e",ea.params )
params<-as.numeric(fit4$params)
x1<-seq(0, 1, by = 0.00001)
x2<-params[1]*exp(params[2]*x1)
plot(xvals,dat,xlab="Dose",ylab="Response")
lines(x1,x2)



# no model assumption

fit5 = eadrm( dat, xvals,  model="all",ea.params )
params<-as.numeric(fit5$params)
n1<-length(params)
t1<-seq(0, 1, by = 0.00001)

if(n1==5) t2<-params[2]+(params[1]-params[2])/((1+exp(params[4]*(log(x1)-log(params[3]))))^params[5])
if(n1==4) t2<-params[1] - (params[1]-params[2])/((1+(t1/params[3])^(params[4])))
if(n1==3) t2<-(params[1])/(1+exp(params[3]*(log(t1)-log(params[2]))))
if(n1==2) t2<-params[1]*exp(params[2]*t1)

plot(xvals,dat,xlab="Dose",ylab="Response")
lines(t1,t2)




