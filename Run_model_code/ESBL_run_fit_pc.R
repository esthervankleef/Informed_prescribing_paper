rm(list=ls())

library(deSolve)
library(manipulate)
setwd("./Model_code")

system("R CMD SHLIB ESBL.c")
dyn.load("ESBL.so")  #use this one for unix

init.state<-c(
  S=58,
  Cs=1,
  Ca=0,
  Cb=0,
  Is_a=1,
  Is_o=0,
  Ia_u=0,
  Ia_a=0,
  Ib_u=0,
  Ib_a=0,
  Inc_s=0,
  Inc_a=0,
  Inc_b=0,
  admis = 0)



#  Do 1000 random draws of p_c and store in a results matrix
#  Low transmission
set.seed(100)
results<-matrix(NA,nrow=200,ncol=365*2+1)
for(i in 1:200){
  parameters1<-c(
    beta=0.056,
    d = 1/9, 
    d_i =1/9*(1-0.14), 
    a_s = 1-(0.05+0.015),
    a_i = 0.015,
    a_c = 0.05, 
    p_s = 1-(0.4+0.08),
    p_A = 0.4,
    p_B = 0.08,  
    c1 = 1, 
    c2 = 1,
    s_1 = 10^-3,
    s_2 = 10^-3, 
    f_cA = 0.18, 
    f_cB = 0.17, 
    f_clast = 0.038, 
    f_isA = 0.82, 
    f_isB = 1-0.82,
    f_iaA = 0.56, 
    f_iaB = 1-0.56, 
    f_ib_last = 0, 
    p_c = runif(1,0,1), 
    r = 1/6, 
    r_delay = 1/(3+6), 
    cycl_period = 42*4, 
    cycling = 0, 
    rdt = 0,
    ts = 3,
    mixing=0)
  
  run1<-ode(y=init.state,times=seq(1,(2*365),by=1), func="derivs",parms= parameters1,dllname="ESBL", initfunc="initmod",nout=15,
            outnames=c("Sum","prop_o","prop_u","prev_s","prev_a","prev_b","fCs","fCa","fCb","prevI_s","prevI_a","prevI_b","fIs","fIa","fIb"),method=rkMethod(method = 'rk4',hini=0.001))
  results[i,]<-c(rowSums(run1[,c(6:11)])/60,parameters1[22])
} 
# 
# #  Med transmission
# results2<-matrix(NA,nrow=100,ncol=365*2+1)
# for(i in 1:100){
#   parameters1<-c(
#     beta=0.1667,
#     d = 1/9, 
#     d_i =1/9*(1-0.14), 
#     a_s = 1-(0.05+0.015),
#     a_i = 0.015,
#     a_c = 0.05, 
#     p_s = 1-(0.2+0.08),
#     p_A = 0.2,
#     p_B = 0.08,  
#     c1 = 1, 
#     c2 = 1,
#     s_1 = 10^-3,
#     s_2 = 10^-3, 
#     f_cA = 0.5, 
#     f_cB = 0.5, 
#     f_clast = 0.03, 
#     f_isA = 0.82, 
#     f_isB = 1-0.82,
#     f_iaA = 0.56, 
#     f_iaB = 1-0.56, 
#     f_ib_last = 0, 
#     p_c = runif(1,0,1), 
#     r = 1/6, 
#     r_delay = 1/(3+6), 
#     cycl_period = 42*4, 
#     cycling = 0, 
#     rdt = 0,
#     ts = 3,
#     mixing=0)
#   
#   run2<-ode(y=init.state,times=seq(1,(2*365),by=1), func="derivs",parms=parameters1,dllname="ESBL", initfunc="initmod",nout=1,outnames="Sum",method=rkMethod(method = 'rk4',hini=0.001))
#   results2[i,]<-c(rowSums(run2[,c(6:11)])/60,parameters1[22])
# }

#  High transmission
results3<-matrix(NA,nrow=200,ncol=365*2+1)
for(i in 1:200){
  parameters1<-c(
    beta=0.222,
    d = 1/9, 
    d_i =1/9*(1-0.14), 
    a_s = 1-(0.05+0.015),
    a_i = 0.015,
    a_c = 0.05, 
    p_s = 1-(0.4+0.08),
    p_A = 0.4,
    p_B = 0.08,  
    c1 = 1, 
    c2 = 1,
    s_1 = 10^-3,
    s_2 = 10^-3, 
    f_cA = 0.17, 
    f_cB = 0.18, 
    f_clast = 0.038, 
    f_isA = 0.82, 
    f_isB = 1-0.82,
    f_iaA = 0.56, 
    f_iaB = 1-0.56, 
    f_ib_last = 0, 
    p_c = runif(1,0,1), 
    r = 1/6, 
    r_delay = 1/(3+6), 
    cycl_period = 42*4, 
    cycling = 0, 
    rdt = 0,
    ts = 3,
    mixing=0)
  
  run3<-ode(y=init.state,times=seq(1,(2*365),by=1), func="derivs",parms= parameters1,dllname="ESBL", initfunc="initmod",nout=15,
                  outnames=c("Sum","prop_o","prop_u","prev_s","prev_a","prev_b","fCs","fCa","fCb","prevI_s","prevI_a","prevI_b","fIs","fIa","fIb"),method=rkMethod(method = 'rk4',hini=0.001))
  
  results3[i,]<-c(rowSums(run3[,c(6:11)])/60,parameters1[22])
}

  
  pdf("~/Dropbox/RGNOSIS/Methods/Empiric prescribing model/Model/Output_C/Figures/Paper/Progression_rate.pdf",
      width=15,height=5)
  par(mfrow=c(1,3))
  
  # cols<-rainbow(100)
  # plot(results2[1,1:(2*365)],type='l',ylim=c(0,0.4), ylab="Fraction infected",main="Rn ~ 1.5")
  # for(i in 2:100) lines(results2[i,1:(2*365)],col=cols[i])
  # lines(c(1:(2*365)),rep(0.07,(2*365)),lty=2)
  
  cols<-rainbow(200)
  plot(results3[1,1:(2*365)],type='l',ylim=c(0,0.4), ylab="Fraction infected",main="Rn ~ 2")
  for(i in 2:200) lines(results3[i,1:(2*365)],col=cols[i])
  lines(c(1:(2*365)),rep(0.1,(2*365)),lty=2)
  
  # best_fit_m = which(abs(results2[,365*2]-0.07)==min(abs(results2[,365*2]-0.07)))
  # best_fit_m
  # plot(results2[best_fit_m,1:(2*365)],type='l',ylim=c(0,0.4), ylab="Fraction infected",main=paste("Rn ~ 1.5, best fit p_c = ",round(results2[best_fit_m,365*2+1],2)))
  # lines(c(1:(2*365)),rep(0.07,(2*365)),lty=2)
  
  best_fit_h = which(abs(results3[,2*365]-0.1)==min(abs(results3[,2*365]-0.1)))
  best_fit_h
  plot(results3[best_fit_h,1:(2*365)],type='l',ylim=c(0,0.4),xlab="Time", ylab="Fraction infected",main=paste("Rn ~ 2, best fit p_c = ",round(results3[best_fit_h,365*2+1],2)))
  lines(c(1:(2*365)),rep(0.1,(2*365)),lty=2)
  
  low = which(abs(round(results[,(365*2)+1],3)-0.060) == min(abs(round(results[,(365*2)+1],3)-0.060)))
  results[low,(365*2)+1]
  
  plot(results[low,1:(2*365)],type='l',ylim=c(0,0.4), xlab="Time", ylab="Fraction infected",main=paste("Rn ~ 0.5, if p_c =",round(results[low,365*2+1],2)))
  lines(c(1:(2*365)),rep(0.035,(2*365)),lty=2)
  
  # cols<-rainbow(200)
  # plot(results[1,1:(2*365)],type='l',ylim=c(0,0.4), ylab="Fraction infected",main="Rn ~ 0.5")
  # for(i in 2:200) lines(results[i,1:(2*365)],col=cols[i])
  # lines(c(1:(2*365)),rep(0.035,(2*365)),lty=2)
  # 
  # best_fit = which(abs(results[,2*365]-0.035)==min(abs(results[,2*365]-0.035)))
  # best_fit
  # plot(results[best_fit,1:(2*365)],type='l',ylim=c(0,0.4), ylab="Fraction infected",main=paste("Rn ~ 0.5, best fit p_c = ",round(results[best_fit,365*2+1],2)))
  # lines(c(1:(2*365)),rep(0.035,(2*365)),lty=2)
  # 
  
  dev.off()
