rm(list=ls())
library(deSolve)
library(manipulate)

setwd("./Model_code")

system("R CMD SHLIB ESBL.c")
dyn.load("ESBL.so")  #use this one for unix


#######################
# SENSITIVITY ANALYIS
######################
scen = 3 # chance to right importation scene
switch = 1 # switch = 1; no switch = 2
sel = 1 # s = 1: selection 10^-3; s=2: no selection

name = paste0((ifelse(switch == 1 & sel == 1, "f_cost_",
                      ifelse(switch == 2 & sel == 1, "f_cost_noswitch_",
                             ifelse(switch == 1 & sel == 2, "f_cost_nosel_", "f_cost_noswitch_nosel_")))),"scen",scen)

scens_A = c(0.1, 0.25, 0.4)
scens_B = c(0, 0.025, 0.1)
scens_ts = c(3,500) # scenario no switch --> ts = 500
scens_s = c(10^-3, 10^-100)

scenA = scens_A[scen] # Importation scenario 1: pA = 0.1; scenario 2: pA = 0.25; scenario 3: pA = 0.4
scenB = scens_B[scen] # Importation scenario 1: pB = 0; scenario 2: pB = 0.025; scenario 3: pB = 0.1
ts = scens_ts[switch]
s = scens_s[sel]

theta=c(beta=0.056, 
        d=1/9, 
        d_i=1/9*(1-0.14),
        a_c = 0.05, 
        a_i= 0.015,
        p_A= scenA, 
        p_B=scenB,
        c1=0.95,
        c2=0.8, 
        s_1=s,
        s_2=s,
        f_cA=0.18,
        f_cB=0.17,
        f_clast = 0.038, 
        f_isA=0.82,
        f_iaA=0.44,
        p_c=0.06,
        r=1/6, 
        r_delay=1/(ts+6), 
        cycl_period=42,
        cycling=0,
        rdt=0,
        ts=ts,
        mixing=0,
        a_s=1-(0.05+0.015),
        p_s=1-(scenA+scenB),
        f_isB=1-0.82,
        f_iaB=1-0.44, 
        f_ib_last=1-0.44)

years = 5 

theta_all = matrix(theta, nrow=29,ncol=100)
theta_all = data.frame(t(theta_all))
names(theta_all) = names(theta)

beta = c(0.023,0.056,0.144,0.222,0.28)
d = c(1/4,1/9,1/12,1/16,1/20)
d_i =c(0.076,0.085,0.096,0.111,0.123)
a_c = c(0, 0.05,0.1,0.15, 0.2)
a_i =c(0,0.015,0.03, 0.045,0.06)
p_A = c(0,0.2,0.4,0.6,0.8)
p_B = c(0,0.2,0.4,0.6,0.8)
c1 = c(0.2,0.4,0.6,0.8,1)
c2 = c(0.2,0.4,0.6,0.8,1)
s_1 =c(0, 0.001,0.01,0.1,0.2)
s_2 =c(0, 0.001,0.01,0.1,0.2)
f_cA =c(0,0.1,0.2,0.3,0.4)
f_cB =c(0,0.1,0.2,0.3,0.4)
f_clast = c(0,0.1,0.2,0.3,0.4)
f_isA =c(0,0.25,0.5,0.75,1)
f_iaA =c(0,0.25,0.5,0.75,1)
p_c =c(0.01, 0.03, 0.06,0.09,0.12)
r = c(1/3,1/6,1/9,1/12,1/15)
r_delay =c(1/(0+6),1/(3+6),1/(10+6),1/(50+6),1/(500+6))



var_value = c(beta,d,d_i,a_c,a_i,p_A,p_B,c1,c2,s_1,s_2,f_cA,f_cB,f_clast,f_isA,f_iaA,p_c,r,r_delay) #m,m_i,

params = c("beta","d","d_i","a_c","a_i","p_A","p_B","c1","c2","s_1","s_2",
           "f_cA","f_cB","f_clast","f_isA","f_iaA","p_c","r","r_delay") 

for(i in unique(params)){
  print(i)
  print(which(names(theta_all)==i)-1)
  theta_all[c(1,2,3,4,5)+5*(which(names(theta_all)==i)-1),which(names(theta_all)==i)] = var_value[c(1,2,3,4,5)+5*(which(names(theta_all)==i)-1)]
  print(theta_all[c(1,2,3,4,5)+5*(which(names(theta_all)==i)-1),which(names(theta_all)==i)])
  
}

theta_all$a_s = 1-(theta_all$a_c+theta_all$a_i)
theta_all$p_A = ifelse(theta_all$p_B==0.8,0.2,theta_all$p_A)
theta_all$p_s = 1-(theta_all$p_A+theta_all$p_B)
theta_all$f_isB = 1-theta_all$f_isA
theta_all$f_iaB = 1-theta_all$f_iaA
theta_all$f_ib_last = 1-theta_all$f_iaA

times = 5
var_change = c(rep("beta",times), rep("d",times), rep("d_i",times),rep("a_c",times), #rep("m",times), rep("m_i",times),
               rep("a_i",times), rep("p_A",times), rep("p_B",times), rep("c1",times), rep("c2",times),rep("s_1",times),rep("s_2",times),
               rep("f_cA",times), rep("f_cB", times), rep("f_clast",times),rep("f_is_A",times),rep("f_ia_A",times),rep("p_c",times),
               rep("r", times), rep("r_delay",times),rep("cycl_period",5))

sens = data.frame(cbind(var_change,var_value=c(var_value,30,42,90,180,365)))

low = matrix(NA,length(rownames(theta_all)), 6) 
colnames(low) = c("Standard", "Informed", "Cycling", "Mixing", "var_change","var_value")
rownames(low) = c(1:100)
low[,"var_change"] = as.character(sens$var_change)
low[,"var_value"] = rep(c(1:5),20)

low = list(prop_o = low, prop_u = low, prev_s = low, prev_a = low, prev_b= low, inc_s = low, inc_a = low, inc_b = low, domI=low,domP=low,
             prevI_s = low,prevI_a = low,prevI_b = low, fIs = low,fIa = low,fIb = low)
data_s = list(l=low,m=low,h=low)

theta_all = theta_all[,c(1,2,3,25,5,4,26,6,7,8,9,10,11,12,13,14,15,27,16,28,29,17,18,19,20,
                         21,22,23,24)]
parameters_c = theta_all
parameters_r = theta_all 
parameters_m = theta_all 

parameters_c["cycling"] = 1;parameters_c["f_ib_last"] = 0;
parameters_c$cycl_period = c(rep(42,95),30,42,90,180,365)
parameters_r["rdt"] = 1; parameters_r["f_isA"] = 0.85;parameters_r["f_isB"] =1-parameters_r["f_isA"];
parameters_r["f_iaA"] = 0.05;parameters_r["f_iaB"] = 1-parameters_r["f_iaA"];parameters_r["f_ib_last"] = 0.95;
parameters_m["mixing"] = 1; parameters_m["f_ib_last"] = 0

par_list = list(theta_all,parameters_r,parameters_c,parameters_m)
par_list = c(par_list,par_list,par_list)

for(d in 5:8){
  par_list[[d]]["beta"] =c(0.056,0.144,0.17,0.222,0.28,rep(0.144,95))
  print(par_list[[d]]["beta"])
}

for(d in 9:12){
  par_list[[d]]["beta"] = c(0.056,0.144,0.17,0.222,0.28,rep(0.222,95))
  print(par_list[[d]]["beta"])
}

# RUN SENSITIVITY ANALYSIS
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


for(d in 1:12){
  par_set=data.frame(par_list[d])
  print(paste("d =",d))
  for(i in 1:length(sens$var_change)){
    print(i)
    parameters1=unlist(par_set[i,])
    print(parameters1)
    run1<-ode(y=init.state,times=seq(1,(15*365),by=1), func="derivs",parms= parameters1,dllname="ESBL", initfunc="initmod",nout=15,
              outnames=c("Sum","prop_o","prop_u","prev_s","prev_a","prev_b","fCs","fCa","fCb","prevI_s","prevI_a","prevI_b","fIs","fIa","fIb"),method=rkMethod(method = 'rk4',hini=0.001))
    
    cols = ifelse(d%in%c(1,5,9),1,ifelse(d%in%c(2,6,10),2, ifelse(d%in%c(3,7,11),3,4)))
    rows = i
    trans = ifelse(d%in%c(1:4),1,ifelse(d%in%c(5:8),2,3))
    run_y = round(sapply(data.frame(run1[(365*(years-1)):(365*years),]),FUN=mean),7)
    admis_y = run1[(365*years),"admis"]-run1[(365*(years-1)),"admis"]
    data_s[[trans]]$prop_o[rows,cols] = run_y["prop_o"]
    data_s[[trans]]$prop_u[rows,cols] = run_y["prop_u"]
    data_s[[trans]]$prev_s[rows,cols] = run_y["prev_s"]
    data_s[[trans]]$prev_a[rows,cols] = run_y["prev_a"]
    data_s[[trans]]$prev_b[rows,cols] = run_y["prev_b"]
    data_s[[trans]]$inc_s[rows,cols] = (run1[(365*years),"Inc_s"]-run1[(365*(years-1)),"Inc_s"])/admis_y*1000
    data_s[[trans]]$inc_a[rows,cols] = (run1[(365*years),"Inc_a"]-run1[(365*(years-1)),"Inc_a"])/admis_y*1000
    data_s[[trans]]$inc_b[rows,cols] = (run1[(365*years),"Inc_b"]-run1[(365*(years-1)),"Inc_b"])/admis_y*1000
    data_s[[trans]]$prevI_s[rows,cols] = run_y["prevI_s"]
    data_s[[trans]]$prevI_a[rows,cols] = run_y["prevI_a"]
    data_s[[trans]]$prevI_b[rows,cols] = run_y["prevI_b"]
    data_s[[trans]]$fIs[rows,cols] = run_y["fIs"]
    data_s[[trans]]$fIa[rows,cols] = run_y["fIa"]
    data_s[[trans]]$fIb[rows,cols] = run_y["fIb"]
  }
}

# WHICH STRATEGY DOMINANT ###
for(d in 1:length(data_s)){
  dom_o=NULL
  dom_u=NULL
  dom_fIa=NULL
  dom_fIb=NULL
  #paste("d=",print(d))
  for(i in 1:length(sens$var_change)){
  #print(paste("i=",i))
  o = which(data_s[[d]]$prop_o[i,c(1:4)]==min(data_s[[d]]$prop_o[i,c(1:4)]))
  o = ifelse(length(o)==4,0,o)
  u = which(data_s[[d]]$prop_u[i,c(1:4)]==min(data_s[[d]]$prop_u[i,c(1:4)]))
  u = ifelse(length(u)==4,0,u)
  fa = which(data_s[[d]]$fIa[i,c(1:4)]==min(data_s[[d]]$fIa[i,c(1:4)]))
  fa = ifelse(length(fa)==4,0,fa)
  fb = which(data_s[[d]]$fIb[i,c(1:4)]==min(data_s[[d]]$fIb[i,c(1:4)]))
  fb = ifelse(length(fb)==4,0,fb)
  dom_o = c(dom_o,o)
  dom_u = c(dom_u,u)
  dom_fIa = c(dom_fIa,fa)
  dom_fIb = c(dom_fIb,fb)
  }
  data_s[[d]][["dom"]] = cbind(dom_o,dom_u,dom_fIa, dom_fIb, data_s[[d]]$fIb[,5],data_s[[d]]$fIb[,6]) 
}

#######################################
# Run model under baseline parameters

theta_l<-c(
  beta = 0.056, 
  d = 1/9, 
  d_i =1/9*(1-0.14), 
  a_s = 1-(0.05+0.015),
  a_i = 0.015,
  a_c = 0.05, 
  p_s = 1-(scenA+scenB),
  p_A = scenA,
  p_B = scenB,  
  c1 = 0.95, 
  c2 = 0.8,
  s_1 = s,
  s_2 = s, 
  f_cA = 0.18, 
  f_cB = 0.17, 
  f_clast = 0.038, 
  f_isA = 0.82, 
  f_isB = 1-0.82,
  f_iaA = 0.44, 
  f_iaB = 1-0.44, 
  f_ib_last = 1-0.44, 
  p_c = 0.06, 
  r = 1/6, 
  r_delay = 1/(ts+6), 
  cycl_period = 42, 
  cycling =0, 
  rdt = 0,
  ts = ts,
  mixing=0)


parameters_c = theta_l
parameters_r = theta_l 
parameters_m = theta_l 

parameters_c["cycling"] = 1;parameters_c["f_ib_last"] = 0;
parameters_r["rdt"] = 1; parameters_r["f_isA"] = 0.85;parameters_r["f_isB"] =1-parameters_r["f_isA"];
parameters_r["f_iaA"] = 0.05;parameters_r["f_iaB"] = 1-parameters_r["f_iaA"];parameters_r["f_ib_last"] = 0.95;
parameters_m["mixing"] = 1; parameters_m["f_ib_last"] = 0

par_list = list(theta_l, parameters_r, parameters_c, parameters_m)
par_list = c(par_list,par_list,par_list)

for(d in 5:8){
  par_list[[d]]["beta"] =0.144
  print(par_list[[d]]["beta"])
}

for(d in 9:12){
  par_list[[d]]["beta"] = 0.222
  print(par_list[[d]]["beta"])
}


low = matrix(NA,1, 4) 
colnames(low) = c("Standard care", "Informed prescribing", "Cycling", "Mixing")

mean=list(prop_o = low, prop_u = low, prev_s = low, prev_a = low, prev_b= low, inc_s = low, inc_a = low, inc_b = low, domI=low,domP=low,
     prevI_s = low,prevI_a = low,prevI_b = low, fIs = low,fIa = low,fIb = low)

data_mean = list(l=mean,m=mean,h=mean)

for(d in 1:12){
  parameters1=unlist(par_list[d])
  print(paste("d =",d))
    print(parameters1)
    run1<-ode(y=init.state,times=seq(1,(years*365),by=1), func="derivs",parms= parameters1,dllname="ESBL", initfunc="initmod",nout=15,
              outnames=c("Sum","prop_o","prop_u","prev_s","prev_a","prev_b","fCs","fCa","fCb","prevI_s","prevI_a","prevI_b","fIs","fIa","fIb"),method=rkMethod(method = 'rk4',hini=0.001))
    cols = ifelse(d%in%c(1,5,9),1,ifelse(d%in%c(2,6,10),2, ifelse(d%in%c(3,7,11),3,4)))
    rows = 1
    trans = ifelse(d%in%c(1:4),1,ifelse(d%in%c(5:8),2,3))
    run_y = round(sapply(data.frame(run1[(365*(years-1)):(365*years),]),FUN=mean),7)
    admis_y = run1[(365*years),"admis"]-run1[(365*(years-1)),"admis"]
    data_mean[[trans]]$prop_o[rows,cols] = run_y["prop_o"]
    data_mean[[trans]]$prop_u[rows,cols] = run_y["prop_u"]
    data_mean[[trans]]$prev_s[rows,cols] = run_y["prev_s"]
    data_mean[[trans]]$prev_a[rows,cols] = run_y["prev_a"]
    data_mean[[trans]]$prev_b[rows,cols] = run_y["prev_b"]
    data_mean[[trans]]$inc_s[rows,cols] = (run1[(365*years),"Inc_s"]-run1[(365*(years-1)),"Inc_s"])/admis_y*1000
    data_mean[[trans]]$inc_a[rows,cols] = (run1[(365*years),"Inc_a"]-run1[(365*(years-1)),"Inc_a"])/admis_y*1000
    data_mean[[trans]]$inc_b[rows,cols] = (run1[(365*years),"Inc_b"]-run1[(365*(years-1)),"Inc_b"])/admis_y*1000
    data_mean[[trans]]$prevI_s[rows,cols] = run_y["prevI_s"]
    data_mean[[trans]]$prevI_a[rows,cols] = run_y["prevI_a"]
    data_mean[[trans]]$prevI_b[rows,cols] = run_y["prevI_b"]
    data_mean[[trans]]$fIs[rows,cols] = run_y["fIs"]
    data_mean[[trans]]$fIa[rows,cols] = run_y["fIa"]
    data_mean[[trans]]$fIb[rows,cols] = run_y["fIb"]
}


#


setwd("~/Dropbox/RGNOSIS/Methods/Empiric prescribing model/Model/Output_C/")
save(file=paste0("sens_all",name,".RData"),data_s)
save(file=paste0("mean_all",name,".RData"),data_mean)
