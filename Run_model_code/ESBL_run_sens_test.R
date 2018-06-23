########################################
# SENSITIVITY TEST PERFORMANCE
########################################
rm(list=ls())

library(deSolve)

setwd("./Model_code") # Change to appropriate directory

system("R CMD SHLIB ESBL.c")
dyn.load("ESBL.so")  #use this one for unix


#  Create matrices for sensitivity performance
spec = c(0,0 + 0.02*c(1:50))
sens = c(0,0 + 0.02*c(1:50))

low = matrix(NA,length(spec), length(sens)) # spec = rows, sens = columns
colnames(low) = sens
rownames(low) = spec
low = list(prop_o = low, prop_u = low, prev_s = low, prev_a = low, prev_b= low,inc_s = low, inc_a = low, inc_b = low,
           prevI_s = low,prevI_a = low,prevI_b = low, fIs = low,fIa = low,fIb = low)
data_test = list(RDT_l=low,RDT_m=low,RDT_h=low)

p_A = c(0 + 0.025*c(0:40))
p_B = as.numeric(c(0 + 0.025*c(0:40)))

scen = 3# chance to right importation scene
switch = 2 # switch = 1; no switch = 2
sel = 2 # s = 1: selection 10^-3; s=2: no selection

name = paste0((ifelse(switch == 1 & sel == 1, "f_cost_",
                      ifelse(switch == 2 & sel == 1, "f_cost_noswitch_",
                             ifelse(switch == 1 & sel == 2, "f_cost_nosel_", "f_cost_noswitch_nosel_")))),"scen",scen)

#name = paste0((ifelse(switch == 1, "f_cost_","f_cost_noswitch_")),"scen",s)

scens_A = c(0.1, 0.25, 0.4)
scens_B = c(0, 0.025, 0.08)
scens_ts = c(3,500) # scenario no switch --> ts = 500
scens_s = c(10^-3, 0)

scenA = scens_A[scen] # Importation scenario 1: pA = 0.1; scenario 2: pA = 0.25; scenario 3: pA = 0.4
scenB = scens_B[scen] # Importation scenario 1: pB = 0; scenario 2: pB = 0.025; scenario 3: pB = 0.1
scenAn = which(p_A == scenA)
scenBn = which(p_B == scenB)
ts = scens_ts[switch]
s = scens_s[sel]
  
  
parameters1<-c(
  beta = 0.056, 
  d = 1/9, 
  d_i =1/9*(1-0.14), 
  a_s = 1-(0.05+0.015),
  a_i = 0.015,
  a_c = 0.05, 
  p_s = 1-(scenA+scenB),
  p_A = scenA,
  p_B = scenB,  
  c1 = 0.95, # scenario of fitness cost: 0.95; see ref Sandegren et al. JAC 2012, 67(1), 74-83
  c2 = 0.8, # scenario of fitness cost: 0.8; see Adler et al. JAC 2013, 68(1), 51-59; Tsai et al 2011 JAC; Garcia-Sureda et al 2011 AAC. The first ref is e. coli, but only ref that quantifies fitness cost
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
  r_delay = 1/(6+ts), 
  cycl_period = 42, 
  cycling =0, 
  rdt = 1,
  ts = ts,
  mixing=0)

years = 5
trans = c("m","h")
for(i in unique(c(0.144,0.222))){
    p = parameters1 
    p["beta"] = i
    assign(value=p,x= paste0("parameters1_",trans[which(unique(c(0.144,0.222))==i)]))
}

par_list = list(parameters1,parameters1_m,parameters1_h)


init.state<-c(
  S=28,
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


for(d in 1:length(data_test)){
  parameters1=unlist(par_list[d])
  print(paste("i=",i))
  for(i in 1:length(sens)){
    parameters1["f_iaA"] = 1-sens[i]
    parameters1["f_iaB"] = sens[i]
    parameters1["f_ib_last"] = sens[i]
      for(b in 1:length(spec)){
        parameters1["f_isA"] = spec[b]
        parameters1["f_isB"] = 1-spec[b]
    # print(paste("i=",i))
        run1<-ode(y=init.state,times=seq(1,(years*365),by=1), func="derivs",parms= parameters1,dllname="ESBL", initfunc="initmod",nout=15,
                  outnames=c("Sum","prop_o","prop_u","prev_s","prev_a","prev_b","fCs","fCa","fCb","prevI_s","prevI_a","prevI_b","fIs","fIa","fIb"),method=rkMethod(method = 'rk4',hini=0.001))
        print(paste("b=",b))
    run_y = sapply(data.frame(run1[(365*(years-1)):(365*years),]),FUN=mean)
    admis_y = run1[(365*years),"admis"]-run1[(365*(years-1)),"admis"]
    data_test[[d]]$prop_o[b,i] = run_y["prop_o"]
    data_test[[d]]$prop_u[b,i] = run_y["prop_u"]
    data_test[[d]]$prev_s[b,i] = run_y["prev_s"]
    data_test[[d]]$prev_a[b,i] = run_y["prev_a"]
    data_test[[d]]$prev_b[b,i] = run_y["prev_b"]
    data_test[[d]]$inc_s[b,i] = (run1[(365*years),"Inc_s"]-run1[(365*(years-1)),"Inc_s"])/admis_y*1000
    data_test[[d]]$inc_a[b,i] = (run1[(365*years),"Inc_a"]-run1[(365*(years-1)),"Inc_a"])/admis_y*1000
    data_test[[d]]$inc_b[b,i] = (run1[(365*years),"Inc_b"]-run1[(365*(years-1)),"Inc_b"])/admis_y*1000
    data_test[[d]]$prevI_s[b,i] = run_y["prevI_s"]
    data_test[[d]]$prevI_a[b,i] = run_y["prevI_a"]
    data_test[[d]]$prevI_b[b,i] = run_y["prevI_b"]
    data_test[[d]]$fIs[b,i] = run_y["fIs"]
    data_test[[d]]$fIa[b,i] = run_y["fIa"]
    data_test[[d]]$fIb[b,i] = run_y["fIb"]
      }
    }
  }


setwd("~/Dropbox/RGNOSIS/Methods/Empiric prescribing model/Model/Output_C/")
save(file=paste0("sens_test_", name, ".RData"),data_test)
 
