########################################
# SENSITIVITY COMMUNITY IMPORTATION
########################################
rm(list=ls())

library(deSolve)
library(manipulate)

setwd("./Model_code") # Change to appropriate directory

system("R CMD SHLIB ESBL.c")
dyn.load("ESBL.so")  #use this one for unix

imp = read.csv("~/Scenarios_parameters/cont_Cycl_low.csv") # Change this to appropriate scene 

imp_pA = imp$p_A
imp_pB = imp$p_B

#  Create matrices for sensitivity importation
p_A = c(0 + 0.025*c(0:40))
p_B = as.numeric(c(0 + 0.025*c(0:40)))

low = matrix(NA,length(p_A), length(p_A)) # p_A = rows, p_B = columns
colnames(low) = p_B
rownames(low) = p_A
low = list(prop_o = low, prop_u = low, prev_s = low, prev_a = low, prev_b= low, inc_s = low, inc_a = low, inc_b = low, domI=low,domP=low,
           prevI_s = low,prevI_a = low,prevI_b = low, fIs = low,fIa = low,fIb = low, prop_ia = low)
data_com = list(base_l=low,RDT_l=low,Cycling_l=low,Mixing_l=low,
            base_m=low,RDT_m=low,Cycling_m=low,Mixing_m=low,
            base_h=low,RDT_h=low,Cycling_h=low,Mixing_h=low)

# Importation scenario 1: pA = 0.1; scenario 2: pA = 0.25; scenario 3: pA = 0.4
# Importation scenario 1: pB = 0; scenario 2: pB = 0.025; scenario 3: pB = 0.1

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
  rdt = 0,
  ts = ts,
  mixing=0)

years = 5

parameters_c = parameters1 
parameters_r = parameters1 
parameters_m = parameters1 

parameters_c["cycling"] = 1;parameters_c["cycl_period"] = 42; parameters_c["f_ib_last"] = 0
parameters_r["rdt"] = 1; parameters_r["f_isA"] = 0.85;parameters_r["f_isB"] =1-parameters_r["f_isA"] ;
parameters_r["f_iaA"] = 0.05;parameters_r["f_iaB"] = 1-parameters_r["f_iaA"];parameters_r["f_ib_last"] = 0.95;
parameters_m["mixing"] = 1; parameters_m["f_ib_last"] = 0

par_list = list(parameters1,parameters_r,parameters_c,parameters_m)
par_list = c(par_list,par_list,par_list)

for(d in 5:8){
  par_list[[d]]["beta"] = 0.144
  print(par_list[[d]])
}

for(d in 9:12){
  par_list[[d]]["beta"] = 0.222
  print(par_list[[d]])
}

par_list

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


for(d in 1:12){
  parameters1=unlist(par_list[d])
  print(paste("d =",d))
  print(parameters1)
  for(i in 1:length(imp_pA)){
    parameters1["p_A"] = imp_pA[i]
    parameters1["p_B"] = imp_pB[i]
    parameters1["p_s"] = 1-(parameters1["p_A"]+parameters1["p_B"])
    #print(paste("i=",i))
    #print(parameters1)
    run1<-ode(y=init.state,times=seq(1,(years*365),by=1), func="derivs",parms= parameters1,dllname="ESBL", initfunc="initmod",nout=15,
              outnames=c("Sum","prop_o","prop_u","prev_s","prev_a","prev_b","fCs","fCa","fCb","prevI_s","prevI_a","prevI_b","fIs","fIa","fIb"),method=rkMethod(method = 'rk4',hini=0.001))
    cols = which(colnames(data_com[[d]]$prop_o)==imp_pB[i])
    rows = which(rownames(data_com[[d]]$prop_o)==imp_pA[i])
    run_y = sapply(data.frame(run1[(365*(years-1)):(365*years),]),FUN=mean)
    admis_y = run1[(365*years),"admis"]-run1[(365*(years-1)),"admis"]
        data_com[[d]]$prop_o[rows,cols] = ifelse(run_y["prop_o"]<0.0001,0,run_y["prop_o"])
        data_com[[d]]$prop_u[rows,cols] = ifelse(run_y["prop_u"]<0.0001,0,run_y["prop_u"])
        data_com[[d]]$prev_s[rows,cols] = ifelse(run_y["prev_s"]<0.0001,0,run_y["prev_s"])
        data_com[[d]]$prev_a[rows,cols] =  ifelse(run_y["prev_a"]<0.0001,0,run_y["prev_a"])
        data_com[[d]]$prev_b[rows,cols] =  ifelse(run_y["prev_b"]<0.0001,0,run_y["prev_b"])
        data_com[[d]]$inc_s[rows,cols] = (run1[(365*years),"Inc_s"]-run1[(365*(years-1)),"Inc_s"])/admis_y*1000
        data_com[[d]]$inc_a[rows,cols] = (run1[(365*years),"Inc_a"]-run1[(365*(years-1)),"Inc_a"])/admis_y*1000
        data_com[[d]]$inc_b[rows,cols] = (run1[(365*years),"Inc_b"]-run1[(365*(years-1)),"Inc_b"])/admis_y*1000
        data_com[[d]]$prevI_s[rows,cols] = ifelse(run_y["prevI_s"]<0.0001,0,run_y["prevI_s"])
        data_com[[d]]$prevI_a[rows,cols] = ifelse(run_y["prevI_a"]<0.0001,0,run_y["prevI_a"])
        data_com[[d]]$prevI_b[rows,cols] =  ifelse(run_y["prevI_b"]<0.0001,0,run_y["prevI_b"])
        data_com[[d]]$fIs[rows,cols] = run_y["fIs"]
        data_com[[d]]$fIa[rows,cols] = run_y["fIa"]
        data_com[[d]]$fIb[rows,cols] = run_y["fIb"]
        data_com[[d]]$prop_ia[rows,cols] = run_y["prop_o"]+run_y["prop_u"]
        
    }
}

#Relative difference ###
data_dif = data_com


for(d in 1:length(data_com)){
  if(d %in% c(1:4)){
    for(v in c(1:8,12:17)){
      data_dif[[d]][[v]] = round((data_dif[[d]][[v]]-data_com[[1]][[v]])/data_com[[1]][[v]]*100,2) 
      data_dif[[d]][[v]] = ifelse(data_dif[[d]][[v]]=="NaN", 0, data_dif[[d]][[v]])
    }
  }
  else if(d %in% c(5:8)){
    for(v in c(1:8,12:17)){
      data_dif[[d]][[v]] = round((data_dif[[d]][[v]]-data_com[[5]][[v]])/data_com[[5]][[v]]*100,2)
      data_dif[[d]][[v]] = ifelse(data_dif[[d]][[v]]=="NaN", 0, data_dif[[d]][[v]])
      
    }
  }
  else{
    for(v in c(1:8,12:17)){
      data_dif[[d]][[v]] = round((data_dif[[d]][[v]]-data_com[[9]][[v]])/data_com[[9]][[v]]*100,2)
      data_dif[[d]][[v]] = ifelse(data_dif[[d]][[v]]=="NaN", 0, data_dif[[d]][[v]])
    }
  }
}

# WHICH STRAIN DOMINANT ###
for(d in 1:length(data_com)){
  data_com[[d]][["domI"]] = ifelse(data_com[[d]][["inc_s"]]>data_com[[d]][["inc_a"]] & data_com[[d]][["inc_s"]]>data_com[[d]][["inc_b"]], 1, 
                     ifelse(data_com[[d]][["inc_a"]] >data_com[[d]][["inc_s"]] &data_com[[d]][["inc_a"]] > data_com[[d]][["inc_b"]], 2, 3))
  data_com[[d]][["domP"]] = ifelse(data_com[[d]][["prev_s"]]>data_com[[d]][["prev_a"]] & data_com[[d]][["prev_s"]]>data_com[[d]][["prev_b"]], 1, 
                     ifelse(data_com[[d]][["prev_a"]] > data_com[[d]][["prev_s"]] & data_com[[d]][["prev_a"]] > data_com[[d]][["prev_b"]], 2, 3))
}

# RESULTS UNDER BASELINE WITH VARIOUS BETAs - ABSOLUTE VALUES
data_base_abs = data.frame(cbind(scen = c(1:length(data_dif)), trans = rep(NA,length(data_dif)), intv =  rep(NA,length(data_dif)), prop_o =  rep(NA,length(data_dif)),
                             prop_u =  rep(NA,length(data_dif)),inc_s =  rep(NA,length(data_dif)),inc_a =  rep(NA,length(data_dif)),
                             inc_b =  rep(NA,length(data_dif)), prev_s = rep(NA,length(data_dif)),prev_a =  rep(NA,length(data_dif)),
                             prev_b = rep(NA,length(data_dif)),prevI_s = rep(NA,length(data_dif)),prevI_a =  rep(NA,length(data_dif)),
                             prevI_b = rep(NA,length(data_dif)),fIs = rep(NA,length(data_dif)),fIa =  rep(NA,length(data_dif)),
                             fIb = rep(NA,length(data_dif)),prop_ia = rep(NA,length(data_dif))))

data_base = data_base_abs


for(d in 1:length(data_dif)){
  data_base_abs$trans[d] = ifelse(d %in% c(1:4), "Low", ifelse(d %in% c(5:8), "Medium", "High"))
  data_base_abs$intv[d] = ifelse(d %in% c(1,5,9), "Standard care", ifelse(d %in% c(2,6,10), "Informed prescribing", 
                                                                      ifelse(d %in% c(3,7,11), "Cycling", "Mixing")))
  
  data_base_abs[d,c(4:18)] = c(data_com[[d]]$prop_o[scenAn,scenBn],data_com[[d]]$prop_u[scenAn,scenBn],data_com[[d]]$inc_s[scenAn,scenBn],
                           data_com[[d]]$inc_a[scenAn,scenBn],data_com[[d]]$inc_b[scenAn,scenBn],data_com[[d]]$prev_s[scenAn,scenBn],
                           data_com[[d]]$prev_a[scenAn,scenBn],data_com[[d]]$prev_b[scenAn,scenBn],data_com[[d]]$prevI_s[scenAn,scenBn],
                           data_com[[d]]$prevI_a[scenAn,scenBn],data_com[[d]]$prevI_b[scenAn,scenBn],data_com[[d]]$fIs[scenAn,scenBn],
                           data_com[[d]]$fIa[scenAn,scenBn],data_com[[d]]$fIb[scenAn,scenBn],data_com[[d]]$prop_ia[scenAn,scenBn]) 
}
# RESULTS UNDER BASELINE WITH VARIOUS BETAs - RELATIVE DIFFERENCE
for(d in 1:length(data_dif)){
  data_base$trans[d] = ifelse(d %in% c(1:4), "Low", ifelse(d %in% c(5:8), "Medium", "High"))
  data_base$intv[d] = ifelse(d %in% c(1,5,9), "Standard care", ifelse(d %in% c(2,6,10), "Informed prescribing", 
                                                                      ifelse(d %in% c(3,7,11), "Cycling", "Mixing")))
  
  data_base[d,c(4:18)] = c(data_dif[[d]]$prop_o[scenAn,scenBn],data_dif[[d]]$prop_u[scenAn,scenBn],data_dif[[d]]$inc_s[scenAn,scenBn],
                           data_dif[[d]]$inc_a[scenAn,scenBn],data_dif[[d]]$inc_b[scenAn,scenBn],data_dif[[d]]$prev_s[scenAn,scenBn],
                           data_dif[[d]]$prev_a[scenAn,scenBn],data_dif[[d]]$prev_b[scenAn,scenBn],data_dif[[d]]$prevI_s[scenAn,scenBn],
                           data_dif[[d]]$prevI_a[scenAn,scenBn],data_dif[[d]]$prevI_b[scenAn,scenBn],data_dif[[d]]$fIs[scenAn,scenBn],
                           data_dif[[d]]$fIa[scenAn,scenBn],data_dif[[d]]$fIb[scenAn,scenBn],data_dif[[d]]$prop_ia[scenAn,scenBn]) 
}

data_base$intv = factor(data_base$intv, levels=c("Standard care", "Informed prescribing", "Cycling", "Mixing"))
data_base_abs$intv = factor(data_base_abs$intv, levels=c("Standard care", "Informed prescribing", "Cycling", "Mixing"))

data_base$trans = factor(data_base$trans, levels=c("Low", "Medium", "High"))
data_base_abs$trans = factor(data_base_abs$trans, levels=c("Low", "Medium", "High"))

setwd("~/Dropbox/RGNOSIS/Methods/Empiric prescribing model/Model/Output_C/")
save(file=paste0("sens_com_",name,".RData"),data_com)
save(file=paste0("sens_com_reldif_",name,".RData"),data_dif)
save(file=paste0("sens_trans_",name,".Rdata"), data_base_abs)
save(file=paste0("sens_trans_reldif_",name,".Rdata"), data_base)

