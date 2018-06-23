#!bin/env Rscript
###########################################
# ANALYSIS EMPIRIC PRESCRIBING MODEL
###########################################
rm(list = ls())
source("./Run_model_code/Run_scenarios/ESBL_amr_model_ODE_FINAL.R")
theta_all = read.csv("./Scenarios_parameters/cycl_all_final_low.csv")
sens = read.csv("./Scenarios_parameters/sens_final_low.csv")

theta_all$p_A =ifelse(theta_all$p_A==0.399,0.33,theta_all$p_A)
theta_all$p_B =ifelse(theta_all$p_B==0.08,0.05,theta_all$p_B)
theta_all$p_c =ifelse(theta_all$p_c==0.01,0.1,theta_all$p_c)
theta_all$a_c =ifelse(theta_all$a_c==0.05,0.2,theta_all$a_c)
theta_all$a_i =ifelse(theta_all$a_i==0.002,0.012,theta_all$a_i)

#########################
# RUN MODEL
#########################
#set  = Sys.getenv("SLURM_ARRAY_TASK_ID")
for(i in 1:length(theta_all[,1])){
  print(i)
 state.init=SIR_initialiseState(theta_all[i,])

 traj <- SCI_simulateDeterministic(theta_all[i,],state.init, times=0:(15*365))
 output = traj[c(((length(traj$time)-1)-359):(length(traj$time)-1)),]
 traj$Intervention = "Cycling"
 
 # Incidence per 1000 at equilibrium
 Inc_s_p1000 = ((output$Inc_s[length(output$admis)]-output$Inc_s[1])/(output$admis[length(output$admis)]-output$admis[1]))*1000
 Inc_a_p1000 =((output$Inc_a[length(output$admis)]-output$Inc_a[1])/(output$admis[length(output$admis)]-output$admis[1]))*1000
 Inc_b_p1000 =((output$Inc_b[length(output$admis)]-output$Inc_b[1])/(output$admis[length(output$admis)]-output$admis[1]))*1000 # New
 
 output$Inc_s_p1000 = Inc_s_p1000
 output$Inc_a_p1000 = Inc_a_p1000
 output$Inc_b_p1000 = Inc_b_p1000
 
 output = as.data.frame(t(sapply(output, FUN=mean, na.rm=T)))
 output$var_change = sens$var_change[i]
 output$var_value = sens$var_value[i]
# write.table(output, "./Output/Final/Sens/Output_empiric_base_antib_cycl_low_new.csv", append=T, sep = ",",row.names=T,col.names=ifelse(i==1,NA,F))
# write.csv(traj, paste0("./Output/Sensitivity analysis/Trajectories/Traj_empiric_base_antib_cycl_low_new_",i,".csv"))
}