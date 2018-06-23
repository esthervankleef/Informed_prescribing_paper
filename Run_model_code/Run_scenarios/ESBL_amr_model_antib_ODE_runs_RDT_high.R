#!bin/env Rscript
###########################################
# ANALYSIS EMPIRIC PRESCRIBING MODEL
###########################################
rm(list = ls())
source("./Run_model_code/Run_scenarios/ESBL_amr_model_ODE_FINAL.R")
theta_all = read.csv("./Scenarios_parameters/RDT_all_final_high.csv")
sens = read.csv("./Scenarios_parameters/sens_final_high.csv")

#########################
# RUN MODEL
#########################
#set  = Sys.getenv("SLURM_ARRAY_TASK_ID")
for(i in 1:length(theta_all[,1])){
  print(i)
 state.init=SIR_initialiseState(theta_all[i,])

 traj <- SCI_simulateDeterministic(theta_all[i,],state.init, times=0:(15*365))

 output = traj[c(((length(traj$time)-1)-359):(length(traj$time)-1)),]
 traj$Intervention = "Informed prescribing"
 
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
 output$Intervention == "Informed prescribing"
 
 #write.table(output, "~/Dropbox/RGNOSIS/Methods/Empiric prescribing model/Model/Output/Final/Sens/Output_empiric_base_antib_RDT_high.csv", append=T, sep = ",",row.names=T,col.names=ifelse(i==1,NA,F))
 #write.csv(traj, paste0("~/Documents/Esther/RGNOSIS/Methods/Empiric_prescribing_model/Output/Sensitivity analysis/Trajectories/Traj_empiric_base_antib_RDT_high_",i,".csv"))
}