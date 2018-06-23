#STANDARD CARE MED
theta=c(beta=0.111, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.82,f_iaA=0.56, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=1,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#STANDARD CARE LOW
theta=c(beta=0.027, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.82,f_iaA=0.56, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=1,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#STANDARD CARE HIGH
theta=c(beta=0.24, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.82,f_iaA=0.56, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=1,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#STANDARD CARE MED - NO CLEARANCE
theta=c(beta=0.111, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.82,f_iaA=0.56, r=6, ts=365,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=1,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)


#RDT MED
theta=c(beta=0.111, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.95,f_iaA=0.15, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0.85, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=1,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

# RDT LOW
theta=c(beta=0.027, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.95,f_iaA=0.15, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0.85, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=1,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

# RDT HIGH
theta=c(beta=0.24, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.95,f_iaA=0.15, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0.85, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=1,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#Cycl MED
theta=c(beta=0.111, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0,f_iaA=0, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=30, alt_cycl_period=30,cycling=1,rdt=0,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)
#Cycl LOW
theta=c(beta=0.027, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0,f_iaA=0, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=30, alt_cycl_period=30,cycling=1,rdt=0,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#Cycl HIGH
theta=c(beta=0.24, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0,f_iaA=0, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=30, alt_cycl_period=30,cycling=1,rdt=0,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#Cycl MED - NO CLEARANCE
theta=c(beta=0.11, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0,f_iaA=0, r=6, ts=365,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=30, alt_cycl_period=30,cycling=1,rdt=0,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#Mix MED
theta=c(beta=0.111, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.5,f_iaA=0.5, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=0,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#Mix LOW
theta=c(beta=0.027, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.5,f_iaA=0.5, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=0,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#Mix HIGH
theta=c(beta=0.24, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.5,f_iaA=0.5, r=6, ts=3,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=0,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)

#Mix MED - NO CLEARANCE
theta=c(beta=0.111, d=9, d_i=0.86, #m=0.01,m_i=0.01,
        a_c = 0.05, a_i= 0.002, p_A= 0.399, p_B=0.075,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_clast = 0.03, f_isA=0.5,f_iaA=0.5, r=6, ts=365,
        s_1=0.001,s_2=0.001, p_c=0.01, f_ib_last = 0, cycl_period=0, alt_cycl_period=0,cycling=0,rdt=0,
        pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0,N=60,t=365*15)


state.init=SIR_initialiseState(theta)

traj <- SCI_simulateDeterministic(theta,state.init, times=0:(365*15))

output = traj[c(((length(traj$time)-1)-359):(length(traj$time)-1)),]
#traj$Intervention = "Cycling"

# Incidence per 1000 at equilibrium
Inc_s_p1000 = ((output$Inc_s[length(output$admis)]-output$Inc_s[1])/(output$admis[length(output$admis)]-output$admis[1]))*1000
Inc_a_p1000 =((output$Inc_a[length(output$admis)]-output$Inc_a[1])/(output$admis[length(output$admis)]-output$admis[1]))*1000
Inc_b_p1000 =((output$Inc_b[length(output$admis)]-output$Inc_b[1])/(output$admis[length(output$admis)]-output$admis[1]))*1000 # New

output$Inc_s_p1000 = Inc_s_p1000
output$Inc_a_p1000 = Inc_a_p1000
output$Inc_b_p1000 = Inc_b_p1000

output = as.data.frame(t(sapply(output, FUN=mean, na.rm=T)))
output
