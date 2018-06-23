    #!bin/env Rscript
########################################################
# EMPERIC PRESCRIBING MODEL                            #
########################################################

# Author: E van Kleef
# Date: September 2015

#rm(list=ls())
require(deSolve)
require(manipulate)
#require(reshape2); require(ggplot2)

###########################
# INITIALISE STATES
###########################
#!bin/env Rscript
SIR_initialiseState <- function(theta) {
  
  # constant pop size
  N <- theta[["N"]]
  
  # number of infected and immune
  Cs <- round(theta[["pCs0"]]*N)
  Ca <- round(theta[["pCa0"]]*N)
  Cb <- round(theta[["pCb0"]]*N)
  Is <- round(theta[["pIs0"]]*N)
  Ia <- round(theta[["pIa0"]]*N)
  Ib <- round(theta[["pIb0"]]*N)
  
  if(Cs+Ca+Cb+Is+Ia+Ib>N){
    stop("Initial conditions not valid")
  }
  
  return(c(S=N-Cs-Ca-Cb-Is-Ia-Ib,Cs=Cs,Ca=Ca, Cb=Cb,Is_a=Is*theta[["f_isA"]],Is_o=Is*(1-theta[["f_isA"]]),
           Ia_u=Ia*theta[["f_iaA"]],Ia_a=Ia*(1-theta[["f_iaA"]]), Ib_u=Ib*(1-theta[["f_ib_last"]]),Ib_a=Ib*theta[["f_ib_last"]],
           Total=N, Inc_s=0, Inc_a=0, Inc_b=0,admis=0, foi_s = 0,foi_a = 0,foi_b = 0)) # death=0, 
}

################################################
# MODEL TRAJECTORY
################################################

# Conventional
SCI_simulateDeterministic <- function(theta,state.init,times) {
  # The following function compute the derivative of the ODE system
  SCI_ode <- function(time,state,theta) {
    # parms
    beta = theta[["beta"]]
    d = 1/theta[["d"]]
    d_i = (1/theta[["d"]])*theta[["d_i"]]
    N = theta[["N"]]
    #m=theta[["m"]]
    #m_i=theta[["m"]]+theta[["m_i"]]    
    a_s = 1-(theta[["a_c"]]+theta[["a_i"]]) # Fraction of total susceptible patients on admission
    a_c = theta[["a_c"]]                  # Fraction of total colonised patients on admission
    a_i = theta[["a_i"]]                # Fraction of total infected patients on admission
    p_s = 1-(theta[["p_A"]]+theta[["p_B"]]) # Among the colonised/infected, the fraction that is infected with a susceptible strain
    p_A = theta[["p_A"]]           # Among the colonised/infected, the fraction that is infected with an ESBL strain
    p_B = theta[["p_B"]]            # Among the colonised/infected, the fraction that is infected with a carbapenem resistant strain
    c1 = theta[["c1"]]
    c2 = theta[["c2"]]
    s_1 = theta[["s_1"]] # Fraction of patients that carry an ESBL carrying KPN undetectable
    s_2 = theta[["s_2"]] # Fraction of patients that carry a carb resistant carrying KPN undetectable
    f_cA =theta[["f_cA"]]
    f_cB =theta[["f_cB"]]
    f_clast =theta[["f_clast"]]
    f_isA = theta[["f_isA"]]
    f_isB= 1-theta[["f_isA"]]
    f_iaA = theta[["f_iaA"]]
    f_iaB= 1-theta[["f_iaA"]]
    f_ib_last = theta[["f_ib_last"]]
    p_c= theta[["p_c"]]
    r=1/theta[["r"]]
    r_delay=ifelse(theta[["ts"]]<0, 0, 1/(theta[["r"]]+theta[["ts"]]))
    cycl_period = theta["cycl_period"]
    alt_cycl_period=theta["alt_cycl_period"]
    t=theta[["t"]]
    antib_cycling_scheme = data.frame(cbind(Time=c(1:t), reg=rep(c(rep(1,cycl_period), rep(0,alt_cycl_period)),10000)[1:t]))
    cycl_reg = antib_cycling_scheme$reg[time+1]
    cycling=theta[["cycling"]]
    rdt = theta[["rdt"]]
    parms = c(beta, d, d_i, N, a_s, a_c, a_i,p_s,p_A,p_B, c1, c2,
              s_1,s_2, f_cA, f_cB, f_clast, f_isA, f_isB,f_iaA, f_iaB, f_ib_last, p_c, r, r_delay, cycl_reg, cycling, rdt)
    #m, m_i,
    # Init states
    S = state[["S"]]
    Cs = state[["Cs"]]
    Ca = state[["Ca"]]
    Cb = state[["Cb"]]
   # Is = state[["Is"]]
    Is_a = state[["Is_a"]]
    Is_o = state[["Is_o"]]
  #  Ia = state[["Ia"]]
    Ia_u = state[["Ia_u"]]
    Ia_a = state[["Ia_a"]]
   # Ib = state[["Ib"]]
    Ib_u = state[["Ib_u"]]
    Ib_a = state[["Ib_a"]]
    total=state[["Total"]]
    Inc_s = state[["Inc_s"]]
    Inc_a = state[["Inc_a"]]
    Inc_b = state[["Inc_b"]]
    #death = state[["death"]]
    admis = state[["admis"]]
    foi_s = state[["foi_s"]]
    foi_a = state[["foi_a"]]
    foi_b = state[["foi_b"]]
    if(cycling==1){
      f_isA=cycl_reg
      f_iaA=cycl_reg
      f_isB = 1-f_isA
      f_iaB = 1-f_iaA
    }
    
    # Functions
    dDS = (d)*S #+m
    dDCs = (d)*Cs #+m
    dDCa = (d)*Ca #+m
    dDCb = (d)*Cb #+m
    dDIs_o = (d_i)*Is_o #+m_i
    dDIs_a = (d_i)*Is_a #+m_i
    dDIa_u = (d_i)*Ia_u #+m_i
    dDIa_a = (d_i)*Ia_a #+m_i
    dDIb_u = (d_i)*Ib_u #+m_i
    dDIb_a = (d_i)*Ib_a #+m_i
    dis = dDS+dDCs+dDCa+dDCb+dDIs_o+dDIs_a+dDIa_u+dDIa_a+dDIb_u+dDIb_a

      dS  = a_s*dis-((beta*S)*(Cs+Is_a+Is_o))/N-((beta*c1*S)*(Ca+Ia_u+Ia_a))/N-((beta*c2*S)*(Cb+Ib_u+Ib_a))/N+ 
            r*(1-s_1)*(f_cA*Cs + Is_a) + r*(1-s_2)*(f_cB*(Cs+Ca) + Is_o + Ia_a) + r_delay*(1-s_2)*Ia_u+
            r_delay*Ib_u+r*Ib_a + r*f_clast*Cb - dDS
    
      dCs = a_c*p_s*dis+(beta*S*(Cs+Is_a+Is_o))/N-p_c*Cs - r*(f_cA+f_cB)*Cs - dDCs
      dCa = a_c*p_A*dis+(beta*c1*S*(Ca+Ia_u+Ia_a))/N-p_c*Ca - r*f_cB*Ca + r*s_1*f_cA*Cs - dDCa 
      dCb = a_c*p_B*dis+((beta*c2*S)*(Cb+Ib_u+Ib_a))/N-p_c*Cb - r*f_clast*Cb + r*s_2*f_cB*(Cs+Ca) - dDCb
     # dIs = a_i*p_s*dis+p_c*Cs - Is 
      dIs_o =  (a_i*p_s*dis)*f_isB + p_c*Cs*f_isB-r*Is_o-dDIs_o
      dIs_a =  (a_i*p_s*dis)*f_isA + p_c*Cs*f_isA-r*Is_a-dDIs_a
     # dIa = a_i*p_A*dis+p_c*Ca + r*s_1*Is_a - Ia 
      dIa_u = (a_i*p_A*dis)*f_iaA + p_c*Ca*f_iaA + r*s_1*Is_a*f_iaA - r_delay*Ia_u - dDIa_u  
      dIa_a = (a_i*p_A*dis)*f_iaB + p_c*Ca*f_iaB + r*s_1*Is_a*f_iaB - r*Ia_a - dDIa_a
      #dIb = a_i*p_B*dis+p_c*Cb + r*s_2*(Is_o+Ia_a)+r_delay*s_2*Ia_u-Ib
      dIb_u = (a_i*p_B*dis)*(1-f_ib_last)+p_c*Cb*(1-f_ib_last) + r*s_2*(Is_o+Ia_a)*(1-f_ib_last)+r_delay*s_2*Ia_u*(1-f_ib_last)- r_delay*Ib_u-dDIb_u  
      dIb_a = (a_i*p_B*dis)*f_ib_last+p_c*Cb*f_ib_last + r*s_2*(Is_o+Ia_a)*f_ib_last+r_delay*s_2*Ia_u*f_ib_last- r*Ib_a - dDIb_a
      dTotal = dS+dCs+dCa+dCb+dIs_o+dIs_a+dIa_u+dIa_a+dIb_u+dIb_a

    dInc_s = beta*S*(Cs+Is_a+Is_o)/N#+alpha*S
    dInc_a = beta*c1*S*(Ca+Ia_u+Ia_a)/N#+alpha*S
    dInc_b = beta*S*c2*(Cb+Ib_u+Ib_a)/N#+alpha*S
    #death = m*(S+Cs+Ca+Cb)+m_i*(Is_a+Is_o+Ia_u+Ia_a+Ib_u+Ib_a) 
    admis = dis
    dfoi_s = beta*((Cs+Is_a+Is_o)/N)
    dfoi_a = beta*c1*((Ca+Ia_u+Ia_a)/N)
    dfoi_b = beta*c2*((Cb+Ib_u+Ib_a)/N)
    return(list(c(dS,dCs,dCa,dCb,dIs_a,dIs_o,dIa_u,dIa_a,dIb_u,dIb_a,dTotal, 
                  dInc_s,dInc_a,dInc_b, admis, dfoi_s, dfoi_a, dfoi_b))) #, death
  }

  # simulate and return a data.frame
  trajectory <- data.frame(ode(y=state.init,times=times,func=SCI_ode,parms=theta))
  trajectory$Is          = trajectory$Is_a+trajectory$Is_o
  trajectory$Ia          = trajectory$Ia_a+trajectory$Ia_u
  trajectory$Ib          = trajectory$Ib_a+trajectory$Ib_u
  trajectory$Inc_s_p1000 = (trajectory$Inc_s/trajectory$admis)*1000
  trajectory$Inc_a_p1000 = (trajectory$Inc_a/trajectory$admis)*1000
  trajectory$Inc_b_p1000 = (trajectory$Inc_b/trajectory$admis)*1000
  #trajectory$death_p1000 = (trajectory$death/trajectory$admis)*1000
  trajectory$prev_s      = (trajectory$Cs+trajectory$Is_a+trajectory$Is_o)/trajectory$Total[1]
  trajectory$prev_a      = (trajectory$Ca+trajectory$Ia_u+trajectory$Ia_a)/trajectory$Total[1]
  trajectory$prev_b      = (trajectory$Cb+trajectory$Ib_u+trajectory$Ib_a)/trajectory$Total[1]
  trajectory$prop_s      = trajectory$S/trajectory$Total[1]
  trajectory$prop_a      = (trajectory$Is_a + trajectory$Ia_a + trajectory$Ib_a)/(trajectory$Is_a+trajectory$Is_o + trajectory$Ia_u + trajectory$Ia_a + trajectory$Ib_u + trajectory$Ib_a)
  trajectory$prop_u      = (trajectory$Ia_u + trajectory$Ib_u)/(trajectory$Is_a+trajectory$Is_o + trajectory$Ia_u + trajectory$Ia_a + trajectory$Ib_u + trajectory$Ib_a)
  trajectory$prop_o      = trajectory$Is_o/(trajectory$Is_a+trajectory$Is_o + trajectory$Ia_u + trajectory$Ia_a + trajectory$Ib_u + trajectory$Ib_a)                
  return(trajectory)
}

theta=c(beta=0.06,N=100,t=365,
        d=7, d_i=7*2, a_c = 0.044, a_i= 0.004, p_A= 0.297, p_B=0.04,c1=1,c2=1, f_cA=0.1,f_cB=0.1,
        f_isA=0.82,f_iaA=0.56, f_ib_last = 0, r=5, ts=3,f_clast=0.03,
        s_1=0.001,s_2=0.001, cycl_period=365*2, alt_cycl_period=365*2,cycling=0,rdt=0,
        p_c=0.01, pS0=0.95, pCs0=0.05,pCa0=0,pCb0=0, pIs0=0, pIa0=0,pIb0=0)

# state.init=SIR_initialiseState(theta)
# times = c(0:365)
# #traj <- SCI_simulateDeterministic(theta,state.init, times=0:365)
# 
# 
# ## for plotting the ode results
# plot.ode<-function(theta,state.init,times, beta, a_c, a_i, c1, c2, s_1,s_2,p_A,p_B,f_cA,f_cB){
#   #solve equations
#   theta[["beta"]] = beta
#   theta[["a_c"]] = a_c
#   theta[["a_i"]] = a_i
#   theta[["c1"]] = c1
#   theta[["c2"]] = c2
#   theta[["s_1"]] = s_1
#   theta[["s_2"]] = s_2
#   theta[["p_A"]] = p_A
#   theta[["p_B"]] = p_B
#   theta[["f_cA"]] = f_cA
#   theta[["f_cB"]] = f_cB
#   out<-SCI_simulateDeterministic(theta,state.init, times=times)
#   #plot the results
#   par(mfrow=c(2,2))
#   for(i in c("Inc_a_p1000","Inc_b_p1000", "prop_u","prop_o")){
#   plot(times,out[,which(names(out)%in%i)], main=i, xlab="Time", ylab=i, type="l") 
#   }
# }
# 
# ## show the results with sliders
# manipulate(plot.ode(theta,state.init, times, beta, a_c, a_i, c1, c2, s_1,s_2,p_A,p_B,f_cA,f_cB),
#            beta=slider(0,1,step=0.01,label="Transmission rate"),
#            a_c=slider(0,1,step=0.1,label="Fraction colonised on admission"),
#            a_i=slider(0,1,step=0.1,label="Fraction infected on admission"),
#            c1=slider(0,1,step=0.1,label="ESBL fitness cost"),
#            c2=slider(0,1,step=0.1,label="Carb-R fitness cost"),
#            s_1=slider(0,1,step=0.1,label="pre-existing ESBL"),         
#            s_2=slider(0,1,step=0.1,label="pre-existing Carb-R"),          
#            p_A=slider(0,1,step=0.1,label="Fraction ESBL among colonised on admission") ,         
#            p_B=slider(0,1,step=0.1,label="Fraction Carb-R among colonised on admission") ,         
#            f_cA=slider(0,1,step=0.1,label="Fraction on background A antibiotic")  ,        
#            f_cB=slider(0,1,step=0.1,label="Fraction on background B antibiotic")          
# )
#            
