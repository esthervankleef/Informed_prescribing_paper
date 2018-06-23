
//  seir model
//
//  Created by Esther van Kleef
// To compile from R: system("R CMD SHLIB ESBL.c")
// To load from R (in Unix): dyn.load("ESBL.so")
// To call from R use deSolve package and see details in help on using compiled code.
// Make sure to set the work directory to the .c file

#include <R.h>
#include <math.h>
#include <stdio.h>


static double parms[29];

#define beta parms[0]
#define d parms[1]
#define d_i parms[2]
#define a_s parms[3]
#define a_c parms[4]
#define a_i parms[5]
#define p_s parms[6]
#define p_A parms[7]
#define p_B parms[8]
#define c1 parms[9]
#define c2 parms[10]
#define s_1 parms[11]
#define s_2 parms[12]
#define f_cA parms[13]
#define f_cB parms[14]
#define f_clast parms[15]
#define f_isA parms[16]
#define f_isB parms[17]
#define f_iaA parms[18]
#define f_iaB parms[19]
#define f_ib_last parms[20]
#define p_c parms[21]
#define r parms[22]
#define r_delay parms[23]
#define cycl_period parms[24]
#define cycling parms[25]
#define rdt parms[26]
#define ts parms[27]
#define mixing parms[28]



/* initializer */
void initmod(void (* odeparms)(int * , double *))
{
    int N=29;
    odeparms(&N, parms);
}

/*Function to change parameters for cycling */


/* derivatives and one output variable 
y[0] is S (susceptibles)
y[1] is Cs 
y[2] is Ca
y[3] is Cb
y[4] is Is_a 
y[5] is Is_o 
y[6] is Ia_u 
y[7] is Ia_a
y[8] is Ib_u 
y[9] is Ib_a
y[10] is Inc_s
y[11] is Inc_a
y[12] is Inc_b
y[13] is admis
y[14] is fraction Cs
y[15] is fraction Ca
y[16] is fraction Cb
*/

void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if(ip[0]<1) error("nout should be at least 1") ;
    if(y[1]<1 && *t==1) error("Cs should be at least 1") ;
    if(y[4]<1 && *t==1) error("Is_a should be at least 1") ;
    double N;
    N= y[0]+y[1]+y[2]+y[3]+y[4]+y[5]+y[6]+y[7]+y[8]+y[9];

    double dDS;
    dDS = d*y[0] ;
    double dDCs;
    dDCs = d*y[1];
    double dDCa;
    dDCa = d*y[2];
    double dDCb;
    dDCb = d*y[3];
    double dDIs_o;
    dDIs_o = d_i*y[5];
    double dDIs_a;
    dDIs_a = d_i*y[4];
    double dDIa_u;
    dDIa_u = d_i*y[6];
    double dDIa_a;
    dDIa_a = d_i*y[7];
    double dDIb_u;
    dDIb_u = d_i*y[8];
    double dDIb_a;
    dDIb_a = d_i*y[9];
    double dis ;
    dis = dDS+dDCs+dDCa+dDCb+dDIs_o+dDIs_a+dDIa_u+dDIa_a+dDIb_u+dDIb_a;
    double time = *t;
    
    double weeks ;
    weeks = ceil(time/cycl_period) ;
  
    double nf_isA = f_isA;
    double nf_iaA = f_iaA;
    double nf_isB = f_isB;
    double nf_iaB = f_iaB;
    double nf_csA = f_cA ;
    double nf_csB = f_cB ;
    double nf_caB = f_cB ;
    double nf_cslast = f_clast;
    double nf_calast = f_clast;
    double nf_cblast = f_clast;
    double nf_cA = f_cA;
    double nf_cB = f_cB;
    double nf_clast = f_clast;
    
    if (cycling == 1){ 
      if ((int) weeks % 2 == 0){
      nf_isA= 0 ;
      nf_iaA= 0 ;
      nf_isB = 1 ;
      nf_iaB = 1 ;
      nf_csA = 0 ;
      nf_csB = f_cA+f_cB;
      nf_caB = f_cA+f_cB;

      }
    else{
      nf_isA= 1 ;
      nf_iaA= 1 ;
      nf_isB = 0 ;
      nf_iaB = 0 ;
      nf_csA = f_cA+f_cB;
      nf_csB = 0;
      nf_caB = 0;
      }
    }
    
    if(mixing == 1){
      nf_isA= 0.5 ;
      nf_iaA= 0.5 ;
      nf_isB = 0.5 ;
      nf_iaB = 0.5 ;
      nf_csA = (f_cA+f_cB)/2;
      nf_csB = (f_cA+f_cB)/2;
      nf_caB = (f_cA+f_cB)/2;
    }
    
    if(rdt == 1 || (cycling == 0 && mixing == 0 && rdt == 0)){
      nf_cA = (yout[12]*f_isA+yout[13]*f_iaA+(yout[14]*(1-f_ib_last)))*(f_cA+f_cB+f_clast);
      nf_cB = (yout[12]*f_isB+yout[13]*f_iaB)*(f_cA+f_cB+f_clast);
      nf_clast = (yout[14]*f_ib_last)*(f_cA+f_cB+f_clast);
      nf_csA = nf_cA ;
      nf_csB = nf_cB ;
      nf_caB = nf_cB ;
      nf_cslast = nf_clast;
      nf_calast = nf_clast;
      nf_cblast = nf_clast;
    }
    /*printf("f_cA %f\n",nf_cA);
    printf("f_cB %f\n",nf_cB);
    printf("f_clast %f\n",nf_clast);

    printf("f_isA %f\n",nf_csA);
    printf("%f\n",weeks);
    printf("time %f\n",time);
    printf("t %f\n",*t);
    printf("t-1 %f\n",*t-1);
*/  ydot[0] = a_s*dis-((beta*y[0])*(y[1]+y[4]+y[5]))/N-((beta*c1*y[0])*(y[2]+y[6]+y[7]))/N-((beta*c2*y[0])*(y[3]+y[8]+y[9]))/N+ 
      r*(1-s_1)*(nf_csA*y[1] + y[4]) + r*(1-s_2)*(nf_csB*y[1]+nf_caB*y[2] + y[5] + y[7]) + r_delay*(1-s_2)*y[6]+
      r_delay*y[8]+r*y[9] + r*(nf_cslast*y[1]+nf_calast*y[2]+nf_cblast*y[3]) - dDS; // S
    ydot[1] =  a_c*p_s*dis+(beta*y[0]*(y[1]+y[4]+y[5]))/N-p_c*y[1] - r*(nf_csA+nf_csB+nf_cslast)*y[1] - dDCs;  // Cs
    ydot[2] =  a_c*p_A*dis+(beta*c1*y[0]*(y[2]+y[6]+y[7]))/N-p_c*y[2] - r*(nf_caB+nf_calast)*y[2] + r*s_1*nf_csA*y[1] - dDCa ;  //Ca
    ydot[3] =  a_c*p_B*dis+((beta*c2*y[0])*(y[3]+y[8]+y[9]))/N-p_c*y[3] - r*nf_cblast*y[3] + r*s_2*(nf_csB*y[1]+nf_caB*y[2]) - dDCb; //Cb
    ydot[4] = (a_i*p_s*dis)*nf_isA + p_c*y[1]*nf_isA-r*y[4]-dDIs_a; // Is with appropriate treatment
    ydot[5] = (a_i*p_s*dis)*nf_isB + p_c*y[1]*nf_isB-r*y[5]-dDIs_o; // Is with overtreatment
    ydot[6] = (a_i*p_A*dis)*nf_iaA + p_c*y[2]*nf_iaA + r*s_1*y[4]*nf_iaA - r_delay*y[6] - dDIa_u ; // Ia with undertreatment
    ydot[7] = (a_i*p_A*dis)*nf_iaB + p_c*y[2]*nf_iaB + r*s_1*y[4]*nf_iaB - r*y[7] - dDIa_a; // Ia with appropriate treatment
    ydot[8] = (a_i*p_B*dis)*(1-f_ib_last)+p_c*y[3]*(1-f_ib_last) + r*s_2*(y[5]+y[7])*(1-f_ib_last)+r_delay*s_2*y[6]*(1-f_ib_last)- r_delay*y[8]-dDIb_u  ; // Ib with undertreatment
    ydot[9] =(a_i*p_B*dis)*f_ib_last+p_c*y[3]*f_ib_last + r*s_2*(y[5]+y[7])*f_ib_last+r_delay*s_2*y[6]*f_ib_last- r*y[9] - dDIb_a; // Ib with overtreatment
    ydot[10] = beta*y[0]*(y[1]+y[4]+y[5])/N ;
    ydot[11] = beta*c1*y[0]*(y[2]+y[6]+y[7])/N ;
    ydot[12] = beta*c2*y[0]*(y[3]+y[8]+y[9])/N ; 
    ydot[13] = d*(y[0]+y[1]+y[2]+y[3])+d_i*(y[4]+y[5]+y[6]+y[7]+y[8]+y[9]); /*admis*/
   
      
    yout[0]=y[0]+y[1]+y[2]+y[3]+y[4]+y[5]+y[6]+y[7]+y[8]+y[9];
    yout[1] = y[5]/(y[4]+y[5]+y[6]+y[7]+y[8]+y[9]); /*prop_o*/
    yout[2] = (y[6]+y[8])/(y[4]+y[5]+y[6]+y[7]+y[8]+y[9]); /*prop_u*/
    yout[3] = (y[1]+y[4]+y[5])/N; /*prev_s*/
    yout[4] = (y[2]+y[6]+y[7])/N; /*prev_a*/
    yout[5] = (y[3]+y[8]+y[9])/N; /*prev_b*/
    yout[6] = y[1]/(y[1]+y[2]+y[3]) ;  /*f_Cs*/
    yout[7] = y[2]/(y[1]+y[2]+y[3]) ;/*f_Ca*/
    yout[8] = y[3]/(y[1]+y[2]+y[3]) ;/*f_Cb*/
    yout[9] = (y[4]+y[5])/N; /*prevI_s*/
    yout[10] = (y[6]+y[7])/N; /*prevI_a*/
    yout[11] = (y[8]+y[9])/N; /*prevI_b*/
    yout[12] = (y[4]+y[5])/(y[4]+y[5]+y[6]+y[7]+y[8]+y[9]); /*f_Is*/
    yout[13] = (y[6]+y[7])/(y[4]+y[5]+y[6]+y[7]+y[8]+y[9]); /*f_Ia*/
    yout[14] = (y[8]+y[9])/(y[4]+y[5]+y[6]+y[7]+y[8]+y[9]); /*f_Ib*/
    
    
}

