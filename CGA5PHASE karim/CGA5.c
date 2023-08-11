/*
 * sfuntmpl_basic.c: Basic 'C' template for a level 2 S-function.
 *
 * Copyright 1990-2013 The MathWorks, Inc.
 */


/*
 * You must specify the S_FUNCTION_NAME as the name of your S-function
 * (i.e. replace sfuntmpl_basic with the name of your S-function).
 */

#define S_FUNCTION_NAME  CGA5
#define S_FUNCTION_LEVEL 2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h" 
#include "mex.h" 
#include "math.h" 
#include "stdlib.h"
#include "stdio.h"

#define  h  1e-6
#define  PI 3.14159265358979323846264338327950288419716939937510582
#define  NB_PARAM 4

#define  STATOR_PARAM(S)     ssGetSFcnParam(S,0) 
#define  ROTOR_PARAM(S)      ssGetSFcnParam(S,1)  
#define  Lm_PARAM(S)         ssGetSFcnParam(S,2) 
#define  JFP_PARAM(S)        ssGetSFcnParam(S,3) 

real_T Rs,Rr,Ls,Lr,Lm,J,F,P   ,Tr;
real_T t,k,n1,n2  ,ew,ve,n   ,nref ,nn1,nn2,Phir_d_estim,Phir_q_estim,Te,Phir_estim;
real_T NB1,NM1,NS1,Z1,PS1,PM1,PB1;
real_T faper1[7],fapder1[7],faper2[7],fapder2[7];

real_T Is_d_max,Is_q_max,e0,ei,Is_d_ref,Is_q_ref,Phir_ref ,alpha1,alpha2,Phir_Err,veph ;


real_T Wsg,Ws,n1est,n2est,Theta_s,Id,Iq;




//========================================================================================
int_T sign(real_T m)
{
if ( m == 0 ) { return 0  ; }
if ( m >  0 ) { return 1  ; }
if ( m <  0 ) { return -1 ; }
}


//====================== PARK ==========================================================================
void PARK( real_T Va , real_T Vb , real_T Vc , real_T Vf, real_T Ve, real_T Angle ,real_T* Vd, real_T* Vq)
{   
   real_T K0 = sqrt(2.0/5) ; // Constant de PARK
   *Vd =  K0 * ( Va * cos( Angle ) + Vb * cos( Angle + 2 * PI / 5 ) + Vc * cos( Angle + 4 * PI / 5 ) + Vf * cos( Angle + 6 * PI / 5 )+ Ve * cos( Angle + 8 * PI / 5 )) ; 
   *Vq = -K0 * ( Va * sin( Angle ) + Vb * sin( Angle + 2 * PI / 5 ) + Vc * sin( Angle + 4 * PI / 5 ) + Vf * sin( Angle + 6 * PI / 5 )+ Ve * sin( Angle + 8 * PI / 5 )) ;
}
//=======================================================================================================

//====================== PARK_INV =======================================================================
void PARK_INV( real_T Vd , real_T Vq , real_T Angle ,real_T* Va, real_T* Vb , real_T* Vc , real_T* Vf, real_T* Ve)
{   
    *Va = Vd   * cos( Angle )              - Vq   * sin( Angle ) ;                // Va 
    *Vb = Vd   * cos( Angle + 2 * PI / 5 ) - Vq   * sin( Angle + 2 * PI / 5 ) ;   // Vb   
    *Vc = Vd   * cos( Angle + 4 * PI / 5 ) - Vq   * sin( Angle + 4 * PI / 5 ) ;   // Vc   
    *Vf = Vd   * cos( Angle + 6 * PI / 5 ) - Vq   * sin( Angle + 6 * PI / 5 ) ;   // Vf  
    *Ve = Vd   * cos( Angle + 8 * PI / 5 ) - Vq   * sin( Angle + 8 * PI / 5 ) ;   // Ve  
}
//=======================================================================================================


real_T max1(real_T x, real_T y){ if (x > y) return x;else return y;}
real_T min1(real_T x, real_T y){ if (x < y) return x;else return y;}

void FUZ( real_T x  , real_T* NB , real_T* NM, real_T* NS, real_T* Z, real_T* PS, real_T* PM, real_T* PB )

{ real_T a1,a2,a3 ,nbp,nbm ,nmp,nmm ,nsp,nsm,zp,zm,psp,psm,pmp,pmm,pbp,pbm;  
 
 a1 = 0.1875 ;
 a2 = 0.3750 ;
 a3 = 0.7500 ;

//==================== NB ========================
  nbp =  1;
  nbm = -( x + a2 ) / ( a3 - a2 ) ;
  *NB =  max1( min1( nbp , nbm ) , 0 ) ;
//==================== NM ======================== 
  nmp =  ( x + a3 ) / ( a3 - a2 ) ;
  nmm = -( x + a1 ) / ( a2 - a1 ) ;
  *NM =  max1( min1( nmp , nmm ) , 0 ) ;
//==================== NS ======================== 
  nsp =  ( x + a2 ) / ( a2 - a1 ) ;
  nsm = -x / a1 ;
  *NS =  max1( min1( nsp , nsm ) , 0 ) ;
//==================== Z ======================== 
  zp  =  1 + x / a1 ;
  zm  =  1 - x / a1 ;
  *Z  =  max1( min1( zp , zm ) , 0 ) ;
//==================== PS ======================== 
  psp =  x  /a1 ;
  psm =  ( a2 - x ) / ( a2 - a1 ) ;
  *PS =  max1( min1( psp , psm ) , 0 ) ;
//==================== PM ========================
  pmp =  ( x - a1 ) / ( a2 - a1 ) ;
  pmm =  ( a3 - x ) / ( a3 - a2 ) ;
  *PM =  max1( min1( pmp , pmm ) , 0 ) ;
//==================== PB ======================== 
  pbp =  ( x - a2 ) / ( a3 - a2 ) ;
  pbm =  1;
  *PB =  max1( min1( pbp , pbm ) , 0 ) ;
//================================================ 

  
}
//========================= GAINE =========================================================
real_T GAINE(real_T fapder[7],real_T faper[7],real_T gain)
{
  int_T s1,s2,s3,s4,s5,r;
  int_T  a[7] ;
  int_T numa1,numa2,numa3,numa4,numa5,numa6,numa7;
  int_T dena1,dena2,dena3,dena4,dena5,dena6,dena7;
  int i,j;
  
  
   s1 = 0.2 * gain ;
   s2 = 0.4 * gain ;
   s3 = 0.6 * gain ;
   s4 = 0.8 * gain ;
   s5 = 1.0 * gain ;
   
  
   for ( i = 0 ; i <= 6 ; i++ )
 {
       
     for ( j = 0 ; j <= 6 ; j++ ) { a[j] = 0 ; }
      
       
       
  if (faper[i] != 0 )
 {       
   if ( fapder[0] != 0 ) { a[0] = faper[i] * fapder[0] ; }
   if ( fapder[1] != 0 ) { a[1] = faper[i] * fapder[1] ; }
   if ( fapder[2] != 0 ) { a[2] = faper[i] * fapder[2] ; }
   if ( fapder[3] != 0 ) { a[3] = faper[i] * fapder[3] ; }
   if ( fapder[4] != 0 ) { a[4] = faper[i] * fapder[4] ; }
   if ( fapder[5] != 0 ) { a[5] = faper[i] * fapder[5] ; }
   if ( fapder[6] != 0 ) { a[6] = faper[i] * fapder[6] ; }
 }
  
  
  if ( i==0 ) {numa1= a[0] * s5 + a[1] * s5 + a[2] * s4 + a[3] * s3 + a[4] * s2 + a[5] * s1 + a[6] * s1;
              dena1=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];} 
  
  if ( i==1 ) {numa2=a[0]*s5+a[1]*s4+a[2]*s3+a[3]*s2+a[4]*s1+a[5]*s1+a[6]*s1;
               dena2=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==2 ) {numa3=a[0]*s4+a[1]*s3+a[2]*s2+a[3]*s1+a[4]*s1+a[5]*s1+a[6]*s2;
               dena3=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==3 ) {numa4=a[0]*s3+a[1]*s2+a[2]*s1+a[3]*s1+a[4]*s1+a[5]*s2+a[6]*s3;
               dena4=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==4 ) {numa5=a[0]*s3+a[1]*s1+a[2]*s1+a[3]*s1+a[4]*s2+a[5]*s3+a[6]*s4;
               dena5=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==5 ) {numa6=a[0]*s1+a[1]*s1+a[2]*s1+a[3]*s2+a[4]*s3+a[5]*s4+a[6]*s5;
               dena6=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==6 ) {numa7=a[0]*s1+a[1]*s1+a[2]*s2+a[3]*s3+a[4]*s4+a[5]*s5+a[6]*s5;
               dena7=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
 }
  if ((dena1+dena2+dena3+dena4+dena5+dena6+dena7)==0) {r=0;}
  else { r=(numa1+numa2+numa3+numa4+numa5+numa6+numa7)/(dena1+dena2+dena3+dena4+dena5+dena6+dena7) ;}
   return r ;
}        

//========================= FLO =========================================================
real_T FLO(real_T Is_q_ref1, real_T fapder[7],real_T faper[7],real_T alpha ,real_T gain)
{
  int_T s1,s2,s3,s4,s5,s6,s7,iq,Is_q_ref00;
  int_T  a[7] ;
  int_T num1,num2,num3,num4,num5,num6,num7;
  int_T den1,den2,den3,den4,den5,den6,den7;
  int i,j;

  
   s1 = -1.00 * gain ;
   s2 = -0.66 * gain ;
   s3 = -0.33 * gain ;
   s4=   0.00 * gain ;
   s5=   0.33 * gain ;
   s6=   0.66 * gain ;
   s7=   1.00 * gain ;
   
  
   for ( i = 0 ; i <= 6 ; i++ )
 {
       
     for ( j = 0 ; j <= 6 ; j++ ) { a[j] = 0 ; }
      
       
       
  if (faper[i] != 0 )
 {       
   if ( fapder[0] != 0 ) { a[0] = faper[i] * fapder[0] ; }
   if ( fapder[1] != 0 ) { a[1] = faper[i] * fapder[1] ; }
   if ( fapder[2] != 0 ) { a[2] = faper[i] * fapder[2] ; }
   if ( fapder[3] != 0 ) { a[3] = faper[i] * fapder[3] ; }
   if ( fapder[4] != 0 ) { a[4] = faper[i] * fapder[4] ; }
   if ( fapder[5] != 0 ) { a[5] = faper[i] * fapder[5] ; }
   if ( fapder[6] != 0 ) { a[6] = faper[i] * fapder[6] ; }
 }
  
  
  if ( i==0 ) {num1=    a[0] * s1 + a[1] * s1 + a[2] * s1 + a[3] * s1 + a[4] * s2 + a[5] * s3 + a[6] * s4 ;
               den1=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];} 
  
  if ( i==1 ) {num2=a[0]*s1+a[1]*s1+a[2]*s1+a[3]*s2+a[4]*s3+a[5]*s4+a[6]*s5;
               den2=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==2 ) {num3=a[0]*s1+a[1]*s1+a[2]*s2+a[3]*s3+a[4]*s4+a[5]*s5+a[6]*s6;
               den3=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==3 ) {num4=a[0]*s1+a[1]*s2+a[2]*s3+a[3]*s4+a[4]*s5+a[5]*s6+a[6]*s7;
               den4=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==4 ) {num5=a[0]*s2+a[1]*s3+a[2]*s4+a[3]*s5+a[4]*s6+a[5]*s7+a[6]*s7;
               den5=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==5 ) {num6=a[0]*s3+a[1]*s4+a[2]*s5+a[3]*s6+a[4]*s7+a[5]*s7+a[6]*s7;
               den6=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
  if ( i==6 ) {num7=a[0]*s4+a[1]*s5+a[2]*s6+a[3]*s7+a[4]*s7+a[5]*s7+a[6]*s7;
               den7=a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6];}
  
 }
if ((den1+den2+den3+den4+den5+den6+den7)==0) {iq =0 ;}
else { iq=(num1+num2+num3+num4+num5+num6+num7)/(den1+den2+den3+den4+den5+den6+den7);}
Is_q_ref00=Is_q_ref1+alpha*iq;
   return Is_q_ref00 ;
}        

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, NB_PARAM);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 4)) return;
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortWidth(S, 2, 1);
    ssSetInputPortWidth(S, 3, 5);
    

    ssSetInputPortRequiredContiguous(S, 0, true); /*direct input signal access*/
    ssSetInputPortRequiredContiguous(S, 1, true); /*direct input signal access*/
    ssSetInputPortRequiredContiguous(S, 2, true); /*direct input signal access*/
    ssSetInputPortRequiredContiguous(S, 3, true); /*direct input signal access*/

    /*
     * Set direct feedthrough flag (1=yes, 0=no).
     * A port has direct feedthrough if the input is used in either
     * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
     */
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortDirectFeedThrough(S, 3, 1);

  if (!ssSetNumOutputPorts(S, 5)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortWidth(S, 3, 5);
    ssSetOutputPortWidth(S, 4, 1);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Specify the sim state compliance to be same as a built-in block */
    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

    ssSetOptions(S, 0);
}



//=======================================================================================================

/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, h);
    ssSetOffsetTime(S, 0, 0.0);

}



#define MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)


  static void mdlInitializeConditions(SimStruct *S)
  {const real_T       *STATOR  = mxGetPr(STATOR_PARAM(S)); 
     const real_T       *ROTOR   = mxGetPr(ROTOR_PARAM(S));      
     const real_T       *Lm1     = mxGetPr(Lm_PARAM(S));        
     const real_T       *JFP     = mxGetPr(JFP_PARAM(S));        
     
    
//===============================================================
Rs = STATOR[0];     //% Résistance du stator
Ls = STATOR[1];     //% H Inductance du stator

Rr = ROTOR[0];      //% Résistance du rotor
Lr = ROTOR[1];      //% H Inductance du rotor

Lm = *Lm1;          //% H Inductance Mutuelle

J =  JFP[0];        //% Kg.m2 Moment d’inertie
F =  JFP[1];        //% SI Coefficient de frottement
P =  JFP[2];        //% Nombre de paire de pôle.



Ls = Ls + Lm ;
Lr = Lr + Lm ;

Tr = Lr / Rr ;
//====================================================
      
      k  = 0 ;
      t  = 0 ;
      n1 = 0 ;
      n2 = 0;
      Is_d_ref = 0;
      
      nn1 = 0 ;
      nn2 = 0;
      
Is_d_max = 5  ;
Is_q_max = 13 ;

e0 = 0 ;
ei = 0 ;

Is_d_ref = 0 ;
Is_q_ref = 0 ;

Phir_estim = 0.0001 ;



Phir_d_estim = 0.1 ;
Phir_q_estim = 0.1 ;

Ws      = 2 * PI * 50 ;
n1est   = 0 ; n2est = 0 ;
Theta_s = 0 ;
  }
#endif /* MDL_INITIALIZE_CONDITIONS */



#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
  /* Function: mdlStart =======================================================
   * Abstract:
   *    This function is called once at start of model execution. If you
   *    have states that should be initialized once, this is the place
   *    to do it.
   */
  static void mdlStart(SimStruct *S)
  {
  }
#endif /*  MDL_START */



/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block.
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    const real_T     *Wm_Ref1        = ( const real_T* )  ssGetInputPortSignal(S,0);
    const real_T     *Phir_d_ref     = ( const real_T* )  ssGetInputPortSignal(S,1);
    const real_T     *Wm1            = ( const real_T* )  ssGetInputPortSignal(S,2);
    const real_T     *Iabc1          = ( const real_T* )  ssGetInputPortSignal(S,3);
   
    real_T       *Vs_d1               = ssGetOutputPortSignal(S,0);
    real_T       *Vs_q1               = ssGetOutputPortSignal(S,1);
    real_T       *Ws1                 = ssGetOutputPortSignal(S,2);
    real_T       *Iabc                = ssGetOutputPortSignal(S,3);
    real_T       *Err                 = ssGetOutputPortSignal(S,4);
 //==================================================================================   
  
    PARK( Iabc1[0] , Iabc1[1] , Iabc1[2] ,Iabc1[3] ,Iabc1[4] , Theta_s ,&Id, &Iq) ;
    
    
    nref     = 60 * Wm_Ref1[0] / ( 2 * P * PI ) ;
    Phir_ref = Phir_d_ref[0] ;
    n        = 60 * Wm1[0] / ( 2 * P * PI ) ;
    //n    =  Wm1[0];
   //====================================   
    k       = k + 1 ; 
    t       = t + h;
    Te      = 2 * h;
    n1      = k / 20 ;
    n2      = round( k / 20 ) ;
    
    nn1     = k / 2 ;
    nn2     = round( k / 2 ) ;
    
    
    
//====================================  
    if ( n1est == n2est )
    {   
   Phir_d_estim = Phir_d_estim + Te * ( Lm * Id + Tr * Wsg * Phir_q_estim - Phir_d_estim ) / Tr ;
   Phir_q_estim = Phir_q_estim + Te * ( Lm * Iq - Tr * Wsg * Phir_d_estim - Phir_q_estim ) / Tr ;
   Phir_estim  = sqrt( Phir_d_estim * Phir_d_estim + Phir_q_estim * Phir_q_estim ) ;
   //Phir_estim0 = Phir_estim;
    }
    
    Wsg     =( Tr * Phir_estim ) / ( Lm * Iq )     ; 
    Ws      = Wsg-Wm1[0] ;
    
    
    if ( Ws =0 ){Ws      = 1 ;}
    Theta_s = Theta_s + h * Wm_Ref1[0];
//===========================================================    
    
    if ( n1 == n2 )
    {   
     ew     = ( nref - n ) / fabs( nref ) ;
     ve     = ( ew - e0 ) * 2.75 ;
      FUZ( ew ,&NB1,&NM1,&NS1,&Z1,&PS1,&PM1,&PB1) ;
      
      
      faper1[0] = NB1 ;
      faper1[1] = NM1 ;
      faper1[2] = NS1 ;
      faper1[3] = Z1 ;
      faper1[4] = PS1 ;
      faper1[4] = PM1 ;
      faper1[6] = PB1 ;
           
      FUZ( ve ,&NB1,&NM1,&NS1,&Z1,&PS1,&PM1,&PB1) ;
      
      fapder1[0] = NB1 ;
      fapder1[1] = NM1 ;
      fapder1[2] = NS1 ;
      fapder1[3] = Z1 ;
      fapder1[4] = PS1 ;
      fapder1[4] = PM1 ;
      fapder1[6] = PB1 ;
     e0       = ew ;
    }
    
   alpha1 = GAINE( fapder1 , faper1 , 1 ) ;
   
  if ( fabs( Is_q_ref ) > Is_q_max )
        { Is_q_ref = Is_q_max * sign( Is_q_ref ) ;}
   else { Is_q_ref = FLO( Is_q_ref , faper1 , fapder1 , alpha1 , 20 ) ;} 
//=========================================================================================== 
   
  
   
 if ( nn1 == nn2 )
    {   
     Phir_Err     = ( Phir_ref - Phir_estim ) ;
     veph     = ( Phir_Err - ei ) *18 ;
      FUZ( Phir_Err ,&NB1,&NM1,&NS1,&Z1,&PS1,&PM1,&PB1) ;
      
      
      faper2[0] = NB1 ;
      faper2[1] = NM1 ;
      faper2[2] = NS1 ;
      faper2[3] = Z1 ;
      faper2[4] = PS1 ;
      faper2[4] = PM1 ;
      faper2[6] = PB1 ;
           
      FUZ( veph ,&NB1,&NM1,&NS1,&Z1,&PS1,&PM1,&PB1) ;
      
      fapder2[0] = NB1 ;
      fapder2[1] = NM1 ;
      fapder2[2] = NS1 ;
      fapder2[3] = Z1  ;
      fapder2[4] = PS1 ;
      fapder2[4] = PM1 ;
      fapder2[6] = PB1 ;
     ei      = Phir_Err ;
    }
    
   alpha2 = GAINE( fapder2 , faper2 , 8 ) ;  
   
   if ( fabs( Is_d_ref ) > Is_d_max )
              { Is_d_ref  = Is_d_max * sign( Is_d_ref ) ;}
       else   { Is_d_ref  = FLO( Is_d_ref , faper2 , fapder2 , alpha2 , 8 ) ; }
 //=================================================================================
 PARK_INV( Is_d_ref , Is_q_ref , Theta_s ,&Iabc[0], &Iabc[1] , &Iabc[2] , &Iabc[3] , &Iabc[4]);
 
                               //Theta_s
 
 //Iabc[0] = (650/14)* Iabc[0];
 //Iabc[1] = (650/14)* Iabc[1];
 //Iabc[2] = (650/14)* Iabc[2];
 Vs_d1[0] = Is_d_ref ;  
 Vs_q1[0] = Is_q_ref ; 
 Ws1[0]   = Ws ;
 Err[0]   = ew ;
}



#define MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
  /* Function: mdlUpdate ======================================================
   * Abstract:
   *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
  }
#endif /* MDL_UPDATE */



#define MDL_DERIVATIVES  /* Change to #undef to remove function */
#if defined(MDL_DERIVATIVES)
  /* Function: mdlDerivatives =================================================
   * Abstract:
   *    In this function, you compute the S-function block's derivatives.
   *    The derivatives are placed in the derivative vector, ssGetdX(S).
   */
  static void mdlDerivatives(SimStruct *S)
  {
  }
#endif /* MDL_DERIVATIVES */



/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
}


/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
