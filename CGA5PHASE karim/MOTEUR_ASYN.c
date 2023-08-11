#define  S_FUNCTION_NAME  MOTEUR_ASYN
#define  S_FUNCTION_LEVEL 2

#include "simstruc.h" 
#include "mex.h" 
#include "math.h" 
#include "stdlib.h"
#include "stdio.h"

#define  h  1e-6
#define  PI 3.14159265358979323846264338327950288419716939937510582

#define  NB_SORTIE 10 

#define  NB_PARAM 4 // [Rs Ls],[Rr Lr],Lm,[J F P] 

#define  STATOR_PARAM(S)     ssGetSFcnParam(S,0) // [Rs Ls]
#define  ROTOR_PARAM(S)      ssGetSFcnParam(S,1) // [Rr Lr] 
#define  Lm_PARAM(S)         ssGetSFcnParam(S,2) // Lm
#define  JFP_PARAM(S)        ssGetSFcnParam(S,3) // [J F P]
//========================================================== VARIABLE GLOBAL =================================================================================
real_T  Rs,Rr,Ls,Lr,Lm,J,F,P,fs ,Teta_s,
        Vs_d,Vs_q,Is_d,Is_q,Ir_d,Ir_q,D_Phis_d,D_Phis_q,D_Phir_d,D_Phir_q,Phis_d,Phis_q,Phir_d,Phir_q,
        D_Wm,Ce,Cr1,Ws,Wm,Wr,W,
        Sg,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4;
//============================================================================================================================================================

//====================== PARK ==========================================================================
void PARK( real_T Va , real_T Vb , real_T Vc, real_T Vf, real_T Ve , real_T Angle ,real_T* Vd, real_T* Vq)
{   
   real_T K0 = (2.0/5) ; // Constant de PARK
   
   *Vd =  K0 * ( Va * cos( Angle ) + Vb * cos( Angle + 2 * PI / 5) + Vc * cos( Angle + 4 * PI / 5 )+ Vf * cos( Angle + 6 * PI / 5 )+ Ve * cos( Angle + 8 * PI / 5 ) ) ; 
   *Vq = -K0 * ( Va * cos( Angle ) + Vb * cos( Angle + 2 * PI / 5) + Vc * cos( Angle + 4 * PI / 5 )+ Vf * cos( Angle + 6 * PI / 5 )+ Ve * cos( Angle + 8 * PI / 5 ) ) ; 
}
//=======================================================================================================

//====================== PARK_INV =======================================================================
void PARK_INV( real_T Vd , real_T Vq , real_T Angle ,real_T* Va, real_T* Vb , real_T* Vc , real_T* Vf , real_T* Ve )
{   
    *Va = Vd   * cos( Angle )              - Vq   * sin( Angle ) ;                // Va 
    *Vb = Vd   * cos( Angle + 2 * PI / 5 ) - Vq   * sin( Angle + 2 * PI / 5 ) ;   // Vb   
    *Vc = Vd   * cos( Angle + 4 * PI / 5 ) - Vq   * sin( Angle + 4 * PI / 5 ) ;   // Vc   
    *Vf = Vd   * cos( Angle + 6 * PI / 5 ) - Vq   * sin( Angle + 6 * PI / 5 ) ;   // Vf   
    *Ve = Vd   * cos( Angle + 8 * PI / 5 ) - Vq   * sin( Angle + 8 * PI / 5 ) ;   // Ve   
}
//=======================================================================================================
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, NB_PARAM);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 3)) return;
    ssSetInputPortWidth(S, 0, 5);
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortWidth(S, 2, 1);
    ssSetInputPortRequiredContiguous(S, 0, true); /*direct input signal access*/
    ssSetInputPortRequiredContiguous(S, 1, true); /*direct input signal access*/
    ssSetInputPortRequiredContiguous(S, 2, true); /*direct input signal access*/
    
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 1);

    if (!ssSetNumOutputPorts(S, NB_SORTIE)) return;
    ssSetOutputPortWidth(S, 0, 2);   // Vs_d , Vs_q
    ssSetOutputPortWidth(S, 1, 2);   // Is_d , Is_q
    ssSetOutputPortWidth(S, 2, 2);   // Ir_d , Ir_q
    ssSetOutputPortWidth(S, 3, 2);   // Phis_d , Phis_d
    ssSetOutputPortWidth(S, 4, 2);   // Phir_d , Phir_d
    ssSetOutputPortWidth(S, 5, 5);   // Is_a , Is_b , Is_c , Is_f , Is_e
    ssSetOutputPortWidth(S, 6, 5);   // Ir_a , Ir_b , Ir_c , Ir_f , Ir_e
    ssSetOutputPortWidth(S, 7, 5);   // Phis_a , Phis_b , Phis_c  , Phis_f, Phis_e 
    ssSetOutputPortWidth(S, 8, 5);   // Phir_a , Phir_b , Phir_c  , Phir_f , Phir_e
    ssSetOutputPortWidth(S, 9, 4);   // Ce , Wm , Wr ,fs
    
    
 
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


static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, h);
    ssSetOffsetTime(S, 0, 0.0);

}



#define MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)
 
static void mdlInitializeConditions(SimStruct *S)
  {
     const real_T       *STATOR  = mxGetPr(STATOR_PARAM(S)); 
     const real_T       *ROTOR   = mxGetPr(ROTOR_PARAM(S));      
     const real_T       *Lm1     = mxGetPr(Lm_PARAM(S));        
     const real_T       *JFP     = mxGetPr(JFP_PARAM(S));  
//===============================================================
Rs = STATOR[0];     //% Résistance du stator
Ls = STATOR[1];     //% H Inductance du stator

Rr = ROTOR[0];      //% Résistance du rotor
Lr = ROTOR[1] ;      //% H Inductance du rotor

Lm = *Lm1;          //% H Inductance Mutuelle

J =  JFP[0];        //% Kg.m2 Moment d’inertie
F =  JFP[1];        //% SI Coefficient de frottement
P =  JFP[2];        //% Nombre de paire de pôle.

Ls = Ls + Lm ;
Lr = Lr + Lm ;



Is_d  = 0  ;
Is_q  = 0  ;

Ir_d  = 0  ;
Ir_q  = 0  ;

Phis_d = 0 ;
Phis_q = 0 ;

Phir_d = 0 ;
Phir_q = 0 ;
D_Wm  = 0 ;

Teta_s = 0 ;
Wm = 1500 * 2 * PI / 60 ; // Demarage Avec 1500 tr/min ;


  }
#endif /* MDL_INITIALIZE_CONDITIONS */



#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
 
  static void mdlStart(SimStruct *S)
  {
  }
#endif /*  MDL_START */

  static void mdlOutputs(SimStruct *S, int_T tid)
{
    
//============================================================================================================================================================    
    const real_T     *Vs            = ( const real_T* )  ssGetInputPortSignal(S,0) ; // Vs_abcfe
    const real_T     *Ws1           = ( const real_T* )  ssGetInputPortSignal(S,1) ; // Ws
    const real_T     *Cr            = ( const real_T* )  ssGetInputPortSignal(S,2) ; // Cr
       
    real_T           *Y0            = ssGetOutputPortSignal(S,0) ;                   // Vs_d , Vs_q
    real_T           *Y1            = ssGetOutputPortSignal(S,1) ;                   // Is_d , Is_q
    real_T           *Y2            = ssGetOutputPortSignal(S,2) ;                   // Ir_d , Ir_q
    real_T           *Y3            = ssGetOutputPortSignal(S,3) ;                   // Phis_d , Phis_d
    real_T           *Y4            = ssGetOutputPortSignal(S,4) ;                   // Phir_d , Phir_d
    real_T           *Y5            = ssGetOutputPortSignal(S,5) ;                   // Is_a , Is_b , Is_c, Is_f, Is_e
    real_T           *Y6            = ssGetOutputPortSignal(S,6) ;                   // Ir_a , Ir_b , Ir_c, Ir_f, Ir_e
    real_T           *Y7            = ssGetOutputPortSignal(S,7) ;                   // Phis_a , Phis_b , Phis_c, Phis_f, Phis_e
    real_T           *Y8            = ssGetOutputPortSignal(S,8) ;                   // Phir_a , Phir_b , Phir_c ,Phir_f , Phir_e    
    real_T           *Y9            = ssGetOutputPortSignal(S,9) ;                   // Ce , Wm , Wr ,fs
//============================================================================================================================================================     
    Cr1 = Cr[0] ; 
    //fs  = 50 ;
    //Ws =  2 * PI * fs ; 
    Ws = Ws1[0] ;
    Teta_s = Teta_s + h * Ws ;
    //Teta_s = PI/2 ;
//=============================================== PARK =======================================================================================================        
    PARK( Vs[0] , Vs[1] , Vs[2], Vs[3], Vs[4] , Teta_s , &Vs_d , &Vs_q )  ;
//========================================== Eq Electrique  ==================================================================================================
   Sg = ( Ls * Lr ) - ( Lm * Lm ) ;       
   W  = Ws - Wr ;
   
   A1 = -Rs * Lr / Sg    ;A2 =  Ws                ;A3 =  Rs * Lm / Sg     ;A4 = 0 ;
   B1 = -Ws              ;B2 = -Rs * Lr / Sg      ;B3 =  0                ;B4 =  Rs * Lm / Sg ;
   C1 =  Rr * Lm / Sg    ;C2 =  0                 ;C3 = -Rr * Ls / Sg     ;C4 =  W     ;
   D1 =  0               ;D2 =  Rr * Lm / Sg      ;D3 = -W                ;D4 = -Rr * Ls / Sg ;
 //============================================================================================================================================================     
   D_Phis_d  = A1 * Phis_d + A2 * Phis_q + A3 * Phir_d + A4 * Phir_q + Vs_d ;          
   D_Phis_q  = B1 * Phis_d + B2 * Phis_q + B3 * Phir_d + B4 * Phir_q + Vs_q ;
   D_Phir_d  = C1 * Phis_d + C2 * Phis_q + C3 * Phir_d + C4 * Phir_q ;
   D_Phir_q  = D1 * Phis_d + D2 * Phis_q + D3 * Phir_d + D4 * Phir_q ;

   Phis_d = Phis_d + h * D_Phis_d ;                                                     // integral D_Phis_d
   Phis_q = Phis_q + h * D_Phis_q ;                                                     // integral D_Phis_q
   Phir_d = Phir_d + h * D_Phir_d ;                                                     // integral D_Phir_d
   Phir_q = Phir_q + h * D_Phir_q ;                                                     // integral D_Phir_q
  
   Is_d   =  ( Lr / Sg ) * Phis_d - ( Lm / Sg ) * Phir_d ;
   Is_q   =  ( Lr / Sg ) * Phis_q - ( Lm / Sg ) * Phir_q ;
   Ir_d   =  ( Ls / Sg ) * Phir_d - ( Lm / Sg ) * Phis_d ;
   Ir_q   =  ( Ls / Sg ) * Phir_q - ( Lm / Sg ) * Phis_q ;
 
//================================ Eq Mecanique ===============================================================================================================
     Ce     = ( P / 2 ) * ( Lm / Lr ) * ( Phir_d * Is_q - Phir_q * Is_d ) ;              // Electrical torque Ce 
     D_Wm  =  ( Ce - Cr1 -F * Wm ) / J ; 
     Wm    =  Wm + h * D_Wm ;                      //Wm * 60 / ( 2 * PI )  ;             // integral D_Wm 
     Wr     =  P  * Wm ;
     
//================================= SORTIE =====================================================================================================================    
    Y0[0] = Vs_d ;                                                                        // Vs_d
    Y0[1] = Vs_q ;                                                                        // Vs_q
    
    Y1[0] = Is_d ;                                                                        // Is_d
    Y1[1] = Is_q ;                                                                        // Is_q
      
    Y2[0] = Ir_d ;                                                                        // Ir_d
    Y2[1] = Ir_q ;                                                                        // Ir_q
    
    Y3[0] = Phis_d ;                                                                      // Phis_d
    Y3[1] = Phis_q ;                                                                      // Phis_q
        
    Y4[0] = Phir_d ;                                                                      // Phir_d
    Y4[1] = Phir_q ;                                                                      // Phir_q
     
    PARK_INV( Is_d   , Is_q   , Teta_s , &Y5[0] , &Y5[1] , &Y5[2], &Y5[3], &Y5[4] );                      // Is_abcfe
    PARK_INV( Ir_d   , Ir_q   , Teta_s , &Y6[0] , &Y6[1] , &Y6[2], &Y6[3], &Y6[4] );                      // Ir_abcfe
    PARK_INV( Phis_d , Phis_q , Teta_s , &Y7[0] , &Y7[1] , &Y7[2], &Y7[3], &Y7[4] );                      // Phis_abcfe
    PARK_INV( Phir_d , Phir_q , Teta_s , &Y8[0] , &Y8[1] , &Y8[2], &Y8[3], &Y8[4] );                      // Phir_abcfe 
        
    Y9[0] = Ce  ;                                                                         // Ce
    Y9[1] = Wm ;                                                                          // Wm
    Y9[2] = Wr  ;                                                                         // Wr
    Y9[3] = Ws/(2* PI)  ;                                                                 // fs
//================================================================================================================================================================     
}



#define MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
   static void mdlUpdate(SimStruct *S, int_T tid)
  {
        
       
  }
#endif /* MDL_UPDATE */



#define MDL_DERIVATIVES  /* Change to #undef to remove function */
#if defined(MDL_DERIVATIVES)
  static void mdlDerivatives(SimStruct *S)
  {
  }
#endif /* MDL_DERIVATIVES */

static void mdlTerminate(SimStruct *S)
{
}


#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
