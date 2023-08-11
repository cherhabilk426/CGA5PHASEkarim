/*
 * sfuntmpl_basic.c: Basic 'C' template for a level 2 S-function.
 *
 * Copyright 1990-2013 The MathWorks, Inc.
 */


/*
 * You must specify the S_FUNCTION_NAME as the name of your S-function
 * (i.e. replace sfuntmpl_basic with the name of your S-function).
 */

#define S_FUNCTION_NAME  CVV
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






real_T  Rs,Rr,Ls,Lr,Lm,J,F,P   ,Ts,Tr,Sg,Wm,Wm_ref,
        Wr_ref,Phir_d,Wr,S_1,P_I,S_1,Ce,Def_Wm  ,Is_d,Is_q,Ws,Vs_d,Vs_q,Teta_s;

real_T Psi,Wn,Ki,Kp;


//====================== PARK ==========================================================================
void PARK( real_T Va , real_T Vb , real_T Vc, real_T Vf, real_T Ve , real_T Angle , real_T* Vd , real_T* Vq)
{   
   real_T K0 = (2.0/5) ; // Constant de PARK
   *Vd =  K0 * ( Va * cos( Angle ) + Vb * cos( Angle + 2 * PI / 5 ) + Vc * cos( Angle + 4 * PI / 5 )+ Vf * cos( Angle + 6  * PI / 5 )+ Ve * cos( Angle + 8 * PI / 5 ) ) ; 
   *Vq = -K0 * ( Va * sin( Angle ) + Vb * sin( Angle + 2 * PI / 5 ) + Vc * sin( Angle + 4 * PI / 5 )+ Vf * sin( Angle + 6  * PI / 5 )+ Ve * sin( Angle + 8 * PI / 5 ) ) ;  
}
//=======================================================================================================

//====================== PARK_INV =======================================================================
void PARK_INV( real_T Vd , real_T Vq , real_T Angle ,real_T* Va, real_T* Vb , real_T* Vc, real_T* Vf, real_T* Ve)
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
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortWidth(S, 2, 1);
    

    ssSetInputPortRequiredContiguous(S, 0, true); /*direct input signal access*/
    ssSetInputPortRequiredContiguous(S, 1, true); /*direct input signal access*/
    ssSetInputPortRequiredContiguous(S, 2, true); /*direct input signal access*/
   

    /*
     * Set direct feedthrough flag (1=yes, 0=no).
     * A port has direct feedthrough if the input is used in either
     * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
     */
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    


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
  /* Function: mdlInitializeConditions ========================================
   * Abstract:
   *    In this function, you should initialize the continuous and discrete
   *    states for your S-function block.  The initial states are placed
   *    in the state vector, ssGetContStates(S) or ssGetRealDiscStates(S).
   *    You can also perform any other initialization activities that your
   *    S-function may require. Note, this routine will be called at the
   *    start of simulation and if it is present in an enabled subsystem
   *    configured to reset states, it will be call when the enabled subsystem
   *    restarts execution to reset the states.
   */
  static void mdlInitializeConditions(SimStruct *S)
  {  const real_T       *STATOR  = mxGetPr(STATOR_PARAM(S)); 
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


/*
   Psi    Wn
   0.4    7.70
   0.5    5.30
   0.6    5.20
   0.7    3.00
   1.0    4.75
*/
Psi = 1.00 ; Wn = 4.75;

Ki =  J * Wn * Wn ;
Kp = (2 * Psi * Ki / Wn ) - F ;

S_1 = 0 ;

//=====================
Sg = 1 - ( Lm * Lm / ( Ls * Lr ) );
Tr = Lr / Rr ; 
Ts = Ls / Rs ;
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
    const real_T     *Phir_d_ref      = ( const real_T* )  ssGetInputPortSignal(S,1);
    const real_T     *Wm1            = ( const real_T* )  ssGetInputPortSignal(S,2);
    
    real_T       *Vs_d1               = ssGetOutputPortSignal(S,0);
    real_T       *Vs_q1               = ssGetOutputPortSignal(S,1);
    real_T       *Ws1                 = ssGetOutputPortSignal(S,2);
    real_T       *Vabcfe              = ssGetOutputPortSignal(S,3);
    
    real_T       *Err                 = ssGetOutputPortSignal(S,4);
    
//==================================================================================           
    Wm_ref      = Wm_Ref1[0] ;                                    //Wm Ref
    Phir_d      = Phir_d_ref[0] ;                                 //Phir_d Ref
    Wm          = Wm1[0] ;                                        //Wm Sortie moteur
//=======================regulation PI = kp + ki/s ==============================
    Def_Wm = Wm_ref - Wm ;                                  //Wr_ref - Wr
    S_1 = S_1 + h * Def_Wm;
    Ce    = Kp * Def_Wm + Ki * S_1 ;   
//=================================================================================
   
    
    
   Is_d =  Phir_d / Lm   ;  
   Is_q =  Ce * Lr / ( P * Lm * Phir_d ) ; 
   
   Wr = Wm_ref * P ;
   Ws   = Lm * Is_q / ( Tr * Phir_d ) + Wr;  
   
   Vs_d = Rs * Is_d - Ws * Sg * Ls * Is_q;  
   Vs_q = Rs * Is_q + Ws * Sg * Ls * Is_d + Ws * Lm * Phir_d / Lr;  
     
      
   Teta_s = Teta_s + h * Ws ;
   
    Vs_d1[0]  = Vs_d ;
    Vs_q1[0]  = Vs_q ;
    
    Ws1[0]    = Ws   ;
    
    Err[0]    = Def_Wm ;
    
    PARK_INV( Vs_d , Vs_q , Teta_s ,&Vabcfe[0], &Vabcfe[1] , &Vabcfe[2], &Vabcfe[3], &Vabcfe[4] );

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
