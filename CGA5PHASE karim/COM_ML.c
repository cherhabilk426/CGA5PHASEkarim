/*
 * sfuntmpl_basic.c: Basic 'C' template for a level 2 S-function.
 *
 * Copyright 1990-2013 The MathWorks, Inc.
 */


/*
 * You must specify the S_FUNCTION_NAME as the name of your S-function
 * (i.e. replace sfuntmpl_basic with the name of your S-function).
 */

#define S_FUNCTION_NAME  COM_ML
#define S_FUNCTION_LEVEL 2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */

#include "stdio.h"
#include "string.h"
//#include "stdbool.h"  //bool

#include "stdlib.h"
#include "simstruc.h"
#include "math.h" 

#define PI 3.14159265358979323846264338327950288419716939937510582 
#define h  1e-6
real_T t , def_1 , def_2 , def_3, def_4, def_5,Sin_a,Sin_b,Sin_c,Sin_f,Sin_e,tri_1,f_tri; 
/* Error handling
 * --------------
 *
 * You should use the following technique to report errors encountered within
 * an S-function:
 *
 *       ssSetErrorStatus(S,"Error encountered due to ...");
 *       return;
 *
 * Note that the 2nd argument to ssSetErrorStatus must be persistent memory.
 * It cannot be a local variable. For example the following will cause
 * unpredictable errors:
 *
 *      mdlOutputs()
 *      {
 *         char msg[256];         {ILLEGAL: to fix use "static char msg[256];"}
 *         sprintf(msg,"Error due to %s", string);
 *         ssSetErrorStatus(S,msg);
 *         return;
 *      }
 *
 */

  
/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 2)) return;
    ssSetInputPortWidth(S, 0, 5);
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortRequiredContiguous(S, 0, true); /*direct input signal access*/
    ssSetInputPortRequiredContiguous(S, 1, true); /*direct input signal access*/
    /*
     * Set direct feedthrough flag (1=yes, 0=no).
     * A port has direct feedthrough if the input is used in either
     * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
     */
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);

    if (!ssSetNumOutputPorts(S, 16)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortWidth(S, 3, 1);
    ssSetOutputPortWidth(S, 4, 1);
    ssSetOutputPortWidth(S, 5, 1);
    ssSetOutputPortWidth(S, 6, 1);
    ssSetOutputPortWidth(S, 7, 1);
    ssSetOutputPortWidth(S, 8, 1);
    ssSetOutputPortWidth(S, 9, 1);
    ssSetOutputPortWidth(S, 10, 1);
    ssSetOutputPortWidth(S, 11, 1);
    ssSetOutputPortWidth(S, 12, 1);
    ssSetOutputPortWidth(S, 13, 1);
    ssSetOutputPortWidth(S, 14, 1);
    ssSetOutputPortWidth(S, 15, 1);

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
  {
      t=0;
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
   
    //entree
    const real_T *Sin_abcfe    = (const real_T*)    ssGetInputPortSignal(S,0);  // Vs_abcfe
    const real_T *Ws1          = (const real_T*)    ssGetInputPortSignal(S,1);  // Ws
    //sortie
    real_T       *sin_1    = ssGetOutputPortSignal(S,0);
    real_T       *sin_2    = ssGetOutputPortSignal(S,1);
    real_T       *sin_3    = ssGetOutputPortSignal(S,2);
    real_T       *sin_4    = ssGetOutputPortSignal(S,3);
    real_T       *sin_5    = ssGetOutputPortSignal(S,4);
    
    real_T       *tri_0    = ssGetOutputPortSignal(S,5);
    
    real_T       *T1       = ssGetOutputPortSignal(S,6);
    real_T       *T2       = ssGetOutputPortSignal(S,7);
    real_T       *T3       = ssGetOutputPortSignal(S,8);
    real_T       *T4       = ssGetOutputPortSignal(S,9);
    real_T       *T5       = ssGetOutputPortSignal(S,10);
    real_T       *T6       = ssGetOutputPortSignal(S,11);
    real_T       *T7       = ssGetOutputPortSignal(S,12);
    real_T       *T8       = ssGetOutputPortSignal(S,13);
    real_T       *T9       = ssGetOutputPortSignal(S,14);
    real_T       *T10      = ssGetOutputPortSignal(S,15);
  
    t = t + h ;
    
    f_tri = 15000;  // 15 K frequance tri
    
    Sin_a = 0.8 * Sin_abcfe[0];
    Sin_b = 0.8 * Sin_abcfe[1];
    Sin_c = 0.8 * Sin_abcfe[2];
    Sin_f = 0.8 * Sin_abcfe[3];
    Sin_e = 0.8 * Sin_abcfe[4];
    
    sin_1[0] = Sin_a ;
    sin_2[0] = Sin_b ;
    sin_3[0] = Sin_c ;
    sin_4[0] = Sin_f ;
    sin_5[0] = Sin_e ;
       
    tri_1    = 10 * -2.0 * asin( cos( 2 * PI * f_tri * t ) ) / PI ;
    tri_0[0] = tri_1 ;
       
    def_1 = Sin_a - tri_1 ; if ( def_1 >= 0 ) { T1[0] = 1 ; T2[0]  = 0 ; } else { T1[0] = 0 ; T2[0]  = 1 ; } ;
    def_2 = Sin_b - tri_1 ; if ( def_2 >= 0 ) { T3[0] = 1 ; T4[0]  = 0 ; } else { T3[0] = 0 ; T4[0]  = 1 ; } ;
    def_3 = Sin_c - tri_1 ; if ( def_3 >= 0 ) { T5[0] = 1 ; T6[0]  = 0 ; } else { T5[0] = 0 ; T6[0]  = 1 ; } ;
    def_4 = Sin_f - tri_1 ; if ( def_4 >= 0 ) { T7[0] = 1 ; T8[0]  = 0 ; } else { T7[0] = 0 ; T8[0]  = 1 ; } ;
    def_5 = Sin_e - tri_1 ; if ( def_5 >= 0 ) { T9[0] = 1 ; T10[0] = 0 ;} else { T9[0] = 0 ; T10[0]  = 1 ; } ;
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
