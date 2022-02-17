/*********************************************************************************
 * Transition matrix for defaults
 *
 ********************************************************************************/

#ifndef _SMILE_H_
#define _SMILE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
    int                 totalName;        /* total num of names                 */
    int                 numStrike;        /* num of strikes                     */
    double              prob;             /* protfolio prob of default          */
    
    double              strike[MAX_NB];   /* strike                             */
    double              beta[MAX_NB];     /* base corr                          */
    
    double              normStrike[MAX_NB]; /* strike                           */    
    double              A[MAX_NB];        /* A coefs                            */
    double              C;                /* normalization                      */
    double              p;                /* beta parameter                     */
    double              rho;              /* rho */
    double              cdf[MAX_NB];      /* cdf                                */
    
    int                 calib;            /* 1 = calibrated, 0 = not calib'ed   */
    
    double              K1;               /* aux parms                          */
    double              K2;               /* aux parms                          */
}BCSmile;


/*********************************************************************************
 * Base Corr Smile Calibration
 *
 ********************************************************************************/
int      BCSmileCalib(
    BCSmile          *bcSmile);           /* (I/O) x                            */
    
/*********************************************************************************
 * Base Corr Smile Interoplation
 *
 ********************************************************************************/
int      BCSmileInterp(
    double           *res,                /* (O) interpolated base corr         */
    double           strike,              /* (I) strike                         */
    BCSmile          *bc,                 /* (I) smile structure                */
    char             flag);               /* (I) 'B': return beta, 'C':  CDF
                                             'P': PDF                           */

/*********************************************************************************
 * Test Func
 *
 ********************************************************************************/
int     CMTest(
    double          *res,                 /* (O) res                            */
    int             n,                    /* (I) n                              */
    int             N,                    /* (I) N                              */    
    double          rho,                  /* (I) rho                            */
    double          prob,                 /* (I) prob                           */
    char            flag);                /* (I) flag                           */    
    
#ifdef __cplusplus
}
#endif

#endif
    
