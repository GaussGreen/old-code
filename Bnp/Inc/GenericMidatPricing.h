#ifndef GENERICMIDATPRICING
#define GENERICMIDATPRICING

#include "GenericMidatUtil.h"
#include "MCEBOptimisation.h"
#include "srt_h_all.h"

typedef struct
{
    int iPricingMethod;

    /* PDE */
    int  iDiscType;
    long lNbTime;
    long lNbX;

    double dNbStd;
    double dTheta;

    int iJumpQuadra;

    int iDiffMethod;
    int iAdaptGrid;

    /* MC */
    long lNbPaths;

} genmidat_pdeparams, *GENMIDAT_PDEPAMS;

void genmidat_set_pdeparams_default(GENMIDAT_PDEPAMS sParams);

typedef struct
{
    int iDiffMethod;

    double dCoefX1;
    double dCoefX2;
    double dCoefY1;
    double dCoefY2;

    double dLastCoefX1;
    double dLastCoefX2;
    double dLastCoefY1;
    double dLastCoefY2;

} genmidat_recons, *GENMIDAT_RECONS;

Err GenericMidat_OneFactor_PDE(/* Time Discretisation */
                               long    lNbTime,
                               double* dTimes,

                               /* Model */
                               GENMIDAT_MODEL sModel,

                               /*	Product data */
                               void** sFuncParmTab,
                               int*   iEvalEvt,

                               /*	Payoff Function */
                               Err (*payoff_func)(/* Event */
                                                  double dTime,

                                                  /* Parameters */
                                                  void* sFuncParm,

                                                  /* Gride data	*/
                                                  long            lStartIndex,
                                                  long            lEndIndex,
                                                  double*         dXt,
                                                  GENMIDAT_RECONS sReconsParams,

                                                  /* Vector of results to be updated */
                                                  int      iNbProd,
                                                  double** dProdVal),

                               /* PDE Parameters */
                               GENMIDAT_PDEPAMS sParams,

                               /* Result */
                               int     iNbProd,
                               double* dProdVal);

Err GenericMidat_TwoFactor_PDE(/* Time Discretisation */
                               long    lNbTime,
                               double* dTimes,

                               /* Model */
                               GENMIDAT_MODEL sModel,

                               /*	Product data */
                               void** sFuncParmTab,
                               int*   iEvalEvt,

                               /*	Payoff Function */
                               Err (*payoff_func)(/* Event */
                                                  double dTime,

                                                  /* Parameters */
                                                  void* sFuncParm,

                                                  /* Gride data	*/
                                                  long lStartIndex1,
                                                  long lEndIndex1,
                                                  long lStartIndex2,
                                                  long lEndIndex2,

                                                  double*         dXt,
                                                  double*         dZt,
                                                  GENMIDAT_RECONS sReconsParams,

                                                  /* Vector of results to be updated */
                                                  int       iNbProd,
                                                  double*** dProdVal),

                               /* PDE Parameters */
                               GENMIDAT_PDEPAMS sParams,

                               /* Result */
                               int     iNbProd,
                               double* dProdVal);

Err GenericMidat_OneFactor_MC(/* Time Discretisation */
                              int     iNbEvent,
                              double* dTimes,

                              /* Model */
                              GENMIDAT_MODEL sModel,

                              /*	Product data */
                              void** sFuncParmTab,

                              /*	Payoff Function */
                              Err (*payoff_func)(/* Event */
                                                 double dTime,

                                                 /* Parameters */
                                                 void* sFuncParm,

                                                 /* Gride data	*/
                                                 long            lNumPaths,
                                                 double*         dXt,
                                                 GENMIDAT_RECONS sReconsParams,

                                                 /* Vector of results to be updated */
                                                 int      iNbProd,
                                                 double** dProdVal),

                              /* MC Parameters */
                              GENMIDAT_PDEPAMS sParams,

                              /* for Optimisation of exercise boundary */
                              int        iDoOptim,
                              int*       iOptimise,
                              MCEBPARAMS sMCEBParams,

                              /* Result */
                              int      iNbProd,
                              double** dProdVal);

Err GenericMidat_TwoFactor_MC(/* Time Discretisation */
                              int     iNbEvent,
                              double* dTimes,

                              /* Model */
                              GENMIDAT_MODEL sModel,

                              /*	Product data */
                              void** sFuncParmTab,

                              /*	Payoff Function */
                              Err (*payoff_func)(/* Event */
                                                 double dTime,

                                                 /* Parameters */
                                                 void* sFuncParm,

                                                 /* Gride data	*/
                                                 long            lNumPaths,
                                                 double*         dXt,
                                                 double*         dYt,
                                                 GENMIDAT_RECONS sReconsParams,

                                                 /* Vector of results to be updated */
                                                 int      iNbProd,
                                                 double** dProdVal),

                              /* MC Parameters */
                              GENMIDAT_PDEPAMS sParams,

                              /* for Optimisation of exercise boundary */
                              int        iDoOptim,
                              int*       iOptimise,
                              MCEBPARAMS sMCEBParams,

                              /* Result */
                              int      iNbProd,
                              double** dProdVal);

#endif