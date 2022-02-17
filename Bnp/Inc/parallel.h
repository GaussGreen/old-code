/* ===========================================================================

         FILENAME:     parallel.h

     PURPOSE:      Contains structure and function declarations for use in the
                       parallel version of the library.

     C. Godart: 16/06/99

   ===========================================================================
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif
#include "srt_h_grfnparams.h"
#include "srt_h_types.h"
#include "swp_h_swapdp.h"

/* ----------------------------------------------------------------------- */
struct prl_context {
  char bLocalComp;
  char *szGroupName;
  int iMsgTag;
  int iNumCPUs;   /* Number of machines in the VM */
  int *piCPUs;    /* array of TID for the CPUs as set by the master */
  int bAggregate; /* For use in mrqcof function */
};
extern struct prl_context PRL_CONTEXT;

typedef struct INSTRUMENTS_STRUCT INSTRUMENTS_STRUCT;
struct INSTRUMENTS_STRUCT {
  long lClcnDate;
  long lNumInstruments;
  SwapDP *psSwapDp;
  double *pdStrike;
  double *pdBondStrike;
  SrtReceiverType *peRecPay;
  StructType *peProductType;
  String *pszRefRateCode;
  long lNumTenors;
  double *pdFraMaturities;
  SrtGrfnParam *psGrfnParams;
};

typedef struct PARAMS_STRUCT PARAMS_STRUCT;
struct PARAMS_STRUCT {
  long lNumData;
  double *pdMarketTargets; /* From [1] to [lNumData] */
  double *pdMarketWeights; /* From [1] to [lNumData] */
  long lNumParams;
  double **ppdParamBounds; /* [1] is minimum ; [2] is maximum */

  double **ppdSigmaValues;
  long lNumSigmas;
  long lNumSigmaCols;
  double **ppdTauValues;
  long lNumTaus;
  long lNumTauCols;

  SRT_Boolean bFreezeTau;
  SRT_Boolean bOneTau;
  double *dFixedTau;
  SRT_Boolean bFreezeBeta;
  SRT_Boolean bOneBeta;
  double dFixedBeta;
  SRT_Boolean bFreezeOmega;
  SRT_Boolean bOneOmega;
  double dFixedOmega;

  double dFixedAlpha;
  double dFixedGamma;
  double dFixedRho;

  SrtCalibType eCalibType;

  SrtMdlDim eModelDim;
  SrtModelType eModelType;
  SRT_Boolean bSmoothSigma; /* use smooth sigma criterion */

  char szUndName[64];
};

/* Initially mrqcof.c stuff */
void SendCalibData(const double *x, int ndata, const double *a, int ma);
Err CheckRcvCalibData(int ma, int *pnb_local_data, int **prcvd, double **y,
                      double ***derivs, int **local_indexes, double *pmax_time,
                      int reset);
void TagTimes(double max_time);

/* ======================================================================= */
/* srt_f_calibcore.c
 * Changed : Cyril Godart : 16/08/98
 * From static to extern...
 * Needed in the parallel version
 * set_static_params_for_calib
 * set_static_instr_for_calib
 *
 ========================================================================== */

Err set_static_params_for_calib(
    SrtCalibType eCalibType, long lNumData, double *pdMarketWeights,
    double *pdMarketTargets, long lNumParams, long lNumSigmas, long lNumTaus,
    double **ppdParamBounds, SRT_Boolean bFreezeTau, SRT_Boolean bOneTau,
    double **dFixedTau, SRT_Boolean bFreezeBeta, SRT_Boolean bOneBeta,
    double dFixedBeta, SRT_Boolean bFreezeOmega, SRT_Boolean bOneOmega,
    double dFixedOmega, double dFixedAlpha, double dFixedGamma,
    double dFixedRho, SrtModelType eModelType, SrtMdlDim eModelDim,
    SRT_Boolean bSmoothSigma, String szUndName);

Err set_static_instr_for_calib(long lClcnDate, long lNumInstruments,
                               SwapDP *psSwapDp, double *pdStrike,
                               double *pdBondStrike, SrtReceiverType *peRecPay,
                               StructType *peProductType,
                               String *pszRefRateCode, long lNumTenors,
                               double *pdFraMaturities,
                               SrtGrfnParam *psGrfnParams);
Err allocate_space_for_static_sigtau_for_calib(SrtMdlDim eModelDim,
                                               long lNumSigmas,
                                               double *pdSigmaDates,
                                               long lNumTaus,
                                               double *pdTauDates);

Err levenberg_calib_funcs(double instr_index, double pdOptParams[],
                          double *value, double deriv[], int lNumParams);

Err free_calib_memory(SrtMdlDim eModelDim, long lNumSigmas, long lNumTaus);

Err change_und_name(char *the_und);

Err aggregate_price_value(double *pdOptParams, double *pdMktPrice);
Err set_static_sprdsht_for_calib(String und_name, long today);
Err set_parallel_context(int bLocalComp, char *szGroupName, int iNumCPUs,
                         int *piCPUs, int Version);
void SendParams(SrtCalibType eCalibType, long lNumData, double *pdMarketWeights,
                double *pdMarketTargets, long lNumParams, long lNumSigmas,
                long lNumTaus, double **ppdParamBounds, SRT_Boolean bFreezeTau,
                SRT_Boolean bOneTau, double **dFixedTau,
                SRT_Boolean bFreezeBeta, SRT_Boolean bOneBeta,
                double dFixedBeta, SRT_Boolean bFreezeOmega,
                SRT_Boolean bOneOmega, double dFixedOmega, double dFixedAlpha,
                double dFixedGamma, double dFixedRho, SrtModelType eModelType,
                SrtMdlDim eModelDim, SRT_Boolean bSmoothSigma,
                String szUndName);

void SendInstruments(long lClcnDate, long lNumInstruments, SwapDP *psSwapDp,
                     double *pdStrike, double *pdBondStrike,
                     SrtReceiverType *peRecPay, StructType *peProductType,
                     String *pszRefRateCode, long lNumTenors,
                     double *pdFraMaturities, SrtGrfnParam *psGrfnParams);

void SendParamForAllocation(SrtMdlDim eModelDim, long lNumSigmas,
                            double *pdSigmaDates, long lNumTaus,
                            double *pdTauDates);

void SendParamToDeallocate(SrtMdlDim eModelDim, long lNumSigmas, long lNumTaus);

///////////////////////////////////////////////////////////////
//                   New Version
//                These functions are called in
//                 levenberg_aggregate_functions
//////////////////////////////////////////////////////////////////
void StoreOptParams(double *pdOptParams, int size, int tag);
Err ParallelComputation(int lNumData, int lNumParams, double *pdOptParams);
Err RemoteErrorRecovery();
void SendSprdsht(char *und_name, long today);
void GetCPUs(int *pnb_cpus, int **pcpus);

/* srtcalibrateall */
void RequestChangeUndName(char *the_und);

/* srtinitirund.c */
char *SendInitIRUndData(char *undName, char *ycName, char *model,
                        int volCrvRows, int volCrvCols, double **volCrvVals,
                        int tauCrvRows, int tauCrvCols, double **tauCrvVals,
                        double beta, double alpha, double gamma, double rho,
                        double vovol, double etaOrBeta2);

void shift_param(double *param, int sign);

#ifdef _DEBUG
void MyTRACE(char *);
void MyTRACE1(char *, double d);
void MyTRACE1n(char *, double d);
void MyTRACE2(char *, double d1, double d2);
void MyTRACE3(char *, double d1, double d2, double d3);
#endif
#ifdef __cplusplus
}
#endif
