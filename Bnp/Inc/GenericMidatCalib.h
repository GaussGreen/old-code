#ifndef GENERICMIDATCALIB
#define GENERICMIDATCALIB

#include "DiagCalibGen.h"
#include "GenericMidatUtil.h"
#include "srt_h_all.h"

typedef struct
{
    int    iCalibLambda;
    int    iShortIsCap;
    double dMinTime;
    int    iSkipLast;

    double dBetaStep;
    int    iContinueOnFail;
    int    iMaxBetaIterations;
    int    iMaxBadIter;
    int    iMaxReallyBadIter;
    int    iMaxSameDir;

    double dMinChangeVol;
    double dMaxChangeVol;

    double dLastBeta2;
    double dLastBeta2Tolerance;

    int iInterpFailedVol;
    int iInterpVolType; /* 0: linear in Vol, 1: Linear in Vol / (TEnd - TStart) */

} genmidat_calibparams, *GENMIDAT_CALIBPARAMS;

void genmidat_set_calibparams_default(GENMIDAT_CALIBPARAMS params);

typedef struct
{
    int     iNbExe;
    double* dLongVols;
    double* dShortVols;

    int    iFailedOpt;
    double dFailedSigma2;

} genmidat_calibinfos, *GENMIDAT_CALIBINFOS;

Err genmidat_alloc_calibinfos(int iNbExe, GENMIDAT_CALIBINFOS sCalibInfos);

void genmidat_free_calibinfos(GENMIDAT_CALIBINFOS sCalibInfos);

Err generic_midat_calib(/* Input */
                        int     iNbExe,
                        double* dLongOption,
                        double* dShortOption,

                        /* Model */
                        GENMIDAT_MODEL sModel,

                        /* Parameters */
                        GENMIDAT_CALIBPARAMS sParams,

                        /* Optional Output */
                        GENMIDAT_CALIBINFOS sCalibInfos);

Err generic_midat_calib_one_factor(/* Input */
                                   int     iNbExe,
                                   int*    iCalibLong,
                                   double* dLongVols,
                                   int*    iCalibShort,
                                   double* dShortVols,

                                   /* Model */
                                   GENMIDAT_MODEL sModel,

                                   /* Parameters */
                                   GENMIDAT_CALIBPARAMS sParams,

                                   /* Optional Infos */
                                   GENMIDAT_CALIBINFOS sCalibInfos);

Err generic_midat_calib_two_factor_total(/* Input */
                                         int     iNbExe,
                                         int*    iCalibLong,
                                         double* dLongVols,
                                         int*    iCalibShort,
                                         double* dShortVols,

                                         /* Model */
                                         GENMIDAT_MODEL sModel,

                                         /* Parameters */
                                         GENMIDAT_CALIBPARAMS sParams,

                                         /* Optional Infos */
                                         GENMIDAT_CALIBINFOS sCalibInfos);

Err generic_midat_calib_two_factor_new(/* Input */
                                       int     iNbExe,
                                       int*    iCalibLong,
                                       double* dLongVols,
                                       int*    iCalibShort,
                                       double* dShortVols,

                                       /* Model */
                                       GENMIDAT_MODEL sModel,

                                       /* Parameters */
                                       GENMIDAT_CALIBPARAMS sParams,

                                       /* Optional Infos */
                                       GENMIDAT_CALIBINFOS sCalibInfos);

Err generic_midat_calib_two_factor(/* Input */
                                   int     iNbExe,
                                   int*    iCalibLong,
                                   double* dLongVols,
                                   int*    iCalibShort,
                                   double* dShortVols,

                                   /* Model */
                                   GENMIDAT_MODEL sModel,

                                   /* Parameters */
                                   GENMIDAT_CALIBPARAMS sParams,

                                   /* Optional Infos */
                                   GENMIDAT_CALIBINFOS sCalibInfos);

typedef struct
{
    GENMIDAT_MODEL sModel;
    double         dSigma2;
    double         dLambda;
    double         dLastBeta;
    double         dBeta;

    double dLastVarX;
    double dLastVarY;
    double dLastCoVar;

    double dVarX;
    double dVarY;
    double dCoVar;

} genmidat_calib_model, *GenMidatCalib_MODEL;

struct genmidat_calib_inst;
typedef struct genmidat_calib_inst genmidat_calib_inst;
typedef genmidat_calib_inst*       GenMidatCalib_INST;

struct genmidat_calib_inst
{
    int    iIndex;
    double dTime;
    double dDt;

    double dLongVariance;
    double dShortVol;
    double dExpFact1;     /* exp(-gamma Ti) */
    double dExpFact2;     /* exp(-2gamma Ti) */
    double dLastExpFact1; /* exp(-gamma Ti) */
    double dLastExpFact2; /* exp(-2gamma Ti) */

    double dFact0;
    double dFact1;
    double dFact2;

    GenMidatCalib_INST sLastInst;
};

typedef struct
{
    GenMidatCalib_INST* AllInst;

} genmidat_calib_params, *GenMidatCalib_PARAMS;

#endif