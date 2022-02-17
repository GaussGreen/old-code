#ifndef GENERICMIDATAUTOCAL
#define GENERICMIDATAUTOCAL

#include "GenericMidatCalib.h"
#include "GenericMidatPricing.h"
#include "GenericMidatUtil.h"
#include "srt_h_all.h"

typedef struct {
  int iOneTimeCallIndex;
  int iOneTimeShortLong;
  int iCalcExeBound;

  long lFindBestOptim;

} genmidat_autocalparams, *GENMIDAT_AUTOCALPARAMS;

void genmidat_set_autocalparams_default(GENMIDAT_AUTOCALPARAMS sParams);

typedef struct {
  GENMIDAT_MODEL sModel;
  GENMIDAT_CALIBINFOS sCalibInfos;

  int iCalcExeBound;
  int iNbExe;
  double *dExeBound;
  double **dCoefLin;

} genmidat_autocalinfos, *GENMIDAT_AUTOCALINFO;

void genmidat_init_autocalinfos(GENMIDAT_AUTOCALINFO sInfos);

Err genmidat_allocate_autocalinfos(int iNbExe, GENMIDAT_AUTOCALPARAMS sParams,
                                   GENMIDAT_AUTOCALINFO sInfos);

void genmidat_free_autocalinfos(GENMIDAT_AUTOCALINFO sInfos);

typedef struct {
  int iIndex;
  int iEventIndex;
  double dForward;
  double dForward2;
  double dBeta1;
  double dBeta2;

  GENMIDAT_PDEPAMS sPdeParams;
  GENMIDAT_AUTOCALINFO sInfos;

} genmidat_payoffparams, *GENMIDAT_PAYOFFPARAMS;

Err GenericMidatAutocal(int iNbExe, double *dLongOptions, double *dShortOptions,
                        int *iIsCallDate,

                        /* Model */
                        GENMIDAT_MODEL sModel,

                        /* Parameters */
                        GENMIDAT_AUTOCALPARAMS sAutocalParams,
                        GENMIDAT_CALIBPARAMS sCalibParams,
                        GENMIDAT_PDEPAMS sPDEParams,

                        /* Outputs */
                        double *MultiCall, GENMIDAT_AUTOCALINFO sInfos);

#endif