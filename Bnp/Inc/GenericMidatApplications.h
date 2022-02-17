#ifndef GENERICMIDATAPPLICATIONS
#define GENERICMIDATAPPLICATIONS

#include "GenericMidatCalib.h"
#include "GenericMidatPricing.h"
#include "GenericMidatUtil.h"
#include "srt_h_all.h"

typedef struct {
  int iCalibATM;
  int iCalibFwdVol;
  char cCorrelName[100];

  int iNbShortVolShift;
  double *dShortVolShiftTime;
  double *dShortVolShift;

  double dVolShift;
  int iVolShiftType;

  int iNewAlgo;

} GMA_MidatParams, *GMA_MIDATPARAMS;

Err GMA_Midat_allocate_params(int iNbShortVolShift, double *dShortVolShiftTime,
                              double *dShortVolShift, GMA_MIDATPARAMS sParams);

void GMA_MIDAT_set_default_params(GMA_MIDATPARAMS sParams);

void GMA_MIDAT_free_params(GMA_MIDATPARAMS sParams);

Err GMA_MidatAutocal(/* Cash Flows description */
                     long lToday, char *cYieldCurve, char *cVolCurve,

                     /*		Generic Info */
                     long lStartDate, long lTheoEndDate,

                     /* Fixed Leg */
                     int iNbFixedCoupon, double *dCoupon,
                     long *lCouponStartDate, long *lCouponEndDate,
                     long *lCouponPayDate, char *cCouponBasis,

                     /* Floating Leg */
                     int iNbFixedFunding, char *cRefRate,
                     double *dFundingMargin, long *lFundingStartDate,
                     long *lFundingEndDate, long *lFundingPayDate,
                     char *cFundingBasis,

                     /* Exercise Dates */
                     int iNbExercise, long *lExerciseDate,
                     long *lSettlementDate, double *lExerciseFee,

                     int iPayRec,

                     /* Model parameters */
                     GENMIDAT_MODEL sModel,

                     /* Parameters */
                     GENMIDAT_AUTOCALPARAMS sAutocalParams,
                     GENMIDAT_CALIBPARAMS sCalibParams,
                     GENMIDAT_PDEPAMS sPDEParams,

                     /* Extra Parameters */
                     GMA_MIDATPARAMS sMidatParams,

                     /* Outputs */
                     double *dIV, double *dCall, GENMIDAT_AUTOCALINFO sInfos);

Err GMA_MidatAutocalNew(/* Cash Flows description */
                        long lToday, char *cYieldCurve, char *cVolCurve,
                        int iEODFixFlag, int iEODPayFlag, int iEODExeFlag,

                        /*		Generic Info */
                        long lStartDate, long lTheoEndDate,

                        /* Floating Leg */
                        int iNbFunding, char *cRefRate, long *lFundingFixDate,
                        long *lFundingStartDate, long *lFundingPayDate,
                        double *dFundingCoverage, double *dFundingMargin,
                        double *dPastFixings,

                        /* Fixed Leg */
                        int iNbCoupon, long *lCouponStartDate,
                        long *lCouponPayDate, double *dCouponCoverage,
                        double *dCoupon,

                        /* Exercise Dates */
                        int iNbExercise, long *lExerciseDate,
                        long *lSettlementDate, double *lExerciseFee,

                        double dPayRec,

                        /* Model parameters */
                        GENMIDAT_MODEL sModel,

                        /* Parameters */
                        GENMIDAT_AUTOCALPARAMS sAutocalParams,
                        GENMIDAT_CALIBPARAMS sCalibParams,
                        GENMIDAT_PDEPAMS sPDEParams,

                        /* Extra Parameters */
                        GMA_MIDATPARAMS sMidatParams,

                        /* Outputs */
                        double *dIV, double *dCall,
                        GENMIDAT_AUTOCALINFO sInfos);

#endif