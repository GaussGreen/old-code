/* ===================================================================================
   FILENAME:      swp_h_swaption.h

   PURPOSE:       Compute a fwd swaption with the CMS methos compatible with the
   smile
   ===================================================================================
 */

#ifndef SWP_H_AMORTSWAPTION_H
#define SWP_H_AMORTSWAPTION_H

#include "swp_h_all.h"

Err AmortizedSwaptionShiftedLog(char *cYCname, char *cVCname, char *cRefRname,

                                long xlStartDate, long xlEndDate,

                                SrtCompounding srtFreq, SrtBasisCode srtBasis,

                                double exer_fee,

                                long lNFixNot, double *dFixNotionals,
                                double *dFixRates,

                                long lNFloatNot, double *dFloatNotionals,
                                double *dMargins,

                                SrtCallPutType srtCallPut,

                                int xlStudDegree, int xlNumSim,
                                double xlUpBound, int xlnClasses,

                                int UseVol,

                                double **Correl, double *dPrice);

Err AmortizedSwaptionShiftedLogForMAD_NotOpt(
    char *cYCname, char *cVCname, char *cRefRname,

    SrtCompounding srtFreq, SrtBasisCode srtBasis,

    double exer_fee,

    long *lFixStartDates, long *lFixEndDates, long lNFixDates,
    double *dFixNotionals, double *dFixRates,

    long *lFloatStartDates, long *lFloatEndDates, long lNFloatDates,
    double *dFloatNotionals, double *dFloatMargins, double *dFloatSpreads,

    SrtCallPutType srtCallPut,

    int xlNumSim,

    int iCalibShift, double Shift, int UseVol,

    double **Correl, double *dPrice);

Err AmortizedSwaptionShiftedLogForMAD(
    char *cYCname, char *cVCname, char *cRefRname,

    SrtCompounding srtFreq, SrtBasisCode srtBasis,

    double exer_fee,

    long *lFixStartDates, long *lFixEndDates, long lNFixDates,
    double *dFixNotionals, double *dFixRates,

    long *lFloatStartDates, long *lFloatEndDates, long lNFloatDates,
    double *dFloatNotionals, double *dFloatMargins, double *dFloatSpreads,

    SrtCallPutType srtCallPut,

    int xlNumSim,

    int iCalibShift, double Shift, int UseVol,

    double **Correl, double *dPrice);

Err AmortizedSwaptionShiftedLogForMAD_(
    char *cYCname, char *cVCname, char *cRefRname,

    SrtCompounding srtFreq, SrtBasisCode srtBasis,

    double exer_fee,

    long *lFixStartDates, long *lFixEndDates, long lNFixDates,
    double *dFixNotionals, double *dFixRates,

    long *lFloatStartDates, long *lFloatEndDates, long lNFloatDates,
    double *dFloatNotionals, double *dFloatMargins, double *dFloatSpreads,

    SrtCallPutType srtCallPut,

    int xlNumSim,

    int iCalibShift, double Shift, int UseVol,

    double **Correl, double *dPrice);

//--------------------------------------------------------------------
//---------------I don't like to delete code... ----------------------
//--------------------------------------------------------------------

/*

Err		AmortizedSwaptionHeston(
                                                 char *cYCname  ,
                                                 char *cVCname  ,
                                                 char *cRefRname  ,

                                                 long xlStartDate  ,
                                                 long xlEndDate  ,

                                                 SrtCompounding srtFreq  ,
                                                 SrtBasisCode srtBasis  ,

                                                 double exer_fee  ,

                                                 long	lNFixNot  ,
                                                 double *dFixNotionals  ,
                                                 double *dFixRates  ,

                                                 long	lNFloatNot  ,
                                                 double *dFloatNotionals  ,
                                                 double *dMargins  ,

                                                 SrtCallPutType srtCallPut  ,

                                                 double xlnStdDev  ,
                                                 int xlStudDegree  ,
                                                 int xlNumSim  ,
                                                 int xlNPts  ,
                                                 double xlUpBound  ,
                                                 int xlnClasses  ,

                                                 int UseVol  ,

                                                 double **Correl  ,
                                                 double *dPrice
                                    );

*/

/*
Err		PriceAmortizedSwaptionHeston(
                                                 char *cYCname  ,
                                                 char *cVCname  ,
                                                 char *cRefRname  ,
                                                 long xlStartDate  ,
                                                 long xlEndDate  ,
                                                 SrtCompounding srtFreq  ,
                                                 SrtBasisCode srtBasis  ,
                                                 long	lNFixNot  ,
                                                 double *dFixNotionals  ,
                                                 long	lNFloatNot  ,
                                                 double *dFloatNotionals  ,
                                                 double *dFixRates  ,
                                                 SrtCallPutType srtCallPut  ,
                                                 double xlnStdDev  ,
                                                 int xlStudDegree  ,
                                                 int xlNumSim  ,
                                                 char *xlTypeVol  ,
                                                 int xlNPts  ,
                                                 double xlUpBound  ,
                                                 int xlnClasses  ,
                                                 double **Correl  ,
                                                 double *dPrice  ,
                                                 int iPayoffMethod  ,
                                                 int iDistribMethod
                                    );
*/

Err PriceAmortizedSwaptionSABR(char *cYCname, char *cVCname, char *cRefRname,

                               long xlStartDate, long xlEndDate,

                               SrtCompounding srtFreq, SrtBasisCode srtBasis,

                               long lNFixNot, double *dFixNotionals,
                               long lNFloatNot, double *dFloatNotionals,
                               double *dFixRates,

                               SrtCallPutType srtCallPut,

                               int xlStudDegree, int xlNumSim, double xlUpBound,
                               int xlnClasses, double **Correl, double *dPrice,
                               int iPayoffMethod, int iDistribMethod);

Err PriceAmortizedSwaptionShiftedLog(
    char *cYCname, char *cVCname, char *cRefRname,

    long xlStartDate, long xlEndDate,

    SrtCompounding srtFreq, SrtBasisCode srtBasis,

    long lNFixNot, double *dFixNotionals, long lNFloatNot,
    double *dFloatNotionals, double *dFixRates,

    SrtCallPutType srtCallPut,

    int xlStudDegree, int xlNumSim, double xlUpBound, int xlnClasses,
    double **Correl, double *dPrice, int iPayoffMethod, int iDistribMethod);

#endif