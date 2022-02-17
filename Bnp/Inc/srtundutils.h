/* ===============================================================================

   FILENAME	:	SrtUndUtils.h

   PURPOSE:     Utility Functions to deal with underlying for the outside world

   ===============================================================================
 */

#ifndef SRTUNDUTILS_H
#define SRTUNDUTILS_H

SRT_Boolean SrtIsUnderlyingDefined(char *und_name);

SRT_Boolean SrtIsUnderlyingInterestRate(char *und_name);

SRT_Boolean SrtIsUnderlyingIrOneFactor(char *und_name);

SRT_Boolean SrtIsUnderlyingIrTwoFactor(char *und_name);

void SrtGetUnderlyingModelName(char *und_name, char *mdl_name);

void SrtGetUnderlyingYieldCurveName(char *und_name, char *yc_name);

void SrtGetUnderlyingCurrency(char *und_name, char *ccy_name);

long SrtGetUnderlyingTicker(char *und_name);

/* Allocate memory for arrays where the underlying TS is displayed */
Err SrtDisplayUndTermStruct(char *UndName,

                            double ***SigmaCurve, long *NumSigmaRows,
                            long *NumSigmaCols,

                            double ***TauCurve, long *NumTauRows,
                            long *NumTauCols,

                            double *Alpha, double *Gamma, double *Rho,

                            double *Beta, double *Omega,

                            double *VoVol);
#endif