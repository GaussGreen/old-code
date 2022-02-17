/*	-------------------------------------------- */
/*	---- Continuous Volatility Model (CVM) ----- */
/*	-------------------------------------------- */

#ifndef __CVM_CALIB_H
#define __CVM_CALIB_H

#include "CVM.h"

char *
SrtCalibCVMUnd(char *undName,
               // Market ID
               char *ycname, char *vcname,
               // Swaption Details
               SrtCompounding srtFreq, SrtBasisCode srtBasis, char *refRate,

               // Interpollation Type
               double shift, int converging_sliding, int interpollation_mode,

               // Pillars
               int nbOfPillars, char **pillars_tenors,

               // Pillars Correlations
               double ***pillars_correls,

               // Sensitivities
               int sensitype, double **lambdas,

               int nbOfSensiTenors, char **sensi_tenors, double **sensitivities,

               // Numerical Parameter
               int nQuadLegendre,

               // Calibration Parameters
               int nbOptionTenors, char **option_tenors, int nbUnderlyingTenors,
               char **underlying_tenors, double **market_vols,
               int allow_negativeVols);

#endif