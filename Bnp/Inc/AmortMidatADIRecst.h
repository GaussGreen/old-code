// prevent multiple inclusions
#pragma once

//////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                // in the debug information"
// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

//#include "AmortMidatADIUtils.h"
#include "AmortMidatADIPrice.h"

double _DF_nojump_(double dX1, double dX3, double dDF_0_t, double dDF_0_T,
                   double dLAM1_0_t, double dLAM2_0_t, double dLAM1_0_T,
                   double dLAM2_0_T, double dEILam1_t, double dEILam2_t,
                   double dPHI1_t, double dPHI2_t, double dPHICross_t,
                   const _PCQ_Seq *pNumeraire, double dRho, double dAlphaRho,
                   double dAlpha_Sqrt1mRhoSqrd);

void _df_constant(double dDF_0_t, double dDF_0_T, double dLAM1_0_t,
                  double dLAM2_0_t, double dLAM1_0_T, double dLAM2_0_T,
                  double dEILam1_t, double dEILam2_t, double dPHI1_t,
                  double dPHI2_t, double dPHICross_t, double dExp_R1,
                  double dExp_R2, double dRho,
                  // results
                  double *pdConstant, double *pdCoeff1, double *pdCoeff2);

double _DF_(double dX1, double dX3, double dDF_0_t, double dDF_0_T,
            double dLAM1_0_t, double dLAM2_0_t, double dLAM1_0_T,
            double dLAM2_0_T, double dEILam1_t, double dEILam2_t,
            double dPHI1_t, double dPHI2_t, double dPHICross_t, double dExp_R1,
            double dExp_R2, double dRho, double dAlphaRho,
            double dAlpha_Sqrt1mRhoSqrd);

void _fill_ExpR(const double *pdLAM1_0_T_Begin, const double *pdLAM1_0_T_End,
                const double *pdLAM2_0_T_Begin, const double *pdEILam1_T_Begin,
                const double *pdEILam2_T_Begin, const double *pdPHI1_T_Begin,
                const double *pdPHI2_T_Begin, const double *pdPHICross_T_Begin,
                double dRho, const _PCQ_Seq *pNumeraire,
                // results
                double *pdExp_R1_Begin, double *pdExp_R2_Begin);

void _df_singleslice(
    // grid info
    const double *pdX1_Begin, const double *pdX1_End, const double *pdX3_Begin,
    const double *pdX3_End, double dDF_0_t, double dLAM1_0_t, double dLAM2_0_t,
    double dExpIntglam1_t, double dExpIntglam2_t, double dPHI1_t,
    double dPHI2_t, double dPHICross_t, double dDF_0_T, double dLAM1_0_T,
    double dLAM2_0_T, double dExp_R1, double dExp_R2, double dRho,
    double dAlphaRho, double dAlpha_Sqrt1mRhoSqrd,
    /// result
    double **ppdDF, int nR_Begin, int nC_Begin);
