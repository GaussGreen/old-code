// prevent multiple inclusions
#pragma once

//////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                // in the debug information"
// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "AmortMidatADIUtils.h"

typedef struct _model_specs {
  const double *m_pdSigT_Begin;
  const double *m_pdSigT_End;
  const double *m_pdSig_Begin;

  const double *m_pdLamT_Begin;
  const double *m_pdLamT_End;
  const double *m_pdLam_Begin;

  const int m_nNumFactor;
  const double m_dAlpha;
  const double m_dGamma;
  const double m_dRho;
  const double m_dAlphaRho;
  const double m_dAlpha_Sqrt1mRhoSqrd;
  const double m_dGammaRho_Over_Sqrt1mRhoSqrd;

  // storage
  const double *const m_pdLocalT_Begin;
  const double *const m_pdLocalT_End;
  const double *const m_pdLocalLam_Begin;
  const double *const m_pdLocalSig_Begin;

} _Model;

//-------------------------------------------------------------------------------------------------------
//	Description:
//
//  returns
//  (1) LAMBDA(0      ,pdLamT_Begin[i]) in pdResult
//  (2) exp( - integration of lambda term structure(u)) from u = 0 to u =
//  pdLamT_Begin[i] in  pdEIL_Begin
//-------------------------------------------------------------------------------------------------------
void _LAMBDA_Zero(const double *pdLamT_Begin, const double *pdLamT_End,
                  const double *pdLam_Begin, double dGamma,
                  double *pdResult_Begin,
                  double *pdEIL_Begin // EIL a.k.a exp(integration of lambda)
);

void _fill_PHI(const double *pdT_Begin, // exercise time
               const double *pdT_End, const double *pdEIL1_Begin,
               const double *pdEIL2_Begin, const double *pdLamT_Begin,
               const double *pdLamT_End, const double *pdLam_Begin,
               const double *pdSigT_Begin, const double *pdSigT_End,
               const double *pdSig_Begin, double dAlpha, double dGamma,
               double *pdPHI1_Begin, double *pdPHI2_Begin,
               double *pdPHI12_Begin);

void _fill_LAMBDA(const double *pdT_Begin, // exercise time
                  const double *pdT_End, const double *pdLamT_Begin,
                  const double *pdLamT_End, const double *pdLam_Begin,
                  double dGamma,
                  double *pdResult_Begin, /// of size exercise time only!
                  double *pdEIL_Begin);

void _fill_LAMBDA(const double *pdT_Begin, // exercise time
                  const double *pdT_End, const double *pdLamT_Begin,
                  const double *pdLamT_End, const double *pdLam_Begin,
                  double dGamma,
                  double *pdResult_Begin, /// of size exercise time only!
                  double *pdEIL_Begin);

const char *_alloc_Model(int nSizeSig, int nSizeLam, _Model *pModel);

void _init_Model(const double *pdSigT_Begin, const double *pdSigT_End,
                 const double *pdSig_Begin, const double *pdLamT_Begin,
                 const double *pdLamT_End, const double *pdLam_Begin,
                 const double *pdAlpha, const double *pdGamma,
                 const double *pdRho, _Model *pModel);

void _free_Model(_Model *pModel);
