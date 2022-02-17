#include "AmortMidatADIContract.h"
#include "AmortMidatADIGrid.h"
#include "AmortMidatADIModel.h"

//-------------------------------------------------------------------------------------------------------
//	Description: computes the local vol of X1 in LGM2F
//
//	NB:	local vol is a function of time only;sgnature required for the
//PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _SIGMA1_Sqrd_(const double *pdX1, const double *pdX3,
                               const double *pdT, void *pArg) {
  // cast
  const _Model *pModel = (const _Model *)pArg;
  const int n = pdT - pModel->m_pdLocalT_Begin;
  _ASSERTE(n >= 0 && n < (pModel->m_pdLocalT_End - pModel->m_pdLocalT_Begin));
  return pModel->m_pdLocalSig_Begin[n] * pModel->m_pdLocalSig_Begin[n];
  // unreferenced
  pdX1, pdX3;
}

//-------------------------------------------------------------------------------------------------------
//	Description: computes the local vol of X3 in LGM2F
//
//	NB:	local vol is a function of time only;sgnature required for the
//PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _SIGMA3_Sqrd_(const double *pdX1, const double *pdX3,
                               const double *pdT, void *pArg) {
  return _SIGMA1_Sqrd_(pdX1, pdX3, pdT, pArg);
}

//-------------------------------------------------------------------------------------------------------
//	Description: computes the drift of X1 in LGM2F in the backward equation
//
//	NB:	local vol is a function of time only;sgnature required for the
//PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _MU1_BACKWARD_(const double *pdX1, const double *pdX3,
                                const double *pdT, void *pArg) {
  const _Model *pModel = (const _Model *)pArg;
  const int n = pdT - pModel->m_pdLocalT_Begin;
  _ASSERTE(n >= 0 && n < (pModel->m_pdLocalT_End - pModel->m_pdLocalT_Begin));
  return -pModel->m_pdLocalLam_Begin[n] * (*pdX1);
  // not referenced
  pdX3;
}

//-------------------------------------------------------------------------------------------------------
//	Description: computes the drift of X1 in LGM2F in the forward equation
//
//	NB:	local vol is a function of time only;sgnature required for the
//PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _MU1_FORWARD_(const double *pdX1, const double *pdX3,
                               const double *pdT, void *pArg) {
  return -1. * _MU1_BACKWARD_(pdX1, pdX3, pdT, pArg);
}

//-------------------------------------------------------------------------------------------------------
//	Description: computes the drift of X3 in LGM2F in the backward equation
//
//	NB:	local vol is a function of time only;sgnature required for the
//PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _MU3_BACKWARD_(const double *pdX3, const double *pdX1,
                                const double *pdT, void *pArg) {
  // cast
  _Model *pModel = (_Model *)pArg;
  const int n = pdT - pModel->m_pdLocalT_Begin;
  double dLambda2 = 0.;

  _ASSERTE(n >= 0 && n < (pModel->m_pdLocalT_End - pModel->m_pdLocalT_Begin));

  // const double dLambda2 = pModel->m_pdLocalLam_Begin[n] + pModel->m_dGamma;
  dLambda2 = pModel->m_pdLocalLam_Begin[n] + pModel->m_dGamma;

  return -dLambda2 * (*pdX3) - (*pdX1) * pModel->m_dGammaRho_Over_Sqrt1mRhoSqrd;
}

//-------------------------------------------------------------------------------------------------------
//	Description: computes the drift of X3 in LGM2F in the forward equation
//
//	NB:	local vol is a function of time only;sgnature required for the
//PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _MU3_FORWARD_(const double *pdX3, const double *pdX1,
                               const double *pdT, void *pArg) {
  return -1. * _MU3_BACKWARD_(pdX3, pdX1, pdT, pArg);
}

//-------------------------------------------------------------------------------------------------------
//	Description: local helper for preparingthe common parameters
//						for the LGM drift and vol
//functors
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
static void _Prepare_CommPDE(const double *pdT_Begin, const double *pdT_End,
                             const double *pdLocalLam_Begin,
                             const double *pdLocalSig_Begin, void *pComm_PDE) {
  // unpack PDE parameters
  _Model *pComm = (_Model *)pComm_PDE;
  const int nNumT = pdT_End - pdT_Begin;

  *(const double **)(&pComm->m_pdLocalT_Begin) = pdT_Begin;
  *(const double **)(&pComm->m_pdLocalT_End) = pdT_End;
  *(const double **)(&pComm->m_pdLocalLam_Begin) =
      pdLocalLam_Begin; //_alloca(nNumT * sizeof(double));
  *(const double **)(&pComm->m_pdLocalSig_Begin) =
      pdLocalSig_Begin; //_alloca(nNumT * sizeof(double));

  _bmatch_ex(pComm->m_pdLamT_Begin, sizeof(*pComm->m_pdLamT_Begin),
             pComm->m_pdLam_Begin, sizeof(*pComm->m_pdLam_Begin),
             pComm->m_pdLamT_End - pComm->m_pdLamT_Begin,
             pComm->m_pdLocalT_Begin, nNumT,
             (double *)pComm->m_pdLocalLam_Begin);

  _bmatch_ex(pComm->m_pdSigT_Begin, sizeof(*pComm->m_pdSigT_Begin),
             pComm->m_pdSig_Begin, sizeof(*pComm->m_pdSig_Begin),
             pComm->m_pdSigT_End - pComm->m_pdSigT_Begin,
             pComm->m_pdLocalT_Begin, nNumT,
             (double *)pComm->m_pdLocalSig_Begin);
}

//-------------------------------------------------------------------------------------------------------
//	Description:  Solve for the "contract" U(t  ,x)
//						U solves U(t  ,x1  ,x3)=E[Integrand(T
//,X(T))|X1(t)=x  ,X3(t)=x3]
//
//						NB: U traverses from *pdT_Begin to
//*(pdT_End-1) 						To solve for the backward equation  , pdT should be decreasing
//and for the forward equation  , increasing.
//
//
//	Returns : U(0  ,0)
//-------------------------------------------------------------------------------------------------------
const char *Price_Contract(
    /// time direction  , n = pdTKnot_End- pdTKnot_Begin-1 intervals in total  ,
    /// i.e. pdTKnot_Begin[0]-pdTKnot_Begin[1]  , ...
    /// pdTKnot_Begin[n-1]-pdTKnot_Begin[n]
    const double *pdTKnot_Begin, const double *pdTKnot_End,
    // Grid
    _Grid *pGrid,
    // size of the Grid corresponding to the intervals
    const int *pnSizeX1_Begin, const int *pnSizeX3_Begin,
    // Functor for integrand  ,i.e. "payoff" function(t  ,X1(t)  ,X3(t))
    // IMPORTANT : see .h for specs of the payoff function
    PFN_INTEGRAND _INTEGRAND_,
    // Auxillary arguments to be passed to _INTEGRAND_
    void *pComm_INTD,
    // functor to update vol of X1
    // functor must accept the argument seqence (X1  ,X3  ,T)
    PFN_COEFF_3DPDE _S1_Sqrd_,
    // functor to update vol of X3
    // functor must accept the argument seqence (X3  ,X1  ,T)
    PFN_COEFF_3DPDE _S3_Sqrd_,
    // functor to update drift of X1
    // functor must accept the argument seqence (X1  ,X3  ,T)
    PFN_COEFF_3DPDE _M1_,
    // functor to update drift of X3
    // functor must accept the argument seqence (X3  ,X1  ,T)
    PFN_COEFF_3DPDE _M3_,
    // Auxillary arguments to be passed to _S1_Sqrd_  ,_S3_Sqrd_  ,_M1_  ,_M3_
    void *pComm_PDE,
    // result
    double *pdResult) {
  const int nSize_T = pdTKnot_End - pdTKnot_Begin - 1;

  if (nSize_T >= 1) {
    const char *szErr = "mem allocation failed for Price_Contract(...)!";
    int nX1_Begin, nX1_End, nX3_Begin, nX3_End;
    double **ppdU = 0;        /// Grid
    _ADI_Storage store = {0}; // PDE storage

    // max the size of the grid
    Max_Grid(pnSizeX1_Begin, pnSizeX1_Begin + nSize_T, pnSizeX3_Begin, pGrid,
             &nX1_Begin, &nX1_End, &nX3_Begin, &nX3_End);

    // alloca mem for the Grid
    if (szErr = _alloc_dmatrix(&ppdU, nX1_Begin, nX1_End, nX3_Begin, nX3_End))
      return szErr;

    // alloc mem for PDE solver's storage
    Init_ADIStorage(&store, nX1_End - nX1_Begin, nX3_End - nX3_Begin);

    if (store.m_pdata == 0 || store.m_pm == 0) {
      // free storage
      Free_ADIStorage(&store);
      // free matrix
      free_dmatrix(ppdU, nX1_Begin, nX1_End - 1, nX3_Begin,
                   nX3_End - 1); // temporarily
      return "Price_Contract(...): mem allocation failed!";
    }

    // Check ...
    _ASSERTE(*pdTKnot_Begin >= 0.);

    // Solve from on excise date to the next
    for (; pdTKnot_Begin < (pdTKnot_End - 1);
         ++pdTKnot_Begin, ++pnSizeX1_Begin, ++pnSizeX3_Begin) {
      const double dT_From = pdTKnot_Begin[0];
      const double dT_To = pdTKnot_Begin[1];

      // actual # of t intervals and delta_t used
      const int nNumTInterval = NumInterval(dT_From, dT_To, pGrid->m_dDeltaT);

      // NB: n intervals but n+1 points
      double *pdLocalT_Begin = _alloca((nNumTInterval + 1) * sizeof(double));
      const double *pdLocalT_End = pdLocalT_Begin + nNumTInterval + 1;
      double *pdLocalLam_Begin = _alloca((nNumTInterval + 1) * sizeof(double));
      double *pdLocalSig_Begin = _alloca((nNumTInterval + 1) * sizeof(double));

      // Populate the time direction
      // Populate(dT_From  ,dT_To  ,pdLocalT_Begin  ,pdLocalT_End);
      Discretize_1D(dT_From, dT_To, nNumTInterval, pdLocalT_Begin);

      // Step 1. Shift the Grid
      _Shift_X1_(*pnSizeX1_Begin, pGrid);
      _Shift_X3_(*pnSizeX3_Begin, pGrid);

      // Step 2. Update integrand at dT_From
      _INTEGRAND_(pdTKnot_Begin, pGrid->m_pdX1_Begin, pGrid->m_pdX1_End,
                  pGrid->m_pdX3_Begin, pGrid->m_pdX3_End, ppdU,
                  pGrid->m_pdX1_Begin - pGrid->m_pdX1_Begin_,
                  pGrid->m_pdX3_Begin - pGrid->m_pdX3_Begin_, pComm_INTD);

      // Step 3. Solve the expectation

      // prepare the PDE parameters for this iteration
      _Prepare_CommPDE(pdLocalT_Begin, pdLocalT_End, pdLocalLam_Begin,
                       pdLocalSig_Begin, pComm_PDE);

      // delegate ...
      /// dominant in X1
      _ADI_MultipleSlices(
          pdLocalT_Begin, pdLocalT_End, pGrid->m_pdX1_Begin, pGrid->m_pdX1_End,
          pGrid->m_pdX3_Begin, pGrid->m_pdX3_End, pComm_PDE, _S1_Sqrd_, _M1_,
          _S3_Sqrd_, _M3_, ppdU, pGrid->m_pdX1_Begin - pGrid->m_pdX1_Begin_,
          pGrid->m_pdX3_Begin - pGrid->m_pdX3_Begin_, &store);
    }

    _ASSERTE(pGrid->m_pdX3_Begin_[pGrid->m_nX3_Zero] == 0.);
    _ASSERTE(pGrid->m_pdX1_Begin_[pGrid->m_nX1_Zero] == 0.);
    _ASSERTE(pGrid->m_nX3_Zero ==
             (pGrid->m_pdX3_Begin - pGrid->m_pdX3_Begin_ +
              (pGrid->m_pdX3_End - pGrid->m_pdX3_Begin) / 2));
    _ASSERTE(pGrid->m_nX1_Zero ==
             (pGrid->m_pdX1_Begin - pGrid->m_pdX1_Begin_ +
              (pGrid->m_pdX1_End - pGrid->m_pdX1_Begin) / 2));

    *pdResult = ppdU[pGrid->m_nX1_Zero][pGrid->m_nX3_Zero];

    // free storage
    Free_ADIStorage(&store);
    // free matrix
    free_dmatrix(ppdU, nX1_Begin, nX1_End - 1, nX3_Begin,
                 nX3_End - 1); // temporarily
  }

  return 0;
}
