#include "AmortMidatADIExProb.h"
#include "AmortMidatADI.h"
#include "AmortMidatADIContract.h"

//////////////////////////////////////////////////////////////////////
//	 Common parameters for exercise probability
//   (to be supplied to Price_Contract)
typedef struct _GenMidAt_ExProb_Common_ {
  // for backward method
  // m_nCompute_ExProb = 0 to compute probability of not exercising at all
  int m_nCompute_ExProb;

  const double *m_pdT_Begin;
  const double *m_pdT_End;
  const double *m_pdVarX1_Begin;
  const double *m_pdVarX3_Begin;
  const double ***m_pppdExIndicator_Begin;

  // result
  // exercise probabilities
  double *m_pdExProb_Begin;
  // total probabilities
  double *m_pdTotalProb_Begin;

  const int *m_pnExIndicator_RBegin;
  const int *m_pnExIndicator_REnd;
  const int *m_pnExIndicator_CBegin;
  const int *m_pnExIndicator_CEnd;

} _GenMidAt_ExProb_Comm;

//-------------------------------------------------------------------------------------------------------
//	Description: initializes the common parameters
//						needed for updating exercise
// probability
//
//	Returns :
//-------------------------------------------------------------------------------------------------------
static void _init_ExProb_Comm(const double *pdEx_Begin, const double *pdEx_End,
                              const double ***pppdExIndicator,
                              const int *pnExIndicator_RBegin,
                              const int *pnExIndicator_REnd,
                              const int *pnExIndicator_CBegin,
                              const int *pnExIndicator_CEnd,
                              int nCompute_ExProb,
                              _GenMidAt_ExProb_Comm *pComm) {
  pComm->m_nCompute_ExProb = nCompute_ExProb;
  pComm->m_pdT_Begin = pdEx_Begin;
  pComm->m_pdT_End = pdEx_End;
  pComm->m_pppdExIndicator_Begin = pppdExIndicator;

  pComm->m_pnExIndicator_RBegin = pnExIndicator_RBegin;
  pComm->m_pnExIndicator_REnd = pnExIndicator_REnd;
  pComm->m_pnExIndicator_CBegin = pnExIndicator_CBegin;
  pComm->m_pnExIndicator_CEnd = pnExIndicator_CEnd;
}

static void _Init_Density(const double *pdX1_Begin, const double *pdX1_End,
                          const double *pdX3_Begin, const double *pdX3_End,
                          const double dVarX1, const double dVarX3,
                          double **ppdDensity, int nX1_Begin, int nX3_Begin) {
  const double dPi = 3.1415926535897932;
  const double dC = 1. / sqrt(2. * dPi);
  const double dDeltaX1 = pdX1_Begin[1] - pdX1_Begin[0];
  const double dDeltaX3 =
      (pdX3_End - pdX3_Begin) > 1 ? pdX3_Begin[1] - pdX3_Begin[0] : 1.;

  const double dDeltaArea =
      dDeltaX1 * dDeltaX3 == 0. ? dDeltaX1 : dDeltaX1 * dDeltaX3;
  const int nX1_End = nX1_Begin + pdX1_End - pdX1_Begin;
  const int nX3_End = nX3_Begin + pdX3_End - pdX3_Begin;
  int nX1, nX3;

  double dSumDensity = 0.;

  _ASSERTE(dDeltaArea > 1.e-15);
  _ASSERTE(dDeltaX1 > 1.e-15);

  for (nX1 = nX1_Begin; nX1 < nX1_End; ++nX1) {
    const double dX1 = pdX1_Begin[nX1 - nX1_Begin];
    for (nX3 = nX3_Begin; nX3 < nX3_End; ++nX3) {
      const double dX3 = pdX3_Begin[nX3 - nX3_Begin];
      ppdDensity[nX1][nX3] =
          _NormDensity(dC, dX1, 0., dVarX1) * _NormDensity(dC, dX3, 0., dVarX3);
      dSumDensity += ppdDensity[nX1][nX3];
    }
  }

  dSumDensity *= dDeltaArea;
  _ASSERTE(dSumDensity != 0.);
  dSumDensity = 1. / dSumDensity;

  // resample to sum to probabilitities 1.
  for (nX1 = nX1_Begin; nX1 < nX1_End; ++nX1)
    for (nX3 = nX3_Begin; nX3 < nX3_End; ++nX3)
      ppdDensity[nX1][nX3] *= dSumDensity;

#ifdef _DEBUG
  {
    const double dSum_Prob =
        dDeltaArea *
        _Sum_dMatrix(ppdDensity, nX1_Begin, nX1_End, nX3_Begin, nX3_End);
    _ASSERTE(fabs(dSum_Prob - 1.) < 1.e-14);
  }
#endif
}

static void _EXPROB_BACKWARD_(const double *pdEx, // exercise time
                              const double *pdX1_Begin, const double *pdX1_End,
                              const double *pdX3_Begin, const double *pdX3_End,
                              double **ppdU, int nX1_Begin, int nX3_Begin,
                              void *pAuxArg) {
  //// cast
  _GenMidAt_ExProb_Comm *pComm = (_GenMidAt_ExProb_Comm *)pAuxArg;

  // alias
  const double ***pppdExIndicator = pComm->m_pppdExIndicator_Begin;
  const double *pdEx_Begin = pComm->m_pdT_Begin;
  const double *pdEx_End = pComm->m_pdT_End;
  const int nNumEx = pdEx_End - pdEx_Begin;

  const int nEx_I = (const double *)bsearch(pdEx, pdEx_Begin, nNumEx,
                                            sizeof(*pdEx), _dcompare) -
                    pdEx_Begin;
  const double **ppdExIndicator = pppdExIndicator[nEx_I];

  const int nX3_End = nX3_Begin + pdX3_End - pdX3_Begin;
  const int nX1_End = nX1_Begin + pdX1_End - pdX1_Begin;
  int nX1, nX3;

  // if it is the last exercise and we are
  // interested in the exercise probabilities (as opposed to not exercise at
  // all)
  if (nEx_I == nNumEx - 1) {
    for (nX3 = nX3_Begin; nX3 < nX3_End; ++nX3) {
      for (nX1 = nX1_Begin; nX1 < nX1_End; ++nX1) {
        _ASSERTE(ppdExIndicator[nX1][nX3] == 0. ||
                 ppdExIndicator[nX1][nX3] == 1.);
        ppdU[nX1][nX3] = pComm->m_nCompute_ExProb
                             ? ppdExIndicator[nX1][nX3]
                             : 1. - ppdExIndicator[nX1][nX3];
      }
    }
  }

  else {
    for (nX3 = nX3_Begin; nX3 < nX3_End; ++nX3) {
      for (nX1 = nX1_Begin; nX1 < nX1_End; ++nX1) {
        _ASSERTE(ppdExIndicator[nX1][nX3] == 0. ||
                 ppdExIndicator[nX1][nX3] == 1.);
        ppdU[nX1][nX3] *= (1. - ppdExIndicator[nX1][nX3]);
      }
    }
  }
}

static void _GENMIDAT_EXPROB_FORWARD_(
    // exercise time
    const double *pdT, const double *pdX1_Begin, const double *pdX1_End,
    const double *pdX3_Begin, const double *pdX3_End, double **ppdU,
    int nX1_Begin, int nX3_Begin, void *pAuxArg) {
  //// cast
  const _GenMidAt_ExProb_Comm *pComm = (_GenMidAt_ExProb_Comm *)pAuxArg;

  // #th of exercise
  const int nSize = pComm->m_pdT_End - pComm->m_pdT_Begin;
  const double dT = *pdT;
  const int nI = (const double *)bsearch(&dT, pComm->m_pdT_Begin, nSize,
                                         sizeof(dT), _dcompare) -
                 pComm->m_pdT_Begin;

  // exercise indicator function at ith exercise
  const double **ppdIndicator = pComm->m_pppdExIndicator_Begin[nI];
  // variance
  const double dVar_X1 = pComm->m_pdVarX1_Begin[nI];
  const double dVar_X3 = pComm->m_pdVarX3_Begin[nI];

  // results
  double *pdExProb = pComm->m_pdExProb_Begin + nI;

  const int nX1_End = nX1_Begin + (pdX1_End - pdX1_Begin);
  const int nX3_End = nX3_Begin + (pdX3_End - pdX3_Begin);
  int nX1, nX3;

  // if first exercise         , then populate the grid with initial
  // probabilities
  if (nI == 0)
    _Init_Density(pdX1_Begin, pdX1_End, pdX3_Begin, pdX3_End, dVar_X1, dVar_X3,
                  ppdU, nX1_Begin, nX3_Begin);

  *pdExProb = 0.;

  _ASSERTE(pComm->m_pnExIndicator_RBegin[nI] >= nX1_Begin &&
           nX1_End >= pComm->m_pnExIndicator_REnd[nI]);
  _ASSERTE(pComm->m_pnExIndicator_CBegin[nI] >= nX3_Begin &&
           nX3_End >= pComm->m_pnExIndicator_CEnd[nI]);

  for (nX1 = nX1_Begin; nX1 < nX1_End; ++nX1) {
    for (nX3 = nX3_Begin; nX3 < nX3_End; ++nX3) {
      const double dIndiator = ppdIndicator[nX1][nX3];
      const double dExDensity = ppdU[nX1][nX3];

      _ASSERTE(dIndiator == 0. || dIndiator == 1.);
      _ASSERTE(nX1 < pComm->m_pnExIndicator_RBegin[nI] ||
                       nX1 >= pComm->m_pnExIndicator_REnd[nI]
                   ? dIndiator == 0.
                   : 1);
      _ASSERTE(nX3 < pComm->m_pnExIndicator_CBegin[nI] ||
                       nX3 >= pComm->m_pnExIndicator_CEnd[nI]
                   ? dIndiator == 0.
                   : 1);

      // Sum up the probabilities - result
      *pdExProb += (dIndiator == 1.) ? dExDensity : 0.;

      // Update
      ppdU[nX1][nX3] = dExDensity * (1. - dIndiator);
    }
  }
}

static void _GENMIDAT_TOTALPROB_FORWARD_(
    // exercise time
    const double *pdT, const double *pdX1_Begin, const double *pdX1_End,
    const double *pdX3_Begin, const double *pdX3_End, double **ppdU,
    int nX1_Begin, int nX3_Begin, void *pAuxArg) {
  //// cast
  const _GenMidAt_ExProb_Comm *pComm = (_GenMidAt_ExProb_Comm *)pAuxArg;

  // #th of exercise
  const int nSize = pComm->m_pdT_End - pComm->m_pdT_Begin;
  const double dT = *pdT;
  const int nI = (const double *)bsearch(&dT, pComm->m_pdT_Begin, nSize,
                                         sizeof(dT), _dcompare) -
                 pComm->m_pdT_Begin;

  // results
  double *pdTotalProb = pComm->m_pdTotalProb_Begin + nI;

  const int nX1_End = nX1_Begin + (pdX1_End - pdX1_Begin);
  const int nX3_End = nX3_Begin + (pdX3_End - pdX3_Begin);

  // if first exercise         , then populate the grid with initial
  // probabilities
  if (nI == 0) {
    // variance
    const double dVar_X1 = pComm->m_pdVarX1_Begin[nI];
    const double dVar_X3 = pComm->m_pdVarX3_Begin[nI];
    _Init_Density(pdX1_Begin, pdX1_End, pdX3_Begin, pdX3_End, dVar_X1, dVar_X3,
                  ppdU, nX1_Begin, nX3_Begin);
  }

  *pdTotalProb = _Sum_dMatrix(ppdU, nX1_Begin, nX1_End, nX3_Begin, nX3_End);
}

const char *Compute_ExProb_Forward_(
    const double *pdT_Begin, const double *pdT_End, _Grid *pGrid,
    const int *pnSizeX1_Begin, const int *pnSizeX3_Begin,
    const double *pdVarX1_Begin, const double *pdVarX3_Begin,
    const double ***pppdExIndicator_Begin, const int *pnExIndicator_RBegin,
    const int *pnExIndicator_REnd, const int *pnExIndicator_CBegin,
    const int *pnExIndicator_CEnd, const _Model *pModel,
    // results
    double *pdExProb_Begin) {
  const int nSize = (pdT_End - pdT_Begin);
  double *pdTKnot_Begin = memcpy(_alloca((nSize + 1) * sizeof(double)),
                                 pdT_Begin, nSize * sizeof(double));
  int *pnSizeX1_Begin_ = memcpy(_alloca((nSize + 1) * sizeof(int)),
                                pnSizeX1_Begin, nSize * sizeof(int));
  int *pnSizeX3_Begin_ = memcpy(_alloca((nSize + 1) * sizeof(int)),
                                pnSizeX3_Begin, nSize * sizeof(int));
  double *pdTotalProb = _alloca(nSize * sizeof(double)),
         *pdExProb = pdExProb_Begin;
  const char *szErr = 0;

  _GenMidAt_ExProb_Comm Comm = {0,
                                pdT_Begin,
                                pdT_End,
                                pdVarX1_Begin,
                                pdVarX3_Begin,
                                pppdExIndicator_Begin,
                                pdExProb_Begin,
                                pdTotalProb,
                                pnExIndicator_RBegin,
                                pnExIndicator_REnd,
                                pnExIndicator_CBegin,
                                pnExIndicator_CEnd};

  pdTKnot_Begin[nSize] = pdT_Begin[nSize - 1];
  pnSizeX1_Begin_[nSize] = pnSizeX1_Begin[nSize - 1];
  pnSizeX3_Begin_[nSize] = pnSizeX3_Begin[nSize - 1];

  /// diffuse the exercise probabilities
  szErr = Price_Contract(
      // T knots        , excluding today(t=0) because we know
      // analytically the joint densities of X1 and X3
      pdTKnot_Begin, pdTKnot_Begin + nSize + 1, pGrid,
      // NB: plus 1 because moving forward between Ti and Ti+1        ,
      // we want to use the size determined by the variance at Ti+1
      pnSizeX1_Begin_ + 1, pnSizeX3_Begin_ + 1, _GENMIDAT_EXPROB_FORWARD_,
      &Comm, _SIGMA1_Sqrd_, _SIGMA3_Sqrd_, _MU1_FORWARD_, _MU3_FORWARD_,
      (void *)pModel,
      // "Dummy" needed for the call to Price_Contract to proceed
      (double *)_alloca(sizeof(double)));

  if (szErr)
    return szErr;

  /// diffuse the total probabilities (for resampling purpose)
  szErr = Price_Contract(
      // T knots        , excluding today(t=0) because we know
      // analytically the joint densities of X1 and X3
      pdTKnot_Begin, pdTKnot_Begin + nSize + 1, pGrid,
      // NB: plus 1 because moving forward between Ti and Ti+1        ,
      // we want to use the size determined by the variance at Ti+1
      pnSizeX1_Begin_ + 1, pnSizeX3_Begin_ + 1, _GENMIDAT_TOTALPROB_FORWARD_,
      &Comm, _SIGMA1_Sqrd_, _SIGMA3_Sqrd_, _MU1_FORWARD_, _MU3_FORWARD_,
      (void *)pModel,
      // "Dummy" needed for the call to Price_Contract to proceed
      (double *)_alloca(sizeof(double)));

  if (szErr)
    return szErr;

  // re-sample
  for (; pdT_Begin < pdT_End; ++pdT_Begin, ++pdExProb, ++pdTotalProb)
    *pdExProb /= *pdTotalProb;

  // probability of not exercising
  *pdExProb = 1. - _Sum_dVec(pdExProb_Begin, pdExProb_Begin + nSize);

  return 0;
}

static const char *_Compute_ExProb_Backward(
    const double *pdEx_Begin, const double *pdEx_End, _Grid *pGrid,
    const int *pnSizeX3_Begin, const int *pnSizeX1_Begin,
    const double ***pppdExIndicator_Begin, const int *pnExIndicator_RBegin,
    const int *pnExIndicator_REnd, const int *pnExIndicator_CBegin,
    const int *pnExIndicator_CEnd, const _Model *pModel, int nCompute_ExProb,
    // results
    double *pdExProb) {
  // initialize common parameter
  _GenMidAt_ExProb_Comm Comm = {0};

  // sizes of the grid at each exercise point
  const int nNumEx = pdEx_End - pdEx_Begin;
  double *pdT_Begin = memset(_alloca((nNumEx + 1) * sizeof(double)), 0,
                             (nNumEx + 1) * sizeof(double));
  const char *szErr = 0;

  _reverse((int *)pnSizeX3_Begin, nNumEx, sizeof(int));
  _reverse((int *)pnSizeX1_Begin, nNumEx, sizeof(int));

  memcpy(pdT_Begin + 1, pdEx_Begin, nNumEx * sizeof(double));
  _reverse(pdT_Begin, nNumEx + 1, sizeof(double));

  _init_ExProb_Comm(pdEx_Begin, pdEx_End, pppdExIndicator_Begin,
                    pnExIndicator_RBegin, pnExIndicator_REnd,
                    pnExIndicator_CBegin, pnExIndicator_CEnd, nCompute_ExProb,
                    &Comm);

  // result
  szErr = Price_Contract(pdT_Begin, pdT_Begin + nNumEx + 1, pGrid,
                         pnSizeX1_Begin, pnSizeX3_Begin, _EXPROB_BACKWARD_,
                         &Comm, _SIGMA1_Sqrd_, _SIGMA3_Sqrd_, _MU1_BACKWARD_,
                         _MU3_BACKWARD_, (void *)pModel, pdExProb);

  // restore
  _reverse((int *)pnSizeX3_Begin, nNumEx, sizeof(int));
  _reverse((int *)pnSizeX1_Begin, nNumEx, sizeof(int));

  return szErr;
}

const char *Compute_ExProb_Backward(
    const double *pdEx_Begin, const double *pdEx_End, _Grid *pGrid,
    const int *pnSizeX3_Begin, const int *pnSizeX1_Begin,
    const double ***pppdExIndicator_Begin, const int *pnExIndicator_RBegin,
    const int *pnExIndicator_REnd, const int *pnExIndicator_CBegin,
    const int *pnExIndicator_CEnd, const _Model *pModel,
    // results
    double *pdExProb) {
  const double *const pdT_Begin = pdEx_Begin;
  const double *pdT_End = pdEx_End;
  const int nNumT = pdT_End - pdT_Begin;
  double *pdProb = pdExProb + nNumT - 1;
  const char *szErr = 0;

  // compute the probability of not exercise at all
  if (szErr = _Compute_ExProb_Backward(
          pdT_Begin, pdT_End, pGrid, pnSizeX3_Begin, pnSizeX1_Begin,
          pppdExIndicator_Begin, pnExIndicator_RBegin, pnExIndicator_REnd,
          pnExIndicator_CBegin, pnExIndicator_CEnd, pModel,
          0, // Compute the probability of not exercising at all
          pdExProb + nNumT))
    return szErr;

  for (; pdT_End > pdT_Begin; --pdT_End, --pdProb) {
    // compute probability
    if (szErr = _Compute_ExProb_Backward(
            pdT_Begin, pdT_End, pGrid, pnSizeX3_Begin, pnSizeX1_Begin,
            pppdExIndicator_Begin, pnExIndicator_RBegin, pnExIndicator_REnd,
            pnExIndicator_CBegin, pnExIndicator_CEnd, pModel,
            1, // Compute exercise probability
            pdProb))
      return szErr;
  }

  return 0;
}

static void _Compute_ExProb_Forward(const _Grid *pGrid,
                                    const double **ppdExIndicator,
                                    // results
                                    double *pdExProb, double **ppdExDensity) {
  const int nR_Begin = pGrid->m_pdX3_Begin - pGrid->m_pdX3_Begin_;
  const int nR_End = pGrid->m_pdX3_End - pGrid->m_pdX3_Begin_;
  const int nC_Begin = pGrid->m_pdX1_Begin - pGrid->m_pdX1_Begin_;
  const int nC_End = pGrid->m_pdX1_End - pGrid->m_pdX1_Begin_;
  int nR, nC;

  *pdExProb = 0.;

  for (nR = nR_Begin; nR < nR_End; ++nR) {
    for (nC = nC_Begin; nC < nC_End; ++nC) {
      const double dExIndicator = ppdExIndicator[nR][nC];
      const double dExDensity = ppdExDensity[nR][nC];
      // const double dDensity = ppdDensity[nR][nC];

      _ASSERTE(dExIndicator == 0. || dExIndicator == 1.);

      // sum up the probability density for this exercise
      *pdExProb += (dExIndicator == 1.) ? dExDensity : 0.;

      // prepare for the next iteration
      ppdExDensity[nR][nC] = dExDensity * (1. - dExIndicator);
    }
  }
}

static void _Init_ExDensity(const double **ppdDensity,
                            const double **ppdIndicator, const _Grid *pGrid,
                            double **ppdExDensity) {
  const int nR_Begin = pGrid->m_pdX3_Begin - pGrid->m_pdX3_Begin_;
  const int nR_End = pGrid->m_pdX3_End - pGrid->m_pdX3_Begin_;
  const int nC_Begin = pGrid->m_pdX1_Begin - pGrid->m_pdX1_Begin_;
  const int nC_End = pGrid->m_pdX1_End - pGrid->m_pdX1_Begin_;
  int nR, nC;
  for (nR = nR_Begin; nR < nR_End; ++nR) {
    for (nC = nC_Begin; nC < nC_End; ++nC) {
      const double dIndicator = ppdIndicator[nR][nC];
      _ASSERTE(dIndicator == 1. || dIndicator == 0.);
      ppdExDensity[nR][nC] = (1. - dIndicator) * ppdDensity[nR][nC];
    }
  }
}

const char *
Compute_ExProb(const double *pdEx_Begin, const double *pdEx_End, _Grid *pGrid,
               const int *pnSizeX3_Begin, const int *pnSizeX1_Begin,
               const double *pdVarX3_Begin, const double *pdVarX1_Begin,
               const _Model *pModel, const double ***pppdExIndicator_Begin,
               const int *pnExIndicator_RBegin, const int *pnExIndicator_REnd,
               const int *pnExIndicator_CBegin, const int *pnExIndicator_CEnd,
               // results
               double *pdExProb, int nUse_Backward) {
  const char *szErr = 0;

  // use the backward methodology
  // to compute exercise probabilities
  if (nUse_Backward) {
    szErr = Compute_ExProb_Backward(
        pdEx_Begin, pdEx_End, pGrid, pnSizeX3_Begin, pnSizeX1_Begin,
        pppdExIndicator_Begin, pnExIndicator_RBegin, pnExIndicator_REnd,
        pnExIndicator_CBegin, pnExIndicator_CEnd, pModel, pdExProb);
  }
  // use the forward methodology
  // to compute exercise probabilities
  else {
    szErr = Compute_ExProb_Forward_(
        pdEx_Begin, pdEx_End, pGrid, pnSizeX1_Begin, pnSizeX3_Begin,
        pdVarX1_Begin, pdVarX3_Begin, pppdExIndicator_Begin,
        pnExIndicator_RBegin, pnExIndicator_REnd, pnExIndicator_CBegin,
        pnExIndicator_CEnd, pModel, pdExProb);
  }

  if (szErr)
    return szErr;

  return 0;
}
