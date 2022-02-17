#include "AmortMidatADIUtils.h"
#include "stdlib.h"
#include "utconst.h"

void _transform_vector(
    const void *pIn,
    unsigned nElem,    // # elements of input
    unsigned nSize_In, // size of elements in input (in bytes)
    unsigned
        nSize_Out, // size of element (in bytes) returned by the unary functor
    PFN_U pfnU,    // unary functor
    void *pComm,   // common function parameters
    void *pOut) {
  const void *pIn_end = (char *)pIn + nElem * nSize_In;
  for (; pIn != pIn_end; (char *)pIn += nSize_In, (char *)pOut += nSize_Out)
    pfnU(pIn, pComm, pOut);
}

//-------------------------------------------------------------------------------------------------------
//	Description: transform the matrix pOut by calling the binary function
//pfnBI
//
//	Returns : if(n1R)
//						pOut[n_RBegin+i][nC_Begin+j]=pfnBI(pIn1[i]
//,pIn2[j]  ,pComm) 					else 						pOut[n_RBegin+j][nC_Begin+i]=pfnBI(pIn1[i]  ,pIn2[j]
//,pComm)
//
//	NB:  pfnBI must take the order of the arguments pIn1 and pIn2
//-------------------------------------------------------------------------------------------------------
void _transform_matrix(
    const void *pIn1,
    unsigned nElem_In1, // # elements of input 1
    unsigned nSize_In1, // size of elements in input1(in bytes)
    const void *pIn2,
    unsigned nElem_In2, // # elements of input 2
    unsigned nSize_In2, // size of elements in input2(in bytes)
    unsigned
        nSize_Out, // size of elements (in bytes) returned by the binary functor
    PFN_BI pfnBI,  // binary functor
    void *pComm,   // common function parameters
    void **pOut,   // output-	pOut[nR_Begin  ,...][nC_Begin  ,...]
    unsigned nR_Begin, unsigned nC_Begin,
    int n1R // align 1 in dimension R
) {
  unsigned n1, n2;

  for (n1 = 0; n1 < nElem_In1; ++n1) {
    // input 1  ,2
    const void *p1 = (char *)pIn1 + n1 * nSize_In1;
    const char *p2 = pIn2;
    for (n2 = 0; n2 < nElem_In2; ++n2, p2 += nSize_In2) {
      // output
      const unsigned nR = n1R ? n1 + nR_Begin : n2 + nR_Begin;
      const unsigned nC = n1R ? n2 + nC_Begin : n1 + nC_Begin;
      pfnBI(p1, p2, pComm, (char *)pOut[nR] + nC * nSize_Out);
    }
  }
}

//-------------------------------------------------------------------------------------------------------
//	Description: binary functor that computes the product of 2 doubles
//
//	Returns : product
//-------------------------------------------------------------------------------------------------------
void __stdcall _MULTD_PFN_BI(const double *pIn1, const double *pIn2,
                             const void *pComm, double *pProduct) {
  *pProduct = (*pIn1) * (*pIn2);
  // not referenced
  pComm;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// function retuns 1 if dX equals machine precision  , 0 otherwise
int is_zero(double dX) {
  const double dZero = 1.e-15;
  return fabs(dX) > dZero ? 0 : 1;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
const void *lower_bound(const void *pvKey,   // pointer to key
                        const void *pvBegin, // pointer to base element
                        unsigned nElems,     // # elements
                        unsigned nSize,      // size of elements (in bytes)
                        PFN_COMPARE_BINARY pfnCompare // binary compare function
) {
  int nN = nElems;
  int nN2 = (nN >> 1);
  const void *pvMid = 0;
  const void *pvEnd = (char *)pvBegin + nN * nSize;
  for (; nN > 0; nN2 = (nN >> 1)) {
    pvMid = (char *)pvBegin + nN2 * nSize;
    if (pfnCompare(pvMid, pvKey) < 0) {
      pvBegin = (char *)pvMid + nSize;
      nN -= nN2 + 1;
    } else
      nN = nN2;
  }
  return pvBegin;
}

const void *upper_bound(const void *pvKey,   // pointer to key
                        const void *pvBegin, // pointer to base element
                        unsigned nElems,     // # elements
                        unsigned nSize,      // size of elements (in bytes)
                        PFN_COMPARE_BINARY pfnCompare // binary compare function
) {
  int nN = nElems;
  int nN2 = (nN >> 1);
  const void *pvMid = 0;
  const void *pvEnd = (char *)pvBegin + nN * nSize;
  for (; nN > 0; nN2 = (nN >> 1)) {
    pvMid = (char *)pvBegin + nN2 * nSize;
    if (pfnCompare(pvMid, pvKey) <= 0) {
      pvBegin = (char *)pvMid + nSize;
      nN -= nN2 + 1;
    } else
      nN = nN2;
  }
  return pvBegin;
}

///// returns min(n : pvBegin[n]  > *pvKey  , for 0 <= n <nElems )
///// if no such an index can be found  , return pvBegin+nElems*nSize
const void *
find_min_greater_than(const void *pvKey,   // pointer to key
                      const void *pvBegin, // pointer to base element
                      unsigned nElems,     // # elements
                      unsigned nSize,      // size of elements (in bytes)
                      PFN_COMPARE_BINARY pComp) {
  const void *pvResult = upper_bound(pvKey, pvBegin, nElems, nSize, pComp);

  // Check ...
#ifdef _DEBUG
  if (nElems > 0 && nSize == sizeof(double)) {
    _ASSERTE(_non_descending_(pvBegin, nElems, nSize));
    _ASSERTE(pvResult == (char *)pvBegin + nElems * nSize
                 ? *(double *)((char *)pvBegin + (nElems - 1) * nSize) <=
                       *(double *)pvKey
                 : 1);

    // if found
    if (pvResult != (char *)pvBegin + nElems * nSize) {
      _ASSERTE(*(double *)pvResult > *(double *)pvKey);
      if (pvResult != (char *)pvBegin)
        _ASSERTE(*(double *)((char *)pvResult - nSize) <= *(double *)pvKey);
    }
  }
#endif

  return pvResult;
}

static unsigned __stdcall _sort_interpolate(
    const void *pcBegin, // linear coordinates
    unsigned ncSize,     // size of coordinates (in bytes)
    const void *pvBegin, // # elements in values
    unsigned nvSize,     // size of values (in bytes)
    unsigned nElems,     // # elements in the linear coordinates
    const void *pc,      // coordinate to be interpolated
    void *pv             /// results
) {
  // find the coordinate to be interpolated in pcBegin
  unsigned n =
      ((char *)find_min_greater_than(pc, pcBegin, nElems, ncSize, dless) -
       (char *)pcBegin) /
      ncSize;
  // if found  , interpolate; else extrapolate
  n = n < nElems ? n : (nElems - 1);
  *(double *)pv = *(double *)((char *)pvBegin + n * nvSize);
  return n;
}

void sort_interpolate(const void *pcBegin, // linear coordinates
                      unsigned ncSize,     // size of coordinates (in bytes)
                      const void *pvBegin, // # elements in values
                      unsigned nvSize,     // size of values (in bytes)
                      unsigned nElems, // # elements in the linear coordinates
                      const void *pcInBegin, // coordinates to be interpolated
                      unsigned nInElems, // # elements in the linear coordinates
                      void *pvOut        /// results
) {
  const char *pcInEnd = (char *)pcInBegin + nInElems * ncSize;
  const void *pcEnd = (char *)pcBegin + nElems * ncSize;
  const char *pcin = pcInBegin, *pc = pcBegin, *pv = pvBegin;
  char *pvout = pvOut;

  unsigned n = 0;

  _ASSERTE(_non_descending_(pcBegin, nElems, ncSize));
  _ASSERTE(_non_descending_(pcInBegin, nInElems, ncSize));

  // NB: increment pcBegin improves performance
  for (; pcin < pcInEnd; pcin += ncSize, pc += n * ncSize, pv += n * nvSize,
                         pvout += nvSize, nElems -= n)
    n = _sort_interpolate(pc, ncSize, pv, nvSize, nElems, pcin, pvout);
}

void memswap(void *pm1, void *pm2, unsigned nSize) {
  void *p = _alloca(nSize);
  memcpy(p, pm1, nSize);
  memcpy(pm1, pm2, nSize);
  memcpy(pm2, p, nSize);
}

const void *_reverse(void *pvBegin,   // pointer to base element
                     unsigned nElems, // # elements
                     unsigned nSize   // size of elements (in bytes)
) {
  void *pvEnd = (char *)pvBegin + (nElems - 1) * nSize;
  const void *preturn = pvBegin;
  for (; (char *)pvBegin <= (char *)pvEnd - nSize;
       (char *)pvBegin += nSize, (char *)pvEnd -= nSize)
    memswap(pvBegin, pvEnd, nSize);
  return preturn;
}

#ifdef _DEBUG
int _non_descending_(const void *pvBegin, // pointer to base element
                     unsigned nElems,     // # elements
                     unsigned nSize       // size of elements (in bytes)
) {
  const void *pvEnd = (const char *)pvBegin + nElems * nSize;
  const char *pv = (const char *)pvBegin + nSize;
  for (; pv < (const char *)pvEnd; pv += nSize) {
    // bad  , will come back later ....
    const double dL = nSize == sizeof(double) ? *((double *)pv) : *((int *)pv);
    const double dR = nSize == sizeof(double) ? *((double *)(pv - nSize))
                                              : *((int *)(pv - nSize));
    if (dL < dR)
      return 0;
  }
  return 1;
}
#endif

int __stdcall dless(const double *pdLeft, const double *pdRight) {
  return (*pdLeft < *pdRight) ? -1 : (*pdLeft > *pdRight);
}

int __stdcall lless(const long *pdLeft, const long *pdRight) {
  return (*pdLeft < *pdRight) ? -1 : (*pdLeft > *pdRight);
}

int __stdcall dgreater(const double *pdLeft, const double *pdRight) {
  //	const int n = dless(pdLeft  ,pdRight);
  //	return n?-n:0;
  return (*pdLeft > *pdRight) ? -1 : (*pdLeft < *pdRight);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
Err check_null_ptr(void *ptr) {
  if (!ptr) {
    return "Null pointer detected!";
  }
  return 0;
}

void **_alloca_matrix(void **pm, void *pdata, int nR_Begin, int nR_End,
                      int nC_Begin, int nC_End, unsigned nSize) {
  const int nSize_R = nR_End - nR_Begin;
  const int nSize_C = nC_End - nC_Begin;
  const void **p = pm, **pm_end = pm + nSize_R;
  const char *data = (char *)pdata - nC_Begin * nSize; // offset by nC_Begin
  const int nSize_Inc = nSize_C * nSize;

  // set to 0.
  memset(pdata, 0, nSize_R * nSize_C * nSize);
  memset(pm, 0, nSize_R);

  for (; p < pm_end && (*p = data); ++p, data += nSize_Inc)
    ;

  return pm - nR_Begin; // offset by nR_Begin
}

/*
void _alloca_multiple_matrices(
           void **pm  ,
           void *pdata  ,
           int nR_Begin  ,
           int nR_End  ,
           int nC_Begin  ,
           int nC_End  ,
           unsigned nSize
           )
{

}

*/
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
Err _alloc_dvector_(int nBegin, int nEnd, double **ppd) {
  char *szErrMessage = 0;
  *ppd = dvector(nBegin, nEnd - 1);
  if (szErrMessage = check_null_ptr(*ppd))
    return szErrMessage;
  return 0;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
void _free_dvector(double *pd_Begin, int nSize) {
  if (pd_Begin)
    free_dvector(pd_Begin, 0, nSize - 1);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
void _sort2(const void *pv1Begin, // base 1
            unsigned n1Elems,     // # elements in pv1Begin
            const void *pv2Begin, // base 2
            unsigned n2Elems,     // # elements pv2Begin
            unsigned nvSize,      // size of v1/v2 (in bytes)
            void *pvOut) {
  int _dcompare(const void *, const void *);
  memcpy(pvOut, pv1Begin, n1Elems * nvSize);
  memcpy((char *)pvOut + nvSize * n1Elems, pv2Begin, n2Elems * nvSize);
  qsort(pvOut, n1Elems + n2Elems, nvSize, _dcompare);
}

void _sort3(const void *pv1Begin, // base 1
            unsigned n1Elems,     // # elements in pv1Begin
            const void *pv2Begin, // base 2
            unsigned n2Elems,     // # elements pv2Begin
            const void *pv3Begin, // base 3
            unsigned n3Elems,     // # elements pv2Begin
            unsigned nvSize,      // size of v1/v2/v3 (in bytes)
            void *pvOut) {
  int _dcompare(const void *, const void *);
  memcpy(pvOut, pv1Begin, n1Elems * nvSize);
  memcpy((char *)pvOut + nvSize * n1Elems, pv2Begin, n2Elems * nvSize);
  memcpy((char *)pvOut + nvSize * n1Elems + nvSize * n2Elems, pv3Begin,
         n3Elems * nvSize);

  qsort(pvOut, n1Elems + n2Elems + n3Elems, nvSize, _dcompare);
}

int _dcompare(const void *pvL, const void *pvR) {
  const double *pdL = (const double *)pvL;
  const double *pdR = (const double *)pvR;
  if (*pdL == *pdR)
    return 0;
  return (*pdL < *pdR) ? -1 : 1;
}

int _dequal(const void *pvL, const void *pvR) {
  const double *pdL = (const double *)pvL;
  const double *pdR = (const double *)pvR;
  return (*pdL == *pdR) ? 0 : 1;
}

int _nequal(const void *pvL, const void *pvR) {
  const int *pnL = (const int *)pvL;
  const int *pnR = (const int *)pvR;
  return (*pnL == *pnR) ? 0 : 1;
}

void _bmatch_ex(const void *pcBegin,   // linear coordinates
                unsigned ncSize,       // size of coordinates (in bytes)
                const void *pvBegin,   // # elements in values
                unsigned nvSize,       // size of values (in bytes)
                unsigned nElems,       // # elements in the linear coordinates
                const void *pcInBegin, // sub linear coordinates
                unsigned nInElems, // # elements in the sub linear coordinates
                void *pInvOut) {
  const char *pInc = pcInBegin,
             *pIncEnd = (char *)pcInBegin + nInElems * ncSize;
  char *pInv = pInvOut;

  _ASSERTE(_non_descending_(pcBegin, nElems, ncSize));
  //_ASSERTE(_non_descending_(pcInBegin  ,nInElems  ,ncSize));

  for (; pInc < pIncEnd; pInc += ncSize, pInv += nvSize) {
    const char *p = bsearch(pInc, pcBegin, nElems, ncSize, _dcompare);

    // if found  , use the match
    if (p) {
      const int n = (p - (char *)pcBegin) / ncSize;
      *(double *)pInv = *(double *)((char *)pvBegin + n * nvSize);
    } else // interpolate
      _sort_interpolate(pcBegin, ncSize, pvBegin, nvSize, nElems, pInc, pInv);
  }
}

void _bmatch(const void *pcBegin,    // linear coordinates
             unsigned ncSize,        // size of coordinates (in bytes)
             const void *pvBegin,    // # elements in values
             unsigned nvSize,        // size of values (in bytes)
             unsigned nElems,        // # elements in the linear coordinates
             const void *psubcBegin, // sub linear coordinates
             unsigned nsubElems,     // # elements in the sub linear coordinates
             void *psubvOut) {
  const char *psubc = psubcBegin,
             *psubcEnd = (char *)psubcBegin + nsubElems * ncSize;
  const char *pc = pcBegin;
  const char *pv = pvBegin;
  char *psubv = psubvOut;
  int n = 0;

  _ASSERTE(nsubElems <= nElems);
  _ASSERTE(_non_descending_(pcBegin, nElems, ncSize));
  _ASSERTE(_non_descending_(psubcBegin, nsubElems, ncSize));

  for (; psubc < psubcEnd; psubc += ncSize, psubv += nvSize, pc += n * ncSize,
                           pv += n * nvSize, nElems -= n) {
    n = ((char *)bsearch(psubc, pc, nElems, ncSize, _dcompare) - pc) / ncSize;
    *(double *)psubv = *(double *)(pv + n * nvSize);
  }
}

/// return 1 is vector is ascending (X[n+1] >= X[n]   , for nBegin =<n< nEnd)
//// return 0 otherwise
int _is_ascending_(const double *pdX_Begin, const double *pdX_End) {
  const double *pd = pdX_Begin;
  for (; pd < pdX_End - 1; ++pd) {
    if (pd[1] < (*pd))
      return 0;
  }
  return 1;
}

void decomp_sysmetric_tridiagonal(
    /// matrix to be decomposed
    const double *pdDiag, int nDiag_Begin, int nDiag_End,
    const double *pdSubDiag, int nSubDiag_Begin,
    /// result
    double *pdSubDiagL, int nSubDiagL_Begin, double *pdPivot,
    int nPivot_Begin) {
  int nI;
  double dProd = 0.;
  const int nMatrixSize = nDiag_End - nDiag_Begin;

  for (nI = 0; nI < nMatrixSize - 1;
       ++nI) // NB: pdSubDiagL is one element less than nSize
  {
    pdPivot[nPivot_Begin + nI] = pdDiag[nDiag_Begin + nI] - dProd;
    pdSubDiagL[nSubDiagL_Begin + nI] =
        pdSubDiag[nSubDiag_Begin + nI] / pdPivot[nPivot_Begin + nI];
    dProd = pdSubDiag[nSubDiag_Begin + nI] * pdSubDiagL[nSubDiagL_Begin + nI];
  }

  // dont forget the last pivot
  pdPivot[nPivot_Begin + nMatrixSize - 1] = pdDiag[nDiag_End - 1] - dProd;
}

const char *_fill_swapcashrate_diag(const char *szYC, const long *plStart_Begin,
                                    const long *plStart_End, long lEndDate,
                                    const char *szFreq, const char *szBasis,
                                    const char *szRefRate,
                                    double *pdSwapCashRate_Begin) {
  SrtCompounding eFreq = SRT_SIMPLE;
  SrtBasisCode eBasis = BASIS_ACT_ACT;

  const char *szErr = interp_basis(szBasis, &eBasis);
  if (szErr)
    return szErr;

  interp_compounding(szFreq, &eFreq);

  for (; plStart_Begin < plStart_End; ++pdSwapCashRate_Begin, ++plStart_Begin) {
    _ASSERTE(*plStart_Begin <= lEndDate);
    *pdSwapCashRate_Begin =
        swp_f_swapcashrate(*plStart_Begin, lEndDate, eBasis, eFreq,
                           (char *)szYC, (char *)szRefRate);
  }

  return 0;
}

Err _fill_df(const char *szYieldCurve, long lToday, const long *plDate_Begin,
             const long *plDate_End,
             // output
             double *pdDF_Begin) {
  int nI;
  double dDF = 0.;
  const int nSize = (plDate_End - plDate_Begin);
  _ASSERTE(nSize > 0);

  for (nI = 0; nI < nSize; ++nI) {
    if (lToday > plDate_Begin[nI])
      return "_fill_df(...) : invalid dates!";
    dDF = swp_f_df(lToday, plDate_Begin[nI], szYieldCurve);
    if (dDF > 1. || dDF < 0.)
      return "_fill_df(...) : invalid discount factor!";
    pdDF_Begin[nI] = dDF;
  }
  return 0;
}

const double *_fill_tenor(long lBase, const long *plDate_Begin,
                          const long *plDate_End,
                          // output
                          double *pdTenor_Begin) {
  int nI;
  const int nSize = (plDate_End - plDate_Begin);
  _ASSERTE(nSize >= 0);
  for (nI = 0; nI < nSize; ++nI)
    pdTenor_Begin[nI] = (plDate_Begin[nI] - lBase) * YEARS_IN_DAY;

  return pdTenor_Begin;
}

const long *_fill_date(long lBase, const double *pdTenor_Begin,
                       const double *pdTenor_End,
                       // output
                       long *plDate_Begin) {
  const double *pdTenor = pdTenor_Begin;
  long *plDate = plDate_Begin;
  for (; pdTenor < pdTenor_End; ++pdTenor, ++plDate)
    *plDate = lBase + (int)((*pdTenor) / YEARS_IN_DAY + 0.5);

#ifdef _DEBUG
  {
    const int nSize = (pdTenor_End - pdTenor_Begin);
    double *pdTenorBack = _alloca(nSize * sizeof(double));
    _fill_tenor(lBase, plDate_Begin, plDate_Begin + nSize, pdTenorBack);
    for (; pdTenor_Begin < pdTenor_End; ++pdTenor_Begin, ++pdTenorBack)
      _ASSERTE(fabs(*pdTenor_Begin - *pdTenorBack) < 1. / 360.);
  }
#endif

  return plDate_Begin;
}

const void *_memset(const void *pKey,
                    void *pOut_Begin, // pointer to base element
                    unsigned nElems,  // # elements
                    unsigned nSize    // size of Key (in bytes)
) {
  const void *pOut_End = (char *)pOut_Begin + nElems * nSize;
  void *pv = pOut_Begin;

  while (pv < pOut_End) {
    pv = (char *)memcpy(pv, pKey, nSize) + nSize;
  };

  return pOut_Begin;
}

const double *_multiply_vector2(const double *pdIn_Begin,
                                const double *pdIn_End, double dConstant,
                                double *pdOut_Begin) {
  double *pd = pdOut_Begin;
  for (; pdIn_Begin < pdIn_End; ++pdIn_Begin, ++pd)
    *pd = (*pdIn_Begin) * dConstant;
  return pdOut_Begin;
}

const double *_multiply_vector(const double *pdIn1_Begin,
                               const double *pdIn1_End,
                               const double *pdIn2_Begin, double *pdOut_Begin) {
  int nI = 0;
  const int nSize = pdIn1_End - pdIn1_Begin;
  for (nI = 0; nI < nSize; ++nI)
    pdOut_Begin[nI] = pdIn1_Begin[nI] * pdIn2_Begin[nI];
  return pdOut_Begin;
}

const double *_max_dEle(const double *pdBegin, const double *pdEnd) {
  const double *pdMaxEle = pdBegin;
  for (; pdBegin < pdEnd; ++pdBegin)
    pdMaxEle = (*pdBegin) > (*pdMaxEle) ? pdBegin : pdMaxEle;
  return pdMaxEle;
}

const int *_max_nEle(const int *pnBegin, const int *pnEnd) {
  const int *pnMaxEle = pnBegin;
  for (; pnBegin < pnEnd; ++pnBegin)
    pnMaxEle = (*pnBegin) > (*pnMaxEle) ? pnBegin : pnMaxEle;
  return pnMaxEle;
}

const double *_add_vector(const double *pdIn1_Begin, const double *pdIn1_End,
                          const double *pdIn2_Begin, double *pdOut_Begin) {
  int nI = 0;
  const int nSize = pdIn1_End - pdIn1_Begin;
  for (nI = 0; nI < nSize; ++nI)
    pdOut_Begin[nI] = pdIn1_Begin[nI] + pdIn2_Begin[nI];
  return pdOut_Begin;
}

const char *_alloc_dvector(const double **ppd_Begin, int nSize) {
  *ppd_Begin = 0;
  *ppd_Begin = dvector(0, nSize - 1);
  if (!(*ppd_Begin))
    return "_alloc_dvector(.): mem allocation failed!";
  return 0;
}

const char *_alloc_dmatrix(const double ***pppd_Begin, int nR_Begin, int nR_End,
                           int nC_Begin, int nC_End) {
  *pppd_Begin = 0;
  *pppd_Begin = dmatrix(nR_Begin, nR_End - 1, nC_Begin, nC_End - 1);

  if (!(*pppd_Begin))
    return "_alloc_dmatrix(...): mem allocation failed!";
  return 0;
}

double _NormDensity(double dOne_Over_Sqr2Pi, double dX, double dMu,
                    double dVar) {
  if (dVar == 0.)
    return 1.;
  else {
    const double dX_ = (dX - dMu);
    _ASSERTE(dVar > 0.);
    return dOne_Over_Sqr2Pi * exp(-0.5 * dX_ * dX_ / dVar) / sqrt(dVar);
  }
}

const char *_free2(const char *szErr, void *p1, void *p2) {
  if (p1)
    free(p1);
  if (p2)
    free(p2);
  return szErr;
}

void _free_and_zero(const void *pPointer) {
  if (pPointer)
    free((void *)pPointer);
  pPointer = 0;
}

double _Sum_dVec(const double *pdIn_Begin, const double *pdIn_End) {
  double dResult = 0.;
  for (; pdIn_Begin < pdIn_End; ++pdIn_Begin)
    dResult += *pdIn_Begin;
  return dResult;
}

double _Sum_dMatrix(const double **ppdMatrix, int nR_Begin, int nR_End,
                    int nC_Begin, int nC_End) {
  double dResult = 0.;
  int nR, nC;
  for (nR = nR_Begin; nR < nR_End; ++nR) {
    for (nC = nC_Begin; nC < nC_End; ++nC)
      dResult += ppdMatrix[nR][nC];
  }

  return dResult;
}

long _ToDate(long lBase, double dT) {
  long lDate;
  _fill_date(lBase, &dT, 1 + &dT, &lDate);
  return lDate;
}

#ifdef _DEBUG
void _print_vector(const double *pdBegin, int nSize) {
  int nI;
  _RPT0(_CRT_WARN, "\n Vector \n");
  for (nI = 0; nI < nSize; ++nI) {
    _RPT1(_CRT_WARN, "\t%.10f", ((double)pdBegin[nI] + 0.));
  }
  _RPT0(_CRT_WARN, "\n");
}

void _print_nvector(const int *pnBegin, int nSize) {
  int nI;
  _RPT0(_CRT_WARN, "\n Vector \n");
  for (nI = 0; nI < nSize; ++nI) {
    _RPT1(_CRT_WARN, "\t%.0f", (double)pnBegin[nI] + 0.);
  }
  _RPT0(_CRT_WARN, "\n");
}

void _print_matrix(const double **ppdMatrix, int nR_Begin, int nR_End,
                   int nC_Begin, int nC_End) {
  int nR, nC;
  _RPT0(_CRT_WARN, "\n Matrix \n");
  for (nR = nR_Begin; nR < nR_End; ++nR) {
    for (nC = nC_Begin; nC < nC_End; ++nC) {
      _RPT1(_CRT_WARN, "\t%.6f", ppdMatrix[nR][nC]);
    }
    _RPT0(_CRT_WARN, "\n");
  }
  _RPT0(_CRT_WARN, "\n");
}
#endif //_DEBUG

const char *Decode_RefRate(const char *szRefRate, char **pszFeq,
                           char **pszBasis) {
  const char *szErr = 0;
  SrtCompounding eFreq;
  SrtBasisCode eBasis;

  szErr = swp_f_get_ref_rate_details((char *)szRefRate, &eBasis, &eFreq);
  if (szErr)
    return szErr;

  szErr = translate_compounding(pszFeq, eFreq);
  if (szErr)
    return szErr;

  return translate_basis(pszBasis, eBasis);
}

double _PayRec(const char *szPayRec) {
  // return szPayRec=="REC" ? 1.:-1.;
  return _stricmp(szPayRec, "REC") ? -1 : 1;
}
