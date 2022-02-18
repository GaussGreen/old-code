#include "AmortMidatADIGrid.h"

#define _MAX_NUM_STD_ 7.

//-------------------------------------------------------------------------------------------------------
//	Description: computes # of intervals between dHB and dLB according to dDelta_Base
//
//	Returns :  indicated, rounded up
//
//-------------------------------------------------------------------------------------------------------
int NumInterval(
    double dB1,  // bound 1
    double dB2,  // bound 2
    double dDelta_Base)
{
    // establish the base
    const double dNumInterval = fabs(dB1 - dB2) / dDelta_Base;
    const int    nNumInterval =
        (int)(dNumInterval + 0.5) > 1 ? (int)(dNumInterval + 0.5) : 1;  /// truncate
    //_ASSERTE(dT_HB>=dT_LB && nNumTInterval>=0);
    _ASSERTE(dDelta_Base != 0.);
    _ASSERTE(nNumInterval >= 0);
    return nNumInterval;
}

//-------------------------------------------------------------------------------------------------------
//	Description: divides up [dFront,dBack] in nNumInterval idential intevals
//						and fills pdBegin with the result of the
// discretization
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
void Discretize_1D(double dFront, double dBack, int nNumInterval, double* pdBegin)
{
    // const int nNumInterval = pdEnd-pdBegin-1;
    const double dDelta = (dBack - dFront) / ((double)nNumInterval);
    double *     pd = pdBegin, *pdEnd = pdBegin + nNumInterval + 1;
    //_ASSERTE(dT_HB>=dT_LB);
    _ASSERTE(nNumInterval != 0.);

    // populate the grid in the time direction
    for (; pd < pdEnd; ++pd)
        *pd = dFront + dDelta * (pd - pdBegin);

    _ASSERTE(fabs(*(pdEnd - 1) - dBack) < 1.e15);
}

//-------------------------------------------------------------------------------------------------------
//	Description: re-instate the 2D Grid to its max size
//
//	Returns : begin and end indeces of the Grid
//-------------------------------------------------------------------------------------------------------
void Max_Grid(
    const int* pnSizeX1_Begin,
    const int* pnSizeX1_End,
    const int* pnSizeX3_Begin,
    _Grid*     pGrid,
    int*       pnX1_Begin,
    int*       pnX1_End,
    int*       pnX3_Begin,
    int*       pnX3_End)
{
    const int nMaxSize_X1 = *_max_nEle(pnSizeX1_Begin, pnSizeX1_End);
    const int nMaxSize_X3 =
        *_max_nEle(pnSizeX3_Begin, pnSizeX3_Begin + (pnSizeX1_End - pnSizeX1_Begin));
    _Shift_X1(nMaxSize_X1, pGrid, pnX1_Begin, pnX1_End);
    _Shift_X3(nMaxSize_X3, pGrid, pnX3_Begin, pnX3_End);
    _ASSERTE(*pnX1_Begin >= 0 && *pnX3_Begin >= 0);
}

//-------------------------------------------------------------------------------------------------------
//	Description: shifts *ppdX_Begin and *ppdX_End equal distance to reach size nSize_After
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
static int __stdcall _Shift_dVec_(
    int nSize_After, const double** ppdX_Begin, const double** ppdX_End)
{
    // bounds prior to adjustment
    const int nSize_Before = *ppdX_End - *ppdX_Begin;
    const int nShift       = (nSize_After - nSize_Before) / 2;

    // NB: if Grid shrinks, nShift_X <0; else nShift_X>0
    *ppdX_Begin -= nShift;
    *ppdX_End += nShift;

    _ASSERTE(nSize_Before % 2 == 1);
    _ASSERTE(nSize_After % 2 == 1);

    return nShift;
}

//-------------------------------------------------------------------------------------------------------
//	Description: shifts X1 dimensino of the grid to nSize_After
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
void __stdcall _Shift_X1_(int nSize_After, _Grid* pGrid)
{
    _Shift_dVec_(nSize_After, &pGrid->m_pdX1_Begin, &pGrid->m_pdX1_End);
    _ASSERTE(pGrid->m_pdX1_Begin >= pGrid->m_pdX1_Begin_);
    _ASSERTE(pGrid->m_pdX1_End <= pGrid->m_pdX1_Begin_ + 2 * pGrid->m_nX1_Zero + 1);
}

//-------------------------------------------------------------------------------------------------------
//	Description: shifts X3 dimension of the grid to nSize_After
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
void __stdcall _Shift_X3_(int nSize_After, _Grid* pGrid)
{
    _Shift_dVec_(nSize_After, &pGrid->m_pdX3_Begin, &pGrid->m_pdX3_End);
    _ASSERTE(pGrid->m_pdX3_Begin >= pGrid->m_pdX3_Begin_);
    _ASSERTE(pGrid->m_pdX3_End <= pGrid->m_pdX3_Begin_ + 2 * pGrid->m_nX3_Zero + 1);
}

//-------------------------------------------------------------------------------------------------------
//	Description: shifts *ppdX_Begin and *ppdX_End equal distance to reach size nSize_After
//
//	Returns : updated nBegin and nEnd
//-------------------------------------------------------------------------------------------------------
static void _Shift_dVec(
    int            nSize_After,
    const double** ppdX_Begin,
    const double** ppdX_End,
    int*           pnBegin,  // on entry, current begin index; output: updated begin index
    int*           pnEnd     // on entry, current end index; output: updated end index
)
{
    const int nShift = _Shift_dVec_(nSize_After, ppdX_Begin, ppdX_End);
    *pnBegin -= nShift;
    *pnEnd += nShift;
    _ASSERTE((*ppdX_End - *ppdX_Begin) == (*pnEnd - *pnBegin));
}

//-------------------------------------------------------------------------------------------------------
//	Description: shifts X1 dimension of the grid to nSize_After
//
//	Returns : nX1_Begin and nX1_End after being resized
//-------------------------------------------------------------------------------------------------------
void __stdcall _Shift_X1(
    int    nSize_After,
    _Grid* pGrid,
    int*   pnX1_Begin,  // result
    int*   pnX1_End     // result
)
{
    *pnX1_Begin = pGrid->m_pdX1_Begin - pGrid->m_pdX1_Begin_;
    *pnX1_End   = pGrid->m_pdX1_End - pGrid->m_pdX1_Begin_;
    _Shift_dVec(nSize_After, &pGrid->m_pdX1_Begin, &pGrid->m_pdX1_End, pnX1_Begin, pnX1_End);
}

//-------------------------------------------------------------------------------------------------------
//	Description: shifts X3 dimension of the grid to nSize_After
//
//	Returns : nX3_Begin and nX3_End after being resized
//-------------------------------------------------------------------------------------------------------
void __stdcall _Shift_X3(int nSize_After, _Grid* pGrid, int* pnX3_Begin, int* pnX3_End)
{
    *pnX3_Begin = pGrid->m_pdX3_Begin - pGrid->m_pdX3_Begin_;
    *pnX3_End   = pGrid->m_pdX3_End - pGrid->m_pdX3_Begin_;
    _Shift_dVec(nSize_After, &pGrid->m_pdX3_Begin, &pGrid->m_pdX3_End, pnX3_Begin, pnX3_End);
}

//-------------------------------------------------------------------------------------------------------
//	Description: local helper for computing the variance of X3
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
static double _Variance_X3(double dPHI1, double dPHI2, double dPHI12, double dRho, double dAlpha)
{
    // variance of r1
    const double dVar_X1    = dPHI1;
    const double dAlphaRho  = dRho * dAlpha;
    const double dVar_X2    = dPHI2;
    const double dCoVar_X12 = dRho * dPHI12;

    // if 1 Factor
    double dVar_R3 = dVar_X2;
    dVar_R3 += dAlphaRho * dAlphaRho * dVar_X1;
    dVar_R3 -= 2. * dAlphaRho * dCoVar_X12;
    dVar_R3 /= (dAlpha * dAlpha);
    dVar_R3 /= (1. - dRho * dRho);

    return dVar_R3;
}

static void _preprocess_grid_(
    const double* pdVarX_Begin,
    const double* pdVarX_End,
    int           nBase_NumX,
    int*          pnSizeX_Begin,
    double*       pdDeltaX)
{
    const double dMaxVarX = *_max_dEle(pdVarX_Begin, pdVarX_End);

    int *pnSizeX = pnSizeX_Begin, *pnSizeX_End = pnSizeX_Begin + (pdVarX_End - pdVarX_Begin);

    // delta x is determined by the max variance
    if (nBase_NumX == 1)  // 1 factor
    {
        _ASSERTE(dMaxVarX == 0.);
        (*pdDeltaX) = 0.;
    }

    else  // 2 factor
    {
        _ASSERTE(dMaxVarX != 0.);
        (*pdDeltaX) = 2. * _MAX_NUM_STD_ * sqrt(dMaxVarX) / (nBase_NumX - 1);
        _ASSERTE(nBase_NumX != 1);
    }

    for (; pdVarX_Begin < pdVarX_End; ++pdVarX_Begin, ++pnSizeX)
    {
        *pnSizeX = (*pdDeltaX) == 0.
                       ? 1
                       : 1 + (int)(2. * _MAX_NUM_STD_ * sqrt(*pdVarX_Begin) / (*pdDeltaX));
        *pnSizeX = _validate_gridpoints(*pnSizeX);
    }

    _ASSERTE((*pdDeltaX) >= 0.);
    _ASSERTE(nBase_NumX == _validate_gridpoints(nBase_NumX));
    _ASSERTE(*_max_nEle(pnSizeX_Begin, pnSizeX_End) == nBase_NumX);

    // in case the grid is not monotonically decreasing in size, do sth
    {
        int* pnSize = pnSizeX_Begin + 1;
        for (; pnSize < pnSizeX_End; ++pnSize)
            *pnSize = *pnSize < *(pnSize - 1) ? *(pnSize - 1) : *pnSize;
    }
}
//-------------------------------------------------------------------------------------------------------
//	Description: computes the variance of X3
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
void Variance_X3(
    const double* pdPHI1_Begin,
    const double* pdPHI1_End,
    const double* pdPHI2_Begin,
    const double* pdPHI12_Begin,
    double        dRho,
    double        dAlpha,
    double*       pdVarX3_Begin)
{
    for (; pdPHI1_Begin < pdPHI1_End;
         ++pdPHI1_Begin, ++pdPHI2_Begin, ++pdPHI12_Begin, ++pdVarX3_Begin)
    {
        *pdVarX3_Begin = _Variance_X3(*pdPHI1_Begin, *pdPHI2_Begin, *pdPHI12_Begin, dRho, dAlpha);
    }
}

const char* _init_amort_grid(
    const double* pdVarX1_Begin,  // =_alloca(nNumEx*sizeof(double));
    const double* pdVarX1_End,    // =_alloca(nNumEx*sizeof(double));
    const double* pdVarX3_Begin,  // =_alloca(nNumEx*sizeof(double));
    int           nBase_NumX1,
    int           nBase_NumX3,
    double        dEx_Back,
    int           nNumT,
    /// results
    _Grid* pGrid,
    int*   pnSizeX3_Begin,
    int*   pnSizeX1_Begin)
{
    const char* szErr = 0;
    // alias
    const int nNumEx = pdVarX1_End - pdVarX1_Begin;
    double    dDeltaX3, dDeltaX1;

    nBase_NumX1 = _validate_gridpoints(nBase_NumX1);
    nBase_NumX3 = _validate_gridpoints(nBase_NumX3);

    _preprocess_grid_(
        pdVarX1_Begin, pdVarX1_Begin + nNumEx, nBase_NumX1, pnSizeX1_Begin, &dDeltaX1);
    _preprocess_grid_(
        pdVarX3_Begin, pdVarX3_Begin + nNumEx, nBase_NumX3, pnSizeX3_Begin, &dDeltaX3);

    if (szErr = _init_Grid(
            *_max_nEle(pnSizeX1_Begin, pnSizeX1_Begin + nNumEx),
            dDeltaX1,
            *_max_nEle(pnSizeX3_Begin, pnSizeX3_Begin + nNumEx),
            dDeltaX3,
            nNumT,
            dEx_Back,
            pGrid))
    {
        return szErr;
    }

    return 0;
}

int _DefaultNumT(const double* pdEx_Begin, const double* pdEx_End, int nNumPointT_In)
{
    const double dEx_Back = *(pdEx_End - 1);
    const long   nNumEx   = pdEx_End - pdEx_Begin;
    int          nDefault = (int)(12. * dEx_Back);  // number of months
    nDefault              = nDefault > nNumEx ? nDefault : nNumEx;
    nDefault              = nDefault > 128 ? nDefault : 128;  // lower bound 128
    return nNumPointT_In > nDefault ? nNumPointT_In : nDefault;
}

int _DefaultNumX(int nNumPointX_In)
{
    const int nDefaultNumPointX = 128;
    return nNumPointX_In > nDefaultNumPointX ? nNumPointX_In : nDefaultNumPointX;
}

void _free_Grid(_Grid* pGrid)
{
    if (pGrid)
    {
        const int nX1_Begin   = 0;
        const int nX3_Begin   = 0;
        const int nNumPointX1 = 2 * pGrid->m_nX1_Zero + 1;
        const int nNumPointX3 = 2 * pGrid->m_nX3_Zero + 1;

        _free_dvector((double*)pGrid->m_pdX1_Begin_, nNumPointX1);
        _free_dvector((double*)pGrid->m_pdX3_Begin_, nNumPointX3);
    }
}

int _validate_gridpoints(int nInitPoint)
{
    return (nInitPoint = (nInitPoint % 2 ? nInitPoint : nInitPoint + 1)) > 1 ? nInitPoint : 1;
}

const char* _init_Grid(
    int    nInitNumPointX1,
    double dDeltaX1,
    int    nInitNumPointX3,
    double dDeltaX3,
    int    nNumPointT,
    double dT,
    _Grid* pGrid)
{
    const char* szErr = 0;
    const char* _space_discretize_(int, double, double**, int*);

    // number of space discretization points must be
    // an odd number no less than 3
    const int nNumPointX1 = _validate_gridpoints(nInitNumPointX1);
    const int nNumPointX3 = _validate_gridpoints(nInitNumPointX3);

    *(double*)(&pGrid->m_dDeltaT) = dT / nNumPointT;

    /// space direction
    if (szErr = _space_discretize_(
            nNumPointX1, dDeltaX1, (double**)(&pGrid->m_pdX1_Begin), (int*)(&pGrid->m_nX1_Zero)))
        return szErr;

    if (szErr = _space_discretize_(
            nNumPointX3, dDeltaX3, (double**)(&pGrid->m_pdX3_Begin), (int*)(&pGrid->m_nX3_Zero)))
        return szErr;

    *(const double**)(&pGrid->m_pdX1_Begin_) = pGrid->m_pdX1_Begin;
    *(const double**)(&pGrid->m_pdX3_Begin_) = pGrid->m_pdX3_Begin;

    // initialize m_pdX1_End
    pGrid->m_pdX1_End = pGrid->m_pdX1_Begin + nNumPointX1;

    // initialize m_pdX3_End
    pGrid->m_pdX3_End = pGrid->m_pdX3_Begin + nNumPointX3;

    return 0;
}

static const char* _space_discretize_(
    /// inputs
    int    nNumPoint,
    double dDeltaX,
    /// outputs
    double** ppd_Begin,
    int*     pn_Zero)
{
    int         nI;
    const char* szErr     = 0;
    double      dhalfline = 0., *pd_Begin = 0;

    _ASSERTE(nNumPoint % 2 == 1);

    *pn_Zero  = (nNumPoint - 1) / 2;
    dhalfline = (*pn_Zero) * dDeltaX;
    // fixme if(szErr = _alloc_dvector(ppd_Begin,nNumPoint)) return szErr;
    pd_Begin = (*ppd_Begin);
    for (nI = 0; nI < *pn_Zero; ++nI)
    {
        pd_Begin[nI]                 = -dhalfline + dDeltaX * nI;
        pd_Begin[nNumPoint - 1 - nI] = -pd_Begin[nI];
    }
    pd_Begin[*pn_Zero] = 0.;

    return 0;
}
