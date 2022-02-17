// prevent multiple inclusions
#pragma once

////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                // in the debug information"

// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "AmortMidatADIUtils.h"
#include "AmortmidatADIPrice.h"

typedef struct _amortizingmidat_grid {
  const double *const m_pdX1_Begin_;
  const double *m_pdX1_Begin; /// discretization in space X1 direction
  const double *m_pdX1_End;
  const int m_nX1_Zero;

  const double *const m_pdX3_Begin_;
  const double *m_pdX3_Begin; /// discretization in space X1 direction
  const double *m_pdX3_End;
  const int m_nX3_Zero;

  const double m_dDeltaT;

  ////////////
  // storage
  // double **m_ppd_3L1;
  double **m_ppd_3DInv1;
  double **m_ppd_3U1;
  double **m_ppd_1L3;
  double **m_ppd_1DInv3;
  double **m_ppd_1U3;

} _Grid;

int _validate_gridpoints(int nInitPoint);

void _free_Grid(_Grid *pGrid);

void Discretize_1D(double dFront, double dBack, int nNumInterval,
                   double *pdBegin);

int NumInterval(double dFront, // bound 1
                double dBack,  // bound 2
                double dDelta_Base);

void Max_Grid(const int *pnSizeX1_Begin, const int *pnSizeX1_End,
              const int *pnSizeX3_Begin, _Grid *pGrid, int *pnX1_Begin,
              int *pnX1_End, int *pnX3_Begin, int *pnX3_End);

void __stdcall _Shift_X1(int nSize_After, _Grid *pGrid, int *pnR_Begin,
                         int *pnR_End);
void __stdcall _Shift_X3(int nSize_After, _Grid *pGrid, int *pnR_Begin,
                         int *pnR_End);
void __stdcall _Shift_X1_(int, _Grid *);
void __stdcall _Shift_X3_(int, _Grid *);

void Variance_X3(const double *pdPHI1_Begin, const double *pdPHI1_End,
                 const double *pdPHI2_Begin, const double *pdPHI12_Begin,
                 double dRho, double dAlpha, double *pdVarX3_Begin);

int _DefaultNumT(const double *pdEx_Begin, const double *pdEx_End,
                 int nNumPointT_In);

int _DefaultNumX(int nNumPointX_In);

const char *_init_Grid(int nInitNumPointX1, double dDeltaX1,
                       int nInitNumPointX3, double dDeltaX3, int nNumPointT,
                       double dT, _Grid *pGrid);

const char *_init_amort_grid(
    const double *pdVarX1_Begin, // =_alloca(nNumEx*sizeof(double));
    const double *pdVarX1_End,   // =_alloca(nNumEx*sizeof(double));
    const double *pdVarX3_Begin, // =_alloca(nNumEx*sizeof(double));
    int nBase_NumX1, int nBase_NumX3, double dEx_Back, int nNumT,
    /// results
    _Grid *pGrid, int *pnSizeX3_Begin, int *pnSizeX1_Begin);