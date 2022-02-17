// prevent multiple inclusions
#pragma once

////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                // in the debug information"

// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "AmortMidatADI.h"
#include "AmortMidatADIGrid.h"
#include "AmortMidatADIUtils.h"

/// signature of the integrand needed by Price_Contract
typedef void (*PFN_INTEGRAND)(
    const double *, /// dT
    const double *, // pdX1_Begin
    const double *, // pdX1_End
    const double *, // pdX3_Begin
    const double *, // pdX3_End
    double **,      // Grid - Row index [nX1_Begin      ,...)
                    //     , Grid-Col index [nX3_Begin      ,...)
    int,            // nX1_Begin      ,
    int,            // nX3_Begin      ,
    void *);

//-------------------------------------------------------------------------------------------------------
//	Description: computes the local vol of X1 in LGM2F
//
//	NB:	local vol is a function of time only;sgnature required for the
// PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _SIGMA1_Sqrd_(const double *pdX1, const double *pdX3,
                               const double *pdT, void *pArg);
//-------------------------------------------------------------------------------------------------------
//	Description: computes the local vol of X3 in LGM2F
//
//	NB:	local vol is a function of time only;sgnature required for the
// PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _SIGMA3_Sqrd_(const double *pdX1, const double *pdX3,
                               const double *pdT, void *pArg);

//-------------------------------------------------------------------------------------------------------
//	Description: computes the drift of X1 in LGM2F in the backward equation
//
//	NB:	local vol is a function of time only;sgnature required for the
// PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _MU1_BACKWARD_(const double *pdX1, const double *pdX3,
                                const double *pdT, void *pArg);
//-------------------------------------------------------------------------------------------------------
//	Description: computes the drift of X3 in LGM2F in the backward equation
//
//	NB:	local vol is a function of time only;sgnature required for the
// PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _MU3_BACKWARD_(const double *pdX1, const double *pdX3,
                                const double *pdT, void *pArg);
//-------------------------------------------------------------------------------------------------------
//	Description: computes the drift of X1 in LGM2F in the forward equation
//
//	NB:	local vol is a function of time only;sgnature required for the
// PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _MU1_FORWARD_(const double *pdX1, const double *pdX3,
                               const double *pdT, void *pArg);
//-------------------------------------------------------------------------------------------------------
//	Description: computes the drift of X3 in LGM2F in the forward equation
//
//	NB:	local vol is a function of time only;sgnature required for the
// PDE engine
//
//	Returns : as indicated
//-------------------------------------------------------------------------------------------------------
double __stdcall _MU3_FORWARD_(const double *pdX1, const double *pdX3,
                               const double *pdT, void *pArg);

//-------------------------------------------------------------------------------------------------------
//	Description:  Solve for the "contract" U(t      ,x)
//						U solves U(t      ,x1
//,x3)=E[Integrand(T
//    ,X(T))|X1(t)=x      ,X3(t)=x3]
//
//						NB: U traverses from *pdT_Begin
//to
//*(pdT_End-1) 						To solve for the backward
//equation      , pdT should be decreasing and for the forward equation      ,
// increasing.
//
//
//	Returns : U(0      ,0)
//-------------------------------------------------------------------------------------------------------
const char *Price_Contract(
    /// time direction      , n = pdTKnot_End- pdTKnot_Begin-1 intervals in
    /// total      , i.e. pdTKnot_Begin[0]-pdTKnot_Begin[1]      , ...
    /// pdTKnot_Begin[n-1]-pdTKnot_Begin[n]
    const double *pdTKnot_Begin, const double *pdTKnot_End,
    // Grid
    _Grid *pGrid,
    // size of the Grid corresponding to the intervals
    const int *pnSizeX1_Begin, const int *pnSizeX3_Begin,
    // Functor for integrand      ,i.e. "payoff" function(t      ,X1(t) ,X3(t))
    // IMPORTANT : see .h for specs of the payoff function
    PFN_INTEGRAND _INTEGRAND_,
    // Auxillary arguments to be passed to _INTEGRAND_
    void *pComm_INTD,
    // functor to update vol of X1
    // functor must accept the argument seqence (X1      ,X3      ,T)
    PFN_COEFF_3DPDE _S1_Sqrd_,
    // functor to update vol of X3
    // functor must accept the argument seqence (X3      ,X1      ,T)
    PFN_COEFF_3DPDE _S3_Sqrd_,
    // functor to update drift of X1
    // functor must accept the argument seqence (X1      ,X3      ,T)
    PFN_COEFF_3DPDE _M1_,
    // functor to update drift of X3
    // functor must accept the argument seqence (X3      ,X1      ,T)
    PFN_COEFF_3DPDE _M3_,
    // Auxillary arguments to be passed to _S1_Sqrd_      ,_S3_Sqrd_      ,_M1_
    // ,_M3_
    void *pComm_PDE,
    // result
    double *pdResult);
