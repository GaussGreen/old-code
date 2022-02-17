/*===================================================================
** PROGRAM:  ZFDSolver1D.h
**
** AUTHOR:   Milan Kovacevic April 04, 2001.
**
** CONTAINS: Class for FD 1D solver
**=================================================================*/

#ifndef ZFDSOLVER1D_H
#define ZFDSOLVER1D_H

//#include <ospace/std/string>
//#include <zcurve.h>
//#include <zstrgvec.h>
//#include <global.h> // for jdate
//#include <matrix.h>

#include "edginc/config.hpp"
#include "edginc/FDUtils.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////
// ZFDSolver1D //
//////////////////

class TREE_DLL ZFDSolver1D
{
protected:
    double              dx;
    int                 m_iMaxMem;      // maximum memory allocated
    int                 m_iMax;         
    double*             J1;
    double*             m_rhs;
/*
//      double**        V;
double**        Vnew;
*/
    double*             spot1;
    double*             left1;
    double*             diag1;
    double*             right1;
    double*             gam;
    int                 m_nonNegativeValues;
    
    // Precompute for performance
    double*             JdriftTerm;
    double*             JvolTermDown;
    double*             JvolTermUp;
    
public:
    ZFDSolver1D();
    
    
    ZFDSolver1D(                                                                // COPY constructor
        const ZFDSolver1D&);                                    // (I) Original object
    
    int
    allocateMem(int iMax);
    
    void
    resetMem();
    
    void
    deleteMem();
    
    int
    getGridSpots(
        double          S1Min,
        double          S1Max,
        int                     iMax,
        double*         spots,  // (O) spot vector to be populated
        double*         dx,             // (O) physical step size of the grid
        int                     coordinateType=0,
        double          K=1.);
    
    int
    createGrid(
        double          S1Min,
        double          S1Max,
        int                     iMax,
        int                     coordinateType=0,
        double          K=1.);
    
    inline double*
    getSpotsPtr() { return spot1;}
    
    int
    solve(
        double          upBarrier1,
        double          upBC1,
        double          downBarrier1,
        double          downBC1,
        double          divYld1,
        int             iMax,
        int             iMin,
//              double*         spot1,  
        double          r,
        double*         sigma1,
        double*         cs,
        double          coupon,         // Continuous coupon
        double          dt,
//              double          s1Min,
//              double          s1Max,
        double*         V,
        double*         Vnew,
        int                     american=0,
        double*         payoffMin=0,
        double          toleranceSOR=CHEM_EPSILON9,
        double          payoffMax=0,
        double          theta=0.5,              // Discretization method: 0 - explicit, 0.5 - Crank - Nicholson, 1 - implicit, 2 - three time level
        double*         Vold=0,              // Old option prices for the three time level method
        int         varMethod=0
        );
    
    int
    solveBarrier(
        int                     m,
        double          r,
        double          divYld1,
        double*         sigma1,
        double*         cs,
        double          coupon,         // Continuous coupon
        double          dt,
        double*         V,
        double*         Vnew,
        double          upBarrierNew=-1,
        double          upPayoutNew=0,
        double          upPayoutDeltaNew=0,
        double          upBarrierOld=-1,
        double          upPayoutOld=0,
        double          upPayoutDeltaOld=0,
        double          downBarrierNew=-1,
        double          downPayoutNew=0,
        double          downPayoutDeltaNew=0,
        double          downBarrierOld=-1,
        double          downPayoutOld=0,
        double          downPayoutDeltaOld=0,
        int                     american=0,
        double*         payoffMin=0,
        double          toleranceSOR=CHEM_EPSILON9,
        double          payoffMax=0,
        double          theta=0.5,              // Discretization method: 0 - explicit, 0.5 - Crank - Nicholson, 1 - implicit, 2 - three time level
        double*         Vold=0,         // Old option prices for the three time level method
        int         varMethod=0
        );
    
    int triDiag1D(
        double*         a, 
        double*         b, 
        double*         c, 
        double*         r, 
        double*         u,
        int                     n);
    
    int
    sor(
        double*         a,
        double*         b,
        double*         c,
        double*         d,
        double*         payoff,
        double*         u,
        int                     m,
        double          tolerance=CHEM_EPSILON9
        );
    
    int sorMinMax(
        double*         a,
        double*         b,
        double*         c,
        double*         d,
        double*         payoffMin,
        double          payoffMax,
        double*         u,
        int                     n,
        double          tolerance=CHEM_EPSILON9);
    
    int isNonNegative() const {return m_nonNegativeValues; } 
    void forceNonNegative(int   nonNegative=TRUE) {m_nonNegativeValues = (nonNegative ? TRUE : FALSE);}
    
    // Destructor
///     virtual 
    ~ZFDSolver1D();
    
/*
  // Assignment
  int
  setFXData(
  int                           payoutZCID,
  int                           zeroCurveID,
  int                           fxVolID,
  double                        fxVol,
  const char*           fxCorrID,
  double                        fxCorr,
  const ZMarket*        mkt=NULL);
*/
    
    //ZFDSolver1D & operator=(const ZFDSolver1D &);
};

DRLIB_END_NAMESPACE
#endif
