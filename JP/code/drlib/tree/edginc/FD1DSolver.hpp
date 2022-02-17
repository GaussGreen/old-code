#ifndef QLIB_FD1D_SOLVER_HPP
#define QLIB_FD1D_SOLVER_HPP

#include "edginc/FD1D.hpp"

DRLIB_BEGIN_NAMESPACE

class FD1D;

class TREE_DLL FD1DSolver : virtual public FDModel::IFDSolver
{
public:
    FD1DSolver( const FD1D & model );

private: 
    // FDModel::IFDSolver implementation
    virtual void roll();

    void rollOneStep();

    void preCalc();

    void calcCoeffPde();
    void calcCoeffPdeFix();
    void calcCoeffPdeVar();

    void calcSource();

    void solveOneStep(
        double pv, bool nonNegative,
        const TreeSliceEQ & prevSlice,
        TreeSliceEQ & currSlice );

    /**-------------------------------------------------------------------------
        setup the coefficients of the PIDE needed to be solved
        Assume the most general PDE forms

        PDE 1 factor is
        U_t + a*U_x + c*U_xx + f*U + Jump(x) = g;
        coeffs are functions of t, x

        discrete FD equation
        aX0_(i) * U_(i-1)^(n)   + bX0_(i) * U_(i)^(n)   + cX0_(i) * U_(i+1)^(n) =
        aX1_(i) * U_(i-1)^(n+1) + bX1_(i) * U_(i)^(n+1) + cX1_(i) * U_(i+1)^(n+1)

        -------------------------------------------------------------

        PDE 2 factors
        U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = g;
        coeffs are functions of t, x, y

    --------------------------------------------------------------------------*/

    // coefficients retrieved from the model at each step
    TreeSliceEQRef a[ 2 ];
    TreeSliceEQRef c[ 2 ];
    TreeSliceEQRef f[ 2 ];
    TreeSliceEQRef g[ 2 ];

    // solver scheme coefficients
    TreeSliceEQRef aX[ 2 ];
    TreeSliceEQRef bX[ 2 ];
    TreeSliceEQRef cX[ 2 ];

    // solver source slice
    TreeSliceEQRef xSrc;

    // contains one slice for each combination of dimensions
    vector< TreeSliceEQRef > sliceCache;

    // solver scheme theta between 0 and 1: 0 - fully explicit, 1 - fully implicit
    double theta;

    // current step
    int step;

    // current solver range, for convinience
    int bot1, top1;

    const FD1D & model;

    // products retrieved from the model
    const FDProductArray & products;

    // variable slice retrieved from the model
    const TreeSliceEQ & xVar;
};

DRLIB_END_NAMESPACE

#endif
