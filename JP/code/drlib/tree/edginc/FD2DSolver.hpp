//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DSolver.hpp
//
//   Description : two factor finite difference algorithm
//
//   Author      : Ning Shen
//               : Xiaolan Zhang
//
//   Date        : November 29, 2004
//
//----------------------------------------------------------------------------

#ifndef FD2DSOLVER_HPP
#define FD2DSOLVER_HPP

#include "edginc/FD2D.hpp"
#include "edginc/FD2DSVCJ.hpp"
#include "edginc/FD2DSolverJumps.hpp"

DRLIB_BEGIN_NAMESPACE

//----------------------------------------------------------------------------
//class FD2DSolverJumps;

/** FD2DSolver algorithm class that supports FD2DSolver */
class TREE_DLL FD2DSolver : public virtual FDModel::IFDSolver{
public:

    FD2DSolver(FD2D* engine);
    virtual ~FD2DSolver();

    /** solve backward or forward induction through all steps. returns price*/
    virtual void roll();

    //-----------------------------------------------------------------------------------

private: 
    //  ---------------------------  local data -------------------------                         

    // temporary slices to perform DEV into
    vector< vector< TreeSliceSP > > tempSlices;

    /**--------------------------------------------------
       coeff:
       dimension 1 size =numCoeff of PDE (cont. form)
       dimension 2 size = xNum  
       dimension 3 size = yNum 
       --------------------------------------------------*/
    
    int numCoeff;
    double*** coeff; 
    
    /**FD params*/
    int xNum;
    int yNum;
    
    /** -----------------------------------------------------------------------------
    space steps, copied from engine->v_dxM, v_dyM
        if const grid, FD scheme uses v_dx[1], v_dy[1]
        -----------------------------------------------------------------------------*/
    
    vector<double> v_dx;
    vector<double> v_dy;
    
    /**----------------------------------------------------------
       1 for LOD, no cross term
       2 for operator splitting (called ADI advanced in the codes)
       -----------------------------------------------------------*/
    
    FD2D::TFdSolveType solveMethod; 
    
    //int solveMethod; 
    
    /**-------------------------------------------------------------------
       dim1: nb of products
       dim2: nb of values
       dim3: xNum*yNum (X direction is vertical, Y direction is horizontal)
       -------------------------------------------------------------------*/

    double***  mSourceX;
    double***  mSourceY;
    
    /**---------
       PDE coef
       ---------*/
    
    /**----------------------------------------------------------------
       setup the coefficients of the PIDE needed to be solved
       Assume the most general PDE forms
       PDE 2 facrors
       U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
       coeffs are functions of t, x, y
       ----------------------------------------------------------------*/
    
    vector<vector<vector<double> > > a;   
    vector<vector<vector<double> > > b;
    vector<vector<vector<double> > > c;
    vector<vector<vector<double> > > d;
    vector<vector<vector<double> > > e;
    vector<vector<vector<double> > > f;
    
    //  ---------------------------  factor 1--------------------                         
    
    /**--------------------------------------------------------------------
       discrete FD equation at direction X (similar for Y)
       alphaX_i * U_(i-1)^(n) + betaX_i * U_i^(n) + gammaX_i * U_(i+1)^(n) =
       aX * U_(i-1)^(n+1) + bX_i * U_i^(n+1) + cX_i * U_(i+1)^(n+1) + source
       --------------------------------------------------------------------*/
    
    vector<double> alphaX;
    vector<double> betaX;
    vector<double> gammaX;
    vector<double> aX;
    vector<double> bX;
    vector<double> cX;

    vector<double> alphaY;
    vector<double> betaY;
    vector<double> gammaY;
    vector<double> aY;
    vector<double> bY;
    vector<double> cY;

    //  ---------------------------  local ft --------------------                         
    
    void rollOneStep(int step, int rollDirection); 
    
    void nextStepFdBoundary(int step,
                            int* low1, int* top1, int* low2, int* top2,
                            int pStart, int pEnd, double pv_x,
                            const vector< TreeSliceSP > & mcurrP,
                            const vector< TreeSliceSP > & mlastP);
    
    //lastIndexTop will be the same for these 4 vectors low1, top1, low2, top2
    void updateFdBoundary(int step,
                          int* low1, int* top1, int* low2, int* top2,
                          int pStart, int pEnd,
                          const vector< TreeSliceSP > & mcurrP,
                          const vector< TreeSliceSP > & mlastP);
    
    /**------------------------------------------------------
       update solver from model, FD Grid boundaries..  
       be careful, if we decide to chg dx, dy at it,
       FD grid between two steps might be diff.
       so need to update solver,
       for now, model chg dx, dy only in the case of barriers
       ------------------------------------------------------*/
    //need rollDirection when adding barriers special treatment, so, keep it here
    void updateSolverInfo(int step, int idx, int rollDirection);
    
    /** keep a,b,c,d,e,f, to make code easy to read*/
    void passCoeff(int step, int idx, int low1, int top1, 
                   int low2, int top2);
    
    /** changing time line segment of a different density */
    bool switchSegment();
    
    //*  useful methods /moved to FD2D
    int getfdSeg(int step);
    
    //-----------------------------------------------------------------------------
    
    /**----------------
       ADI or Splitting
       ----------------*/
    
    /**Dim1: nb of price, Dim2: X (vertical direction), Dim 3: Y*/
    void euroOneStepwithSource2D(int lastIndexTop, int* low1, int* top1,
                                 int* low2, int* top2,
                                 int pStart, int pEnd, int idx,
                                 double** mSrcX,
                                 double** mSrcY,
                                 const vector< TreeSliceSP > & mcurrP,
                                 const vector< TreeSliceSP > & mlastP);
    
    
    void calcCoeffPdeAll_X(const FD2D::TFdSolveType solveMethod, int low, int top, 
                           int idx, int index_y );
    
    void calcCoeffPdeAll_Y(const FD2D::TFdSolveType solveMethod, int low, int top, 
                           int idx, int index_x );
    
    void calcSourceAll(const FD2D::TFdSolveType solveMethod, int low1, int top1, 
                       int low2, int top2, int idx, 
                       double* vSrcX, double* vSrcY, double* M_price);
    
    //  ---------------------------  schema 1 --------------------                         
    
    /**solveMethod = 1, ADI without cross term*/
    double theta;
    void calcCoeffPdeADI_X(int low, int top,int idx, int index_y );
    void calcCoeffPdeADI_Y(int low, int top, int idx, int index_x );
    
    void calcSource(int low1, int top1, int low2, int top2,
                    int idx, double* vSrc, double* M_price);
    
    //  ---------------------------  schema 2 --------------------                         
    
    /**solveMethod =2, advanced ADI with cross term*/
    void calcCoeffPdeAdvancedADI_X(int low, int top, int idx, int index_y);
    void calcCoeffPdeAdvancedADI_Y(int low, int top, int idx,  int index_x);
    
    void calcSourceAdvancedADI_X(int low1, int top1, int low2, int top2, 
                                 int idx, double* vSrc, double* M_price);
    void calcSourceAdvancedADI_Y(int low1, int top1, int low2, int top2, 
                                 int idx, double* vSrc, double* M_price);      
    
//----------------------------------------------------------------------------

    FD2D* engine;
    //FD2DSolverJumpsSP fd2fsolverjumps;    
    
    //allocating mem
    void init_mem();
};

DRLIB_END_NAMESPACE
#endif
