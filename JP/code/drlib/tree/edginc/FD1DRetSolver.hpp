//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRetSolver.hpp
//
//   Description : one factor finite difference algorithm
//
//   Author      : Xiaolan Zhang
//               : Ning Shen
//   Date        : Apr, 2005
//
//----------------------------------------------------------------------------

#ifndef FD1D_RET_SOLVER_HPP
#define FD1D_RET_SOLVER_HPP
#include "edginc/FD1DRet.hpp"

DRLIB_BEGIN_NAMESPACE

//----------------------------------------------------------------------------

/** FD1DRetSolver algorithm class that supports FD1DRetSolver */

class TREE_DLL FD1DRetSolver : public virtual FDModel::IFDSolver{


public:

    FD1DRetSolver(FD1DRet* engine);
    virtual ~FD1DRetSolver();

    friend class FD1DRetStochGarfSolver;

    /** solve backward or forward induction through all steps. returns price*/
    virtual void roll();

    //-----------------------------------------------------------------------------------

private: 
    //  ---------------------------  local data -------------------------                         

    // contains one slice for each combination of dimensions
    vector< TreeSliceEQSP > sliceCache;

    /**--------------------------------------------------
        coeff:
        dimension 1 size =numCoeff of PDE (cont. form)
        dimension 2 size = xNum  
    --------------------------------------------------*/

    int numCoeff;
    double** coeff; 

    /**FD params copying from model*/
    int xNum;

    /** -----------------------------------------------------------------------------
        space steps, copied from engine->v_dxM
        if const grid, FD scheme uses v_dx[1]
     -----------------------------------------------------------------------------*/

    vector<double> v_dx;

    /**-------------------------------------------
        for 2 factor solver
        1 for LOD, no cross term
        2 for operator splitting (called ADI advanced in the codes)
    
        for 1 factor Solver
        solveMethod = DEFAULT, only 1 methode is proposed
    --------------------------------------------*/

    FD1DRet::TFdSolveType solveMethod; 

    /**-------------------------------------------------------------------
        dim1: nb of products
        dim2: nb of values
        dim3: xNum*yNum with yNUm = 1 in 1D
    -------------------------------------------------------------------*/

    double*** mSourceX;

    /**---------
        PDE coef
    ---------*/

    /**----------------------------------------------------------------
        setup the coefficients of the PIDE needed to be solved
        Assume the most general PDE forms
        PDE 2 factors
        U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = g;
        coeffs are functions of t, x, y
    ----------------------------------------------------------------*/

    /**----------------------------------------------------------------
        PDE 1 factor is
        U_t + a*U_x + c*U_xx + f*U + Jump(x) = g;
        coeffs are functions of t, x
    ----------------------------------------------------------------*/

    vector<vector<double> > a;   
    vector<vector<double> > c;
    vector<vector<double> > f;
    vector<vector<double> > g;

    //  ---------------------------  FD Scheme Coeff--------------------                         

    /**--------------------------------------------------------------------
        discrete FD equation
        alphaX_i * U_(i-1)^(n) + betaX_i * U_i^(n) + gammaX_i * U_(i+1)^(n) =
        aX * U_(i-1)^(n+1) + bX_i * U_i^(n+1) + cX_i * U_(i+1)^(n+1) 
    --------------------------------------------------------------------*/
    vector<double> alphaX;
    vector<double> betaX;
    vector<double> gammaX;
    vector<double> aX;
    vector<double> bX;
    vector<double> cX;

    //  ---------------------------  local ft --------------------                         

    void rollOneStep(int step, int rollDirection);

    void nextStepFdBoundary(
        int step,
        int* bot1, int* top1,
        int k, double pv_x,
        const TreeSliceEQ & mcurrP,
        const TreeSliceEQ & mlastP );

    void updateFdBoundary(
        int step,
        int* bot1, int* top1,
        int pStart, int pEnd,
        int k, double pv_x,
        const TreeSliceEQ & mcurrP,
        const TreeSliceEQ & mlastP );

    void interpBoundaryUp(int t, double* p);

    void interpBoundaryLow(int l, double* p);

    void forceNonNegative( int k, const TreeSliceEQ & price );


    /**--------------------------------------------------------------------
        update solver from model, FD Grid boundaries..  
        be careful, if we decide to let solver chg dx, dy at it,
        FD grid between two steps might be diff.
        so need to update solver,
        for now, solver chg dx only in the case of VAR_GRID in barriers cases
    ---------------------------------------------------------------------*/
    void updateSolverInfo(int step, int idx, int rollDirection);

    void updatePayoffIndex(int step, FDProductSP payoffIndex);

    /**--------------------------------------------------------
        only used in VAR_GRID for now,
        copy back some info. 
    --------------------------------------------------------*/

    void postUpdateSolverInfo(int step);

    /**keep a,b,c,d,e,f, to make code easy to read*/
    void passCoeff(int idx, int bot1, int top1);

    int getfdSeg(int step);

    /** changing time line segment of a different density */
    bool switchSegment(int step, int idx, int currSeg, int RollDirection);
    
    //-----------------------------------------------------------------------------
    
    /**-----------
        FD scheme
    -----------*/
    /**Dim1: nb of price, Dim2: X+Y, big vector: (xNum )*(yNum)*/

    void euroOneStepwithSource1D(
        int lastIndexTop,
        int* bot1, int* top1,
        int k, int idx,
        double ** mSrcX,
        const TreeSliceEQ & mcurrP,
        const TreeSliceEQ & mlastP);

    void calcCoeffPdeAll_X(const FD1DRet::TFdSolveType solveMethod, int low, int top, int idx );

    void calcSourceAll(const FD1DRet::TFdSolveType solveMethod, int bot1, int top1, 
                                int idx, 
                                double* vSrc, 
                                double* vlastP);

    
    //  ---------------------------  schema 1 --------------------                         
    
    double theta;
    void calcCoeffPdeConstdx_X(int low, int top,int idx);

    /**calc variable FD scheme's coeff*/
    void calcCoeffPdeVardx_X(int low, int top, int idx);

    void calcSource(int bot1, int top1, 
                    int idx, double* vSrc, double* vlastP);


    //---------------------------------------------------------------------------

    //  ---------------------------  special barriers ---------------------------                         

    //---------------------------------------------------------------------------

    
    //  ---------------------------  loacl data for barriers --------------------                         


    /**for barrier*/
    /**fixed grid or non fixed grid*/
    FD1DRet::TFdBarrierMethodType whichBarrierMethod; 

    /**------------------------------------
        created in order to call payoffBCDs
    ------------------------------------*/

    /** insNodePrice will be obtained by linear interpolation*/
    int                  numOfInsertNode; // size = 2* number values
    TreeSliceGeneralSP          insNode;  //size = numOfInsertNode, contain inserted nodes
    /** price relating to inserted points    */
    TreeSliceGeneralContSP insNodePrice;     //nbOfP, numOfInsertNode

    vector<int > insNode_FdIndex;  //size = numOfInsertNode, 
                                   //the corresponding index of insNode in FD grid
    vector<bool > insNodeActive;  //default =false, set to true if there is a barrier, 

    /**----------------------------------------------------------------------
        used to store the inserted old price,
        copy them back after "euroOneStepwithSource1D"
        to insure that at any time, the price vectors reflect 
        the right value at initial Grid
    ----------------------------------------------------------------------*/

    vector<vector<double > > insNodePriceBak;  //0: down node, 1: upNode


    /**------------------------------------------------------------
        set up the value for solver BotDim, TopDim
        for barriers: get barrier levels and looking for the cloest index
        it's called at updateSolverInfo 
    ------------------------------------------------------------*/

    void setFDSolverDims(int step);

    /** calc index cloest to crit. pts 
        and convert bar. based on the choose of chg of variable*/
    void setFDSolverDimsSpecial(int step);

    /**recalc dx based on new info. from prod.*/
    void resetdx(vector<double>& vdx, int kkProd, int pStart, int pEnd  );

    /** get info. of barrier and put in a vector*/
    void setInsNodeInfo(int pStart, int pEnd);

    /** calc coeffs at inserted pts at step*/
    void calcBCDsBefore(
        int step,
        int* bot1, int* top1,
        int pStart, int pEnd, int k,
        const TreeSliceEQ & mlastP );

    void ajustLastPatInsNodes(
        int* bot1, int* top1,
        int pStart, int pEnd, int k,
        vector<vector<double> >& pBackup,
        const TreeSliceEQ & mlastP );


    /** update boundary conditions */
    void nextStepFdBoundaryBar(
        int step,
        int* bot1, int* top1,
        int k, double pv_x,
        const TreeSliceEQ & mcurrP,
        const TreeSliceEQ & mlastP );

    void updateFdBoundaryBar(
        int step,
        int* bot1, int* top1,
        int pStart, int pEnd,
        int k, double pv_x,
        const TreeSliceEQ & mcurrP,
        const TreeSliceEQ & mlastP );

    /**update the current value at bar. for const grid
        here, bar aren't on the grids */
    void calcBCDsAfterFixGrid(
        int step,
        int* bot1, int* top1,
        int pStart, int pEnd, int k,
        const TreeSliceEQ & mcurrP,
        const TreeSliceEQ & mlastP );

    void calcBCDsAfterVarGrid(
        int step,
        int* bot1, int* top1,
        int pStart, int pEnd, int k,
        const TreeSliceEQ & mcurrP,
        const TreeSliceEQ & mlastP );

    
    /** for a VAR_GRID, 
        borrower some index for variable grid, 
        so need to copy back the correct price */
    void copyBackLastPatInsNodes( 
        int* bot1, int* top1,
        int pStart, int pEnd, int k,
        vector<vector<double> >& pBackup,
        const TreeSliceEQ & mlastP );


    /**allocating mem*/
    void init_mem();

    FD1DRet* engine;
};

DRLIB_END_NAMESPACE
#endif
