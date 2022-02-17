#ifndef QLIB_FD1D_HPP
#define QLIB_FD1D_HPP

#include "edginc/LatticeModelEQ.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1D : public LatticeModelEQ
{
public:
    static const CClassConstSP TYPE;
    static void load( CClassSP & type );

protected:
    FD1D( const CClassConstSP & type );

    virtual void validatePop2Object();

    virtual void retrieveFactor();
    virtual void initModel();
    virtual void finaliseModel( CControl * control );

    /** insert nodes helper functions */
    bool insertNode( int step, int index, double x, double x1, double x2 );
    void insertNodesBelow( int count, int step, int index, double x, double x1, double x2 );
    void insertNodesAbove( int count, int step, int index, double x, double x1, double x2 );

    /** product should call this function in 'preCalc' to
        add critical 'level' of 'kind' to 'target' slice based on 'base' slice,
        optionally specifying 'value' of 'target' at the 'level' */
    virtual void addCriticalLevel(
        int step,
        const TreeSlice & base,
        double level,
        TreeSlice & target,
        LevelKind kind = LEVEL_GENERIC,
        double value = DBL_MAX );

private:

    /** accept or reject factors */
    virtual bool acceptFactor( IMarketFactor * factor );

    /** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
    virtual double getPrice0( const TreeSlice & price ) const; 

    /** override FDModel::makeProduct to retrieve spot index product */
    virtual FDProductSP makeProduct( const IProdCreatorSP & creator );

    /** create a solver */
    virtual IFDSolverSP createSolver();
    friend class FD1DSolver;

    /**-----------------------------------------------------
        calculate the stuff which specifies the grid 
        set up the FD boundaries, and calc space step, can be const or var.
        shouldn't be called if SameGrid = true
    /-----------------------------------------------------*/
    virtual void initGrid();

    /** initializes value based on var */
    virtual void initValue( int step ) const;

    /** updates value based on var */
    virtual void updateValue( int step ) const;

    /**-----------------------------------------------------
        each derived model can set diff. boundaries 
        based on the dynamics of underlying,
        alpha is the input truncation1D (nb of std) 
        outLowB, outUpB are low and up boundaries for FD 
    -----------------------------------------------------*/
    virtual void setFdBounds(
        double alpha,
        double & outLowB,
        double & outUpB ) const = 0;

    /**----------------------------------------------------------
        calculate pde operator coefficients 
        extract coeff computed from the market data for the PIDE.
        be aware that the structure and coefficients for forward 
        and backward equation can be very different.
        forward  equation not available yet
    ---------------------------------------------------------*/
    /**--------------------------------------------------------------
        setup the coefficients of the PIDE needed to be solved
        Assume the most general PDE for 2 factors  is 
        U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
        this ft is called at each time step, so they're time depd.
    ----------------------------------------------------------------*/
    /**---------------------------------------------------------------
        PDE 1 factor is
        U_t + a*U_x + c*U_xx + f*U + Jump(x) = g;
    ----------------------------------------------------------------*/
    virtual void pdeCoeff(
        int step,
        TreeSliceEQ & a,
        TreeSliceEQ & c,
        TreeSliceEQ & f,
        TreeSliceEQ & g ) const = 0;

    /** convert input strings into enum variables */
    void convertInputStrings();

    /** diff. methods for set dx */
    void setVarSpaceSteps( int bot1, int top1 );

protected:
    CAssetConstSP underlying;

    /** fwd stock at each step */    
    DoubleArray stepForward;

    typedef enum{ NONE, STD } VarStepSpacing;

    /** ------------------------------------------------------------------
        truncation1D : num of stdev for truncation for each dimension
        dim1 : original dimension of Und1 size, need to be odd nb
        stepsPerYear: num of steps per year, actual steps 
                    created in TimePts may be different
    -------------------------------------------------------------------*/
    double truncation1D;
    int    dim1;
    int    stepsPerYear;

    bool isVariableGrid;

    // chached solver
    IFDSolverSP solver;

    /** registration */
    ChangeOfVar changeOfVar;
    string changeOfVarStr;

    ///** -------------------------------------------------------------
    //    control which method to be used for setting var space step 
    //    NONE: const dx
    //    STD: only based on std
    //    3: more barriers such as in ladder... to be added 
    //-------------------------------------------------------------*/
    VarStepSpacing varStepSpacing;
    string varStepSpacingStr;

    int nodesPerCritLevel;
    string DEBUG_DumpToFile;
};

DRLIB_END_NAMESPACE

#endif
