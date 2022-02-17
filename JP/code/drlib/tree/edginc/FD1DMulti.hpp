//----------------------------------------------------------------------------
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DMulti.hpp
//
//   Description : 1-factor finite difference engine to solve multiple PDEs, which have the form of 
//
//					U(i,t,x)_t + a(i)*U(i,t,x)_x + c(i)*U(i,t,x)_xx + f(i)*U(i,t,x) 
//						
//						+ Sum_j{q(i,j)*U(j,t,x + jump(i,j))} = g(i)
//
//				   for i = 1, ..., I
//
//   Author      : Zhijiang Huang
//
//   Date        : Aug 08, 2006
//----------------------------------------------------------------------------
#ifndef FD1DMulti_HPP
#define FD1DMulti_HPP

#include "edginc/LatticeModelEQ.hpp"
#include "edginc/TreeSliceEQ.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1DMulti : public LatticeModelEQ {
public:
    static const CClassConstSP TYPE;
    static void load( CClassSP & type );

	//temporary, to make barrier option work, to review
    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end) = 0;

	double getTruncationStd();

	/** get date at step on the time line */
    double getTrdYrFrac( int step ) const;

protected:
    FD1DMulti( const CClassConstSP & type );

    virtual void validatePop2Object();

    virtual void retrieveFactor();
    virtual void initModel();
    virtual void finaliseModel( CControl * control );

    /** product should call this function in 'preCalc' to
        add critical level at 'value' of 'slice' */
    virtual void addCriticalLevel( int step, const TreeSlice & slice, double value );

private:
    /** accept or reject factors */
    virtual bool acceptFactor( IMarketFactor * factor );

    /** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
    virtual double getPrice0( const TreeSlice & price ) const; 

    /** override FDModel::makeProduct to retrieve spot index product */
    virtual FDProductSP makeProduct( const IProdCreatorSP & creator );

    /** create a solver */
    virtual IFDSolverSP createSolver();
    friend class FD1DMultiSolver;

    /**-----------------------------------------------------
        calculate the stuff which specifies the grid 
        set up the FD boundaries, and calc space step, can be const or var.
        shouldn't be called if SameGrid = true
    /-----------------------------------------------------*/
    virtual void initVarGrid();

    /**-----------------------------------------------------
        initializes value based on var
    /-----------------------------------------------------*/
    virtual void initValue( int step ) const;

    /**-----------------------------------------------------
        updates value based on var
    /-----------------------------------------------------*/
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
        setup the coefficients of the PIDEs needed to be solved
        Assume the most general PDE is 
		U(i,t,x)_t + a(i)*U(i,t,x)_x + c(i)*U(i,t,x)_xx + f(i)*U(i,t,x)	+ Sum_j{q(i,j)*U(j,t,x + jump(i,j))} = g(i)
		for i = 1, ..., I
    ----------------------------------------------------------------*/
	virtual void pdeCoeff(
        int step,
        DoubleArray & a,
        DoubleArray & c,
        DoubleArray & f,
        DoubleArray & g,
		DoubleArray & q,
		DoubleArray & jump) const = 0;

    /** convert input strings into enum variables */
    void convertInputStrings();

    /** diff. methods for set dx */
    void setVarSpaceSteps();

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
		nPdes: num of PDEs
    -------------------------------------------------------------------*/
    double	truncation;
    int		dim;
    int		stepsPerYear;
	int		nPdes;
	int		iState;

    bool isVariableGrid;

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

    string DEBUG_DumpToFile;
};

DRLIB_END_NAMESPACE

#endif
