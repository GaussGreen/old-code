//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CountingTree.hpp
//
//   Author      : Charles Morcom, July 31, 2006
//
//   Description : An n-dimensional tree where each factor can have a number
//                 of counted states as well as a rate/spread. The counting
//                 states are suitable for representing, e.g., the number of
//                 defaults as a state-variable.
//                 Features include:
//                 - n factors (limited by memory and speed)
//                 - counts discrete states (e.g. default/no-default states)
//                 - 'real' factors are mapped O-U, or CEV
//                 - diffusion variables can have simple local jumps
//                 - permanent contagion handled by transition prob. mapping function
//
//                 This is not yet complete. I have written all major sections
//                 except the clipping 
//                 - computeInnerLimits()
//                 and the slice resizing during the forward induction 
//                 - resizeGlobalLimits()
//                 the remaining code has all been written, but not yet 
//                 tested.
//                 The only major feature that has not been implemented at all
//                 yet is transient contagion (jumps in diffusion variables
//                 caused by transitions in counting-state variables). You'll
//                 find some sections of code commented out, but some more
//                 thinking is needed.
//----------------------------------------------------------------------------

/*============================================================================
 * TODO:
 * - write clipping code
 * - resize slices for change in limits (resizeGlobalLimits())
 * - initializers and constructors
 *==========================================================================*/

#ifndef QR_COUNTINGTREE_HPP
#define QR_COUNTINGTREE_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/FDModel.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL CountingTree : public FDModel {
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    /** we re-use same slice for each time step in all trees, so this index is always 0 */
    virtual int getSliceIndex( int step ) const  { return 0;}

	// override FDModel::initModel() so can precompute various things
	virtual void initModel();

    virtual ~CountingTree();

protected:

    CountingTree(const CClassConstSP &type = TYPE);
    
    /*=========================================================================
     * METHODS TO ADD/MODIFY MARKET FACTORS
     *=======================================================================*/
    /** Add a CEV market factor: the 'x' variable modeled is x = y^(1-beta)/(1-beta),
        where y is CEV with parameter beta, given mean-reversion, and time-dependent
        mean-reversion level determined during calibration. */
    void addCEVMarketFactor(double meanReversionRate, double cevBetaParameter);

    /** Add a single dimension mapped OU market factor (usu. 2Q) */
    void addOUMarketFactor(double meanReversionRate);

    /** Add a multi-dimensional mapped OU market factor */
    void addOUMarketFactor(const vector<double>& meanReversionRates);

private:
    /* This is called internally only */
    void addMarketFactor(bool addF, double mr, bool isCEV, double beta, int f, int n);
protected:

    /** Adds a non-trivial counting-state to the nth dimension driver whose transitions
        depend on the value of market-factor f. The states can vary from 0 (intially) up
        to maxCState. */
    void setCountingState(int n, int f, int maxCState);

    /** Adds local jumps to market factor f (where f has a multi-dimensional driver,
        the jumps are added to the _last_ driver). */
    void setLocalJumps(int f, double jumpHazard, double jumpMeanSize, double jumpWidth);

    /** Set the local diffusion correlations for the drivers (x). These are constant
        for all time. */
    void setDriverCorrelation(const DoubleMatrix& corr);

    /** Set the (time-dependent) driver diffusion volatilities. vol[t][n] is the
        volatility of driver n from t to (t+1). */
    void setFactorVolatility(const DoubleMatrix& vol);

private:
    /** Checks to see that market factors, drivers, counting-states, local jumps
        and dependencies are a AOK. */
    void checkFactorSetUpOK();
protected:

    /** Get slice value at middle of range (for d-dimensions). For O-U d-states, this will be 0.
        For c-dimensions this is the zero state (so not mid-range)*/
    virtual double getPrice0( const TreeSlice& price ) const;

    /*=========================================================================
     * MAPPINGS FROM (c,d)-SPACE TO MARKET FACTOR f-SPACE
     * The functions factorMapping and inverseFactorMapping are virtual. 
     * The default implementations are trivial, and it is not a very good
     * idea to use them.
     * As an implementing class you must make sure 
     * that you implement these functions consistently with whether the mapping
     * is CEV or (usually) 2Q. I cannot implement them at this parent class
     * level in a way which is efficient, since the function needs not just the
     * contagion function alpha(c), but also, in the case of a 2Q mapping, the
     * forwards of the factor in question. Using function pointers or virtual
     * calls to these, while you might think it's more elegant, would be
     * too slow, I suspect.
     * The general idea of the mappings is the following:
     * f = alpha(c).[(1-beta).v]^{1/(1-beta)} for a CEV f, and
     * f = alpha(c).f_bar.2Q(v + hDriftAdjustment) for a 2Q process.
     * The drift adjustment for the CEV process is used inside the transition
     * probability and shift computations, which is why it doesn't appear
     * in the mapping function for the CEV process.
     *=========================================================================
    /* Maps a state (c,v=sum(x)) to the value of the fth factor. To get v, you should
       use vMapping(f, d). For safety, this is pure virtual. */
    virtual double factorMapping(int f, double v, const vector<int>& c) const;
    /* Returns the v = sum(x) value that corresponds to a given market factorVal
       and a given c-state. For safety, this is pure virtual. */
    virtual double inverseFactorMapping(int f, double factorVal, const vector<int>& c) const;
    /* Given v_n = sum(x), what single d-value in dimension dDimension maps to this? 
       This is used to compute the change in d required to produce a required
       change in v. */
    double singleDInverseMapping(int f, double v, int dDimension) const;
    /* Maps a diffusion grid vector, d, to v_f = sum(x), where x = Dd */
    double vMapping(int f, const vector<int>& d) const;

    /*=========================================================================
     * INTERNAL METHODS FOR CHECKING AND COMPUTATION
     *=======================================================================*/
    /* Compute the transition probabilities and associated values between tIdx and
       the next time point */
    void computeTransition();
    /* Compute the transition probabilities from t to (t+1), and also compute
       the tree limits (use during forward induction) */
    void computeTransitionAndLimits();

	/**(D)EV of a valueSlice to the time point tIdx from tIdx+1. The operation is EV if the 
	   discountSlice pointer is null, else DEV. */
	void sliceEV(
        TreeSliceNDimCounting& valueSlice, 
        TreeSliceNDimCounting* discountSlice) const;

    /** Compute the D (and DInverse) matrices which relate d to x, once you know the factor vols,
        time-steps and correlations.*/
    void generateDMatrices();

    /** Check that global limit are big enough to contain inner, mid, and outer limits
        up to the current time. If not, increase, and resize all internal slices to
        reflect changes. */
    void resizeGlobalLimits();
    
    /** Use state-probabilities and outer limits to compute inner limits at current
        time. */
    void computeInnerLimits();
    
    /*=========================================================================
     * PRIMITIVE INPUT VARIABLES
     *=======================================================================*/
    /* MARKET FACTORS AND X->F MAPPINGS 
       Each factor f (there are F of them) is a mapping of a sum of driver 
       x processes.
       f_i = H_i(sum_{k=xForFStart[i]}^{xForFStart[i+1]-1} x_k), c, h) */
    int F; // number of market factors
    vector<int> xForFStart; // what is the start index of the xs that drive f?

    /* DIFFUSION AND LOCAL JUMP x/d DRIVERS */
    int N; // diffusion dimension of tree (N>=F)
    vector<double> localJumpHazard;   // one rate for each dimension: p(jump) = 1-exp(-hazard*dt)
    vector<double> localJumpMeanSize; // p(u)JU + p(d)JD: mean proportional size. This refers to the size of the jump in the associated market factor.
    vector<double> localJumpWidth;    // JU/JD: if up jump is twice down jump this is 2.0
    DoubleMatrix dfCorrelation; // correlations between x-diffusions
    DoubleMatrix dfSpotVolatility; // dfSpotVolatility[t][n] is the spot volatility of x_n between t and t+1.

    
    /* COUNTING STATES */
    vector<bool> hasCState; // true if factor n has a counting process attached
    vector<int> cMax; // the maximum number of c-state transitions in factor n (only relevant if hasCState[n], else foribly set to 0.
    double cTailProbCutOff; // in each time-step, cap c'-c once this cumulative probability has been reached

    /*=========================================================================
     * VARIABLES THAT ARE PRECOMPUTED ONCE ONLY AT INITIALIZATION
     *=======================================================================*/
    int T; // This is the index of the last date in the timeline so tIdx = 0,...,T $unregistered
    /* length contains the same information as timeLine->TradeYrFrac, but is from t to t+1 rather than
       from t-1 to t+1. This is to make it less confusing. This way 'update' at tIdx means looking
       at/computing changes between tIdx and tIdx+1. */
    vector<double> length;  // $unregistered length[tIdx] is the time from tIdx to tIdx+1
    vector<double> lengthD; // $unregistered lengthD[tIdx] is the asset-space 'effective' length from tIdx to tIdx+1 used to compute the mesh-spacing
    vector<int> fForX; // $unregistered market factor that x_i is a driver of.
    DoubleMatrix sqrtDFCorrelation;    // Choleski decomp of diffusionFactorCorrelation $unregistered
    vector<DoubleMatrix> D; // D[tIdx] is the time tIdx matrix that determines the transition from correlated (x) Gaussian drivers to uncorrelated (d) $unregistered
    vector<DoubleMatrix> DInverse; // d = DInverse.x

    // MARKET FACTOR SPECIFICATIONS
    
    vector< vector<double> > marketFactorValues; // slices with values of market factors
    
    // x DRIVER VARIABLE SPECIFICATIONS: OU or CEV with (time) constant mean-reversion
    vector<bool> xIsOU; // is x_i O-U? If not, it is a CEV driver.
    vector<double> xMeanReversion; // mean reversion rate of x_i - for CEV, this is 'k'
    vector<double> xCEVParameter;  // value of 1/2 = CIR, 1 is lognormal, 0 is Gaussian

    /* h[t] is the time-dependent drift adjustment for the mean-reversion level. 
       It is calibrated during forward induction, but can only be done by a concrete
       tree that knows the meaning of the market factors. hDriftAdjustment[f][t] is the
       adjustment that applies to market-factor f at time t. It is needed here
       for the Euler approximation to the CEV SDE's dt increment. */
    vector< vector<double> > hDriftAdjustment;  

    vector<int> transitionHazardFactorIndex; // which factor determines the transition hazard?
    //vector<int> transientContagionDriver; // which c-state change drives transient contagion for this dimension d-state?
    // LOCAL JUMPS: THEY ARE ALL PROPORTIONAL, NOT ABSOLUTE.
    bool hasLocalJumps;            // true if any factors have local jumps
    vector<int> localJumpFactorIndex; // which market (f) factor is affected by the jump?
    
    vector< vector<double> > localJumpUProb;   // p(u) at each time-step, in each dimension
    vector< vector<double> > localJumpDProb;   // p(d) at each time-step, in each dimension
    vector< vector<double> > localJumpNJProb;  // 1-p(u)-p(d) at each time-step, in each dimension
    vector<double> localJumpUFactor;  // JU
    vector<double> localJumpDFactor;  // JD

    /*=========================================================================
     * VARIABLES CALIBRATED DURING FORWARD INDUCTION
     *=======================================================================*/
    /* This defines the inner, mid, outer and global limits for the tree, and
       the current time step. Structured this way so that the slice can access
       this information. */
    TreeSliceNDimCounting::Range range;
    /* Slice of previous time state-probabilities */
    TreeSliceNDimCountingSP lastProbabilities;

    /* Current time-step c-state transition probability slices
       The probability of an increment of dc in dimension i is
       cTransitionProbs[i][dc][offset], where offset is the slice-location
       of the (c,d) state node */
    vector< vector< vector<double> > > cTransitionProbs;
    /* d-state transition probabilities. indices are [i][DTRANS_TYPE][offset] */
    vector< vector< vector<double> > > dTransitionProbs;
    /* Change in d-state induced by changes in c-states: this you need to rework when
       you implement transient contagion. */
    //int*** dStarMapping;

    /* shift[i][offset] is the ith d-coordinate shift from the t-slice to the t+1-slice 
       for the dimension i diffusion: the possible node mappings are shift[i][offset] 
       +1, -1, +0 depending on the trinomial branch. shift[i] is an i-dimensional
       slice (i.e. it depends on all the state-variables up to dimension i. */
    vector< vector<int> > shift;
    /* Offsets to the local jump destinations. As above, localJumpXXDestination[i]
       is an i-dimensional slice. */
    vector< vector<int> > localJumpUDestination;
    vector< vector<int> > localJumpDDestination;

    /* Returns true if (c,d) is within the inner limits at range.tIdx */
    static bool withinInnerLimits(const TreeSliceNDimCounting::Range& range, int sliceDim, 
        const vector<int>& c, const vector<int>& d);

    /*=========================================================================
     * METHODS TO LOOP THROUGH POSSIBLE C-STATE TRANSITIONS AT A GIVEN NODE
     * AND COMPUTE TRANSITION PROBABILITIES
     *=======================================================================*/
    /* The start point for this is the zero-change state, so cPrime=c and dStar=d.
       the probability of this is computed.*/
    void cPrimeDStarIterBegin(
        int             sliceDim,   // slice dimension
        const vector<int>& c,          // input c-state vector
        const vector<int>& d,          // input d-state vector
        int             cdOffset,   // cdOffset of (c,d) into sliceDim slice
        vector<int>&       cPrime,     // output c' c-state
        vector<int>&       dStar,      // output d* d-state w. transient contagion
        double*         transProb   // output p(c->c')
        ) const;
    /* Increments to next possible (c',d*) given (c,d) and computes transition
       probability. When gets to end, sets c'[0]=-1. For the moment, I have
       forbidden transient contagion, to d*=d always. If you try to create
       a tree with transient contagion, you will get an exception. */
    void cPrimeDStarIterNext(
        int             sliceDim, 
        const           vector<int>& c, 
        const           vector<int>& d, 
        int             cdOffset,
        vector<int>&       cPrime, 
        vector<int>&       dStar, 
        double*         transProb
        ) const;

    /*=========================================================================
     * METHODS TO LOOP THROUGH POSSIBLE d*->d' TRANSITIONS AT A GIVEN NODE
     *=======================================================================*/
    /* Given (cPrime, dStar), find the starting point of the accessible set 
       of dPrime in the t+1 slice. Also compute the transition probability and
       the shiftedOffset into the t+1 slice. */
    void dPrimeIterBegin(
        int sliceDim,
        const vector<int>& cPrime, 
        const vector<int>& dStar, 
        int offset,
        vector<int>& dTransState, // have values DTRANS_DU etc.
        vector<int>& dPrime, // d-values in t+1 slice
        double* p, // probability of transition
        int* shiftedOffset // offset into t+1 value slice
        ) const;
    /* Advance to the next dTransState/dPrime possible node in the (t+1) 
       slice, compute probability of transition and the shifted offset into
       the (t+1) slice. When get to end, set dPrime[0]=-1. */
    void dPrimeIterNext(
        int sliceDim,
        const vector<int>& cPrime, 
        const vector<int>& dStar, 
        int offset, // offset into this (t) slice DO YOU NEED THIS?
        vector<int>& dTransState, // have values DTRANS_DU etc.
        vector<int>& dPrime, // d-values in (t+1) slice
        double* p, // probability of transition
        int* shiftedOffset // offset into (t+1) value slice
        ) const;

    /* Computes the slice value array offset, given a set of (c,d) coordinates
       and a slice dimension. c and d coordinates are interleaved into
       p = (d0, c0, d1, c1, ... ) and each are converted into [0, Si] values
       where Si is upperlimit - lowerlimit + 1, by subtracting the lower limits. 
       Note that slice array offsets are controlled by the _global_ limits, not
       those for this particular time-point. */
    static int sliceOffset(
        const TreeSliceNDimCounting::Range& range, int sliceDim, 
        const vector<int>& c, const vector<int>& d);

    static void nodeIterBegin(int limitType, const TreeSliceNDimCounting::Range& range, int sliceDim, vector<int>& c, vector<int>& d);
    static bool nodeIterEnd(int limitType, const TreeSliceNDimCounting::Range& range, int sliceDim, vector<int>& c, vector<int>& d) {return c[0]<0;};
    static void nodeIterNext(int limitType, const TreeSliceNDimCounting::Range& range, int sliceDim, vector<int>& c, vector<int>& d);

    friend class TreeSliceNDimCounting;

};

DRLIB_END_NAMESPACE

#endif
