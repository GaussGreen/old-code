//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1F.hpp
//
//   Description : one factor finite difference engine
//                  U_t + a*U_x + c*U_xx + f*U + Jump(x) = g;
//
//   Author      : Xiaolan Zhang
//
//   Date        : Apr, 2005
//
//----------------------------------------------------------------------------

#ifndef FD1D_RET_HPP
#define FD1D_RET_HPP

#include "edginc/FDUtils.hpp"
#include "edginc/LatticeModelEDR.hpp"
#include "edginc/IndexSpecEDR.hpp"

DRLIB_BEGIN_NAMESPACE

/** initialisation data for setting up the model (2 part)*/
/** add more data */

class TREE_DLL FD1DRet : public LatticeModelEDR
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    
    friend class FD1DRetSolver;
    friend class FD1DRetStochGarfSolver;

    virtual ~FD1DRet();
    FD1DRet(CClassConstSP type);

    /** calc the underlying value */
    virtual void computeUndLevel1(int bot, int top, int step, double*s) const;

    //temporary, to make barrier option work, to review
    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end) = 0;

    /** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
    virtual double getPrice0( const TreeSlice & price ) const; 

    /** add data for fd initialisation, timeline setup */
    //temporary, to chg
    void initSegments(const DateTimeArray& seg, 
                            const IntArray&   dens, 
                            DoubleArray& critSpacePts, //critical pts at space dimension
                            IntArray* isAddedSeg);


    void setNbOfProd(int nProdIn);
    void setMaxNumOfValue(int maxNumOfValueIn);
    void setNumOfInsertNode(int    numOfInsertNodeIn); // size = 2*number of prices in Price array
    int getNumOfInsertNode();

    /** override FDModel::makeProduct to retrieve spot index product */
    virtual FDProductSP makeProduct( const IProdCreatorSP & creator );

    /** create a solver */
    virtual IFDSolverSP createSolver();

    // temp
    bool DEBUG_UseCtrlVar; // $unregistered

    //--------------------------  local methods   ------------------------------------------
protected:
    /** Invoked after instrument has got its market data. */
    virtual void getMarket(const MarketData* market, IInstrumentCollectionSP instruments);

    /** accept or reject factors */
    virtual bool acceptFactor( IMarketFactor * factor );

    /** retrieving market data from MDF */    
    virtual void retrieveFactor();

    /**-----------------------------------------------------
        model initialisation (first part)
        get vol 
        call timeline builder using initData 
        set the dates and FD params at time direction 
        called before initProd 
    -----------------------------------------------------*/

    /**-----------------------------------------------------
        if more than 2 products we restrict:
        only has 1 barrier product and we assume it's the first product 
    -----------------------------------------------------*/

    virtual void initModel();

    /** model parameter validation*/
    /**-----------------------------------------------------
        model initialisation  (second part)
        most memory init should be here
        finalise model initialisation once after initModel() and product init() 
    -----------------------------------------------------*/

    virtual void finaliseModel(CControl*    control);

    /**-----------------------------------------------------
        calculate the stuff which specifies the grid 
        set up the FD boundaries, and calc space step, can be const or var.
        shouldn't be called if SameGrid = true
    /-----------------------------------------------------*/

    virtual void SetSpaceGridData();

    /**-----------------------------------------------------
        each derived model can set diff. boundaries 
        based on the dynamics of underlying,
        alpha is the input truncation1D (nb of std) 
        outLowB, outUpB are low and up boundaries for FD 
    -----------------------------------------------------*/

    virtual     void setFdBounds(double& volForBound,  
                                double alpha, 
                                double& outLowB, 
                                double& outUpB){};


    /** -----------------------------
        one factors level set up 
     -----------------------------*/

    /** s : gridLevel */
    virtual void computeFdGridLevel(int bot, int top, vector<double>& vdx, 
                double lBoundary, double uBoundary, vector<double>& s);

    /** get segment index of a fd*/
    //int getSegSegment(int step);

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

    virtual void pdeCoeff(int step, double** coeff, 
                            int bot1, int top1) = 0;

    //--------------------------  local data   --------------------------

    CAssetConstSP  underlying; // $unregistered

    FDProductSP    payoffIndex; // $unregistered
    FDProductSP    payoffIndexOrig; // $unregistered

    /**general PDE 1F has 4 coeffs, jumps model will be more.*/
    int numCoeffPDE;   // $unregistered

    /** fwd stock at each step */
    //CDoubleArray stepForward1;
    /** fwd stock at each step */    
    vector<DoubleArray>   stepForward; //should replace stepForward1 and stepForward2 $unregistered

    vector<double> gridLevel1; // 1-d grid levels $unregistered
    
    /**-----------------------------------------------------------------------------------
        the following two var tell us FD will be done between which index 
        the last element in the topDim1 contains the biggest index amonge the numOfValues
        the last element in the botDim1 contains the smallest index amonge the numOfValues
    ------------------------------------------------------------------------------------*/
    int** topDim1; // top grid index for Und1 $unregistered
    int** botDim1; // bottom grid index for Und1 $unregistered

    double  lBound1; //Low boundary for factor1 $unregistered
    double  uBound1; //Upper boundary for factor1 $unregistered

    /** -----------------------------------------------------------------------------
        v_dxM is the working vector, may change based on barriers level at step
        v_dxMOrg is the reference vector, don't change with barrier level at step
        initially, v_dxM = v_dxMOrg
        M refer to Model
     -----------------------------------------------------------------------------*/

    /** -----------------------------------------------------------------------------
        space steps
        if const grid, FD scheme uses v_dxM[1]
    -----------------------------------------------------------------------------*/

    vector<double> v_dxM; // $unregistered
    /** a copy of original v_dx */
    vector<double> v_dxMOrg; //only used in Barriers $unregistered

    /** ----------------------------------------------------------------------------
        model save it, otherwise, has some trouble at delete 
        payoff maybe deleted before model destructor
    -----------------------------------------------------------------------------*/

    int nProd; // $unregistered
    int maxNumOfValue; // $unregistered

    /**--------
        flags
    --------*/
    /** -----------------------------------------------------------------------------
        hasJumps : false, solver will not call jumpsolver such as FFT...
    -------------------------------------------------------------------------------*/

    bool hasJumps; // $unregistered
    bool needRecalcUndValue1;  //true, undValue depends on t. If chg var. = log(fwd), this should be set to true $unregistered
                                //or some barrier cases

    /**------------------------------------------------------------------
        User Input
    ------------------------------------------------------------------*/

    /** Different type of FD scheme ADI or ADVANCE_ADI */
    typedef  enum{DEFAULT} TFdSolveType;    
    typedef  enum{FIX_GRID, VAR_GRID} TFdBarrierMethodType;
    typedef  enum{NONE, STD, STD_BAR} TFdSetSpaceStepType;


    ChangeOfVar  whichChgVarDim1;    // $unregistered
    bool hasDiscTerminPDE;   //if yes, PDE includes term -ru; if no, chg var. such that -ru term disappear in PDE $unregistered

    /**------------------------
        It's a solver param, put it here only for temporary input purpose
        olny 1 scheme for FD1DRetsolver is proposed at this moment
    ------------------------*/

    TFdSolveType solveMethod;  // $unregistered

    /** ------------------------------------------------------------------
        StepsPerYear: num of steps per year, actual steps 
                    created in TimePts may be different
        Truncation1D : num of stdev for truncation for each dimension
        Dim1 : original dimension of Und1 size, need to be odd nb
        Dim2 : original dimension of Und2 size
    -------------------------------------------------------------------*/
    double truncation1D;
    int    dim1;
    int    stepsPerYear;

private:
    int    GetStepsPerYear() const;
    double ConvertStoGridVar(int seg, double s);


    //------------------------------------------------------------------------------

    //--------------------------  Barriers -----------------------------------------

    //------------------------------------------------------------------------------


public: 
    /** for barrier, to see where to put it.*/

    FDUtilsBarrierSP*  barrier;   //it's a vector of size = NumOfPrice $unregistered

    /** copy barrier level from the prods*/
    bool getDownBarrier(int iP, double s);
    bool getUpBarrier(int iP, double s);
    bool hasBarrier(int step, int pStart, int pEnd);
//    int     numOfInsertNodeM; // size = 2*number of prices in Price array

protected :

    /** --------------------------------------------------------------------------
        only mv the inserted pts, the remaining pts in the grids stay the same 
        after each time step, all pts are back to original grids 
        a copy of original spot vector,
    --------------------------------------------------------------------------*/

    /** --------------------------------------------------------------------------
        control if we need to do a variable grid calc, including coeffs, recalc S    
        isVariableGrid is set to true if either 
        whichBarrierMethodM == VAR_GRID
        whichSetSpaceStepM = STD or STD_BAR 
    --------------------------------------------------------------------------*/

    bool isVariableGrid;     // $unregistered
    int     numOfInsertNodeM; // size = 2*number of prices in Price array $unregistered

    /** for barriers */
    TreeSliceGeneral::RangeSP insRangeM; // $unregistered
    TreeSliceGeneralSP insNodeM;  //size = nb of node added, only to be consistent with tree, fd don't need 2 slice $unregistered
    /** price relating to inserted points    */
    TreeSliceGeneralContSP insNodePriceM;     //nbOfP, sizeOf insNode, only to be consistent with tree, fd don't need 2 slice $unregistered

    /** --------------------------------------------------------------------------
        convert the barGrid consistent with chg of variable in FD grid
        need to add the new chg of variable if there is in the derived class
     --------------------------------------------------------------------------*/

    virtual  void convertToGridBar(int step, int kValue);
    
private:

    /** containing crit. pts at space which need to do special treatment*/
    DoubleArray critSpacePts; //critical pts at space dimension $unregistered
    IntArray isAddedSeg ;  // $unregistered

    /** registration */
    string  whichChgVarDim1String; // registered string chgVarType
    string  whichBarrierMethodMString; // registered string chgVarType
    string  whichSetSpaceStepMString; // registered string chgVarType
    string  solveMethodString; // registered string TFdSolveType

    /** convert the registered variables to ChangeOfVar*/
    void convertRegisteredString();

    /** diff. methods for set dx */
    void setVarSpaceSteps(double volForBound);
    void setVarSpaceStepsBasedonStd(double volForBound, TFdSetSpaceStepType whichSetSpaceStep, int n);
    
    /** -------------------------------------------------------------
        control which method to be used for setting var space step 
        NONE: const dx
        STD: only based on std
        STD_BAR: 1+ more pts around barriers, only allow 2 barriers at these moments
        3: more barriers such as in ladder... to be added 
    -------------------------------------------------------------*/

    TFdSetSpaceStepType whichSetSpaceStepM; // $unregistered

    /** --------------------------------------------------------------------------
        flag to distinguish normal FD and special treatment FD such as barrier, 
        It's a input 
    --------------------------------------------------------------------------*/

    bool needSpecialFD; 

    /** -------------------------------------------
        diff barrier methods
        const dx or variable dx
        FIX_GRID:  Fixed FD grid,
        VAR_GRID: FD Grid can be changed over time
    --------------------------------------------*/

    TFdBarrierMethodType whichBarrierMethodM; // $unregistered
    
    /** min time step for the added segment around critical dates such as barrier*/
    int minStepInAddedSegM; 

    /**--------------------------------------------
        false: special for no more than 2 barriers 
        true: general case for more than 2 barriers
    --------------------------------------------*/

    bool DEBUG_NODE_GEN; 

    bool DEBUG_forceNonNegative; // $unregistered

    string DEBUG_DumpToFile;

protected:
    TreeSliceEQ::RangeSP range;

    virtual TreeSliceSP createSlice( int dimBits ) const
    {
        return TreeSliceEQ::create( *range, dimBits ? dimBits : ( 1 << range->nDim ) - 1 );
    }

    /** state variable support for spot
        this basic spot support can be overriden in child engines */
    class TREE_DLL Spot : public IndexSpecEDR
    {
    public:
        Spot( const string & factorName ) :
            IndexSpecEDR( factorName )
        {}

        class TREE_DLL Product : public IndexSpecEDR::Product
        {
        public:
            Product( FDModel * model, const string & factorName ) :
                IndexSpecEDR::Product( model, factorName )
            {}

        protected:
            virtual void update( int & step, FDProduct::UpdateType type );

        public:
            void update( int & step, int & bot, int & top, FDProduct::UpdateType type );
        };

        virtual FDProductSP createProduct( FDModel * model ) const
        {
            return FDProductSP( new Product( model, factorName ) );
        }
    };
};

DRLIB_END_NAMESPACE
#endif
