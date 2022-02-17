//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2D.hpp
//
//   Description : two factor finite difference engine
//
//   Author      : Ning Shen
//               : Xiaolan Zhang
//
//   Date        : November 29, 2004
//
//----------------------------------------------------------------------------

#ifndef FD2D_HPP
#define FD2D_HPP

#include "edginc/LatticeModelEDR.hpp"
#include "edginc/IndexSpecEDR.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD2D : public LatticeModelEDR {
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );

    friend class FD2DSpot;
    friend class FD2DSpot2;
    friend class FD2DSolver;
    friend class FD2DSolverJumps;

    virtual ~FD2D();
    FD2D(CClassConstSP type);

    //--------------------------  public methods   ------------------------------------------

    /** validate some inputs */
    virtual void validatePop2Object();
    
    /** calc the underlying value */
    virtual void computeUndLevel1(int bot, int top, int step, double*s) const;
    virtual void computeUndLevel2(int bot, int top, int step, double*s) const;
    
protected:
    /** Invoked after instrument has got its market data. */
    virtual void getMarket(const MarketData* market, IInstrumentCollectionSP instruments);

    /** accept or reject factors */
    virtual bool acceptFactor( IMarketFactor * factor );

    /**-----------------------------------------------------*/
    /** model initialisation (first part)                   */
    /** get vol                                             */
    /** call timeline builder using initData                */
    /** set the dates and FD params at time direction       */
    /** called before initProd                              */
    /**-----------------------------------------------------*/

    virtual void initModel();

    /** model parameter validation*/

    /**-----------------------------------------------------*/
    /** model initialisation  (second part)                 */
    /** most memory init should be here                     */
    /** finalise model initialisation once after initModel()*/
    /** and product init()                                  */
    /** called by validate                                  */
    /**-----------------------------------------------------*/

    virtual void finaliseModel(CControl*    control);

    /**-----------------------------------------------------*/    
    /** calculate the stuff which specifies the grid        */
    /** set up the FD boundaries, and calc space step,      */
    /**  can be const or var.                               */
    /** shouldn't be called if SameGrid = true              */
    /**-----------------------------------------------------*/
    
    virtual void SetSpaceGridData();  //to be reviewed
    
    
    /**-----------------------------------------------------*/
    /** each derived model can set diff. boundaries         */
    /** based on the dynamics of underlying,                */
    /** alpha is the input truncation1D (nb of std)         */
    /** outLowB, outUpB are low and up boundaries for FD    */
    /**-----------------------------------------------------*/
    
    virtual  void setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1,
                              double& volForBound2, double alpha2, double& outLowB2, double& outUpB2){};
    
    
    /**-----------------------------*/
    /** fd two factors level set up */
    /** ----------------------------*/
    
    void computeFdGridLevel (int bot, int top, vector<double>& vdx, 
                             double lBoundary, double uBoundary, vector<double>& s);


    /**-----------------------------*/
    /** update payoffIndex          */
    /** ----------------------------*/
   
    virtual void payoffIndexUpdate (int& step, FDProduct::UpdateType);

    /**----------------------------------------------------------*/
    /** calculate pde operator coefficients                      */
    /** extract coeff computed from the market data for the PIDE.*/
    /** be aware that the structure and coefficients for forward */
    /** and backward equation can be very different.             */
    /** forward  equation not available yet                      */
    /**---------------------------------------------------------*/
    
    /**----------------------------------------------------------------*/
    /**  setup the coefficients of the PIDE needed to be solved        */
    /**  Assume the most general PDE is                                */
    /**  U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;*/
    /**  this ft is called at each time step, so they're time depd.    */
    /**----------------------------------------------------------------*/
    virtual  void pdeCoeff(int step, double*** coeff, 
                           int bot1, int top1, int bot2, int top2 ) = 0;

    /** get segment index of a fd*/
    //int getSegSegment(int step);
    
    /** override FDModel::makeProduct to retrieve spot index product */
    virtual FDProductSP makeProduct( const IProdCreatorSP & creator );

    /** create a solver */
    virtual IFDSolverSP createSolver();

    // ------------------------------ local data ----------------------------------

    CAssetConstSP  underlying; // $unregistered

    FDProductSP    payoffIndex;     // $unregistered

    /**general PDE 2F has 6 coeffs, jumps model will be more.*/
    int numCoeffPDE;    // $unregistered
    
    /** fwd stock at each step */    
    vector<DoubleArray>   stepForward; //should replace stepForward1 and stepForward2 $unregistered
    
    vector<double> gridLevel1; // 2-d grid levels $unregistered
    vector<double> gridLevel2; // 2-d grid levels $unregistered
    
    /**-----------------------------------------------------------------------------------*/
    /** the following two var tell us FD will be done between which index                 */
    /** the last element in the topDim1 contains the biggest index amonge the numOfValues */
    /** the last element in the botDim1 contains the smallest index amonge the numOfValues*/
    /**------------------------------------------------------------------------------------*/
    
    // not good : to be reviewed !!!
    int** topDim1; // top grid index for Und1 $unregistered
    int** botDim1; // bottom grid index for Und1 $unregistered
    int** topDim2; // top grid index for Und2 $unregistered
    int** botDim2; // bottom grid index for Und2 $unregistered
    

    double  lBound1; //Low boundary for factor1 $unregistered
    double  uBound1; //Upper boundary for factor1 $unregistered
    double  lBound2; //Low boundary for factor2 $unregistered
    double  uBound2; //Upper boundary for factor2 $unregistered

    /** -----------------------------------------------------------------------------*/
    /** v_dxM is the working vector, may change based on barriers level at step      */
    /** v_dxMOrg is the reference vector, don't change with barrier level at step    */
    /** initially, v_dxM = v_dxMOrg                                                  */
    /** M refer to IModel                                                             */
    /** -----------------------------------------------------------------------------*/

    /** -----------------------------------------------------------------------------*/
    /** space steps                                                                  */
    /** if const grid, FD scheme uses v_dxM[1], v_dyM[1]                             */
    /** -----------------------------------------------------------------------------*/
    
    vector<double> v_dxM;  //M refer to Model $unregistered
    vector<double> v_dyM;  //M refer to Model $unregistered

    /** ----------------------------------------------------------------------------*/
    /** model save it, otherwise, has some trouble at delete                        */
    /** payoff maybe deleted before model destructor                                */
    /** -----------------------------------------------------------------------------*/
    
    int nProd;  //nb of Products $unregistered
    int maxNumOfValue; //max nb of values for all Products $unregistered
    
    /**--------*/
    /** flags  */
    /**--------*/
    /** -----------------------------------------------------------------------------*/
    /** hasJumps : false, solver will not call jumpsolver such as FFT...             */
    /**------------------------------------------------------------------------------*/
    bool hasJumps; // $unregistered

    bool needRecalcUndValue1;  //true, undValue depends on t. If chg var. = log(fwd), this should be set to true $unregistered
    bool needRecalcUndValue2; // $unregistered

    //temporary
    //for test
    bool isVariableGrid;         
    
    /**------------------------------------------------------------------*/
    /** User Input                                                       */
    /**------------------------------------------------------------------*/
    
    /** Different type of FD scheme ADI or ADVANCE_ADI */
    typedef  enum{ADI, ADVANCE_ADI} TFdSolveType;
    
    /** convert the registered variables to ChangeOfVar*/
    void convertRegisteredString();
    
    ChangeOfVar    whichChgVarDim1; // not registered string chgVarType $unregistered
    ChangeOfVar    whichChgVarDim2; // $unregistered
    
    /** ------------------------*/
    /** It's a param for solver, put it here only for temporary input purpose*/
    /** ADI: traditional ADI with no cross term*/
    /** ADVANCE_ADI : advanced ADI with cross term*/
    /** ------------------------*/
    TFdSolveType solveMethod; // not registered string TFdSolveType $unregistered

    /**-------------------------------------------------------------------*/
    /** StepsPerYear: num of steps per year, actual steps                 */
    /**               created in TimePts may be different                 */
    /** Truncation1D : num of stdev for truncation for each dimension     */
    /** Dim1 : original dimension of Und1 size, need to be odd nb         */
    /** Dim2 : original dimension of Und2 size                            */
    /**-------------------------------------------------------------------*/
    double truncation1D;
    double truncation2D;
    int    dim1;
    int    dim2;
    int    stepsPerYear;
    
    // ------------------------------ For Barriers ----------------------------------
    
    /**Barriers*/
    bool needSpecialFD; // $unregistered
    
private:
    string  whichChgVarDim1String; // registered string chgVarType
    string  whichChgVarDim2String; // registered string chgVarType
    string  solveMethodString; // registered string TFdSolveType

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
        };

        virtual FDProductSP createProduct( FDModel * model ) const
        {
            return FDProductSP( new Product( model, factorName ) );
        }
    };
    class TREE_DLL Spot2 : public IndexSpecEDR
    {
    public:
        Spot2( const string & factorName ) :
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
        };

        virtual FDProductSP createProduct( FDModel * model ) const
        {
            return FDProductSP( new Product( model, factorName ) );
        }
    };
};

typedef smartPtr<FD2D> FD2DSP;

DRLIB_END_NAMESPACE
#endif
