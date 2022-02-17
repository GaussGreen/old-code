//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FDModel.hpp
//
//   Description : FD/tree model base class.
//
//   Author      : Ning Shen
//
//   Date        : February 8, 2005
//
//----------------------------------------------------------------------------

#ifndef FDMODEL_HPP
#define FDMODEL_HPP

#include "edginc/Model.hpp"
#include "edginc/Format.hpp"
#include "edginc/FDProduct.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/TimeLine.hpp"
#include "edginc/LastSensDate.hpp"

DRLIB_BEGIN_NAMESPACE

// *****************************
// fd/tree model parent class
class TREE_DLL FDModel :public CModel,
    virtual public LastProductSensDate
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );

    /** fd solver interfcace */
    class TREE_DLL IFDSolver
    {
    public:
        virtual ~IFDSolver() {}

        virtual void roll() = 0;
    };
    typedef refCountPtr<IFDSolver> IFDSolverSP;

    /** interface for instrument to create product */
    class TREE_DLL IIntoProduct :
        virtual public CModel::IModelIntoProduct,
        virtual public IProdCreator
    {
    public:
        static CClassConstSP const TYPE;
        static void load(CClassSP& );

    private:
        friend class FDModel;

        virtual FDProductSP createProduct( FDModel * model ) const = 0;
    };

    /** interface for products modifiers (concrete types can be defined within products) */
    class TREE_DLL IProdModifier :
        virtual public IObject
    {
    public:
        static CClassConstSP const TYPE;
        virtual ~IProdModifier() {}
    private:
        static void load(CClassSP& );
    };
    DECLARE(IProdModifier);

    FDModel(const CClassConstSP &type);

    virtual ~FDModel() {}

    /** Implementation of LastProductSensDate - when to stop tweaking */
    virtual DateTime endDate(const CInstrument *instrument,
                             const Sensitivity *sensControl) const;

    //////////////  access methods ////////////
    /** get products */
    const FDProductArray & getProducts() const;
    /** get time metric */
    virtual TimeMetricConstSP getTimeMetric(){return timeMetric;}
    /** get last step (total num of steps) on time line */
    virtual int getLastStep() const { return timeLine->NumOfStep; }
    /** get date at step on the time line */
    virtual DateTime getDate(int step) const { return timeLine->StepDates[step]; }
    /** get date at step on the time line */
    virtual DateTimeArray getDates() const { return timeLine->StepDates; }
    //////////////  access methods ////////////

    /** creates new slice */
    virtual TreeSliceSP createSlice(
        const string & curveToDEV = "",
        const string & factorName1 = "",
        const string & factorName2 = "",
        const string & factorName3 = "") const = 0;

    /** override a control shift (eg for delta on tree/fd) */
    virtual SensControl* AlterControl(const SensControl* currSensControl) const;

    /** get the initial conditions and corresponding indexes in the slices
        for forward induction */
    virtual void getInitialConditions(
        IntArray& initialIndex, 
        DoubleArray& initialValue) const{ /*do nothing by default*/ }

     /** ask the model whether it is doing forward or backward induction */
    virtual bool IsFwdInduction() const;

    struct InductionType {
        enum Enum {  FWD, BACK, EXPRESSBACK };
    };

    virtual InductionType::Enum getInductionType() const ;

    /**  populate the list with dates for state prices, default implementation : do nothing */
    virtual void storeStatePriceDates(const DateTimeArray& statePriceDates ) {};

    /** asked the underlying models if they have a preference going 
        forward or backward and setting the FDModel such to satisfy those demands */
    virtual void decideDirection();

    enum LevelKind
    {
        LEVEL_GENERIC,
        LEVEL_BARRIER_UP,
        LEVEL_BARRIER_DOWN,
    };

    /** product should call this function in 'preCalc' to
        add critical 'level' of 'kind' to 'target' slice based on 'base' slice,
        optionally specifying 'value' of 'target' at the 'level' */
    virtual void addCriticalLevel(
        int step,
        const TreeSlice & base,
        double level,
        TreeSlice & target,
        LevelKind kind = LEVEL_GENERIC,
        double value = DBL_MAX )
    { /*do nothing by default*/ }

    /** This is meant to be a helper function for createSlice().
        It converts a simple slice to a 
        TreeSliceLayer of TreeSliceLayer of ... of the slice
        with respect to the state-variable supports registered
        with registerStateSupport */
    TreeSliceSP createLayer( const TreeSliceSP & slice ) const;

    /** add critical dates for timeline setup 
        default to one segment timeline */
    void addCritDates(const DateTimeArray& critDates);
    void addCritDate(const DateTime critDate);

    /** return a slive of statePrices (compute Arrow-Debreu prices) */
    virtual const TreeSlice & getStatePriceSlice(const int t) const
            { throw ModelException(getClass()->getName(), 
                            "getStatePriceSlice: not implemented yet"); };

    /** register the handling of a new sate variable */
    void registerStateSupport( TreeSliceLayer::StateSupport * owner );

    /** collects all the 'layers' of a given 'owner'
        if 'expandedOnly' then skip layers that are not expanded in state variable dimension  */
    void collectLayers(
        const TreeSliceLayer::StateSupport * owner,
        vector< TreeSliceLayer * > & layers,
        bool expandedOnly = true ) const
    {
        for( size_t k = 0; k < products.size(); ++k ) 
        {
            if( products[ k ]->isCalcOff() )
                continue;

            TreeSliceLayer::collectLayers(
                owner, products[ k ]->getSlicesToDEV(), layers );
        }
    }

    /** add data for timeline setup */
    void initSegments(
        const DateTimeArray& segDates, 
        const IntArray&      density);

    /** this is the main entry point for a model */
    virtual void Price(
        CInstrument* instrument, 
        CControl*    control, 
        CResults*    results);

    /** register with the model required zero/discounters */
    virtual void registerZero(
        DateTime useDate, 
        DateTime matDate, 
        string curveName) {
            throw ModelException(__FUNCTION__, 
                "Function not supported by model " + getClass()->getName());
    }

    /** Get slice representing zero/discounter between useDate and matDate 
     *  If "slice" is null it creates a new slice (no ownership kept), 
     *  otherwise it uses the slices passed to return the zero */
    virtual void getZero(
        DateTime useDate, 
        DateTime matDate, 
        string curveName, 
        TreeSliceSP &slice) {
            throw ModelException(__FUNCTION__, 
                "Function not supported by model " + getClass()->getName());
    }

    /** rates tree products distinguish between the today date and the value date -
        they may or may not be the same date.  The getToday returns the 
        reference date in the market that defines today.  If there is any possibility
        that valueDate != today, derived FDModel should override CModel::getValueDate
        to return the appropriate value date */
    virtual DateTime getToday() const {
        throw ModelException(__FUNCTION__, 
            "Function not supported by model " + getClass()->getName());
    }

    /** get discount/pricing yield curve and currency from the model */
    YieldCurveConstSP  getDiscountYC() const { return discYC; };

    ////////////////////////////////////////
    //*******  public interface methods **********

    /** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
    virtual double getPrice0( const TreeSlice & price ) const = 0;
    ////////////////////////////////////////

    /** accept or reject factors */
    virtual bool acceptFactor( IMarketFactor * factor ) { return true; }

    /** product creation request is always done through this function
        it takes topmost base class (IProdCreator not FDModel::IIntoProduct)
        to avoid down-casting before each call */
    FDProductSP createProduct( const IProdCreatorSP & creator );

    IProdModifierSP getProdModifier(CClassConstSP const &type);

    template <class C>
    smartPtr<C> getProdModifier() {
        return smartPtr<C>::dynamicCast(getProdModifier(C::TYPE));
    }

protected:
    /** creates actual products
        can be overridden by inherited models */
    virtual FDProductSP makeProduct(const IProdCreatorSP & creator);

    ////////////////////////////////////////
    //*******  model interface methods **********
    /** collect model initialization data, set up timeline */
    virtual void initModel() = 0;

    /** create a solver  */
    virtual IFDSolverSP createSolver() = 0;

    /** finalize model initialization */
    virtual void finaliseModel(CControl*    control) = 0;
    ////////////////////////////////////////

    /**  helper function, aggregate segDates, critDates */
    void arrangeDates();

    /** model to record any additional information requested */
    virtual void recordOutput(Control* ctrl, Results* results) const {};

    /**  rebuild grid or not  */
    virtual bool isRebuilt(const CControl* control);

    /** Invoked after instrument has got its market data. */
    virtual void getMarket(const MarketData* market, IInstrumentCollectionSP instruments);

    /** further selecting market data */
    virtual void retrieveFactor() {}

    /** chg of variable for underlying X, LOG_X, LOG_FWDX */
    typedef enum{ X, LOG_X, LOG_FWDX } ChangeOfVar;

    // to be reviewed
    int     stepsPerYearFD;
    bool    sameGridTweak;
    bool    isFwdInduction; // instead of isFwdInduction as this is currently used in some derived classes
    InductionType::Enum inductionType;  // normally backwards, but can be forward or express backward as well

    double sameGridDeltaShiftSize;

    // unregistered

    /** hold main product, can have products within prod */
    FDProductSP prod; // $unregistered

    YieldCurveConstSP  discYC;  // discount curve, could also be a factor
    IMarketFactorArray factors; // factors
    CDoubleMatrixSP    correls; // factors' correlations

    TimeMetricConstSP  timeMetric; // $unregistered
    TimeLineSP         timeLine; // $unregistered

    DateTimeArray   segDates; // timeline segment dates
    IntArray        density;  // relative points density for segments on timeline
    DateTimeArray   critDates; // critical dates that must be on timeline
    bool            sliceCreationMode;

private:
    /** this is the main sequence of function calls called by Price() */
    void price(const IProdCreatorSP & creator,
               CControl* control, CResults* results, bool rebuild);

    /** call init() on all products */
    void initCall(CControl* control);

    /** call initProd() on all on all products */
    void initProdCall();

    /** print to stream any debug/other information each FDProduct chooses to do */
    void printProductInfo(ostream& outputStream);

    /** call recordOutput() on all products, decorating results with instrument/component/indexSpec name:: */
    void recordOutputCall(Control* control, YieldCurveConstSP disc, Results* results);

    /** Helper method to answer date of last step for endDate processing. */
    DateTime lastStepDate(IIntoProduct &creator, const Sensitivity *sensControl);

    /** keep all the products in the right order */
    FDProductArray products; // $unregistered

    /** cache of all prodCreators used to create products */
    vector< IProdCreatorSP > creators; // $unregistered
    /** indecies of products corresponding to creators */
    vector< int > prodIdxs; // $unregistered

    /** for now, filename of file stream to printProductInfo to */
    string productInfoFileName;

    /** list of component modifiers */
    IProdModifierArray prodModifiers;

    /** list of all unique StateSupports */
    mutable vector< TreeSliceLayer::StateSupport * > stateSupportCache;
};

typedef smartPtr<FDModel> FDModelSP;


DRLIB_END_NAMESPACE
#endif
