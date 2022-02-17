//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : RateTree.hpp
//
//   Description : top level model class for rates trees
//
//----------------------------------------------------------------------------

#ifndef RATE_TREE_HPP
#define RATE_TREE_HPP

#include "edginc/FDModel.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/IRModelConfig.hpp"
#include "edginc/IndexSpec.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/IRExoticParam.hpp"
#include "edginc/MarketTable.hpp"
#include "edginc/MarketObjectQualifiers.hpp"
#include "esl.h"

DRLIB_BEGIN_NAMESPACE

/******************************* RateTree ******************************/
// forward declarations
class ZeroBond;
class IndexSpecEQ;
class IndexSpecFX;
class IndexSpecIR;

struct ZeroInterpStyle {
    enum Enum { LINEAR, FLAT_FWD };
};
typedef BoxedEnum<ZeroInterpStyle::Enum> ZeroInterpStyleBoxedEnum;

struct ZeroBankMode {
    enum Enum { ZEROBANK, CLAIMBANK };
};
typedef BoxedEnum<ZeroBankMode::Enum> ZeroBankModeBoxedEnum;


class TREE_DLL RateTree : public FDModel,
                 virtual public FDModel::IFDSolver,  // must implement solver
                 // Subclasses must provide grid point sensitivity/exposure
                 // info for smart vega tweaking and exposure reporting.
                 virtual public IRVegaPointwise::ISensitivePoints
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    //enum ZeroBankMode {CLAIMBANK,  ZEROBANK};
    enum Currency {FOREIGN=0, DOMESTIC=1};

    // simple structure used as a container for zero bank relevant fields
    struct TREE_DLL NamedZeroBank
    {
        NamedZeroBank(const string &name="")
            : name(name), nbZeros(0), nbCurrentZeros(0), zeroBankDim(0) {}

        // fields populated by client
        string name;       // associated indexSpecIR name
        string curveName;  // name of curve index is reseting/fixing on
        int    nbZeros;    // number of zeros required at any one time
        DateTimeArray zeroMatDateList;
        string currency;

        // ??? internal parameters used by tree - move these to structure
        // inside tree later on - just keep here for now
        int nbCurrentZeros;
        int zeroBankDim;
        int devMode;
        int discCurve;
        vector<IRDate> zeroCritDates;
        vector<double*>  zeroBank;
        Currency currencyEnum;
        virtual ~NamedZeroBank();
    };

    // simple structure used as a container for the optimized/non optimized
    // claim bank
    struct NamedClaimBank {

        NamedClaimBank() : dimension(-1), foreignDiscountedAsDomestic(false),
            critDatesIdx(-1) {
            Esl_CbkInit(&zeroBank);
        };

        // used only for hyb based tree models (not used in single currency)
        int dimension;  // slice dimension (currency dependent)
        // used for the domestic equivalent discounting of foreign denominated
        // discount factor/zero calculated/managed within the hybrids tree.
        bool foreignDiscountedAsDomestic;

        string curveName;
        int critDatesIdx;
        CLAIM_BANK  zeroBank;
        vector<IRDate> critZeroMatDates;
        vector<IRDate> critZeroUseDates;
        vector<IRDate> optZeroMatDates;
        vector<IRDate> optZeroUseDates;
        // combined critical and optional date list
        vector<IRDate> combZeroMatDates;
        vector<IRDate> combZeroUseDates;
        // unsorted/pre duplicate removal data used to print out
        // additional debug information
        vector<IRDate> critZeroMatDatesUnsorted;
        vector<IRDate> critZeroUseDatesUnsorted;
        vector<string> critZeroLabel;  // to improve debug info

        void saveUnsorted(); // save for debug info print
    };

    // container for interest rate model parameters
    struct TREE_DLL CcyIRParams : public CObject {

        int nbFactors;
        string volCalibIndex;
        IRCalibWrapper irCalib;
        string smileSet;
        string modelSet;
        YieldCurveWrapper curveToDiffuse;
        YieldCurveWrapper curveToDiscount;
        IRExoticParamTableWrapper smileTable;     // smile map collection
        IRExoticParamTableWrapper modelTable;     // model map collection
        IRVol::VolType          calibVolType;
        IRVol::CalibType        calibType;
        ExpirySP                calibTenor;
        CDoubleSP               calibVolOverride;


        /************* methods *************/
        CcyIRParams(const CClassConstSP &type=TYPE) : CObject(type),nbFactors(1)
            ,calibVolType(IRVol::BASE_VOL),calibType(IRVol::CAP){}
        void getData(CModel *model, const MarketData*  market);
        void validate(string const&ccy, string const& discName) const;

        /************* type *************/
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);
        static IObject* defaultConstructor() { return new CcyIRParams(); }
    };
    typedef smartPtr<CcyIRParams> CcyIRParamsSP;
    typedef array<CcyIRParamsSP, CcyIRParams> CcyIRParamsArray;

    // simple structure to store the IR model parameters
    struct TREE_DLL IRModelParams
    {
        int     nbFactors;       // number of factors (1..3)
        int     PPY;             // Points per year
        double  FactorMR[3];     // Mean Reversion
        double  FactorVol[3];    // Factor Vol (Weighting)
        double  FactorCorr[3];   // Correlation [1&2] [1&3] [2&3]

        /* smile parameters */
        double  QLeft;
        double  QRight;
        double  FwdShift;
        int     CetNbIter;

        /* others */
        int    nbStdDevs;
        double backBone;
        int    nbStateVariables;
        int    stateVarStdDevs;
        IRModelParams();
    };

    // Recurse the supplied instruments to retrieve the first YieldCurve object
    // which name is "name"
    static YieldCurveConstSP collectYC(string const& ycName, IObjectConstSP obj);

    // recurse the instrument components to find and store all the IMarketFactor types found
    // in the exported fields. Ex: these are yield curve types for fix3
    static void collectFactors(FDModel * model,                // IN
                               IObjectSP instruments,          // IN
                               IMarketFactorArray & factors);  // OUT

    /** implement creates new slice */
    virtual TreeSliceSP createSlice(
        const string & curveToDEV = "",
        const string & factorName1 = "",
        const string & factorName2 = "",
        const string & factorName3 = "" ) const;

    virtual TreeSliceSP createSimpleSlice(
        const string & curveToDEV) const;

    /****** implement solver interface *******/
    /** rate tree specific roll function (contains the time loop) */
    virtual void roll(void);

    /** Forces the tree to be rebuilt for each greek computation */
    virtual bool isRebuilt(const CControl* control) {return true;}

    /** create the rate tree specific solver */
    virtual IFDSolverSP createSolver() { 
        NullDeleter deleter;
        IFDSolverSP p(this, deleter);
        return p; 
    }

    /****** declare Rate Tree interface methods ******/

    /** insert into the tree a product configured zeroBank if tree is running in zerobank mode */
    virtual void insertNamedZeroBank(NamedZeroBank& zeroBank) { 
        throw ModelException("insertNamedZeroBank", "method not supported by rate model " + TYPE->getName());
    }

    /** register with the tree that the ParYield defined by rateSpec will be required */
    virtual void insertIRIndex(const IndexSpecIR& rateSpec, DateTime date) {};
    virtual void insertIRIndex(const IndexSpecIR& rateSpec, DateTime resetDate, DateTime lateResetDate, bool isCrit);

    /** populate the slice with the ParYield value associated with the rateSpec */
    // ??? depricate
    virtual void getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, TreeSlice& slice) {
        throw ModelException("getIRindex(rateSpec, currentDate, treeSlice) not supported by model");
    }
    /** populate the slice with ParYield value with resetDate/currentDate adjustments */
    virtual void getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, DateTime resetDate, 
                            TreeSlice& treeSlice) {
        throw ModelException("getIRindex(rateSpec, currentDate, resetDate, treeSlice) not supported by model");
    }

    /** populate the slice with the current FX spot rate from the tree */
    // ??? depricate
    virtual void getFXIndex(TreeSlice& slice) {
        throw ModelException(__FUNCTION__, "model " + getClass()->getName()+" does not support FX");
    }

    /** populate slice with FX spot rate with resetDate/currentDate adjustments */
    virtual void getFXIndex(const IndexSpecFX& fxSpec, DateTime currentDate, DateTime resetDate,
                            TreeSlice& treeSlice) {
        throw ModelException(__FUNCTION__, "model " + getClass()->getName()+" does not support FX");
    }

    /** populate the slice with the current EQ spot rate from the tree */
    virtual void getEQIndex(TreeSlice& slice) {
        throw ModelException(__FUNCTION__, getClass()->getName()+" does not support Equity");
    }

    virtual void getEQIndex(const IndexSpecEQ& eqSpec, DateTime currentDate, DateTime resetDate,
                            TreeSlice& treeSlice) {
        throw ModelException(__FUNCTION__, "model " + getClass()->getName()+" does not support EQ");
    }

    /** default assumes trees running in claimBank mode */
    virtual ZeroBankMode::Enum getZeroBankMode() const { return ZeroBankMode::CLAIMBANK; }

    /** insert zerobank dates/details into the tree */
    virtual void insertZero(const ZeroBond& zeroBond) {
        throw ModelException("RateTree::insertZero", getClass()->getName() + " model does not support this function.");
    }

    /** get from tree the slice corresponding to zeroIndexSpec */
    virtual void getZero(const ZeroBond& zeroBond, int step, TreeSlice& slice) {
        throw ModelException("RateTree::getZero", getClass()->getName() + " model does not support this function.");
    }

    /** we re-use same slice for each time step in all trees, so this index is always 0 */
    virtual int getSliceIndex( int step ) const  { return 0;}

    /** returns an ExposureHighlighter - a "model" that does everything except
        actually price, so you get to see what market data it uses
        Default implementation supplied */
    virtual ExposureHighlighter* exposureHighlighter();

    /** utility function to flag event as true if stepDate is found in the array of dates */
    static bool happensNow(const DateTimeArray& array, DateTime stepDate, int *arrayPos=0);

    /** adjust a stochastic indexSpec index slice value to account for the differences in the
        date at which is was calculated by the tree and the actual reset date (essentially
        first moment matching) - usually used to simplify pricing algorithms for stubs and
        forward underlyings */
    virtual void dateAdjustIRIndex(const IndexSpecIR& indexSpec, DateTime currentDate, 
        DateTime resetDate, TreeSlice& treeSlice) {
        throw ModelException("RateTree::dateAdjustIRIndex", "Function not suppoted by model");
    }

    /**  populate the list with dates for state prices, default implementation : do nothing */
    virtual void storeStatePriceDates(const DateTimeArray& inpDates ) { statePriceDates = inpDates; };


    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant --- not sure if this is always right in
     * the long run.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** override FDModel::makeProduct to retrieve IndexSpec product */
    virtual FDProductSP makeProduct(const IProdCreatorSP & creator);

    // Class to select volatility info from QLib IRVol objects and populate
    // Rates MKTVOL_DATA structure.  Also provides vol grid point exposure
    // recording.
    class IRVolSelector {
    public:
        IRVolSelector(
            IRVolRawSP     volSP,
            const T_CURVE& diffusionCurve,
            const string&  volCalibIndex,
            bool           smoothing);

        IRVolSelector(
            IRVolRawSP     volSP,
            const T_CURVE& diffusionCurve,
            const string&  volCalibIndex,
            bool           smoothing,
            CcyIRParamsSP  ccyIRParams);

        // Select volatility info from QLib IRVol objects and populate
        // Rates MKTVOL_DATA structure.
        void getMktVolData(MKTVOL_DATA &mktVolData)
        {
            selectVols(mktVolData, IRGridPointAbsArraySP(), OutputNameConstSP());
        }

        // Populate vol exposures with a list of grid points for the
        // specified IRVol.
        void getVolExposures(IRGridPointAbsArraySP volExposuresSP, OutputNameConstSP irVolName);

    private:
        void selectVols(MKTVOL_DATA &mktVolData,
                        IRGridPointAbsArraySP volExposuresSP,
                        OutputNameConstSP irVolNameSP);

        void selectVolsMQ(MKTVOL_DATA &mktVolData,
                          IRGridPointAbsArraySP volExposuresSP,
                          OutputNameConstSP irVolNameSP);

        void validate();

        void addBaseVolExposures(IRGridPointAbsArraySP volExposuresSP);
        void addSwapVolExposures(IRGridPointAbsArraySP volExposuresSP);

        IRVolRawSP     mVolSP;
        const T_CURVE& mDiffusionCurve;
        string         mVolCalibIndex;

        IRVol::VolType   mCalibVolType;
        IRVol::CalibType mCalibType;
        ExpirySP         mCalibTenor;
        CDoubleSP        mCalibVolOverride;

        BASEVOL_EXPOSURE_DATA mSelectedBV;
        SWAPVOL_EXPOSURE_DATA mSelectedSV;
        ExpiryArray           mBaseVolTenors;
        ExpiryArray           mBaseVolExpiries;
        ExpiryArray           mSwapVolTenors;
        ExpiryArray           mSwapVolExpiries;
        DateTime              mBaseVolBaseDate;
        bool                  mSmoothing;
    };

    static IRVolRawSP getIRVolRaw(RateTree::CcyIRParamsSP modelIRParams);

protected:

    RateTree(const CClassConstSP &type = TYPE);
    virtual ~RateTree() {}

    /****** implement solver interface *******/
    /** rate tree specific roll function (contains the time loop) */
    virtual void rollFront(void); // implementation for forward induction
    virtual void rollBack(void); // implementation for backward induction
    virtual void produceStatePrices(void); // rolls forward and stores state prices

    /** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
    virtual double getPrice0( const TreeSlice & price ) const;

    /** Update tree variables at the "top" of the time loop */
    virtual void update(int t, FDProduct::UpdateType type) = 0;

    /** Implement any updating logic that must occur at the "bottom" of the time loop after
        all product updates */
    virtual void updateFinalize(int t, FDProduct::UpdateType type) = 0;

    /** Update state prices in the tree (compute Arrow-Debreu prices) */
    virtual void updateStatePrice(int t, FDProduct::UpdateType type) 
            { throw ModelException(getClass()->getName(),"updateStatePrice: not implemented yet"); };

    /** Perform one backward induction step from t+1 to t by calculating expected
        value and discounting on the provided IR curve */
    virtual void sliceDev(TreeSlice& slice, int curveIdx) const = 0;

    /** Perform one backward induction step from t+1 to t by calculating expected
        value only. */
    virtual void sliceEv(TreeSlice& slice) const = 0;

    /** Expand the dimension of the slice so that it can be DEV'd. */
    // ??? remove
    virtual void sliceExpand(TreeSlice& slice, const string& curveName) const {
        throw ModelException("RateTree::sliceExpand", getClass()->getName() +
                             " model does not support function");
    }

    /** return the value date of the tree, or for the particular IR curve name */
    virtual DateTime getCurveValueDate(string curveName) const = 0;

    // ??? CModel::getValueDate  here for legacy tree modes - remove and make each
    // model override this function when removing legacy trees
    virtual DateTime getValueDate() const { return getCurveValueDate(""); }

    /** print debug information for the tree - nodes, discount factors, vols etc. */
    virtual void print() = 0;

    /** return curev index by name */
    virtual int getCurveIdx( const string & curveName ) const = 0;

    /** Extend the zero curve to the specified new start date */
    static void ExtendTreeCurve(T_CURVE *tc, const DateTime &newStartDate);

    RatesSliceRangeSP range; // $unregistered
    int         mTpIdx;                   ///< current time point = currStep in FDModel

    IndexSpecArray     indexSpecs;  // array of index specs if required by model
    IMarketFactorArray fxFactors;   // array of FXAsset type market factors
    IMarketFactorArray eqFactors;   // array of EQAsset type market factors
    string             engineSet;   // engine key name
    IRModelConfigTableWrapper engineTable;

    // temporary claim bank models legacy pricing mode flag - false by default and should only
    // be exported by CB models as it will only effect them.  It has to be in rate tree as the
    // insertIRIndex claim bank mode is defined at this level.
    bool CBLegacyPricingMode;

    // max number of fwd starting swap days to permit determinstic ratio adjustment
    static const int maxFwdPeriodIRAdj = 10;  
    // max number of days between reset and observation date to permit ratio adjustmnet
    static const int maxIRDateOffsetAdj = 180;   

protected:
    map<string, NamedClaimBank> claimBank;

    // information of state prices
    DateTimeArray  statePriceDates; 
    map<int, TreeSliceSP > statePrices;

private:
    // helper to loop slice layers
    void loopDev(const vector< TreeSliceSP >& slices, int curveIdx = -1);
    void loopExpand(TreeSlice &slice, const string &curveToDev);

};
typedef smartPtr<RateTree> RateTreeSP;


// ************ ZeroBond, is temp here ************
    // Collection of associated zero bond start/end dates
    // ??? for the meantime in hyb3 assume this is always 2D
class TREE_DLL ZeroBond : public CObject, virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;

    ZeroBond(const DateTimeArray& startDates=DateTimeArray(),
             const DateTimeArray& matDates=DateTimeArray(),
             const string& curveName="");
    ~ZeroBond();

    /** returns critical dates */
    virtual DateTimeArraySP getCritDates(void) const;

    const string& getName(void) const { return name; }

    DateTimeArray startDates; // $unregistered
    DateTimeArray matDates; // $unregistered
    string curveName;  // each discount bond may be different curve $unregistered
    string name; // unique name $unregistered

    // ??? internal parameters used by tree - move these to structure
    // inside tree later on - just keep here for now
    // ??? internal slices used by the tree should be tree slice, but just use
    // double* for now as there are various template/constructor complaints
    string currency; //!!! $unregistered
    double *zeroSlice;   // ??? should be a stack to support overlapping zero dates $unregistered
    RateTree::Currency currencyEnum; // $unregistered
    int sliceDim; // $unregistered
    int devMode; // $unregistered

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new ZeroBond(); }

    virtual FDProductSP createProduct( FDModel * model ) const;
};
typedef smartPtr<ZeroBond> ZeroBondSP;
typedef smartConstPtr<ZeroBond> ZeroBondConstSP;

/**********************************
 *
 *  Elementary Product for RateTree
 *  Zero/Discounter
 *
 **********************************
?? zero bond/discounter elementary product -
currently only takes list of dates and assumes zeros don't
overlap - generalise to array of 1 per coupon payment date to
allow for delayed payments
 */
class TREE_DLL ZeroBondProd : public FDProduct {
public:
    ZeroBondProd(const ZeroBondConstSP &inst, FDModel* model);

    virtual void init(Control*) const; // model init
    virtual void initProd(void);       // product init
    const TreeSlice & getValue(int step, DateTime matDate) const;

    virtual const TreeSlice & getValue(int step) const { return getValue(step, DateTime()); }
    virtual void     update(int & step, FDProduct::UpdateType) {}
    virtual DateTime getStartDate(void) const { return model->getDate(0);}
    virtual bool     isElementary(void) const {return true;}
    virtual void     recordOutput(Control*, YieldCurveConstSP, Results*) {}

private:
    TreeSliceSP mainSlice;
    ZeroBondConstSP inst;
    RateTree *rateTree;
};
typedef refCountPtr<ZeroBondProd> ZeroBondProdSP;

/****************************** RateTreeMDF ********************************/
class TREE_DLL RateTreeMDF : public MarketDataFetcher {
public:
    RateTreeMDF(int num, const MarketObjectQualifiersSP quaifiers = MarketObjectQualifiersSP());

    /** Resolve from multiple market data of same name and type */
    MarketObjectConstSP resolveMultiData(const MarketObjectArray& marketObjs,
                                         CClassConstSP            type) const;
private:
    RateTreeMDF(const RateTreeMDF& rhs);
    RateTreeMDF& operator=(const RateTreeMDF& rhs);

    MarketObjectQualifiersSP qualifiers; 
};

DRLIB_END_NAMESPACE

#endif
