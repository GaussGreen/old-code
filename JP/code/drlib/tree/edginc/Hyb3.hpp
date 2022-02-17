//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3.hpp
//
//   Description : hyb3 FD model class
//
//----------------------------------------------------------------------------

#ifndef HYB3_RATE_TREE_HPP
#define HYB3_RATE_TREE_HPP

#include "edginc/RateTree.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "cupsmodl.h"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL Hyb3 : public RateTree {
public:

    /********* FDModel interface ****************/

    /**  collect model initialisation data, set up timeline  */
    virtual void initModel(void);

    /**  finalise model initialisation, create payoffIndex  */
    virtual void finaliseModel(CControl* /*control*/) {}

    /** initialise tree */
    void initTreeData(void);

public:

    /********* RateTree interface ***************/

    /** insert into the tree a product configured zeroBank if tree is running in zerobank mode */
    virtual void insertNamedZeroBank(NamedZeroBank& zeroBank);

    /** register with the tree that the ParYield defined by rateSpec will be required */
    virtual void insertIRIndex(const IndexSpecIR& rate, DateTime date);

    /** populate the slice with the ParYield value associated with the rateSpec */
    virtual void getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, TreeSlice& slice);

    /** populate the slice with the current FX spot rate from the tree */
    virtual void getFXIndex(TreeSlice& slice);

    /** set tree mode to be either zerobank or claim bank */
    virtual void setZeroBankMode(ZeroBankMode::Enum mode) { zeroBankMode = mode; }

    /** retrieve zeroBank Mode for tree */
    virtual ZeroBankMode::Enum getZeroBankMode() const { return zeroBankMode; }

    /** insert zerobank dates/details into the tree */
    virtual void insertZero(const ZeroBond& zeroBond);
    
    void registerZero(
            DateTime obsDate, 
            DateTime matDate, 
            string curveName);

    void getZero(
            DateTime obsDate, 
            DateTime matDate, 
            string curveName, 
            TreeSliceSP &slice);

    /** get from tree the slice corresponding to zeroIndexSpec */
    virtual void getZero(const ZeroBond& zero, int step, TreeSlice& slice);

    virtual int getStep() const { return range->treeStep; }

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    virtual MarketDataFetcherSP createMDF() const;

protected:

    virtual void update(int t, FDProduct::UpdateType type);

    virtual void updateFinalize(int t, FDProduct::UpdateType type);

    virtual void updateStatePrice(int t, FDProduct::UpdateType type);

    virtual const TreeSlice & getStatePriceSlice(const int t) const ;

    virtual void sliceDev(TreeSlice& slice, int curveIdx) const;

    virtual void sliceEv(TreeSlice& slice) const;

    virtual void sliceExpand(TreeSlice& treeSlice, const string& curveName) const;

    virtual DateTime  getToday() const;

    virtual DateTime getDate(int step) const;

    virtual DateTimeArray getDates() const;

    virtual int getLastStep() const;

    DateTime getCurveValueDate(string curveName) const;

    virtual void print();   

protected:

    void populateIRParams(MKTVOL_DATA* mktvol, HYB3_TREE_DATA* tree, RateTree::IRModelParams* mp);

    virtual int getCurveIdx( const string & curveName ) const;
    int getCrvIdx(const string& curveName, int currIdx) const;
    int getCurrIdx(const string& curveName) const;

    HYB3_TREE_DATA   mTreeData;                ///< tree data $unregistered
    MKTVOL_DATA      mMktVolData[2];           ///< market volatility data $unregistered
    HYB3_DEV_DATA    mDevData;                 ///< dev data $unregistered

    DateTime    mToday;                        ///< today's date
    bool   cetSmoothing;

    T_CURVE     mTCurves[2][3];                ///< zero curves $unregistered

    FX_DATA     mFxData; // $unregistered
    EQ_DATA     mEqData; // $unregistered

    int         nbFactorsDom;
    int         nbFactorsFgn;

    // internal functions
    void init();
    void clear();


    bool                  mSmoothFlag;                  ///< tree node smoothing flag $unregistered
    ZeroBankMode::Enum    zeroBankMode;
    ZeroInterpStyle::Enum zeroInterpStyle;    /**< zero curve interpolation type */

    // StatePrice Information
    TreeSliceRatesSP    statePriceSlice;       /**< state Price at time point t */
    TreeSliceRatesSP    tempStatePriceSlice;   // temporary slice for state price calculation

    //--------------------------------------------------------------------------
    // ZEROBANK mode fields
    //--------------------------------------------------------------------------
    map<string, NamedZeroBank> mZeroBank;  // $unregistered

    typedef map<string, ZeroBondProdSP> ZeroProdMap;
    ZeroProdMap zeroProdList;

    //--------------------------------------
    // named zero bond start/maturity dates
    //--------------------------------------
    // generalised zero bond framework that allows for overlapping zeros, and 
    // news/frees slice internally when required
    map<string, ZeroBond> mZeroBond;         // $unregistered

    // ??? temp to stop deleting memory for unbuilt tree
    // until I improve the clear function
    bool mTreeBuilt; // $unregistered

    string mCurveName[2][3]; // $unregistered
    string mCurrency[2]; // $unregistered
    string mFXName;  // FX factor name $unregistered
    string mEQName;  // EQ factor name $unregistered

    bool   momentMatching;  // hy3b moment matching on/off $unregistered

    string treeDataDebugFile; // $unregistered
    // ZeroInterp mZeroInterpStyle;  // nicer here than anonymous global variable

    vector<IRDate> sortedCritDates;

public:

    /********************** reflection **************************/
    static void load(CClassSP& clazz);
    static CClassConstSP const TYPE;

    Hyb3(const CClassConstSP &type = TYPE);
    ~Hyb3();

};
typedef smartPtr<Hyb3> Hyb3SP;

/* For simplicity, there are 3 hyb3 models - one for each generalized 
   mode of the tree */


/*********************************** Hyb3 - FX Mode *******************************/

class TREE_DLL Hyb3FX : public Hyb3 {
public:

    /********************** reflection **************************/
    static IObject* defaultConstructor(void) { return new Hyb3FX(); }
    static void load(CClassSP& clazz);
    static CClassConstSP const TYPE;

    Hyb3FX() : Hyb3(TYPE), corrIR(0), corrForIRFX(0), corrDomIRFX(0), 
        fxIndexSpecSupplied(false), FXSmileDriftAdjustment(true) {}
    virtual ~Hyb3FX() {};

    virtual void initModel(void);
    virtual void getMarket(const MarketData*  market, IInstrumentCollectionSP instruments);
    virtual void retrieveFactor();
    virtual void initTreeData(void);
    virtual void recordOutput(Control* ctrl, Results* results) const;

    /****** implement IRVegaPointise::ISensitivePoints interface *******/
    /** returns all the points on the ir vol surface to which this
    model for the supplied instrument is sensitive.
    A null return value is ok - it
    will cause an exception to be thrown (rather than crash) */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst) const;

    // correlations for model
    double corrIR;
    double corrForIRFX;
    double corrDomIRFX;

    bool fxIndexSpecSupplied;

    /********************** exported fields *********************/
    int ppy;                   /**< minimum tree nodes/year */
    int maxStdDeviations;      /**< nb stdDevs to cut/trim tree */
    CcyIRParamsSP forIRParams; /**< Collection of IR parameters */
    CcyIRParamsSP domIRParams; /**< Collection of IR parameters */
    string FXSmileParams;
    string FXIRCorrelations;
    bool FXSmileDriftAdjustment;   /**< FX cups drift adjustment on or off */    
    CDoubleSP FXSpotVolOverride; 

private:
    void populateTreeIRParams(Currency curr);

    inline void getForeignDiffusionTCurve(T_CURVE &diffusionCurve) const;
    inline void getDomesticDiffusionTCurve(T_CURVE &diffusionCurve) const;
    inline IRVolRawSP getForeignIRVolRaw() const;
    inline IRVolRawSP getDomesticIRVolRaw() const;
};


DRLIB_END_NAMESPACE
#endif
