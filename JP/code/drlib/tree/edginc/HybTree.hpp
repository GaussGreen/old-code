//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : HybTree.hpp
//
//   Description : HybTree
//
//----------------------------------------------------------------------------

#ifndef HYBTREE_HPP
#define HYBTREE_HPP

#include "edginc/RateTree.hpp"
#include "edginc/DateTime.hpp"

#include "edginc/IRCalib.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/DividendList.hpp"
#include "cupslib.h"
#include "hyb4_lib.h"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL Hyb4temp : public RateTree {
public:

    /********* CModel interface *****************/
    /** value date */
    virtual DateTime getValueDate() const { return valueDate; }

    /********* FDModel interface ****************/

    /** finalise model initialisation - nothing to do for hyb3 */
    virtual void finaliseModel(CControl* control) {};  // ??? make blank function in rateTree

    /** register with the model required zero/discounters */
    virtual void registerZero(DateTime useDate, DateTime matDate, 
                              string curveName);

    /** request from model a zero/discounter slice */
    virtual void getZero(DateTime useDate, DateTime matDate, string curveName, 
                         TreeSliceSP &slice);

    /** today date defined as market reference/today date, including time (supports SOD or EOD)*/
    virtual DateTime getToday() const { return today; }

public:

    /********* RateTree interface ***************/

    /** populate the slice with the ParYield value associated with the rateSpec */
    virtual void getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, 
                            DateTime resetDate, TreeSlice& treeSlice);

    /** set tree mode to be either zerobank or claim bank */
    virtual ZeroBankMode::Enum getZeroBankMode() const { return ZeroBankMode::CLAIMBANK; }

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    virtual MarketDataFetcherSP createMDF() const;

    /****** implement IRVegaPointise::ISensitivePoints interface *******/
    /** returns all the points on the ir vol surface to which this
    model for the supplied instrument is sensitive.
    A null return value is ok - it
    will cause an exception to be thrown (rather than crash) */
    //virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
    //    OutputNameConstSP  outputName,
    //    const CInstrument* inst) const;

    /** adjust the value of a stochastic IR index to account for currentDate != resetDate */
    void dateAdjustIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, 
                                   DateTime resetDate, TreeSlice& treeSlice);

    virtual TreeSliceSP createSimpleSlice( const string & curveToDEV ) const;

protected:

    virtual void update(int t, FDProduct::UpdateType type);

    virtual void updateFinalize(int t, FDProduct::UpdateType type);

    virtual void sliceExpand(TreeSlice& treeSlice, const string& curveName) const;

    virtual void updateStatePrice(int t, FDProduct::UpdateType type);

    virtual const TreeSlice & getStatePriceSlice(const int t) const;

    // virtual void sliceDev(TreeSlice& slice, const string& curveName) const;

    // virtual void sliceEv(TreeSlice& slice) const;

    virtual DateTime getDate(int step) const;

    virtual DateTimeArray getDates() const;

    virtual int getLastStep() const;

    virtual DateTime getCurveValueDate(string curveName) const;

protected:  
    
    // these functions are accessed by all the derived hyb3 models

    void clear();

    size_t getSliceSize(size_t nbDim) const;

    inline void getDiffusionTCurve(T_CURVE& diffusionCurve) const;

    void configureClaimBank(const CRIT_DATE* critDate, int nbCritDate);

    DateTime getCurrentDate() const {
        return DateTime::fromIrDate(treeData3.TPDate[mTpIdx]);
    }

    void populateIRParams(MKTVOL_DATA* mktvol, HYB3_TREE_DATA* tree, 
                          RateTree::IRModelParams* mp);

    /** RateTree interface */
    virtual int getCurveIdx( const string & curveName ) const;

    /** return internal T_CURVE index number from the available list */
    virtual int getCrvIdx(const string& curveName, int currIdx) const;

    /** return internal array index number for the currency of the curve */
    virtual int getCurrIdx(const string& curveName) const = 0;

    /** return the internal slice dimension of the curve Name - tree mode dependent */
    virtual int getCrvDim(const string& curveName) const = 0;

    /** initialise all tree structures */
    void init();


    /***************************** variables ********************************/

    enum {NBCRV = 3};             ///< number of curves

    DateTime today;
    DateTime currencyValueDate[2];  // foreign and domestic may be different
    DateTime valueDate;  // defined as EOD today date
    string curveName[2][NBCRV];
    string currency[2];
    string mFXName;  // FX factor name
    string mEQName;  // FX factor name

    vector<IRDate> sortedCritDates;  // only used in debug print function

    // ???? AK: should be part of the underlying curve object
    double today2ValDateZeros[2][NBCRV];   ///< today to value date discount factors

    bool treeBuilt;  // ??? remove when decent clear function in Hyb3

    // all the internal tree structures for the various hyb3 derived modes
    HYB3_TREE_DATA   treeData3;
    HYB4_TREE_DATA   treeData4;
    LATTICE_PROG     latticeProg;
    MKTVOL_DATA      mktVolData[2];
    HYB4_DEV_DATA    devData;
    T_CURVE     treeCurves[2][NBCRV];
    FX_DATA     fxData;
    EQ_DATA     eqData;
    ASSET_IR    assetFor;
    ASSET_IR    assetDom;
    ASSET_FX    assetFx;
    ASSET_EQ    assetEq;
    bool        cetSmoothing;

    // StatePrice Information
    TreeSliceRatesCompactSP    statePriceSlice;       /**< state Price at time point t */
    TreeSliceRatesCompactSP    tempStatePriceSlice;   // temporary slice for state price calculation


    Hyb4temp(const CClassConstSP &type);
    virtual ~Hyb4temp(void);

public:
    /********************** reflection **************************/
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    /********* exported variables **********/
protected:
    // used by all derived models
    string treeDataDebugFile; /**< print debug info for the tree if supplied */
    ZeroInterpStyle::Enum zeroInterpStyle;   /**< zero curve interpolation type */
    int nbFactorsDom;         /**< number of domestic rate factors (default 1)*/ 
    int nbFactorsFgn;         /**< number of foreign rate factors (default 1)*/ 
};

typedef smartPtr<Hyb4temp> Hyb4tempSP;

/***********************************
*
* Hyb4 FX tree mode
* Tree runs in either 2D mode or 3D
* 1D - foreign currency
* 2D - domestic currency
* 3D - FX (if tree running in 3D
*
* Based on instrument configuration, tree 
* decides if running in CUPS mode (2D)
* or full 3D mode
*
***********************************/

class TREE_DLL Hyb4CB : public Hyb4temp {
public:

    /********************** reflection **************************/
    static IObject* defaultConstructor(void) { return new Hyb4CB(); }
    static void load(CClassSP& clazz);
    static CClassConstSP const TYPE;

    Hyb4CB() : Hyb4temp(TYPE), legacyHyb4FXEQMode(false),
        corrIR(0.), corrForIRFX(0.), corrDomIRFX(0.),
        corrDomIREQ(0.), corrForIREQ(0.), corrFXEQ(0.) {}
    virtual ~Hyb4CB();

    /********* CModel interface *****************/
    virtual void getMarket(const MarketData*  market, IInstrumentCollectionSP instruments);
    virtual void retrieveFactor();

    /********* FDModel interface ****************/
    /** register with the model required currency aware zero/discounters */
    virtual void registerZero(DateTime useDate, DateTime matDate, 
                              string curveName);

    // request from model a zero/discounter slice
    virtual void getZero(DateTime useDate, DateTime matDate, string curveName, 
                         TreeSliceSP &slice);

    virtual void initModel(void);

    /********* RateTree interface ***************/
    virtual void initTreeData(void) {};
    virtual void sliceDev(TreeSlice& slice, int curveIdx) const;
    virtual void sliceEv(TreeSlice& slice) const {};
    virtual void sliceExpand(TreeSlice& treeSlice, const string& curveName) const;
    virtual void update(int t, FDProduct::UpdateType type);
    
    virtual DateTime getCurveValueDate(string curveName) const;
    virtual void print();
    virtual void recordOutput(Control* ctrl, Results* results) const;

    virtual void getFXIndex(TreeSlice& slice);
    /** populate slice with FX spot rate with resetDate/currentDate adjustments */
    virtual void getFXIndex(const IndexSpecFX& fxSpec, DateTime currentDate, DateTime resetDate,
                            TreeSlice& treeSlice);

    virtual void getEQIndex(TreeSlice& slice);
    /** populate slice with FX spot rate with resetDate/currentDate adjustments */
    virtual void getEQIndex(const IndexSpecEQ& eqSpec, DateTime currentDate, DateTime resetDate,
                            TreeSlice& treeSlice);

    virtual double getPrice0(const TreeSlice& price) const;

    /****** implement IRVegaPointise::ISensitivePoints interface *******/
    /** returns all the points on the ir vol surface to which this
    model for the supplied instrument is sensitive.
    A null return value is ok - it
    will cause an exception to be thrown (rather than crash) */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
            OutputNameConstSP  outputName,
            const CInstrument* inst) const;

protected:

    /** return internal array index number for the currency of the curve */
    virtual int getCurrIdx(const string& curveName) const;

    /** return the internal slice dimension of the curve Name - tree mode dependent */
    virtual int getCrvDim(const string& curveName) const;

    /********* local functions ****************/
    string foreignZeroName(string curveName);

    /********************** exported fields *********************/
    int ppy;                     /**< minimum tree nodes/year */
    int maxStdDeviations;        /**< nb stdDevs to cut/trim tree */
    CcyIRParamsSP forIRParams; /**< Collection of IR parameters */
    CcyIRParamsSP domIRParams; /**< Collection of IR parameters */
    CDoubleSP FXSpotVolOverride; 
    CDoubleSP EQSpotVolOverride; 
    DividendListConstSP dividendList;
    bool legacyHyb4FXEQMode;

private:
    void clear();

    // correlations for model
    double corrIR;
    double corrForIRFX;
    double corrDomIRFX;
    double corrForIREQ;
    double corrDomIREQ;
    double corrFXEQ;

    void populateTreeIRParams(Currency curr);

    /** adjust the value of a stochastic FX index to account for current date != resetDate */
    void dateAdjustFXIndex(const IndexSpecFX& fxSpec, DateTime currentDate, 
                           DateTime resetDate, TreeSlice& treeSlice);

    void dateAdjustEQIndex(const IndexSpecEQ& eqSpec, DateTime currentDate, 
                           DateTime resetDate, TreeSlice& treeSlice);

};


DRLIB_END_NAMESPACE

#endif
