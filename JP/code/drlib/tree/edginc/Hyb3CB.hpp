//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3CB.hpp
//
//   Description : temporary version of fix3 to work on the claim/bank and
//                 tree starting today without cluttering up the existing Fix3
//
//----------------------------------------------------------------------------

#ifndef _HYB3CB_HPP
#define _HYB3CB_HPP

#include "edginc/RateTree.hpp"
#include "edginc/DateTime.hpp"
#include "cupslib.h"


DRLIB_BEGIN_NAMESPACE


class TREE_DLL Hyb3CB : public RateTree {
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
        return DateTime::fromIrDate(treeData.TPDate[mTpIdx]);
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
    virtual void reset() = 0;


    /***************************** variables ********************************/

    enum {NBCRV = 3};             ///< number of curves

    DateTime today;
    DateTime currencyValueDate[2];  // foreign and domestic may be different
    DateTime valueDate;  // defined as EOD today date
    string curveName[2][NBCRV];
    string currency[2];
    string mFXName;  // FX factor name

    vector<IRDate> sortedCritDates;  // only used in debug print function

    // ???? AK: should be part of the underlying curve object
    double today2ValDateZeros[2][NBCRV];   ///< today to value date discount factors

    bool treeBuilt;  // ??? remove when decent clear function in Hyb3

    // all the internal tree structures for the various hyb3 derived modes
    HYB3_TREE_DATA   treeData;
    HYB3_DEV_DATA    devData;
    MKTVOL_DATA      mktVolData[2]; // ??? AK: should be vector ..
    T_CURVE     treeCurves[2][NBCRV]; // ??? AK: to do...
    FX_DATA     fxData;
    EQ_DATA     eqData;
    bool        cetSmoothing;

    // StatePrice Information
    TreeSliceRatesSP    statePriceSlice;       /**< state Price at time point t */
    TreeSliceRatesSP    tempStatePriceSlice;   // temporary slice for state price calculation


    Hyb3CB(const CClassConstSP &type);
    virtual ~Hyb3CB(void);

public:
    /********************** reflection **************************/
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    /********* exported variables **********/
protected:
    // used by all derived models
    string treeDataDebugFile; /**< print debug info for the tree if supplied */
    ZeroInterpStyle::Enum zeroInterpStyle;   /**< zero curve interpolation type */
};

typedef smartPtr<Hyb3CB> Hyb3CBSP;

DRLIB_END_NAMESPACE

#endif

