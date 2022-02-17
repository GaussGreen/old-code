//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Fix3CB.hpp
//
//   Description : temporary version of fix3 to work on the claim/bank and
//                 tree starting today without cluttering up the existing Fix3
//
//----------------------------------------------------------------------------

#ifndef _FIX3CB_HPP
#define _FIX3CB_HPP

#include "edginc/RateTree.hpp"
#include "edginc/DateTime.hpp"
#include "fix123head.h"


DRLIB_BEGIN_NAMESPACE

class TREE_DLL Fix3CB : public RateTree {
public:

    /********* CModel interface *****************/

    virtual void getMarket(const MarketData*  market, IInstrumentCollectionSP instruments);
    /** value date as defined by the model */
    virtual DateTime getValueDate() const { return valueDate; }

    /********* FDModel interface ****************/

    /** collect model initialisation data, set up timeline  */
    virtual void initModel(void);

    /** finalise model initialisation - nothing to do*/
    virtual void finaliseModel(CControl* control) {};

    /** initialise tree */
    void initTreeData(void);

    /** register with the model required zero/discounters */
    virtual void registerZero(DateTime useDate, DateTime matDate, 
                              string curveName);

    /** request from model a zero/discounter slice */
    virtual void getZero(DateTime useDate, DateTime matDate, string curveName, 
                         TreeSliceSP &slice);

    /** today date defined as market reference/today date, including time (supports SOD or EOD)*/
    virtual DateTime getToday() const;

    /** this forwards the price back from today to the valueDate as is the rates convention */
    virtual double getPrice0(const TreeSlice& price) const;

    /** model to record any additional information requested */
    virtual void recordOutput(Control* ctrl, Results* results) const;

protected:

    /** retrieve market data into model */
    virtual void retrieveFactor();

public:

    /********* RateTree interface ***************/

    /** populate the slice with the ParYield value associated with the rateSpec */
    virtual void getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, 
                            DateTime resetDate, TreeSlice& treeSlice);
        
    /** populate the slice with the "current" FX spot, which will always take value 1.0 across the slice */
    virtual void getFXIndex(TreeSlice& slice);

    /** set tree mode to be either zerobank or claim bank */
    virtual ZeroBankMode::Enum getZeroBankMode() const { return ZeroBankMode::CLAIMBANK; }

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    virtual MarketDataFetcherSP createMDF() const;

    /****** implement IRVegaPointise::ISensitivePoints interface *******/
    /** returns all the points on the ir vol surface to which this
    model for the supplied instrument is sensitive.
    A null return value is ok - it
    will cause an exception to be thrown (rather than crash) */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst) const;

    /** adjust the value of a stochastic IR index to account for currentDate != resetDate */
    virtual void dateAdjustIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, 
                                   DateTime resetDate, TreeSlice& treeSlice);

protected:

    virtual void update(int t, FDProduct::UpdateType type);

    virtual void updateFinalize(int t, FDProduct::UpdateType type);

    virtual void sliceExpand(TreeSlice& treeSlice, const string& curveName) const;

    virtual void sliceDev(TreeSlice& slice, int curveIdx) const;

    virtual void sliceEv(TreeSlice& slice) const;

    virtual DateTime getDate(int step) const;

    virtual DateTimeArray getDates() const;

    virtual int getLastStep() const;

    virtual DateTime getCurveValueDate(string curveName) const;

    virtual int getCurveIdx( const string & curveName ) const;

    virtual void print();

    /********************** local definitions  **************************/
public:
    static CClassConstSP const TYPE;
    virtual ~Fix3CB(void);

private:


    static IObject* defaultConstructor(void) { return new Fix3CB(); }
    static void load(CClassSP& clazz); /* declare what is being exposed */

private:

    void clear();

    int getCrvIdx(const string& curveId) const;

    void configureClaimBank(const CRIT_DATE* critDate, int nbCritDate);

    /***************************** variables ********************************/

    enum {NBCRV = 3};             ///< number of curves

    DateTime today;           // date the tree starts diffusing from
    DateTime valueDate;       // date trade is discounted/priced (currently = EOD today)
    DateTime curveValueDate;  // spot date of the curve
    string curveName[NBCRV];
    string currency;

    vector<IRDate> sortedCritDates;  // only used in debug print function

    double today2ValDateZeros[NBCRV];   ///< today to value date discount factors

    bool treeBuilt;  // ??? remove when decent clear function in fix3

    FIX3_TREE_DATA treeData;           ///< tree data
    FIX3_DEV_DATA  devData;            ///< dev data
    T_CURVE        treeCurves[NBCRV];  ///< zero curves
    MKTVOL_DATA    mktVolData;         ///< market volatility data

    Fix3CB(const CClassConstSP &type = TYPE);

    /********* exported fields **********/
protected:
    int nbFactors;            /**< Number of factors */
    CcyIRParamsSP IRParams;   /**< Collection of IR parameters */
    string treeDataDebugFile; /**< print debug info for the tree if supplied */
    ZeroInterpStyle::Enum zeroInterpStyle;   /**< zero curve interpolation type */
    bool cetSmoothing;
};
typedef smartPtr<Fix3CB> Fix3CBSP;

DRLIB_END_NAMESPACE

#endif

