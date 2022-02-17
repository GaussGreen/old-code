//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Fix3.hpp
//
//   Description : Fix3 tree
//
//----------------------------------------------------------------------------

#ifndef _FIX3_HPP
#define _FIX3_HPP

#include "edginc/RateTree.hpp"
#include "edginc/DateTime.hpp"
#include "fix123head.h"


DRLIB_BEGIN_NAMESPACE

class TREE_DLL Fix3 : public RateTree {
public:

    /********* CModel interface *****************/

    virtual void getMarket(const MarketData*  market, IInstrumentCollectionSP instruments);

    /********* FDModel interface ****************/

    /** collect model initialisation data, set up timeline  */
    virtual void initModel(void);

    /** finalise model initialisation, create payoffIndex  */
    virtual void finaliseModel(CControl* control);

    /** initialise tree */
    void initTreeData(void);

    /** register with the model required zero/discounters */
    virtual void registerZero(
        DateTime obsDate, 
        DateTime matDate, 
        string curveName);

    void getZero(
        DateTime obsDate, 
        DateTime matDate, 
        string curveName, 
        TreeSliceSP &slice);

protected:

    /** retrieve market data into model */
    virtual void retrieveFactor();

public:

    /********* RateTree interface ***************/

    /** insert into the tree a product configured zeroBank if tree is running in zerobank mode */
    virtual void insertNamedZeroBank(NamedZeroBank& zeroBank);

    /** register with the tree that the ParYield defined by rateSpec will be required */
    virtual void insertIRIndex(const IndexSpecIR& rate, DateTime date);

    /** populate the slice with the ParYield value associated with the rateSpec */
    virtual void getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, TreeSlice& slice);

    /** populate the slice with the "current" FX spot, which will always take value 1.0 across the slice */
    virtual void getFXIndex(TreeSlice& slice);

    /** set tree mode to be either zerobank or claim bank */
    virtual ZeroBankMode::Enum getZeroBankMode() const { return zeroBankMode; }

    /** insert zerobank dates/details into the tree */
    virtual void insertZero(const ZeroBond& zeroBond);

    /** get from tree the slice corresponding to zeroIndexSpec */
    virtual void getZero(const ZeroBond& zero, int step, TreeSlice& slice);

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

protected:

    virtual void update(int t, FDProduct::UpdateType type);

    virtual void updateFinalize(int t, FDProduct::UpdateType type);

    virtual void sliceDev(TreeSlice& slice, int curveIdx) const;

    virtual void sliceEv(TreeSlice& slice) const;

    virtual void sliceExpand(TreeSlice& slice, const string& curveName) const;

    virtual DateTime getDate(int step) const;

    virtual DateTimeArray getDates() const;

    virtual int getLastStep() const;

    virtual DateTime getCurveValueDate(string curveName) const;

    virtual DateTime getToday() const;

    virtual void print();

    virtual int getCurveIdx( const string & curveName ) const;

    /********************** local definitions  **************************/
public:
    static CClassConstSP const TYPE;
    virtual ~Fix3(void);

private:
    void clear();
    int getCrvIdx(const string& curveId) const;
    void getDiffusionTCurve(T_CURVE& diffusionCurve) const;

    Fix3(const CClassConstSP &type = TYPE);
    static IObject* defaultConstructor(void) { return new Fix3(); }
    static void load(CClassSP& clazz); /* declare what is being exposed */

    /***************************** variables ********************************/
private:
    enum {NBCRV = 3};             ///< number of curves

    IRDate    mToday;
    IRDate    mValueDate;

    MKTVOL_DATA    mMktVolData;           ///< market volatility data
    IR_SIM         mIRSim;                ///< IR_SIM data

    string mCurveName[NBCRV];
    string mCurrency;   // pricing/domestic currency

    typedef map<string, ZeroBondProdSP> ZeroProdMap;
    ZeroProdMap zeroProdList;

    // ZEROBANK mode fields
    map<string, NamedZeroBank> mZeroBank;
    map<string, ZeroBond> mZeroBond;

    // ??? temp to stop deleting memory for unbuilt tree
    // until I improve the clear function
    bool mTreeBuilt;

protected:
    FIX3_TREE_DATA mTreeData;             ///< tree data
    FIX3_DEV_DATA  mDevData;              ///< dev data
    T_CURVE     mTCurves[NBCRV];         ///< zero curves

    /********* exported variables **********/
protected:
    int nbFactors;            /**< Number of factors */
    CcyIRParamsSP IRParams;   /**< Collection of IR parameters */
    string treeDataDebugFile; /**< print debug info for the tree if supplied */
    ZeroInterpStyle::Enum zeroInterpStyle;   /**< zero curve interpolation type */
    ZeroBankMode::Enum   zeroBankMode;        ///< zero mode - "traditional" zerobank or claimbank
    bool cetSmoothing;
};
typedef smartPtr<Fix3> Fix3SP;

DRLIB_END_NAMESPACE

#endif
