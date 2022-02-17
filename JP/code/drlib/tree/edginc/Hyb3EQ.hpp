//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3Eq.hpp
//
//   Description : hyb3 domestic equity mode FD model class
//
//----------------------------------------------------------------------------

#ifndef HYB3EQ_HPP
#define HYB3EQ_HPP

#include "edginc/Hyb3.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/DividendList.hpp"
#include "cupsmodl.h"


DRLIB_BEGIN_NAMESPACE


/***********************************
*
* Hyb3 Domestic Equity tree mode
* Equity is in domestic currency
* Tree runs in 2D mode 
* 1D - domestic currency
* 2D - equity
*
***********************************/

class TREE_DLL Hyb3EQ : public Hyb3 {
public:

    /********************** reflection **************************/
    static IObject* defaultConstructor(void) { return new Hyb3EQ(); }
    static void load(CClassSP& clazz);
    static CClassConstSP const TYPE;

    Hyb3EQ();
    ~Hyb3EQ();

    virtual void initModel(void);
    virtual void getMarket(const MarketData*  market, IInstrumentCollectionSP instruments);
    virtual void retrieveFactor();
    virtual void initTreeData(void) {};
    virtual void sliceDev(TreeSlice& slice, int curveIdx) const;
    virtual void sliceExpand(TreeSlice& treeSlice, const string& curveName) const;

    virtual DateTime getCurveValueDate(string curveName) const;
    virtual void print();

    virtual void getEQIndex(TreeSlice& slice);

    /** retrieve zeroBank Mode for tree */
    virtual ZeroBankMode::Enum getZeroBankMode() const { return zeroBankMode; }


    /****** implement IRVegaPointise::ISensitivePoints interface *******/
    /** returns all the points on the ir vol surface to which this
    model for the supplied instrument is sensitive.
    A null return value is ok - it
    will cause an exception to be thrown (rather than crash) */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
            OutputNameConstSP  outputName,
            const CInstrument* inst) const {
        return IRGridPointAbsArraySP(   );
    }

    /********************** exported fields *********************/
    int ppy;                     /**< minimum tree nodes/year */
    int maxStdDeviations;        /**< nb stdDevs to cut/trim tree */
    CcyIRParamsSP IRParams;      /**< Collection of IR parameters */
    string EQSmileParams;
    string IREQCorrelation;      /**< correlation override */
    double EQCalibrationStrike;  /**< strike to calibrate to for equity vols */
    string zeroInterpStyle;      /**< zero curve interpolation type */
    DateTime lastDividendDate;   /**< stop inserting dividend dates after this date. Typically the product maturity */
    double *eqSpotVolOverride;

protected:
    DividendListConstSP dividendList;

private:
    void populateTreeIRParams();

    double corrIREQ;  // correlation for model
    string mEQName;
};


/*************************************
*
* Hyb3 FX Equity tree mode
* covers equity modes:
*   1) domestic equity, 2 (cups) IRs
*   2) foreign equity, 2 (cups) IRs
*   3) composite equity, 2 (cups) IRs
* Tree runs in 3D mode for all cases
* 1D - foreign currency
* 2D - domestic currency
* 3D - equity
*
***********************************/
#if 0
class TREE_DLL Hyb3FXEQ : public Hyb3 {
public:

    /********************** reflection **************************/
    static IObject* defaultConstructor(void) { return new Hyb3FXEQ(); }
    static void load(CClassSP& clazz);
    static CClassConstSP const TYPE;

    Hyb3FXEQ() : Hyb3(TYPE) {}
    virtual ~Hyb3FXEQ() {};

    virtual void initModel(void);
    virtual void getMarket(const MarketData*  market, IInstrumentCollectionSP instruments);
    virtual void retrieveFactor();
    virtual void initTreeData(void);

    /****** implement IRVegaPointise::ISensitivePoints interface *******/
    /** returns all the points on the ir vol surface to which this
    model for the supplied instrument is sensitive.
    A null return value is ok - it
    will cause an exception to be thrown (rather than crash) */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
            OutputNameConstSP  outputName,
            const CInstrument* inst) const {
        return IRGridPointAbsArraySP(   );
    }

    // correlations for model
    double corrIR;
    double corrForIRFX;
    double corrDomIRFX;
    double corrForIREQ;
    double corrDomIREQ;
    double corrFXEQ;
    
    /********************** exported fields *********************/
    int ppy;                   /**< minimum tree nodes/year */
    double maxStdDeviations;   /**< nb stdDevs to cut/trim tree */
    CcyIRParamsSP forIRParams; /**< Collection of IR parameters */
    CcyIRParamsSP domIRParams; /**< Collection of IR parameters */
    string EQSmileParams;
    string IRFXEQCorrelations; /**< correlation override */
    string zeroInterpStyle;    /**< zero curve interpolation type */

private:
    void populateTreeIRParams(Currency curr);
    
};
#endif

DRLIB_END_NAMESPACE

#endif
