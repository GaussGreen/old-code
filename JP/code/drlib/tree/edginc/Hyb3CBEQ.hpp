//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3Eq.hpp
//
//   Description : hyb3 domestic equity mode FD model class
//
//----------------------------------------------------------------------------

#ifndef HYB3CBEQ_HPP
#define HYB3CBEQ_HPP

#include "edginc/Hyb3CB.hpp"
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

class TREE_DLL Hyb3CBEQ : public Hyb3CB {
public:

    /********************** reflection **************************/
    static IObject* defaultConstructor(void) { return new Hyb3CBEQ(); }
    static void load(CClassSP& clazz);
    static CClassConstSP const TYPE;

    Hyb3CBEQ();
    virtual ~Hyb3CBEQ();

    virtual void initModel(void);
    virtual void getMarket(const MarketData*  market, IInstrumentCollectionSP instruments);
    virtual void retrieveFactor();
    virtual void initTreeData(void) {};
    virtual void sliceDev(TreeSlice& slice, int curveIdx) const;
    virtual void sliceExpand(TreeSlice& treeSlice, const string& curveName) const;

    //virtual void update(int t, FDProduct::UpdateType type);

    virtual void sliceEv(TreeSlice& slice) const {};

    virtual DateTime getCurveValueDate(string curveName) const;
    virtual void print();
    virtual void recordOutput(Control* ctrl, Results* results) const;

    virtual void getEQIndex(TreeSlice& slice);
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
    
    /** initialise all tree structures */
    virtual void reset() ;

    /********************** exported fields *********************/
    int ppy;                     /**< minimum tree nodes/year */
    int maxStdDeviations;        /**< nb stdDevs to cut/trim tree */
    CcyIRParamsSP IRParams;      /**< Collection of IR parameters */
    double *EQCalibrationStrike; /**< strike to calibrate to for equity vols */
    DateTime lastDividendDate;   /**< stop inserting dividend dates after this date. Typically the product maturity */

    // ??? it is not clear if these should remain for audit reasons
    double *eqSpotVolOverride;
    double *IREQCorrelationOverride;

protected:
    DividendListConstSP dividendList;

private:
    void clear();
    void populateTreeIRParams();

    double corrIREQ;  // correlation for model
    string eqName;
};


DRLIB_END_NAMESPACE

#endif
