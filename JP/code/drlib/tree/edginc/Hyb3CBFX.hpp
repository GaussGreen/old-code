//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3FX.hpp
//
//   Description : hyb3 FX mode (Foreign and Domestic currency + FX)
//
//----------------------------------------------------------------------------

#ifndef HYB3CBFX_HPP
#define HYB3CBFX_HPP

#include "edginc/Hyb3CB.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "cupsmodl.h"


DRLIB_BEGIN_NAMESPACE


/***********************************
*
* Hyb3 FX tree mode
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

class TREE_DLL Hyb3CBFX : public Hyb3CB {
public:

    /********************** reflection **************************/
    static IObject* defaultConstructor(void) { return new Hyb3CBFX(); }
    static void load(CClassSP& clazz);
    static CClassConstSP const TYPE;

    Hyb3CBFX() : Hyb3CB(TYPE), legacyHyb3FXMode(false), fxIndexSpecSupplied(false), 
        corrIR(0.), corrForIRFX(0.), corrDomIRFX(0.) {}
    virtual ~Hyb3CBFX();

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
    /** initialise all tree structures */
    virtual void reset() ;

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
    bool legacyHyb3FXMode;
    
private:
    void clear();

    bool fxIndexSpecSupplied;

    // correlations for model
    double corrIR;
    double corrForIRFX;
    double corrDomIRFX;

    void populateTreeIRParams(Currency curr);

    /** adjust the value of a stochastic FX index to account for current date != resetDate */
    void dateAdjustFXIndex(const IndexSpecFX& fxSpec, DateTime currentDate, 
                           DateTime resetDate, TreeSlice& treeSlice);

};


DRLIB_END_NAMESPACE

#endif
