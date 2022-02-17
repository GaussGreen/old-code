//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FastQuoteEnv.hpp
//
//   Description : pricing environment for market making and exchange linked trading.
//                 performance is the key requirement to be close to real time..
//                 the pricing environment performs optimisation and contains cached data.
//
//   Author      : Ning Shen
//
//   Date        : 24 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef FAST_QUOTE_ENV_H
#define FAST_QUOTE_ENV_H

#include "edginc/config.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/EquityCache.hpp"

DRLIB_BEGIN_NAMESPACE

/** input for FastQuote */
class PRODUCTS_DLL FastQuoteInput : public CObject
{
public:

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    static IObject* defaultFastQuoteInput(){
        return new FastQuoteInput();
    }

    /** return control flag string to calc */
    string needCalc();

    /** throw if unsupported greeks are reqested */
    void validatePop2Object();

private:
    FastQuoteInput():  CObject(TYPE){isAmerican = true;}
 
    friend class FastQuoteEnv;

    string      optionType; //CALL , PUT, STRADDLE
    double      strike;
    DateTime    maturity;
    double      volatility;
    double      marketPrice;
    // "[I][P][D][V][R]" : I=imp vol, P=price, D=delat/gamma, V=vega, R=rho
    string      resultFlag;
    bool        isAmerican; // this is to replace the flag in FastQuoteEnv, use this one ! 
};

typedef smartPtr<FastQuoteInput> FastQuoteInputSP;
typedef array<FastQuoteInputSP, FastQuoteInput> FastQuoteInputArray;
typedef smartPtr<FastQuoteInputArray> FastQuoteInputArraySP;

///////////////////////////////////////////////////
// FastQuoteSpline class
///////////////////////////////////////////////////
class PRODUCTS_DLL FastQuoteSpline
{
public:
    string resultFlag;
    string resultCcy;
    double volCached;
    double spotBot;
    double spotTop;

    StringArray          packetNames;
    OutputNameArray      outNames;

    DoubleArrayArray y2;
    DoubleArray spotGrid;
    DoubleArrayArray valueGrid;

    /** interpolate a price  */
    bool canInterpolate(const string& req, double spot, double vol);

    /** store and prepare spline interpolation */
    void prepareSpline(const string& inFlags, double vol, const CResultsArray& in);

    /** interpolate a price  */
    void interpolate(double spot, CResults& out);

    FastQuoteSpline()
    {
        resultFlag = "";
        spotBot = 0.0;
        spotTop = 0.0;
        volCached = 0.0;
    }

};

///////////////////////////////////////////////////
// FastImpVol  class
///////////////////////////////////////////////////
class PRODUCTS_DLL FastImpVol
{
public:
    double oldVol;
    double oldPrice;
    double oldVega;
    double oldSpot;

    /** compute implied vol. returns true if result converged, false if not */
    bool FastImplied(double spot, double price, double volGuess, CInstrument* inst, IModel* model,
                     CAsset* asset, CResults& out);

    FastImpVol()
    {
        oldVol = 0.3;
        oldPrice = -1.0;
        oldVega = -1.0;
        oldSpot = -1.0;
    }

    FastImpVol(double vol, double price, double vega, double spot)
    {
        oldVol = vol;
        oldPrice = price;
        oldVega = vega;
        oldSpot = spot;
    }
};

/////////////////////////////////////////////////////////////
// FastQuoteEnv class
/////////////////////////////////////////////////////////////
class PRODUCTS_DLL FastQuoteEnv: public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    static IObject* defaultFastQuoteEnv(){
        return new FastQuoteEnv();
    }
    
    void validatePop2Object();

    /** call model for each input */
    void Price(FastQuoteInput& in, CResults& out);

    /**  initialise pricing environment  */
    void init(double spot);

    /** compute spline and result */
    void computeSpline(IModelSP& model, CControlSP& ctrl, double spot, double vol, const string& resultFlag, 
                       CResults& out, FastQuoteSpline& store);

protected:
    //** set input to the required array element */
    void setInput(const FastQuoteInput& input);

private:
    FastQuoteEnv(const FastQuoteEnv& rhs);
    FastQuoteEnv& operator=(const FastQuoteEnv& rhs);

    CMarketDataSP   market;
    IModelSP        model;
    CVanillaSP      baseInst; // base instrument
    CDateTime       lastMaturity; // longest maturty for building data cache
    double          spotRef; // spotRef the pricing environment is set up for
    bool            isAmerican;  // to be removed
    /** false=use flat vol, true=use term structure and offset */
    bool            isFlatVol;

    /** false=use input vol for price and sensitivities, 
        true(default)=use implied vol if requested, input vol if IV not requested */
    bool IVpricing;

    /** false(default)=do not use spline interp
        true = cache and use spline interpolation when possible */
    bool useInterp;

    bool useFastImplied; // true =used fat implied method implemented here, false(default) = use standard object

    /** for reflection */
    FastQuoteEnv():  CObject(TYPE){
        isFlatVol = true;
        spotRef=-1;
        needInit = true;
        IVpricing = true;
        useInterp = false;
        useFastImplied = false;
        closedForm = IModelSP(new CClosedFormLN("VolPreferred"));
        isAmerican = true; // to be removed
    }

    /** unregistered */

    // if init needed
    bool needInit; // $unregistered

    // to do:
    //vector<DoubleArray> discFactors;
    //vector<DoubleArray> settlePV; // pv factor for settlement

    // only used internally
    // this is the equity that we can change spot and vol directly and using cached fwds
    EquityCacheSP       assetCache; // $unregistered

    CVanillaSP          baseInstCopy; // copy used in init() so that original can be serialised correctly. $unregistered

    // unregistered
    map<string, FastQuoteSpline> resultCache; // for spline interp price and greeks $unregistered
    map<string, FastImpVol> impVolCache; // for fast implied vol $unregistered

    // a closed form model for european
    IModelSP closedForm; // $unregistered
};

typedef smartPtr<FastQuoteEnv> FastQuoteEnvSP;

DRLIB_END_NAMESPACE

#endif
