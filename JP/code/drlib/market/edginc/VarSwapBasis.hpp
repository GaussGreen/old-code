//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VarSwapBasis.hpp
//
//   Description : 
//
//   Author      : Zhijiang Huang
//
//   Date        : 16 August 2004
//
//
//----------------------------------------------------------------------------

#ifndef VARSWAPBASIS_HPP
#define VARSWAPBASIS_HPP

#include "edginc/config.hpp"
#include "edginc/Asset.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/VolRequest.hpp"
#include "edginc/YieldCurve.hpp"


DRLIB_BEGIN_NAMESPACE

class VarSwapBasis;
typedef smartPtr<VarSwapBasis> VarSwapBasisSP;

/** New market data for Variance Swaps.
    Contains a delta cutoff level and a term structure of volatility 
    basis points to be added to the Variance Swap Fair Volatility */
class MARKET_DLL VarSwapBasis: public MarketObject {
public:
    /** Fetches a VarSwapBasis from the market data. If it does not exist
        a NULL is returned */
    static VarSwapBasisSP fetchBasis(const IModel*     model,
                                     const MarketData* market,
                                     const string&     name);
    
    static CClassConstSP const TYPE;

    class VarSwapBasisRequest;

    ////////////////////////////////////////////////////////////////////////////////////////
    
    /** Vol Processed for VarSwapBasis Interpolant. 
        Interpolates Vol Basis points and Delta strikes between benchmark dates */
    class MARKET_DLL VarSwapBasisProc: public CObject,
                            public virtual IVolProcessed {
    public:
        
        static CClassConstSP const TYPE;
        
        /** Full constructor */
        VarSwapBasisProc(const VarSwapBasis* basis,
                         const CAsset* asset, 
                         const VolSurface* backbone,
                         const VarSwapBasisRequest* req,
                         const string& tenorTime);

        /** Computes by interpolation the basis for a particular maturity */
        double interpVolBasis(double time) const;

        /** Computes by interpolation the strike for a given delta for a particular maturity */
        double interpPriceCutoff(double time) const;

        /** Computes by interpolation the strike for a given delta for a particular maturity */
        double interpSkewCutoff(double time) const;

        /** Records delta cutoff and basis in Results */
        void recordRequests(double tradYear, Control* control, CResults* results);


        // IVolProcessed Methods
        
        /** Vol name */
        virtual string getName() const;
        
        /** calculates the trading time between two dates */
        virtual double calcTradingTime(const DateTime &date1, 
                                       const DateTime &date2) const;
        /** retieve time measure for the vol */
        virtual TimeMetricConstSP GetTimeMetric()const;

    private:

        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);

        static IObject* defaultVarSwapBasisProc();

        /** Default constructor */
        VarSwapBasisProc();
        
        bool                            tweakDeltaCutoff;       //!< Whether to recompute the cutoff in tweaks
        int                             methodId;               //!< Cutoff method
        
        LinearInterpolantNonVirtualSP   strikeCutoffInterp;     //!< Interpolant of absolute strike cutoffs
        LinearInterpolantNonVirtualSP   volBasisInterp;         //!< Interpolant of input vol basis points along the time line
        TimeMetricSP                    metric;                 //!< Time metric $unregistered
    };

    typedef smartPtr<VarSwapBasisProc> VarSwapBasisProcSP;

    
    ////////////////////////////////////////////////////////////////////////////////////////
    
    /** Vol Processed when computation failed */
    class MARKET_DLL VarSwapBasisProcError: public CObject,
                                 public virtual IVolProcessed {
    public:
        static CClassConstSP const TYPE;

        /** Vol name */
        virtual string getName() const;
        
        /** calculates the trading time between two dates */
        virtual double calcTradingTime(const DateTime &date1, 
                                       const DateTime &date2) const;
        /** retieve time measure for the vol */
        virtual TimeMetricConstSP GetTimeMetric()const;

        /** Default constructor */
        VarSwapBasisProcError(const ModelException& e);

        /** Allows access to exception */
        const ModelException& getException() const;

    private:
        ModelException e;           //!< Error $unregistered
        
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);
    };

    typedef smartPtr<VarSwapBasisProcError> VarSwapBasisProcErrorSP;

    
    ////////////////////////////////////////////////////////////////////////////////////////


    /** Vol request for basis VarSwapBasis Interpolant */
    class MARKET_DLL VarSwapBasisRequest: public CVolRequest {
    public:
        static CClassConstSP const TYPE;
        
        /** Full Constructor */
        VarSwapBasisRequest(const DateTimeArray& dates, YieldCurveConstSP yc);

        /** Gives access to the dates */
        const DateTimeArray& getDates() const;

        /** Gives access to YC */
        YieldCurveConstSP getAssetYC() const;

    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);

        DateTimeArray       dates;      //!< Relevant dates
        YieldCurveConstSP   yc;         //!< Asset YC for domestic currency
    };

    typedef smartPtr<VarSwapBasisRequest> VarSwapBasisRequestSP;


    ////////////////////////////////////////////////////////////////////////////////////////


    friend class VarSwapBasisProc;

    /** Vol name */
    virtual string getName() const;

    /** Computes by interpolation the basis for a particular maturity */
    double interpVolBasis(double time) const;
    
    /** Computes by interpolation the strike for a given delta for a particular maturity */
    double interpPriceCutoff(double time) const;

    /** Computes by interpolation the strike for a given delta for a particular maturity */
    double interpSkewCutoff(double time) const;
    
    /** Sets up the interpolant given the asset information */
    VarSwapBasisProc* setup(const CAsset* asset, 
                            const VolSurface* backbone, 
                            const VarSwapBasisRequest* req, 
                            const string& tenorTime);

    /** Validation and population of transient fields */
    void validatePop2Object();

    virtual ~VarSwapBasis();

    IObject* clone() const;

private:

    // Methodology
    static const string NONE;
    static const string DEFAULT;
    static const string SKEW_CUTOFF;
    static const string PRICE_CUTOFF;
    static const string BASIS_ONLY;

    // Methodology Ids for efficiency
    static const int NONE_ID;
    static const int SKEW_CUTOFF_ID;
    static const int PRICE_CUTOFF_ID;
    static const int BASIS_ONLY_ID;
    
    /** Performs conversion of delta denominated strike to absolute strike */
    class DeltaImpliedStrike;
    
    // OPTIONAL FIELDS FROM PYRAMID
    string              name;
    
    bool                tweakDeltaCutoff;       //!< Whether to recompute the cutoff in tweaks
    double              deltaCutoff;            //!< Delta cutoff level
    string              method;                 //!< Type of cutoff i.e. NONE, DEFAULT, SKEW_CUTOFF, PRICE_CUTOFF
    StringArray         volBasisBMs;            //!< Basis points benchmarks
    DoubleArray         volBasis;               //!< Vol Basis points at each benchmark

    // TRANSIENT FIELDS
    int                 methodId;               //!< Integer representation of method for speed
    VarSwapBasisProcSP  basisProc;              //!< Interpolant for strikes & vols

    /** Default constructor */
    VarSwapBasis();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static IObject* defaultVarSwapBasis();
};

DRLIB_END_NAMESPACE

#endif
