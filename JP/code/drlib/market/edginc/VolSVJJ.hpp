//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSVJJ.hpp
//
//   Description : Heston + StockOnlyJump (Merton) + VolOnlyJump + CommonStockVolJump
//
//   Date        : 20 Nov 2002
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLSVJJ_HPP
#define EDR_VOLSVJJ_HPP

#include "edginc/IDynamicsParameter.hpp"
#include "edginc/VolBaseParam.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Heston.hpp"
#include "edginc/MarketDataConvert.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VolSV)
FORWARD_DECLARE(VolSVJ)
FORWARD_DECLARE(VolSVCJ)

/** A parameterised CVolBase (does not support BS and DVF views of the world for now) */
class MARKET_DLL VolSVJJ: public VolBaseParam,
               public virtual Calibrator::IAdjustable,
               public virtual IDynamicsParameter,
               public virtual MarketDataConvert::IConvert<VolSV>,
               public virtual MarketDataConvert::IConvert<VolSVJ>,
               public virtual MarketDataConvert::IConvert<VolSVCJ>
/*               virtual public IVolatilityBS,
               virtual public IVolatilityDVF */ {
public:
    /** Full constructor */
    VolSVJJ(const string& name,
            double initialVol,
            double meanVol,
            double meanReversRate,
            double volVol,
            double correlation,
            double stockCrashRate,
            double stockCrashSizeMean,
            double stockCrashSizeUncertainty,
            double volCrashRate,
            double volCrashSizeMean,
            double commonCrashRate,
            double commonStockCrashSizeMean,
            double commonStockCrashSizeUncertainty,
            double commonVolCrashSizeMean,
            double stockVolCrashSizeCorrelation);

    /** Downgrade to VolSV if jump parameters are zero*/
    virtual VolSVSP convert(VolSV* p = 0) const;

    /** Downgrade to VolSVJ if vol jump parameters are zero*/
    virtual VolSVJSP convert(VolSVJ* p = 0) const;

    /** Downgrade to VolSVCJ if idiosyncratic jump parameters are zero*/
    virtual VolSVCJSP convert(VolSVCJ* p = 0) const;

private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class MARKET_DLL SVJJVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        SVJJVolParam(): CVolParam(TYPE){}

        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const;

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface* spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const;

    private:
        static void load(CClassSP& clazz){
            REGISTER(SVJJVolParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };

    /* Methods called by SVJJVolParam. Implemented on VolSVJJ to avoid dereferencing */ 
    void ComputeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const;

    VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray&   strikes) const;

    /** Calculates the components alpha and beta that appear in 
        the exponent of the time-t joint Bilateral Laplace 
        of X_T = (Y_T, V_T, I_T) where Y_T = ln(S_T / F(0, T)) is the dimension-less
        log spot at time T, V_T is the instantaneous variance at time T
        and I_T is the instantaneous variance from 0 to time T
        The Bilateral Laplace transform is of the form
            exp(alpha(tau, u1, u2, u3) 
                + u1 * Y_t 
                + beta(tau, u1, u2, u3) * V_t 
                + u3 * I_t)
        where tau = T - t and the complex numbers u1, u2, u3 are the frequencies 
        wrt Y_T, V_T and I_T, respectively. */
    void calcJointLapAlphaBeta(double         tau,    // tau = T - t (in years)
                               const Complex& u1,
                               const Complex& u2,
                               const Complex& u3,
                               Complex&       alpha,
                               Complex&       beta) const;

    // registered fields
    double initialVol;
    double correlation;
    double volVol;
    double meanVol;
    double meanReversRate;
    /** The jump process is driven by 3 Poisson processes. The first Poisson drives 
       jumps in the stock only. The second Poisson drives jumps in the vol only. 
       The third Poisson drives common jumps to the stock and the vol.
       The stock-only jump size is lognormally-distributed. The vol-only jump size is 
       exponentially distributed. The common vol jump size is exponentially distributed
       and conditional on that jump the common stock jump size is lognormally ditributed. */
    double stockCrashRate;                  // Stock (only) jump rate
    double stockCrashSizeMean;              // Stock (only) jump size mean
    double stockCrashSizeUncertainty;       // Stock (only) jump size std dev
    double volCrashRate;                    // Stock (only) jump rate
    double volCrashSizeMean;                // Stock (only) jump size mean
    double commonCrashRate;                 // Common jump rate
    double commonStockCrashSizeMean;        // Common stock jump size (conditional/unconditional) mean 
                                            // @ zero vol jump size mean and/or zero crash correlation
    double commonStockCrashSizeUncertainty; // Common stock jump size conditional std dev
    double commonVolCrashSizeMean;          // Vol jump size (unconditional) mean
    double stockVolCrashSizeCorrelation;    // 'Correlation' parameter between stock and vol common jump sizes
                                            // Not correlation exactly, but drives it. Note: lives in the real line

    // transient fields (won't appear in dd interface)
    DateTime      baseDate;
    HestonSP              diffusion;
    MertonCrashSP         stockCrash;
    Heston::VolCrashSP    volCrash;
    Heston::CommonCrashSP commonCrash;
        
    /* Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VolSVJJ();

    static IObject* defaultCtor();

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache();

public:
    friend class SVJJVolParam;

    static CClassConstSP const TYPE;

    void validatePop2Object();

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const;

    /* Methods called by FourierProcessSVJJ. Implemented on VolSVJJ to
       avoid dereferencing */
    Complex scalelessCumulant(const StFourierProcessLogRtn& process,
                              const StFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;
    Complex scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                              const FwdStFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;
    Complex cumulant(const StFourierProcessIntVar& process,
                     const StFourierProductIntVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;
    Complex cumulant(const FwdStFourierProcessIntVar& process,
                     const FwdStFourierProductIntVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const;

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields);
private:
    void update();
};
typedef smartPtr<VolSVJJ> VolSVJJSP;

typedef MarketDataConvert::IConvert<VolSVJJ> VolSVJJConvert;

DRLIB_END_NAMESPACE
#endif
