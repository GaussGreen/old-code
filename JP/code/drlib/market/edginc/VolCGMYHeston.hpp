//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolCGMYHeston.hpp
//
//   Description : 
//
//   Date        : 11 November 2002
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLCGMYHESTON_HPP
#define EDR_VOLCGMYHESTON_HPP

#include "edginc/IDynamicsParameter.hpp"
#include "edginc/VolBaseParam.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Range.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Heston.hpp"
DRLIB_BEGIN_NAMESPACE

// helper class
class MARKET_DLL CGMYe: public CObject{
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object();

    Complex calcLaplaceExponent(double         tau,
                                const Complex& u) const;

    CGMYe(double  cNegative,
          double  cPositive,
          double  g,
          double  m,
          double  yNegative,
          double  yPositive,
          double  gaussianVol);

private:
    CGMYe(const CGMYe& rhs);
    CGMYe& operator=(const CGMYe& rhs);

    CGMYe();

    static IObject* defaultCtor();

    static void load(CClassSP& clazz);

    // registered fields
    double cNegative;
    double cPositive;
    double g;
    double m;
    double yNegative;
    double yPositive;
    double gaussianVol;
    // transient field
    double cNegativeTimesGammaYNegative;
    double cPositiveTimesGammaYPositive;
    double mToThePowerYPositive;
    double gToThePowerYNegative;
    double sqGaussianVol;
};
typedef smartPtr<CGMYe> CGMYeSP;

/** A parameterised CVolBase (does not support BS and DVF views of the world for now) */
class MARKET_DLL VolCGMYHeston: public VolBaseParam,
                     public virtual Calibrator::IAdjustable,
                     public virtual IDynamicsParameter
                     /*,
                     virtual public IVolatilityBS,
                     virtual public IVolatilityDVF */ {
private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class MARKET_DLL VolCGMYHestonParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        VolCGMYHestonParam(): CVolParam(TYPE){}

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
            REGISTER(VolCGMYHestonParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };

    /* Methods called by VolCGMYHestonParam. Implemented on VolCGMYHeston to avoid dereferencing */ 
    void ComputeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const;

    VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray&   strikes) const;

    // registered fields
    double initialVol;
    double volVol;
    double meanVol;
    double meanReversRate;
    
    double cNegative;
    double cPositive;
    double g;
    double m;
    double yNegative;
    double yPositive;
    double gaussianVol;

    // transient fields (won't appear in dd interface)
    DateTime     baseDate;
    CGMYeSP      cgmye;
    HestonSP     heston;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VolCGMYHeston();

    static IObject* defaultCtor();

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache();

public:
    friend class VolCGMYHestonParam;
    friend class CGMYHestonModelAddin;

    static CClassConstSP const TYPE;

    virtual ~VolCGMYHeston();

    static void calcCumulantComponents(double          tau,    // tau = T - t (in years)
                                       const Complex&  u,
                                       const CGMYe&    cgmye,
                                       const Heston&   heston,
                                       Complex&        alpha,
                                       Complex&        beta);
    
    /* Methods called by FourierProcessSVJ. Implemented on VolCGMYHeston to
       avoid dereferencing */
    Complex scalelessCumulant(const StFourierProcessLogRtn& process,
                              const StFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;
    Complex scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                              const FwdStFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;

    void validatePop2Object();

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const;
    
    /** calculate a simple expected variance */
    double FutureVariance(double mat) const;

    /** adjust mean vol and reversion rate for risk premium */
    void AdjustRiskPremium(double volRiskPrice, double& mVar, double& mReversion) const;

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const;

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields);
private:
    /** Called after adjustments have been made to fields (eg calibrator) */
    void update();
};
typedef smartPtr<VolCGMYHeston> VolCGMYHestonSP;
DRLIB_END_NAMESPACE
#endif
