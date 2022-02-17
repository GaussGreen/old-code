//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierEngine.hpp
//
//   Description : Fourier Engine Integrator
//
//   Date        : February 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FOURIERENGINE_HPP
#define EDR_FOURIERENGINE_HPP
#include "edginc/Class.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Asset.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/Model.hpp"
#include "edginc/Integrator.hpp"
#include "edginc/Fourier.hpp"

DRLIB_BEGIN_NAMESPACE

class InstrumentSettlement;

/** Fourier Engine */
class FOURIER_DLL FourierEngine : public CModel {
public:
    static CClassConstSP const TYPE;
    friend class FourierEngineHelper;
    friend class FourierProduct;
    friend class PDFFourier;

    struct FOURIER_DLL IntegratorType{
        enum {
            INTEGRAL = 0,
            FFT,
            NB_ENUMS
        };
        enum {
            defaultIndex = INTEGRAL
        };
    };

    /** Instrument Specific Algorithm Parameters */
    class FOURIER_DLL ISAP: virtual public IObject{
    public:
        static CClassConstSP const TYPE;
    private:
        static void load(CClassSP& clazz);
    };
   
    typedef smartPtr<ISAP> ISAPSP;
    typedef smartConstPtr<ISAP> ISAPConstSP;

    /** interface that the instrument must implement */
    class FOURIER_DLL IIntoProduct: virtual public CModel::IModelIntoProduct{
    public:
        static CClassConstSP const TYPE;

        /** Creates an instance of an FourierProduct */
        virtual FourierProduct* createProduct(const FourierEngine * model) const = 0;
    };
    
    /** Implementation of method in IModel */
    virtual void Price(CInstrument*  instrument,
                       Control*      control, 
                       CResults*     results);

    /** Override default createMDF in order to set the MDF 
     * returned by the process */
    virtual MarketDataFetcherSP createMDF() const;

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns FourierProcess::wantsRiskMapping() from the process.  See
     * IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    const FourierProcess& getProcess() const;

    const FourierEngine::ISAP& getISAP() const;

    void validatePop2Object();

    FourierEngine(const FourierProcessSP&  p,
                  const string&            integratorType,
                  const Integrator1DSP&    integrator,
                  const FFTIntegrator1DSP& fftIntegrator);

    FourierEngine(const FourierProcessSP&  p);

protected:
    FourierEngine(CClassConstSP clazz);

private:
    FourierEngine();

    ISAPSP            isap;
    FourierProcessSP  process;
    Integrator1DSP    integrator;
    FFTIntegrator1DSP fftIntegrator;    // this is weak. FFTIntegrator1D really should be derived from Integrator1D
    string            integratorType;   // so only one integrator field would be required.

    // transient
    int              integratorTypeIndex;
};


/** Fourier Product Interface that needs to be implemented by products 
    that support (default) 1D integrators */
class FOURIER_DLL FourierProductIntegrator1D {
public:
    friend class FourierProduct;

    typedef double Integral;
    typedef CDoubleArray IntegralArray;

protected:
    /** Returns an array of integrands that will be used by the Integrator */
    virtual Function1DDoubleArrayConstSP Integrand(const FourierEngine* model,
                                                   const Integrator1D*  integrator) = 0;

    /** Does post price processing and returns results */
    virtual void postResults(const FourierEngine* model,
                             const Integrator1D*  integrator,
                             const FourierProductIntegrator1D::IntegralArray& integrals,
                             CControl*            control,
                             CResults*            results) = 0;

    /** Gives error reporting when failing to price the i-th payoff */
    virtual string errorHandling(int                  iIntegrand,
                                 const Integrator1D*  integrator) const;
};

/** Fourier Product Interface that needs to be implemented by products 
    that support FFT 1D integrators */
class FOURIER_DLL FourierProductFFTIntegrator1D {
public:
    friend class FourierProduct;

    typedef FFTIntegrator1D::IntegralConstSP Integral;
    typedef vector<FFTIntegrator1D::IntegralConstSP> IntegralArray;

protected:
    /** Returns an array of integrands that will be used by the Integrator */
    virtual Function1DComplexArrayConstSP Integrand(const FourierEngine*    model,
                                                    const FFTIntegrator1D*  integrator) = 0;

    /** Does post price processing and returns results */
    virtual void postResults(const FourierEngine*    model,
                             const FFTIntegrator1D*  integrator,
                             const FourierProductFFTIntegrator1D::IntegralArray& integrals,
                             CControl*               control,
                             CResults*               results) = 0;

    /** Gives error reporting when failing to price the i-th payoff */
    virtual string errorHandling(int                    iIntegrand,
                                 const FFTIntegrator1D* integrator) const;
};

/** Fourier Product */
class FOURIER_DLL FourierProduct {
public:
    friend class FourierEngine;

    // Returns the maturity date
    const DateTime& getToday() const {return today;}

    virtual ~FourierProduct(){};

    static double intersectRanges(double l1,
                                  double u1,
                                  double l2,
                                  double u2,
                                  double p);

    // loop on integration
    virtual void price(const FourierEngine* model,
                       Control*             control, 
                       Results*             results);

    const IMultiFactors& getMultiAsset() const;

    virtual void validate(FourierProcess* process) const;

protected:

    void recordFwdAtMat(Control*             control,
                        const DateTime&      matDate,
                        Results*             results);

    /** 'Full' constructor for single factor payoffs */
    FourierProduct(const CAsset*               asset,        // single factor
                   const DateTime&             today,        // value date
                   const YieldCurve*           discount,     // for pv'ing
                   const InstrumentSettlement* instSettle);  // used for pay date

    /** 'Full' constructor for N factor payoffs */
    FourierProduct(const IMultiFactors*        mFactors,     // N factors
                   const DateTime&             today,        // value date
                   const YieldCurve*           discount,     // for pv'ing
                   const InstrumentSettlement* instSettle);  // used for pay date

private:

    // turns Asset into IMultiFactors
    static IMultiFactorsConstSP convertToMultiFactors(const Asset* asset);

protected:
    IMultiFactorsConstSP    mAsset;     // This is the FourierProduct's view of assets
    DateTime                today;      // value date
    const YieldCurve*       discount;   // discount curve
    const InstrumentSettlement* instSettle;   // instrument settlement
    int                     NbAssets;   // number of assets (not factors)
};
typedef smartPtr<FourierEngine> FourierEngineSP;


class FOURIER_DLL EmptyISAP: public CObject, 
                 public FourierEngine::ISAP {
public:
    static CClassConstSP const TYPE;
    friend class EmptyISAPHelper;

private:
    EmptyISAP(): CObject(TYPE){}
};

typedef smartPtr<EmptyISAP> EmptyISAPSP;
typedef smartConstPtr<EmptyISAP> EmptyISAPConstSP;

DRLIB_END_NAMESPACE

#endif
