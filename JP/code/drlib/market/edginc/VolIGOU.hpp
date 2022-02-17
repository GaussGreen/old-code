//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolIGOU.hpp
//
//   Description : 
//
//   Date        : 07 November 2002
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLIGOU_HPP
#define EDR_VOLIGOU_HPP

#include "edginc/IDynamicsParameter.hpp"
#include "edginc/VolBaseParam.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Range.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/VolAJDSuper.hpp"
#include "edginc/VolOUHelper.hpp"

DRLIB_BEGIN_NAMESPACE

/** A parameterised CVolBase (does not support BS and DVF views of the world for now) */
class MARKET_DLL VolIGOU: public VolBaseParam,
               public virtual IDynamicsParameter,
               public virtual Calibrator::IAdjustable,
               public virtual VolAJDSuper::ISuperposable {
              // virtual public IVolatilityBS,
              // virtual public IVolatilityDVF
private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class MARKET_DLL VolIGOUParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        VolIGOUParam(): CVolParam(TYPE){}

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
            REGISTER(VolIGOUParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };

    /* Methods called by SVJVolParam. Implemented on VolIGOU to avoid dereferencing */ 
    void ComputeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const;

    VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray&   strikes) const;

    /** Calculate the components alpha and beta that appear in the exponent 
        of the time-t joint charcteristic function of ln(S_T / F(0, T)) in the 
        IGOU model.
        The Bilateral Laplace transform is of the form
             exp(alpha(tau, u1, u2) + u1 * Y_t + beta(tau, u1, u2) * V_t)
        where tau = T - t and the complex u1, u2 are the frequencies wrt Y_T and 
        V_T, respectively. */
    void calcCFComponents(double         tau,    // tau = T - t (in years)
                          const Complex& u,
                          Complex&       comp1,
                          Complex&       comp2) const;

    // registered fields
    double initialVol;
    double correlation;
    double volVol;
    double meanVol;
    double meanReversRate;

    /** Helper class for IG OU-type process - does the actual work */
    class MARKET_DLL Helper: public VolOUHelper{
    public:
        static CClassConstSP const TYPE;

        Helper(double meanVol,
               double meanReversRate,
               double volVol,                                  
               double correlation);

    #if 0
        static void mapFromVarsToVols(const double& gamma,
                                      const double& delta,
                                      double&       meanVol,
                                      double&       volVol);
    #endif

    protected:
        /** Computes lambda * \int_0^t k(f(s)) where
            k(x) = delta * x / sqrt(gamma^2 - 2 * x), gamma^2 - 2 * x > 0 and    
            f(s) = c1 + c2 exp(-lambda(t-s)) */
        virtual Complex integral(double         tau,
                                 const Complex& c1,
                                 const Complex& c2);

        /** Returns log E_0 exp(u Z(1)) = nu * x / (alpha - x) for alpha > x*/
        virtual Complex cumulantBDLP(const Complex& u);

    private:
        static Complex atan(Complex x);

        static Complex atanh(Complex x);

        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);

        Helper();

        static IObject* defaultCtor();

        // vars
        double gamma;
        double delta;
    };
    typedef smartPtr<Helper> HelperSP;
    // transient fields (won't appear in dd interface)
    DateTime     baseDate;
    HelperSP     helper;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VolIGOU();

    static IObject* defaultCtor();

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache();

public:
    friend class VolIGOUParam;

    static CClassConstSP const TYPE;

    virtual ~VolIGOU();

    void validatePop2Object();

    /* Methods called by FourierProcessIGOU. Implemented on VolIGOU to
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

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const;
    
    /** calculate a simple expected variance */
    double FutureVariance(double mat) const;

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const;

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields);

    /** Calculates the components alpha and betas that appear in 
        the exponent of the time-t joint Bilateral Laplace
        of X_T = (Y_T,..., X^i_T,...) where Y_T = ln(S_T / F(0, T)) 
        is the 'weighted' dimension-less log spot at time T and 
        the X^i_T's are any additional factors.
        The Bilateral Laplace transform is of the form
            exp(alpha(t, T, u) + sum_i betas^i(t, T, u) * X^i_t) */
    virtual void calcCumulantComponents(const DateTime&     fromDate,  // time t
                                        const DateTime&     toDate,    // time T
                                        const ComplexArray& u,         // frequencies
                                        double              weight,
                                        Complex&            alpha,
                                        ComplexArray&       betas) const;

    /** Returns the initial values Y_0,..., X^i_0,...
        Naturally, Y_0 = 0. */
    virtual double getInitialValue(int i) const;

    /** Returns the nber of factors. The first factor is understood to be 
        the log spot. */
    virtual int getNbFactors() const;
    
    virtual void getMarket(const IModel*     model, 
                           const MarketData* market, 
                           const string&     name);
private:
    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void update();
};
typedef smartPtr<VolIGOU> VolIGOUSP;

DRLIB_END_NAMESPACE
#endif
