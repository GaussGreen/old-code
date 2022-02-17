//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolAJDSuper.hpp
//
//   Description : Superposition of Affine Jump-Diffusions
//
//   Date        : 10 April 2003
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLAJDSUPER_HPP
#define EDR_VOLAJDSUPER_HPP

#include "edginc/VolBaseParam.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** A parameterised CVolBase (does not support BS and DVF views of the
    world for now) */
class MARKET_DLL VolAJDSuper: public VolBaseParam,
                   public virtual Calibrator::IAdjustable
/* virtual public IVolatilityBS,
   virtual public IVolatilityDVF */ {
public:
    /** Interface that needs to be implemented by a class for it 
        to support superposition */
    class MARKET_DLL ISuperposable: public virtual IObject{
    public:
        static CClassConstSP const TYPE;

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
                                            ComplexArray&       betas) const = 0;

        /** Returns the initial values Y_0,..., X^i_0,...
            Naturally, Y_0 = 0. */
        virtual double getInitialValue(int i) const = 0;

        /** Returns the nber of factors. The first factor is understood to be 
            the log spot. */
        virtual int getNbFactors() const = 0;

        /** This is needed because the vol that implements ISuperposable
            cannot rely on its own name to get its market (eg the backbone
            vol surface, the time metrics). This is because the name of the
            superposable vol has got to be different from the name of the 
            actual vol that it is going to be superposed into. */
        virtual void getMarket(const IModel*     model, 
                               const MarketData* market, 
                               const string&     name) = 0;

    private:
        static void load(CClassSP& clazz);
    };
    typedef ISuperposable Superposable;
    typedef smartPtr<Superposable> SuperposableSP;
    typedef array<SuperposableSP, Superposable> SuperposableArray;

private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class MARKET_DLL AJDSuperVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        AJDSuperVolParam(): CVolParam(TYPE){}

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
            REGISTER(AJDSuperVolParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };

    /* Methods called by AJDSuperVolParam. Implemented on VolAJDSuper to avoid dereferencing */ 
    void ComputeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const;

    VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray&   strikes) const;

    /** Returns the nber of factors. The first factor is understood to be 
        the log spot. */
    virtual int getNbFactors() const;

    /** Returns the initial values Y_0,..., X^i_0,...
        Naturally, Y_0 = 0. */
    virtual void getInitialValues(DoubleArray& initVals) const;

    /** Calculates the components alpha and betas that appear in 
        the exponent of the time-t joint Bilateral Laplace
        of X_T = (Y_T,..., X^i_T,...) where Y_T = ln(S_T / F(0, T)) 
        is the dimension-less log spot at time T and 
        the X^i_T's are any additional factors.
        The Bilateral Laplace transform is of the form
            exp(alpha(t, T, u) + sum_i betas^i(t, T, u) * X^i_t) */
    virtual void calcCumulantComponents(const DateTime&     fromDate,  // time t
                                        const DateTime&     toDate,    // time T
                                        const ComplexArray& u,         // frequencies
                                        Complex&            alpha,
                                        ComplexArray&       betas) const;

    // registered fields
    DoubleArray   sqWeights;
    StringArray   vols;      // array of superposable vol types
    StringArray   names;     // array of market names

    // transient fields (won't appear in dd interface)
    DateTime          baseDate;
    DoubleArray       usedSqWeights;
    DoubleArray       usedWeights;
    SuperposableArray supervols;    // array of superposable vols

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VolAJDSuper();

    static IObject* defaultCtor();

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache();

public:
    friend class AJDSuperVolParam;

    static CClassConstSP const TYPE;

    void validatePop2Object();

    virtual void getMarket(const IModel* model, 
                           const MarketData* market);

    /* Methods called by FourierProcessAJDSuper. Implemented on VolAJDSuper to
       avoid dereferencing */
    Complex scalelessCumulant(const StFourierProcessLogRtn& process,
                              const StFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;
    Complex scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                              const FwdStFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;
#if 0
    Complex cumulant(const StFourierProcessIntVar& process,
                     const StFourierProductIntVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;
    Complex cumulant(const FwdStFourierProcessIntVar& process,
                     const FwdStFourierProductIntVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;
    Complex cumulant(const StFourierProcessQuadVar& process,
                     const StFourierProductQuadVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;
    double expectation(const StFourierProcessQuadVar& process,
                       const StFourierProductQuadVar& product, 
                       const DateTime& matDate) const;
#endif

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const;

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const;

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields);
private:    
    /** Called after adjustments have been made to fields */
    void update();
};
typedef smartPtr<VolAJDSuper> VolAJDSuperSP;

DRLIB_END_NAMESPACE
#endif
