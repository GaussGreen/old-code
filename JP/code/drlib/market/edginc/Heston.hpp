//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Heston.hpp
//
//   Description : Building blocks for Heston-type vol objects 
//                 and their Laplaces
//
//   Date        : 27 Nov 2002
//
//
//----------------------------------------------------------------------------
#ifndef EDR_HESTON_HPP
#define EDR_HESTON_HPP

#include "edginc/Object.hpp"
#include "edginc/Class.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Complex.hpp"
#include "edginc/Range.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL MertonCrash: public CObject{
public:
    static CClassConstSP const TYPE;

    /** DefaultVal or suggested parameter values */
    struct MARKET_DLL DefaultVal{
        static const double crashRate;
        static const double crashSizeMean;
        static const double crashSizeUncertainty;
    };

    /** Ranges of definition */
    struct MARKET_DLL RangeDef{
        static const Range crashRate;
        static const Range crashSizeMean;
        static const Range crashSizeUncertainty;
    };

    MertonCrash(double crashRate,
                double crashSizeMean,
                double crashSizeUncertainty);

    virtual void validatePop2Object();

    /** Calculates the contribution of the jump component to the component 
        alpha that appears in the exponent of the time-t Bilateral Laplace 
        of Y_T in the Merton model where Y_T = ln(S_T / F(0, T)) is the 
        dimension-less log-return at time T.
        The Bilateral Laplace transform is of the form
            exp(alpha(tau, u1) + u1 * Y_t)
        where tau = T - t and the complex u1 is the frequencies wrt Y_T. */
    Complex calcAlpha(double         tau,    // tau = T - t (in years)
                      const Complex& u1) const;

    /** Calculates the cumulant of the quadratic variation Sum logJump^2
        of the jump process in the Merton model */
    Complex calcQuadVarCumulant(double         tau,    // tau = T - t (in years)
                                const Complex& z) const;

    /** Calculates the annualized expected quadratic variation of the
        jump process in Merton */
    double calcAnnualQuadVar(double tau) const;

protected:
    MertonCrash(const CClassConstSP& clazz);

private:
    MertonCrash(const MertonCrash& rhs);
    MertonCrash& operator=(const MertonCrash& rhs);

    MertonCrash();

    static IObject* defaultCtor();
    static void load(CClassSP& clazz);

    // registered fields
    double crashRate;
    double crashSizeMean;
    double crashSizeUncertainty;
    // transient field
    bool   zeroCrashRate;
    double crashGamma;
    double sqCrashSizeUncertainty;
};
typedef smartPtr<MertonCrash> MertonCrashSP;

class Heston;
typedef smartPtr<Heston> HestonSP;
/** Repository class for the parameters of the Hestom model and its associated
    Laplace transform(s). Nested in the Heston class are also repository classes 
    for the parameters of the jump components of the vol and the stock when
    the Heston model is augmented with jumps in the vol. The latter classes are 
    nested in the Heston class in order to emphasize the fact that the Laplaces 
    associated with the jump components are only applicable in the Heston model */
class MARKET_DLL Heston: public CObject{
public:
    class VolCrash;
    friend class VolCrash;
    class CommonCrash;
    friend class CommonCrash;

    static CClassConstSP const TYPE;

    /** DefaultVal or suggested parameter values */
    struct MARKET_DLL DefaultVal{
        static const double initialVol;
        static const double meanVol;
        static const double meanReversRate;
        static const double volVol;
        static const double correlation;
    };

    /** Ranges of definition */
    struct MARKET_DLL RangeDef{
        static const Range initialVol;
        static const Range meanVol;
        static const Range meanReversRate;
        static const Range volVol;
        static const Range correlation;
    };

    Heston(double initialVol,
           double meanVol,
           double meanReversRate,
           double volVol,
           double correlation);

    virtual void validatePop2Object();

    /** Calculates the components alpha and beta that appear in the exponent 
        of the time-t joint Bilateral Laplace of X_T = (Y_T, V_T, I_T) in the 
        Heston model where Y_T = ln(S_T / F(0, T)) is the dimension-less
        log spot at time T, V_T is the instantaneous variance at time T 
        and I_T is the instantaneous variance from 0 to time T.
        The Bilateral Laplace transform is of the form
            exp(alpha(tau, u1, u2, u3) 
                + u1 * Y_t 
                + beta(tau, u1, u2, u3) * V_t 
                + u3 * I_t)
        where tau = T - t and the complex numbers u1, u2, u3 are the frequencies 
        wrt Y_T, V_T and I_T, respectively. */
    void calcAlphaBeta(double         tau,    // tau = T - t (in years)
                       const Complex& u1,
                       const Complex& u2,
                       const Complex& u3,
                       Complex&       alpha,
                       Complex&       beta) const;

    /** Calculates the annualized expected quadratic variation in Heston,
        i.e. the annualized expected integrated variance */
    double calcAnnualQuadVar(double tau) const;

    /** Calculates the forward starting annualized expected quadratic variation in Heston,
        i.e. the annualized expected integrated variance, tau = T - t (in years) */
    double calcAnnualQuadVar(double t,
                                     double tau) const;

    /** Repository class for the parameters of the vol-only jump component of the 
        vol when the Heston model is augmented with jumps in the vol. The jump in the 
        vols are assumed to be exponentially distributed here (see Duffie et al 1999) */
    class MARKET_DLL VolCrash: public CObject{
    public:
        static CClassConstSP const TYPE;

        /** DefaultVal or suggested parameter values */
        struct MARKET_DLL DefaultVal{
            static const double crashRate;
            static const double crashSizeMean;
        };

        /** Ranges of definition */
        struct MARKET_DLL RangeDef{
            static const Range crashRate;
            static const Range crashSizeMean;
        };

        VolCrash(double          crashRate,
                 double          crashSizeMean,
                 const HestonSP& heston);

        virtual void validatePop2Object();

        /** Calculates the contribution of the vol jump component to the component 
            alpha that appears in the exponent of the time-t joint Bilateral Laplace 
            of X_T = (Y_T, V_T, I_T) where Y_T = ln(S_T / F(0, T)) is the dimension-less
            log spot at time T, V_T is the instantaneous variance at time T
            and I_T is the instantaneous variance from 0 to time T
            when (i) the diffusion component of V is given by the Heston model and 
            (ii) V is augmented with a vol jump component which is exponentially 
            distributed.
            The Bilateral Laplace transform is of the form
                exp(alpha(tau, u1, u2, u3) 
                    + u1 * Y_t 
                    + beta(tau, u1, u2, u3) * V_t 
                    + u3 * I_t)
            where tau = T - t and the complex numbers u1, u2, u3 are the frequencies 
            wrt Y_T, V_T and I_T, respectively. */
        Complex calcAlpha(double         tau,    // tau = T - t (in years)
                          const Complex& u1,
                          const Complex& u2,
                          const Complex& u3) const;

    protected:
        VolCrash(const CClassConstSP& clazz);

    private:
        VolCrash(const VolCrash& rhs);
        VolCrash& operator=(const VolCrash& rhs);

        VolCrash();

        static IObject* defaultCtor();
        static void load(CClassSP& clazz);

        // registered fields
        double crashRate;
        double crashSizeMean;
        double crashSizeUncertainty; // $unregistered
        HestonSP heston;
        // transient field
        double meanReversRate;
        double volVol;
        double correlation;
    };
    typedef smartPtr<VolCrash> VolCrashSP;

    /** Repository class for the parameters of the jump components that are common
        to the stock and the vol when the Heston model is augmented with jumps in 
        the vol. The jump in the vols are assumed to be exponentially distributed 
        here whereas the jump in the stock are assumed to be normally distributed
        conditional on the size of the vol jump (see Duffie et al 1999) */
    class MARKET_DLL CommonCrash: public CObject{
    public:
        static CClassConstSP const TYPE;

        /** DefaultVal or suggested parameter values */
        struct MARKET_DLL DefaultVal{
            static const double crashRate;
            static const double stockCrashSizeMean;
            static const double stockCrashSizeUncertainty;
            static const double volCrashSizeMean;
            static const double crashCorrelation;
        };

        /** Ranges of definition */
        struct MARKET_DLL RangeDef{
            static const Range crashRate;
            static const Range stockCrashSizeMean;
            static const Range stockCrashSizeUncertainty;
            static const Range volCrashSizeMean;
            static const Range crashCorrelation;
        };

        CommonCrash(double crashRate,
                    double stockCrashSizeMean,
                    double stockCrashSizeUncertainty,
                    double volCrashSizeMean,
                    double crashCorrelation,
                    const HestonSP& heston);

        virtual void validatePop2Object();

        /** Calculates the contribution of the jump components that are common to 
            the stock and to the vol, to the component alpha that appears in 
            the exponent of the time-t joint Bilateral Laplace 
            of X_T = (Y_T, V_T, I_T) where Y_T = ln(S_T / F(0, T)) is the dimension-less
            log spot at time T, V_T is the instantaneous variance at time T
            and I_T is the instantaneous variance from 0 to time T
            when (i) the diffusion component of V is given by the Heston model 
            and (ii) V is augmented with a vol jump 
            component which is exponentially distributed and Y is augmented with a 
            simultaneous stock jump component which is normally distributed conditional 
            on the size of the vol jump.
            The Bilateral Laplace transform is of the form
                exp(alpha(tau, u1, u2, u3) 
                    + u1 * Y_t 
                    + beta(tau, u1, u2, u3) * V_t 
                    + u3 * I_t)
            where tau = T - t and the complex numbers u1, u2, u3 are the frequencies 
            wrt Y_T, V_T and I_T, respectively. */
        Complex calcAlpha(double         tau,    // tau = T - t (in years)
                          const Complex& u1,
                          const Complex& u2,
                          const Complex& u3) const;

    protected:
        CommonCrash(const CClassConstSP& clazz);

    private:
        CommonCrash(const CommonCrash& rhs);
        CommonCrash& operator=(const CommonCrash& rhs);

        CommonCrash();

        static IObject* defaultCtor();
        static void load(CClassSP& clazz);

        // registered fields
        double crashRate;                 // Jump rate
        double stockCrashSizeMean;        // Stock jump size (conditional/unconditional) mean 
                                          // @ zero vol jump size mean and/or zero crash correlation
        double stockCrashSizeUncertainty; // Stock jump size conditional std dev
        double volCrashSizeMean;          // Vol jump size (unconditional) mean
        double crashCorrelation;          // 'Correlation' parameter between stock and vol common jump sizes
                                          // Not correlation exactly, but drives it. Note: lives in the real line
        HestonSP heston;
        // transient field
        double meanReversRate;
        double volVol;
        double correlation;
    };
    typedef smartPtr<CommonCrash> CommonCrashSP;

protected:
    Heston(const CClassConstSP& clazz);

private:
    Heston(const Heston& rhs);
    Heston& operator=(const Heston& rhs);

    Heston();

    static IObject* defaultCtor();
    static void load(CClassSP& clazz);

    // registered fields
    double initialVol;
    double meanVol;
    double meanReversRate;
    double volVol;
    double correlation;
    // transient field
    bool   zeroVolVolZeroMeanReversRate;
    bool   zeroVolVolNonZeroMeanReversRate;
    double meanReversRateSqMeanVol;
    double correlationVolVol;
    double sqVolVol;
};

DRLIB_END_NAMESPACE
#endif
