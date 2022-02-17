//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolMerton.hpp
//
//   Description : 
//
//   Author      : Oliver Brockhaus
//
//   Date        : 03 April 2003
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLMERTON_HPP
#define EDR_VOLMERTON_HPP

#include "edginc/IDynamicsParameter.hpp"
#include "edginc/VolBaseParam.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Heston.hpp"
#include "edginc/VolRequestMerton.hpp"
#include "edginc/VolProcessed.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////////
// VolMerton
//
//  this file contains also classes
//      - FourierProcessMerton 
//      - VolProcessedMerton
//
///////////////////////////////////////////////////////////////////////////////////
/** A parameterised CVolBase (does not support BS and DVF views of the
    world for now) */
class MARKET_DLL VolMerton: public VolBaseParam,
                 public virtual Calibrator::IAdjustable,
                 public virtual IDynamicsParameter
/* virtual public IVolatilityBS,
   virtual public IVolatilityDVF */ {
private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class MARKET_DLL MertonVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        MertonVolParam(): CVolParam(TYPE){}

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
            REGISTER(MertonVolParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };

    /* Methods called by MertonVolParam. Implemented on VolMerton to avoid dereferencing */ 
    void ComputeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const;

    VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray&   strikes) const;

    /** Calculate the components alpha and beta that appear in the exponent 
        of the time-t joint Laplace transform of X_T = (Y_T, V_T) in the 
        Heston Merton model where Y_T = ln(S_T / F(0, T))
        The Bilateral Laplace transform is of the form
        exp(alpha(tau, u1, u2) + u1 * Y_t + beta(tau, u1, u2) * V_t)
        where tau = T - t and the complex u1, u2 are the frequencies
        wrt Y_T and V_T, respectively. */
    void calcJointLapAlphaBeta(double         tau,    // tau = T - t (in years)
                               const Complex& u1,
                               const Complex& u2,
                               Complex&       alpha,
                               Complex&       beta) const;

    /** Calculate the components alpha and beta that appear in the exponent 
        of the time-t joint Laplace transform of X_T = (V_T, Y_T) in the 
        Heston Merton model where Y_T = Integrated Variance [t,T]
        The Bilateral Laplace transform is of the form
        exp(alpha(tau, u1, u2) + beta(tau, u1, u2) * V_t)
        where tau = T - t and the complex u1, u2 are the frequencies 
        wrt V_T and Y_T, respectively. */
    /*
    void calcIntVarLapAlphaBeta(double         tau,    // tau = T - t (in years)
                                const Complex& u1,     // Variance 
                                const Complex& u2,     // Integrated Variance
                                Complex&       alpha,
                                Complex&       beta) const;
    */

    // registered fields
    double      initialVol;
    double      meanVol;
    double      meanReversRate;
    double      crashRate;            // Poisson rate of jumps. Same as Merton's p313 jumpLambda.
    double      crashSizeMean;        // k=E[dS/S] for one jump as Merton p 313/321 (k=exp(jumpGamma)-1).
    double      crashSizeUncertainty; // Std dev of log(spot) under one jump (Merton's delta p 321).
    double      volRiskPrice;         // coefficient of market price of risk for stochastic vol
    double      beta;                 // holds the beta wrt the "market" jumps

    // transient fields (won't appear in dd interface)
    DateTime      baseDate;
    HestonSP      heston;
    MertonCrashSP mertonCrash;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VolMerton();

    static IObject* defaultCtor();

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache();

public:
    friend class MertonVolParam;
    friend class VolProcessedMerton;

    static CClassConstSP const TYPE;

    void validatePop2Object();

    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

    CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const;

    /* Methods called by FourierProcessMerton. Implemented on VolMerton to
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

    Complex cumulant(const StFourierProcessQuadVar& process,
                     const StFourierProductQuadVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;
    double expectation(const StFourierProcessQuadVar& process,
                       const StFourierProductQuadVar& product, 
                       const DateTime& matDate) const;

    /** calculate a simple expected variance */
    double FutureVariance(double mat) const;

    /** adjust mean vol and reversion rate for risk premium */
    void AdjustRiskPremium(double volRiskPrice, double& mVar, 
                           double& mReversion) const;

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const;

    /** Needed for IAdjustable interface. Returns base date for vol */
    virtual DateTime getBaseDate() const;
    double getBeta() const { return beta; };

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields);

    class MARKET_DLL Request {
    private:
        smartPtr<VolMerton> theVol;
    };

    class MARKET_DLL Processed {
    private:
        smartPtr<VolMerton> theVol;
    };

    friend class Processed;
private:
    /** Called after adjustments have been made to fields */
    void update();
};
typedef smartPtr<VolMerton> VolMertonSP;
typedef smartConstPtr<VolMerton> VolMertonConstSP;

///////////////////////////////////////////////////////////////////////////////////
// VolProcessedMerton
///////////////////////////////////////////////////////////////////////////////////

class CVolBase;
class CVolRequest;

/** Defines what a Black-Scholes processed volatility can do. Can be thought
    of as a Vol Curve */
class MARKET_DLL VolProcessedMerton: public CObject,
                          public virtual IVolProcessed{
public:
    static CClassConstSP const TYPE;

    /** Used to indicate the type of calculation should be performed when
        an array of dates is passed to CalcVar/CalcVol */
    typedef enum _TCalcType{
        fromFirst = 0, // calculate vol/variance from first date
        forward,       // calculate vol/variance between successive dates
        toLast         // calculate vol/variance to last date
    } TCalcType;

    // constructor
    VolProcessedMerton(const VolMertonConstSP& vol);

    /** calculates the trading time between two dates */
    virtual double calcTradingTime( const DateTime &date1, 
                                    const DateTime &date2) const;

    /** identifies the market data name of the volatility */
    virtual string getName() const;

    /** identifies the base date of the volatility */
    virtual DateTime getBaseDate() const;

    /** retrieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric() const;

    /** Calculates variance between 2 dates */
    virtual double CalcVar(const DateTime &date1,
                           const DateTime &date2) const;
    
    /** Calculates variance from 0 to time (in years) */
    virtual double CalcVar( const double time) const;

    /** Calculates variance between a series of dates. If the dateList
        has n dates in it, n-1 variances will be calculated. */
    virtual void CalcVar(const DateTimeArray& dateList,
                         TCalcType            calcType, 
                         CDoubleArray&        vars) const;

    /** Calculates variance beginning at dateFrom. If the dateList
        has n dates in it, n variances will be calculated. */
    virtual void CalcVar(const DateTime&      dateFrom,
                         const DateTimeArray& datesTo,
                         TCalcType            calcType, 
                         CDoubleArray&        vars) const;

    /** Calculates volatility between dates. */
    virtual double CalcVol(const DateTime& date1, 
                           const DateTime& date2) const;
    virtual void CalcVol(const DateTimeArray& dateList, 
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;
    virtual void CalcVol(const DateTime&      dateFrom,
                         const DateTimeArray& datesTo,
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;

    /** Calculates drift factor between dates. */
    virtual double CalcDrift(const DateTime& date1, 
                             const DateTime& date2) const;
    virtual void CalcDrift(const DateTimeArray& dateList, 
                           TCalcType            calcType, 
                           CDoubleArray&        factors) const;
    virtual void CalcDrift(const DateTime&      dateFrom,
                           const DateTimeArray& datesTo,
                           TCalcType            calcType, 
                           CDoubleArray&        factors) const;

    /** Calculates the number of jumps n such that 
        the probability of more than n jumps between dates date1 and date2 
        is less than epsilon. */
    void Quantile(
        const DateTime &date1,
        const DateTime &date2,
        double         epsilon,
        int            maxJumps,
        double*        quantile,
        int*           nbJumps ) const;

    /** Calculates jump factor between dates. */
    virtual double CalcJump(double noise,
                            int    numJumps) const;

    /** Preprocesses fields required for generateCumulatives() */
    void setupCumulatives(
        const DateTime&      dateFrom,
        const DateTimeArray& datesTo,
        const double*        driverSpotVars
    );

    /** Samples from copula of VolProcessedMerton */
    void generateCumulatives(
        const double*       random,
        const double*       randomJump,
        const double*       numJumps,
        double*             path) const;

    /**  */
    double getBeta() const { return volMerton->beta; };
    double getCrashRate() const { return volMerton->crashRate; };
    double getCrashSizeUncertainty() const { return volMerton->crashSizeUncertainty; };

protected:
    VolProcessedMerton(const CClassConstSP& clazz);
private:
    VolProcessedMerton(const VolProcessedMerton &rhs);
    VolProcessedMerton& operator=(const VolProcessedMerton& rhs);

    // transient fields
    VolMertonConstSP volMerton; // $unregistered
    string          name; // $unregistered
    DateTime        baseDate; // $unregistered

    /** optional */
    double          volRiskPrice;         // coefficient of market price of risk for stochastic vol $unregistered

    /** for performance in CalcJump */
    double          lnOnePlusMeanMinusHalfVar; // $unregistered

    /** for cumulatives */
    CDoubleArraySP  volCont; // $unregistered
    CDoubleArraySP  volJump; // $unregistered
    CDoubleArraySP  meanJump; // $unregistered
    CDoubleMatrixSP factor; // $unregistered
    CDoubleMatrixSP proba; // $unregistered
    CIntArraySP     maxNumJumps; // $unregistered
    int             numSteps; // $unregistered
    int             maxMaxNumJumps; // $unregistered
    CDoubleArraySP  driverFwdVol; // $unregistered
    CDoubleArraySP  driverSpotVol; // $unregistered

};
typedef smartConstPtr<VolProcessedMerton> VolProcessedMertonConstSP;
typedef smartPtr<VolProcessedMerton> VolProcessedMertonSP;

DRLIB_END_NAMESPACE
#endif
