//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSVCJ.hpp
//
//   Description : stock diffusion, vol diffusion,  commonStockVolJump
//
//   Date        : 20 Apr 2004
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLSVCJ_HPP
#define EDR_VOLSVCJ_HPP

#include "edginc/VolBaseParam.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/IDynamicsParameter.hpp"
#include "edginc/MarketDataConvert.hpp"
#include <list>

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VolSV)
FORWARD_DECLARE(VolSVJ)
FORWARD_DECLARE(VolSVJJ)


/** A parameterised CVolBase (does not support BS and DVF views of the world for now) */
class MARKET_DLL VolSVCJ: public VolBaseParam,
               public virtual Calibrator::IAdjustable,
               public virtual IDynamicsParameter,
               public virtual MarketDataConvert::IConvert<VolSV>,
               public virtual MarketDataConvert::IConvert<VolSVJ>,
               public virtual MarketDataConvert::IConvert<VolSVJJ>{
public:
    /** Full constructors */
    VolSVCJ(const string& name,
            double initialVol,
            double correlation,
            double volVol,
            double meanVol,
            double meanReversRate,
            double commonCrashRate,
            double commonStockCrashSizeMean,
            double commonStockCrashSizeUncertainty,
            double commonVolCrashSizeMean,
            double stockVolCrashSizeCorrelation,
            bool randInitVol);

    /** Downgrade to VolSV if jump parameters are zero*/
    virtual VolSVSP convert(VolSV* p = 0) const;

    /** Downgrade to VolSVJ if vol jump parameters are zero*/
    virtual VolSVJSP convert(VolSVJ* p = 0) const;

    /** Upgrade to VolSVJJ with zero idiosyncratic jump parameters */
    virtual VolSVJJSP convert(VolSVJJ* p = 0) const;
    
    typedef list<double>::iterator LIST_ITER;
    typedef enum{START_DATE, SAMPLE_DATE, DIFF_DATE, JUMP_DATE} STEP_TYPE;
    typedef list<STEP_TYPE>::iterator STEP_TYPE_ITER;
protected:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class MARKET_DLL SVCJVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        SVCJVolParam(): CVolParam(TYPE){}

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
            REGISTER(SVCJVolParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };

    /* Methods called by SVCJVolParam. Implemented on VolSVCJ to avoid dereferencing */ 
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
    
	//GAD 19/02/2006

    //Laplace transform of the quadratic variation
    void calcJointLapAlphaBetaCommon(double         tau,    //tau = T - t (in years)
                                     const Complex& u1, //quad var
                                     const Complex& u2, //spot var
                                     Complex&       alpha,
                                     Complex&       beta) const;

	//end of GAD 19/02/2006

    // registered fields
    double initialVol;
    double correlation;
    double volVol;
    double meanVol;
    double meanReversRate;
    /** Commom jump, The common vol jump size is a constant shift (binomial to be considered)
       and conditional on that jump the common stock jump size is lognormally ditributed. */
    double commonCrashRate;                 // Common jump rate
    double commonStockCrashSizeMean;        // Common stock jump size (conditional/unconditional) mean 
                                            // @ zero vol jump size mean and/or zero crash correlation
    double commonStockCrashSizeUncertainty; // Common stock jump size conditional std dev
    double commonVolCrashSizeMean;          // Vol jump size (unconditional) mean
    double stockVolCrashSizeCorrelation;    // 'Correlation' parameter between stock and vol common jump sizes
                                            // Not correlation exactly, but drives it. Note: lives in the real line
    bool   randInitVol;                     // if true, use a random initial volatility (gamma distribution)

    // transient fields (won't appear in dd interface)
    DateTime      baseDate;

    /* Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VolSVCJ(CClassConstSP clazz = TYPE);

    static IObject* defaultCtor(){
        return new VolSVCJ();
    }

    /*** Build the parameterised vol and cache any values **/
    virtual void getMarket(const IModel* model, const MarketData* market);

public:
    friend class SVCJVolParam;
    friend class FD2DSVCJ;
    friend class CalcAlpha; // local helper functions
    friend class CommonAlphaHelper; // local helper functions for the computation pof alpha common for integrated variance
    friend class FD2DSolverJumps;

    static CClassConstSP const TYPE;

    //specific method implemented in order to be used in MCPathConfigSVJ for volRequest (vol + trading time) for CCYP options
    //Taken form SRM Should be careful that it cannot be called for other products

    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

    void validatePop2Object();

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
        return new SVCJVolParam();
    }

    /* Methods called by FourierProcessSVCJ. Implemented on VolSVCJ to
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
    
	//GAD 19/02/2006

    Complex cumulant(const StFourierProcessQuadVar&    process,
                     const StFourierProductQuadVar&    product, 
                     const Complex&	                   z, 
                     const DateTime&				   matDate) const;
    
    Complex cumulant(const FwdStFourierProcessQuadVar& process,
                     const FwdStFourierProductQuadVar& product, 
                     const Complex&                    z, 
                     const DateTime&                   matDate) const;

    double expectation(const StFourierProcessQuadVar& process,
                       const StFourierProductQuadVar& product, 
                       const DateTime& matDate) const;

    double expectation(const FwdStFourierProcessQuadVar& process,
                       const FwdStFourierProductQuadVar& product, 
                       const DateTime& matDate) const;

	//end of GAD 19/02/2006

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const;

    enum EnumParam{INITIAL_VOL = 0, CORRELATION, VOL_VOL, MEAN_VOL, MEAN_REVERS_RATE ,
                    COMMON_CRASH_RATE, COMMON_STOCK_CRASH_SIZE_MEAN, COMMON_STOCK_CRASH_SIZE_UNCERTAINTY, 
                    COMMON_VOL_CRASH_SIZE_MEAN, STOCK_VOL_CRASH_SIZE_CORRELATION,
                    STATIONARY_VARIANCE};
    double getSVCJParam(const EnumParam param) const;
    
    bool isRandomInitialVolatility(){
        return randInitVol;
    }

    // *********** MC support for simulating arrival times*********************/
    /** compute shape param for gamma distribution */
    double gammaShape(){
        return (initialVol*initialVol*meanReversRate*2.0/volVol/volVol);
    }

    /** compute jump contributions for var and log spots */
    void addJumps(list<double>& tradYrs,
                    const double* jumpTimes,
                    int           numJumps,
                    BoolArray     jumpParticipation,
                    const double* jumpDeviatesVar,
                    const double* jumpDeviatesSpot,
                    list<STEP_TYPE>& stepTypes,
                    list<double>& instVars,
                    list<double>& integratedVars,
                    list<double>& logSpots,
                    vector<double>& varJumps,
                    vector<STEP_TYPE_ITER>& jumpTypesAdded,
                    vector<vector<LIST_ITER> >&jumpsAdded,
                    double          randomInitVol) const;

    /** Given an array of gaussian deviates, populate the arrays 'instVars' 
        and 'integratedVars' with simulated values using an Euler
        discretization scheme */
    virtual void generateVarPathsEuler(LIST_ITER       tradYears,
                                       int             nbSteps,
                                       const double*   deviates,
                                       LIST_ITER       instVars,
                                       LIST_ITER       integratedVars) const;

    /** Given an array of gaussian deviates, populate the arrays 'instVars' 
        and 'integratedVars' with simulated values using a variable transform
        scheme together with an Euler discretization scheme */
    virtual void generateVarPathsTransEuler(LIST_ITER       tradYears,
                                            int             nbSteps,
                                            const double*   deviates,
                                            LIST_ITER       instVars,
                                            LIST_ITER       integratedVars) const;

    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
        values using an Euler discretization scheme for the spot */
    virtual void generateSpotPathsEuler(LIST_ITER          tradYears,
                                        int                nbSteps,
                                        STEP_TYPE_ITER     stepTypes,
                                        const double*      spotDeviates,
                                        const double*      varDeviates,
                                        LIST_ITER          logSpots,
                                        LIST_ITER          instVars,
                                        LIST_ITER          integratedVars) const;

    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
        values. The log spots are simulated exactly */
    virtual void generateSpotPathsExact(LIST_ITER          tradYears,
                                        int                nbSteps,
                                        STEP_TYPE_ITER     stepTypes,
                                        const double*      spotDeviates,
                                        vector<double>&    varJumps,
                                        LIST_ITER          logSpots,
                                        LIST_ITER          instVars,
                                        LIST_ITER          integratedVars) const;
};

typedef smartPtr<VolSVCJ> VolSVCJSP;
typedef smartConstPtr<VolSVCJ> VolSVCJConstSP;
typedef array<VolSVCJSP, VolSVCJ>  VolSVCJArray;



/************************/
/*** SVCJ + LV scheme ***/
/************************/

class MARKET_DLL VolSVCJLV: public VolSVCJ,
                 public virtual Calibrator::IAdjustable {
private:
    /** return a copy of LV expiries */
    static ExpiryArraySP getLocVolExpiries(const IObject* obj);

    /** compute local volatility function for SVCJ + LV scheme */
    double calcLocVol(const double logSpot, const double paramA, const double paramB, const double paramC) const;

    // SVCJ + LV scheme
    ExpiryArraySP locVolExpiries;  // expiries for local volatility adjustment in SVCJ + LV scheme
    DoubleArraySP locVolParamsA;   // parameters A for local volatility adjustment in SVCJ + LV scheme
    DoubleArraySP locVolParamsB;   // parameters B for local volatility adjustment in SVCJ + LV scheme
    DoubleArraySP locVolParamsC;   // parameters C for local volatility adjustment in SVCJ + LV scheme

    /* Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VolSVCJLV();

    static IObject* defaultCtor(){
        return new VolSVCJLV();
    }

public:
    static CClassConstSP const TYPE;

    void validatePop2Object();

    /** extend local volatility parameters arrays
        in order to get array with size equal to the number of steps in MC */
    void extendLocVolParams(LIST_ITER tradYears,
                            int nbSteps,
                            DoubleArraySP locVolParamsAExt,
                            DoubleArraySP locVolParamsBExt,
                            DoubleArraySP locVolParamsCExt) const;

    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
        values using an Euler discretization scheme for the spot */
    virtual void generateSpotPathsEuler(LIST_ITER           tradYears,
                                        int                nbSteps,
                                        STEP_TYPE_ITER     stepTypes,
                                        const double*      spotDeviates,
                                        const double*      varDeviates,
                                        LIST_ITER          logSpots,
                                        LIST_ITER          instVars,
                                        LIST_ITER          integratedVars) const;

    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
        values. The log spots are simulated exactly */
    virtual void generateSpotPathsExact(LIST_ITER           tradYears,
                                        int                nbSteps,
                                        STEP_TYPE_ITER     stepTypes,
                                        const double*      spotDeviates,
                                        vector<double>&    varJumps,
                                        LIST_ITER          logSpots,
                                        LIST_ITER          instVars,
                                        LIST_ITER          integratedVars) const;
};

typedef smartPtr<VolSVCJLV> VolSVCJLVSP;
typedef smartConstPtr<VolSVCJLV> VolSVCJLVConstSP;

typedef MarketDataConvert::IConvert<VolSVCJ> VolSVCJConvert;

DRLIB_END_NAMESPACE
#endif
