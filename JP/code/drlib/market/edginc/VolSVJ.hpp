//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSVJ.hpp
//
//   Description : 
//
//   Date        : 20 May 2002
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLSVJ_HPP
#define EDR_VOLSVJ_HPP


#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Heston.hpp"
#include "edginc/VolAJDSuper.hpp"
#include "edginc/Enum2StringListHelper.hpp"
#include "edginc/Random.hpp"
#include "edginc/VolSV.hpp"
#include "edginc/IDynamicsParameter.hpp"
#include "edginc/MarketDataConvert.hpp"


DRLIB_BEGIN_NAMESPACE

//Quanto structure used for CCYProtected products
class QuantoParamSVJ{
public:
    bool            isQuanto;
    double          eqFXCorr;
    DoubleArray     fxSqrtVar;

    QuantoParamSVJ(bool quanto,double FXcorr,DoubleArray fxVar):isQuanto(quanto),eqFXCorr(FXcorr),fxSqrtVar(fxVar){};

private:
        QuantoParamSVJ();
};
typedef refCountPtr<QuantoParamSVJ> QuantoParamSVJSP;
typedef array<QuantoParamSVJSP, QuantoParamSVJ> QuantoParamSVJArray;

FORWARD_DECLARE(VolSV)
FORWARD_DECLARE(VolSVJJ)
FORWARD_DECLARE(VolSVCJ)

/** A parameterised CVolBase (does not support BS and DVF views of the
    world for now) */
class VolSVJ;
typedef smartPtr<VolSVJ> VolSVJSP;
class VolSVJ: public VolBaseParam,
              public virtual Calibrator::IAdjustable,
              public virtual VolAJDSuper::ISuperposable,
              public virtual IDynamicsParameter,
              public virtual MarketDataConvert::IConvert<VolSV>,
              public virtual MarketDataConvert::IConvert<VolSVJJ>,
              public virtual MarketDataConvert::IConvert<VolSVCJ>
/* virtual public IVolatilityBS,
   virtual public IVolatilityDVF */ {
public:
    /** Full constructor */
    MARKET_DLL VolSVJ(const string& name,
                      double initialVol,
                      double meanVol,
                      double meanReversRate,
                      double volVol,
                      double correlation,
                      double crashRate,
                      double crashSizeMean,
                      double crashSizeUncertainty,
                      double volRiskPrice);
    
    /** Downgrade to VolSV if jump parameters are zero*/
    virtual VolSVSP convert(VolSV* p = 0) const;

    /** Upgrade to VolSVJJ with zero vol jump parameters */
    virtual VolSVJJSP convert(VolSVJJ* p = 0) const;

    /** Upgrade to VolSVCJ with zero vol jump parameters */
    virtual VolSVCJSP convert(VolSVCJ* p = 0) const;

private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class MARKET_DLL SVJVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        SVJVolParam(): CVolParam(TYPE){}

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
            REGISTER(SVJVolParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };

    /* Methods called by SVJVolParam. Implemented on VolSVJ to avoid dereferencing */ 
    MARKET_DLL void ComputeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const;

    MARKET_DLL VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray&   strikes) const;

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
    MARKET_DLL void calcJointLapAlphaBeta(double         tau,    // tau = T - t (in years)
                               const Complex& u1,     // log spot
                               const Complex& u2,     // instantaneous variance
                               const Complex& u3,     // integrated variance
                               Complex&       alpha,
                               Complex&       beta) const;

    // registered fields
    double initialVol;
    double meanVol;
    double meanReversRate;
    double volVol;
    double correlation;
    double crashRate;            // Poisson rate of jumps. Same as Merton's p313 jumpLambda.
    double crashSizeMean;        // k=E[dS/S] for one jump as Merton p 313/321 (k=exp(jumpGamma)-1).
    double crashSizeUncertainty; // Std dev of log(spot) under one jump (Merton's delta p 321).
    double volRiskPrice;         // coefficient of market price of risk for stochastic vol

    // transient fields (won't appear in dd interface)
    DateTime      baseDate;
    HestonSP      heston;
    MertonCrashSP mertonCrash;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VolSVJ();

    static IObject* defaultCtor();

protected:
    /*** Build the parameterised vol and cache any values **/
    MARKET_DLL void buildCache();

public:
    friend class SVJVolParam;

    MARKET_DLL static CClassConstSP const TYPE;

    //specific method implemented in order to be used in MCPathConfigSVJ for volRequest (vol + trading time) for CCYP options
    //Taken form SRM Should be careful that it cannot be called for other products
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

    struct MARKET_DLL MCParams{
        // types of spot discretization schemes
        struct MARKET_DLL SpotDiscreteSchemeType{
            enum {
                EXACT = 0,
                EULER,
                NB_ENUMS
            };
            enum {
                defaultIndex = EXACT
            };
        };
        typedef SpotDiscreteSchemeType SpotDSType;
        typedef Enum2StringListHelper<SpotDSType> SpotDSTypeHelper;

        // types of var discretization schemes
        struct MARKET_DLL VarDiscreteSchemeType{
            enum {
                EULER = 0,
                VAR_TRANSFORM_EULER,
                NB_ENUMS
            };
            enum {
                defaultIndex = EULER
            };
        };
        typedef VarDiscreteSchemeType VarDSType;
        typedef Enum2StringListHelper<VarDSType> VarDSTypeHelper;
    };

    MARKET_DLL void validatePop2Object();

    /** method that builds a CVolParam. */
    MARKET_DLL virtual CVolParam* createVolParam() const;

    /* Methods called by FourierProcessSVJ. Implemented on VolSVJ to
       avoid dereferencing */
    MARKET_DLL Complex scalelessCumulant(const StFourierProcessLogRtn& process,
                              const StFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;
    MARKET_DLL Complex scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                              const FwdStFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;
	/**ARNAUD  Expected Quad Var**/ 
    MARKET_DLL Complex cumulant(const FwdStFourierProcessExpQuadVar& process,              
                     const FwdStFourierProductExpQuadVar& product, 
                     const Complex& z, 
                     const DateTime& matDate) const;

    MARKET_DLL Complex cumulant(const StFourierProcessIntVar& process,
                     const StFourierProductIntVar& product, 
                     const Complex& z, 
                     const DateTime& matDate) const;
    MARKET_DLL Complex cumulant(const FwdStFourierProcessIntVar& process,
                     const FwdStFourierProductIntVar& product, 
                     const Complex& z, 
                     const DateTime& matDate) const;
    MARKET_DLL Complex cumulant(const StFourierProcessQuadVar& process,
                     const StFourierProductQuadVar& product, 
                     const Complex& z, 
                     const DateTime& matDate) const;
    MARKET_DLL Complex cumulant(const FwdStFourierProcessQuadVar& process,
                     const FwdStFourierProductQuadVar& product, 
                     const Complex& z, 
                     const DateTime& matDate) const;
    MARKET_DLL double expectation(const StFourierProcessQuadVar& process,
                       const StFourierProductQuadVar& product, 
                       const DateTime& matDate) const;
    MARKET_DLL double expectation(const FwdStFourierProcessQuadVar& process,
                       const FwdStFourierProductQuadVar& product, 
                       const DateTime& matDate) const;

    /** Needed for IAdjustable interface. Returns market data name for vol */
    MARKET_DLL virtual string getName() const;

    /** Calculates the components alpha and betas that appear in 
        the exponent of the time-t joint Bilateral Laplace
        of X_T = (Y_T,..., X^i_T,...) where Y_T = ln(S_T / F(0, T)) 
        is the 'weighted' dimension-less log spot at time T and 
        the X^i_T's are any additional factors.
        The Bilateral Laplace transform is of the form
            exp(alpha(t, T, u) + sum_i betas^i(t, T, u) * X^i_t) */
    MARKET_DLL virtual void calcCumulantComponents(const DateTime&     fromDate,  // time t
                                        const DateTime&     toDate,    // time T
                                        const ComplexArray& u,         // frequencies
                                        double              weight,
                                        Complex&            alpha,
                                        ComplexArray&       betas) const;

    /** Returns the initial values Y_0,..., X^i_0,...
        Naturally, Y_0 = 0. */
    MARKET_DLL virtual double getInitialValue(int i) const;

    /** Price a PutCliquet*/
    MARKET_DLL virtual double calcPricePutCliquet(double				strike,
                                                  double				coupon,
									   DateTime				valueDate,
									   DateTimeArray		maturities,
									   DateTimeArray		qMaturities,
									   YieldCurveWrapper	discount,
									   CAssetWrapper		asset) const;

    /** Returns the nber of factors. The first factor is understood to be 
        the log spot. */
    MARKET_DLL virtual int getNbFactors() const;
    
    MARKET_DLL virtual void getMarket(const IModel*     model, 
                           const MarketData* market, 
                           const string&     name);


    /** Called after adjustments have been made to fields (eg calibrator) */
    MARKET_DLL virtual void fieldsUpdated(const CFieldArray& fields);

    /** Interface that must be implemented by the many possible
        path generators of VolSVJ */
    class MARKET_DLL IPathGen{
    public:
        virtual void generatePath(
              const DoubleArray&         tradYears,
              const IntArray&            spotOffsets,
              const IMCRandNormal::Path& spotDeviates,
              const IMCRandNormal::Path& varDeviates,
              MCRandPoisson&             jumpNbDeviates,
              MCRandNormal&              jumpSizeDeviates,
              QuantoParamSVJSP           quanto,
              DoubleArray&               logSpots,
              DoubleArray&               instVars,
              DoubleArray&               quadVars) const = 0;
        virtual ~IPathGen(){}
    };
    typedef refCountPtr<IPathGen> IPathGenSP;
    typedef vector<IPathGenSP> IPathGenArray;

    /** IPathGen's factory */
    /*  The reason why a virtual constructor is needed is the following.
        The implementation of an IPathGen must be done in the mcarlo directory, 
        since it requires -- among other things -- the evaluation of MCRandom-type
        objects. Yet an IPathGen must be created here, in the market directory,
        since only VolSVJ can populate the IPathGen with its parameters. */
    class MARKET_DLL IPathGenMaker{
    public:
        virtual IPathGen* make(
            double initialVol,
            double meanVol,
            double meanReversRate,
            double volVol,
            double correlation,
            double crashRate,
            double crashSizeMean,
            double crashSizeUncertainty,
            double volRiskPrice) const = 0;
        virtual ~IPathGenMaker(){}
    };
    typedef refCountPtr<IPathGenMaker> IPathGenMakerSP;
    typedef vector<IPathGenMakerSP> IPathGenMakerArray;

    MARKET_DLL IPathGen* createPathGen(const IPathGenMaker& maker) const;

    /** Given an array of gaussian deviates, populate the arrays 'instVars' 
        and 'integratedVars' with simulated values using an Euler
        discretization scheme */
    MARKET_DLL void MCGenerateVarPaths(Int2Type<MCParams::VarDSType::EULER>,
                                       const DoubleArray&  tradYears,
                                       const double*       deviates,
                                       DoubleArray&        instVars,
                                       DoubleArray&        integratedVars) const;

    /** Given an array of gaussian deviates, populate the arrays 'instVars' 
        and 'integratedVars' with simulated values using a variable transform
        scheme together with an Euler discretization scheme */
    MARKET_DLL void MCGenerateVarPaths(Int2Type<MCParams::VarDSType::VAR_TRANSFORM_EULER>,
                                       const DoubleArray&  tradYears,
                                       const double*       deviates,
                                       DoubleArray&        instVars,
                                       DoubleArray&        integratedVars) const;
#if 0

    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
        values using a var and a spot discretization scheme respectively */
    template <int SpotSchemeType, int VarSchemeType>
    void MCGenerateExJumpPaths(Int2Type<SpotSchemeType>,
                               Int2Type<VarSchemeType>,
                               const DoubleArray& tradYears,
                               const IntArray&    spotOffsets,
                               const double*      spotDeviates,
                               const double*      varDeviates,
                               QuantoParamSVJSP   quanto,
                               DoubleArray&       spots,
                               DoubleArray&       instVars,
                               DoubleArray&       integratedVars) const;
#endif

    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
        values using an Euler discretization scheme for the spot */
    template <int VarSchemeType>
    void MCGenerateExJumpPaths(Int2Type<MCParams::SpotDSType::EULER>,
                               Int2Type<VarSchemeType>,
                               const DoubleArray&     tradYears,
                               const IntArray&        spotOffsets,
                               const double*          spotDeviates,
                               const double*          varDeviates,
                               QuantoParamSVJSP       quanto,
                               DoubleArray&           exJumpLogSpots,
                               DoubleArray&           instVars,
                               DoubleArray&           integratedVars) const{
        // first of, generate the var paths
        MCGenerateVarPaths(Int2Type<VarSchemeType>(),
                           tradYears,
                           varDeviates,
                           instVars,
                           integratedVars);
        // given the vars, generate the log spots
#ifdef VOLSV_DEBUG
        ASSERT(exJumpLogSpots.size() == spotOffsets.size() + 1);
#endif
        double beta1 = correlation;
        double beta2 = sqrt(1.0 - beta1 * beta1);
        int nbSteps = tradYears.size();
        double logSpot = 0.0;
        exJumpLogSpots[0] = logSpot;
        int iOffset = 0, nextSpotStep = spotOffsets[0];
        // get quanto params
        bool isQuanto = quanto->isQuanto;
        double eqFXCorr = quanto->eqFXCorr;
        const DoubleArray& fxSqrtVar = quanto->fxSqrtVar;
        //loop
        for (int iStep = 1, iLastStep = 0,
             iCoarseStep = 1;
             iStep < nbSteps; iLastStep = iStep++) {
#ifdef VOLSV_DEBUG
            ASSERT(nextSpotStep < tradYears.size());
#endif
            double drift = -0.5 * (integratedVars[iStep] - integratedVars[iLastStep]);
            double dt = tradYears[iStep] - tradYears[iLastStep];
            double instVarPlus = Maths::max(0.0, instVars[iLastStep]);
            double diffusion = sqrt(instVarPlus * dt) 
                               * (beta1 * varDeviates[iLastStep]
                                  + beta2 * spotDeviates[iLastStep]);
            //if isQuanto, compute drift adjustment otherwise leave at zero
            double quantoAdjustment = 0.0;
            if(isQuanto){
                quantoAdjustment = - sqrt(Maths::max(instVars[iLastStep],0.0)) * eqFXCorr * fxSqrtVar[iLastStep];
            }
            logSpot += drift + diffusion + quantoAdjustment; //quantoAdjustment is added here
            if (nextSpotStep == iStep){
                exJumpLogSpots[iCoarseStep] = logSpot;
                iOffset = iCoarseStep++;
                if (iOffset < spotOffsets.size()){
                    nextSpotStep += spotOffsets[iOffset];
                }
            }
        }
    }

    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
        values. The log spots are simulated exactly */
    template <int VarSchemeType>
    void MCGenerateExJumpPaths(Int2Type<MCParams::SpotDSType::EXACT>,
                               Int2Type<VarSchemeType>,
                               const DoubleArray&   tradYears,
                               const IntArray&      spotOffsets,
                               const double*        spotDeviates,
                               const double*        varDeviates,
                               QuantoParamSVJSP     quanto,
                               DoubleArray&         exJumpLogSpots,
                               DoubleArray&         instVars,
                               DoubleArray&         integratedVars) const{
        // first of, generate the var paths
        MCGenerateVarPaths(Int2Type<VarSchemeType>(),
                           tradYears,
                           varDeviates,
                           instVars,
                           integratedVars);
        // given the vars, generate the log spots
#ifdef VOLSV_DEBUG
        ASSERT(exJumpLogSpots.size() == spotOffsets.size() + 1);
#endif
        double beta1 = correlation;
        double beta2 = sqrt(1.0 - beta1 * beta1);
        double sqMeanVol = Maths::square(meanVol);
        int nbCoarseSteps = spotOffsets.size() + 1;
        exJumpLogSpots[0] = 0.0;
        // get quanto params
        bool isQuanto = quanto->isQuanto;
        double eqFXCorr = quanto->eqFXCorr;
        const DoubleArray& fxSqrtVar = quanto->fxSqrtVar;
        //loop over
        int iLastFineStep = 0;
        for (int iCoarseStep = 1, iOffset = 0, iFineStep = 0; 
             iCoarseStep < nbCoarseSteps; 
             iOffset = iCoarseStep++) {
            iFineStep += spotOffsets[iOffset];
            int iLastCoarseStep = iOffset;
            // special treatment for the deterministic case
            if(Maths::isZero(volVol)){
                double integratedVarDiff = integratedVars[iFineStep] 
                                           - integratedVars[iLastFineStep];
                double mean = -0.5 * integratedVarDiff;
                double sqrtVar = sqrt(integratedVarDiff);
                //if isQuanto, compute drift adjustment otherwise leave at zero
                double quantoAdjustment = 0.0;
                if (isQuanto) {
                    for (int iQuanto = iLastFineStep; iQuanto < iFineStep - 1; iQuanto ++) {
                        quantoAdjustment -= sqrt(Maths::max(instVars[iQuanto],0.0)) * eqFXCorr * fxSqrtVar[iQuanto];
                    }
                }
                exJumpLogSpots[iCoarseStep] = exJumpLogSpots[iLastCoarseStep] 
                                        + mean + sqrtVar * spotDeviates[iOffset]
                                        +quantoAdjustment;//quanto drift adjustment is added here
            }
            // general case
            else{
                double dt = tradYears[iFineStep] - tradYears[iLastFineStep];
                double instVarDiff = instVars[iFineStep] 
                                     - instVars[iLastFineStep];
                double integratedVarDiff = integratedVars[iFineStep] 
                                           - integratedVars[iLastFineStep];
                double meanAdjTerm = (instVarDiff 
                                      - meanReversRate * (sqMeanVol * dt 
                                                          - integratedVarDiff))
                                     / volVol;
                double mean = -0.5 * integratedVarDiff + beta1 * meanAdjTerm;
                double sqrtVar = beta2 * sqrt(integratedVarDiff);
                //if isQuanto, compute drift adjustment otherwise leave at zero
                double quantoAdjustment = 0.0;
                if (isQuanto) {
                    for (int iQuanto = iLastFineStep; iQuanto < iFineStep - 1; iQuanto ++) {
                        quantoAdjustment -= sqrt(Maths::max(instVars[iQuanto],0.0)) * eqFXCorr * fxSqrtVar[iQuanto];
                    }
                }
                exJumpLogSpots[iCoarseStep] = exJumpLogSpots[iLastCoarseStep] 
                                        + mean + sqrtVar * spotDeviates[iOffset]
                                        +quantoAdjustment; //quanto drift adjustment added here
            }
            iLastFineStep = iFineStep;
        }
#ifdef VOLSV_DEBUG
        ASSERT(iLastFineStep < tradYears.size());
#endif
    }

    /** Given a uniform deviate, simulate a Poisson deviate with a given mean
        using the transform method.
        Based on algo 3.9 p. 128 Glasserman Monte Carlo Method in Finance */
    MARKET_DLL static int generatePoisson(double mean, double uni){
        double p = exp(-mean);
        double F = p;
        int N = 0;
        while (Maths::isPositive(uni - F)){
            ++N;
            p *= mean / N;
            F += p;
        }
        return N;
    }


    /** Given two arrays of (independent) gaussian deviates and one further, independent array 
        of uniform deviates, populate the arrays 'instVars', 'integratedVars' and 'logSpots' 
        with simulated values. */
    template <class SpotSchemeType, class VarSchemeType>
    void MCGeneratePaths(SpotSchemeType     spotScheme,
                         VarSchemeType      varScheme,
                         const DoubleArray& tradYears,
                         const IntArray&    spotOffsets,
                         const double*      spotDeviates,
                         const double*      varDeviates,
                         const double*      jumpUniDeviates,
                         const double*      jumpNormalDeviates,
                         QuantoParamSVJSP   quanto,
                         DoubleArray&       logSpots,
                         DoubleArray&       instVars,
                         DoubleArray&       integratedVars) const{
        // first of, generate the ex-jump path
        MCGenerateExJumpPaths(spotScheme,
                              varScheme,
                              tradYears,
                              spotOffsets,
                              spotDeviates,
                              varDeviates,
                              quanto,
                              logSpots,
                              instVars,
                              integratedVars); 

        // given the ex-jump log spots, generate the jumps and add on top of ex-jump
        // log spots 
#ifdef VOLSV_DEBUG
        ASSERT(logSpots.size() == spotOffsets.size() + 1);
#endif
        if (Maths::isPositive(crashRate)){
            double crashGamma = log(1.0 + crashSizeMean);
            double sqCrashSizeUncertainty = Maths::square(crashSizeUncertainty);
            double logJumpMean = crashGamma - 0.5 * sqCrashSizeUncertainty;
            double dJump = 0.0;
            int nbCoarseSteps = spotOffsets.size() + 1;
            int iLastFineStep = 0;
            for (int iCoarseStep = 1, iOffset = 0, iFineStep = 0; 
                 iCoarseStep < nbCoarseSteps; 
                 iOffset = iCoarseStep++) {
                iFineStep += spotOffsets[iOffset];
                double dt = tradYears[iFineStep] - tradYears[iLastFineStep];
                if (int nbJumps = generatePoisson(crashRate * dt, jumpUniDeviates[iOffset])){
                    double mean = nbJumps * logJumpMean;
                    double vol = sqrt(static_cast<double>(nbJumps)) * crashSizeUncertainty;
                    dJump += mean + vol * jumpNormalDeviates[iOffset];
                }
                dJump -= crashRate * crashSizeMean * dt;
                logSpots[iCoarseStep] += dJump;
                iLastFineStep = iFineStep;
            }
#ifdef VOLSV_DEBUG
            ASSERT(iLastFineStep < tradYears.size());
#endif
        }
    }

#if 0
    /** Given two arrays of (independent) gaussian deviates and one further, independent array 
        of uniform deviates, populate the arrays 'instVars', 'integratedVars' and 'logSpots' 
        with simulated values. */
    template <class SpotSchemeType, class VarSchemeType>
    void MCGeneratePaths2(SpotSchemeType     spotScheme,
                          VarSchemeType      varScheme,
                          const DoubleArray& tradYears,
                          const IntArray&    spotOffsets,
                          const double*      spotDeviates,
                          const double*      varDeviates,
                          MCRandPoisson&     poissonDeviates,
                          const double*      jumpNormalDeviates,
                          QuantoParamSVJSP   quanto,
                          DoubleArray&       logSpots,
                          DoubleArray&       instVars,
                          DoubleArray&       integratedVars) const{
        // first of, generate the ex-jump path
        MCGenerateExJumpPaths(spotScheme,
                              varScheme,
                              tradYears,
                              spotOffsets,
                              spotDeviates,
                              varDeviates,
                              quanto,
                              logSpots,
                              instVars,
                              integratedVars); 

        // given the ex-jump log spots, generate the jumps and add on top of ex-jump
        // log spots 
#ifdef VOLSV_DEBUG
        ASSERT(logSpots.size() == spotOffsets.size() + 1);
#endif
        if (Maths::isPositive(crashRate)){
            double crashGamma = log(1.0 + crashSizeMean);
            double sqCrashSizeUncertainty = Maths::square(crashSizeUncertainty);
            double logJumpMean = crashGamma - 0.5 * sqCrashSizeUncertainty;
            double dJump = 0.0;
            int nbCoarseSteps = spotOffsets.size() + 1;
            int iLastFineStep = 0;
            for (int iCoarseStep = 1, iOffset = 0, iFineStep = 0;
                iCoarseStep < nbCoarseSteps; 
                iOffset = iCoarseStep++;) {
                iFineStep = spotOffsets[0]; 
                double dt = tradYears[iFineStep] - tradYears[iLastFineStep];
                int nbJumps = static_cast<int>(poissonDeviates.draw(crashRate * dt));
                if (nbJumps){
                    double mean = nbJumps * logJumpMean;
                    double vol = sqrt(static_cast<double>(nbJumps)) * crashSizeUncertainty;
                    dJump += mean + vol * jumpNormalDeviates[iOffset];
                }
                dJump -= crashRate * crashSizeMean * dt;
                logSpots[iCoarseStep] += dJump;
                iLastFineStep = iFineStep;
            }
#ifdef VOLSV_DEBUG
            ASSERT(iLastFineStep < tradYears.size());
#endif
        }
    }
#endif
private:
    /** Called after adjustments have been made to fields */
    MARKET_DLL void update();

};
typedef array<VolSVJSP, VolSVJ>  VolSVJArray;
typedef VolSVJ::MCParams::SpotDSType VolSVJ_SpotDSType;
template<> MARKET_DLL string nameForType<VolSVJ_SpotDSType>(VolSVJ_SpotDSType*);
typedef VolSVJ::MCParams::SpotDSTypeHelper VolSVJ_SpotDSTypeHelper;
typedef VolSVJ::MCParams::VarDSType VolSVJ_VarDSType;
template<> MARKET_DLL string nameForType<VolSVJ_VarDSType>(VolSVJ_VarDSType*);
typedef VolSVJ::MCParams::VarDSTypeHelper VolSVJ_VarDSTypeHelper;

typedef MarketDataConvert::IConvert<VolSVJ> VolSVJConvert;

// For windows dlls you must EXTERN/INSTANTIATE the template for 
// Enum2StringListHelper if you want to use it in another dll
#ifndef QLIB_VOLSVJ_CPP
EXTERN_TEMPLATE(struct MARKET_DLL Enum2StringListHelper<VolSVJ::MCParams::SpotDSType>);
EXTERN_TEMPLATE(struct MARKET_DLL Enum2StringListHelper<VolSVJ::MCParams::VarDSType>);
#else
INSTANTIATE_TEMPLATE(struct MARKET_DLL Enum2StringListHelper<VolSVJ::MCParams::SpotDSType>);
INSTANTIATE_TEMPLATE(struct MARKET_DLL Enum2StringListHelper<VolSVJ::MCParams::VarDSType>);
#endif

DRLIB_END_NAMESPACE
#endif
