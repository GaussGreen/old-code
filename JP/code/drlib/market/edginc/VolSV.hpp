//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSV.hpp
//
//   Description : 
//
//   Date        : 25 April 03
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLSV_HPP
#define EDR_VOLSV_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Heston.hpp"
#include "edginc/VolAJDSuper.hpp"
#include "edginc/Enum2StringListHelper.hpp"
#include "edginc/Maths.hpp"
#include "edginc/IDynamicsParameter.hpp"
#include "edginc/MarketDataConvert.hpp"

DRLIB_BEGIN_NAMESPACE

//Quanto structure used for CCYProtected products
class QuantoParamSV{
public:
    bool            isQuanto;
    double          eqFXCorr;
    DoubleArray     fxSqrtVar;
    QuantoParamSV(bool quanto,double FXCorr,DoubleArray fxVar):isQuanto(quanto),eqFXCorr(FXCorr),fxSqrtVar(fxVar){};
    
private:
    QuantoParamSV();
};

typedef refCountPtr<QuantoParamSV> QuantoParamSVSP;
typedef array<QuantoParamSVSP,QuantoParamSV >  QuantoParamSVArray;


#define VOLSV_DEBUG 

template<int v>
class Int2Type{
    enum{
        value = v
    };
};

FORWARD_DECLARE(VolSVJ)
FORWARD_DECLARE(VolSVJJ)
FORWARD_DECLARE(VolSVCJ)

/** A parameterised CVolBase (does not support BS and DVF views of the
    world for now) */
class VolSV: public VolBaseParam,
             public virtual Calibrator::IAdjustable,
             public virtual VolAJDSuper::ISuperposable,
             public virtual IDynamicsParameter,
             public virtual MarketDataConvert::IConvert<VolSVJ>,
             public virtual MarketDataConvert::IConvert<VolSVJJ>,
             public virtual MarketDataConvert::IConvert<VolSVCJ>
/* virtual public IVolatilityBS,
   virtual public IVolatilityDVF */ {
public:
    /** Full constructor */
    MARKET_DLL VolSV(const string& name,
                     double initialVol,
                     double meanVol,
                     double meanReversRate,
                     double volVol,
                     double correlation);
    
    /** Upgrade to VolSVJ with zero jump parameters */
    virtual VolSVJSP convert(VolSVJ* p = 0) const;

    /** Upgrade to VolSVJJ with zero jump parameters */
    virtual VolSVJJSP convert(VolSVJJ* p = 0) const;

    /** Upgrade to VolSVCJ with zero jump parameters */
    virtual VolSVCJSP convert(VolSVCJ* p = 0) const;
    
private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class MARKET_DLL SVVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        SVVolParam(): CVolParam(TYPE){}

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
            REGISTER(SVVolParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };

    /* Methods called by SVVolParam. Implemented on VolSV to avoid dereferencing */ 
    MARKET_DLL void ComputeImpVol(const CLatticeDouble&      strikes,
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

    // transient fields (won't appear in dd interface)
    DateTime      baseDate;
    HestonSP      heston;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    MARKET_DLL VolSV();

    static IObject* defaultCtor();

protected:
    /*** Build the parameterised vol and cache any values **/
    MARKET_DLL void buildCache();

public:
    friend class SVVolParam;

    friend class FD2DSV;

    static MARKET_DLL CClassConstSP const TYPE;

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

    MARKET_DLL double getVolOfVol() const { return volVol; }
    MARKET_DLL double getMeanReversion() const { return meanReversRate; }

    MARKET_DLL void validatePop2Object();

    /** method that builds a CVolParam. */
    virtual MARKET_DLL  CVolParam* createVolParam() const;

    /** adjust mean vol and reversion rate for risk premium */
    MARKET_DLL void AdjustRiskPremium(double volRiskPrice, double& mVar, 
                                      double& mReversion) const;

    /* Methods called by FourierProcessSV. Implemented on VolSV to
       avoid dereferencing. Made public to avoid friend issue. */
    MARKET_DLL Complex scalelessCumulant(const StFourierProcessLogRtn& process,
                              const StFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;
    MARKET_DLL Complex scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                              const FwdStFourierProductLogRtn& product, 
                              const Complex& z, 
                              const DateTime& matDate) const;
    MARKET_DLL Complex cumulant(const StFourierProcessIntVar& process,
                     const StFourierProductIntVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;
    MARKET_DLL Complex cumulant(const FwdStFourierProcessIntVar& process,
                     const FwdStFourierProductIntVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;
    MARKET_DLL Complex cumulant(const StFourierProcessQuadVar& process,
                     const StFourierProductQuadVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;
    MARKET_DLL Complex cumulant(const FwdStFourierProcessQuadVar& process,
                     const FwdStFourierProductQuadVar& product, 
                     const Complex z, 
                     const DateTime& matDate) const;
    MARKET_DLL double expectation(const StFourierProcessQuadVar& process,
                       const StFourierProductQuadVar& product, 
                       const DateTime& matDate) const;
    MARKET_DLL double expectation(const FwdStFourierProcessQuadVar& process,
                       const FwdStFourierProductQuadVar& product, 
                       const DateTime& matDate) const;

    MARKET_DLL Complex cumulant(const FwdStFourierProcessExpQuadVar& process,              
                     const FwdStFourierProductExpQuadVar& product, 
                     const Complex& z, 
                     const DateTime& matDate) const;


        /* Used by VDax pricer */
    MARKET_DLL double varSwap(const DateTime& valueDate, const DateTime& maturity) const;
    MARKET_DLL double impliedSqVol(const DateTime& valueDate, const DateTime& maturity1, const DateTime& maturity2) const;

    /** Needed for IAdjustable interface. Returns market data name for vol */
    MARKET_DLL virtual string getName() const;

    /** Called after adjustments have been made to fields (eg calibrator) */
    MARKET_DLL virtual void fieldsUpdated(const CFieldArray& fields);

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

    /** Returns the nber of factors. The first factor is understood to be 
        the log spot. */
    MARKET_DLL virtual int getNbFactors() const;
    
    MARKET_DLL virtual void getMarket(const IModel*     model, 
                           const MarketData* market, 
                           const string&     name);

    /** Given an array of gaussian deviates, populate the arrays 'instVars' 
        and 'integratedVars' with simulated values using an Euler
        discretization scheme */
    MARKET_DLL void MCGenerateVarPaths(Int2Type<MCParams::VarDSType::EULER>,
                            const DoubleArray& tradYears,
                            const double*      deviates,
                            DoubleArray&       instVars,
                            DoubleArray&       integratedVars) const;

    /** Given an array of gaussian deviates, populate the arrays 'instVars' 
        and 'integratedVars' with simulated values using a variable transform
        scheme together with an Euler discretization scheme */
    MARKET_DLL void MCGenerateVarPaths(Int2Type<MCParams::VarDSType::VAR_TRANSFORM_EULER>,
                            const DoubleArray& tradYears,
                            const double*      deviates,
                            DoubleArray&       instVars,
                            DoubleArray&       integratedVars) const;
#if 0
    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'logSpots' with simulated 
        values using a var and a spot discretization scheme respectively */
    template <int SpotSchemeType, int VarSchemeType>
    void MCGeneratePaths(Int2Type<SpotSchemeType>,
                         Int2Type<VarSchemeType>,
                         const DoubleArray& tradYears,
                         const IntArray&    spotOffsets,
                         const double*      spotDeviates,
                         const double*      varDeviates,
                         QuantoParamSVSP        quanto,
                         DoubleArray&       spots,
                         DoubleArray&       instVars,
                         DoubleArray&       integratedVars) const;
#endif

    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'logSpots' with simulated 
        values using an Euler discretization scheme for the spot */
    template <int VarSchemeType>
    void MCGeneratePaths(Int2Type<MCParams::SpotDSType::EULER>,
                         Int2Type<VarSchemeType>,
                         const DoubleArray& tradYears,
                         const IntArray&    spotOffsets,
                         const double*      spotDeviates,
                         const double*      varDeviates,
                         QuantoParamSVSP        quanto,
                         DoubleArray&       logSpots,
                         DoubleArray&       instVars,
                         DoubleArray&       integratedVars) const{
        // first of, generate the var paths
        MCGenerateVarPaths(Int2Type<VarSchemeType>(),
                           tradYears,
                           varDeviates,
                           instVars,
                           integratedVars);
        // given the vars, generate the log spots
#ifdef VOLSV_DEBUG
        ASSERT(logSpots.size() == spotOffsets.size() + 1);
#endif
        double beta1 = correlation;
        double beta2 = sqrt(1.0 - beta1 * beta1);
        int nbSteps = tradYears.size();
        double logSpot = 0.0;
        logSpots[0] = logSpot;
        int iOffset = 0, nextSpotStep = spotOffsets[0];
        // get quanto params
        bool isQuanto = quanto->isQuanto;
        double eqFXCorr = quanto->eqFXCorr;
        DoubleArray fxSqrtVar = quanto->fxSqrtVar;
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
            //if isQuanto, compute drift adjustment
            double quantoAdjustment = 0.0;
            if (isQuanto) {
                quantoAdjustment -= sqrt(Maths::max(instVars[iLastStep],0.0)) * eqFXCorr * fxSqrtVar[iLastStep];
            }
            logSpot += drift + diffusion  + quantoAdjustment;
            if (nextSpotStep == iStep){
                logSpots[iCoarseStep] = logSpot;
                iOffset = iCoarseStep++;
                if (iOffset < spotOffsets.size()){
                    nextSpotStep += spotOffsets[iOffset];
                }
            }
        }
    }

    /** Given two arrays of (independent) gaussian deviates, populate the 
        arrays 'instVars', 'integratedVars' and 'logSpots' with simulated 
        values. The log spots are simulated exactly */
    template <int VarSchemeType>
    void MCGeneratePaths(Int2Type<MCParams::SpotDSType::EXACT>,
                         Int2Type<VarSchemeType>,
                         const DoubleArray& tradYears,
                         const IntArray&    spotOffsets,
                         const double*      spotDeviates,
                         const double*      varDeviates,
                         QuantoParamSVSP        quanto,
                         DoubleArray&       logSpots,
                         DoubleArray&       instVars,
                         DoubleArray&       integratedVars) const{
        // first of, generate the var paths
        MCGenerateVarPaths(Int2Type<VarSchemeType>(),
                           tradYears,
                           varDeviates,
                           instVars,
                           integratedVars);
        // given the vars, generate the log spots
#ifdef VOLSV_DEBUG
        ASSERT(logSpots.size() == spotOffsets.size() + 1);
#endif
        double beta1 = correlation;
        double beta2 = sqrt(1.0 - beta1 * beta1);
        double sqMeanVol = Maths::square(meanVol);
        int nbCoarseSteps = spotOffsets.size() + 1;
        logSpots[0] = 0.0;
        // get quanto params
        bool isQuanto = quanto->isQuanto;
        double eqFXCorr = quanto->eqFXCorr;
        DoubleArray fxSqrtVar = quanto->fxSqrtVar;
        int iLastFineStep = 0;
		if (nbCoarseSteps > 1){
			for (int iCoarseStep = 1, iOffset = 0,
				iFineStep = spotOffsets[0]; 
				iCoarseStep < nbCoarseSteps;) {
				int iLastCoarseStep = iOffset;
				// special treatment for the deterministic case
				if(Maths::isZero(volVol)){
					double integratedVarDiff = integratedVars[iFineStep] 
											- integratedVars[iLastFineStep];
					double mean = -0.5 * integratedVarDiff;
					double sqrtVar = sqrt(integratedVarDiff);
                    //if isQuanto, compute drift adjustment
                    double quantoAdjustment = 0.0;
                    if (isQuanto) {
                        for (int iQuanto = iLastFineStep; iQuanto < iFineStep - 1; iQuanto ++) {
                            quantoAdjustment -= sqrt(Maths::max(instVars[iQuanto],0.0)) * eqFXCorr * fxSqrtVar[iQuanto];
                        }
                    }
					logSpots[iCoarseStep] = logSpots[iLastCoarseStep] 
											+ mean + sqrtVar * spotDeviates[iOffset] + quantoAdjustment;
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
                    //if isQuanto, compute drift adjustment
                    double quantoAdjustment = 0.0;
                    if (isQuanto) {
                        for (int iQuanto = iLastFineStep; iQuanto < iFineStep - 1; iQuanto ++) {
                            quantoAdjustment -= sqrt(Maths::max(instVars[iQuanto],0.0)) * eqFXCorr * fxSqrtVar[iQuanto];
                        }
                    }
					logSpots[iCoarseStep] = logSpots[iLastCoarseStep] 
											+ mean + sqrtVar * spotDeviates[iOffset] + quantoAdjustment;
				}
				iLastFineStep = iFineStep;
				iOffset = iCoarseStep++;
				if (iOffset < spotOffsets.size()){
					iFineStep += spotOffsets[iOffset];
				}
			}
		}
#ifdef VOLSV_DEBUG
        ASSERT(iLastFineStep < tradYears.size());
#endif
    }
private:
    /** Called after adjustments have been made to fields */
    MARKET_DLL void update();
};
typedef smartPtr<VolSV> VolSVSP;
typedef smartConstPtr<VolSV> VolSVConstSP;
typedef array<VolSVSP, VolSV>  VolSVArray;
typedef VolSV::MCParams::SpotDSType VolSV_SpotDSType;
template<> MARKET_DLL string nameForType<VolSV_SpotDSType>(VolSV_SpotDSType*);
typedef VolSV::MCParams::SpotDSTypeHelper VolSV_SpotDSTypeHelper;
typedef VolSV::MCParams::VarDSType VolSV_VarDSType;
template<> MARKET_DLL string nameForType<VolSV_VarDSType>(VolSV_VarDSType*);
typedef VolSV::MCParams::VarDSTypeHelper VolSV_VarDSTypeHelper;

typedef MarketDataConvert::IConvert<VolSV> VolSVConvert;

// For windows dlls you must EXTERN/INSTANTIATE the template for 
// Enum2StringListHelper if you want to use it in another dll
#ifndef QLIB_VOLSV_CPP
EXTERN_TEMPLATE(struct MARKET_DLL Enum2StringListHelper<VolSV::MCParams::SpotDSType>);
EXTERN_TEMPLATE(struct MARKET_DLL Enum2StringListHelper<VolSV::MCParams::VarDSType>);
#else
INSTANTIATE_TEMPLATE(struct MARKET_DLL Enum2StringListHelper<VolSV::MCParams::SpotDSType>);
INSTANTIATE_TEMPLATE(struct MARKET_DLL Enum2StringListHelper<VolSV::MCParams::VarDSType>);
#endif

DRLIB_END_NAMESPACE
#endif
