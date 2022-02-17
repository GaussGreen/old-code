//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSVJ.cpp
//
//   Description : Monte Carlo path generator for VolSVJ
//
//   Date        : 20 Sep 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/MCPathConfigParametric.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/VolSVJ.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/StruckAsset.hpp"
#include "edginc/DependenceGauss.hpp"

#define STATE_VARIABLES

#ifdef STATE_VARIABLES
#include "edginc/mathlib.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#endif

DRLIB_BEGIN_NAMESPACE

#ifdef STATE_VARIABLES
///////////////////////////////////////////////////////////////////////////////////
// VolSVJ::IPathGen's children
///////////////////////////////////////////////////////////////////////////////////
struct VarPathGenEuler{
    template<class NotUsed>
    class In{
    protected:
        /** Given an array of gaussian deviates, populate the arrays 'instVars' 
            and 'integratedVars' with simulated values using an Euler
            discretization scheme */
        void generatePath(const DoubleArray&        tradYears,
                          const IMCRandNormal::Path& deviates,
                          DoubleArray&              instVars,
                          DoubleArray&              integratedVars) const{ 
            ASSERT(instVars.size() == tradYears.size());
            ASSERT(instVars.size() == integratedVars.size());
            int nbSteps = tradYears.size();
            double meanVar = Maths::square(meanVol);
            instVars[0] = Maths::square(initialVol);
            integratedVars[0] = 0.0;
            for (int iStep = 1, iLastStep = 0; iStep < nbSteps; 
                 iLastStep = iStep++) {
                double dt = tradYears[iStep] - tradYears[iLastStep];
                // drift contribution
                double drift = meanReversRate * (meanVar - instVars[iLastStep]) * dt;
                // diffusion contribution
                double instVarPlus = Maths::max(0.0, instVars[iLastStep]);
                double diffusion = volVol * sqrt(instVarPlus * dt) * deviates[iLastStep];
                // compute inst variance
                instVars[iStep] = instVars[iLastStep] + drift + diffusion;
                // compute integrated var using an Euler scheme too
                integratedVars[iStep] = integratedVars[iLastStep] + instVarPlus * dt;
            }
        }

        // seems to be needed by gcc !
        In(){}

        In(double initialVol,
           double meanVol,
           double meanReversRate,
           double volVol,
           double volRiskPrice):
        initialVol(initialVol),
        meanVol(meanVol),
        meanReversRate(meanReversRate),
        volVol(volVol),
        volRiskPrice(volRiskPrice){}

        // fields
        double initialVol;
        double meanVol;
        double meanReversRate;
        double volVol;
        double volRiskPrice;
    };
};

struct VarPathGenTransformEuler{
    template<class NotUsed>
    class In{
    protected:
        /** Given an array of gaussian deviates, populate the arrays 'instVars' 
            and 'integratedVars' with simulated values using a variable transform
            scheme together with an Euler discretization scheme */
        void generatePath(const DoubleArray&        tradYears,
                          const IMCRandNormal::Path& deviates,
                          DoubleArray&              instVars,
                          DoubleArray&              integratedVars) const{
            ASSERT(instVars.size() == tradYears.size());
            ASSERT(instVars.size() == integratedVars.size());
            int nbSteps = tradYears.size();
            double alpha = meanReversRate * Maths::square(meanVol) 
                           - Maths::square(volVol) / 4.0;
            double vol = initialVol;
            instVars[0] = Maths::square(vol);
            integratedVars[0] = 0.0;
            for (int iStep = 1, iLastStep = 0; iStep < nbSteps; 
                 iLastStep = iStep++) {
                double dt = tradYears[iStep] - tradYears[iLastStep];
                // integrated var is function of last step values only
                integratedVars[iStep] = integratedVars[iLastStep] + instVars[iLastStep] * dt;
                // Euler approx for the vol
                vol = sqrt(instVars[iLastStep]);
                ASSERT(!Maths::isZero(vol));
                double drift = 0.5 * (alpha / vol - meanReversRate * vol) * dt;
                double diffusion = 0.5 * volVol * sqrt(dt) * deviates[iLastStep];
                vol += drift + diffusion;
                instVars[iStep] = Maths::square(vol);
            }
        }

        // seems to be needed by gcc !
        In(){}

        In(double initialVol,
           double meanVol,
           double meanReversRate,
           double volVol,
           double volRiskPrice):
        initialVol(initialVol),
        meanVol(meanVol),
        meanReversRate(meanReversRate),
        volVol(volVol),
        volRiskPrice(volRiskPrice){}

        // fields
        double initialVol;
        double meanVol;
        double meanReversRate;
        double volVol;
        double volRiskPrice;
    };
};

class NullType{};

struct ExJumpSpotPathGenEuler{
    template<class VarPathGen>
    class In: public MSVC6Helper::Apply1<VarPathGen, NullType>::Base{
        typedef typename MSVC6Helper::Apply1<VarPathGen, NullType>::Base Base;
    protected:
        /** Given two arrays of (independent) gaussian deviates, populate the 
            arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
            values using an Euler discretization scheme for the spot */
        void generatePath(const DoubleArray&         tradYears,
                          const IntArray&            spotOffsets,
                          const IMCRandNormal::Path& spotDeviates,
                          const IMCRandNormal::Path& varDeviates,
                          QuantoParamSVJSP           quanto,
                          DoubleArray&               exJumpLogSpots,
                          DoubleArray&               instVars,
                          DoubleArray&               integratedVars) const{
            // first of, generate the var paths
            Base::generatePath(tradYears,
                               varDeviates,
                               instVars,
                               integratedVars);
            // given the vars, generate the log spots
            ASSERT(exJumpLogSpots.size() == spotOffsets.size() + 1);
            double beta1 = correlation;
            double beta2 = sqrt(1.0 - beta1 * beta1);
            int nbSteps = tradYears.size();
            double logSpot = 0.0;
            exJumpLogSpots[0] = logSpot;
            int iOffset = 0, nextSpotStep = spotOffsets[0];
            // get quanto params
            bool isQuanto = quanto->isQuanto;
            double eqFXCorr = quanto->eqFXCorr;
            DoubleArray & fxSqrtVar = quanto->fxSqrtVar;
            for (int iStep = 1, iLastStep = 0,
                 iCoarseStep = 1;
                 iStep < nbSteps; iLastStep = iStep++) {
                ASSERT(nextSpotStep < tradYears.size());
                double drift = -0.5 * (integratedVars[iStep] - integratedVars[iLastStep]);
                double dt = tradYears[iStep] - tradYears[iLastStep];
                double instVarPlus = Maths::max(0.0, instVars[iLastStep]);
                double diffusion = sqrt(instVarPlus * dt) 
                                   * (beta1 * varDeviates[iLastStep]
                                      + beta2 * spotDeviates[iLastStep]);
                //if isQuanto, compute drift adjustment
                double quantoAdjustment = 0.0;
                if (isQuanto) {
                    quantoAdjustment -= sqrt(instVars[iLastStep]) * eqFXCorr * fxSqrtVar[iLastStep];
                }
                logSpot += drift + diffusion + quantoAdjustment;
                if (nextSpotStep == iStep){
                    exJumpLogSpots[iCoarseStep] = logSpot;
                    iOffset = iCoarseStep++;
                    if (iOffset < spotOffsets.size()){
                        nextSpotStep += spotOffsets[iOffset];
                    }
                }
            }
        }

        // seems to be needed by gcc !
        In(){}

        In(double initialVol,
           double meanVol,
           double meanReversRate,
           double volVol,
           double correlation,
           double volRiskPrice):
        Base(initialVol,
             meanVol,
             meanReversRate,
             volVol,
             volRiskPrice),
        correlation(correlation){}

        // fields
        double correlation;
    };
};

struct ExJumpSpotPathGenExact{
    template<class VarPathGen>
    class In: public MSVC6Helper::Apply1<VarPathGen, NullType>::Base{
        typedef typename MSVC6Helper::Apply1<VarPathGen, NullType>::Base Base;

    protected:
        /** Given two arrays of (independent) gaussian deviates, populate the 
            arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
            values. The log spots are simulated exactly */
        void generatePath(const DoubleArray&         tradYears,
                          const IntArray&            spotOffsets,
                          const IMCRandNormal::Path& spotDeviates,
                          const IMCRandNormal::Path& varDeviates,
                          QuantoParamSVJSP           quanto,
                          DoubleArray&               exJumpLogSpots,
                          DoubleArray&               instVars,
                          DoubleArray&               integratedVars) const{
            // first of, generate the var paths
            Base::generatePath(tradYears,
                               varDeviates,
                               instVars,
                               integratedVars);
            // given the vars, generate the log spots
            ASSERT(exJumpLogSpots.size() == spotOffsets.size() + 1);
            double beta1 = correlation;
            double beta2 = sqrt(1.0 - beta1 * beta1);
            double sqMeanVol = Maths::square(this->meanVol);
            int nbCoarseSteps = spotOffsets.size() + 1;
            exJumpLogSpots[0] = 0.0;
            int iLastFineStep = 0;
            // get quanto params
            bool isQuanto = quanto->isQuanto;
            double eqFXCorr = quanto->eqFXCorr;
            DoubleArray & fxSqrtVar = quanto->fxSqrtVar;
            for (int iCoarseStep = 1, iOffset = 0,
                 iFineStep = spotOffsets[0]; 
                 iCoarseStep < nbCoarseSteps; 
                 iFineStep += spotOffsets[iOffset]) {
                int iLastCoarseStep = iOffset;
                // special treatment for the deterministic case
                if(Maths::isZero(this->volVol)){
                    double integratedVarDiff = integratedVars[iFineStep] 
                                               - integratedVars[iLastFineStep];
                    double mean = -0.5 * integratedVarDiff;
                    double sqrtVar = sqrt(integratedVarDiff);
                    double quantoAdjustment = 0.0;
                    if (isQuanto) {
                        for (int iQuanto = iLastFineStep; iQuanto < iFineStep - 1; iQuanto++) {
                        quantoAdjustment -= sqrt(instVars[iQuanto]) * eqFXCorr * fxSqrtVar[iQuanto];
                        }
                    }
                    exJumpLogSpots[iCoarseStep] = exJumpLogSpots[iLastCoarseStep] 
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
                                          - this->meanReversRate * 
                                          (sqMeanVol * dt 
                                                              - integratedVarDiff))
                                         / this->volVol;
                    double mean = -0.5 * integratedVarDiff + beta1 * meanAdjTerm;
                    double sqrtVar = beta2 * sqrt(integratedVarDiff);
                    double quantoAdjustment = 0.0;
                    if (isQuanto) {
                        for (int iQuanto = iLastFineStep; iQuanto < iFineStep - 1; iQuanto++) {
                        quantoAdjustment -= sqrt(instVars[iQuanto]) * eqFXCorr * fxSqrtVar[iQuanto];
                        }
                    }
                    exJumpLogSpots[iCoarseStep] = exJumpLogSpots[iLastCoarseStep] 
                                            + mean + sqrtVar * spotDeviates[iOffset] + quantoAdjustment;
                }
                iLastFineStep = iFineStep;
                iOffset = iCoarseStep++;
            }
            ASSERT(iLastFineStep < tradYears.size());
        }

        // seems to be needed by gcc !
        In(){}

        In(double initialVol,
           double meanVol,
           double meanReversRate,
           double volVol,
           double correlation,
           double volRiskPrice):
        Base(initialVol,
             meanVol,
             meanReversRate,
             volVol,
             volRiskPrice),
        correlation(correlation){}

        // fields
        double correlation;
    };
};

template<class VarPathGen, class ExJumpSpotPathGen>
class SpotPathGenExact: public MSVC6Helper::Apply1<ExJumpSpotPathGen, VarPathGen>::Base{
    typedef typename MSVC6Helper::Apply1<ExJumpSpotPathGen, VarPathGen>::Base Base;

protected:
    /** Given two arrays of (independent) gaussian deviates and one further, independent array 
        of uniform deviates, populate the arrays 'instVars', 'integratedVars' and 'logSpots' 
        with simulated values. */
    void generatePath(const DoubleArray&         tradYears,
                      const IntArray&            spotOffsets,
                      const IMCRandNormal::Path& spotDeviates,
                      const IMCRandNormal::Path& varDeviates,
                      MCRandPoisson&             jumpNbDeviates,
                      const IMCRandNormal::Path& jumpSizeDeviates,
                      QuantoParamSVJSP           quanto,
                      DoubleArray&               logSpots,
                      DoubleArray&               instVars,
                      DoubleArray&               integratedVars) const{
        // first of, generate the ex-jump path
        Base::generatePath(tradYears,
                           spotOffsets,
                           spotDeviates,
                           varDeviates,
                           quanto,
                           logSpots,
                           instVars,
                           integratedVars); 
        // given the ex-jump log spots, generate the jumps and add on top of ex-jump
        // log spots 
        ASSERT(logSpots.size() == spotOffsets.size() + 1);
        if (Maths::isPositive(crashRate)){
            double crashGamma = log(1.0 + crashSizeMean);
            double sqCrashSizeUncertainty = Maths::square(crashSizeUncertainty);
            double logJumpMean = crashGamma - 0.5 * sqCrashSizeUncertainty;
            double dJump = 0.0;
            int nbCoarseSteps = spotOffsets.size() + 1;
            int iLastFineStep = 0;
            for (int iCoarseStep = 1, iOffset = 0,
                 iFineStep = spotOffsets[0]; 
                 iCoarseStep < nbCoarseSteps; 
                 iFineStep += spotOffsets[iOffset]) {
//                int iLastCoarseStep = iOffset;
                double dt = tradYears[iFineStep] - tradYears[iLastFineStep];
                int nbJumps = static_cast<int>(jumpNbDeviates.draw(crashRate * dt));
                if (nbJumps){
                    double mean = nbJumps * logJumpMean;
                    double vol = sqrt(static_cast<double>(nbJumps)) * crashSizeUncertainty;
                    dJump += mean + vol * jumpSizeDeviates[iOffset];
                }
                dJump -= crashRate * crashSizeMean * dt;
                logSpots[iCoarseStep] += dJump;
                iLastFineStep = iFineStep;
                iOffset = iCoarseStep++;
            }
            ASSERT(iLastFineStep < tradYears.size());
        }
    }

    SpotPathGenExact(double initialVol,
                     double meanVol,
                     double meanReversRate,
                     double volVol,
                     double correlation,
                     double crashRate,
                     double crashSizeMean,
                     double crashSizeUncertainty,
                     double volRiskPrice):
    Base(initialVol,
         meanVol,
         meanReversRate,
         volVol,
         correlation,
         volRiskPrice),
    crashRate(crashRate),
    crashSizeMean(crashSizeMean),
    crashSizeUncertainty(crashSizeUncertainty){}

    // fields
    double crashRate;
    double crashSizeMean;
    double crashSizeUncertainty;
};

template<class VarPathGen, class ExJumpSpotPathGen>
class SpotQuadVarPathGenExact: public MSVC6Helper::Apply1<ExJumpSpotPathGen, VarPathGen>::Base{
    typedef typename MSVC6Helper::Apply1<ExJumpSpotPathGen, VarPathGen>::Base Base;

protected:
    /** Given two arrays of (independent) gaussian deviates and one further, independent array 
        of uniform deviates, populate the arrays 'instVars', 'integratedVars' and 'logSpots' 
        with simulated values. */
    void generatePath(const DoubleArray&         tradYears,
                      const IntArray&            spotOffsets,
                      const IMCRandNormal::Path& spotDeviates,
                      const IMCRandNormal::Path& varDeviates,
                      MCRandPoisson&             jumpNbDeviates,
                      MCRandNormal&              jumpSizeDeviates,
                      QuantoParamSVJSP           quanto,
                      DoubleArray&               logSpots,
                      DoubleArray&               instVars,
                      DoubleArray&               quadVars) const{
        // first of, generate the ex-jump path
        DoubleArray& integratedVars = quadVars; // just an alias
        Base::generatePath(tradYears,
                           spotOffsets,
                           spotDeviates,
                           varDeviates,
                           quanto,
                           logSpots,
                           instVars,
                           integratedVars); 
        // given the ex-jump log spots, generate the jumps and add on top of ex-jump
        // log spots 
        ASSERT(logSpots.size() == spotOffsets.size() + 1);
        if (Maths::isPositive(crashRate)){
            double crashGamma = log(1.0 + crashSizeMean);
            double sqCrashSizeUncertainty = Maths::square(crashSizeUncertainty);
            double logJumpMean = crashGamma - 0.5 * sqCrashSizeUncertainty;
            double dJump = 0.0;
            double dSqJump = 0.0;
            int nbCoarseSteps = spotOffsets.size() + 1;
            int iLastFineStep = 0;
            for (int iCoarseStep = 1, iOffset = 0,
                 iFineStep = spotOffsets[0]; 
                 iCoarseStep < nbCoarseSteps; 
                 iFineStep += spotOffsets[iOffset]) {
//                int iLastCoarseStep = iOffset;
                double dt = tradYears[iFineStep] - tradYears[iLastFineStep];
                int nbJumps = static_cast<int>(jumpNbDeviates.draw(crashRate * dt));
                while (nbJumps){
                    double jump = logJumpMean + crashSizeUncertainty * jumpSizeDeviates.draw();
                    dJump += jump;
                    dSqJump += jump * jump;
                    --nbJumps;
                }
                dJump -= crashRate * crashSizeMean * dt;
                logSpots[iCoarseStep] += dJump;
                quadVars[iFineStep] += dSqJump;
                iLastFineStep = iFineStep;
                iOffset = iCoarseStep++;
            }
            ASSERT(iLastFineStep < tradYears.size());
        }
    }

    SpotQuadVarPathGenExact(double initialVol,
                     double meanVol,
                     double meanReversRate,
                     double volVol,
                     double correlation,
                     double crashRate,
                     double crashSizeMean,
                     double crashSizeUncertainty,
                     double volRiskPrice):
    Base(initialVol,
         meanVol,
         meanReversRate,
         volVol,
         correlation,
         volRiskPrice),
    crashRate(crashRate),
    crashSizeMean(crashSizeMean),
    crashSizeUncertainty(crashSizeUncertainty){}

    // fields
    double crashRate;
    double crashSizeMean;
    double crashSizeUncertainty;
};
#endif

///////////////////////////////////////////////////////////////////////////////////
// MCPathConfig
///////////////////////////////////////////////////////////////////////////////////
class MCPathConfigSVJ: public MCPathConfigParametric {
public: 
    static CClassConstSP const TYPE;

    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
        for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const {
        return MarketDataFetcherSP(new MDFAssetVol(VolSVJ::TYPE->getName(), 
                                                   VolSurface::TYPE->getName()));
    }

    class PathGenFineProcess;
    class Generator;

    class Random: public RandomImpl{
    public:
        friend class PathGenFineProcess;
        friend class Generator;

        static CClassConstSP const TYPE;

        Random(int seed):
        RandomImpl(TYPE), 
        uniRand(new RandUniformDefault(seed)){
            validatePop2Object();
        }

        Random(const RandUniformDefaultSP& uniRand):
        RandomImpl(TYPE), uniRand(uniRand){
            validatePop2Object();
        }

        void validatePop2Object(){
            init();
        }

        /** Initializes the genarator */
        virtual void init() {
            normal = RandNormalDefaultSP(new RandNormalDefault(uniRand));
            normal->init();
            poisson = RandPoissonDefaultSP(new RandPoissonDefault(uniRand));
            poisson->init();
        }

        /** returns a State object capturing the state of the IRandom */
        virtual IRandom::State* getState() const{
            State* state = new State();
            state->poisson = IRandom::StateSP(poisson->getState());
            state->normal = IRandom::StateSP(normal->getState());
            return state;
        }

        /** restores the state of an IRandom */
        virtual void setState(const IRandom::State* state){
            const State& myState = dynamic_cast<const State&>(*state);
            poisson->setState(myState.poisson.get());
            normal->setState(myState.normal.get());
        }

    private:
        Random():
        RandomImpl(TYPE){}

        class State: public IRandom::State{
        friend class Random;
        private:
            IRandom::StateSP  normal;
            IRandom::StateSP  poisson;
        };
        typedef refCountPtr<State> StateSP;

        static IRandom* fromRandom(int seed){
            return new Random(seed);
        }

        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            REGISTER(Random, clazz);
            SUPERCLASS(RandomImpl);
            clazz->setPublic(); // make visible to EAS/spreadsheet
            EMPTY_SHELL_METHOD(defaultCtor);
            FIELD(uniRand, "Uniform random generator");
            FIELD(normal, "");
            FIELD_MAKE_TRANSIENT(normal);
            FIELD(poisson, "");
            FIELD_MAKE_TRANSIENT(poisson);
            NAMESPACE::Random::registerIRandomFromRandomMethod(TYPE, fromRandom);
        }

        static IObject* defaultCtor(){
            return new Random();
        }

        // registered field
        RandUniformDefaultSP uniRand;
        // transient fields
        RandNormalDefaultSP  normal;
        RandPoissonDefaultSP poisson;
    };

#ifdef STATE_VARIABLES
    template<class VarPathGen, class ExJumpSpotPathGen>
    class PathGen: public SpotQuadVarPathGenExact<VarPathGen, ExJumpSpotPathGen>,
                   public VolSVJ::IPathGen{
        typedef SpotQuadVarPathGenExact<VarPathGen, ExJumpSpotPathGen> Base;

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
              DoubleArray&               quadVars) const{
            Base::generatePath(tradYears,
                               spotOffsets,
                               spotDeviates,
                               varDeviates,
                               jumpNbDeviates,
                               jumpSizeDeviates,
                               quanto,
                               logSpots,
                               instVars,
                               quadVars);
        }

        PathGen(double initialVol,
                double meanVol,
                double meanReversRate,
                double volVol,
                double correlation,
                double crashRate,
                double crashSizeMean,
                double crashSizeUncertainty,
                double volRiskPrice):
        Base(initialVol,
             meanVol,
             meanReversRate,
             volVol,
             correlation,
             crashRate,
             crashSizeMean,
             crashSizeUncertainty,
             volRiskPrice){}
    };

    template<class VarPathGen, class ExJumpSpotPathGen>
    class PathGenMaker: public VolSVJ::IPathGenMaker{
    public:
        virtual VolSVJ::IPathGen* make(
                double initialVol,
                double meanVol,
                double meanReversRate,
                double volVol,
                double correlation,
                double crashRate,
                double crashSizeMean,
                double crashSizeUncertainty,
                double volRiskPrice) const{
            return new PathGen<VarPathGen, ExJumpSpotPathGen>(
                            initialVol,
                            meanVol,
                            meanReversRate,
                            volVol,
                            correlation,
                            crashRate,
                            crashSizeMean,
                            crashSizeUncertainty,
                            volRiskPrice);
        }

        PathGenMaker(){}
    };
#endif

private: 
#ifdef STATE_VARIABLES
    // cheap dispatcher method
    VolSVJ::IPathGenMaker* createPathGenMaker(){
        static const string method("MCPathConfigSVJ::createPathGenMaker");
        // XXX this list of nested if statements is not pretty...
        if (varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::EULER){
            if (spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EULER){
                return new PathGenMaker<VarPathGenEuler, ExJumpSpotPathGenEuler>();
            }
            else if(spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EXACT){
                return new PathGenMaker<VarPathGenEuler, ExJumpSpotPathGenExact>();
            }
            else{ 
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_SpotDSTypeHelper::getName(spotDiscreteSchemeIndex));
            }
        }
        else if (varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
            if (spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EULER){
                return new PathGenMaker<VarPathGenTransformEuler, ExJumpSpotPathGenEuler>();
            }
            else if(spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EXACT){
                return new PathGenMaker<VarPathGenTransformEuler, ExJumpSpotPathGenExact>();
            }
            else{
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_SpotDSTypeHelper::getName(spotDiscreteSchemeIndex));
            }
        }
        else{
            throw ModelException(method, 
                                 "unexpected variance discretization scheme of type "
                                 + VolSVJ_VarDSTypeHelper::getName(varDiscreteSchemeIndex));
        }
    }
#endif

    virtual void validatePop2Object(){
        static const string routine("MCPathConfigSVJ::validatePop2Object");
        try {
            // converts the 'rand' to the appropriate MCPathConfigSVJ::Random type if needed
            // typically, this means that MCPathConfigSVJ can be built with an old world
            // 'Random' type rand and the rand will be converted on the fly to an MCPathConfigSVJ::Random
            rand = IRandomSP(rand->convert(Random::TYPE));
            // further validation
            Maths::checkPositive(nbStepsPerYear, "nbStepsPerYear");
            spotDiscreteSchemeIndex = VolSVJ_SpotDSTypeHelper::getIndex(spotDiscreteScheme);
            varDiscreteSchemeIndex = VolSVJ_VarDSTypeHelper::getIndex(varDiscreteScheme);
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Returns the number of bytes used for random number storage per
        path. Do not invoke if there are no sim dates in future */
    virtual int randomStoragePerPath(IMCProduct* product) const{
        // Otherwise have to worry about how many random numbers we use
        // This is a bit of a pain as currently there is no easy way to get hold
        // of the number of dates unless we build the entire path generator
        smartPtr<MCPathConfigSVJ> pathConfig(copy(this)); // copy to avoid const problems
        MCPathGeneratorSP pastPathGenerator(pathConfig->pastPathGenerator(product));
        SensitivityArrayConstSP    sens;
        OutputRequestArrayConstSP  request;
        Control control(sens, request, false, "");
        Results results;
        DateTimeArray simDates;
        MCPathGeneratorSP pathGen(pathConfig->futurePathGenerator(0, // no caching
                                                                  2, // num paths 
                                                                  pastPathGenerator,
                                                                  product,
                                                                  &control,
                                                                  &results,
                                                                  simDates ));
        if (!pathConfig->simDates){
            throw ModelException("MCPathConfigSVJ::storagePerPath", "Internal "
                                 "error - no simulation dates");
        }
        // we store both correlated and uncorrelated numbers but only for
        // every other path. Hence no times by 2.
        return (sizeof(double) * pathConfig->simDates->size() * 
                product->getNumAssets());
    }

    virtual bool vegaMatrixSupported() const {
        return true;
    }

    virtual bool carefulRandoms() const {
        return isCarefulRandoms;
    }

    /** Throws away cached sim dates for Theta-type tweaks */
    virtual bool sensShift(Theta* shift) {
        MCPathConfig::sensShift(shift); // call parents method
        simDates.reset();
        return true;
    }

    /** Creates a future path generator */
    MCPathGeneratorSP makePathGenerator(
        bool                               cachingRequested,    
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                  prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     simDates );

    /** Creates a past path generator */
    MCPathGeneratorSP pastPathGenerator(const IMCProduct* prod);

    /** Helper function used by both Generator and PathGenFineProcess.
        Creates a time line that is uniformally spaced between each
        fixing date aka 'futureDate' */
    static DateTimeArraySP genSimDates(const DateTimeArray& futureDates,
                                       const TimeMetric&    timeMetric,
                                       int                  nbStepsPerYear,
                                       IntArray&            fixDateOffsets) {
        static const string method("MCPathConfigSVJ::Generator::genSimDates");
        try{
            // calculate approx number of sim dates
            int nbFutureDates = futureDates.size();
            ASSERT(nbFutureDates > 1);
            DateTime startDate = futureDates.front();
            DateTime horizonDate = futureDates.back();
            double horizonYearFrac = timeMetric.yearFrac(startDate, horizonDate);
            int approxNbSteps = static_cast<int>(ceil(horizonYearFrac * nbStepsPerYear));
            DateTimeArraySP simDates(new DateTimeArray());
            simDates->reserve(approxNbSteps + 1);
            // loop over future dates and insert additional points as needed
            // calculate fixing date offsets simultaneously 
            fixDateOffsets.resize(nbFutureDates - 1);
            DateTime currFutureDate = startDate;
            simDates->push_back(currFutureDate);
            for (int iDate = 1, iOffset = 0; iDate < nbFutureDates; iOffset = iDate++){
                // calculate mesh size
                DateTime nextFutureDate = futureDates[iDate];
                double yearFrac = timeMetric.yearFrac(currFutureDate, nextFutureDate);
                int nbSteps = static_cast<int>(ceil(yearFrac * nbStepsPerYear));
                double delta = yearFrac / nbSteps;
                // create sim dates
                DateTime currDate = currFutureDate; 
                fixDateOffsets[iOffset] = 0;
                for (int iStep = 1; iStep < nbSteps; ++iStep){
                    double notused;
                    currDate = timeMetric.impliedTime(currDate,
                                                      delta,
                                                      notused);
                    simDates->push_back(currDate);
                    fixDateOffsets[iOffset]++;
                }
                simDates->push_back(nextFutureDate);
                fixDateOffsets[iOffset]++;
                currFutureDate = nextFutureDate;
            }
            return simDates;
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(MCPathConfigSVJ, clazz);
        SUPERCLASS(MCPathConfigParametric);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        EMPTY_SHELL_METHOD(defaultMCPathConfigSVJ);
        FIELD(dependenceType, "not used");
        FIELD_MAKE_OPTIONAL(dependenceType);
        FIELD(nbStepsPerYear, "Nb of steps per year");
        FIELD_MAKE_OPTIONAL(nbStepsPerYear);
        FIELD(spotDiscreteScheme, "Spot discretization scheme (" 
                                         + VolSVJ_SpotDSTypeHelper::getNameList()
                                         + ")");
        FIELD_MAKE_OPTIONAL(spotDiscreteScheme);
        FIELD(spotDiscreteSchemeIndex, "");
        FIELD_MAKE_TRANSIENT(spotDiscreteSchemeIndex);
        FIELD(varDiscreteScheme, "Variance discretization scheme (" 
                                         + VolSVJ_VarDSTypeHelper::getNameList()
                                         + ")");
        FIELD_MAKE_OPTIONAL(varDiscreteScheme);
        FIELD(varDiscreteSchemeIndex, "");
        FIELD_MAKE_TRANSIENT(varDiscreteSchemeIndex);
        FIELD_NO_DESC(simDates);
        FIELD_MAKE_TRANSIENT(simDates); 
        FIELD(isCarefulRandoms, "Use careful randoms if true; don't otherwise");
        FIELD_MAKE_OPTIONAL(isCarefulRandoms);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMCPathConfigSVJ(){
        return new MCPathConfigSVJ();
    }

    // for reflection
    MCPathConfigSVJ():
    MCPathConfigParametric(TYPE),
    isCarefulRandoms(false),
    dependenceType("not used"),
    nbStepsPerYear(12),
    spotDiscreteScheme(VolSVJ_SpotDSTypeHelper::getDefaultName()),
    varDiscreteScheme(VolSVJ_VarDSTypeHelper::getDefaultName()),
    spotDiscreteSchemeIndex(-1),
    varDiscreteSchemeIndex(-1){}

//    class Generator;
    friend class Generator;

#ifdef STATE_VARIABLES
    // Referee class
    class Gen;
    friend class Gen;
    // PathGen for spot
    // class PathGenFineProcess;
    friend class PathGenFineProcess;
    typedef refCountPtr<PathGenFineProcess> PathGenFineProcessSP;
#endif

    // visible fields
    bool            isCarefulRandoms;
    string          dependenceType;
    int             nbStepsPerYear;
    string          spotDiscreteScheme;
    string          varDiscreteScheme;

    // transient fields
    int                spotDiscreteSchemeIndex;
    int                varDiscreteSchemeIndex;
    DateTimeArraySP    simDates; // cached between tweaks
};

CClassConstSP const MCPathConfigSVJ::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigSVJ", typeid(MCPathConfigSVJ), load);

CClassConstSP const MCPathConfigSVJ::Random::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigSVJ::Random", typeid(MCPathConfigSVJ::Random), load);

///////////////////////////////////////////////////////////////////////////////////
// Generator
///////////////////////////////////////////////////////////////////////////////////
class MCPathConfigSVJ::Generator: public MCPathBase,
                                  virtual public IMCRandom::Callbacks,
                                  public DependenceMakerGauss::Support {
public:
    Generator(int                      numSimPaths,
              MCPathConfigSVJ*         pathConfig, 
              const MCPathGeneratorSP& pastPathGenerator,
              const IMCProduct*         prod):
    pathConfig(pathConfig),
    mAsset(prod->getMultiFactors()), 
    nbAssets(prod->getNumAssets()),
    simDates(pathConfig->simDates),
    fwds(nbAssets),
    quanto(nbAssets),
    vols(nbAssets),
    timeMetrics(nbAssets),
    tradYears(nbAssets),
    productPaths(nbAssets),
    instVars(nbAssets),
    integratedVars(nbAssets),
    maxDrifts(nbAssets, 1.0),
    isCarefulRandoms(pathConfig->carefulRandoms()) {
        static const string method("MCPathConfigSVJ::Generator::Generator");
        try{
            // Obtain product timeline
            timeline = getProductTimeline(prod, pastPathGenerator);

            // Some consistency checks
            if (prod->getToday() != prod->getEffectiveSimStartDate()){
                throw ModelException(method, "fwd start style not supported");
            }
            ASSERT(timeline->futureDates[0] == prod->getEffectiveSimStartDate());

            // Obtain market data
            IntArray nbPaths(nbAssets, 1);
            refData = getRefData(timeline, nbPaths, mAsset, 
                                 prod, pastPathGenerator);

            // Initialize vols and time metrics
            int iAsset;
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                // copy the vols that are in the multiasset
                vols[iAsset] = VolSVJSP::dynamicCast(
                    VolRequestRaw::copyVolBase(*mAsset, iAsset));
                // copy time metric
                timeMetrics[iAsset] = TimeMetricSP(copy(&vols[iAsset]->getTimeMetric()));
            }

            // reuse sim dates between tweaks for numerical stability
            if (simDates->empty()){
                simDates = genSimDates(timeline->futureDates,
                                       *timeMetrics[0],       // pick a time metric !
                                       pathConfig->nbStepsPerYear,
                                       fixDateOffsets);
            }
        
            // Initialize paths, i.e. logSpot, inst var and integrated var paths
            // Also, initialize productPaths. 
            // NB: the latter contains the past (if any), the fomer don't
            logSpots.resize(timeline->futureDates.size());
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                productPaths[iAsset].resize(timeline->totalNumSteps);
                instVars[iAsset].resize(simDates->size()) ;
                integratedVars[iAsset].resize(simDates->size()) ;
            }

            // compute trading years
            DateTime startDate = (*simDates)[0];
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                tradYears[iAsset].resize(simDates->size());
                for (int iDate = 0; iDate < simDates->size(); ++iDate) {
                    tradYears[iAsset][iDate] 
                        = timeMetrics[iAsset]->yearFrac(startDate,
                                                        (*simDates)[iDate]);
                }
            }

            // pre compute fwds

            preComputeFwds();

		    // set up Dependence
            dependence = pathConfig->dependenceMaker->createDependence(this);

            // Initialize random number generator for spots
            // NB the size of the random deviates depends on the scheme 
            // that is requested
            int nbSpotDeviates;
            if (pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EULER){
                nbSpotDeviates = simDates->size() - 1;
            }
            else if(pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EXACT){
                nbSpotDeviates = logSpots.size() - 1;
            }
            else{
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
            }
            int nbVarDeviates;
            if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else{
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_VarDSTypeHelper::getName(pathConfig->varDiscreteSchemeIndex));
            }
            IRandomSP rand(pathConfig->getRandomGenerator());
            MCPathConfigSVJ::Random* mcsvjrand = dynamic_cast<MCPathConfigSVJ::Random*>(rand.get());
            if (!mcsvjrand){
                throw ModelException(method, 
                                     "the random generator must be of type MCPathConfigSVJ::Random");
            }
            randomGenSpot = MCRandNormalNoCacheSP(new MCRandNormalNoCache(
                0,
                dependence,
                mcsvjrand->normal,
                pathConfig->carefulRandoms(),
                nbSpotDeviates,
                nbAssets,
                timeline->totalNumPastDates));
            randomGenVar = MCRandNormalNoCacheSP(new MCRandNormalNoCache(
                0,
                dependence,
                mcsvjrand->normal,
                pathConfig->carefulRandoms(),
                nbVarDeviates,
                nbAssets,
                0));
            randomGenJumpNb = MCRandNormalNoCacheSP(new MCRandNormalNoCache(
                0,
                dependence,
                mcsvjrand->normal,
                pathConfig->carefulRandoms(),
                nbSpotDeviates,
                nbAssets,
                0));
            randomGenJumpSize = MCRandNormalNoCacheSP(new MCRandNormalNoCache(
                this,
                dependence,
                mcsvjrand->normal,
                pathConfig->carefulRandoms(),
                nbSpotDeviates,
                nbAssets,
                0));
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
        the drift between simulation date i and simulation date i+1.
        (With d_0 = simulation start date). Range of dates is restricted to
        those in the future so the past path generator will
        return 1.0 */
    double maxDriftProduct(int iAsset) const{
        return maxDrifts[iAsset];
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigSVJ::Generator::getGaussData");
        try {
            return mAsset->factorsCorrelationMatrix();
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

private:

    /** Configures pathGen for antithetics */
    void configureAntithetics() {
    };

    /** Configures pathGen for nonAntithetics. We need a better name */
    void configureNonAntithetics() {
    };

    /** Draws random numbers */
    virtual void drawRandomNumbers(int pathIdx){
        // Draw random numbers for paths and jump sizes
        randomGenSpot->generate(pathIdx);
        randomGenVar->generate(pathIdx);
        randomGenJumpNb->generate(pathIdx);
        randomGenJumpSize->generate(pathIdx); 
    } 
    
    virtual void generatePath(int pathIdx, int iAsset, int iPath) {
        static const string method = "MCPathConfigSVJ::Generator::generatePath";
        try{
            const DoubleMatrix& randomSpot = randomGenSpot->getRandomNumbers();
            const DoubleMatrix& randomVol = randomGenVar->getRandomNumbers();
            const DoubleMatrix& randomUniJumpMatrix = randomGenJumpNb->getRandomNumbers();
            const DoubleMatrix& randomNormalJump = randomGenJumpSize->getRandomNumbers();
            DoubleArray randomUniJump(randomUniJumpMatrix.numRows());
            for (int iUniJump = 0; iUniJump < randomUniJump.size(); ++iUniJump){
                randomUniJump[iUniJump] = N1(randomUniJumpMatrix[iAsset][iUniJump]);
            }
            // XXX this list of nested if statements is obviously not viable
            if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::EULER){
                if (pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EULER){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSVJ_SpotDSType::EULER>(),
                                                  Int2Type<VolSVJ_VarDSType::EULER>(),
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset], 
                                                  &randomUniJump[0],
                                                  randomNormalJump[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]); 
                }
                else if(pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EXACT){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSVJ_SpotDSType::EXACT>(),
                                                  Int2Type<VolSVJ_VarDSType::EULER>(), 
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  &randomUniJump[0],
                                                  randomNormalJump[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else{
                    throw ModelException(method, 
                                         "unexpected spot discretization scheme of type "
                                         + VolSVJ_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
                }
            }
            else if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
                if (pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EULER){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSVJ_SpotDSType::EULER>(),
                                                  Int2Type<VolSVJ_VarDSType::VAR_TRANSFORM_EULER>(),
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  &randomUniJump[0],
                                                  randomNormalJump[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else if(pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EXACT){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSVJ_SpotDSType::EXACT>(),
                                                  Int2Type<VolSVJ_VarDSType::VAR_TRANSFORM_EULER>(), 
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  &randomUniJump[0],
                                                  randomNormalJump[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else{ 
                    throw ModelException(method, 
                                         "unexpected spot discretization scheme of type "
                                         + VolSVJ_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
                }
            }
            else{
                throw ModelException(method, 
                                     "unexpected variance discretization scheme of type "
                                     + VolSVJ_VarDSTypeHelper::getName(pathConfig->varDiscreteSchemeIndex));
            }
            // then, populate spot path
            const DoubleArray& fwds = this->fwds[iAsset];
            DoubleArray& spots = productPaths[iAsset];
            int nbFutFixDates = logSpots.size();
            // NB: spots[timeline->numPastDates] corresponds to futureDates[1], i.e.
            // the first future product/fixing date -- not today
            
            for (int iFutFixDates = 1, iFixDates = timeline->numPastDates;
                 iFutFixDates < nbFutFixDates; ++iFutFixDates, ++iFixDates) {
                    spots[iFixDates] = fwds[iFutFixDates] * exp(logSpots[iFutFixDates]);
                    //quantoAdjustment had been taken out and put in log-spot path
	        }
        }
        catch (exception& e){ 
            throw ModelException(e, method);
        }
    }

    /** Obtains timeline object from base */
    MCProductTimelineConstSP getTimeline() const {
        return timeline;
    }

    /** Returns number of assets */
    int NbSimAssets() const {
        return nbAssets;
    }

    /** Returns the reference level for iAsset, iPath */
    const double& refLevel(int iAsset, int iPath) const {
        return refData->refLevels[iAsset][0];
    }
    
    /** Returns the reference level for iAsset, iPath */
    double& refLevel(int iAsset, int iPath) {
        return refData->refLevels[iAsset][0];
    }

    /** Returns the reflevel path */
    IRefLevel::IMCPathSP& refLevelPath() const {
        return refData->refLevelPath;
    }

    /** Returns the path for iAsset, iPath */
    const double* Path(int iAsset, int iPath) const {
        return &productPaths[iAsset][0];
    }

    /** Returns the path for iAsset, iPath */
    double* Path(int iAsset, int iPath)  {
        return &productPaths[iAsset][0];
    }

    /** Returns the number of paths per asset */
    int nbPaths(int iAsset) const {
        return 1;
    }

    /** Returns if it is a single path generator or not */
    bool isSinglePath() const {
        return true;
    }

    /** pre-compute fwds arrays, get quanto data here for now */
    // Reminder about timelines:
    // No state variables: coarse path (timeline->futureDates); fine path (simDates);
    // State variables: coarse path (simTimeline); fine path (simDates);
    void preComputeFwds(){

        static const string routine = "MCPathConfigSVJ::preComputeFwds";

        try {

            // for each asset
            for (int iAsset = 0; iAsset < nbAssets; iAsset++){
                // need to cache the forward at all simulation dates
                fwds[iAsset].resize(timeline->numFutSteps + 1);
                
                
                CAssetConstSP asset = CAssetConstSP::attachToRef(&mAsset->getAsset(iAsset));

                if (CAsset::IStruck::TYPE->isInstance(asset)) {
                    throw ModelException(routine, "Struck assets (" + asset->getTrueName() + 
                                        ") are not supported yet");
                } 

                //using quanto of type quantoParamSVJ
                bool isQuanto = Asset::IQuanto::TYPE->isInstance(asset);
                DoubleArray  fxSqrtVar = DoubleArray(simDates->size() - 1);
                double eqFXCorr = 0.0;
                if(isQuanto) {
                    const Asset::IQuanto& prot = dynamic_cast<const Asset::IQuanto&>(*asset.get());
                    
                    // Get the unprotected forward values
                    prot.unadjustedFwdValue(timeline->futureDates, fwds[iAsset]);
                    
                    // get fx corr
                    eqFXCorr = prot.getCorrelation()->getCorrelation();
                    
                    // get fx vols atm
                    ATMVolRequestSP fxVolRequest(new ATMVolRequest());
                    CVolProcessedSP volFX(prot.getProcessedFXVol( fxVolRequest.get() ));
                    // cast to the type of vol we're expecting
                    CVolProcessedBSSP volFXBS = CVolProcessedBSSP::dynamicCast(volFX);
                    volFXBS->CalcVol(*simDates, CVolProcessedBS::forward, fxSqrtVar);
                    
                    // calculate fx vars
                    for (int iStep = 0; iStep < fxSqrtVar.size(); iStep++) {
                        double dt = tradYears[iAsset][iStep+1]-tradYears[iAsset][iStep];
                        fxSqrtVar[iStep] *= dt; // using asset dt
                    }
                    quanto[iAsset] = QuantoParamSVJSP(new QuantoParamSVJ(isQuanto,eqFXCorr,fxSqrtVar));
                } else{
                    quanto[iAsset] = QuantoParamSVJSP(new QuantoParamSVJ(isQuanto,eqFXCorr,fxSqrtVar));
                    mAsset->factorFwdValues(iAsset, timeline->futureDates, fwds[iAsset]);
                }
            }
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // fields
    MCPathConfigSVJ*        pathConfig;
    const IMultiFactors*    mAsset;            // multi asset interface
    int                     nbAssets;          // number of assets
    DateTimeArraySP         simDates;          // fine grid of simulation dates
    vector<DoubleArray>     fwds;              // fwds [iAsset][iStep]
    VolSVJArray             vols;              // array of VolSVJs
    TimeMetricArray         timeMetrics;       // array of time metrics
    vector<DoubleArray>     tradYears;         // fine grid of trad years
    vector<DoubleArray>     productPaths;      // spot path [iAsset][iStep]
    DoubleArray             logSpots;          // spot path [iAsset][iStep]
    vector<DoubleArray>     instVars;          // instantaneous variance path [iAsset][iStep]
    vector<DoubleArray>     integratedVars;    // integrated variance path [iAsset][iStep]
    IntArray                fixDateOffsets;    // offsets for fixing dates
    MCRandNormalNoCacheSP   randomGenSpot;     // Random number generator for spots
    MCRandNormalNoCacheSP   randomGenVar;      // Random number generator for vars
    MCRandNormalNoCacheSP   randomGenJumpNb;   // Random number generator for jump uniform deviates
    MCRandNormalNoCacheSP   randomGenJumpSize; // Random number generator for jump normal deviates
    DependenceSP            dependence;        // Dependence object
    MCProductTimelineSP     timeline;          // Product timeline
    RefLevelDataSP          refData;           // RefLevel data
    DoubleArray             maxDrifts;         // XXX product of MAX(1, drifts)
    bool                    isCarefulRandoms;  // XXX

    // currency protection
    vector<QuantoParamSVJSP>     quanto;            // quanto parameters
};

#ifdef STATE_VARIABLES
/** Spot path generator using state variable approach */
class MCPathConfigSVJ::PathGenFineProcess: virtual public MCPathGen,
                                           virtual public IMCRandom::Callbacks,
                                           public DependenceMakerGauss::Support {
public:

    virtual ~PathGenFineProcess(){}

    PathGenFineProcess(int                             numSimPaths, 
                       MCPathConfigSVJ*                pathConfig, 
                       const PastPathGenSpotSP&        pastPathGenSpot,
                       const MCProductClient*          prodClient,
                       const SVGenSpotArray&              spotGenArray,
                       const MCSqrtAnnualQuadVarArray& sqrtAnnualQuadVarGenArray,
                       StateVarDBase&                  svDBase):
    pathConfig(pathConfig),
    mAsset(prodClient->getMultiFactors()), 
    nbAssets(prodClient->getNumAssets()),
    simDates(pathConfig->simDates), 
    fwds(nbAssets),
    quanto(nbAssets),
    vols(nbAssets),
    timeMetrics(nbAssets),
    pathGenMakers(nbAssets),
    pathGens(nbAssets),
    tradYears(nbAssets),
    instVars(nbAssets),
    quadVars(nbAssets),
    randomPathSpot(nbAssets),
    randomPathVar(nbAssets),
    maxDrifts(nbAssets, 1.0),
    simHasPast(pastPathGenSpot->hasPast()),
    spotProdPaths(nbAssets),
    spotProdOffset(nbAssets, 0),
    sqrtAnnualQuadVarProdPaths(nbAssets),
    sqrtAnnualQuadVarProdOffset(nbAssets, 0) {

        static const string method("MCPathConfigSVJ::Generator::Generator");
        try{
            const DateTime& today = prodClient->getToday();

            // Create simulation timeline
            DateTimeArraySP spotGenDates = MCPath::getAllDates(spotGenArray);
            DateTimeArraySP sqrtAnnualQuadVarGenDates = MCPath::getAllDates(sqrtAnnualQuadVarGenArray);
            DateTimeArray allDates = DateTime::merge(*sqrtAnnualQuadVarGenDates, *spotGenDates);
            simTimeline = today.getFutureDates(allDates);
            numFutSteps = simTimeline.size();
            if (!numFutSteps){
                throw ModelException(method,
                                     "there is no future date to simulate");
            }
            simTimeline.insert(simTimeline.begin(), today);

            // Obtain market data
            IntArray nbPaths(nbAssets, 1);

            // Initialize vols and time metrics
            int iAsset;
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                // copy the vols that are in the multiasset
                vols[iAsset] = VolSVJSP::dynamicCast(
                    VolRequestRaw::copyVolBase(*mAsset, iAsset));
                // copy time metric
                timeMetrics[iAsset] = TimeMetricSP(copy(&vols[iAsset]->getTimeMetric()));
                // create the pathgen makers
                pathGenMakers[iAsset] = VolSVJ::IPathGenMakerSP(
                    pathConfig->createPathGenMaker());
                pathGens[iAsset] = VolSVJ::IPathGenSP(
                    vols[iAsset]->createPathGen(*pathGenMakers[iAsset]));
            }

            // reuse sim dates between tweaks for numerical stability
            if (simDates->empty()){
                simDates = genSimDates(simTimeline,
                                       *timeMetrics[0],       // pick a time metric !
                                       pathConfig->nbStepsPerYear,
                                       fixDateOffsets);
            }
        
            // Initialize paths, i.e. logSpot, inst var and integrated var paths
            logSpots.resize(simTimeline.size());
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                instVars[iAsset].resize(simDates->size()) ;
                quadVars[iAsset].resize(simDates->size()) ;
            }

            if (spotGenArray.size()){
                // Create spot paths per asset
                DateTimeArrayArray spotDatesPerAsset(nbAssets);
                vector<const double*> spotPtrs(nbAssets);
                vector<int> spotBeginInd(nbAssets);
                vector<int> spotEndInd(nbAssets);

                const IPastValues* pastValues = prodClient->getMCPastValues();
                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    // Populate spot past values
                    DateTimeArraySP spotProdDates = MCPath::getAllDates(spotGenArray, iAsset);
                    spotProdPaths[iAsset] = DoubleArray(spotProdDates->size());
                    DateTimeArray spotPastProdDates = today.getPastDates(*spotProdDates);
                    DoubleArray spotPastValues = pastValues->getPastValues(spotPastProdDates, iAsset, today);
                    spotProdOffset[iAsset] = spotPastValues.size();
                    int iStep;
                    for(iStep = 0; iStep < spotPastValues.size(); iStep++) {
                        spotProdPaths[iAsset][iStep] = spotPastValues[iStep];
                    }
                
                    // Create spot mappings
                    spotDatesPerAsset[iAsset] = *spotProdDates;
                    spotPtrs[iAsset] = &spotProdPaths[iAsset][0];
                    spotBeginInd[iAsset] = spotPastProdDates.size();
                    spotEndInd[iAsset]   = spotProdDates->size();
                }

                // Create SVGenSpot::IStateVars and put them in database
                MCPath::IStateVarArray spotSVJArray(MCPath::createPaths(
                    false,
                    spotGenArray,
                    spotDatesPerAsset,
                    spotBeginInd,
                    spotEndInd,
                    spotPtrs,
                    maxDrifts));

                for(unsigned int iVar = 0; iVar < spotGenArray.size(); iVar++) {
                    svDBase.append(spotGenArray[iVar], spotSVJArray[iVar]);
                }
            }

            if (sqrtAnnualQuadVarGenArray.size()){
                // Create quad var paths per asset
                DateTimeArrayArray sqrtAnnualQuadVarDatesPerAsset(nbAssets);
                vector<const double*> sqrtAnnualQuadVarPtrs(nbAssets);
                vector<int> sqrtAnnualQuadVarBeginInd(nbAssets);
                vector<int> sqrtAnnualQuadVarEndInd(nbAssets);

                // some historical values for the past
                pastAnnualQuadVars.resize(nbAssets);
                pastWeights.resize(nbAssets);

                const IPastValues* pastValues = prodClient->getMCSqrtAnnualQuadVarPastValues();
                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    // populate quad var past values
                    DateTimeArraySP sqrtAnnualQuadVarProdDates = MCPath::getAllDates(sqrtAnnualQuadVarGenArray, iAsset);
                    sqrtAnnualQuadVarProdPaths[iAsset] = DoubleArray(sqrtAnnualQuadVarProdDates->size());
                    DateTimeArray sqrtAnnualQuadVarPastProdDates = today.getPastDates(*sqrtAnnualQuadVarProdDates);
                    DoubleArray sqrtAnnualQuadVarPastValues 
                        = pastValues->getPastValues(sqrtAnnualQuadVarPastProdDates, iAsset, today);
                    sqrtAnnualQuadVarProdOffset[iAsset] = sqrtAnnualQuadVarPastValues.size();
                    int iStep;
                    for(iStep = 0; iStep < sqrtAnnualQuadVarPastValues.size(); iStep++) {
                        sqrtAnnualQuadVarProdPaths[iAsset][iStep] = sqrtAnnualQuadVarPastValues[iStep];
                    }

                    // some precomputation need to be done to evaluate total variance
                    // when there is a past
                    if (sqrtAnnualQuadVarProdOffset[iAsset] > 0){
                        // get the realized var at valueDate
                        DateTimeArray valueDate(1, (*simDates)[0]);
                        DoubleArray pastVols = pastValues->getPastValues(valueDate, iAsset, valueDate[0]);
                        pastAnnualQuadVars[iAsset] = Maths::square(pastVols[0]);
                        // get nb of past values up to and including valueDate                    
                        double nbPastValues = pastValues->getNbPastValues(valueDate[0], iAsset);
                        // get total nb of past values
                        double nbValues = pastValues->getNbPastValues(iAsset);
                        pastWeights[iAsset] = (nbPastValues - 1.0) / (nbValues - 1.0);
                    }

                    // Create sqrtAnnualQuadVar mappings
                    sqrtAnnualQuadVarDatesPerAsset[iAsset] = *sqrtAnnualQuadVarProdDates;
                    sqrtAnnualQuadVarPtrs[iAsset] = &sqrtAnnualQuadVarProdPaths[iAsset][0];
                    sqrtAnnualQuadVarBeginInd[iAsset] = sqrtAnnualQuadVarPastProdDates.size();
                    sqrtAnnualQuadVarEndInd[iAsset]   = sqrtAnnualQuadVarProdDates->size();
                }

                // Create MCQuadVar::IStateVars and put them in database
                MCPath::IStateVarArray sqrtAnnualQuadVarSVJArray(MCPath::createPaths(
                    false,
                    sqrtAnnualQuadVarGenArray,
                    sqrtAnnualQuadVarDatesPerAsset,
                    sqrtAnnualQuadVarBeginInd,
                    sqrtAnnualQuadVarEndInd,
                    sqrtAnnualQuadVarPtrs,
                    maxDrifts));

                for(unsigned int iVar = 0; iVar < sqrtAnnualQuadVarGenArray.size(); iVar++) {
                    svDBase.append(sqrtAnnualQuadVarGenArray[iVar], sqrtAnnualQuadVarSVJArray[iVar]);
                }
            }

            // compute trading years
            DateTime startDate = (*simDates)[0];
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                tradYears[iAsset].resize(simDates->size());
                for (int iDate = 0; iDate < simDates->size(); ++iDate) {
                    tradYears[iAsset][iDate] 
                        = timeMetrics[iAsset]->yearFrac(startDate,
                                                        (*simDates)[iDate]);
                }
            }

            // pre compute fwds
            preComputeFwds();

            // set up Dependence
            dependence = pathConfig->dependenceMaker->createDependence(this);

            // Initialize random number generator for spots
            // NB the size of the random deviates depends on the scheme 
            // that is requested
            int nbSpotDeviates;
            if (pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EULER){
                nbSpotDeviates = simDates->size() - 1;
            }
            else if(pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EXACT){
                nbSpotDeviates = logSpots.size() - 1;
            }
            else{
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
            }
            int nbVarDeviates;
            if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else{
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_VarDSTypeHelper::getName(pathConfig->varDiscreteSchemeIndex));
            }
            IRandomSP rand(pathConfig->getRandomGenerator());
            MCPathConfigSVJ::Random* mcsvjrand = dynamic_cast<MCPathConfigSVJ::Random*>(rand.get());
            if (!mcsvjrand){
                throw ModelException(method, 
                                     "the random generator must be of type MCPathConfigSVJ::Random");
            }
            randomGenSpot = MCRandNormalNoCacheSP(new MCRandNormalNoCache(
                0,
                dependence,
                mcsvjrand->normal,
                pathConfig->carefulRandoms(),
                nbSpotDeviates,
                nbAssets,
                spotGenDates->size() - numFutSteps));
            randomGenVar = MCRandNormalNoCacheSP(new MCRandNormalNoCache(
                0,
                dependence,
                mcsvjrand->normal,
                pathConfig->carefulRandoms(),
                nbVarDeviates,
                nbAssets,
                0));
            randomGenJumpNb = MCRandPoissonSP(new MCRandPoisson(
                mcsvjrand->poisson));
            randomGenJumpSize = MCRandNormalSP(new MCRandNormal(
                mcsvjrand->normal));
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                randomPathSpot[iAsset] = IMCRandNormal::PathSP(randomGenSpot->getPath(iAsset));
                randomPathVar[iAsset] = IMCRandNormal::PathSP(randomGenVar->getPath(iAsset));
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
        the drift between simulation date i and simulation date i+1.
        (With d_0 = simulation start date). Range of dates is restricted to
        those in the future so the past path generator will
        return 1.0 */
    double maxDriftProduct(int iAsset) const{
        return maxDrifts[iAsset];
    }

    bool hasPast() const {
        return simHasPast;
    }

    /** MCPathGen method */
    void generatePath(int pathIdx) {
        // Draw random numbers
        drawRandomNumbers(pathIdx);
        // loop over assets
        for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
            generatePath(pathIdx, iAsset);
        }
    }

    bool doingPast() const {
        return false;
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigSVJ::PathGenFineProcess::getGaussData");
        try {
            return mAsset->factorsCorrelationMatrix();
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

private:

    /** Configures pathGen for antithetics */
    void configureAntithetics() {
    };

    /** Configures pathGen for nonAntithetics. We need a better name */
    void configureNonAntithetics() {
    };

    /** Draws random numbers */
    virtual void drawRandomNumbers(int pathIdx){
        // Draw random numbers for paths and jump sizes
        randomGenSpot->generate(pathIdx);
        randomGenVar->generate(pathIdx);
//        randomGenJumpNb->generate(pathIdx);
//        randomGenJumpSize->generate(pathIdx);
    }
    
    virtual void generatePath(int pathIdx, int iAsset) {
        static const string method = "MCPathConfigSVJ::Generator::generatePath";
        try{
            // generate the "core" path
            pathGens[iAsset]->generatePath(tradYears[iAsset],
                                           fixDateOffsets,
                                           *randomPathSpot[iAsset],
                                           *randomPathVar[iAsset], 
                                           *randomGenJumpNb,
                                           *randomGenJumpSize,
                                           quanto[iAsset],
                                           logSpots,
                                           instVars[iAsset],
                                           quadVars[iAsset]);
            // then, populate spot path
            const DoubleArray& fwds = this->fwds[iAsset];
            int nbFutFixDates = logSpots.size();
            DoubleArray& spots = spotProdPaths[iAsset];
            if (spots.size()){
                for (int iFutFixDate = 1, iSpotProdDate = spotProdOffset[iAsset];
                     iFutFixDate < nbFutFixDates; ++iFutFixDate, ++iSpotProdDate) {
                        spots[iSpotProdDate] = fwds[iFutFixDate] * exp(logSpots[iFutFixDate]);
	            }
            }
            // and sqrt annualized quad var path
            DoubleArray& sqrtAnnualQuadVars = sqrtAnnualQuadVarProdPaths[iAsset];
            if (sqrtAnnualQuadVars.size()){
                bool hasPast = sqrtAnnualQuadVarProdOffset[iAsset] > 0;
                for (int iFutFixDate = 1, iSqrtAnnualQuadVarProdDate = sqrtAnnualQuadVarProdOffset[iAsset],
                     iIntegratedVar = fixDateOffsets[0];
                     iFutFixDate < nbFutFixDates; 
                     ++iFutFixDate,
                     ++iSqrtAnnualQuadVarProdDate) {
                    if (hasPast){
                        // future annualized var from today till current sim date
						double tradYear = tradYears[iAsset][iIntegratedVar];
                        double futAnnualQuadVar = Maths::isZero(tradYear) ? 
                                                  instVars[iAsset][iIntegratedVar] :
												  quadVars[iAsset][iIntegratedVar] / tradYear;
                        // total annualized var from first sampling date till current sim date
                        double annualQuadVar = pastWeights[iAsset] * pastAnnualQuadVars[iAsset]
                                               + (1.0 - pastWeights[iAsset]) * futAnnualQuadVar;
                        // total annualized vol
                        sqrtAnnualQuadVars[iSqrtAnnualQuadVarProdDate] 
                            = sqrt(annualQuadVar);
                    }
                    else{
                        int integratedVarStartIdx = fixDateOffsets[0];
                        if (iIntegratedVar == integratedVarStartIdx){
                            sqrtAnnualQuadVars[iSqrtAnnualQuadVarProdDate] = 0.0;
                        }
                        else{
                            // quad var increment between first sampling date and current sim date
                            double quadVarInc = quadVars[iAsset][iIntegratedVar]
                                                - quadVars[iAsset][integratedVarStartIdx];
                            // time elapsed since first sampling date
                            double timeDiff = tradYears[iAsset][iIntegratedVar]
                                              - tradYears[iAsset][integratedVarStartIdx];
                            // total annualized vol
                            sqrtAnnualQuadVars[iSqrtAnnualQuadVarProdDate] 
                                = sqrt(quadVarInc / timeDiff);
                        }
                    }
                    if (iFutFixDate < fixDateOffsets.size()){
                        iIntegratedVar += fixDateOffsets[iFutFixDate];
                    }
	            }
            }
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** pre-compute fwds arrays, get quanto data here for now */
    // Reminder about timelines:
    // No state variables: coarse path (timeline->futureDates); fine path (simDates);
    // State variables: coarse path (simTimeline); fine path (simDates);
    void preComputeFwds() {

        static const string routine = "MCPathConfigSVJ::preComputeFwds";

        try { 
            // for each asset
            for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
                // need to cache the forward at all simulation dates
                fwds[iAsset].resize(simTimeline.size());

                CAssetConstSP asset = CAssetConstSP::attachToRef(&mAsset->getAsset(iAsset));

                if (CAsset::IStruck::TYPE->isInstance(asset)) {
                    throw ModelException(routine, "Struck assets (" + asset->getTrueName() + 
                                        ") are not supported yet");
                } 

                //using quanto of type quantoParamSVJ
                bool isQuanto = Asset::IQuanto::TYPE->isInstance(asset);
                double eqFXCorr = 0.0;
                DoubleArray fxSqrtVar = DoubleArray(simDates->size() - 1);
          
                if(isQuanto){
                    const Asset::IQuanto& prot = dynamic_cast<const Asset::IQuanto&>(*asset.get());
                    
                    // Get the unprotected forward values
                    prot.unadjustedFwdValue(simTimeline, fwds[iAsset]);
                    
                    // get fx corr
                    eqFXCorr = prot.getCorrelation()->getCorrelation();
                    
                    // get fx vols atm
                    ATMVolRequestSP fxVolRequest(new ATMVolRequest());
                    CVolProcessedSP volFX(prot.getProcessedFXVol( fxVolRequest.get() ));
                    // cast to the type of vol we're expecting
                    CVolProcessedBSSP volFXBS = CVolProcessedBSSP::dynamicCast(volFX);
                    volFXBS->CalcVol(*simDates, CVolProcessedBS::forward, fxSqrtVar);
                    
                    // calculate fx vars
                    for (int iStep = 0; iStep < fxSqrtVar.size(); iStep++) {
                        double dt = tradYears[iAsset][iStep+1]-tradYears[iAsset][iStep];
                        fxSqrtVar[iStep] *= dt; // using asset dt
                    }
                    quanto[iAsset] = QuantoParamSVJSP(new QuantoParamSVJ(isQuanto,eqFXCorr,fxSqrtVar));
		        }
                else{
                    quanto[iAsset] = QuantoParamSVJSP(new QuantoParamSVJ(isQuanto,eqFXCorr,fxSqrtVar)); 
                    mAsset->factorFwdValues(iAsset, simTimeline, fwds[iAsset]);
                }
            }
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // fields
    MCPathConfigSVJ*        pathConfig;
    const IMultiFactors*    mAsset;            // multi asset interface
    int                     nbAssets;          // number of assets
    DateTimeArraySP         simDates;          // fine grid of simulation dates
    vector<DoubleArray>     fwds;              // fwds [iAsset][iStep]
    VolSVJArray             vols;              // array of VolSVJs
    TimeMetricArray         timeMetrics;       // array of time metrics
    VolSVJ::IPathGenMakerArray pathGenMakers;  // array of path gen makers         
    VolSVJ::IPathGenArray   pathGens;          // array of path gens
    vector<DoubleArray>     tradYears;         // fine grid of trad years
    DoubleArray             logSpots;          // spot path [iAsset][iStep]
    vector<DoubleArray>     instVars;          // instantaneous variance path [iAsset][iStep]
    vector<DoubleArray>     quadVars;          // quadratic variation path [iAsset][iStep]
    IntArray                fixDateOffsets;    // offsets for fixing dates
    MCRandNormalNoCacheSP   randomGenSpot;     // Random number generator for spots
    IMCRandNormal::PathArray randomPathSpot;   // Random path for spots  
    MCRandNormalNoCacheSP   randomGenVar;      // Random number generator for vars
    IMCRandNormal::PathArray randomPathVar;    // Random path for vars
    MCRandPoissonSP         randomGenJumpNb;   // Random number generator for jump uniform deviates
    MCRandNormalSP          randomGenJumpSize; // Random number generator for jump normal deviates
    DependenceSP            dependence;        // Dependence object
    vector<double>          pastAnnualQuadVars;// past realized vars at value date [iAsset]
    vector<double>          pastWeights;       // past weights at value date [iAsset]
    vector<double>          maxDrifts;         // XXX product of MAX(1, drifts)
    
    vector<QuantoParamSVJSP>     quanto; //quanto parameters

    // Timeline fields
    DateTimeArray           simTimeline;     //!< Today + strictly future merged dates
    int                     numFutSteps;     //!< Number of strictly future merged dates
    bool                    simHasPast;      //!< Whether product has past

    vector<DoubleArray>     spotProdPaths;   // spot path [iAsset][iStep]
    vector<int>             spotProdOffset;  //!< Starting point for future path per asset

    vector<DoubleArray>     sqrtAnnualQuadVarProdPaths;      // sqrtAnnualQuadVar path [iAsset][iStep]
    vector<int>             sqrtAnnualQuadVarProdOffset;  //!< Starting point for future path per asset
};

/** Referee class that distributes simulation to components
    e.g. spot, hitTime etc. */
class MCPathConfigSVJ::Gen: virtual public MCPathGenerator, // For backward compatibility
                            virtual public MCPathGen,
                            virtual public IStateVariableGen::IStateGen {
public:
    /** Constructor */
    Gen(int                      numSimPaths,
        MCPathConfigSVJ*         pathConfig, 
        const MCPathGeneratorSP& pastPathGenerator,
        const MCProductClient*   prodClient):
    nowPathIdx(0) {
        static const string routine = "MCPathConfigSVJ::Gen::Gen";
        
        try {
            // Collect state variables from product and categorize them
            StateVariableCollectorSP svCollector(new StateVariableCollector());
            prodClient->collectStateVars(svCollector);
            IElemStateVariableGenArray stateVarGenArray = svCollector->getElemStateVarGens();
            
            // Spot requests
            SVGenSpotArray spotGenArray = filterStateVars<SVGenSpot>(stateVarGenArray);
            
            // SqrtAnnualQuadVar requests
            MCSqrtAnnualQuadVarArray sqrtAnnualQuadVarGenArray = filterStateVars<MCSqrtAnnualQuadVar>(stateVarGenArray);

            // Create discount factor state variables
            SVGenDiscFactorArray discFactors = filterStateVars<SVGenDiscFactor>(stateVarGenArray);
            for(unsigned int iVar = 0; iVar < discFactors.size(); iVar++) {
                IStateVariableSP sv(discFactors[iVar]->determinsticSV(false /* not doing past */));
                svDBase.append(discFactors[iVar], sv);
            }
            // Nothing else is supported
            if(stateVarGenArray.size()) {
                throw ModelException("Unable to recognize all state variable types.");
            }
            
            // XXX should pass in a PastPathGen in the first place, really
            PastPathGen* pastPathGen = dynamic_cast<PastPathGen*>(pastPathGenerator.get());
            if(!pastPathGen) {
                throw ModelException("Past path generator is not of PastPathGen type.");
            }

            // Create a (spot, quad var, sqrt annualized quad var) path generator
            pathGenFineProcess = MCPathConfigSVJ::PathGenFineProcessSP(
                new MCPathConfigSVJ::PathGenFineProcess(numSimPaths, 
                                                       pathConfig, 
                                                       pastPathGen->getPathGenSpot(), 
                                                       prodClient, 
                                                       spotGenArray,
                                                       sqrtAnnualQuadVarGenArray,
                                                       svDBase));
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** MCpathGenerator methods */
    // Deprecated methods
    int NbSimAssets() const {
        static const string routine = "MCPathConfigSVJ::Gen::NbSimAssets";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    const double* Path(int iAsset, int iPath) const {
        static const string routine = "MCPathConfigSVJ::Gen::Path";
        throw ModelException(routine, "Method is retired for StateVars");
    }; 

    double refLevel(int iAsset, int iPath) const {
        static const string routine = "MCPathConfigSVJ::Gen::refLevel";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    double maxDriftProduct(int iAsset) const {
        static const string routine = "MCPathConfigSVJ::Gen::maxDriftProduct";
        try {
            return pathGenFineProcess->maxDriftProduct(iAsset);
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    int begin(int iAsset) const {
        static const string routine = "MCPathConfigSVJ::Gen::begin";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    int end(int iAsset) const{
        static const string routine = "MCPathConfigSVJ::Gen::end";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    // Live methods
    bool hasPast() const {
        return pathGenFineProcess->hasPast();
    }

    bool doingPast() const {
        return false;
    }

    void generatePath(int pathIdx) {
        nowPathIdx = pathIdx;
        pathGenFineProcess->generatePath(pathIdx);
    }

    int getPathIndex() const {
        return nowPathIdx;
    }

    /** Returns the state variable corresponding to generator.
        Part of the IStateVariableGen::IStateGen IFace */
    IStateVariableSP create(const IStateVariableGen* svGen) {
        static const string routine = "MCPathConfigSVJ::Gen::create";

        try {
            return svDBase.find(svGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }
  
private:
    StateVarDBase                         svDBase;               //!< Collection of Generators + statevars
    MCPathConfigSVJ::PathGenFineProcessSP pathGenFineProcess;    //!< Spot generator
    int                                   nowPathIdx;            //!< Current path idx
};
#endif

/** Creates a future path generator */
MCPathGeneratorSP MCPathConfigSVJ::makePathGenerator(
        bool                               cachingRequested,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     mainSimDates ){
    static const string routine = "MCPathConfigSVJ::makePathGenerator";

    if(cachingRequested) {
        throw ModelException(routine, "Paths caching is not supported in MCPathConfigSVJ");
    }

    // create empty sim dates array if simDates is null or on new pricing run
    if (!simDates || control->isPricing()){
        simDates = DateTimeArraySP(new DateTimeArray());
    }

    const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
    if(!prodClient){
        refCountPtr<MCPathConfigSVJ::Generator> futurePathGenBase(
            new MCPathConfigSVJ::Generator(numPaths,
                                          this,
                                          pastPathGenerator,
                                          prod));
    
        // store num of sim dates
        if (control && control->isPricing()){
            int nbSimDates = simDates->size();
            results->storeScalarGreek(nbSimDates - 1, 
                                      Results::DEBUG_PACKET, 
                                      OutputNameSP(
                                          new OutputName("MertonLV_SIM_STEPS_USED")));
        }

        // Construct the MCPathGenerator
        return MCPathBase::createPathGenerator(pastPathGenerator,
                                               this,
                                               prod,
                                               futurePathGenBase);
    } else {
        // State variables approach
        return MCPathGeneratorSP(new MCPathConfigSVJ::Gen(
            numPaths,
            this, 
            pastPathGenerator,
            prodClient));
    }
}

/** Creates a past path generator */
MCPathGeneratorSP MCPathConfigSVJ::pastPathGenerator(const IMCProduct* prod) {
    static const string method = "MCPathConfigSVJ::pastPathGenerator";
    try{
        const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
        if(!prodClient){
            // Old approach
            return MCPathConfig::pastPathGenerator(prod);
        } else {
            // State variables approach
            return MCPathGeneratorSP(new PastPathGen(prodClient));
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** MC LV model */
class MonteCarloSVJDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;
    
private:
    MonteCarloSVJDefault():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloSVJDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloSVJ);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloSVJ(){
        return new MonteCarloSVJDefault();
    }
};

CClassConstSP const MonteCarloSVJDefault::TYPE = CClass::registerClassLoadMethod(
    "MonteCarloSVJDefault", typeid(MonteCarloSVJDefault), load);

// external symbol to allow class to be forced to be linked in
bool MCPathConfigSVJLoad(){
    return (MCPathConfigSVJ::TYPE != 0 && MonteCarloSVJDefault::TYPE !=0);
}

DRLIB_END_NAMESPACE
