/**
 * @file CCMDBasketBetaDSpread.cpp
 *
 *
 */

#include "edginc/config.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/GenericAllNameScalarShift.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/ParSpreadPropShift.hpp"
#include "edginc/CCMAbsoluteBetaTweak.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Sensitivity to absolute change in all CCM betas
 */

class CCMAbsoluteBetaSens:
    public GenericAllNameScalarShift<CCMAbsoluteBetaTweak> {

public:

    CCMAbsoluteBetaSens(double spreadScaleShift = 0.01):
        GenericAllNameScalarShift<CCMAbsoluteBetaTweak>(
            TYPE, "CCM_ABSOLUTE_BETA_SENS", spreadScaleShift) {}

    /**
     * Reflection type of this class
     */

    static CClassConstSP const TYPE;

    /**
     * Output name for first derivative
     */

    const static string NAME;

    /**
     * Shift size to use if none provided
     */

    const static double DEFAULT_SHIFT;

private:

    /**
     * Invoked when Class is 'loaded'
     */

    static void load(CClassSP& clazz){
        REGISTER(CCMAbsoluteBetaSens, clazz);
        SUPERCLASS(GenericAllNameScalarShift<CCMAbsoluteBetaTweak>);
    }
};

CClassConstSP const CCMAbsoluteBetaSens::TYPE =
CClass::registerClassLoadMethod("CCMAbsoluteBetaSens",
                                typeid(CCMAbsoluteBetaSens),
                                CCMAbsoluteBetaSens::load);

template<> CClassConstSP const GenericAllNameScalarShift<CCMAbsoluteBetaTweak>::TYPE =
    CClass::registerClassLoadMethod("GenericAllNameScalarShift<CCMAbsoluteBetaTweak>",
                                    typeid(GenericAllNameScalarShift<CCMAbsoluteBetaTweak>),
                                    GenericAllNameScalarShift<CCMAbsoluteBetaTweak>::load);

/**
 * Sensitivity of CCM "whole-basket beta-sensitivity" to spreads
 *
 * See CCMDBetaDSpread in CCMDBetaDSpread.cpp: this class calculates
 * CCM_D_BASKET_BETA_D_PAR_SPREAD_RHO_PROP.
 */

class CCMDBasketBetaDSpread: public Sensitivity {

public:

    static CClassConstSP const TYPE;
    static const string NAME;

    // Just for GenericSensitivityFactory:
    static const double DEFAULT_SHIFT;

    /**
     * @param spreadScaleShift  The relative amount to tweak the spread curves
     *                          (via ParSpreadPropShift)
     * @param betaShift         The amount to tweak the betas
     *                          (via CCMAbsoluteBetaSens)
     */

    CCMDBasketBetaDSpread(double spreadScaleShift = DEFAULT_SHIFT,
                          double betaShift = 0.01):
        Sensitivity(TYPE),
        spreadScaleShift(spreadScaleShift),
        betaShift(betaShift)
    {}

    const string& getPacketName() const {
        return Results::INSTRUMENT_PACKET;
    }

    /** identifies the name used for storing associated results in the output*/

    const string& getSensOutputName() const {
        return NAME;
    }

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one (return value: false) */

    bool discreteShift() const {
        return false;
    }

    /** calculates sensitivity - invoked by calculateSens */

    void calculate(TweakGroup*  tweakGroup,
                   CResults*    results) {
        string packetName = getPacketName();
        OutputNameConstSP outputName(new OutputName(getSensOutputName()));
        try {
            if (results->exists(packetName, outputName)) {
                // Done this before: don't overwrite anything
                return;
            }

            CCMAbsoluteBetaSens vByBeta(getBetaShift());
            OutputNameConstSP vByBetaOutputName(
                new OutputName(vByBeta.getSensOutputName()));

            Results baseResultsTmp;
            Results *baseResults;
            if (getControl()->sensitivityRequested(
                                  vByBeta.getClass()).get() != NULL) {
                baseResults = results;
            }
            else {
                baseResultsTmp.storePrice(results->retrievePrice(),
                                        results->getCcyName());
                baseResults = &baseResultsTmp;
            }
            
            vByBeta.calculateSens(tweakGroup->getModel(),
                                   tweakGroup->getInstrument(),
                                   getControl(), baseResults);

            if (baseResults->isNotApplicable(&vByBeta)) {
                results->storeNotApplicable(this);
                return;
            }

            // DON'T use retrieveScalarGreek(..., vByBeta)
            double vByBetaBase =
                baseResults->retrieveScalarGreek(vByBeta.getPacketName(),
                                                 vByBetaOutputName);
            double sprDownShift = getSpreadScaleShift();
            if (Maths::isZero(sprDownShift)){
                throw ModelException("CCMDBasketBetaDSpread::calculate",
                                     "spreadScaleShift is zero");
            }
            ParSpreadPropShift sprDown(sprDownShift);
            SensMgrOpt sensMgr(tweakGroup);

            if (sensMgr.allNames(&sprDown)->empty()) {
                results->storeNotApplicable(this);
                return;
            }

            sensMgr.shift(&sprDown);

            Results downResults;
            try {
                tweakGroup->getModel()->Price(tweakGroup->getInstrument(),
                                              getControl(), &downResults);
                vByBeta.calculateSens(tweakGroup->getModel(),
                                       tweakGroup->getInstrument(),
                                       getControl(), &downResults);

                sensMgr.restore();
            }
            catch (exception&) {
                sensMgr.restore();
                throw;
            }

            if (downResults.isNotApplicable(&vByBeta)) {
                results->storeNotApplicable(this);
                return;
            }

            // DON'T use retrieveScalarGreek(..., vByBeta)
            double vByBetaSprDown =
                downResults.retrieveScalarGreek(vByBeta.getPacketName(),
                                                vByBetaOutputName);

            // DON'T use storeScalarGreek(..., this)
            results->storeScalarGreek(
                (vByBetaSprDown - vByBetaBase) / sprDownShift,
                packetName, outputName);
        }
        catch (exception& e) {
            // DON'T use storeGreek(..., this)
            results->storeGreek(IObjectSP(new Untweakable(e)),
                                packetName, outputName);
        }
    }

    /**
     * The relative amount the spread curves are tweaked (via ParSpreadPropShift)
     */

    double getSpreadScaleShift() const {
        return spreadScaleShift;
    }

    /**
     * The amount the betas are tweaked (via CCMAbsoluteBetaSens)
     */

    double getBetaShift() const {
        return betaShift;
    }

private:

    static IObject *defaultCCMDBasketBetaDSpread() {
        return new CCMDBasketBetaDSpread();
    }

    /** Invoked when this class is 'loaded' */

    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CCMDBasketBetaDSpread, clazz);
        SUPERCLASS(Sensitivity);
        FIELD(spreadScaleShift, "Relative amount to shift spreads by");
        FIELD(betaShift, "Amount to shift beta by");
        FIELD_MAKE_OPTIONAL(betaShift);
        EMPTY_SHELL_METHOD(defaultCCMDBasketBetaDSpread);

        SensitivityFactory::addSens(
            NAME, new GenericSensitivityFactory<CCMDBasketBetaDSpread>(),
            new CCMDBasketBetaDSpread(),
            CCMDBasketBetaDSpread::TYPE);
    }

    double spreadScaleShift;
    double betaShift;
};

const string CCMDBasketBetaDSpread::NAME = "CCM_D_BASKET_BETA_D_PAR_SPREAD_RHO_PROP";

const double CCMDBasketBetaDSpread::DEFAULT_SHIFT = 0.01;

CClassConstSP const CCMDBasketBetaDSpread::TYPE = CClass::registerClassLoadMethod(
    "CCMDBasketBetaDSpread", typeid(CCMDBasketBetaDSpread), CCMDBasketBetaDSpread::load);

/**
 * Included in RiskMgrLib::linkInClasses() to force CCMDBasketBetaDSpreadLinkIn
 * to get linked into the Windows exe.
 */

bool CCMDBasketBetaDSpreadLinkIn() {
    return CCMDBasketBetaDSpread::TYPE != NULL &&
           CCMAbsoluteBetaSens::TYPE != NULL &&
           GenericAllNameScalarShift<CCMAbsoluteBetaTweak>::TYPE != NULL;
}

DRLIB_END_NAMESPACE

