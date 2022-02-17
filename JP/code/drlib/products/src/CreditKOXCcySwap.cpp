//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditKOXCcySwap.cpp
//
//   Description : Credit knockout cross-currency swap
//
//   Author      : Bruno O Melka
//
//   Date        : 23 Feb 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/SVGenIRSwap.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

/** Credit KO XCcy Swap intrument */
class CreditKOXCcySwap: public CInstrument,
                        public virtual LastSensDate,
                        public virtual Theta::Shift,
                        public virtual IMCIntoProduct {
public:
    static CClassConstSP const TYPE;

    /** validates object */
    virtual void validatePop2Object() {
        static const string method = "CreditKOXCcySwap::validatePop2Object";
    }

    /** parses knockInOrOut (KO = true, KI = false) */
    bool getKnockInOrOutBool() const {
        static const string method = "CreditKOXCcySwap::getKnockInOrOutBool";
        bool knockInOrOutBool;

        if (CString::equalsIgnoreCase(knockInOrOut, "O") ||
            CString::equalsIgnoreCase(knockInOrOut, "Out") ||
            CString::equalsIgnoreCase(knockInOrOut, "KO")) {
                knockInOrOutBool = true;
        } else if (CString::equalsIgnoreCase(knockInOrOut, "I") ||          
            CString::equalsIgnoreCase(knockInOrOut, "In") ||
            CString::equalsIgnoreCase(knockInOrOut, "KI")) {
                knockInOrOutBool = false;
        } else {
            throw ModelException(method, "Unsupported knockInOrOut (" 
                + knockInOrOut + "). Valid values are I and O.");
        }

        return knockInOrOutBool;
    }

    /** performs extra validation, after GetMarket is called */
    virtual void Validate() {
        static const string method = "CreditKOXCcySwap::Validate";
        try {
            if (payDatesDom[payDatesDom.size() - 1].isLess(resetDatesDom[0])) {
                throw ModelException(method,
                    "swapStartDate ("+resetDatesDom[0].toString()
                    + ") is > or = to swapMatDate ("
                    + payDatesDom[payDatesDom.size() - 1].toString() + ")");
            }
            if (payDatesFor[payDatesFor.size() - 1].isLess(resetDatesFor[0])) {
                throw ModelException(method,
                    "swapStartDate ("+resetDatesFor[0].toString()
                    + ") is > or = to swapMatDate ("
                    + payDatesFor[payDatesFor.size() - 1].toString() + ")");
            }

            // check knockInOrOut
            getKnockInOrOutBool(); 
        }    
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** maturity date */
    DateTime maturityDate() const {
        DateTime mat = payDatesDom.back().max(payDatesFor.back());
        return mat;
    }

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const {
        return maturityDate();
    }

    /** the value date, implementing CInstrument */
    virtual DateTime getValueDate() const {
        return valueDate;
    }

    /** roll through time, implementing Theta::Shift */
    virtual bool sensShift(Theta* theta) {
        valueDate = theta->rollDate(valueDate);
        return true;
    }
    
    /** gets market data */
    virtual void GetMarket(const IModel* model, const CMarketDataSP market)
    {
        market->GetReferenceDate(valueDate);
        discountDom.getData(model, market);
        discountFor.getData(model, market);

        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, discountDom, fx);
        ICDSParSpreads::getMarketData(model, market.get(), discountDom.getName(), credit);

        couponDom = YieldCurveSP(discountDom.getSP().clone());
        couponDom->setProjectionCurve(); // use 'growth' zc

        couponFor = YieldCurveSP(discountFor.getSP().clone());
        couponFor->setProjectionCurve(); // use 'growth' zc
    }

    /** returns the name of the instrument's discount currency */
    virtual string discountYieldCurveName() const {
        return discountDom.getName();
    }

    /** creates product, implementing IMCIntoProduct
        see below for implementation*/
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;

private:

    friend class CreditKOXCcySwapMC;

    /** invoked when class is 'loaded' (for reflection) */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditKOXCcySwap, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultCreditKOXCcySwap);

        FIELD(accStartDatesDom, "domestic accrual start dates");
        FIELD(accEndDatesDom, "domestic accrual end dates");
        FIELD(resetDatesDom, "domestic reset dates");
        FIELD(payDatesDom, "domestic payment dates");
        FIELD(dccDom, "domestic day count convention");
        FIELD(indexFreqDom, "domestic index frequency");
        FIELD(indexMaturityDom, "domestic index maturity");
        FIELD(indexDCCDom, "domestic index day count convention");
        FIELD(indexStubTypeDom, "domestic index stub type");
        FIELD(indexStubAtEndDom, "only matters if indexStubType != N");
        FIELD(indexAccBadDayConvDom, "domestic index accrual bad day convention");
        FIELD(indexPayBadDayConvDom, "domestic index payment bad day convention");
        FIELD(indexHolsDom, "domestic index holidays");
        FIELD(weightsDom, "domestic weights");
        FIELD(spreadsDom, "domestic spreads");
        FIELD(capsDom, "domestic caps");
        FIELD(floorsDom, "domestic floors");
        FIELD(notionalsDom, "domestic notionals");
        FIELD(recoveryNotionalsDom, "domestic notionals for recovery");
        FIELD(payInitialNotionalDom, "TRUE = pay initial notional on domestic leg");
        FIELD(payFinalNotionalDom, "TRUE = pay final notional on domestic leg");
        FIELD(payAccOnDefaultDom, "TRUE = pay accrued interest on default on domestic leg");
        FIELD(recRateOnDefaultDom, "recovery rate on domestic notional on default");
        FIELD(useRecRateOnDefaultDom, "FALSE = use par recovery rate");
        FIELD(discountDom, "domestic discount curve");
        FIELD_NO_DESC(couponDom);
        FIELD_MAKE_TRANSIENT(couponDom);

        FIELD(accStartDatesFor, "foreign accrual start dates");
        FIELD(accEndDatesFor, "foreign accrual end dates");
        FIELD(resetDatesFor, "foreign reset dates");
        FIELD(payDatesFor, "foreign payment dates");
        FIELD(dccFor, "foreign day count convention");
        FIELD(indexFreqFor, "foreign index frequency");
        FIELD(indexMaturityFor, "foreign index maturity");
        FIELD(indexDCCFor, "foreign index day count convention");
        FIELD(indexStubTypeFor, "foreign index stub type");
        FIELD(indexStubAtEndFor, "only matters if indexStubType != N");
        FIELD(indexAccBadDayConvFor, "foreign index accrual bad day convention");
        FIELD(indexPayBadDayConvFor, "foreign index payment bad day convention");
        FIELD(indexHolsFor, "foreign index holidays");
        FIELD(weightsFor, "foreign weights");
        FIELD(spreadsFor, "foreign spreads");
        FIELD(capsFor, "foreign caps");
        FIELD(floorsFor, "foreign floors");
        FIELD(notionalsFor, "foreign notionals");
        FIELD(recoveryNotionalsFor, "foreign notionals for recovery");
        FIELD(payInitialNotionalFor, "TRUE = pay initial notional on foreign leg");
        FIELD(payFinalNotionalFor, "TRUE = pay final notional on foreign leg");
        FIELD(payAccOnDefaultFor, "TRUE = pay accrued interest on default on foreign leg");
        FIELD(recRateOnDefaultFor, "recovery rate on foreign notional on default");
        FIELD(useRecRateOnDefaultFor, "FALSE = use par recovery rate");
        FIELD(discountFor, "foreign discount curve");
        FIELD_NO_DESC(couponFor);
        FIELD_MAKE_TRANSIENT(couponFor);

        FIELD(valueDate, "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(receiveDomLeg, "TRUE = receive domestic leg and pay foreign leg");
        FIELD(knockInOrOut, "I = knock-in, O = knock-out");
        FIELD(ccyTreatment, "currency Treatment");
        FIELD(fx, "fx");
        FIELD(credit, "credit");
    }

    /** default constructor (for reflection) */
    static IObject* defaultCreditKOXCcySwap() {
        return new CreditKOXCcySwap();
    }

    /** constructor */
    CreditKOXCcySwap():CInstrument(TYPE) {}; 

    /** copy constructor (not implemented) */
    CreditKOXCcySwap(const CreditKOXCcySwap& rhs);

    /** assignment (not implemented) */
    CreditKOXCcySwap& operator=(const CreditKOXCcySwap& rhs);

    // domestic leg fields
    DateTimeArray             accStartDatesDom;          // accrual start dates
    DateTimeArray             accEndDatesDom;            // accrual end dates
    DateTimeArray             resetDatesDom;             // reset dates
    DateTimeArray             payDatesDom;               // payment dates
    DayCountConventionConstSP dccDom;                    // day count convention
    MaturityPeriodConstSP     indexFreqDom;              // index frequency
    MaturityPeriodConstSP     indexMaturityDom;          // index maturity
    DayCountConventionConstSP indexDCCDom;               // index day count convention
    string                    indexStubTypeDom;          // index stub type
    bool                      indexStubAtEndDom;         // only matters if indexStubType != NONE
    BadDayConventionConstSP   indexAccBadDayConvDom;     // index accrual bad day convention
    BadDayConventionConstSP   indexPayBadDayConvDom;     // index payment bad day convention
    HolidaySP                 indexHolsDom;              // index holidays
    DoubleArray               weightsDom;                // weights
    DoubleArray               spreadsDom;                // spreads
    DoubleArray               capsDom;                   // caps
    DoubleArray               floorsDom;                 // floors
    DoubleArray               notionalsDom;              // notionals
    DoubleArray               recoveryNotionalsDom;      // notionals to be used for recovery
    bool                      payInitialNotionalDom;     // TRUE = pay initial notional
    bool                      payFinalNotionalDom;       // TRUE = pay final notional
    bool                      payAccOnDefaultDom;        // TRUE = pay accrued interest on default
    double                    recRateOnDefaultDom;       // recovery rate on notional on default
    bool                      useRecRateOnDefaultDom;    // FALSE = use par recovery rate
    YieldCurveWrapper         discountDom;               // discount curve
    YieldCurveSP              couponDom;                 // transient - set to 'projection' curve
    
    // foreign leg fields
    DateTimeArray             accStartDatesFor;          // accrual start dates
    DateTimeArray             accEndDatesFor;            // accrual end dates
    DateTimeArray             resetDatesFor;             // reset dates
    DateTimeArray             payDatesFor;               // payment dates
    DayCountConventionConstSP dccFor;                    // day count convention
    MaturityPeriodConstSP     indexFreqFor;              // index frequency
    MaturityPeriodConstSP     indexMaturityFor;          // index maturity
    DayCountConventionConstSP indexDCCFor;               // index day count convention
    string                    indexStubTypeFor;          // index stub type
    bool                      indexStubAtEndFor;         // only matters if indexStubType != NONE
    BadDayConventionConstSP   indexAccBadDayConvFor;     // index accrual bad day convention
    BadDayConventionConstSP   indexPayBadDayConvFor;     // index payment bad day convention
    HolidaySP                 indexHolsFor;              // index holidays
    DoubleArray               weightsFor;                // weights
    DoubleArray               spreadsFor;                // spreads
    DoubleArray               capsFor;                   // caps
    DoubleArray               floorsFor;                 // floors
    DoubleArray               notionalsFor;              // notionals
    DoubleArray               recoveryNotionalsFor;      // notionals to be used for recovery
    bool                      payInitialNotionalFor;     // TRUE = pay initial notional
    bool                      payFinalNotionalFor;       // TRUE = pay final notional
    bool                      payAccOnDefaultFor;        // TRUE = pay accrued interest on default
    double                    recRateOnDefaultFor;       // recovery rate on notional on default
    bool                      useRecRateOnDefaultFor;    // FALSE = use par recovery rate
    YieldCurveWrapper         discountFor;               // discount curve
    YieldCurveSP              couponFor;                 // transient - set to 'projection' curve

    // other instrument fields
    DateTime                  valueDate;
    bool                      receiveDomLeg;             // TRUE = receive domestic leg and pay foreign leg
    string                    knockInOrOut;              // I = knock-in, O = knock-out
    string                    ccyTreatment;              // none(vanilla), protected or struck
    CAssetWrapper             fx;                        // fx
    ICDSParSpreadsWrapper     credit;                    // credit
};

CClassConstSP const CreditKOXCcySwap::TYPE = CClass::registerClassLoadMethod(
    "CreditKOXCcySwap", typeid(CreditKOXCcySwap), CreditKOXCcySwap::load);

/** MC product class for CreditKOXCcySwap */
class CreditKOXCcySwapMC : public MCProductClient
{

private:

    // state variables
    SVDiscFactorSP                 discFacPayDatesSVDom;  // disc facs on (all) domestic pay dates  (numCoups+1)
    SVSurvivalDiscFactorSP         survProbPayDatesSVDom; // surv probs on (all) domestic pay dates (numCoups+1)
    SVGenIRSwap::StateVarArraySP   indexSVDom;            // domestic indices                       (numCoups)
    SVDiscFactorSP                 discFacPayDatesSVFor;  // disc facs on (all) foreign pay dates   (numCoups+1)
    SVSurvivalDiscFactorSP         survProbPayDatesSVFor; // surv probs on (all) foreign pay dates  (numCoups+1)
    SVGenIRSwap::StateVarArraySP   indexSVFor;            // foreign indices                        (numCoups)
    SVGenSpot::IStateVarSP         fxPayDatesSV;          // fx rate on (all) foreign pay dates     (numCoups+1)

    // state variable generators
    SVGenDiscFactorSP              discFacPayDatesSVGenDom;
    SVGenSurvivalDiscFactorSP      survProbPayDatesSVGenDom;
    SVGenIRSwapArraySP		       indexSVGenDom;
    SVGenDiscFactorSP              discFacPayDatesSVGenFor;
    SVGenSurvivalDiscFactorSP      survProbPayDatesSVGenFor;
    SVGenIRSwapArraySP		       indexSVGenFor;
    SVGenSpotSP                    fxPayDatesSVGen;

    // from instrument
    DoubleArray                    weightsDom;
    DoubleArray                    spreadsDom;
    DoubleArray                    capsDom;
    DoubleArray                    floorsDom;
    DoubleArray                    notionalsDom;
    DoubleArray                    recoveryNotionalsDom;
    bool                           payInitialNotionalDom;
    bool                           payFinalNotionalDom;
    bool                           payAccOnDefaultDom;
    double                         recRateOnDefaultDom;
    bool                           useRecRateOnDefaultDom;
    DoubleArray                    weightsFor;
    DoubleArray                    spreadsFor;
    DoubleArray                    capsFor;
    DoubleArray                    floorsFor;
    DoubleArray                    notionalsFor;
    DoubleArray                    recoveryNotionalsFor;
    bool                           payInitialNotionalFor;
    bool                           payFinalNotionalFor;
    bool                           payAccOnDefaultFor;
    double                         recRateOnDefaultFor;
    bool                           useRecRateOnDefaultFor;
    bool                           receiveDomLeg;
  
    // other
    int                            numCoupsDom;             // number of domestic coupon periods
    DoubleArray                    yearFracsDom;            // domestic day count fracs
    double                         actRecRateOnDefaultDom;  // actual domestic rec rate to use
    bool                           payRecOnDefaultDom;      // TRUE if recovery is non-zero
    int                            numCoupsFor;             // number of foreign coupon periods
    DoubleArray                    yearFracsFor;            // foreign day count fracs
    double                         actRecRateOnDefaultFor;  // actual foreign rec rate to use
    bool                           payRecOnDefaultFor;      // TRUE if recovery is non-zero
    bool                           knockInOrOutBool;        // true = KO, false = KI

protected:

    /** update state variables */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) {
        discFacPayDatesSVDom = discFacPayDatesSVGenDom->getSVDiscFactor(newPathGen);
        survProbPayDatesSVDom = survProbPayDatesSVGenDom->getSVSurvivalDiscFactor(newPathGen);
        for (int i = 0; i < numCoupsDom; i++) {
            (*indexSVDom)[i] = SVGenIRSwap::IStateVarSP(((*indexSVGenDom)[i])->getIRSwapSV(newPathGen));
        }
        discFacPayDatesSVFor = discFacPayDatesSVGenFor->getSVDiscFactor(newPathGen);
        survProbPayDatesSVFor = survProbPayDatesSVGenFor->getSVSurvivalDiscFactor(newPathGen);
        for (int i = 0; i < numCoupsFor; i++) {
            (*indexSVFor)[i] = SVGenIRSwap::IStateVarSP(((*indexSVGenFor)[i])->getIRSwapSV(newPathGen));
        }
        fxPayDatesSV = fxPayDatesSVGen->getSpotSV(newPathGen);
    }

public:

    /** requests state variables */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const {
        svCollector->append(discFacPayDatesSVGenDom.get());
        svCollector->append(survProbPayDatesSVGenDom.get());
        for (int i = 0; i < numCoupsDom; i++) {
            svCollector->append((*indexSVGenDom)[i].get());
        }
        svCollector->append(discFacPayDatesSVGenFor.get());
        svCollector->append(survProbPayDatesSVGenFor.get());
        for (int i = 0; i < numCoupsFor; i++) {
            svCollector->append((*indexSVGenFor)[i].get());
        }
        svCollector->append(fxPayDatesSVGen.get()); 
    }

    /** constructor */
    CreditKOXCcySwapMC (
        const CreditKOXCcySwap* inst,
        const SimSeriesSP& simSeries,
        InstrumentSettlementSP instSettle):
        MCProductClient(
            inst->fx.get(),
            inst->valueDate,
            inst->discountDom.get(),
            IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)), // fix!
            simSeries, // fix!
            IPastValuesSP(IPastValues::Util::makeTrivial(inst->valueDate, 0.0)), // fix
            instSettle.get(),
            inst->maturityDate()),
        weightsDom(inst->weightsDom),
        spreadsDom(inst->spreadsDom),
        capsDom(inst->capsDom),
        floorsDom(inst->floorsDom),
        notionalsDom(inst->notionalsDom),
        recoveryNotionalsDom(inst->notionalsDom),
        payInitialNotionalDom(inst->payInitialNotionalDom),
        payFinalNotionalDom(inst->payFinalNotionalDom),
        payAccOnDefaultDom(inst->payAccOnDefaultDom),
        recRateOnDefaultDom(inst->recRateOnDefaultDom),
        useRecRateOnDefaultDom(inst->useRecRateOnDefaultDom),
        weightsFor(inst->weightsFor),
        spreadsFor(inst->spreadsFor),
        capsFor(inst->capsFor),
        floorsFor(inst->floorsFor),
        notionalsFor(inst->notionalsFor),
        recoveryNotionalsFor(inst->notionalsFor),
        payInitialNotionalFor(inst->payInitialNotionalFor),
        payFinalNotionalFor(inst->payFinalNotionalFor),
        payAccOnDefaultFor(inst->payAccOnDefaultFor),
        recRateOnDefaultFor(inst->recRateOnDefaultFor),
        useRecRateOnDefaultFor(inst->useRecRateOnDefaultFor),
        receiveDomLeg(inst->receiveDomLeg)
    {
        numCoupsDom = inst->notionalsDom.size();

        // gather all domestic pay dates i.e. initial notional date & coupon dates
        DateTime initialNotionalDateDom = inst->accStartDatesDom[0];
        DateTimeArray allPayDatesDom = inst->payDatesDom;
        allPayDatesDom.insert(allPayDatesDom.begin(), initialNotionalDateDom);

        discFacPayDatesSVGenDom = SVGenDiscFactorSP(new SVGenDiscFactor(
            inst->valueDate,
            inst->discountDom.getSP(),
            allPayDatesDom));

        survProbPayDatesSVGenDom = SVGenSurvivalDiscFactorSP(new SVGenSurvivalDiscFactor(
            inst->valueDate,
            inst->credit.getSP(),
            allPayDatesDom));

        indexSVDom = SVGenIRSwap::StateVarArraySP(new SVGenIRSwap::StateVarArray(numCoupsDom));
        indexSVGenDom = SVGenIRSwapArraySP(new SVGenIRSwapArray(numCoupsDom));
        for (int i = 0; i < numCoupsDom; i++) {
            (*indexSVGenDom)[i] = SVGenIRSwapSP(new SVGenIRSwap(
                inst->couponDom,
                inst->discountDom.getSP(),
                inst->resetDatesDom[i],
                inst->indexAccBadDayConvDom->adjust(inst->indexMaturityDom->toDate(inst->resetDatesDom[i]), inst->indexHolsDom.get()),
                inst->indexFreqDom->toString(),
                inst->indexDCCDom->toString(),
                inst->indexStubTypeDom,
                inst->indexStubAtEndDom,
                inst->indexAccBadDayConvDom->toString(),
                inst->indexPayBadDayConvDom->toString(),
                inst->indexHolsDom,
                true)); // not used?
        }

        yearFracsDom.resize(numCoupsDom);
        for (int i = 0; i < numCoupsDom; i++) {
            yearFracsDom[i] = inst->dccDom->years(inst->accStartDatesDom[i], inst->accEndDatesDom[i]);
        }

        actRecRateOnDefaultDom = useRecRateOnDefaultDom ? recRateOnDefaultDom : inst->credit->getRecovery();
        payRecOnDefaultDom = Maths::isPositive(actRecRateOnDefaultDom);

        numCoupsFor = inst->notionalsFor.size();

        // gather all foreign pay dates i.e initial notional date & coupon dates
        DateTime initialNotionalDateFor = inst->accStartDatesFor[0];
        DateTimeArray allPayDatesFor = inst->payDatesFor;
        allPayDatesFor.insert(allPayDatesFor.begin(), initialNotionalDateFor);

        discFacPayDatesSVGenFor = SVGenDiscFactorSP(new SVGenDiscFactor(
            inst->valueDate,
            inst->discountFor.getSP(),
            allPayDatesFor));

        survProbPayDatesSVGenFor = SVGenSurvivalDiscFactorSP(new SVGenSurvivalDiscFactor(
            inst->valueDate,
            inst->credit.getSP(),
            allPayDatesFor));

        indexSVFor = SVGenIRSwap::StateVarArraySP(new SVGenIRSwap::StateVarArray(numCoupsFor));
        indexSVGenFor = SVGenIRSwapArraySP(new SVGenIRSwapArray(numCoupsFor));
        for (int i = 0; i < numCoupsFor; i++) {
            (*indexSVGenFor)[i] = SVGenIRSwapSP(new SVGenIRSwap(
                inst->couponFor,
                inst->discountFor.getSP(),
                inst->resetDatesFor[i],
                inst->indexAccBadDayConvFor->adjust(inst->indexMaturityFor->toDate(inst->resetDatesFor[i]), inst->indexHolsFor.get()),
                inst->indexFreqFor->toString(),
                inst->indexDCCFor->toString(),
                inst->indexStubTypeFor,
                inst->indexStubAtEndFor,
                inst->indexAccBadDayConvFor->toString(),
                inst->indexPayBadDayConvFor->toString(),
                inst->indexHolsFor,
                true)); // not used?
        }

        fxPayDatesSVGen = SVGenSpotSP(new SVGenSpot(1, allPayDatesFor));

        yearFracsFor.resize(numCoupsFor);
        for (int i = 0; i < numCoupsFor; i++) {
            yearFracsFor[i] = inst->dccFor->years(inst->accStartDatesFor[i], inst->accEndDatesFor[i]);
        }

        actRecRateOnDefaultFor = useRecRateOnDefaultFor ? recRateOnDefaultFor : inst->credit->getRecovery();
        payRecOnDefaultFor = Maths::isPositive(actRecRateOnDefaultFor);

        knockInOrOutBool = inst->getKnockInOrOutBool();
    }

    /** returns payoff (called within the simulation loop) */
    virtual void payoff(const IPathGenerator* pathGen, IMCPrices& prices) {

        // DOMESTIC LEG

        // disc facs
        const SVPath& pathDiscFacPayDatesDom = discFacPayDatesSVDom->path();
        
        // surv probs
        const SVPath& pathSurvProbPayDatesDom = survProbPayDatesSVDom->path();

        // index rates
        double parYieldDom, annuityDom;
        DoubleArraySP indexRatesDom(new DoubleArray(numCoupsDom));
        for (int i = 0; i < numCoupsDom; i++) {
            (*indexSVDom)[i].get()->parYield(parYieldDom, annuityDom);
            (*indexRatesDom)[i] = parYieldDom;
        }

        double payoffULDom = 0.0; // payoff of underlying leg
        double payoffKODom = 0.0; // payoff of knock-out leg
        double payoffKIDom = 0.0; // payoff of knock-in leg

        // initial notional
        if (payInitialNotionalDom) {
            double initialNotionalFlowDom = - notionalsDom[0];
            double initialNotionalRiskFreePayoffDom = initialNotionalFlowDom * pathDiscFacPayDatesDom[0];
            payoffULDom += initialNotionalRiskFreePayoffDom;
            payoffKODom += initialNotionalRiskFreePayoffDom * pathSurvProbPayDatesDom[0];
        }

        for (int j = 0; j < numCoupsDom; j++) {

            // final notionals
            if (payFinalNotionalDom) {
                double finalNotionalFlowDom = ( j < numCoupsDom - 1 ?  notionalsDom[j] - notionalsDom[j+1] : notionalsDom[j] );
                double finalNotionalRiskFreePayoffDom = finalNotionalFlowDom * pathDiscFacPayDatesDom[j+1];
                payoffULDom += finalNotionalRiskFreePayoffDom;
                payoffKODom += finalNotionalRiskFreePayoffDom * pathSurvProbPayDatesDom[j+1];
            }

            // recovery on default
            if (payRecOnDefaultDom) {
                double defaultProb = pathSurvProbPayDatesDom[j] - pathSurvProbPayDatesDom[j+1];
                payoffKODom += actRecRateOnDefaultDom * recoveryNotionalsDom[j] * pathDiscFacPayDatesDom[j+1] * defaultProb;
            }

            // coupons
            double couponRate = Maths::collar((*indexRatesDom)[j] * weightsDom[j] + spreadsDom[j], capsDom[j], floorsDom[j]);
            double couponFlowDom = notionalsDom[j] * couponRate * yearFracsDom[j];
            double couponRiskFreePayoffDom = couponFlowDom * pathDiscFacPayDatesDom[j+1];
            payoffULDom += couponRiskFreePayoffDom;
            payoffKODom += couponRiskFreePayoffDom * pathSurvProbPayDatesDom[j+1];

            // accrued interest on default
            if (payAccOnDefaultDom) {
                double defaultProb = pathSurvProbPayDatesDom[j] - pathSurvProbPayDatesDom[j+1];
                payoffKODom += 0.5 * notionalsDom[j] * couponRate * yearFracsDom[j] * pathDiscFacPayDatesDom[j+1] * defaultProb;
            }
        }

        payoffKIDom = payoffULDom - payoffKODom;

        double payoffDom = (knockInOrOutBool) ? payoffKODom : payoffKIDom;

        // FOREIGN LEG

        // disc facs
        const SVPath& pathDiscFacPayDatesFor = discFacPayDatesSVFor->path();

        // surv probs
        const SVPath& pathSurvProbPayDatesFor = survProbPayDatesSVFor->path();

        // index rates
        double parYieldFor, annuityFor;
        DoubleArraySP indexRatesFor(new DoubleArray(numCoupsFor));
        for (int i = 0; i < numCoupsFor; i++) {
            (*indexSVFor)[i].get()->parYield(parYieldFor, annuityFor);
            (*indexRatesFor)[i] = parYieldFor;
        }

        // fx rates
        const SVPath& pathFxPayDates = fxPayDatesSV->path(0);

        double payoffULFor = 0.0; // payoff of underlying leg
        double payoffKOFor = 0.0; // payoff of knock-out leg
        double payoffKIFor = 0.0; // payoff of knock-in leg

        // initial notional
        if (payInitialNotionalFor) {
            double initialNotionalFlowFor = - pathFxPayDates[0] *notionalsFor[0];
            double initialNotionalRiskFreePayoffFor = initialNotionalFlowFor * pathDiscFacPayDatesFor[0];
            payoffULFor += initialNotionalRiskFreePayoffFor;
            payoffKOFor += initialNotionalRiskFreePayoffFor * pathSurvProbPayDatesFor[0];
        }

        for (int j = 0; j < numCoupsFor; j++) {

            // final notionals
            if (payFinalNotionalFor) {
                double finalNotionalFlowFor = pathFxPayDates[j+1] * ( j < numCoupsFor - 1 ?  notionalsFor[j] - notionalsFor[j+1] : notionalsFor[j] );
                double finalNotionalRiskFreePayoffFor = finalNotionalFlowFor * pathDiscFacPayDatesFor[j+1];
                payoffULFor += finalNotionalRiskFreePayoffFor;
                payoffKOFor += finalNotionalRiskFreePayoffFor * pathSurvProbPayDatesFor[j+1];
            }

            // recovery on default
            if (payRecOnDefaultFor) {
                double defaultProb = pathSurvProbPayDatesFor[j] - pathSurvProbPayDatesFor[j+1];
                payoffKOFor += actRecRateOnDefaultFor * pathFxPayDates[j+1] * recoveryNotionalsFor[j] * pathDiscFacPayDatesFor[j+1] * defaultProb;
            }

            // coupons
            double couponRate = Maths::collar((*indexRatesFor)[j] * weightsFor[j] + spreadsFor[j], capsFor[j], floorsFor[j]);
            double couponFlowFor = pathFxPayDates[j+1] * notionalsFor[j] * couponRate * yearFracsFor[j];
            double couponRiskFreePayoffFor = couponFlowFor * pathDiscFacPayDatesFor[j+1];
            payoffULFor += couponRiskFreePayoffFor;
            payoffKOFor += couponRiskFreePayoffFor * pathSurvProbPayDatesFor[j+1];

            // accrued interest on default
            if (payAccOnDefaultFor) {
                double defaultProb = pathSurvProbPayDatesFor[j] - pathSurvProbPayDatesFor[j+1];
                payoffKOFor += 0.5 * pathFxPayDates[j+1] * notionalsFor[j] * couponRate * yearFracsFor[j] * pathDiscFacPayDatesFor[j+1] * defaultProb;
            }
        }

        payoffKIFor = payoffULFor - payoffKOFor;

        double payoffFor = (knockInOrOutBool) ? payoffKOFor : payoffKIFor;

        // SWAP

        double payoff = payoffDom - payoffFor;

        if (!receiveDomLeg) {
            payoff *= -1.0;
        }

        prices.add(payoff);
    }

    /** returns pv factor, overriding default method (pv factor already applied in payoff) */
    virtual double pvFromPaymentDate() const {
        return 1.0;
    }

};

/** creates product, implementing IMCIntoProduct interface */
IMCProduct* CreditKOXCcySwap::createProduct(const MonteCarlo* model) const {
    // the simSeries passed to the IMCProduct is redundant
    // well not quite. The CreditKOXCcySwapMC wants to know the last sim date
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(DateTimeArray(1, maturityDate()));
    InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
    return new CreditKOXCcySwapMC(this, simSeries, instSettle);
}

/** for class loading */
bool CreditKOXCcySwapLoad() {
    return (CreditKOXCcySwap::TYPE != 0);
}

//declare smartPtr versions
DECLARE(CreditKOXCcySwap);

DRLIB_END_NAMESPACE
