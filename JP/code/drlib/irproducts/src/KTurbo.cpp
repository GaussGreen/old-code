//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KTurb.cpp
//
//   Description : turbo leg component
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KTurbo.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/RateTree.hpp"

DRLIB_BEGIN_NAMESPACE

class KTurboTree : public FDProduct {
public:
    /******************** variables ********************/
    KTurboConstSP inst;
    int           nbCoupons;
    ZeroBondProdSP zeroProd;
    FDProductArray legFXprods;
    FDProductArray legISprods; // legs IndexSpec products
    CouponSchedDates& sched;


    /******************** methods ********************/
    KTurboTree(const KTurboConstSP &inst, FDModel* model);

    /** if this is elementary product - model supplied state varibles */
    virtual void init(Control*) const;
    virtual void initProd(void);
    virtual void update(int & step, FDProduct::UpdateType);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual const TreeSlice & getCashFlow(int step) const { return *cashFlowSlice; }

    virtual string getOutputName(void) const { return inst->outputName; }

    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
        string ccy;
        if (inst->discount.get()) ccy = inst->discount->getCcy();
        recordSliceToOutputName(ctrl, results, model,
            inst->isTopComponent(), inst->outputName, ccy, getValue(0, model->getDate(0)));
        inst->recordExtraOutput(ctrl, results);
    }
private:
    TreeSliceSP mainSlice;
    TreeSliceSP cashFlowSlice;
    TreeSliceSP getValueSlice;
};

/********************************** KTurboTree *******************************/

KTurboTree::KTurboTree(const KTurboConstSP &inst, FDModel* model) :
    FDProduct(model), inst(inst), nbCoupons(inst->sched->nbCoupons), sched(*inst->sched)
{
    try {
        // Create the underlying KIndexSpecs products
        legFXprods.resize(inst->legs.size());
        legISprods.resize(inst->legs.size());
        for (int legIdx=0; legIdx < inst->legs.size(); ++legIdx) {
            KTurbo::Leg &leg = *inst->legs[legIdx];
            // FX IndexSpecs
            legFXprods[legIdx] = model->createProduct(leg.fx);
            legFXprods[legIdx]->addModelResetDates(sched.resetEff);
            // Fixing indexes
            if (leg.indexes.size()) {
                legISprods[legIdx] = model->createProduct(leg.indexes[0]);
            /* [0]: should handle more than one in the future */
                DateTimeArray fixingDates;
                FlexDates::filterOnWeight(sched.resetEff, *leg.weights, fixingDates);
                legISprods[legIdx]->addModelResetDates(fixingDates);
            }
        }
        zeroProd = DYNAMIC_POINTER_CAST<ZeroBondProd>(
            model->createProduct(IProdCreatorSP(new ZeroBond(
            sched.resetEff, sched.pay, inst->discountYieldCurveName()
            ))));
    }
    catch (exception& e) {
        throw ModelException(e, "KTurboTree::KTurboTree");
    }
}

void KTurboTree::init(Control*) const{
    try {
        model->addCritDates(sched.resetEff);
        model->addCritDates(sched.pay);
    }
    catch (exception& e) { throw ModelException(e, "KTurboTree::init"); }
}

void KTurboTree::initProd(void) {
    try {
        mainSlice = model->createSlice(inst->discount->getName());
        mainSlice->name = inst->outputName + "_mainSlice";
        *mainSlice = 0.;

        cashFlowSlice = model->createSlice(inst->discount->getName());
        cashFlowSlice->name = inst->outputName + "_cashFlowSlice";
        *cashFlowSlice = 0.;

        getValueSlice = model->createSlice();
        getValueSlice->name = inst->outputName + "_getValueSlice";
    }
    catch (exception& e) { throw ModelException(e, "KTurboTree::initProd"); }
}

const TreeSlice & KTurboTree::getValue(int step, DateTime eventDate) const {
try
{
    DateTime stepDate = model->getDate(step);
    UnadjustedConvention unadjustedConvention;
    ObservationExact observationExact;

    if (eventDate != stepDate) {
        throw ModelException("Cannot be valued at a date != currentDate");
    }

    *getValueSlice = *mainSlice;
    if (step) 
        return *getValueSlice;

    for (int cpnI=0; cpnI<nbCoupons; ++cpnI) {
        DateTime payDate = sched.pay[cpnI];
        DateTime resetDate = sched.resetEff[cpnI];
        double legRate;
        double legCpn;
        double coupon = 0;

        if ((payDate > stepDate) &&  (resetDate <= model->getValueDate()))
        {
            for (int legI=0; legI < inst->legs.size(); ++legI) {
                KTurbo::Leg &leg = *(inst->legs[legI]);
                double notional = (*leg.notionals)[cpnI];
                double weight   = (*leg.weights)[cpnI];
                double dcf      = leg.dcfs[cpnI];

                if (Maths::isZero(weight)) legRate = 0.;
                else {
                    IMarketObservable *mobs = dynamic_cast<IMarketObservable*>(leg.indexes[0].get());
                    if (!mobs) throw ModelException("leg["+Format::toString(legI)+"].index does not implement IMarketObservable");

                    legRate = weight * mobs->pastValue(resetDate, &observationExact,
                        IMarketObservable::getDefaultObsSource().get(),
                        0, 0, &unadjustedConvention);
                }

                legRate += (*leg.spreads)[cpnI];

                if (leg.rateType==RateType::SIMPLE) {
                    legCpn = dcf * legRate;
                } else {
                    legCpn = ::pow (1. + legRate, dcf) - 1.;
                }

                legCpn *= notional * leg.fx->pastValue(resetDate, &observationExact,
                    IMarketObservable::getDefaultObsSource().get(),
                    0, 0, &unadjustedConvention);

                coupon += legCpn;
            }
            if (inst->notional.get()) {
                if (inst->floor.get()) {
                    double floorVal = (*inst->floor)[cpnI] * ::abs((*inst->notional)[cpnI]) * inst->dcfs[cpnI];
                    coupon = ::max(floorVal, coupon);
                }
                if (inst->cap.get()) {
                    double capValue = (*inst->cap)[cpnI] * ::abs((*inst->notional)[cpnI]) * inst->dcfs[cpnI];
                    coupon = ::min(coupon, capValue);
                }
            }
            *getValueSlice += coupon * zeroProd->getValue(step, sched.pay[cpnI]);
        }
    }
    return *getValueSlice;
}
catch (exception& e) { throw ModelException(e, "KTurboTree::getValue"); }
}



void KTurboTree::update(int & step, FDProduct::UpdateType) {
    try {
        DateTime stepDate = model->getDate(step);
        int cpnI; // coupon index
        int legI; // leg index
        *cashFlowSlice = 0.;

        for (cpnI=0; cpnI<nbCoupons; ++cpnI) {
            if (sched.resetEff[cpnI] != stepDate) continue;

            TreeSliceSP legRateSlice = cashFlowSlice->clone(false); 
            legRateSlice->name = inst->outputName + "_legRateSlice";

            TreeSliceSP legCpnSlice = cashFlowSlice->clone(false); 
            legCpnSlice->name = inst->outputName + "_legCpnSlice";

            TreeSliceSP couponSlice = cashFlowSlice->clone(false);
            couponSlice->name = inst->outputName + "_couponSlice";
            *couponSlice = 0;

            for (legI=0; legI < inst->legs.size(); ++legI) {
                KTurbo::Leg &leg = *(inst->legs[legI]);
                double notional = (*leg.notionals)[cpnI];
                double weight   = (*leg.weights)[cpnI];
                double dcf      = leg.dcfs[cpnI];

                if (Maths::isZero(weight)) *legRateSlice = 0.;
                else {
                    *legRateSlice = legISprods[legI]->getValue(step, stepDate);
                    *legRateSlice *= weight;
                }

                *legRateSlice += (*leg.spreads)[cpnI];

				if (leg.rateType==RateType::SIMPLE) *legCpnSlice = dcf * (*legRateSlice);
                else                                    *legCpnSlice = pow (1. + (*legRateSlice), dcf) - 1.;
                *legCpnSlice *= notional * legFXprods[legI]->getValue(step, stepDate);
                *couponSlice += *legCpnSlice;
            }
            if (inst->notional.get()) {
                if (inst->floor.get()) {
                    double floorSlice = (*inst->floor)[cpnI] * ::abs((*inst->notional)[cpnI]) * inst->dcfs[cpnI];
                    *couponSlice = smax(floorSlice, *couponSlice);
                }
                if (inst->cap.get()) {
                    double capSlice = (*inst->cap)[cpnI] * ::abs((*inst->notional)[cpnI]) * inst->dcfs[cpnI];
                    *couponSlice = smin(*couponSlice, capSlice);
                }
            }
            if (sched.pay[cpnI] != stepDate) {
                *legCpnSlice = zeroProd->getValue(step, sched.pay[cpnI]);
                *couponSlice *= *legCpnSlice;
            }
            *cashFlowSlice += *couponSlice;
            startDEV(mainSlice);
        }
        // payInitial (principal)
        if (stepDate==sched.accStart[0]) {
            for (legI=0; legI < inst->legs.size(); ++legI) {
                if (inst->legs[legI]->payInitial) {
                    *cashFlowSlice -= (*inst->legs[legI]->notionals)[0] * legFXprods[legI]->getValue(step, stepDate);
                    startDEV(mainSlice);
                }
            }
        }
        // payPrincipal
        for (cpnI=0; cpnI<nbCoupons; ++cpnI) {
            if (sched.pay[cpnI] != stepDate) continue;
            for (legI=0; legI < inst->legs.size(); ++legI) {
                if (inst->legs[legI]->payPrincipal) {
                    *cashFlowSlice += (*inst->legs[legI]->notionals)[cpnI] - (*inst->legs[legI]->notionals)[cpnI+1];
                    startDEV(mainSlice);
                }
            }
        }
        if (!cashFlowSlice->isZero()) 
            *mainSlice += *cashFlowSlice;
    }
    catch (exception& e) { throw ModelException(e, "KTurboTree::update"); }
}

/************************************ KTurbo *********************************/

FDProductSP KTurbo::createProduct(FDModel * model) const {
    if (!setupCalled) throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KTurboTree(KTurboConstSP(this), model));
}

static void checkBadSize(const string &name, int a, int b) {
    if (a==b) return;
    throw ModelException(name+" size ("+Format::toString(a)
        +") should equal number of coupons ("+Format::toString(b)+")");
}

void KTurbo::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        int nbCoupons = sched->nbCoupons;

        if (notional.get()) {
            checkBadSize("notional", notional->size(), nbCoupons);
            if (floor.get()) checkBadSize("floor", floor->size(), nbCoupons);
            if (cap.get()) checkBadSize("cap", cap->size(), nbCoupons);
            if (!dcc.get()) throw ModelException("dcc needed if notional provided");
        } else {
            if (floor.get() || cap.get()) throw ModelException(
                "Notional must be provided if cap/floor are provided");
        }
        if (legs.size()<1) throw ModelException("legs.size()=0 but at least one leg is needed.");
        int legI;
        for (legI=0; legI<legs.size(); ++legI) {
            if (!legs[legI].get() || legs[legI]->indexes.size()>1)
                throw ModelException("legs["+Format::toString(legI)
                +"].indexes.size() must be 1 (temporary limitation of this component)");
            legs[legI]->dcfs.resize(nbCoupons);
        }
        dcfs.resize(nbCoupons);

        // calc dcf variables
        sched->setup(model, market);

        for (legI=0; legI<legs.size(); ++legI) {
            sched->calcDcfs(*legs[legI]->dcc, legs[legI]->dcfs);
        }
        if (dcc.get())
            sched->calcDcfs(*dcc, dcfs);

        if (notional.get() && notional->size() == nbCoupons) notional->push_back(0);
        discount.getData(model, market);
        for (legI=0; legI<legs.size(); ++legI) {
            Leg &leg = *(legs[legI]);
            // add final notional 0
            if (leg.notionals->size() == nbCoupons) leg.notionals->push_back(0);

            // creates fx KIndexSpecs
            leg.fx.reset(new IndexSpecFX(leg.discount->getName(), discount->getName()));
            leg.fx->addResetDates(sched->resetEff);
            leg.fx->setup(model, market);
            // send reset dates to the underlying indexes
            if (leg.indexes.size()) {
                DateTimeArray fixingDates;
                FlexDates::filterOnWeight(sched->resetEff, *leg.weights, fixingDates);
                leg.indexes[0]->addResetDates(fixingDates);
                leg.indexes[0]->setup(model, market);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, "KTurbo::setup");
    }
}

DateTime KTurbo::getLastDate(void) const {
    return sched->pay.back();
}

void KTurbo::reportEvents(const KnownCashflows*, IModel* model,
    const DateTime& eDate, EventResults* events) const
{
    try {
        CashFlowArraySP cashFlows(new CashFlowArray(0));
        CashflowInfoArraySP info(new CashflowInfoArray(0));
        DateTime valueDate = model->getValueDate();
        ObservationExact     observationExact;
        UnadjustedConvention unadjustedConvention;

        // principal payments
        for (int legI=0; legI<legs.size(); ++legI) {
            Leg &leg = *legs[legI];
            // payInitial (principal)
            if (leg.payInitial) {
                CashFlow cf;
                CashflowInfoSP cfi(new CashflowInfo);
                cfi->componentName = outputName;
                cfi->cfType = CashflowInfo::CfType::PRINCIPAL;

                cf.date = sched->accStart[0];
                cf.amount = (*leg.notionals)[0];
                
                cashFlows->push_back(cf);
                info->push_back(cfi);
            }
            
            // payPrincipal
            if (leg.payPrincipal) {
                for (int cpnI = 0; cpnI < sched->pay.size(); cpnI++) {
                    CashFlow cf;
                    CashflowInfoSP cfi(new CashflowInfo);
                    cfi->componentName = outputName;
                    cfi->cfType = CashflowInfo::CfType::PRINCIPAL;

                    cf.date = sched->pay[cpnI];
                    cf.amount = (*leg.notionals)[cpnI] - (*leg.notionals)[cpnI+1];
                    
                    cashFlows->push_back(cf);
                    info->push_back(cfi);  
                }
            }
        }

        // coupon payments
        for (int cpnI = 0; cpnI < sched->pay.size(); cpnI++) 
        {
            CashflowInfoSP cfi(new CashflowInfo);
            cfi->componentName = outputName;
            cfi->cfType = CashflowInfo::CfType::COUPON;
            CashFlow cf;
            cf.date = sched->pay[cpnI];
            
            if (model->getValueDate() <= sched->pay[cpnI]) {
                info->push_back(cfi);
                cashFlows->push_back(cf);
                continue;
            }        

            string keyPrefix = outputName+".";

            double coupon=0;
            cfi->amountType = CashflowInfo::AmountType::KNOWN;

            for (int legI=0; legI<legs.size(); ++legI) {
                Leg &leg = *legs[legI];
                string legPrefix = keyPrefix+"leg["+Format::toString(legI)+"].";
                double pmt=0;

                // compute leg payment
                if (!Maths::isZero((*leg.weights)[cpnI])) {
                    pmt = leg.indexes[0]->getValue(sched->resetEff[cpnI], *cfi);
                    pmt *= (*leg.weights)[cpnI];
                    cfi->push(legPrefix+"weight", (*leg.weights)[cpnI]);
                }
                pmt += (*leg.spreads)[cpnI];
                cfi->push(legPrefix+"spread", (*leg.spreads)[cpnI]);

                pmt *= (*leg.notionals)[cpnI] * leg.dcfs[cpnI];
                cfi->push(legPrefix+"notional", (*leg.notionals)[cpnI]);
                cfi->push(keyPrefix+"accStart", sched->accStart[cpnI]);
                cfi->push(keyPrefix+"accEnd", sched->accEnd[cpnI]);
                cfi->push(legPrefix+"dcf", leg.dcfs[cpnI]);

                if (leg.discount->getCcy() != discount->getCcy()) {
                    pmt *= leg.fx->getValue(sched->resetEff[cpnI], *cfi);
                }
                cfi->push(legPrefix+"amount", pmt);
                coupon += pmt;
            }
            if (cfi->amountType == CashflowInfo::AmountType::UNKNOWN)
                continue;

            // compute turbo leg payment
            if (notional.get()) {
                double capFloorBase = fabs((*notional)[cpnI]) * dcfs[cpnI];

                coupon = MAX(MIN(coupon, capFloorBase * (*cap)[cpnI]),
                                        capFloorBase * (*floor)[cpnI]);

                // store the results into requested structures
                cfi->push("capFloorNotional",(*notional)[cpnI]);
                cfi->push("cap",capFloorBase * (*cap)[cpnI]);
                cfi->push("floor",capFloorBase * (*floor)[cpnI]);
            }
            cfi->push(keyPrefix+"ccy",discount->getCcy());
            info->push_back(cfi);

            cf.amount = coupon;
            cashFlows->push_back(cf);
        }
        ASSERT(info->size() == cashFlows->size());
        if (cashFlows->size()) events->addEvent(new KnownCashflows(eDate, cashFlows, discount->getCcy(), info));
    }
    catch (exception& e) {
        throw ModelException(e, "KTurbo::getEvents");
    }
}

void KTurbo::getAccEndDates(DateTimeArray &accEndDates) const {
    accEndDates.insert(accEndDates.end(), sched->accEnd.begin(), sched->accEnd.end());
}

void KTurbo::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(KTurbo, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(legs, ""); FIELD_MAKE_OPTIONAL(legs);

    FIELD(sched, "");

    FIELD(notional, ""); FIELD_MAKE_OPTIONAL(notional);
    FIELD(floor, ""); FIELD_MAKE_OPTIONAL(floor);
    FIELD(cap, ""); FIELD_MAKE_OPTIONAL(cap);
    FIELD(dcc, ""); FIELD_MAKE_OPTIONAL(dcc);

    FIELD(dcfs, ""); FIELD_MAKE_TRANSIENT(dcfs);
    Addin::registerConstructor(Addin::UTILITIES, KTurbo::TYPE);
 }

CClassConstSP const KTurbo::TYPE = CClass::registerClassLoadMethod(
    "KTurbo", typeid(KTurbo), KTurbo::load);

DEFINE_TEMPLATE_TYPE(KTurboArray);

/************************************ KTurbo::Leg *********************************/

void KTurbo::Leg::validatePop2Object(void) {
    try {
        int nbCoupons = notionals->size();

        if (nbCoupons<1) throw ModelException("notional: One row per coupon, at least one coupon expected");
        checkBadSize("weights", weights->size(), nbCoupons);
        if (!spreads.get()) spreads.reset(new DoubleArray(weights->size(), 0));
        checkBadSize("spreads", spreads->size(), nbCoupons);

        int nbIS = indexes.size();
        if (nbIS==0) {
            for (int i=0; i<nbCoupons; ++i) {
                if (!Maths::isZero((*weights)[i])) throw ModelException(
                    "indexes empty but weights are not all zeros");
            }
        }
        else if (nbIS != 1 && nbIS != nbCoupons) throw ModelException(
            "indexes: Just one row or one row per coupon expected");
    }
    catch (exception& e) {
        throw ModelException(e, "KTurbo::Leg::validatePop2Object");
    }
}

void KTurbo::Leg::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(KTurbo::Leg, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(notionals, "");
    FIELD(spreads, "Fixed spread/coupon"); FIELD_MAKE_OPTIONAL(spreads);
    FIELD(weights, "Weight/coupon for floating index");
    FIELD(discount, "");
    FIELD(payInitial, "");
    FIELD(payPrincipal, "");
    FIELD(payExercise, ""); FIELD_MAKE_OPTIONAL(payExercise);
    FIELD(rateType, "");
    FIELD(dcc, "");
    FIELD(indexes, ""); FIELD_MAKE_OPTIONAL(indexes);
    FIELD(dcfs, ""); FIELD_MAKE_TRANSIENT(dcfs);
    FIELD(fx, ""); FIELD_MAKE_TRANSIENT(fx);
    Addin::registerConstructor("KTurbo_Leg", Addin::UTILITIES, "Creates a handle to KTurbo_Leg", KTurbo::Leg::TYPE);
}
CClassConstSP const KTurbo::Leg::TYPE = CClass::registerClassLoadMethod(
    "KTurbo::Leg", typeid(KTurbo::Leg), KTurbo::Leg::load);

typedef KTurbo::LegArray KTurboLegArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("KTurbo::LegArray", KTurboLegArray);

bool KTurboLoad() { return KTurbo::TYPE != 0; }

DRLIB_END_NAMESPACE
