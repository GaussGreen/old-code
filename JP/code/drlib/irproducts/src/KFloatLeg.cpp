//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KFloatLeg.cpp
//
//   Description : floater component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KFloatLeg.hpp"
#include "edginc/RateTree.hpp"
#include "edginc/KFloatLegMC.hpp"

DRLIB_BEGIN_NAMESPACE

/////////// KFloatLegTree class /////////////
class KFloatLegTree : public FDProduct, public IGetCfDate {
public:

    /******************** variables ********************/
    KFloatLegConstSP inst;
    ZeroBondProdSP   zeroProd;
    FDProductSP      indexProd;
    int              couponPaid; // index (in inst->resetDates[]) of the last coupon paid

    // convenience copies from inst
    int rateType;
    CouponSchedDates &sched;

    /******************** methods ********************/
    // constructor
    KFloatLegTree(const KFloatLegConstSP &inst, FDModel* model);

    /** initialisation, called ONCE only after InitState() for each new model instance */
    virtual void init(Control*) const;

    /** initialising and setting product variables */
    // this is called per pricing call before each pricing
    virtual void initProd(void);

    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int & step, FDProduct::UpdateType);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual const TreeSlice & getCashFlow(int step) const { 
        //DateTime stepDate = model->getDate(step);
        //if (!RatesUtils::happensNow(inst->sched->resetEff, stepDate, 0)) {
            //*cashFlowSlice =0.0;
            //throw ModelException("KFloatLeg::getCashFlow", "Cannot get cash flow of "
            //+inst->outputName+" on "+stepDate.toString());
        //}
        return *cashFlowSlice; 
    }

    virtual DateTime getStartDate(void) const { return model->getDate(0);}
    virtual string getOutputName(void) const { return inst->outputName; }

    void printInfo(ostream& outputStream) const;

    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
        string ccy;
        if (inst->discount.get()) ccy = inst->discount->getCcy();
        recordSliceToOutputName(ctrl, results, model,
            inst->isTopComponent(), inst->outputName, ccy, getValue(0, model->getDate(0)));
        inst->recordExtraOutput(ctrl, results);
    }
    virtual void getCfDate(DateTimeArray &cfDates) { 
        cfDates.insert(cfDates.end(), inst->sched->resetEff.begin(), inst->sched->resetEff.end());
    }
private:
    TreeSliceSP mainSlice, getValueSlice, cashFlowSlice;
    // to help debugging and know which product instance this is, take copy of the instrument's
    // outputName field as the debuggers are unable to get from instrumentSP in the watch windows
    string instOutputName;
};

/******************************* KFloatLegTree ******************************/

KFloatLegTree::KFloatLegTree(const KFloatLegConstSP &inst, FDModel* model) :
    FDProduct(model), inst(inst), couponPaid(-1),
    // convenience aliases & copies from inst
    rateType(inst->rateType), sched(*inst->sched)
{
    if (inst->index.get()) {
        // if any weights are not zero, register index
        int i;
        for (i = 0; i < inst->weights->size(); i++) {
            if (!Maths::isZero((*inst->weights)[i])) {
                indexProd = model->createProduct(inst->index);
                indexProd->addModelResetDates(inst->sched->resetEff);
                break;
            }
        }
    }
    string discYCName = inst->discountYieldCurveName();
    if (discYCName.empty())
        throw ModelException("KFloatLegTree::KFloatLegTree",
                             "discount YieldCurve name of component must be set");
    zeroProd = DYNAMIC_POINTER_CAST<ZeroBondProd>(model->createProduct(IProdCreatorSP(new ZeroBond(
        sched.resetEff, sched.pay, discYCName))));
}

void KFloatLegTree::init(Control*) const{
    model->addCritDates(sched.resetEff);
    model->addCritDates(sched.pay);
    model->addCritDates(inst->principalDates);
}

void KFloatLegTree::initProd(void){
    string discCurve = inst->discount->getName();

    mainSlice = model->createSlice(discCurve);
    mainSlice->name = inst->outputName + "_mainSlice";
    *mainSlice = 0.;
    startDEV(mainSlice);

    getValueSlice = model->createSlice();
    getValueSlice->name = inst->outputName + "_getValueSlice";

    cashFlowSlice = model->createSlice();
    cashFlowSlice->name = inst->outputName + "_cashFlowSlice";
    *cashFlowSlice = 0.;

    instOutputName = inst->outputName;  // to help debugging, store local copy of instrument name
}


const TreeSlice & KFloatLegTree::getValue(int step, DateTime eventDate) const {
    
    DateTime stepDate = model->getDate(step);
    int cpnIdx; // coupon index, loop counter
    int nbCoupons = sched.pay.size();
        
    if (eventDate != stepDate) {
        throw ModelException(__FUNCTION__, "Cannot be valued at a date != currentDate");
    }

    if (inst->sched->pay.back() < stepDate)
        throw ModelException("KFloatLegTree::getValue", "Unable to getValue for leg on " +
                             stepDate.toString() + " as last leg payment date is " +
                             inst->sched->pay.back().toString());

    // if getValue is not a reset date or a <= accrual start date, stub is required and
    // we don't yet support stubs
    bool resetDateFound = false;
    for (cpnIdx=0; cpnIdx<nbCoupons; ++cpnIdx) {
        DateTime resetDate = inst->sched->resetEff[cpnIdx];
        if (stepDate == resetDate) {
            resetDateFound = true;
            break;
        }
    }
    if (step != 0 && (stepDate > inst->sched->accStart.front()) && (resetDateFound == false) &&
        (stepDate < inst->sched->pay.back()))
        throw ModelException("KFloatLegTree::getValue", "Leg does not support stubs, getValue "
                             "date = " + stepDate.toString() + " at time step = " +
                             Format::toString(step));
    
    *getValueSlice = *mainSlice;

    if (step != 0)
        return *getValueSlice; 

    // ??? keep logic separate for each type to simplify removal of 
    // includePricingValDateEvents later on, even though they're very similar
    
    if (model->getToday().getTime() == DateTime::START_OF_DAY_TIME)
    {
        for (cpnIdx=0; cpnIdx<nbCoupons; ++cpnIdx) {
            DateTime payDate = inst->sched->pay[cpnIdx];
            DateTime resetDate = inst->sched->resetEff[cpnIdx];

            if ((payDate >= stepDate) &&  (resetDate < model->getValueDate())) 
            {
                // past fixing required
                double fixing=0.;
                double weight = (*inst->weights)[cpnIdx];
                if (!Maths::isZero(weight)) 
                {
                    IMarketObservable *mobs = dynamic_cast<IMarketObservable*>(inst->index.get());
                    if (!mobs) throw ModelException("index does not implement IMarketObservable");

                    UnadjustedConvention unadjustedConvention;
                    ObservationExact observationExact;
                    fixing = weight * mobs->pastValue(resetDate, &observationExact, 
                                            IMarketObservable::getDefaultObsSource().get(), 
                                            0, 0, &unadjustedConvention);
                }
                double spread = (*inst->spreads)[cpnIdx];
                double dcf = (*inst->dcfs)[cpnIdx];
                double notional = (*inst->notionals)[cpnIdx];

                const TreeSlice &zeroSlice  = zeroProd->getValue(step, payDate);

                if (rateType==RateType::SIMPLE) {
                    *getValueSlice += notional * zeroSlice * dcf * (fixing + spread);
                } else {
                    *getValueSlice += notional * zeroSlice * (::pow (1. + (fixing +  spread), dcf) - 1.);
                }
            }
        }
    }
    else
    {
        // drop events on the valueDate - reset = pastReset and thus deterministic payment on this date
        // and if payment on the valuedate, drop it completely
        for (cpnIdx=0; cpnIdx<nbCoupons; ++cpnIdx) {
            DateTime payDate = inst->sched->pay[cpnIdx];
            DateTime resetDate = inst->sched->resetEff[cpnIdx];

            if ((payDate > stepDate) &&  (resetDate <= model->getValueDate())) 
            {
                double fixing=0.;
                double weight = (*inst->weights)[cpnIdx];
                if (!Maths::isZero(weight)) 
                {
                    IMarketObservable *mobs = dynamic_cast<IMarketObservable*>(inst->index.get());
                    if (!mobs) throw ModelException("index does not implement IMarketObservable");

                    UnadjustedConvention unadjustedConvention;
                    ObservationExact observationExact;
                    fixing = weight * mobs->pastValue(resetDate, &observationExact, 
                                            IMarketObservable::getDefaultObsSource().get(), 
                                            0, 0, &unadjustedConvention);
                }
                double spread = (*inst->spreads)[cpnIdx];
                double dcf = (*inst->dcfs)[cpnIdx];
                double notional = (*inst->notionals)[cpnIdx];

                const TreeSlice &zeroSlice  = zeroProd->getValue(step, payDate);

                if (rateType==RateType::SIMPLE) {
                    *getValueSlice += notional * zeroSlice * dcf * (fixing + spread);
                } else {
                    *getValueSlice += notional * zeroSlice * (::pow (1. + (fixing +  spread), dcf) - 1.);
                }
            }
        }
    }
    return *getValueSlice; 
}

void KFloatLegTree::update(int & step, FDProduct::UpdateType update) {
    try {
        DateTime stepDate = model->getDate(step);
        int cpnIdx; // coupon index, loop counter
        int i;
        int nbCoupons = inst->sched->pay.size();
        TreeSlice &cashFlowSliceL = *cashFlowSlice;
        cashFlowSliceL = 0.;

        // ??? if valueDate (for now step==0), if ignoring any events on this day
        // there is nothing to do and can simply return
        // ??? remove this option later on when we've checked enough wrappers vs qlib
        // and assume events are always dropped on value date
        if (stepDate < model->getToday())
            return;

        // principal repayment, pays initial, final coupon payments and any intermediate payments
        for (i=0; i<inst->principalDates.size(); ++i) {
            DateTime princDate = inst->principalDates[i];
            if (princDate == stepDate) {
                double principalPmt = inst->principalPayments[i];
                cashFlowSliceL += principalPmt;
                break;
            }
        }

        // coupon payment
        for (cpnIdx=0; cpnIdx<nbCoupons; ++cpnIdx) {
            if (sched.resetEff[cpnIdx] != stepDate) continue;


            couponPaid = cpnIdx;
            double spread   = (*inst->spreads)[cpnIdx];
            double dcf      = (*inst->dcfs)[cpnIdx];
            double notional = (*inst->notionals)[cpnIdx];
            double weight   = (*inst->weights)[cpnIdx];


            const TreeSlice & zeroSlice = zeroProd->getValue(step, sched.pay[cpnIdx]);

            if (!Maths::isZero(weight)) {
                const TreeSlice &indexSlice = indexProd->getValue(step, stepDate);
                if (rateType==RateType::SIMPLE) {
                    cashFlowSliceL += (notional * dcf * weight) * (indexSlice + spread/weight) * zeroSlice;
                } else {
                    cashFlowSliceL += notional * (pow (1. + (weight * indexSlice + spread), dcf) - 1.) * zeroSlice;
                }
            } else {
                if (rateType==RateType::SIMPLE) {
                    cashFlowSliceL += (notional * dcf * spread) * zeroSlice;
                } else {
                    cashFlowSliceL += notional * (::pow (1. + spread, dcf) - 1.) * zeroSlice;
                }
            }
        }
        if (!cashFlowSlice->isZero()) 
            *mainSlice += *cashFlowSlice;
    }
    catch (exception& e) {
        throw ModelException(e, "KFloatLegTree::update");
    }
}

void KFloatLegTree::printInfo(ostream& outputStream) const {
    string prodType = "FLoateLegTree";
    string name = getOutputName();

    outputStream << "Product type: " << prodType << " name: " << name << endl << endl;

    int i;
    // ??? tidy up spacings/formatting
    outputStream << "AccrualStart   AccrualEnd   Payment   Notional   IndexWeight   Spread   dcf" << endl;
    for (i = 0; i < (sched.accStart).size(); i++) {
        outputStream << sched.accStart[i].toString() << "   " <<
                        sched.accEnd[i].toString() << "   " <<
                        sched.pay[i].toString() << "   " <<
                        (*inst->notionals)[i] << "   " <<
                        (*inst->weights)[i] << "   " <<
                        (*inst->spreads)[i] << "   " <<
                        (*inst->dcfs)[i] << endl;
    }

    outputStream << endl;
    outputStream << "Pay Initial Principal: " << (inst->payInitialPrincipal ? "TRUE" : "FALSE") << endl;
    outputStream << "Pay Principal(s): " << (inst->payPrincipal ? "TRUE" : "FALSE") << endl;
    for (i = 0; i < inst->principalDates.size(); i++) {
        outputStream << (inst->principalDates)[i].toString() << "   " <<
                        (inst->principalPayments)[i] << endl;
    }

    outputStream << endl;
    if (!(!inst->index))
        outputStream << "Underlying index type : " << typeid(inst->index).name() << endl;
    else
        outputStream << "Reset/Floating index not supplied" << endl;

    // internal slices and slice dimensions info
    outputStream << endl;
    outputStream << "Internal product slice details:" << endl;
    outputStream << "instanceName   dimension   devCurveName   descriptiveName" << endl;
    outputStream << "mainSlice   " << DYNAMIC_POINTER_CAST<TreeSliceRates>(mainSlice)->getDim() << "   " <<
                     mainSlice->getCurveToDEV() << endl; // "   " << mainSlice->getSliceName() << endl;
    outputStream << "getValueSlice   " << DYNAMIC_POINTER_CAST<TreeSliceRates>(getValueSlice)->getDim() << "   " <<
                     getValueSlice->getCurveToDEV() << endl; // "   " << getValueSlice->getSliceName() << endl;

}

/******************************* KFloatLeg ******************************/

static void checkBadSize(const string outputName, const string &name, int a, int b) {
    if (a==b) return;
    throw ModelException(name+" size ("+Format::toString(a)
        +") should equal payDates size (number of coupons = "+Format::toString(b)+")");
}

static void missing(const string &name) {
    throw ModelException("Missing "+name+" array");
}

void KFloatLeg::validatePop2Object(void) {
    try {
        if (!sched.get()) missing("sched");
        if (!notionals.get()) missing("notionals");
        //if (!weights.get()) missing("weights");
        //if (!spreads.get()) missing("spreads");
    }
    catch (exception& e) {
        throw ModelException(e, "KFloatLeg::validatePop2Object");
    }
}

void KFloatLeg::reportCashFlows(CashflowInfoArray &cashflowInfos, bool amountsNeeded ) const
{
    try {
        DateTime today = getToday();
        ObservationExact     observationExact;
        UnadjustedConvention unadjustedConvention;

        // principal payments
        for (int i=0; i<principalPayments.size(); ++i) 
        {
            CashflowInfoSP cfi(new CashflowInfo);
            cfi->componentName = outputName;
            cfi->cfType = CashflowInfo::CfType::PRINCIPAL;
            cfi->amount = principalPayments[i];
            cfi->date = principalDates[i];
            cfi->updateAmountType(CashflowInfo::AmountType::KNOWN);
            cashflowInfos.push_back(cfi);
        }     

        // coupon payments
        for (int cpnI = 0; cpnI < sched->pay.size(); cpnI++) {
            CashflowInfoSP cfi(new CashflowInfo);
            cashflowInfos.push_back(cfi);
            cfi->componentName = outputName;
            cfi->cfType = CashflowInfo::CfType::COUPON;
            cfi->date = sched->pay[cpnI];

            if (!amountsNeeded || sched->pay[cpnI] > today)
                continue;

            string keyPrefix = outputName+".";

            double pmt=0;
            cfi->updateAmountType(CashflowInfo::AmountType::KNOWN);

            // compute leg payment
            if (!Maths::isZero((*weights)[cpnI])) {
                pmt = index->getValue(sched->resetEff[cpnI], *cfi);
                if (cfi->amountType == CashflowInfo::AmountType::UNKNOWN)
                    continue;

                pmt *= (*weights)[cpnI];
                cfi->push(keyPrefix+"weight", (*weights)[cpnI]);
            }

            pmt += (*spreads)[cpnI];
            cfi->push(keyPrefix+"spread", (*spreads)[cpnI]);

            if (rateType==RateType::SIMPLE)
                pmt = (*dcfs)[cpnI] * pmt;
            else
                pmt = (::pow(1+pmt, (*dcfs)[cpnI]) - 1);

            cfi->push(keyPrefix+"accStart", sched->accStart[cpnI]);
            cfi->push(keyPrefix+"accEnd", sched->accEnd[cpnI]);
            cfi->push(keyPrefix+"dcf", (*dcfs)[cpnI]);

            pmt *= (*notionals)[cpnI];
            cfi->push(keyPrefix+"notional", (*notionals)[cpnI]);

            cfi->push(keyPrefix+"estimated", sched->resetEff[cpnI] > today);
            cfi->push(keyPrefix+"amount", pmt);
            cfi->push(keyPrefix+"ccy",discount->getCcy());
            cfi->amount = pmt;
        }
    }
    catch (exception& e) {
        throw ModelException(e, "KFloatLeg::getEvents");
    }
}

void KFloatLeg::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        sched->setup(model, market);

        int i;
        int nbCoupons = sched->nbCoupons;

        // If the underlying index is a constant, we set index==NULL,
        // weights = 0, spreads = constant index value.
        bool noIndex = (index.get()==NULL);

        if (weights.get() == NULL ) {
            if (noIndex)
                weights = DoubleArraySP(new DoubleArray(nbCoupons,0.0));
            else // If the underlying index is a KRib, weights and spreads are processed by KRib.
                weights = DoubleArraySP(new DoubleArray(nbCoupons,1.0));;
        }
        if ( spreads.get() == NULL )
            spreads = DoubleArraySP(new DoubleArray(nbCoupons,0.0));
        if (weights->size()!=nbCoupons)
            throw ModelException("weights.size="+Format::toString(weights->size())+
                                 " != nbCoupons="+Format::toString(nbCoupons));

        // complement the notionals schedule
        checkBadSize(outputName, "notionals", notionals->size(), nbCoupons);

        if (dcfs.get()) {
            checkBadSize(outputName, "dcfs", dcfs->size(), nbCoupons);
        }
        else {
            if (!dcc.get())
                throw ModelException("Both dcc and dcfs are missing");

            dcfs = DoubleArraySP(new DoubleArray(nbCoupons));
            sched->calcDcfs(*dcc, *dcfs);

            for (i = 0; i < nbCoupons; i++) {
                if (noIndex && !Maths::isZero((*weights)[i])) {
                    throw ModelException("Non null weights["
                    +Format::toString(i)+"] and no index specified");
                }
            }
        }

        // precalculate principal payments to simplify the update
        if (payInitialPrincipal) {
            // initial payment
            principalDates.push_back(sched->accStart[0]);
            principalPayments.push_back(-(*notionals)[0]);

        }
        if (payPrincipal) {
            for (i = 0; i < nbCoupons-1; i++) {
                // for each coupon
                double principalPmt = (*notionals)[i] - (*notionals)[i+1];
                if (!Maths::isZero(principalPmt)) {
                    principalDates.push_back(sched->pay[i]);
                    principalPayments.push_back(principalPmt);
                }
            }
            // add final payment
            double lastPmt = (*notionals).back();
            principalDates.push_back(sched->pay.back());
            principalPayments.push_back(lastPmt);
       }
        if (index.get()) {
            index->addResetDates(sched->resetEff);
            index->setup(model, market);
        }
    }
    catch (exception& e) {
        throw ModelException(e, "KFloatLeg::setup, " + outputName);
    }
}

void KFloatLeg::getAccEndDates(DateTimeArray &accEndDates) const {
    accEndDates.insert(accEndDates.end(), sched->accEnd.begin(), sched->accEnd.end());
}

DateTime KFloatLeg::getLastDate() const{
    return sched->pay.back();
}

FDProductSP KFloatLeg::createProduct( FDModel * model ) const{
    if (!setupCalled) throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KFloatLegTree(KFloatLegConstSP(this), model));
}

IMCProduct* KFloatLeg::createProduct(const MonteCarlo* model) const {
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    DateTimeArray resetDates(sched->resetEff.getDates());
    simSeries->addDates(DateTimeArray(resetDates.begin(), resetDates.end())); //, resetDates));
    InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
    return new KFloatLegMC(simSeries, instSettle,
                           getValueDate(),
                           sched,
                           discount.getSP(),
                           index,
                           notionals,
                           weights,
                           spreads,
                           dcfs,
                           rateType,
                           principalDates,
                           principalPayments
                          );
}

void KFloatLeg::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Wrapper style Simple IR swap interface component");
    REGISTER(KFloatLeg, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    IMPLEMENTS(IMCIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(sched,"");
    FIELD(notionals, "Notional value for each coupon");
    FIELD(spreads, "Spread as value (ie. not in basis points or percentage)");
    FIELD_MAKE_OPTIONAL(spreads);
    FIELD(weights, "Weight to apply to each floating index for each coupon");
    FIELD_MAKE_OPTIONAL(weights);
    FIELD(dcfs,"");
    FIELD_MAKE_OPTIONAL(dcfs);
    FIELD(rateType, "");
    FIELD_MAKE_OPTIONAL(rateType);
    FIELD(payInitialPrincipal,
        "Pay Initial notional value as principal payment on accrual start date - true by default");
    FIELD_MAKE_OPTIONAL(payInitialPrincipal);
    FIELD(payPrincipal,
        "Paid when notional amortizes and pays final coupon notional on final coupon payment date - true by default");
    FIELD_MAKE_OPTIONAL(payPrincipal);
    FIELD(dcc, "");
    FIELD_MAKE_OPTIONAL(dcc);
    FIELD(index, "floating/reset index");
    FIELD_MAKE_OPTIONAL(index);
    FIELD(principalDates,"");    FIELD_MAKE_TRANSIENT(principalDates);
    FIELD(principalPayments,""); FIELD_MAKE_TRANSIENT(principalPayments);
    Addin::registerConstructor(Addin::UTILITIES, KFloatLeg::TYPE);
}

CClassConstSP const KFloatLeg::TYPE = CClass::registerClassLoadMethod(
    "KFloatLeg", typeid(KFloatLeg), KFloatLeg::load);

bool KFloatLegLoad(){ return (KFloatLeg::TYPE != 0); }

DRLIB_END_NAMESPACE
