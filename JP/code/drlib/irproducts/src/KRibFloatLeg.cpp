//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KRibFloatLeg.cpp
//
//   Description : floater component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KRibFloatLeg.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/RateTree.hpp"

DRLIB_BEGIN_NAMESPACE

/////////// KRibFloatLegTree class /////////////
class KRibFloatLegTree : public FDProduct, public IGetCfDate {
public:

    /******************** variables ********************/
    KRibFloatLegConstSP inst;
    ZeroBondProdSP   zeroProd;
    FDProductSP      indexProd;
    int              couponPaid; // index (in inst->resetDates[]) of the last coupon paid

    // convenience copies from inst
    int rateType;
    CouponSchedDates &sched;

    /******************** methods ********************/
    // constructor
    KRibFloatLegTree(const KRibFloatLegConstSP &inst, FDModel* model);

    /** initialisation, called ONCE only after InitState() for each new model instance */
    virtual void init(Control*) const;

    /** initialising and setting product variables */
    // this is called per pricing call before each pricing
    virtual void initProd(void);
    virtual void update(int & step, FDProduct::UpdateType);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual const TreeSlice & getCashFlow(int step) const { return *cashFlowSlice; }

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
    TreeSliceSP mainSlice;
    TreeSliceSP cashFlowSlice;
    TreeSliceSP ribCpn;
    TreeSliceSP ribIndex;
    TreeSliceSP obsIndex;

    bool resetInArrears;  // ??? temporary field

    // observation underlying is a composite Knockout at this stage that
    // returns true/false weights for the event observations
    FDProductSP  ribObsInd;
    int          nbObs;

};

/******************************* KRibFloatLegTree ******************************/

static bool happensNow(const DateTimeArray& array, DateTime stepDate, int *pos=0) 
{
    for (int i=0; i < array.size(); ++i) {
        if (array[i] == stepDate) {
            if (pos) *pos = i;
            return true;
        }
    }
    return false;
}

KRibFloatLegTree::KRibFloatLegTree(const KRibFloatLegConstSP &inst, FDModel* model) :
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
                indexProd->addModelResetDates(sched.resetEff);
                break;
            }
        }
    }
    string discYCName = inst->discountYieldCurveName();
    if (discYCName.empty())
        throw ModelException("KRibFloatLegTree::KRibFloatLegTree",
                             "discount YieldCurve name of component must be set");
    
    // if in rib mode, no need for zero discounter as the discounting is done by DEV
    // in the leg itself
    if (!inst->rib.get()) {
        zeroProd = DYNAMIC_POINTER_CAST<ZeroBondProd>(model->createProduct(IProdCreatorSP(new ZeroBond(
                   sched.resetEff, sched.pay, discYCName))));
    }

    // ??? for now, check that the last rib observation date is < swap leg end date
    if (inst->rib.get()) {
        if (inst->rib->obsUnd.get()) {
            ribObsInd = model->createProduct(inst->rib->obsUnd);

            if (!ribObsInd.get()) 
                throw ModelException("Missing ribObsInd");

            ribObsInd->addModelResetDates(inst->rib->obsDates);
            nbObs = 0;
        }
        DateTime swapEndDate = inst->sched->accEnd.back();
        DateTime lastObsDate = inst->rib->obsEffDates.back();
        if (lastObsDate >= swapEndDate)
            throw ModelException("Final RIB Observation date " + 
                                lastObsDate.toString() +
                                " must be before final coupon accrual end date " +
                                swapEndDate.toString());
    }

}

void KRibFloatLegTree::init(Control*) const{
    model->addCritDates(sched.resetEff);
    model->addCritDates(sched.pay);
    model->addCritDates(inst->principalDates);

    // add RIB critical dates
    if (inst->rib.get()) {
        if (inst->rib->obsEffDates.size() != inst->rib->obsDates.size())
            throw ModelException("KRibFloatLegTree::init", 
                                 "Rates tree model requires RIB observation effective/spot dates as "
                                 "well as actual observation dates");

        model->addCritDates(inst->rib->obsEffDates);
    }
}

void KRibFloatLegTree::initProd(void){
    string discCurve = inst->discount->getName();
    mainSlice = model->createSlice(inst->discount->getName());
    mainSlice->name = getOutputName() + "_mainSlice";
    *mainSlice = 0.;
    startDEV(mainSlice);

    cashFlowSlice = model->createSlice(inst->discount->getName());
    cashFlowSlice->name = getOutputName() + "_cashFlowSlice";
    *cashFlowSlice = 0.;

    // ??? these are created in potentially wrong dimension for now
    // RIB slices
    ribCpn = model->createSlice(inst->discount->getName());
    ribCpn->name = getOutputName() + "_ribCpn";
    *ribCpn=0.;
    startDEV(ribCpn);

    ribIndex = model->createSlice(inst->discount->getName());
    ribIndex->name = getOutputName() + "_ribIndex";
    *ribIndex=1.;

    obsIndex = model->createSlice(inst->discount->getName());
    obsIndex->name = getOutputName() + "_obsIndex";

    nbObs = 0;

    // for now, just check if in arrears or advance from the first coupon
    if (sched.accStart[0] == sched.reset[0])
        resetInArrears = false;
    else if (sched.accEnd[0] == sched.reset[0])
        resetInArrears = true;
    else
        throw ModelException("void KRibFloatLegTree::initProd",
                             "RibFloatLeg must explicitly fix in either arrears or advance");

}

const TreeSlice & KRibFloatLegTree::getValue(int step, DateTime eventDate) const {
    DateTime stepDate = model->getDate(step);
    DateTime lastPayDate = inst->sched->pay.back();

    if (eventDate != stepDate) {
        throw ModelException(__FUNCTION__, "Cannot be valued at a date != currentDate");
    }

    if (stepDate > lastPayDate)
        throw ModelException("KRibFloatLegTree::getValue", 
                             "Unable to calculate getValue function as request at step " +
                             Format::toString(step) + " (date =  " +
                             stepDate.toString() +
                             ") is after the last payment date defined in the leg " +
                             lastPayDate.toString());

    // ??? if stepDate = 0, return the value of the leg - not sure if this is completely
    // safe to do, and will have to review logic
    if (step == 0)
        return *mainSlice;

    // ??? for now RIB does not support stubs, so this function can only return the
    // value of the FloatingLeg on a reset date

    if (!happensNow(inst->rib->accStartDates, stepDate)
        && stepDate < inst->sched->accEnd.back()
        && stepDate > inst->sched->accStart.front()) 
    {
        throw ModelException("KRibFloatLegTree::getValue", 
                             "KRibFloatLeg does not support stubs ("+stepDate.toString()+")");
    }
    return *mainSlice;
}


void KRibFloatLegTree::update(int & step, FDProduct::UpdateType update) {
    try {
        DateTime stepDate = model->getDate(step);
        DateTime today = model->getToday();
        int accrualPos;
        int i;
        // local slice references to make debugging easier
        TreeSlice& cashFlowSliceL = *cashFlowSlice;
        TreeSlice& ribCpnL = *ribCpn;
        TreeSlice& ribIndexL = *ribIndex;
        cashFlowSliceL = 0.;

        // principal repayment, pays initial, final coupon payments and any intermediate payments
        for (i = 0; i < inst->principalDates.size(); ++i) {
            DateTime princDate = inst->principalDates[i];
            if (princDate == stepDate) {
                double principalPmt = inst->principalPayments[i];

                if (today <= stepDate) {
                    // don't add payment if on the value date
                    cashFlowSliceL += principalPmt;
                }
                break;
            }
        }

        // if step == 0, assume observation events on value date are past fixings
        // and are handled below, otherwise
        // if rib observation date, add rib obs coupon slice value
        if (inst->rib.get() && happensNow(inst->rib->obsDates, stepDate) 
            && (today <= stepDate))
        {
            const TreeSlice &indSlice = ribObsInd->getValue(step, stepDate);
            ribCpnL += ribIndexL * indSlice;
            nbObs++;
        }

        if (happensNow(inst->sched->accStart, stepDate, &accrualPos)) {
            // on the accrual start date, calculate the floating reset + spread, or if there is
            // no floating reset reset = spread.  Multiply this reset slice by the unit "Ribbed"
            // slice to calculate the actual RIB reset slice for this current coupon
            double notional = (*inst->notionals)[accrualPos];
            double dcf = (*inst->dcfs)[accrualPos];

            if (resetInArrears == false) {
                double weight = (*inst->weights)[accrualPos];
                double spread = (*inst->spreads)[accrualPos];
                if (!Maths::isZero(weight)) {
                    const TreeSlice &indexSlice = indexProd->getValue(step, stepDate);
                    ribIndexL = (indexSlice * weight) + spread;
                }
                else
                    ribIndexL = spread;
                
                // apply cap and floor
                if (inst->capRates.get())
                    ribIndexL = smin(ribIndexL, (*inst->capRates)[accrualPos]); 
                if (inst->floorRates.get())
                    ribIndexL = smax(ribIndexL, (*inst->floorRates)[accrualPos]); 

                ribCpnL *= ribIndexL;
            }

            // Calculate the current RIB coupon and add it to the floating leg's mainSlice
            // on the accrual start date as at this point the absolute value of the 
            // ribIndex is known 

            if (inst->rateType!=RateType::SIMPLE)
                throw ModelException("Analytics only support simple rates- not compounding");

            // ??? later on relax this so that we can have 0 rib Obs dates for coupon and
            // thus more flexible RIB date schedules
            if (nbObs == 0) {
                // non rib leg
                cashFlowSliceL += (notional * dcf ) * ribIndexL;
            }
            else {
                cashFlowSliceL += (notional * dcf / (double)nbObs) * ribCpnL;
            }
            nbObs = 0;
            ribCpnL = 0.;
        }

        if (happensNow(inst->sched->accEnd, stepDate, &accrualPos)) {
            if (resetInArrears == true) {
                double weight = (*inst->weights)[accrualPos];
                double spread = (*inst->spreads)[accrualPos];
                if (indexProd.get()) {
                    const TreeSlice &indexSlice = indexProd->getValue(step, stepDate);
                    ribIndexL = (indexSlice * weight) + spread;
                }
                else
                    ribIndexL = spread;
                
                // apply cap and floor
                if (inst->capRates.get())
                    ribIndexL = smin(ribIndexL, (*inst->capRates)[accrualPos]); 
                if (inst->floorRates.get())
                    ribIndexL = smax(ribIndexL, (*inst->floorRates)[accrualPos]); 
            }
            else
                ribIndexL = 1.;

            startDEV(ribIndex);
        }

        if (step==0) {
            CashflowInfo cfi;
            cfi.amountType = CashflowInfo::AmountType::KNOWN;
            DateTime valueDate = model->getValueDate();
            cfi.cfType = CashflowInfo::CfType::COUPON;

            for (int cpnI=0; cpnI<sched.nbCoupons; ++cpnI) {
                DateTime payDate  = sched.pay[cpnI];
                DateTime accStart = sched.accStart[cpnI];
                DateTime accEnd   = sched.accEnd[cpnI];
                ASSERT(payDate==accEnd);

                if (!((accStart < valueDate) && (valueDate <= payDate)))
                    continue;
                if (inst->rib.get()) {
                    DateTimeArray const &obsDates = inst->rib->obsDates.getDates();

                    double pastObsProfile = 0;
                    for (int i=0; i<obsDates.size(); ++i) {
                        DateTime obsDate = obsDates[i];

                        // note that if obsDate == valueDate the value has already 
                        // been added above (as a slice)
                        if (accStart <= obsDate && obsDate < accEnd && obsDate < valueDate) {
                            pastObsProfile += inst->rib->obsUnd->getValue(obsDates[i], cfi);
                            nbObs++;
                        }
                    }
                    ribCpnL += ribIndexL * pastObsProfile;
                }
                

                if (!resetInArrears) {
                    double weight = (*inst->weights)[cpnI];
                    double spread = (*inst->spreads)[cpnI];

                    if (!Maths::isZero(weight)) {
                        double idxVal = inst->index->getValue(accStart, cfi);
                        ribIndexL = idxVal * weight + spread;
                    }
                    else {
                        ribIndexL = spread;
                    }
                    // apply cap and floor
                    if (inst->floorRates.get())
                        ribIndexL = smax(ribIndexL, (*inst->floorRates)[cpnI]); 
                    if (inst->capRates.get())
                        ribIndexL = smin(ribIndexL, (*inst->capRates)[cpnI]); 
                    
                    ribCpnL *= ribIndexL;
                }
                double notional = (*inst->notionals)[cpnI];
                double dcf = (*inst->dcfs)[cpnI];
                if (nbObs == 0) {
                    cashFlowSliceL += (notional * dcf ) * ribIndexL;
                }
                else {
                    cashFlowSliceL += (notional * dcf / (double)nbObs) * ribCpnL;
                }
            }
            if (cfi.amountType == CashflowInfo::AmountType::UNKNOWN)
                throw ModelException("Error calculating past fixing logic for RIB");
        }
        if (!cashFlowSlice->isZero()) 
            *mainSlice += *cashFlowSlice;
    }
    catch (exception& e) {
        throw ModelException(e, "KRibFloatLegTree::update");
    }
}

void KRibFloatLegTree::printInfo(ostream& outputStream) const {
    string prodType = "FRibLoateLegTree";
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

    // print out Rib Details
    outputStream << endl << "RIB observation dates" << endl;
    for (i = 0; i < inst->rib->obsDates.size(); i++)
        outputStream << inst->rib->obsDates[i].toString() << endl;


    // internal slices and slice dimensions info
    outputStream << endl;
    outputStream << "Internal product slice details:" << endl;
    outputStream << "instanceName   dimension   devCurveName   descriptiveName" << endl;
    outputStream << "mainSlice   " << DYNAMIC_POINTER_CAST<TreeSliceRates>(mainSlice)->getDim() << "   " <<
                     mainSlice->getCurveToDEV() << endl; // "   " << mainSlice->getSliceName() << endl;

}

/******************************* KRibFloatLeg ******************************/

static void checkBadSize(const string outputName, const string &name, int a, int b) {
    if (a==b) return;
    throw ModelException(name+" size ("+Format::toString(a)
        +") should equal payDates size (number of coupons = "+Format::toString(b)+")");
}

static void missing(const string &name) {
    throw ModelException("Missing "+name+" array");
}

void KRibFloatLeg::validatePop2Object(void) {
    try {
        if (!sched.get()) 
            missing("sched");
        if (!notionals.get()) 
            missing("notionals");
    }
    catch (exception& e) {
        throw ModelException(e, "KRibFloatLeg::validatePop2Object");
    }
}

void KRibFloatLeg::getAccEndDates(DateTimeArray &accEndDates) const {
    accEndDates.insert(accEndDates.end(), sched->accEnd.begin(), sched->accEnd.end());
}

// ??? must add cap/floor to known cashflows
void KRibFloatLeg::reportEvents(const KnownCashflows*, IModel* model,
    const DateTime& eDate, EventResults* events) const
{
    try {
        ObservationExact observationExact;
        UnadjustedConvention unadjustedConvention;
        CashFlowArraySP cashFlows(new CashFlowArray(0));
        CashflowInfoArraySP info(new CashflowInfoArray(0));

        // principal payments
        for (int i=0; i<principalPayments.size(); ++i) 
        {
            CashflowInfoSP cfi(new CashflowInfo);
            cfi->componentName = outputName;
            cfi->cfType = CashflowInfo::CfType::PRINCIPAL;
            info->push_back(cfi);

            CashFlow cf;
            cf.amount = principalPayments[i];
            cf.date = principalDates[i];
            cashFlows->push_back(cf);
        }     
           
        // coupon payments
        for (int cpnI = 0; (cpnI < sched->pay.size()) 
                           && (sched->pay[cpnI] <= model->getValueDate()); cpnI++) 
        {
            CashflowInfoSP cfi(new CashflowInfo);
            CashFlow cf;
            cfi->componentName = outputName;
            cfi->cfType = CashflowInfo::CfType::COUPON;
            cf.date = sched->pay[cpnI];

            if (model->getValueDate() <= sched->pay[cpnI]) {
                info->push_back(cfi);
                cashFlows->push_back(cf);
                continue;
            }            

            string keyPrefix = outputName+".";

            double pmt=0;
            cfi->amountType = CashflowInfo::AmountType::KNOWN;

            // compute leg payment
            if (!Maths::isZero((*weights)[cpnI]))
            {
                pmt = index->getValue(sched->resetEff[cpnI], *cfi);
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

            pmt *= rib->getValue(sched->accStart[cpnI], *cfi);
            ASSERT(cfi->amountType != CashflowInfo::AmountType::UNKNOWN);

            cfi->push(keyPrefix+"ccy",discount->getCcy());
            info->push_back(cfi);

            cf.amount = pmt;
            cashFlows->push_back(cf);
        }
        ASSERT(info->size() == cashFlows->size());
        if (cashFlows->size()) 
            events->addEvent(new KnownCashflows(eDate, cashFlows, discount->getCcy(), info));   
    }
    catch (exception& e) {
        throw ModelException(e, "KRibFloatLeg::getEvents");
    }
}

void KRibFloatLeg::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        if (rib.get()) {
            rib->setAcsDates(sched->accStart);
            rib->setAceDates(sched->accEnd);
            rib->setup(model, market); // must find a better way to do that
        }
        if (index.get()) {
            index->addResetDates(sched->resetEff);
            index->setup(model, market);
        }

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
        if (spreads->size()!=nbCoupons) 
            throw ModelException("spreads.size="+Format::toString(spreads->size())+
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

        // can supply either/or/both cap/floor rates.
        if (capRates.get())
            checkBadSize(outputName, "capRates", capRates->size(), nbCoupons);

        if (floorRates.get())
            checkBadSize(outputName, "floorRates", floorRates->size(), nbCoupons);

        if (capRates.get() && floorRates.get()) {
            for (i = 0; i < nbCoupons; i++) {
                if ((*capRates)[i] <= (*floorRates)[i])
                    throw ModelException("When both supplied, floorRate must be < capRate. "
                                         "For coupon index " + Format::toString(i) +
                                         ", capRate = " + Format::toString((*capRates)[i]) +
                                         " but floorRate = " + Format::toString((*floorRates)[i]));
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
    }
    catch (exception& e) {
        throw ModelException(e, "KRibFloatLeg::setup, " + outputName);
    }
}

DateTime KRibFloatLeg::getLastDate() const{
    return sched->pay.back();
}

FDProductSP KRibFloatLeg::createProduct( FDModel * model ) const {
    if (!setupCalled) throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KRibFloatLegTree(KRibFloatLegConstSP(this), model));
}

void KRibFloatLeg::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Wrapper style Simple IR swap interface component");
    REGISTER(KRibFloatLeg, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(sched,"");
    FIELD(notionals, "Notional value for each coupon");
    FIELD(spreads, "Spread as value (ie. not in basis points or percentage)");
    FIELD_MAKE_OPTIONAL(spreads);
    FIELD(weights, "Weight to apply to each floating index for each coupon");
    FIELD_MAKE_OPTIONAL(weights);
    FIELD(dcfs,"");
    FIELD_MAKE_OPTIONAL(dcfs);
    FIELD(capRates, "Rate as decimal to calculate coupon payment cap amount");
    FIELD_MAKE_OPTIONAL(capRates);
    FIELD(floorRates, "Rate as decimal to calculate coupon payment floor amount");
    FIELD_MAKE_OPTIONAL(floorRates);
    FIELD(rateType, "");
    FIELD_MAKE_OPTIONAL(rateType);
    FIELD(payInitialPrincipal,
        "Pay Initial notional value as principal payment on accrual start date - false by default");
    FIELD_MAKE_OPTIONAL(payInitialPrincipal);
    FIELD(payPrincipal,
        "Paid when notional amortizes and pays final coupon notional on final coupon payment date - false by default");
    FIELD_MAKE_OPTIONAL(payPrincipal);
    FIELD(dcc, "");
    FIELD_MAKE_OPTIONAL(dcc);
    FIELD(index, "floating/reset index");
    FIELD_MAKE_OPTIONAL(index);
    FIELD(rib, "RIB interface component containing rib schedules and observation configuration");
    FIELD_MAKE_OPTIONAL(rib);

    FIELD(principalDates,"");    FIELD_MAKE_TRANSIENT(principalDates);
    FIELD(principalPayments,""); FIELD_MAKE_TRANSIENT(principalPayments);
    Addin::registerConstructor(Addin::UTILITIES, KRibFloatLeg::TYPE);
}

CClassConstSP const KRibFloatLeg::TYPE = CClass::registerClassLoadMethod(
    "KRibFloatLeg", typeid(KRibFloatLeg), KRibFloatLeg::load);

bool KRibFloatLegLoad(){ return (KRibFloatLeg::TYPE != 0); }

DRLIB_END_NAMESPACE
