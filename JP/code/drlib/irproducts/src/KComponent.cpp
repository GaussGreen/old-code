//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KComponent.cpp
//
//   Description : generic instrument for components
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/RateTree.hpp"

DRLIB_BEGIN_NAMESPACE
/************************** KComponent::CashflowNode ****************************/

void KComponent::CashflowNode::setMultiplier(double multiplier) {
    for (int i=0; i<cashflows->size(); ++i) {
        (*cashflows)[i]->setMultiplier(multiplier);
    }
    for (int i=0; i<underlyings.size(); ++i) {
        underlyings[i]->setMultiplier(multiplier);
    }
}

void KComponent::CashflowNode::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(KComponent::CashflowNode, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(comp,"");
}

CClassConstSP const KComponent::CashflowNode::TYPE = CClass::registerClassLoadMethod(
    "KComponent::CashflowNode", typeid(KComponent::CashflowNode), KComponent::CashflowNode::load);

typedef KComponent::CashflowNodeArray KComponentCashflowNodeArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("KComponent::CashflowNodeArray", KComponentCashflowNodeArray);

/*********************************** KComponent *****************************/

string KComponent::discountYieldCurveName() const { 
    return discountYieldCurveWrapper().getName();
}

DateTime KComponent::getValueDate() const {
    /* This function should be called getToday(), not getValueDate()
     * To avoid confusion, always call:
     * model->getValueDate() for the value date
     * model->getToday() for today
     * nevertheless we keep returning today for compatibility with the tweaking system
     */
    return today; 
    /*throw ModelException(__FUNCTION__,"Not available");*/
}


bool KComponent::priceDeadInstrument(CControl* control, CResults* results) const {
    try {
        const DeadPricer* dp = dynamic_cast<const DeadPricer*>(this);
        if (!dp) 
            return false;

        double price;
        if (!dp->isDead(DateTime(), &price)) 
            return false;
        results->storePrice(price, discountYieldCurveWrapper()->getCcy());
        return true;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KComponent::setup(const IModel* model, const MarketData* market) {
    try {
        setupCalled=true;
        market->GetReferenceDate(today);

        modelDomesticYC.setName(model->getDomesticYCName());
        modelDomesticYC.getData(model, market);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

class MarketWrappersAction: public ObjectIteration::IAction {
public:
    const IModel *model;
    MarketData* market;

    MarketWrappersAction(const IModel *model, MarketData* market)
        : model(model), market(market) {}

    virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
        MarketObjectWrapper* mow = dynamic_cast<MarketObjectWrapper*>(state.getObject().get());
        if (!mow->getName().empty()) {
            mow->getData(model, market, mow->getMOType());
        }
        return false; // do not go into the market object we just fetched
    }
};

class CalcIsPhysical: public ObjectIteration::IAction {
public:
    vector<bool> stack;

    CalcIsPhysical() {
        stack.reserve(20);
        stack.push_back(true);
    }

    virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
        KComponent* comp = dynamic_cast<KComponent*>(state.getObject().get());
        size_t i = state.getDepth();
        if (i>=stack.size()) // some depths migh be skipped
            stack.resize(i+1, stack.back());
        comp->isPhysical = comp->isPhysical && stack[i];
        if (stack.size() < i+2)
            stack.resize(i+2);
        stack[i+1] = comp->isPhysical && comp->undIsPhysical;
        return true; // continue recursion
    }
};

void KComponent::GetMarket(const IModel *model, const CMarketDataSP market) {
    try {
        // setup top compoment specific things
        topComponent = true;
        ExpandableIProdCreator::expandAll(IProdCreatorSP(this));
        if (discount.getName().empty()) 
            throw ModelException("Top component (outputName = " + outputName + 
                                ", type = " + TYPE->getName() + 
                                ") must have a discount curve");

        // populate market wrappers to keep simple cases simple
        MarketWrappersAction action(model, market.get());
        ObjectIteration iter(MarketObjectWrapper::TYPE);
        iter.recurse(action, IObjectSP(this));

        // start the recursive setup() pyramid call
        setup(model, market.get());

        // spread includeNonPricingEvents to children
        CalcIsPhysical action2;
        ObjectIteration iter2(KComponent::TYPE);
        iter2.recurse(action2, IObjectSP(this));
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

YieldCurveConstSP KComponent::getActualDiscount() const {
    string modelCcy = getModelDomesticYC()->getCcy();
    string discCcy = discount->getCcy();
    if (modelCcy == discCcy)
        return discount.getSP();
    return getModelDomesticYC();
}

DateTime KComponent::getLastInstrDate(void) const
{
    struct Action: public ObjectIteration::IActionConst{
        DateTime lastDate;

        virtual bool invoke(const ObjectIteration::State& state, IObjectConstSP obj) {
            DateTime d = dynamic_cast<const KComponent*>(obj.get())->getLastDate();

            if (!d.empty() && (lastDate.empty() || lastDate.isLess(d))) 
                lastDate = d;

            return true;
        }
    };
    try {
        Action action;
        ObjectIteration iter(KComponent::TYPE);
        iter.recurse(action, IObjectConstSP(this));
        return action.lastDate;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


void KComponent::recordExtraOutput(Control* ctrl, Results* results) const {
    try {
        if (!isTopComponent()) 
            return;

        OutputRequest* request;

        request = ctrl->requestsOutput("LAST_INSTR_DATE");
        if (request) results->storeRequestResult(request, IObjectSP(new DateTime(getLastInstrDate())));

        request = ctrl->requestsOutput("PAYMENT_DATES");
        if (request) {
            DateTimeArray cfDate;
            getCfDates(cfDate);
            OutputRequestUtil::recordPaymentDates(ctrl, results, &cfDate); 
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

static void packCfDates(KComponent::CashflowNodeSP cashflowRoot, DateTimeArray &cfDate) {
    if (cashflowRoot->cashflows->size()) {
        CashflowInfoArray &cfi = *cashflowRoot->cashflows;
        for (int i=0; i<cfi.size(); ++i) {
            cfDate.push_back(cfi[i]->date);
        }
    }
    for (int i=0; i<cashflowRoot->underlyings.size(); ++i) {
        packCfDates(cashflowRoot->underlyings[i], cfDate);
    }
}

static void packKnownCashflows(KComponent::CashflowNodeSP cashflowRoot, EventResults* eventResults, DateTime eDate) {
    if (cashflowRoot->cashflows->size()) {
        CashflowInfoArray &cfi = *cashflowRoot->cashflows;
        CashFlowArraySP cashflows(new CashFlowArray);
        CashflowInfoArraySP cashflowInfos(new CashflowInfoArray);
        for (int i=0; i<cfi.size(); ++i) {
            if (cfi[i]->amountType == CashflowInfo::AmountType::KNOWN) {
                cashflows->push_back(cfi[i]->getCashFlow());
                cashflowInfos->push_back(cfi[i]);
            }
        }
        eventResults->addEvent(new KnownCashflows(eDate, cashflows, cashflowRoot->comp->discount->getCcy()));
    }
    for (int i=0; i<cashflowRoot->underlyings.size(); ++i) {
        packKnownCashflows(cashflowRoot->underlyings[i], eventResults, eDate);
    }
}

static void packAllCashflows(KComponent::CashflowNodeSP cashflowRoot, EventResults* eventResults, DateTime eDate) {
    if (cashflowRoot->cashflows->size()) {
        eventResults->addEvent(new AllCashflows(eDate, cashflowRoot->cashflows, cashflowRoot->comp->discount->getCcy()));
    }
    for (int i=0; i<cashflowRoot->underlyings.size(); ++i) {
        packAllCashflows(cashflowRoot->underlyings[i], eventResults, eDate);
    }
}

void KComponent::getEvents(const KnownCashflows*, IModel* model,const DateTime& eDate, EventResults* eventResults ) const
{
    try {
        if (!isTopComponent()) return;
        CashflowNodeSP cashflowRoot = reportCashflowsTree(true);
        packKnownCashflows(cashflowRoot, eventResults, eDate);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}    

void KComponent::getEvents(const AllCashflows*, IModel* model,const DateTime& eDate, EventResults* eventResults ) const
{
    try {
        if (!isTopComponent()) return;
        CashflowNodeSP cashflowRoot = reportCashflowsTree(true);
        packAllCashflows(cashflowRoot, eventResults, eDate);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}    

void KComponent::getCfDates(DateTimeArray &cfDate) const {
    try {
        CashflowNodeSP cashflowRoot = reportCashflowsTree(true);
        packCfDates(cashflowRoot, cfDate);
        DateTime::doSortUniq(cfDate);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

ModelException KComponent::makeException(exception &e, const string &functionName) const {
    return ModelException(e, functionName+"(), "+outputName);
}

KComponent::CashflowNodeSP KComponent::reportCashflowsTree(bool amountsNeeded) const {
    CashflowNodeSP cfn(new CashflowNode);
    cfn->comp = KComponentConstSP(this);
    reportCashFlows(*cfn->cashflows, amountsNeeded);
    if (cfn->cashflows->size()==0) 
        return CashflowNodeSP();
    return cfn;
}

void KComponent::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("QLIB product representation of existing swaption_t wrapper product");
    REGISTER(KComponent, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(IProdCreator);
    IMPLEMENTS(KnownCashflows::IEventHandler);
    IMPLEMENTS(AllCashflows::IEventHandler);
    IMPLEMENTS(LastSensDate);

    FIELD(outputName,"Name, to display outputRequests and debug messages");
    FIELD_MAKE_OPTIONAL(outputName);
    // comment too long for Excel to display 
    /*FIELD(recordOutputName, "flag governing whether the outputRequest for the component is called. "
                 "default Value = True.  This flag is added to essentially overcome an issue with the tree "
                 "where we don't want to call ParYield if there is no zeroDate, but we want to supply "
                 "outputNames to components in order to assist debugging in general.  This field will be "
                 "reviewed in due course, particularly when we formalize today/valueDate behaviour");*/
    FIELD(recordOutputName, "flag governing whether the outputRequest for the component is called. "
                 "default Value = True.  The use of this flag to be reviewed.");
    FIELD_MAKE_OPTIONAL(recordOutputName);
    FIELD(discount,"Discount curve");
    FIELD_MAKE_OPTIONAL(discount);

    FIELD(undIsPhysical,""); FIELD_MAKE_TRANSIENT(undIsPhysical);
    FIELD(isPhysical,""); FIELD_MAKE_TRANSIENT(isPhysical);
    FIELD(setupCalled,""); FIELD_MAKE_TRANSIENT(setupCalled);
    FIELD(marketRetrieved,""); FIELD_MAKE_TRANSIENT(marketRetrieved);
    FIELD(topComponent,""); FIELD_MAKE_TRANSIENT(topComponent);
    FIELD(modelDomesticYC,""); FIELD_MAKE_TRANSIENT(modelDomesticYC);
}

CClassConstSP const KComponent::TYPE = CClass::registerClassLoadMethod(
    "KComponent", typeid(KComponent), KComponent::load);

/************************ ExpandableIProdCreator ***************************/

// expands IProdCreators
class ExpandIProdCreatorAction: public ObjectIteration::IAction {
public:
    virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
        // needs to be changed so that it browses only ExpandableIProdCreator
        // it does not work because ExpandableIProdCreator is an interface
        if (!dynamic_cast<const ExpandableIProdCreator*>(obj.get())) return true;
        state.setObject(dynamic_cast<ExpandableIProdCreator*>(state.getObject().get())->expandIProdCreator());
        return true; // continue iteration
    }
};

IProdCreatorSP ExpandableIProdCreator::expandAll(const IProdCreatorSP &root) {
    try {
        // expands IProdCreator
        ExpandIProdCreatorAction action;
        ObjectIteration iter(CObject::TYPE);
        return IProdCreatorSP::dynamicCast(iter.recurse(action, root));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);     
    }
}

void ExpandableIProdCreator::load(CClassSP& clazz) {
    REGISTER_INTERFACE(ExpandableIProdCreator, clazz);
    EXTENDS(IProdCreator);
}

CClassConstSP const ExpandableIProdCreator::TYPE = CClass::registerClassLoadMethod(
    "ExpandableIProdCreator", typeid(ExpandableIProdCreator), ExpandableIProdCreator::load);


/******************************************************************************/

/** When to stop tweaking - delegate to Model (within Sensitivity) */
DateTime KComponent::endDate(const Sensitivity* sensControl) const
{
    return productEndDate(this, sensControl);
}

bool KComponentLoad(){
    return KComponent::TYPE != 0;
}

DRLIB_END_NAMESPACE

