//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIUtil.hpp
//
//   Description : Synthetic Portfolio Insurance utilities
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPIUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// a base class to handle the delayed validation in Pyramid
// so we trap after construction
SPIInterfaceIMS::SPIInterfaceIMS(): CObject(TYPE), 
        isOK(true), err("") {} // for reflection

SPIInterfaceIMS::SPIInterfaceIMS(CClassConstSP clazz): CObject(clazz), 
        isOK(true), err("") {} // for reflection

bool SPIInterfaceIMS::isValid() {
    return isOK;
}

string SPIInterfaceIMS::errString() {
    if (!isOK) {
        return err;
    }
    return "";
}

IObject* SPIInterfaceIMS::defaultSPIInterfaceIMS(){
    return new SPIInterfaceIMS();
}

/** Invoked when Class is 'loaded' */
void SPIInterfaceIMS::load(CClassSP& clazz){
    REGISTER(SPIInterfaceIMS, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultSPIInterfaceIMS);
    FIELD(isOK,            "isOK");
    FIELD_MAKE_TRANSIENT(isOK);
    FIELD(err,                "error string");
    FIELD_MAKE_TRANSIENT(err);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPIInterfaceIMS::TYPE = CClass::registerClassLoadMethod(
    "SPIInterfaceIMS", typeid(SPIInterfaceIMS), SPIInterfaceIMS::load);

/*****************************************************************************/

// for the output requests PAYMENT_DATES & KNOWN_CASHFLOWS

// knowing 'today' allows any known but not yet past (i.e. unpaid) 
// flows to be calculated and provided as additional value
SPIKnownFlows::SPIKnownFlows(const DateTime&   today,
                             const YieldCurve* disc,
                             const DateTime*   whenValued) :  // optional - default today
    today(today), knownFlows(0), disc(disc), whenValued(whenValued) {}
    
// MC handles its own PV so will be most convenient to get value
// at maturity date
// Actually we don't really want to do this inside the simulation
// - better to add as a constant afterwards. Though SV requires
// that it IS done within the sim.
double SPIKnownFlows::getUnpaidValue() {
    double unpaid = 0.0;
    double pvWhenValued = whenValued?disc->pv(*whenValued):1.0;
    // performance best if read backwards ...
    for(int i= (knownFlows.get()?knownFlows->size() : 0) -1;
        i>=0; i--) {
        if ((*knownFlows)[i].date<=today) {
            break;
        }
        double pvFlowDate = disc->pv((*knownFlows)[i].date);
        unpaid += (*knownFlows)[i].amount * pvFlowDate / pvWhenValued;
    }
    return unpaid;
}

// XXX Might need to have a pair of dates - the payment("notification") date and
// XXX the corresponding settlement date. For the moment (I think)
// XXX this implementation assumes a flow is added only when it
// XXX has been notified.
void SPIKnownFlows::addFlow(const DateTime& when,
                            double          howMuch) {
    // take care to maintain date order
    CashFlow cf(when, howMuch);
    CashFlowArraySP cfa(new CashFlowArray(0));
    cfa->push_back(cf);
    // merge needs 2 existing CFAs - so start carefully
    if (knownFlows.get()) {
        knownFlows = CashFlow::merge(cfa, knownFlows);
    } else {
        knownFlows = cfa;
    }
}

const CashFlowArray* SPIKnownFlows::getFlows() {
    return knownFlows.get();
}
DRLIB_END_NAMESPACE
