//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPILoanCost.cpp
//
//   Description : Loan Cost interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPILoanCost.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/
// this is the external interface for LoanCost so that users can bolt in 
// any LoanCost type they want - note this includes the SPILoanCostWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface ILoanCostSPI as soon as possible
void ILoanCostSPIInterface::load(CClassSP& clazz) {
    REGISTER_INTERFACE(ILoanCostSPIInterface, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}
CClassConstSP const ILoanCostSPIInterface::TYPE = CClass::registerInterfaceLoadMethod(
    "ILoanCostSPIInterface", typeid(ILoanCostSPIInterface), ILoanCostSPIInterface::load);

/*****************************************************************************/

ILoanCostSPI::~ILoanCostSPI() {};

void SPILoanCostStd::init(const DateTime&      today,
                          const YieldCurve*    ycForFixings,
                          double               Basis,
                          const DateTimeArray& rebalDates) {
    this->today = today;
    this->ycForFixings = ycForFixings;
    this->Basis = Basis;
    this->numDays = vector<int>(rebalDates.size(), 0);
    for(int i=1; i<rebalDates.size(); i++) {
        this->numDays[i] = rebalDates[i].daysDiff(rebalDates[i-1]);
    }
}

ILoanCostSPISP SPILoanCostStd::getLoanCostSPI() {
    static const string routine = "SPILoanCostStd::getLoanCostSPI";
    // essentially just returns itself
    ILoanCostSPISP theLoanCost = ILoanCostSPISP(this, NullDeleter()); // FIXME can be a problem
    return theLoanCost;
}


void SPILoanCostStd::refresh(const vector<SVExpectedDiscFactorSP>* dfSV) {
    this->dfSV = dfSV;
}

// This interface says it all - this is a very poorly defined class. 
// It should contain at least Basis and probably YieldCurve.
void SPILoanCostStd::roll(const Theta::Util&   thetaUtil,
                          const YieldCurve*    yc,
                          double               Basis) {
    // if no pastMMRates then nothing to do - note this means instruments will need to 
    // have past value dates supplied if about to roll over a forward start date.
    // Same behaviour as asset PastValues
    // note also nothing to do if we're not using Libor (i.e. spreadOnly)
    if (!pastMMRates.get() || spreadOnly) {
        return;
    }
    const DateTime& origDate = thetaUtil.getOriginalValueDate();
    const DateTime& newDate = thetaUtil.getNewValueDate();
    bool            pastRollDate = false;
    
    for (int i = 0; !pastRollDate && i < pastMMRates->size() - 1; i++){
        const DateTime& aDateTime = (*pastMMRates)[i].date;
        if ((aDateTime.getDate()>origDate.getDate() && aDateTime.getDate()<=newDate.getDate()) ||
            (aDateTime.getDate()==newDate.getDate() && Maths::isZero((*pastMMRates)[i].amount))) {
            // well, ...! XXX
            // Need to get something from market data, and an interval rate is reasonable.
            // A shame the endDate is not well defined, but theta then is arguably of little value.
            // With introduction of ycForFixings should probably use that - but of little value and more hassle.
            DateTime endDate = (i<pastMMRates->size()) ? (*pastMMRates)[i+1].date:
                aDateTime.rollDate(1);
            double ndays = endDate.daysDiff(aDateTime);
            double pv = yc->pv(aDateTime, endDate);
            (*pastMMRates)[i].amount = Basis / ndays * (1.0 / pv - 1.0);
        }
        if (aDateTime>newDate) {
            pastRollDate = true;
        }
    }
}

//helper for libor fees object
double SPILoanCostStd::getLiborRate(const YieldCurve* yc,
                                    const DateTime&   startDate, 
                                    const DateTime&   endDate,
                                    double            Basis) {
    static const string method = "SPILoanCostStd::getLiborRate";
    int                 numPast = pastMMRates.get()? pastMMRates->size() : 0;
    double              ndays = endDate.daysDiff(startDate);
    int                 iPast;
    double              yield;

    // Careful with past - and rolling both "to today" and for theta 
    for(iPast=0 ; iPast < numPast && (*pastMMRates)[iPast].date.getDate() < startDate.getDate(); iPast++) {
        // empty - just finding index
    }
    // now iPast refers to date on or after startDate

    if (startDate < today) {
        // If past rely on user to supply the yield - locate
        // based on date to be safe. 
        if (iPast < numPast &&
            (*pastMMRates)[iPast].date.getDate() == startDate.getDate()) {
            if (!Maths::isZero((*pastMMRates)[iPast].amount)) {
                if (iPast+1<numPast && (*pastMMRates)[iPast+1].date.getDate()==endDate.getDate()) {
                    // given past value to use
                    yield = (*pastMMRates)[iPast].amount;
                } else if (iPast+1<numPast && (*pastMMRates)[iPast+1].date.getDate()<endDate.getDate()) {
                    // when skipping dates may have gaps where we have part of a past rate and need
                    // to supplement with an interpolation
                    double pv1 = 1.0 / ( 1.0 +(*pastMMRates)[iPast].amount * ndays / Basis);
                    double pv2 = yc->pv((*pastMMRates)[iPast+1].date, endDate);
                    double pv12 = pv1 * pv2;
                    yield = Basis / ndays * (1.0 / pv12 - 1.0);
                } else {
                    throw ModelException(method, 
                                         "Past value missing for rate from " + startDate.toString() + 
                                         " to " + endDate.toString());
                }
            } else {
                throw ModelException(method, 
                                     "Past value missing for " + startDate.toString());
            }
        } else {
            throw ModelException(method, 
                                 "Past value missing for " + startDate.toString());
        }
    } else {
        if (startDate.getDate() == today.getDate()) {
            // We're using the dates in pastMMRates to define the rate
            // iPast refers to date on or after today
            if (iPast < numPast-1 &&
                (*pastMMRates)[iPast].date.getDate() == today.getDate()) {
                // record the rate which can be used as the next fixing...
                // This should be the same as endDate else we'll get an error
                // when we next price, but doing this is at least explicitly correct
                // (so long as noone changes the pastMMRates dates!)
                DateTime myEndDate = (*pastMMRates)[iPast+1].date;
                double   myPV = ycForFixings->pv(startDate, myEndDate);
                double   myNdays = myEndDate.daysDiff(startDate);
                double   fixingRate = Basis / myNdays * (1.0 / myPV - 1.0);
                newFixingMMRate = CashFlowSP(new CashFlow(startDate, fixingRate));
            }
        }
        double   pv = yc->pv(startDate, endDate);
    
        yield = Basis / ndays * (1.0 / pv - 1.0);
    }
    
    return yield;
}

double SPILoanCostStd::getAccrueFactor(const YieldCurve* yc,
                                       const DateTime&   startDate, 
                                       const DateTime&   endDate, 
                                       double            Basis) {
    
    static const string method = "SPILoanCostStd::getAccrueFactor";
    double              ndays = endDate.daysDiff(startDate);
    double              yield;
    double              spread;

    if (spreadOnly) {
        yield = 0.0;
    } else {
        yield = getLiborRate(yc, startDate, endDate, Basis);
    }
    if (startDate < today) {
        spread = pastSpread;
    } else {
        spread = futureSpread;
    }
    
    double factor = (yield + spread) * ndays / Basis;
    return factor;
}

double SPILoanCostStd::getFutureLiborRate(int               iStep) {
    double pv = (*dfSV)[iStep-1]->path()[0]; // expected DF from iStep to iStep-1.
    double yield = Basis / numDays[iStep] * (1.0 / pv - 1.0);
    return yield;
}

double SPILoanCostStd::getFutureAccrueFactor(int               iStep) {
    double              yield;

    if (spreadOnly) {
        yield = 0.0;
    } else {
        yield = getFutureLiborRate(iStep);
    }
    double factor = (yield + futureSpread) * numDays[iStep] / Basis;
    return factor;
}


CashFlowSP SPILoanCostStd::getLoanCostRateForPyramid() {
    // Only report this on dates in the pastMMRates list
    // And done for each "start date"
    return newFixingMMRate;
}

SPILoanCostStd::SPILoanCostStd(): SPIInterfaceIMS(TYPE), 
    pastSpread(0.0), futureSpread(0.0), spreadOnly(false), 
    ycForFixings(0), newFixingMMRate(0), Basis(0.) {} // for reflection

IObject* SPILoanCostStd::defaultSPILoanCostStd(){
    return new SPILoanCostStd();
}

/** Invoked when Class is 'loaded' */
void SPILoanCostStd::load(CClassSP& clazz){
    REGISTER(SPILoanCostStd, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(ILoanCostSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPILoanCostStd);
    FIELD(pastSpread,           "pastSpread");
    FIELD(futureSpread,         "futureSpread");
    FIELD(spreadOnly,        "spreadOnly");
    FIELD(pastMMRates,             "pastMMRates");
    FIELD(today,            "today");
    FIELD_MAKE_TRANSIENT(today);
    FIELD(Basis,            "Basis");
    FIELD_MAKE_TRANSIENT(Basis);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(pastSpread);
    FIELD_MAKE_OPTIONAL(futureSpread);
    FIELD_MAKE_OPTIONAL(pastMMRates);
    FIELD_MAKE_OPTIONAL(spreadOnly);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPILoanCostStd::TYPE = CClass::registerClassLoadMethod(
    "SPILoanCostStd", typeid(SPILoanCostStd), SPILoanCostStd::load);

/*****************************************************************************/

ILoanCostSPISP SPILoanCostWrapper::getLoanCostSPI() {
    static const string routine = "SPILoanCostWrapper::getLoanCostSPI";
    // possibly some choices later
    if (!loanCostStd->isValid()) {
        throw ModelException(routine, "Invalid SPILoanCost (type " + SPILoanCostType + 
                             ") : " + loanCostStd->errString());
    }
//    ILoanCostSPISP theLoanCost = ILoanCostSPISP::attachToRef(loanCostStd.get());
    ILoanCostSPISP theLoanCost(loanCostStd.get(), NullDeleter()); // FIXME FIXME: we take smartPtr and wrap into refCount!
    return theLoanCost;
}

// validation
void SPILoanCostWrapper::validatePop2Object(){
    static const string routine = "SPILoanCostWrapper::validatePop2Object";

    if (SPILoanCostType.empty()){
        throw ModelException(routine, "Blank SPILoanCostType specified!");
    }
    if (SPILoanCostType==SPI_LOAN_COST_TYPE_STD) {
        if (!loanCostStd.get()) {
            throw ModelException(routine, "Expected loanCostStd but none supplied!");
        }
    } else {
        throw ModelException(routine, "Unrecognised SPILoanCostType " + SPILoanCostType + 
                             ". Expected " + SPI_LOAN_COST_TYPE_STD);
    }
}

    /** Invoked when Class is 'loaded' */
void SPILoanCostWrapper::load(CClassSP& clazz){
    REGISTER(SPILoanCostWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ILoanCostSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPILoanCostWrapper);
    FIELD(SPILoanCostType, "SPILoanCostType");
    FIELD(loanCostStd,  "loanCostStd");
    FIELD_MAKE_OPTIONAL(loanCostStd);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
 
// for reflection
SPILoanCostWrapper::SPILoanCostWrapper(): CObject(TYPE){}

IObject* SPILoanCostWrapper::defaultSPILoanCostWrapper(){
    return new SPILoanCostWrapper();
}
CClassConstSP const SPILoanCostWrapper::TYPE = CClass::registerClassLoadMethod(
    "SPILoanCostWrapper", typeid(SPILoanCostWrapper), SPILoanCostWrapper::load);

DRLIB_END_NAMESPACE
