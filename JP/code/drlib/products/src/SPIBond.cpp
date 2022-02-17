//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIBond.cpp
//
//   Description : Bond interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPIBond.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Model.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/DayCountConventionFactory.hpp"

DRLIB_BEGIN_NAMESPACE

// this is the external interface for abstraction so that users can bolt in 
// any Bond type they want - note this includes the SPIBondWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface IBondSPI as soon as possible
void IBondSPIInterface::load(CClassSP& clazz) {
    REGISTER_INTERFACE(IBondSPIInterface, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}

CClassConstSP const IBondSPIInterface::TYPE = CClass::registerInterfaceLoadMethod(
    "IBondSPIInterface", typeid(IBondSPIInterface), IBondSPIInterface::load);

/*****************************************************************************/

IBondSPI::~IBondSPI() {}

const int SPIBondStd::stubNone = 0;
const int SPIBondStd::stubSwap = 1;
const int SPIBondStd::stubBond = 2;

IBondSPISP SPIBondStd::getBondSPI() {
    static const string routine = "SPIBondStd::getBondSPI";
    // essentially just returns itself
    IBondSPISP theBond = IBondSPISP(this, NullDeleter()); // FIXME can be a problem
    return theBond;
}

void SPIBondStd::validatePop2Object(){
    static const string method = "SPIBondStd::validatePop2Object";

    // To allow IMS to operate we must allow all fields to be optional
    // and construct ok, and fail later. Existence checks are therefore
    // made here but reported later
    // Not sure how to handle yc, datedDate, redemptionAmt, fundingSpread.
    if (stubTypeString.empty()) {
        isOK = false;
        err = "No stubTypeString supplied.";
        return;
    }
    if (!cashFlows.get()) {
        isOK = false;
        err = "No cashFlows supplied.";
        return;
    }
    if (dayCountConvString.empty()) {
        isOK = false;
        err = "No dayCountConvString supplied.";
        return;
    }

    // haven't done this yet
    if (!Maths::isZero(fundingSpread)) {
        isOK = false;
        err = "Funding spread " + Format::toString(fundingSpread)
            + " must be 0.0 - not yet implemented";
        return;
    }

    // convert the day count convention string to an object
    accruedDCC = DayCountConventionSP(
        DayCountConventionFactory::make(dayCountConvString));

    if (cashFlows->size()<1) {
        //throw ModelException(method, "Require at least one cash flow!");
        isOK = false;
        err = "Require at least one cash flow!";
        return;
   }
    // validate that cashFlow dates are increasing
    int i;
    for (i = 1; i < cashFlows->size(); i++){
        if ((*cashFlows)[i-1].date.isGreater((*cashFlows)[i].date)){
/*                throw ModelException(method,
              "CashFlow dates are not increasing : [" +
              Format::toString(i-1) +
              "] " + (*cashFlows)[i-1].date.toString() +
              " > [" + Format::toString(i) + "] " +
              (*cashFlows)[i].date.toString()); */
            isOK = false;
            err = "CashFlow dates are not increasing : [" +
                Format::toString(i-1) +
                "] " + (*cashFlows)[i-1].date.toString() +
                " > [" + Format::toString(i) + "] " +
                (*cashFlows)[i].date.toString();
            return;
        }
    }
    // validate that datedDate is before first cashFlowDate
    if (datedDate.isGreater((*cashFlows)[0].date)) {
        /*throw ModelException(method, "datedDate " + datedDate.toString() +
                             " must not be after first cashFlowDate " +
                             (*cashFlows)[0].date.toString());*/
        isOK = false;
        err = "datedDate " + datedDate.toString() +
            " must not be after first cashFlowDate " +
            (*cashFlows)[0].date.toString();
        return;
    }

    switch (stubTypeString[0])
    {
    case 'n':
    case 'N':
        stubType = stubNone;
        break;
    case 's':
    case 'S':
        stubType = stubSwap;
        break;
    case 'b':
    case 'B':
        stubType = stubBond;
        break;
    default:
/*            throw ModelException(method, "Unrecognised stub type " + stubTypeString +
          "\nExpect one of None, Swap, Bond");*/
        isOK = false;
        err = "Unrecognised stub type " + stubTypeString +
            "\nExpect one of None, Swap, Bond";
        break;
    }
    if (hasLinearBondFloor) {
        if (!linearBondFloor.get()) {
            isOK = false;
            err = "Indicating a linear bond floor but none provided!";
        }
    }

}

DateTimeArray SPIBondStd::getEssentialDates() const {
    // don't bother ... bondDates.push_back(datedDate);
    return CashFlow::dates(*cashFlows.get());
}

// This may return 0 which means it is not (piecewise)linear
// Otherwise the Schedule returned IS the bond floor level
const Schedule* SPIBondStd::getLinearBondFloor() const {
    if (hasLinearBondFloor) {
        return linearBondFloor.get();
    }
    return 0;
}

void SPIBondStd::init(const DateTime&      today,
                      const DateTime&      maturity,
                      const DateTimeArray* simDates) {
    this->today = today;
    this->maturityDate = maturity;
    if (!simDates) {
        throw ModelException("SPIBondStd::init",
                             "No sim dates!");
    }
    this->simDates = simDates;

}

void SPIBondStd::getYieldCurveData(const IModel*          model,
                                   const CMarketDataSP    market) {
    yc.getData(model, market);
}

const YieldCurve* SPIBondStd::getYC() {
    return yc.get();
}

// used when rolling to turn off fixing curve
void SPIBondStd::setYC(YieldCurveWrapper& newYC) {
    YieldCurveSP clone = YieldCurveSP(dynamic_cast<YieldCurve*>(newYC.getMO()->clone()));
    setYC(clone);
}

// used when retrospectively building bond floor history
void SPIBondStd::setYC(YieldCurveSP newYC) {
    yc.setObject(newYC);
}

void SPIBondStd::getBondPrices(const YieldCurve*    disc,
                               DoubleArray&         bondPrices,
                               const DateTimeArray* forDates) {
    static const string routine = "SPIBondStd::getBondPrices";

    // allow local override
    const DateTimeArray* theDates = forDates?forDates:simDates;

    /* At each date will provide the PV of future complete coupons
       Plus possibly stubbed first coupon.
       By working back from maturity to today we can do this in a single pass.
       This should work for any sim date : before datedDate and after last cashFlow date */
    int        iCashFlow = cashFlows->size() - 1;
    double     price;
    double     dirtyPrice = 0;
    DateTime   dirtyDate = (*cashFlows)[iCashFlow].date; // doubles as accrueEndDate
    DateTime   bondMat = dirtyDate;
    bool  isRedemptionIncluded = false;
    int   iPastDate = pastBondValues.get()? pastBondValues->size()-1 : -1;

    if (dirtyDate.getDate() < maturityDate.getDate()) {
        // require final rebal date not after bond maturity
        throw ModelException(routine,
                             "Bond maturity (" +
                             dirtyDate.toString() +
                             ") must not be before final rebalance date ("
                             + maturityDate.toString() + ")");
    }

    int numDates = theDates->size();
    for (int iSimDate = numDates-1; iSimDate >=0 ; iSimDate--)
    {
        const DateTime& valueDate = (*theDates)[iSimDate];

        // Careful with past - and rolling both "to today" and for theta
        // To accomodate the "EURIBOR bond price today" we need to have strict "<" here
        // since else price today will be a past value and not calculated.
        if (valueDate.getDate() < today.getDate()) {
            // If past rely on user to supply the bond price - locate
            // based on date to be safe
            while (iPastDate>=0 && (*pastBondValues)[iPastDate].date.getDate() > valueDate.getDate()) {
                iPastDate--;
            }
            if (iPastDate>=0 &&
                (*pastBondValues)[iPastDate].date.getDate() == valueDate.getDate()) {
                if (Maths::isPositive((*pastBondValues)[iPastDate].amount)) {
                    price = (*pastBondValues)[iPastDate].amount;
                } else {
                    throw ModelException(routine,
                                         "Past Bond Price missing for " + valueDate.toString());
                }
            } else {
                throw ModelException(routine,
                                     "Past Bond Price missing for " + valueDate.toString());
            }
        } else if (valueDate.getDate() > bondMat.getDate()) {
            price = 0.0;
        } else {

            // Redemption amount included (once only) on maturity (final cash flow date)
            if (!isRedemptionIncluded)
            {
                dirtyPrice = redemptionAmt;
                isRedemptionIncluded = true;
            }

            // Capture cash flows future of this valueDate
            for(; iCashFlow>=0 &&
                    (*cashFlows)[iCashFlow].date.getDate() >= valueDate.getDate();
                iCashFlow--)
            {
                dirtyPrice *= disc->pv((*cashFlows)[iCashFlow].date, dirtyDate);
                dirtyDate = (*cashFlows)[iCashFlow].date;
                dirtyPrice += (*cashFlows)[iCashFlow].amount;
            }
            // Treat any stub - iCashFlow now references "next earlier" cash flow
            DateTime accrueStartDate = iCashFlow<0? datedDate : (*cashFlows)[iCashFlow].date;
            double   accrued = 0;
            if (stubType!=stubNone &&
                accrueStartDate < valueDate &&
                dirtyDate.getDate() >= valueDate.getDate()) // accrueEndDate >= valueDate means it's all past so ignore
            {
                // By validating at least one cash flow we are ok referencing [iCashFlow+1]
                accrued = (*cashFlows)[iCashFlow+1].amount *
                    accruedDCC->years(accrueStartDate, valueDate) /
                    accruedDCC->years(accrueStartDate, dirtyDate);

                if (stubType==stubSwap) {
                    price = dirtyPrice - accrued; // swap convention pays accrued at next coupon date
                } else {
                    price = dirtyPrice; // accrued handled later
                }
            }
            else
            {
                price = dirtyPrice;
            }

            // Price at valueDate ...
            price *= disc->pv(valueDate, dirtyDate);

            // Bond convention pays accrued at valueDate
            if (stubType==stubBond) {
                price -= accrued;
            }

            /* Before maturity a zero bond price is an error */
            if (Maths::isZero(price)) {
                throw ModelException(routine,
                                     "Zero bond price at " + (*theDates)[iSimDate].toString());
            }
        }
        bondPrices[iSimDate] = price;
    }
}

void SPIBondStd::getFutureBondPrices(vector<SVExpectedDiscFactorSP>& dfSVBond,
                                     int                             iFirstFutureStep,
                                     DoubleArray&                    bondPrices) {
    static const string routine = "SPIBondStd::getFutureBondPrices";

    // We support only zero coupon bonds 
    if (cashFlows->size()!=1) {
        throw ModelException(routine,
                             "We support a single cash flow only, but " +
                             Format::toString(cashFlows->size()) + 
                             " have been supplied");
    }
    // and stubs make little sense in this case
    if (stubType!=stubNone) {
        throw ModelException(routine,
                             "We support only stubType = stubNone");
    }
    DateTime  bondMat = (*cashFlows)[0].date; // XXX do this with a precomputed index
    double    notl = redemptionAmt + (*cashFlows)[0].amount;
    // only overwrite any that will change - does not include the 1.0 at mat
    // controlled by size of dfSVBond
    for (unsigned int i = iFirstFutureStep; i<dfSVBond.size(); i++) {
        bondPrices[i] = notl * dfSVBond[i]->firstDF(); // expected df at iSimDate from bondMat to iSimDate  
    }
}

void SPIBondStd::roll(const Theta::Util&   thetaUtil,
                      const YieldCurve*    disc) {
    // if no pastBondValues then nothing to do - note this means instruments will need to
    // have past value dates supplied if about to roll over a forward start date.
    // Same behaviour as asset PastValues
    if (!pastBondValues.get()) {
        return;
    }
    // operate on date only (not time)
    const DateTime& origDate = thetaUtil.getOriginalValueDate();
    const DateTime& newDate = thetaUtil.getNewValueDate();
    bool pastRollDate = false;
    DateTimeArray passedDates;
    IntArray passedIdx;
    for (int i = 0; !pastRollDate && i < pastBondValues->size(); i++){
        const DateTime& aDateTime = (*pastBondValues)[i].date;
        if ((aDateTime.getDate()>origDate.getDate() && aDateTime.getDate()<=newDate.getDate()) ||
            (aDateTime.getDate()==newDate.getDate() && Maths::isZero((*pastBondValues)[i].amount))) {
            passedDates.push_back(aDateTime);
            passedIdx.push_back(i);
        }
        if (aDateTime>newDate) {
            pastRollDate = true;
        }
    }
    // For the roll we use the discount curve and the curve in the bond for today
    // This is correct since the idea is for only prices "today" to use the curve in the bond
    // and all others the discount curve.
    if (passedDates.size()>0) {
        // This may well be called before any pricing, so before init() so we
        // patch for this case
        if (today.empty() || maturityDate.empty()) {
            init(origDate,  cashFlows->back().date, &passedDates);
        }

        
        DoubleArray passedPrices(passedDates.size());
        getBondPrices(disc, passedPrices, &passedDates);
        for(int j=0; j<passedIdx.size(); j++) {
            if (origDate.equals(passedDates[j], false)) {
                // note we use curve in bond for fixing values today
                // i.e. if we have a rebalance EOD and we're now SOD
                (*pastBondValues)[passedIdx[j]].amount = getBondPriceToday(origDate);
            } else {
                (*pastBondValues)[passedIdx[j]].amount = passedPrices[j];
            }
        }
    }
}

// This uses the special yield curve specified in the SPIBond
double SPIBondStd::getBondPriceToday(const DateTime& today) {
    try{
        DateTimeArray todayArray(1);
        todayArray[0] = today;
        DoubleArray Z0(1);
        getBondPrices(yc.get(), Z0, &todayArray);
        return Z0[0];
    } catch (exception& e){
        throw ModelException(e, "SPIBondStd::getBondPriceToday");
    }
}

SPIBondStd::SPIBondStd(): SPIInterfaceIMS(TYPE),
    redemptionAmt(0.0), fundingSpread(0.0),
    hasLinearBondFloor(false), linearBondFloor(0),
    today(), maturityDate(), simDates(0),
    stubType(0) {} // for reflection

IObject* SPIBondStd::defaultSPIBondStd(){
    return new SPIBondStd();
}

void SPIBondStd::load(CClassSP& clazz){
    REGISTER(SPIBondStd, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(IBondSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPIBondStd);
    FIELD(yc,                 "yc");
    FIELD(stubTypeString,     "N S or B");
    FIELD(datedDate,          "accrue start date for first cash flow");
    FIELD(cashFlows,                 "contiguous accrual periods");
    FIELD(redemptionAmt,      "included at mat if no accruing");
    FIELD(dayCountConvString, "day count conv for accrued calc");
    FIELD(pastBondValues,            "pastBondValues");
    FIELD_MAKE_OPTIONAL(pastBondValues);
    FIELD(fundingSpread,      "fundingSpread");
    FIELD(hasLinearBondFloor, "true=>override bond floor with supplied schedule");
    FIELD_MAKE_OPTIONAL(hasLinearBondFloor);
    FIELD(linearBondFloor,    "Override bond floor - only if hasLinearBondFloor");
    FIELD_MAKE_OPTIONAL(linearBondFloor);
    FIELD(bondPrices,         "bondPrices");
    FIELD_MAKE_TRANSIENT(bondPrices);
    FIELD(today,              "today");
    FIELD_MAKE_TRANSIENT(today);
    FIELD(maturityDate,       "maturityDate");
    FIELD_MAKE_TRANSIENT(maturityDate);
    FIELD(accruedDCC,                "accruedDCC");
    FIELD_MAKE_TRANSIENT(accruedDCC);
    FIELD(stubType,           "stubType");
    FIELD_MAKE_TRANSIENT(stubType);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(yc);
    FIELD_MAKE_OPTIONAL(stubTypeString);
    FIELD_MAKE_OPTIONAL(datedDate);
    FIELD_MAKE_OPTIONAL(cashFlows);
    FIELD_MAKE_OPTIONAL(redemptionAmt);
    FIELD_MAKE_OPTIONAL(dayCountConvString);
    FIELD_MAKE_OPTIONAL(fundingSpread);
    clazz->setPublic(); // make visible to EAS/spreadsheet
/*        Addin::registerConstructor("SPI_BOND_STD",
                               Addin::RISK,
                               "Creates an SPIBondStd instance",
                               TYPE);*/
}

CClassConstSP const SPIBondStd::TYPE = CClass::registerClassLoadMethod(
    "SPIBondStd", typeid(SPIBondStd), SPIBondStd::load);

/*****************************************************************************/

#define SPI_BOND_TYPE_STD   "Standard"

IBondSPISP SPIBondWrapper::getBondSPI() {
    static const string routine = "SPIBondWrapper::getBondSPI";
    // possibly some choices later
    if (!bondStd->isValid()) {
        throw ModelException(routine, "Invalid SPIBond (type " + SPIBondType +
                             ") : " + bondStd->errString());
    }
    IBondSPISP theBond = IBondSPISP(bondStd.get(), NullDeleter()); // FIXME can be a problem
    return theBond;
}

// validation
void SPIBondWrapper::validatePop2Object(){
    static const string routine = "SPIBondWrapper::validatePop2Object";

    if (SPIBondType.empty()){
        throw ModelException(routine, "Blank SPIBondType specified!");
    }
    if (SPIBondType==SPI_BOND_TYPE_STD) {
        if (!bondStd.get()) {
            throw ModelException(routine, "Expected bondStd but none supplied!");
        }
    } else {
        throw ModelException(routine, "Unrecognised SPIBondType " + SPIBondType +
                             ". Expected " + SPI_BOND_TYPE_STD);
    }
}

/** Invoked when Class is 'loaded' */
void SPIBondWrapper::load(CClassSP& clazz){
    REGISTER(SPIBondWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IBondSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPIBondWrapper);
    FIELD(SPIBondType, "SPIBondType");
    FIELD(bondStd,  "bondStd");
    FIELD_MAKE_OPTIONAL(bondStd);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

// for reflection
SPIBondWrapper::SPIBondWrapper(): CObject(TYPE){}

IObject* SPIBondWrapper::defaultSPIBondWrapper(){
    return new SPIBondWrapper();
}

CClassConstSP const SPIBondWrapper::TYPE = CClass::registerClassLoadMethod(
    "SPIBondWrapper", typeid(SPIBondWrapper), SPIBondWrapper::load);

DRLIB_END_NAMESPACE
