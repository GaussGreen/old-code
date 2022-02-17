//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIFees.cpp
//
//   Description : Fees interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPIFees.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/HolidayCollector.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/
// this is the external interface for Fees so that users can bolt in 
// any Fees type they want - note this includes the SPIFeesWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface IFeesSPI as soon as possible
void IFeesSPIInterface::load(CClassSP& clazz) {
    REGISTER_INTERFACE(IFeesSPIInterface, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}
CClassConstSP const IFeesSPIInterface::TYPE = CClass::registerInterfaceLoadMethod(
    "IFeesSPIInterface", typeid(IFeesSPIInterface), IFeesSPIInterface::load);

/*****************************************************************************/

IFeesSPI::~IFeesSPI() {}

// helper for getting settlement holidays
HolidayConstSP IFeesSPI::getSettleHolidays(const InstrumentSettlement* instSettle) {
    HolidayCollectorSP holidayVisitor = HolidayCollectorSP(
            new HolidayCollector());
    instSettle->accept(holidayVisitor.get());

    if (!holidayVisitor->getHoliday()) {
        return HolidayConstSP(Holiday::weekendsOnly());
    } else {
        return holidayVisitor->getHoliday();
    }
}

/*****************************************************************************/

// Fee based on basket value
SPIFeesStd::SPIFeesStd(): SPIInterfaceIMS(TYPE), feeTimeFactor(0),
        contingentBondFee(0.0), contingentThreshold(0.0), hasContingentFee(false),
        isPayDaily(true), paymentDates(0), realPaymentDates(0){} // for reflection

IFeesSPISP SPIFeesStd::getFeesSPI() {
    static const string routine = "SPIFeesStd::getFeesSPI";
    // essentially just returns itself
    IFeesSPISP theFees = IFeesSPISP(this, NullDeleter()); // FIXME can be a problem
    return theFees;
}

const string SPIFeesStd::feeType() {
    return SPI_FEES_TYPE_STD;
}

void SPIFeesStd::init(int numAlgAssets,
          const DateTimeArray& rebalDates,
          int                  iStepFirstRebal,
          double               Basis,
          ILoanCostSPI*  loanCost,
          const YieldCurve*    yc) {
    // check numAlgAssets agrees with number of equity fees
    if (equityFee.size() != numAlgAssets) {
        throw ModelException("SPIFeesStd::init",
                             "Require same number of equity fees (" +
                             Format::toString(equityFee.size()) +
                             ") as assets (" +
                             Format::toString(numAlgAssets) + ")");
    }
    if (hasContingentFee &&
            contingentEquityFee.size() != numAlgAssets) {
        throw ModelException("SPIFeesStd::init",
                             "Require same number of contingent equity fees (" +
                             Format::toString(contingentEquityFee.size()) +
                             ") as assets (" +
                             Format::toString(numAlgAssets) + ")");
    }
    // build feeTimeFactor
    feeTimeFactor = DoubleArray(rebalDates.size(), 0.0);
    for(int i=iStepFirstRebal+1; i<rebalDates.size(); i++) {
        feeTimeFactor[i] = (rebalDates[i].getDate()-rebalDates[i-1].getDate())/Basis;
    }
}

void SPIFeesStd::refresh(int iFirstFutureStep) {}

double SPIFeesStd::getFeeAmount(int                iStep,
                    double             nZ,
                    double             Z,
                    const DoubleArray& nE,
                    const DoubleArray& E) const {

    double F = nZ * Z * bondFee * feeTimeFactor[iStep];
    for(int iAsset=0; iAsset<equityFee.size(); iAsset++) {
        F += nE[iAsset] * E[iAsset] * equityFee[iAsset] * feeTimeFactor[iStep];
    }
    return F;
}

double SPIFeesStd::getContingentFeeAmount(int                iStep,
                              double             nZ,
                              double             Z,
                              const DoubleArray& nE,
                              const DoubleArray& E) const {

    double F = nZ * Z * contingentBondFee * feeTimeFactor[iStep];
    for(int iAsset=0; iAsset<equityFee.size(); iAsset++) {
        F += nE[iAsset] * E[iAsset] * contingentEquityFee[iAsset] * feeTimeFactor[iStep];
    }
    return F;
}


double SPIFeesStd::getPaidFeeAmount(int                iStep,
                                double             nZ,
                                double             Z,
                                const DoubleArray& nE,
                                const DoubleArray& E) const {
    return 0.0;
}

double SPIFeesStd::getContingentPaidFeeAmount(int                iStep,
                                              double             nZ,
                                              double             Z,
                                              const DoubleArray& nE,
                                              const DoubleArray& E) const {
    return 0.0;
}

bool SPIFeesStd::hasContingentFees() const {
    return hasContingentFee;
}

bool SPIFeesStd::isThresholdBreached(double level) const {
    return Maths::isPositive(level-contingentThreshold);
}

double SPIFeesStd::getFeeAtMin(const DoubleArray& exposureUpperBound,
                           double             equityExposureMin) const {
    double FEmin;
    if (equityFee.size()==1) {
        FEmin = equityFee[0];
    } else { // must be 2
        FEmin = (exposureUpperBound[1] * equityFee[1] +
                 (1. - exposureUpperBound[1]) * equityFee[0]);
    }
    return equityExposureMin*FEmin + (1.-equityExposureMin)*bondFee;
}

const DateTimeArray& SPIFeesStd::getNotificationDates(const DateTimeArray& rebalDates) const {
    return realPaymentDates;
}

DateTimeArray SPIFeesStd::getPaymentDates() const {
    return realPaymentDates;
}

DateTimeArray SPIFeesStd::getEssentialDates() const {
    return realPaymentDates;
}

// yeah this is a bit lame but the DFs for feeDates are basically built
// at every single rebalance date so we just return the date itself
// cince there are no paid fees this won't ever make a difference
// but we've got to return something
const DateTime SPIFeesStd::getFeePayDate(const DateTime& rebalDate) const {
    return rebalDate;
}

// post getMarket validation (handles settlement of fees)
void SPIFeesStd::Validate(const InstrumentSettlement* instSettle) {
    // do nothing
    // no paid fees so no settlement
}

IObject* SPIFeesStd::defaultSPIFeesStd(){
    return new SPIFeesStd();
}

/** Invoked when Class is 'loaded' */
void SPIFeesStd::load(CClassSP& clazz){
    REGISTER(SPIFeesStd, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(IFeesSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPIFeesStd);
    FIELD(bondFee,             "bondFee");
    FIELD(equityFee,           "equityFee");
    FIELD(contingentBondFee,             "Contingent Bond Fee");
    FIELD_MAKE_OPTIONAL(contingentBondFee);
    FIELD(contingentEquityFee,           "Contingent Equity Fee");
    FIELD_MAKE_OPTIONAL(contingentEquityFee);
    FIELD(contingentThreshold,           "Contingent Fee Threshold");
    FIELD_MAKE_OPTIONAL(contingentThreshold);
    FIELD(hasContingentFee,           "Has Contingent Fee?");
    FIELD_MAKE_OPTIONAL(hasContingentFee);
    FIELD(isPayDaily,          "isPayDaily");
    FIELD_MAKE_OPTIONAL(isPayDaily);
    FIELD(paymentDates,        "paymentDates");
    FIELD_MAKE_OPTIONAL(paymentDates);
    FIELD(realPaymentDates,    "realPaymentDates");
    FIELD_MAKE_OPTIONAL(realPaymentDates);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(bondFee);
    FIELD_MAKE_OPTIONAL(equityFee);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

CClassConstSP const SPIFeesStd::TYPE = CClass::registerClassLoadMethod(
    "SPIFeesStd", typeid(SPIFeesStd), SPIFeesStd::load);
/*****************************************************************************/
// Fee based on basket value wit some kept and some paid

SPIFeesKeptAndPaid::SPIFeesKeptAndPaid(): SPIInterfaceIMS(TYPE),
        keptBondFee(0.0), paidBondFee(0.0),
        contingentKeptBondFee(0.0), contingentPaidBondFee(0.0),
        contingentThreshold(0.0), hasContingentFee(false),
        hasSettlePeriod(false), settlePeriod(0),
        hasPaidNotDeductedFee(false) {} // for reflection

IFeesSPISP SPIFeesKeptAndPaid::getFeesSPI() {
    static const string routine = "SPIFeesKeptAndPaid::getFeesSPI";
    // essentially just returns itself
    IFeesSPISP theFees = IFeesSPISP(this, NullDeleter()); // FIXME can be a problem
    return theFees;
}

const string SPIFeesKeptAndPaid::feeType() {
    return SPI_FEES_TYPE_KNP;
}

// validation
void SPIFeesKeptAndPaid::validatePop2Object(){
    static const string routine = "SPIFeesKeptAndPaid::validatePop2Object";
    if (keptEquityFee.size() != paidEquityFee.size()) {
        isOK = false;
        err = "Kept and Paid equity fee arrays should be equal length!";
        return;
    }
    if (hasPaidNotDeductedFee &&
        paidNotDeductedEquityFee.size() != paidEquityFee.size()) {
        isOK = false;
        err = "Paid Not Deducted and Paid equity fee arrays should be equal length!";
        return;
    }
    totalEquityFee = DoubleArray(keptEquityFee.size());
    bool havePaidFee = !Maths::isZero(paidBondFee);
    for(int i=0; i<totalEquityFee.size(); i++) {
        totalEquityFee[i] = keptEquityFee[i] + paidEquityFee[i];
        if (!Maths::isZero(paidEquityFee[i])) {
            havePaidFee = true;
        }
        if (hasPaidNotDeductedFee && !Maths::isZero(paidNotDeductedEquityFee[i])) {
            havePaidFee = true;
        }
    }
    if (hasContingentFee) {
        if (contingentKeptEquityFee.size() != contingentPaidEquityFee.size()) {
            isOK = false;
            err = "Contingent Kept and Paid equity fee arrays should be equal length!";
            return;
        }
        totalContingentEquityFee = DoubleArray(contingentKeptEquityFee.size());
        if(!Maths::isZero(contingentPaidBondFee)) {
            havePaidFee = true;
        }
        for(int i=0; i<totalContingentEquityFee.size(); i++) {
            totalContingentEquityFee[i] = contingentKeptEquityFee[i]
                                            + contingentPaidEquityFee[i];
            if (!Maths::isZero(contingentPaidEquityFee[i])) {
                havePaidFee = true;
            }
        }
    }

    // If any non-zero paid fee require at least one fee payment date
    if (havePaidFee && paidFeePaymentDates.size()<1) {
        isOK = false;
        err = "Non-zero paid fees so require at least one paid fee payment date!";
        return;
    }
}

void SPIFeesKeptAndPaid::init(int                  numAlgAssets,
          const DateTimeArray& rebalDates,
          int                  iStepFirstRebal,
          double               Basis,
          ILoanCostSPI*  loanCost,
          const YieldCurve*    yc) {
    // check numAlgAssets agrees with number of equity fees
    if (keptEquityFee.size() != numAlgAssets) {
        throw ModelException("SPIFeesKeptAndPaid::init",
                             "Require same number of equity fees (" +
                             Format::toString(keptEquityFee.size()) +
                             ") as assets (" +
                             Format::toString(numAlgAssets) + ")");
    }
    if (hasContingentFee &&
            contingentKeptEquityFee.size() != numAlgAssets) {
        throw ModelException("SPIFeesKeptAndPaid::init",
                             "Require same number of contingent equity fees (" +
                             Format::toString(contingentKeptEquityFee.size()) +
                             ") as assets (" +
                             Format::toString(numAlgAssets) + ")");
    }
    // build feeTimeFactor
    feeTimeFactor = DoubleArray(rebalDates.size(), 0.0);
    for(int i=iStepFirstRebal+1; i<rebalDates.size(); i++) {
        feeTimeFactor[i] = (rebalDates[i].getDate()-rebalDates[i-1].getDate())/Basis;
    }
}

void SPIFeesKeptAndPaid::refresh(int iFirstFutureStep) {}

double SPIFeesKeptAndPaid::getFeeAmount(int                iStep,
                    double             nZ,
                    double             Z,
                    const DoubleArray& nE,
                    const DoubleArray& E) const {

    double F = nZ * Z * (keptBondFee+paidBondFee) * feeTimeFactor[iStep];
    for(int iAsset=0; iAsset<totalEquityFee.size(); iAsset++) {
        F += nE[iAsset] * E[iAsset] * totalEquityFee[iAsset] * feeTimeFactor[iStep];
    }
    return F;
}

double SPIFeesKeptAndPaid::getContingentFeeAmount(int                iStep,
                    double             nZ,
                    double             Z,
                    const DoubleArray& nE,
                    const DoubleArray& E) const {

    double F = nZ * Z * (contingentKeptBondFee+contingentPaidBondFee) * feeTimeFactor[iStep];
    for(int iAsset=0; iAsset<totalContingentEquityFee.size(); iAsset++) {
        F += nE[iAsset] * E[iAsset] * totalContingentEquityFee[iAsset] * feeTimeFactor[iStep];
    }
    return F;
}

double SPIFeesKeptAndPaid::getPaidFeeAmount(int                iStep,
                                double             nZ,
                                double             Z,
                                const DoubleArray& nE,
                                const DoubleArray& E) const {
    double F = nZ * Z * paidBondFee * feeTimeFactor[iStep];
    for(int iAsset=0; iAsset<paidEquityFee.size(); iAsset++) {
        double fee = paidEquityFee[iAsset];
        if (hasPaidNotDeductedFee) {
            fee += paidNotDeductedEquityFee[iAsset];
        }
        F += nE[iAsset] * E[iAsset] * fee * feeTimeFactor[iStep];
    }
    return F;
}

double SPIFeesKeptAndPaid::getContingentPaidFeeAmount(int                iStep,
                                double             nZ,
                                double             Z,
                                const DoubleArray& nE,
                                const DoubleArray& E) const {
    double F = nZ * Z * contingentPaidBondFee * feeTimeFactor[iStep];
    for(int iAsset=0; iAsset<contingentPaidEquityFee.size(); iAsset++) {
        F += nE[iAsset] * E[iAsset] * contingentPaidEquityFee[iAsset] * feeTimeFactor[iStep];
    }
    return F;
}

bool SPIFeesKeptAndPaid::isThresholdBreached(double level) const {
    return Maths::isPositive(level-contingentThreshold);
}

bool SPIFeesKeptAndPaid::hasContingentFees() const {
    return hasContingentFee;
}

double SPIFeesKeptAndPaid::getFeeAtMin(const DoubleArray& exposureUpperBound,
                           double             equityExposureMin) const {
    double FEmin;
    if (totalEquityFee.size()==1) {
        FEmin = totalEquityFee[0];
    } else { // must be 2
        FEmin = (exposureUpperBound[1] * totalEquityFee[1] +
                 (1. - exposureUpperBound[1]) * totalEquityFee[0]);
    }
    return equityExposureMin*FEmin + (1.-equityExposureMin)*(keptBondFee+paidBondFee);
}

const DateTimeArray& SPIFeesKeptAndPaid::getNotificationDates(const DateTimeArray& rebalDates) const {
    // check paidFeePaymentDates subset of rebal dates
    if (!DateTime::isSubset(rebalDates, paidFeePaymentDates)) {
        throw ModelException("SPIFeesKeptAndPaid::getNotificationDates",
                             "Fee payment dates must be subset of rebalanceDates");
    }
    return paidFeePaymentDates;
}

DateTimeArray SPIFeesKeptAndPaid::getPaymentDates() const {
    return paidFeeSettlementDates;
}

DateTimeArray SPIFeesKeptAndPaid::getEssentialDates() const {
    return paidFeePaymentDates;
}

const DateTime SPIFeesKeptAndPaid::getFeePayDate(const DateTime& rebalDate) const {
    return feeSettle->settles(rebalDate, 0/*asset*/);
}

// post getMarket validation (handles settlement of fees)
void SPIFeesKeptAndPaid::Validate(const InstrumentSettlement* instSettle) {
    // sort out the settlement
    paidFeeSettlementDates = DateTimeArray(paidFeePaymentDates.size());
    if (hasSettlePeriod) {
        HolidayConstSP settleHols = IFeesSPI::getSettleHolidays(instSettle);
        feeSettle = CashSettlePeriodSP(new CashSettlePeriod(settlePeriod, settleHols.get()));
    } else {
        // this has to work cos we validated in SPI validatePop2Object
        CashSettlePeriod* sett = dynamic_cast<CashSettlePeriod*>(copy(instSettle));
        feeSettle = CashSettlePeriodSP(sett);
    }
    for(int j=0; j<paidFeePaymentDates.size(); j++) {
        paidFeeSettlementDates[j] = feeSettle->settles(paidFeePaymentDates[j],
                                                0/*asset*/);
    }
}

IObject* SPIFeesKeptAndPaid::defaultSPIFeesKeptAndPaid(){
    return new SPIFeesKeptAndPaid();
}

/** Invoked when Class is 'loaded' */
void SPIFeesKeptAndPaid::load(CClassSP& clazz){
    REGISTER(SPIFeesKeptAndPaid, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(IFeesSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPIFeesKeptAndPaid);
    FIELD(keptBondFee,             "keptBondFee");
    FIELD(paidBondFee,             "paidBondFee");
    FIELD(keptEquityFee,           "keptEquityFee");
    FIELD(paidEquityFee,           "paidEquityFee");
    FIELD(totalEquityFee,          "totalEquityFee");
    FIELD_MAKE_TRANSIENT(totalEquityFee);
    FIELD(contingentKeptBondFee,             "contingentKeptBondFee");
    FIELD(contingentPaidBondFee,             "contingentPaidBondFee");
    FIELD(contingentKeptEquityFee,           "contingentKeptEquityFee");
    FIELD(contingentPaidEquityFee,           "contingentPaidEquityFee");
    FIELD(totalContingentEquityFee,          "totalContingentEquityFee");
    FIELD(contingentThreshold,          "contingentThreshold");
    FIELD(hasContingentFee,          "Has Contingent Fee?");
    FIELD(hasPaidNotDeductedFee,     "Has Paid Not Deducted Fee?");
    FIELD(paidNotDeductedEquityFee,  "Paid Not Deducted Equity Fee");
    FIELD_MAKE_TRANSIENT(totalContingentEquityFee);
    FIELD(paidFeePaymentDates,     "paidFeePaymentDates");
    FIELD(paidFeeSettlementDates,     "the actual settlement dates");
    FIELD(feeSettle, "The settlement object for the fees");
    FIELD(hasSettlePeriod, "Do fees have their own settlement");
    FIELD(settlePeriod, "#Days settlement for paid fees");
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(keptBondFee);
    FIELD_MAKE_OPTIONAL(paidBondFee);
    FIELD_MAKE_OPTIONAL(keptEquityFee);
    FIELD_MAKE_OPTIONAL(paidEquityFee);
    FIELD_MAKE_OPTIONAL(contingentKeptBondFee);
    FIELD_MAKE_OPTIONAL(contingentPaidBondFee);
    FIELD_MAKE_OPTIONAL(contingentKeptEquityFee);
    FIELD_MAKE_OPTIONAL(contingentPaidEquityFee);
    FIELD_MAKE_OPTIONAL(contingentThreshold);
    FIELD_MAKE_OPTIONAL(hasContingentFee);
    FIELD_MAKE_OPTIONAL(hasPaidNotDeductedFee);
    FIELD_MAKE_OPTIONAL(paidNotDeductedEquityFee);
    FIELD_MAKE_OPTIONAL(paidFeePaymentDates);
    FIELD_MAKE_OPTIONAL(hasSettlePeriod);
    FIELD_MAKE_OPTIONAL(settlePeriod);
    FIELD_MAKE_TRANSIENT(paidFeeSettlementDates);
    FIELD_MAKE_TRANSIENT(feeSettle);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPIFeesKeptAndPaid::TYPE = CClass::registerClassLoadMethod(
    "SPIFeesKeptAndPaid", typeid(SPIFeesKeptAndPaid), SPIFeesKeptAndPaid::load);

/*****************************************************************************/
// Fee rate is based on notional, not basket.
// So no split across equity/bond at least
SPIFeesPerNotional::SPIFeesPerNotional(): SPIInterfaceIMS(TYPE),
        keptFee(0.0), paidFee(0.0),
        contingentKeptFee(0.0),  contingentPaidFee(0.0),
        contingentThreshold(0.0), hasContingentFee(false),
        hasSettlePeriod(false), settlePeriod(0),
		hasPaidNotDeductedFee(false) {} // for reflection

IFeesSPISP SPIFeesPerNotional::getFeesSPI() {
    static const string routine = "SPIFeesPerNotional::getFeesSPI";
    // essentially just returns itself
    IFeesSPISP theFees = IFeesSPISP(this, NullDeleter()); // FIXME can be a problem
    return theFees;
}

const string SPIFeesPerNotional::feeType() {
    return SPI_FEES_TYPE_PNL;
}

// validation
void SPIFeesPerNotional::validatePop2Object(){
    static const string routine = "SPIFeesPerNotional::validatePop2Object";
    // If any non-zero paid fee require at least one fee payment date
    if (!Maths::isZero(paidFee) && paidFeePaymentDates.size()<1) {
        isOK = false;
        err = "Non-zero paid fees so require at least one paid fee payment date!";
        return;
    }
    if (hasContingentFee &&
            !Maths::isZero(contingentPaidFee) && paidFeePaymentDates.size()<1) {
        isOK = false;
        err = "Non-zero contingent paid fees so require at least one paid fee payment date!";
        return;
    }
	if (hasPaidNotDeductedFee) {
		bool havePaidFee = false;
		for(int i=0; i<paidNotDeductedEquityFee.size(); i++) {
			if (!Maths::isZero(paidNotDeductedEquityFee[i])) {
				havePaidFee = true;
			}
		}
		if (havePaidFee && paidFeePaymentDates.size()<1) {
			isOK = false;
			err = "Non-zero paid not deducted fees so require at least one paid fee payment date!";
		}
	}
}

void SPIFeesPerNotional::init(int                  numAlgAssets,
          const DateTimeArray& rebalDates,
          int                  iStepFirstRebal,
          double               Basis,
          ILoanCostSPI*  loanCost,
          const YieldCurve*    yc) {
    if (hasPaidNotDeductedFee &&
            paidNotDeductedEquityFee.size() != numAlgAssets) {
        throw ModelException("SPIFeesPerNotional::init",
                             "Require same number of paid not deducted equity fees (" +
                             Format::toString(paidNotDeductedEquityFee.size()) +
                             ") as assets (" +
                             Format::toString(numAlgAssets) + ")");
    }
    // build feeTimeFactor
    feeTimeFactor = DoubleArray(rebalDates.size(), 0.0);
    for(int i=iStepFirstRebal+1; i<rebalDates.size(); i++) {
        feeTimeFactor[i] = (rebalDates[i].getDate()-rebalDates[i-1].getDate())/Basis;
    }
}

void SPIFeesPerNotional::refresh(int iFirstFutureStep) {}

double SPIFeesPerNotional::getFeeAmount(int                iStep,
                    double             nZ,
                    double             Z,
                    const DoubleArray& nE,
                    const DoubleArray& E) const {

    double F = (keptFee+paidFee) * feeTimeFactor[iStep];
    return F;
}

// when we know it's fee on notional ... i.e. in bond floor calcs
double SPIFeesPerNotional::getFeeAmount(int iStep) const {
    double F = (keptFee+paidFee) * feeTimeFactor[iStep];
    return F;
}

double SPIFeesPerNotional::getContingentFeeAmount(int                iStep,
                    double             nZ,
                    double             Z,
                    const DoubleArray& nE,
                    const DoubleArray& E) const {

    double F = (contingentKeptFee+contingentPaidFee) * feeTimeFactor[iStep];
    return F;
}

double SPIFeesPerNotional::getPaidFeeAmount(int                iStep,
                                double             nZ,
                                double             Z,
                                const DoubleArray& nE,
                                const DoubleArray& E) const {
    double F = paidFee * feeTimeFactor[iStep];

	if (hasPaidNotDeductedFee) {
		for(int iAsset=0; iAsset<paidNotDeductedEquityFee.size(); iAsset++) {
			F += nE[iAsset] * E[iAsset] * paidNotDeductedEquityFee[iAsset]
							* feeTimeFactor[iStep];
		}
	}
    return F;
}

double SPIFeesPerNotional::getContingentPaidFeeAmount(int                iStep,
                                double             nZ,
                                double             Z,
                                const DoubleArray& nE,
                                const DoubleArray& E) const {
    return (contingentPaidFee * feeTimeFactor[iStep]);
}

bool SPIFeesPerNotional::isThresholdBreached(double level) const {
    return Maths::isPositive(level-contingentThreshold);
}

bool SPIFeesPerNotional::hasContingentFees() const {
    return hasContingentFee;
}

double SPIFeesPerNotional::getFeeAtMin(const DoubleArray& exposureUpperBound,
                           double             equityExposureMin) const {
      throw ModelException("SPIFeesPerNotional::getFeeAtMin",
                            "Function not implemented - bond floor doesn't use fee at min for fees per notional");
//        return (keptFee+paidFee);
}

const DateTimeArray& SPIFeesPerNotional::getNotificationDates(const DateTimeArray& rebalDates) const {
    // check paidFeePaymentDates subset of rebal dates
    if (!DateTime::isSubset(rebalDates, paidFeePaymentDates)) {
        throw ModelException("SPIFeesPerNotional::getNotificationDates",
                             "Fee payment dates must be subset of rebalanceDates");
    }
    return paidFeePaymentDates;
}

DateTimeArray SPIFeesPerNotional::getPaymentDates() const {
    return paidFeeSettlementDates;
}

DateTimeArray SPIFeesPerNotional::getEssentialDates() const {
    return paidFeePaymentDates;
}

const DateTime SPIFeesPerNotional::getFeePayDate(const DateTime& rebalDate) const {
    return feeSettle->settles(rebalDate, 0/*asset*/);
}

// post getMarket validation (handles settlement of fees)
void SPIFeesPerNotional::Validate(const InstrumentSettlement* instSettle) {
    // sort out the settlement
    paidFeeSettlementDates = DateTimeArray(paidFeePaymentDates.size());
    if (hasSettlePeriod) {
        HolidayConstSP settleHols = IFeesSPI::getSettleHolidays(instSettle);
        feeSettle = CashSettlePeriodSP(new CashSettlePeriod(settlePeriod, settleHols.get()));
    } else {
        // this has to work cos we validated in SPI validatePop2Object
        CashSettlePeriod* sett = dynamic_cast<CashSettlePeriod*>(copy(instSettle));
        feeSettle = CashSettlePeriodSP(sett);
    }
    for(int j=0; j<paidFeePaymentDates.size(); j++) {
        paidFeeSettlementDates[j] = feeSettle->settles(paidFeePaymentDates[j],
                                                0/*asset*/);
    }
}

IObject* SPIFeesPerNotional::defaultSPIFeesPerNotional(){
    return new SPIFeesPerNotional();
}

/** Invoked when Class is 'loaded' */
void SPIFeesPerNotional::load(CClassSP& clazz){
    REGISTER(SPIFeesPerNotional, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(IFeesSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPIFeesPerNotional);
    FIELD(keptFee,             "Kept Fee");
    FIELD(paidFee,             "Paid Fee");
    FIELD(contingentKeptFee,   "Contingent Kept Fee");
    FIELD(contingentPaidFee,   "Contingent Paid Fee");
    FIELD(contingentThreshold, "Contingent Threshold");
    FIELD(hasContingentFee, "Contingent Threshold");
    FIELD(paidFeePaymentDates,     "paidFeePaymentDates");
    FIELD(hasPaidNotDeductedFee,     "Has Paid Not Deducted Fee?");
    FIELD(paidNotDeductedEquityFee,  "Paid Not Deducted Equity Fee");
    FIELD(paidFeeSettlementDates,     "the actual settlement dates");
    FIELD(feeSettle, "The settlement object for the fees");
    FIELD(hasSettlePeriod, "Do gfees have their own settlement");
    FIELD(settlePeriod, "#Days settlement for paid fees");
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(keptFee);
    FIELD_MAKE_OPTIONAL(paidFee);
    FIELD_MAKE_OPTIONAL(paidFeePaymentDates);
    FIELD_MAKE_OPTIONAL(contingentKeptFee);
    FIELD_MAKE_OPTIONAL(contingentPaidFee);
    FIELD_MAKE_OPTIONAL(contingentThreshold);
    FIELD_MAKE_OPTIONAL(hasContingentFee);
    FIELD_MAKE_OPTIONAL(hasPaidNotDeductedFee);
    FIELD_MAKE_OPTIONAL(paidNotDeductedEquityFee);
    FIELD_MAKE_OPTIONAL(hasSettlePeriod);
    FIELD_MAKE_OPTIONAL(settlePeriod);
    FIELD_MAKE_TRANSIENT(feeSettle);
    FIELD_MAKE_TRANSIENT(paidFeeSettlementDates);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPIFeesPerNotional::TYPE = CClass::registerClassLoadMethod(
    "SPIFeesPerNotional", typeid(SPIFeesPerNotional), SPIFeesPerNotional::load);

/*****************************************************************************/
// Fee based on basket value and a vriable libor rate with some kept and some paid

SPIFeesLibor::SPIFeesLibor(): SPIInterfaceIMS(TYPE),
        hasSettlePeriod(false), settlePeriod(0) {} // for reflection

IFeesSPISP SPIFeesLibor::getFeesSPI() {
    static const string routine = "SPIFeesLibor::getFeesSPI";
    // essentially just returns itself
    IFeesSPISP theFees = IFeesSPISP(this, NullDeleter()); // FIXME can be a problem
    return theFees;
}

const string SPIFeesLibor::feeType() {
    return SPI_FEES_TYPE_LIB;
}

// validation
void SPIFeesLibor::validatePop2Object(){
    static const string routine = "SPIFeesLibor::validatePop2Object";
    if (keptEquityFeeMult.size() != paidEquityFeeMult.size()) {
        isOK = false;
        err = "Kept/paid equity fee multiplier arrays should be of equal length!";
        return;
    }
    if ((!keptEquityFeeSpread.get() || (*keptEquityFeeSpread).size() == 0)) {
        keptEquityFeeSpread = DoubleArraySP(new DoubleArray(keptEquityFeeMult.size(), 0.0));
    }
    if ((!paidEquityFeeSpread.get() || (*paidEquityFeeSpread).size() == 0)) {
        paidEquityFeeSpread = DoubleArraySP(new DoubleArray(paidEquityFeeMult.size(), 0.0));
    }
    if (keptEquityFeeMult.size() != (*keptEquityFeeSpread).size()) {
        isOK = false;
        err = "Kept fee and spread arrays should be equal length!";
        return;
    }
    if (paidEquityFeeMult.size() != (*paidEquityFeeSpread).size()) {
        isOK = false;
        err = "Paid fee and spread arrays should be equal length!";
        return;
    }

    bool havePaidFee = false;
    bool haveLiborFee = false;
    for(int i=0; i<keptEquityFeeMult.size(); i++) {
        if (Maths::isPositive(keptEquityFeeMult[i])) {
            haveLiborFee = true;
        } else if (Maths::isNegative(keptEquityFeeMult[i])) {
            isOK = false;
            err = "Libor kept fee multiplier for asset " +
                        Format::toString(i+1) + " is negative";
            return;
        }
        if (Maths::isPositive(paidEquityFeeMult[i])) {
            havePaidFee = true;
            haveLiborFee = true;
        } else if (Maths::isNegative(paidEquityFeeMult[i])) {
            isOK = false;
            err = "Libor paid fee multiplier for asset " +
                        Format::toString(i+1) + " is negative";
            return;
        }
    }
    if (!haveLiborFee) {
        isOK = false;
        err = "All libor fee multipliers are zero";
        return;
    }

    // If any non-zero paid fee require at least one fee payment date
    if (havePaidFee && paidFeePaymentDates.size()<1) {
        isOK = false;
        err = "Non-zero paid fees so require at least one paid fee payment date!";
        return;
    }
}

// by using the LoanCost object to load up the libor rates
// we can rely on the loan cost handling the theta roll etc
void SPIFeesLibor::init(int                  numAlgAssets,
                        const DateTimeArray& rebalDates,
                        int                  iStepFirstRebal,
                        double               Basis,
                        ILoanCostSPI*  loanCost,
                        const YieldCurve*    yc) {
    // check numAlgAssets agrees with number of equity fees
    if (keptEquityFeeMult.size() != numAlgAssets) {
        throw ModelException("SPIFeesLibor::init",
                             "Require same number of equity fees (" +
                             Format::toString(keptEquityFeeMult.size()) +
                             ") as assets (" +
                             Format::toString(numAlgAssets) + ")");
    }
    // build feeFactor - per asset, per time point
    // this is Libor rate * multiplier plus spread
    // and includes the year fraction as well
    feeFactors = DoubleMatrix(numAlgAssets, rebalDates.size());
    paidFeeFactors = DoubleMatrix(numAlgAssets, rebalDates.size());
    timeFactors = DoubleArray(rebalDates.size(), 0.0);
    for(int i=iStepFirstRebal+1; i<rebalDates.size(); i++) {
        double liborRate = loanCost->getLiborRate(yc, rebalDates[i-1],
                                                  rebalDates[i], Basis);
        timeFactors[i] = (rebalDates[i].getDate()-rebalDates[i-1].getDate())/Basis;
        for (int iAsset = 0; iAsset < numAlgAssets; ++iAsset) {
            // libor bit
            feeFactors[iAsset][i] = liborRate * paidEquityFeeMult[iAsset];
            // then the spread
            feeFactors[iAsset][i] += (*paidEquityFeeSpread)[iAsset];
            paidFeeFactors[iAsset][i] = feeFactors[iAsset][i] * timeFactors[i]; 

            // now the same for the kept fees
            feeFactors[iAsset][i] += liborRate * keptEquityFeeMult[iAsset];
            feeFactors[iAsset][i] += (*keptEquityFeeSpread)[iAsset];
            feeFactors[iAsset][i] *= timeFactors[i];
            if(Maths::isNegative(feeFactors[iAsset][i])) {
                throw ModelException("SPIFeesLibor::init",
                        "Negative Libor fee generated for asset " +
                        Format::toString(iAsset+1) +
                        " between dates " +
                        rebalDates[i-1].toString() +
                        " and " + rebalDates[i].toString());
            }
        }
    }
    this->numAlgAssets = numAlgAssets;
    this->loanCost = loanCost;
}
    
void SPIFeesLibor::refresh(int iFirstFutureStep) {
    // update each payoff() call for stoch rates
    // future only
    for(int i=iFirstFutureStep+1; i<timeFactors.size(); i++) {
        double liborRate = loanCost->getFutureLiborRate(i);
        for (int iAsset = 0; iAsset < numAlgAssets; ++iAsset) {
            feeFactors[iAsset][i] = liborRate * paidEquityFeeMult[iAsset];
            // then the spread 
            feeFactors[iAsset][i] += (*paidEquityFeeSpread)[iAsset];;
            paidFeeFactors[iAsset][i] = feeFactors[iAsset][i] * timeFactors[i]; 
            
            // now the same for the kept fees
            feeFactors[iAsset][i] += liborRate * keptEquityFeeMult[iAsset];
            feeFactors[iAsset][i] += (*keptEquityFeeSpread)[iAsset];
            feeFactors[iAsset][i] *= timeFactors[i];
            if(Maths::isNegative(feeFactors[iAsset][i])) {
                throw ModelException("SPIFeesLibor::init",
                        "Negative Libor fee generated for asset " +
                        Format::toString(iAsset+1) +
                        " at rebalance date " +
                        Format::toString(i+1));
            }
        }
    }

}

double SPIFeesLibor::getFeeAmount(int                iStep,
                    double             nZ,
                    double             Z,
                    const DoubleArray& nE,
                    const DoubleArray& E) const {

    double F = 0.0;
    for(int iAsset=0; iAsset<keptEquityFeeMult.size(); iAsset++) {
        F += nE[iAsset] * E[iAsset] * feeFactors[iAsset][iStep];
    }
    return F;
}

double SPIFeesLibor::getContingentFeeAmount(int                iStep,
                    double             nZ,
                    double             Z,
                    const DoubleArray& nE,
                    const DoubleArray& E) const {
    return 0.0;
}

double SPIFeesLibor::getPaidFeeAmount(int                iStep,
                                double             nZ,
                                double             Z,
                                const DoubleArray& nE,
                                const DoubleArray& E) const {
    double F = 0.0;
    for(int iAsset=0; iAsset<paidEquityFeeMult.size(); iAsset++) {
        F += nE[iAsset] * E[iAsset] * paidFeeFactors[iAsset][iStep];
    }
    return F;
}

double SPIFeesLibor::getContingentPaidFeeAmount(int                iStep,
                                double             nZ,
                                double             Z,
                                const DoubleArray& nE,
                                const DoubleArray& E) const {
    return 0.0;
}

bool SPIFeesLibor::isThresholdBreached(double level) const {
    return true;
}

bool SPIFeesLibor::hasContingentFees() const {
    return false;
}

// Not sure what to do here in presence of Libor fees.
// do we adjust the bond floor for the fixed spread fees?? disable for now
// as we validate against spreads with min equity exposure > 0
double SPIFeesLibor::getFeeAtMin(const DoubleArray& exposureUpperBound,
                           double             equityExposureMin) const {
//    double FEmin = 0.0;
//    if (totalEquityFee.size()==1) {
//        FEmin = keptEquityFeeSpread[0] + paidEquityFeeSpread[0];
//    } else { // must be 2
//        FEmin = (exposureUpperBound[1] * (keptEquityFeeSpread[1] + paidEquityFeeSpread[1]) +
//                 (1. - exposureUpperBound[1]) * (keptEquityFeeSpread[0] + paidEquityFeeSpread[0]));
//    }
//    return equityExposureMin*FEmin;
    return 0.0;
}

const DateTimeArray& SPIFeesLibor::getNotificationDates(const DateTimeArray& rebalDates) const {
    // check paidFeePaymentDates subset of rebal dates
    if (!DateTime::isSubset(rebalDates, paidFeePaymentDates)) {
        throw ModelException("SPIFeesLibor::getNotificationDates",
                             "Fee payment dates must be subset of rebalanceDates");
    }
    return paidFeePaymentDates;
}

DateTimeArray SPIFeesLibor::getPaymentDates() const {
    return paidFeeSettlementDates;
}

DateTimeArray SPIFeesLibor::getEssentialDates() const {
    return paidFeePaymentDates;
}

const DateTime SPIFeesLibor::getFeePayDate(const DateTime& rebalDate) const {
    return feeSettle->settles(rebalDate, 0/*asset*/);
}

// post getMarket validation (handles settlement of fees)
void SPIFeesLibor::Validate(const InstrumentSettlement* instSettle) {
    // sort out the settlement
    paidFeeSettlementDates = DateTimeArray(paidFeePaymentDates.size());
    if (hasSettlePeriod) {
        HolidayConstSP settleHols = IFeesSPI::getSettleHolidays(instSettle);
        feeSettle = CashSettlePeriodSP(new CashSettlePeriod(settlePeriod, settleHols.get()));
    } else {
        // this has to work cos we validated in SPI validatePop2Object
        CashSettlePeriod* sett = dynamic_cast<CashSettlePeriod*>(copy(instSettle));
        feeSettle = CashSettlePeriodSP(sett);
    }
    for(int j=0; j<paidFeePaymentDates.size(); j++) {
        paidFeeSettlementDates[j] = feeSettle->settles(paidFeePaymentDates[j],
                                                0/*asset*/);
    }
}

// overriden version for Libor fees. Stops spreadd if min equity exposure
void SPIFeesLibor::crossValidate(double equityExposureMin) const {
    if (Maths::isPositive(equityExposureMin)) {
        bool valid = true;
        if (keptEquityFeeSpread.get() && (*keptEquityFeeSpread).size() > 0) {
            for (int i=0; i<(*keptEquityFeeSpread).size(); ++i) {
                if(Maths::isPositive((*keptEquityFeeSpread)[i])) {
                    valid = false;
                }
            }
        }
        if (paidEquityFeeSpread.get() && (*paidEquityFeeSpread).size() > 0) {
            for (int i=0; i<(*paidEquityFeeSpread).size(); ++i) {
                if(Maths::isPositive((*paidEquityFeeSpread)[i])) {
                    valid = false;
                }
            }
        }
        if (!valid) {
            throw ModelException("SPIFeesLibor::crossValidate",
                                "Cannot have spread on Libor fees together with "
                                "a positive minimum equity exposure.");
        }
    }
}


IObject* SPIFeesLibor::defaultSPIFeesLibor(){
    return new SPIFeesLibor();
}
/** Invoked when Class is 'loaded' */
void SPIFeesLibor::load(CClassSP& clazz){
    REGISTER(SPIFeesLibor, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(IFeesSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPIFeesLibor);
    FIELD(keptEquityFeeMult,        "Kept Equity Fee Multiplier");
    FIELD(paidEquityFeeMult,        "Paid Equity Fee Multiplier");
    FIELD(keptEquityFeeSpread,     "keptEquityFeeSpread");
    FIELD(paidEquityFeeSpread,     "paidEquityFeeSpread");
    FIELD(paidFeePaymentDates,     "paidFeePaymentDates");
    FIELD(paidFeeSettlementDates,     "the actual settlement dates");
    FIELD(feeSettle, "The settlement object for the fees");
    FIELD(hasSettlePeriod, "Do gfees have their own settlement");
    FIELD(settlePeriod, "#Days settlement for paid fees");
    FIELD(timeFactors,     "timefactors");
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(keptEquityFeeMult);
    FIELD_MAKE_OPTIONAL(paidEquityFeeMult);
    FIELD_MAKE_OPTIONAL(keptEquityFeeSpread);
    FIELD_MAKE_OPTIONAL(paidEquityFeeSpread);
    FIELD_MAKE_OPTIONAL(paidFeePaymentDates);
    FIELD_MAKE_OPTIONAL(hasSettlePeriod);
    FIELD_MAKE_OPTIONAL(settlePeriod);
    FIELD_MAKE_TRANSIENT(feeSettle);
    FIELD_MAKE_TRANSIENT(paidFeeSettlementDates);
    FIELD_MAKE_TRANSIENT(timeFactors);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPIFeesLibor::TYPE = CClass::registerClassLoadMethod(
    "SPIFeesLibor", typeid(SPIFeesLibor), SPIFeesLibor::load);

/***********************************************************************************/

const string SPIFeesWrapper::feeType() {
    if (SPIFeesType.empty()){
        throw ModelException("SPIFeesWrapper::feeType", "Blank SPIFeesType specified!");
    }
    return SPIFeesType;
}

IFeesSPISP SPIFeesWrapper::getFeesSPI() {
    static const string routine = "SPIFeesWrapper::getFeesSPI";
    IFeesSPISP theFees;

    if (SPIFeesType.empty()){
        throw ModelException(routine, "Blank SPIFeesType specified!");
    }
    if (SPIFeesType==SPI_FEES_TYPE_STD) {
        if (!feesStd.get()) {
            throw ModelException(routine, "Expected feesStd but none supplied!");
        }
        theFees = IFeesSPISP(feesStd.get(), NullDeleter());
    } else if (SPIFeesType==SPI_FEES_TYPE_KNP) {
        if (!feesKeptAndPaid.get()) {
            throw ModelException(routine, "Expected feesKeptAndPaid but none supplied!");
        }
        theFees = IFeesSPISP(feesKeptAndPaid.get(), NullDeleter());
    } else if (SPIFeesType==SPI_FEES_TYPE_PNL) {
        if (!feesPerNotional.get()) {
            throw ModelException(routine, "Expected feesPerNotional but none supplied!");
        }
        theFees = IFeesSPISP(feesPerNotional.get(), NullDeleter());
    } else if (SPIFeesType==SPI_FEES_TYPE_LIB) {
        if (!feesLibor.get()) {
            throw ModelException(routine, "Expected feesLibor but none supplied!");
        }
        theFees = IFeesSPISP(feesLibor.get(), NullDeleter());
    } else {
        throw ModelException(routine, "Unrecognised SPIFeesType " + SPIFeesType +
                             ". Expected " + SPI_FEES_TYPE_STD + ", " +
                             SPI_FEES_TYPE_KNP + ", " +
                             SPI_FEES_TYPE_PNL + " or " +
                             SPI_FEES_TYPE_LIB);
    }

	SPIInterfaceIMS* actualFees = dynamic_cast<SPIInterfaceIMS*>(theFees.get());
	if (!actualFees) {
        throw ModelException(routine, "Internal error in SPI fee validation");
	}
    if (!actualFees->isValid()) {
        throw ModelException(routine, "Invalid SPIFees (type " + SPIFeesType +
                             ") : " + actualFees->errString());
    }
    return theFees;
}

// validation
void SPIFeesWrapper::validatePop2Object(){
    static const string routine = "SPIFeesWrapper::validatePop2Object";

}

/** revert the convert method of SPIFeesWrapperSuper
    needed for PYRAMID and IMS */
IObjectSP SPIFeesWrapper::revert(const string& interfaceType) const {
    static const string method = "SPIFeesWrapper::revert";
    try {
        if (interfaceType != IRevertTypeConvert::PYRAMID) {
            throw ModelException( method,
                                  "Cannot convert a SPIFeesWrapper for"
                                  " the interface " + interfaceType);
        }

        IObjectSP wrapper(new SPIFeesWrapperSuper(SPIFeesType,
                                                  feesStd,
                                                  feesKeptAndPaid,
                                                  feesPerNotional,
                                                  feesLibor));
        return wrapper;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Invoked when Class is 'loaded' */
void SPIFeesWrapper::load(CClassSP& clazz){
    REGISTER(SPIFeesWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IFeesSPIInterface);
    IMPLEMENTS(IRevertTypeConvert);
    EMPTY_SHELL_METHOD(defaultSPIFeesWrapper);
    FIELD(SPIFeesType, "SPIFeesType");
    FIELD(feesStd,  "feesStd");
    FIELD_MAKE_OPTIONAL(feesStd);
    FIELD(feesKeptAndPaid,  "feesKeptAndPaid");
    FIELD_MAKE_OPTIONAL(feesKeptAndPaid);
    FIELD(feesPerNotional,  "feesPerNotional");
    FIELD_MAKE_OPTIONAL(feesPerNotional);
    FIELD(feesLibor,  "feesLibor");
    FIELD_MAKE_OPTIONAL(feesLibor);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

// for reflection
SPIFeesWrapper::SPIFeesWrapper(): CObject(TYPE){}

// constructor for use with the SPIFeesWrapperSuper class
SPIFeesWrapper::SPIFeesWrapper(string type, SPIFeesStdSP std,
                SPIFeesKeptAndPaidSP  keptAndPaid,
                SPIFeesPerNotionalSP  perNotional,
                SPIFeesLiborSP        libor) :
            CObject(TYPE), SPIFeesType(type), feesStd(std),
            feesKeptAndPaid(keptAndPaid),
            feesPerNotional(perNotional),
            feesLibor(libor) {}

IObject* SPIFeesWrapper::defaultSPIFeesWrapper(){
    return new SPIFeesWrapper();
}
CClassConstSP const SPIFeesWrapper::TYPE = CClass::registerClassLoadMethod(
    "SPIFeesWrapper", typeid(SPIFeesWrapper), SPIFeesWrapper::load);

// a wrapper for the SPIFeesWrapper so we can split the IMS interface into 2
SPIFeesWrapperSuper::SPIFeesWrapperSuper():CObject(TYPE){}

/** create a proper SPIFeesWrapper */
void SPIFeesWrapperSuper::convert(IObjectSP&    object,
                                  CClassConstSP requiredType) const {
    static const string method = "SPIFeesWrapperSuper::convert";
    try {
        // note either converting to old style 'abstraction'
        // or the new proper abstraction
        if (requiredType != SPIFeesWrapper::TYPE
            && requiredType != IFeesSPIInterface::TYPE) {
            throw ModelException(method,
                                 "Cannot convert a SPIFeesWrapperSuper into "
                                 "object of type "+requiredType->getName());
        }
        object = IObjectSP(new SPIFeesWrapper(spiFeesType,
                                              feesStd,
                                              feesKeptAndPaid,
                                              feesPerNotional,
                                              feesLibor));
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Invoked when Class is 'loaded' */
void SPIFeesWrapperSuper::load(CClassSP& clazz){
    REGISTER(SPIFeesWrapperSuper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ITypeConvert);
    EMPTY_SHELL_METHOD(defaultSPIFeesWrapperSuper);
    FIELD(spiFeesType, "SPIFeesType");
    FIELD(feesStd,  "feesStd");
    FIELD_MAKE_OPTIONAL(feesStd);
    FIELD(feesKeptAndPaid,  "feesKeptAndPaid");
    FIELD_MAKE_OPTIONAL(feesKeptAndPaid);
    FIELD(feesPerNotional,  "feesPerNotional");
    FIELD_MAKE_OPTIONAL(feesPerNotional);
    FIELD(feesLibor,  "feesLibor");
    FIELD_MAKE_OPTIONAL(feesLibor);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

// constructor for use with the SPIFeesWrapper class
SPIFeesWrapperSuper::SPIFeesWrapperSuper(string type, SPIFeesStdSP std,
                                         SPIFeesKeptAndPaidSP  keptAndPaid,
                                         SPIFeesPerNotionalSP  perNotional,
                                         SPIFeesLiborSP        libor) :
    CObject(TYPE),
    spiFeesType(type),
    feesStd(std),
    feesKeptAndPaid(keptAndPaid),
    feesPerNotional(perNotional),
    feesLibor(libor) {}

IObject* SPIFeesWrapperSuper::defaultSPIFeesWrapperSuper(){
    return new SPIFeesWrapperSuper();
}

CClassConstSP const SPIFeesWrapperSuper::TYPE = CClass::registerClassLoadMethod(
    "SPIFeesWrapperSuper", typeid(SPIFeesWrapperSuper), SPIFeesWrapperSuper::load);

DRLIB_END_NAMESPACE
