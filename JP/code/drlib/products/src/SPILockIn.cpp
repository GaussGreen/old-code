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
#include "edginc/SPILockIn.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/
// this is the external interface for LockIn so that users can bolt in 
// any LockIn type they want - note this includes the SPILockInWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface ILockInSPI as soon as possible
void ILockInSPIInterface::load(CClassSP& clazz) {
    REGISTER_INTERFACE(ILockInSPIInterface, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}
CClassConstSP const ILockInSPIInterface::TYPE = CClass::registerInterfaceLoadMethod(
    "ILockInSPIInterface", typeid(ILockInSPIInterface), ILockInSPIInterface::load);


/*****************************************************************************/

double SPILockInStd::getInitialLockIn() const {
    return initialLockIn;
}

const ILockInSPI* SPILockInStd::getLockInSPINoInit() const {
    return this;
}

ILockInSPI* SPILockInStd::getLockInSPI(const DateTimeArray& rebalanceDates,
                                 const DateTime&      lastRebalDate) {
    // just call the init method
    init(rebalanceDates, lastRebalDate);

    return this;
};

void SPILockInStd::init(const DateTimeArray& sampleDates,
          const DateTime&      lastRebalDate) {
    static const string method = "SPILockInStd::init";
    if (lockInDates.size()>0 &&
        lockInDates.back()>lastRebalDate) {
        throw ModelException(method,
                             "Cannot have lock-in date (" + lockInDates.back().toString() + 
                             ") after final rebalance date (" + lastRebalDate.toString() +
                             ")");
    }

    lockInFlag = BoolArray(sampleDates.size(), false);
    int j=0;
    for(int i=0; i<sampleDates.size() && j<lockInDates.size(); i++) {
        if (sampleDates[i] == lockInDates[j]) {
            lockInFlag[i] = true;
            j++;
        }
    }
    // check all lock-in dates are sample dates
    if (j<lockInDates.size()) {
        throw ModelException(method,
                             "Lock-In dates must be a subset of rebalance dates - check lock-in date #"
                             + Format::toString(j+1) + " = " + lockInDates[j].toString());
        
    }
}

void SPILockInStd::apply(double&  BL,
           double   B,
           double   BF,
           int      iStep) const {
    if (lockInFlag[iStep]) {
        BL = Maths::max(BL, lockInPercentage * B);
    }
}

SPILockInStd::SPILockInStd(): SPIInterfaceIMS(TYPE), 
    initialLockIn(0.0), lockInPercentage(0.0), 
	lockInDates(0) {} // for reflection

void SPILockInStd::validatePop2Object(){
    static const string method = "SPILockInStd::validatePop2Object";
    
    if (Maths::isNegative(initialLockIn)) {
        isOK = false;
        err = "Initial Lock-In cannot be negative : " + Format::toString(initialLockIn);
        return;
    }
    if (Maths::isNegative(lockInPercentage)) {
        isOK = false;
        err = "Lock-In Percentage cannot be negative : " + Format::toString(lockInPercentage);
        return;
    }
    DateTime::ensureIncreasing(lockInDates, "Lock-In dates", false /*failIfEmpty*/);
}

DateTimeArray SPILockInStd::getEssentialDates() const {
    return lockInDates;
}

bool SPILockInStd::doesNothing() const {
    // deemed to do nothing if lock in % is 0 - only allow initial lock in
    return (Maths::isZero(lockInPercentage));
}

IObject* SPILockInStd::defaultSPILockInStd(){
    return new SPILockInStd();
}

/** Invoked when Class is 'loaded' */
void SPILockInStd::load(CClassSP& clazz){
    REGISTER(SPILockInStd, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(ILockInSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPILockInStd);
    FIELD(initialLockIn,  "initialLockIn");
    FIELD(lockInPercentage,  "lockInPercentage");
    FIELD(lockInDates,       "lockInDates");
    FIELD(lockInFlag,        "lockInFlag");
    FIELD_MAKE_TRANSIENT(lockInFlag);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(initialLockIn);
    FIELD_MAKE_OPTIONAL(lockInPercentage);
    FIELD_MAKE_OPTIONAL(lockInDates);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPILockInStd::TYPE = CClass::registerClassLoadMethod(
    "SPILockInStd", typeid(SPILockInStd), SPILockInStd::load);

/*****************************************************************************/
double SPILockInExcess::getInitialLockIn() const {
    return initialLockIn;
}

const ILockInSPI* SPILockInExcess::getLockInSPINoInit() const {
    return this;
}

ILockInSPI* SPILockInExcess::getLockInSPI(const DateTimeArray& rebalanceDates,
                                 const DateTime&      lastRebalDate) {
    // just call the init method
    init(rebalanceDates, lastRebalDate);

    return this;
};

void SPILockInExcess::init(const DateTimeArray& sampleDates,
          const DateTime&      lastRebalDate) {
    static const string method = "SPILockInExcess::init";

    if (lockInDates.size()>0 &&
        lockInDates.back()>lastRebalDate) {
        throw ModelException(method,
                             "Cannot have lock-in date (" + lockInDates.back().toString() + 
                             ") after final rebalance date (" + lastRebalDate.toString() +
                             ")");
    }

    lockInFlag = BoolArray(sampleDates.size(), false);
    int j=0;
    for(int i=0; i<sampleDates.size() && j<lockInDates.size(); i++) {
        if (sampleDates[i] == lockInDates[j]) {
            lockInFlag[i] = true;
            j++;
        }
    }
    // check all lock-in dates are sample dates
    if (j<lockInDates.size()) {
        throw ModelException(method,
                             "Lock-In dates must be a subset of rebalance dates - check lock-in date #"
                             + Format::toString(j+1) + " = " + lockInDates[j].toString());
    
    }
}

void SPILockInExcess::apply(double&  BL,
           double   B,
           double   BF,
           int      iStep) const {
    if (lockInFlag[iStep]) {
        BL = Maths::max(BL, BL + lockInPercentage * (B - BL));
    }
}

SPILockInExcess::SPILockInExcess(): SPIInterfaceIMS(TYPE), 
    initialLockIn(0.0), lockInPercentage(0.0), 
	lockInDates(0) {} // for reflection

DateTimeArray SPILockInExcess::getEssentialDates() const {
    return lockInDates;
}

bool SPILockInExcess::doesNothing() const {
    // all types other than standard are always deemed to be doing something
    return false;
}

void SPILockInExcess::validatePop2Object(){
    static const string method = "SPILockInStd::validatePop2Object";
    
    if (Maths::isNegative(initialLockIn)) {
        isOK = false;
        err = "Initial Lock-In cannot be negative : " + Format::toString(initialLockIn);
        return;
    }
    if (Maths::isNegative(lockInPercentage)) {
        isOK = false;
        err = "Lock-In Percentage cannot be negative : " + Format::toString(lockInPercentage);
        return;
    }
    DateTime::ensureIncreasing(lockInDates, "Lock-In dates", false /*failIfEmpty*/);
}

IObject* SPILockInExcess::defaultSPILockInExcess(){
    return new SPILockInExcess();
}

/** Invoked when Class is 'loaded' */
void SPILockInExcess::load(CClassSP& clazz){
    REGISTER(SPILockInExcess, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(ILockInSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPILockInExcess);
    FIELD(initialLockIn,  "initialLockIn");
    FIELD(lockInPercentage,  "lockInPercentage");
    FIELD(lockInDates,       "lockInDates");
    FIELD(lockInFlag,        "lockInFlag");
    FIELD_MAKE_TRANSIENT(lockInFlag);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(initialLockIn);
    FIELD_MAKE_OPTIONAL(lockInPercentage);
    FIELD_MAKE_OPTIONAL(lockInDates);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPILockInExcess::TYPE = CClass::registerClassLoadMethod(
    "SPILockInExcess", typeid(SPILockInExcess), SPILockInExcess::load);

/*****************************************************************************/

double SPILockInSchedule::getInitialLockIn() const {
    return initialLockIn;
}

const ILockInSPI* SPILockInSchedule::getLockInSPINoInit() const {
    return this;
}

ILockInSPI* SPILockInSchedule::getLockInSPI(const DateTimeArray& rebalanceDates,
                                 const DateTime&      lastRebalDate) {
    // just call the init method
    init(rebalanceDates, lastRebalDate);

    return this;
};

void SPILockInSchedule::init(const DateTimeArray& sampleDates,
          const DateTime&      lastRebalDate) {
    static const string method = "SPILockInSchedule::init";

    if (lockIns.size()>0 &&
        lockIns.back().date>lastRebalDate) {
        throw ModelException(method,
                             "Cannot have lock-in date (" + lockIns.back().date.toString() + 
                             ") after final rebalance date (" + lastRebalDate.toString() +
                             ")");
    }

    if (buffers.size()!=lockIns.size()) {
        throw ModelException(method,
                             "Number of buffer levels (" + Format::toString(buffers.size()) +
                             ") must equal number of lock-ins (" + 
                             Format::toString(lockIns.size()) + ")");
    }
    lockInIdx = IntArray(sampleDates.size(), -1);
    int j=0;
    for(int i=0; i<sampleDates.size() && j<lockIns.size(); i++) {
        if (sampleDates[i] == lockIns[j].date) {
            lockInIdx[i] = j;
            j++;
        }
    }
     // check all lock-in dates are sample dates
    if (j<lockIns.size()) {
        throw ModelException(method,
                             "Lock-In dates must be a subset of rebalance dates - check lock-in date #"
                             + Format::toString(j+1) + " = " + lockIns[j].date.toString());
        
    }
}

void SPILockInSchedule::apply(double&  BL,
           double   B,
           double   BF,
           int      iStep) const {
    int idx = lockInIdx[iStep];
    if (idx>=0) { // possible lock-in
        double L = lockIns[idx].amount;
        double buffer = buffers[idx];
        // protect a "buffer" so always some equity
        if ( L*(1.0 + buffer) < BL/BF ) {
            BL = Maths::max(BL, L * B);
        } else if (doBestLockIn) {
            L = BL/(BF * (1.0 + buffer));
            BL = Maths::max(BL, L * B);
        }
    }
}

SPILockInSchedule::SPILockInSchedule(): SPIInterfaceIMS(TYPE), initialLockIn(0.0),
               lockIns(0), buffers(0), 
               doBestLockIn(false){} // for reflection

DateTimeArray SPILockInSchedule::getEssentialDates() const {
    DateTimeArray lockInDates;
    for(int i=0; i<lockIns.size(); i++) {
        lockInDates.push_back(lockIns[i].date);
    }
    return lockInDates;
}

bool SPILockInSchedule::doesNothing() const {
    // all types other than standard are always deemed to be doing something
    return false;
}

void SPILockInSchedule::validatePop2Object(){
    static const string method = "SPILockInSchedule::validatePop2Object";
    int i;

    if (Maths::isNegative(initialLockIn)) {
        isOK = false;
        err = "Initial Lock-In cannot be negative : " + Format::toString(initialLockIn);
        return;
    }
    for(i=0; i<buffers.size(); i++) {
        if (Maths::isNegative(buffers[i])) {
            isOK = false;
            err = "buffers["+Format::toString(i+1)+
                "] is negative! : " + Format::toString(buffers[i]);
            return;
        }
    }
    CashFlow::ensureDatesIncreasing(lockIns,
                                    "Lock-In Schedule", false/*failIfEmpty*/);
    for(i=0; i<lockIns.size(); i++) {
        if (Maths::isNegative(lockIns[i].amount)) {
            isOK = false;
            err = "Lock-In level [" + Format::toString(i+1) + 
                "] cannot be negative : " + Format::toString(lockIns[i].amount);
            return;
        }
    }
}

IObject* SPILockInSchedule::defaultSPILockInSchedule(){
    return new SPILockInSchedule();
}

/** Invoked when Class is 'loaded' */
void SPILockInSchedule::load(CClassSP& clazz){
    REGISTER(SPILockInSchedule, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(ILockInSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPILockInSchedule);
    FIELD(initialLockIn,  "initialLockIn");
    FIELD(lockIns,            "lockIns");
    FIELD(buffers,            "Array of buffer values");
    FIELD(doBestLockIn,       "doBestLockIn");
    FIELD(lockInIdx,          "lockInIdx");
    FIELD_MAKE_TRANSIENT(lockInIdx);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(initialLockIn);
    FIELD_MAKE_OPTIONAL(lockIns);
    FIELD_MAKE_OPTIONAL(buffers);
    FIELD_MAKE_OPTIONAL(doBestLockIn);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPILockInSchedule::TYPE = CClass::registerClassLoadMethod(
    "SPILockInSchedule", typeid(SPILockInSchedule), SPILockInSchedule::load);

/*****************************************************************************/

double SPILockInGrowth::getInitialLockIn() const {
    return initialLockIn;
}

const ILockInSPI* SPILockInGrowth::getLockInSPINoInit() const {
    return this;
}

ILockInSPI* SPILockInGrowth::getLockInSPI(const DateTimeArray& rebalanceDates,
                                 const DateTime&      lastRebalDate) {
    // just call the init method
    init(rebalanceDates, lastRebalDate);

    return this;
};

void SPILockInGrowth::init(const DateTimeArray& sampleDates,
          const DateTime&      lastRebalDate) {
    static const string method = "SPILockInGrowth::init";
    if (lockInDates.size()>0 &&
        lockInDates.back()>lastRebalDate) {
        throw ModelException(method,
                             "Cannot have lock-in date (" + lockInDates.back().toString() + 
                             ") after final rebalance date (" + lastRebalDate.toString() +
                             ")");
    }

    lockInFlag = BoolArray(sampleDates.size(), false);
    int j=0;
    for(int i=0; i<sampleDates.size() && j<lockInDates.size(); i++) {
        if (sampleDates[i] == lockInDates[j]) {
            lockInFlag[i] = true;
            j++;
        }
    }
    // check all lock-in dates are sample dates
    if (j<lockInDates.size()) {
        throw ModelException(method,
                             "Lock-In dates must be a subset of rebalance dates - check lock-in date #"
                             + Format::toString(j+1) + " = " + lockInDates[j].toString());
        
    }
}

void SPILockInGrowth::apply(double&  BL,
           double   B,
           double   BF,
           int      iStep) const {
    // Very much like "Excess" but based around initialLockIn and not latest BL...
    if (lockInFlag[iStep]) {
        BL = Maths::max(BL, initialLockIn + lockInPercentage * (B - initialLockIn));
    }
}

SPILockInGrowth::SPILockInGrowth(): SPIInterfaceIMS(TYPE), 
    initialLockIn(0.0), lockInPercentage(0.0), 
	lockInDates(0) {} // for reflection

void SPILockInGrowth::validatePop2Object(){
    static const string method = "SPILockInGrowth::validatePop2Object";
    
    if (Maths::isNegative(initialLockIn)) {
        isOK = false;
        err = "Initial Lock-In cannot be negative : " + Format::toString(initialLockIn);
        return;
    }
    if (Maths::isNegative(lockInPercentage)) {
        isOK = false;
        err = "Lock-In Percentage cannot be negative : " + Format::toString(lockInPercentage);
        return;
    }
    DateTime::ensureIncreasing(lockInDates, "Lock-In dates", false /*failIfEmpty*/);
}

DateTimeArray SPILockInGrowth::getEssentialDates() const {
    return lockInDates;
}

bool SPILockInGrowth::doesNothing() const {
    // all types other than standard are always deemed to be doing something
    return false;
}

IObject* SPILockInGrowth::defaultSPILockInGrowth(){
    return new SPILockInGrowth();
}

/** Invoked when Class is 'loaded' */
void SPILockInGrowth::load(CClassSP& clazz){
    REGISTER(SPILockInGrowth, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(ILockInSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPILockInGrowth);
    FIELD(initialLockIn,  "initialLockIn");
    FIELD(lockInPercentage,  "lockInPercentage");
    FIELD(lockInDates,       "lockInDates");
    FIELD(lockInFlag,        "lockInFlag");
    FIELD_MAKE_TRANSIENT(lockInFlag);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(initialLockIn);
    FIELD_MAKE_OPTIONAL(lockInPercentage);
    FIELD_MAKE_OPTIONAL(lockInDates);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPILockInGrowth::TYPE = CClass::registerClassLoadMethod(
    "SPILockInGrowth", typeid(SPILockInGrowth), SPILockInGrowth::load);

/*****************************************************************************/

double SPILockInCappedReserve::getInitialLockIn() const {
    return initialLockIn;
}

const ILockInSPI* SPILockInCappedReserve::getLockInSPINoInit() const {
    return this;
}

ILockInSPI* SPILockInCappedReserve::getLockInSPI(const DateTimeArray& rebalanceDates,
                                 const DateTime&      lastRebalDate) {
    // just call the init method
    init(rebalanceDates, lastRebalDate);

    return this;
};

void SPILockInCappedReserve::init(const DateTimeArray& sampleDates,
          const DateTime&      lastRebalDate) {
    static const string method = "SPILockInCappedReserve::init";
    if (lockInDates.size()>0 &&
        lockInDates.back()>lastRebalDate) {
        throw ModelException(method,
                             "Cannot have lock-in date (" + lockInDates.back().toString() + 
                             ") after final rebalance date (" + lastRebalDate.toString() +
                             ")");
    }

    lockInFlag = BoolArray(sampleDates.size(), false);
    int j=0;
    for(int i=0; i<sampleDates.size() && j<lockInDates.size(); i++) {
        if (sampleDates[i] == lockInDates[j]) {
            lockInFlag[i] = true;
            j++;
        }
    }
    // check all lock-in dates are sample dates
    if (j<lockInDates.size()) {
        throw ModelException(method,
                             "Lock-In dates must be a subset of rebalance dates - check lock-in date #"
                             + Format::toString(j+1) + " = " + lockInDates[j].toString());
        
    }
}

void SPILockInCappedReserve::apply(double&  BL,
           double   B,
           double   BF,
           int      iStep) const {
    if (lockInFlag[iStep]) {
        // note the floor of (B-BF)
        BL = BL + Maths::min(lockInCap, Maths::max(B - BF, 0.0));
    }
}

SPILockInCappedReserve::SPILockInCappedReserve(): SPIInterfaceIMS(TYPE), 
    initialLockIn(0.0), lockInCap(0.0), lockInDates(0) {} // for reflection

void SPILockInCappedReserve::validatePop2Object(){
    static const string method = "SPILockInCappedReserve::validatePop2Object";
    
    if (Maths::isNegative(initialLockIn)) {
        isOK = false;
        err = "Initial Lock-In cannot be negative : " + Format::toString(initialLockIn);
        return;
    }
    if (Maths::isNegative(lockInCap)) {
        isOK = false;
        err = "Lock-In Cap cannot be negative : " + Format::toString(lockInCap);
        return;
    }
    DateTime::ensureIncreasing(lockInDates, "Lock-In dates", false /*failIfEmpty*/);
}

DateTimeArray SPILockInCappedReserve::getEssentialDates() const {
    return lockInDates;
}

bool SPILockInCappedReserve::doesNothing() const {
    // all types other than standard are always deemed to be doing something
    return false;
}

IObject* SPILockInCappedReserve::defaultSPILockInCappedReserve(){
    return new SPILockInCappedReserve();
}

/** Invoked when Class is 'loaded' */
void SPILockInCappedReserve::load(CClassSP& clazz){
    REGISTER(SPILockInCappedReserve, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(ILockInSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPILockInCappedReserve);
    FIELD(initialLockIn,  "Initial Lock-In Level");
    FIELD(lockInCap,      "Cap");
    FIELD(lockInDates,       "Lock-In Dates");
    FIELD(lockInFlag,        "lockInFlag");
    FIELD_MAKE_TRANSIENT(lockInFlag);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(initialLockIn);
    FIELD_MAKE_OPTIONAL(lockInCap);
    FIELD_MAKE_OPTIONAL(lockInDates);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPILockInCappedReserve::TYPE = CClass::registerClassLoadMethod(
    "SPILockInCappedReserve", typeid(SPILockInCappedReserve), SPILockInCappedReserve::load);

#define SPI_LOCK_IN_TYPE_STD                "Standard"
#define SPI_LOCK_IN_TYPE_EXCESS             "Excess"
#define SPI_LOCK_IN_TYPE_SCHEDULE           "Schedule"
#define SPI_LOCK_IN_TYPE_GROWTH             "Growth"
#define SPI_LOCK_IN_TYPE_CAPPED_RESERVE     "CappedReserve"

// Painful - the optional strike field in SPI instrument needs to 
// be able to get lockin info early to support backwards compatability
// so this is needed publicly. That's the only reason - else would be private.
const ILockInSPI* SPILockInWrapper::getLockInSPINoInit() const {
    static const string routine = "SPILockInWrapper::getLockInSPINoInit";
    ILockInSPI* aLockIn = 0;

    if (SPILockInType.empty()){
        throw ModelException(routine, "Blank SPILockInType specified!");
    }
    if (SPILockInType==SPI_LOCK_IN_TYPE_STD) {
        if (!lockInStd.get()) {
            throw ModelException(routine, "Expected LockInStd but none supplied!");
        }
        aLockIn = lockInStd.get();
    } else if (SPILockInType==SPI_LOCK_IN_TYPE_EXCESS) {
        if (!lockInExcess.get()) {
            throw ModelException(routine, "Expected LockInExcess but none supplied!");
        }
        aLockIn = lockInExcess.get();
    } else if (SPILockInType==SPI_LOCK_IN_TYPE_SCHEDULE) {
        if (!lockInSchedule.get()) {
            throw ModelException(routine, "Expected LockInSchedule but none supplied!");
        }
        aLockIn = lockInSchedule.get();
    } else if (SPILockInType==SPI_LOCK_IN_TYPE_GROWTH) {
        if (!lockInGrowth.get()) {
            throw ModelException(routine, "Expected LockInGrowth but none supplied!");
        }
        aLockIn = lockInGrowth.get();
    } else if (SPILockInType==SPI_LOCK_IN_TYPE_CAPPED_RESERVE) {
        if (!lockInCappedReserve.get()) {
            throw ModelException(routine, "Expected LockInCappedReserve but none supplied!");
        }
        aLockIn = lockInCappedReserve.get();
    } else{
        throw ModelException(routine, "Unrecognised SPILockInType " + SPILockInType + 
                             ". Expected " + 
                             SPI_LOCK_IN_TYPE_STD + ", " +
                             SPI_LOCK_IN_TYPE_EXCESS + " or " + 
                             SPI_LOCK_IN_TYPE_SCHEDULE + " or " + 
                             SPI_LOCK_IN_TYPE_GROWTH + " or " + 
                             SPI_LOCK_IN_TYPE_CAPPED_RESERVE);
    }

    return aLockIn;
}

ILockInSPI* SPILockInWrapper::getLockInSPI(const DateTimeArray& rebalanceDates,
                                 const DateTime&      lastRebalDate) {
    static const string routine = "SPILockInWrapper::getLockInSPI";

    // cast away constness 
    theLockIn = const_cast<ILockInSPI*>(getLockInSPINoInit());

    // Nicer alternative is to create another class here, but for the moment we'll stay small (if a bit messy)
    // rebalanceDates are 
    theLockIn->init(rebalanceDates, 
                    lastRebalDate);

	// need to cast for validity check
	SPIInterfaceIMS* actualLockIn = dynamic_cast<SPIInterfaceIMS*>(theLockIn);
	if (!actualLockIn) {
        throw ModelException(routine, "Internal error in SPI lock-in validation");
	}
    if (!actualLockIn->isValid()) {
        throw ModelException(routine, "Invalid SPILockIn (type " + SPILockInType + 
                             ") : " + actualLockIn->errString());
    }
    return theLockIn;
};

// validation
void SPILockInWrapper::validatePop2Object(){
    static const string routine = "SPILockInWrapper::validatePop2Object";

}

/** Invoked when Class is 'loaded' */
void SPILockInWrapper::load(CClassSP& clazz){
    REGISTER(SPILockInWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ILockInSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPILockInWrapper);
    FIELD(SPILockInType, "SPILockInType");
    FIELD(lockInStd,  "lockInStd");
    FIELD_MAKE_OPTIONAL(lockInStd);
    FIELD(lockInExcess,  "lockInExcess");
    FIELD_MAKE_OPTIONAL(lockInExcess);
    FIELD(lockInSchedule,  "lockInSchedule");
    FIELD_MAKE_OPTIONAL(lockInSchedule);
    FIELD(lockInGrowth,  "lockInGrowth");
    FIELD_MAKE_OPTIONAL(lockInGrowth);
    FIELD(lockInCappedReserve,  "lockInCappedReserve");
    FIELD_MAKE_OPTIONAL(lockInCappedReserve);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

SPILockInWrapper::SPILockInWrapper(): CObject(TYPE), theLockIn(0){}

IObject* SPILockInWrapper::defaultSPILockInWrapper(){
    return new SPILockInWrapper();
}

CClassConstSP const SPILockInWrapper::TYPE = CClass::registerClassLoadMethod(
    "SPILockInWrapper", typeid(SPILockInWrapper), SPILockInWrapper::load);


DRLIB_END_NAMESPACE
