#include "edginc/config.hpp"
#include "edginc/CouponSched.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Format.hpp"
#include "esl.h" // just for RESET_ADVANCE, RESET_ARREARS definitions

DRLIB_BEGIN_NAMESPACE

/************************ CouponSchedInit ************************/

void CouponSchedInit::validatePop2Object(void) {
    if (startDate>endDate) throw ModelException("startDate>endDate");
    if (!frequency.get()) throw ModelException("frequency needed");
}


void CouponSchedInit::modify(CouponSchedDates &sched) {
    try {
        DateTimeArray dates;

        RatesUtils::genNewEvenDateList(startDate, endDate, 
            frequency->approxAnnualFrequency(), stubPos==StubPos::STUB_END, dates);

        DateTimeArray &accStart = getDates(sched.accStart);
        accStart.clear();
        accStart.insert(accStart.end(), dates.begin(), dates.end()-1); 
        DateTimeArray &accEnd = getDates(sched.accEnd);
        accEnd.clear();
        accEnd.insert(accEnd.end(), dates.begin()+1, dates.end()); 
        sched.pay = sched.accEnd;
        switch (resetPos) {
            case RESET_ADVANCE: getDates(sched.reset) = accStart; break;
            case RESET_ARREARS: getDates(sched.reset) = accEnd; break;
        }
        getDates(sched.resetEff) = getDates(sched.reset);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void CouponSchedInit::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CouponSchedInit, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(startDate, "");
    FIELD(endDate, "");
    FIELD(frequency, "");
    FIELD(stubPos, "");   FIELD_MAKE_OPTIONAL(stubPos);
    FIELD(resetPos, "");  FIELD_MAKE_OPTIONAL(resetPos);

    Addin::registerConstructor(Addin::UTILITIES, CouponSchedInit::TYPE);
}

CClassConstSP const CouponSchedInit::TYPE = CClass::registerClassLoadMethod(
    "CouponSchedInit", typeid(CouponSchedInit), load);

/************************ CouponSchedDates ************************/

string CouponSchedDates::calcDcfFrac(int cpnIdx, DayCountConvention& dcc) const {
    int days, denom;
    RatesUtils::calcDcfFrac(accStart[cpnIdx],accEnd[cpnIdx], dcc, days, denom);
    return Format::toString(days)+"/"+Format::toString(denom);
}

void CouponSchedDates::calcDcfs(DayCountConvention& dcc, DoubleArray &dcfs) const {
    int nbCoupons = accStart.size();
    dcfs.resize(nbCoupons);
    for (int i = 0; i < nbCoupons; i++) {
        dcfs[i] = RatesUtils::calcDcf(accStart[i],accEnd[i], dcc);
    }
}

static void check(int a, int b, const string& s, bool zeroOk) {
    if (a==b) return;
    if (a==0 && zeroOk) return;
    throw ModelException(s+".size="+Format::toString(a)+" should be "+Format::toString(b));
}

void CouponSchedDates::validatePop2Object() {
    if (initializer.get()) {
        initializer->modify(*this);
        // re-apply the modifiers to adjust the dates
        accStart.modifyDates();
        accEnd.modifyDates();
        reset.modifyDates();
        resetEff.modifyDates();
        pay.modifyDates();
    }

    nbCoupons = accStart.size();
    if (!nbCoupons) throw ModelException("accStart empty");
    check(accEnd.size(), nbCoupons, "accEnd", false);
    check(reset.size(), nbCoupons, "reset", true);
    check(resetEff.size(), nbCoupons, "resetEff", true);
    check(pay.size(), nbCoupons, "pay", true);

    DateTime::ensureIncreasing(pay.getDates(), "pay", true);

    int i;
	for (i=0; i<nbCoupons; ++i) {
        if (i && accEnd[i-1]!=accStart[i]) {
            DateTime accStartL = accStart[i];
            DateTime accEndL = accEnd[i-1];
            throw ModelException("Accrual start of coupon " + Format::toString(i-1) +
                " date: " + accStartL.toString()+ " and accrual end of next coupon " +
                accEndL.toString() + " do not match");
        }
		if (pay[i] < reset[i]) {
			throw ModelException("Coupon #"+Format::toString(i)+": "
				"Pay date is before reset date.");
		}
	}
    // add more checks on the consistency of the dates arrays
}

void CouponSchedDates::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CouponSchedDates, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(initializer, ""); FIELD_MAKE_OPTIONAL(initializer);
    FIELD(accStart, ""); FIELD_MAKE_OPTIONAL(accStart);
    FIELD(accEnd, "");   FIELD_MAKE_OPTIONAL(accEnd);
    FIELD(reset, "");    FIELD_MAKE_OPTIONAL(reset);
    FIELD(resetEff, ""); FIELD_MAKE_OPTIONAL(resetEff);
    FIELD(pay, "");      FIELD_MAKE_OPTIONAL(pay);

    FIELD(nbCoupons, ""); FIELD_MAKE_TRANSIENT(nbCoupons);

    Addin::registerConstructor(Addin::UTILITIES, CouponSchedDates::TYPE);
}

CClassConstSP const CouponSchedDates::TYPE = CClass::registerClassLoadMethod(
    "CouponSchedDates", typeid(CouponSchedDates), load);




/******************************** CouponSchedDisp ********************************/

class CouponSchedDisp : public CObject {
public:
    static CClassConstSP const TYPE;
    enum DatesEnum { ACCSTART, ACCEND, RESET, RESETEFF, PAY };

    IObjectSP addin() {
        DateTimeArray const *a;
        switch (datesToDisplay) {
            case ACCSTART: a = &val->accStart.getDates(); break;
            case ACCEND:   a = &val->accEnd.getDates(); break;
            case RESET:    a = &val->reset.getDates(); break;
            case RESETEFF: a = &val->resetEff.getDates(); break;
            case PAY:      a = &val->pay.getDates(); break;
            default: ASSERT(0);
        }
        return IObjectSP(new DateTimeArray(*a));
    }

private:
    CouponSchedDatesSP val;
    DatesEnum datesToDisplay;

    CouponSchedDisp() : CObject(TYPE)          {}
    static IObject* defaultConstructor(void) {return new CouponSchedDisp();}
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CouponSchedDisp, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);

        FIELD(val, "The CouponSched to display");    
        FIELD(datesToDisplay, "");    
        Addin::registerObjectMethod(
            "CouponSchedDisp", Addin::UTILITIES,
            "Displays the dates of a CouponSched object",
            false, Addin::expandSimple, &CouponSchedDisp::addin);
    }
};

CClassConstSP const CouponSchedDisp::TYPE = CClass::registerClassLoadMethod(
    "CouponSchedDisp", typeid(CouponSchedDisp), load);

START_PUBLIC_ENUM_DEFINITION(CouponSchedDisp::DatesEnum, "Stub position");
ENUM_VALUE_AND_NAME(CouponSchedDisp::ACCSTART, "ACCSTART", "");
ENUM_VALUE_AND_NAME(CouponSchedDisp::ACCEND, "ACCEND", "");
ENUM_VALUE_AND_NAME(CouponSchedDisp::RESET, "RESET", "");
ENUM_VALUE_AND_NAME(CouponSchedDisp::RESETEFF, "RESETEFF", "");
ENUM_VALUE_AND_NAME(CouponSchedDisp::PAY, "PAY", "");
END_ENUM_DEFINITION(CouponSchedDisp::DatesEnum);

START_PUBLIC_ENUM_DEFINITION(CouponSchedInit::ResetPos::Enum, "Reset position");
ENUM_VALUE_AND_NAME(CouponSchedInit::ResetPos::IN_ADVANCE, "IN_ADVANCE", "");
ENUM_VALUE_AND_NAME(CouponSchedInit::ResetPos::IN_ARREARS, "IN_ARREARS", "");
ENUM_VALUE_AND_NAME(CouponSchedInit::ResetPos::RESET_NONE, "RESET_NONE", "");
END_ENUM_DEFINITION(CouponSchedInit::ResetPos::Enum);


bool CouponSchedLoad() {return CouponSchedDates::TYPE !=0;}

DRLIB_END_NAMESPACE
