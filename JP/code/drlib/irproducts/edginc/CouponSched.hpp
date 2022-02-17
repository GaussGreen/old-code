//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KOption.hpp
//
//   Author      : Marc Alzina
//
//   Description : option komponent
//
//----------------------------------------------------------------------------
#ifndef _INSTRSCHED_HPP
#define _INSTRSCHED_HPP

#include "edginc/DayCountConvention.hpp"
#include "edginc/FlexDates.hpp"

DRLIB_BEGIN_NAMESPACE

class CouponSchedDates;

/******************************** CouponSchedInit *******************************/
// This class is used to initialize the dates of the FlexDates in 
// CouponSchedDates (accStart,accEnd,reset,resetEff,pay).
// It overwrites any date provided in any FlexDates::datesExpl. 
// The modifiers of the FlexDates are then applied on top 
// of it (see CouponSchedDates::validatePop2Object) so that
// we can have different adjustments for each of them
// It does not really make sense to have (for instance) startDate.modifier
// of type DatesSched because it would overwrite what CouponSchedInit
// has done.
class CouponSchedInit : public CObject,
                        public ExternalFlexDatesModifier
{

public:
    static CClassConstSP const TYPE;

    struct ResetPos {
        enum Enum { IN_ADVANCE, IN_ARREARS, RESET_NONE };
    };
    typedef BoxedEnum<ResetPos::Enum> ResetPosBoxedEnum;

    /************* exported fields **************/
protected:
    DateTime startDate, endDate;
    MaturityPeriodSP frequency;
    StubPos::Enum stubPos;
    ResetPos::Enum resetPos;

    /************* methods *************/
public:
    virtual void validatePop2Object();
    virtual void modify(CouponSchedDates &sched);

protected:
    CouponSchedInit() : CObject(TYPE), stubPos(StubPos::STUB_NONE), resetPos(ResetPos::RESET_NONE) {}

private:
    static IObject* defaultConstructor(void) {return new CouponSchedInit();}
    static void load(CClassSP& clazz);
};

typedef smartPtr<CouponSchedInit> CouponSchedInitSP;

/******************************** CouponSchedDates *******************************/
/* The class representing the schedule */

class CouponSchedDates : public CObject {
public:
    static CClassConstSP const TYPE;

    /************* exported fields **************/
    CouponSchedInitSP initializer;

    FlexDates accStart;
    FlexDates accEnd;
    FlexDates reset;
    FlexDates resetEff;
    FlexDates pay;


    /************* transient fields **************/
    int nbCoupons;

    /************** methods **************/

    void calcDcfs(DayCountConvention& dcc, DoubleArray &dcfs) const;
    string calcDcfFrac(int cpnIdx, DayCountConvention& dcc) const;

    virtual void validatePop2Object();
    virtual void setup(const IModel* model, const MarketData* market) {validatePop2Object();}

    CouponSchedDates() :  CObject(TYPE) {}
protected:
    CouponSchedDates(CClassConstSP const &type) :  CObject(type) {}
private:
    static IObject* defaultConstructor(void) {return new CouponSchedDates();}
    static void load(CClassSP& clazz);
};

typedef smartPtr<CouponSchedDates> CouponSchedDatesSP;

DRLIB_END_NAMESPACE

#endif


