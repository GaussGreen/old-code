//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : OptionSched.hpp
//
//   Author      : Steve Marks
//
//
//----------------------------------------------------------------------------
#ifndef _OPTIONSCHED_HPP
#define _OPTIONSCHED_HPP

#include "edginc/MaturityPeriod.hpp"
#include "edginc/FlexDates.hpp"

DRLIB_BEGIN_NAMESPACE

class OptionSchedDates;

/******************************** OptionSchedInit *******************************/
class OptionSchedInit : public CObject,
                        public ExternalFlexDatesModifier
{

public:
    static CClassConstSP const TYPE;

    /************* exported fields **************/
    DateTime startDate;
    DateTime endDate;
    MaturityPeriodSP frequency;
    StubPos::Enum stubPos;

    /************* methods *************/
public:
    virtual void validatePop2Object();
    virtual void modify(OptionSchedDates &sched);

protected:
    OptionSchedInit();

private:
    static IObject* defaultConstructor(void) {return new OptionSchedInit();}
    static void load(CClassSP& clazz);
};

typedef smartPtr<OptionSchedInit> OptionSchedInitSP;

/******************************** OptionSchedDates *******************************/
/* The class representing the schedule */

class IRPRODUCTS_DLL OptionSchedDates : public CObject {
public:
    static CClassConstSP const TYPE;

    /************* exported fields **************/
    OptionSchedInitSP initializer;

    FlexDates notifDate;
    FlexDates exerciseDate;

    /************* transient fields **************/
    int nbDates;

    /************** methods **************/

    virtual void validatePop2Object();
    virtual void setup(const IModel* model, const MarketData* market) {validatePop2Object();}

    OptionSchedDates() :  CObject(TYPE) {}
protected:
    OptionSchedDates(CClassConstSP const &type) :  CObject(type) {}
private:
    static IObject* defaultConstructor(void) {return new OptionSchedDates();}
    static void load(CClassSP& clazz);
};

typedef smartPtr<OptionSchedDates> OptionSchedDatesSP;


DRLIB_END_NAMESPACE

#endif


