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
#ifndef _FLEXDATES_HPP
#define _FLEXDATES_HPP

#include "edginc/MaturityPeriod.hpp"
#include "edginc/SharedEnums.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

/******************************* DatesModifier ************************************/
class IDatesModifier : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    virtual void modify(DateTimeArray&)=0;

private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<IDatesModifier> IDatesModifierSP;

/******************************* FlexDates ************************************/

class FlexDates : public CObject {
    friend class IDatesModifier;
    friend class ExternalFlexDatesModifier;
public:
    static CClassConstSP const TYPE;

    /************** exported fields ************/
private:
    DateTimeArraySP datesExpl; // explicit dates (optional)
    IDatesModifierSP modifier; // modifier (optional)

    /************** transient fields ************/
    // This private field can be modified by
    // friends IDatesModifier and ExternalFlexDatesModifier
    mutable DateTimeArray dates; // the computed dates

    /************** methods **************/
public:
    void modifyDates() {
        if (modifier.get()) 
            modifier->modify(dates);
    }

    virtual void validatePop2Object() {
        if (datesExpl.get()) 
            dates = *datesExpl;
        modifyDates();
    }

    virtual void setup(const IModel* model, const MarketData* market) {
        validatePop2Object();
    }

    /* make the object look like a DateTimeArray for convenience */
    DateTime           operator[](int i) const  {return dates[i];}
    DateTimeArray const& getDates() const       {return dates;}
    operator           DateTimeArray&()         {return dates;}
    operator const     DateTimeArray&() const   {return dates;}
    int                size() const             {return dates.size();}
    DateTime           front() const            {return dates.front();}
    DateTime           back()  const            {return dates.back();}
    bool               empty() const            {return dates.empty();}
    DateTimeArray::const_iterator begin() const {return dates.begin();}
    DateTimeArray::const_iterator end()   const {return dates.end();}

    // Returns datesIn[] where the corresponding weights[] is not zero
    static void filterOnWeight(
        DateTimeArray &datesIn,    /* IN  */
        DoubleArray &weights,      /* IN  */
        DateTimeArray &datesOut);  /* OUT */

    FlexDates() : CObject(TYPE) {}
private:
    IObjectSP addinDisplay() {
        return IObjectSP(new DateTimeArray(dates));
    }
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new FlexDates(); }
};

typedef smartPtr<FlexDates> FlexDatesSP;

/************************** ExternalDatesModifier ******************************/
// Any class which is not a IDatesModifier that still want to mess with 
// FlexDates internal "dates" array has to derive from ExternalDatesModifier 
// so that it is clear who does what.
// This is the case for CouponSched and OptionSched for instance.
class ExternalFlexDatesModifier {
    protected:
    inline DateTimeArray& getDates(FlexDates &flexDates) { return flexDates.dates; }
};

/*********************************** DatesSched ********************************/

// This modifier actually overrides the entire date list with dates
// generated from start/end/frequency
class DatesSched : public CObject,
                   virtual public IDatesModifier
{
public:
    static CClassConstSP const TYPE;
    /************* exported fields *************/
protected:
    DateTime startDate, matDate;
    MaturityPeriodSP frequency;
    StubPos::Enum stubPos;

    /************* methods *************/
public:
    virtual void modify(DateTimeArray &dates);

    DatesSched( 
        DateTime startDate, DateTime matDate, 
        MaturityPeriodSP frequency,
        StubPos::Enum stubPos)
    : CObject(TYPE), 
        startDate(startDate), matDate(matDate), 
        frequency(frequency), stubPos(stubPos)
    { validatePop2Object(); }

protected:
    DatesSched() :  CObject(TYPE) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new DatesSched(); }
};

/******************************* DatesAdj ************************************/
// miscellaneous adjustments to the date list
class DatesAdj : public CObject,
                 virtual public IDatesModifier
{
public:
    static CClassConstSP const TYPE;

    /************* exported fields *************/
public:
    IDatesModifierSP modifier; // reference to another modifiers (to build the schedule for instance)

    struct RangeBound {
        enum Enum { START, END, NONE };
    };
    typedef BoxedEnum<RangeBound::Enum> RangeBoundBoxedEnum;

    struct ShiftMode {
        enum Enum { NONE, FWD, BWD };
    };
    typedef BoxedEnum<ShiftMode::Enum> ShiftModeBoxedEnum;

    RangeBound::Enum rangeBound;
    ExpirySP shiftLength;
    ShiftMode::Enum shiftMode;

    /************* methods *************/
public:
    virtual void modify(DateTimeArray &dates);

protected:
    DatesAdj() :  CObject(TYPE), rangeBound(RangeBound::NONE), shiftMode(ShiftMode::NONE) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) {return new DatesAdj();}
};

DRLIB_END_NAMESPACE

#endif


