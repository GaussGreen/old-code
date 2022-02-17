//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZC3Interval.cpp
//
//   Description : Helper class for ALIB zero curve 3 bootstrapping method.
//
//   Author      : Richard Appleton
//
//   Date        : 9th May 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZC3Interval.hpp"
#include "edginc/Maths.hpp"
#include <algorithm>



DRLIB_BEGIN_NAMESPACE

/* static */
// ALIB: szcrates.c#3310
void ZC3Interval::mergeAndSort(ZC3IntervalArray& intervals)
{
    static const string method = "ZC3IntervalArray::mergeAndSort";
    
    /*
     * Perform interval merging. First we sort the intervals by the end points.
     * However we may have two intervals with the same end points, in which 
     * case we will substitute its larger one by its complement to the 
     * smaller, adjusting the rate respectively.
     *
     * This may well produce a situation where the rate intervals are no longer 
     * sorted.
     *
     * Therefore we iterate this situation until there are no more 
     * replacements. In other words all end dates will be in strictly
     * increasing order.
     */
    int numSubs;

    do
    {
        /*
         * We now sort the rate intervals. This array has not been sorted before
         * this stage in the algorithm for the first pass into this loop.
         */
        sort(intervals.begin(), intervals.end(), ZC3Interval::lessThan);

        numSubs = 0; // we continue until there are no more substitutions

        for (int i = 1 ; i < intervals.size() ; i++)
        {
            ZC3Interval& previous = intervals[i-1];
            ZC3Interval& current = intervals[i];

            if (previous.getEndDate() > current.getEndDate())
            {
                string msg = Format::toString(
                    "Program bug - intervals incorrectly sorted [%s line %d].", 
                    __FILE__, __LINE__);
                throw ModelException(method, msg);
            }

            if (previous.getEndDate() == current.getEndDate())
            {
                if (previous.getStartDate() > current.getStartDate())
                {
                    string msg = Format::toString(
                        "Program bug - intervals incorrectly sorted [%s line %d].",
                        __FILE__, __LINE__);
                    throw ModelException(method, msg);
                }

                if (previous.getStartDate() == current.getStartDate())
                {
                    string msg = Format::toString(
                        "Two rates determined for interval between %s and %s",
                        previous.getStartDate().toString().c_str(),
                        previous.getEndDate().toString().c_str());
                    throw ModelException(method, msg);
                }
        
                previous.setEndDate(current.getStartDate());
                previous.setContRate( previous.getContRate() - current.getContRate() );
                ++numSubs;
            }
        }

    } while(numSubs > 0);

    /*
     * If we were paranoid we could test that all the end dates are in strictly
     * increasing order at this point.
     */
}


/* static */
bool ZC3Interval::lessThan(const ZC3Interval& first, const ZC3Interval& second)
{
    if (first.endDate == second.endDate)
    {
        return first.startDate < second.startDate;
    }
    else
    {
        return first.endDate < second.endDate;
    }
}


// for STL use only
ZC3Interval::ZC3Interval() : CObject(TYPE)
{
}


ZC3Interval::ZC3Interval(const CClassConstSP& clazz) : CObject(clazz)
{
}


ZC3Interval::ZC3Interval(const DateTime& start, const DateTime& end, double rate, const CClassConstSP& clazz)
  : CObject(clazz), startDate(start), endDate(end), contRate(rate)
{
}


string ZC3Interval::toString() const
{
    string result = "[";
    result += startDate.toString();
    result += ",";
    result += endDate.toString();
    result += "]";
    return result;
}


const DateTime& ZC3Interval::getStartDate() const
{
    return startDate;
}


const DateTime& ZC3Interval::getEndDate() const
{
    return endDate;
}


double ZC3Interval::getContRate() const
{
    return contRate;
}


void ZC3Interval::setEndDate(const DateTime& date)
{
    endDate = date;
}


void ZC3Interval::setContRate(double rate)
{
    contRate = rate;
}


bool ZC3Interval::isForSamePeriodAs(const ZC3Interval& other) const
{
    return startDate == other.startDate && endDate == other.endDate;
}


bool ZC3Interval::startsWithin(const DateTime& start, const DateTime& end) const
{
    return startDate >= start && startDate < end;
}



// for STL use only
ZC3Adjustment::ZC3Adjustment() : ZC3Interval(ZC3Adjustment::TYPE)
{
}


ZC3Adjustment::ZC3Adjustment(const DateTime& start, const DateTime& end, double rate, double pDiscount)
  : ZC3Interval(start, end, rate, TYPE), discount(pDiscount)
{
}


bool ZC3Adjustment::isAbsolute() const
{
    return !Maths::isZero(discount);
}


double ZC3Adjustment::getDiscount() const
{
    return discount;
}


void ZC3Adjustment::modify(double pContRate)
{
    contRate = pContRate;
}


/*
 * Reflection support
 */

static IObject* defaultInterval()
{
    return new ZC3Interval();
}


/** Invoked when Class is 'loaded' */
void ZC3Interval::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3Interval, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultInterval);

    //fields
    FIELD_NO_DESC(startDate);
    FIELD_NO_DESC(endDate);
    FIELD_NO_DESC(contRate);
}


CClassConstSP const ZC3Interval::TYPE = 
    CClass::registerClassLoadMethod("ZC3Interval", typeid(ZC3Interval), load);


static IObject* defaultAdjustment()
{
    return new ZC3Adjustment();
}


/** Invoked when Class is 'loaded' */
void ZC3Adjustment::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3Adjustment, clazz);
    SUPERCLASS(ZC3Interval);
    EMPTY_SHELL_METHOD(defaultAdjustment);

    //fields
    FIELD_NO_DESC(discount);
}


CClassConstSP const ZC3Adjustment::TYPE = 
    CClass::registerClassLoadMethod("ZC3Adjustment", typeid(ZC3Adjustment), load);


/**
 * Array support.
 */

DEFINE_TEMPLATE_TYPE(ZC3IntervalArray);


/* static */
IObjectConstSP arrayObjectCast<ZC3Interval>::toIObject(const ZC3Interval& value)
{
    IObjectConstSP ptr(IObjectConstSP::attachToRef(&value));
    return ptr;
}


/* static */
IObjectSP arrayObjectCast<ZC3Interval>::toIObject(ZC3Interval& value)
{
    return IObjectSP::attachToRef(&value);
}


/* static */
ZC3Interval arrayObjectCast<ZC3Interval>::fromIObject(IObjectSP& value)
{
    ZC3Interval* ptr = DYNAMIC_CAST(ZC3Interval, value.get());
    if (!ptr)
    {
        throw ModelException("arrayObjectCast::fromIObject", "Object is not a ZC3Interval");
    }

    return *ptr;
}


DEFINE_TEMPLATE_TYPE(ZC3AdjustmentArray);


/* static */
IObjectConstSP arrayObjectCast<ZC3Adjustment>::toIObject(const ZC3Adjustment& value)
{
    IObjectConstSP ptr(IObjectConstSP::attachToRef(&value));
    return ptr;
}


/* static */
IObjectSP arrayObjectCast<ZC3Adjustment>::toIObject(ZC3Adjustment& value)
{
    return IObjectSP::attachToRef(&value);
}


/* static */
ZC3Adjustment arrayObjectCast<ZC3Adjustment>::fromIObject(IObjectSP& value)
{
    ZC3Adjustment* ptr = DYNAMIC_CAST(ZC3Adjustment, value.get());
    if (!ptr)
    {
        throw ModelException("arrayObjectCast::fromIObject", "Object is not a ZC3Adjustment");
    }

    return *ptr;
}


DRLIB_END_NAMESPACE
