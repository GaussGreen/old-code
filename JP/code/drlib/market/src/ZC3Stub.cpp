//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZC3Stub.cpp
//
//   Description : Helper class for ALIB zero curve 3 bootstrapping method.
//
//   Author      : Richard Appleton
//
//   Date        : 2nd May 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZC3Stub.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/Format.hpp"


DRLIB_BEGIN_NAMESPACE


int ZC3Stub::STUB_BASIS = CompoundBasis::SIMPLE;
DayCountConvention* ZC3Stub::STUB_DCC = DayCountConventionFactory::make("Act/360");


// ALIB: szcprvt.c#639 (GtoSZCPrepareStub)
ZC3Stub::ZC3Stub(
     const DateTime&           pStartDate,
     const DateTime&           pEndDate,
     double                    rate,
     const DayCountConvention& dcc,
     double                    basis,
     const DayCountConvention& basisDcc,
     const BadDayConvention&   badDayConv,
     const Holiday&            holidays)
{
    static const string method = "Stub::Stub";

    /*
     * First of all adjust start and end date for bad days. We ignore the fact 
     * that the basis might have a different holiday file.  We cannot create 
     * contiguous stub instruments in the case that the cash flow dates for 
     * the stub are different.
     */
    startDate = badDayConv.adjust(pStartDate, &holidays);
    endDate = badDayConv.adjust(pEndDate, &holidays);

    if (startDate >= endDate)
    {
        string msg = Format::toString(
            "Invalid rate from %s to %s after adjustment",
            startDate.toString().c_str(),
            endDate.toString().c_str());
        throw ModelException(method, msg);
    }

    /*
     * Convert the input rates into the correct stub rate
     */
    stubRate = RateConversion::rateConvert(startDate, endDate, 
        rate, &dcc, CompoundBasis::SIMPLE, STUB_DCC, STUB_BASIS);

    double stubBasis = RateConversion::rateConvert(startDate, endDate, 
        basis, &basisDcc, CompoundBasis::SIMPLE, STUB_DCC, STUB_BASIS);

    stubRate += stubBasis;
}


ZC3Stub::ZC3Stub(const ZC3Stub& other)
  : startDate(other.startDate), endDate(other.endDate), stubRate(other.stubRate)
{
}


ZC3Stub::~ZC3Stub()
{
}


ZC3Stub& ZC3Stub::operator=(const ZC3Stub& other)
{
    if (this != &other)
    {
        startDate = other.startDate;
        endDate = other.endDate;
        stubRate = other.stubRate;
    }

    return *this;
}


// ALIB:: szcprvt.c#762
bool ZC3Stub::operator<(const ZC3Stub& other) const
{ 
    return startDate < other.startDate;
}


string ZC3Stub::toString() const
{
    string result = "[";
    result += startDate.toString();
    result += ",";
    result += endDate.toString();
    result += "]";
    return result;
}


const DateTime& ZC3Stub::getStartDate() const
{
    return startDate;
}


const DateTime& ZC3Stub::getEndDate() const
{
    return endDate;
}


double ZC3Stub::getStubRate() const 
{ 
    return stubRate; 
}


bool ZC3Stub::isContiguousWith(const ZC3Stub& next) const
{
    // end date for this stub must be start date for next stub.
    return endDate == next.startDate;
}


DRLIB_END_NAMESPACE
