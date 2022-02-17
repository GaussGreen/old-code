//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : VolProcessedCalib.cpp
//
//   Description : Simple vol processed object
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Nov-2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VolProcessedCalib_CPP
#include "edginc/VolProcessedCalib.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

void VolProcessedCalib::load(CClassSP& clazz)
{
    REGISTER(VolProcessedCalib, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IVolProcessed);
    FIELD(name, "name");
    FIELD(metric, "time metric");
}

VolProcessedCalib::VolProcessedCalib(
    DateTimeArraySP swapStartDates,
    DateTimeArraySP swapMaturityDates,
    DoubleArraySP   vols,
    ExpiryArraySP   selectedTenors,
    IntArraySP      selectedTenorIndices,
    ExpiryArraySP   selectedExpiries,
    IntArraySP      selectedExpiryIndices,
    DateTimeArraySP selectedSwapMaturityDates)
    :
    CObject(TYPE),
    m_swapStartDates(swapStartDates),
    m_swapMaturityDates(swapMaturityDates),
    m_vols(vols),
    m_selectedTenors(selectedTenors),
    m_selectedTenorIndices(selectedTenorIndices),
    m_selectedExpiries(selectedExpiries),
    m_selectedExpiryIndices(selectedExpiryIndices),
    m_selectedSwapMaturityDates(selectedSwapMaturityDates)
{
}    

VolProcessedCalib::VolProcessedCalib(const CClassConstSP& clazz) : CObject(clazz) {}

VolProcessedCalib::~VolProcessedCalib() {}

CClassConstSP const VolProcessedCalib::TYPE = CClass::registerClassLoadMethod(
    "VolProcessedCalib", typeid(VolProcessedCalib), load);

//-----------------------
// IVolProcessed methods
//-----------------------
string VolProcessedCalib::getName() const
{
    return name;
}

double VolProcessedCalib::calcTradingTime(const DateTime &date1, const DateTime &date2) const
{
    return metric->yearFrac(date1, date2);
}

TimeMetricConstSP VolProcessedCalib::GetTimeMetric() const
{
    return TimeMetricConstSP(metric);
}

DRLIB_END_NAMESPACE
