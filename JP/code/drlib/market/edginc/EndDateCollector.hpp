//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EndDateCollector.hpp
//
//   Description : End date collector class
//
//   Author      : Andrew McCleery
//
//   Date        : 28 May 2004
//
//
//   $Log: $
//
//----------------------------------------------------------------------------

#ifndef ENDDATECOLLECTOR_HPP
#define ENDDATECOLLECTOR_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE
class Sensitivity;

/** A class to get the maximum endDate of components */
class MARKET_DLL EndDateCollector {
public:
    EndDateCollector(const DateTime&    initialEndDate, 
                     const Sensitivity* sensitivity);
    DateTime getMaxEndDate(IObjectConstSP obj) const;

private:
    // Not implemented
    EndDateCollector(const EndDateCollector &rhs);
    EndDateCollector& operator=(const EndDateCollector& rhs);

    DateTime            endDate;
    const Sensitivity*  sensitivity;
};

typedef smartPtr<EndDateCollector> EndDateCollectorSP;
typedef smartPtr<const EndDateCollector> EndDateCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
