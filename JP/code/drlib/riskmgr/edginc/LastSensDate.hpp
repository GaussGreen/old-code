//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LastSensDate.hpp
//
//   Description : Interface for instruments that know when to stop
//                 tweaking pointwise greeks
//
//   Author      : Andrew J Swain
//
//   Date        : 7 June 2001
//
//
//----------------------------------------------------------------------------

#ifndef LASTSENSDATE_HPP
#define LASTSENSDATE_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE
class Sensitivity;
class CInstrument;

/** Interface for instruments that know when to stop
    tweaking pointwise greeks
*/
class RISKMGR_DLL LastSensDate {
public:
    static CClassConstSP const TYPE;

    virtual ~LastSensDate();

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const = 0;

protected:
    LastSensDate();

    /** Default implementation of endDate (above) for product-provided endDate
        (see below). */
    DateTime productEndDate(const CInstrument *instrument,
                            const Sensitivity *sensControl) const;
};

/** Interface for product-provided endDate */
class RISKMGR_DLL LastProductSensDate {
public:
    static CClassConstSP const TYPE;

    virtual ~LastProductSensDate();

    /** When to stop tweaking */
    virtual DateTime endDate(const CInstrument* instrument,
                             const Sensitivity* sensControl) const = 0;
protected:
    LastProductSensDate();
};


DRLIB_END_NAMESPACE
#endif
