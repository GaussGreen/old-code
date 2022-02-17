//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : Interface to an object that provides protection and can  
//                 determine if a date is covered for protection or not
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IPROTECTIONPROVIDER_HPP
#define QLIB_IPROTECTIONPROVIDER_HPP

#include "edginc/DateTime.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/** Interface for objects capable of performing bad day adjustmens */
class MARKET_DLL IProtectionProvider : public virtual IObject {
public:
    static CClassConstSP const TYPE;
    virtual ~IProtectionProvider();

    /* Checks if the input date is covered for protection */
    virtual bool isDateCoveredForProtection(const DateTime& date) const = 0;

protected:
    IProtectionProvider();

private:
    static void load(CClassSP& clazz);
};

DECLARE(IProtectionProvider);

DRLIB_END_NAMESPACE

#endif
