//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : CreditEventOverride.hpp
//
//   Description : Simple class used by credit instruments (e.g., CIS) in order 
//                 to override (at instrument level) names' default-related 
//                 parameters.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CREDITEVENTOVERRIDE_HPP
#define QLIB_CREDITEVENTOVERRIDE_HPP

#include "edginc/Object.hpp"
#include "edginc/ICreditEventOverride.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CreditEventOverrideName);

/** Simple class used by credit instruments (e.g., CIS) in order to
 * override (at instrument level) names' default-related parameters */
class PRODUCTS_DLL CreditEventOverride: public CObject,
                                        public virtual ICreditEventOverride 
{
public:
    static CClassConstSP const TYPE;
    virtual ~CreditEventOverride();

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Returns the override information for the supplied name. If no override
     * information is available returns ICreditEventOverrideNameSP() */
    virtual ICreditEventOverrideNameSP getOverrideForName(
        const string& name) const;

    /** GetMarket implementation*/
    virtual void getMarket(const IModel* model, const MarketData* market);
    
private:
    // For reflection
    static void load (CClassSP& clazz);
    CreditEventOverride();
    static IObject* defaultCreditEventOverride();
    
    // Fields
    ICreditEventOverrideNameArraySP overrideArray; /* Array of override params
                                                    * for the individual names */
};

DECLARE(CreditEventOverride);


DRLIB_END_NAMESPACE

#endif
