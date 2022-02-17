#include "edginc/config.hpp"
#include "edginc/RiskyLogOfDiscFactorKey.hpp"

DRLIB_BEGIN_NAMESPACE

RiskyLogOfDiscFactorKey::RiskyLogOfDiscFactorKey(const IDiscountCurveRisky* idcr)
{
    discKey  = idcr->getDiscountKey();
    riskyKey = idcr->getRiskyKey();
}

/** Returns the log of the risky discount factor between the two dates */
double RiskyLogOfDiscFactorKey::calc(const DateTime&  loDate,
                                     const DateTime&  hiDate)
{
    //sum the two components logs
    return discKey->calc(loDate,hiDate) + riskyKey->calc(loDate,hiDate);
}

DRLIB_END_NAMESPACE

