#ifndef QR_RISKYLOGOFDISCFACTORKEY_HPP
#define QR_RISKYLOGOFDISCFACTORKEY_HPP

#include "edginc/DefaultRates.hpp"
#include "edginc/IDiscountCurveRisky.hpp"

DRLIB_BEGIN_NAMESPACE

//// Class used for optimising repeated 
//// calculations of risky discount factors
class MARKET_DLL RiskyLogOfDiscFactorKey: public IDiscountCurve::IKey
{
public:
    RiskyLogOfDiscFactorKey(const IDiscountCurveRisky* idcr);

    /** Returns the log of the risky discount factor between the two dates */
    virtual double calc(const DateTime&  loDate,
                        const DateTime&  hiDate);

private:
    IDiscountCurve::IKey* discKey;
    DefaultRates::IKey*   riskyKey;
};

DRLIB_END_NAMESPACE
#endif
