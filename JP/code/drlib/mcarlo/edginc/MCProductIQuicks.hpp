#ifndef EDR_MC_PRODUCT_IQUICKS_HPP
#define EDR_MC_PRODUCT_IQUICKS_HPP

#include "edginc/config.hpp"
#include "edginc/ScalarShift.hpp"

DRLIB_BEGIN_NAMESPACE

class IMCPrices;
class MCPathGenerator;
class Sensitivity;

/** Interface that derived instances of IMCProduct implement in
    order to add support for the 'quick greeks' ie the skipping of
    paths whose value will not change under a tweak.  */
class MCARLO_DLL IMCQuickGreeks {
public:
    virtual ~IMCQuickGreeks() {};

    /** Set IMCPrices object for doing greek sens */
    virtual void setPricesForGreek(
        IMCPrices*               untweakedPrices,
        const MCPathGenerator* futurePathGen,
        const Sensitivity*    sens) = 0;

    /** Set IMCPrices object for doing two sided greek sens */
    virtual void setPricesForTwoSidedGreek(
        IMCPrices*                 untweakedPrices,
        const MCPathGenerator*   futurePathGen,
        const ScalarShiftArray& sens) = 0;
};

/** Interface that derived instances of IMCProduct implement in
    order to add support for the 'quick x gamma' ie the skipping of
    paths whose Cross Gamma is zero */
class MCARLO_DLL IMCQuickXGamma{
public:
    virtual ~IMCQuickXGamma() {};

    /** Create IMCPrices object for doing quick x gamma. The IMCPrices
        object needs to be able to skip over paths that have zero
        x gamma. The sens array contains the shifts defining what range
        the cross gamma will be calculated over */
    virtual void setPricesForXGamma(
        IMCPrices*                  untweakedPrices,
        const MCPathGenerator*    futurePathGen,
        const ScalarShiftArray&  sens) = 0;
};

DRLIB_END_NAMESPACE

#endif // EDR_MC_PRODUCT_IQUICKS_HPP
