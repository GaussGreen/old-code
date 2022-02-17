#ifndef EDR_MC_PRODUCT_ENGINE_HOOKS_HPP
#define EDR_MC_PRODUCT_ENGINE_HOOKS_HPP

#include "edginc/config.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/Model.hpp"
#include "edginc/IMCIntoProduct.hpp"

DRLIB_BEGIN_NAMESPACE


// Interface required of product to support LogNormal path generator
class MCARLO_DLL IMCProductLN {
public:

    virtual ~IMCProductLN() {};

    /** Path Generator will call this function for each asset. The
        returned array is [NbPath] - i.e. number of interp levels for
        this asset. Clients may call refLevel() on the supplied path
        generator (but not other funtions in general) */
    virtual CVolRequestLNArray getVolInterp(
        const MCPathGenerator* pathGenerator,
        int                     iAsset) const = 0;

    /** Use this opportunity to do any model driven initialisation
        of the instrument. e.g closed from barrier adjustments */
    virtual void initialiseLN(const MCPathGenerator* pathGenerator)const = 0;

};

// Interface required of product to support Implied path generator
class MCARLO_DLL IMCProductImplied{
public:
     /** Use this opportunity to do any model driven initialisation
        of the instrument. e.g closed from barrier adjustments */
    virtual void initialiseImplied(const MCPathGenerator* pathGenerator)const = 0;

    virtual ~IMCProductImplied(){}
};

DRLIB_END_NAMESPACE

#endif // EDR_MC_PRODUCT_ENGINE_HOOKS_HPP
