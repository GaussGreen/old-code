/**
 * @file IEqVolNamePair.hpp
 */

#ifndef QLIB_IEqVolNamePair_H
#define QLIB_IEqVolNamePair_H

#include <map>
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/OutputName.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IEqVolNamePair)

/**
 * Interface for objects which know how to associate an asset with its
 * volatility
 *
 * This was part of a mechanism constructed for use by DDeltaDVol; I've
 * reused it for Tier.
 */

class RISKMGR_DLL IEqVolNamePair: public virtual IObject {
public:

    virtual ~IEqVolNamePair();
    static CClassConstSP const TYPE;

    /**
     * Place the name of the asset, and the name of its vol object, in the
     * variables provided
     *
     * Returns true if we want you to recurse into our child objects to look
     * for further IEqVolNamePair's.
     */

    virtual bool getNamePairs(string& eqName, string& volName) const = 0;

    /**
     * A table of all the assets in the world and their associated vols
     */

    static map<string, OutputNameSP> namePairs(IObjectConstSP root);
};

DRLIB_END_NAMESPACE

#endif
