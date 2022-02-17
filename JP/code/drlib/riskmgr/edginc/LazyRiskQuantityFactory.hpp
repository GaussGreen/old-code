/**
 * @file LazyRiskQuantityFactory.hpp
 */

#ifndef DRLIB_LazyRiskQuantityFactory_H
#define DRLIB_LazyRiskQuantityFactory_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IHypothesis)
FORWARD_DECLARE(IRiskQuantityFactory)
FORWARD_DECLARE(LazyRiskQuantityFactory)

/**
 * An IRiskQuantityFactory to be invoked in an alternate state of the world.
 *
 * This class is not used at the moment.  See
 * IRiskQuantityFactory::lazies().
 */

class RISKMGR_DLL LazyRiskQuantityFactory: public CObject {

    LazyRiskQuantityFactory(const LazyRiskQuantityFactory& rhs);
    LazyRiskQuantityFactory& operator=(const LazyRiskQuantityFactory& rhs);
    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    LazyRiskQuantityFactory(IHypothesisConstSP hypothesis,
                            IRiskQuantityFactorySP greek);

    ~LazyRiskQuantityFactory();

    IHypothesisConstSP hypothesis;
    IRiskQuantityFactorySP greek;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
