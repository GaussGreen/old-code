/**
 * @file IQualifiedRiskAxis.hpp
 */

#ifndef QLIB_IQualifiedRiskAxis_H
#define QLIB_IQualifiedRiskAxis_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/IRiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

class Void;
FORWARD_DECLARE(Expiry)
FORWARD_DECLARE(ExpiryWindow)
FORWARD_DECLARE(ExpiryPair)
FORWARD_DECLARE(ExpiryAndStrike)

/**
 * A property of the world, with respect to which risk can be estimated, and
 * which may be associated with a "qualifier".
 *
 * IRiskAxis implementations arising from IRiskProperty's naturally have a
 * "qualifier": an ExpiryWindow for term-structured properties, or Void for
 * "scalar" ones (see RiskProperty::axisFor() for an explanation).  This
 * templated interface extends IRiskAxis with a qualifier() method for
 * retrieving it.  It's only used in RiskMappingMatrix.  The chief
 * implementation is PropertyRiskAxis<PROPERTY>.
 */

template <class QUALIFIER>
class IQualifiedRiskAxis: public virtual IRiskAxis {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    virtual smartConstPtr<QUALIFIER> qualifier() const = 0;

    IQualifiedRiskAxis();
    ~IQualifiedRiskAxis() {};
};

#ifndef QLIB_IQUALIFIEDRISKAXIS_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL IQualifiedRiskAxis<Void>);
EXTERN_TEMPLATE(class RISKMGR_DLL IQualifiedRiskAxis<ExpiryWindow>);
#endif

/**
 * A non-term-structured property of the world, with respect to which risk can
 * be estimated.
 */

typedef IQualifiedRiskAxis<Void> IScalarRiskAxis;
DECLARE(IScalarRiskAxis)

/**
 * A term-structured property of the world, with respect to which risk can be
 * estimated.
 */

typedef IQualifiedRiskAxis<ExpiryWindow> IExpiryRiskAxis;
DECLARE(IExpiryRiskAxis)

/**
 * A 2d-term-structured property which a market object (i.e CDS vol objects) can have.
 */

typedef IQualifiedRiskAxis<ExpiryPair> IExpiryPairRiskAxis;
DECLARE(IExpiryPairRiskAxis)

/**
 * A 2d-term-structured property which a market object (i.e vol surfs) can have.
 */

typedef IQualifiedRiskAxis<ExpiryAndStrike> IExpiryAndStrikeRiskAxis;
DECLARE(IExpiryAndStrikeRiskAxis)

DRLIB_END_NAMESPACE

#endif
