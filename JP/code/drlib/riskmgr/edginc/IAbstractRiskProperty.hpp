/**
 * @file IAbstractRiskProperty.hpp
 */

#ifndef DRLIB_IAbstractRiskProperty_H
#define DRLIB_IAbstractRiskProperty_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(OutputName)
FORWARD_DECLARE(IHypothesis)
FORWARD_DECLARE(IRiskAxis)
FORWARD_DECLARE(IAbstractRiskProperty)
FORWARD_DECLARE(RiskMappingMatrix)

/**
 * A property which a market object can have (Spartan base interface for
 * IRiskProperty)
 *
 * An IRiskProperty is a property which a market object can have and with
 * respect to which risk can be estimated; an IAbstractRiskProperty
 * just knows which objects have it, not how to generate the IRiskAxis
 * along which you can estimate risk for the name w.r.t. the property.
 *
 * There's no terribly compelling reason to split this methods off from
 * IRiskProperty, but it's handy not to have to templatize by the QUALIFIER
 * if we just want to call subjectNames().
 */

class RISKMGR_DLL IAbstractRiskProperty: public virtual IObject {
public:

    static CClassConstSP const TYPE;

    IAbstractRiskProperty();
    ~IAbstractRiskProperty();

    /**
     * See RiskProperty<PROPERTY>::discrete()
     */

    virtual bool discrete() const = 0;

    /**
     * See RiskProperty<PROPERTY>::subjectInterface()
     */

    virtual CClassConstSP subjectInterface() const = 0;

    /**
     * See RiskProperty<PROPERTY>::subjectNames()
     */

    virtual OutputNameArrayConstSP subjectNames(IObjectConstSP world) const = 0;

    /**
     * See RiskProperty<PROPERTY>::riskMappingMatrix()
     * Default implementation returns an empty SP
     */

    virtual RiskMappingMatrixConstSP riskMappingMatrix(IObjectConstSP world,
                                                       OutputNameConstSP name) const;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
