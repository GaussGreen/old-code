/**
 * @file FieldRiskAxis.hpp
 */

#ifndef QLIB_FieldRiskAxis_H
#define QLIB_FieldRiskAxis_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskAxis.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(OutputName)

/**
 * Designates a value, stored in a field on some market object, which can be
 * tweaked to construct scenarios or estimate greeks
 *
 * See the implementations ScalarFieldRiskAxis and ExpiryFieldRiskAxis for
 * details.
 *
 * For general information about the "flexible scenarios and greeks" framework,
 * see FlexibleSensitivity.
 */

class RISKMGR_DLL FieldRiskAxis: public CObject,
                                 public virtual IRiskAxis {

    static void load(CClassSP& clazz);
    FieldRiskAxis();
    static IObject* emptyShell();
    void validate();

public:

    static CClassConstSP const TYPE;

    const IAbstractRiskPropertyConstSP property;     // $required
    const OutputNameConstSP name;                    // $optional
    const IFieldTweakConstSP tweak;                  // $required
    const bool absoluteDistance;                     // $required

    /**
     * Constructor
     */

    FieldRiskAxis(IAbstractRiskPropertyConstSP property,
                  OutputNameConstSP name,
                  IFieldTweakConstSP tweak,
                  bool absoluteDistance);

    void validatePop2Object();

    ~FieldRiskAxis();

    /**
     * @name IRiskAxis implementation
     */

    //@{

    /**
     * The IRiskProperty which changes as we move along this risk axis
     * 
     * This is a FieldRiskProperty, specifying a field to tweak and how to
     * tweak it.
     */

    IAbstractRiskPropertyConstSP abstractProperty() const;

    /**
     * The name of the market object we tweak (or NULL to mean "all of the
     * relevant class")
     */

    OutputNameConstSP marketDataName() const;

    IHypothesisConstSP hypothesis(double coeff) const;

    /**
     * Version of this object suitable for easy database storage
     */

    RiskAxisConstSP frozen() const;

    //@}

    string toString() const;
};

DRLIB_END_NAMESPACE

#endif
