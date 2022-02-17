
#ifndef EDG_SCALAR_SHIFT_TWO_SIDED_HPP
#define EDG_SCALAR_SHIFT_TWO_SIDED_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/Additive.hpp"
#include "edginc/TwoSidedDeriv.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/IRiskQuantityFactory.hpp"
#include "edginc/ICompatibilitySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(AbstractPropertyTweakHypothesis)
FORWARD_DECLARE(RiskMapping)

/** Specialised version of ScalarShift which implements a two sided calcuation
    of the derivative */

FORWARD_DECLARE(ScalarShiftTwoSided)

class RISKMGR_DLL ScalarShiftTwoSided: public ScalarShift,
                           public virtual Additive,
                           public virtual ITwoSidedDeriv,
                           public virtual IRiskQuantityFactory,
                           public virtual ICompatibilitySensitivity {
public:    

    static CClassConstSP const TYPE;

    virtual ~ScalarShiftTwoSided();

    /** identifies the name used for storing the second order derivative 
        in the output*/
    virtual const string& getSecondOrderSensOutputName() const = 0;

    /** Overridden to perform two sided tweak */
    virtual void calculate(TweakGroup*  tweakGroup,
                           Results*     results);

    double divisor() const;

    /** Overridden to scale 1st and 2nd order derivatives */
    virtual void scaleResult(Results*     results,
                             double       scaleFactor) const;

    /** Overridden to add 1st and 2nd order derivatives */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;

    /** Returns the shift(s) which have been made for the current pricing
        call */
    virtual ScalarShiftArray getComponentShifts() const;

    virtual CClassConstSP shiftInterface() const;

    virtual CClassConstSP restorableShiftInterface() const;

    virtual OutputNameArrayConstSP names(const IObject* tweakGroup) const;

    virtual bool nameMatches(const OutputName& name, IObjectConstSP obj);

    virtual void appendName(OutputNameArray& namesList, IObjectConstSP obj);

    virtual bool shift(IObjectSP obj);

    virtual void restore(IObjectSP obj);

    virtual void setShiftSize(double);

    virtual void setMarketDataName(OutputNameConstSP);

    virtual NamedRiskQuantityArraySP riskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const;

    virtual LazyRiskQuantityFactoryArraySP lazies(MultiTweakGroupConstSP world) const;

protected:
    /** Note ScalarShiftTwoSided is abstract. Create a scalar shift of
        type clazz, which uses outputName (eg VEGA_PARALLEL) to
        identify results and with given shiftSize */
    ScalarShiftTwoSided(CClassConstSP clazz,
                        bool constantDivisor,
                        double sensitivityUnit,
                        const string& outputName,
                        const IScalarRiskProperty* property,
                        double coefficient);

private:
    static void load(CClassSP& clazz);

    bool _constantDivisor;
    double _sensitivityUnit;

    IScalarRiskPropertyConstSP property;

    // copy-on-write arrangement for currentHypothesis

    AbstractPropertyTweakHypothesisConstSP _currentHypothesisConst;
    AbstractPropertyTweakHypothesisSP _currentHypothesisMutable;
    AbstractPropertyTweakHypothesisSP currentHypothesisMutable();
    AbstractPropertyTweakHypothesisConstSP currentHypothesis() const;

    // called from HypothesisTree in IRiskQuantityFactory.cpp

    void setCurrentHypothesis(AtomicHypothesisArraySP history,
                              AbstractPropertyTweakHypothesisConstSP current,
                              double oldValue);

    ScalarShiftTwoSidedConstSP alteredControl(
        MultiTweakGroupConstSP tweakGroup) const;
};

DRLIB_END_NAMESPACE

#endif
