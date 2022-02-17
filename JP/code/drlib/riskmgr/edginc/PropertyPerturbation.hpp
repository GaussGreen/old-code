/**
 * @file PropertyPerturbation.hpp
 */

#ifndef QLIB_PropertyPerturbation_H
#define QLIB_PropertyPerturbation_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Perturbation.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Presents a RiskProperty for tweaking in a Scenario as an IPerturbation
 *
 * Essentially an adapter class allowing a PropertyTweakHypothesis to be
 * constructed via an IPerturbation (see findAndShift()).  Strictly speaking it
 * isn't necessary, since any IHypothesis is already an IScenarioShift and so
 * can be used in Scenario's directly without the need for a ScenarioShift /
 * IPerturbation combo.  However the slightly different structure would
 * apparently cause the Analytics guys a vast amount of work.
 */

template <class TAG>
class PropertyPerturbation: public CObject,
                            public virtual IPerturbation {

public:

    static CClassConstSP const TYPE;

    typedef TAG Tag;
    typedef smartConstPtr<Tag> TagConstSP;
    
    typedef typename TAG::Qualifier Qualifier;
    DECLARE(Qualifier)

private:

    static IObject* emptyShell() {
        return new PropertyPerturbation();
    }

    static void load(CClassSP& clazz) {
        REGISTER(PropertyPerturbation, clazz);
        clazz->setPublic();
        SUPERCLASS(CObject);
        IMPLEMENTS(IPerturbation);
        EMPTY_SHELL_METHOD(emptyShell);
        FIELD(qualifier, "qualifier");
        FIELD(shiftSize, "shiftSize");
    }

    QualifierConstSP qualifier;
    double shiftSize;

protected:

    PropertyPerturbation(CClassConstSP type,
                         double shiftSize,
                         QualifierConstSP qualifier = QualifierConstSP()):
        CObject(type),
        qualifier(qualifier),
        shiftSize(shiftSize)
    {}

    virtual TagConstSP tag() const {
        return TagConstSP();
    }

public:

    PropertyPerturbation(double shiftSize = 0,
                         QualifierConstSP qualifier = QualifierConstSP()):
        CObject(TYPE),
        qualifier(qualifier),
        shiftSize(shiftSize)
    {}

    ~PropertyPerturbation() {}

    /**
     * Implementation of IPerturbation
     */

    bool findAndShift(IObjectSP object, OutputNameConstSP name) {
        bool changed;
        PropertyTweakHypothesis<TAG>(shiftSize, name, qualifier, tag()).
            applyTo(object, &changed);
        return changed;
    }

private:
    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each array template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class TAG> CClassConstSP const 
PropertyPerturbation<TAG>::TYPE =
CClass::templateRegisterClass(typeid(PropertyPerturbation<TAG>));
#endif

DRLIB_END_NAMESPACE

#endif
