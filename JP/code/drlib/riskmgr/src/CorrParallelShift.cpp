//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrParallelShift.cpp
//
//   Description : Parallel shift to all correlations
//
//   Author      : Andrew J Swain
//
//   Date        : 17 Jan 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ScenarioShift.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/Correl.hpp"

DRLIB_BEGIN_NAMESPACE

class CorrParallelShift: public CObject,
                         public virtual IScenarioShift,
                         public virtual IPerturbation{ /* IPerturbation for 
                                                          backward
                                                          compatibility only */
public:
    static CClassConstSP const TYPE;
    const static double DEFAULT_SHIFT;

    //// shift all correlations together - just like PhiParallel
    virtual bool applyScenario(IObjectSP object){
        bool changed;
        PropertyTweakHypothesis<Correl>(
            shiftSize,
            OutputNameSP(), // all names
            VoidSP(),
            Correl::SP(true, // asset correls please
                       true, // and fx too
                       true, // and others
                       Correl::ABSOLUTE)).
                applyTo(object, &changed);
        return changed;
    }
        
    // Nothing to do before market data is retrieved.
    virtual bool preapplyScenario(IObjectSP object){
        return false;
    }

    /** IPerturbation implementation - for backwards compatibility only.
        Equivalent here to applyScenario(objectToShift) */
    virtual bool findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name){
        return applyScenario(objectToShift);
    }

private:

    // **** registered fields ****
    double shiftSize;

    /** for reflection */
    CorrParallelShift(): CObject(TYPE), shiftSize(DEFAULT_SHIFT){}

    CorrParallelShift(const CorrParallelShift &rhs);
    CorrParallelShift& operator=(const CorrParallelShift& rhs);

    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CorrParallelShift, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IScenarioShift);
        IMPLEMENTS(IPerturbation); // for backwards compat
        EMPTY_SHELL_METHOD(defaultCorrParallelShift);
        FIELD(shiftSize, "How big to make the tweak");
    }
    
    static IObject* defaultCorrParallelShift(){
        return new CorrParallelShift();
    }
};

const double CorrParallelShift::DEFAULT_SHIFT = 0.01;

CClassConstSP const CorrParallelShift::TYPE = CClass::registerClassLoadMethod(
        "CorrParallelShift", typeid(CorrParallelShift), load);


// to force linker to include file (avoid having header file) */
bool CorrParallelShiftLinkIn() {
    return (CorrParallelShift::TYPE != 0);
}

DRLIB_END_NAMESPACE
