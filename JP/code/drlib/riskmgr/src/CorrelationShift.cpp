//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrelationShift.cpp
//
//   Description : Shifts a correlation between a pair of assets (FX-FX, 
//                 FX-equity or equity-equity)
//                 Parallel shifts all correlations if no name is passed.
//
//   Author      : Andrew McCleery
//
//   Date        : 20 Jan 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyPerturbation.hpp"
#include "edginc/Correl.hpp"

DRLIB_BEGIN_NAMESPACE

class CorrelationShift: public PropertyPerturbation<Correl> {
public:
    static CClassConstSP const TYPE;

    CorrelationShift(): PropertyPerturbation<Correl>(TYPE, 0.01) {}

private:

    CorrelConstSP tag() const {
        return Correl::SP(true, // asset correls please
                          true, // and fx ones
                          true, // and others
                          Correl::ABSOLUTE);
    }

    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CorrelationShift, clazz);
        SUPERCLASS(PropertyPerturbation<Correl>);
        EMPTY_SHELL_METHOD(DefaultConstructor<CorrelationShift>::iObject);
    }
};

CClassConstSP const CorrelationShift::TYPE = CClass::registerClassLoadMethod(
        "CorrelationShift", typeid(CorrelationShift), load);

// to force linker to include file (avoid having header file) */
bool CorrelationShiftLinkIn() {
    return (CorrelationShift::TYPE != 0);
}

DRLIB_END_NAMESPACE
