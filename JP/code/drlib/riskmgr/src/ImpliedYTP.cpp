//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ImpliedYTP.cpp
//
//   Description : Implied Yield to First Put
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : April 16, 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ImpliedYTP.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/TweakGroup.hpp"

DRLIB_BEGIN_NAMESPACE

const string ImpliedYTP::NAME = "IMPLIED_YIELD_TO_FIRST_PUT";

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns true */
bool ImpliedYTP::discreteShift() const{
    return true;
}

const string& ImpliedYTP::getSensOutputName() const{
    return ImpliedYTP::NAME;
}

void ImpliedYTP::calculate(TweakGroup*      tweakGroup,
                           CResults*        results)
{
    static const string method("ImpliedYTP::calculate");

    try {
        double ytp;

        // cast to ImpliedYTP::IFaceYTP
        ImpliedYTP::IFaceYTP *YTPInst = 
            dynamic_cast<ImpliedYTP::IFaceYTP *>(tweakGroup->getInstrument());
       
        if (YTPInst == 0) { // the cast failed
            results->storeNotApplicable(this);
        } else { // the instrument implements the interface
            if ( YTPInst->hasPut() ) {
                ytp = YTPInst->yieldToFirstPut(price, clean);
    
                /* Store implied shift */
                results->storeScalarGreek(ytp,
                                          ImpliedYTP::NAME,
                                          resultOutputName);
            } else {
                results->storeNotApplicable(this);
            }
        }
    } catch (exception& e) {
        results->storeGreek( IObjectSP(new Untweakable(e)), 
                                       ImpliedYTP::NAME,
                                       resultOutputName);
    }

}

ImpliedYTP::ImpliedYTP(double        targetPrice, 
                       bool          quotedClean,
                       const string& resultOutputName):
    Sensitivity(TYPE),
    price(targetPrice),
    clean(quotedClean),
    resultOutputName(new OutputName(resultOutputName)) {}

/** for reflection */
ImpliedYTP::ImpliedYTP():Sensitivity(TYPE) {
    clean = true;
}

void ImpliedYTP::validatePop2Object()
{
    static const string method("ImpliedYTP::validatePop2Object");
}

class ImpliedYTPHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ImpliedYTP, clazz);
        SUPERCLASS(Sensitivity);
        EMPTY_SHELL_METHOD(defaultImpliedYTP);
        FIELD(price, "price of bond");
        FIELD(clean, "Is the bond quoted clean? Default is true.");
        FIELD_MAKE_OPTIONAL(clean);
        FIELD(resultOutputName, "name under which the result is reported back");
}

    static IObject* defaultImpliedYTP(){
        return new ImpliedYTP();
    }
};

CClassConstSP const ImpliedYTP::TYPE = CClass::registerClassLoadMethod(
    "ImpliedYTP", typeid(ImpliedYTP), ImpliedYTPHelper::load);

DRLIB_END_NAMESPACE
