//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ImpliedYTM.cpp
//
//   Description : Implied Yield to Maturity
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : April 16, 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ImpliedYTM.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/TweakGroup.hpp"

DRLIB_BEGIN_NAMESPACE

const string ImpliedYTM::NAME = "IMPLIED_YIELD_TO_MATURITY";

const string& ImpliedYTM::getSensOutputName() const{
    return ImpliedYTM::NAME;
}

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns true */
bool ImpliedYTM::discreteShift() const{
    return true;
}

void ImpliedYTM::calculate(TweakGroup*      tweakGroup,
                           CResults*        results)
{
    static const string method("ImpliedYTM::calculate");

    try {
        double ytm;

        // cast to ImpliedYTM::IFaceYTM
        ImpliedYTM::IFaceYTM *YTMInst =
            dynamic_cast<ImpliedYTM::IFaceYTM *>(tweakGroup->getInstrument());
       
        if (YTMInst == 0) { // the cast failed
            results->storeNotApplicable(this);
        } else { // the instrument implements the interface
            ytm = YTMInst->yieldToMaturity(price, clean);
    
            /* Store implied shift */
            results->storeScalarGreek(ytm,
                                      ImpliedYTM::NAME,
                                      resultOutputName);
        }
    } catch (exception& e) {
        results->storeGreek( IObjectSP(new Untweakable(e)), 
                                       ImpliedYTM::NAME,
                                       resultOutputName);
    }

}

ImpliedYTM::ImpliedYTM(double        targetPrice, 
                       bool          quotedClean,
                       const string& resultOutputName):
    Sensitivity(TYPE),
    price(targetPrice),
    clean(quotedClean),
    resultOutputName(new OutputName(resultOutputName)) {}

/** for reflection */
ImpliedYTM::ImpliedYTM():Sensitivity(TYPE) {
    clean = true;
}

void ImpliedYTM::validatePop2Object()
{
    static const string method("ImpliedYTM::validatePop2Object");
}

class ImpliedYTMHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ImpliedYTM, clazz);
        SUPERCLASS(Sensitivity);
        EMPTY_SHELL_METHOD(defaultImpliedYTM);
        FIELD(price, "price of bond");
        FIELD(clean, "Is the bond quoted clean? Default is true.");
        FIELD_MAKE_OPTIONAL(clean);
        FIELD(resultOutputName, "name under which the result is reported back");
}

    static IObject* defaultImpliedYTM(){
        return new ImpliedYTM();
    }
};

CClassConstSP const ImpliedYTM::TYPE = CClass::registerClassLoadMethod(
    "ImpliedYTM", typeid(ImpliedYTM), ImpliedYTMHelper::load);

DRLIB_END_NAMESPACE
