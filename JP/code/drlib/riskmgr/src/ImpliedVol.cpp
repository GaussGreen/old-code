//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ImpliedVol.cpp
//
//   Description : Implied Vol Sensitivity
//
//   Author      : André Segger
//
//   Date        : 23 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ImpliedVol.hpp"
#include "edginc/VolLevel.hpp"
#include "edginc/ImpliedScalarShift.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns true */
bool ImpliedVol::discreteShift() const{
    return true;
}

const string ImpliedVol::NAME = "IMPLIED_VOL";

const string& ImpliedVol::getSensOutputName() const{
    return ImpliedVol::NAME;
}

void ImpliedVol::calculate(TweakGroup*     tweakGroup,
                           Results*        results)
{
    static const string method("ImpliedVol::calculate");
    
    try {
       // if the initialGuess is 0% then we use the indicative volatility of the instrument 
       // as initial guess 
       OutputRequestSP initialGuessRequest(new OutputRequest(OutputRequest::IND_VOL));
             
       // create ImpliedScalarShift sensitivity with VolLevel scenario
       IScalarPerNameShiftSP volScenario(new VolLevel(initialGuess));
       OutputNameSP          impliedName(new OutputName(ImpliedVol::NAME));
       Results               tmpResults;
       ImpliedScalarShiftSP  impliedShift(new ImpliedScalarShift(
                                              volScenario,
                                              targetValue,
                                              rootFinder,
                                              impliedName,
                                              true,                 // isPositive
                                              initialGuessRequest));           
       
       // calculate implied vol
       impliedShift->calculate(tweakGroup, &tmpResults);
       
       IObjectConstSP result =  tmpResults.retrieveGreek(
           impliedShift->getSensOutputName(),
           impliedName);
       
       /* Store implied shift */
       results->storeGreek(IObjectSP(result->clone()),
			   ImpliedVol::NAME,
			   resultOutputName);
       
    } catch (exception& e) {
       results->storeGreek(IObjectSP(new Untweakable(e)), 
			   ImpliedVol::NAME,
			   resultOutputName);
    }
    
}

ImpliedVol::ImpliedVol(double        target, 
                       double        guess, 
                       double        tolerance, 
                       const string& resultOutputName):
    Sensitivity(TYPE),
    targetValue(target),
    initialGuess(guess),
    rootFinder(new ImpliedScalarShift::DefaultRootFinder(tolerance)),
    resultOutputName(new OutputName(resultOutputName)) {}

/** for reflection */
ImpliedVol::ImpliedVol():Sensitivity(TYPE) {}

void ImpliedVol::validatePop2Object()
{
    static const string method("ImpliedVol::validatePop2Object");
}

class ImpliedVolHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ImpliedVol, clazz);
        SUPERCLASS(Sensitivity);
        EMPTY_SHELL_METHOD(defaultImpliedVol);
        FIELD(targetValue, "Target Value");
        FIELD(initialGuess, "Initial Guess");
        FIELD(rootFinder, "Root finder (Use ImpliedScalarShift::DefaultRootFinder by default)");
        FIELD(resultOutputName, "name under which the result is reported back");
    }

    static IObject* defaultImpliedVol(){
        return new ImpliedVol();
    }
};

CClassConstSP const ImpliedVol::TYPE = CClass::registerClassLoadMethod(
    "ImpliedVol", typeid(ImpliedVol), ImpliedVolHelper::load);

DRLIB_END_NAMESPACE
