//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ImpliedScalarShift.cpp
//
//   Description : Implied Scalar Shift Sensitivity
//
//   Author      : regis Guichard
//
//   Date        : 10 Dec 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/ImpliedScalarShift.hpp"
#include "edginc/IRiskyPricer.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** minimal value for the positive shifts */
static const double EPSILON = 0.001;

/** Class creating the objective function */
class ImpliedScalarShift::ObjectiveFunc: public Func1D::NoDeriv{
    TweakGroup*           tweakGroup;
    IScalarPerNameShiftSP sensToShift;
    OutputNameConstSP     nameToShift;
    double                targetValue;
    CControl*             locCtrl;
    OutputNameConstSP     outputName;
    bool                  isPrice;
    string                packetName;
    bool                  isPositive;
public:
    virtual double operator()(double  x) const{
        try{
            /* if x has to be strictly positive then when put x to
               0+epsilon to prevent failure for 0 used mainly for
               volatility */
            if (isPositive) {
                x = max(EPSILON,x);
            }
            
            /* Shift */

            IHypothesis::AlternateWorldSP alt = sensToShift->appliedTo(
                nameToShift, x, IObjectSP::attachToRef(tweakGroup));
            
            /* Compute Price */
            CResultsSP results;
            try{
                TweakGroupSP tg = TweakGroupSP::dynamicCast(alt->world);
                CInstrument* inst = tg->getInstrument();
                results = CResultsSP(tg->getModel()->Run(inst, locCtrl));
            } catch( exception& ){
                alt->undo(); // must restore if failed
                throw;
            }
            /* Un-Shift */
            alt->undo();
            
            /* Grab whatever sensitivity we ask to be computed */
            double value;
            if (isPrice){
                value = results->retrievePrice();
            }
            else{
                value = results->retrieveScalarGreek(packetName,
                                                     outputName);
            }
            
            return value - targetValue;
        }
        catch(exception& e){
            throw ModelException(e, "ImpliedScalarShift::ObjectiveFunc::"
                                 "operator()");
        }
    }

    ObjectiveFunc(TweakGroup*          tweakGroup,
                  IScalarPerNameShift* sensToShift,
                  OutputNameConstSP    nameToShift,
                  double               targetValue,
                  CControl*            locCtrl,
                  OutputNameConstSP    outputName,
                  bool                 isPositive):
        tweakGroup(tweakGroup), 
        sensToShift(copy(sensToShift)), nameToShift(nameToShift),
        targetValue(targetValue),
        locCtrl(locCtrl),
        outputName(outputName),
        isPositive(isPositive) {
        SensitivityArrayConstSP sensArray(locCtrl->getSens());
        OutputRequestArrayConstSP requestArray(locCtrl->getOutputRequests());
        if (sensArray->size() == 1){
            isPrice = false;
            packetName = (*sensArray)[0]->getSensOutputName();
        }
        else if (requestArray->size() == 1){
            isPrice = false;
            packetName = (*requestArray)[0]->getPacketName();
            /* if the packetname does not equal the request name,
               we're using the output name that has been provided
               externally */
            if ( packetName != (*requestArray)[0]->getRequestName()) {
                this->outputName = OutputNameSP(
                    new OutputName((*requestArray)[0]->getRequestName()));
            }
        } else {
            isPrice = true;
        }
    }  
};

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns true */
bool ImpliedScalarShift::discreteShift() const{
    return true;
}

const string ImpliedScalarShift::NAME = "IMPLIED_SHIFT";

const string& ImpliedScalarShift::getSensOutputName() const{
    return ImpliedScalarShift::NAME;
}

void ImpliedScalarShift::addResult(Results*          results,
                                   const Results*    resultsToAdd,
                                   double            scaleFactor) const{
    try{
        results->merge(getSensOutputName(), resultsToAdd);
    } catch (exception& e){
        throw ModelException(e, "ImpliedScalarShift::addResult");
    }
}

/** calculate a range containing the root to find */
pair<double,double> ImpliedScalarShift::calculateRange(
                                          const Func1D::NoDeriv & objFunc,
                                          const double &          initialGuess,
                                          const double &          lBound,
                                          const double &          hBound,
                                          const bool &            isPositive) 
{
    double lShiftSize, hShiftSize;
    static const string method("ImpliedScalarShift::calculateRange");
    try {
        /* if an input bracket has been specified */
        if (!(Maths::isZero(lBound) && Maths::isZero(hBound))) {
            lShiftSize = lBound;
            hShiftSize = hBound;
            if ( (objFunc(lShiftSize)* objFunc(hShiftSize) )>= 0.) {
                throw ModelException(method, 
                                       "The solution is not in the input range ["
                                     + Format::toString(lShiftSize) 
                                     + "," + Format::toString(hShiftSize) +"]");
            }
        } else {
            try {
                /* Bracket roots in R+ */
                if (isPositive) {
                    hShiftSize = 1.1 * initialGuess;
                    lShiftSize = 0.;
                    ZBracPositive_bracket(objFunc,
                                          lShiftSize,
                                          hShiftSize,
                                          true);
                }
                /* Bracket roots in R */
                else {
                    lShiftSize = 0.9 * initialGuess;
                    hShiftSize = 1.1 * initialGuess;
                    ZBrac_bracket(objFunc, 
                                  lShiftSize, 
                                  hShiftSize);
                }
            } catch(exception& ){
                if (isPositive) {
                    /* something that should be done in ZBracPositive_bracket */
                    lShiftSize = (1.1*initialGuess*0.5) / (2.0 - (lShiftSize/(1.1 * initialGuess*0.5)));
                }
                throw ModelException(method, "Can not find a solution in the range ["
                                     + Format::toString(lShiftSize) +","
                                     + Format::toString(hShiftSize) +"]");
            }
        }
    }
    catch(exception& e){
        throw ModelException(e,method);
    }
    return pair<double,double>(lShiftSize,hShiftSize);
}


/* NOTE that this will not work perfectly for implied PHI shifts - as PHI will always shift 
   towards 0.0 correlation, ie. the resulting shift size may have the wrong sign or no root 
   can be found, although there could exist one. There is no easy fix for this. Since 
   IMPLIED_PHI is probably a rarely requested shift, nothing is done here to fix it. */
void ImpliedScalarShift::calculate(TweakGroup*      tweakGroup,
                                   CResults*        results)
{
    static const string method("ImpliedScalarShift::calculate");
    IRiskyPricer* riskyInst = NULL;
    DateTime    effMaturityDate;

    try {
        // get list of names to calculate result for

        ScalarPerturbation* sens = dynamic_cast<ScalarPerturbation *>(sensToShift.get());
        OutputNameArrayConstSP names;
        if (!!sens) {
            names = OutputName::trim(sens->names(tweakGroup));
        } else { 
            names = OutputName::trim(sensToShift->allNames(tweakGroup));
        }

        if (names->size() > 1 ) {
            throw ModelException(
                method, "Unable to calculate implied sensitivity. "
                "There is more \n"
                " than one name which is sensitive to " + 
                sensToShift->getClass()->getName());
        } 
        
        if (names->size() > 0) {
            // If implied sensitivity hasn't already been calculated and stored
            if (!results->exists(getSensOutputName(), 
                                 resultOutputName)){
                /* Create objective function */
                ObjectiveFunc objFunc(tweakGroup,
                                      sensToShift.get(),
                                      names->front(),
                                      targetValue,
                                      locCtrl.get(),
                                      targetControlOutputName,
                                      isPositive);
            
                        
                //Calculate initial guess
                double initialGuess = 0.0;
                                
                if (!!initialGuessRequest) {
                    /* if an initial guess request has been specified */
                    CControlSP controlLocal(Control::makeFromFlags("",0.0));
                    controlLocal->addRequest(initialGuessRequest);
                    CResultsSP results;
                    CInstrument* inst = tweakGroup->getInstrument();
                    results = CResultsSP(tweakGroup->getModel()->Run(inst,controlLocal.get()));
                    OutputNameSP outputName(new OutputName(initialGuessRequest->getRequestName()));
                    initialGuess =  results->retrieveScalarGreek(initialGuessRequest->getPacketName(),
                                                                 outputName);
                } else {
                    /* otherwise we use the initial guess in sensToShift, it has to be at least positive or null */
                    initialGuess = sensToShift->getShiftSize();     
                }

                /* create a range containing the root of the 
                   objective function */
                 pair<double,double> shiftSizes = calculateRange(objFunc,
                                                                 initialGuess,
                                                                 lBound,
                                                                 hBound,
                                                                 isPositive);
                double lShiftSize = shiftSizes.first;
                double hShiftSize = shiftSizes.second;

                /* Use input bracket */
                Instrument * inst = tweakGroup->getInstrument();
                if (IRiskyPricer::TYPE->isInstance(inst))
                {
                    DateTime valueDate = inst->getValueDate();
                    riskyInst = dynamic_cast<IRiskyPricer*>(inst);
                    effMaturityDate = riskyInst->getEffMaturityDate();
                    /* Set effective maturity for building risky curve to present 
                       date to avoid failures in implied spread calculation */
                    riskyInst->setEffMaturityDate(valueDate);
                }
                /* Solve root */
                double shiftSize = rootFinder->solve(objFunc,
                                                     lShiftSize,
                                                     hShiftSize);
                
                /* Store implied shift */
                results->storeScalarGreek(shiftSize,
                                          getSensOutputName(),  // packet
                                          resultOutputName);
                if (riskyInst)
                    riskyInst->setEffMaturityDate(effMaturityDate);
            }
        } else {
            throw ModelException(method, 
                                 "Could not find any objects which implement " + sensToShift->getClass()->getName());
        }
    }
    catch (exception& e) {
        try {
            results->storeGreek( IObjectSP(new Untweakable(e)), 
                                 getSensOutputName(),  // packet
                                 resultOutputName);
            if (riskyInst)
                riskyInst->setEffMaturityDate(effMaturityDate);
        } catch (exception& ) {
        }
    }
}

/** for reflection */
ImpliedScalarShift::ImpliedScalarShift():
    Sensitivity(TYPE),
    lBound(0.),
    hBound(0.),
    isPositive(false){
}

void ImpliedScalarShift::validatePop2Object()
{
    static const string method("ImpliedScalarShift::validatePop2Object");
    
    SensitivityArrayConstSP sensArray(new SensitivityArray());
    OutputRequestArrayConstSP requestArray(new OutputRequestArray());
    
    if (!!targetControl){
        sensArray = targetControl->getSens();
        requestArray = targetControl->getOutputRequests();

        /* At most 1 sensitivity / request. */
        if (sensArray->size() + requestArray->size() > 1){
            throw ModelException(method,
                                 "targetControl should include only one sensitivity / request; got "
                                 + Format::toString(sensArray->size())
                                 + " and "
                                 + Format::toString(requestArray->size())
                                 + ", respectively");
        }    

        /* If control has 1 sensitivity / request, 
           need to provide market name */
        /* AS: I have taken this out for the time being (possibly), as
           it does not allow for any sensitivities in the
           Instrument packet, for. example IND_VOL, which does not
           require an output name. Proper validation here is quite
           tricky and should probably be left to fail downstream
           (failure in extracting results) if (sensArray->size() +
           requestArray->size() != 0){ if
           (!targetControlOutputName){ throw
           ModelException(method, "targetControlOutputName should
           be provided together with targetControl"); } } */
    }

    locCtrl = CControlSP(new Control(sensArray,
                                     requestArray,
                                     false,
                                     ""));
}

/** explicit constructor which is used to calculate ImpliedVol sensitivity.
    Target control is always assumed to be VALUE                            
    the default mapping is the identity function */
ImpliedScalarShift::ImpliedScalarShift(
    IScalarPerNameShiftSP              sensToShift,
    double                             targetValue,
    RootFinder1D::TwoInitValNoDerivSP  rootFinder,
    OutputNameSP                       resultOutputName):
    Sensitivity(TYPE),
    sensToShift(sensToShift), 
    targetValue(targetValue), 
    rootFinder(rootFinder), 
    resultOutputName(resultOutputName),
    lBound(0.),
    hBound(0.),
    isPositive(false)
{
    validatePop2Object();
}

/** constrcutor allowing to specify if the value to shift has to be always positive (i.e. volatility)
    and an request for the initial value (i.e. indicative volatility) */
ImpliedScalarShift::ImpliedScalarShift(
    IScalarPerNameShiftSP              sensToShift,
    double                             targetValue,
    RootFinder1D::TwoInitValNoDerivSP  rootFinder,
    OutputNameSP                       resultOutputName,
    bool                               isPositive,
    OutputRequestSP                    initialGuessRequest):
    Sensitivity(TYPE),
    sensToShift(sensToShift), 
    targetValue(targetValue), 
    rootFinder(rootFinder), 
    resultOutputName(resultOutputName),
    lBound(0.),
    hBound(0.),
    isPositive(isPositive),
    initialGuessRequest(initialGuessRequest)
{
    validatePop2Object();
}

class ImpliedScalarShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ImpliedScalarShift, clazz);
        SUPERCLASS(Sensitivity);
        EMPTY_SHELL_METHOD(defaultImpliedScalarShift);
        FIELD(sensToShift, "Sensitivity whose shift size we seek to solve for (Inputted shift size will be used as initial guess)");
        FIELD(targetValue, "Target Value");
        FIELD(targetControl, "Control that the target value refers to (eg delta). Price, if omitted");
        FIELD_MAKE_OPTIONAL(targetControl);
        FIELD(targetControlOutputName, "Market name associated with target control");
        FIELD_MAKE_OPTIONAL(targetControlOutputName);
        FIELD(rootFinder, "Root finder (Use ImpliedScalarShift::DefaultRootFinder by default)");
        FIELD(resultOutputName, "name under which the result is reported back");
        FIELD(lBound, "initial low bound for result");
        FIELD_MAKE_OPTIONAL(lBound);
        FIELD(hBound, "initial high bound for result");
        FIELD_MAKE_OPTIONAL(hBound);
        FIELD(locCtrl, "locCtrl");
        FIELD_MAKE_TRANSIENT(locCtrl);
        FIELD(isPositive, "has the domain of the result to be positive?");
        FIELD_MAKE_OPTIONAL(isPositive);
        FIELD(initialGuessRequest, "Control that the initial guess value refers to (eg delta).");
        FIELD_MAKE_OPTIONAL(initialGuessRequest);
    }

    static IObject* defaultImpliedScalarShift(){
        return new ImpliedScalarShift();
    }
};

CClassConstSP const ImpliedScalarShift::TYPE = CClass::registerClassLoadMethod(
    "ImpliedScalarShift", typeid(ImpliedScalarShift), ImpliedScalarShiftHelper::load);


ImpliedScalarShift::DefaultRootFinder::DefaultRootFinder(double shiftAccuracy):
    ZBrent(shiftAccuracy){}

ImpliedScalarShift::DefaultRootFinder::DefaultRootFinder():
    ZBrent(TYPE){}

class ImpliedScalarShiftDefaultRootFinderHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ImpliedScalarShift::DefaultRootFinder, clazz);
        SUPERCLASS(ZBrent);
        EMPTY_SHELL_METHOD(defaultCtor);
    }

    static IObject* defaultCtor(){
        return new ImpliedScalarShift::DefaultRootFinder();
    }
};

CClassConstSP const ImpliedScalarShift::DefaultRootFinder::TYPE = CClass::registerClassLoadMethod(
    "ImpliedScalarShift::DefaultRootFinder", typeid(ImpliedScalarShift::DefaultRootFinder), ImpliedScalarShiftDefaultRootFinderHelper::load);

DRLIB_END_NAMESPACE
