//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : ImpliedScalarShiftMulti.cpp
//
//   Description : Like ImpliedScalarShift, but implemented as a "new greek".
//
//   Author      : Linus Thand
//
//   Date        : 09 June 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ImpliedScalarShiftMulti.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/NamedRiskObjectQuantity.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/RiskMapping.hpp"
#include "edginc/IRiskyPricer.hpp"
#include "edginc/Results.hpp"
#include "edginc/ScalarPerturbation.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/TRACE.hpp"


DRLIB_BEGIN_NAMESPACE

/** minimal value for the positive shifts */
static const double EPSILON = 0.001;
static const string atSolutionIdentifierDefault = "RESULTS_AT_SOLUTION";


double sumMultiTweakGroup(MultiTweakGroupConstSP tweakGroup,
                          CControlSP             locCtrl,
                          OutputNameConstSP      outputName,
                          const bool             isPrice,
                          const string&          packetName) {
    TRACE("Pricing");
    IInstrumentCollectionSP insts = tweakGroup->getInstruments();
    CResultsArraySP results = insts->emptyResults(); //// Get a list of empty results.
    insts->Price(tweakGroup->getModel(), locCtrl.get(), results);
    double value = 0;
    for (int i = 0 ; i < results->size();  ++i) {
        if (isPrice){
            value += (*results)[i]->retrievePrice();
        } else {
            value += (*results)[i]->retrieveScalarGreek(packetName, outputName);
        }
    }
    return value;
}

class ImpliedScalarShiftMulti::ObjectiveFunc : public Func1D::NoDeriv {
 private:
    MultiTweakGroupSP     tweakGroup;
    IScalarPerNameShiftSP sensToShift;
    OutputNameConstSP     nameToShift;
    const double          targetValue;
    CControlSP            locCtrl;
    OutputNameConstSP     outputName;
    const bool            isPositive;
    string                packetName;
    bool                  isPrice;
 public:
    virtual double operator()(double x) const{
        try {
            /** if x has to be strictly positive then when put x to
               0+epsilon to prevent failure for 0 used mainly for
               volatility **/
            if (isPositive) { x = max(EPSILON, x); }            
            //// Shift instruments 
            IHypothesis::AlternateWorldSP alt = 
                sensToShift->appliedTo(nameToShift, 
                                       x, 
                                       tweakGroup);
            try {
                
                const double value = 
                    sumMultiTweakGroup(MultiTweakGroupSP::dynamicCast(alt->world),
                                       locCtrl,
                                       outputName,
                                       isPrice,
                                       packetName);
                //// Un-Shift instruments 
                alt->undo();
                return value - targetValue;

           } catch(exception&){
                alt->undo(); //// must restore even if failed
                throw;
           }

        } catch(exception& e){
           throw ModelException(e, "ImpliedScalarShiftMulti::ObjectiveFunc::"
                                   "operator()");
        }
    }

    ObjectiveFunc(MultiTweakGroupSP     tweakGroup,
                  IScalarPerNameShiftSP sensToShift,
                  OutputNameConstSP     nameToShift,
                  const double          targetValue,
                  CControlSP            locCtrl,
                  OutputNameConstSP     outputName,
                  const bool            isPositive):
        tweakGroup(tweakGroup), 
        sensToShift(sensToShift.clone()), 
        nameToShift(nameToShift),
        targetValue(targetValue),
        locCtrl(locCtrl),
        outputName(outputName),
        isPositive(isPositive){
        SensitivityArrayConstSP sensArray(locCtrl->getSens());
        OutputRequestArrayConstSP requestArray(locCtrl->getOutputRequests());

        if (sensArray->size() == 1){
            isPrice = false;
            packetName = (*sensArray)[0]->getSensOutputName();
        }
        else if (requestArray->size() == 1){
            isPrice = false;
            packetName = (*requestArray)[0]->getPacketName();
            /** if the packetname does not equal the request name,
               we're using the output name that has been provided
               externally **/
            if (packetName != (*requestArray)[0]->getRequestName()) {
                outputName = OutputNameSP(
                    new OutputName((*requestArray)[0]->getRequestName()));
            }
        } else {
            isPrice = true;
        }
    }  
};

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns true */
bool ImpliedScalarShiftMulti::discreteShift() const { return true; }

const string ImpliedScalarShiftMulti::NAME = "IMPLIED";

/* NOTE that this will not work perfectly for implied PHI shifts - as PHI will 
   always shift towards 0.0 correlation, ie. the resulting shift size may have 
   the wrong sign or no root can be found, although there could exist one. There 
   is no easy fix for this. Since IMPLIED_PHI is probably a rarely requested 
   shift, nothing is done here to fix it. */

NamedRiskQuantityArraySP ImpliedScalarShiftMulti::nameRiskQuantities(
         MultiTweakGroupConstSP world, 
         RiskMappingConstSP riskMapping) const {
    static const string method("ImpliedScalarShiftMulti::nameRiskQuantities");
    NamedRiskQuantityArraySP rqs(new NamedRiskQuantityArray());
    MultiTweakGroupSP tweakgroup(world.clone());

    /** From ImpliedScalarShiftMulti, 
    don't know how to handle this case so bail out. **/
    for (int i = 0; i < tweakgroup->getInstruments()->size(); ++i) {
        Instrument * inst = (*tweakgroup->getInstruments())[i].get();
        if (IRiskyPricer::TYPE->isInstance(inst)) {
            throw ModelException("ImpliedScalarShiftMulti not implemented"
                                 " for RiskyPricer");
        }
    }

    try {
        //// Get list of names to calculate result for
        ScalarPerturbation* sens = 
            dynamic_cast<ScalarPerturbation*>(sensToShift.get());
        OutputNameArrayConstSP names;
        if (!!sens) {
            names = OutputName::trim(sens->names(tweakgroup.get()));
        } else { 
            names = OutputName::trim(sensToShift->allNames(tweakgroup.get()));
        }

        if (!shiftAllNames && (names->size() > 1)) {
            throw ModelException(
                method, "Unable to calculate implied sensitivity. "
                "There is more \n"
                " than one name which is sensitive to " + 
                sensToShift->getClass()->getName());
        } 

        if (shiftAllNames || (names->size() > 0)) {
            OutputNameConstSP shiftName;
            if (shiftAllNames) {
                shiftName = OutputNameConstSP(0);
            } else {
                shiftName = names->front();
            }
            ObjectiveFunc objFunc(tweakgroup,
                                  sensToShift,
                                  shiftName,
                                  targetValue,
                                  locCtrl,
                                  targetControlOutputName,
                                  isPositive);

            //// Create a range containing the root of the objective function 
            TRACE("Calculating range");
            pair<double, double> shiftSizes = 
                ImpliedScalarShift::calculateRange(objFunc,
                                                   sensToShift->getShiftSize(),
                                                   lBound,
                                                   hBound,
                                                   isPositive);
            const double lShiftSize = shiftSizes.first;
            const double hShiftSize = shiftSizes.second;
            TRACE(" Done.");
            TRACE("Solving");

            //// Solve root 
            double shiftSize = rootFinder->solve(objFunc, lShiftSize, hShiftSize); 
                    
            //// Store implied shift 
            rqs->push_back(NamedRiskQuantity::SP(
                RiskQuantity::constant(shiftSize),
                IResultsIdentifier::SP(getPacketName(), 
                                       resultOutputName->toString()))); 
            TRACE(" Done.");

            //// Calculate outputs at solution point, if requested
            if (!!ctrlAtSolution) 
            { 
                /** Move back to the solution point
                 *  Note that we don't undo this operation. No need since
                 *  we're working on a copy anyway 
                 */ 
                TRACE("Calculating requests at solution point.");
                IHypothesis::AlternateWorldSP alt = 
                    sensToShift->appliedTo(shiftName, 
                                           shiftSize, 
                                           tweakgroup);
                   
                MultiTweakGroupSP localTG =
                    MultiTweakGroupSP::dynamicCast(alt->world);
                
                //// Get a list of empty results.
                CResultsArraySP resultss = localTG->getInstruments()->emptyResults(); 
                
                //// Calculate price, output requests and sensitivities
                ctrlAtSolution->calculateMulti(localTG->getModel(), 
                                          localTG->getInstruments(), 
                                          resultss);
      
                //// Store resultss
                /*
                rqs->push_back(NamedRiskQuantity::SP( 
                    RiskQuantity::constantIObject(IObjectSP(resultss)),
                    IResultsIdentifier::SP(getPacketName(), 
                                           atSolutionIdentifier)));
*/
                rqs->push_back(NamedRiskObjectQuantity::SP( 
                    IObjectSP(resultss),
                    IResultsIdentifier::SP(getPacketName(), 
                                           atSolutionIdentifier)));
                TRACE("Done.");
            }
        } else {
            throw ModelException(method, 
                                 "Could not find any objects which implement " 
                                 + sensToShift->getClass()->getName());
        }
    } catch (exception& e) {
        try {
            rqs->push_back(NamedRiskQuantity::SP( 
                RiskQuantity::untweakable(e),
                IResultsIdentifier::SP(getPacketName(), 
                                       resultOutputName->toString()))); 
        } catch (exception& ) {

        }
    }
    return rqs;
}

void ImpliedScalarShiftMulti::validatePop2Object()
{
    static const string method("ImpliedScalarShiftMulti::validatePop2Object");
    
    SensitivityArrayConstSP sensArray(new SensitivityArray());
    OutputRequestArrayConstSP requestArray(new OutputRequestArray());
    
    if (!!targetControl){
        sensArray = targetControl->getSens();
        requestArray = targetControl->getOutputRequests();

        /** At most 1 sensitivity / request. **/
        if (sensArray->size() + requestArray->size() > 1){
            throw ModelException(method,
                                 "targetControl should include only one"
                                 " sensitivity / request; got "
                                 + Format::toString(sensArray->size())
                                 + " and "
                                 + Format::toString(requestArray->size())
                                 + ", respectively");
        }    

        /** This part is from ImpliedScalarShift: 
           If control has 1 sensitivity / request, 
           need to provide market name 
           AS: I have taken this out for the time being (possibly), as
           it does not allow for any sensitivities in the
           Instrument packet, for. example IND_VOL, which does not
           require an output name. Proper validation here is quite
           tricky and should probably be left to fail downstream
           (failure in extracting results) if (sensArray->size() +
           requestArray->size() != 0){ if
           (!targetControlOutputName){ throw
           ModelException(method, "targetControlOutputName should
           be provided together with targetControl"); } } **/
    }

    locCtrl = CControlSP(new Control(sensArray,
                                     requestArray,
                                     false,
                                     ""));
}

/** for reflection */
ImpliedScalarShiftMulti::ImpliedScalarShiftMulti():
    RiskQuantityFactorySensitivity(TYPE,NAME),
    lBound(0.),
    hBound(0.),
    isPositive(false),
    shiftAllNames(false),
    atSolutionIdentifier(atSolutionIdentifierDefault)
{ } 

ImpliedScalarShiftMulti::~ImpliedScalarShiftMulti() {}

void ImpliedScalarShiftMulti::load(CClassSP& clazz)
{
    clazz->setPublic(); //// make visible to EAS/spreadsheet
    REGISTER(ImpliedScalarShiftMulti, clazz); 
    SUPERCLASS(RiskQuantityFactorySensitivity);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(sensToShift, "Sensitivity whose shift size we seek to solve for "
                       "(Inputted shift size will be used as initial guess)");
    FIELD(targetValue, "Target Value");
    FIELD(targetControl, "Control that the target value refers to"
                         " (eg delta). Price, if omitted");
    FIELD_MAKE_OPTIONAL(targetControl);
    FIELD(targetControlOutputName, "Market name associated with target control");
    FIELD_MAKE_OPTIONAL(targetControlOutputName);
    FIELD(rootFinder, "Root finder (e.g. "
                      "ImpliedScalarShiftMulti::DefaultRootFinder)");
    FIELD(resultOutputName, "name under which the result is reported back");
    FIELD(lBound, "initial low bound for result");
    FIELD_MAKE_OPTIONAL(lBound);
    FIELD(hBound, "initial high bound for result");
    FIELD_MAKE_OPTIONAL(hBound);
    FIELD(locCtrl, "locCtrl");
    FIELD_MAKE_TRANSIENT(locCtrl);
    FIELD(isPositive, "has the domain of the result to be positive?");
    FIELD_MAKE_OPTIONAL(isPositive);
    FIELD(shiftAllNames, "Try to shift all supported objects, regardless of name.");
    FIELD_MAKE_OPTIONAL(shiftAllNames);
    FIELD(ctrlAtSolution, "Control for final calculation");
    FIELD_MAKE_OPTIONAL(ctrlAtSolution);
    FIELD(atSolutionIdentifier, "Identifier for the results at the solution point");
    FIELD_MAKE_OPTIONAL(atSolutionIdentifier);
}

IObject* ImpliedScalarShiftMulti::defaultConstructor() {
    return new ImpliedScalarShiftMulti();
}

CClassConstSP const ImpliedScalarShiftMulti::TYPE = 
    CClass::registerClassLoadMethod("ImpliedScalarShiftMulti", 
                                        typeid(ImpliedScalarShiftMulti), 
                                        load);

bool ImpliedScalarShiftMultiLinkIn() { 
    return ImpliedScalarShiftMulti::TYPE != NULL; 
}

DRLIB_END_NAMESPACE
