
#include "edginc/config.hpp"
#define QLIB_SCALARSHIFT_CPP
#include "edginc/ScalarShift.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/SensMgr.hpp"

DRLIB_BEGIN_NAMESPACE

ScalarShift::~ScalarShift(){}

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns false */
bool ScalarShift::discreteShift() const{
    return false;
}

/** Returns the scalar shift size */
double ScalarShift::getShiftSize() const{
    return shiftSize;
}

/** Sets the scalar shift size */
void ScalarShift::setShiftSize(double shiftSize){
    this->shiftSize = shiftSize;
}

/** implements a one sided scalar derivative for each instance of the
    market data which is sensitive to this SensControl */
void ScalarShift::calculate(TweakGroup*  tweakGroup,
                            CResults*    results){
    try {
        // get list of names to calculate result for. Remove blanks
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        }

        for (int idx = 0; idx < names->size(); idx++){
            // store what we want to shift
            setMarketDataName((*names)[idx]);
            /* skip over where result has been calculated already */
            if (!results->exists(this)){
                try {
                    // calculate sens
                    double firstDeriv = 
                        calcOneSidedFirstDeriv(tweakGroup, results);
                    // and store it
                    results->storeScalarGreek(firstDeriv, this);
                }
                catch (exception& e) {
                    results->storeGreek(IObjectSP(new Untweakable(e)), this);
                }
            }
        }
    } catch (exception& e){
        throw ModelException(&e,  "ScalarShift::calculate");
    }
}

class ScalarShift_AltWorld: public IHypothesis::AlternateWorld {

public:

    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(ScalarShift_AltWorld, clazz);
        SUPERCLASS(IHypothesis::AlternateWorld);
        // no EMPTY_SHELL_METHOD: can't copy sensMgr
    };

    SensMgrOpt sensMgr; // $unregistered
    ScalarShiftSP shift; // $unregistered
    double shiftSize; // $unregistered

    ScalarShift_AltWorld(ScalarShiftSP shift,
                        double shiftSize,
                        OutputNameConstSP overrideName,
                        IObjectSP tweakGroup):
        IHypothesis::AlternateWorld(TYPE, tweakGroup, 0., true),
        sensMgr(tweakGroup.get()),
        shift(shift),
        shiftSize(shiftSize)
    {
        double oldShiftSize = shift->getShiftSize();
        shift->setShiftSize(shiftSize);
        try {
            world = sensMgr.shift(shift.get(), overrideName);
        }
        catch (...) {
            shift->setShiftSize(oldShiftSize);
            throw;
        }
        shift->setShiftSize(oldShiftSize);
    }

    void _undo() {
        double oldShiftSize = shift->getShiftSize();
        shift->setShiftSize(shiftSize);
        try {
            sensMgr.restore();
        }
        catch (...) {
            shift->setShiftSize(oldShiftSize);
            throw;
        }
        shift->setShiftSize(oldShiftSize);
    }
};

CClassConstSP const ScalarShift_AltWorld::TYPE =
    CClass::registerClassLoadMethod(
        "ScalarShift_AltWorld", typeid(ScalarShift_AltWorld),
        load);

/** IScalarPerNameShift implementation, for ImpliedScalarShift  */
IHypothesis::AlternateWorldSP ScalarShift::appliedTo(
        OutputNameConstSP nameToShift,
        double shiftSize,
        IObjectSP tweakGroup) {
    try {
        return IHypothesis::AlternateWorldSP(
            new ScalarShift_AltWorld(ScalarShiftSP::attachToRef(this),
                                     shiftSize, nameToShift, tweakGroup));
    }
    catch (exception& e) {
        throw ModelException(e, "ScalarShift::appliedTo()");
    }
}

OutputNameArrayConstSP ScalarShift::allNames(const IObject* x) const {
    return SensMgrConst(x).allNames(const_cast<ScalarShift*>(this));
}

/** calculates delta/gamma like derivatives via 2 sided tweak */
void ScalarShift::calculateTwoSidedDeriv(
    const string&    secondDerivName, // eg GAMMA
    TweakGroup*      tweakGroup,
    CResults*        results)
{
    try{
        // set up sens mgr
        SensMgr   sensMgr(tweakGroup);

        // get list of names to calculate result for. Remove blanks
        OutputNameArrayConstSP names(this->names(tweakGroup));

        if (hasOverrideNames()){
            // manually cope with any explicit names that aren't present.
            // This should be a virtual method on SensControl/Sensitivity - 
            // also see comments in SensMgr::names(SensControl*, Results*)
            OutputNameArrayConstSP allNames(sensMgr.allNames(this));
            OutputNameArraySP extraNames(OutputName::difference(names,
                                                                allNames));
            for (int i = 0; i < extraNames->size(); i++){
                results->storeGreek(IObjectSP(new NotApplicable()),
                                    getPacketName(),
                                    (*extraNames)[i]);
                results->storeGreek(IObjectSP(new NotApplicable()),
                                    secondDerivName,
                                    (*extraNames)[i]);
            }
        }
            
        if (names->empty() && !results->packetExists(getPacketName())) {
            // if nothing to do AND not already stored something
            results->storeNotApplicable(this);
            results->storeNotApplicable(secondDerivName);
        }
    
        for (int idx = 0; idx < names->size(); idx++){
            // store what we want to shift
            setMarketDataName((*names)[idx]);
            /* skip over result has been calculated already */
            if (!results->exists(this)){
                // calculate sens
                try {
                    ScalarShift*  alteredScalarShift = NULL; 

                    SensControlSP alteredModelControl(
                        tweakGroup->getModel()->AlterControl(this));

                    SensControlSP alteredControl(
                        tweakGroup->getInstrument()->AlterControl(getModel(),
                                                                  this));

                    if (!alteredControl){
                        //instrument takes priority over control
                        alteredControl = alteredModelControl; 
                    }
                    if ( !(!alteredControl) ) {
                        alteredScalarShift = 
                            dynamic_cast<ScalarShift*>(alteredControl.get());
                        // set control/alorithm if not set
                        if (!alteredScalarShift->control){
                            alteredScalarShift->control = control;
                        }
                        if (!alteredScalarShift->algorithm){
                            alteredScalarShift->algorithm = algorithm;
                        }
                    }

                    if (alteredScalarShift)
                    {
                        alteredScalarShift->
                            calcSingleTwoSidedFirstDeriv(secondDerivName, 
                                                         tweakGroup, 
                                                         results);
                    } else {
                        calcSingleTwoSidedFirstDeriv(secondDerivName, 
                                                     tweakGroup, 
                                                     results);
                    }
                }
                catch (exception& e){
                    results->storeGreek(IObjectSP(new Untweakable(e)), this);
                    results->storeGreek(IObjectSP(new Untweakable(e)),
                                        secondDerivName, (*names)[idx]);
                }
            }
        }
    } catch (exception& e){
        throw ModelException(&e,  "ScalarShift::calculateTwoSidedDeriv");
    }
}

/** calculates a single first order two-sided scalar derivative */
void ScalarShift::calcSingleTwoSidedFirstDeriv(
    const string&  secondDerivName, // eg GAMMA
    TweakGroup*    tweakGroup,
    CResults*      results)
{
    try {
        pair<double, double> derivs = twoSidedDerivs(tweakGroup, results);
    
        results->storeScalarGreek(derivs.first, this);
        // see if 2nd order deriv exists already
        OutputNameConstSP name = getMarketDataName();
        if (!results->exists(secondDerivName, name)){
            results->storeScalarGreek(derivs.second, secondDerivName, name);
            // store divisor used (needed for cross derivatives)
            results->storeScalarGreek(getShiftSize(),
                                      getSensOutputName()+
                                      Results::SHIFT_SIZE_POSTFIX,
                                      name);
        }
    }
    catch (exception& e){
        throw ModelException(&e, "ScalarShift::calcTwoSidedFirstDeriv");
    }
}

pair<double, double> ScalarShift::twoSidedDerivs(
    TweakGroup*    tweakGroup,
    CResults*      results)
{
    try {
        double origPrice = results->retrievePrice();

        double shiftedUpPrice = shiftAndPrice(tweakGroup, origPrice);
        // calc divisor before we flip shift size
        double divisor = this->divisor();
        setShiftSize(-getShiftSize());
        double shiftedDownPrice;
        try {
            shiftedDownPrice = shiftAndPrice(tweakGroup, origPrice);
        }
        catch (exception&) {
            setShiftSize(-getShiftSize());
            throw;
        }
        setShiftSize(-getShiftSize());
        return make_pair((shiftedUpPrice - shiftedDownPrice) / (2.0 * divisor),
                         (shiftedUpPrice + shiftedDownPrice - 2.0 * origPrice) /
                             (divisor * divisor));
    }
    catch (exception& e) {
        throw ModelException(&e, "ScalarShift::twoSidedDerivs");
    }
}

/** Calculates cross derivative using the two supplied shifts. 'This'
    shift is used for storing results. Note that the second order
    derivatives must have been calculated already for each of shift1
    and shift2. Additionally, each one must have stored the shift size
    used in that calculation.

    None of the inputs may be null.

    Algorithm used:
 ****************************************************************
 * cross gamma calculator - gamma with respect to multiple stocks 
 *
 *    gamma(i,j) = f(+,+) - f(-,+) - f(+,-) + f(-,-)  = gamma(j,i)
 *               ---------------------------------
 *                     (2*S(i)*dS) * (2*S(j)*dS)
 *                                                              
 * [where f(+/-,+/-) = f( S(i) +/- S(i)*dS, S(j) +/- S(j)*dS )]
 *
 * This requires four reprices, but is equivalent to the following
 * that uses the already computed stock gammas and requires only
 * two additional pricings
 *
 * gamma(i,j)=f(+,+)+f(-,-)-2*f(0,0)-gamma(i)*(S(i)*dS)^2-gamma(j)*(S(j)*dS)^2
 *            -------------------------------------------------------------
 *                                 2 * (S(i)*dS) * (S(j)*dS)
 *
 * [where f(0,0) is the value of the option]
 ****************************************************************/

void ScalarShift::calculateCrossDerivative(
    Sensitivity*  shift,         // used for storing results
    ScalarShift*  shift1,        // how to make the first shift
    const string& secondDeriv1,  // what 'gamma' is called for shift1
    ScalarShift*  shift2,        // how to make the second shift
    const string& secondDeriv2,  // what 'gamma' is called for shift2
    TweakGroup*   tweakGroup,
    CResults*     results)
{
    static const string routine("ScalarShift::calculateCrossDerivative");
    ScalarShift* shifts[2] = {shift1, shift2};
    try{
        // get original price
        // double origPrice = shift->getSensPrice(results);
        double origPrice = results->retrievePrice();

        // get list of names needed to be shifted for shift1
        OutputNameArrayConstSP names1(shift1->names(tweakGroup));
        const OutputNameArray& nm1 = *names1;
        /* are the two sensitivities essentially the same? (eg true for
          cross gamma, false for fx cross gamma) */
        bool sameSens = shift1->shiftInterface() == shift2->shiftInterface();
        /* the second list of names is eg the list of ccy pairs for FX
           cross gamma otherwise eg it's just the same stock names */
        OutputNameArrayConstSP names2 = sameSens? 
            names1: shift2->names(tweakGroup);
        const OutputNameArray& nm2 = *names2;

        bool applicable = false;
#if 0
        if ( names1->empty() || names2->empty() ||
             ( sameSens && names1->size() < 2 )) {
            results->storeNotApplicable(shift);            
        }
#endif

        // start looping over the names
        for (int idx2 = 0; idx2 < nm2.size(); idx2++) {
            // store what we want to shift
            shift2->setMarketDataName(nm2[idx2]);
            /* skip over blank names */
            if (!nm2[idx2]->isEmpty()){
                /* get first gamma */
                double gammas[2];
                double gammaShifts[2];
                double validDelta[2];

                // get first gamma and the shift size used
                validDelta[0] = ( results->isValidScalarGreek(
                                            secondDeriv2,
                                            nm2[idx2]) &&
                                  results->isValidScalarGreek(
                                            shift2->getSensOutputName()+
                                            Results::SHIFT_SIZE_POSTFIX,
                                            nm2[idx2]));

                if ( validDelta[0] ) {
                    gammas[1] = results->retrieveScalarGreek(
                                        secondDeriv2,
                                        nm2[idx2]);
                
                    gammaShifts[1] = results->
                        retrieveScalarGreek(shift2->getSensOutputName()+
                                            Results::SHIFT_SIZE_POSTFIX,
                                            nm2[idx2]);
                } else {
                    gammas[1]      = 0.0;
                    gammaShifts[1] = 0.0;
                }


                for (int idx1 = sameSens? idx2 + 1: 0;
                     idx1 < nm1.size(); idx1++)  {
                    /* skip over blank names and names for the same sens that match */
                    if (!nm1[idx1]->isEmpty() && 
                        (!sameSens || !nm1[idx1]->equals(nm2[idx2].get()))){
                        /* get second gamma and the shift size used */

                        validDelta[1] = (results->isValidScalarGreek(
                                             secondDeriv1,
                                             nm1[idx1]) &&
                                         results->isValidScalarGreek(
                                             shift1->getSensOutputName()+
                                             Results::SHIFT_SIZE_POSTFIX,
                                             nm1[idx1]));
                        
                        if (validDelta[1]) {
                            gammas[0] = results->retrieveScalarGreek(
                                secondDeriv1, 
                                nm1[idx1]);
                            gammaShifts[0] = results-> retrieveScalarGreek(
                                shift1->getSensOutputName()+
                                Results::SHIFT_SIZE_POSTFIX,
                                nm1[idx1]);
                        } else {
                            gammas[0]      = 0.0;
                            gammaShifts[0] = 0.0;
                        }
                        
                        // store what we want to shift
                        shift1->setMarketDataName(nm1[idx1]);
                        // store result under both permutations
                        OutputNameSP outputNames[2] = {
                            OutputNameSP(new OutputName(nm1[idx1].get(), 
                                                        nm2[idx2].get())),
                            OutputNameSP(new OutputName(nm2[idx2].get(), 
                                                        nm1[idx1].get()))};
                        // calculate and store single cross gamma
                        if ( validDelta[0] && validDelta[1] ) {
                            try{
                                applicable = true;
                                singleCrossGamma(shift,
                                                 shifts,
                                                 origPrice,
                                                 gammas,
                                                 gammaShifts,
                                                 tweakGroup,
                                                 outputNames,
                                                 results);
                            } catch (exception& e){
                                string s("Failed for "+nm1[idx1]->toString()+
                                         " and "+nm2[idx2]->toString());
                                ModelException x(e, routine, s);
                                for (int i = 0; i < 2; i++){
                                    results->storeGreek(
                                        IObjectSP(new Untweakable(x)),
                                        shift->getSensOutputName(),
                                        outputNames[i]);
                                }
                            }
                        } else {
                            applicable = true;
                            int i;
                            for (i = 0; i < 2; i++){
                                string s("Failed for "+nm1[idx1]->toString()+
                                         " and "+nm2[idx2]->toString()+
                                         " because of untweakable deltas");
                                ModelException x(routine, s);
                                results->storeGreek(
                                    IObjectSP(new Untweakable(x)),
                                    shift->getSensOutputName(),
                                    outputNames[i]);
                            }
                        }
                    }
                }
            }
        }

        if (!applicable) {
            results->storeNotApplicable(shift);            
        }
    } catch (exception& e){
        throw ModelException(e, routine, "Failed whilst calculating "+
                             shift->getSensOutputName());
    }
}

void ScalarShift::calculateCrossDerivative(
    Sensitivity*  shift,         // used for storing results
    ScalarShift*  shift1,        // how to make the first shift
    const string& secondDeriv1,  // what 'gamma' is called for shift1
    ScalarShift*  shift2,        // how to make the second shift
    const string& secondDeriv2,  // what 'gamma' is called for shift2
    TweakGroup*   tweakGroup,
    CResults*     results,
    MapStringToName* namePair)
{
    static const string routine("ScalarShift::calculateCrossDerivative");
    ScalarShift* shifts[2] = {shift1, shift2};
    try {
        bool applicable = false;
        // get original price
        double origPrice = results->retrievePrice();
        
        // namePair tells us explicitly what we want to compute for
        MapStringToName::iterator iter;

        for (iter = namePair->begin(); iter != namePair->end(); ++iter) {
            OutputNameSP nm1(new OutputName(iter->first));
            OutputNameSP nm2 = iter->second;

            /* get first gamma */
            double gammas[2];
            double gammaShifts[2];
            double validDelta[2];

            // get first gamma and the shift size used
            validDelta[0] = (results->isValidScalarGreek(secondDeriv2, nm2) &&
                             results->isValidScalarGreek(
                                 shift2->getSensOutputName()+
                                 Results::SHIFT_SIZE_POSTFIX,
                                 nm2));

            if (validDelta[0]) {
                gammas[1] = results->retrieveScalarGreek(secondDeriv2, nm2);
                
                gammaShifts[1] = results->
                    retrieveScalarGreek(shift2->getSensOutputName()+
                                        Results::SHIFT_SIZE_POSTFIX,
                                        nm2);
            } else {
                gammas[1]      = 0.0;
                gammaShifts[1] = 0.0;
            }

            // get second gamma and the shift size used 
            validDelta[1] = (results->isValidScalarGreek(secondDeriv1, nm1) &&
                             results->isValidScalarGreek(
                                 shift1->getSensOutputName()+
                                 Results::SHIFT_SIZE_POSTFIX,
                                 nm1));
                        
            if (validDelta[1]) {
                gammas[0] = results->retrieveScalarGreek(secondDeriv1, nm1);
                gammaShifts[0] = results->retrieveScalarGreek(shift1->getSensOutputName()+
                                                              Results::SHIFT_SIZE_POSTFIX,
                                                              nm1);
            } else {
                gammas[0]      = 0.0;
                gammaShifts[0] = 0.0;
            }
                        
            // store what we want to shift
            shift1->setMarketDataName(nm1);
            shift2->setMarketDataName(nm2);
            // store result under both permutations
            OutputNameSP outputNames[2] = {
                OutputNameSP(new OutputName(nm1.get(), nm2.get())),
                OutputNameSP(new OutputName(nm2.get(), nm1.get()))};
            // calculate and store single cross gamma
            if (validDelta[0] && validDelta[1]) {
                try {
                    applicable = true;
                    singleCrossGamma(shift,
                                     shifts,
                                     origPrice,
                                     gammas,
                                     gammaShifts,
                                     tweakGroup,
                                     outputNames,
                                     results);
                } catch (exception& e){
                    string s("Failed for "+nm1->toString()+ " and "+nm2->toString());
                    ModelException x(e, routine, s);
                    for (int i = 0; i < 2; i++){
                        results->storeGreek(
                            IObjectSP(new Untweakable(x)),
                            shift->getSensOutputName(),
                            outputNames[i]);
                    }
                }
            } else {
                applicable = true;
                int i;
                for (i = 0; i < 2; i++){
                    string s("Failed for "+nm1->toString()+
                             " and "+nm2->toString()+
                             " because of untweakable deltas");
                    ModelException x(routine, s);
                    results->storeGreek(
                        IObjectSP(new Untweakable(x)),
                        shift->getSensOutputName(),
                        outputNames[i]);
                }
            }
        }

        if (!applicable) {
            results->storeNotApplicable(shift);            
        }
    } catch (exception& e){
        throw ModelException(e, routine, "Failed whilst calculating "+
                             shift->getSensOutputName());
    }
}

/** useful for checking internally computed x gamma against those generated by
    default using the external tweaker */
#define xDEBUG_INST_XGAMMA
    
/* calculates single cross gamma - no checking for NULLs. Assumes that this
   ScalarShift has the required market data name set in it for the output and
   that shift1 and shift2 have the required market data names needed for
   the tweaking */
void ScalarShift::singleCrossGamma(
    Sensitivity*   shift,             // used for storing results
    ScalarShift*   shifts[2],
    double         origPrice,         /* (I) original price */
    double         gammas[2],         /* (I) single gammas for the two */
    double         gammaShifts[2],    /* (I) shift size used for gammas */
    TweakGroup*    tweakGroup,
    OutputNameSP   outputNames[2],    // (I) result stored under both names
    CResults*      results)
{
    static const string routine("ScalarShift::singleCrossGamma");
    try{
#ifdef DEBUG_INST_XGAMMA
        double origXGamma;
        bool   haveOrig = false;
        if (results->exists(shift->getSensOutputName(), outputNames[0])){
            haveOrig = true;
            origXGamma = 
                results->retrieveScalarGreek(shift->getSensOutputName(), 
                                             outputNames[0]);
        }
#else
        if (!results->exists(shift->getSensOutputName(), outputNames[0])){
#endif
            double   shiftPrice[2];
            double   divisor[2];
            /* retrieve the shiftsizes used for the delta */
            int i; // MSVC broken
            for (i = 0; i < 2; i++) {
                // set the shifts to be negative to begin with
                shifts[i]->setShiftSize(gammaShifts[i]);
            }
            /* do down shift first */
            for (i = 0; i < 2; i++) {
                shiftPrice[i] = shifts[0]->shiftAndPrice(shifts[1], 
                                                         tweakGroup,
                                                         origPrice);
                /* reverse shift/get shift size of the two 'assets' */
                for (int j = 0; j < 2; j++) {
                    if (i == 1) {
                        /* NOTE - both shifts contain [positive] shiftUp */
                        divisor[j] = shifts[j]->divisor();
                        if (Maths::isZero(divisor[j])){
                            throw ModelException(routine, "absolute shift is "
                                                 "zero");
                        } 
                    } else {
                        /* do positive shift 2nd time round */
                        shifts[j]->setShiftSize(- shifts[j]->getShiftSize());
                    }
                }
            }
            /* calculate gamma */
            double gamma = (shiftPrice[1]+shiftPrice[0] -2.0 * origPrice - 
                            (gammas[0] * divisor[0] * divisor[0] + 
                             gammas[1] * divisor[1] * divisor[1])) /
               (2.0 * divisor[0] * divisor[1]);
#ifdef DEBUG_INST_XGAMMA
            if (haveOrig){
                //if (!Maths::equals(gamma, origXGamma)){
                if (!Maths::areEqualWithinTol(gamma, origXGamma, 1.0e-10)){
                    throw ModelException(
                        routine, "Difference between precomputed"
                        " x gamma ("+
                        Format::toString("%.16f", origXGamma) +") and "
                        "externally tweaked one"+
                        Format::toString("%.16f", gamma) +")");
                }
            } else
#endif
            {
                results->storeScalarGreek(gamma, shift->getSensOutputName(),
                                          outputNames[0]);
                if (outputNames[1].get()){
                    results->storeScalarGreek(gamma,
                                              shift->getSensOutputName(),
                                              outputNames[1]);
                }
            }
#ifndef DEBUG_INST_XGAMMA
        }
#endif
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}
    
ScalarShift::ScalarShift(const CClassConstSP& clazz,
                         const string&        outputName,
                         const double&        shiftSize):
    SensControlPerName(clazz, outputName), shiftSize(shiftSize){}

/** for reflection */
ScalarShift::ScalarShift(const CClassConstSP& clazz,
                         const string&        outputName):
    SensControlPerName(clazz, outputName), shiftSize(0.0){}

class ScalarShiftHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ScalarShift, clazz);
        SUPERCLASS(SensControlPerName);
        IMPLEMENTS(IScalarPerNameShift);
        FIELD(shiftSize, "How big to make the tweak");
    }
};

CClassConstSP const ScalarShift::TYPE = CClass::registerClassLoadMethod(
    "ScalarShift", typeid(ScalarShift), ScalarShiftHelper::load);

DEFINE_TEMPLATE_TYPE(ScalarShiftArray);

DRLIB_END_NAMESPACE
