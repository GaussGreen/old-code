//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : BetaSkewPointwiseTweak.cpp
//
//   Description : Tweak of each relevant point inside a skew surface
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//      This class is now deprecated - Use BetaSkewMatrixTweak instead.
//      See BetaSkewMatrixTweak.hpp for more information
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BetaSkewPointwiseTweak.hpp"
#include "edginc/BetaSkewGridResult.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

BetaSkewPointwiseTweak::IShift::~IShift(){} // empty
BetaSkewPointwiseTweak::IRestorableShift::~IRestorableShift(){} // empty
BetaSkewPointwiseTweak::ISensitivePoints::~ISensitivePoints(){} // empty
BetaSkewPointwiseTweak::ISensitivePointsImt::~ISensitivePointsImt(){} // empty

/** Sens Control for vega pointwise */
const string BetaSkewPointwiseTweak::NAME = "BETA_SKEW_POINTWISE";
const double BetaSkewPointwiseTweak::DEFAULT_SHIFT = 0.001;

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double BetaSkewPointwiseTweak::divisor() const{
    static const string method = "BetaSkewPointwiseTweak::divisor";
    // add check for 0 shift size - somewhere 
    double div = 100.0 * getShiftSize();
    if (Maths::isZero(div))
    {
        throw ModelException(method, "Shift size is zero");
    }
    return div;
}

/** returns the interface identifying what an object has to do in order
    to support the tweak that this object represents */
CClassConstSP BetaSkewPointwiseTweak::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to support the tweak that this object represents which is also
    restorable */
CClassConstSP BetaSkewPointwiseTweak::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool BetaSkewPointwiseTweak::nameMatches(
    const OutputName& name,
    IObjectConstSP obj){
    // cast obj to BetaSkewPointwiseTweak::Shift and then invoke name method
    const IShift& skewObj = dynamic_cast<const IShift&>(*obj);
    return name.equals(skewObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void BetaSkewPointwiseTweak::appendName(OutputNameArray& namesList,
                                        IObjectConstSP   obj){
    // cast obj to BetaSkewPointwiseTweak::Shift and then invoke name method
    const IShift& skewObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(skewObj.sensName(this)));
    namesList.push_back(outputName);
}

bool BetaSkewPointwiseTweak::shift(IObjectSP obj) {
    // cast obj to BetaSkewPointwiseTweak::Shift and then invoke shift method
    IShift& skewObj =
        dynamic_cast<IShift&>(*obj);
    return skewObj.sensShift(this);
}

    
void BetaSkewPointwiseTweak::restore(IObjectSP obj) {
    // cast obj to BetaSkewPointwiseTweak::Shift and then invoke restore method
    IRestorableShift& skewObj =
        dynamic_cast<IRestorableShift&>(*obj);
    skewObj.sensRestore(this);
}

void BetaSkewPointwiseTweak::calculate(TweakGroup* tweakGroup,
                                       CResults*   results)
{
    try{
        // get list of names to calculate result for
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        }
    
        for (int idx = 0; idx < names->size(); idx++){
            // store the name of what we want to shift
            setMarketDataName((*names)[idx]);
            /* skip over where result has been calculated already */
            if (!results->exists(this)){
                try {
                    // now get grid points for which we should tweak
                    BetaSkewGridPointArrayConstSP gridPoints = 
                        getGridPoints(tweakGroup->getInstrument(),
                                      (*names)[idx]);

                    if (!gridPoints){
                        // no grid points found: store "not applicable"
                        results->storeGreek(IObjectSP(new NotApplicable()), this);
                    } else {
    
                        // create room for storing the results
                        BetaSkewGridResultArraySP tweaks(
                            new BetaSkewGridResultArray(gridPoints->size()));
    
                        // then loop over the grid points
                        for (int jdx = 0; jdx < gridPoints->size(); jdx++){
                            currentPointToTweak = *(*gridPoints)[jdx];
                            
                            // calculate sens 
                            double firstDeriv = calcOneSidedFirstDeriv(
                                tweakGroup, results);
                                
                            // store result in array
                            (*tweaks)[jdx] = BetaSkewGridResultSP(
                                new BetaSkewGridResult(currentPointToTweak, firstDeriv));
    
                        }
    
                        // and store it
                        results->storeGreek(tweaks, this);
                    }
                }
                catch (exception& e) {
                    results->storeGreek(IObjectSP(new Untweakable(e)), this);
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e,  "BetaSkewPointwiseTweak::calculate");
    }
}

/** Returns the grid points which are to be tweaked */
BetaSkewGridPointArrayConstSP BetaSkewPointwiseTweak::getGridPoints(
    CInstrument* instrument, 
    OutputNameConstSP outputName)
{
    static const string method("BetaSkewPointwiseTweak::getGridPoints");
    const IModel* model = getModel();
    const ISensitivePoints* skewModel =
        dynamic_cast<const ISensitivePoints*>(model);

    try{
        if (!skewModel) {
            throw ModelException(method, "Error: attempted to calculate sensitive grid"
                                 "points for a model which is not derived from "
                                 "BetaSkewPointwiseTweak::ISensitivePoints");
        }
        BetaSkewGridPointArrayConstSP gridPoints =
            skewModel->getSensitiveBetaSkewPoints(outputName, 
                                                  instrument, 
                                                  tweakAll, 
                                                  numberOfNeighbours);
        return gridPoints;
    } catch (exception& e){
        throw ModelException(e, method," For model "+
                             getModel()->getClass()->getName()+" with "
                             "instrument "+instrument->getClass()->getName());
    }
}

/** Returns the current grid point which is to be tweaked */        
BetaSkewGridPoint BetaSkewPointwiseTweak::getCurrentPointToTweak() const {
    return currentPointToTweak;
}

/** constructor with explicit shift size */
BetaSkewPointwiseTweak::BetaSkewPointwiseTweak(const double shiftSize): 
    ScalarShift(TYPE, NAME, shiftSize),
    tweakAll(DEFAULT_VALUE_FOR_TWEAKALL),
    numberOfNeighbours(DEFAULT_VALUE_FOR_NUM_NEIGHBOURS),
    currentPointToTweak(DateTime(), 0.0)
{
    validatePop2Object();
}


/** for reflection */
BetaSkewPointwiseTweak::BetaSkewPointwiseTweak():
    ScalarShift(TYPE, NAME),
    tweakAll(DEFAULT_VALUE_FOR_TWEAKALL),
    numberOfNeighbours(DEFAULT_VALUE_FOR_NUM_NEIGHBOURS),
    currentPointToTweak(DateTime(), 0.0)
{
    validatePop2Object();
}

BetaSkewPointwiseTweak::BetaSkewPointwiseTweak(CClassConstSP clazz, const string name):
    ScalarShift(clazz, name),
    tweakAll(DEFAULT_VALUE_FOR_TWEAKALL),
    numberOfNeighbours(DEFAULT_VALUE_FOR_NUM_NEIGHBOURS),
    currentPointToTweak(DateTime(), 0.0)
{
    validatePop2Object();
}

BetaSkewPointwiseTweak::BetaSkewPointwiseTweak(CClassConstSP clazz, 
                                               const string name, 
                                               const double shiftSize):
    ScalarShift(clazz, name, shiftSize),
    tweakAll(DEFAULT_VALUE_FOR_TWEAKALL),
    numberOfNeighbours(DEFAULT_VALUE_FOR_NUM_NEIGHBOURS),
    currentPointToTweak(DateTime(), 0.0)
{
    validatePop2Object();
}

BetaSkewPointwiseTweak::BetaSkewPointwiseTweak(CClassConstSP clazz, 
                                               const string name, 
                                               const double shiftSize,
                                               const bool tweakAll,
                                               const int numberOfNeighbours):
    ScalarShift(clazz, name, shiftSize),
    tweakAll(tweakAll),
    numberOfNeighbours(numberOfNeighbours),
    currentPointToTweak(DateTime(), 0.0)
{
    validatePop2Object();
}


/** Called immediately after object constructed */
void BetaSkewPointwiseTweak::validatePop2Object() {
    static const string method("BetaSkewPointwiseTweak::validatePop2Object");
    if (numberOfNeighbours < 0) {
        // Note that a value of 0 IS accepted
        throw ModelException(method, "The numberOfNeighbours must be >= 0");
    }
}

/** Invoked when class is 'loaded' */
void BetaSkewPointwiseTweak::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BetaSkewPointwiseTweak, clazz);
    SUPERCLASS(ScalarShift);
    EMPTY_SHELL_METHOD(defaultBetaSkewPointwiseTweak);
    
    // register how to build our sensitivity
    SensitivityFactory::addSens(
        BetaSkewPointwiseTweak::NAME, 
        new GenericSensitivityFactory<BetaSkewPointwiseTweak>(), 
        new BetaSkewPointwiseTweak(BetaSkewPointwiseTweak::DEFAULT_SHIFT),
        BetaSkewPointwiseTweak::IShift::TYPE);
    
    FIELD(tweakAll,
                 "YES = tweaks ALL points, NO = tweaks only relevant points");
    FIELD_MAKE_OPTIONAL(tweakAll);

    FIELD(numberOfNeighbours,
                 "number of neighbour strikes to tweak at each side of "
                 "the used strikes, when tweakAll is false");
    FIELD_MAKE_OPTIONAL(numberOfNeighbours);

    FIELD(currentPointToTweak, "Current point to tweak");
    FIELD_MAKE_TRANSIENT(currentPointToTweak);
}

IObject* BetaSkewPointwiseTweak::defaultBetaSkewPointwiseTweak(){
    return new BetaSkewPointwiseTweak();
}

CClassConstSP const BetaSkewPointwiseTweak::TYPE =
    CClass::registerClassLoadMethod(
        "BetaSkewPointwiseTweak",
         typeid(BetaSkewPointwiseTweak),
         BetaSkewPointwiseTweak::load);

CClassConstSP const BetaSkewPointwiseTweak::IShift::TYPE =
    CClass::registerInterfaceLoadMethod(
        "BetaSkewPointwiseTweak::IShift",
        typeid(BetaSkewPointwiseTweak::IShift),
        0);

static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(BetaSkewPointwiseTweak::IRestorableShift, clazz);
    EXTENDS(BetaSkewPointwiseTweak::IShift);
}
    
CClassConstSP const BetaSkewPointwiseTweak::IRestorableShift::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "BetaSkewPointwiseTweak::IRestorableShift",
        typeid(BetaSkewPointwiseTweak::IRestorableShift), 
        restorableShiftLoad);

DRLIB_END_NAMESPACE
