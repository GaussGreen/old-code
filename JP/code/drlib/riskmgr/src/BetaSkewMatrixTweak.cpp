//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BetaSkewMatrixTweak.cpp
//
//   Description : Tweak of each relevant point inside a skew surface
//                 This is essentially a copy of BetaSkewPointwiseTweak but
//                 changed to store the array of results inside an object
//                 rather than storing the actual array.
//                 BetaSkewPointwiseTweak is now deprecated.
//
//   Author      : Jose Hilera
//
//   Date        : July 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/BetaSkewGrid.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/BetaSkewMatrixTweak.hpp"


DRLIB_BEGIN_NAMESPACE

/** constructor with explicit shift size */
BetaSkewMatrixTweak::BetaSkewMatrixTweak(const double shiftSize): 
    BetaSkewPointwiseTweak(TYPE, NAME, shiftSize)
{}

/** constructor with explicit shift size and tweakAll*/
BetaSkewMatrixTweak::BetaSkewMatrixTweak(const double shiftSize,
                                         const bool   tweakAll,
                                         const int    numberOfNeighbours): 
    BetaSkewPointwiseTweak(TYPE, NAME, shiftSize, tweakAll, numberOfNeighbours)
{}

/** for reflection */
BetaSkewMatrixTweak::BetaSkewMatrixTweak():
    BetaSkewPointwiseTweak(TYPE, NAME) {}

/** Sens Control for vega matrix */
const string BetaSkewMatrixTweak::NAME = "BETA_SKEW_MATRIX";
const double BetaSkewMatrixTweak::DEFAULT_SHIFT = 0.001;

void BetaSkewMatrixTweak::calculate(
    TweakGroup* tweakGroup,
    CResults* results)
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
                        BetaSkewGridSP tweak(
                            new BetaSkewGrid(gridPoints->size()));

                        // then loop over the grid points
                        for (int jdx = 0; jdx < gridPoints->size(); jdx++){
                            currentPointToTweak = *(*gridPoints)[jdx];
                            
                            // calculate sens 
                            double firstDeriv = calcOneSidedFirstDeriv(
                                tweakGroup, results);
                                
                            // store result in array
                            tweak->setPoint(jdx, 
                                            BetaSkewGridResultSP(
                                               new BetaSkewGridResult(currentPointToTweak, 
                                                                      firstDeriv)));
                        }
    
                        // and store it
                        results->storeGreek(tweak, this);
                    }
                }
                catch (exception& e) {
                    results->storeGreek(IObjectSP(new Untweakable(e)), this);
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e,  "BetaSkewMatrixTweak::calculate");
    }
}


/** Invoked when class is 'loaded' */
void BetaSkewMatrixTweak::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BetaSkewMatrixTweak, clazz);
    SUPERCLASS(BetaSkewPointwiseTweak);
    IMPLEMENTS(Additive);
    EMPTY_SHELL_METHOD(defaultBetaSkewMatrixTweak);
    
    // register how to build our sensitivity
    SensitivityFactory::addSens(
        BetaSkewMatrixTweak::NAME, 
        new GenericSensitivityFactory<BetaSkewMatrixTweak>(), 
        new BetaSkewMatrixTweak(BetaSkewMatrixTweak::DEFAULT_SHIFT),
        BetaSkewPointwiseTweak::IShift::TYPE);
}

IObject* BetaSkewMatrixTweak::defaultBetaSkewMatrixTweak(){
    return new BetaSkewMatrixTweak();
}

CClassConstSP const BetaSkewMatrixTweak::TYPE =
    CClass::registerClassLoadMethod(
        "BetaSkewMatrixTweak",
         typeid(BetaSkewMatrixTweak),
         BetaSkewMatrixTweak::load);

bool BetaSkewMatrixTweakLinkIn() {
    return BetaSkewMatrixTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE
