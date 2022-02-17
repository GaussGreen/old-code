//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRGridShift.cpp
//
//   Description : Pointwise type IR vega sensitivities
//
//   Author      : Andrew J Swain
//
//   Date        : 22 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IRGridShift.hpp"
#include "edginc/IRGridResult.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/TweakGroup.hpp"

DRLIB_BEGIN_NAMESPACE

IRGridShift::~IRGridShift(){}

/** Overridden to copy over gridPoint */
IObject* IRGridShift::clone() const{
    IRGridShift* myCopy = DYNAMIC_CAST(IRGridShift, CObject::clone());
    myCopy->gridPoint = gridPoint; // object is const so need to deep copy
    return myCopy;
}

void IRGridShift::calculate(TweakGroup*  tweakGroup,
                            CResults*    results) {
    try{
        // get list of names to calculate result for
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty() || 
            avoid(tweakGroup->getInstrument(), tweakGroup->getModel())) {
            results->storeNotApplicable(this);            
        }
    
        for (int idx = 0; idx < names->size(); idx++){
            // store the name of what we want to shift
            setMarketDataName((*names)[idx]);
            /* skip over where result has been calculated already */
            if (!results->exists(this)){
                try {
                    // now get grid points for which we should tweak
                    IRGridPointAbsArrayConstSP gridPoints = 
                        getGridPoints(tweakGroup->getInstrument(),
                                      (*names)[idx]);
                    if (gridPoints->empty()){
                        results->storeGreek(IObjectSP(new NotApplicable()),
                                            getPacketName(),
                                            (*names)[idx]);
                    } else {
                        // create room for storing the results
                        IRGridResultArraySP tweaks(new IRGridResultArray(
                                                       gridPoints->size()));
                        // then loop over the grid points
                        for (unsigned int jdx = 0; jdx < gridPoints->size(); 
                             jdx++){
                            gridPoint = (*gridPoints)[jdx];
                            // calculate sens 
                            double firstDeriv = 
                                calcOneSidedFirstDeriv(tweakGroup, results);
                            // store result in array
                            (*tweaks)[jdx] = IRGridResult(gridPoint,
                                                          firstDeriv);
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
        throw ModelException(e,  "IRGridShift::calculate");
    }
}

/** what point is being tweaked ? */
ExpiryConstSP IRGridShift::tenor() const {
    if (!gridPoint) {
       throw ModelException("IRGridShift::tenor", "grid point is null");
    }
    return gridPoint->getTenor();
}

ExpiryConstSP IRGridShift::expiry() const {
    if (!gridPoint) {
       throw ModelException("IRGridShift::expiry", "grid point is null");
    }
    return gridPoint->getExpiry();
}

IRGridShift::IRGridShift(const CClassConstSP& clazz,
                         const string&        outputName,
                         const double         shiftSize): 
    ScalarShift(clazz, outputName, shiftSize){}

/** for reflection */
IRGridShift::IRGridShift(const CClassConstSP& clazz,
                         const string&        sensName):
    ScalarShift(clazz, sensName){}

class IRGridShiftHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IRGridShift, clazz);
        SUPERCLASS(ScalarShift);
    }
};

CClassConstSP const IRGridShift::TYPE = CClass::registerClassLoadMethod(
    "IRGridShift", typeid(IRGridShift), IRGridShiftHelper::load);

DRLIB_END_NAMESPACE
