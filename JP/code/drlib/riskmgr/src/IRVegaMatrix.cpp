//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives
//
//   Filename    : IRVegaMatrix.cpp
//
//   Description : IR vega pointwise sensitivity
//                 This is essentially the same of IRVegaPointwise but
//                 changed to store the array of results inside an object
//                 rather than storing the actual array.
//                 IRVegaPointwise is now deprecated.
//
//   Author      : Jose Hilera
//
//   Cautions    : The sensitivity implemented here is derived from 
//                 IRVegaPointwise, which is to be removed. To avoid duplication
//                 IRVegaPointwise::IShift and IRVegaPointwise::IRestorableShift
//                 are still used (and, indirectly, also 
//                 IRVegaPointwise::ISensitivePoints) - they should be moved into
//                 this class when IRVegaPointwise is finally removed, along the 
//                 methods in IRVegaPointwise which this class is currently 
//                 inheriting (NOT calculate, which VegaProxyMatrix overrides!)
//
//   Date        : July 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Model.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/IRGrid.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/IRVegaMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for vega matrix */
const string IRVegaMatrix::NAME = "IRVEGA_MATRIX";
const double IRVegaMatrix::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
IRVegaMatrix::IRVegaMatrix(double shiftSize):
    IRVegaPointwise(TYPE, NAME, shiftSize) {
}

/** for reflection */
IRVegaMatrix::IRVegaMatrix():
    IRVegaPointwise(TYPE, NAME) {
}

void IRVegaMatrix::calculate(TweakGroup*  tweakGroup,
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
                        IRGridSP tweak(new IRGrid(gridPoints->size()));
                        
                        // then loop over the grid points
                        for (unsigned int jdx = 0; jdx < gridPoints->size(); 
                             jdx++){
                            gridPoint = (*gridPoints)[jdx];
                            // calculate sens 
                            double firstDeriv = 
                                calcOneSidedFirstDeriv(tweakGroup, results);
                            // store result in array
                            tweak->setPoint(jdx, IRGridResult(gridPoint, firstDeriv));
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
        throw ModelException(e,  "IRGridShift::calculate");
    }
}

 
class IRVegaMatrixHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new IRVegaMatrix(IRVegaMatrix::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new IRVegaMatrix(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IRVegaMatrix, clazz);
        SUPERCLASS(IRVegaPointwise);
        EMPTY_SHELL_METHOD(defaultIRVegaMatrix);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(IRVegaMatrix::NAME, 
                                    new Factory(), 
                                    new IRVegaMatrix(IRVegaMatrix::DEFAULT_SHIFT),
                                    IRVegaPointwise::IShift::TYPE);
    }

    static IObject* defaultIRVegaMatrix(){
        return new IRVegaMatrix();
    }
};

CClassConstSP const IRVegaMatrix::TYPE = CClass::registerClassLoadMethod(
    "IRVegaMatrix", typeid(IRVegaMatrix), IRVegaMatrixHelper::load);


bool IRVegaMatrixLinkIn() {
    return IRVegaMatrix::TYPE != NULL;
}

DRLIB_END_NAMESPACE
