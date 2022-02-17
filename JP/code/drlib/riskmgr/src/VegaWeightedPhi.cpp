//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaWeightedPhi.cpp
//
//   Description : Vega weighted phi sensitivity
//
//   Author      : Andrew McCleery
//
//   Date        : 14 Oct 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/VegaWeightedPhi.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/Control.hpp"
#include "edginc/PhiParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

const string VegaWeightedPhi::NAME = "VEGA_WEIGHTED_PHI";
const double VegaWeightedPhi::DEFAULT_SHIFT = PhiParallel_defaultShift;

/** constructor with explicit shift size */
VegaWeightedPhi::VegaWeightedPhi(double shiftSize):
    Sensitivity(TYPE), shiftSize(shiftSize){
}

/** for reflection */
VegaWeightedPhi::VegaWeightedPhi():
    Sensitivity(TYPE), shiftSize(DEFAULT_SHIFT){
}

/** identifies the name used for storing associated results in the output*/
const string& VegaWeightedPhi::getSensOutputName() const{
    return NAME;
}

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one (return value: false) */
bool VegaWeightedPhi::discreteShift() const{
    return false;
}

// Calculate vega-weighted phi (parallel phi weighted by vega parallel)
// (questionable theoretical justification)
void VegaWeightedPhi::calculate(TweakGroup*      tweakGroup,
                                Results*         results) {
    static const string method = "VegaWeightedPhi::calculate";
    try {

        // Kick off vega skew parallel if necessary
        VegaSkewParallelSP vegaSkewSens(new VegaSkewParallel(VegaSkewParallel::DEFAULT_SHIFT));

        // Don't report vega skew parallel if not requested.
        if (!control->sensitivityRequested(vegaSkewSens->getClass())){
            control->removePacketAfterCalc(vegaSkewSens->getPacketName());
        }

        // List of everything to  vega-skew-parallel-tweak
        OutputNameArraySP names(vegaSkewSens->names(tweakGroup).clone());

        // Kick off vega parallels if necessary
        VegaParallelSP vegaSens(new VegaParallel(VegaParallel::DEFAULT_SHIFT));
        
      
        // create a copy of result without any eventually calculated vega
        ResultsSP resultVega(copy(results));
        OutputNameSP vegaOutputName(new OutputName(
                                        vegaSens->getSensOutputName()));
        resultVega->removePacket(vegaSens->getPacketName());
        
        // if the set of names is empty then vega should return notApplicable
        // we use a shortcut here as if the names are empty, they are populated automatically
        // ResultsSP baseResults(tweakGroup->getModel()->Run(tweakGroup->getInstrument(), ctrl.get()));
        if (names->size()==0) {
            // nothing to do
            results->storeNotApplicable(getPacketName());
        }
        else {
            // only use the names in vega skew parallel for vega parallel
            vegaSens->storeOverrideNames(names);
            
            vegaSens->calculateSens(tweakGroup->getModel(),
                                    tweakGroup->getInstrument(),
                                    control,
                                    resultVega.get());
        
            // Kick off phi parallel if necessary
            SensitivitySP phiSens(PhiParallel_SP(shiftSize, true));
            
            // Change to calculateOneSidedFirstDeriv() once derive from ScalarShift rather than Sensivitivity
            phiSens->calculateSens(tweakGroup->getModel(),
                                   tweakGroup->getInstrument(),
                                   control,
                                   results);
            
            OutputNameSP phiOutputName(new OutputName(phiSens->getSensOutputName()));
            bool isNotApplicable = false;
            IObjectConstSP phiResultObj = results->retrieveGreek(phiSens->getPacketName(),
                                                                 phiOutputName);
            
            // Now remove PHI_PARALLEL from results set if not requested
            if (!control->sensitivityRequested(phiSens->getClass())){
                results->removeGreek(phiSens->getPacketName(), phiOutputName);
            }
            
            int i=0, j=0, N=0;
            double phi=0.0;
            DoubleArraySP vega(new DoubleArray(0));
            vector<pair<OutputNameConstSP, IObjectConstSP> > vegaResultList;
            
            if (NotApplicable::TYPE->isInstance(phiResultObj)) {
                isNotApplicable = true;
                
            } else if (Untweakable::TYPE->isInstance(phiResultObj)) {
                throw ModelException(method, phiSens->getSensOutputName() +
                                     " sensitivity is Untweakable");
                
            } else {
                
                phi = (CDoubleConstSP::dynamicCast(phiResultObj))->doubleValue();
                
                // Extract vega parallels from results set
                vegaResultList = resultVega->listPacketResults(vegaSens->getPacketName()/*VegaParallel::NAME*/);
                
                N = vegaResultList.size();      // Number of underlyings
                vega->resize(N);
                
                for (j = 0; j < N; j++) {
                    IObjectConstSP vegaResultObj = resultVega->retrieveGreek(vegaSens->getPacketName(),
                                                                             vegaResultList[j].first);
                    if (NotApplicable::TYPE->isInstance(vegaResultObj)) {
                        isNotApplicable = true;
                    } else if (Untweakable::TYPE->isInstance(vegaResultObj)) {
                        throw ModelException(method, "At least one " +
                                             vegaSens->getSensOutputName() +
                                             " sensitivity is Untweakable");
                    } else {
                        // Copy the double into the vega array
                        double tmp = (CDoubleConstSP::dynamicCast(vegaResultObj))->doubleValue();
                        (*vega)[j] = tmp;
                    }
                }
            }
            
            if (isNotApplicable || N==1) {
                
                results->storeNotApplicable(getPacketName());
                
            } else {
                
                double vegaSum = 0.0;
                for (i = 0; i < N; i++) {
                    for (j = i + 1; j < N; j++) {
                        vegaSum += fabs((*vega)[i] * (*vega)[j]);
                    }
                }
                
                if (Maths::isZero(vegaSum)) {
                    throw ModelException(method, "Sum of products of vega parallel "
                                         "sensitivities is zero.");
                }
                
                // Report N(N-1)/2 numbers given N underlyings
                double vegaWeightedPhi=0.0;
                for (i = 0; i < N; i++) {
                    for (j = i + 1; j < N; j++) {
                        OutputNameConstSP outputName(new OutputName(vegaResultList[i].first->toString(),
                                                                    vegaResultList[j].first->toString()));
                        
                        vegaWeightedPhi = (phi * fabs((*vega)[i] * (*vega)[j])) / vegaSum;
                        
                        results->storeScalarGreek(vegaWeightedPhi,
                                                  getPacketName(),
                                                  outputName);
                    }
                }
            }
        }
    } catch (exception& e) {
      results->storeGreek(IObjectSP(new Untweakable(e)),
                          getSensOutputName(), OutputNameSP(new OutputName("")));
    }
}

class VegaWeightedPhiHelper{
   /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VegaWeightedPhi(VegaWeightedPhi::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VegaWeightedPhi(shiftSize);
        }
    };

public:
        /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VegaWeightedPhi, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVegaWeightedPhi);
        FIELD(shiftSize, "How big to make the tweak");

        // register how to build our sensitivity
        SensitivityFactory::addSens(VegaWeightedPhi::NAME,
                                    new Factory(),
                                    new VegaWeightedPhi(VegaWeightedPhi::DEFAULT_SHIFT),
                                    0);
    }

    static IObject* defaultVegaWeightedPhi(){
        return new VegaWeightedPhi();
    }
};

CClassConstSP const VegaWeightedPhi::TYPE = CClass::registerClassLoadMethod(
    "VegaWeightedPhi", typeid(VegaWeightedPhi), VegaWeightedPhiHelper::load);

DRLIB_END_NAMESPACE
