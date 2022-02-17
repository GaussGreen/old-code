//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MegaShuffleInterface.cpp
//
//   Description : Decorates RiskMgr to enable a shuffled mega test.
//
//   Author      : Jon Dee
//
//   Date        : 22 November 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RiskMgrInterface.hpp"
#include "edginc/MegaShuffleInterface.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Results.hpp"
#include "edginc/SpotPrice.hpp"
#include "edginc/SpotPriceProxy.hpp"

DRLIB_BEGIN_NAMESPACE

MegaShuffleInterface::MegaShuffleInterface(IObjectSP _riskMgr, int _seed) : CObject(TYPE), seed(_seed)
{
    riskMgr = RiskMgrInterfaceSP::dynamicCast(_riskMgr);
    QLIB_VERIFY(!!riskMgr, "Object passed into MegaShuffleInterface is not a valid riskMgr object.");
}

MegaShuffleInterface::MegaShuffleInterface() : CObject(TYPE)
{
}

/** Invoked when Class is 'loaded' */
void MegaShuffleInterface::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MegaShuffleInterface, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ClientRunnable);
    EMPTY_SHELL_METHOD(defaultInterface);
    FIELD(riskMgr, "riskMgr");
    FIELD(seed, "seed");
}

IObject* MegaShuffleInterface::defaultInterface(){
    return new MegaShuffleInterface();
}

CClassConstSP const MegaShuffleInterface::TYPE = CClass::registerClassLoadMethod(
    "MegaShuffleInterface", typeid(MegaShuffleInterface), MegaShuffleInterface::load);

IObjectSP MegaShuffleInterface::runTest() const { 
    // don't want to write out an input file
    riskMgr->ctrl->switchOffWriteToFile();
    return go();
}

IObjectSP MegaShuffleInterface::run() { 
    return go();
}


ControlSP MegaShuffleInterface::makeControl(SensitivitySP      sensToUse) const {
    static const string method = "makeControl";
    try {
        CStringArraySP       requests = OutputRequest::allRequests();
        OutputRequestArraySP req(new OutputRequestArray(0));
        for (int i = 0; i < requests->size(); i++) {
            OutputRequestSP request = OutputRequestSP(new OutputRequest((*requests)[i]));
            req->push_back(request);
        }

        SensitivityArraySP sens(new SensitivityArray(0));
       
        CControlSP ctrl(new Control(sens,
            req,
            false,
            "",
            Control::ASSET_THEO));

        ctrl->addSensitivity(sensToUse);
        return ctrl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// 'CALC_TIME' will change every run, so remove it from the list of
// output requests (as well as other similar requests).
void MegaShuffleInterface::removeUnpredictableOutput(CControlSP control) const { 
    control->removeRequest(OutputRequest::CALC_TIME);
    control->removeRequest(OutputRequest::PRICE_TIME);
    control->removeRequest(OutputRequest::COMPUTE_INDEX);
    control->removeRequest(OutputRequest::COMPUTE_ESTIMATE);

    control->removeSensitivity(SpotPrice::TYPE);
    control->removeSensitivity(SpotPriceProxy::TYPE);
}

const char *getStuff(IObjectSP obj) { 
    return obj->getClass()->getName().c_str();
}

IObjectSP MegaShuffleInterface::go() const { 

    RiskMgrInterfaceSP riskOrig(riskMgr.clone());

    // A 'Mega Control' contains (almost) all sensitivities and
    // output requests.
    CControlSP ctrlMega(SensitivityFactory::megaControl(true));
    ctrlMega->shuffleSens(seed);

    removeUnpredictableOutput(ctrlMega);
    SensitivityArrayConstSP megaSens = ctrlMega->getSens();

    riskMgr->ctrl.reset(ctrlMega.get());
    CResultsSP shuffledResults = CResultsSP::dynamicCast(riskMgr->run());

    // run each sensitivity in isolation, 
    // to spot any undesirable interactions between sensitivities
    // (e.g. those caused by caching)
    CResultsSP incrementalResults;
    int nSize = megaSens->size();
    for (int i = 0;i<nSize;i++) {
        SensitivitySP sens(
            SensitivityConstSP::dynamicCast(megaSens->get(i)).clone());        
        RiskMgrInterfaceSP riskSingle(riskOrig.clone());    
        riskSingle->ctrl = makeControl(sens);
        removeUnpredictableOutput(riskSingle->ctrl);
        CResultsSP tempResults = CResultsSP::dynamicCast(riskSingle->run());
        if (i == 0) { 
            incrementalResults.reset(tempResults.clone());    
        }
        else {
            //Some sensitivities may be associated with more than one packet.
            //For example DeltaToCredit adds info to both DELTA_TO_CREDIT
            //and GAMMA_TO_CREDIT.
            incrementalResults->merge(tempResults.get());   
        }
    }

    // The control is used to get the list of sensitivities 
    // to add results for. We want all sensitivities, so use the mega.
    CControlArray controls;
    controls.append(ctrlMega);
    controls.append(ctrlMega);

    DoubleArray weights;
    weights.append(CDoubleSP(CDouble::create(1.0)));
    weights.append(CDoubleSP(CDouble::create(-1.0)));

    CResultsArray resultsArray;
    resultsArray.append(shuffledResults);
    resultsArray.append(incrementalResults);

    // Add one result with weight 1.0, the other with -1.0.
    // We should then get zeros in every (applicable) packet, 
    // if the results are the same, as expected.
    return CResults::combineResults(controls, weights, resultsArray);
}

int MegaShuffleInterface::createShuffleSeed()
{
    const int day = 60 * 60 * 24; 
    time_t seed;
    time(&seed);
    seed = seed / day;
    return (int)seed;
}
DRLIB_END_NAMESPACE

