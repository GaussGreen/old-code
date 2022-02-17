//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CompositeInstrument.hpp
//
//   Description : Defines a composite instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 1 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CompositeInstrument.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DRWrapper.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/ScaleOutputs.hpp"
#include "edginc/XMLWriter.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE


//// note no copies are made of the inputs
CompositeInstrument::CompositeInstrument(ScenarioSP    scenario,
                                         ObjectArraySP model, 
                                         ObjectArraySP inst, 
                                         ObjectArraySP ctrl,
                                         DoubleArraySP multiplier,
                                         DoubleArraySP weight,
                                         CMarketDataSP market) :
    CObject(TYPE), scenario(scenario), model(model), inst(inst), 
    ctrl(ctrl), multiplier(multiplier), weight(weight), market(market) {
    // empty
}

// for reflection
CompositeInstrument::CompositeInstrument(): CObject(TYPE){
    // empty
}

/** Runs 'regression test' */
IObjectSP CompositeInstrument::runTest() const{
    // leave any "write to file" on - allows ability to pick out individual
    // instruments
    return IObjectSP(run());
}

// EdrAction 
IObjectSP CompositeInstrument::run(){
    const CompositeInstrument* ci = this;
    return IObjectSP(ci->run());
}

/** Calculates price and sensitivities for given instrument and
    context together with supplied market data + any scenario shifts.
    Results are returned in either a ResultSet object. ( >1 instrument) 
    or a Results object (1 instrument) */
IObject* CompositeInstrument::run() const{
    static const string method = "CompositeInstrument::run";
    try {
        int numInsts = inst->size();
        if (numInsts == 0){
            throw ModelException(method, "Zero instruments supplied");
        }
        if (numInsts != ctrl->size()  ||
            numInsts != model->size() ||
            numInsts != multiplier->size() ||
            numInsts != weight->size()) {
            throw ModelException(method, 
                                 "instrument, model, control, multiplier"
                                 " & weight arrays must be the same length");
        }
        // create local variable in case market is null
        CMarketDataSP   marketData(market);
        if (!marketData){
            // review in conjunction with avoiding writing all market data out
            marketData = CMarketDataSP(new MarketData());
        }

        CResultsArraySP results(new CResultsArray(numInsts));
        CControlArray  controls(numInsts);

        if (CommandLineParams::hasParameter(CommandLineParams::Wrapper)) {
            // create a de-wrappered version of input file
            fileConvert();
        }

        // check that we don't have too many scaleable instruments and that
        // the asset-price-source flag is consistent
        if ( numInsts > 1 ) {
            string assetPriceSource;
            for (int i = 0; i < numInsts; i++) {
                CInstrumentSP instSP(CInstrumentSP::dynamicCast(
                    DRWrapper::drWrapperToObject((*inst)[i], CInstrument::TYPE)));
                CControlSP ctrlSP(CControlSP::dynamicCast(
                    DRWrapper::drWrapperToObject((*ctrl)[i], CControl::TYPE)));

                if ( ctrlSP->scaleOutputs() && 
                     IScaleOutputs::TYPE->isInstance(instSP.get()) ) {
                        throw ModelException(method,
                            "The scaleResults flag must be false for "
                                             "composite components"
                            " in composites with more than one component.");
                }
                if (i == 0){
                    assetPriceSource = ctrlSP->getAssetPriceSource();
                } else if (!CString::equalsIgnoreCase(
                    assetPriceSource, ctrlSP->getAssetPriceSource())){
                    throw ModelException(method, "AssetPriceSource setting "
                                         "must be the same in all controls");
                }
            }
        }

        // loop over each instrument in the composite
        // run it and stash the results
        for (int i = 0; i < numInsts; i++) {
            IModelSP      modelSP(IModelSP::dynamicCast(
                DRWrapper::drWrapperToObject((*model)[i], IModel::TYPE)));
            CInstrumentSP instSP(CInstrumentSP::dynamicCast(
                DRWrapper::drWrapperToObject((*inst)[i], CInstrument::TYPE)));
            controls[i] = CControlSP(CControlSP::dynamicCast(
                DRWrapper::drWrapperToObject((*ctrl)[i], CControl::TYPE)));

            // see if anyone wants to add in a greek using its default
            // constructor or go mega
            if (CommandLineParams::hasParameter(CommandLineParams::Mega)) {
                controls[i] = // default everything possible
                    CControlSP(SensitivityFactory::megaControl()); 
            }
            if (CommandLineParams::hasParameter(
                CommandLineParams::AddSens)) {
                IObjectSP param = CommandLineParams::
                    getParameterArgs(CommandLineParams::AddSens);
                const string& s =
                    (dynamic_cast<CString&>(*param)).stringValue();
                SensitivitySP sens(SensitivityFactory::defaultSensitivity(s));
                controls[i]->addSensitivity(sens);
            }
            if (CommandLineParams::hasParameter(CommandLineParams::AddReq)) {
                IObjectSP param = CommandLineParams::
                    getParameterArgs(CommandLineParams::AddReq);
                const string& s = 
                    (dynamic_cast<CString&>(*param)).stringValue();
                OutputRequestSP req(new OutputRequest(s));
                controls[i]->addRequest(req);
            }
            
            // either run scenario or risk mgr directly
            (*results)[i] = modelSP->go(instSP, scenario,
                                        controls[i], marketData);
            // scale the result by multiplier
            double thisMultiplier = numInsts == 1? 
                ((*multiplier)[0] * (*weight)[0]): (*multiplier)[i];
            if (!Maths::equals(thisMultiplier, 1.0)){
                (*results)[i]->scale(controls[i], thisMultiplier, false);
            }
        }
            
        IObject* resultsSet;
        if (numInsts == 1){ // different format of results for n = 1
            resultsSet = (*results)[0].release();
        } else {
            // now we 'just' combine all the results together somehow
            ResultsSP combined = Results::combineResults(controls,
                                                         *weight, *results);
            resultsSet = new ResultsSet(combined.release(), results.release());
        }
        return resultsSet;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

}

// create a de-wrappered version of input file
void CompositeInstrument::fileConvert() const {
    static const string method = "CompositeInstrument::fileConvert";
    try {
        IObjectSP param = CommandLineParams::
            getParameterArgs(CommandLineParams::Wrapper);
        const string& s = 
            (dynamic_cast<CString&>(*param)).stringValue();

        ObjectArraySP m(new ObjectArray(0));
        ObjectArraySP i(new ObjectArray(0));
        ObjectArraySP c(new ObjectArray(0));

        for (int j = 0; j < inst->size(); j++) {
            m->push_back(DRWrapper::drWrapperToObject((*model)[j],
                                                        Model::TYPE));
            i->push_back(DRWrapper::drWrapperToObject((*inst)[j],
                                                        Instrument::TYPE));
            c->push_back(DRWrapper::drWrapperToObject((*ctrl)[j],
                                                        Control::TYPE));
        }

        CompositeInstrument ci(scenario,m,i,c,multiplier,weight,market);
        
        XMLWriter xml(s);
        ci.write("CONVERTED", &xml);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** addin function wrapper for CompositeInstrument*/
static IObjectSP addinComposite(CompositeInstrument* params){
    return IObjectSP(params->run());
}

class CompositeInstrumentHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CompositeInstrument, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultInterface);
        FIELD(scenario, "scenario");
        FIELD_MAKE_OPTIONAL(scenario);
        FIELD(model, "model(s)");
        FIELD(inst, "instrument(s)");
        FIELD(ctrl, "ctrl(s)");
        FIELD(multiplier, "multiplier");
        FIELD(weight, "weighting");
        FIELD(market, "market");
        FIELD_MAKE_OPTIONAL(market);
        // registration for addin function
        Addin::registerInstanceObjectMethod("COMPOSITE_INSTRUMENT",
                                            Addin::RISK,
                                            "Calculates price and greeks "
                                            "for a composite",
                                            CompositeInstrument::TYPE,
                                            true,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*)addinComposite);
        
    }

    static IObject* defaultInterface(){
        return new CompositeInstrument();
    }
};

CClassConstSP const CompositeInstrument::TYPE = 
CClass::registerClassLoadMethod(
    "CompositeInstrument", typeid(CompositeInstrument), 
    CompositeInstrumentHelper::load);

/** Tests aggregation of results */
class CompInstTest: public CObject{
public:
    static CClassConstSP const TYPE;

    // four addin parameters
    CControlSP      control;     // the control
    CResultsArraySP results;     // component results
    CDoubleArray    multipliers; // applied to each result
    CDoubleArray    weights;     // used when combining scaled results

    /** pop data dictionary to an object */
    static IObjectSP combine(CompInstTest* params){
        int numResults = params->results->size();
        if (numResults != params->multipliers.size() ||
            numResults != params->weights.size()){
            throw ModelException("CompInstTest::combine",
                                 "Different number of results,"
                                 " weights and multipliers");
        }
        CResultsArraySP scaledResults(new CResultsArray(numResults));
        CControlArray   ctrls(numResults);
        for (int i = 0; i < numResults; i++){
            ctrls[i] = params->control; // can just use references
            (*scaledResults)[i] = CResultsSP((*params->results)[i].clone());
            (*scaledResults)[i]->scale(ctrls[i], params->multipliers[i],
                                       false);
        }
        // now we 'just' combine all the results together somehow
        ResultsSP combined = Results::combineResults(ctrls, params->weights, 
                                                     *scaledResults);
        return IObjectSP(new ResultsSet(combined.release(),
                                        scaledResults.release()));
    }

    CompInstTest(CClassConstSP clazz): CObject(clazz){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CompInstTest, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCompInstTest);
        FIELD(control, "Control object used to generate outputs");
        FIELD(results, "Results for each instrument");
        FIELD(multipliers, "Scaling factor for each result");
        FIELD(weights, "Used to combine scaled results");
        Addin::registerClassObjectMethod("AGGREGATE_RESULTS",
                                         Addin::RISK,
                                         "Combines several sets of results",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)combine);
    }

    static IObject* defaultCompInstTest(){
        return new CompInstTest(TYPE);
    }
};

CClassConstSP const CompInstTest::TYPE = CClass::registerClassLoadMethod(
    "CompInstTest", typeid(CompInstTest), load);

/** Class to allow EAS to call combine method. Ought to make CompInstTest
    implement the ClientRunnable interface but it's got a poor name. This
    approach is less work */
class CompositeAggregation: public CompInstTest, 
                            virtual public ClientRunnable {
    static CClassConstSP const TYPE;

    // EdrAction 
    IObjectSP run() {
        return combine(this);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CompositeAggregation, clazz);
        SUPERCLASS(CompInstTest);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultCompositeAggregation);
    }
    static IObject* defaultCompositeAggregation(){
        return new CompositeAggregation();
    }
    CompositeAggregation(): CompInstTest(TYPE){}
};

CClassConstSP const CompositeAggregation::TYPE = 
CClass::registerClassLoadMethod("CompositeAggregation", 
                                typeid(CompositeAggregation), load);

    

DRLIB_END_NAMESPACE

