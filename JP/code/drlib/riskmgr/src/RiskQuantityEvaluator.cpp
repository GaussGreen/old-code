/**
 * @file RiskQuantityEvaluator.cpp
 */

// Note: code marked /**/ is temporary, for compatibility.  See notes on
// RiskQuantityEvaluator::storeResults().

#include "edginc/config.hpp"
#include ext_hash_set
#include "edginc/TRACE.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/NotApplicableException.hpp"
#include "edginc/HypotheticalQuantity.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/LazyRiskQuantityFactory.hpp"
#include "edginc/RiskMapping.hpp"
#include "edginc/RiskMappingMatrix.hpp"
#include "edginc/IRiskQuantityFactory.hpp"
#include "edginc/RiskQuantityEvaluator.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"

#include "edginc/ICompatibilitySensitivity.hpp" /**/
#include "edginc/AbstractPropertyTweakHypothesis.hpp"       /**/
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE



// 
// ******************************
//  RiskQuantityEvaluator::Value
// ******************************
// 

RiskQuantityEvaluator::Value::Value(NamedRiskQuantitySP riskQuantity,
                                    int instrument):
    CObject(TYPE),
    _value(0.),
    _valueValid(false),
    riskQuantity(riskQuantity),
    instrument(instrument),
    vals(new DoubleArray(riskQuantity->riskQuantity->parameters()->size())),
    dsts(new DoubleArray(riskQuantity->riskQuantity->parameters()->size())),
    hypothesesApplicable(true),
    instrumentApplicable(true),
    usedRiskMapping(false)
{}

RiskQuantityEvaluator::ValueSP RiskQuantityEvaluator::Value::SP(
        NamedRiskQuantitySP riskQuantity, int instrument) {
    return ValueSP(new Value(riskQuantity, instrument));
}

RiskQuantityEvaluator::Value::~Value() {}

double RiskQuantityEvaluator::Value::value() {
    if (!!oops) {
        throw ModelException("RiskQuantityEvaluator::Value::value()",
                             oops->getMessage());
    }
    else if (!(hypothesesApplicable && instrumentApplicable)) {
        throw NotApplicableException();
    }
    else if (_valueValid) {
        return _value;
    }
    else {
        try {
            _value = riskQuantity->riskQuantity->value(*vals, *dsts);
        }
        catch (NotApplicableException&) {
            hypothesesApplicable = false;
            throw;
        }
        catch (exception& e) {
            oops.reset(new Untweakable(e));
            throw ModelException("RiskQuantityEvaluator::Value::value()",
                                 oops->getMessage());
        }

        _valueValid = true;

        return _value;
    }
}

static bool maybeOK(const RiskQuantityEvaluator::Value& v) {
    return v.hypothesesApplicable && v.instrumentApplicable && !v.oops;
}

bool RiskQuantityEvaluator::Value::ok() {
    if (!hypothesesApplicable || !instrumentApplicable || !!oops) return false;

    try {
        value();
        return true;
    }
    catch (exception&) {
        return false;
    }
}

void RiskQuantityEvaluator::Value::storeResult(CResultsArraySP rs) {
    if (!!oops) {
        riskQuantity->resultsName->storeUntweakable((*rs)[instrument], oops);
    }
    else if (!instrumentApplicable) {
        riskQuantity->resultsName->storeNotApplicableToInstrument(
            (*rs)[instrument]);
    }
    else if (!hypothesesApplicable) {
        riskQuantity->resultsName->storeNotApplicableToName((*rs)[instrument]);
    }
    else {
        riskQuantity->storeResult(*vals, *dsts, (*rs)[instrument]);
    }
}

static IObject* defaultValue() {
    return new RiskQuantityEvaluator::Value(NamedRiskQuantitySP(), 0);
}

void RiskQuantityEvaluator::Value::load(CClassSP& clazz) {
    REGISTER(RiskQuantityEvaluator::Value, clazz);
    EMPTY_SHELL_METHOD(defaultValue);
    SUPERCLASS(CObject);
    FIELD(riskQuantity, "riskQuantity");
    FIELD(instrument, "instrument");
    FIELD(vals, "vals");
    FIELD(dsts, "dsts");
    FIELD(oops, "oops");
    FIELD(instrumentApplicable, "instrumentApplicable");
    FIELD(hypothesesApplicable, "hypothesesApplicable");
    FIELD(usedRiskMapping, "usedRiskMapping");
}

typedef RiskQuantityEvaluator::Value RiskQuantityEvaluator_Value;

CClassConstSP const RiskQuantityEvaluator_Value::TYPE = CClass::registerClassLoadMethod(
    "RiskQuantityEvaluator::Value", typeid(RiskQuantityEvaluator_Value), load);

typedef RiskQuantityEvaluator::ValueArray RiskQuantityEvaluator_ValueArray;

DEFINE_TEMPLATE_TYPE_WITH_NAME("RiskQuantityEvaluator::ValueArray", RiskQuantityEvaluator_ValueArray);



// 
// ****************
//  HypothesisTree
// ****************
// 

// util func used in TRACE output
static string historyString(AtomicHypothesisArrayConstSP history) {
    if (history->empty()) return "no changes, i.e. base case";
    string it;
    for (int h = 0; h < history->size(); ++h) {
        if (h) it += " & ";
        it += (*history)[h]->toString();
    }
    return it;
}

FORWARD_DECLARE(HypothesisTree)

/**
 * Tree of all the states we need to put the world into, by applying various
 * combinations of IHypothesis (=~ tweak), and the HypotheticalQuantity's we
 * want to evaluate in those alternate worlds, in order to compute some
 * greeks.
 *
 * The only "entry point" is evaluate().
 */

class HypothesisTree: public CObject {

public:

    static CClassConstSP const TYPE;

private:

    static void load(CClassSP& clazz) {
        REGISTER(HypothesisTree, clazz);
        SUPERCLASS(CObject);
    }

    /**
     * Where to store the values of the HypotheticalQuantity's as we compute
     * them
     */

    RiskQuantityEvaluator::ValueArraySP values; // $unregistered

    /**
     * Quantities we need to compute at this node to evaluate the required
     * HypotheticalQuantity's, and where to put their values within 'values'
     */

    //@{

    struct QuantityToCompute {
        IResultsFunctionConstSP quantity;
        int instrument;
        vector<pair<RiskQuantityEvaluator::Value *, int> > buffers;
        IHypothesis::IDistanceMetricConstSP distanceMetric;

        QuantityToCompute(IResultsFunctionConstSP quantity, int instrument,
                          IHypothesis::IDistanceMetricConstSP distanceMetric):
            quantity(quantity),
            instrument(instrument),
            distanceMetric(distanceMetric)
        {}

        void append(RiskQuantityEvaluator::Value* buffer, int paramIndex) {
            buffers.push_back(make_pair(buffer, paramIndex));
        }

        void compute(CResultsArrayConstSP results,
                     CDoubleArrayConstSP distances, bool found) const {
            if (buffers.size() > 0) {
                TRACE_METHOD;

                double it;
                try {
                    it = (*quantity)((*results)[instrument].get());
                }
                catch (NotApplicableException&) {
                    for (int b = 0; b < int(buffers.size()); ++b) {
                        buffers[b].first->instrumentApplicable = false;
                    }
                    return;
                }
                catch (ModelException& e) {
                    noteOops(UntweakableConstSP(new Untweakable(e)));
                    return;
                }

                double distance = (*distanceMetric)(*distances);

                for (int b = 0; b < int(buffers.size()); ++b) {
                    TRACE("Store that number, for future calculation of " <<
                          "instrument #" << instrument << "'s " <<
                          *buffers[b].first->riskQuantity->resultsName);

                    (*buffers[b].first->vals)[buffers[b].second] = it;
                    (*buffers[b].first->dsts)[buffers[b].second] = distance;
                    if (!found) {
                        buffers[b].first->hypothesesApplicable = false;
                    }
                }
            }
        }

        void noteOops(UntweakableConstSP oops) const {
            TRACE_METHOD;
            for (int b = 0; b < int(buffers.size()); ++b) {
                TRACE("Noting that " <<
                      "instrument #" << instrument << "'s " <<
                      *buffers[b].first->riskQuantity->resultsName <<
                      " will be Untweakable");
                if (!buffers[b].first->oops) buffers[b].first->oops = oops;
            }
        }

        bool maybeRelevant() const {
            for (int b = 0; b < int(buffers.size()); ++b) {
                if (maybeOK(*buffers[b].first)) {
                    return true;
                }
            }

            TRACE("We don't need to compute " << *quantity << " since "
                  "all of the values it's needed for are already Untweakable "
                  "or NotApplicable");

            return false;
        }

        static bool anyMaybeRelevant(const vector<QuantityToCompute*>& qs) {
            for (size_t q = 0; q < qs.size(); ++q)
                if (qs[q]->maybeRelevant()) return true;

            return false;
        }

        // duplication, groan ...

        static bool anyMaybeRelevant(const vector<QuantityToCompute>& qs) {
            for (size_t q = 0; q < qs.size(); ++q)
                if (qs[q].maybeRelevant()) return true;

            return false;
        }
    };

    vector<QuantityToCompute> quantitiesToCompute; // $unregistered

    void addQuantity(IResultsFunctionConstSP quantity,
                     IHypothesis::IDistanceMetricConstSP distanceMetric,
                     int instrument,
                     RiskQuantityEvaluator::ValueSP buffer, int paramIndex) {
        for (int q = 0; q < int(quantitiesToCompute.size()); ++q) {
            QuantityToCompute& qc = quantitiesToCompute[q];
            if (qc.quantity->equalTo(quantity.get()) &&
                qc.instrument == instrument &&
                qc.distanceMetric->equalTo(distanceMetric.get()))
            {
                qc.append(buffer.get(), paramIndex);
                return;
            }
        }

        quantitiesToCompute.push_back(QuantityToCompute(quantity, instrument,
                                                        distanceMetric));
        quantitiesToCompute.back().append(buffer.get(), paramIndex);
    }

    //@}

    /**
     * RiskQuantityFactory's to be expanded at this node.  In the root node
     * this is the greeks we're actually asked for; in branches and leaves this
     * is almost always empty unless some of the greeks return nonempty
     * IRiskQuantityFactory::lazies().
     */

    IRiskQuantityFactoryArraySP greeks; // $unregistered

    /**
     * Sub-branches of this tree node
     */

    //@{

    struct Branch {
        AtomicHypothesisConstSP hyp;
        HypothesisTreeSP subtree;

        Branch(AtomicHypothesisConstSP hyp, HypothesisTree *subtree):
            hyp(hyp),
            subtree(subtree)
        {}

        Branch() {}
    };

    vector<Branch> branches; // $unregistered

    /**
     * Return (or if it doesn't exist create) the node in the tree representing
     * the state of the world obtained by applying @a hyp to the state of the
     * world at this node.
     *
     * @param i   only nonzero when called recursively
     */

    HypothesisTree *subtree(IHypothesisConstSP hyp, int i = 0) {
        ASSERT(i <= hyp->numAtomics());

        if (hyp->numAtomics() == i) {
            return this;
        }
        else {
            int b;
            for (b = 0;
                 b < int(branches.size()) &&
                     !branches[b].hyp->equalTo(hyp->atomic(i).get());
                 ++b);

            if (b == int(branches.size())) {
                branches.push_back(Branch(hyp->atomic(i),
                                          new HypothesisTree(values)));
            }

            return branches[b].subtree->subtree(hyp, i + 1);
        }
    }

    //@}

    /**
     * Mark all quantities to be computed at this node and those below it as
     * "untweakable"
     */

    void noteOops(UntweakableConstSP oops) {
        for (int q = 0; q < int(quantitiesToCompute.size()); ++q) {
            quantitiesToCompute[q].noteOops(oops);
        }

        for (int b = 0; b < int(branches.size()); ++b) {
            branches[b].subtree->noteOops(oops);
        }

        // FIXME what about lazies ... ?
    }

    HypothesisTree(RiskQuantityEvaluator::ValueArraySP values):
        CObject(TYPE),
        values(values),
        greeks(new IRiskQuantityFactoryArray())
    {}

    static bool allRequestsExist(OutputRequestArrayConstSP requests,
                                 CResultsArrayConstSP results) {
        for (int r = 0; r < requests->size(); ++r) {
            for (int i = 0; i < results->size(); ++i) {
                if (!(*results)[i]->exists((*requests)[r].get())) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Compute some OutputRequest's if they aren't already present in a
     * previously computed set of Results
     *
     * In evaluating some greeks taken with respect to output requests, like
     * NakedBondRhoParallel, we need to ensure that the request in question
     * is actually computed, which it may not be either because
     *
     *    -  we are in a tweaked world, with Control::isPricing() == false:
     *       computation of output requests is conventionally suppressed in
     *       this case; or
     *
     *    -  isPricing == true but we are using previously computed results
     *       [this is actually what always happens at the moment] and the
     *       output wasn't explicitly requested.
     *
     * If it's not there, we need to build a fresh Control, pretending we're in
     * the untweaked world, and reprice.  We need to make sure to use a cloned
     * Model because e.g. rebuilding of tree structure is keyed off the same
     * isPricing flag.  For the same reason the actual price we compute here
     * may be different than in the precomputed results, so for consistency we
     * overwrite the former with the latter.  [We don't go as far as doing two
     * pricings, one with isPricing == true and the other with isPricing ==
     * false, in the case that 'precomputed' is absent.]
     */

    static CResultsArrayConstSP requestResults(
            MultiTweakGroupSP world, CControlSP control,
            OutputRequestArrayConstSP requests,
            CResultsArrayConstSP precomputed) {
        if (!!precomputed && allRequestsExist(requests, precomputed)) {
            return precomputed;
        }
        else {
            TRACE_METHOD;
            TRACE("Re-price to obtain missing output requests");

            CResultsArraySP rs;
            if (!precomputed) rs = world->getInstruments()->emptyResults();
            else rs.reset(copy(precomputed.get()));

            CControlSP c2(new Control(
                SensitivityArray::SP(), OutputRequestArray::SP(),
                false, "", control->skipShifts()));

            for (int r = 0; r < requests->size(); ++r) {
                OutputRequestSP req = (*requests)[r];

                for (int i = 0; i < rs->size(); ++i) {
                    if (!(*rs)[i]->exists(req.get()) &&
                            !c2->requestsOutput(req->getRequestName())) {
                        c2->addRequest(req);
                    }
                }
            }

            c2->calculateMulti(IModelSP(copy(world->getModel())).get(),
                               world->getInstruments(), rs);

            if (!!precomputed) {
                for (int i = 0; i < rs->size(); ++i) {
                    if ((*precomputed)[i]->priceExists()) {
                        (*rs)[i]->storePrice(
                            (*precomputed)[i]->retrievePrice(),
                            (*precomputed)[i]->getCcyName());
                    }
                }
            }

            return rs;
        }
    }

    bool maybeRelevant() const {
        return !greeks->empty() ||
               !branches.empty() ||
               QuantityToCompute::anyMaybeRelevant(quantitiesToCompute);
    }

public:

    HypothesisTree(IRiskQuantityFactoryArrayConstSP greeks,
                   RiskQuantityEvaluator::ValueArraySP values):
        CObject(TYPE),
        values(values),
        greeks(!greeks ? new IRiskQuantityFactoryArray() :
                         new IRiskQuantityFactoryArray(*greeks))
    {}

    void evaluate(MultiTweakGroupSP world, CControlSP control,
                  RiskMappingConstSP riskMapping,
                  // Would like to make these 2 optional but VC71 gives totally
                  // spurious errors
                  CResultsArraySP resultss /* = CResultsArraySP() */,
                  AtomicHypothesisArraySP history /* = AtomicHypothesisArray::SP() */,
                  DoubleArraySP distances = DoubleArray::SP(),
                  bool found = false) {
        try {
            TRACE_METHOD;
            TRACE("We're in a world subjected to " << historyString(history));

            if (!greeks->empty()) {
                TRACE_BLOCK("The user has asked for some greeks:");

                for (int l = 0; l < greeks->size(); ++l) {
                    TRACE_BLOCK((*greeks)[l]->getClass()->getName());

                    NamedRiskQuantityArraySP rqs =
                        (*greeks)[l]->riskQuantities(world, riskMapping);

                    if (!!rqs) for (int q = 0; q < rqs->size(); ++q) {
                        HypotheticalQuantityArrayConstSP params =
                            (*rqs)[q]->riskQuantity->parameters();

                        TRACE_BLOCK("For \"" << *(*rqs)[q]->resultsName << "\":");

                        for (int i = 0; i < world->getInstruments()->size(); ++i) {
                            RiskQuantityEvaluator::ValueSP b =
                                RiskQuantityEvaluator::Value::SP((*rqs)[q], i);
                            values->push_back(b);
                            for (int p = 0; p < params->size(); ++p) {
                                IHypothesisConstSP hyp = (*params)[p]->hypothesis();

                                TRACE("It wants us to evaluate " <<
                                      *(*params)[p]->quantity() << ", "
                                      "assuming " << *hyp);

                                IHypothesisConstSP mappedHyp =
                                    !riskMapping ? hyp : riskMapping->mapped(hyp);

                                if (mappedHyp.get() != hyp.get()) {
                                    TRACE("RiskMapping says instead to do " <<
                                          *mappedHyp);
                                    b->usedRiskMapping = true;
                                }

                                TRACE("... adding that to the tree for future "
                                      "evaluation");

                                subtree(mappedHyp)->addQuantity(
                                    (*params)[p]->quantity(),
                                    mappedHyp->distanceMetric(), i, b, p);
                            }
                        }
                    }

                    // We don't use this feature at the moment

                    LazyRiskQuantityFactoryArraySP lgs = (*greeks)[l]->lazies(world);

                    if (!!lgs) for (int g = 0; g < lgs->size(); ++g) {
                        subtree((*lgs)[g]->hypothesis)->greeks->push_back(
                            (*lgs)[g]->greek);
                    }
                }
            }

            // Split quantities into those which involve OutputRequest's and
            // those which don't, and collate the output requests required

            OutputRequestArraySP outputReqs(new OutputRequestArray());
            vector<QuantityToCompute*> priceQtys, outputReqQtys;

            for (int q = 0; q < int(quantitiesToCompute.size()); ++q) {
                OutputRequestArrayConstSP reqs =
                    quantitiesToCompute[q].quantity->outputRequests();
                if (!reqs || reqs->empty()) {
                    priceQtys.push_back(&quantitiesToCompute[q]);
                }
                else {
                    outputReqQtys.push_back(&quantitiesToCompute[q]);
                    for (int r = 0; r < reqs->size(); ++r) {
                        outputReqs->push_back((*reqs)[r]);
                    }
                }
            }

            // Compute qtys depending on price only (i.e. not on output reqs)

            if (QuantityToCompute::anyMaybeRelevant(priceQtys)) {
                if (!resultss) {
                    TRACE("Now we price the instrument(s) in this world state");

                    resultss = world->getInstruments()->emptyResults();
                    world->getInstruments()->Price(
                        world->getModel(), control.get(), resultss);
                }
                else {
                    TRACE("We don't have to do any pricing: can use the "
                          "Results passed in from the Control");
                }

                TRACE_BLOCK("and store the prices in buffers");

                for (int q = 0; q < int(priceQtys.size()); ++q) {
                    priceQtys[q]->compute(
                        resultss, distances, history->size() == 0 || found);
                }
            }

            // Compute qtys depending on output requests

            if (QuantityToCompute::anyMaybeRelevant(outputReqQtys)) {
                CResultsArrayConstSP rs =
                    requestResults(world, control, outputReqs, resultss);

                TRACE_BLOCK("Also store OutputRequest values in buffers");

                for (int q = 0; q < int(outputReqQtys.size()); ++q) {
                    outputReqQtys[q]->compute(
                        rs, distances, history->size() == 0 || found);
                }
            }

            // Recursively evaluate sub-branches

            for (int b = 0; b < int(branches.size()); ++b) if (branches[b].subtree->maybeRelevant()) {
                IHypothesis::AlternateWorldSP alt;
                UntweakableConstSP oops;

                try {
                    TRACE("Applying hypothesis \"" << *branches[b].hyp << "\"");

                    alt = branches[b].hyp->appliedTo(world);

                    AtomicHypothesisArraySP historyPlus(
                        new AtomicHypothesisArray(*history));
                    historyPlus->push_back(
                        AtomicHypothesisSP::constCast(branches[b].hyp));
                    DoubleArraySP distancesPlus(new DoubleArray(*distances));
                    distancesPlus->push_back(alt->distance);

/**/                const AbstractPropertyTweakHypothesis* th =
/**/                        dynamic_cast<const AbstractPropertyTweakHypothesis*>(branches[b].hyp.get());
/**/
/**/                SensitivitySP oldCurrentSens;
/**/                if (th && !!th->_originatingSensitivity) {
/**/                    th->_originatingSensitivity->setCurrentHypothesis(
/**/                        history, AbstractPropertyTweakHypothesisConstSP(th), alt->oldValue);
/**/                    oldCurrentSens = control->currentSens;
/**/                    control->currentSens = SensitivitySP::dynamicCast(th->_originatingSensitivity);
/**/                }
/**/
                    bool isPricing = control->isPricingRun;
                    control->isPricingRun = false;
                    TRACE("Recursing into hypothetical world ...");
                    try {
                        branches[b].subtree->evaluate(
                            MultiTweakGroupSP::dynamicCast(alt->world), control,
                            riskMapping, CResultsArraySP(),
                            historyPlus, distancesPlus, found || alt->found);
                    }
                    catch (...) {
                        control->isPricingRun = isPricing;
/**/                    if (th) th->_originatingSensitivity.reset();
                        throw;
                    }
                    control->isPricingRun = isPricing;
/**/                if (th) th->_originatingSensitivity.reset();
/**/
/**/                if (!!oldCurrentSens) control->currentSens = oldCurrentSens;
                }
                catch (ModelException& e) {
                    oops.reset(new Untweakable(e));
                }
                catch (exception& e) {
                    oops.reset(new Untweakable(string(e.what())));
                }

                if (!!alt) alt->undo();

                if (!!oops) {
                    branches[b].subtree->noteOops(oops);
                }
            }

        }
        catch (ModelException& e) {
            throw ModelException(e, "HypothesisTree::evaluate()");
        }
    }
};

CClassConstSP const HypothesisTree::TYPE = CClass::registerClassLoadMethod(
    "HypothesisTree", typeid(HypothesisTree), load);



// 
// ***********************
//  RiskQuantityEvaluator
// ***********************
// 

RiskQuantityEvaluator::RiskQuantityEvaluator(bool riskMappingEnabled):
    CObject(TYPE),
    riskMappingEnabled(riskMappingEnabled),
    riskMapping(new RiskMapping())
{}

RiskQuantityEvaluator::RiskQuantityEvaluator(RiskMappingSP riskMapping):
    CObject(TYPE),
    riskMappingEnabled(true),
    riskMapping(riskMapping)
{}

RiskQuantityEvaluator::~RiskQuantityEvaluator() {}

void RiskQuantityEvaluator::getMarket(
        IModelConstSP model, MarketDataConstSP market,
        IInstrumentCollectionConstSP instruments) {
    try {
        riskMapping = riskMappingEnabled ?
                          RiskMapping::fromMarket(model, instruments, market) :
                          RiskMapping::SP();
    }
    catch (exception& e) {
        throw ModelException(e,
            "RiskQuantityEvaluator::getMarket()",
            "Initialising for a " + model->getClass()->getName() + " model");
    }
}

static void setupSensitivities(IRiskQuantityFactoryArrayConstSP greeks,
                               Control* control, IModel* model) {
    for (int g = 0; g < greeks->size(); ++g) {
        Sensitivity* sens = dynamic_cast<Sensitivity*>((*greeks)[g].get());
        if (sens) {
            sens->setControl(control);
            sens->setAlgorithm(model);
        }
    }
}

RiskQuantityEvaluator::ValueArraySP RiskQuantityEvaluator::values(
        IRiskQuantityFactoryArrayConstSP greeks,
        MultiTweakGroupSP world) const {
    try {
        CControlSP control(new Control(
            SensitivityArray::SP(), OutputRequestArray::SP(), false, ""));
            
/**/    setupSensitivities(greeks, control.get(), world->getModel());

        // Buffers for storing the values needed for evaluating the required
        // RiskQuantity's

        ValueArraySP values(new ValueArray());

        // Tree of all the different states we want to put the "world" into

        HypothesisTree tree(greeks, values);

        // Run through all those states evaluating the HypotheticalQuantity's
        // needed for the RiskQuantity's

        tree.evaluate(world, control, riskMapping,
                      CResultsArraySP(), AtomicHypothesisArray::SP());

        return values;
    }
    catch (ModelException& e){
/**/    setupSensitivities(greeks, NULL, NULL);
        throw ModelException(e, "RiskQuantityEvaluator::values()");
    }
/**/setupSensitivities(greeks, NULL, NULL);
}

/**
 * Entry point for evaluating some greeks
 *
 * <H3>Compatibility: provision for currentSens and getControl enquiries</H3>
 *
 * There are several mechanisms in the code, such as MC "quick Greeks", which
 * use Control::getCurrentSensitivity() to find out which Sensitivity is currently
 * being evaluated, and others which use Sensitivity::getControl() and
 * Sensitivity::getModel() to find out which Control is associated with a given
 * Sensitivity.
 *
 * Those usages can be expressed (more simply and extensibly) in terms of
 * looking at the differences between the base state of the world and the
 * tweaked state being priced, as a CompoundHypothesis.  However, I don't want
 * to have to do everything at once, so for now I've put in a compatibility
 * mechanism.
 *
 * (The compatibility-related code is marked / * * /.)
 *
 * Each AbstractPropertyTweakHypothesis knows which Sensitivity it plays a part
 * in evaluating, via its field
 * AbstractPropertyTweakHypothesis::_originatingSensitivity of type
 * ICompatibilitySensitivity.  The field is set on all hypothesis extracted
 * from the IRiskQuantityFactory's by
 * HypothesisTree::setOriginatingSensitivity().  It's read in
 * HypothesisTree::evaluate(), so that
 *
 *    -  the Sensitivity can be caused to be returned from
 *       Control::getCurrentSensitivity() for the duration of the pricing, and
 *    -  its shift size, market data name and whatever else look as
 *       expected to anyone who asks, via
 *       ICompatibilitySensitivity::setCurrentHypothesis().
 *
 * We also make set the Control and the Algorithm on the Sensitivity in
 * case anyone asks.
 */

void RiskQuantityEvaluator::storeResults(
        IRiskQuantityFactoryArrayConstSP greeks,
        MultiTweakGroupSP world,
        CControlSP control,
        CResultsArraySP resultss) const {
    try {
        TRACE_METHOD;

        ASSERT(resultss->size() == world->getInstruments()->size() ||
               resultss->size() == 0);

/**/    setupSensitivities(greeks, control.get(), world->getModel());

        // Buffers for storing the values needed for evaluating the required
        // RiskQuantity's

        ValueArraySP values(new ValueArray());

        // Tree of all the different states we want to put the "world" into

        HypothesisTree tree(greeks, values);

        // Run through all those states evaluating the HypotheticalQuantity's
        // needed for the RiskQuantity's; reuse previously computed results for
        // base case

        tree.evaluate(world, control, riskMapping,
                      resultss->size() == 0 ? CResultsArraySP() : resultss,
                      AtomicHypothesisArray::SP());

        {
            TRACE_BLOCK("Transfer values from buffers to Results, "
                        "for external consumption");

            if (resultss->size() == 0) {
                *resultss = *world->getInstruments()->emptyResults();
            }

            for (int v = 0; v < values->size(); ++v) {
                (*values)[v]->storeResult(resultss);
            }
        }

        OutputRequest* riskmapped = control->requestsOutput(OutputRequest::RISK_MAPPED_SENS);
        if (riskmapped && !riskmapped->getHasFinished()) {
            StringArrayArraySP rms(new StringArrayArray);
            rms->push_back(StringArraySP(new StringArray));  // greeks
            rms->push_back(StringArraySP(new StringArray));  // market data names

            hash_set<string, Hashtable::StringHash> seen;

            for (int v = 0; v < values->size(); ++v) {
                if ((*values)[v]->usedRiskMapping && (*values)[v]->ok()) {
                    const IResultsIdentifier* id =
                        (*values)[v]->riskQuantity->resultsName.get();
                    string k = id->packet() + " " +
                                   (!id->entry() ? id->packet() :
                                                   id->entry()->toString());
                    if (seen.find(k) == seen.end()) {
                        seen.insert(k);
                        (*rms)[0]->push_back(id->packet());
                        (*rms)[1]->push_back(!id->entry() ? id->packet() :
                                                   id->entry()->toString());
                    }
                }
            }

            for (int r = 0; r < resultss->size(); ++r) {
                (*resultss)[r]->storeRequestResult(riskmapped, rms);
            }
        }

        // Hack: harmless but annoying test diffs arise the fact that we now
        // include DELTA_SHIFT_SIZE even if the corresponding DELTA is
        // Untweakable.  Maybe I should fix these in the tests but for now ...

        for (int i = 0; i < resultss->size(); ++i) {
            ResultsSP results = (*resultss)[i];
            if (results->packetExists("DELTA_SHIFT_SIZE")) {
                OutputNameArraySP names =
                    results->packetContents("DELTA_SHIFT_SIZE");
                for (int n = 0; n < names->size(); ++n) {
                    try {
                        results->retrieveScalarGreek("DELTA", (*names)[n]);
                    }
                    catch (exception&) {
                        results->removeGreek("DELTA_SHIFT_SIZE", (*names)[n]);
                    }
                }
            }
        }

        // Bit of extra info for use in test cases

        if (!!riskMapping && resultss->size() > 0) {
            riskMapping->storeRequestResults(control, (*resultss)[0]);
        }
    }
    catch (ModelException& e){
/**/    setupSensitivities(greeks, NULL, NULL);
        throw ModelException(e, "IRiskQuantityFactory::evaluate");
    }

/**/setupSensitivities(greeks, NULL, NULL);
}

static void evaluateExposures(
                IRiskQuantityFactoryArrayConstSP greeks,
                MultiTweakGroupSP world,
                RiskMappingConstSP riskMapping,
                CResultsArraySP resultss) {
    try {
        TRACE_METHOD;

        if (!greeks->empty()) {
            TRACE_BLOCK("The user has asked for some exposure characteristics:");

            // Extract the quantities we're being asked to compute in the base
            // world state and store dummy sensitivities to indicate exposure

            IRiskQuantityFactoryArraySP subgreeks(new IRiskQuantityFactoryArray());
            double indicator;

            for (int l = 0; l < greeks->size(); ++l) {
                TRACE_BLOCK((*greeks)[l]->getClass()->getName());

                NamedRiskQuantityArraySP rqs = (*greeks)[l]->riskQuantities(
                                                                world, riskMapping);

                if (!!rqs) for (int q = 0; q < rqs->size(); ++q) {
                    TRACE_BLOCK("For \"" << *(*rqs)[q]->resultsName << "\":");

                    RiskQuantityConstSP riskQuantity((*rqs)[q]->riskQuantity);
                    if (riskQuantity->isExceptional()) {
						ResultsSP results;
						try {
							// Pretty ugly, but our RiskQuantity will throw an exception
							// and we'll get its "value" from that.
							riskQuantity->value(CDoubleArray(), CDoubleArray());
						}
						catch (NotApplicableException&)
						{
							TRACE_BLOCK("Store "
								<< *(*rqs)[q]->resultsName
								<< " NotApplicable");
							for (int i = 0; i < world->getInstruments()->size(); ++i)
								(*rqs)[q]->resultsName->storeNotApplicableToName(
															(*resultss)[i]);
						}
						catch (ModelException& me)
						{
							TRACE_BLOCK("Store "
								<< *(*rqs)[q]->resultsName
								<< " Untweakable");
							for (int i = 0; i < world->getInstruments()->size(); ++i)
								(*rqs)[q]->resultsName->storeUntweakable(
									(*resultss)[i],
									UntweakableConstSP(new Untweakable(me)));
						}
                    }
                    else {
                        // We assume that any non-constant RiskQuantity indicates
                        // exposure to the associated resultsName and indicate that
                        // exposure with a non-zero value.
                        indicator = riskQuantity->isConstant() ? 0 : 1;
                        TRACE_BLOCK((indicator ? "" : "no ") << "exposure");

                        for (int i = 0; i < world->getInstruments()->size(); ++i)
                            (*rqs)[q]->resultsName->store((*resultss)[i], indicator);
                    }
                }

                // Record any lazy greeks for subsequent processing
                LazyRiskQuantityFactoryArraySP lgs = (*greeks)[l]->lazies(world);
                if (!!lgs)
                    for (int g = 0; g < lgs->size(); ++g)
                        subgreeks->push_back((*lgs)[g]->greek);   
            }

            // Recursively process exposures for any (lazy) subgreeks
            if (subgreeks->size() > 0)
                evaluateExposures(subgreeks, world, riskMapping, resultss);
        }
    }
    catch (ModelException& me) {
        throw ModelException(me, "HypothesisTree::evaluateExposures()");
    }
}

void RiskQuantityEvaluator::storeExposureResults(
        IRiskQuantityFactoryArrayConstSP greeks,
        MultiTweakGroupSP world,
        CControlSP control,
        CResultsArraySP resultss) const {
    try {
        ASSERT(resultss->size() == world->getInstruments()->size() ||
               resultss->size() == 0);

        /**/    setupSensitivities(greeks, control.get(), world->getModel());

        if (resultss->size() == 0) {
            *resultss = *world->getInstruments()->emptyResults();
        }

        evaluateExposures(greeks, world, riskMapping, resultss);
    }
    catch (ModelException& e){
        /**/    setupSensitivities(greeks, NULL, NULL);
        throw ModelException(e, "IRiskQuantityFactory::evaluate");
    }

    /**/setupSensitivities(greeks, NULL, NULL);
}

static IObject* defaultRiskQuantityEvaluator() {
    return new RiskQuantityEvaluator();
}

void RiskQuantityEvaluator::load(CClassSP& clazz) {
    REGISTER(RiskQuantityEvaluator, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultRiskQuantityEvaluator);
    FIELD(riskMappingEnabled, "riskMappingEnabled");
    FIELD(riskMapping, "riskMapping");
    FIELD_MAKE_TRANSIENT(riskMapping);
}

CClassConstSP const RiskQuantityEvaluator::TYPE = CClass::registerClassLoadMethod(
    "RiskQuantityEvaluator", typeid(RiskQuantityEvaluator), load);

DRLIB_END_NAMESPACE
