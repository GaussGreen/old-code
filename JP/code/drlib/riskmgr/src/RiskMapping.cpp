/**
 * @file RiskMapping.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/Format.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/ExpiryAndStrike.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/Model.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IDynamicsParameter.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/IHypothesis.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/RiskMappingMatrix.hpp"
#include "edginc/RiskMapping.hpp"

DRLIB_BEGIN_NAMESPACE

// 
// ************************
//  RiskMappingMatrixCache
// ************************
// 

/**
 * MarketData "plugin" for indexed retrieval of all the RiskMappingMatrix's for
 * a particular parametric vol
 *
 * This object is created by RiskMapping::fromMarket() via the
 * MarketData::plugin() method.  It scans the market and creates an index of
 * the RiskMappingMatrix's contained therein.  Because it's attached to the
 * MarketData, it persists across pricings, and the (possibly expensive) scan
 * doesn't have to be repeated every time.
 *
 * <H3>Entry points</H3>
 *
 * The rmmsForParameters(ObjectArrayConstSP) method returns all the RMM's in
 * the market containing a row labelled with an IRiskAxis which parameterises a
 * tweak to one of the supplied objects.  Called by RiskMapping::fromMarket().
 *
 * In order to support that the constructor SP(const MarketData&) sets up an
 * index of all the RiskMappingMatrix's in the market, keyed by the
 * subjectInterface and marketDataName of the IRiskAxis's in their
 * _outputBasis's.
 *
 *     rmmCache (clazz, name) = {
 *         rmm in MarketData
 *         where exists axis in rmm._outputBasis
 *                      where name = axis.marketDataName and
 *                            axis.abstractProperty.subjectInterface = clazz
 *     }
 *
 * Called via MarketData::plugin() from RiskMapping::fromMarket().
 */

class RiskMappingMatrixCache: public CObject {

    RiskMappingMatrixCache(): CObject(TYPE) {}
    static IObject* emptyShell() { return new RiskMappingMatrixCache(); }

    static void load(CClassSP& clazz) {
        REGISTER(RiskMappingMatrixCache, clazz);
        SUPERCLASS(CObject);
        FIELD(rmmCache, "rmmCache");
        EMPTY_SHELL_METHOD(emptyShell);
    }

    HashtableSP rmmCache; // string -> RiskMappingMatrixArraySP; $required

    static string rmmCacheKey(CClassConstSP clazz, string name) {
        return clazz->getName() + "" + name;
    }

    RiskMappingMatrixCache(const MarketData& market):
        CObject(TYPE),
        rmmCache(new Hashtable())
    {
        try {
            TRACE_METHOD;

            MarketObjectArraySP available = market.GetAllDataWithType(
                RiskMappingMatrix::TYPE);

            for (int m = 0; m < available->size(); ++m) {
                RiskMappingMatrixSP rmm =
                    RiskMappingMatrixSP::dynamicCast((*available)[m]);

                TRACE_BLOCK("Indexing RMM " << rmm->getName());

                for (int o = 0; o < rmm->_outputBasis->size(); ++o) {
                    IRiskAxisConstSP axis = (*rmm->_outputBasis)[o];
                    IAbstractRiskPropertyConstSP prop = axis->abstractProperty();
                    OutputNameConstSP mdn = axis->marketDataName();

                    try {
                        if (!prop) {
                            throw ModelException(__FUNCTION__,
                                "Axis won't say what property it parameterises");
                        }

                        if (!mdn) {
                            throw ModelException(__FUNCTION__,
                                "Axis won't say what market name it's for");
                        }
                    }
                    catch (exception& e) {
                        throw ModelException(e, __FUNCTION__,
                            "outputBasis entry #" + Format::toString(o) + " (" +
                            axis->getClass()->getName() + " \"" +
                            axis->toString() + "\"");
                    }

                    string key = rmmCacheKey(prop->subjectInterface(),
                                             mdn->idGet(0));

                    TRACE("Adding " << key << " -> " << rmm->getName() <<
                          " to index");

                    RiskMappingMatrixArraySP rmms =
                        RiskMappingMatrixArraySP::dynamicCast(
                            rmmCache->get(key));
                    if (!rmms) {
                        rmms.reset(new RiskMappingMatrixArray());
                        rmmCache->put(key, rmms);
                    }

                    int n = rmms->size();
                    while (--n >= 0 && !((*rmms)[n] == rmm));
                    if (n < 0) {
                        rmms->push_back(rmm);
                    }
                    else {
                        TRACE("(but it's been added already)");
                    }
                }
            }
        }
        catch (exception& e) {
            throw ModelException(
                e, "RiskMappingMatrixCache::RiskMappingMatrixCache()");
        }
    }

    /**
     * them <- union(them, RMM's affecting parameters of given clazz/market name)
     *
     * That's
     *
     *     union { rmmCache (c, name)
     *             for c in [ clazz ] + clazz.ancestors
     *                 where c issubclass MarketObject }
     */

    void addRmmsForParameter(CClassConstSP clazz, const string& name,
                             RiskMappingMatrixArraySP them) const {
        TRACE_METHOD;

        if (MarketObject::TYPE->isAssignableFrom(clazz)) {
            TRACE_SHOW(clazz->getName());
            TRACE_SHOW(name);

            RiskMappingMatrixArraySP rmms =
                RiskMappingMatrixArraySP::dynamicCast(
                    rmmCache->get(rmmCacheKey(clazz, name)));

            if (!!rmms) {
                for (int m = 0; m < rmms->size(); ++m) {
                    RiskMappingMatrixSP rmm = (*rmms)[m];
                    TRACE("Found " << rmm->getName());
                    int n = them->size();
                    while (--n >= 0 && !((*them)[n] == rmm));
                    if (n < 0) {
                        them->push_back(RiskMappingMatrixSP(rmm.clone()));
                    }
                    else {
                        TRACE("(but it would be a duplicate)");
                    }
                }
            }

#           if MarketObjectCanBeAnInterface
                TRACE_BLOCK("Considering interfaces");

                CClassVec ifaces = clazz->getInterfaces();
                for (size_t i = 0; i < ifaces.size(); ++i) {
                    addRmmsForParameter(ifaces[i], name, them);
                }

#           endif

            TRACE_BLOCK("Considering superclass");
            if (clazz->getSuperClass()) {
                addRmmsForParameter(clazz->getSuperClass(), name, them);
            }
        }
        else {
            TRACE(clazz->getName() << " is not a MarketObject");
        }
    }

public:

    static const CClassConstSP TYPE;

    /**
     * Construct a RiskMappingMatrixCache for a given market data cache
     */

    static IObjectConstSP SP(const MarketData& market) {
        return IObjectConstSP(new RiskMappingMatrixCache(market));
    }

    /**
     * The RiskMappingMatrix's in the market which might affect the given model
     * parameters (i.e. stoch vol objects)
     *
     * @param params    A list of model parameter objects, like VolSVJ,
     *                  SRMEQ::Vol, ...
     *
     * @return  All the RMM's in the market containing a row labelled with
     *          an IRiskAxis which parameterises a tweak to one of the @a
     *          params.  Typically the axes in question will be FieldRiskAxis's
     *          with explicitly FieldRiskAxis::className = param->getClass()
     *          and FieldRiskAxis::_marketDataName = param->getName().
     */

    RiskMappingMatrixArraySP rmmsForParameters(
            ObjectArrayConstSP params) const {
        try {
            TRACE_METHOD;

            RiskMappingMatrixArraySP them(new RiskMappingMatrixArray());

            for (int p = 0; p < params->size(); ++p) {
                MarketObjectSP pa = MarketObjectSP::dynamicCast((*params)[p]);
                TRACE_BLOCK("Getting RMMs for object " << pa->getClass()->getName()
                            << " \"" << pa->getName() << "\"");
                addRmmsForParameter(pa->getClass(), pa->getName(), them);
            }

            return them;
        }
        catch (exception& e) {
            throw ModelException(e, "RiskMappingMatrixCache::rmmsForParameters");
        }
    }
};

CClassConstSP const RiskMappingMatrixCache::TYPE = CClass::registerClassLoadMethod(
  "RiskMappingMatrixCache", typeid(RiskMappingMatrixCache), RiskMappingMatrixCache::load);

// 
// =============
//  RiskMapping
// =============
// 

RiskMapping::RiskMapping(RiskMappingMatrixArrayConstSP matrices):
    CObject(TYPE),
    matrices(matrices)
{}

RiskMapping::RiskMapping():
    CObject(TYPE),
    matrices(new RiskMappingMatrixArray())
{}

IObject* RiskMapping::defaultOne() {
    return new RiskMapping(RiskMappingMatrixArrayConstSP());
}

RiskMappingSP RiskMapping::SP() {
    return RiskMappingSP(new RiskMapping());
}

RiskMappingSP RiskMapping::fromMarket(IModelConstSP model,
                                      IInstrumentCollectionConstSP instruments,
                                      MarketDataConstSP market) {
    TRACE_METHOD;
    try {
        MultiTweakGroupSP tg(new MultiTweakGroup(
            IInstrumentCollectionSP::constCast(instruments),
            IModelSP::constCast(model)));

        ObjectArrayConstSP params = SensMgrConst(tg.get()).all(
                                        IDynamicsParameter::TYPE,
                                        (ITweakNameResolver*)0);

        if (!params) {
            TRACE("Didn't find any IDynamicsParameter's: no RiskMapping needed");
            return RiskMappingSP();
        }
        else {
            TRACE("Found at least one IDynamicsParameter (" <<
                  (*params)[0]->getClass()->getName() << ") "
                  "so can do RiskMapping");

            return RiskMappingSP(new RiskMapping(
                market->plugin<RiskMappingMatrixCache>()->
                    rmmsForParameters(params)));
        }
    }
    catch (exception& e) {
        throw ModelException(e, "RiskMapping::fromMarket()");
    }
}

RiskMapping::~RiskMapping() {}

IHypothesisConstSP RiskMapping::mapped(IHypothesisConstSP hypothesis) const {
    TRACE_METHOD;

    try {
        TRACE_SHOW(*hypothesis);

        for (int m = 0; m < matrices->size(); ++m) {
            IHypothesisConstSP mh = (*matrices)[m]->mapped(hypothesis);

            // We assume that the matrices are disjoint, so that the first one
            // which maps hypothesis to a different object is the only one
            // we need to consider

            if (mh.get() != hypothesis.get()) {
                TRACE("RiskMappingMatrix " << (*matrices)[m]->getName() <<
                      " applied, returning transformed hypothesis");
                return mh;
            }
        }

        TRACE("None of our RiskMappingMatrix's applied");
        return hypothesis;
    }
    catch (exception& e) {
        throw ModelException(e, "RiskMapping::mapped()");
    }
}

OutputNameArrayConstSP RiskMapping::subjectNames(
        IAbstractRiskPropertyConstSP property,
        IObjectConstSP world) const {
    TRACE_METHOD;
    TRACE("Figuring out the names we can tweak w.r.t. " << *property);

    OutputNameArraySP names(new OutputNameArray());

    for (int m = 0; m < matrices->size(); ++m) {
        (*matrices)[m]->appendMappedNames(property, names);
    }

    if (names->empty()) {
        TRACE("None of our RiskMappingMatrix's applied; "
              "using names which are directly sensitive to " << *property);
    }
    else {
        TRACE("Using names obtained from our RiskMappingMatrix's");
    }

    return names->empty() ? property->subjectNames(world) : names;
}

template <>
VoidArrayConstSP RiskMapping::subjectQualifiers(
        smartConstPtr<IRiskProperty<Void> > property,
        OutputNameConstSP name,
        IObjectConstSP world) const {
    TRACE_METHOD;

    TRACE("Figuring out what qualifiers we can tweak for property " <<
          *property << " of name " << *name);

    //provide the facility for a property to return a risk mapping matrix
    //on the fly...and use its qualifiers instead
    TRACE("Testing " << *property << " for its own RiskMappingMatrix ...");
    RiskMappingMatrixConstSP rmm = property->riskMappingMatrix(world, name);
    if (!!rmm) {
        //keep hold of the new matrix
        TRACE("Yes, it does...add to the previous set of RiskMappings");
        extendMatrices(rmm);
    }

    //standard case for Void
    TRACE("--- there's only one possible qualifier in this case, i.e. 'Void'");
    return VoidArray::SP(1, VoidSP());
}

template <class Q>
smartConstPtr<array<smartPtr<Q>, Q> > RiskMapping::subjectQualifiers(
        smartConstPtr<IRiskProperty<Q> > property,
        OutputNameConstSP name,
        IObjectConstSP world) const {
    TRACE_METHOD;

    TRACE("Figuring out what qualifiers we can tweak for property " <<
          *property << " of name " << *name);

    smartPtr<array<smartPtr<Q>, Q> > them(new array<smartPtr<Q>, Q>());

    bool any = false;

    for (int m = 0; m < matrices->size(); ++m) {
        any = any ||
            (*matrices)[m]->appendMappedQualifiers(property, name, them);
    }

    //provide the facility for a property to return a risk mapping matrix
    //on the fly...and use its qualifiers instead
    if (!any) {
        TRACE("None of our RiskMappingMatrix's applied; "
              "test " << *property << " for its own instead");
        RiskMappingMatrixConstSP rmm = property->riskMappingMatrix(world, name);
        if (!!rmm) {
            //keep hold of the new matrix
            TRACE("Yes, it does...add to the previous set of RiskMappings");
            extendMatrices(rmm);
            TRACE("Try applying it now");
            any = any ||
                rmm->appendMappedQualifiers(property, name, them);
        }
    }

    if (any) {
        TRACE("Using qualifiers from our RiskMappingMatrix's");
    }
    else {
        TRACE("None of our RiskMappingMatrix's applied; "
              "using qualifiers from the actual property");
    }

    return any ? them : property->subjectQualifiers(world, name);
}

void RiskMapping::extendMatrices(RiskMappingMatrixConstSP mtx) const
{
    //firstly the new mtx needs to be non-const
    RiskMappingMatrixSP ncMtx = RiskMappingMatrixSP::constCast(mtx);
    //as does the array
    RiskMappingMatrixArraySP ncMatrices = RiskMappingMatrixArraySP::constCast(matrices);
    //now we can finally add to the matrices
    ncMatrices->push_back(ncMtx);
}

template ExpiryWindowArrayConstSP RiskMapping::subjectQualifiers(
    smartConstPtr<IRiskProperty<ExpiryWindow> > property,
    OutputNameConstSP name,
    IObjectConstSP world) const;

template ExpiryPairArrayConstSP RiskMapping::subjectQualifiers(
    smartConstPtr<IRiskProperty<ExpiryPair> > property,
    OutputNameConstSP name,
    IObjectConstSP world) const;

template ExpiryAndStrikeArrayConstSP RiskMapping::subjectQualifiers(
    smartConstPtr<IRiskProperty<ExpiryAndStrike> > property,
    OutputNameConstSP name,
    IObjectConstSP world) const;

template BoxedIntArrayConstSP RiskMapping::subjectQualifiers(
    smartConstPtr<IRiskProperty<BoxedInt> > property,
    OutputNameConstSP name,
    IObjectConstSP world) const;

void RiskMapping::storeRequestResults(CControlSP control, ResultsSP results)
        const {
    if (control->requestsOutput(
            OutputRequest::DEBUG_RELEVANT_RISK_MAPPING_MATRICES)) {

        StringArraySP s(new StringArray());
        for (int m = 0; m < matrices->size(); ++m) {
            s->push_back((*matrices)[m]->getName());
        }
        OutputRequest outputRequest(OutputRequest::DEBUG_RELEVANT_RISK_MAPPING_MATRICES);
        results->storeRequestResult(
            &outputRequest,
            s);
    }
}

void RiskMapping::load(CClassSP& clazz) {
    REGISTER(RiskMapping, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultOne);
    FIELD(matrices, "matrices");
}

CClassConstSP const RiskMapping::TYPE = CClass::registerClassLoadMethod(
    "RiskMapping", typeid(RiskMapping), load);

DRLIB_END_NAMESPACE
