//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AggregateInstrument.hpp
//
//   Description : Defines the idea of an 'aggregate instrument'
//
//   Author      : Mark A Robson
//
//   Date        : 10 Nov 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RegressionTest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/ScaleOutputs.hpp"
#include "edginc/AggregateModel.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/ResultsSet.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE
/** A class that represents the idea of pricing an array of instruments in one
    go using one model. The inputs are similar to CompositeInstrument but there
    is only one model and one control */
class AggregateInstrument: public CObject,
                           virtual public IRegressionTest,
                           virtual public ClientRunnable{
    // fields
    ScenarioSP         scenario; // applied before any calcs done
    IAggregateModelSP  model;    // special type of model
    CInstrumentArraySP inst;     // which can price array of instruments
    CControlSP         ctrl;     // same ctrl for all instruments
    DoubleArraySP      multiplier;
    DoubleArraySP      weight;
    CMarketDataSP      market;

    AggregateInstrument();
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
    class ModelAdapter;
    class InstrumentAdapter;
public:
    static CClassConstSP const TYPE;

    // EdrAction
    virtual IObjectSP run();

    /** Runs 'regression test' */
    virtual IObjectSP runTest() const;
};

/** Make one instrument look like several - all with the same behaviour.
    Note that the lack of support for proxies makes this cumbersome. We also
    have to implement/wrap all possible methods to ensure we get the same
    behaviour */
class AggregateInstrument::InstrumentAdapter: public CInstrument,
                           public virtual ISensitiveStrikes,
                           public virtual LastSensDate{
public:
    /** Class to manipulate the results returned from the AggregateModel */
    class PropertiesCache{
        friend class InstrumentAdapter;
        vector<bool>            validateOK;
        vector<int>             validateIdx;
        // there's probably a nice way of doing this with templates or
        // treating everything as an object etc. In fact if the language had
        // support for proxies then it would be very succinct
        DateTimeArray           valueDates;
        vector<bool>            valueDatesOK;
        vector<int>             valueDatesIdx;

        CBoolArray              deadInsts;
        vector<bool>            deadInstsOK;
        vector<int>             deadInstsIdx;

        DateTimeArray           endDates;
        vector<bool>            endDatesOK;
        vector<int>             endDatesIdx;

        CBoolArray              avoidVegaMatrix;
        vector<bool>            avoidVegaMatrixOK;
        vector<int>             avoidVegaMatrixIdx;

        DoubleArrayArray        sensStrikes;
        vector<bool>            sensStrikesOK;
        vector<int>             sensStrikesIdx;

        SensControlArray        sensControl;
        vector<bool>            sensControlOK;
        vector<int>             sensControlIdx;

        StringArray             discountYieldCurveNames;
        vector<bool>            discountYieldCurveNamesOK;
        vector<int>             discountYieldCurveNamesIdx;

        bool                    failed;
        string                  errorMsg;

        // record error and throw exception;
        void internalError(const string& method, const string& m){
            if (!failed){
                failed = true;
                errorMsg = method + ": Internal error. "+m;
            }
            throw ModelException(errorMsg);
        }

    public:
        // confirms internal state is ok. Call before doing any operations
        void check(){
            if (failed){
                throw ModelException(errorMsg);
            }
        }

        void checkIdx(int index, int length){
            if (index >= length){
                internalError("PropertiesCache::checkIdx",
                              "Internal error");
            }
        }
        /** Checks that all objects in the cache have been used */
        void validateAllResultsUsed(){
            static const string method("AggregateInstrument::InstrumentAdapter"
                                       "::validateAllResultsUsed");
            for (unsigned int i = 0; i < valueDatesIdx.size(); i++){
                if (validateIdx[i] != (int) validateOK.size()){
                    internalError(method, "Not all value dates used");
                }
                if (valueDatesIdx[i] != valueDates.size()){
                    internalError(method, "Not all value dates used");
                }
                if (deadInstsIdx[i] != deadInsts.size()){
                    internalError(method, "Not all deadInsts used");
                }
                if (endDatesIdx[i] != (int) endDates.size()){
                    internalError(method, "Not all endDates used");
                }
                if (avoidVegaMatrixIdx[i] != avoidVegaMatrix.size()){
                    internalError(method, "Not all avoidVegaMatrix used");
                }
                if (sensStrikesIdx[i] != sensStrikes.size()){
                    internalError(method, "Not all sensStrikesIdx used");
                }
                if (sensControlIdx[i] != (int) sensControl.size()){
                    internalError(method, "Not all sensControl used");
                }
                if (discountYieldCurveNamesIdx[i] != (int) discountYieldCurveNames.size()){
                    internalError(method, "Not all discountYieldCurveName used");
                }

            }
        }

        PropertiesCache(int numInst):
            validateIdx(numInst), valueDatesIdx(numInst),
            deadInstsIdx(numInst), endDatesIdx(numInst),
            avoidVegaMatrixIdx(numInst), sensStrikesIdx(numInst),
            sensControlIdx(numInst),
            discountYieldCurveNamesIdx(numInst),
            failed(false)
        {}

    };

private:
    InstrumentSP                 theInst;
    refCountPtr<PropertiesCache> cache; /* shared across InstrumentAdapters $unregistered */
    int                          index; // which one we're on
    friend class  ModelAdapter;

    InstrumentAdapter(const InstrumentSP&                 theInst,
                      const refCountPtr<PropertiesCache>& cache,
                      int index):
        CInstrument(TYPE), theInst(theInst), cache(cache), index(index){}

    InstrumentAdapter(): CInstrument(TYPE),cache(){}

    static IObject* defaultConstructor(){
        return new InstrumentAdapter();
    }
    static void load(CClassSP& clazz){
        // not public
        REGISTER(InstrumentAdapter, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(theInst, "Wrapped Instrument");
        FIELD(index, "Index of corresponding instrument");
    }

    void throwError(const string& m) const {
        throw ModelException(m, "Original call on pricing run failed");
    }

public:
    static CClassConstSP const TYPE;

    //// overridden to ensure cache is preserved
    IObject* clone() const{
        IObject* copy = Instrument::clone();
        dynamic_cast<InstrumentAdapter&>(*copy).cache = cache;
        return copy;
    }

     static CInstrumentArraySP create(
        const InstrumentSP& inst,
        int   numComponents){
        refCountPtr<PropertiesCache> cache(new PropertiesCache(numComponents));
        CInstrumentArraySP insts(new CInstrumentArray(numComponents));
        for (int i = 0; i < numComponents; i++){
            (*insts)[i] = CInstrumentSP(new InstrumentAdapter(inst, cache, i));
        }
        return insts;
    }

    /** Checks that all objects in the cache have been used */
    void validateAllResultsUsed(){
        cache->check();
        cache->validateAllResultsUsed();
    }

    // implementation below
    virtual void GetMarket(const IModel* model, const CMarketDataSP market);

    /** Called once before the initial pricing */
    virtual void Validate(){
        bool success = true;
        cache->check();
        if (index == 0){
            try{
                theInst->Validate();
            } catch (exception&){
                success = false;
            }
            cache->validateOK.push_back(success);
        } else {
            cache->checkIdx(cache->validateIdx[index],
                            cache->validateOK.size());
            success = cache->validateOK[cache->validateIdx[index]];
        }
        cache->validateIdx[index]++;
        if (!success){
            throwError("AggregateInstrument::InstrumentAdapter::Validate");
        }
    }

    /** override a control shift (eg for delta on trees) - may return
        null to use original. Default implementation returns null. */
    virtual CSensControl* AlterControl(
        const IModel*          model,
        const CSensControl*    sensControl) const;
    // implementation below

    /** Returns the value date (aka today) the instrument is currently
        pricing for */
    virtual DateTime getValueDate() const{
        cache->check();
        DateTime myValDate;
        bool success = true;
        if (index == 0){
            try{
                myValDate = theInst->getValueDate();
            } catch(exception&){
                success = false;
            }
            cache->valueDatesOK.push_back(success);
            cache->valueDates.push_back(myValDate); // store regardless
        } else {
            cache->checkIdx(cache->valueDatesIdx[index],
                            cache->valueDates.size());
            success = cache->valueDatesOK[cache->valueDatesIdx[index]];
            myValDate = cache->valueDates[cache->valueDatesIdx[index]];
        }
        cache->valueDatesIdx[index]++;
        if (!success){
            throwError("AggregateInstrument::InstrumentAdapter::getValueDate");
        }
        return myValDate;
    }

    /** price a dead instrument until settlement - exercised, expired,
        knocked out etc.  returns true if it is dead (and priced), false
        if it is not dead */
    virtual bool priceDeadInstrument(CControl* control,
                                     CResults* results) const{
        cache->check();
        bool success = true;
        bool dead = false;
        if (index == 0){
            try{
                dead = theInst->priceDeadInstrument(control, results);
            } catch(exception&){
                success = false;
            }
            cache->deadInstsOK.push_back(success);
            cache->deadInsts.push_back(dead);
        } else {
            cache->checkIdx(cache->deadInstsIdx[index],
                            cache->deadInsts.size());
            success = cache->valueDatesOK[cache->deadInstsIdx[index]];
            dead = cache->deadInsts[cache->deadInstsIdx[index]];
        }
        cache->deadInstsIdx[index]++;
        if (!success){
            throwError("AggregateInstrument::InstrumentAdapter::"
                       "priceDeadInstrument");
        }
        return dead;
    }

    /** when to stop tweaking */
    DateTime endDate(const Sensitivity* sensControl) const {
        cache->check();
        DateTime myEndDate;
        bool success = true;
        if (index == 0){
            LastSensDate* lsd = dynamic_cast<LastSensDate*>(theInst.get());
            try{
                if (lsd) {
                    myEndDate =  lsd->endDate(sensControl);
                } else {
                    // ugh !
                    MaturityPeriod ages("50Y");
                    myEndDate = ages.toDate(theInst->getValueDate());
                }
            } catch (exception&){
                success = false;
            }
            cache->endDatesOK.push_back(success);
            cache->endDates.push_back(myEndDate);
        } else {
            cache->checkIdx(cache->endDatesIdx[index],
                            cache->endDates.size());
            success = cache->endDatesOK[cache->endDatesIdx[index]];
            myEndDate = cache->endDates[cache->endDatesIdx[index]];
        }
        cache->endDatesIdx[index]++;
        if (!success){
            throwError("AggregateInstrument::InstrumentAdapter::endDate");
        }
        return myEndDate;
    }

    /** Returns the name of the instrument's discount currency */
    virtual string discountYieldCurveName() const;

    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    virtual bool avoidVegaMatrix(const IModel*model); // defined below

    /** Returns all strikes the ComputeEstimator is sensitve to  */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);
    // defined below
};
CClassConstSP const AggregateInstrument::InstrumentAdapter::TYPE =
CClass::registerClassLoadMethod(
    "AggregateInstrument::InstrumentAdapter", typeid(InstrumentAdapter), load);


/** Class to wrap a single AggregateModel to give N models which can then
    be used to price a composite. It makes the AggregateModel look like
    it is pricing the components of a composite one at a time. */
class AggregateInstrument::ModelAdapter: public Model{
public:
    /** Class to manipulate the results returned from the AggregateModel */
    class ResultsCache{
        friend class ModelAdapter;
        vector<CResultsArraySP> results;
        vector<ModelException>  exceptions;
        vector<unsigned int>    exceptionsIdx;
        vector<vector<bool> >   outputRequestsDone;
        vector<unsigned int>    outputRequestsDoneIdx;
        SensControlArray        sensControl;
        vector<bool>            sensControlOK;
        vector<int>             sensControlIdx;
        bool                    failed;
        string                  errorMsg;
    private:
        // record error and throw exception;
        void internalError(const string& method, const string& m){
            if (!failed){
                failed = true;
                errorMsg = method + ": Internal error. "+m;
            }
            throw ModelException(errorMsg);
        }
    public:
        // confirms internal state is ok. Call before doing any operations
        void check(){
            if (failed){
                throw ModelException(errorMsg);
            }
        }

        void checkIdx(int index, int length){
            if (index >= length){
                internalError("AggregateInstrument::ModelAdapter::checkIdx",
                              "Internal error");
            }
        }
        //// Creates an empty cache for a composite with numInst components.
        ResultsCache(int numInst): results(numInst), exceptionsIdx(numInst),
                                   outputRequestsDoneIdx(numInst),
                                   sensControlIdx(numInst),
                                   failed(false){
            for (int i = 0; i < numInst; i++){
                results[i] =  CResultsArraySP(new CResultsArray(0));
            }
        }

        //// Creates an array of results ready to call AggregateModel::price
        CResultsArraySP newResultsSet(){
            CResultsArraySP newResults(new CResultsArray(results.size()));
            for (unsigned int i = 0; i < results.size(); i++){
                results[i]->push_back(CResultsSP(new Results()));
                (*newResults)[i] =  results[i]->back();
            }
            return newResults;
        }

        // record exception
        void recordException(const ModelException& e){
            exceptions.push_back(e);
            // then clear out all the results for the current pricing run
            for (unsigned int i = 0; i < results.size(); i++){
                results[i]->back() = CResultsSP();
            }
        }

        /** Removes the Results object at the front of the cache for the
            specified component. The object removed is returned */
        CResultsSP deleteResultsFromFront(int index){
            static const string m("AggregateInstrument::ModelAdapter::"
                                  "deleteResultsFromFront");
            static const string s("Run out of results");
            CResultsArraySP& instResult = results[index];
            if (instResult->empty()){
                internalError(m, s);
            }
            CResultsSP frontResult((instResult->front()));
            instResult->erase(instResult->begin());
            if (!frontResult){
                // there was an exception
                if (exceptionsIdx[index] >= exceptions.size()){
                    internalError(m, s);
                }
                exceptionsIdx[index]++;
                throw exceptions[exceptionsIdx[index]-1];
            }
            return frontResult;
        }
        /** record which output requests have been calcualted */
        void recordOutputRequestsStatus(
            const OutputRequestArrayConstSP& requests){
            outputRequestsDone.push_back(vector<bool>(requests->size()));
            vector<bool>&   outputRequestsStatus = outputRequestsDone.back();
            // record which output requests have been calculated
            for (unsigned int j = 0; j < outputRequestsStatus.size(); j++){
                outputRequestsStatus[j] = (*requests)[j]->getHasFinished();
            }
        }

        /** Update the supplied requests to flag which output requests
            have been calculated */
        void getOutputRequestsStatus(const OutputRequestArrayConstSP& requests,
                                     int                              index){
            static const string m("AggregateInstrument::ModelAdapter::"
                                  "getOutputRequestsStatus");
            unsigned int& requestIdx = outputRequestsDoneIdx[index];
            if (requestIdx >= outputRequestsDone.size()){
                internalError(m, "Ran out of output requests status");
            }
            vector<bool>&   outputRequestsStatus =
                outputRequestsDone[requestIdx];
            if ((int) outputRequestsStatus.size() != requests->size()){
                internalError(m, "Wrong number of output requests");
            }
            for (unsigned int j = 0; j < outputRequestsStatus.size(); j++){
                (*requests)[j]->setHasFinished(outputRequestsStatus[j]);
            }
            requestIdx++;
        }
        /** Checks that all Results objects in the cache have been used */
        void validateAllResultsUsed(){
            static const string method("AggregateInstrument::ModelAdapter::"
                                       "validateAllResultsUsed");
            for (unsigned int i = 0; i < results.size(); i++){
                if (!results[i]->empty()){
                    internalError(method, "Not all results used");
                }
                if (exceptionsIdx[i] != exceptions.size()){
                    internalError(method, "Not all exceptions used");
                }
                if (outputRequestsDoneIdx[i] != outputRequestsDone.size()){
                    internalError(method, "Not all output requests used");
                }
                if (sensControlIdx[i] != (int) sensControl.size()){
                    internalError(method, "Not all sensControl used");
                }
            }
        }
    };

private:
    IAggregateModelSP         aggModel;
    int                       index;  /* different index for each instrument
                                         in composite */
    refCountPtr<ResultsCache> resultsCache; // shared across ModelAdapters $unregistered
    friend class  InstrumentAdapter;
public:
    static CClassConstSP const TYPE;

    /** Creates an array of ModelAdapters which can be used to price
        a composite in the regular way */
    static CModelArraySP create(const IAggregateModelSP& aggModel,
                                int                      numComponents){
        CModelArraySP models(new CModelArray(numComponents));
        refCountPtr<ResultsCache> cache(new ResultsCache(numComponents));
        for (int i = 0; i < numComponents; i++){
            (*models)[i] = IModelSP(new ModelAdapter(aggModel, i, cache));
        }
        return models;
    }

    /** Checks that all Results objects in the cache have been used */
    void validateAllResultsUsed(){
        resultsCache->check();
        resultsCache->validateAllResultsUsed();
    }

    //// overridden to ensure cache is preserved
    IObject* clone() const{
        IObject* copy = Model::clone();
        dynamic_cast<ModelAdapter&>(*copy).resultsCache = resultsCache;
        return copy;
    }

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument,
                       CControl*     control,
                       CResults*     results){
        resultsCache->check();
        OutputRequestArrayConstSP requests(control->getOutputRequests());
        if (index == 0){
            // price all instruments in composite now
            CResultsArraySP resultSet(resultsCache->newResultsSet());
            try{
                // our InstrumentAdapter must be used in conjunction with the
                // ModelAdapter class
                InstrumentAdapter& instAdapter =
                    dynamic_cast<InstrumentAdapter&>(*instrument);
                aggModel->price(instAdapter.theInst.get(), control, *resultSet);
                resultsCache->recordOutputRequestsStatus(requests);
            } catch (exception& e){
                resultsCache->recordException(ModelException(e));
                // we will throw this exception below
            }
        } else {
            // don't do anything
        }
        // this will throw the original exception if that was what was
        // returned originally
        ResultsSP resultsToUse(resultsCache->deleteResultsFromFront(index));
        // copy results over
        results->overwrite(resultsToUse.get());
        // and update the output requests in the control
        resultsCache->getOutputRequestsStatus(requests, index);
    }

    /** Returns a [deep] copy of the market data with supplied name
        and type from the given market data cache. This gives the
        model a chance to choose a specific type of market data rather
        than just a general instance. For example, the method could
        request a Black-Scholes Vol rather than just any old vol. The
        default implementation provided by CModel just asks the market
        data for the object of the given type */
    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const{
        resultsCache->check();
        return aggModel->GetMarket(market, name, type);
    }

    /** Invoked for each piece of market data (whether already inline in
        instrument or pulled from cache). The default implementation just
        returns mo. Derived classes can use this to replace market data
        objects with other instances. It compliments GetMarket in that this
        method works for instruments which have the market data inside them
        already */
    virtual MarketObjectSP modifyMarketData(
        const MarketData*     market,
        const CClassConstSP&  clazz,  // what type was originally requested
        const MarketObjectSP& mo) const{ /* what GetMarket returned or what was
                                            "inline" already */
        resultsCache->check();
        return aggModel->modifyMarketData(market, clazz, mo);
    }


    /** Invoked after instrument has got its market data. Allows model to
        get any extra data required. Default implementation does nothing */
     virtual void getMarket(const MarketData*  market,
                            IInstrumentCollectionSP instruments){
        resultsCache->check();
        // our InstrumentAdapter must be used in conjunction with the
        // ModelAdapter class

        //# warning FIXME this is sticking plaster
        ASSERT(instruments->size() > 0);
        InstrumentAdapter* instAdapter = dynamic_cast<InstrumentAdapter *>((*instruments)[0].get());
        ASSERT(instAdapter);

        aggModel->getMarket(market,
                            IInstrumentCollection::singleton(instAdapter->theInst.get()));
     }

    /** override a control shift (eg for delta on trees)
        returns true if new control is constructed else returns 0 */
    virtual SensControl* AlterControl(
        const SensControl* currSensControl) const{
        resultsCache->check();
        //return aggModel->AlterControl(currSensControl);
        SensControl* newControl = 0;
        bool success = true;
        if (index == 0){
            try{
                const IModel& mdl = dynamic_cast<const IModel&>(*aggModel);
                newControl = mdl.AlterControl(currSensControl);
            } catch(exception&){
                success = false;
            }
            resultsCache->sensControlOK.push_back(success);
            resultsCache->sensControl.push_back(
                SensControlSP(newControl? copy(newControl): 0));
        } else {
            resultsCache->checkIdx(resultsCache->sensControlIdx[index],
                                   resultsCache->sensControl.size());
            success = resultsCache->sensControlOK[resultsCache->
                                                  sensControlIdx[index]];
            SensControl* ctrl =
                resultsCache->sensControl[resultsCache->
                                          sensControlIdx[index]].get();
            newControl = ctrl? copy(ctrl): 0;
        }
        resultsCache->sensControlIdx[index]++;
        if (!success){
            throw ModelException("AggregateInstrument::ModelAdapter::"
                                 "AlterControl",
                                 "Original call on pricing run failed");
        }
        return newControl;
    }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * @see IModel::wantsRiskMapping()
     */

    WantsRiskMapping wantsRiskMapping() const {
        return aggModel->wantsRiskMapping();
    }

    /** called after a series of pricings to indicate that the object
        will not be used again to calculate any
        sensitivities. */
    virtual void flush(){
        resultsCache->check();
        if (index == 0){
            aggModel->flush();
        }
    }
private:
    ModelAdapter(const IAggregateModelSP&   aggModel, // shared between adapters
                 int                        index,
                 refCountPtr<ResultsCache>& resultsCache): /* shared between
                                                            adapters */
        Model(TYPE), aggModel(aggModel), index(index),
        resultsCache(resultsCache){
        if (!Model::TYPE->isInstance(aggModel.get())){
            throw ModelException("ModelAdapter", "Aggregate model must be "
                                 "derived from CModel");
        }
    }

    ModelAdapter(): CModel(TYPE) {}

    static IObject* defaultConstructor(){
        return new ModelAdapter();
    }

    static void load(CClassSP& clazz){
        // not public
        REGISTER(ModelAdapter, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(aggModel, "Wrapped AggregateModel");
        FIELD(index, "Index of corresponding instrument");
    }
};
CClassConstSP const AggregateInstrument::ModelAdapter::TYPE =
CClass::registerClassLoadMethod(
    "AggregateInstrument::ModelAdapter", typeid(ModelAdapter), load);

void AggregateInstrument::InstrumentAdapter::GetMarket(
    const IModel*       model,
    const CMarketDataSP market){
    cache->check();
    // our InstrumentAdapter must be used in conjunction with the
    // ModelAdapter class
    const ModelAdapter& mdlAdp = dynamic_cast<const ModelAdapter&>(*model);
    mdlAdp.aggModel->getInstrumentAndModelMarket(market.get(), theInst.get());
}


/** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
bool AggregateInstrument::InstrumentAdapter::avoidVegaMatrix(
    const IModel *model) {
    bool avoid = true;
    cache->check();
    bool success = true;
    if (index == 0){
        ISensitiveStrikes* vm = dynamic_cast<ISensitiveStrikes*>(theInst.get());
        if (vm){
            try{
                // our InstrumentAdapter must be used in conjunction with the
                // ModelAdapter class
                const ModelAdapter& mdlAdp =
                    dynamic_cast<const ModelAdapter&>(*model);
                avoid = vm->avoidVegaMatrix(mdlAdp.aggModel.get());
            } catch(exception&){
                success = false;
            }
        }
        cache->avoidVegaMatrixOK.push_back(success);
        cache->avoidVegaMatrix.push_back(avoid);
    } else {
        cache->checkIdx(cache->avoidVegaMatrixIdx[index],
                        cache->avoidVegaMatrix.size());
        avoid = cache->avoidVegaMatrix[cache->avoidVegaMatrixIdx[index]];
        success = cache->avoidVegaMatrixOK[cache->avoidVegaMatrixIdx[index]];
    }
    cache->avoidVegaMatrixIdx[index]++;
    if (!success){
        throwError("AggregateInstrument::InstrumentAdapter::"
                   "avoidVegaMatrix");
    }
    return avoid;
}

/** Returns all strikes the ComputeEstimator is sensitve to  */
DoubleArraySP AggregateInstrument::InstrumentAdapter::getSensitiveStrikes(
    OutputNameConstSP outputName,
    const IModel*      model){
    cache->check();
    DoubleArraySP strikes;
    bool success = true;
    if (index == 0){
        ISensitiveStrikes* vm = dynamic_cast<ISensitiveStrikes*>(theInst.get());
        try{
            if (!vm) {
                cache->internalError("InstrumentAdapter::getSensitiveStrikes",
                                     "Internal error");
            }
            const ModelAdapter& mdlAdp =
                dynamic_cast<const ModelAdapter&>(*model);
            strikes = vm->getSensitiveStrikes(outputName, mdlAdp.aggModel.get());
        } catch(exception&){
            success = false;
        }
        cache->sensStrikesOK.push_back(success);
        cache->sensStrikes.push_back(*strikes);
    } else {
        cache->checkIdx(cache->sensStrikesIdx[index],
                        cache->sensStrikes.size());
        success = cache->sensStrikesOK[cache->sensStrikesIdx[index]];
        strikes = DoubleArraySP(new DoubleArray(
                                    cache->sensStrikes[cache->
                                                       sensStrikesIdx[index]]));
    }
    cache->sensStrikesIdx[index]++;
    if (!success){
        throwError("AggregateInstrument::InstrumentAdapter::"
                   "getSensitiveStrikes");
    }
    return strikes;
}

/** override a control shift (eg for delta on trees) - may return
    null to use original. Default implementation returns null. */
SensControl* AggregateInstrument::InstrumentAdapter::AlterControl(
    const IModel*          model,
    const SensControl*    sensControl) const{
    SensControl* newControl = 0;
    cache->check();
    bool success = true;
    if (index == 0){
        try{
            // our InstrumentAdapter must be used in conjunction with the
            // ModelAdapter class
            const ModelAdapter& mdlAdp =
                dynamic_cast<const ModelAdapter&>(*model);
            newControl = theInst->AlterControl(mdlAdp.aggModel.get(), sensControl);
        } catch(exception&){
            success = false;
        }
        cache->sensControlOK.push_back(success);
        cache->sensControl.push_back(SensControlSP(newControl?
                                                   copy(newControl): 0));
    } else {
        cache->checkIdx(cache->sensControlIdx[index],
                        cache->sensControl.size());
        success = cache->sensControlOK[cache->sensControlIdx[index]];
        SensControl* ctrl = cache->sensControl[cache->
                                               sensControlIdx[index]].get();
        newControl = ctrl? copy(ctrl): 0;
    }
    cache->sensControlIdx[index]++;
    if (!success){
        throwError("AggregateInstrument::InstrumentAdapter::AlterControl");
    }
    return newControl;
}



/** override a control shift (eg for delta on trees) - may return
    null to use original. Default implementation returns null. */
string AggregateInstrument::InstrumentAdapter::discountYieldCurveName() const {
    cache->check();
    string ycName;
    bool success = true;
    if (index == 0) {
        try {
            ycName = theInst->discountYieldCurveName();
        } catch(exception&){
            success = false;
        }
        cache->discountYieldCurveNamesOK.push_back(success);
        cache->discountYieldCurveNames.push_back(ycName);    // regardless
    } else {
        cache->checkIdx(cache->discountYieldCurveNamesIdx[index],
                        cache->discountYieldCurveNames.size());
        success = cache->discountYieldCurveNamesOK[cache->discountYieldCurveNamesIdx[index]];
        ycName = cache->discountYieldCurveNames[cache->discountYieldCurveNamesIdx[index]];
    }
    cache->discountYieldCurveNamesIdx[index]++;
    if (!success){
        throwError("AggregateInstrument::InstrumentAdapter::discountYieldCurveName");
    }
    return ycName;
}


// for reflection
AggregateInstrument::AggregateInstrument(): CObject(TYPE){
    // empty
}

/** Runs 'regression test' */
IObjectSP AggregateInstrument::run(){
    ctrl->switchOffWriteToFile();
    return runTest();
}

/** Calculates price and sensitivities for given instrument and
    context together with supplied market data + any scenario shifts.
    Results are returned in either a ResultSet object. ( >1 instrument)
    or a Results object (1 instrument) */
IObjectSP AggregateInstrument::runTest() const{
    static const string method = "AggregateInstrument::run";
    try {
        int numInsts = inst->size();
        if (numInsts == 0){
            throw ModelException(method, "Zero instruments supplied");
        }
        if (numInsts != multiplier->size() ||
            numInsts != weight->size()) {
            throw ModelException(method,
                                 "instrument, multiplier"
                                 " & weight arrays must be the same length");
        }
        // create local variable in case market is null
        CMarketDataSP   marketData(market);
        if (!marketData){
            // review in conjunction with avoiding writing all market data out
            marketData = CMarketDataSP(new MarketData());
        }

        CResultsArraySP  results(new CResultsArray(numInsts));
        CControlSP ctrl(this->ctrl); // hide field
        // check for command line options
        if (CommandLineParams::hasParameter(CommandLineParams::Mega)){
            ctrl = // default everything possible
                CControlSP(SensitivityFactory::megaControl());
        }
        if (CommandLineParams::hasParameter(CommandLineParams::AddSens)) {
            IObjectSP param(CommandLineParams::
                            getParameterArgs(CommandLineParams::AddSens));
            const string& s = (dynamic_cast<CString&>(*param)).stringValue();
            SensitivitySP extraSens(SensitivityFactory::defaultSensitivity(s));
            ctrl->addSensitivity(extraSens);
        }
        if (CommandLineParams::hasParameter(CommandLineParams::AddReq)) {
            IObjectSP param(CommandLineParams::
                            getParameterArgs(CommandLineParams::AddReq));
            const string& s = (dynamic_cast<CString&>(*param)).stringValue();
            OutputRequestSP extraRequest(new OutputRequest(s));
            ctrl->addRequest(extraRequest);
        }

        // check that we don't have any scaleable instruments
        for (int j = 0; j < numInsts; j++) {
            if (ctrl->scaleOutputs() &&
                IScaleOutputs::TYPE->isInstance((*inst)[j].get())) {
                throw ModelException(method,
                                     "The scaleResults flag must be false "
                                     "for aggregate instrument components");
            }
        }

        // get the IAggregateModel so create a single instrument from the
        // array of instruments
        DoubleArray  pos(*weight);
        for (int q = 0; q < pos.size(); q++){
            pos[q] *= (*multiplier)[q];
        }
        InstrumentSP  aggInst(model->createAggregateInstrument(*inst, pos));
        if (!aggInst){ // ought to be unnecessary
            throw ModelException(method, "AggregateModel of type "+
                                 model->getClass()->getName()+
                                 " failed to build an aggregate instrument");
        }
        // then set up our model adapters to make the single IAggregateModel
        // look like an array of regular models
        CModelArraySP      aggAdapters(ModelAdapter::create(model, numInsts));
        // then set up our instrument adapters to make the aggInst
        // look like an array of regular instruments
        CInstrumentArraySP instAdapters(InstrumentAdapter::create(aggInst,
                                                                  numInsts));
        // loop over each instrument in the composite
        // run it and stash the results
        for (int i = 0; i < numInsts; i++) {
            CInstrumentSP inst = (*instAdapters)[i];
            IModelSP      model = (*aggAdapters)[i];
            // either run scenario or risk mgr directly
            (*results)[i] = model->go(inst, scenario, ctrl, marketData);
            // scale the result by multiplier
            double thisMultiplier = numInsts == 1?
                ((*multiplier)[0] * (*weight)[0]): (*multiplier)[i];
            if (!Maths::equals(thisMultiplier, 1.0)){
                (*results)[i]->scale(ctrl, thisMultiplier, false);
            }
        }
        // extra validation on adapters
        dynamic_cast<ModelAdapter&>(*aggAdapters->front()).
            validateAllResultsUsed();
        dynamic_cast<InstrumentAdapter&>(*instAdapters->front()).
            validateAllResultsUsed();


        IObject* resultsSet;
        if (numInsts == 1){ // different format of results for n = 1
            resultsSet = (*results)[0].release();
        } else {
            // now we 'just' combine all the results together somehow
            CControlArray controls(numInsts, ctrl);
            ResultsSP combined = Results::combineResults(controls,
                                                         *weight, *results);
            resultsSet = new ResultsSet(combined.release(), results.release());
        }
        return IObjectSP(resultsSet);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Invoked when Class is 'loaded' */
void AggregateInstrument::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(AggregateInstrument, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ClientRunnable);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(scenario, "scenario");
    FIELD_MAKE_OPTIONAL(scenario);
    FIELD(model, "model");
    FIELD(inst, "instrument(s)");
    FIELD(ctrl, "ctrls");
    FIELD(multiplier, "multiplier");
    FIELD(weight, "weighting");
    FIELD(market, "market");
    FIELD_MAKE_OPTIONAL(market);
}

IObject* AggregateInstrument::defaultConstructor(){
    return new AggregateInstrument();
}

CClassConstSP const AggregateInstrument::TYPE =
CClass::registerClassLoadMethod("AggregateInstrument",
                                typeid(AggregateInstrument),
                                AggregateInstrument::load);

// symbol (referenced by AddinLib.cpp) to ensure file gets linked in
bool AggregateInstrumentLinked = true;

DRLIB_END_NAMESPACE

