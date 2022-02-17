//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AggregateModel.cpp
//
//   Description : Way of pricing composite instruments in one go
//
//   Author      : Mark A Robson
//
//   Date        : 30 Oct 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AggregateModel.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE
IAggregateModel::~IAggregateModel(){}

void IAggregateModel::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IAggregateModel, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IAggregateModel::TYPE = CClass::registerInterfaceLoadMethod(
    "IAggregateModel", typeid(IAggregateModel), load);

class IAggregateModel::MyAggInst: public Instrument,
                       public virtual ISensitiveStrikes,
                       public virtual LastSensDate{
    CInstrumentArray  inst;
    friend class MyAggregateModel;

    MyAggInst(): Instrument(TYPE) {}

    static IObject* defaultConstructor(){
        return new MyAggInst();
    }
    // type information
    static void load(CClassSP& clazz){
        // private class
        REGISTER(MyAggInst, clazz);
        SUPERCLASS(Instrument);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(inst, "The component instruments");
    }
public:
    static CClassConstSP const TYPE;
    MyAggInst(const CInstrumentArray& inst): Instrument(TYPE), inst(inst){}

    // implementation below
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP   mkt);

    /** Called once before the initial pricing */
    virtual void Validate(){
        for (int i = 0; i < inst.size(); i++){
            inst[i]->Validate();
        }
    }
        
    virtual CSensControl* AlterControl(
        const IModel*          modelParams,
        const CSensControl*    sensControl) const; // defined below

    /** Returns the value date (aka today) the instrument is currently
        pricing for */
    virtual DateTime getValueDate() const{
        return inst.front()->getValueDate();
    }


    /** Returns the name of the discount currency of the instrument currently
     * pricing for */
    string discountYieldCurveName() const {
        return inst.front()->discountYieldCurveName();
    }

        
    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const{
        DateTime theEnd;
        bool     endIsSet = false;
        for (int i = 0; i < inst.size(); i++){
            LastSensDate* lsd = dynamic_cast<LastSensDate*>(inst[i].get());
            if (lsd) {
                DateTime newEnd(lsd->endDate(sensControl));
                if (!endIsSet || newEnd.isGreater(theEnd)){
                    theEnd = newEnd;
                    endIsSet = true;
                }
            }
        }
        if (endIsSet){
            return theEnd;
        }
        // ugh !
        MaturityPeriod ages("50Y");
        return ages.toDate(inst.front()->getValueDate());
    }

    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    virtual bool avoidVegaMatrix(const IModel*model); // defined below

    /** Returns all strikes the ComputeEstimator is sensitve to  */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model); 
    // defined below
};

class IAggregateModel::MyAggregateModel: public Model,
                       public virtual IAggregateModel{
    CModelArray models;
    mutable int currentInstIndex; // when doing getMarket etc
    friend class MyAggInst;

    MyAggregateModel(): CModel(TYPE), currentInstIndex(-1){}

    static IObject* defaultConstructor(){
        return new MyAggregateModel();
    }
    // type information
    static void load(CClassSP& clazz){
        // private class
        REGISTER(MyAggregateModel, clazz);
        SUPERCLASS(Model);
        IMPLEMENTS(IAggregateModel);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(models, "The component models");
        FIELD_NO_DESC(currentInstIndex);
        FIELD_MAKE_TRANSIENT(currentInstIndex);
    }
public:
    static CClassConstSP const TYPE;
    // simple constructor - no copying
    MyAggregateModel(const CModelArray& models):
        // note getMarket is also used for the control/model so at times the
        // current instrument index is not clear
        Model(TYPE), models(models), currentInstIndex(0){}

    /** Whether to enable RiskMapping: see IModel::wantsRiskMapping() */
    virtual WantsRiskMapping wantsRiskMapping() const {
        {for (int m = 0; m < models.size(); ++m) {
            if (models[m]->wantsRiskMapping() == riskMappingDisallowed) {
                return riskMappingDisallowed;
            }
        }}

        {for (int m = 0; m < models.size(); ++m) {
            if (models[m]->wantsRiskMapping() == riskMappingAllowed) {
                return riskMappingAllowed;
            }
        }}

        return riskMappingIrrelevant;
    }

    /** Returned object will be passed to price method */
    virtual Instrument* createAggregateInstrument(
        const CInstrumentArray& instruments,
        const DoubleArray&      weights){
        // we don't use the weights
        return new MyAggInst(instruments);
    }

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results){
        // we should never be here anyway
        if (models.size() != 1){
            throw ModelException("MyAggregateModel::Price",
                                 "Internal error");
        }
        CResultsArray resultsArray(1, CResultsSP::attachToRef(results));
        price(instrument, control, resultsArray);
    }

    /** calculate n prices and store result in CResult */
    virtual void price(Instrument*          compositeInstrument, 
                       Control*             control, 
                       CResultsArray&       results){
        // cast to our instrument
        MyAggInst& myCompInst = 
            dynamic_cast<MyAggInst&>(*compositeInstrument);
        for (int i = 0; i < results.size(); i++){
            models[i]->Price(myCompInst.inst[i].get(), control,
                             results[i].get());
        }
    }
    /** Instruments call this to select correct piece of market data from
        cache */
    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const{
        // need to use right model
        return models[currentInstIndex]->GetMarket(market, name, type);
    }

    virtual MarketObjectSP modifyMarketData(
        const MarketData*     market,
        const CClassConstSP&  clazz,  // what type was originally requested
        const MarketObjectSP& mo) const{ /* what GetMarket returned or what
                                            was "inline" already */
        // need to use right model
        return models[currentInstIndex]->modifyMarketData(market, 
                                                          clazz, mo);
    }
    /** Invoked after instrument has got its market data. Allows model to
        get any extra data required. */
    virtual void getMarket(const MarketData*    market,
                           IInstrumentCollectionSP instruments){
        // cast to our instrument
        //# warning FIXME this is sticking plaster
        ASSERT(instruments->size() > 0);
        MyAggInst* myCompInst = dynamic_cast<MyAggInst *>((*instruments)[0].get());
        ASSERT(myCompInst);

        for (int i = 0; i < models.size(); i++){
            models[i]->getMarket(market, IInstrumentCollection::singleton(myCompInst->inst[i]));
        }
    }

    /** override a control shift (eg for delta on trees)
        returns true if new control is constructed else returns 0 */
    virtual SensControl* AlterControl(
        const SensControl* currSensControl) const{
        // problematic because potentially each model might alter it to
        // something different
        for (int i = 0; i < models.size(); i++){
            SensControl* ctrl = models[i]->AlterControl(currSensControl);
            if (ctrl){
                return ctrl;
            }
        }
        return 0;
    }
    virtual void flush(){
        for (int i = 0; i < models.size(); i++){
            models[i]->flush();
        }
    }
};


void IAggregateModel::MyAggInst::GetMarket(
    const IModel*                model, 
    const CMarketDataSP         mkt){
    const MyAggregateModel& myAggModel =
        dynamic_cast<const MyAggregateModel&>(*model);
    for (int i = 0; i < inst.size(); i++){
        myAggModel.currentInstIndex = i;
        myAggModel.models[i]->getInstrumentAndModelMarket(mkt.get(), inst[i].get());
        myAggModel.currentInstIndex = 0; // set back to default value
    }
}

SensControl* IAggregateModel::MyAggInst::AlterControl(
    const IModel*          model,
    const CSensControl*    sensControl) const{
    // problematic because potentially each instrument might alter it to
    // something different
    const MyAggregateModel& myAggModel =
        dynamic_cast<const MyAggregateModel&>(*model);
    for (int i = 0; i < inst.size(); i++){
        SensControlSP sensCtrl(
            inst[i]->AlterControl(myAggModel.models[i].get(), sensControl));
        if (sensCtrl.get()){
            return sensCtrl.release();
        }
    }
    return 0;
}

/** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
bool IAggregateModel::MyAggInst::avoidVegaMatrix(const IModel*model){
    const MyAggregateModel& myAggModel =
        dynamic_cast<const MyAggregateModel&>(*model);
    for (int i = 0; i < inst.size(); i++){
        ISensitiveStrikes* sensStrikes = 
            dynamic_cast<ISensitiveStrikes*>(inst[i].get());
        if (!sensStrikes || 
            sensStrikes->avoidVegaMatrix(myAggModel.models[i].get())){
            return true;
        }
    }
    return false;
}

/** Returns all strikes the ComputeEstimator is sensitve to  */
DoubleArraySP IAggregateModel::MyAggInst::getSensitiveStrikes(
    OutputNameConstSP outputName,
    const IModel*      model){
    DoubleArraySP allStrikes(new DoubleArray());
    const MyAggregateModel& myAggModel =
        dynamic_cast<const MyAggregateModel&>(*model);
    for (int i = 0; i < inst.size(); i++){
        ISensitiveStrikes& sensStrikes =
            dynamic_cast<ISensitiveStrikes&>(*inst[i]);
        DoubleArraySP strikes(sensStrikes.
                              getSensitiveStrikes(outputName, 
                                                  myAggModel.models[i].get()));
        allStrikes->insert(allStrikes->begin(),
                           strikes->begin(), strikes->end());
    }
    return allStrikes;
}

 
CClassConstSP const IAggregateModel::MyAggInst::TYPE = 
CClass::registerClassLoadMethod(
    "IAggregateModel::MyAggInst", typeid(MyAggInst), load);

CClassConstSP const IAggregateModel::MyAggregateModel::TYPE = 
CClass::registerClassLoadMethod(
    "IAggregateModel::MyAggregateModel", typeid(MyAggregateModel),
    load);


DRLIB_END_NAMESPACE
