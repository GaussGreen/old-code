//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DerivativeAsset.cpp
//
//   Description : interface class for a derivative asset (e.g. AssetCVB)
//
//   Author      : Jay Blumenstein
//
//   Date        : 18 Sep 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DerivativeAsset.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

/** Invoked when Class is 'loaded' - needed to use Addin feature in load */
static void loadDerivativeAsset(CClassSP& clazz) {
    REGISTER_INTERFACE(IObject, clazz);
    EXTENDS(IObject);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

/** Identifier for MTM results when we are doing ASSET_THEO_MTM. The
    results live in the instrument packet using the supplied identifier */
string DerivativeAsset::ASSET_MTM_ID = "ASSET_MTM";

CClassConstSP const DerivativeAsset::TYPE = 
CClass::registerInterfaceLoadMethod(
    "DerivativeAsset", typeid(DerivativeAsset), loadDerivativeAsset);


/** This class is used to drive the updating of any DerivativeAssets with
    the current control. Used once before each pricing */
class DerivativeAsset::UpdateAssetWithControlAndModel: 
    virtual public ObjectIteration::IAction {
private:    
    CControlSP       ctrl;
    HashtableSP      assetModels;
public:
    UpdateAssetWithControlAndModel(CControl*          ctrl, 
                                   const HashtableSP& assetModels):
        // copy control as asset might overwrite (eg outputRequest status) as
        // it calculates its theo price. Could probably restrict the copy
        // to the case when isPricing() is true but safer this way.
        ctrl(copy(ctrl)), assetModels(assetModels){}
    
    virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj){
        IObjectSP object(state.getObject()); // get non const access
        DerivativeAssetSP asset(DerivativeAssetSP::dynamicCast(object));
        const string& name = asset->getName(); // get name
        IModelSP model(IModelSP::dynamicCast(assetModels->get(name)));
        asset->setControlAndModel(ctrl, model);
        return true;
    }
};

/** This class is used to drive the updating of any DerivativeAssets with
    the value of the useTheoAssetPrice flag. This is done once before each
    set of price+greek calculations */
class DerivativeAsset::UpdateAssetWithUseTheoFlag: 
    virtual public ObjectIteration::IAction {
private:
    bool useTheoAssetPrice;
public:
    UpdateAssetWithUseTheoFlag(bool             useTheoAssetPrice):
        useTheoAssetPrice(useTheoAssetPrice) {};
    
    /** ObjectIteration calls this each time an object of type DerivativeAsset
        is found */
    virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj){
        IObjectSP object(state.getObject()); // get non const access
        DerivativeAssetSP::dynamicCast(object)->
            setUseTheoAssetPrice(useTheoAssetPrice);
        return true;
    }
};

/** This class is used to drive the updating of any DerivativeAssets with
    the value of the useTheoAssetPrice flag. This is done once before each
    set of price+greek calculations */
class DerivativeAsset::FindDerivAssets: 
    virtual public ObjectIteration::IActionConst {
public:
    FindDerivAssets(): foundDerivAsset(false){}
    
    /** ObjectIteration calls this each time an object of type DerivativeAsset
        is found */
    virtual bool invoke(const ObjectIteration::State& state,IObjectConstSP obj){
        foundDerivAsset = true;
        state.quitRecursion(true); // get out now
        return false;
    }
    //// were any derivative assets found
    bool derivAssetFound(){
        return foundDerivAsset;
    }
private:
    bool foundDerivAsset;
};

/** This is a 'private' model used by the infrastructure to ensure
    that any 'derivative assets' are updated with what the current
    control is before each pricing. This is necessary since some
    sensitivites may create their own control while the model used by
    the derivative asset may need to know what the control is (eg for
    better greeks). It also manages the assets model (used to
    calculate asset's theo price). */
class DerivativeAsset::Model: public CModel,
                       virtual public ObjectIteration::IActionConst {
private:
    IModelSP       realModel; // how to price the insts
    HashtableSP    assetModels; // how to calc theo price of asset
public:
    static CClassConstSP const TYPE;

    virtual ~Model(){}

    /** ObjectIteration calls this each time an object of type DerivativeAsset
        is found */
    virtual bool invoke(const ObjectIteration::State& state,IObjectConstSP obj){
        IObjectSP object(IObjectSP::constCast(obj)); // we respect constness
        DerivativeAssetSP asset(DerivativeAssetSP::dynamicCast(object));
        const string& name = asset->getName(); // get name
        assetModels->put(name, asset->getModel()); // store in hash
        return true;
    }
    /** Constructs our model - stores reference to each assets model. This
        means that when it is cloned so will the models */
    Model(const IModelSP&      model, 
          IInstrumentCollectionSP insts): CModel(TYPE), 
        realModel(model), assetModels(new Hashtable()){
        try{
            // get hold of all the assetModels
            ObjectIteration iter(DerivativeAsset::TYPE);
            iter.recurse(*this, insts);
        } catch (exception& e){
            throw ModelException(e, "DerivativeAsset::Model");
        }
    }

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results) {
        PriceMulti(IInstrumentCollection::singleton(instrument),
                   CControlSP::attachToRef(control),
                   CResultsArraySP(new CResultsArray(1, CResultsSP::attachToRef(results))));
    }

    /** calculate single price and store result in CResult */
    virtual void PriceMulti(IInstrumentCollectionSP instruments,
                            CControlSP control, 
                            CResultsArraySP results){
        // update any DerivativeAsset with the new control
        ObjectIteration iter(DerivativeAsset::TYPE); // create iteration
        UpdateAssetWithControlAndModel action(control.get(), assetModels);
        iter.recurse(action, instruments);
        // then do normal price
        realModel->PriceMulti(instruments, control.get(), results);
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
        return realModel->GetMarket(market, name, type);
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
        return realModel->modifyMarketData(market, clazz, mo);
    }

    /** Invoked after instrument has got its market data. Allows model to
        get any extra data required. Default implementation does nothing */
    virtual void getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instruments){
        realModel->getMarket(market, instruments);
    }

    /** override a control shift (eg for delta on trees)
        returns true if new control is constructed else returns 0 */
    virtual SensControl* AlterControl(
        const SensControl* currSensControl) const{
        return realModel->AlterControl(currSensControl);
    }

    /** called after a series of pricings to indicate that the object
        will not be used again to calculate any
        sensitivities. Basically, state information can be stored
        inside the Model - some of which might be expensive in terms
        of memory. Since the model is constructed by clients, we do
        not control when the Model object is freed. This gives a
        chance for Models to free any expensive caches. The default
        implementation does nothing */
    virtual void flush(){
        realModel->flush();
    }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const {
        return realModel->wantsRiskMapping();
    }

private:
    Model(): CModel(TYPE){}
    
    static IObject* defaultModel(){
        return new Model();
    }
    
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        // note: not visible to EAS/spreadsheet
        REGISTER(Model, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultModel);
        FIELD(realModel, "Model actually used for pricing");
        FIELD(assetModels, "Models used by DerivativeAssets "
              "to calc theo price");
    }
};        

CClassConstSP const DerivativeAsset::Model::TYPE = 
CClass::registerClassLoadMethod(
    "DerivativeAsset::Model", typeid(DerivativeAsset::Model), load);

/** Creates a model that should be used to price a product which has
    a DerivativeAsset in it. Note copy of model is not made */
IModel* DerivativeAsset::createModel(const IModelSP&      model, 
                                     IInstrumentCollectionSP insts){
    return new Model(model, insts);
}

/** Returns true if the instrument contains derivatives assets */
bool DerivativeAsset::derivativeAssetsExist(IInstrumentCollectionSP insts){
    ObjectIteration iter(DerivativeAsset::TYPE);
    FindDerivAssets action;
    iter.recurse(action, insts);
    return action.derivAssetFound();
}

/** Sets all DerivativeAssets in inst to use supplied value for
    whether to use theoretical price or mtm. Returns false if no
    DerivativeAssets were found */
void DerivativeAsset::setUseTheoAssetPrice(IInstrumentCollectionSP insts,
                                           bool useTheoAssetPrice){
    ObjectIteration iter(DerivativeAsset::TYPE);
    UpdateAssetWithUseTheoFlag action(useTheoAssetPrice);
    iter.recurse(action, insts);
}

    


DRLIB_END_NAMESPACE
