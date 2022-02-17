//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLObjViewAddins.cpp
//
//   Author      : André Segger
//
//   Date        : 21 Sep 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvert.hpp"
#include "edginc/XCB.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/Fund.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/RiskMgrInterface.hpp"
#include "edginc/ScenarioInterface.hpp"
#include "edginc/CompositeInstrument.hpp"
#include "edginc/YieldNameCollector.hpp"
#include "edginc/IXLInterfaceMap.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/DeltaToCredit.hpp"
#include "edginc/AssetVegaParallel.hpp"
#include "edginc/ClosedFormFA.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/PseudoSimpleEquity.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/ClientRunnable.hpp"

#include <list>

DRLIB_BEGIN_NAMESPACE

/* non static variable to force linkage of this file. Outside of class to
   avoid necessity of header file */
bool XLObjViewAddinsRegistered = true;

/** Addin to map one object into another. For situations where
    the spreadsheet interface is not the same as the internal C++
    data model */
class XLMapAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP  object;

    /** the 'addin function' */
    static IObjectSP XLmapper(XLMapAddin* params){
        static const string routine = "XLMapAddin::XLmapper";
        
        if (IXLInterfaceMap::TYPE->isInstance(params->object))
        {
            // The object requires converting
            IXLInterfaceMap* objectToMap = dynamic_cast<IXLInterfaceMap*>(params->object.get());
            
            IObjectSP mappedObject(objectToMap->map());
            return mappedObject;
        }
        else
        {
            // no conversion necessary, pass back the original object
            return params->object;
        }
    }

    /** for reflection */
    XLMapAddin():  CObject(TYPE){}

    
 /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLMapAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLMapAddin);
        FIELD(object,"object to be mapped");
        Addin::registerClassObjectMethod(
            "XL_MAP",
            Addin::UTILITIES,
            "maps objects between C++ and spread sheet representations",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)XLmapper);
    }

    static IObject* defaultXLMapAddin() {
        return new XLMapAddin();
    }
};

CClassConstSP const XLMapAddin::TYPE = CClass::registerClassLoadMethod(
    "XLMapAddin", typeid(XLMapAddin), load);

/** Addin to invoke the EDRAction on an object */
class VolParamSetAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP  vol;
    string     fieldName;
    IObjectSP  fieldValue;

    /** the 'addin function' */
    static IObjectSP volParamSet(VolParamSetAddin* params){
        static const string routine = "VolParamSetAddin::VolParamSet";
        try
        {
            IObjectSP volclone(params->vol->clone());
            if (CVolBaseParamSurface::TYPE->isInstance(params->vol.get()))
            {
                // we want to return a copy
                CClassConstSP   c  = volclone.get()->getClass();
                const CFieldArray& fields = c->getDeclaredFields();
                
                for (unsigned int i = 0; i < fields.size(); i++) 
                {
                    string name = fields[i]->getName();
                    int modifiers = fields[i]->getModifiers();
                    if (CString::equalsIgnoreCase(name,params->fieldName) &&
                        !(Modifier::isTransient(modifiers)))
                    {
                        fields[i]->set(volclone, params->fieldValue);
                    }
                }
            }
            return volclone;
        }
        catch (exception& e) 
        {
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    VolParamSetAddin():  CObject(TYPE){}

    
 /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VolParamSetAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVolParamSetAddin);
        FIELD(vol,"invoke action on this object");
        FIELD(fieldName, "field name");
        FIELD(fieldValue," field value");
        Addin::registerClassObjectMethod(
            "VOL_PARAM_SET",
            Addin::UTILITIES,
            "sets the value of a vol parameter",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)volParamSet);
    }

    static IObject* defaultVolParamSetAddin() {
        return new VolParamSetAddin();
    }
};

CClassConstSP const VolParamSetAddin::TYPE = CClass::registerClassLoadMethod(
    "VolParamSetAddin", typeid(VolParamSetAddin), load);


/** Addin to change the parameters of a cloned object 
    calls the validation function on the new object
    the function is not recursive */

class ObjectParamSetAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP       object;
    CStringArraySP  fieldNames;
    ObjectArraySP   fieldValues;

    /** the 'addin function' */
    static IObjectSP ObjectParamSet(ObjectParamSetAddin* params){
        static const string routine = "ObjectParamSetAddin::ObjectParamSet";
        try
        {
            if (params->fieldNames->size()!=params->fieldValues->size())
            {
                throw ModelException(routine, "The  number of fields should be equal to the number of values"); 
            }
            IObjectSP objectClone(params->object->clone());
            
            // we want to return a copy
            CClassConstSP   c  = objectClone.get()->getClass();
            const CFieldArray& fields = c->getDeclaredFields();
            
            for (unsigned int i = 0; i < fields.size(); i++) 
            {
                string name = fields[i]->getName();
                int modifiers = fields[i]->getModifiers();
                
                for (int j = 0; j < (params->fieldNames)->size(); j++) 
                    if (CString::equalsIgnoreCase(name,(*params->fieldNames)[j]) &&
                        !(Modifier::isTransient(modifiers)))
                    {
                        fields[i]->set(objectClone, (*params->fieldValues)[j]);
                    }
            }
            
            objectClone->validatePop2Object();
            
            return objectClone;
        }
        catch (exception& e) 
        {
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    ObjectParamSetAddin():  CObject(TYPE){}

    
 /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ObjectParamSetAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObjectParamSetAddin);
        FIELD(object,"invoke action on this object");
        FIELD(fieldNames, "field names");
        FIELD(fieldValues," field values");
        Addin::registerClassObjectMethod(
            "OBJECT_PARAM_SET",
            Addin::UTILITIES,
            "Change the values of the specified parameters of an object",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)ObjectParamSet);
    }

    static IObject* defaultObjectParamSetAddin() {
        return new ObjectParamSetAddin();
    }
};

CClassConstSP const ObjectParamSetAddin::TYPE = CClass::registerClassLoadMethod(
    "ObjectParamSetAddin", typeid(ObjectParamSetAddin), load);


/** Addin to determine run-time type of a handle */
class XLGetAssetNameAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    AssetSP  asset;

    /** the 'addin function' - construct object from components */
    static string getAssetName(XLGetAssetNameAddin* params){
        static const string routine = "XLGetAssetNameAddin::getAssetName";

        return (params->asset->getName());
    }

    /** for reflection */
    XLGetAssetNameAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetAssetNameAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetAssetNameAddin);
        FIELD(asset, "asset handle");
        Addin::registerClassStringMethod("GET_ASSET_NAME",
                                         Addin::UTILITIES,
                                         "Returns the name of an asset",
                                         TYPE,
                                         (Addin::StringMethod*)getAssetName);
    }

    static IObject* defaultGetAssetNameAddin() {
        return new XLGetAssetNameAddin();
    }
};

CClassConstSP const XLGetAssetNameAddin::TYPE = CClass::registerClassLoadMethod(
        "XLGetAssetNameAddin", typeid(XLGetAssetNameAddin), XLGetAssetNameAddin::load);

/** Addin to determine run-time type of a handle */
class XLIsWrapperObjectAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP  object;

    /** the 'addin function' - construct object from components */
    static bool isWrapper(XLIsWrapperObjectAddin* params){
        static const string routine = "XLIsWrapperObjectAddin::isWrapper";

        return (MarketObjectWrapper::TYPE->isInstance(params->object.get()));
    }

    /** for reflection */
    XLIsWrapperObjectAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLIsWrapperObjectAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultWrapperObjectAddin);
        FIELD(object, "object handle");
        Addin::registerClassBoolMethod("IS_WRAPPER",
                                       Addin::MARKET,
                                       "determines whether the handle represents"
                                       " a market wrapper object",
                                       TYPE,
                                       (Addin::BoolMethod*)isWrapper);
    }

    static IObject* defaultWrapperObjectAddin() {
        return new XLIsWrapperObjectAddin();
    }
};

CClassConstSP const XLIsWrapperObjectAddin::TYPE = CClass::registerClassLoadMethod(
    "XLIsWrapperObjectAddin", typeid(XLIsWrapperObjectAddin), XLIsWrapperObjectAddin::load);

/** Addin to determine run-time type of a handle */
class XLGetCorrelationMatrixAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP      basket;
    CMarketDataSP  market;

    /** the 'addin function' - construct object from components */
    static IObjectSP getMatrix(XLGetCorrelationMatrixAddin* params){
        static const string routine = "XLGetCorrelationMatrixAddin::getMatrix";

        CStringArray          corrNames;
        CAssetWrapperArraySP  assets;
        CDoubleMatrixSP       correlations;

        if (XCB::TYPE->isInstance(params->basket.get()))        {
            XCB* xcb = dynamic_cast<XCB*>(params->basket.get());
            assets       = CAssetWrapperArraySP::attachToRef(&xcb->assets);
            correlations = CDoubleMatrixSP::attachToRef(&xcb->correlations);
        } else if (Fund::TYPE->isInstance(params->basket.get())) {
            Fund* fund = dynamic_cast<Fund*>(params->basket.get());
            assets       = CAssetWrapperArraySP::attachToRef(&fund->components);
            correlations = CDoubleMatrixSP::attachToRef(&fund->correlations);
        }

        CClosedFormLN model("VolSurface");

        if ( correlations->empty() ) {
            CDoubleMatrixSP       mdCorrelations;

            // generate double matrix of correlations
            // First allocate space
            int numAssets = assets->size();
            corrNames.resize((numAssets * numAssets - numAssets)/2);
            mdCorrelations = CDoubleMatrixSP(new DoubleMatrix(numAssets, numAssets));
            // then look up names and values (stored in DoubleMatrix)
            int pos = 0;
            for (int i = 0; i < numAssets; i++) {
                (*mdCorrelations)[i][i] = 1.0;
                for (int j = i + 1; j < numAssets; j++, pos++) {
                    // look up name in cache
                    corrNames[pos] = 
                        params->market->getCorrelationName((*assets)[i].getName(),
                                                           (*assets)[j].getName());
                    // then get hold of the actual object
                    MarketObjectSP corr(params->market->GetData(&model, corrNames[pos], 
                                                                Correlation::TYPE));

                    IObject* obj = dynamic_cast<IObject*>(corr.get());
                    Correlation& correlation =
                        dynamic_cast<Correlation&>(*obj);
                    (*mdCorrelations)[i][j] = correlation.getCorrelation();
                    (*mdCorrelations)[j][i] = (*mdCorrelations)[i][j];
                }
            }
            return mdCorrelations;

        } else {
            return correlations;
        }
    }

    /** for reflection */
    XLGetCorrelationMatrixAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetCorrelationMatrixAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCorrelationMatrixAddin);
        FIELD(basket, "XCB/Fund handle");
        FIELD(market, "Market handle");
        Addin::registerClassObjectMethod("GET_CORRELATION_MATRIX",
                                         Addin::MARKET,
                                         "retrieves the correlation matrix for"
                                         " a basket from the market data cache",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getMatrix);
    }

    static IObject* defaultCorrelationMatrixAddin() {
        return new XLGetCorrelationMatrixAddin();
    }
};

CClassConstSP const XLGetCorrelationMatrixAddin::TYPE = CClass::registerClassLoadMethod(
    "XLGetCorrelationMatrixAddin", typeid(XLGetCorrelationMatrixAddin), 
    XLGetCorrelationMatrixAddin::load);

/** Addin to determine run-time type of a handle */
class XLGetFXRatesAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    CStringArraySP  fxNames;
    CMarketDataSP  market;

    /** the 'addin function' - construct object from components */
    static IObjectSP getFXMatrix(XLGetFXRatesAddin* params){
        static const string routine = "XLGetFXRatesAddin::getFXMatrix";

        CStringArraySP  fxNames = params->fxNames;
        string           fxAssetName;
        CClosedFormLN model("FlatFXVol");
        double spot;

        int numFXAssets = params->fxNames->size();
        CDoubleMatrixSP  fxRates(new DoubleMatrix(numFXAssets, numFXAssets));

        for (int i = 0; i < numFXAssets; i++) {
            (*fxRates)[i][i] = 1.0;
            for (int j = i + 1; j < numFXAssets; j++) {
                if ( params->market->hasFXData((*fxNames)[j], (*fxNames)[i])) {
                    fxAssetName = params->market->getFXName((*fxNames)[j],
                                                            (*fxNames)[i]);
                    MarketObjectSP fxAssetMO(params->market->GetData(
                                                        &model,
                                                        fxAssetName,
                                                         FXAsset::TYPE));
                    IObject* obj = dynamic_cast<IObject*>(fxAssetMO.get());
                    FXAsset* fxAsset = dynamic_cast<FXAsset*>(obj);
                    spot = fxAsset->getSpot();
                    (*fxRates)[j][i] = spot;
                    (*fxRates)[i][j] = Maths::isZero(spot)?0.0:(1.0/spot);

                } else if ( params->market->hasFXData((*fxNames)[i], 
                                                      (*fxNames)[j])) {
                    fxAssetName = params->market->getFXName((*fxNames)[i],
                                                            (*fxNames)[j]);
                    MarketObjectSP fxAssetMO(params->market->GetData(
                                                        &model,
                                                        fxAssetName,
                                                        FXAsset::TYPE));
                    IObject* obj = dynamic_cast<IObject*>(fxAssetMO.get());
                    FXAsset* fxAsset = dynamic_cast<FXAsset*>(obj);
                    spot = fxAsset->getSpot();
                    (*fxRates)[i][j] = Maths::isZero(spot)?0.0:(1.0/spot);
                    (*fxRates)[j][i] = spot;

                } else {
                    // couldn't find corresponding fx asset
                    (*fxRates)[i][j] = 1.0;
                    (*fxRates)[j][i] = 1.0;
                }
            }
        }

        return fxRates;
    }

    /** for reflection */
    XLGetFXRatesAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetFXRatesAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetFXRatesAddin);
        FIELD(fxNames, "currency names");
        FIELD(market, "Market handle");
        Addin::registerClassObjectMethod("GET_FX_RATES",
                                         Addin::MARKET,
                                         "retrieves the fx rates matrix for a"
                                         " set of currencies",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getFXMatrix);
    }

    static IObject* defaultGetFXRatesAddin() {
        return new XLGetFXRatesAddin();
    }
};

CClassConstSP const XLGetFXRatesAddin::TYPE = CClass::registerClassLoadMethod(
    "XLGetFXRatesAddin", typeid(XLGetFXRatesAddin), XLGetFXRatesAddin::load);

/** Addin to determine run-time type of a handle */
class XLGetFXVolsAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    CStringArraySP  fxNames;
    CMarketDataSP   market;

    /** the 'addin function' - construct object from components */
    static IObjectSP getFXVolMatrix(XLGetFXVolsAddin* params){
        static const string routine = "XLGetFXVolsAddin::getFXVolMatrix";

        CStringArraySP  fxNames = params->fxNames;
        string           fxAssetName;
        CClosedFormLN model("FlatFXVol");
        string fxVolName, fxVolNameReverse;
        double fxVol;
        

        int numFXAssets = params->fxNames->size();
        CDoubleMatrixSP  fxVols(new DoubleMatrix(numFXAssets, numFXAssets));

        for (int i = 0; i < numFXAssets; i++) {
            (*fxVols)[i][i] = 0.0;
            for (int j = i + 1; j < numFXAssets; j++) {
                fxVolName = (*fxNames)[i] + "_" + (*fxNames)[j];
                fxVolNameReverse = (*fxNames)[j] + "_" + (*fxNames)[i];
                if ( params->market->hasFXData((*fxNames)[i], (*fxNames)[j])) {
                    fxAssetName = params->market->getFXName((*fxNames)[i],
                                                            (*fxNames)[j]);
                    MarketObjectSP fxAssetMO(params->market->GetData(&model,
                                                                     fxAssetName,
                                                                     FXAsset::TYPE));
                    IObject* obj = dynamic_cast<IObject*>(fxAssetMO.get());
                    FXAsset* fxAsset = dynamic_cast<FXAsset*>(obj);
                    CVolBaseWrapper fxVolWrapper(fxAsset->getVolName());
                    fxVolWrapper.getData(&model, params->market);

                    ATMVolRequestSP   fxVolRequest(new ATMVolRequest());
                    CVolProcessedSP   fxProcVol(fxVolWrapper->getProcessedVol(fxVolRequest.get(), NULL));
                    CVolProcessedBSSP fxVolBS = CVolProcessedBSSP::dynamicCast(fxProcVol);
                    fxVol = fxVolBS->CalcVol(DateTime(0,0), DateTime(0,0));
                    (*fxVols)[i][j] = fxVol;
                    (*fxVols)[j][i] = fxVol;
                }// if we have a ProtEquity with no conventional cache their will be an FXVol but no FXAssets in the addtional cache
                else if (params->market->hasData(fxVolName, FXVolBase::TYPE))
                {
                    CVolBaseWrapper fxVolWrapper(fxVolName);
                    fxVolWrapper.getData(&model, params->market);

                    ATMVolRequestSP   fxVolRequest(new ATMVolRequest());
                    CVolProcessedSP   fxProcVol(fxVolWrapper->getProcessedVol(fxVolRequest.get(), NULL));
                    CVolProcessedBSSP fxVolBS = CVolProcessedBSSP::dynamicCast(fxProcVol);
                    fxVol = fxVolBS->CalcVol(DateTime(0,0), DateTime(0,0));
                    (*fxVols)[i][j] = fxVol;
                    (*fxVols)[j][i] = fxVol;
                }
                else if (params->market->hasData(fxVolNameReverse, FXVolBase::TYPE))
                {
                    CVolBaseWrapper fxVolWrapper(fxVolNameReverse);
                    fxVolWrapper.getData(&model, params->market);

                    ATMVolRequestSP   fxVolRequest(new ATMVolRequest());
                    CVolProcessedSP   fxProcVol(fxVolWrapper->getProcessedVol(fxVolRequest.get(), NULL));
                    CVolProcessedBSSP fxVolBS = CVolProcessedBSSP::dynamicCast(fxProcVol);
                    fxVol = fxVolBS->CalcVol(DateTime(0,0), DateTime(0,0));
                    (*fxVols)[i][j] = fxVol;
                    (*fxVols)[j][i] = fxVol;
                }
                else if ( params->market->hasFXData((*fxNames)[j], 
                                                      (*fxNames)[i])) 
                {
                    fxAssetName = params->market->getFXName((*fxNames)[j],
                                                            (*fxNames)[i]);
                    MarketObjectSP fxAssetMO(params->market->GetData(&model,
                                                                     fxAssetName,
                                                                     FXAsset::TYPE));
                    IObject* obj = dynamic_cast<IObject*>(fxAssetMO.get());
                    FXAsset* fxAsset = dynamic_cast<FXAsset*>(obj);
                    CVolBaseWrapper fxVolWrapper(fxAsset->getVolName());
                    fxVolWrapper.getData(&model, params->market);

                    ATMVolRequestSP   fxVolRequest(new ATMVolRequest());
                    CVolProcessedSP   fxProcVol(fxVolWrapper->getProcessedVol(fxVolRequest.get(), NULL));
                    CVolProcessedBSSP fxVolBS = CVolProcessedBSSP::dynamicCast(fxProcVol);
                    fxVol = fxVolBS->CalcVol(DateTime(0,0), DateTime(0,0));
                    (*fxVols)[i][j] = fxVol;
                    (*fxVols)[j][i] = fxVol;
                } else {
                    // couldn't find corresponding fx asset
                    (*fxVols)[i][j] = 0.0;
                    (*fxVols)[j][i] = 0.0;
                }
            }
        }

        return fxVols;
    }

    /** for reflection */
    XLGetFXVolsAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetFXVolsAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetFXVolsAddin);
        FIELD(fxNames, "currency names");
        FIELD(market, "Market handle");
        Addin::registerClassObjectMethod("GET_FX_VOLS",
                                         Addin::MARKET,
                                         "retrieves the fx vol matrix for a"
                                         " set of currencies",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getFXVolMatrix);
    }

    static IObject* defaultGetFXVolsAddin() {
        return new XLGetFXVolsAddin();
    }
};

CClassConstSP const XLGetFXVolsAddin::TYPE = CClass::registerClassLoadMethod(
    "XLGetFXVolsAddin", typeid(XLGetFXVolsAddin), XLGetFXVolsAddin::load);

// addin support
class GetRelatedObjectsAddin: public CObject{
    static CClassConstSP const TYPE;
    
    CMarketDataSP    marketDataCache;
    IObjectSP        object;
    
    /** create a market data cache */
    static IObjectSP getObjects(GetRelatedObjectsAddin* params){
        static const string routine("GetRelatedObjectsAddin::getObjects");

        typedef multimap<const string, MarketObjectSP> CDataMap;

        MarketObjectArraySP marketObjArray(new MarketObjectArray(0));

        if ( SimpleEquity::TYPE->isInstance(params->object) ||
             ProtEquity::TYPE->isInstance(params->object) ||
             StruckEquity::TYPE->isInstance(params->object) )
        {
            if ( SimpleEquity::TYPE->isInstance(params->object))
            {
                // cast to simple equity object
                IObject* obj         = dynamic_cast<IObject*>(params->object.get());
                SimpleEquity* equity = dynamic_cast<SimpleEquity*>(obj);
                            
            
                marketObjArray = params->marketDataCache->getEquityFXCorrelations(equity->getName(),
                                                                              FXAsset::TYPE,
                                                                              Correlation::TYPE);

                // get all vols for asset with same volName
                string volName  = equity->getVolName();
                
                try
                {
                    MarketObjectArraySP vols = params->marketDataCache->GetAllDataWithName(volName);
                    // The vol we are really interested in might actually be in the equity
                    MarketObjectSP vol = equity->getVol().getMO();

                    if ( vols->size() == 0 && !vol) 
                    {
                        // It's not in the cache or the Equity - oh dear!
                        throw ModelException(routine,
                                             "Vol with name " + volName +
                                             " is not in the cache or in the equity asset");
                    }
                    else
                    {
                        for (int i = 0; i < vols->size(); i++)
                        {
                            marketObjArray->push_back((*vols)[i]);
                        }
                        if (!(!vol))
                        {
                            marketObjArray->push_back(vol);
                        }
                    }
                }
                catch (exception& e)
                {
                    throw ModelException(&e, routine);
                }
                
            }
            else if ( ProtEquity::TYPE->isInstance(params->object))
            {
                // cast to simple equity object
                IObject* obj         = dynamic_cast<IObject*>(params->object.get());
                ProtEquity* equity = dynamic_cast<ProtEquity*>(obj);
                            
                if (!equity->useCorrFromCache())
                {
                    MarketObjectSP corr(copy(equity->getCorrelation()));
                    marketObjArray->push_back(corr);
                    params->marketDataCache->AddData(corr);
                }
                else
                {
                    marketObjArray = params->marketDataCache->getEquityFXCorrelations(equity->getName(),
                                                                                      FXAsset::TYPE,
                                                                                      Correlation::TYPE);
                }

                FXVolBaseWrapper fxvolWrapper = equity->getFXVol();
                if (!fxvolWrapper.usingCache())
                {
                    MarketObjectSP fxvol = fxvolWrapper.getMO();
                    marketObjArray->push_back(fxvol);
                    params->marketDataCache->AddData(fxvol);
                }

                // get all vols for asset with same volName
                string volName  = equity->getVolName();
                
                try
                {
                    MarketObjectArraySP vols = params->marketDataCache->GetAllDataWithName(volName);
                    // The vol we are really interested in might actually be in the equity
                    MarketObjectSP vol = equity->getVol().getMO();

                    if ( vols->size() == 0 && !vol) 
                    {
                        // It's not in the cache or the Equity - oh dear!
                        throw ModelException(routine,
                                             "Vol with name " + volName +
                                             " is not in the cache or in the equity asset");
                    }
                    else
                    {
                        for (int i = 0; i < vols->size(); i++)
                        {
                            marketObjArray->push_back((*vols)[i]);
                        }
                        if (!(!vol))
                        {
                            marketObjArray->push_back(vol);
                                }
                    }
                }
                catch (exception& e)
                {
                    throw ModelException(&e, routine);
                }
                
                
            }
            else if ( StruckEquity::TYPE->isInstance(params->object))
            {
                // cast to simple equity object
                IObject* obj         = dynamic_cast<IObject*>(params->object.get());
                StruckEquity* equity = dynamic_cast<StruckEquity*>(obj);
                            
                
                if (!equity->useCorrFromCache())
                {
                    MarketObjectSP corr(copy(equity->getCorrelation()));
                    marketObjArray->push_back(corr);
                    params->marketDataCache->AddData(corr);
                }
                else
                {
                    marketObjArray = params->marketDataCache->getEquityFXCorrelations(equity->getName(),
                                                                                      FXAsset::TYPE,
                                                                                      Correlation::TYPE);
                }

                FXAssetWrapper fxAssetWrapper = equity->getFXAsset();
                if (!fxAssetWrapper.usingCache())
                {
                    MarketObjectSP fxAsset = fxAssetWrapper.getMO();
                    marketObjArray->push_back(fxAsset);
                    params->marketDataCache->AddData(fxAsset);
                    IObject* iobj         = dynamic_cast<IObject*>(fxAsset.get());
                    FXAsset* fxa = dynamic_cast<FXAsset*>(iobj);
                    FXVolBaseWrapper fxVolWrapper = fxa->getFXVol();
                    MarketObjectSP fxaMO = fxVolWrapper.getMO();
                    marketObjArray->push_back(fxaMO);
                    params->marketDataCache->AddData(fxaMO);
                }

                // get all vols for asset with same volName
                string volName  = equity->getVolName();
                
                try
                {
                    MarketObjectArraySP vols = params->marketDataCache->GetAllDataWithName(volName);
                    // The vol we are really interested in might actually be in the equity
                    MarketObjectSP vol = equity->getVol().getMO();

                    if ( vols->size() == 0 && !vol) 
                    {
                        // It's not in the cache or the Equity - oh dear!
                        throw ModelException(routine,
                                             "Vol with name " + volName +
                                             " is not in the cache or in the equity asset");
                    }
                    else
                    {
                        for (int i = 0; i < vols->size(); i++)
                        {
                            marketObjArray->push_back((*vols)[i]);
                        }
                        if (!(!vol))
                        {
                            marketObjArray->push_back(vol);
                        }
                    }
                }
                catch (exception& e)
                {
                    throw ModelException(&e, routine);
                }
            }
        }
        else if ( Fund::TYPE->isInstance(params->object) ) 
        {
            // cast to simple equity object
            IObject* obj  = dynamic_cast<IObject*>(params->object.get());
            Fund*    fund = dynamic_cast<Fund*>(obj);
                
            marketObjArray = params->marketDataCache->getEquityFXCorrelations(
                fund->getName(),
                FXAsset::TYPE,
                Correlation::TYPE);
        } 
        else 
        {
            throw ModelException("Expecting object of type SimpleEquity, ProtEquity, StruckEquity or Fund - input handle\n"
                                 "is of type " + params->object->getClass()->getName(), 
                                 routine);
        }
        
        return marketObjArray;
    }
    
    GetRelatedObjectsAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetRelatedObjectsAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetRelatedObjectsAddin);
        FIELD(object,          "Handle to SimpleEquity object");
        FIELD(marketDataCache, "Handle to Market Data Cache");
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "GET_MARKET_DATA_FOR_EQUITY",
            Addin::MARKET,
            "Retrieves all market data objects related to an equity asset",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)getObjects);
    }

    static IObject* defaultGetRelatedObjectsAddin(){
        return new GetRelatedObjectsAddin();
    }
};

CClassConstSP const GetRelatedObjectsAddin::TYPE= CClass::registerClassLoadMethod(
    "GetRelatedObjectsAddin", typeid(GetRelatedObjectsAddin), GetRelatedObjectsAddin::load);

class GetObjectFromMarket: public CObject {
    static CClassConstSP const TYPE;

    CMarketDataSP    marketDataCache;
    string           objectName;
    string           objectType;
    IModelSP         model;
    bool             getInternalMarket;  // i.e populate any wrappers within this object

    /** create a market data cache */
    static IObjectSP getObject(GetObjectFromMarket* params){
        static const string routine("GetObjectFromMarket::getObject");

        CClassConstSP clazz = CClass::forName(params->objectType);

        if (!clazz) {
            throw ModelException("Invalid object type " + params->objectType, routine);
        }

        IModel* model = params->model.get();
        MarketObjectSP marketObj;
        if (model)
        {
            marketObj = params->marketDataCache->GetData(params->model.get(),
                                                                        params->objectName,
                                                                        clazz);
        }
        else
        {
            marketObj = params->marketDataCache->GetData(params->objectName,
                                                                        clazz);        
        }

        if (params->getInternalMarket)
        {
            IModelSP cfln(new CClosedFormLN("VolPreferred"));
            IModelSP mdl = params->model.get() ? params->model : cfln;
            marketObj->getMarket(mdl.get(), params->marketDataCache.get());
        }
        return marketObj;
    }
 
    GetObjectFromMarket(): CObject(TYPE), getInternalMarket(true){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GetObjectFromMarket, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetObjectFromMarket);
        FIELD(marketDataCache,      "Handle to Market Data Cache");
        FIELD(objectName,    "name of market data object");
        FIELD(objectType,    "type of market data object");
        FIELD(model,                "(optional) model");
        FIELD_MAKE_OPTIONAL(model);
        FIELD(getInternalMarket, "get internal market?");
         FIELD_MAKE_OPTIONAL(getInternalMarket);

        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "GET_OBJECT_FROM_MARKET",
            Addin::MARKET,
            "Retrieves an object from the market",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)getObject);
    }

    static IObject* defaultGetObjectFromMarket() {
        return new GetObjectFromMarket();
    }
};

CClassConstSP const GetObjectFromMarket::TYPE= CClass::registerClassLoadMethod(
    "GetObjectFromMarket", typeid(GetObjectFromMarket), GetObjectFromMarket::load);


class GetPVFromMarket: public CObject {
    static CClassConstSP const TYPE;

    CMarketDataSP    marketDataCache;
    string           objectName;
    MarketObjectSP   marketObj; // $unregistered
    DateTimeArraySP dates;     

    /** create a market data cache */
    static IObjectSP getObject(GetPVFromMarket* params){
        static const string routine("GetPVFromMarket::getObject");

        DoubleArraySP output(new DoubleArray(params->dates->size()));
        
        CClassConstSP clazz = CClass::forName("YieldCurve");
        if (!clazz) {
            throw ModelException("Invalid object type", routine);
        }

        MarketObjectSP marketObj = params->marketDataCache->GetData(params->objectName,
                                                                    clazz);
        if ( YieldCurve::TYPE->isInstance(marketObj.get()))
        {
            YieldCurve* yc = dynamic_cast<YieldCurve*>(marketObj.get());
            CClosedFormLN model("VolSurface");
            yc->getMarket(&model, params->marketDataCache.get());
            for (int i = 0; i < params->dates->size(); i++) {
                (*output)[i] = yc->pv((*params->dates)[i]);
           }
        }
       
        return IObjectSP(output);
    }
 
    GetPVFromMarket(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GetPVFromMarket, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetPVFromMarket);
        FIELD(marketDataCache,      "Handle to Market Data Cache");
        FIELD(objectName,    "name of market datat object");
        FIELD(dates, "Date(s) for pv factor");
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "PV_FACTOR_FROM_MARKET",
            Addin::MARKET,
            "Returns pv factor given yc name and market",
            TYPE,
            false,
            Addin::expandSimple,
            (Addin::ObjMethod*)getObject);
    }

    static IObject* defaultGetPVFromMarket() {
        return new GetPVFromMarket();
    }
};

CClassConstSP const GetPVFromMarket::TYPE= CClass::registerClassLoadMethod(
    "GetPVFromMarket", typeid(GetPVFromMarket), GetPVFromMarket::load);

class GetLastSensDateAddin: public CObject {
    static CClassConstSP const TYPE;
    
    IModelSP      model;
    CInstrumentSP inst;
    CMarketDataSP market;

    static IObjectSP getLastSensDateAddinFunc(GetLastSensDateAddin* params){
        static const string method("GetLastSensDateAddin::getObject");

        DateTimeSP output(0);
        if ( LastSensDate::TYPE->isInstance(params->inst.get()))
        {
            CInstrumentSP inst = CInstrumentSP(
                dynamic_cast<CInstrument*>(params->inst->clone())); // copy so we don't modify orig when we get market
            // get market data; some instruments require before end date (e.g Vanilla)
            params->model->getInstrumentAndModelMarket(params->market.get(), inst.get());
            LastSensDate* lsd = dynamic_cast<LastSensDate*>(inst.get());
            RhoParallel aShift(0.0001); 
            DateTime tmp = lsd->endDate(&aShift);
            output = DateTimeSP(new DateTime(tmp));
        }
        else
        {
            throw ModelException(method, "Only instrument types supporting lastSensDate are supported by this function!");
        }
        return IObjectSP(output);
    }
 
    GetLastSensDateAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GetLastSensDateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetLastSensDateAddin);
        FIELD(model, "model");
        FIELD(inst, "instrument");
        FIELD(market, "market");
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "CREDIT_LAST_EXPOSURE_DATE",
            Addin::UTILITIES,
            "computes last sensitivity date of an instrument",
            TYPE,
            false,
            Addin::expandSimple,
            (Addin::ObjMethod*)getLastSensDateAddinFunc);
    }

    static IObject* defaultGetLastSensDateAddin() {
        return new GetLastSensDateAddin();
    }
};

CClassConstSP const GetLastSensDateAddin::TYPE= CClass::registerClassLoadMethod(
    "GetLastSensDateAddin", typeid(GetLastSensDateAddin), GetLastSensDateAddin::load);

class GetFwdFromMarket: public CObject,
                        public ClientRunnable, 
                        public IRegressionTest  {
    static CClassConstSP const TYPE;

    CMarketDataSP    marketDataCache;
    string           objectName; // $unregistered
    MarketObjectSP   marketObj; // $unregistered
    DateTimeArraySP  dates;     
    string            ccyTreatment;
    YieldCurveWrapper discountCcy;
    mutable CAssetWrapper asset; // mutable so that it can be populated with the market data

    /** create a market data cache */
    static IObjectSP getFwd(const GetFwdFromMarket* params){
        static const string method("GetFwdFromMarket::getfwd");

        DoubleArraySP output(new DoubleArray(params->dates->size()));
        
        if (!(params->ccyTreatment == "N" || params->ccyTreatment == "S" || params->ccyTreatment == "P" )) {
            throw ModelException(method, "only ccyTreatment types N, S, and P currently supported");
        }
         
        if (params->ccyTreatment == "S" || params->ccyTreatment == "P") {
            if (params->discountCcy.getName() == "") {
                throw ModelException(method, "discountCcy must be set when ccyTreatment is S or P");
            }
        }

        CClosedFormLN model("VolSurface");
        CAsset::getAssetMarketData(&model, params->marketDataCache.get(), params->ccyTreatment, 
                                   params->discountCcy, params->asset);
    
        for (int i = 0; i < params->dates->size(); i++) {
            (*output)[i] = params->asset->fwdValue((*params->dates)[i]);
        }
        
       
        return IObjectSP(output);
    }

    IObjectSP run()
    {
        return (IObjectSP) getFwd(this);
    }

    /** for regression run */
    IObjectSP runTest() const
    {
        return (IObjectSP) getFwd(this);
    }
 
    GetFwdFromMarket(): CObject(TYPE){
        ccyTreatment = "N";}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GetFwdFromMarket, clazz);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetFwdFromMarket);
        FIELD(marketDataCache,      "Handle to Market Data Cache");
        FIELD(asset,    "name of equity");
        FIELD(dates, "Date(s) for Fwd factor");
        FIELD(ccyTreatment,    "N(one), S(truck) or P(rotected)");
        FIELD_MAKE_OPTIONAL(ccyTreatment);
        FIELD(discountCcy,    "discount ccy when struck or protected");
        FIELD_MAKE_OPTIONAL(discountCcy);
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "FORWARD_FROM_MARKET",
            Addin::MARKET,
            "Returns equity forward given stock name and market",
            TYPE,
            false,
            Addin::expandSimple,
            (Addin::ObjMethod*)getFwd);
    }

    static IObject* defaultGetFwdFromMarket() {
        return new GetFwdFromMarket();
    }
};

CClassConstSP const GetFwdFromMarket::TYPE= CClass::registerClassLoadMethod(
    "GetFwdFromMarket", typeid(GetFwdFromMarket), GetFwdFromMarket::load);


class BusDaysFwdFromMarket: public CObject {
    static CClassConstSP const TYPE;

    CMarketDataSP     marketDataCache;
    HolidayWrapper    hols;
    DateTimeArraySP   dates;
    CIntArraySP       numBusDays;

    static IObjectSP busDaysFwd(BusDaysFwdFromMarket* params){
        if (params->dates->size() != params->numBusDays->size()){
            throw ModelException("BusDaysFwdFromMarket", "Both arrays"
                                 " must be the same length");
        }

        CClosedFormLN model("VolSurface");
        params->hols.getData(&model, params->marketDataCache);

        DateTimeArraySP results(new DateTimeArray(params->dates->size()));
        for (int i = 0; i < results->size(); i++){
            (*results)[i] = 
                params->hols->addBusinessDays((*params->dates)[i], 
                                              (*params->numBusDays)[i]);
        }
        return results;
    }

    BusDaysFwdFromMarket(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BusDaysFwdFromMarket, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBusDaysFwdFromMarket);
        FIELD(marketDataCache,      "Handle to Market Data Cache");
        FIELD(hols, "Holiday name");
        FIELD(dates, "dates");
        FIELD(numBusDays, "number of business days to add");
        Addin::registerInstanceObjectMethod("BUS_DAYS_FWD_FROM_MARKET",
                                            Addin::MARKET,
                                            "Calculates new dates using offset"
                                            " of given number of business "
                                            "dates",
                                            TYPE,
                                            false,
                                            Addin::expandSimple,
                                            (Addin::ObjMethod*)busDaysFwd);
    }
    
    static IObject* defaultBusDaysFwdFromMarket(){
        return new BusDaysFwdFromMarket();
    }
};

CClassConstSP const BusDaysFwdFromMarket::TYPE= CClass::registerClassLoadMethod(
    "BusDaysFwdFromMarket", typeid(BusDaysFwdFromMarket), BusDaysFwdFromMarket::load);



class GetObjectPairFromMarket: public CObject {
    static CClassConstSP const TYPE;

    CMarketDataSP    marketDataCache;
    string           objectType;
    string           objectName1;
    string           objectName2;
    MarketObjectSP   marketObj; // $unregistered

    /** create a market data cache */
    static IObjectSP getObject(GetObjectPairFromMarket* params) {
        static const string routine("GetObjectPairFromMarket::getObject");

        CClassConstSP clazz = CClass::forName(params->objectType);

        if (!clazz) {
            throw ModelException("Invalid object type " + params->objectType, routine);
        }

        MarketObjectSP marketObj;
        CClosedFormLN model("VolSurface");

        if ( params->objectType == "Correlation" ) {
            // get correlation data from the market
            string correlationName = params->marketDataCache->getCorrelationName(
                                                   params->objectName1,params->objectName2);

            // then get hold of the actual object
            marketObj = params->marketDataCache->GetData(&model, correlationName, Correlation::TYPE);
        } else if ( params->objectType == "FXAsset" ) {
            // get an fx asset from the market
            if ( params->marketDataCache->hasFXData(params->objectName1, 
                                                    params->objectName2)) {
                string fxAssetName = params->marketDataCache->getFXName(params->objectName1, 
                                                                        params->objectName2);
                marketObj = params->marketDataCache->GetData(&model,fxAssetName, FXAsset::TYPE);
            }
        } else {
            throw ModelException("Invalid object type " + params->objectType, routine);
        }

        return marketObj;
    }
 
    GetObjectPairFromMarket(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GetObjectPairFromMarket, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetObjectPairFromMarket);
        FIELD(marketDataCache,      "Handle to Market Data Cache");
        FIELD(objectName1,   "name (risk ccy for FXAsset)");
        FIELD(objectName2,   "name (base ccy for FXAsset)");
        FIELD(objectType,    "type of market datat object");

        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "GET_OBJECT_PAIR_FROM_MARKET",
            Addin::MARKET,
            "Retrieves an object from the market given a pair of names",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)getObject);
    }

    static IObject* defaultGetObjectPairFromMarket() {
        return new GetObjectPairFromMarket();
    }
};

CClassConstSP const GetObjectPairFromMarket::TYPE= CClass::registerClassLoadMethod(
    "GetObjectPairFromMarket", typeid(GetObjectPairFromMarket), GetObjectPairFromMarket::load);

class GetEquityNames: public CObject {
    static CClassConstSP const TYPE;

    CInstrumentSP    instrument;
    CMarketDataSP    marketDataCache;

    /** create a market data cache */
    static IObjectSP getNames(GetEquityNames* params){
        static const string routine("GetEquityNames::getNames");

        int i,j;
        bool found;

        DeltaSP         delta(new Delta(0.005, NULL, NULL));

        CClosedFormLN   model("VolSurface");

        CInstrumentSP   inst(copy(params->instrument.get()));

        model.getInstrumentAndModelMarket(params->marketDataCache.get(), inst.get());

        // get list of equity names
        OutputNameArrayConstSP names(delta->names(inst.get()));

        CStringArraySP equityNames(new StringArray(0));
        for(i=0;i<names->size();++i) {
           equityNames->push_back((*names)[i]->toString());
        }

        try {
            DeltaToCreditSP deltaToCredit(new DeltaToCredit(0.005, NULL, NULL));
            AssetVegaParallelSP assetVega(new AssetVegaParallel(0.005));

            ClosedFormFA   firmAssetModel;
            firmAssetModel.getInstrumentAndModelMarket(params->marketDataCache.get(), inst.get());
            OutputNameArrayConstSP namesE2C(deltaToCredit->names(inst.get()));
            OutputNameArrayConstSP assetNames(assetVega->names(inst.get()));

            // append deltaToCredit names 
            for(i=0;i<namesE2C->size();++i) {
                found = false;
                for(j=0;j<equityNames->size();++j) {
                    if ( (*namesE2C)[i]->toString() == (*equityNames)[j] ) {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    equityNames->push_back((*namesE2C)[i]->toString());
                }
            }

            // append firma asset names 
            for(i=0;i<assetNames->size();++i) {
                found = false;
                for(j=0;j<equityNames->size();++j) {
                    if ( (*assetNames)[i]->toString() == (*equityNames)[j] ) {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    equityNames->push_back((*assetNames)[i]->toString());
                }
            }
        } catch(exception&) {
        }

        return equityNames;
    }
 
    GetEquityNames(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GetEquityNames, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetEquityNames);
        FIELD(instrument,      "Handle to the instrument");
        FIELD(marketDataCache, "Handle to the Market Data");

        // register addin for adding objects to cache
        Addin::registerClassObjectMethod("GET_EQUITY_NAMES",
                                         Addin::UTILITIES,
                                         "Retrieves the names of all equities",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)getNames);
    }

    static IObject* defaultGetEquityNames() {
        return new GetEquityNames();
    }
};

CClassConstSP const GetEquityNames::TYPE= CClass::registerClassLoadMethod(
    "GetEquityNames", typeid(GetEquityNames), GetEquityNames::load);

/** Addin to determine run-time type of a handle */
class XLGetAssetCurrencyAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    CMarketDataSP    market;
    string           assetName;

    /** the 'addin function' - construct object from components */
    static string getAssetCcy(XLGetAssetCurrencyAddin* params) {
        static const string routine = "XLGetAssetNameAddin::getAssetCcy";

        CClosedFormLN model("VolSurface");
        MarketObjectSP assetMO(params->market->GetData(&model,
                                                       params->assetName,
                                                       SimpleEquity::TYPE));
        if (!(!assetMO)) {
            IObject*      obj     = dynamic_cast<IObject*>(assetMO.get());
            SimpleEquity* equity  = dynamic_cast<SimpleEquity*>(obj);

            return (equity->getYCName());
        } else {
            return ("");
        }
    }

    /** for reflection */
    XLGetAssetCurrencyAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetAssetCurrencyAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLGetAssetCurrencyAddin);
        FIELD(assetName, "The asset's name");
        FIELD(market, "Handle to the Market Data");
        Addin::registerClassStringMethod("GET_ASSET_CURRENCY",
                                         Addin::MARKET,
                                         "Returns the currency in which the asset is denominated",
                                         TYPE,
                                         (Addin::StringMethod*)getAssetCcy);
    }

    static IObject* defaultXLGetAssetCurrencyAddin() {
        return new XLGetAssetCurrencyAddin();
    }
};

CClassConstSP const XLGetAssetCurrencyAddin::TYPE= CClass::registerClassLoadMethod(
    "XLGetAssetCurrencyAddin", typeid(XLGetAssetCurrencyAddin), XLGetAssetCurrencyAddin::load);


/** Addin to determine run-time type of a handle */
class XLGetSpotFXAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    string           baseCcy;
    string           riskCcy;
    CMarketDataSP    market;

    /** the 'addin function' - construct object from components */
    static double getSpotFX(XLGetSpotFXAddin* params) {
        static const string routine = "XLGetSpotFXAddin::getSpotFX";

        double spotFX = 1.0;

        CClosedFormLN model("FlatFXVol");
        string        baseCcy = params->baseCcy;
        string        riskCcy = params->riskCcy;

        if ( params->market->hasFXData(riskCcy, baseCcy)) {
            string fxAssetName = params->market->getFXName(riskCcy, baseCcy);
            MarketObjectSP fxAssetMO(params->market->GetData(
                                                &model,
                                                fxAssetName,
                                                 FXAsset::TYPE));
            IObject* obj = dynamic_cast<IObject*>(fxAssetMO.get());
            FXAsset* fxAsset = dynamic_cast<FXAsset*>(obj);

            spotFX = fxAsset->getSpot();
        } else {
            throw ModelException("Could not find fx asset with risk currency " + riskCcy + " \n and"
                                 " base currency " + baseCcy, routine);
        }
        return spotFX;
    }

    /** for reflection */
    XLGetSpotFXAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetSpotFXAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLGetSpotFXAddin);
        FIELD(baseCcy, "The fx rate's base currency");
        FIELD(riskCcy, "The fx rate's risk currency");
        FIELD(market, "Handle to the Market Data");
        Addin::registerClassDoubleMethod("GET_SPOT_FX",
                                         Addin::MARKET,
                                         "Returns the a spot fx rate",
                                         TYPE,
                                         (Addin::DoubleMethod*)getSpotFX);
    }

    static IObject* defaultXLGetSpotFXAddin() {
        return new XLGetSpotFXAddin();
    }
};

CClassConstSP const XLGetSpotFXAddin::TYPE= CClass::registerClassLoadMethod(
    "XLGetSpotFXAddin", typeid(XLGetSpotFXAddin), XLGetSpotFXAddin::load);

class GetCurrencyNames: public CObject {
    static CClassConstSP const TYPE;

    CInstrumentSP    instrument;
    CMarketDataSP    marketDataCache;

    /** create a market data cache */
    static IObjectSP getNames(GetCurrencyNames* params){
        static const string routine("GetCurrencyNames::getNames");

        // could possibly better use some kind of collector

        CClosedFormLN   model("VolSurface");
        CInstrumentSP   inst(copy(params->instrument.get()));

        model.getInstrumentAndModelMarket(params->marketDataCache.get(), inst.get());

        // get list of currency names
        OutputNameArrayConstSP names = RiskProperty<RateParallel>().subjectNames(inst);

        CStringArraySP currencyNames(new StringArray(0));
        for(int i=0;i<names->size();++i) {
            // check whether currency already exists
            bool foundCcy = false;
            for(int j=0;j<currencyNames->size();++j) {
                if ( (*names)[i]->toString() == (*currencyNames)[j]) {
                    foundCcy = true;
                }
            }
            if (!foundCcy) {
                currencyNames->push_back((*names)[i]->toString());
            }
        }
        return currencyNames;
    }
 
    GetCurrencyNames(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GetCurrencyNames, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetCurrencyNames);
        FIELD(instrument,      "Handle to the instrument");
        FIELD(marketDataCache, "Handle to the Market Data");

        // register addin for adding objects to cache
        Addin::registerClassObjectMethod("GET_CURRENCY_NAMES",
                                         Addin::UTILITIES,
                                         "Retrieves the names of all currencies",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)getNames);
    }

    static IObject* defaultGetCurrencyNames() {
        return new GetCurrencyNames();
    }
};

CClassConstSP const GetCurrencyNames::TYPE= CClass::registerClassLoadMethod(
    "GetCurrencyNames", typeid(GetCurrencyNames), GetCurrencyNames::load);

/** Addin in order to retrieve all asset names */
class GetAssetNames: public CObject {
    static CClassConstSP const TYPE;

    CInstrumentSP    instrument;
    CMarketDataSP    market;
    bool             viewAllCorrs;

       /** create a market data cache */
    static IObjectSP getNames(GetAssetNames* params){       
        static const string method("GetAssetNames::getNames");
        try {
            int i,j;

            /** dummy model */
            CClosedFormLN   model("VolSurface");
            CInstrumentSP   inst(copy(params->instrument.get()));
            model.getInstrumentAndModelMarket(params->market.get(), inst.get());
            CStringArraySP assetNames(new StringArray(0));
            CStringArraySP assetNamesHelper(new StringArray(0));

            /** only EQEQ -- hardly anything todo */
            if(!(params->viewAllCorrs)) {
                DeltaSP delta(new Delta(0.005, NULL, NULL));
                OutputNameArrayConstSP eqNames(delta->names(inst.get()));
                for (i = 0; i < eqNames->size(); i++) {
                    assetNames->push_back((*eqNames)[i]->toString());
                    assetNamesHelper->push_back("EQ");
                }
                ObjectArraySP result(new ObjectArray(2));
                (*result)[0] = assetNames;
                (*result)[1] = assetNamesHelper;
                return result;
            }

            /** retrieve base ccy */
            string baseCcy = inst->discountYieldCurveName(); // eg "currency-EUR-MM"
            MarketObjectSP mo(params->market->GetData(&model, baseCcy, YieldCurve::TYPE)); 
            IObject*      obj   = dynamic_cast<IObject*>(mo.get());
            YieldCurve*   yc    = dynamic_cast<YieldCurve*>(obj);
            string baseCcyName = yc->getCcy(); // eg "EUR"

            /** retrieve IR names, which are "EUR" rather than "currency-EUR-MM" */ 
            OutputNameArrayConstSP ccyNames = RiskProperty<RateParallel>().subjectNames(inst);
            for(i = 0; i < ccyNames->size(); i++) {
                bool found = false; // check whether currency already exists
                for(j = 0; j < assetNames->size(); j++) {
                    if ( (*ccyNames)[i]->toString() == (*assetNames)[j]) {
                        found = true;
                    }
                }
                if (!found) {
                    assetNames->push_back((*ccyNames)[i]->toString());
                    assetNamesHelper->push_back("IR");
                }
            }

            /** retrieve FX names, which are in basket */
            FXDeltaSP fxDelta(new FXDelta(0.005, NULL, NULL));
            OutputNameArrayConstSP fxNames(fxDelta->names(inst.get()));
            for(i = 0; i < fxNames->size(); i++) {
                string thisName = (*fxNames)[i]->toString();
                if (params->market->hasData(thisName, FXAsset::TYPE)) {
                    MarketObjectSP moFx(params->market->GetData(&model, thisName, FXAsset::TYPE)); 
                    IObject*    objFx  = dynamic_cast<IObject*>(moFx.get());
                    FXAsset*    forex  = dynamic_cast<FXAsset*>(objFx);                    
                    bool found = false; // check whether fx already exists
                    for(j = 0; j < assetNames->size(); j++) {
                        if ( (*fxNames)[i]->toString() == (*assetNames)[j]) {
                            found = true;
                        }
                    }
                    if (!found) {
                        assetNames->push_back((*fxNames)[i]->toString());
                        assetNamesHelper->push_back(baseCcyName + "-" + forex->getRiskCcyIsoCode());
                    }
                }
            }

            /** retrieve EQ names (put them only at the end into assetNames) */
            DeltaSP delta(new Delta(0.005, NULL, NULL));
            OutputNameArrayConstSP eqNames(delta->names(inst.get()));

            /** retrieve potentially further FX names (which are not in basket) */
            StringArraySP furtherCcys(new StringArray(0));
            StringArraySP furtherCcysNames(new StringArray(0));
            for(i = 0; i < eqNames->size(); i++) {
                string thisName = (*eqNames)[i]->toString();
                if (params->market->hasData(thisName, SimpleEquity::TYPE)) {
                    MarketObjectSP moEq(params->market->GetData(&model, thisName, SimpleEquity::TYPE)); 
                    IObject*      objEq   = dynamic_cast<IObject*>(moEq.get());
                    SimpleEquity* equity  = dynamic_cast<SimpleEquity*>(objEq);
                    string ccy = equity->getYCName();
                    MarketObjectSP moIr(params->market->GetData(&model, ccy, YieldCurve::TYPE)); 
                    IObject*    objIr   = dynamic_cast<IObject*>(moIr.get());
                    YieldCurve* yc      = dynamic_cast<YieldCurve*>(objIr);
                    string ccyName = yc->getCcy();
                    bool found = false;
                    for (j = 0; j < furtherCcys->size(); j++) {
                        if ( ccy == (*furtherCcys)[j] ) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        furtherCcys->push_back(ccy);
                        furtherCcysNames->push_back(ccyName);
                    }
                } 
            }

            /** now potentially further FX assets */
            for (i = 0; i < furtherCcys->size(); i++) {
                string riskCcy = (*furtherCcys)[i];
                if (! (baseCcy == riskCcy) ){
                    string fxAssetName = params->market->getFXName(riskCcy, baseCcy); 
                    bool found = false; // check whether fx already exists
                    for(j = 0; j < assetNames->size(); j++) {
                        if ( fxAssetName == (*assetNames)[j]) {
                            found = true;
                        }
                    }
                    if (!found) {
                        assetNames->push_back(fxAssetName);
                        assetNamesHelper->push_back(baseCcyName + "-" + (*furtherCcysNames)[i]);
                    }
                }
            }

            /** append EQ assets only at the end */
            for (i = 0; i < eqNames->size(); i++) {
                assetNames->push_back((*eqNames)[i]->toString());
                assetNamesHelper->push_back("EQ");
            }
            ObjectArraySP result(new ObjectArray(2));
            (*result)[0] = assetNames;
            (*result)[1] = assetNamesHelper;
            return result;
        } catch (exception& e) {
            throw ModelException(&e, method);
        }
    }
 
    GetAssetNames(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GetAssetNames, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetAssetNames);
        FIELD(instrument,      "Handle to the instrument");
        FIELD(market, "Handle to the Market Data");
        FIELD(viewAllCorrs, "whether or not all corrs");

        // register addin for adding objects to cache
        Addin::registerClassObjectMethod("GET_ASSET_NAMES",
                                         Addin::UTILITIES,
                                         "Retrieves the names of all currencies",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)getNames);
    }

    static IObject* defaultGetAssetNames() {
        return new GetAssetNames();
    }
};

CClassConstSP const GetAssetNames::TYPE= CClass::registerClassLoadMethod(
    "GetAssetNames", typeid(GetAssetNames), GetAssetNames::load);

/** Addin to determine run-time type of a handle */
class XLGetNFactorCorrelationMatrixAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    CStringArray     equityNames;
    CMarketDataSP    market;

    /** the 'addin function' - construct object from components */
    static IObjectSP getMatrix(XLGetNFactorCorrelationMatrixAddin* params){
        static const string routine = "XLGetNFactorCorrelationMatrixAddin::getMatrix";

        CDoubleMatrixSP  correlations;
        string           corrName;

        CClosedFormLN model("VolSurface");

        // generate double matrix of correlations
        // First allocate space
        int numAssets = params->equityNames.size();
        correlations = CDoubleMatrixSP(new DoubleMatrix(numAssets, numAssets));
        // then look up names and values (stored in DoubleMatrix)
        int pos = 0;
        for (int i = 0; i < numAssets; i++) {
            (*correlations)[i][i] = 1.0;
            for (int j = i + 1; j < numAssets; j++, pos++) {
                // look up name in cache
                corrName = params->market->getCorrelationName(params->equityNames[i],
                                                              params->equityNames[j]);
                // then get hold of the actual object
                MarketObjectSP corr(params->market->GetData(&model, corrName, 
                                                            Correlation::TYPE));

                IObject* obj = dynamic_cast<IObject*>(corr.get());
                Correlation& correlation =
                    dynamic_cast<Correlation&>(*obj);
                (*correlations)[i][j] = correlation.getCorrelation();
                (*correlations)[j][i] = (*correlations)[i][j];
            }
        }
        return correlations;
    }

    /** for reflection */
    XLGetNFactorCorrelationMatrixAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetNFactorCorrelationMatrixAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCorrelationMatrixAddin);
        FIELD(equityNames, "Equity underlyings");
        FIELD(market,      "Market handle");
        Addin::registerClassObjectMethod("GET_NFACTOR_CORRELATION_MATRIX",
                                         Addin::MARKET,
                                         "retrieves the correlation matrix for"
                                         " a set of equities from the market data cache",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getMatrix);
    }

    static IObject* defaultCorrelationMatrixAddin() {
        return new XLGetNFactorCorrelationMatrixAddin();
    }
};

CClassConstSP const XLGetNFactorCorrelationMatrixAddin::TYPE = CClass::registerClassLoadMethod(
    "XLGetNFactorCorrelationMatrixAddin", typeid(XLGetNFactorCorrelationMatrixAddin), 
    XLGetNFactorCorrelationMatrixAddin::load);


/** Addin to determine run-time type of a handle */
class GetMarketObjectName: public CObject,
                           virtual public ClientRunnable {
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    MarketObjectSP marketObject;

    // EdrAction version of addin
    IObjectSP run() {
        return IObjectSP(CString::create(getName(this)));
    }

    /** the 'addin function' - construct object from components */
    static string getName(GetMarketObjectName* params){
        return params->marketObject->getName();
    }

    /** for reflection */
    GetMarketObjectName():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(GetMarketObjectName, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultGetMarketObjectName);
        FIELD(marketObject, "handle to market object");
        Addin::registerClassStringMethod("GET_MARKET_OBJECT_NAME",
                                         Addin::MARKET,
                                         "Returns the name of a market object",
                                         TYPE,
                                         (Addin::StringMethod*)getName);
    }

    static IObject* defaultGetMarketObjectName() {
        return new GetMarketObjectName();
    }
};

CClassConstSP const GetMarketObjectName::TYPE = CClass::registerClassLoadMethod(
        "GetMarketObjectName", typeid(GetMarketObjectName), GetMarketObjectName::load);


class XLGetObjectFromWrapper: public CObject {
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP      marketWrapperObject;

    /** the 'addin function' - construct object from components */
    static IObjectSP getObject(XLGetObjectFromWrapper* params){
        static const string routine = "XLGetObjectFromWrapper::getObject";

        if ( MarketObjectWrapper::TYPE->isInstance(params->marketWrapperObject) ) {
            MarketObjectWrapper* mow = 
                    dynamic_cast<MarketObjectWrapper*>(params->marketWrapperObject.get());
            return mow->getMO();

        } else {
            throw ModelException("Expecting object of type MarketObjectWrapper\n", 
                                 routine);
        }

        return IObjectSP(NULL);
    }

    /** for reflection */
    XLGetObjectFromWrapper():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetObjectFromWrapper, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCorrelationMatrixAddin);
        FIELD(marketWrapperObject,      "Market wrapper object");
        
        Addin::registerClassObjectMethod("GET_OBJECT_FROM_WRAPPER",
                                         Addin::MARKET,
                                         "converts a wrapper object into a market object",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getObject);
    }

    static IObject* defaultCorrelationMatrixAddin() {
        return new XLGetObjectFromWrapper();
    }
};

CClassConstSP const XLGetObjectFromWrapper::TYPE = CClass::registerClassLoadMethod(
    "XLGetObjectFromWrapper", typeid(XLGetObjectFromWrapper), 
    XLGetObjectFromWrapper::load);

class XLUsingCacheAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP      marketWrapperObject;

    /** the 'addin function' - construct object from components */
    static bool usingCache(XLUsingCacheAddin* params){
        static const string routine = "XLUsingCacheAddin::usingCache";

        bool usingCache = false;

        if ( MarketObjectWrapper::TYPE->isInstance(params->marketWrapperObject) ) {
            MarketObjectWrapper* mow = 
                    dynamic_cast<MarketObjectWrapper*>(params->marketWrapperObject.get());
            return mow->usingCache();

        } else {
            throw ModelException("Expecting object of type MarketObjectWrapper\n", 
                                 routine);
        }
        return usingCache;
    }

    /** for reflection */
    XLUsingCacheAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLUsingCacheAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLUsingCacheAddin);
        FIELD(marketWrapperObject, "Handle to market object wrapper");
        Addin::registerClassBoolMethod("WRAPPER_IS_USING_CACHE",
                                       Addin::MARKET,
                                       "adds an object to the market data cache",
                                       TYPE,
                                       (Addin::BoolMethod*)usingCache);
    }

    static IObject* defaultXLUsingCacheAddin() {
        return new XLUsingCacheAddin();
    }
};


CClassConstSP const XLUsingCacheAddin::TYPE = CClass::registerClassLoadMethod(
    "XLUsingCacheAddin", typeid(XLUsingCacheAddin), XLUsingCacheAddin::load);


class XLGetMarket: public CObject {
    static CClassConstSP const TYPE;

    IObjectSP      fileHandle;

    /** the 'addin function' - construct object from components */
    static IObjectSP getMarket(XLGetMarket* params){
        static const string routine = "XLGetMarket::getObject";

        CMarketDataSP market;

        if ( RiskMgrInterface::TYPE->isInstance(params->fileHandle) ) {
            // risk mgr interface
            RiskMgrInterfaceSP riskMgr = RiskMgrInterfaceSP::dynamicCast(params->fileHandle);
            market = riskMgr->market;
        } else if (ScenarioInterface::TYPE->isInstance(params->fileHandle) ) {
            // scenario interface
            ScenarioInterfaceSP scenario = ScenarioInterfaceSP::dynamicCast(params->fileHandle);
            market = scenario->market;
        } else if ( CompositeInstrument::TYPE->isInstance(params->fileHandle) ) {
            // composite instruments
            CompositeInstrumentSP composite = CompositeInstrumentSP::dynamicCast(params->fileHandle);
            market = composite->market;
        } else {
            throw ModelException("Expecting object of type CompositeInstrument or RiskMgrInterface\n", 
                                 routine);
        }

        // return the market
        return market;
    }

    /** for reflection */
    XLGetMarket():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetMarket, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLGetMarket);
        FIELD(fileHandle,      "file handle");
        
        Addin::registerClassObjectMethod("GET_MARKET_DATA",
                                         Addin::MARKET,
                                         "retrieves the market data cache from a regression file handle",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getMarket);
    }

    static IObject* defaultXLGetMarket() {
        return new XLGetMarket();
    }
};

CClassConstSP const XLGetMarket::TYPE = CClass::registerClassLoadMethod(
    "XLGetMarket", typeid(XLGetMarket), XLGetMarket::load);

class XLGetWrappers: public CObject {
    static CClassConstSP const TYPE;

    IObjectSP      fileHandle;

    /** the 'addin function' - construct object from components */
    static IObjectSP getMarket(XLGetWrappers* params){
        static const string routine = "XLGetWrappers::getObject";

        IModelSP      model;
        CInstrumentSP inst;
        CControlSP    ctrl;
        CMarketDataSP market;


        if ( RiskMgrInterface::TYPE->isInstance(params->fileHandle) ) {
            // risk mgr interface
            RiskMgrInterfaceSP riskMgr = RiskMgrInterfaceSP::dynamicCast(params->fileHandle);
            market = riskMgr->market;
            model  = riskMgr->model;
            inst   = riskMgr->inst;
            ctrl   = riskMgr->ctrl;

        } else if (ScenarioInterface::TYPE->isInstance(params->fileHandle) ) {
            // scenario interface
            ScenarioInterfaceSP scenario = ScenarioInterfaceSP::dynamicCast(params->fileHandle);
            market = scenario->market;
            model  = scenario->model;
            inst   = scenario->inst;
            ctrl   = scenario->ctrl;

        } else if ( CompositeInstrument::TYPE->isInstance(params->fileHandle) ) {
            // composite instruments
            CompositeInstrumentSP composite = CompositeInstrumentSP::dynamicCast(params->fileHandle);
            market = composite->market;
            model  = IModelSP::dynamicCast((*composite->model)[0]);
            inst   = CInstrumentSP::dynamicCast((*composite->inst)[0]);
            ctrl   = CControlSP::dynamicCast((*composite->ctrl)[0]);

        } else {
            throw ModelException("Expecting object of type CompositeInstrument or RiskMgrInterface\n", 
                                 routine);
        }

        // create clones of the instrument and control
        CInstrumentSP cloneImnt(copy(inst.get()));
        CControlSP    cloneContrl(copy(ctrl.get()));

        if (!market) {
            market = CMarketDataSP(new MarketData());
        }

        model->getInstrumentAndModelMarket(market.get(), cloneImnt.get());
        cloneContrl->getMarket(
            model, market, IInstrumentCollection::singleton(cloneImnt));

        WrapperNameCollector wrapperNameColl;

        CStringArraySP  imntWrapperNames = wrapperNameColl.getWrapperNames(cloneImnt);
        CStringArraySP  ctrlWrapperNames = wrapperNameColl.getWrapperNames(cloneContrl);

        for (int i=0 ; i<ctrlWrapperNames->size() ; ++i) {
            imntWrapperNames->push_back((*ctrlWrapperNames)[i]);
        }

        // return the market
        return imntWrapperNames;
    }

    /** for reflection */
    XLGetWrappers():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetWrappers, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLGetWrappers);
        FIELD(fileHandle,      "file handle");
        
        Addin::registerClassObjectMethod("GET_WRAPPER_NAMES",
                                         Addin::MARKET,
                                         "retrieves all market object wrapper names",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)getMarket);
    }

    static IObject* defaultXLGetWrappers() {
        return new XLGetWrappers();
    }
};

CClassConstSP const XLGetWrappers::TYPE = CClass::registerClassLoadMethod(
    "XLGetWrappers", typeid(XLGetWrappers), XLGetWrappers::load);

class XLGetYieldCurveNames: public CObject {
    static CClassConstSP const TYPE;

    IObjectSP      fileHandle;

    /** the 'addin function' - construct object from components */
    static IObjectSP getYieldNames(XLGetYieldCurveNames* params){
        static const string routine = "XLGetYieldCurveNames::getYieldNames";

        IModelSP      model;
        CInstrumentSP inst;
        CControlSP    ctrl;
        CMarketDataSP market;


        if ( RiskMgrInterface::TYPE->isInstance(params->fileHandle) ) {
            // risk mgr interface
            RiskMgrInterfaceSP riskMgr = RiskMgrInterfaceSP::dynamicCast(params->fileHandle);
            market = riskMgr->market;
            model  = riskMgr->model;
            inst   = riskMgr->inst;
            ctrl   = riskMgr->ctrl;
            
        } else if (ScenarioInterface::TYPE->isInstance(params->fileHandle) ) {
            // scenario interface
            ScenarioInterfaceSP scenario = ScenarioInterfaceSP::dynamicCast(params->fileHandle);
            market = scenario->market;
            model  = scenario->model;
            inst   = scenario->inst;
            ctrl   = scenario->ctrl;

        } else if ( CompositeInstrument::TYPE->isInstance(params->fileHandle) ) {
            // composite instruments
            CompositeInstrumentSP composite = CompositeInstrumentSP::dynamicCast(params->fileHandle);
            market = composite->market;
            model  = IModelSP::dynamicCast((*composite->model)[0]);
            inst   = CInstrumentSP::dynamicCast((*composite->inst)[0]);
            ctrl   = CControlSP::dynamicCast((*composite->ctrl)[0]);

        } else {
            throw ModelException("Expecting object of type CompositeInstrument or RiskMgrInterface\n", 
                                 routine);
        }

        // create clones of the instrument and control
        CInstrumentSP cloneImnt(copy(inst.get()));
        CControlSP    cloneContrl(copy(ctrl.get()));

        if (!market) {
            market = CMarketDataSP(new MarketData());
        }
        else
        {
            //CMarketSP cloneMkt(copy(market.get()));
        }

        model->getInstrumentAndModelMarket(market.get(), cloneImnt.get());
        cloneContrl->getMarket(model, market,
                               IInstrumentCollection::singleton(cloneImnt));

        YieldNameCollector yieldNameColl;

        CStringArraySP  imntYieldNames = yieldNameColl.getYieldNames(cloneImnt);
        CStringArraySP  ctrlYieldNames = yieldNameColl.getYieldNames(cloneContrl);

        for (int i=0 ; i<ctrlYieldNames->size() ; ++i) {
            imntYieldNames->push_back((*ctrlYieldNames)[i]);
        }

        // sort and make unique
        list<string> yieldList;
        for (int j = 0; j < imntYieldNames->size(); j++)
        {
            yieldList.push_back((*imntYieldNames)[j]);
        }

        yieldList.sort();
        yieldList.unique();

        // copy them back to a string array
        CStringArraySP uniqueNames(new CStringArray(0));

        while (!yieldList.empty())
        {
            const string& theString = yieldList.front();
            uniqueNames->push_back(theString);
            yieldList.pop_front();
        }

        // return the unique yield names
        return uniqueNames;
    }

    /** for reflection */
    XLGetYieldCurveNames():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetYieldCurveNames, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLGetYieldCurveNames);
        FIELD(fileHandle,      "file handle");
        
        Addin::registerClassObjectMethod("GET_YIELD_NAMES",
                                         Addin::MARKET,
                                         "retrieves all yield curve names",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)getYieldNames);
    }

    static IObject* defaultXLGetYieldCurveNames() {
        return new XLGetYieldCurveNames();
    }
};

CClassConstSP const XLGetYieldCurveNames::TYPE = CClass::registerClassLoadMethod(
    "XLGetYieldCurveNames", typeid(XLGetYieldCurveNames), XLGetYieldCurveNames::load);


/** Addin to determine run-time type of a handle */
class FileInstTypeAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP      fileHandle;

    /** the 'addin function' - construct object from components */
    static string getTypeKey(FileInstTypeAddin* params){
        static const string routine = "FileInstTypeAddin::getTypeKey";

        CInstrumentSP inst;

        if ( RiskMgrInterface::TYPE->isInstance(params->fileHandle) ) {
            // risk mgr interface
            RiskMgrInterfaceSP riskMgr = RiskMgrInterfaceSP::dynamicCast(params->fileHandle);
            inst   = riskMgr->inst;

        } else if (ScenarioInterface::TYPE->isInstance(params->fileHandle) ) {
            // scenario interface
            ScenarioInterfaceSP scenario = ScenarioInterfaceSP::dynamicCast(params->fileHandle);
            inst   = scenario->inst;
            
        } else if ( CompositeInstrument::TYPE->isInstance(params->fileHandle) ) {
            // composite instruments
            CompositeInstrumentSP composite = CompositeInstrumentSP::dynamicCast(params->fileHandle);
            inst   = CInstrumentSP::dynamicCast((*composite->inst)[0]);
        } else {
            throw ModelException("Expecting object of type CompositeInstrument or RiskMgrInterface\n", 
                                 routine);
        }
        return inst->getClass()->getName();
    }

    /** for reflection */
    FileInstTypeAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(FileInstTypeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultTypeKeyAddin);
        FIELD(fileHandle, "file handle");
        Addin::registerClassStringMethod("GET_INSTRUMENT_TYPE",
                                         Addin::UTILITIES,
                                         "determines the run-time type of an "
                                         " object handle",
                                         TYPE,
                                         (Addin::StringMethod*)getTypeKey);
    }

    static IObject* defaultTypeKeyAddin() {
        return new FileInstTypeAddin();
    }
};

CClassConstSP const FileInstTypeAddin::TYPE = CClass::registerClassLoadMethod(
    "FileInstTypeAddin", typeid(FileInstTypeAddin), load);

class ExpandObjects: public CObject {
    static CClassConstSP const TYPE;

    CInstrumentSP    instrument;
    CMarketDataSP    marketDataCache;
    IModelSP         model;

    /** create a market data cache */
    static IObjectSP expand(ExpandObjects* params){
        static const string routine("ExpandObjects::expand");

        CInstrumentSP   inst(copy(params->instrument.get()));
        params->model->getInstrumentAndModelMarket(params->marketDataCache.get(), inst.get());

        return inst;
    }
 
    ExpandObjects(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(ExpandObjects, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultExpandObjects);
        FIELD(instrument,      "Handle to the instrument");
        FIELD(marketDataCache, "Handle to the Market Data");
        FIELD(model,           "Handle to the Model");

        // register addin for adding objects to cache
        Addin::registerClassObjectMethod("EXPAND_OBJECT",
                                         Addin::UTILITIES,
                                         "Retrieves all wrapper objects from the market",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)expand);
    }

    static IObject* defaultExpandObjects() {
        return new ExpandObjects();
    }
};

CClassConstSP const ExpandObjects::TYPE= CClass::registerClassLoadMethod(
    "ExpandObjects", typeid(ExpandObjects), ExpandObjects::load);

/** Addin to determine run-time type of a handle */
class N2Addin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    double           a;
    double           b;
    double           r;



    /** the 'addin function' - construct object from components */
    static double getN2(N2Addin* params) {
        static const string routine = "N2Addin::getN2";
        return N2(params->a, params->b, params->r);
    }

    /** for reflection */
    N2Addin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(N2Addin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultN2Addin);
        FIELD(a, "a");
        FIELD(b, "b");
        FIELD(r, "r");
        Addin::registerClassDoubleMethod("N2",
                                         Addin::UTILITIES,
                                         "Returns the cumulative bivariate normal",
                                         TYPE,
                                         (Addin::DoubleMethod*)getN2);
    }

    static IObject* defaultN2Addin() {
        return new N2Addin();
    }
};

CClassConstSP const N2Addin::TYPE= CClass::registerClassLoadMethod(
    "N2Addin", typeid(N2Addin), N2Addin::load);


// Addin to determine run-time type of a handle
class PseudoEquityVolAddin: public CObject{
    static CClassConstSP const TYPE;

    // parameters that our addin takes
    AssetSP         equity;
    DateTime        valueDate;
    DateTime        maturity;
    double          strike;

    // the 'addin function' - construct object from components
    static double getVol(PseudoEquityVolAddin* params){
        static const string routine = "PseudoEquityVolAddin::getVol";

        DateTimeArray divCritDates;
        PseudoSimpleEquitySP pseudoAsset(PseudoSimpleEquity::create(params->equity.get(),
                                                                      divCritDates,
                                                                      true, 0));

        CVolRequestSP volRequest(new LinearStrikeVolRequest(params->strike,
                                                            params->valueDate,
                                                            params->maturity,
                                                            false));

        CVolProcessed*    vol   = pseudoAsset->getProcessedVol((CVolRequest*)volRequest.get());
        CVolProcessedBSSP volLN = CVolProcessedBSSP(dynamic_cast<CVolProcessedBS*>((IObject*)vol));
        if (!volLN) {
            throw ModelException("CAsset::getProcessedVol", "Failed to get LN"
                                 " vol from request of type "+
                                 volRequest->getClass()->getName());
        }


        // calculate the pseudo spot price
        pseudoAsset->getSpot();

        // calculate the volatility
        double pseudoVol = volLN->CalcVol(params->valueDate, params->maturity);

        return pseudoVol;
    }

    // for reflection
    PseudoEquityVolAddin():  CObject(TYPE){}

    // Invoked when Class is 'loaded'
    static void load(CClassSP& clazz) {
        REGISTER(PseudoEquityVolAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPseudoEquityVolAddin);
        FIELD(equity,            "asset handle");
        FIELD(valueDate,  "today")
        FIELD(maturity,   "option maturity")
        FIELD(strike,     "strike")

        Addin::registerClassDoubleMethod("GET_PSEUDO_EQUITY_VOL",
                                         Addin::MARKET,
                                         "Returns the pseudo equity volatility",
                                         TYPE,
                                         (Addin::DoubleMethod*)getVol);
    }

    static IObject* defaultPseudoEquityVolAddin() {
        return new PseudoEquityVolAddin();
    }
};

CClassConstSP const PseudoEquityVolAddin::TYPE= CClass::registerClassLoadMethod(
    "PseudoEquityVolAddin", typeid(PseudoEquityVolAddin), PseudoEquityVolAddin::load);


// Addin to determine run-time type of a handle
class StockFloorAddin: public CObject { 
    static CClassConstSP const TYPE;

    // parameters that our addin takes
    AssetSP         equity;

    // the 'addin function' - construct object from components
    static double getStockFloor(StockFloorAddin* params) {
        static const string routine = "StockFloorAddin::getStockFloor";

        DateTimeArray divCritDates;
        PseudoSimpleEquitySP pseudoAsset(PseudoSimpleEquity::create(params->equity.get(),
                                                                      divCritDates,
                                                                      true, 0));
        double stockFloor = pseudoAsset->getSpot();
        return stockFloor;
    }

    // for reflection
    StockFloorAddin():  CObject(TYPE){}

    // Invoked when Class is 'loaded'
    static void load(CClassSP& clazz) {
        REGISTER(StockFloorAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultStockFloorAddin);
        FIELD(equity,            "asset handle");

        Addin::registerClassDoubleMethod("GET_STOCK_FLOOR",
                                         Addin::MARKET,
                                         "Returns the stock floor of an equity",
                                         TYPE,
                                         (Addin::DoubleMethod*)getStockFloor);
    }

    static IObject* defaultStockFloorAddin() {
        return new StockFloorAddin();
    }
};

CClassConstSP const StockFloorAddin::TYPE= CClass::registerClassLoadMethod(
    "StockFloorAddin", typeid(StockFloorAddin), StockFloorAddin::load);

DRLIB_END_NAMESPACE
