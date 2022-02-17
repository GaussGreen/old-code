//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolPreferred.cpp
//
//   Description : Redirects to a real vol, on an asset basis function
//
//   Author      : JNJ
//
//   Date        : 09 Nov 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/NextStrike.hpp"
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE

/** A CVolBase which wraps another class of CVolBase, keyed from its
    name.  Need to review whether we should implement the
    IVolatilityBS and IVolatilityDVF interfaces - since the vol
    contained might not. (The alternative would be to alter the market
    data retrieval somehow so that the wrapped vol would be obtained
    and not the VolPreferred.) Also having this object in the cache which
    is set to ask for a generic black-scholes vol will probably fail because
    this object [currently] is a black-scholes vol - as well as the object
    being asked for. If this class no longer implements these two interfaces
    then we need to alter ModelLN.cpp so when it checks to see if the
    specified class is a IVolatilityBS it has special code for VolPreferred */
class VolPreferred: public CVolBase,
                    virtual public INextStrike,
                    virtual public IVolatilityBS,
                    virtual public IVolatilityDVF,
                    virtual public IPDFCalculator,
                    virtual public VolRelativeShift::IShift
{
public:
    static CClassConstSP const TYPE;

    /** Returns name of vol */
    virtual string getName() const{
        return name;
    }

    virtual CClassConstSP proxyType(const MarketData* marketData,
                                    const IModel*     model) const{
        static const string method = "VolPreferred::proxyType";
        if (model) {
            // specialise via map (engine, volToUse)
            const IModelFamily* modelFamily = dynamic_cast<const IModelFamily*>(model);
            // getFamilyName() can return 0
            const string* modelFamilyName = modelFamily?modelFamily->getFamilyName() : 0;
            if (modelFamilyName) {  
                int engineId = 0;
                for(; engineId<engine.size(); engineId++) {
                    if (engine[engineId] == *modelFamilyName) {
                        break;
                    }
                }
                if (engineId>=engine.size()) {
                    // not found - so ....
                    throw ModelException(method, 
                                         "Failed to locate engine family " + 
                                         *modelFamilyName + 
                                         " for " + getName());
                }
                return CClass::forName(volToUse[engineId]);
            }
            // Note the fall-through which means any models not implementing
            // IModelFamily simply proceed as before - and use volClass.
        }
        return volClass;
    }

    /** Creates processed vol */
    CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                   const CAsset*      asset) const{
        if (!myRealVolBase){
            throw ModelException("VolPreferred::getProcessedVol",
                                 "Market data missing");
        }
        return myRealVolBase->getProcessedVol(volRequest, asset);
    }

    /** Creates Struck processed vol */
    CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const{
        if (!myRealVolBase){
            throw ModelException("VolPreferred::getProcessedVol",
                                 "Market data missing");
        }
        return myRealVolBase->getProcessedVol(volRequest,
                                              eqAsset,
                                              fxAsset,
                                              eqFXCorr);
    }

    void validatePop2Object(){
        static const char routine[] = "VolPreferred::validatePop2Object";
        try {
            volClass = CClass::forName(volType);
        } catch (exception& ) {
            throw ModelException(routine, 
                                 volType + " is not an existing Vol Type");
        }
        if (!CVolBase::TYPE->isAssignableFrom(volClass)) {
            throw ModelException(routine,
                                 volType + " is not an existing Vol Type");
        }
        if (volType == TYPE->getName()) {
            throw ModelException(routine,
                                 volType + " is not a valid choice of VolPreferred");           
        }
        if (engine.size() != volToUse.size()) {
            throw ModelException(routine, "engine and volToUse arrays are of "
                                 "different lengths");
        }            
    }

    /** override clone method to copy over volClass field */
    virtual IObject* clone() const{
        VolPreferred& volPref = 
            dynamic_cast<VolPreferred&>(*CVolBase::clone());
        volPref.volClass = volClass;
        return &volPref;
    }

    //// pull out vol with same name as ourselves but of specified type 
    void getMarket(const IModel* model, const MarketData* market) {
        CVolBase::getMarket(model, market);
        validatePop2Object();
        try{
            myRealVolBase = CVolBaseSP::dynamicCast(
                (IObjectSP)(market->GetData(getName(), volClass)));
            /* we now have to manually call getMarket on our realVolBase here.
               Normally this would be done for us, but we are bypassing the
               normal route here (normally routed through Model but that
               would result in us going in circles) */
            myRealVolBase->getMarket(model, market);
        } catch (exception& e){
            throw ModelException(e, "CVolPreferred::getMarket");
        }
    }

    PDFCalculator* getPDFCalculator(
        const PDFRequest* request,
        const CAsset*     asset) const {
        static const string method("VolPreferred::getPDFCalculator");
        try {
            if (IPDFCalculator::TYPE->isInstance(myRealVolBase.get())) {
                const IPDFCalculator* pdf = dynamic_cast< const IPDFCalculator*>(myRealVolBase.get());
                return pdf->getPDFCalculator(request, asset);
            }
            throw ModelException(method,
                                 "vol of type (" + myRealVolBase->getClass()->getName() +
                                 ") has no pdf calculator");
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    double getNextStrike(const double& strike,
                         bool          isUp,
                         bool&         offSurface) const
    {
        double volStrike;
        const INextStrike* nextStrike = dynamic_cast<const INextStrike*>(myRealVolBase.get());
        if ( nextStrike )
        {
            volStrike = nextStrike->getNextStrike(strike,
                                                  isUp,
                                                  offSurface);
        }
        else
        {
            volStrike  = 0.0;
            offSurface = true;
        }
        return volStrike;
    }

    /** VolRelativeShift Interface */
    // a bit nasty, but it's all very LN vol based and don't want (say) CEVJ
    // to fail, so asset needs to check if vol will work, so we need to pass
    // down to the actual vol in use
    string sensName(VolRelativeShift* shift) const{
        return VolRelativeShift::IShift::TYPE->isInstance(myRealVolBase.get())
            ? getName(): "";
    }
    bool sensShift(VolRelativeShift* shift){
        return true;  // tweak the real vol
    }

private:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolPreferred, clazz);
        SUPERCLASS(CVolBase);
        IMPLEMENTS(INextStrike);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(IPDFCalculator);
        IMPLEMENTS(VolRelativeShift::IShift);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(name, "Vol identifier");
        FIELD(volType, "the Preferred Vol Type");
        FIELD(engine, "engines with preferred over-rides");
        FIELD_MAKE_OPTIONAL(engine);
        FIELD(volToUse, "vol to use for given engine");
        FIELD_MAKE_OPTIONAL(volToUse);

        FIELD(myRealVolBase, "Internal");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(myRealVolBase);
    }

    /** default constructor. NB important to set volClass to null since it's
        really a pointer */
    VolPreferred(): CVolBase(TYPE), volType("VolSurface"), volClass(0) { }

    static IObject* defaultCtor(){ 
        return new VolPreferred();
    }


    // registered fields
    string          name;    // name of the vol
    string          volType;
    // volType is default, allow to over-ride by engine family
    StringArray     engine;    
    StringArray     volToUse;

    CVolBaseSP      myRealVolBase;  // transient
    // cached field
    CClassConstSP   volClass;       // cached in validatePop2Object $unregistered
};

CClassConstSP const VolPreferred::TYPE =
CClass::registerClassLoadMethod("VolPreferred", typeid(VolPreferred), load);

// external symbol to allow class to be forced to be linked in
bool VolPreferredLinkIn(){
    return (VolPreferred::TYPE != 0);
}

DRLIB_END_NAMESPACE
