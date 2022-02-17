//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PyramidLegacyFFX.cpp
//
//   Description : An FX Forward that looks and acts just like the one in the
//                 EDG library to allow migration onto QLib without having
//                 to answer those tricky questions like "why can't it be a
//                 forward on an FX asset?"
//
//   Author      : Andrew J Swain
//
//   Date        : 26 June 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/Theta.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Results.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/Spot.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/RhoParallel.hpp"
#include "edginc/RhoPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL BogusFFXLeg: public CInstrument, 
                                virtual public ClosedForm::IIntoProduct,
                                virtual public IRestorableWithRespectTo<Spot>,
                                virtual public LastSensDate,                  
                                virtual public Theta::Shift {
public:
    static CClassConstSP const TYPE;
    friend class BogusFFXLegClosedForm;

    BogusFFXLeg(const string&     name,
                const DateTime&   valueDate,
                const DateTime&   maturity,
                double            notional,
                double            spotFX,
                double            scale,
                const YieldCurve* legCcy,
                const YieldCurve* usdCcy): CInstrument(TYPE), 
        name(name), valueDate(valueDate), maturity(maturity), notional(notional),
        spotFX(spotFX), scale(scale), legCcy(copy(legCcy)), usdCcy(copy(usdCcy)) {}


    // nothing to see here
    virtual void GetMarket(const IModel* model, const CMarketDataSP market) {}

    /** instrument validation */
    virtual void Validate() {
        // think of something later
    }
        
    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;
    
    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift) {
        valueDate = shift->rollDate(valueDate);
        return true;  //carry on recursing
    }

    /** what's today ? */
    virtual DateTime getValueDate() const {
        return valueDate;
    }
    
    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const {
        return maturity;
    }

    // value a "leg"
    void price(Control* control, CResults* results) const {
	static const string method("BogusFFXLeg::price");
        try {
            double value = notional*legCcy->pv(maturity);
            results->storePrice(value/scale, legCcy->getCcy());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
      
    // handle any output requests
    void requests(Control* control, CResults* results) const {
	static const string method("BogusFFXLeg::requests");
        try {
            double fwdFX = spotFX*usdCcy->pv(maturity)/legCcy->pv(maturity);
            double value = notional*legCcy->pv(maturity);
            string name  = legCcy->getName();

            if (control && control->isPricing()) {

                OutputRequest* request = control->requestsOutput(OutputRequest::FWD_FX_RATE);
                if (request) {
                    results->storeRequestResult(request, fwdFX, name);
                }
                request = control->requestsOutput(OutputRequest::FFX_PV);
                if (request) {
                    results->storeRequestResult(request, value, name);
                }
                request = control->requestsOutput(OutputRequest::FFX_FUTURE_VALUE);
                if (request) {
                    results->storeRequestResult(request, notional/fwdFX, name);
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const {
        return usdCcy->getName(); 
    }

    /** Returns name identifying this object for Delta */
    virtual string sensName(const Spot*) const {
        if (usdCcy->getName() == legCcy->getName()) {
            return ""; // no point doing USD/USD FX delta now is there?
        }
        return name; 
    }

    /** Shifts the object using given shift (see Delta::Shift)*/
    virtual TweakOutcome sensShift(const PropertyTweak<Spot>& shift) {
        double orig = spotFX;
        spotFX *= (1.0 + shift.coefficient);
        // Arg1 is the initial spot to be stored on the delta shift
        // Arg2-Arg1 is the required divisor
        return TweakOutcome(orig, spotFX, false);   // none of our components has a delta type sensitivity        
    }

    /** Restores the object to its original form */
    void sensRestore(const PropertyTweak<Spot>& shift) {
        spotFX /= (1.0 + shift.coefficient);
    }
        
    double spot() const {
        return spotFX;
    }

private:
    BogusFFXLeg(): CInstrument(TYPE), notional(0.0), spotFX(0.0), scale(0.0) {}

    BogusFFXLeg(const BogusFFXLeg& rhs);
    BogusFFXLeg& operator=(const BogusFFXLeg& rhs);

    // fields
    string              name;
    DateTime            valueDate;
    DateTime            maturity;
    double              notional;
    double              spotFX;
    double              scale;
    YieldCurveSP        usdCcy;
    YieldCurveSP        legCcy;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BogusFFXLeg, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedForm::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(IRestorableWithRespectTo<Spot>);
        IMPLEMENTS(Theta::IShift);
        EMPTY_SHELL_METHOD(defaultBogusFFXLeg);
        FIELD(name, "name");
        FIELD(valueDate, "valueDate");
        FIELD(maturity, "maturity");
        FIELD(notional, "notional");
        FIELD(spotFX, "spotFX");
        FIELD(scale, "scale");
        FIELD(usdCcy, "usdCcy");
        FIELD(legCcy, "baseCcy");
    }

    static IObject* defaultBogusFFXLeg(){
        return new BogusFFXLeg();
    }

};

typedef smartPtr<BogusFFXLeg> BogusFFXLegSP;

class DuffFFXLeg: public ClosedForm::IProduct{
private:
    const BogusFFXLeg* ffx; // a reference

public:
    DuffFFXLeg(const BogusFFXLeg* ffx): ffx(ffx){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const{
        ffx->price(control, results);
    }
};


/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* BogusFFXLeg::createProduct(ClosedForm* model) const {
    return new DuffFFXLeg(this);
}
 

CClassConstSP const BogusFFXLeg::TYPE = CClass::registerClassLoadMethod(
    "BogusFFXLeg", typeid(BogusFFXLeg), BogusFFXLeg::load);



class PRODUCTS_DLL PyramidLegacyFFX: public CInstrument, 
                   virtual public ClosedForm::IIntoProduct,
                   virtual public LastSensDate,                  
                   virtual public Theta::Shift {
public:
    static CClassConstSP const TYPE;
    friend class PyramidLegacyFFXClosedForm;

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel* model, const CMarketDataSP market) {
	static const string method("PyramidLegacyFFX::GetMarket");
        try {
            market->GetReferenceDate(valueDate);
            usdCcy.getData(model, market);
            baseCcy.getData(model, market);
            riskCcy.getData(model, market);

            if (baseCcy->getName() != usdCcy->getName()) {
                baseFX = market->getFXName(usdCcy->getName(), baseCcy->getName());
            }
            if (riskCcy->getName() != usdCcy->getName()) {
                riskFX = market->getFXName(usdCcy->getName(), riskCcy->getName());
            }

            double baseNotl = notional;
            double riskNotl = -notional*dealFXRisk/dealFXBase;

            baseLeg = BogusFFXLegSP(new BogusFFXLeg(baseFX,
                                                    valueDate,
                                                    maturity,
                                                    baseNotl,
                                                    spotFXBase,
                                                    notional,
                                                    baseCcy.get(),
                                                    usdCcy.get()));

            riskLeg = BogusFFXLegSP(new BogusFFXLeg(riskFX,
                                                    valueDate,
                                                    maturity,
                                                    riskNotl,
                                                    spotFXRisk,
                                                    notional,
                                                    riskCcy.get(),
                                                    usdCcy.get()));
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** instrument validation */
    virtual void Validate() {
        // think of something later
    }
        
    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;
    
    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift) {
        valueDate = shift->rollDate(valueDate);
        return true;  //carry on recursing
    }

    /** what's today ? */
    virtual DateTime getValueDate() const {
        return valueDate;
    }
    
    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const {
        return maturity;
    }
       
    // create a control to handle the bizarre leg-specific tweaking
    Control* fxCtrl(Control* ctrl) const {
        CControlSP rhoCtrl(new Control());
        rhoCtrl->validatePop2Object();

        SensitivitySP parallel(ctrl->sensitivityRequested(RhoParallel::TYPE));
        SensitivitySP pwise(ctrl->sensitivityRequested(RhoPointwise::TYPE));

        if (parallel.get()) {
            rhoCtrl->addSensitivity(parallel);
        }
        if (pwise.get()) {
            rhoCtrl->addSensitivity(pwise);
        }

        return rhoCtrl.release();
    }

    // the pricing itself
    void price(ClosedForm* model, Control* control, CResults* results) const  {
	static const string method("PyramidLegacyFFX::price");
        try {
            double value = 0.0;

            if (valueDate < maturity) {
 
                ResultsSP  baseResults(new Results());
                ResultsSP  riskResults(new Results());
                CControlSP baseCtrl(fxCtrl(control));
                CControlSP riskCtrl(fxCtrl(control));

                baseCtrl->calculate(model, baseLeg.get(), baseResults.get());
                double base = baseResults->retrievePrice();
                riskCtrl->calculate(model, riskLeg.get(), riskResults.get());
                double risk = riskResults->retrievePrice();

                const CControlSP addmebase(baseCtrl);
                const CControlSP addmerisk(riskCtrl);

                results->add(baseResults.get(),
                             addmebase,
                             1.0,
                             true);

                results->add(riskResults.get(),
                             addmerisk,
                             1.0,
                             false);

                value = (base/baseLeg->spot() + risk/riskLeg->spot())*baseLeg->spot();
           
                if (control) {

                    OutputRequest* request = control->requestsOutput(OutputRequest::USD_DISC_RATE);
                    if (request) {
                        double usdZero = usdCcy->zero(maturity);
                        results->storeRequestResult(request, usdZero);
                    }

                    request = control->requestsOutput(OutputRequest::FFX_UFV);
                    if (request) {
                        results->storeRequestResult(request, value/spotFXBase);
                    } 

                    baseLeg->requests(control, results);
                    riskLeg->requests(control, results);
                }
            }

            results->storePrice(value, baseCcy->getCcy());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const {
        return usdCcy.getName(); 
    }


private:
    PyramidLegacyFFX(): CInstrument(TYPE), notional(0.0), spotFXBase(0.0),
                        spotFXRisk(0.0), dealFXBase(0.0), dealFXRisk(0.0) {}

    PyramidLegacyFFX(const PyramidLegacyFFX& rhs);
    PyramidLegacyFFX& operator=(const PyramidLegacyFFX& rhs);

    // fields
    DateTime            valueDate;
    DateTime            maturity;
    double              notional;
    YieldCurveWrapper   usdCcy;
    YieldCurveWrapper   baseCcy;
    YieldCurveWrapper   riskCcy;
    double              spotFXBase;
    double              spotFXRisk;
    double              dealFXBase;
    double              dealFXRisk;

    // transients
    string              baseFX;
    string              riskFX;
    BogusFFXLegSP       baseLeg;    
    BogusFFXLegSP       riskLeg;    

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("EDG lib style FX forward - don't use this unless"
                              " you really know what you are doing and why");
        REGISTER(PyramidLegacyFFX, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedForm::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::IShift);
        EMPTY_SHELL_METHOD(defaultPyramidLegacyFFX);
        FIELD(valueDate, "valueDate");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(maturity, "maturity");
        FIELD(notional, "notional");
        FIELD(usdCcy, "usdCcy");
        FIELD(baseCcy, "baseCcy");
        FIELD(riskCcy, "riskCcy");
        FIELD(spotFXBase, "spotFXBase");
        FIELD(spotFXRisk, "spotFXRisk");
        FIELD(dealFXBase, "dealFXBase");
        FIELD(dealFXRisk, "dealFXRisk");
        FIELD(baseFX, "");
        FIELD_MAKE_TRANSIENT(baseFX);
        FIELD(riskFX, "");
        FIELD_MAKE_TRANSIENT(riskFX);
        FIELD(baseLeg, "");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(baseLeg);
        FIELD(riskLeg, "");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(riskLeg);
    }

    static IObject* defaultPyramidLegacyFFX(){
        return new PyramidLegacyFFX();
    }

};

class DuffFFX: public ClosedForm::IProduct{
private:
    const PyramidLegacyFFX* ffx; // a reference

public:
    DuffFFX(const PyramidLegacyFFX* ffx): ffx(ffx){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const {
        ffx->price(model, control, results);
    }
};


/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* PyramidLegacyFFX::createProduct(ClosedForm* model) const {
    return new DuffFFX(this);
}
 

CClassConstSP const PyramidLegacyFFX::TYPE = CClass::registerClassLoadMethod(
    "PyramidLegacyFFX", typeid(PyramidLegacyFFX), 
    PyramidLegacyFFX::load);


/* for class loading */
bool PyramidLegacyFFXLoad() {
    return (PyramidLegacyFFX::TYPE != 0);
}

DRLIB_END_NAMESPACE
