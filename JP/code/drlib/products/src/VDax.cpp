//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VDax.cpp
//
//   Description : VDax do not call variance index forward
//
//   Author      : 
//
//   Date        : 
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/VolSV.hpp"
#include "edginc/VolVarSwap.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/AssetUtil.hpp"

static const double MAX_TENOR = 0.3; // hard coded, however, is just a rule of thumb and pretty arbitrary

DRLIB_BEGIN_NAMESPACE

class VDaxModel : public CModel {
public:
    static CClassConstSP const TYPE;
    friend class VDax;
    friend class VDaxPricer;

// **************************************************
    class VDaxBasisApprox : public CObject { 
// **************************************************
// purpose: VDax is priced as approx + basis, where the basis is computed using the Heston Model
//          this class is used in order to compute the basis using the Heston Model (see VolSV.cpp for details)
    public:
        static CClassConstSP const TYPE;

        void validatePop2Object(){
            static const string method = "VDaxBasisApprox::validatePop2Object";
            try {
                if (!(volType==VolSV::TYPE->getName())) {
                    throw ModelException("VDaxBasisApprox::validatePop2Object",
                        "VDaxBasisApprox only supports VolSV, but " + volType + " is chosen.");
                }
            }
            catch (exception& e){
                throw ModelException(e, method);
            }
        };

        void getData(const IModel* model, const MarketData* market, const string& name) const{
            static const string method = "VDaxBasisApprox::getData";
            try {
                static CClassConstSP volClass(CClass::forName("VolSV"));
                if (!theVol.get()) {
                    MarketObjectSP data(market->GetData(name, volClass));
                    if (data.get()){
                        VolSV* temp = dynamic_cast<VolSV*>(data.get());
                        theVol = VolSVSP(copy(temp));
                    }
                    if (!theVol){
                        throw ModelException("Failed to get VolSV");
                    }
                    // must absolutely getMarket for copied vol object
                    theVol->getMarket(model, market, name);
                }
            }
            catch (exception& e){
                throw ModelException(e, method);
            }
        }

        double varSwap(const DateTime& valueDate, const DateTime& maturity) const {
            static const string method = "VDaxBasisApprox::varSwap";
            try {
                return theVol->varSwap(valueDate, maturity);
            }
            catch (exception& e){
                throw ModelException(e, method);
            }
        }
        double impliedSqVol(const DateTime& valueDate, 
                            const DateTime& maturityFwd, 
                            const DateTime& maturityEnd) const {
            static const string method = "VDaxBasisApprox::impliedSqVol";
            try {
                return theVol->impliedSqVol(valueDate, maturityFwd, maturityEnd);
            }
            catch (exception& e){
                throw ModelException(e, method);
            }
        }
        double calcTradingTime(const DateTime& maturity1, const DateTime& maturity2) const {
            static const string method = "VDaxBasisApprox::impliedSqVol";
            try {
                const TimeMetric& metric = theVol->getTimeMetric(); 
                return metric.yearFrac(maturity1, maturity2);
            }
            catch (exception& e){
                throw ModelException(e, method);
            }
        }

        // registration, invoked when class is 'loaded' 
        static void load(CClassSP& clazz){
            clazz->setPublic(); 
            REGISTER(VDaxBasisApprox, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultVDaxBasisApprox);
            FIELD(volType, "Type of vol to use");
            FIELD(theVol, "");
            FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(theVol);
        }
    
        static IObject* defaultVDaxBasisApprox(){
            return new VDaxBasisApprox();
        }
    
    private:
        // constructor
        VDaxBasisApprox():CObject(TYPE) {}
        
        string volType;
        mutable VolSVSP theVol;
    }; 
    typedef smartPtr<VDaxBasisApprox> VDaxBasisApproxSP;
    // **************************************************

    virtual void validatePop2Object() {
        static const string method("VDaxModel::validatePop2Object");
        try {
            if (!(CString::equalsIgnoreCase(methodChoice,"VAR_SWAP_REL")) && 
                !(CString::equalsIgnoreCase(methodChoice,"VAR_SWAP_ABS")) && 
                !(CString::equalsIgnoreCase(methodChoice,"ATM_VAR_REL")) &&
                !(CString::equalsIgnoreCase(methodChoice,"ATM_VAR_ABS")) ) {
                throw ModelException("Method Choice has to be either VAR_SWAP_ABS or VAR_SWAP_REL or ATM_VAR_ABS or ATM_VAR_REL, but it is " + methodChoice );
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
        
    /** the class that the product must be able to create */
    class IProduct{
    public:
        virtual void price(VDaxModel* model,
                           Control*   control, 
                           Results*   results) const = 0;
        virtual ~IProduct(){};
    };

    /** interface that the instrument must implement */
    class IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(VDaxModel* model) const = 0;
    };

    // inherited from CModel - mandatory to implement
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results) {
        static const string method = "VDaxModel::Price";
        if (!IIntoProduct::TYPE->isInstance(instrument)){
            throw ModelException("Instrument of type "+
                                 instrument->getClass()->getName() +
                                 " does not support VDaxModel::IntoProduct");
        }
        try {
            if (instrument->priceDeadInstrument(control, results)) {
                return; // done for a dead instrument
            }
            // cast to VDaxModel::IIntoProduct
            IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
            // create the product
            auto_ptr<IProduct> product(intoProd.createProduct(this));
            product->price(this, control, results);
        } 
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
     
    MarketObjectSP GetMarket(const MarketData*    market,
                             const string&        name,
                             const CClassConstSP& type) const{
        static const string method = "VDaxModel::GetMarket";
        try {
            basisApproxModel->getData(this, market,name);
            return modelForVarSwap->GetMarket(market, name, type);
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const {
        return modelForVarSwap->wantsRiskMapping();
    }

    MarketObjectSP modifyMarketData(
        const MarketData*     market,
        const CClassConstSP&  clazz,  // what type was originally requested
        const MarketObjectSP& mo) const   /* what GetMarket returned or what was
                                         "inline" already */
    {
        return modelForVarSwap->modifyMarketData(market, clazz, mo);
    }


    // for VDaxModel::IIntoProduct
    static void IntoProduct_load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(Model::IModelIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    // registration, invoked when class is 'loaded'
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(VDaxModel, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultVDaxModel);
        FIELD(modelForVarSwap,"model");
        FIELD(basisApproxModel, "basisApproxModel");
        FIELD(methodChoice, "either VAR_SWAP_REL, VAR_SWAP_ABS, ATM_VAR_REL or ATM_VAR_ABS");
        FIELD_MAKE_OPTIONAL(methodChoice); // default value is VAR_SWAP
    }

    // for VDaxModel::IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(VDaxModel::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultVDaxModel(){
        return new VDaxModel();
    }
    // constructor 
    VDaxModel():CModel(TYPE), methodChoice("VAR_SWAP_ABS") {}
private:
    
    // registered fields
    IModelSP                modelForVarSwap;    // either ImpliedIntegration of ClosedFormIntegrateLN
    VDaxBasisApproxSP       basisApproxModel;   // VDaxBasisApprox
    string                  methodChoice;       // two possibilities: via VAR_SWAP or via ATM_VAR
};
typedef smartPtr<VDaxModel> VDaxModelSP;

class VDax : public Generic1Factor, 
             virtual public VDaxModel::IIntoProduct,
             virtual public ISensitiveStrikes,
             virtual public LastSensDate {
public:
    static CClassConstSP const TYPE; 
    
    friend class VDaxModel;

    // do some validation at construction -- not mandatory
    virtual void validatePop2Object() {
        static const string method("VDax::validatePop2Object");
        try {
            if (instSettle->isPhysical() || instSettle->isMargin()) {
                throw ModelException("Only cash settlement is allowed.");
            }
            if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
                throw ModelException("Swap can't be ccy struck.");
            }
            if (ccyTreatment == CAsset::CCY_TREATMENT_PROTECTED) {
                throw ModelException("Swap can't be ccy protected.");
            }
            if (FXAsset::TYPE->isInstance(asset.get())){
                throw ModelException("FX underlying not allowed");
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
        
    // do some asset specific validation -- mandatory since VDax : Generic1Factor : CInstrument
    virtual void Validate() {
        static const string method("VDax::validate");
        try {
            AssetUtil::assetCrossValidate(asset.get(),
                               false,
                               valueDate,
                               valueDate,
                               discount,
                               this);
        }    
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // indicates whether VEGA_MATRIX is sensible for this instrument
    virtual bool avoidVegaMatrix(const IModel* model) {
        return false;
    }
        
    // returns all strikes on the vol surface to which this instrument is sensitive
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model) {
        static const string method("VolVarShell::getSensitiveStrikes");
        try {
            if (avoidVegaMatrix(model)) {
                throw ModelException(method, "VEGA_MATRIX is not valid for this instrument.");
            }
            const VDaxModel* vdx = dynamic_cast<const VDaxModel*>(model);

            if (!vdx) {
                throw ModelException("Only VDaxModel allowed for VDax.");
            }
            
            const DateTime& maturityEnd = endDate();
            DoubleArraySP sensStrikes = VarianceSwapUtil::getSensitiveStrikes(outputName,
                                                                              asset.get(),
                                                                              valueDate,
                                                                              maturityEnd,
                                                                              vdx->modelForVarSwap.get());
            return sensStrikes;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    DateTime endDate() const{
        MaturityPeriod periodRef(tenorPeriod);            
        DateTime maturityEnd = DateTime(periodRef.toDate(maturityFwd).getDate(), maturityFwd.getTime());
        return maturityEnd;
    }

    DateTime endDate(const Sensitivity* sensControl) const {
        return endDate();
    }

    // price a dead instrument until settlement
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const {
        static const string method = "VDax::priceDeadInstrument";
        try {
            bool deadinstrument = false;
            // compare maturityFwd + settlement Date vs valueDate 
            DateTime settleDate  = instSettle->settles(maturityFwd, asset.get());

            if (maturityFwd.isGreater(valueDate)) {
                // nthg happens, i.e. normal price
            } else if (valueDate.isGreaterOrEqual(settleDate)) {
                // zero value
                results->storePrice(0.0, discount->getCcy());
                
                deadinstrument = true;   
                
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, 0.0);
                }
            } else if (valueDate.isGreaterOrEqual(maturityFwd)&&valueDate.isLess(settleDate)) {
                // check that there is user input
                if (Maths::isNegative(priceAtFwd)) {
                    throw ModelException("Pls enter a valid price at the maturity of the foward");
                }
                double value = priceAtFwd * discount->pv(settleDate);
                results->storePrice(value, discount->getCcy());
                
                deadinstrument = true;
                
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, value);
                }
            }
            return deadinstrument;
        }
        catch (exception& e){
                throw ModelException(e, method);
        }
    }

    // extra output request
    void recordOutputRequests(Control* control, 
                              Results* results, 
                              double   value) const {
        static const string method = "VDax::priceDeadInstrument";
        try {
            OutputRequest* request = NULL;
        
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTime paymentDate = instSettle->settles(maturityFwd, asset.get());
                DateTimeArray date(1, paymentDate);
                OutputRequestUtil::recordPaymentDates(control,results,&date); 
            }

            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                if (valueDate.isGreaterOrEqual(maturityFwd)) {
                    DateTime paymentDate = instSettle->settles(maturityFwd, asset.get());
                    CashFlow cf(paymentDate, value);
                    CashFlowArray cfl(1, cf);
                    OutputRequestUtil::recordKnownCashflows(control,
                                                            results,
                                                            discount->getCcy(),
                                                            &cfl);  
                }
            }
            
            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           maturityFwd,
                                           valueDate,
                                           asset.get());
        }
        catch (exception& e){
                throw ModelException(e, method);
        }
    }
    
    // registration, invoked when class is 'loaded'
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VDax, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(VDaxModel::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultVDax);
        FIELD(maturityFwd, "maturity of forward contract");
        FIELD(tenorPeriod, "length of tenor");
        FIELD(priceAtFwd,"price when the forward contract expires");
        FIELD_MAKE_OPTIONAL(priceAtFwd);
    }
    
    static IObject* defaultVDax(){
        return new VDax();
    }

    VDaxModel::IProduct* createProduct(VDaxModel* model) const;
    
private:
    friend class VDaxPricer;
    // registered fields
    DateTime                maturityFwd;        // maturity of forward contract (ie date of cash settlement)
    string                  tenorPeriod;
    double                  priceAtFwd;
    
    // constructor
    VDax(): Generic1Factor(TYPE), priceAtFwd(0.0) {}
  
}; // end of class VDax


// need a pricer class, since pricing is done here (neither mc nor tree do it) 
class VDaxPricer: virtual public VDaxModel::IProduct {
private:
    const VDax* inst; // a reference
public:
    VDaxPricer(const VDax* instrument): inst(instrument){}
    
    void price(VDaxModel*       model,
               Control*         control, 
               CResults*        results) const{
        static const string method = "VDaxPricer::price";
        try {
            MaturityPeriod periodRef(inst->tenorPeriod);            
            DateTime maturityEnd = periodRef.toDate(inst->maturityFwd);
                        
            // computation of time as in VolVarSwap.cpp, so that consistency
            double T = model->basisApproxModel->calcTradingTime(inst->valueDate,inst->maturityFwd);
            double delta = model->basisApproxModel->calcTradingTime(inst->maturityFwd,maturityEnd);
            if (Maths::isZero(delta)) {
                throw ModelException("The tenor period must not be zero!");
            }
            if (delta>MAX_TENOR){
                throw ModelException("The tenor period must not be longer than 3 months, but is "
                    + inst->tenorPeriod + ".");
            }

            // Compute VarSwap from MatFwd to MatEnd in Heston Model
            double HestonVarSwapToMatEnd = model->basisApproxModel->varSwap(inst->valueDate, maturityEnd);
            double HestonVarSwapToMatFwd = model->basisApproxModel->varSwap(inst->valueDate, inst->maturityFwd);

            double HestonVarSwapBetwMatFwdMatEnd = 
                ((T+delta) * HestonVarSwapToMatEnd - T * HestonVarSwapToMatFwd) / delta;
 
            /* final value of product = price(VDax) = approxValue + basis, where basis and approxValue is different 
               for the different methodologies => two possibilities: absolute basis and relative basis, ie 
               price(VDax) = approxValue + basis or price(VDax) + approxValue * (1+betaBasis) */
            double value;

            double VDaxHeston = 
                model->basisApproxModel->impliedSqVol(inst->valueDate, inst->maturityFwd, maturityEnd);

            if (CString::equalsIgnoreCase(model->methodChoice,"VAR_SWAP_ABS") ||
                    CString::equalsIgnoreCase(model->methodChoice,"VAR_SWAP_REL")){
                // from VolVarSwap.cpp, we get variance as output                
               
                /** if we are doing vega matrix, then we use in any case ClosedFormIntegrateLN as model, 
                    even in the case where ImpliedIntegration is the "main" model; 
                    hence, pass control in order to find out whether we are doing vega matrix or not */
                    
                double approxValue = VarianceSwapUtil::futureFwdVar(inst->asset.get(),
                                                             inst->discount.get(),
                                                             inst->valueDate,
                                                             inst->maturityFwd,
                                                             maturityEnd,
                                                             model->modelForVarSwap.get(),
                                                             control);
                approxValue /= delta;
                
                // deduct this value obtained from ImpliedATMVol
                double basis = VDaxHeston - HestonVarSwapBetwMatFwdMatEnd;
                
                if (CString::equalsIgnoreCase(model->methodChoice,"VAR_SWAP_REL")) {
                    double betaBasis = basis / HestonVarSwapBetwMatFwdMatEnd;
                    value = approxValue * (1 + betaBasis);
                } else {
                    value = approxValue + basis;
                }
                
            } else if (CString::equalsIgnoreCase(model->methodChoice,"ATM_VAR_ABS") || 
                           CString::equalsIgnoreCase(model->methodChoice,"ATM_VAR_REL")) {
                // get implied variance for intervals (0,T) and (0,T+delta), i.e. sigma(t)^2*t
                ATMVolRequest volRequest;
                CVolProcessedBSSP volBS(inst->asset->getProcessedVol(&volRequest));
                double varToFwd = volBS->CalcVar(inst->valueDate, inst->maturityFwd);
                double varToEnd = volBS->CalcVar(inst->valueDate, maturityEnd);

                double approxValue = (varToEnd-varToFwd) / delta;

                double impliedATMVolMatEnd = 
                    model->basisApproxModel->impliedSqVol(inst->valueDate, inst->valueDate, maturityEnd);
                double impliedATMVolMatFwd = 
                    model->basisApproxModel->impliedSqVol(inst->valueDate, inst->valueDate, inst->maturityFwd);
                
                double basis = VDaxHeston - ( impliedATMVolMatEnd*(T+delta) - impliedATMVolMatFwd*T ) / delta;
                
                if (CString::equalsIgnoreCase(model->methodChoice,"ATM_VAR_REL")) {
                    double betaBasis = basis / HestonVarSwapBetwMatFwdMatEnd; 
                    value = approxValue * (1 + betaBasis);
                } else {
                    value = approxValue + basis;
                }
            } 

            double pv = inst->instSettle->pv(inst->valueDate,
                                             inst->maturityFwd,
                                             inst->discount.get(), 
                                             inst->asset.get());
            value *= inst->notional * pv;
            
            results->storePrice(value, inst->discount->getCcy());
            
            if (control && control->isPricing()) {
                inst->recordOutputRequests(control, results, value);

                OutputRequest* request = control->requestsOutput(OutputRequest::DBG);
                if (request) {
                    results->storeRequestResult(request, VDaxHeston, "VDAX_HESTON");
                }
            }
        }
        catch (exception& e){
                throw ModelException(e, method);
        }
    }    
};

CClassConstSP const VDaxModel::TYPE = 
    CClass::registerClassLoadMethod("VDaxModel", typeid(VDaxModel),
    VDaxModel::load);

CClassConstSP const VDaxModel::VDaxBasisApprox::TYPE = 
    CClass::registerClassLoadMethod("VDaxModel::VDaxBasisApprox", typeid(VDaxModel::VDaxBasisApprox), 
    VDaxModel::VDaxBasisApprox::load);

CClassConstSP const VDax::TYPE = 
    CClass::registerClassLoadMethod("VDax", typeid(VDax),
    VDax::load);

CClassConstSP const VDaxModel::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("VDaxModel::IIntoProduct",
                                    typeid(VDaxModel::IIntoProduct), 
                                    VDaxModel::loadIntoProduct);

// at the end due to use of VDaxPricer
VDaxModel::IProduct* VDax::createProduct(VDaxModel* model) const
{
    return new VDaxPricer(this);
}

bool VDaxLoad() {
    return (VDax::TYPE != 0);
}

/** "bread & butter" that can be captured in Pyramid using current IMS */
class VDaxII: public VDaxModel {
public:
    static CClassConstSP const TYPE;
    
private:
    VDaxII():VDaxModel() {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(VDaxII, clazz);
        SUPERCLASS(VDaxModel);
        EMPTY_SHELL_METHOD(defaultVDaxII);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultVDaxII(){
        return new VDaxII();
    }
};

CClassConstSP const VDaxII::TYPE = 
CClass::registerClassLoadMethod("VDaxII", typeid(VDaxII), load);

/** "bread & butter" that can be captured in Pyramid using current IMS */
class VDaxCF: public VDaxModel {
public:
    static CClassConstSP const TYPE;
    
private:
    VDaxCF():VDaxModel() {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(VDaxCF, clazz);
        SUPERCLASS(VDaxModel);
        EMPTY_SHELL_METHOD(defaultVDaxCF);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultVDaxCF(){
        return new VDaxCF();
    }
};

CClassConstSP const VDaxCF::TYPE = 
CClass::registerClassLoadMethod("VDaxCF", typeid(VDaxCF), load);


DRLIB_END_NAMESPACE
