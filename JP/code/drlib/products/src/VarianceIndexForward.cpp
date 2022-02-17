//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VarianceIndexForward.cpp
//
//   Description : VarianceIndexForward contract: 
//
//   Author      : Jay Blumenstein
//
//   Date        : 19 Apr 2004
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/VolSV.hpp"
#include "edginc/VolVarSwap.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/AssetUtil.hpp"


DRLIB_BEGIN_NAMESPACE

class VarSwapImpliedIntegration : public CModel {
public:
    static CClassConstSP const TYPE;
    
    friend class VarianceIndexForwardPricer;

    virtual void validatePop2Object() {
        static const string method("VarSwapImpliedIntegration::validatePop2Object");

    }
        
    /** the class that the product must be able to create */
    class IProduct{
    public:
        virtual void price(VarSwapImpliedIntegration* model,
                           Control*   control, 
                           Results*   results) const = 0;
        virtual ~IProduct(){};
    };

    /** interface that the instrument must implement */
    class IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(VarSwapImpliedIntegration* model) const = 0;
    };

    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results) {
        static const string method = "VarSwapImpliedIntegration::Price";
        if (!IIntoProduct::TYPE->isInstance(instrument)){
            throw ModelException(method, "Instrument of type "+
                                 instrument->getClass()->getName() +
                                 " does not support VarSwapImpliedIntegration::IntoProduct");
        }
        IProduct* product = 0;
        try {
            if (instrument->priceDeadInstrument(control, results)) {
                return; // done for a dead instrument
            }

            // cast to VarSwapImpliedIntegration::IIntoProduct
            IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
            // create the product
            product = intoProd.createProduct(this);
            product->price(this, control, results);
        } 
        catch (exception& e) {
            delete product;
            throw ModelException(e, method);
        }
        delete product;
    }
     
    MarketObjectSP GetMarket(const MarketData*    market,
                             const string&        name,
                             const CClassConstSP& type) const{
        static const string method = "VarSwapImpliedIntegration::GetMarket";
        try {
            return ImpIntModel->GetMarket(market, name, type);
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
        return riskMappingIrrelevant;
    }


    // for VarSwapImpliedIntegration::IIntoProduct
    static void IntoProduct_load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(Model::IModelIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    // registration, invoked when class is 'loaded'
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(VarSwapImpliedIntegration, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultVarSwapImpliedIntegration);
        FIELD(ImpIntModel,"Implied Integration");
    }

    // for VarSwapImpliedIntegration::IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(VarSwapImpliedIntegration::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultVarSwapImpliedIntegration(){
        return new VarSwapImpliedIntegration();
    }


private:
    // constructor
    VarSwapImpliedIntegration():CModel(TYPE) {}

    // registered fields
    ImpliedIntegrationSP    ImpIntModel;        // ImpliedIntegration
};

/*********************************************************************/

class VarianceIndexForward : public Generic1Factor, 
             virtual public VarSwapImpliedIntegration::IIntoProduct,
             virtual public ISensitiveStrikes,
             virtual public LastSensDate {
public:
    static CClassConstSP const TYPE; 
    
    friend class VarSwapImpliedIntegration;

    static const DateTime offsetDateByTenor(const string& tenorPeriod, const DateTime& aDate)
    {
        MaturityPeriod periodRef(tenorPeriod);            
        return DateTime(periodRef.toDate(aDate).getDate(), aDate.getTime());          
    }

    // do some validation at construction -- not mandatory
    virtual void validatePop2Object() {
        static const string method("VarianceIndexForward::validatePop2Object");
        try {
            if (instSettle->isPhysical() || instSettle->isMargin()) {
                throw ModelException(method, "Only cash settlement is allowed");
            }
            if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
                throw ModelException(method, "forward variance contract can't be ccy struck");
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
        
    // do some asset specific validation -- mandatory since VarianceIndexForward : Generic1Factor : CInstrument
    virtual void Validate() {
        static const string method("VarianceIndexForward::validate");
        try {
            validate(); // validation of Generic1Factor
            if (FXAsset::TYPE->isInstance(asset.get())){
                throw ModelException(method, "FX underlying not allowed");
            }
            if (frontMaturity.isGreaterOrEqual(backMaturity)) {
                throw ModelException(method, "front maturity " + frontMaturity.toString() + 
                    " cannot be after back maturity " + backMaturity.toString() + ".");
            }
            if (forwardDate.isGreaterOrEqual(frontMaturity)) {
                throw ModelException(method, "variance index forward date " + forwardDate.toString() + 
                    " cannot be after front maturity " + frontMaturity.toString() + ".");
            }

            DateTime f_plus_tenor = offsetDateByTenor(tenorPeriod, forwardDate);

            if (f_plus_tenor.isGreaterOrEqual(backMaturity)) {
                throw ModelException(method, "variance index forward date + tenor = " + f_plus_tenor.toString() + 
                    " cannot be after \nbackMaturity = " + backMaturity.toString() + ".");
            }
             
        }    
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // indicates whether VEGA_MATRIX is sensible for this instrument
    virtual bool avoidVegaMatrix(const IModel* model) {
        // similar stuff as in or delegated to VolVarSwap.cpp?
        // todo
        bool temp=0;
        return temp;
    }
        
    // returns all strikes on the vol surface to which this instrument is sensitive
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model) {
        // similar stuff as in or delegated to VolVarSwap.cpp?
        // todo
        DoubleArraySP temp;
        return temp;
    }

    DateTime endDate(const Sensitivity* sensControl) const {
        // todo
        DateTime temp;
        return temp;
    }

    // price a dead instrument until settlement
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const {
        static const string method = "VarianceIndexForward::priceDeadInstrument";
        try {
            bool deadinstrument = false;
            // compare forwardDate + settlement Date vs valueDate 
            DateTime settleDate  = instSettle->settles(forwardDate, asset.get());

            if (forwardDate.isGreater(valueDate)) {
                // nthg happens, i.e. normal price
            } else if (valueDate.isGreaterOrEqual(settleDate)) {
                // zero value
                results->storePrice(0.0, discount->getCcy());
                
                deadinstrument = true;   
                
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, 0.0);
                }
            } else if (valueDate.isGreaterOrEqual(forwardDate)&&valueDate.isLess(settleDate)) {
                // check that there is user input
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
        static const string method = "VarianceIndexForward::recordOutputRequests";
        try {
            OutputRequest* request = NULL;
        
                request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
                if (request) {
                    DateTime paymentDate = instSettle->settles(forwardDate, asset.get());
                    DateTimeArray date(1, paymentDate);
                    OutputRequestUtil::recordPaymentDates(control,results,&date); 
                }

                request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
                if (request) {
                    if (valueDate.isGreaterOrEqual(forwardDate)) {
                        DateTime paymentDate = instSettle->settles(forwardDate, asset.get());
                        CashFlow cf(paymentDate, value);
                        CashFlowArray cfl(1, cf);
                        OutputRequestUtil::recordKnownCashflows(control,
                                                                results,
                                                                discount->getCcy(),
                                                                &cfl);  
                    }
                }
        }
        catch (exception& e){
                throw ModelException(e, method);
        }
    }
    
    // registration, invoked when class is 'loaded'
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VarianceIndexForward, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(VarSwapImpliedIntegration::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultVarianceIndexForward);
        FIELD(forwardDate, "maturity of forward contract");
        FIELD(frontMaturity, "first maturity date which is on or after forwardDate");
        FIELD(backMaturity, "first date which is on or after forwardDate + tenorPeriod");

        FIELD(tenorPeriod, "length of tenor");
        FIELD_MAKE_OPTIONAL(tenorPeriod);
        FIELD(priceAtFwd,"price of variance index when forward contract expires");
        FIELD_MAKE_OPTIONAL(priceAtFwd);
        FIELD(noDivAdj, "if true, no adjustment in realized variance computation for dividends");
        FIELD(observationsPerYear, "number of sample points per year");
        FIELD(numCalendarDaysPerYear, "number of calendar days per year");
        FIELD(strikeVol, "strike");
        FIELD(dontScaleByStrike, "TRUE = don't divide by 2K");
        FIELD_MAKE_OPTIONAL(dontScaleByStrike);
        

        
    }
    
    static IObject* defaultVarianceIndexForward(){
        return new VarianceIndexForward();
    }
    
    VarSwapImpliedIntegration::IProduct* createProduct(VarSwapImpliedIntegration* model) const;
    
private:
    friend class VarianceIndexForwardPricer;
    
    // registered fields
    DateTime        forwardDate;        // maturity of forward contract 
    DateTime        frontMaturity;        // first maturity date which is on or after forwardDate
    DateTime        backMaturity;         // first date which is on or after forwardDate + tenorPeriod
    
    double          priceAtFwd;
    string          tenorPeriod;
    bool            noDivAdj;
    
    double          strikeVol;
    bool            dontScaleByStrike;

    int observationsPerYear;
    int numCalendarDaysPerYear;          // in termsheet, variance is converted to vol using calendar time

    // constructor
    VarianceIndexForward(): Generic1Factor(TYPE), 
        tenorPeriod("30D"), dontScaleByStrike(false) {}  
    
}; // end of class VarianceIndexForward


/*********************************************************************/

class VarianceIndexForwardPricer: virtual public VarSwapImpliedIntegration::IProduct{
private:
    const VarianceIndexForward* inst; // a reference
public:
    VarianceIndexForwardPricer(const VarianceIndexForward* instrument): inst(instrument){}
    
    
    void price(VarSwapImpliedIntegration*       model,
        Control*         control, 
        CResults*        results) const{
        static const string method = "VarianceIndexForwardPricer::price";
        try {
            // MaturityTimePeriod        
            //int tenorTimeInt = DateTime::timeConvert(inst->forwardDate.getTime());            
            int tenorTimeInt = inst->forwardDate.getTime();            
            MaturityPeriod periodRef(inst->tenorPeriod);            
            DateTime maturityEnd = DateTime(periodRef.toDate(inst->forwardDate).getDate(), tenorTimeInt);
            
            // ---v------f-------T1----f + 30---------T2-----
            
            // total variance from f to T1
            HolidayConstSP hols(AssetUtil::getHoliday(inst->asset.get()));

            int busdays_f_T1 = hols->businessDaysDiff(inst->forwardDate, inst->frontMaturity);
            int busdays_f_T2 = hols->businessDaysDiff(inst->forwardDate, inst->backMaturity);

            ATMVolRequest volReq;
            CVolProcessedBSSP vol(inst->asset->getProcessedVol(&volReq));
            
            double yrs_f_T1 = vol->calcTradingTime(inst->forwardDate, inst->frontMaturity);
            double yrs_f_T2 = vol->calcTradingTime(inst->forwardDate, inst->backMaturity);

            double sigmaSquaredT1 = VarianceSwapUtil::priceFwdStartingVarSwap(
                inst->asset.get(), 
                inst->discount.get(), 
                1.0, /* notional */
                1.0, /* scale factor */
                true, /* don't scale by strike */
                false, /* no div adj */
                inst->valueDate, 
                inst->forwardDate, 
                inst->frontMaturity, /* maturity */
                0.0, /* strike vol */
                inst->observationsPerYear, 
                busdays_f_T1,
                model->ImpIntModel.get(),
                control); 
            
            double varT1 = sigmaSquaredT1 * yrs_f_T1;
                
            double sigmaSquaredT2 = VarianceSwapUtil::priceFwdStartingVarSwap(
                inst->asset.get(), 
                inst->discount.get(), 
                1.0, /* notional */
                1.0, /* scale factor */
                true, /* don't scale by strike */
                false, /* no div adj */
                inst->valueDate, 
                inst->forwardDate, 
                inst->frontMaturity, /* maturity */
                0.0, /* strike vol */
                inst->observationsPerYear, 
                busdays_f_T2,
                model->ImpIntModel.get(),
                control); 
            
            double varT2 = sigmaSquaredT2 * yrs_f_T2;
            double T1 = (double) inst->frontMaturity.getDate();
            double T2 = (double) inst->backMaturity.getDate();
            double Tm = (double) inst->offsetDateByTenor(inst->tenorPeriod, inst->forwardDate).getDate();
            double numDaysInBenchmark = Tm - (double)inst->forwardDate.getDate();
            
            double varBenchmark = (T2 - Tm)/(T2 - T1) * varT1 + (Tm - T1)/(T2 - T1) * varT2;
            double volSquaredBenchmark = varBenchmark * inst->numCalendarDaysPerYear / numDaysInBenchmark;
            
            double pv = inst->instSettle->pv(inst->valueDate,
                inst->forwardDate,
                inst->discount.get(), 
                inst->asset.get());
            
            // note the VarianceIndexForward contract is structured as Var Swaps
            double forwardContractValue = VarianceSwapUtil::priceVarSwapSimple(
                        volSquaredBenchmark, inst->strikeVol, inst->notional, 100.0, inst->dontScaleByStrike);

            double value = pv * forwardContractValue;
            
            results->storePrice(value, inst->discount->getCcy());
            
            if (control && control->isPricing()) {
                inst->recordOutputRequests(control, results, value);

                // output to debug packet
                OutputNameConstSP out1(new OutputName("varT1"));
                results->storeGreek(CDoubleSP(CDouble::create(varT1)), Results::DEBUG_PACKET, out1);

                OutputNameConstSP out2(new OutputName("varT2"));
                results->storeGreek(CDoubleSP(CDouble::create(varT2)), Results::DEBUG_PACKET, out2);

                OutputNameConstSP out3(new OutputName("varBenchmark"));
                results->storeGreek(CDoubleSP(CDouble::create(varBenchmark)), Results::DEBUG_PACKET, out3);
            }
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }    
};

/*********************************************************************/

// additional stuff -- however, it also compiles if this stuff is not here ... 
CClassConstSP const VarSwapImpliedIntegration::TYPE = 
CClass::registerClassLoadMethod("VarSwapImpliedIntegration", typeid(VarSwapImpliedIntegration),
                                VarSwapImpliedIntegration::load);

CClassConstSP const VarianceIndexForward::TYPE = 
CClass::registerClassLoadMethod("VarianceIndexForward", typeid(VarianceIndexForward),
                                VarianceIndexForward::load);

CClassConstSP const VarSwapImpliedIntegration::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("VarSwapImpliedIntegration::IIntoProduct",
                                    typeid(VarSwapImpliedIntegration::IIntoProduct), 
                                    VarSwapImpliedIntegration::loadIntoProduct);

// at the end due to use of VarianceIndexForwardPricer
VarSwapImpliedIntegration::IProduct* VarianceIndexForward::createProduct(VarSwapImpliedIntegration* model) const
{
    return new VarianceIndexForwardPricer(this);
}

bool VarianceIndexForwardLoad() {
    return (VarianceIndexForward::TYPE != 0);
}

DRLIB_END_NAMESPACE
