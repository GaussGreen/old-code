#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Bond.hpp"
#include "edginc/FD2DSolver.hpp"
#include "edginc/FD2DLN.hpp"
#include "edginc/FDModel.hpp"
#include "edginc/FDProduct.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/PhysicalSettlement.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Maths.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/GeneralAsset.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Asset.hpp"
#include "edginc/CreditCurve.hpp"
#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/TreeSliceOperExt.hpp"

DRLIB_BEGIN_NAMESPACE


static CFieldConstSP creditCurveField;

/**
 * Opera instrument
 *
 * <H3>Description</H3>
 *
 * The Opera contract is a convertible bond that can be converted in two different shares.
 * 
 * The holder owns a bond that can be converted on certain dates into a number a shares of
 * two different stocks .
 *
 * <H3>Features</H3>
 * The following features are available:
 * - conversion
 * - call
 * - put
 */

/** The instrument class */
class Opera : public CInstrument,
              public virtual FDModel::IIntoProduct,
              public virtual Theta::IShift,
              public virtual LastSensDate {

protected:
    Opera(CClassConstSP clazz);
    
    /**------------------------------------------------------------------------*/
    /**--------------------------------- Fields -------------------------------*/
    /**------------------------------------------------------------------------*/
        
    /*--------------------- Required fields -----------------------*/

    /** settlement for the put/call cash flows */
    InstrumentSettlementSP instSettle;
    
    /** Underlying bond */
    // the bond
    BondSP                 bond;
    // discount curve of the bond currency
    YieldCurveWrapper      discount;
   
    /** The assets in which the bond can be converted */
    IMultiMarketFactorsSP  assets;
    
    /*---------------------- Optional fields ---------------------*/

    /** the credit spread for the bond */
    CreditCurveWrapper     creditSpread;  
    
    /** Conversion dates and ratios */
    ScheduleSP convSchedule1; /** first asset  */
    ScheduleSP convSchedule2; /** second asset */

    /** Call */
    ScheduleSP callSchedule1; /** first asset  */
    ScheduleSP callSchedule2; /** second asset */

    /** Put */
    ScheduleSP putSchedule;

    /*---------------------- Internal interface ------------------*/

    /** value date */
    DateTime               valueDate;

    /** settlement for the conversion: physical by default */
    InstrumentSettlementSP convSettle; 
    
    /** risky discount curve for the bond (contains the credit spread) */
    YieldCurveConstSP           riskyDiscount;

    /** Arrays of Schedules for the conversion and the call schedules:
       - prevents to duplicate for each asset (1 and 2)
       - the code can be generalized to any number of 
       underlying stocks */
    ScheduleArraySP        convSchedules;
    ScheduleArraySP        callSchedules;
    
public:
    static CClassConstSP const TYPE; // defined later
    friend class OperaFDProd ;       // the Finite Difference product
    
    //-----------------------------------------------------------------------------------------------------------

    /** Validate schedule:
        - all the values are positive
        - the last date is strictly before the bond maturity date */
    void validateSchedule(ScheduleSP schedule, string name) {
        static const string method("Opera::validateSchedule");
        try {
            // only if the schedule is not empty
            if (schedule->length() > 0) {
                // only positive values
                if (!schedule->isNonNegative()) {
                    throw ModelException(method, name +" schedule must be positive");
                }
                
                // all the dates are before the bond maturity
                DateTime bondEnd     = bond->getMaturityDate();
                DateTime scheduleEnd = schedule->lastDate();
                if (scheduleEnd >= bondEnd) {
                    throw ModelException(method, name + " schedule must end before the bond maturity date.\n"
                                         "Last " + name + " date:" + scheduleEnd.toString() +
                                         "\nBond maturity date:" + bondEnd.toString());
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        } 
    }

    //-----------------------------------------------------------------------------------------------------------

    /** Validation with market data */
    /** just check the things that aren't/cannot be checked in validatePop2Object() */
    virtual void Validate(){
        static const string method = "Opera::Validate";
        try {
            int iAsset;
            
            /** a bond can only be converted into a IGeneralAssets */
            for(iAsset=0; iAsset<assets->numFactors(); iAsset++) {
                IMarketFactorConstSP factor(assets->getFactor(iAsset));
                if (!CAsset::TYPE->isInstance(factor)){
                    throw ModelException(method, "Only underlying of type CAsset are supported:"
                                         " the type of asset " + factor->getName() + " is " 
                                         + factor->getClass()->getName());
                }
            }
            
            /** validate the currency
                the asset must be vanilla */
            for(iAsset=0; iAsset<assets->numFactors(); iAsset++) {
                IMarketFactorConstSP factor(assets->getFactor(iAsset));
                CAssetConstSP        asset(CAssetConstSP::dynamicCast(assets->getFactor(iAsset)));
                string               ccyTreatment = asset->getCcyTreatment();
                
                if (ccyTreatment != CAsset::CCY_TREATMENT_NONE && ccyTreatment != CAsset::CCY_TREATMENT_VANILLA) {
                    throw ModelException(method, "The currency treatment of asset " + asset->getTrueName() +
                                         " (" + ccyTreatment +") must be (N)one");
                }
            }
           
            /** validate the market data for the assets */
            assets->crossValidate(valueDate,
                                  valueDate,
                                  discount.get(),
                                  this);

            /** Conversion, Call and Put schedules validation  
                we need the bond maturity date that can be access only 
                once the bond has been initialized by getMarket*/
            // conversion
            for(iAsset=0; iAsset<convSchedules->size();iAsset++) {
                IMarketFactorConstSP factor(assets->getFactor(iAsset));
                CAssetConstSP        asset(CAssetConstSP::dynamicCast(assets->getFactor(iAsset)));
                validateSchedule((*convSchedules)[iAsset],"Asset" + asset->getTrueName() + "'s conversion");
            }
            // call
            for(iAsset=0; iAsset<callSchedules->size();iAsset++) {
                IMarketFactorConstSP factor(assets->getFactor(iAsset));
                CAssetConstSP        asset(CAssetConstSP::dynamicCast(assets->getFactor(iAsset)));
                validateSchedule((*callSchedules)[iAsset],"Asset" + asset->getTrueName() + "'s call");
            }
            // put
            validateSchedule(putSchedule,"Put");

            // initialize the risky discount curve for the bond
            if (creditSpread.get() != 0) {
                DateTime maturityDate = getMaturityDate();
                riskyDiscount = YieldCurveConstSP::dynamicCast(creditSpread.get()->makeRiskyCurve(*discount.get()));
            }
            else {
                riskyDiscount = discount.getSP();
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        } 
    }
    
    //-----------------------------------------------------------------------------------------------------------

    /** Initializes the internal interface and validates the instrument */
    virtual void validatePop2Object(){
        static const string method("Opera::validatePop2Object");
        try {
            /** exactly 2 underlyings are required */
            if (assets->numFactors() != 2) {
                throw ModelException(method, "The Opera product needs exactly 2 underlyings whereas "+
                                     Format::toString(assets->numFactors()) + "underlyings have been provided");
            }
            
            /** Initialized the internal interface */
            // the conversion settlement is physical
            convSettle = InstrumentSettlementSP(new PhysicalSettlement());

            // conversion schedules
            convSchedules = ScheduleArraySP(new ScheduleArray(2));
            (*convSchedules)[0] = convSchedule1;
            (*convSchedules)[1] = convSchedule2;
            // call schedules
            callSchedules = ScheduleArraySP(new ScheduleArray(2));
            (*callSchedules)[0] = callSchedule1;
            (*callSchedules)[1] = callSchedule2;
        
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    
    void fieldsUpdatedconst(CFieldArray& fields) {
        // initialize the risky discount curve for the bond
        if (creditSpread.get() != 0) {
            DateTime maturityDate = getMaturityDate();
            riskyDiscount = YieldCurveConstSP::dynamicCast(creditSpread.get()->makeRiskyCurve(*discount.get()));
        }
        else {
            riskyDiscount = discount.getSP();
        }
    }

    //-----------------------------------------------------------------------------------------------------------   

    /** returns the clean price of the bond with/without the credit spread curve */
    double getCleanBondValue(const DateTime& date, bool isRisky) const {
        static const string method = "Opera::getCleanBondValue";
        try {
            double bondValue;
         
            // the value is 0 on and after the last cash flow
            if (date>=bond->getMaturityDate()) {
                return 0.0;
            }
            else {
                /*  clean price = PV of future cash flows
                    (i.e. coupons + redemption) */
                YieldCurveConstSP disc=isRisky? riskyDiscount: discount.getSP();
                bondValue = bond->presentValue(date, disc);
                // the line below generates a leak on gcc 3.2
                //bondValue = bond->presentValue(date,
                //               isRisky?riskyDiscount:discount.getSP());
                return bondValue;
            }
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }
    
    //-----------------------------------------------------------------------------------------------------------

    /** returns a double array with the conversion ratios for the assets
        take in account the settlements
        if the bond can not be converted in one asset at that date then the conversion ratio is 0.0 */
    DoubleArray getConversionFactors(const DateTime& convDate) const {
        static const string method = "Opera::getConversionFactors";
        try {
            int iAsset;
            DoubleArray ratios(assets->numFactors());
            DoubleArray settles = getSettlements(convDate);
            
            for(iAsset=0;iAsset<assets->numFactors();iAsset++) {
                if ((*convSchedules)[iAsset]->length() >0) {
                    try {
                        ratios[iAsset]=(*convSchedules)[iAsset]->interpolate(convDate);
                    }
                    catch(exception&) {
                        // in case of failure, the conversion ratio is 0.0
                        ratios[iAsset]=0.0;
                    }
                }
                ratios[iAsset] *= settles[iAsset];
            }
            return ratios;
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }


    //-----------------------------------------------------------------------------------------------------------

    /** returns the put level value taking in account the settlements
        if the bond can not beput om this date then the put value is 0.0 */
    double getPutValue(const DateTime& putDate) const {
        static const string method = "Opera::getputValue";
        try {
            double value;
            
            if (putSchedule->length()>0 && putSchedule->coversDateRange(putDate,putDate,false)) {
                value = putSchedule->interpolate(putDate);
            }
            else {
                // in case of failure, the put value is 0.0
                value = 0.0;
            }
            
            // dummy asset not used for cash settlement
            CAssetConstSP asset(CAssetConstSP::dynamicCast(assets->getFactor(0))); 
            DateTime settle  = instSettle->settles(putDate,asset.get());
            value *= discount->pv(putDate,settle);
            return value;
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    //-----------------------------------------------------------------------------------------------------------

    /** returns the call trigger level on the call date
        if not defined, then strictly negative */
    DoubleArray getCallTriggers(const DateTime& callDate) const {
        static const string method = "Opera::getCallTriggers";
        try {
            int iAsset;
            DoubleArray triggers(assets->numFactors());

            for(iAsset=0;iAsset<assets->numFactors();iAsset++) {
                if (((*callSchedules)[iAsset]->length() >0)&&
                    ((*callSchedules)[iAsset]->coversDateRange(callDate,callDate,false))) {
                    triggers[iAsset]=(*callSchedules)[iAsset]->interpolate(callDate);
                }
                else {
                    // in case of failure, the trigger level is -1.0
                    triggers[iAsset]=-1.0;
                }
            }
            return triggers;
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    
    /** return the pv factor due to settlements on the settle date for each asset */
    DoubleArray getSettlements(const DateTime& settleDate) const {
        static const string method = "Opera::getSettlements";
        try {
            int iAsset;
            DoubleArray settles(assets->numFactors());
            
            for(iAsset=0;iAsset<assets->numFactors();iAsset++) {
                CAssetConstSP  asset(CAssetConstSP::dynamicCast(assets->getFactor(iAsset)));

                DateTime settle = instSettle->settles(settleDate, asset.get());
                settles[iAsset] = discount->pv(settleDate,settle);
            }
            
            return settles;
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }
        
    //-----------------------------------------------------------------------------------------------------------

    /** return the maximum value of the product of 2 arrays */
    static void push_back(DateTimeArraySP array1,
                   const  DateTimeArray& array2) {
        static const string method = "Opera::push_back";
        try {
            int i;
            for(i=0;i<array2.size();i++) {
                array1->push_back(array2[i]);
            }
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }
    

    //-----------------------------------------------------------------------------------------------------------

    /** price a dead instrument until settlement
        returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const{
        static const string method = "Opera::priceDeadInstrument";
        try {
            DateTime maturityDate(getMaturityDate());
            
            // after the bond maturity date
            if (valueDate >= maturityDate) {
                double price = 0.0;
                results->storePrice(price, discount->getCcy());
                recordOutputRequests(control, results, price);
                
                return true;
            }
            else {
                return false;
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    } 

    //-----------------------------------------------------------------------------------------------------------
    
    /** output request */
    void recordOutputRequests(Control* control, 
                              Results* results, 
                              double   optionValue) const {
        static const string method = " recordOutputRequests";
        try {
            if (control->isPricing()) {
                OutputRequest* request = NULL;
                
                // Option price
                if (control->requestsOutput(OutputRequest::OPTION_PRICE, request)) {
                    results->storeRequestResult(request, optionValue);
                }
                
                // Bond price
                if ( control->requestsOutput(OutputRequest::NAKED_BOND_PRICE, request)) {
                    double bondValue = getCleanBondValue(valueDate,true);
                    results->storeRequestResult(request, bondValue);
                }
                
                // Known cashflows i.e. bond cash flows
                if (control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS, request)) {
                    CashFlowArrayConstSP cashFlows = bond->getKnownCashFlows();
                    
                    OutputRequestUtil::recordKnownCashflows(control, results, discount->getCcy(), cashFlows.get());
                }

                // Payment dates i.e. bond cash flow payment dates
                if (control->requestsOutput(OutputRequest::PAYMENT_DATES, request)) {
                  DateTimeArraySP dates = bond->getPaymentDates();
                  
                  OutputRequestUtil::recordPaymentDates(control,results,dates.get());
                }

            }
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    //-----------------------------------------------------------------------------------------------------------

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    /** defined later in the code */
    virtual FDProductSP createProduct(FDModel* model) const;
    
    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel* model, const CMarketDataSP market){
        // populates the value date
        market->GetReferenceDate(valueDate);

        // assets
        assets->getMarket(model,market.get());

        // bond
        bond->getMarket(model, market.get());

        // credit spread
        if (!creditSpread.getName().empty()) {
            creditSpread.getData(model, market);
        }

        // discount curve
        discount.getData(model, market);

        // settlements
        instSettle->getMarket(model, market.get());
    }

    //-----------------------------------------------------------------------------------------------------------

    /** what's today ? */
    virtual DateTime getValueDate() const{
        return valueDate;
    }

    //-----------------------------------------------------------------------------------------------------------

    /** what's the maturity date ? */
    virtual DateTime getMaturityDate() const{
        return bond->getMaturityDate();
    }

    //-----------------------------------------------------------------------------------------------------------
    
    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const{
        // the convertible bond dies when the underlying bond dies
        return getMaturityDate();
    }

    //-----------------------------------------------------------------------------------------------------------

    /** how to do a theta shift */
    virtual bool sensShift(Theta* shift){
        valueDate = shift->rollDate(valueDate);
        
        return true;
    }

    //-----------------------------------------------------------------------------------------------------------

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const {
        return discount.getName();
    }

    //-----------------------------------------------------------------------------------------------------------

private:
    /** for reflection */
    Opera(): CInstrument(TYPE),
             convSchedule1(new Schedule(DateTimeArray(0),DoubleArray(0),"N")),
             convSchedule2(new Schedule(DateTimeArray(0),DoubleArray(0),"N")),
             callSchedule1(new Schedule(DateTimeArray(0),DoubleArray(0),"N")),
             callSchedule2(new Schedule(DateTimeArray(0),DoubleArray(0),"N")),
             putSchedule  (new Schedule(DateTimeArray(0),DoubleArray(0),"N")){
    }

    //-----------------------------------------------------------------------------------------------------------

    Opera(const Opera& rhs);            // not implemented
    Opera& operator=(const Opera& rhs); // not implemented

    //-----------------------------------------------------------------------------------------------------------

    /** default constructor */
    static IObject* defaultOpera() {
    return new Opera();
    }

    //-----------------------------------------------------------------------------------------------------------

    /** object registration */
    static void load(CClassSP& clazz){
        REGISTER(Opera, clazz);
        SUPERCLASS(Instrument);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultOpera);
        // public required fields
        FIELD(discount,     "Discount curve of the instrument currency");
        FIELD(instSettle,          "Instrument settlement when conversion occurs");
        FIELD(assets,              "Assets into which the bond can be converted");
        FIELD(bond,                "Underlying bond of the convertible bond");
        // public optional fields
        FIELD(creditSpread, "Credit spread for the bond");
        FIELD_MAKE_OPTIONAL(creditSpread);              
        FIELD(convSchedule1,       "Conversion dates and ratios for the first asset");
        FIELD_MAKE_OPTIONAL(convSchedule1);        
        FIELD(convSchedule2,       "Conversion dates and ratios for the second asset");
        FIELD_MAKE_OPTIONAL(convSchedule2);
        FIELD(callSchedule1,       "Call dates and ratios for the first asset");
        FIELD_MAKE_OPTIONAL(callSchedule1);
        FIELD(callSchedule2,       "Call dates and ratios for the second asset");
        FIELD_MAKE_OPTIONAL(callSchedule2);
        FIELD(putSchedule,         "Put dates and ratios for the second asset");
        FIELD_MAKE_OPTIONAL(putSchedule);
        // internal interface
        FIELD(valueDate,    "Value date");
        FIELD_MAKE_TRANSIENT(valueDate);
        FIELD(riskyDiscount,"Risky bond discount curve");
        FIELD_MAKE_TRANSIENT(riskyDiscount);
        FIELD(convSettle,          "Conversion settlement: physical");
        FIELD_MAKE_TRANSIENT(convSettle);
        FIELD(convSchedules,       "Array of conversion schedules, one schedule per asset");
        FIELD_MAKE_TRANSIENT(convSchedules);
        FIELD(callSchedules,       "Array of call schedules, one schedule per asset");
        FIELD_MAKE_TRANSIENT(callSchedules);

        // make visible to EAS/spreadsheet
        clazz->setPublic(); 

        // look up field for use on recurse
        creditCurveField = clazz->getDeclaredField("creditSpread");
    }
};

//-----------------------------------------------------------------------------------------------------------

/** for class loading (avoid having header file) */
bool OperaLoad() {
    return (Opera::TYPE != 0);
}

/** smart pointer and const smart pointer classes */
typedef smartConstPtr<Opera> OperaConstSP;
typedef smartPtr<Opera> OperaSP;

//-----------------------------------------------------------------------------------------------------------

/**-----------------------------------------------------------------------**/
/**------------------private class for all tree/FD product----------------**/
/**-----------------------------------------------------------------------**/
class OperaFDProd : public FDProduct,
                    public IFDProductLN{
private:
    OperaConstSP     opera;
    FDProductSP      payoffIndex1;
    FDProductSP      payoffIndex2;

    TreeSliceSP      value;

    /** for the payoff (calculation efficiency) */

    // is conversion, call or put happening on a step date ?
    BoolArray isEvent; 
    // number of assets
    int       nbAssets;
    
public:    

    //-----------------------------------------------------------------------------------------------------------

    OperaFDProd(OperaConstSP opera, FDModelSP model) : 
        FDProduct(model.get()),
        opera(opera)
    {
        static const string method = "OperaFDProd::OperaFDProd";
        try {
            // only the FD2DLn model is supported
            if( ! dynamic_cast< FD2DLN * >( model.get() ) )
                throw ModelException(method, "Only the FD2DLN model is supported by the Opera FD Product");

            // second: create spot payoffs
            payoffIndex1 = model->createProduct( IProdCreatorSP( new
                IndexSpecEQ(
                    opera->assets->getName(0),
                    CAssetWrapper( dynamic_cast< CAsset * >(
                        const_cast< IMarketFactor * >( opera->assets->getFactor(0).get() ) ) ),
                    "" ) ) );
            payoffIndex2 = model->createProduct( IProdCreatorSP( new
                IndexSpecEQ(
                    opera->assets->getName(1),
                    CAssetWrapper( dynamic_cast< CAsset * >(
                        const_cast< IMarketFactor * >( opera->assets->getFactor(1).get() ) ) ),
                    "" ) ) );
        }         
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
   
    //-----------------------------------------------------------------------------------------------------------

    /** start date */
    virtual DateTime getStartDate() const{
        return opera->valueDate;
    }
    
    //-----------------------------------------------------------------------------------------------------------

    /** output prices and any additional results */
    virtual void recordOutput(Control* control, YieldCurveConstSP discount, Results* results){
        // get prices at t=0
        // save price
        double price = model->getPrice0( *value );
        double bondValue = opera->getCleanBondValue(opera->valueDate,true);
        results->storePrice(price + bondValue, discount->getCcy()); 

        opera->recordOutputRequests(control, results, price);
    }

    //-----------------------------------------------------------------------------------------------------------

    /**---------------------------------------------------------------------------**/
    /**---------------------------FD-2 factors engine-----------------------------**/
    /**---------------------------------------------------------------------------**/

    /**************************** initialisation functions *************************/

    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(Control*  control) const {
        static const string method = "OperaFDProd::init";
        try {
            /** customize tree parameters here and set up the FD grid */
            DateTimeArray segDates(2);
            segDates[0] = getStartDate();
            segDates[1] = opera->getMaturityDate(); 

            // create the list of the critical dates
            int i;
            DateTimeArraySP critDates(new DateTimeArray(0));
            
            // adds bond payment dates
            CashFlowArraySP bondCFs = opera->bond->getExAdjCashFlows(opera->valueDate, opera->discount.getSP());
            for (i=0; i<bondCFs->size(); i++) {
                critDates->push_back((*bondCFs)[i].date);
            }

            // adds conversion dates
            for (i=0;i<opera->convSchedules->size();i++) {
                if ((*opera->convSchedules)[i]->length()>0) {
                    Opera::push_back(critDates,(*opera->convSchedules)[i]->getDates());
                }
            }
            // adds call dates
            for (i=0;i<opera->callSchedules->size();i++) {
                if ((*opera->callSchedules)[i]->length()>0) {
                     Opera::push_back(critDates,(*opera->callSchedules)[i]->getDates());
                }
            }
            // adds put dates
            if (opera->putSchedule->length()>0) {
                Opera::push_back(critDates,opera->putSchedule->getDates());
            }
            
            // dividend dates ??????? to do
            //DateTimeArraySP divDates(new DateTimeArray(0));

            // add critical dates
            model->addCritDates( *critDates );

            IntArray density(1,1);
            // prepare timeline set up
            model->initSegments( segDates, density );
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    //-----------------------------------------------------------------------------------------------------------

    /** initialising and setting product variables 
        this is called per pricing call before each pricing */
    virtual void initProd() {
        static const string method = "OperaFDProd::initProd()";
        try {
            // number of assets
            nbAssets = opera->assets->numFactors();           
            
            startDEV( value = model->createSlice() ); 

            // call, put or conversion on a step date?
            int nbSteps = model->getLastStep()+1;
            isEvent = BoolArray(nbSteps,false);
            
            // set up event dates i.e. when there is a possible conversion, call or put
            int iStep, iAsset, iEvent;
            
            IntArray      idxEvent(2*nbAssets+1,0);
            ScheduleArray eventList(2*nbAssets+1);
            
            // list of all the possible events dates
            eventList[0] = opera->putSchedule;
            for (iAsset=0; iAsset<nbAssets; iAsset++) {
                eventList[iAsset+1]            = (*opera->callSchedules)[iAsset];
                eventList[iAsset+1 + nbAssets] = (*opera->convSchedules)[iAsset];
            }

            //to allow go into payoff ft at maturity.
            isEvent[nbSteps -1] = true;
                        
            // loops through the simulation dates to identify the event dates
            const DateTimeArray fdDates =  model->getDates();
            for (iStep=0; iStep<nbSteps; iStep++) {
                for(iEvent=0;iEvent<eventList.size();iEvent++) {
                
                    ScheduleSP event = eventList[iEvent];
                    DateTimeArray eventDates = event->getDates();

                    if (eventDates.size()>0) {
                        // discret set of dates
                        if (event->getInterp() == Schedule::INTERP_NONE) {
                            if (eventDates[idxEvent[iEvent]] ==  fdDates[iStep]) {
                                idxEvent[iEvent] = Maths::min(idxEvent[iEvent]++,eventDates.size()-1);
                                isEvent[iStep] = true;
                            }
                        }
                        // continuous interval
                        else {
                            isEvent[iStep] = isEvent[iStep]|| 
                                ((fdDates[iStep]>= event->firstDate()) && 
                                 (fdDates[iStep]<= event->lastDate()));
                        }
                    }
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    
    //-----------------------------------------------------------------------------------------------------------

    /**************************** payoff functions *********************************************/
    
    /** payoff for backward induction*/
    void payoff( int step, FDProduct::UpdateType type )
    {
        // value date
        DateTime stepDate(model->getDate(step));

        // conversion ratio for each asset
        DoubleArray convFactors = opera->getConversionFactors(stepDate);

        // call triggers for each asset
        DoubleArray callTriggers = opera->getCallTriggers(stepDate);

        // value of the bond: in the option part, the bond is riskless
        double bondValue = opera->getCleanBondValue(stepDate,false);

        // put value
        double putValue = opera->getPutValue(stepDate);

        // first asset
        const TreeSlice & s1 = payoffIndex1->getValue(step);
        // second asset
        const TreeSlice & s2 = payoffIndex2->getValue(step);

        // at maturity then the reference price is 0.0
        if( type == FDProduct::BWD_T )
            *value = 0.;

        // is trigger condition verified?
        // remember the call trigger is -1.0 if not defined
        if( callTriggers[0] >= 0. )
        {
            // then the convertible maybe be called
            *value =
                IF( callTriggers[0] < s1 * convFactors[0] )
                    smin( *value, s1 * convFactors[0] - bondValue )
                ELSE
                    *value
                ENDIF ;
        }
        if( callTriggers[1] >= 0. )
        {
            // then the convertible maybe be called
            *value =
                IF( callTriggers[1] < s2 * convFactors[1] )
                    smin( *value, s2 * convFactors[1] - bondValue )
                ELSE
                    *value
                ENDIF ;
        }

        // final value
        *value =
            smax(
                *value,
                smax(
                    smax( s1 * convFactors[ 0 ], s2 * convFactors[ 1 ] ), // maximum conversion value
                    putValue ) - bondValue );
    }

    //-----------------------------------------------------------------------------------------------------------

    /** payoff at T for backward or value at t=0 for fwd induction
        payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type){
        static const string method = "OperaFDProd::ypdate()";
        try {
            /** generic obscure par related to the FD2D interface 
                more comments are needed to explain the use of each function 
                this function is calculating the spot levels fir the current step
                and call the payoff function */

            if (isEvent[step])
            {
                if( type == FDProduct::BWD_T || type == FDProduct::BWD )
                {
                    // payoff for backward induction
                    payoff( step, type );
                }
                else
                {
                    // payoff for forward induction
                    throw ModelException(method, "Only backward induction is supported by the Opera FD2DLN Product");
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    //-----------------------------------------------------------------------------------------------------------

    /**-----------------------------------------------------------------------------------**/
    /**---------------------- Log normal model -------------------------------------------**/
    /**-----------------------------------------------------------------------------------**/

    /** Volatility level interpolation */
    CVolRequestLNSP getVolInterp(int iAsset) const {
        
        //CVolRequestLNArray getVolInterp(int iAsset) const {
        
        // get strike and maturity date from instrument
        DateTime matDate = opera->getMaturityDate();
        
        double volStrike;
        ScheduleSP convSchedule = (*opera->convSchedules)[iAsset];

        // find the first convertion ratio not null
        bool     found = false;
        DateTime convDate;
        double   ratio ;
        for (int iDate=convSchedule->length()-1; (iDate>=0)&(!found); iDate--) {
            convDate = (convSchedule->getDates())[iDate];
            ratio = convSchedule->interpolate(convDate);
            found = !Maths::isZero(ratio);
        }      
        
        if (found) {
            // for the option part, the bond is riskless
            double bondValue = opera->getCleanBondValue(convDate,false);
            volStrike = bondValue/ratio;
        }
        // if there is no conversion ratio then the vol level has no impact on the price
        // therefore any strike can be chosen
        else {
            CAssetConstSP asset(CAssetConstSP::dynamicCast(opera->assets->getFactor(iAsset)));
            volStrike = asset->getSpot();
        }
        DateTime imntStartDate = opera->getValueDate();
        CVolRequestLNSP  reqarr;
        reqarr = CVolRequestLNSP(new LinearStrikeVolRequest(volStrike, imntStartDate, 
                                                            matDate, false));        
        return reqarr;
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN()const{
        // empty
    }

};

//-----------------------------------------------------------------------------------------------------------

/**************************************************************************/
/************** Last methods of the class Opera ***************************/
/**************************************************************************/    

/** create a fd payoff product */
FDProductSP Opera::createProduct(FDModel* model) const{
    return FDProductSP(new OperaFDProd(OperaConstSP(this), FDModelSP(model)));
};

/*** Register the type */
CClassConstSP const Opera::TYPE = 
    CClass::registerClassLoadMethod("Opera", typeid(Opera),Opera::load);

DRLIB_END_NAMESPACE
