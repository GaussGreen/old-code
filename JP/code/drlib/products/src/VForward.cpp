//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VForward.hpp
//
//   Description : Implement Volatility Index (VIX) Forward with VAsset
//
//   Author      : xiaolan zhang
//
//   Date        : 22 Aug 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/VAsset.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Schedule.hpp"

DRLIB_BEGIN_NAMESPACE

/** Instrument */
class VForward:   public Generic1Factor, 
                    //public CInstrument,
                    virtual public LastSensDate,
                    virtual public ISensitiveStrikes,
                    virtual public VIXFModel::IIntoProduct{  //in order to use var swap pricer for BS
public:
    static CClassConstSP const TYPE;
    static IObject* defaultVForward();

    /** instrument validation */
    virtual void Validate();

    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** Rolls the value date and sets initial spot if rolling over start date */
    virtual bool sensShift(Theta* shift);

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    virtual DateTime getValueDate() const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    /** The following fts are due the addition of the second model to pricer Vanilla Var swap in Close Form*/
    /** Implementation of VIXFModel::IntoProduct interface */
    virtual VIXFModel::IProduct* createProduct(const VIXFModel* model) const;

protected:
    friend class VixForwardProd;

    VForward();
    VForward(const VForward& rhs);
    VForward& operator=(const VForward& rhs);

    static void load(CClassSP& clazz);

private:

    void requests(Control* control, CResults* results) const;

    // just for pricing on maturity date, transiend fields
    double spotAtMat; 

    // **** exported fields ****
    ScheduleSP             exerciseSchedule; //maturity and strike of forward contract
    string                 assetHistorySource;  //!< Designates asset history source

};

typedef smartPtr<VForward> VForwardSP;

void VForward::Validate() {
    static const string method = "VForward::Validate";
    try {
        VAsset* vAsset = dynamic_cast<VAsset*>(asset.get());
        if (!vAsset){
            throw ModelException(method, "failed to get VAsset!");
        }

        // can't get exercise schedule from Market - fail if it is NULL
        if (!exerciseSchedule) {
            throw ModelException(method, "Exercise schedule is NULL");
        }

        // check that we have at least one entry in the exercise schedule
        int numDates = exerciseSchedule->length();
        if (numDates != 1) {
            throw ModelException(method, "Exercise schedule must have one date and strike");
        }

        if (!instSettle) {
            throw ModelException(method, "Instrument settlement is NULL");
        }

        // check that settlement is cash
        if (instSettle->isPhysical() || instSettle->isMargin()) {
            throw ModelException(method,
                                "Only cash settlement is allowed");
        }

    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** Instrument:: Get the asset and discount market data */
void VForward::GetMarket(const IModel*          model, 
                               const CMarketDataSP    market)
{
    CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                               discount, asset);

    const VIXFModel* vsModel = dynamic_cast<const VIXFModel*>(model);

    if (!vsModel) {
        throw ModelException("Model must be of type VIXFModel for Volatility Index Forward.");
    }

    discount.getData(vsModel->varSwapModel.get(), market);

    if (!!instSettle)
        instSettle->getMarket(model, market.get());

    //populate transient fields
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    DateTime maturity = exerciseSchedule->lastDate();

    if (getValueDate().isGreater(maturity) || getValueDate().equals(maturity, false) ){
        spotAtMat = vAsset->getFixing (maturity, assetHistorySource); 

        //if it's the same date, and if the assethistory is zero, then, use spot
        if ((getValueDate().equals(maturity, false)) && Maths::isZero(spotAtMat)){
            spotAtMat = vAsset->getMarketSpot();
        }
    }else{
        spotAtMat = 0.0; 
    }

}

/** when to stop tweaking */
DateTime VForward::endDate(const Sensitivity* sensControl) const {
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    DateTime maturity = exerciseSchedule->lastDate();

    return vAsset->getTenor()->toDate(maturity);
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool VForward::sensShift(Theta* shift)
{    
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    DateTime valueDate = getValueDate(); 
    DateTime maturity = exerciseSchedule->lastDate();
    DateTime aDate = shift->rollDate(valueDate);

    if ((aDate >= maturity && valueDate < maturity ) ||
        (valueDate == maturity && Maths::isZero(spotAtMat))) {
        spotAtMat = asset->getThetaSpotOnDate(shift, maturity);
    }
  
    return true;
};

/** price a dead instrument until settlement*/
bool VForward::priceDeadInstrument(CControl* control, CResults* results) const {
    static const string method = "VForward::priceDeadInstrument";
    try {
        //to review settle????????
        const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());

        bool deadinstrument = false;
        double value;
        // compare maturity + settlement Date vs valueDate 
        DateTime valueDate = getValueDate();
        DateTime maturity = exerciseSchedule->lastDate();
        DateTime settleDate  = instSettle->settles(maturity, asset.get());

        if (maturity.isGreater(valueDate)) { //valueDate < maturity 
            // nthg happens, i.e. normal price
        } else {
            if (valueDate.isGreater(settleDate)) {
                // zero value
                value = 0.0;                
            } else{ // maturityFut <= valueDate <= settleDate
                // VIX as input has already in vol term ( already X 100.0, ex: 20 means 20% )
                value = spotAtMat; 
                if (!Maths::isPositive(value)){
                    throw ModelException("Historic Sample on Asset (" + vAsset->getName() + ") at maturity " + maturity.toString() 
                        + " is missing.");
                }
                double strike = exerciseSchedule->lastValue();
                value -= strike;
                // pv from settlement to today
                value *= discount->pv(valueDate, settleDate)*notional;
            }
            deadinstrument = true;               
            results->storePrice(value, discount->getCcy());  
            if (control && control->isPricing()) {
                requests(control, results);
            }
        }                       
        return deadinstrument;
    }
    catch (exception& e){
            throw ModelException(e, method);
    }
}

void VForward::requests(Control* control, CResults* results) const {
    static const string method = "VixForwardProd::requests";
    try {
        DateTime maturity = exerciseSchedule->lastDate();
        double strike = exerciseSchedule->lastValue();
        DateTime settleDate  = instSettle->settles(maturity, asset.get());
        const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());

        OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(new DateTimeArray(1, settleDate));
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
        }

        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            if (getValueDate().isGreaterOrEqual(maturity)) {

                CashFlowArraySP payoff(new CashFlowArray(1));
                (*payoff)[0].date   = maturity;
                (*payoff)[0].amount = notional * (spotAtMat - strike);
//                (*payoff)[0].amount = (vAsset->getFixing (maturity) - strike);

                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        payoff.get()); 
            }
        }

        if (!(getValueDate().isGreater(maturity))) {
   
            // FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                results,
                maturity,
                getValueDate(),
                vAsset);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//the following 2 ft are copied from VolVarSwap.cpp
/** indicates whether VEGA_MATRIX is sensible for this instrument*/
//?? to ask, is this good?? VIXFModel->VarSwapModel
bool VForward::avoidVegaMatrix(const IModel* model) {
    return false;   
}

/** returns all strikes on the vol surface to which this instrument is sensitive */
DoubleArraySP VForward::getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model) {
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    const VIXFModel* m = dynamic_cast<const VIXFModel*>(model);

    if (!m) {
        throw ModelException("Model must be of type VIXFModel for Volatility Index Forward.");
    }
    return vAsset->sensitiveStrikes(outputName, m);
}

DateTime VForward::getValueDate() const
{
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    DateTime valueDate = vAsset->getAlgorithm()->model->getValueDate();
    return valueDate;
}

/** Returns the name of the instrument's discount currency. */
string VForward::discountYieldCurveName() const {
    return discount.getName();
}

VForward::VForward(): Generic1Factor(TYPE),
                assetHistorySource(IMarketObservable::DEFAULT_SOURCE){
    instSettle = InstrumentSettlementSP(new CashSettlePeriod(0));
};

/** Invoked when Class is 'loaded' */
void VForward::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VForward, clazz);
    SUPERCLASS(Generic1Factor);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(ISensitiveStrikes);

    IMPLEMENTS(VIXFModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultVForward);

    // transient fields
    FIELD(spotAtMat, "spot at maturity");
    FIELD_MAKE_TRANSIENT(spotAtMat);

    // input Fields
    //maturity and strike for forward contract
    FIELD(exerciseSchedule, "Matuirty and strike, one element only");

    FIELD(assetHistorySource, "Asset history source");
    FIELD_MAKE_OPTIONAL(assetHistorySource);
}

IObject* VForward::defaultVForward(){
    return new VForward();
}

/** Type registration*/
CClassConstSP const VForward::TYPE = CClass::registerClassLoadMethod(
    "VForward", typeid(VForward), VForward::load);


//////////////////////////////////////////////////////////////////////////////////

/** VIXFModel : two models in order to call Var Swap BS and SV
*/
class VixForwardProd: virtual public VIXFModel::IProduct {
public:
    
    VixForwardProd(const VForward* inst): inst(inst){}
    
    /** This is the method responsible for calling vAsset->fwdValue
        and output additional components (BS, SV) 
    */
    void price(VIXFModel*       model,
               Control*      control, 
               CResults*     results) const;

private:
    const VForward* inst; // a reference

};

/** This is the method responsible for calling vAsset->fwdValue
    and output additional components (BS, SV) 
*/
void VixForwardProd::price(VIXFModel*       model,               
						  Control*          control, 
                          CResults*         results) const {
    static const string method = "VixForwardProd::price";
    try {

       /** BS_SV : BS + convexity + basis
        SV: E_t[sqrt(Fwd Var(T,T+ tau))] + basis
        BS_SV_NO_BASIS: BS + convexity
        SV_NO_BASIS: E_t[sqrt(Fwd Var(T,T+ tau))]
        VSW_BS: Var Swap using BS
        VSW_SV: Var Swap using SVJ
        */

        // there can be only one ...
        DateTime maturity = inst->exerciseSchedule->lastDate();
        double   strike  = inst->exerciseSchedule->lastValue();
        DateTime valueDate = inst->getValueDate();

        const VAsset* vAsset = dynamic_cast<const VAsset*>((inst->asset.get()));

        CControlSP ctrl(copy(control));
        vAsset->getAlgorithm()->setControl(ctrl);

        double totalValue = 100.0* vAsset->fwdValue(maturity);

        // calculate the discount factor back to value date
        double discFactor = inst->discount->pv(valueDate, maturity);

        // discounting 
        double settlementPV = inst->instSettle->pvAdjust(maturity, 
                                                      inst->discount.get(), 
                                                      inst->asset.get());
        // premium for fwd contract
        totalValue = settlementPV*discFactor* inst->notional*(totalValue - strike);

        results->storePrice(totalValue, inst->discount->getCcy());

        /** vector "out" contains BS, SV components, only has values after calling vAsset->fwdValue*/
        /** record "out" into results based on the control */
        vAsset->recordExtraOutput(control, results);

        if (control && control->isPricing()) {
            inst->requests(control, results);
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
};

//////////////////////////////////////////////////////////////////////////////////

/** Implementation of VIXFModel::IntoProduct interface */
VIXFModel::IProduct* VForward::createProduct(const VIXFModel* model) const {
    return new VixForwardProd(this);
}

bool VForwardLoad()
{
    return (VForward::TYPE != 0);
}

DRLIB_END_NAMESPACE
