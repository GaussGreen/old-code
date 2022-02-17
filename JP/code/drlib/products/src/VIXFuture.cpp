//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VIXFuture.hpp
//
//   Description : Implement VIX Future with VAsset
//
//   Author      : xiaolan zhang
//
//   Date        : 28 February 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/VAsset.hpp"
#include "edginc/IndexSpecEQ.hpp"
#include "edginc/InstrumentUtil.hpp"

DRLIB_BEGIN_NAMESPACE

/** Instrument */
class VIXFuture:    public CInstrument,
                    // we'll use below two after removing CInstrument
                    //public IndexSpecEQ,
                    //virtual public IInstrument,
                    virtual public LastSensDate,
                    virtual public ISensitiveStrikes,
                    virtual public VIXFModel::IIntoProduct{  //in order to use var swap pricer for BS
public:
    static CClassConstSP const TYPE;
    static IObject* defaultVIXFuture();

    /** instrument validation */
    virtual void Validate();

    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    virtual DateTime getValueDate() const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    /** Implementation of VIXFModel::IntoProduct interface */
    virtual VIXFModel::IProduct* createProduct(const VIXFModel* model) const;

protected:
    friend class VIXFutureProd;

    VIXFuture();
    VIXFuture(const VIXFuture& rhs);
    VIXFuture& operator=(const VIXFuture& rhs);

    static void load(CClassSP& clazz);

private:

    void requests(Control* control, CResults* results) const;


  	/**!!!!************** exported fields in IndexSpecEQ ************************/
    string                  name;
    CAssetWrapper           asset;

    CVolBaseWrapper         vol; // $unregistered
    BorrowCurveSP           borrowCurve; // $unregistered
    InstrumentSettlementSP  settle;

    string                  assetHistorySource;  //!< Designates asset history source

    /****************** transiend fields ************/
    //will be used only when VIXFuture can be used as an underlying
    //AssetHistoryContainerSP  history;

    // to be initialised by child classes when needed
    //string                  ccyMode;
    //YieldCurveWrapper       yc;
  	/*!!!!!************* end exported fields in IndexSpecEQ ************************/

    // **** exported fields ****
    ExpirySP                 maturityFut;	        // maturity of future contract (ie date of cash settlement)
    double                   notional;
    string                   ccyTreatment;
    YieldCurveWrapper        discount;

    // just for pricing on maturity date, transiend fields
    double spotAtMat; 
};

typedef smartPtr<VIXFuture> VIXFutureSP;

void VIXFuture::Validate() {
    static const string routine = "VIXFuture::Validate";
    try {
        VAsset* vAsset = dynamic_cast<VAsset*>(asset.get());
        if (!vAsset){
            throw ModelException(routine, "failed to get VAsset!");
        }

        if (!settle) {
            throw ModelException(routine, "Instrument settlement is NULL");
        }

        // check that settlement is cash
        if (settle->isPhysical() || settle->isMargin()) {
            throw ModelException(routine,
                                "Only cash settlement is allowed");
        }

//        DateTime valueDate = getValueDate();
//        DateTime maturity =  maturityFut->toDate(valueDate);
//        // Validate that firstDate and lastDate are not weekends
//        if (maturity.isWeekend()) {
//            throw ModelException(routine, "Future Maturity (" + maturity.toString() +
//                                         ") can't be on a weekend");
//        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** Instrument:: Get the asset and discount market data */
void VIXFuture::GetMarket(const IModel*          model, 
                               const CMarketDataSP    market)
{
    CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                               discount, asset);

    const VIXFModel* vsModel = dynamic_cast<const VIXFModel*>(model);

    if (!vsModel) {
        throw ModelException("Model must be of type VIXFModel for Volatility Index Future.");
    }

    discount.getData(vsModel->varSwapModel.get(), market);

// we'll remove these when inheriting from IndexSpecEQ
    if (!!settle)
        settle->getMarket(model, market.get());

   // now try and retrieve the asset histories
//    MarketObjectArraySP objs = 
//        market->getAllObservableHistories(model, name, AssetHistory::TYPE);
//    history = AssetHistoryContainerSP(new AssetHistoryContainer(objs)); 
// we'll remove these when inheriting from IndexSpecEQ
// then, will need to change "name" from optional to mandatory

    //populate transient fields
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    DateTime maturity = maturityFut->toDate(getValueDate());

    if (getValueDate().isGreater(maturity) || getValueDate().equals(maturity, false) ){
        spotAtMat = vAsset->getFixing (maturity, assetHistorySource); 

        //if it's the same date, and if the assethistory is zero, then, use spot
        if ((getValueDate().equals(maturity, false)) && Maths::isZero(spotAtMat)){
            spotAtMat = vAsset->getMarketSpot();
        }
    }else{
        spotAtMat = 0.0; 
    }

    // populate transient fields
//    yc = discount;
//    ccyMode = ccyTreatment;
}

DateTime VIXFuture::endDate(const Sensitivity* sensControl) const {
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    DateTime valueDate = getValueDate(); 
    DateTime maturity = maturityFut->toDate(valueDate);

    return vAsset->getTenor()->toDate(maturity);
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool VIXFuture::sensShift(Theta* shift)
{    
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    DateTime valueDate = getValueDate(); 
    DateTime maturity = maturityFut->toDate(valueDate);
    DateTime aDate = shift->rollDate(valueDate);

    if ((aDate >= maturity && valueDate < maturity ) ||
        (valueDate == maturity && Maths::isZero(spotAtMat))) {
        spotAtMat = asset->getThetaSpotOnDate(shift, maturity);
    }

    return true;
};

/** price a dead instrument until settlement*/
bool VIXFuture::priceDeadInstrument(CControl* control, CResults* results) const {
    static const string method = "VIXFuture::priceDeadInstrument";
    try {
        //to review settle????????
        const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());

        bool deadinstrument = false;
        double value;
        // compare maturityFut + settlement Date vs valueDate 
        DateTime valueDate = getValueDate();
        DateTime maturity =  maturityFut->toDate(valueDate);
        DateTime settleDate  = settle->settles(maturity, asset.get());

        if (maturity.isGreater(valueDate)) { //valueDate < maturity 
            // nthg happens, i.e. normal price
        } else {
            if (valueDate.isGreater(settleDate)) {
                // zero value
                value = 0.0;                
            } else{ // maturityFut <= valueDate <= settleDate
                // VIX as input has already in vol term ( already X 100.0, ex: 20 means 20% )
                if (!Maths::isPositive(spotAtMat)){
                    throw ModelException("Historic Sample on Future (" + name + ") maturity " + maturity.toString() 
                        + " is missing.");
                }
                value = spotAtMat *notional;
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

void VIXFuture::requests(Control* control, CResults* results) const {
    static const string method = "VIXFutureProd::requests";
    try {
        DateTime maturity =  maturityFut->toDate(getValueDate());
        DateTime settleDate  = settle->settles(maturity, asset.get());
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
                (*payoff)[0].amount = notional * spotAtMat;

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
bool VIXFuture::avoidVegaMatrix(const IModel* model) {
    return false;   
}

/** returns all strikes on the vol surface to which this instrument is sensitive */
DoubleArraySP VIXFuture::getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model) {
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    const VIXFModel* m = dynamic_cast<const VIXFModel*>(model);

    if (!m) {
        throw ModelException("Model must be of type VIXFModel for Volatility Index Future.");
    }
    return vAsset->sensitiveStrikes(outputName, m);
}

DateTime VIXFuture::getValueDate() const
{
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());
    VAssetAlgorithmSP       dataTemp= vAsset->getAlgorithm();
    VIXFModelConstSP   modelTemp  = dataTemp->model;

    DateTime valueDate = modelTemp->getValueDate();

    //DateTime valueDate = vAsset->getAlgorithm()->model->getValueDate();
    return valueDate;
}


/** Returns the name of the instrument's discount currency. */
string VIXFuture::discountYieldCurveName() const {
    return discount.getName();
}

VIXFuture::VIXFuture(): CInstrument(TYPE), //IndexSpecEQ(TYPE), 
                notional(1.0),
                ccyTreatment("V"),
                name("") ,
                assetHistorySource(IMarketObservable::DEFAULT_SOURCE){
    settle = InstrumentSettlementSP(new CashSettlePeriod(0));
};

/** Invoked when Class is 'loaded' */
void VIXFuture::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VIXFuture, clazz);
    SUPERCLASS(CInstrument);
    //SUPERCLASS(IndexSpecEQ);
    //IMPLEMENTS(IInstrument);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(ISensitiveStrikes);

    IMPLEMENTS(VIXFModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultVIXFuture);

    // !!!!! fields in IndexSpecEQ, to be removed once inherit from IndexSpecEQ
    FIELD(name,"IndexSpec name");
    FIELD_MAKE_OPTIONAL(name); //put optional for now, will need to change if we use IndexSpecEQ

    FIELD(asset, "underlying asset");

//    FIELD(vol, "Volatility");
//    FIELD_MAKE_OPTIONAL(vol);
//    FIELD(borrowCurve, "Borrow Curve");
//    FIELD_MAKE_OPTIONAL(borrowCurve);
    FIELD(settle, "settlement for VIX Future");
    FIELD_MAKE_OPTIONAL(settle);

    FIELD(assetHistorySource, "Asset history source");
    FIELD_MAKE_OPTIONAL(assetHistorySource);

    // transient fields
//    FIELD(history, "Historical values");
//    FIELD_MAKE_TRANSIENT(history);
//    FIELD(ccyMode,"");
//    FIELD_MAKE_TRANSIENT(ccyMode);
//    FIELD(yc,"");
//    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(yc);
    // !!!!! end fields in IndexSpecEQ, to be removed once inherit from IndexSpecEQ

    // transient fields
    FIELD(spotAtMat, "spot at maturity");
    FIELD_MAKE_TRANSIENT(spotAtMat);

    // input Fields
    FIELD(maturityFut, "maturity of future contract");
    FIELD(discount, " ");
    FIELD(notional, "contract notional");
    FIELD_MAKE_OPTIONAL(notional);
    FIELD(ccyTreatment, "ccy treatment");
    FIELD_MAKE_OPTIONAL(ccyTreatment);    
}

IObject* VIXFuture::defaultVIXFuture(){
    return new VIXFuture();
}

/** Type registration*/
CClassConstSP const VIXFuture::TYPE = CClass::registerClassLoadMethod(
    "VIXFuture", typeid(VIXFuture), VIXFuture::load);

//////////////////////////////////////////////////////////////////////////////////

/** two models in order to call Var Swap BS and SV
*/
class VIXFutureProd: virtual public VIXFModel::IProduct {
public:
    
    VIXFutureProd(const VIXFuture* inst): inst(inst){}
    
    /** This is the method responsible for calling vAsset->fwdValue
        and output additional components (BS, SV) 
    */
    void price(VIXFModel*       model,
               Control*      control, 
               CResults*     results) const;

private:
    const VIXFuture* inst; // a reference

};

/** This is the method responsible for calling vAsset->fwdValue
    and output additional components (BS, SV) 
*/
void VIXFutureProd::price(VIXFModel*       model,               
						  Control*          control, 
                          CResults*         results) const {
    static const string method = "VIXFutureProd::price";
    try {
        const VAsset* vAsset = dynamic_cast<const VAsset*>((inst->asset.get()));

        CControlSP ctrl(copy(control));
        vAsset->getAlgorithm()->setControl(ctrl);

        DateTime valueDate = inst->getValueDate();
        DateTime maturity =  inst->maturityFut->toDate(valueDate);

        double totalValue = 100.0* inst->notional * vAsset->fwdValue(maturity);

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
VIXFModel::IProduct* VIXFuture::createProduct(const VIXFModel* model) const {
    return new VIXFutureProd(this);
}

bool VIXFutureLoad()
{
    return (VIXFuture::TYPE != 0);
}

DRLIB_END_NAMESPACE
