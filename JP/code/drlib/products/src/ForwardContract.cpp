//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ForwardContract.cpp
//
//   Description   forward contract
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ForwardContract.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/PDFRequestLNStrike.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/PhysicalDelivery.hpp"


DRLIB_BEGIN_NAMESPACE

// helpers
void CForwardContract::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spread sheet
    REGISTER(CForwardContract, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    IMPLEMENTS(NumericalIntegrationLN::IIntoProduct);
    IMPLEMENTS(VIXFModel::IIntoProduct);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(Theta::IShift);
    EMPTY_SHELL_METHOD(defaultForwardContract);
    FIELD(valueDate,        "valuation Date");
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD(startDate,        "Option start date");
    FIELD_MAKE_OPTIONAL(startDate);
    FIELD(premiumSettle, "Settlement of fwd premium");
    FIELD_MAKE_OPTIONAL(premiumSettle);
    FIELD(fwdStarting, "Is it a fwd starting option");
    FIELD(exerciseSchedule, "Matuirty and strike, one element only");
    FIELD(oneContract,      "Calc price for 1 contract");
    FIELD(notional,         "Option notional");
    FIELD_MAKE_OPTIONAL(notional);
    FIELD(initialSpot,      "Initial spot price");
    FIELD_MAKE_OPTIONAL(initialSpot);
    FIELD(spotAtMaturity,   "underlying spot at maturity date");
    FIELD_MAKE_OPTIONAL(spotAtMaturity);
    FIELD(asset,            "Underlying of option");
    FIELD(discount,         "Discount curve");
    FIELD(ccyTreatment,     "Currency Treatment");
    FIELD(instSettle, "Instrument settlement at maturity");

    // old dividend re-invest (zeroing dividend) flag is back
    FIELD(divReinvest, "use dividend re-invest mode");
    FIELD_MAKE_OPTIONAL(divReinvest);
}

CClassConstSP const CForwardContract::TYPE = CClass::registerClassLoadMethod(
    "ForwardContract", typeid(CForwardContract), load);

bool  CForwardContractLoad() {
    return (CForwardContract::TYPE != 0);
	}

void CForwardContract::Validate()
{
    static const string method = "ForwardContract::Validate";
    try {
        DateTime            matDate;
        int                 numDates;

        if (divReinvest) {
            throw ModelException(method, 
                                 "Forward has reinvested dividends. "
                                 "Check you really want that functionality, "
                                 "if so check with EDG Quant Research");
        }

        if (!oneContract && !fwdStarting) {
            if (!Maths::isPositive(initialSpot)) {
                throw ModelException(method,
                                     "initial spot (" + 
                                     Format::toString(initialSpot) + 
                                     ") <= 0.0");
            }
        }

        if (fwdStarting && oneContract) {
            throw ModelException(method,
                                 "fwd starting contracts must be notional based");

        }
        // can't get exercise schedule from Market - fail if it is NULL
        if (!exerciseSchedule) {
            throw ModelException(method, "Exercise schedule is NULL");
        }

        if (!instSettle) {
            throw ModelException(method, "Instrument settlement is NULL");
        }

        // check that we have at least one entry in the exercise schedule
        numDates = exerciseSchedule->length();
        if (numDates != 1) {
            throw ModelException(method, "Exercise schedule must have one date and strike");
        }

        // asset and discount curve could come from the market - ie. will
        // be NULL after pop2obj. Do not cross validate if either of them
        // is NULL.
        if ( !(!asset) && !(!discount) ) {
            AssetUtil::assetCrossValidate(asset.get(),
                                          fwdStarting,
                                          startDate,
                                          valueDate,
                                          discount,
                                          this);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// initiate GetMarket 
void CForwardContract::GetMarket(const IModel* model, const CMarketDataSP market)
{
    static const string method = "CForwardContract::GetMarket";
    try {
        market->GetReferenceDate(valueDate);
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, discount, asset);

        const VIXFModel* vsModel = dynamic_cast<const VIXFModel*>(model);
        if (vsModel) {
            discount.getData(vsModel->varSwapModel.get(), market);            
        }else{
            discount.getData(model, market);
        }

        instSettle->getMarket(model, market.get());

        if (premiumSettle.get()) {
            premiumSettle->getMarket(model, market.get());
        }

        //validations for volatility INdex forward 
        if (vsModel) {//VIX Forward
            // check that settlement is cash
            if (instSettle->isPhysical() || instSettle->isMargin()) {
                throw ModelException(method,
                                    "Only cash settlement is allowed for Volatility Index Forward");
            }
            if (fwdStarting)  {
                throw ModelException(method,
                                    "fwd starting volatility index contracts aren't allowed.");
            }
            if (!oneContract) {
                throw ModelException(method,
                                    "oneContract has to be TRUE for volatility index contracts!");
            }
        }

        // since the assets have changed, we have to validate the instrument
        validatePop2Object();

        // asset and discount curve must not be null 
        if (!asset) {
            throw ModelException(method, "Asset is NULL.");
        }

        if (!discount) {
            throw ModelException(method, "Discount curve is NULL.");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
/** false for VAsset, True otherwise */
bool CForwardContract::avoidVegaMatrix(const IModel* model){
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());

    if (vAsset) {
        return false;
    }else{
        return true;
    }
}

/** returns all strikes on the vol surface to which this instrument is sensitive */
/** should be called only if asset = VAsset */
DoubleArraySP CForwardContract::getSensitiveStrikes(OutputNameConstSP outputName,
                                            const IModel*      model){

    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());

    if (vAsset) {
        const VIXFModel* m = dynamic_cast<const VIXFModel*>(model);
        if (!m) {
            throw ModelException("Model must be of type VIXFModel for Volatility Index Forward.");
        }
        return vAsset->sensitiveStrikes(outputName, m);
    }else{ //shouldn't be called
        DoubleArraySP temp(new DoubleArray(asset->getSpot()));
        return temp;
    }
}

// constructor
CForwardContract::CForwardContract(): CInstrument(TYPE), notional(0.0),
                                      initialSpot(0.0), spotAtMaturity(0.0),
                                      divReinvest(false) {
    // empty
}

/** private class */
class CForwardContractClosedFormProd: public CClosedFormLN::IProduct{
private:
    const CForwardContract*  instrFwd; // a reference

public:
    CForwardContractClosedFormProd(const CForwardContract* instr): instrFwd(instr){}

    void price(CClosedFormLN*  model,
               Control*        control, 
               CResults*       results) const;
};

void CForwardContractClosedFormProd::price(CClosedFormLN*   model,
                                           Control*        control, 
                                           CResults*       results) const
{
    static const string method = "CForwardContractClosedForm::price";
    try {
        double         discFactor;       
        double         fwdAtStart = 0.0;      // forward price at start date
        DateTime       settlementDate;  // instrument settlement date
        double         premium;         // the fair value 
        double         fwdAtMat;        // forward at maturity
        double         settlementPV;    // df between maturity and settlement

        if (instrFwd->priceDeadInstrument(control, results)){
            return; // dead instrument priced
        }

        // there can be only one ...
        DateTime matDate = instrFwd->exerciseSchedule->lastDate();
        double   strike  = instrFwd->exerciseSchedule->lastValue();

        // get settlement date
        settlementDate = instrFwd->instSettle->settles(matDate, 
                                                       instrFwd->asset.get());

        CAsset::FwdValueAlgorithm algo(instrFwd->divReinvest);
        if (instrFwd->fwdStarting) {
            DateTimeArray fwdDates(1, instrFwd->startDate);
            DoubleArray   fwds(1);
            instrFwd->asset->fwdValue(fwdDates, algo, fwds);
            fwdAtStart = fwds[0];
            strike *= fwdAtStart;
        }

        DateTimeArray fwdDates(1, matDate);
        DoubleArray   fwds(1);
        instrFwd->asset->fwdValue(fwdDates, algo, fwds);
        fwdAtMat = fwds[0];


        // calculate the discount factor back to value date
        discFactor = instrFwd->discount->pv(instrFwd->valueDate, matDate);

        // discounting 
        settlementPV = instrFwd->instSettle->pvAdjust(matDate, 
                                                      instrFwd->discount.get(), 
                                                      instrFwd->asset.get());
        // premium for fwd contract
        premium = settlementPV*discFactor*(fwdAtMat - strike);

        if (!instrFwd->oneContract) {
            // handle fixed notional 
            if (instrFwd->fwdStarting) {
                if (Maths::isZero(fwdAtStart)) {
                    throw ModelException(method, 
                                         "Forward at start is 0.0. Infinite premium.");
                }
                premium *= instrFwd->notional/fwdAtStart;
            }
            else {
                if (Maths::isZero(instrFwd->initialSpot)) {
                    throw ModelException(method, 
                                         "initial Spot is 0.0. Infinite number of contracts.");
                }
                // handle position 
                premium *= instrFwd->notional/instrFwd->initialSpot;
            }
        }
        

        results->storePrice(premium, instrFwd->discount->getCcy());

        // take care of additional outputs
        if (control && control->isPricing()) {
            instrFwd->recordRequests(control, results);

            // DELAY_PRICE
            InstrumentUtil::delayPriceHelper(control,
                                             results,
                                             premium,
                                             instrFwd->valueDate,
                                             instrFwd->discount.get(),
                                             instrFwd->asset.get(),
                                             instrFwd->premiumSettle.get());
        } 
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}
    
/** Rolls the value date and sets initial spot if rolling over start date */
bool CForwardContract::sensShift(Theta* shift)
{    
    DateTime aDate = shift->rollDate(valueDate);

    DateTime matDate = exerciseSchedule->lastDate();

    if ( ( aDate >= matDate && valueDate < matDate ) ||
         ( valueDate == matDate && Maths::isZero(spotAtMaturity)))
        spotAtMaturity = asset->getThetaSpotOnDate(shift, matDate);

    if (fwdStarting && aDate.isGreaterOrEqual(startDate) &&
        startDate.isGreaterOrEqual(valueDate))
    {
        fwdStarting = false;
        initialSpot = asset->getThetaSpotOnDate(shift, startDate);
        exerciseSchedule->scale(initialSpot);
    }

    valueDate = aDate;

    return true;
};

DateTime CForwardContract::getValueDate() const 
{
    return valueDate;
}

/** when to stop tweaking */
DateTime CForwardContract::endDate(const Sensitivity* sensControl) const {
    DateTime end;

    DateTime maturity = exerciseSchedule->lastDate();
    DateTime instEnd  = instSettle->settles(maturity, asset.get());
    DateTime assetEnd = asset->settleDate(maturity); // this is for spot asset

    if (const LastSensDate* ptr = dynamic_cast<const LastSensDate*>(asset.get())){
        end = ptr->endDate(sensControl); // this is more general
        int daysDiff = end.daysDiff(valueDate);
        end = maturity.rollDate(daysDiff);
        assetEnd = assetEnd > end? assetEnd : end;
    }
    end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;

    return end;
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool CForwardContract::priceDeadInstrument(CControl* control, CResults* results) const
{
    double value = 0.0;

    DateTime matDate = exerciseSchedule->lastDate();
    if (valueDate < matDate) {
        return false; // not dead yet
    }
    
    DateTime settlementDate = instSettle->settles(matDate, asset.get());
    if (valueDate >= settlementDate) {
        value = 0.0; // all done
    }
    else {
        // not yet settled
        double strike = exerciseSchedule->lastValue();
        value = spotAtMaturity - strike;
        // pv from settlement to today
        value *= discount->pv(valueDate, settlementDate);
        if (!oneContract) {
            value *= notional/initialSpot;
        }    

        if (instSettle->isPhysical()) {
            DateTime dropDate = PhysicalDelivery::exclusionDate(matDate, 
                                                                asset.get());

            if (valueDate > dropDate) {
                value = 0.0;
            }
        }
    }

    // store results
    results->storePrice(value, discount->getCcy());
    if (control && control->isPricing()) {
        recordRequests(control, results);
    }
    return true;
}

CreditSupportSP CForwardContract::createCreditSupport(CMarketDataSP market){
    return CreditSupportSP(new ForwardContractCreditSupport(this, market));
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* CForwardContract::createProduct(
    CClosedFormLN* model) const{
    return new CForwardContractClosedFormProd(this);
}


// handle output requests
void CForwardContract::recordRequests(Control* control, Results* results) const {
    try {
        if (control && control->isPricing()) {
            OutputRequest* request = NULL;

            DateTime matDate = exerciseSchedule->lastDate();
            DateTime settles = instSettle->settles(matDate, asset.get());
            double   strike  = exerciseSchedule->lastValue();

            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArray paydates(1, settles);
                OutputRequestUtil::recordPaymentDates(control,results,&paydates); 
            }

            if (instSettle->isPhysical()) {
                request = control->requestsOutput(OutputRequest::PHYSICAL_DELIVERY);
                if (request && valueDate >= matDate) {
                    double shares = oneContract ? 1.0 : notional/initialSpot;

                    PhysicalDelivery::recordPhysicalDelivery(shares,
                                                             strike,
                                                             matDate,
                                                             asset.get(),
                                                             control, 
                                                             results);
                }
            }
            else {
                request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
                if (request && valueDate >= matDate) {
                    CashFlowArray cfl;
                    double payoff = spotAtMaturity - strike;
                    if (!oneContract) {
                        payoff *= notional/initialSpot;
                    }
                    CashFlow cf(settles, payoff);
                    cfl.push_back(cf);
                    
                    OutputRequestUtil::recordKnownCashflows(control,
                                                            results,
                                                            discount->getCcy(),
                                                            &cfl);   
                }       

            }

            if (valueDate <= matDate) {
                InstrumentUtil::recordFwdAtMat(control,
                                               results,
                                               matDate,
                                               valueDate,
                                               asset.get());
            }                
        }
    }
    catch (exception&) {
        // don't die for this
    }
}


/** Returns the name of the instrument's discount currency. */
string CForwardContract::discountYieldCurveName() const {
    return discount.getName();
}


class ForwardContractNumerical: public NumericalIntegrationLN::IProduct{
private:
    const CForwardContract* fwd; // a reference
    
    DateTime maturity;
    DateTime settles;
    double   strike;

public:
    ForwardContractNumerical(const CForwardContract* fwd): fwd(fwd) {
        maturity = fwd->exerciseSchedule->lastDate();
        strike   = fwd->exerciseSchedule->lastValue();
        settles  = fwd->instSettle->settles(maturity, fwd->asset.get());
    }

    void price(NumericalIntegrationLN* model,
               Control*                control, 
               CResults*               results) {
        static const string method = "ForwardContractNumerical::price";
        try {
            double value = model->integrate(this);

            results->storePrice(value, fwd->discount->getCcy());

            // take care of additional outputs
            if (control && control->isPricing()) {
                // DELAY_PRICE
                InstrumentUtil::delayPriceHelper(control,
                                                 results,
                                                 value,
                                                 fwd->valueDate,
                                                 fwd->discount.get(),
                                                 fwd->asset.get(),
                                                 fwd->premiumSettle.get());
                // FWD_AT_MAT
                InstrumentUtil::recordFwdAtMat(control,
                                               results,
                                               time(),
                                               fwd->valueDate,
                                               fwd->asset.get());
            } 
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    PDFCalculator* pdfCalculator() const {
        static const string method = "ForwardContractNumerical::pdfCalculator";
        try {
            LinearStrikeVolRequest volRequest(0.0,
                                              fwd->fwdStarting?fwd->startDate:
                                                               fwd->valueDate,
                                              maturity,
                                              fwd->fwdStarting);

            PDFRequestLNStrike pdfRequest(&volRequest);
            
            return fwd->asset->pdfCalculator(&pdfRequest);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    double payoff(double spot) const {
        static const string method = "ForwardContractNumerical::payoff";
        try {
            double k = strike;
            if (fwd->valueDate.isGreaterOrEqual(settles)) {
                return 0.0;
            }
            else {
                CAsset::FwdValueAlgorithm algo(fwd->divReinvest);
                double fwdAtStart;
                if (fwd->fwdStarting) {
                    DateTimeArray fwdDates(1, fwd->startDate);
                    DoubleArray   fwds(1);
                    fwd->asset->fwdValue(fwdDates, algo, fwds);
                    k *= fwds[0];
                    spot *= fwds[0];
                    fwdAtStart = fwds[0];
                }

                // taking care of beyond maturity
                if (fwd->valueDate.isGreaterOrEqual(maturity)) {
                    spot = fwd->spotAtMaturity;
                }

                double pv = fwd->discount->pv(settles);

                double value = pv * (spot - k);

                if (!fwd->oneContract) {
                    // handle fixed notional 
                    if (fwd->fwdStarting) {
                        value *= fwd->notional;
                    }
                    else {
                        // handle position 
                        value *= fwd->notional/fwd->initialSpot;
                    }
                }     

                return value;
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    double centre(const DateTime& date) const {
        if (fwd->fwdStarting) {
            return 1.0;
        }
        return fwd->asset->fwdValue(date);
    }

    double variance(double strike, const DateTime& date) const {
        static const string method = "ForwardContractNumerical::variance";
        try {
            LinearStrikeVolRequest volRequest(strike,
                                              fwd->fwdStarting?fwd->startDate:
                                                               fwd->valueDate,
                                              maturity,
                                              fwd->fwdStarting);

            CVolProcessedBSSP volBS(fwd->asset->getProcessedVol(&volRequest));

            return volBS->CalcVar(fwd->valueDate, date);                                   
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }   
    }

    DateTime time() const {
        return maturity;
    }
};

/** Implementation of NumericalIntegrationLN::IntoProduct interface */
NumericalIntegrationLN::IProduct* CForwardContract::createProduct(
    NumericalIntegrationLN* model) const{
    return new ForwardContractNumerical(this);
}

////////////////////////////Vol Index Forward////////////////////////////

/** two models in order to call Var Swap BS and SV
*/
class VForwardProd: virtual public VIXFModel::IProduct {
public:
    
    VForwardProd(const CForwardContract* inst): inst(inst){}
    
    /** This is the method responsible for calling vAsset->fwdValue
        and output additional components (BS, SV) 
    */
    void price(VIXFModel*       model,
               Control*      control, 
               CResults*     results) const;

private:
    const CForwardContract* inst; // a reference

};

/** This is the method responsible for calling vAsset->fwdValue
    and output additional components (BS, SV) 
*/
void VForwardProd::price(VIXFModel*       model,               
						  Control*          control, 
                          CResults*         results) const {
    static const string method = "VForwardProd::price";
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
        totalValue = settlementPV*discFactor*(totalValue - strike);

        /** vector "out" contains BS, SV components, only has values after calling vAsset->fwdValue*/
        /** record "out" into results based on the control */
        vAsset->recordExtraOutput(control, results);

        results->storePrice(totalValue, inst->discount->getCcy());

        if (control && control->isPricing()) {            
            // take care of additional outputs
            inst->recordRequests(control, results);
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
};

//////////////////////////////////////////////////////////////////////////////////

/** Implementation of VIXFModel::IntoProduct interface */
VIXFModel::IProduct* CForwardContract::createProduct(const VIXFModel* model) const {
    return new VForwardProd(this);
}

DRLIB_END_NAMESPACE
