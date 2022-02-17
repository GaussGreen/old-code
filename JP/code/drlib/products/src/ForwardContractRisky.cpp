//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ForwardContractRisky.cpp
//
//   Description   forward contract
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ForwardContractRisky.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/CanBeRisky.hpp"
#include "edginc/RiskyCDSCurve.hpp"
#include "edginc/BadDayConventionFactory.hpp"


DRLIB_BEGIN_NAMESPACE

class CForwardContractRiskyHelper {
public:

    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spread sheet
        REGISTER(CForwardContractRisky, clazz);
        SUPERCLASS(Generic1FactorCredit);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(LastSensDate);   
        EMPTY_SHELL_METHOD(defaultCForwardContractRisky);
        FIELD(startDate,        "Option start date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(fwdStarting, "Is it a fwd starting option");
        FIELD(exerciseSchedule, "Maturity and strike, one element only");
        FIELD(initialSpot,      "Initial spot price");
        FIELD_MAKE_OPTIONAL(initialSpot);
        FIELD(spotAtMaturity,   "underlying spot at maturity date");
        FIELD_MAKE_OPTIONAL(spotAtMaturity);   
        FIELD(zeroDivs, "compute forward price with zeroed divs.");
        FIELD(zeroBorrow, "compute forward price with zeroed borrow.");
    }
        
    static IObject* defaultCForwardContractRisky() {
            return new CForwardContractRisky();
    }
};

CClassConstSP const CForwardContractRisky::TYPE = 
					CClass::registerClassLoadMethod("ForwardContractRisky", 
					typeid(CForwardContractRisky), 
					CForwardContractRiskyHelper::load);
bool   CForwardContractRiskyLoad() {
    return (CForwardContractRisky::TYPE != 0);
}



void CForwardContractRisky::Validate()
{
    static const string method = "ForwardContractRisky::Validate";
    try {
        DateTime            matDate;
        int                 numDates;

        if (zeroBorrow) {
            throw ModelException(method, 
                                 "Zeroing borrow cost is not yet supported");
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

        if (fwdStarting)
        {
            throw ModelException(method,
                                 "risky fwd starting cannot be fwdStarting");

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
void CForwardContractRisky::GetMarket(const IModel* model, const CMarketDataSP market)
{
    static const string method = "CForwardContractRisky::GetMarket";
    try 
    {
        market->GetReferenceDate(valueDate);
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, discount, asset);
        discount.getData(model, market);
        instSettle->getMarket(model, market.get());
        if (premiumSettle.get()) {
            premiumSettle->getMarket(model, market.get());
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

        Generic1FactorCredit::GetMarket(model, market);


    }

    catch (exception& e) {
        throw ModelException(e, method);
    }
    
    #ifdef CDS_BACKWARD_COMPATIBILITY
    cdsParSpreads->setHolidays(HolidaySP(Holiday::weekendsOnly()));
    cdsParSpreads->setBadDayConvention(
        BadDayConventionSP(BadDayConventionFactory::make("None")));
    #endif
}

// constructor
CForwardContractRisky::CForwardContractRisky(): Generic1FactorCredit(TYPE), 
                                      initialSpot(0.0), spotAtMaturity(0.0),
                                      zeroDivs(false) {
    // empty
}

/** private class */
class CForwardContractRiskyClosedFormProd: public CClosedFormLN::IProduct{
private:
    const CForwardContractRisky*  instrFwd; // a reference

public:
    CForwardContractRiskyClosedFormProd(const CForwardContractRisky* instr): instrFwd(instr){}

    void price(CClosedFormLN*  model,
               Control*        control, 
               CResults*       results) const;
};

void CForwardContractRiskyClosedFormProd::price(CClosedFormLN*   model,
                                                Control*        control, 
                                                CResults*       results) const
{
    static const string method = "CForwardContractRiskyClosedForm::price";
    try {
        double         riskyDiscFactor;       
        double         settlementPV;
        double         fwdAtStart = 0.0;      // forward price at start date
        DateTime       settlementDate;  // instrument settlement date
        double         premium;         // the fair value 
        double         fwdAtMat;        // forward at maturity

        // create a risky curve
        IYieldCurveSP tmpCSC = instrFwd->getICDSParSpreads()->makeRiskyCurve(*instrFwd->discount.get());
        YieldCurveSP riskyCurve(dynamic_cast<YieldCurve*>(tmpCSC.get()));

        // there can be only one ...
        DateTime matDate = instrFwd->exerciseSchedule->lastDate();
        double   strike  = instrFwd->exerciseSchedule->lastValue();

        // get settlement date
        settlementDate = instrFwd->instSettle->settles(matDate, 
                                                       instrFwd->asset.get());

        if (instrFwd->valueDate >= settlementDate) {
            premium = 0.0;
        }
        else {
            CAsset::FwdValueAlgorithm algo(instrFwd->zeroDivs);
            if (instrFwd->fwdStarting) {
                throw; // never get here until fwd starting is supported 
            }

            // typical case: before maturity
            if (!(instrFwd->valueDate >= matDate)) {
                DateTimeArray fwdDates(1, matDate);
                DoubleArray   fwds(1);
                                // check to see that we can create a risky asset
                if (ICanBeRisky::TYPE->isInstance(instrFwd->asset.get()))
                {
                    // cast ICanBeRisky pointer
                    CAssetSP riskyAsset = CAssetSP(copy(instrFwd->asset.get()));
                    IObject* obj            = dynamic_cast<IObject*>(riskyAsset.get());
                    smartPtr<ICanBeRisky> riskyAssetI = smartPtr<ICanBeRisky>(
                        dynamic_cast<ICanBeRisky*>(obj));
                    //CreditCurveWrapper* cs = makeRiskyCurve(parSpreadCurve);
                    ICreditCurveSP cdsParSpreads2(copy(instrFwd->cdsParSpreads.get()));
                    const DateTime endDate( matDate.isGreater(settlementDate)? matDate : settlementDate) ;
                    riskyAssetI->makeRisky(cdsParSpreads2,
                                           &endDate);
                    // now we have a risky asset, can call fwdValue to get risky 
                    riskyAsset->fwdValue(fwdDates, algo, fwds);
                    fwdAtMat = fwds[0];
                }
            }
            else { // taking care of beyond maturity
                fwdAtMat = instrFwd->spotAtMaturity;
            }

            // calculate the discount factor back to value date
            riskyDiscFactor = riskyCurve->pv(instrFwd->valueDate, matDate);

            // discounting for settlement: 
            settlementPV = instrFwd->instSettle->pvAdjust(matDate, 
                                                          riskyCurve.get(),
                                                          instrFwd->asset.get());


            // premium for fwd contract
            premium = riskyDiscFactor*settlementPV*(fwdAtMat - strike);

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
        }

        results->storePrice(premium, instrFwd->discount->getCcy());

        // take care of additional outputs
        if (control && control->isPricing()) {
            // DELAY_PRICE
            InstrumentUtil::delayPriceHelper(control,
                                             results,
                                             premium,
                                             instrFwd->valueDate,
                                             instrFwd->discount.get(),
                                             instrFwd->asset.get(),
                                             instrFwd->premiumSettle.get());
            // FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           matDate,
                                           instrFwd->valueDate,
                                           instrFwd->asset.get());
        } 
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}
    

/** Rolls the value date and sets initial spot if rolling over start date */
bool CForwardContractRisky::sensShift(Theta* shift)
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

DateTime CForwardContractRisky::getValueDate() const 
{
    return valueDate;
}

/** when to stop tweaking */
DateTime CForwardContractRisky::endDate(const Sensitivity* sensControl) const {
    DateTime end;

    DateTime maturity = exerciseSchedule->lastDate();
    DateTime instEnd  = instSettle->settles(maturity, asset.get());
    DateTime assetEnd = asset->settleDate(maturity);
    end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;

    return end;
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool CForwardContractRisky::priceDeadInstrument(CControl* control, CResults* results) const
{
    double value = 0.0;

    DateTime matDate = exerciseSchedule->lastDate();
    if (valueDate < matDate)
        return false; // not dead yet
    
    DateTime settlementDate = instSettle->settles(matDate, asset.get());
    if (valueDate < settlementDate)
    {// not yet settled
        double strike        = exerciseSchedule->lastValue();
        value = GetIntrinsic(spotAtMaturity,
                             strike,
                             true, /*isCall*/
                             false /* isOption */);
        
        // pv from settlement to today
        value *= discount->pv(valueDate, settlementDate);
    }

    // store results
    results->storePrice(value, discount->getCcy());

    return true;
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* CForwardContractRisky::createProduct(
    CClosedFormLN* model) const{
    return new CForwardContractRiskyClosedFormProd(this);
}

DRLIB_END_NAMESPACE
