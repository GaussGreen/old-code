//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Vanilla.cpp
//
//   Description : Vanilla instrument
//
//   Author      : Andre X Segger
//
//   Date        : 23 April 2001
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/Black.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/XCB.hpp"
#include "edginc/FD1FDDE.hpp"
#include "edginc/FD1FDDE.hpp"
#include "edginc/RootFinder.hpp"
#include "edginc/FunctionWrapper.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/VegaMatrixLite.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/EndDateCollector.hpp"

#include "edginc/TreeSliceOper.hpp"
#include "edginc/VarSwapUtilities.hpp"
#include "edginc/SimpleEquity.hpp"

DRLIB_BEGIN_NAMESPACE

void CVanilla::Validate()
{
    static const string method = "Vanilla::Validate";
    // just check the things that aren't/cannot be checked in
    // validatePop2Object
    if (!asset){
        throw ModelException(method, "Asset is null");
    }
    if (!discount){
        throw ModelException(method, "Discount YC is null");
    }

    AssetUtil::assetCrossValidate(asset.get(),
                                  fwdStarting,
                                  startDate,
                                  valueDate,
                                  discount,
                                  this);

    if (noExerciseWindow < 0) {
        throw ModelException(method,
                             "noExerciseWindow must not be negative");
    }
}

void CVanilla::validatePop2Object()
{
    static const string method("Vanilla::validatePop2Object");
    int                 numDates;

    // can't get exercise schedule from Market - fail if it is NULL
    if( !exerciseSchedule )
    {
        throw ModelException(method, "Exercise schedule is NULL");
    }

    // can't get instrument settlement from Market - fail if it is NULL
    if( !instSettle )
    {
        throw ModelException(method, "Instrument settlement is NULL");
    }

    // if fwd starting, can't be one contract
    if ( fwdStarting && oneContract)
    {
        throw ModelException(method, "Can't be forward starting and "
                             "one contract");
    }

    // check that we have at least one entry in the exercise schedule
    numDates = exerciseSchedule->length();
    if ( numDates < 1)
    {
        throw ModelException(method, "Exercise schedule is empty");
    }

    // not sure whether we should override any of these value ...
    if (oneContract)
    {
        /* one contract is pricing for a fixed notional where
           notional/initial spot = 1.0 */
        notional    = 1.0;
        initialSpot = 1.0;
    }

    if (!canExerciseEarly && noExerciseWindow != 0){
        throw ModelException(method, "When the option is not American, noExerciseWindow should be == 0.");
    }
}

void CVanilla::GetMarket(const IModel*          model,
                         const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);

    CAsset::getAssetMarketData(model, market.get(), ccyTreatment,
                               discount, asset);

    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
    if (premiumSettle.get())
    {
        premiumSettle->getMarket(model, market.get());
    }
}

DateTime CVanilla::getValueDate() const
{
  return valueDate;
}

/** when to stop tweaking */
DateTime CVanilla::endDate(const Sensitivity* sensControl) const {
    DateTime matDate = exerciseSchedule->lastDate();
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;

    // Check any components of the asset implementing the LastSensDate interface
    DateTime maxDate;
    EndDateCollector collector(maxDate, sensControl);
    maxDate = collector.getMaxEndDate(IObjectSP(asset.getMO()));

    end = end.isGreater(maxDate) ? end : maxDate;

    return end;
}

void CVanilla::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime matDate = exerciseSchedule->lastDate();
        // DELAY_PRICE
        InstrumentUtil::delayPriceHelper(control,
                                         results,
                                         fairValue,
                                         valueDate,
                                         discount.get(),
                                         asset.get(),
                                         premiumSettle.get());
        // IND_VOL
        OutputRequest* request = NULL;
        if (control->requestsOutput(OutputRequest::IND_VOL, request))
        {
             if (isCall)
             {
                 double indVolThreshold = 0.01;
                 double strike = exerciseSchedule->lastValue();
                 if (!fwdStarting)
                 {
                     double spot = asset->getSpot();
                     indVolThreshold *= spot;
                 }
                 // If strike is less than 1% of the spot store NotApplicable.
                 // This is required by the result combining code later
                 if (strike < (indVolThreshold))
                 {
                     results->storeNotApplicable(request);
                 }
                 else
                 {
                     results->storeRequestResult(request, indVol);
                 }
             }
             else
             {
                 results->storeRequestResult(request, indVol);
             }
        }

        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       matDate,
                                       valueDate,
                                       asset.get());
    }
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool CVanilla::priceDeadInstrument(CControl* control, CResults* results) const
{
    double    strike        = 0.0;
    double    value         = 0.0;
    bool      foundExerDate = false;
    DateTime  exerDate;

    static string method = "Vanilla::priceDeadInstrument";

    DateTime matDate = exerciseSchedule->lastDate();

    bool expired = (valueDate >= matDate || (isExercised && canExerciseEarly));
    if (!expired)
        return false; // not dead yet

    DateTime settlementDate = instSettle->settles(matDate, asset.get());

    if (valueDate >= settlementDate)
    {// settled already
        results->storePrice(0.0, discount->getCcy());
        addOutputRequests(control, results, 0.0, 0.0);
        return true;
    }

    // we may simply consider isExercised flag
    // but here we give the instrinsic value at maturity regardless
    if (!isExercised)
    {
        // maturity instrinsic value
        exerDate      = matDate;
        strike        = exerciseSchedule->lastValue();
        foundExerDate = true;
    }
    else
    {
        // may early exercise
        // check that exercise date is not in the Future
        // soft check on date as in Pyramid EOD still runs as of SOD
        if (dateExercised.getDate() > valueDate.getDate()) {
            throw ModelException(method,
                    "Option exercised on " +
                    dateExercised.toString() + ". " +
                    "This date is after the current value date " +
                    valueDate.toString());
        }

        // check whether exercise date is valid and find corresponding strike
        if ( !canExerciseEarly ) {
            // multi-european case
            DateTimeArray exerciseDates = exerciseSchedule->getDates();

            for (int i=0 ; i<exerciseDates.size() ; ++i)
            {
                if ( dateExercised.equals(exerciseDates[i])) {
                    foundExerDate = true;
                    exerDate      = dateExercised;
                    strike        = exerciseSchedule->interpolate(
                                                    dateExercised);
                    break;
                }
                else
                {
                    throw ModelException(method,
                            "Option has been exercised on " +
                            dateExercised.toString() + ". " +
                            "This date could not be found in the exercise schedule.");
                }
            }
        } else {
            if (  dateExercised.isGreaterOrEqual(
                        exerciseSchedule->firstDate()) &&
                 !dateExercised.isGreater(
                        exerciseSchedule->lastDate())  ) {
                // american case - note: not checking for weekend yet
                foundExerDate = true;
                exerDate      = dateExercised;
                strike        = exerciseSchedule->interpolate(
                                                dateExercised);
            }
            else
            {
                throw ModelException(method,
                        "Option has been exercised on " +
                        dateExercised.toString() + " but valid " +
                        "exercise range is from " +
                        exerciseSchedule->firstDate().toString() +
                        " to " +
                        exerciseSchedule->lastDate().toString());
            }
        }
    }

    if ( foundExerDate )
    {
        if (!Maths::isPositive(spotAtMaturity)){
            throw ModelException(method,
                                 "Option has exercised on "+exerDate.toString() +
                                 ", but spot at maturity for asset "
                                 + asset->getTrueName() + " is missing");
        }

        value = GetIntrinsic(spotAtMaturity,
                                 strike,
                                 isCall,
                                 true /* isOption */);

        settlementDate = instSettle->settles(exerDate, asset.get());
        // pv from settlement to today
        value *= discount->pv(valueDate, settlementDate);
    }

    double scalingFactor = InstrumentUtil::scalePremium(
        oneContract,
        false,
        notional,
        0.0,
        initialSpot);

    value *= scalingFactor;

    // store results
    results->storePrice(value, discount->getCcy());
    addOutputRequests(control, results, value, 0.0);

    return true;
}

/** for ITaxableInst::Basic */
const DateTime CVanilla::getFinalPaymentDate() const {
    // can't support anything which may have early cash flows
    if (canExerciseEarly || isExercised) {
        throw ModelException("Vanilla::getFinalPaymentDate",
                             "Tax is not supported if early exercise is possible");
    }
    DateTime matDate = exerciseSchedule->lastDate();
    return instSettle->settles(matDate, asset.get());
}

/** private class */
class CVanillaClosedForm: public CClosedFormLN::IProduct{
private:
    const CVanilla*  vanilla; // a reference

public:
    CVanillaClosedForm(const CVanilla* vanilla): vanilla(vanilla){}

    void price(CClosedFormLN*   model,
               Control*        control,
               CResults*       results) const;
};

void CVanillaClosedForm::price(CClosedFormLN* model,
                               Control*       control,
                               CResults*      results) const
{
    static const string method = "VanillaClosedForm::price";
    try {
        double         discFactor;
        double         variance;        // variance between start date and mat.
        double         fwdAtStart;      // forward price at start date
        double         premium;         // the fair value
        double         fwdPrice;        // forward at maturity
        double         ivol;            // indicative vol

        DateTime        matDate   = vanilla->exerciseSchedule->lastDate();
        double          strike    = vanilla->exerciseSchedule->lastValue();
        double          absStrike = strike;

        // valueDate >= matDate is taken care of here
        if(vanilla->priceDeadInstrument(control, results)){
            return; // dead instrument priced
        }

        /* if it's forward starting, convert payoffToProcFreqBoundWeight strikes
           to absolute values based on spot at start date */
        DateTime imntStartDate = vanilla->fwdStarting?
            vanilla->startDate: vanilla->valueDate;

        if (vanilla->fwdStarting)
        {
            fwdAtStart = vanilla->asset->fwdValue(vanilla->startDate);
            absStrike *= fwdAtStart;
        }

        // price live instrument
        fwdPrice = vanilla->asset->fwdValue(matDate);

        // choose how to interpolate the vol - go for traditional route for now
        LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(
                                                   strike,
                                                   imntStartDate,
                                                   matDate,
                                                   vanilla->fwdStarting));

        // interrogate the model to see if neg fwd vars are allowed
        volRequest->allowNegativeFwdVar(model->negativeFwdVarAllowed());

        // interpolate the vol using our LN request
        CVolProcessedBSSP volBS(vanilla->asset->
                                getProcessedVol(volRequest.get()));
        // calculate the variance
        variance = volBS->CalcVar(imntStartDate, matDate);

        // calculate the indicative vol
        try {
            ivol = volBS->CalcVol(imntStartDate, matDate);
        }
        catch (exception& ) {
            ivol = 0.0;
        }

        // call Black model without discounting
        premium = Black::price(vanilla->isCall, fwdPrice,
                               absStrike, 1.0, variance);
        // and discount
        discFactor = vanilla->instSettle->pv(vanilla->valueDate,
                                             matDate,
                                             vanilla->discount.get(),
                                             vanilla->asset.get());

        premium *= discFactor;

        double scalingFactor = InstrumentUtil::scalePremium(
                                        vanilla->oneContract,
                                        vanilla->fwdStarting,
                                        vanilla->notional,
                                        fwdAtStart,
                                        vanilla->initialSpot);

        premium *= scalingFactor;

        // VEGA_MATRIX_LITE
        if (control)
        {
            SensitivitySP sens(control->getCurrentSensitivity());
            VegaMatrixLiteSP vml(dynamic_cast<VegaMatrixLite *>(sens.get()));
            if (vml.get()) {
                // Record options
                VanillaContractsRecorderSP recorder = VanillaContractsRecorder::createVanillaOptionRecorder(control);
                double yearFrac = volBS->calcTradingTime(imntStartDate, matDate);
                recorder->recordContract(matDate, discFactor * scalingFactor, vanilla->isCall, 
                    fwdPrice, absStrike, ivol, yearFrac);

                // Save results
                VanillaInfo::storeVegaMatrix(
                    vml,
                    const_cast<CVanilla *>(vanilla),
                    model,
                    recorder,
                    CAssetSP::constCast(vanilla->asset.getSP()),
                    vanilla->valueDate,
                    matDate,
                    results);
                
                return;
            }
        }

        results->storePrice(premium, vanilla->discount->getCcy());

        vanilla->addOutputRequests(control,
                                   results,
                                   premium,
                                   ivol);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/////////////////////////////////////////////////////////
//     private class for all tree/FD product
//     new state variable interface for any num of factors
/////////////////////////////////////////////////////////
class VanillaFDProd : public LatticeProdEDR, virtual public IFDProductLN
{
private:
    const CVanilla*     instrVan;
    vector<double>      stepStrike;
    vector<bool>        stepCanExercise;

    // to allow switching between original and using slice operators update
    typedef void ( VanillaFDProd::*prod_FUNC )(
        int step, const TreeSlice & spot, const vector< TreeSliceSP > & price );
    prod_FUNC prod_BWD_T, prod_BWD;

public:
    VanillaFDProd(const CVanilla* van, FDModel* m) :
        LatticeProdEDR(m), instrVan(van)
    {
        if( tree1f )
        {
            if (tree1f->DivsTreatedAsAbsolute() && instrVan->fwdStarting){
                tree1f->SetDivAmountTreatment(false);
            }
        }

        // first: set discount curve
        if( tree1f )
            tree1f->setDiscountCurve( instrVan->discount.getSP() );

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ(
                instrVan->asset.getName(),
                instrVan->asset,
                instrVan->ccyTreatment ) ) );
    }

    void AdjustDeltaShift(CControl* control) const;

    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const {return instrVan->ccyTreatment;}

    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(Control*  control) const;

    /** initialising and setting product variables */
    // this is called per pricing call before each pricing
    virtual void initProd();

    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type);

    /** return price after postPrice process */
    double postPrice(double euroPrice, double exer_premium, YieldCurveConstSP disc);
    // scale for notional
    double scalePremium(const double& fairValue);

    /** output prices and any additional results */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results);

    virtual DateTime getStartDate() const {
        return instrVan->fwdStarting ? instrVan->startDate : instrVan->valueDate;
    }

    // for the LogNormal model
    CVolRequestLNSP getVolInterp(int iAsset) const {

    //CVolRequestLNArray getVolInterp(int iAsset) const {

        // get strike and maturity date from instrument
        DateTime        matDate = instrVan->exerciseSchedule->lastDate();

        double volStrike  = instrVan->exerciseSchedule->lastValue();

        DateTime imntStartDate = instrVan->fwdStarting?
                             instrVan->startDate: instrVan->valueDate;

        CVolRequestLNSP   reqarr;

        reqarr = CVolRequestLNSP(
                            new LinearStrikeVolRequest(volStrike, imntStartDate,
                            matDate, instrVan->fwdStarting));
        return reqarr;
    }

private:
    //local methods

    class PenultSmooth : public SliceMarker< PenultSmooth >
    {
    public:
        PenultSmooth(
            const TreeSlice & spot,
            const TreeSlice & price,
            bool isCall,
            vector< double > & vol_arr,
            vector< double > & drift_arr,
            double dt,
            bool fastRoll,
            double df,
            double stepStrike,
            double truncationStd,
            double stockFloor0,
            double stockFloor1 )
            :
            spot( spot ),
            spotValues( spot.getValues() ),
            price( price ),
            isCall( isCall ),
            vol_arr( vol_arr ),
            drift_arr( drift_arr ),
            dt( dt ),
            fastRoll( fastRoll ),
            df( df ),
            stepStrike( stepStrike ),
            truncationStd( truncationStd ),
            stockFloor0( stockFloor0 ),
            stockFloor1( stockFloor1 )
        {
            variance = vol_arr[0]*vol_arr[0]*dt;
            drift = drift_arr[0];
        }

        // TreeSlice "expression template" primitives
        static const int sliceCount = 2;
        template< typename S >
        const S** listSlices(const S** list) const
        {
            return price.listSlices( spot.listSlices( list ) );
        }
        inline double calc() const
        {
            return apply( spot.iter - spotValues, spot.calc(), price.calc() );
        }
        void printDebug(char *s) const
        {
            strcat(s, "(PenultSmooth)");
        }

    private:
        double apply( int i, double spot, double price ) const
        {
            if (!fastRoll)
            {
                if (vol_arr.size() > 1) // Bug fix RG
                    variance = vol_arr[i]*vol_arr[i]*dt; // needs one vol per node - local vol

                if (drift_arr.size() > 1) // Bug fix RG
                    drift = drift_arr[i]; // needs drift per node
            }

            if (spot>0.0 && fabs(log(spot/stepStrike)) < truncationStd*sqrt(variance))
            {
                /* If dollar div treatment (in which case stock floor is <> 0.0), the price of the option
                is that of an option with same maturity written on the pseudo asset (i.e., S - StockFloor)
                and with strike adjusted by the stock floor at maturity (i.e., K - StockFloor).
                NB Black::price takes care of potentially negative fwd and strike values. */
                return Black::price( isCall,
                    drift * (spot - stockFloor0),  // pseudo asset's (step+1)-maturity fwd
                    // as viewed from t = step
                    stepStrike - stockFloor1,  // adjusted strike
                    df,
                    variance );
            }
            else
                return price;
        }

        const TreeSlice & spot;
        const double * spotValues;
        const TreeSlice & price;
        const bool isCall;
        const vector< double > & vol_arr;
        const vector< double > & drift_arr;
        const double dt;
        const bool fastRoll;
        const double df;
        const double stepStrike;
        const double truncationStd;
        const double stockFloor0;
        const double stockFloor1;

        mutable double variance;
        mutable double drift;
    };

    /** product payoff method at maturity */
    void prod_BWD_T_orig( int step, const TreeSlice & spot, const vector< TreeSliceSP > & price )
    {
        int bot, top;
        spot.getCalcRange( bot, top );
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        double settlementPV = instrVan->instSettle->pvAdjust(instrVan->exerciseSchedule->lastDate(),
                                                             instrVan->discount.get(),
                                                             instrVan->asset.get());

        int nPrice = price.size();
        for( int i = 0; i < nPrice; ++i )
        {
            for( int j = bot; j <= top; ++j )
            {
                p[i][j] = settlementPV *
                    GetIntrinsic( s[j], stepStrike[step], instrVan->isCall, true );
            }
        }
    }

    /** product payoff method at maturity */
    void prod_BWD_T_oper( int step, const TreeSlice & spot, const vector< TreeSliceSP > & price )
    {
        double settlementPV = instrVan->instSettle->pvAdjust(
            instrVan->exerciseSchedule->lastDate(),
            instrVan->discount.get(),
            instrVan->asset.get() );

        int nPrice = price.size();
        for( int i = 0; i < nPrice; ++i )
        {
            // Using template slice operators
            *price[i] = smax( 0., ( spot - stepStrike[step] ) * ( instrVan->isCall ? settlementPV : -settlementPV ) );
        }
    }

    /** product payoff method at steps earlier than maturity */
    void prod_BWD_orig( int step, const TreeSlice & spot, const vector< TreeSliceSP > & price )
    {
        static const string method = "VanillaFDProd::prod_BWD";
        try {
            int bot, top;
            spot.getCalcRange( bot, top );
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            if (tree1f && step == model->getLastStep()-1 && payoffIndex->isElementary())
            {
                vector <double> vol_arr, drift_arr;
                double dt = tree1f->getTrdYrFrac(step+1);

                bool fastRoll = tree1f->CalcStepDriftAndVar(s, bot, top, vol_arr, &drift_arr); // just to get drift
                // use GetStepVol() to get vol, do not use GetStepVar() which may not have simple BS conversion
                tree1f->GetStepVol(step, vol_arr, s, bot, top);

                double variance = vol_arr[0]*vol_arr[0]*dt;
                double drift = drift_arr[0];
                double df = instrVan->discount->pv(model->getDate(step),
                                                   model->getDate(step+1));

                /* If dollar dividend treatment, we must adjust fwd and strike by stock floor.
                   If not dollar dividend treatment, stock floor == 0.0 */
                double stockFloor0 = tree1f->GetStockFloor(step);
                double stockFloor1 = tree1f->GetStockFloor(step + 1);

                for( int j = bot; j <= top; ++j )
                {
                    if (!fastRoll)
                    {
                        if (vol_arr.size() > 1){    // Bug fix RG
                            variance = vol_arr[j-bot]*vol_arr[j-bot]*dt; // needs one vol per node - local vol
                        }
                        if (drift_arr.size() > 1){    // Bug fix RG
                            drift = drift_arr[j-bot]; // needs drift per node
                        }
                    }

                    if (s[j]>0.0 && fabs(log(s[j]/stepStrike[step+1])) < tree1f->TruncationStd*sqrt(variance))
                    {
                        /* If dollar div treatment (in which case stock floor is <> 0.0), the price of the option
                            is that of an option with same maturity written on the pseudo asset (i.e., S - StockFloor)
                            and with strike adjusted by the stock floor at maturity (i.e., K - StockFloor).
                            NB Black::price takes care of potentially negative fwd and strike values. */
                        p[0][j] = p[1][j] =
                            Black::price(
                                instrVan->isCall,
                                drift * (s[j] - stockFloor0),  // pseudo asset's (step+1)-maturity fwd
                                // as viewed from t = step
                                stepStrike[step+1] - stockFloor1,  // adjusted strike
                                df,
                                variance);
                    }
                }
            }

            if(!stepCanExercise[step])
                return;

            double settlementPV = instrVan->instSettle->pvAdjust(model->getDate(step),
                                                                 instrVan->discount.get(),
                                                                 instrVan->asset.get());

            double callput = instrVan->isCall ? settlementPV : -settlementPV;

            for( int j = bot; j <= top; ++j )
            {
                double intrinsic = callput * ( s[j] - stepStrike[step] );
                if( p[0][j] < intrinsic )
                    p[0][j] = intrinsic; // American
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** product payoff method at steps earlier than maturity */
    void prod_BWD_oper( int step, const TreeSlice & spot, const vector< TreeSliceSP > & price )
    {
        static const string method = "VanillaFDProd::prod_BWD";
        try {
            if (tree1f && step == model->getLastStep()-1 && payoffIndex->isElementary())
            {
                vector <double> vol_arr, drift_arr;
                double dt = tree1f->getTrdYrFrac(step+1);

                int bot, top;
                spot.getCalcRange( bot, top );
                double * s = spot.getValues();

                bool fastRoll = tree1f->CalcStepDriftAndVar(s, bot, top, vol_arr, &drift_arr); // just to get drift
                // use GetStepVol() to get vol, do not use GetStepVar() which may not have simple BS conversion
                tree1f->GetStepVol(step, vol_arr, s, bot, top);

                double df = instrVan->discount->pv(model->getDate(step),
                                                   model->getDate(step+1));

                /* If dollar dividend treatment, we must adjust fwd and strike by stock floor.
                   If not dollar dividend treatment, stock floor == 0.0 */
                double stockFloor0 = tree1f->GetStockFloor(step);
                double stockFloor1 = tree1f->GetStockFloor(step + 1);

                *price[1] = PenultSmooth(
                    spot,
                    *price[0],
                    instrVan->isCall,
                    vol_arr,
                    drift_arr,
                    dt,
                    fastRoll,
                    df,
                    stepStrike[step+1],
                    tree1f->TruncationStd,
                    stockFloor0,
                    stockFloor1 );

                *price[0] = *price[1];
            }

            if(!stepCanExercise[step])
                return;

            double settlementPV = instrVan->instSettle->pvAdjust(
                model->getDate(step),
                instrVan->discount.get(),
                instrVan->asset.get() );

            // Using template slice operators
            *price[0] = smax( *price[0], ( spot - stepStrike[step] ) * ( instrVan->isCall ? settlementPV : -settlementPV ) );
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
};

// adjust delta
void VanillaFDProd::AdjustDeltaShift(CControl* control) const {
    // only do this for tree1f for now
    if (tree1f) {

        double deltaShift = control->getDeltaShiftSize();

        if (Maths::isPositive(deltaShift)){
            DividendListSP dollarDivs =
                AssetUtil::getDollarDivsBetweenDates(instrVan->asset.get(),
                                                     instrVan->valueDate,
                                                     instrVan->exerciseSchedule->lastDate());
            bool isFwdStart = instrVan->fwdStarting && instrVan->startDate>instrVan->valueDate;

            tree1f->DeltaSizeAdjust(deltaShift,
                                    instrVan->valueDate,
                                    instrVan->exerciseSchedule->firstDate(),
                                    instrVan->exerciseSchedule->lastDate(),
                                    instrVan->exerciseSchedule->length()>1 && !instrVan->canExerciseEarly,
                                    dollarDivs->getDivAmounts()->size()>0,
                                    isFwdStart);
        }
    }
}

//need to replace this function by getVolInterp, to review
/** returns a vol request for log-normal vol */
//CVolRequestConstSP VanillaFDProd::GetLNRequest() const
//{
//    // get strike and maturity date from instrument
//    DateTime matDate = instrVan->exerciseSchedule->lastDate();
//    double volStrike = instrVan->exerciseSchedule->lastValue();
//
//    CVolRequestConstSP volRequest(
//        new LinearStrikeVolRequest(volStrike, getStartDate(),
//                                   matDate, instrVan->fwdStarting));
//    return volRequest;
//}

/** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
    isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
void VanillaFDProd::update(int& step, FDProduct::UpdateType type)
{
    // we assume just need one und level for spot here
    const TreeSlice & spot = payoffIndex->getValue( step );

    if (type == FDProduct::BWD_T)
        ( this->*prod_BWD_T )( step, spot, slices );
    else if(type == FDProduct::BWD)
        ( this->*prod_BWD )( step, spot, slices );
    else{
        // to do fwd induction
    }
};

/** initialise tree1f - allow product customisation
    must not init product variables here, use initProd() instead */
void VanillaFDProd::init(CControl* control) const{
    static const string method = "VanillaFDProd::init()";
    try {
        /** customize tree parameters here and set up the tree */
        DateTimeArray segDates;
        segDates.resize(2);
        // this needs change if fwd start tree has to start today !!!
        if (instrVan->fwdStarting && instrVan->startDate>instrVan->valueDate) {
            segDates[0] = instrVan->startDate;
            if (tree1f){
                tree1f->controlSameGridFwdStart(instrVan->ccyTreatment);
            }
        }
        else {
            segDates[0] = instrVan->valueDate;
        }

        segDates[1] = instrVan->exerciseSchedule->lastDate();
        IntArray density(1,1);
        int numOfPriceArray = 2;

        // all exercise dates are copied to critical dates
        DateTimeArray critDates = instrVan->exerciseSchedule->getDates();

        // add div event dates if needed
        EventAssetMove divEvent;
        DateTimeArraySP divCritDates;
        if (tree1f) {
            int numOfInsertNode = (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ? 1 : 0); // use one inserted node for exercise boundary

            if (instrVan->canExerciseEarly || tree1f->DivsTreatedAsAbsolute()) {
                // American exercise or dollar div interp treatment
                DateTime start(segDates[0]);
                DateTime end(segDates[1]);
                if (tree1f->DivsTreatedAsAbsolute()) {
                    start = min(start, instrVan->valueDate);    // fwd starting not supported anyway when dol divs
                    end = max(end, Equity::calcDivTransPeriodEndDate(instrVan->valueDate));
                }
                int numDivs = (int)(4*start.yearFrac(end))+1; // 4 divs per year selected as critical dates
                if (AssetUtil::getJumpEvents(instrVan->asset.get(),
                                             start,
                                             end,
                                             numDivs,
                                             divEvent)) {

                    // calculate critical dates
                    divCritDates = divEvent.getCritDate(instrVan->noExerciseWindow, instrVan->isCall);
                    /* If dol divs, filter past critical dividend dates out as they are used as bench
                       dates by pseudo asset. Have to do this here (as opposed to in Equity), since
                       if there are no critical div dates left after filtering I revert to
                       DivAmountTreatment == false immediately below */
                    if (tree1f->DivsTreatedAsAbsolute() && divCritDates->size() > 0){
                        DateTimeArray temp(segDates[0].getFutureDates(*divCritDates));
                        divCritDates = DateTimeArraySP(copy(&temp));
                    }
                }
            }

            /* If no dividend critical dates, no point to bother with dol divs */
            if (tree1f->DivsTreatedAsAbsolute() && (!divCritDates || divCritDates->size() == 0)){
                tree1f->SetDivAmountTreatment(false);
            }
            /* If dollar div,
               - kill off control var if call
               - switch on control var if put. */
            if (tree1f->DivsTreatedAsAbsolute()){
                if (instrVan->isCall){
                    tree1f->DEBUG_UseCtrlVar = false;
                }
                else{
                    tree1f->DEBUG_UseCtrlVar = true;
                }
            }

            // use simple delta size adjustment if needed, note that TreeDeltaShift is stored in base tree only.
            if (!tree1f->DEBUG_SameGridDelta && control->isPricing()) {
                AdjustDeltaShift(control);
            }

            // temp solution
            tree1f->NumOfPrice = numOfPriceArray;
            tree1f->NumOfInsertNode = numOfInsertNode;

            tree1f->isCall = instrVan->isCall;
            tree1f->noExerciseWindow = instrVan->noExerciseWindow;

            if( divCritDates.get() )
                tree1f->addDivCritDates( *divCritDates );
        }
        else {
            // finite difference
            if (instrVan->canExerciseEarly) {
                // American exercise or dollar div interp treatment
                // create div events, this call is very expensive for basket with lots of div dates
                // should only need once for pricing call but need to think how to store/copy for tweaks
                int numDivs = (int)(4*segDates[0].yearFrac(segDates[1]))+1; // 4 divs per year selected as critical dates
                if (AssetUtil::getJumpEvents(instrVan->asset.get(),
                                             segDates[0],
                                             segDates[1],
                                             numDivs,
                                             divEvent)){

                    // calculate critical dates
                    divCritDates = divEvent.getCritDate(instrVan->noExerciseWindow, instrVan->isCall);
                }
            }
            if( divCritDates.get() )
                model->addCritDates( *divCritDates );
        }

        // add critical dates
        model->addCritDates( critDates );

        // prepare timeline set up
        model->initSegments( segDates, density );
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** initialising and setting product variables
    this is called per pricing call before tree sweep call (after InitTree)  */
void VanillaFDProd::initProd() {
    static const string method = "VanillaFDProd::initProd()";
    try {
        int i;
        int lastStep = model->getLastStep();
        double fwdAtStart = 1.0;

        initSlices( 2, instrVan->discount->getName() );

        // if tree1f engine then use original way (for performance) otherwise slice operators
        if( tree1f )
        {
            prod_BWD_T = &VanillaFDProd::prod_BWD_T_orig;
            prod_BWD = &VanillaFDProd::prod_BWD_orig;
        }
        else
        {
            prod_BWD_T = &VanillaFDProd::prod_BWD_T_oper;
            prod_BWD = &VanillaFDProd::prod_BWD_oper;
        }

        stepStrike.resize(lastStep+1);
        stepCanExercise.resize(lastStep+1);
        // decide first about steps that can exercise
        AssetUtil::setStepExercise(stepCanExercise,
                                 model->getDates(),
                                 instrVan->exerciseSchedule,
                                 instrVan->canExerciseEarly,
                                 instrVan->asset.getSP());

        /* If no exercise before ex-div date applies */
        if (instrVan->canExerciseEarly && instrVan->noExerciseWindow > 0) {
            /* Get dividend list */
            DateTime lastDate = AssetUtil::getHoliday(instrVan->asset.get())->addBusinessDays(
                model->getDate(lastStep),
                instrVan->noExerciseWindow);

            DividendListSP divList(AssetUtil::getAllDivsBetweenDates(instrVan->asset,
                                                                     model->getDate(0),
                                                                     lastDate));

            /* If dividend list not empty */
            if (divList->getArray().size() > 0){
                HolidayConstSP holiday(
                    AssetUtil::getHoliday(instrVan->asset.get()));
                DividendConstSP div;    // null
                DateTime nextDivDate(0, 0);   // past date
                int iStep = 0;
                for(; iStep <= lastStep; ++iStep){
                    const DateTime& currStepDate = model->getDate(iStep);
                    /* If need to get next ex div date */
                    if (currStepDate.isGreaterOrEqual(nextDivDate)){
                        /* attempt to get next div */
                        div = divList->getNextDivFromDate(currStepDate);
                        /* if there is no next div, break */
                        if (!div){
                            break;
                        }
                        /* otherwise, get ex div date */
                        nextDivDate = div->getExDate();
                    }
                    int daysToNextDiv;
                    /* compute nb of bus days from current time step to next ex date */
                    daysToNextDiv = holiday->businessDaysDiff(currStepDate,
                                                              nextDivDate);
                    /* if current step is exercise but falls within noExerciseWindow days,
                       override exercise to no exercise */
                    if (stepCanExercise[iStep] == true && daysToNextDiv <= instrVan->noExerciseWindow) {
                        stepCanExercise[iStep] = false;
                    }
                }
            }
        }

        bool canInterpExSched = (instrVan->exerciseSchedule->length() > 1);
        // get spot at start if needed
        if (instrVan->fwdStarting) {
            fwdAtStart = instrVan->asset->fwdValue(instrVan->startDate);
        }

        CIntArray  canEx;
        canEx.resize(lastStep+1); // for debug
        canEx[lastStep] = 1;

        // get strikes
        stepStrike[lastStep] = fwdAtStart*instrVan->exerciseSchedule->lastValue();

        for (i=0; i<lastStep; i++) {
            if (stepCanExercise[i]) {
                canEx[i] = 1;
                if (canInterpExSched) {
                    stepStrike[i] = fwdAtStart*instrVan->exerciseSchedule->interpolate(model->getDate(i));
                }
                else {
                    stepStrike[i] = stepStrike[lastStep]; // strike level is constant for american with no schedule
                }
            }
            else {
                canEx[i] = 0;
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// postPrice process
double VanillaFDProd::postPrice(double euroPrice, double exer_premium, YieldCurveConstSP disc)
{
    double price;
    double fwdStartDF = 1.0;

    if (instrVan->fwdStarting && model->getDate(0).getDate() > instrVan->valueDate.getDate())
        fwdStartDF = disc->pv(instrVan->valueDate, model->getDate(0));

    if ((tree1f && tree1f->DEBUG_UseCtrlVar) || (fd1dRet && fd1dRet->DEBUG_UseCtrlVar)) {
        CVanillaClosedForm closedForm(instrVan);
        CClosedFormLN model;
        Control ctrl;
        CResults result;
        closedForm.price(&model, &ctrl, &result);
        double closedFormPrice = result.retrievePrice();

        price = scalePremium(exer_premium*fwdStartDF) + closedFormPrice;
    }
    else {
        price = fwdStartDF*scalePremium(euroPrice+exer_premium);
    }
    return price;
}

/** premium scaling */
double VanillaFDProd::scalePremium(const double& fairValue)
{
    double fwdAtStart = 0.0;
    if (instrVan->fwdStarting)
    {
        fwdAtStart = instrVan->asset->fwdValue(instrVan->startDate);
    }

    double scalingFactor = InstrumentUtil::scalePremium(
                                    instrVan->oneContract,
                                    instrVan->fwdStarting,
                                    instrVan->notional,
                                    fwdAtStart,
                                    instrVan->initialSpot);

    return fairValue * scalingFactor;
}

//output results
void VanillaFDProd::recordOutput(Control* control, YieldCurveConstSP disc, Results* results)
{
    // get prices at t=0
    vector<double> price0(2);
    price0[0] = model->getPrice0( *slices[0] );
    price0[1] = model->getPrice0( *slices[1] );
    // save price
    double price = postPrice(price0[1], price0[0] - price0[1], disc);
    results->storePrice(price, disc->getCcy());

    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime       matDate = instrVan->exerciseSchedule->lastDate();
        double         indVol;
        // calculate indicative vol
        if ( matDate.isGreater(instrVan->valueDate) )
        {
            DateTime imntStartDate = instrVan->fwdStarting?
                             instrVan->startDate: instrVan->valueDate;

            // get vol request
            double volStrike  = instrVan->exerciseSchedule->lastValue();
            LinearStrikeVolRequest volRequest(volStrike, imntStartDate,
                                           matDate, instrVan->fwdStarting);

            try{
                // interpolate the vol
                CVolProcessedSP  vol(instrVan->asset->getProcessedVol(&volRequest));
                // cast to the BS vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);

                // this should never happen if our get market data has worked properly
                if (!volBS){
                    throw ModelException("VanillaFDProd::recordOutput", 
                                         "No Black Scholes Vol");
                }
                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
            }
            catch (exception& ) {
                indVol = 0.0;
            }
        }
        else
        {
            indVol = 0.0;
        }
        instrVan->addOutputRequests(control,
                                        results,
                                        price,
                                        indVol);
    }
}

/** create a fd payoff product */
FDProductSP Vanilla::createProduct(FDModel* model) const
{
    return FDProductSP(new VanillaFDProd(this, model));
}

/// ------------  end of new FD product -----------


/////////////////////////////////////////////////////////
//           private class for tree1f/FD product
/////////////////////////////////////////////////////////
/** vanilla product payoff for a tree/FD */
// this old interface is left here just to keep DDE working
class Vanilla1fProd: virtual public FD1F::IProduct,
virtual public FD1FGeneric::IProduct, virtual public IDDEInitiator
{
public:
    friend class CVanilla;

    Vanilla1fProd(const CVanilla* vanilla):instrVan(vanilla),
        useInsertNode(false), isStaticSpread(true) {}

    virtual CAssetConstSP GetAssetRef();

    virtual YieldCurveConstSP GetDiscCurveRef()
    {
        return YieldCurveConstSP::dynamicCast(
            (IObjectConstSP)instrVan->discount.getSP());
    }

    virtual bool Positive() {
        // option payoff
        return true;
    }

    /** called before PayoffAtMat, PayoffBeforeMat and tree roll() */
    virtual void preCalc(int step, int idx) {}

    virtual CVolRequestConstSP GetLNRequest();

    /** initialise model1f - allow product customisation */
    virtual void Init(CControl* control);

    /** initialising and setting product variables */
    // this is called per pricing call before tree sweep call (after InitTree)
    virtual void InitProd();

    /** product payoff method at maturity */
    virtual void PayoffAtMat(const double * s, int step, int bot,
                                 int top, int pStart, int pEnd,
                                 double * const * price);

    /** product payoff method at steps earlier than maturity */
    virtual void PayoffBeforeMat(const double * s, int step, int bot,
                                     int top, int pStart, int pEnd,
                                     double * const * price);

    virtual string getCcyTreatment(){return instrVan->ccyTreatment;}

    /** premium scaling */
    virtual double scalePremium(const double& fairValue);

    /** extra output requests */
    virtual void recordOutputRequests(Control* control, Results* results,
                                      double fairValue);

    /** control variate */
    virtual double RefinePrice(double basePrice, double fwdStartDF,
                               bool useCtrlVar);

    // functions for FD1FGeneric and spot dependent spread
    virtual double getCoupon(int step, const double* s, int start, int end)
    { return 0.0; }

    virtual bool hasEquityLayer()
    { return false; }

    // for DDE
    DateTime maxMaturity() const;

    void sensitiveDates(  DateTimeArray    &dates) const;

    void sensitiveStrikes(  const DateTimeArray     dates,
                            DoubleArray             &strikes,   // same dimension as dates
                            bool                    &strikeIsPct) const;

    // temp here
    void AdjustDeltaShift(CControl* control);

    virtual bool GetFwdStartLV() {
        return instrVan->fwdStarting;
    }

    virtual DateTime GetFwdStartDateLV()
    {
        if (instrVan->fwdStarting)
            return instrVan->startDate;
        else
            return instrVan->valueDate;
    }

private:
    const CVanilla* instrVan;
    bool            useInsertNode;
    vector<double>  stepStrike;
    vector<bool>    stepCanExercise;
    bool            isStaticSpread;
};

/** only LN model can use protected asset. Assuming all other models need plain asset  */
CAssetConstSP Vanilla1fProd::GetAssetRef()
{
    return CAssetConstSP::dynamicCast(
        (IObjectConstSP)instrVan->asset.getSP());
}

/** returns a vol request for log-normal vol */
CVolRequestConstSP Vanilla1fProd::GetLNRequest()
{
    // get strike and maturity date from instrument
    DateTime        matDate = instrVan->exerciseSchedule->lastDate();

    double volStrike  = instrVan->exerciseSchedule->lastValue();

    DateTime imntStartDate = instrVan->fwdStarting?
                         instrVan->startDate: instrVan->valueDate;

    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(volStrike, imntStartDate,
                                   matDate, instrVan->fwdStarting));
    return volRequest;
}

/** initialise tree1f - allow product customisation */
void Vanilla1fProd::Init(CControl* control)
{
    static const string method = "Vanilla1fProd::Init";
    try {
        /** customize tree parameters here and set up the tree */
        DateTimeArray segDates;
        segDates.resize(2);
        // this needs change if fwd start tree has to start today !!!
        if (instrVan->fwdStarting && instrVan->startDate>instrVan->valueDate) {
            segDates[0] = instrVan->startDate;
        }
        else {
            segDates[0] = instrVan->valueDate;
        }

        segDates[1] = instrVan->exerciseSchedule->lastDate();
        vector<int> density;
        density.resize(1);
        density[0] = 1;
        double minGap = 0; // insert all points
        bool useEqualTime = false; // false = equal variance steps
        int numOfPriceArray = 2;
        // if not static spread and is put, need to add slice for default prob
        if( !isStaticSpread && !instrVan->isCall ) numOfPriceArray++;

        // all exercise dates are copied to critical dates
        DateTimeArray critDates = instrVan->exerciseSchedule->getDates();
        // remove exercise date from crit date
        critDates.erase(critDates.end()-1);
        // add div event dates if needed
        EventAssetMove divEvent;
        DateTimeArraySP divCritDates;
        // finite difference
        if (instrVan->canExerciseEarly) {
            // American exercise or dollar div interp treatment
            // create div events, this call is very expensive for basket with lots of div dates
            // should only need once for pricing call but need to think how to store/copy for tweaks
            int numDivs = (int)(4*segDates[0].yearFrac(segDates[1]))+1; // 4 divs per year selected as critical dates
            if (AssetUtil::getJumpEvents(instrVan->asset.get(),
                                            segDates[0],
                                            segDates[1],
                                            numDivs,
                                            divEvent)){

                // calculate critical dates
                divCritDates = divEvent.getCritDate(instrVan->noExerciseWindow, instrVan->isCall);

                for (int i = 0; i < divCritDates->size(); i++)
                    critDates.push_back((*divCritDates)[i]);
            }
        }
        // call tree step set up routine
        if( fdModel )
            fdModel->Setup(instrVan->valueDate, segDates, density, &critDates,
                        minGap, useEqualTime, numOfPriceArray);
        else
            genericFDModel->Setup(instrVan->valueDate, segDates, density, &critDates,
                        minGap, useEqualTime, numOfPriceArray);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** initialising and setting product variables
    this is called per pricing call before tree sweep call (after InitTree)  */
void Vanilla1fProd::InitProd() {
    static const string method = "Vanilla1fProd::InitProd";
    try {
        int i;
        int num = model1F->TimePts.NumOfStep;
        double fwdAtStart = 1.0;

        stepStrike.resize(num+1);
        stepCanExercise.resize(model1F->TimePts.NumOfStep+1);
        // ask tree to decide first about steps that can exercise
        AssetUtil::setStepExercise(stepCanExercise,
                                 model1F->TimePts.StepDates,
                                 instrVan->exerciseSchedule,
                                 instrVan->canExerciseEarly,
                                 instrVan->asset.getSP());

        /* If no exercise before ex-div date applies */
        if (instrVan->canExerciseEarly && instrVan->noExerciseWindow > 0) {

            int numSteps = model1F->TimePts.StepDates.size();

            /* Get dividend list */
            DateTime lastDate = AssetUtil::getHoliday(instrVan->asset.get())->addBusinessDays(
                model1F->TimePts.StepDates[numSteps - 1],
                instrVan->noExerciseWindow);

            DividendListSP divList(AssetUtil::getAllDivsBetweenDates(instrVan->asset,
                                                                     model1F->TimePts.StepDates[0],
                                                                     lastDate));

            /* If dividend list not empty */
            if (divList->getArray().size() > 0){
                HolidayConstSP holiday(
                    AssetUtil::getHoliday(instrVan->asset.get()));
                DividendConstSP div;    // null
                DateTime nextDivDate(0, 0);   // past date
                int iStep = 0;
                for(; iStep < numSteps; ++iStep){
                    const DateTime& currStepDate = model1F->TimePts.StepDates[iStep];
                    /* If need to get next ex div date */
                    if (currStepDate.isGreaterOrEqual(nextDivDate)){
                        /* attempt to get next div */
                        div = divList->getNextDivFromDate(currStepDate);
                        /* if there is no next div, break */
                        if (!div){
                            break;
                        }
                        /* otherwise, get ex div date */
                        nextDivDate = div->getExDate();
                    }
                    int daysToNextDiv;
                    /* compute nb of bus days from current time step to next ex date */
                    daysToNextDiv = holiday->businessDaysDiff(currStepDate,
                                                              nextDivDate);
                    /* if current step is exercise but falls within noExerciseWindow days,
                       override exercise to no exercise */
                    if (stepCanExercise[iStep] == true && daysToNextDiv <= instrVan->noExerciseWindow) {
                        stepCanExercise[iStep] = false;
                    }
                }
            }
        }

        bool canInterpExSched = (instrVan->exerciseSchedule->length() > 1);
        // get spot at start if needed
        if (instrVan->fwdStarting) {
            fwdAtStart = instrVan->asset->fwdValue(instrVan->startDate);
        }

        CIntArray  canEx;
        canEx.resize(num+1); // for debug
        canEx[num] = 1;

        // get strikes
        stepStrike[num] = fwdAtStart*instrVan->exerciseSchedule->lastValue();

        for (i=0; i<num; i++) {
            if (stepCanExercise[i]) {
                canEx[i] = 1;
                if (canInterpExSched) {
                    stepStrike[i] = fwdAtStart*instrVan->exerciseSchedule->interpolate(model1F->TimePts.StepDates[i]);
                }
                else {
                    stepStrike[i] = stepStrike[num]; // strike level is constant for american with no schedule
                }
            }
            else {
                canEx[i] = 0;
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** product payoff method at maturity */
void Vanilla1fProd::PayoffAtMat(const double * s, int step, int bot,
                                int top, int pStart, int pEnd,
                                double * const * price) {
    int i, j;

    double settlementPV = instrVan->instSettle->pvAdjust(instrVan->exerciseSchedule->lastDate(),
                                                         instrVan->discount.get(),
                                                         instrVan->asset.get());
    for (i=pStart; i<=pEnd; i++)
    {
        for (j=-bot; j<=top; j++)
        {
            price[i][j] = GetIntrinsic(s[j], stepStrike[step], instrVan->isCall, true/* this allow fwd */);
            price[i][j] *= settlementPV;
        }
    }

    if( !isStaticSpread && !instrVan->isCall )
    {
        FD1FDDE *fdDDE = dynamic_cast<FD1FDDE *>(genericFDModel);
        fdDDE->setDefaultProb(step, price[2], -bot, top);
    }
}

// SmoothMax routine
void VanillaSMax(int bot, int top, double strike, int callPutMult,
                 double settlePV, double* option, const double* under)
{
    const double h = 0.5; // can discuss all day what's best for this one
    double diffUp = 0;
    double valueDown, diffDown, hMaxDiff;
    double valueMid = option[top] - under[top]*callPutMult;
    for (int idx=top; idx > bot; idx--)
    {
        valueDown = option[idx-1] - under[idx-1]*callPutMult;
        diffDown = fabs(valueMid - valueDown);
        hMaxDiff = h*Maths::max(diffUp, diffDown);
        option[idx] = SMax((under[idx]-strike)*callPutMult, option[idx], hMaxDiff);
        diffUp = diffDown;
        valueMid = valueDown;
    }
    hMaxDiff = diffUp;
    option[bot] = settlePV*SMax((under[bot]-strike)*callPutMult, option[bot], hMaxDiff);
}


/** product payoff method at steps earlier than maturity */
// just to support old FD for now
void Vanilla1fProd::PayoffBeforeMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                                      double * const * price)
{
    static const string method = "Vanilla1fProd::PayoffBeforeMatFD";
    try {
        int j;

        ASSERT(pStart==0 && pEnd==(1+(!isStaticSpread && !instrVan->isCall)));

        if (step == model1F->TimePts.NumOfStep-1 && isStaticSpread ) {
            // do penultimate smoothing
            vector <double> vol_arr, drift_arr;
            double dt = model1F->TimePts.TradeYrFrac[step+1];
            model1F->GetStepVol(step, vol_arr, s, -bot, top);

            double variance = vol_arr[0]*vol_arr[0]*dt;
            double drift = instrVan->asset->fwdValue(model1F->TimePts.StepDates[step + 1])/instrVan->asset->fwdValue(model1F->TimePts.StepDates[step]);
            double df = instrVan->discount->pv(model1F->TimePts.StepDates[step],
                                               model1F->TimePts.StepDates[step+1]);
            for (j=-bot; j<=top; j++)
            {

                if (s[j]>0.0 && fabs(log(s[j]/stepStrike[step+1])) < model1F->TruncationStd*sqrt(variance))
                {

                    price[0][j] = price[1][j] = Black::price(instrVan->isCall, drift*s[j],
                                                             stepStrike[step+1], df, variance);
                }
            }
        }

        if( !isStaticSpread && !instrVan->isCall )
        {
            FD1FDDE *fdDDE = dynamic_cast<FD1FDDE *>(genericFDModel);
            // for american, find next exercisable point to get the def payoff
            int strkIdx=step+1;
            while( strkIdx<model1F->TimePts.NumOfStep && !stepCanExercise[strkIdx] ) strkIdx++;
            double pv = instrVan->discount->pv(model1F->TimePts.StepDates[step + 1],
                                               model1F->TimePts.StepDates[strkIdx]);
            fdDDE->adjustPriceByDefPO(step+1, stepStrike[strkIdx] * pv, price[2], price[0], -bot, top);

            // for european, the default value is IR discounted strike from mat
            strkIdx = model1F->TimePts.NumOfStep;
            pv = instrVan->discount->pv(model1F->TimePts.StepDates[step + 1],
                                        model1F->TimePts.StepDates[strkIdx]);
            fdDDE->adjustPriceByDefPO(step+1, stepStrike[strkIdx] * pv, price[2], price[1], -bot, top);

            // prep for the next time step
            if( step ) fdDDE->setDefaultProb(step, price[2], -bot, top);
        }

        if (!stepCanExercise[step]) {
            return;
        }

        double settlementPV = instrVan->instSettle->pvAdjust(model1F->TimePts.StepDates[step],
                                                             instrVan->discount.get(),
                                                             instrVan->asset.get());


        double callput = settlementPV*(instrVan->isCall? 1.0:-1.0);

        double intrinsic = 0.0;

        for (j=-bot; j<=top; j++) {
            intrinsic = callput*(s[j] - stepStrike[step]);

            if (price[0][j] < intrinsic) {
                price[0][j] = intrinsic; // American
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** control variate */
double Vanilla1fProd::RefinePrice(double basePrice, double fwdStartDF, bool useCtrlVar)
{
    double price;

    if (useCtrlVar) {
        CVanillaClosedForm closedForm(instrVan);
        CClosedFormLN model;
        Control ctrl;
        CResults result;
        closedForm.price(&model, &ctrl, &result);
        double closedFormPrice = result.retrievePrice();
        double exer_premium = model1F->PriceEnd[0] - model1F->PriceEnd[1];

        price = scalePremium(exer_premium*fwdStartDF) + closedFormPrice;
    }
    else {
        price = fwdStartDF*scalePremium(basePrice);
    }
    return price;
}

/** premium scaling */
double Vanilla1fProd::scalePremium(const double& fairValue)
{
    double fwdAtStart = 0.0;
    if (instrVan->fwdStarting)
    {
        fwdAtStart = instrVan->asset->fwdValue(instrVan->startDate);
    }

    double scalingFactor = InstrumentUtil::scalePremium(
                                    instrVan->oneContract,
                                    instrVan->fwdStarting,
                                    instrVan->notional,
                                    fwdAtStart,
                                    instrVan->initialSpot);

    return fairValue * scalingFactor;
}


void Vanilla1fProd::recordOutputRequests(Control* control, Results* results, double fairValue)
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime       matDate = instrVan->exerciseSchedule->lastDate();
        double         indVol;
        // calculate indicative vol
        if ( matDate.isGreater(instrVan->valueDate) )
        {
            DateTime imntStartDate = instrVan->fwdStarting?
                             instrVan->startDate: instrVan->valueDate;

            // get vol request
            CVolRequestConstSP lnVolRequest = GetLNRequest();

            try{
                // interpolate the vol
                CVolProcessedSP  vol(instrVan->asset->getProcessedVol(lnVolRequest.get()));
                // cast to the BS vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);

                // this should never happen if our get market data has worked properly
                if (!volBS){
                    throw ModelException("Vanilla1fProd::recordOutputRequests",
                                         "No Black Scholes Vol");
                }
                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
            }
            catch (exception& ) {
                indVol = 0.0;
            }
        }
        else
        {
            indVol = 0.0;
        }

        instrVan->addOutputRequests(control,
                                        results,
                                        fairValue,
                                        indVol);
    }
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* CVanilla::createProduct(
    CClosedFormLN* model) const{

    // there is some additional validation which is dependent
    // on the model type here
    if ( exerciseSchedule->length() > 1 )
    {
        throw ModelException("Vanilla::createProduct",
         "Closed form vanilla model can only handle single exercise dates.");
    }

    if ( canExerciseEarly )
    {
        throw ModelException("Vanilla::createProduct",
         "Closed form vanilla model can not handle early exercise feature.");
    }

    return new CVanillaClosedForm(this);
}

/** create a fd payoff product */
FD1F::IProduct* CVanilla::createProduct(FD1F* model) const
{

    Vanilla1fProd* fdProd = new Vanilla1fProd(this);
    fdProd->fdModel = model;
    fdProd->model1F = model;

    return fdProd;
}

/** Implementation of DDEInitiator interface, built on FD1FGeneric
 */
    DateTime Vanilla1fProd::maxMaturity() const { return instrVan->exerciseSchedule->lastDate(); }

    void Vanilla1fProd::sensitiveDates(  DateTimeArray    &dates) const
    {
        dates = instrVan->exerciseSchedule->getDates();
    }

    void Vanilla1fProd::sensitiveStrikes(  const DateTimeArray     dates,
                            DoubleArray             &strikes,   // same dimension as dates
                            bool                    &strikeIsPct) const
    {
        if( dates.size() != strikes.size() )
            throw ModelException("Vanilla1fProdDDE::sensitiveStrikes", "Dates and strikes dimension mismatch");

        /************ doesn't support fwd start ***************/
        if ( instrVan->fwdStarting )
            throw ModelException("Vanilla1fProdDDE::sensitiveStrikes", "Forward starting option not supported under DDE");

        double strike = instrVan->exerciseSchedule->lastValue();
        for(int i=dates.size()-1; i>=0; i--)
            strikes[i] = strike;

        strikeIsPct = false;
    }

/** create a fd payoff product */
FD1FGeneric::IProduct* CVanilla::createProduct(FD1FGeneric* model) const
{
    if(FD1FDDE::TYPE->isInstance(model))
    {
        // don't support fwd starting
        if( fwdStarting )
            throw ModelException("Vanilla::createProduct", "Do not support fwdStarting under DDE");
        if( ccyTreatment == CAsset::CCY_TREATMENT_STRUCK || ccyTreatment == CAsset::CCY_TREATMENT_PROTECTED )
            throw ModelException("Vanilla::createProduct", "Do not support ccy feature under DDE");
    }
    else
        throw ModelException("Vanilla::createProduct", "Do not support non-DDE FD1FGeneric");

    Vanilla1fProd* fdProd = new Vanilla1fProd(this);
    fdProd->genericFDModel = model;
    fdProd->model1F = model;
    if(FD1FDDE::TYPE->isInstance(model)) fdProd->isStaticSpread = false;

    return fdProd;
}


/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool CVanilla::avoidVegaMatrix(const IModel* model)
{
    /* this should possibly be false for local vol models etc. */
    return false;
}

/** returns all strikes on the vol surface to which
    this instrument is sensitive */
DoubleArraySP CVanilla::getSensitiveStrikes(OutputNameConstSP outputName,
                                            const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("Vanilla::getSensitiveStrikes",
                             "VEGA_MATRIX is not valid for this instrument");
    }

    // get start date for vol interpolation
    DateTime imntStartDate = fwdStarting ? startDate:valueDate;
    // get last exercise date in exercise schedule
    DateTime maturityDate = exerciseSchedule->lastDate();
    // get last strike in exercise schedule
    double   strike       = exerciseSchedule->lastValue();

    // create a vol request object to be passed on
    LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(strike,
                                                                   imntStartDate,
                                                                   maturityDate,
                                                                   fwdStarting));

    SensitiveStrikeDescriptor sensStrikeDesc;
    sensStrikeDesc.forwardOnly = false;

    asset->getSensitiveStrikes(volRequest.get(), outputName,
                               sensStrikeDesc, sensStrikes);

    return sensStrikes;
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool CVanilla::sensShift(Theta* shift)
{
    DateTime newDate = shift->rollDate(valueDate);

    DateTime matDate = exerciseSchedule->lastDate();

    if ( ( newDate >= matDate && valueDate < matDate ) ||
         ( valueDate == matDate && Maths::isZero(spotAtMaturity)))
        spotAtMaturity = asset->getThetaSpotOnDate(shift, matDate);

    if (fwdStarting && newDate.isGreaterOrEqual(startDate) &&
        startDate.isGreaterOrEqual(valueDate))
    {
        fwdStarting = false;
        initialSpot = asset->getThetaSpotOnDate(shift, startDate);
        exerciseSchedule->scale(initialSpot);
    }

    // roll today
    valueDate = newDate;

    return true;
};

double CVanilla::priceSpread(const DateTime& valueDate,
                             const DateTime& startDate,
                             const DateTime& matDate,
                             bool isCall,
                             bool fwdStarting,
                             bool oneContract,
                             double notional,
                             double initialSpot,
                             double lowStrike,
                             double highStrike,
                             const InstrumentSettlement* instSettle,
                             const Asset* asset,
                             const YieldCurve* discount)
{
    static const string method = "Vanilla::priceSpread";

    /* if it's forward starting, convert payoffToProcFreqBoundWeight strikes
       to absolute values based on spot at start date */
    DateTime imntStartDate = fwdStarting? startDate: valueDate;

    double fwdAtStart;
    double lowAbsStrike = lowStrike;
    double highAbsStrike = highStrike;
    if (fwdStarting)
    {
        fwdAtStart = asset->fwdValue(startDate);
        lowAbsStrike *= fwdAtStart;
        highAbsStrike *= fwdAtStart;
    }

    // get the fwd price - only need to do this once
    double fwdPrice = asset->fwdValue(matDate);

    // choose how to interpolate the vol - go for traditional route for now
    LinearStrikeVolRequestSP lowVolRequest(new LinearStrikeVolRequest(
                                                lowStrike,
                                                imntStartDate,
                                                matDate,
                                                fwdStarting));

    LinearStrikeVolRequestSP highVolRequest(new LinearStrikeVolRequest(
                                                highStrike,
                                                imntStartDate,
                                                matDate,
                                                fwdStarting));

    // interpolate the vol for each strike
    CVolProcessedSP lowVol(asset->getProcessedVol(lowVolRequest.get()));
    CVolProcessedSP highVol(asset->getProcessedVol(highVolRequest.get()));

    // cast to the type of vol we're expecting
    CVolProcessedBSSP lowVolBS = CVolProcessedBSSP::dynamicCast(lowVol);
    CVolProcessedBSSP highVolBS = CVolProcessedBSSP::dynamicCast(highVol);

    // this should never happen if our get market data has worked properly
    if (!lowVol || !highVol)
    {
        throw ModelException(method, "No Black Scholes Vol");
    }

    // calculate the variance
    double lowVariance = lowVolBS->CalcVar(imntStartDate, matDate);
    double highVariance = highVolBS->CalcVar(imntStartDate, matDate);

    // calculate the discount factor back to today
    double discFactor = instSettle->pv(valueDate,
                                       matDate,
                                       discount,
                                       asset);

    // finally call Black model for each strike
    double lowStrikePremium = Black::price(isCall, fwdPrice,
                                           lowAbsStrike, discFactor, lowVariance);

    double highStrikePremium = Black::price(isCall, fwdPrice,
                                            highAbsStrike, discFactor, highVariance);

    double scalingFactor = InstrumentUtil::scalePremium(oneContract,
                                                        fwdStarting,
                                                        notional,
                                                        fwdAtStart,
                                                        initialSpot);


    lowStrikePremium *= scalingFactor;
    highStrikePremium *= scalingFactor;


    // Now combine the premiums
    double spreadPrice = isCall?
        lowStrikePremium - highStrikePremium : /* call spread */
        highStrikePremium - lowStrikePremium;  /* put spread */

    return spreadPrice;
}


void CVanilla::avoidVegaMatrixLite(const IModel* model) {
    static const string method = "CVanilla::avoidVegaMatrixLite";
    
    if (avoidVegaMatrix(model)) {
        // Check basic VEGA_MATRIX
        throw ModelException(method, "Instrument does not support VEGA_MATRIX and hence not VEGA_MATRIX_LITE");
    } else if (!SimpleEquity::TYPE->isInstance(asset.get())) {
        // Allow LITE only for SimpleEquity
        throw ModelException(method, "Only SimpleEquity underlyings supported");
    } else if(fwdStarting) {
        // Cannot reproduce Vega for fwd starting options
        throw ModelException(method, "Fwd starting options not supported");
    } else if(canExerciseEarly) {
        // Cannot reproduce Vega for american options
        throw ModelException(method, "Early exercise options not supported");
    }
}



// Just like the call spread routine except it prices only one option.
// Useful for other products which include groups of vanilla options (e.g. KI Fwd)
// Perhaps the spread routine should be made to call this one twice. Would have
// some redundant calcs but that's not too big a price to pay for consistancy.
double CVanilla::priceBS(const DateTime& valueDate,
                         const DateTime& startDate,
                         const DateTime& matDate,
                         bool isCall,
                         bool fwdStarting,
                         bool oneContract,
                         double notional,
                         double initialSpot,
                         double strike,
                         const InstrumentSettlement* instSettle,
                         const Asset* asset,
                         const YieldCurve* discount)
{
    static const string method = "Vanilla::priceBS";

    /* if it's forward starting, convert payoffToProcFreqBoundWeight strikes
       to absolute values based on spot at start date */
    DateTime imntStartDate = fwdStarting? startDate: valueDate;

    double fwdAtStart;
    double absStrike = strike;
    if (fwdStarting)
    {
        fwdAtStart = asset->fwdValue(startDate);
        absStrike *= fwdAtStart;
    }

    // get the fwd price - only need to do this once
    double fwdPrice = asset->fwdValue(matDate);

    // choose how to interpolate the vol - go for traditional route for now
    LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(
                                        strike,
                                        imntStartDate,
                                        matDate,
                                        fwdStarting));

    // interpolate the vol for each strike
    CVolProcessedSP vol(asset->getProcessedVol(volRequest.get()));

    // cast to the type of vol we're expecting
    CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);

    // this should never happen if our get market data has worked properly
    if (!volBS)
    {
        throw ModelException(method, "No Black Scholes Vol");
    }

    // calculate the variance
    double variance = volBS->CalcVar(imntStartDate, matDate);

    // calculate the discount factor back to today
    double discFactor = instSettle->pv(valueDate,
                                       matDate,
                                       discount,
                                       asset);

    // finally call Black model
    double strikePremium = Black::price(isCall, fwdPrice,
                                        absStrike, discFactor, variance);

    double scalingFactor = InstrumentUtil::scalePremium(oneContract,
                                                        fwdStarting,
                                                        notional,
                                                        fwdAtStart,
                                                        initialSpot);

    strikePremium *= scalingFactor;

    return strikePremium;
}


CSensControl* CVanilla::AlterControl(const IModel*       modelParams,
                                     const CSensControl* sensControl) const
{
    SensControlPerName* alteredControl = NULL;
    if (Delta::TYPE->isInstance(sensControl)         &&
        CClosedFormLN::TYPE->isInstance(modelParams) )
    {
       const Delta* delta =
            dynamic_cast<const Delta*>((IObject*)sensControl);
        double  strike  = exerciseSchedule->lastValue();
        ShiftSizeCollectorSP shiftSizeVisitor(new ShiftSizeCollector(
            delta,
            strike,
            fwdStarting?
            ShiftSizeCollector::FWD_START_ADJUSTMENT:
            ShiftSizeCollector::SPOT_START_ADJUSTMENT));

        asset->accept(shiftSizeVisitor.get());

        if ( Maths::isPositive(shiftSizeVisitor->getShiftSize()) )
        {
            alteredControl = new Delta(shiftSizeVisitor->getShiftSize());
            alteredControl->
                setMarketDataName(sensControl->getMarketDataName());
        }
    }

    return alteredControl;
}

// for reflection
CVanilla::CVanilla():
CInstrument(TYPE),isCall(false),
fwdStarting(false),
canExerciseEarly(false),
oneContract(true),
notional(0.0),
initialSpot(0.0),
spotAtMaturity(0.0),
isExercised(false),
noExerciseWindow(0){}

// Monte Carlo
class VanillaMC : public IMCProduct, virtual public IMCProductLN {
private:
    const CVanilla*     inst;
    double              strike;
    double              mult;
    bool                useRatio;
public:

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    // equivalent to InstIntoMCProduct
    VanillaMC(const CVanilla*           inst,
              const IRefLevelConstSP&   refLevel,    // how to 'avg in'
              const SimSeriesConstSP&   simSeries,   // simulation dates
              const IPastValuesConstSP& mcPastValues, // historic values
              const DateTime&           matDate):    // maturity date
    IMCProduct(inst->asset.get(),
              inst->valueDate,
              inst->discount.get(),
              refLevel,
              simSeries,
              mcPastValues, // historic values
              inst->instSettle.get(),
              matDate),
    inst(inst) {
        strike = inst->exerciseSchedule->lastValue();
        if (inst->fwdStarting) {
            // forward starting => notional based, and % strike
            useRatio = true;
            mult = inst->notional;
        } else {
            useRatio = false;
            mult = inst->oneContract? 1.0: inst->notional / inst->initialSpot;
        }
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        double inValue = pathGen->refLevel(0,0);
        double outValue = 0.0;
        for(int iStep = pathGen->begin(0); iStep < pathGen->end(0); iStep++) {
            outValue += pathGen->Path(0,0)[iStep];
        }

        double price = useRatio ? (outValue / inValue - strike): (outValue - strike);
        if (!inst->isCall){
            price = -price;
        }
        prices.add(price < 0.0? 0.0: (price * mult));
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {

        DateTime imntStartDate = inst->fwdStarting?
            inst->startDate: inst->valueDate;
        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here

        reqarr[0] = CVolRequestLNSP(new LinearStrikeVolRequest(
            inst->exerciseSchedule->lastValue(),
            imntStartDate,
            inst->exerciseSchedule->lastDate(),
            inst->fwdStarting));

        return reqarr;
    }

};

// State var compliant Monte Carlo product
class VanillaMCSV : public MCProductClient,
                    virtual public IMCProductLN {
private:
    const CVanilla*     inst;
    double              strike;
    double              mult;
    bool                useRatio;
    // state vars
    SVGenSpotSP                  spotGen;      //!< Generator for spot
    SVGenSpot::IStateVarSP       spotSV;       //!< Spot state variable
    IRefLevel::IStateVarGenSP refLevelGen;  //!< Generator for ref level
    IRefLevel::IStateVarSP    refLevelSV;   //!< Ref level state variable
    SVGenDiscFactorSP            dfGen;        //!< Generator for discount factors
    SVDiscFactorSP dfSV;         //!< Df state variable

    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        static const string routine = "VanillaMCSV::pathGenUpdated";

        try {
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };

public:
    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        // ask for a reference level State Variable
        svCollector->append(refLevelGen.get());
        svCollector->append(spotGen.get());
        svCollector->append(dfGen.get());
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    // equivalent to InstIntoMCProduct
    VanillaMCSV(const CVanilla*           inst,
                const IRefLevelConstSP&   refLevel,     // how to 'avg in'
                const SimSeriesConstSP&   simSeries,    // simulation dates
                const IPastValuesConstSP& mcPastValues, // historic values
                const DateTime&           matDate):     // maturity date
    MCProductClient(inst->asset.get(),
                    inst->valueDate,
                    inst->discount.get(),
                    refLevel,
                    simSeries,
                    mcPastValues, // historic values
                    inst->instSettle.get(),
                    matDate),
    inst(inst),
    spotGen(new SVGenSpot(simSeries)),
    refLevelGen(refLevel->createStateVarGen(getMultiFactors(), inst->valueDate)),
    dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                           inst->instSettle, simSeries->getLastDate())) {
        strike = inst->exerciseSchedule->lastValue();
        if (inst->fwdStarting) {
            // forward starting => notional based, and % strike
            useRatio = true;
            mult = inst->notional;
        } else {
            useRatio = false;
            mult = inst->oneContract? 1.0: inst->notional / inst->initialSpot;
        }
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        int iAsset = 0; // only 1 asset
        double inValue = refLevelSV->refLevel(iAsset);
        const SVPath& path = spotSV->path(iAsset);
        double outValue = 0.0;
        // NB: might be the case that the payoff method is called for the past
        // even though 'path' is empty. That is why we need this test even though
        // path contains 1 date at most
        for (int iStep = path.begin(); iStep < path.end(); ++iStep){
            outValue = path[iStep];
        }
        double price = useRatio ? (outValue / inValue - strike): (outValue - strike);
        if (!inst->isCall){
            price = -price;
        }

        // Discount
        price *= dfSV->firstDF();

        prices.add(price < 0.0? 0.0: (price * mult));
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {
        DateTime imntStartDate = inst->fwdStarting ?
            inst->startDate: inst->valueDate;
        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here

        reqarr[0] = CVolRequestLNSP(new LinearStrikeVolRequest(
            inst->exerciseSchedule->lastValue(),
            imntStartDate,
            inst->exerciseSchedule->lastDate(),
            inst->fwdStarting));

        return reqarr;
    }

};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* CVanilla::createProduct(const MonteCarlo* model) const {
    static const string routine("Vanilla::createProduct");
    if (exerciseSchedule->length() > 1) {
        throw ModelException(routine, "MC vanilla model can only handle"
                             " single exercise dates.");
    }

    if (canExerciseEarly) {
        throw ModelException(routine, "MC vanilla model can not"
                             " handle early exercise feature.");
    }

    // v simple simSeries
    SimSeriesSP simSeries(new SimSeries(1));
    const DateTime& matDate = exerciseSchedule->lastDate();
    simSeries->addDates(DateTimeArray(1, matDate));
    // v simple RefLevel
    DateTime refDate = fwdStarting? startDate : valueDate;
    IRefLevelSP refLevel(IRefLevel::Util::makeFwdStart(refDate));
    // need to add start date and maturity to past values
    DateTimeArray  allDates(1, refDate);
    allDates.push_back(matDate);
    DoubleArray    allValues(1, initialSpot);
    allValues.push_back(spotAtMaturity);
    DoubleMatrix allValuesAsMatrix(allValues);
    IPastValuesSP pastValues(
        IPastValues::Util::makeSimple(allDates, allValuesAsMatrix));
    // if state vars requested
    if (model->stateVarUsed()){
        return new VanillaMCSV(this, refLevel, simSeries,
                               pastValues, matDate);
    }
    // otherwise, use old methodology
    return new VanillaMC(this, refLevel, simSeries,
                        pastValues, matDate);
}

class VanillaISAPHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VanillaISAP, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FourierEngine::ISAP);
        EMPTY_SHELL_METHOD(defaultVanillaISAP);
        FIELD(useOneIntegral, "Use 1-integral approach if true; use 2-integral approach otherwise");
        FIELD_MAKE_OPTIONAL(useOneIntegral);
        FIELD(payoffToProcFreqBoundWeight, "Determines the imaginary line of integration. "
                                                  "The closer to 1.0, the closer to the payoff's relevant frequency boundary");
        FIELD_MAKE_OPTIONAL(payoffToProcFreqBoundWeight);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultVanillaISAP(){
        return new VanillaISAP();
    }
};

CClassConstSP const VanillaISAP::TYPE = CClass::registerClassLoadMethod(
    "VanillaISAP", typeid(VanillaISAP), VanillaISAPHelper::load);

/** Single integrand - Default integrator case */
template <class Process, class Product>
class VanillaSingleIntegrand: public Function1DDouble {
public:
    VanillaSingleIntegrand(const Process&   process,
                    const Product&   product,
                    double           strike,
                    double           fwd,
                    const DateTime&  maturity,
                    double           omega):
    Function1DDouble(Range(OpenBoundary(0.0), Infinity(Infinity::Plus))),
    logMoneyness(log(strike/fwd)),
    omega(omega),
    matdate(maturity),
    process(process),
    product(product){}

    virtual double operator()(double  u) const {    // u == frequency
        Complex  z(omega, u);
        Complex  Laplace = exp(process.scalelessCumulant(product, z, matdate) - z * logMoneyness) / (z * (z - 1.0));
        return Laplace.real();
    }

private:
    double logMoneyness;
    double omega;
    const DateTime& matdate;
    const Process& process;
    const Product& product;
};

/** Single integrand - FFT case */
template <class Process, class Product>
class VanillaSingleIntegrandFFT: public Function1DComplex {
public:
    VanillaSingleIntegrandFFT(const Process&   process,
                       const Product&   product,
                       const DateTime&  maturity,
                       double           omega):
    // Infinite range by default
    omega(omega),
    matdate(maturity),
    process(process),
    product(product){}

    const Range& getInterval() const {return this->interval;}

    virtual Complex operator()(double  u) const {    // u == frequency
        Complex  z(omega, u);
        Complex  Laplace = exp(process.scalelessCumulant(product, z, matdate)) / (z * (z - 1.0));
        return Laplace;
    }

private:
    double omega;
    const DateTime& matdate;
    const Process& process;
    const Product& product;
};

/** Double integrand - Default integrator case */
template <class Process, class Product>
class VanillaDoubleIntegrand: public Function1DDouble {
public:
    VanillaDoubleIntegrand(const Process&   process,
                    const Product&   product,
                    double           strike,
                    double           fwd,
                    const DateTime&  maturity,
                    int              j):    // 0 or 1
    Function1DDouble(Range(OpenBoundary(0.0), Infinity(Infinity::Plus))),
    logMoneyness(log(strike/fwd)),
    matdate(maturity),
    j(j),
    process(process),
    product(product){}

    virtual double operator()(double  u) const {    // u == frequency
        Complex  z(0.0, u); // = i * u
        Complex  Laplace = exp(process.scalelessCumulant(product, z + j, matdate) - z * logMoneyness) / z;
        return Laplace.real();
    }

private:
    double logMoneyness;
    const DateTime& matdate;
    double j;
    const Process& process;
    const Product& product;
};

/** Double integrand - FFT integrator case */
template <class Process, class Product>
class VanillaDoubleIntegrandFFT: public Function1DComplex {
public:
    VanillaDoubleIntegrandFFT(const Process&   process,
                       const Product&   product,
                       const DateTime&  maturity,
                       int              j):    // 0 or 1
    // Infinite range by default
    matdate(maturity),
    j(j),
    process(process),
    product(product){}

    virtual Complex operator()(double  u) const {    // u == frequency
        if (Maths::isZero(u)){
            throw ModelException("VanillaDoubleIntegrandFFT::operator(double)",
                                 "Zero frequency is not supported yet");
        }
        Complex  z(0.0, u); // = i * u
        Complex  Laplace = exp(process.scalelessCumulant(product, z + j, matdate)) / z;
        return Laplace;
    }

private:
    const DateTime& matdate;
    double j;
    const Process& process;
    const Product& product;
};

/** Fourier Product */
class VanillaFP: public FourierProduct,
                 public FourierProductIntegrator1D,     // supports default 1D integrators
                 public FourierProductFFTIntegrator1D,  // supports FFT 1D integrators
                 public StFourierProductLogRtn,         // requires a StFourierProcessLogRtn
                 public FwdStFourierProductLogRtn{      // requires a FwdStFourierProcessLogRtn
public:
    // equivalent to InstIntoFourierProduct
    VanillaFP(const CVanilla*          inst,
              const DateTime&          matDate):    // maturity date
    FourierProduct(inst->asset.get(),
                   inst->valueDate,
                   inst->discount.get(),
                   inst->instSettle.get()),
    inst(inst),
    maturity(matDate),
    paymentDate(instSettle->settles(matDate, inst->asset.get())),
    strike(inst->exerciseSchedule->lastValue()){
        if (inst->fwdStarting) {
            // forward starting => notional based, and % strike
            mult = inst->notional;
        }
        else {
            mult = inst->oneContract ? 1.0: inst->notional / inst->initialSpot;
        }
    }

    virtual void price(const FourierEngine* model,
                       Control*             control,
                       Results*             results){
        recordFwdAtMat(control,
                       maturity,
                       results);
        FourierProduct::price(model,
                              control,
                              results);
        if (control && control->isPricing() && useOneIntegral){
            OutputNameConstSP omegaOutput(new OutputName("omega"));
            results->storeScalarGreek(omega, Results::DEBUG_PACKET, omegaOutput);
        }

    }

    /** Constructs default 1D integrands for Vanilla */
    virtual Function1DDoubleArrayConstSP Integrand(const FourierEngine* model,
                                                   const Integrator1D*  integrator) {
        static const string method = "VanillaFP::Integrand";
        try{
            const VanillaISAP* isap = dynamic_cast<const VanillaISAP*>(&model->getISAP());
            if(!isap) {
                throw ModelException(method, "FourierEngine needs to be provided with a VanillaISAP." );
            }
            useOneIntegral = isap->useOneIntegral;
            isCall = inst->isCall || !useOneIntegral;

            if (inst->fwdStarting) {
                DateTimeArray dates(2);
                dates[0] =  getStartDate(); dates[1] = maturity;
                CDoubleArray fwds(2);
                mAsset->assetFwdValue(0, dates, fwds);
                fwd = fwds[1] / fwds[0];
            }
            else {
                fwd = mAsset->assetFwdValue(0, maturity);
            }

            int nbIntegrals = useOneIntegral ? 1 : 2;
            Function1DDoubleArraySP functions(new Function1DDoubleArray(nbIntegrals));

            if (useOneIntegral) {
                if (inst->fwdStarting) {
                    const FwdStFourierProductLogRtn& thisProd = *this;
                    const FwdStFourierProcessLogRtn* thisProc = dynamic_cast<const FwdStFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support FwdStFourierProcessLogRtn interface");
                    }

                    omega = FourierProduct::intersectRanges(isCall ? callLowerBound : putLowerBound,
                                                            isCall ? callUpperBound : putUpperBound,
                                                            thisProc->lowerRealBound(thisProd, maturity),
                                                            thisProc->upperRealBound(thisProd, maturity),
                                                            isCall ? isap->payoffToProcFreqBoundWeight : 1.0 - isap->payoffToProcFreqBoundWeight);

                    (*functions)[0] = Function1DDoubleSP(new VanillaSingleIntegrand<FwdStFourierProcessLogRtn, FwdStFourierProductLogRtn>
                                                                            (*thisProc,
                                                                             thisProd,
                                                                             strike,
                                                                             fwd,
                                                                             maturity,
                                                                             omega));
                }
                else {  // started
                    const StFourierProductLogRtn& thisProd = *this;
                    const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support StFourierProcessLogRtn interface");
                    }

                    omega = FourierProduct::intersectRanges(isCall ? callLowerBound : putLowerBound,
                                                            isCall ? callUpperBound : putUpperBound,
                                                            thisProc->lowerRealBound(thisProd, maturity),
                                                            thisProc->upperRealBound(thisProd, maturity),
                                                            isCall ? isap->payoffToProcFreqBoundWeight : 1.0 - isap->payoffToProcFreqBoundWeight);

                    (*functions)[0] = Function1DDoubleSP(new VanillaSingleIntegrand<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                            (*thisProc,
                                                                             thisProd,
                                                                             strike,
                                                                             fwd,
                                                                             maturity,
                                                                             omega));
                }
            }
            else {  // !useOneIntegral
                if (inst->fwdStarting) {
                    const FwdStFourierProductLogRtn& thisProd = *this;
                    const FwdStFourierProcessLogRtn* thisProc = dynamic_cast<const FwdStFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support FwdStFourierProcessLogRtn interface");
                    }

                    int iIntegral = 0;
                    for (; iIntegral < nbIntegrals; ++iIntegral){
                        (*functions)[iIntegral] = Function1DDoubleSP(new VanillaDoubleIntegrand<FwdStFourierProcessLogRtn, FwdStFourierProductLogRtn>
                                                                                (*thisProc,
                                                                                 thisProd,
                                                                                 strike,
                                                                                 fwd,
                                                                                 maturity,
                                                                                 iIntegral));
                    }
                }
                else {  // started
                    const StFourierProductLogRtn& thisProd = *this;
                    const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support StFourierProcessLogRtn interface.");
                    }

                    int iIntegral = 0;
                    for (; iIntegral < nbIntegrals; ++iIntegral){
                        (*functions)[iIntegral] = Function1DDoubleSP(new VanillaDoubleIntegrand<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                                (*thisProc,
                                                                                 thisProd,
                                                                                 strike,
                                                                                 fwd,
                                                                                 maturity,
                                                                                 iIntegral));
                    }
                }
            }
            return functions;
        }
        catch (exception& e){
            throw ModelException(e, method, "Failed to construct integrand(s)");
        }
    }

    /** Post process method for default integrator */
    virtual void postResults(const FourierEngine* model,
                             const Integrator1D*  integrator,
                             const FourierProductIntegrator1D::IntegralArray& integrals,
                             CControl*            control,
                             CResults*            results) {
        static const string method = "VanillaFP::postResults";
        try {
            double price = 0.0;

            // calculate the discount factor back to today
            double pv = inst->instSettle->pv(today,
                                             paymentDate,
                                             inst->discount.get(),
                                             inst->asset.get());

            if (useOneIntegral){
                price = integrals[0] * strike / Maths::PI;
            }
            else{   // !useOneIntegral
                price = (fwd * integrals[1] - strike * integrals[0]) / Maths::PI;
                if (inst->isCall) {
                    price += 0.5 * (fwd - strike);
                }
                else {  // put
                    price -= 0.5 * (fwd - strike);
                }
            }
            price *= mult * pv;
            results->storePrice(price, discount->getCcy());
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Constructs FFT 1D integrands for Vanilla */
    virtual Function1DComplexArrayConstSP Integrand(const FourierEngine*    model,
                                                    const FFTIntegrator1D*  integrator) {
        static const string method = "VanillaFP::Integrand";
        try{
            const VanillaISAP* isap = dynamic_cast<const VanillaISAP*>(&model->getISAP());
            if(!isap) {
                throw ModelException(method, "FourierEngine needs to be provided with a VanillaISAP." );
            }
            useOneIntegral = isap->useOneIntegral;
            isCall = inst->isCall || !useOneIntegral;

            if (inst->fwdStarting) {
                DateTimeArray dates(2);
                dates[0] =  getStartDate(); dates[1] = maturity;
                CDoubleArray fwds(2);
                mAsset->assetFwdValue(0, dates, fwds);
                fwd = fwds[1] / fwds[0];
            }
            else {
                fwd = mAsset->assetFwdValue(0, maturity);
            }

            int nbIntegrals = useOneIntegral ? 1 : 2;
            Function1DComplexArraySP functions(new Function1DComplexArray(nbIntegrals));

            if (useOneIntegral) {
                if (inst->fwdStarting) {
                    const FwdStFourierProductLogRtn& thisProd = *this;
                    const FwdStFourierProcessLogRtn* thisProc = dynamic_cast<const FwdStFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support FwdStFourierProcessLogRtn interface");
                    }

                    omega = FourierProduct::intersectRanges(isCall ? callLowerBound : putLowerBound,
                                                            isCall ? callUpperBound : putUpperBound,
                                                            thisProc->lowerRealBound(thisProd, maturity),
                                                            thisProc->upperRealBound(thisProd, maturity),
                                                            isCall ? isap->payoffToProcFreqBoundWeight : 1.0 - isap->payoffToProcFreqBoundWeight);

                    (*functions)[0] = Function1DComplexSP(new VanillaSingleIntegrandFFT<FwdStFourierProcessLogRtn, FwdStFourierProductLogRtn>
                                                                            (*thisProc,
                                                                             thisProd,
                                                                             maturity,
                                                                             omega));
                }
                else {  // started
                    const StFourierProductLogRtn& thisProd = *this;
                    const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support StFourierProcessLogRtn interface");
                    }

                    omega = FourierProduct::intersectRanges(isCall ? callLowerBound : putLowerBound,
                                                            isCall ? callUpperBound : putUpperBound,
                                                            thisProc->lowerRealBound(thisProd, maturity),
                                                            thisProc->upperRealBound(thisProd, maturity),
                                                            isCall ? isap->payoffToProcFreqBoundWeight : 1.0 - isap->payoffToProcFreqBoundWeight);

                    (*functions)[0] = Function1DComplexSP(new VanillaSingleIntegrandFFT<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                            (*thisProc,
                                                                             thisProd,
                                                                             maturity,
                                                                             omega));
                }
            }
            else {  // !useOneIntegral
                if (inst->fwdStarting) {
                    const FwdStFourierProductLogRtn& thisProd = *this;
                    const FwdStFourierProcessLogRtn* thisProc = dynamic_cast<const FwdStFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support FwdStFourierProcessLogRtn interface");
                    }

                    int iIntegral = 0;
                    for (; iIntegral < nbIntegrals; ++iIntegral){
                        (*functions)[iIntegral] = Function1DComplexSP(new VanillaDoubleIntegrandFFT<FwdStFourierProcessLogRtn, FwdStFourierProductLogRtn>
                                                                                (*thisProc,
                                                                                 thisProd,
                                                                                 maturity,
                                                                                 iIntegral));
                    }
                }
                else {  // started
                    const StFourierProductLogRtn& thisProd = *this;
                    const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support StFourierProcessLogRtn interface.");
                    }

                    int iIntegral = 0;
                    for (; iIntegral < nbIntegrals; ++iIntegral){
                        (*functions)[iIntegral] = Function1DComplexSP(new VanillaDoubleIntegrandFFT<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                                (*thisProc,
                                                                                 thisProd,
                                                                                 maturity,
                                                                                 iIntegral));
                    }
                }
            }
            return functions;
        }
        catch (exception& e){
            throw ModelException(e, method, "Failed to construct integrand(s)");
        }
    }

    /** Post process method for FFT integrator */
    virtual void postResults(const FourierEngine* model,
                             const FFTIntegrator1D*  integrator,
                             const FourierProductFFTIntegrator1D::IntegralArray& integrals,
                             CControl*            control,
                             CResults*            results) {
        static const string method = "VanillaFP::postResults";
        try {
            double price = 0.0;

            // calculate the discount factor back to today
            double pv = inst->instSettle->pv(today,
                                             paymentDate,
                                             inst->discount.get(),
                                             inst->asset.get());

            double logMoneyness = log(strike / fwd);
            if (useOneIntegral){
                price = exp(-omega * logMoneyness) * integrals[0]->getValue(logMoneyness) * strike;
            }
            else{   // !useOneIntegral
                price = (fwd * integrals[1]->getValue(logMoneyness) - strike * integrals[0]->getValue(logMoneyness)) / Maths::PI;
                if (inst->isCall) {
                    price += 0.5 * (fwd - strike);
                }
                else {  // put
                    price -= 0.5 * (fwd - strike);
                }
            }
            price *= mult * pv;
            results->storePrice(price, discount->getCcy());
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    const DateTime& getStartDate() const {return inst->startDate;}

private:
    const CVanilla* inst;
    DateTime        maturity;
    DateTime        paymentDate;
    double          strike;
    double          mult;
    bool            isCall;
    double          fwd;
    bool            useOneIntegral;
    double          omega;

    static const double callLowerBound;
    static const double callUpperBound;
    static const double putLowerBound;
    static const double putUpperBound;
};

const double VanillaFP::callLowerBound = 1.0;
const double VanillaFP::callUpperBound = 100.0;     // infinity
const double VanillaFP::putLowerBound = -100.0;     // infinity
const double VanillaFP::putUpperBound = 0.0;

/** Implementation of FourierEngine::IntoProduct interface */
FourierProduct* CVanilla::createProduct(const FourierEngine* model) const {
    static const string routine("Vanilla::createProduct");
    if (exerciseSchedule->length() > 1) {
        throw ModelException(routine, "Fourier vanilla model can only handle"
                             " single exercise dates.");
    }

    if (canExerciseEarly) {
        throw ModelException(routine, "Fourier vanilla model can not"
                             " handle early exercise feature.");
    }

    if(!Maths::isPositive(exerciseSchedule->lastValue())) {
        throw ModelException(routine, "Strike must be strictly positive.");
    }

    return new VanillaFP(this, exerciseSchedule->lastDate());
}

/** make a simple started vanilla ready for pricing */
CVanilla* CVanilla::make(const DateTime&             valueDate,
                         bool                        isCall,
                         bool                        american,
                         const Schedule*             exerciseSchedule,
                         const CAsset*               asset,
                         const YieldCurve*           discount,
                         const InstrumentSettlement* settle,
                         int                         noExerciseWindow) {
    static const string routine = "Vanilla::make";
    try {
        CVanillaSP vanilla(new CVanilla());

        vanilla->valueDate        = valueDate;
        vanilla->isCall           = isCall;
        vanilla->canExerciseEarly = american;
        vanilla->exerciseSchedule = ScheduleSP(copy(exerciseSchedule));

        vanilla->oneContract      = true;
        vanilla->notional         = 1.0;
        vanilla->initialSpot      = 1.0;
        vanilla->asset            = CAssetWrapper(copy(asset));
        vanilla->ccyTreatment     = CAsset::CCY_TREATMENT_NONE;
        vanilla->discount         = YieldCurveWrapper(copy(discount));
        vanilla->instSettle       = InstrumentSettlementSP(copy(settle));
        vanilla->noExerciseWindow = noExerciseWindow;

        return vanilla.release();
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** make a simple started vanilla ready for pricing from market wrappers */
CVanilla* CVanilla::make(const DateTime&             valueDate,
                         bool                        isCall,
                         bool                        american,
                         bool                        oneContract,
                         double                      notional,
                         double                      initialSpot,
                         const Schedule*             exerciseSchedule,
                         const CAssetWrapper&        asset,
                         const YieldCurveWrapper&    discount,
                         const InstrumentSettlement* settle,
                         int                         noExerciseWindow) {
    static const string routine = "Vanilla::make";
    try {
        CVanillaSP vanilla(new CVanilla());

        vanilla->valueDate        = valueDate;
        vanilla->isCall           = isCall;
        vanilla->canExerciseEarly = american;
        vanilla->exerciseSchedule = ScheduleSP(copy(exerciseSchedule));

        vanilla->oneContract      = oneContract;
        vanilla->notional         = notional;
        vanilla->initialSpot      = initialSpot;
        vanilla->ccyTreatment     = CAsset::CCY_TREATMENT_NONE;
        vanilla->instSettle       = InstrumentSettlementSP(copy(settle));
        vanilla->noExerciseWindow = noExerciseWindow;

        // Deep copy of the wrapper objects
        smartPtr<CAssetWrapper> assetCopy(copy(&asset));
        vanilla->asset = *assetCopy;

        smartPtr<YieldCurveWrapper> discountCopy(copy(&discount));
        vanilla->discount = *discountCopy;

        return vanilla.release();
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** make a vanilla ready for pricing */
CVanilla* CVanilla::make(const DateTime&             valueDate,
                         bool                        isCall,
                         bool                        american,
                         const Schedule*             exerciseSchedule,
                         bool                        fwdStarting,
                         const DateTime&             startDate,
                         bool                        oneContract,
                         double                      notional,
                         double                      initialSpot,
                         const CAsset*               asset,
                         const string&               ccyTreatment,
                         const YieldCurve*           discount,
                         const InstrumentSettlement* settle) {
    static const string routine = "Vanilla::make";
    try {
        CVanillaSP vanilla(new CVanilla());

        vanilla->valueDate        = valueDate;
        vanilla->isCall           = isCall;
        vanilla->canExerciseEarly = american;
        vanilla->exerciseSchedule = ScheduleSP(copy(exerciseSchedule));

        vanilla->fwdStarting      = fwdStarting;
        vanilla->startDate        = startDate;

        vanilla->oneContract      = oneContract;
        vanilla->notional         = notional;
        vanilla->initialSpot      = initialSpot;
        vanilla->asset            = CAssetWrapper(copy(asset));
        vanilla->ccyTreatment     = CAsset::CCY_TREATMENT_NONE;
        vanilla->discount         = YieldCurveWrapper(copy(discount));
        vanilla->instSettle       = InstrumentSettlementSP(copy(settle));

        return vanilla.release();
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

class CVanillaHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CVanilla, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ISupportVegaMatrixLite);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(FourierEngine::IIntoProduct);
        IMPLEMENTS(FD1F::IIntoProduct);
        IMPLEMENTS(FD1FGeneric::IIntoProduct);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultVanilla);
        FIELD(valueDate,        "valuation Date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(startDate,        "Option start date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(premiumSettle, "Settlement of option premium");
        FIELD_MAKE_OPTIONAL(premiumSettle);
        FIELD(isCall,           "Is it a call option");
        FIELD(fwdStarting,      "Is it a fwd starting option");
        FIELD(exerciseSchedule,        "Exercise Schedule");
        FIELD(canExerciseEarly, "Can option be exercised early");
        FIELD(oneContract,      "Calc price for 1 contract");
        FIELD(notional,         "Option notional");
        FIELD_MAKE_OPTIONAL(notional);
        FIELD(initialSpot,      "Initial spot price");
        FIELD_MAKE_OPTIONAL(initialSpot);
        FIELD(spotAtMaturity,   "underlying spot level when exercised");
        FIELD_MAKE_OPTIONAL(spotAtMaturity);
        FIELD(asset,            "Underlying of option");
        FIELD(discount,         "Discount curve");
        FIELD(ccyTreatment,     "Currency Treatment");
        FIELD(instSettle, "Instrument settlement at maturity");
        FIELD(isExercised, "Indicates whether option has been exercised");
        FIELD_MAKE_OPTIONAL(isExercised);
        FIELD(dateExercised,        "Date on which option has been exercised");
        FIELD_MAKE_OPTIONAL(dateExercised);
        FIELD(noExerciseWindow, "Number of buisness days prior to ex-dividend date for "
                                       "which early exercise is disallowed");
        FIELD_MAKE_OPTIONAL(noExerciseWindow);
    }

    static IObject* defaultVanilla(){
        return new CVanilla();
    }
};

CClassConstSP const CVanilla::TYPE = CClass::registerClassLoadMethod(
    "Vanilla", typeid(CVanilla), CVanillaHelper::load);
bool  CVanillaLoad() {
    return (CVanilla::TYPE != 0);
   }

/** Class that creates a DeltaImpledStrike calculator */
CVanilla::DeltaImpliedStrikeMakerLN::DeltaImpliedStrikeMakerLN():
CObject(TYPE) {}

IObject* CVanilla::DeltaImpliedStrikeMakerLN::defaultDeltaImpliedStrikeMakerLN() {
    return new DeltaImpliedStrikeMakerLN();
}

void CVanilla::DeltaImpliedStrikeMakerLN::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CVanilla::DeltaImpliedStrikeMakerLN, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IDeltaToStrikeMaker);
    EMPTY_SHELL_METHOD(defaultDeltaImpliedStrikeMakerLN);
}

const IDeltaToStrikeMaker::IDeltaToStrike* CVanilla::DeltaImpliedStrikeMakerLN::make(
    const DateTime&             valueDate,
    const DateTime&             maturityDate,
    bool                        isCall,
    const CAsset*               asset,
    const YieldCurve*           discount,
    const InstrumentSettlement* settle,
    double                      deltaShiftSize,
    double                      tgtDelta,
    const string&               volType,
    bool                        allowNegativeFwdVar) {

    model = IModelSP(new CClosedFormLN(volType,allowNegativeFwdVar));

    IDeltaToStrikeMaker::IDeltaToStrike* ptr = new CVanilla::DeltaImpliedStrike(
        *model,
        valueDate,
        maturityDate,
        isCall,
        asset,
        discount,
        settle,
        deltaShiftSize,
        tgtDelta);

    return ptr;
}

CClassConstSP const CVanilla::DeltaImpliedStrikeMakerLN::TYPE =
CClass::registerClassLoadMethod("CVanilla::DeltaImpliedStrikeMakerLN", typeid(CVanilla::DeltaImpliedStrikeMakerLN), load);

/** Implementation of Vanilla::DeltaImpliedStrike hidden here */
class CVanilla::DeltaImpliedStrike::Imp{
    friend class CVanilla::DeltaImpliedStrike;

private:
    Imp(const IModel&                model,
        const DateTime&             valueDate,
        const DateTime&             maturityDate,
        bool                        isCall,
        const CAsset*               asset,
        const YieldCurve*           discount,
        const InstrumentSettlement* settle,
        double                      deltaShiftSize,
        double                      tgtDelta):
    model(copy(&model)),
    maturityDate(maturityDate),
    tgtDelta(tgtDelta),
    isCall(isCall) {
        double dummyStrike = asset->getSpot();
        Schedule exoSchedule(DateTimeArray(1, maturityDate),
                             DoubleArray(1, dummyStrike),
                             "N");      // no interpolation
        vanilla = CVanillaSP(
            CVanilla::make(valueDate,
                           isCall,
                           false,   // european
                           &exoSchedule,
                           asset,
                           discount,
                           settle,
                           0));     // noExerciseWindow (not used)
        SensitivityArraySP sens(new SensitivityArray(1));
        delta = DeltaSP(new Delta(deltaShiftSize));
        (*sens)[0] = SensitivitySP(new Delta(deltaShiftSize));
        OutputRequestArraySP outReqs(new OutputRequestArray(0));
        control = CControlSP(new Control(sens,
                                         outReqs,
                                         false,
                                         ""));
    }

    static double resetVanillaStrike(CVanilla* vanilla,
                                     double strike){
        double oldStrike = vanilla->exerciseSchedule->lastValue();
        vanilla->exerciseSchedule = ScheduleSP(
            new Schedule(*vanilla->exerciseSchedule,
                         strike));
        return oldStrike;
    }

    // calc delta
    double calcDelta(double logStrike) const{
        static const string method = "Vanilla::DeltaImpliedStrike::Imp::calcDelta";
        try {
            double strike = exp(logStrike);
            // double oldStrike = resetVanillaStrike(vanilla.get(), strike);
            resetVanillaStrike(vanilla.get(), strike);
            CResultsSP results(model->Run(vanilla.get(), control.get()));
            // resetVanillaStrike(vanilla.get(), oldStrike);
            OutputNameArraySP names(results->packetContents(Delta::NAME));
            ASSERT(names->size() == 1);
            IObjectConstSP deltaobj(
                results->retrieveGreek(Delta::NAME, (*names)[0]));

            if(Untweakable::TYPE->isInstance(deltaobj)) {
                const Untweakable* tmp = dynamic_cast<const Untweakable*>(deltaobj.get());
                throw ModelException(method, tmp->getMessage());
            }

            CDoubleConstSP deltaval(CDoubleConstSP::dynamicCast(deltaobj));
            return deltaval->doubleValue();
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // needed by bracketer + solver
    typedef const Imp* Ptr;
    typedef double (Imp::* Func1DConstPtr)(double) const;
    typedef MemFuncWrapper<Ptr, Func1DConstPtr> DeltaFunc;
    typedef FuncShiftWrapper<DeltaFunc> DeltaDiffFunc;

    // calc strike implied by delta
    double calcStrike(double& lowerStrike,
                      double& upperStrike,
                      double  strikeAbsAcc) const{
        static const string method = "Vanilla::DeltaImpliedStrike::Imp::calcStrike";
        try {
            // validate 0 <= lowerStrike < upperStrike
            if (Maths::isNegative(lowerStrike)
                || !Maths::isPositive(upperStrike - lowerStrike)){
                throw ModelException(method, "we should have 0.0 <= lowerStrike < upperStrike; got "
                                              + Format::toString(lowerStrike)
                                              + " and "
                                              + Format::toString(upperStrike)
                                              + ", respectively");
            }

            if(Maths::isZero(tgtDelta) && !isCall) {
                return 0.0;
            }

            // bracket the root
            double lowerLogStrike = log(lowerStrike);
            double upperLogStrike = log(upperStrike);
            DeltaFunc deltaFunc(this, &Imp::calcDelta);
            DeltaDiffFunc deltaDiffFunc(deltaFunc, -tgtDelta);
            try{
                ZBrac_bracket(deltaDiffFunc,
                              lowerLogStrike,
                              upperLogStrike);
            }
            catch (exception& e) {
                throw ModelException::addTextToException(e,
                                                         "Failed to bracket root at maturity "
                                                         + maturityDate.toString()
                                                         + " and delta level "
                                                         + Format::toString(tgtDelta));
            }
            lowerStrike = exp(lowerLogStrike);
            upperStrike = exp(upperLogStrike);
            // find the root using ZBrent
            double avgStrike = 0.5 * (lowerStrike + upperStrike);
            double logStrikeAbsAcc = strikeAbsAcc / avgStrike;
            double logStrike;
            try{
                logStrike = ZBrent_solve(deltaDiffFunc,
                                         lowerLogStrike,
                                         upperLogStrike,
                                         logStrikeAbsAcc);
            }
            catch (exception& e) {
                throw ModelException::addTextToException(e,
                                                         "Failed to solve for root at maturity "
                                                         + maturityDate.toString()
                                                         + " and delta level "
                                                         + Format::toString(tgtDelta));
            }
            return exp(logStrike);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // fields
    IModelSP           model;
    DateTime           maturityDate;
    double             tgtDelta;
    double             isCall;
    CVanillaSP         vanilla;
    DeltaSP            delta;
    CControlSP         control;
};

CVanilla::DeltaImpliedStrike::DeltaImpliedStrike(
        const IModel&                model,
        const DateTime&             valueDate,
        const DateTime&             maturityDate,
        bool                        isCall,
        const CAsset*               asset,
        const YieldCurve*           discount,
        const InstrumentSettlement* settle,
        double                      deltaShiftSize,
        double                      tgtDelta):
me(new Imp(model,
           valueDate,
           maturityDate,
           isCall,
           asset,
           discount,
           settle,
           deltaShiftSize,
           tgtDelta)) { }

double CVanilla::DeltaImpliedStrike::calcStrike(
        double             lowerStrike,
        double             upperStrike,
        double             strikeAbsAcc) const{
    return me->calcStrike(lowerStrike,
                          upperStrike,
                          strikeAbsAcc);
}

CVanilla::DeltaImpliedStrike::~DeltaImpliedStrike(){}

void CVanilla::DeltaImpliedStrike::Helper::getMarket(
        const MarketData*        market,
        const IModel*             model,
        CAssetWrapper&           asset,
        YieldCurveWrapper&       discount,
        InstrumentSettlement*    settle){
    CAsset::getAssetMarketData(model,
                               market,
                               CAsset::CCY_TREATMENT_NONE,
                               discount,
                               asset);
    discount.getData(model, market);
    settle->getMarket(model, market);
}


/** Returns the name of the instrument's discount currency. */
string CVanilla::discountYieldCurveName() const {
    return discount.getName();
}

/*********************************************************************/

/* The ultimate wrapping of FourierEngine::ISAP, mainly for use in Pyramid
 */
#define ISAP_TYPE_NONE              "none"
#define ISAP_TYPE_VANILLA           "VanillaISAP"
#define ISAP_TYPE_VANILLAGRID       "VanillaGridISAP"

class ISAPWrapper : public CObject,
                    virtual public ITypeConvert {
public: // how can I have this protected or private?
    string                  isapType;     // Vanilla or VanillaGrid
    EmptyISAPSP             emptyISAP;
    VanillaISAPSP           vanillaISAP;
    VanillaGridISAPSP       vanillaGridISAP;

private:
    FourierEngine::ISAPSP   realISAP;

public:
    static CClassConstSP const TYPE;

    // validation
    void validatePop2Object(){
        static const string routine = "ISAPWrapper::validatePop2Object";
        try{
            if (isapType.empty()){
                throw ModelException(routine,
                                     "Blank FourierEngine ISAP specified!");
            }
            if (CString::equalsIgnoreCase(isapType,ISAP_TYPE_NONE)) {
//                realISAP = FourierEngine::ISAPSP(); // null ISAP
                if (emptyISAP.get()) {
                    realISAP = emptyISAP;
                } else {
                    throw ModelException(routine, "Expected EmptyISAP "
                                         "but none supplied!");
                }
            } else if (CString::equalsIgnoreCase(isapType,ISAP_TYPE_VANILLA)) {
                if (vanillaISAP.get()) {
                    realISAP = vanillaISAP;
                } else {
                    throw ModelException(routine, "Expected VanillaISAP "
                                         "but none supplied!");
                }
            } else if (CString::equalsIgnoreCase(isapType,ISAP_TYPE_VANILLAGRID)) {
                if (vanillaGridISAP.get()) {
                    realISAP = vanillaGridISAP;
                } else {
                    throw ModelException(routine, "Expected VanillaGridISAP "
                                         "but none supplied!");
                }
            } else {
                throw ModelException(routine, "Unrecognised Intergrator1D "
                                     + isapType + ". Expected "
                                     + ISAP_TYPE_VANILLA + " or "
                                     + ISAP_TYPE_VANILLAGRID);
            }
        }  catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** create a proper ISAP */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const {
        static const string method = "ISAPWrapper::convert";
        try {
            if (requiredType != FourierEngine::ISAP::TYPE) {
                throw ModelException(method,
                                     "Cannot convert a ISAPWrapper into "
                                     "object of type "+requiredType->getName());
            }
            object = realISAP;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ISAPWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ITypeConvert);
        EMPTY_SHELL_METHOD(defaultISAPWrapper);
        FIELD(isapType, "VanillaISAP or VanillaGridISAP");
        FIELD(emptyISAP,  "EmptyISAP");
        FIELD_MAKE_OPTIONAL(emptyISAP);
        FIELD(vanillaISAP,  "VanillaISAP");
        FIELD_MAKE_OPTIONAL(vanillaISAP);
        FIELD(vanillaGridISAP,  "VanillaGridISAP");
        FIELD_MAKE_OPTIONAL(vanillaGridISAP);
        FIELD(realISAP, "real ISAP");
        FIELD_MAKE_TRANSIENT(realISAP);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    // for reflection
    ISAPWrapper(): CObject(TYPE){}

    static IObject* defaultISAPWrapper(){
        return new ISAPWrapper();
    }
};

typedef smartPtr<ISAPWrapper> ISAPWrapperSP;

CClassConstSP const ISAPWrapper::TYPE =
CClass::registerClassLoadMethod("ISAPWrapper",
                                typeid(ISAPWrapper), load);


DRLIB_END_NAMESPACE

