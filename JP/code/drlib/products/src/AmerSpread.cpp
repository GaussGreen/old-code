//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AmerSpread.cpp
//
//   Description : American spread vanilla option
//
//   Author      : Ning shen
//
//   Date        : 06 Jan 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Black.hpp"
#include "edginc/AmerSpread.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/XCB.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VegaParallel.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

void CAmerSpread::Validate()
{
    static const string method = "CAmerSpread::Validate";
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
    
    if ( FXAsset::TYPE->isInstance(asset.get()) )
    {
        throw ModelException(method,
                             "Options on FX assets are not allowed yet");
    }
    
    if (noExerciseWindow < 0) {
        throw ModelException(method,
                             "noExerciseWindow must not be negative");
    }
}

void CAmerSpread::validatePop2Object()
{
    static const string method("CAmerSpread::validatePop2Object");
    int                 numDates;

    // can't get exercise schedule from Market - fail if it is NULL
    if( !exerciseScheduleHi || !exerciseScheduleLo)
    {
        throw ModelException(method, "either hi or lo exercise schedule is NULL");
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
    numDates = exerciseScheduleHi->length();
    if ( numDates < 1 || numDates != exerciseScheduleLo->length())
    {
        throw ModelException(method, "one of exercise schedules is empty or they are not equal in size");
    }

    for (int i=0; i<numDates; i++)
    {
        if (exerciseScheduleHi->getDates()[i] != exerciseScheduleLo->getDates()[i])
            throw ModelException(method, "exercise date in Hi schedule must match that in Lo schedule");
        if (exerciseScheduleHi->getValues()[i] < exerciseScheduleLo->getValues()[i])
            throw ModelException(method, "hi strike cannot be less than lo strike");
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

DateTime CAmerSpread::getValueDate() const
{
  return valueDate;
}

/** when to stop tweaking */
DateTime CAmerSpread::endDate(const Sensitivity* sensControl) const {
    DateTime matDate = exerciseScheduleHi->lastDate();
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

void CAmerSpread::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime matDate = exerciseScheduleHi->lastDate();
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
             results->storeNotApplicable(request); // LV does not return this
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
bool CAmerSpread::priceDeadInstrument(CControl* control, CResults* results) const
{
    double    strikeHi        = 0.0;
    double    strikeLo        = 0.0;
    double    value         = 0.0;
    bool      foundExerDate = false;
    DateTime  exerDate;

    static string method = "CAmerSpread::priceDeadInstrument";

    DateTime matDate = exerciseScheduleHi->lastDate();

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
    {// maturity instrinsic value
        exerDate      = matDate;
        strikeHi        = exerciseScheduleHi->lastValue();
        strikeLo        = exerciseScheduleLo->lastValue();
        foundExerDate = true;
    }
    else
    {// may early exercise
        // check that exercise date is not in the Future
        if ( dateExercised.isGreater(valueDate))
        {
            throw ModelException(method,
                    "Option exercised on " + 
                    dateExercised.toString() + ". " +
                    "This date is after the current value date (" + 
                    valueDate.toString());
        }

        // check whether exercise date is valid and find corresponding strike
        if ( !canExerciseEarly ) {
            // multi-european case
            DateTimeArray exerciseDates = exerciseScheduleHi->getDates();

            for (int i=0 ; i<exerciseDates.size() ; ++i)
            {
                if ( dateExercised.equals(exerciseDates[i])) {
                    foundExerDate = true;
                    exerDate      = dateExercised;
                    strikeHi        = exerciseScheduleHi->interpolate(
                                                    dateExercised);
                    strikeLo        = exerciseScheduleLo->interpolate(
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
                        exerciseScheduleHi->firstDate()) &&
                 !dateExercised.isGreater(
                        exerciseScheduleHi->lastDate())  ) {
                // american case - note: not checking for weekend yet
                foundExerDate = true;
                exerDate      = dateExercised;
                strikeHi        = exerciseScheduleHi->interpolate(
                                                dateExercised);
                strikeLo        = exerciseScheduleLo->interpolate(
                                                dateExercised);
            }
            else 
            {
                throw ModelException(method, 
                        "Option has been exercised on " + 
                        dateExercised.toString() + " but valid " +
                        "exercise range is from " + 
                        exerciseScheduleHi->firstDate().toString() +
                        " to " +
                        exerciseScheduleHi->lastDate().toString());
            }
        }
    }

    if ( foundExerDate )
    {
        double valueHi = GetIntrinsic(spotAtMaturity,
                                 strikeHi,
                                 isCall, 
                                 true /* isOption */);
        double valueLo = GetIntrinsic(spotAtMaturity,
                                 strikeLo,
                                 isCall, 
                                 true /* isOption */);

        settlementDate = instSettle->settles(exerDate, asset.get());
        // pv from settlement to today
        value = fabs(valueHi - valueLo)*discount->pv(valueDate, settlementDate);
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

class AmerSpreadFDProd: public LatticeProdEDRIns
{
public:
    AmerSpreadFDProd( const CAmerSpread * amsp, FDModel * mdl ) :
        LatticeProdEDRIns( mdl, 1, 1 ),
        inst( amsp )
    {
        if( ! tree1f )
        {
            throw ModelException( "AmerSpreadFDProd::AmerSpreadFDProd", "Instrument of type "+
                                 inst->getClass()->getName() +
                                 " can be priced by CTree1f only" );
        }

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );
    }

    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const
    {
        return inst->ccyTreatment;
    }

    /** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

    /** this sets up the timeline */
    virtual void init( CControl * control ) const
    {
        static const string method = "AmerSpreadFDProd::Init";
        try {
            // default to NODE_INSERTION smoothing
            if (tree1f->GetSmoothMethod() == CTree1f::DEFAULT) {
                tree1f->SetSmoothMethod(CTree1f::NODE_INSERTION);
            }
            tree1f->NumOfPrice = numPrices;
            tree1f->NumOfInsertNode =
                ( tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ? 1 : 0 );

            if( inst->fwdStarting )
                tree1f->controlSameGridFwdStart(inst->ccyTreatment);            
            
            // all exercise dates are copied to critical dates
            DateTimeArray critDates = inst->exerciseScheduleHi->getDates();
            // remove exercise date from crit date
            critDates.erase(critDates.end()-1);

            // don't bother with dollar divs
            tree1f->SetDivAmountTreatment(false);

            // add div event dates if needed
            EventAssetMove divEvent;
            DateTimeArraySP divCritDates;

            // define start and end points of tree
            DateTimeArray segDates(2);
            segDates[0] = getStartDate();
            segDates[1] = inst->exerciseScheduleHi->lastDate();

            if( inst->canExerciseEarly )
            {
                int numDivs = (int)(4*segDates[0].yearFrac(segDates[1]))+1; // 4 divs per year selected as critical dates
                if (AssetUtil::getJumpEvents(inst->asset.get(), 
                                             segDates[0], 
                                             segDates[1],
                                             numDivs,
                                             divEvent)) {
                    // calculate critical dates
                    divCritDates = divEvent.getCritDate(inst->noExerciseWindow, inst->isCall);
                }
            }

            // use simple delta size adjustment if needed, note that TreeDeltaShift is stored in base tree only.
            if (!tree1f->DEBUG_SameGridDelta && control->isPricing()) {
                double deltaShift = control->getDeltaShiftSize();

                if (Maths::isPositive(deltaShift)){
                    DividendListSP dollarDivs = 
                        AssetUtil::getDollarDivsBetweenDates(inst->asset.get(),
                                                             inst->valueDate,
                                                             inst->exerciseScheduleHi->lastDate());
                    tree1f->DeltaSizeAdjust(deltaShift,
                                            inst->valueDate,
                                            inst->exerciseScheduleHi->firstDate(),
                                            inst->exerciseScheduleHi->lastDate(), 
                                            inst->exerciseScheduleHi->length()>1 && !inst->canExerciseEarly,
                                            dollarDivs->getDivAmounts()->size()>0,
                                            inst->fwdStarting);
                }
            }

            // add critical dates
            model->addCritDates( critDates );
            if( divCritDates.get() )
                tree1f->addDivCritDates( *divCritDates );

            // timeline configuration
            // 'density factor' for timeline
            IntArray density( 1, 1 );

            // prepare timeline set up
            model->initSegments( segDates, density );
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** initialising and setting product variables */
    // this is called per pricing call before tree sweep call (after InitTree)
    virtual void initProd()
    {
        static const string method = "AmerSpreadFDProd::InitProd";
        try {
            int i;
            // timeline exists at this point
            int lastStep = model->getLastStep();

            initSlices( numPrices );
            initInsertNode();

            CAssetConstSP asset = inst->asset.getSP();
            double fwdAtStart = 1.0;

            stepStrikeHi.resize(lastStep+1);
            stepStrikeLo.resize(lastStep+1);
            stepCanExercise.resize(lastStep+1);
            // ask tree to decide first about steps that can exercise
            AssetUtil::setStepExercise(stepCanExercise,
                                    model->getDates(),
                                    inst->exerciseScheduleHi,
                                    inst->canExerciseEarly,
                                    asset);

            /* If no exercise before ex-div date applies */
            if (inst->canExerciseEarly && inst->noExerciseWindow > 0) {
                /* Get dividend list */
                DateTime lastDate = AssetUtil::getHoliday(asset.get())->addBusinessDays(
                    model->getDate(lastStep), 
                    inst->noExerciseWindow);

                DividendListSP divList(AssetUtil::getAllDivsBetweenDates(inst->asset,
                                                                         model->getDate(0),
                                                                         lastDate));

                /* If dividend list not empty */
                if (divList->getArray().size() > 0){
                    HolidayConstSP holiday(
                        AssetUtil::getHoliday(inst->asset.get()));
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
                        if (stepCanExercise[iStep] == true && daysToNextDiv <= inst->noExerciseWindow) {
                            stepCanExercise[iStep] = false;
                        }
                    }
                }
            }

            bool canInterpExSched = (inst->exerciseScheduleHi->length() > 1);
            // get spot at start if needed
            if (inst->fwdStarting) {
                fwdAtStart = inst->asset->fwdValue(inst->startDate);
            }

            CIntArray  canEx;
            canEx.resize(lastStep+1); // for debug
            canEx[lastStep] = 1;

            // get strikes
            stepStrikeHi[lastStep] = fwdAtStart*inst->exerciseScheduleHi->lastValue();
            stepStrikeLo[lastStep] = fwdAtStart*inst->exerciseScheduleLo->lastValue();

            for (i=0; i<= lastStep; i++) {
                if (stepCanExercise[i]) {
                    canEx[i] = 1;
                    if (canInterpExSched) {
                        stepStrikeHi[i] = fwdAtStart*inst->exerciseScheduleHi->interpolate(model->getDate(i));
                        stepStrikeLo[i] = fwdAtStart*inst->exerciseScheduleLo->interpolate(model->getDate(i));
                    }
                    else {
                        stepStrikeHi[i] = stepStrikeHi[lastStep]; // strike level is constant for american with no schedule
                        stepStrikeLo[i] = stepStrikeLo[lastStep]; // strike level is constant for american with no schedule
                    }
                }
                else {
                    canEx[i] = 0;
                }
            }

#ifdef TREE_DEBUG
            tree1f->stepCanExercise = CIntArraySP(dynamic_cast<CIntArray*>(canEx.clone()));
#endif
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // here just take care of scaling by notional and any additional 
    // discounting for fwd starting case
    double scalePremium(const double& fairValue, YieldCurveConstSP disc)
    {
        double fwdAtStart = 0.0;
        if (inst->fwdStarting)
        {
            fwdAtStart = inst->asset->fwdValue(inst->startDate);
        }

        double scalingFactor = InstrumentUtil::scalePremium(
                                        inst->oneContract,
                                        inst->fwdStarting,
                                        inst->notional,
                                        fwdAtStart,
                                        inst->initialSpot);

        return (fairValue*scalingFactor);
    }

    /** extra output requests */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results)
    {
        // save price
        double price = scalePremium(model->getPrice0( *slices[0] ), disc);
        results->storePrice(price, disc->getCcy());

        // take care of additional outputs
        if ( control && control->isPricing() )
        {
            double indVol = 0.0; // LV cannot return indicative vol

            inst->addOutputRequests(control,
                                    results,
                                    price,
                                    indVol);
        }
    }
    
    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type)
    {
        // we assume just need one und level for spot here
        const TreeSlice & s = payoffIndex->getValue( step );
        int bot, top;
        s.getCalcRange( bot, top );

        const vector< TreeSliceSP > & price = slices;
        int pStart = 0, pEnd = price.size() - 1;

        if (type == FDProduct::BWD_T){
            prod_BWD_T(   s,
                          step,
                          bot,
                          top,
                          pStart, 
                          pEnd,
                          price);

            //insert nodes
            if (tree1f && tree1f->NumOfInsertNode>0){
                prod_BWD_T(   *insNodes,
                              step,
                              0,
                              tree1f->NumOfInsertNode-1,
                              pStart, 
                              pEnd,
                              *insPrices);
            }

        }
        else if(type == FDProduct::BWD){
            prod_BWD( s,
                      step,
                      bot,
                      top,
                      pStart, 
                      pEnd,
                      price);

            //insert nodes
            if (tree1f && tree1f->NumOfInsertNode>0)
            {
                prod_BWD(   *insNodes,
                            step,
                            0,
                            tree1f->NumOfInsertNode-1,
                            pStart, 
                            pEnd,
                            *insPrices);
            }
        }
    }

    /** product payoff method at maturity */
    void prod_BWD_T(
        const TreeSlice & spot,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price)
    {
        static const string method("AmerSpreadFDProd::prod_BWD_T");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int i, j;
            double v1, v2;
    
            double settlementPV = inst->instSettle->pvAdjust(inst->exerciseScheduleHi->lastDate(),
                                                             inst->discount.get(), 
                                                             inst->asset.get());
            for (i=pStart; i<=pEnd; i++)
            {
                for (j=bot; j<=top; j++)
                {
                    v1 = GetIntrinsic(s[j], stepStrikeHi[step], inst->isCall, true/* this allow fwd */);
                    v2 = GetIntrinsic(s[j], stepStrikeLo[step], inst->isCall, true/* this allow fwd */);
                    (p[i])[j] = fabs(v1 - v2)*settlementPV;
                }
            }
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    // this is payoff boundary condition, for KO, early exercise etc.
    void prod_BWD(
        const TreeSlice & spot,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price)
    {
        static const string method("AmerSpreadFDProd::prod_BWD");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int j;
            double v1, v2;
   
            if (tree1f && step == model->getLastStep()-1)
            {// do penultimate smoothing
                vector <double> vol_arr, drift_arr;
                double dt = tree1f->getTrdYrFrac(step+1);
                bool fastRoll = tree1f->CalcStepDriftAndVar(s, bot, top, vol_arr, &drift_arr); // just to get drift
                // use GetStepVol() to get vol, do not use GetStepVar() which may not have simple BS conversion
                tree1f->GetStepVol(step, vol_arr, s, bot, top);

                double variance = vol_arr[0]*vol_arr[0]*dt;
                double drift = drift_arr[0];
                double df = inst->discount->pv(model->getDate(step), model->getDate(step+1));

                /* If dollar dividend treatment, we must adjust fwd and strike by stock floor.
                   If not dollar dividend treatment, stock floor == 0.0 */
                double StockFloor0 = 0.0;
                double StockFloor1 = 0.0;
     
                for (j=bot; j<=top; j++)
                {
                    if (!fastRoll)
                    {
                        if (vol_arr.size() > 1){
                            variance = vol_arr[j-bot]*vol_arr[j-bot]*dt; // needs one vol per node - local vol
                        }
                        if (drift_arr.size() > 1){
                            drift = drift_arr[j-bot]; // needs drift per node
                        }
                    }

                    if (s[j]>0.0)
                    {
                        if (fabs(log(s[j]/stepStrikeHi[step+1])) < tree1f->TruncationStd*sqrt(variance)
                            ||fabs(log(s[j]/stepStrikeLo[step+1])) > tree1f->TruncationStd*sqrt(variance))
                        {
                            /* If dollar div treatment (in which case stock floor is <> 0.0), the price of the option
                               is that of an option with same maturity written on the pseudo asset (i.e., S - StockFloor)
                               and with strike adjusted by the stock floor at maturity (i.e., K - StockFloor).
                               NB Black::price takes care of potentially negative fwd and strike values. */
                            v1 = Black::price(inst->isCall, 
                                                 drift * (s[j] - StockFloor0),  // pseudo asset's (step+1)-maturity fwd
                                                 // as viewed from t = step
                                                 stepStrikeHi[step+1] - StockFloor1,  // adjusted strike
                                                 df,
                                                 variance);
                            v2 = Black::price(inst->isCall, 
                                                 drift * (s[j] - StockFloor0),  // pseudo asset's (step+1)-maturity fwd
                                                 // as viewed from t = step
                                                 stepStrikeLo[step+1] - StockFloor1,  // adjusted strike
                                                 df,
                                                 variance);
                            (p[0])[j] = fabs(v1 - v2);
                        }
                    }
                }
            }

            if (!stepCanExercise[step]) {
                return;
            }

            double settlementPV = inst->instSettle->pvAdjust(model->getDate(step),
                                                                 inst->discount.get(), 
                                                                 inst->asset.get());
   
            double callput = settlementPV*(inst->isCall? 1.0:-1.0);
            // for node insertion use
            int jIns = bot;
            bool needInsert = useInsertNode;

            double intrinsic, priceLast = 0.0, priceCurr = 0.0;

            for (j=bot; j<=top; j++)
            {
                v1 = Maths::max(0.0, callput*(s[j] - stepStrikeHi[step]));
                v2 = Maths::max(0.0, callput*(s[j] - stepStrikeLo[step]));
                intrinsic = callput*(v2 - v1);
                if (needInsert)
                {
                    priceLast = priceCurr;
                    priceCurr = (p[0])[j]; // keep it
                }
                if ((p[0])[j] < intrinsic)
                {
                    (p[0])[j] = intrinsic; // American
                    if (needInsert && inst->isCall) // use node insertion by simple linear interpolation
                    {
                        jIns = j;
                        needInsert = false; // done for call
                    }
                }
                else if (needInsert && !inst->isCall)
                {
                    jIns = j;
                    needInsert = false; // done for put
                }
            }
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }
    
private:
    const CAmerSpread* inst;
    bool              useInsertNode;
    vector<double>    stepStrikeHi;
    vector<double>    stepStrikeLo;
    vector<bool>      stepCanExercise;
};

/** create a fd payoff tree - new for all fd/tree state variable interface */
FDProductSP CAmerSpread::createProduct(FDModel* model) const
{
    return FDProductSP( new AmerSpreadFDProd(this, model) );
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool CAmerSpread::avoidVegaMatrix(const IModel*)
{
    return true;
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool CAmerSpread::sensShift(Theta* shift)
{    
    DateTime newDate = shift->rollDate(valueDate);

    DateTime matDate = exerciseScheduleHi->lastDate();

    if ( ( newDate >= matDate && valueDate < matDate ) ||
         ( valueDate == matDate && Maths::isZero(spotAtMaturity)))
        spotAtMaturity = asset->getThetaSpotOnDate(shift, matDate);

    if (fwdStarting && newDate.isGreaterOrEqual(startDate) &&
        startDate.isGreaterOrEqual(valueDate))
    {
        fwdStarting = false;
        initialSpot = asset->getThetaSpotOnDate(shift, startDate);
        exerciseScheduleHi->scale(initialSpot);
        exerciseScheduleLo->scale(initialSpot);
    }

    // roll today 
    valueDate = newDate;
    
    return true;
};

// for reflection
CAmerSpread::CAmerSpread(): 
Generic1Factor(TYPE),isCall(false),
canExerciseEarly(false),
spotAtMaturity(0.0),
isExercised(false),
noExerciseWindow(0){}

class CAmerSpreadHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CAmerSpread, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultAmerSpread);
        FIELD(isCall,           "Is it a call option");
        FIELD(exerciseScheduleHi,        "High strike exercise Schedule");
        FIELD(exerciseScheduleLo,        "Lo strike exercise Schedule");
        FIELD(canExerciseEarly, "Can option be exercised early");
        FIELD(spotAtMaturity,   "underlying spot level when exercised");
        FIELD_MAKE_OPTIONAL(spotAtMaturity);
        FIELD(isExercised, "Indicates whether option has been exercised");
        FIELD_MAKE_OPTIONAL(isExercised);
        FIELD(dateExercised,        "Date on which option has been exercised");
        FIELD_MAKE_OPTIONAL(dateExercised);
        FIELD(noExerciseWindow, "Number of buisness days prior to ex-dividend date for "
                                       "which early exercise is disallowed");
        FIELD_MAKE_OPTIONAL(noExerciseWindow);
    }

    static IObject* defaultAmerSpread(){
        return new CAmerSpread();
    }
};

CClassConstSP const CAmerSpread::TYPE = CClass::registerClassLoadMethod(
    "AmerSpread", typeid(CAmerSpread), CAmerSpreadHelper::load);

bool AmerSpreadLoad()
{
    return (true && CAmerSpread::TYPE);
}

DRLIB_END_NAMESPACE

