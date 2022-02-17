//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Installment.cpp
//
//   Author      : Ning Shen
//
//   Description   Installment option
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/FD1F.hpp"
#include "edginc/Maths.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Black.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/XCB.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Tree1fLN.hpp"
#include "edginc/VegaParallel.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////// class declaration
class CInstallment : public Generic1Factor, 
                     public FDModel::IIntoProduct,
                     public LastSensDate
{
public:
    static CClassConstSP const TYPE;

    // override base implementation if required
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);
    
    virtual void Validate();
 
    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    // below are copied from Vanilla
    virtual DateTime getValueDate() const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;
    
    bool sensShift(Theta* shift);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel *     model);

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    CashFlowArraySP getKnownCashFlows() const;

    DateTimeArraySP getPaymentDates() const;
 
private:
    friend class CInstallmentHelper;
    friend class CInstallmentFDProd;

protected:
    static void load(CClassSP& clazz);
    CInstallment();

    // this block is the same as Vanilla
    bool                    isCall;
    bool                    canExerciseEarly;
    ScheduleSP              exerciseSchedule;

    double                  spotAtMaturity;

    DateTime                dateExercised;
    bool                    isExercised;

    CashFlowArraySP         instCoupons;
};
typedef smartPtr<CInstallment> CInstallmentSP;

// helpers
class CInstallmentHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(CInstallment, clazz);
        SUPERCLASS(Generic1Factor);
        EMPTY_SHELL_METHOD(defaultInstallment);
        // same as vanilla
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        FIELD(isCall,        "true = call, false = put");
        FIELD(exerciseSchedule,        "Exercise Schedule");
        FIELD(canExerciseEarly, "Can option be exercised early");
        FIELD(spotAtMaturity,   "underlying spot level when exercised");
        FIELD_MAKE_OPTIONAL(spotAtMaturity);
        FIELD(isExercised, "Indicates whether option has been exercised");
        FIELD_MAKE_OPTIONAL(isExercised);
        FIELD(dateExercised,        "Date on which option has been exercised");
        FIELD_MAKE_OPTIONAL(dateExercised);
        // barrier parts
        FIELD(instCoupons, "installment payment dates and amounts");
    }

    static IObject* defaultInstallment(){
        return new CInstallment();
    }
};


CClassConstSP const CInstallment::TYPE = CClass::registerClassLoadMethod(
    "Installment", typeid(CInstallment), CInstallmentHelper::load);


bool CInstallment::avoidVegaMatrix(const IModel* model)
{
    return false;
}

// constructor
CInstallment::CInstallment(): Generic1Factor(TYPE), isCall(true)
{
    canExerciseEarly = false;
    spotAtMaturity           = 0.0;
    isExercised              = false;
}

void CInstallment::Validate()
{
    static const string method = "CInstallment::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object
    if (!asset){
        throw ModelException(method, "Asset is null");
    }
    if (!discount){
        throw ModelException(method, "Discount YC is null");
    }

    // validate against struck until thoroughly tested and agreed methodology
    if (ccyTreatment == Asset::CCY_TREATMENT_STRUCK) {
        throw ModelException(method, "ccy struck type not supported.");
    }

    // american is turned off for now, need precise time of day definition for coupon
    if (canExerciseEarly){
        throw ModelException(method, "American style not yet allowed.");
    }
    // no need to know if it is exercised or not, as the product is european for the moment.
    /*if (isExercised){
        throw ModelException(method, "European product is exercised at maturity only.");
    }*/

    if (exerciseSchedule->getDates().size() != 1) {
        throw ModelException(method, "European option: only 1 strike allowed in schedule.");
    }

    if ( instCoupons->size() > 0 &&
        (*instCoupons)[instCoupons->size() - 1].date > exerciseSchedule->getDates()[0]) {
        throw ModelException(method, "All Installment coupons must lie on or before maturity.");
    }

    // validation against one coupon and forward start
    if ( fwdStarting && oneContract ) {
        throw ModelException(method, "Cannot be One Contract and Forward Starting");
    }
    

    AssetUtil::assetCrossValidate(asset.get(),
                                  fwdStarting,
                                  startDate,
                                  valueDate,
                                  discount,
                                  this);  
}

// initiate GetMarket 
void CInstallment::GetMarket(const IModel*         model, 
                            const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);

    if (asset.usingCache() || !Tree1fLN::TYPE->isInstance(model))
    {// should always come in - just to cope with old regression convertion
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                                   discount, asset);
    }

    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
    if (premiumSettle.get())
    {
        premiumSettle->getMarket(model, market.get());
    }
}

/** returns the current value date */
DateTime CInstallment::getValueDate() const
{
   return valueDate;
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool CInstallment::priceDeadInstrument(CControl* control, CResults* results) const
{
    double    strike        = 0.0;
    double    value         = 0.0;
    bool      foundExerDate = false;
    DateTime  exerDate;
    bool priceDead = false;

    static string method = "CInstallment::priceDeadInstrument";

    DateTime matDate = exerciseSchedule->lastDate();

    DateTime settlementDate = instSettle->settles(matDate, asset.get());

    if (valueDate >= settlementDate)
    {    // settled already
        results->storePrice(0.0, discount->getCcy());
        addOutputRequests(control, results, 0.0, 0.0);
        priceDead = true;
    }
       
    if (isExercised)
    {// may early exercise
        // check that exercise date is not in the Future
        if ( dateExercised.isGreater(valueDate))
        {
            throw ModelException(method,
                    "Installment terminated on " + 
                    dateExercised.toString() + ". " +
                    "This date is after the current value date (" + 
                    valueDate.toString());
        }

        // check whether exercise date is on a valid date
        for (int i=0 ; i<instCoupons->size() ; ++i)
        {
            if ( dateExercised.equals((*instCoupons)[i].date)) 
            {
                foundExerDate = true;
                exerDate      = dateExercised;
                break;
            }
        }
        if (!foundExerDate)
        {  
            throw ModelException(method, 
                    "Installment cannot be terminated on " + 
                    dateExercised.toString() + ". " +
                    "This date is not in the coupon schedule.");        
        }
    }

    if (foundExerDate)
    {
        priceDead = true;
        if (exerDate.equals(matDate))
        {
            value = GetIntrinsic(spotAtMaturity,
                                 strike,
                                 isCall, 
                                 true) - (*instCoupons)[instCoupons->size()-1].amount;
        }

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

    // record KNOWN_CASHFLOWS
    if ( control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS) ) {
        CashFlowArraySP allCFs = getKnownCashFlows();
        
        OutputRequestUtil::recordKnownCashflows(control,
            results,
            discount->getCcy(),
            allCFs.get()); 
    }

    // PAYMENT_DATES
    if ( control->requestsOutput(OutputRequest::PAYMENT_DATES) ) {
        DateTimeArraySP dates = getPaymentDates();
        
        OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
    }    

    return priceDead;
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool CInstallment::sensShift(Theta* shift)
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

/** when to stop tweaking */
DateTime CInstallment::endDate(const Sensitivity* sensControl) const {
    DateTime matDate = exerciseSchedule->lastDate();
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

void CInstallment::addOutputRequests(Control* control,
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
        InstrumentUtil::recordIndicativeVol(control,results,indVol);
        // FWD_AT_MAT
        try{
            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           matDate,
                                           valueDate,
                                           asset.get());
        }
        catch(exception&)
        {// continue if fwd failed - this hapens now for quanto asset with CEVj vol
        }
    }
}


CashFlowArraySP CInstallment::getKnownCashFlows() const
{
    double scalingFactor,fwdStrt;
    DateTime exerDate;    
    CashFlowArraySP cfl(new CashFlowArray(0));
    DateTime matDate = exerciseSchedule->lastDate();
                                                            
    // computing the scaling factor
    fwdStrt = 0.0;
    if (fwdStarting)
    {
        fwdStrt = asset->fwdValue(startDate);
    }
    scalingFactor = InstrumentUtil::scalePremium(oneContract,fwdStarting,notional,fwdStrt,initialSpot);

    if (isExercised)
    {
        exerDate = dateExercised;
    }
    else
    {
        exerDate = matDate;
    }

    // selecting the dates that are prior to valueDate
    int j = 0;
    int maxSize = instCoupons->size();
    while (j<maxSize && ((*instCoupons)[j].date < valueDate) && ((*instCoupons)[j].date < exerDate) )
    {
        cfl->push_back(CashFlow((*instCoupons)[j].date, (*instCoupons)[j].amount));
        j++;
    }
        
    return cfl;
}

DateTimeArraySP CInstallment::getPaymentDates() const
{
    int maxSize = instCoupons->size();
    DateTimeArraySP payDates(new DateTimeArray(maxSize));
    int j;

    for(j=0; j<maxSize; j++)
    {
        (*payDates)[j] = instSettle->settles((*instCoupons)[j].date, asset.get());
    }

    return payDates;
}


/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP CInstallment::getSensitiveStrikes(OutputNameConstSP outputName,
                                               const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("CInstallment::getSensitiveStrikes", 
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

// product class
class CInstallmentFDProd: public LatticeProdEDRIns
{
public:
    CInstallmentFDProd( const CInstallment * installment, FDModel * mdl ) :
        LatticeProdEDRIns( mdl, 1, 2 ),
        inst( installment )
    {
        if( ! tree1f )
        {
            throw ModelException( "CInstallmentFDProd::CInstallmentFDProd", "Instrument of type "+
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
        static const string method = "CInstallmentFDProd::Init";
        try {
            int i;
            // define start and end points of tree
            DateTimeArray segDates(2);
            segDates[0] = getStartDate();
            segDates[1] = inst->exerciseSchedule->lastDate();

            // timeline configuration
            // 'density factor' for timeline
            IntArray density( 1, 1 );

            // all exercise dates are copied to critical dates
            DateTimeArraySP critDates(new DateTimeArray(inst->exerciseSchedule->getDates()));
            // remove exercise date from crit date
            critDates->erase(critDates->end()-1);

            int num = inst->instCoupons->size();
            for (i=0; i<num; i++) {          
                critDates->push_back((*(inst->instCoupons))[i].date);
            }

            // add div event dates if needed
            EventAssetMove divEvent;
            DateTimeArraySP divCritDates;

            if( tree1f )
            {
                tree1f->NumOfPrice = numPrices;
                tree1f->NumOfInsertNode =
                    ( tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ? 1 : 0); // use one inserted node for exercise boundary

                if( inst->fwdStarting )
                    tree1f->controlSameGridFwdStart(inst->ccyTreatment);                

                if (inst->canExerciseEarly || tree1f->DivsTreatedAsAbsolute()) {
                    // American exercise or dollar div interp treatment
                    DateTime start(segDates[0]);
                    DateTime end(segDates[1]);
                    if (tree1f->DivsTreatedAsAbsolute()) {
                        start = min(start, inst->valueDate);    // fwd starting not supported anyway when dol divs
                        end = max(end, Equity::calcDivTransPeriodEndDate(inst->valueDate));
                    }
                    int numDivs = (int)(4*start.yearFrac(end))+1; // 4 divs per year selected as critical dates
                    if (AssetUtil::getJumpEvents(inst->asset.get(), 
                                                 start, 
                                                 end,
                                                 numDivs,
                                                 divEvent)) {

                        // calculate critical dates
                        divCritDates = divEvent.getCritDate(0, inst->isCall);
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
                    if (inst->isCall){
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

                if( divCritDates.get() )
                    tree1f->addDivCritDates( *divCritDates );
            }

            // add critical dates
            model->addCritDates( *critDates );

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
        static const string method = "CInstallmentFDProd::InitProd";
        try {
            int i;
            // timeline exists at this point
            int lastStep = model->getLastStep();

            initSlices( numPrices );
            initInsertNode();

            CAssetConstSP asset = inst->asset.getSP();
            double fwdAtStart = 1.0;

            stepStrike.resize(lastStep+1);
            stepCanExercise.resize(lastStep+1);
            // ask tree to decide first about steps that can exercise
            AssetUtil::setStepExercise(stepCanExercise,
                                    model->getDates(),
                                    inst->exerciseSchedule,
                                    inst->canExerciseEarly,
                                    asset);

            bool canInterpExSched = (inst->exerciseSchedule->length() > 1);
            // get spot at start if needed
            if (inst->fwdStarting) {
                fwdAtStart = inst->asset->fwdValue(inst->startDate);
            }

            CIntArray  canEx;
            canEx.resize(lastStep+1); // for debug
            canEx[lastStep] = 1;

            // past coupons are not taken into account
            int cfIdx = 0;
            while ( (*(inst->instCoupons))[cfIdx].date.getDate() < model->getDate(0).getDate() ){
                cfIdx ++;
            }

            // validation against coupon dates occurring before start date
            if ( (*(inst->instCoupons))[0].date.getDate() < inst->startDate.getDate() ) {
                throw ModelException(method, "coupon dates prior to start date are not allowed");
            }

            
            // compute installment amount at each step  
            stepCoupon.resize(lastStep+1);
            for (i=0; i<=lastStep; i++) {
                if ( cfIdx < inst->instCoupons->size() && 
                    model->getDate(i).getDate() == (*(inst->instCoupons))[cfIdx].date.getDate() )
                {
                    stepCoupon[i] = (*(inst->instCoupons))[cfIdx].amount;
                    cfIdx++;

                    // validation against negative coupons
                    if (stepCoupon[i]<0.0) {
                        throw ModelException(method, "negative coupons are not allowed");
                    }
                
                    // validation against non-ordered coupon dates
                    if ( (cfIdx < inst->instCoupons->size()) && ((*(inst->instCoupons))[cfIdx-1].date.getDate() > (*(inst->instCoupons))[cfIdx].date.getDate())){
                        throw ModelException(method, "coupon dates must be in increasing order");
                    }

                } else {
                    stepCoupon[i] = 0.0;
                }
            }

            // get strikes
            stepStrike[lastStep] = fwdAtStart*inst->exerciseSchedule->lastValue();

            for (i=0; i<lastStep; i++) {
                if (stepCanExercise[i]) {
                    canEx[i] = 1;
                    if (canInterpExSched) {
                        stepStrike[i] = fwdAtStart*inst->exerciseSchedule->interpolate(model->getDate(i));
                    }
                    else {
                        stepStrike[i] = stepStrike[lastStep]; // strike level is constant for american with no schedule
                    }
                }
                else {
                    canEx[i] = 0;
                }
            }
            if(tree1f) {
                useInsertNode = (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION); 
    #ifdef TREE_DEBUG
                tree1f->stepCanExercise = CIntArraySP(dynamic_cast<CIntArray*>(canEx.clone()));
    #endif
            }
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
        double fwdStartDF = 1.0;
        if (inst->fwdStarting)
        {
            fwdAtStart = inst->asset->fwdValue(inst->startDate);
            fwdStartDF = disc->pv(inst->valueDate, model->getDate(0));
        }

        double scalingFactor = InstrumentUtil::scalePremium(
            inst->oneContract,
            inst->fwdStarting,
            inst->notional,
            fwdAtStart,
            inst->initialSpot);

        return fairValue * scalingFactor * fwdStartDF;
    }

    /** extra output requests */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results)
    {
        // get prices at t=0
        // save price
        double price = scalePremium(model->getPrice0( *slices[0] ), disc);
        results->storePrice(price, disc->getCcy());

        // take care of additional outputs
        if ( control && control->isPricing() )
        {
            DateTime       matDate = inst->exerciseSchedule->lastDate();
            double         indVol;
            // calculate indicative vol
            try{
                if ( matDate.isGreater(inst->valueDate) )
                {
                    // get vol request
                    CVolRequestConstSP lnVolRequest = GetLNRequest();

                    // interpolate the vol
                    CVolProcessedSP  vol(inst->asset->getProcessedVol(lnVolRequest.get()));
                    // cast to the type of vol we're expecting
                    CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                    // this should never happen if our get market data has worked properly
                    if (!vol){
                        throw ModelException("CInstallmentFDProd::recordOutput", 
                                             "No Black Scholes Vol");
                    }

                    // calculate the indicative vol
                    indVol = volBS->CalcVol(getStartDate(), matDate);
                }
                else
                {
                    indVol = 0.0;
                }
            }
            catch(exception&)
            {// continue if indicative vol fail - non BS vol
                indVol = 0.0;
            }

            inst->addOutputRequests(control,
                                    results,
                                    price,
                                    indVol);


            OutputNameConstSP lastStep(new OutputName("LastStepDebug"));
            results->storeGreek(d_prices, Results::DEBUG_PACKET, lastStep);

        }
    }
    
    /** returns a vol request for log-normal vol */
    //for set up timeline only,  to be reviewed/removed */
    virtual CVolRequestConstSP GetLNRequest() const
    {
        // get strike and maturity date from instrument
        DateTime matDate = inst->exerciseSchedule->lastDate();
        double volStrike  = inst->exerciseSchedule->lastValue();

        CVolRequestConstSP volRequest(
            new LinearStrikeVolRequest(volStrike, getStartDate(),
                                       matDate, inst->fwdStarting));
        return volRequest;
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
            prod_BWD_T(  s,
                          step,
                          bot,
                          top,
                          pStart, 
                          pEnd,
                          price);

            //insert nodes
            if (tree1f && tree1f->NumOfInsertNode>0)
            {
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
                prod_BWD(     *insNodes,
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
        static const string method("CInstallmentFDProd::prod_BWD_T");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int i, j;
    
            double settlementPV = inst->instSettle->pvAdjust(inst->exerciseSchedule->lastDate(),
                                                                 inst->discount.get(), 
                                                                 inst->asset.get());
    
            double strike = stepStrike[step] + (inst->isCall? stepCoupon[step] : -stepCoupon[step]) ;
            for (i=pStart; i<=pEnd; i++)
            {
                for (j=bot; j<=top; j++)
                {
                    (p[i])[j] = GetIntrinsic(s[j], strike, inst->isCall, true/* this allow fwd */);
                    (p[i])[j] *= settlementPV;
                }
            }
            if (useInsertNode && tree1f)
                tree1f->SetInsertNodeAndPrice(tree1f->GetSliceIdx(), 0, stepStrike[step], 0, pStart, pEnd, 0.0);

            d_prices = DoubleArraySP(new DoubleArray(top - bot + 1));
            for (int d_i = bot; d_i <= top; d_i++)
            {
                (*d_prices)[d_i-bot] = (p[pStart])[d_i];
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
        static const string method("CInstallmentFDProd::prod_BWD");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int j;
            ASSERT(pStart==0 && pEnd==1);

            if (tree1f && step == model->getLastStep()-1)
            {// do penultimate smoothing
                vector <double> vol_arr, drift_arr;
                double dt = tree1f->getTrdYrFrac(step+1);
                bool fastRoll = tree1f->CalcStepDriftAndVar(s, bot, top, vol_arr, &drift_arr); // just to get drift
                // use GetStepVol() to get vol, do not use GetStepVar() which may not have simple BS conversion
                tree1f->GetStepVol(step, vol_arr, s, bot, top);

                double variance = vol_arr[0]*vol_arr[0]*dt;
                double drift = drift_arr[0];
                double df = inst->discount->pv(model->getDate(step),
                                                   model->getDate(step+1));

                // If dollar dividend treatment, we must adjust fwd and strike by stock floor.
                //   If not dollar dividend treatment, stock floor == 0.0 
                double StockFloor0 = tree1f->GetStockFloor(step);
                double StockFloor1 = tree1f->GetStockFloor(step + 1);
     
                for (j=bot; j<=top; j++)
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
                        // If dollar div treatment (in which case stock floor is <> 0.0), the price of the option
                        //  is that of an option with same maturity written on the pseudo asset (i.e., S - StockFloor)
                        //  and with strike adjusted by the stock floor at maturity (i.e., K - StockFloor).
                        //  NB Black::price takes care of potentially negative fwd and strike values. 
                        (p[0])[j] = (p[1])[j] = Black::price(inst->isCall, 
                                                                 drift * (s[j] - StockFloor0),  // pseudo asset's (step+1)-maturity fwd
                                                                 // as viewed from t = step
                                                                 stepStrike[step+1] - StockFloor1,  // adjusted strike
                                                                 df,
                                                                 variance);
                    }
                }

            }

            // it is just european at the moment, american is easily allowed 
            // but be careful about time of day issue on coupon date !
            for (j=bot; j<=top; j++)
            {
                (p[0])[j] = (p[1])[j] = Maths::max((p[0])[j] - stepCoupon[step], 0.0);
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
                intrinsic = callput*(s[j] - stepStrike[step]);
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
            if (useInsertNode && jIns > bot+1 && jIns < top-1)
                setExerciseBoundary(s, (p[1]), jIns, priceLast, priceCurr, callput, step);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    
    void AdjustDeltaShift(CControl* control) const
    {
        // only do this for tree1f for now
        if (tree1f) {

            double deltaShift = control->getDeltaShiftSize();

            if (Maths::isPositive(deltaShift)){
                DividendListSP dollarDivs = 
                    AssetUtil::getDollarDivsBetweenDates(inst->asset.get(),
                                                         inst->valueDate,
                                                         inst->exerciseSchedule->lastDate());
                bool isFwdStart = inst->fwdStarting && inst->startDate>inst->valueDate;

                tree1f->DeltaSizeAdjust(deltaShift,
                                        inst->valueDate,
                                        inst->exerciseSchedule->firstDate(),
                                        inst->exerciseSchedule->lastDate(), 
                                        inst->exerciseSchedule->length()>1 && !inst->canExerciseEarly,
                                        dollarDivs->getDivAmounts()->size()>0,
                                        isFwdStart);
            }
        }
    }

    void setExerciseBoundary(const double* s, const double* euro, int jIns, 
                             double p1, double p2, 
                             double cp, int step)
    {
        int idx = tree1f->GetSliceIdx();

        // locate node position to insert
        double a = (p2-p1)/(s[jIns]-s[jIns-1]);
        double insNode =s[jIns-1]+ (cp*(s[jIns-1]-stepStrike[step]) - p1)/(a-cp);
        double gap = (s[jIns]-s[jIns-1])/100.0; // within 1% to either side do not need to be treated.
        if (fabs(insNode-s[jIns-1]) > gap && fabs(insNode-s[jIns]) > gap && fabs(p1-p2) > gap)
        {
            // interpolate for european
            double insEuro = euro[jIns-1]+(insNode-s[jIns-1])*(euro[jIns]-euro[jIns-1])/(s[jIns]-s[jIns-1]);
            // set the inserted node and price here
            tree1f->SetInsertNodeAndPrice(idx, 0, insNode, 0, 1, 1, insEuro);
            // intrinhsic for american
            double insPrice = cp*(insNode-stepStrike[step]);
            // set the inserted node and price here
            ASSERT((insPrice >= p2 && insPrice <= p1) || (insPrice <= p2 && insPrice >= p1));
            tree1f->SetInsertNodeAndPrice(idx, 0, insNode, 0, 0, 0, insPrice);
        }
    }

    // --- debug ---
    static DoubleArraySP d_prices;

private:
    const CInstallment* inst;
    bool            useInsertNode;
    vector<double>  stepStrike;
    vector<bool>    stepCanExercise;
    // installment amount at each step
    vector<double>  stepCoupon;
};

DoubleArraySP CInstallmentFDProd::d_prices = DoubleArraySP(   );

/** create a fd payoff tree - new for all fd/tree state variable interface */
FDProductSP CInstallment::createProduct(FDModel* model) const
{
    return FDProductSP( new CInstallmentFDProd(this, model) );
}

extern bool InstallmentLoad()
{
    return true && CInstallment::TYPE;
}

DRLIB_END_NAMESPACE
