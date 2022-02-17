//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AmericanLadder.cpp
//
//   Description : American Ladder  + rebate at maturity / exedate.
// 
//   Author      : Keiji Kitazawa
//
//   Date        : 02 Oct 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/Asset.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/Tree1fLV.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/VegaParallel.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

class AmericanLadder : public Generic1Factor,
                       virtual public FDModel::IIntoProduct,
                       virtual public LastSensDate
{
public:
    static CClassConstSP const TYPE;  

    virtual void Validate(){
        static const string method = "AmericanLadder::Validate";
        int i;
        DateTime maturity = exerciseSchedule->lastDate();
        DoubleArray strikes = exerciseSchedule->getValues();
        double tmpV;

        // Call Generic1Factor validate()
        validate();

        if (fwdStarting && startDate <= valueDate) {
            throw ModelException(method, 
                                 "instrument is fwd starting but start date ("+
                                 startDate.toString() + ") is <= today ("+
                                 valueDate.toString());
        }

        if (fwdStarting && startDate >= maturity) {
            throw ModelException(method, 
                                 "instrument is fwd starting but start date ("+
                                 startDate.toString() + 
                                 ") is on or after maturity ("+
                                 maturity.toString());
        }

        if (fwdStarting && oneContract) {
            throw ModelException(method, 
                                 "Can't be forward starting and one contract");
        }

        // let's not go there shall we
        if (instSettle->isMargin()) {
            throw ModelException(method, "margin settlement not supported");
        }
        if (instSettle->isPhysical()) {
            throw ModelException(method, "physical settlement not supported");
        }

        if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
            throw ModelException(method, "ccy struck not supported");
        }
         
        if (rungs.empty()) {
            throw ModelException(method, "need at least one rung");
        }
            
        if (Maths::isNegative(maxSoFar)) {
            throw ModelException(method, "max so far (" + 
                                 Format::toString(maxSoFar) + 
                                 ") is negative");
        }

        if (isCapped && Maths::isNegative(cap)) {
            throw ModelException(method, "cap (" + Format::toString(cap) + 
                                 ") is negative");
        }

        for (i = 0; i < strikes.size(); i++){
            if (Maths::isNegative(strikes[i])) {
                throw ModelException(method, "strike[" + Format::toString(i) + "] ("
                                    +Format::toString(strikes[i]) + ") is negative");
            }
        }

        if (isCapped && cap < rungs[rungs.size()-1]) {
            throw ModelException(method, "cap (" + Format::toString(cap) + 
                                 ") is below the top rung (" + 
                                 Format::toString(rungs[rungs.size()-1])+")");
        }

        for (i = 0; i < strikes.size(); i++){
            if (rungs[0] < strikes[i]) {
                throw ModelException(method, "first rung (" + 
                                     Format::toString(rungs[0]) + 
                                     ") is below the strike[" + Format::toString(i) + "]  (" + 
                                     Format::toString(strikes[i]) + ")");
            }            
        }

        if (rungs.size() != ladderScales.size()) {
                throw ModelException(method, "the number of rung (" + 
                                     Format::toString(rungs.size()) + 
                                     ") is not sames as the number of ladderScales(" + 
                                     Format::toString(ladderScales.size()) + ")");
        }

        for (i = 1; i < rungs.size(); i++) {
            if (rungs[i] < rungs[i-1]) {
                throw ModelException(method, 
                                     "rungs  must be in ascending order - " 
                                     "rung " + Format::toString(i-1) + 
                                     "(" + Format::toString(rungs[i-1]) + 
                                     ") is higher than " + 
                                     "rung " + Format::toString(i) + 
                                     "(" + Format::toString(rungs[i]) +")");
            }
            if (Maths::isNegative(ladderScales[i-1])){
                throw ModelException(method, 
                                     "ladderScales should be positive number - " 
                                     "ladderScales[" + Format::toString(i-1) + 
                                     "] = " + Format::toString(ladderScales[i-1]) + " ");
            }                
        }

        // validate if isExercised case
        if (isExercised){
            // check exe date is specified.
/*            if(!!dateExercised) 
                throw ModelException(method,"It's already exercised, but no info about exercised date." );
*/
            // check that exercise date is not in the Future
            if ( dateExercised.isGreater(valueDate)) {
                throw ModelException(method,
                                     "Option exercised on " + 
                                     dateExercised.toString() + ". " +
                                     "This date is after the current value date (" + 
                                     valueDate.toString());
            }

            // check whether exercise date is valid and find corresponding strike
            // note: not checking for weekend yet
            try{
                tmpV = exerciseSchedule->interpolate(dateExercised);
            }
            catch (ModelException& e) {
                e = ModelException::addTextToException(e, method + ": Invalid dateExercised");
                throw ModelException(e, method);
            }
        }

        // avoid multiple exedate if it's not american.
        if (!canExerciseEarly) {
            if (exerciseSchedule->length() > 1){
                throw ModelException(method, 
                                     "your input for 'canExerciseEarly' = false, but you set multiple exercise schedule.");
            }
        }

        // Check the schedule doesn't contains holiday, for the case interp = "N".
        DateTimeArray exedates = exerciseSchedule->getDates();
        for (i = 0; i<exedates .size(); i++){
            if (AssetUtil::getHoliday(asset.get())->isBusinessDay(exedates[i])==false) {
                throw ModelException(method, "You can not set holiday for exerciseSchedule.  ("
                    +Format::toString(i+1) + "-th dates in exercise schedule is holiday).");
            }
        }

        // if rebate is not completed input, make it null.
        if (!!rebate) {
            if (rebate->length() == 0)
                rebate = ScheduleSP(   );
            else if (rebate->length() > 0) {
                DateTimeArray rebDates = rebate->getDates();
                try{
                    double tmpValue;
                    for (i = 0; i<rebDates.size(); i++){
                        tmpValue = exerciseSchedule->interpolate(rebDates[i]);
                    }
                }
                catch (exception& ) {
                    throw ModelException(method, "Your rebate on [" + rebDates[i].toString() +
                                                "] is not found in exercise schedule. Check dates or InterpType.");
                }
            }
        }
                                        

    }

private:
    friend class AmericanLadderHelper;
    friend class AmericanLadderProd;
    friend class AmericanLadderFDProd;

    AmericanLadder():Generic1Factor(TYPE), maxSoFar(0.0),
             isCapped(false), cap(0.0), spotAtMaturity(0.0) {
        earlyExePayType = "DEFAULT";
        isExercised = false;
        rebate = ScheduleSP(   );    
        RebateNotScaled = false;
        IntraDayMonitor = false;
    }; 

    AmericanLadder(const AmericanLadder& rhs);
    AmericanLadder& operator=(const AmericanLadder& rhs);
   
    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    bool avoidVegaMatrix(const IModel*model){
        if (CTree1fLV::TYPE->isInstance(model)) {
            return true; // do pointwise instead
        }
        return false;
    }

    // make a LN vol request - not in real use as prefer LV tree
    CVolRequestLN* makeVolRequest() const {
        static const string method("AmericanLadder::makeVolRequest");
        try {
            DateTime imntStartDate = fwdStarting ? startDate : valueDate;
            DateTime maturity = exerciseSchedule->lastDate();
            double strike = exerciseSchedule->lastValue();
            if (fwdStarting == false)
                strike *= initialSpot;

            return new LinearStrikeVolRequest(strike,
                                              imntStartDate,
                                              maturity,
                                              fwdStarting);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Returns all strikes the AmericanLadder is sensitive to  */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model){
        static const string method("AmericanLadder::getSensitiveStrikes");
        try {
            DoubleArraySP sensStrikes(new DoubleArray(0));
            if (avoidVegaMatrix(model)) {
                throw ModelException(method, 
                                     "VEGA_MATRIX is not valid for this "
                                     "instrument");
            }
            CVolRequestConstSP volRequest(makeVolRequest());

            SensitiveStrikeDescriptor sensStrikeDesc;
            sensStrikeDesc.forwardOnly = false;
            asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                       sensStrikeDesc, sensStrikes);

            return sensStrikes;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** when to stop tweaking */
    DateTime endDate(const Sensitivity* sensControl) const{
        DateTime maturity = exerciseSchedule->lastDate();
        DateTime instEnd  = instSettle->settles(maturity, asset.get());
        DateTime assetEnd = asset->settleDate(maturity);
        DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;    
    }

    bool sensShift(Theta* shift) {
        static const string method = "AmericanLadder::sensShift";
        try  {
            DateTime newDate = shift->rollDate(valueDate);
            DateTime maturity = exerciseSchedule->lastDate();
            
            if ((newDate >= maturity && valueDate < maturity) ||
                (valueDate == maturity && Maths::isZero(spotAtMaturity))) {
                spotAtMaturity = asset->getThetaSpotOnDate(shift, maturity);
                maxSoFar = Maths::max(maxSoFar, spotAtMaturity);
            }

            if (fwdStarting && newDate.isGreaterOrEqual(startDate) &&
                startDate.isGreaterOrEqual(valueDate))
            {
                fwdStarting = false;
                initialSpot = asset->getThetaSpotOnDate(shift, startDate);
/*                exerciseSchedule->scale(initialSpot);        
                cap    *= initialSpot;
                for (int i = 0; i < rungs.size(); i++) {
                    rungs[i] *= initialSpot;
                }
                if (!RebateNotScaled)
                {
                    if (!!rebate)
                        rebate->scale(initialSpot);
                }*/
            }

            // roll the parent (updates value date etc)
            Generic1Factor::sensShift(shift);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }    
        return true; // our components have theta type sensitivity
    }
      
    // what's the final payoff (excluding settlement etc)
    // Here, early exercise is also handled.    
    double exePayoff() const {
        double spotScale = fwdStarting ? asset->fwdValue(startDate) : initialSpot;
        double strike = spotScale * exerciseSchedule->lastValue();
        double s_level = spotAtMaturity;
        DateTime exeDate = isExercised ? dateExercised : exerciseSchedule->lastDate();

        // find highest hit rung
        double maxSpot = Maths::max(maxSoFar, s_level);    //re-write the maxSoFar.
        double highestRung = strike;
        for (int i = 0; i < rungs.size(); i++) {
            if (maxSpot >= spotScale * rungs[i]) {
                highestRung = spotScale * rungs[i];
            }
        }

        if (isExercised) {
            strike = spotScale * exerciseSchedule->interpolate(dateExercised);            
            if (earlyExePayType == "NO_INTRINSIC") {
                s_level = highestRung;
            }
        }
        
        // to make payoff equation simpler
        double capLevel = isCapped ? spotScale * cap: s_level;

        double payoff = Maths::max(highestRung-strike, 0.0) + 
            Maths::max(0.0, Maths::min(capLevel, s_level) - highestRung);

        if (!!rebate){
            try{
                if (RebateNotScaled)
                    payoff += rebate->interpolate(exeDate);
                else
                    payoff += spotScale * rebate->interpolate(exeDate);
            }
            catch (exception& ) {
                // if it's fail, just pass it.
            }
        }


        double scalingFactor = InstrumentUtil::scalePremium(
            oneContract,
            false,
            notional,
            0.0,
            initialSpot);

        payoff *= scalingFactor;

        return payoff;
    }


    bool priceDeadInstrument(CControl* control, CResults* results) const{
        static const string method = "AmericanLadder::priceDeadInstrument";
        try  {
            double   value;
            bool     deadInstrument = false;
            DateTime maturity = exerciseSchedule->lastDate();

            DateTime end = instSettle->settles(maturity,asset.get());

            if (valueDate >= end) {
                // all over, worth zero
                results->storePrice(0.0, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results);
                }
                              
                deadInstrument = true;  
            }
            else if (valueDate >= maturity || isExercised) {
                if(isExercised){
                    end = instSettle->settles(dateExercised, asset.get());
                }
                value = exePayoff()*discount->pv(end);
                                
                results->storePrice(value, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results);
                }
                              
                deadInstrument = true;  
            } 
            return deadInstrument;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }            
    }

    /** extra output requests */
    void recordOutputRequests(Control* control, 
                              Results* results) const {
        try {
            DateTime maturity = exerciseSchedule->lastDate();
            DateTime exeDate = isExercised ? dateExercised : maturity;
            
            // FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           maturity,
                                           valueDate,
                                           asset.get());

            OutputRequest* request = NULL;
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArray paydates(1,instSettle->settles(exeDate,asset.get()));
                OutputRequestUtil::recordPaymentDates(control,results,&paydates); 
            }
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                CashFlowArray cfl;
                if (valueDate >= exeDate ) {
                    double value = exePayoff();
                    DateTime pays = instSettle->settles(exeDate, asset.get());
                    CashFlow cf(pays, value);
                    cfl.push_back(cf);
                }
       
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        &cfl);   
            }  
            request = control->requestsOutput(OutputRequest::IND_VOL);
            if (request) {
                DateTime imntStartDate = fwdStarting ? startDate : valueDate;
                CVolRequestLNSP volRequest(makeVolRequest());
                CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
                double ivol = volBS->CalcVol(imntStartDate, maturity);
                results->storeRequestResult(request, ivol); 
            }

            request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
            if (request) {
                // report ladder level today 
                double spotScale = fwdStarting ? asset->fwdValue(startDate) : initialSpot;
                BarrierLevelArray levels;
                double barlevel = 0.0;
                
                // store all ladder level as upper barrier
                for (int i=0; i<rungs.size(); i++){
                    barlevel = rungs[i] * spotScale; 
                    if (maxSoFar < rungs[i] * spotScale){
                        break;
                    }
                }

                levels.push_back(BarrierLevel(true, valueDate, barlevel,IntraDayMonitor));
                OutputRequestUtil::recordBarrierLevels(control,
                                                       results,
                                                       asset->getTrueName(),
                                                       &levels);
            } 

        }
        catch (exception&) {
            // don't die if any of these fail
        }
    }  

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

private:
    //DateTime    maturity;
    ScheduleSP  exerciseSchedule;
    bool        canExerciseEarly;
    string      earlyExePayType;

    bool        isExercised;
    DateTime    dateExercised;

    //double      strike;
    DoubleArray rungs;
    double      maxSoFar;
    bool        isCapped;
    double      cap;
    double      spotAtMaturity;
    DoubleArray ladderScales;
    ScheduleSP  rebate;
    bool        RebateNotScaled;
    bool        IntraDayMonitor;

};

class AmericanLadderFDProd: public LatticeProdEDRIns
{
public:
    AmericanLadderFDProd( const AmericanLadder * lad, FDModel * mdl ) :
        LatticeProdEDRIns( mdl, lad->rungs.size(), lad->rungs.size() + 1 ),
        inst( lad )
    {
        if( ! tree1f )
        {
            throw ModelException( "AmericanLadderFDProd::AmericanLadderFDProd", "Instrument of type "+
                                 inst->getClass()->getName() +
                                 " can be priced by CTree1f only" );
        }

        // first: set discount curve
        if( tree1f )
            tree1f->setDiscountCurve( inst->discount.getSP() );

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );

        koStrike.resize(numIns);
        koBarrier.resize(numIns);
        koRebate.resize(numIns);

        koStrike[0]  = inst->exerciseSchedule->lastValue();
        koBarrier[0] = inst->rungs[0];
        koRebate[0]  = (inst->rungs[0] - inst->exerciseSchedule->lastValue())*inst->ladderScales[0];
  
        int i;
        for (i = 1; i < numIns; i++) {
            koStrike[i]  = inst->rungs[i-1];
            koBarrier[i] = inst->rungs[i];
            koRebate[i]  = koRebate[i-1] + (inst->rungs[i] - inst->rungs[i-1])*inst->ladderScales[i];
        }
      
        // and the last vanilla
        //vanStrike = inst->rungs[numRungs-1];
        // vanilla strike.
        vanStrike = koStrike[0];
        cap = inst->isCapped ? inst->cap : 0;

        rebScale = 1.0; // rebate Scale is 1.0 as default.
        if (inst->fwdStarting)
            spotScale = inst->asset->fwdValue(inst->startDate);
        else
            spotScale = inst->initialSpot;

        for (i = 0; i < numIns; i++) {
            koStrike[i]  *= spotScale;
            koBarrier[i] *= spotScale;
            koRebate[i]  *= spotScale;
        }            
        vanStrike *= spotScale;
        if (inst->isCapped) {
            cap *= spotScale;
        }
        if (!!inst->rebate){
            if (!inst->RebateNotScaled){
                rebScale = spotScale;                    
            }
        }

        maxSoFar = inst->maxSoFar;
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
        static const string method = "AmericanLadderFDProd::Init";
        try {
            // default to NODE_INSERTION smoothing
            if (tree1f->GetSmoothMethod() == CTree1f::DEFAULT) {
                tree1f->SetSmoothMethod(CTree1f::NODE_INSERTION);
            }

            tree1f->NumOfPrice = numPrices;
            tree1f->NumOfInsertNode =
                ( tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ? numIns : 0 );

            if( inst->fwdStarting )
                tree1f->controlSameGridFwdStart(inst->ccyTreatment);
            
            // all exe dates are critical dates
            DateTimeArray criticalDates = inst->exerciseSchedule->getDates();

            // don't bother with dollar divs
            tree1f->SetDivAmountTreatment(false);

            // kill off control var - unlikely to be very useful here
            tree1f->DEBUG_UseCtrlVar = false;

            // add critical dates
            model->addCritDates( criticalDates );

            // define start and end points of tree
            DateTimeArray segDates(2);
            segDates[0] = getStartDate();
            segDates[1] = inst->exerciseSchedule->lastDate();
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
        static const string method = "AmericanLadderFDProd::InitProd";
        try {
            int i;
            int lastStep = model->getLastStep();

            initSlices( numPrices );
            initInsertNode();

            // set up adjBarrier
            adjBarrier.resize( numIns );

            // set up an array of flags indicating if time step is a KO date
            stepIsKO.resize(lastStep+1);

            // disallow rung hit at weekends/holidays
            CAssetConstSP asset = inst->asset.getSP();

            stepIsKO[lastStep] = true; // by definition
            stepIsKO[0] = AssetUtil::getHoliday(asset.get())->isBusinessDay(model->getDate(0));

            for (i = 1; i < lastStep; i++) {
                stepIsKO[i] = AssetUtil::hasTradingTime(asset,
                                                        model->getDate(i-1),
                                                        model->getDate(i));
            } 

            bool isAmerican = inst->canExerciseEarly;
            if ( inst->exerciseSchedule->getInterp() == "N" ){
                isAmerican = false;
            }

            // set up an array of flags indicating if time step is a early exercise date
            stepCanExercise.resize(lastStep+1);
            // ask tree to decide first about steps that can exercise
            AssetUtil::setStepExercise(stepCanExercise,
                                    model->getDates(),
                                    inst->exerciseSchedule,
                                    isAmerican,
                                    asset);
            // set up rebate, too.
            stepIsRebate.resize(lastStep+1);
            if (!!inst->rebate) {
                AssetUtil::setStepExercise(stepIsRebate,
                                        model->getDates(),
                                        inst->rebate,
                                        isAmerican,
                                        asset);               
            }
            else {
                for (i=0;i<=lastStep;i++){
                    stepIsRebate[i] = false;
                }
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

        return (fairValue*scalingFactor*fwdStartDF);
    }

    /** extra output requests */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results)
    {
        // get prices at t=0
        // save price
        double price = scalePremium(model->getPrice0( *slices[0] ), disc);
        results->storePrice(price, disc->getCcy());

        // throw this back to the instrument itself
        inst->recordOutputRequests(control, results);
    }
    
    /** returns a vol request for log-normal vol */
    //for set up timeline only,  to be reviewed/removed */
    virtual CVolRequestConstSP GetLNRequest() const
    {
        return CVolRequestConstSP( inst->makeVolRequest() );
    }

    /** called before each update() */
    virtual void preCalc(int step)
    {
        static const string method("AmericanLadderFDProd::preCalc");
        try {
            int i;
            vector<double> vol;
            int idx = tree1f->getSliceIndex(step);

            for (i = 0; i < koBarrier.size(); i++) {
                adjBarrier[i] = koBarrier[i];

                if (!inst->IntraDayMonitor){
                    // taking care of discrete monitoring, treat as daily monitor
                    // get vol at barrier
                    double level = adjBarrier[i];
                    tree1f->GetStepVol(step, vol, &level, 0, 0);

                    // please fix me .......
                    Barrier::BarrierAdjustment(vol[0], true, adjBarrier[i]); 
                }

                if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION) { 
                    // apparently last param = 1 means KO
                    tree1f->SetInsertNode(idx, i, adjBarrier[i], 1); 
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }


    /** Activate only next ladder level. */
    // NODE_INSERTION is not good at vanilla option pricing, 
    // because it gives some bias and it's necessary extremely high steps to be disappered.
    // so Ladder is going to activate only next Ladder Level.
    virtual bool moveInsertNode(int currStep, int iPrice)
    {
        if (tree1f->GetSmoothMethod() != CTree1f::NODE_INSERTION)
            return false;
        int idx = tree1f->getSliceIndex(currStep);

        // if the iPrice is top, it's plain call (or capped call).
        // no need insert node, so shouldn't go into the loop.  
        // activation is controlled by priority rather than flag.  Flag is currently working entirely.
        // by giving negative priority, we can switch off on the specific level of insStock.
        // need to operate on "1-idx" because RollTree looks 1-idx.
        for (int i = 0; i < (int)adjBarrier.size(); i++) {
            tree1f->SetInsertNode(1-idx, i, adjBarrier[i], (iPrice <= i ? 1 : -1));     
        }       
        return true;
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
            prod_BWD_T( s,
                        step,
                        bot,
                        top,
                        pStart, 
                        pEnd,
                        price);

            //insert nodes
            if (tree1f && tree1f->NumOfInsertNode>0){
                prod_BWD_T( *insNodes,
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
                // at pEnd (i.e. vanilla option), the inserted node price is calculated by RollTree from 1-CurrIdx.
                // this could be not good and generate vial, because the higher moment is not consistent.
                // so, replaced by interpolated value from CurrIdx.
                const double *s_ins = insNodes->getValues();
                int idx = tree1f->getSliceIndex(step);
                for (int j=0; j<tree1f->NumOfInsertNode; j++){
                    double interpPrice = tree1f->TreeInterp(s_ins[j], true, pEnd, 0);
                    //Use priority = 1, but it would be turned off at moveInsertNode.
                    tree1f->SetInsertNodeAndPrice(idx, j, s_ins[j], 1, pEnd, pEnd, interpPrice);    
                }
                prod_BWD( *insNodes,
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
        static const string method("AmericanLadderFDProd::prod_BWD_T");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int numRungs = koStrike.size();
            int i, j;
            DateTime matDate = inst->exerciseSchedule->lastDate();
            double settlePV = inst->instSettle->pvAdjust(matDate,
                                                        inst->discount.get(), 
                                                        inst->asset.get());

            // calculate step rebate for early exercise
            double reb = 0.0;
            if (!!inst->rebate){
                if(stepIsRebate[step])
                    reb = rebScale * inst->rebate->interpolate(matDate);
            }

            // we've a vanilla (main) + 1 slice per rung.
            // cap is attached for each slice.
            
            // the main slice is no ladder achived. i.e. vanilla
            for (j=bot; j<=top; j++) {
                (p[0])[j] = GetIntrinsic(s[j], vanStrike, true, true);
                if (inst->isCapped) {
                    (p[0])[j] -= GetIntrinsic(s[j],cap,true,true);
                }
                (p[0])[j] += reb;
                (p[0])[j] *= settlePV;
            }

            // +  various state variable w.r.t. achived ladder
            double value  =0.0;
            for (i = 1; i <= numRungs; i++) {
                for (j=bot; j<=top; j++) {
                    // price all the KO's
                    value   = GetIntrinsic(s[j], koStrike[0], true, true);

                    // no barrier Adjustemet is necessary for Mat, with closing sampling.
                    // Also, not judget the S is above or not ladder, but just floored by koRebate,
                    // becaure if it's scaled, should be protected by koRebate.
                    value = Maths::max(value, koRebate[i-1]);

                    if (inst->isCapped) {
                        value -= GetIntrinsic(s[j],cap,true,true);
                    }                    
                    (p[i])[j] = (value + reb ) * settlePV;
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
        static const string method("AmericanLadderFDProd::prod_BWD");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int numRungs = koStrike.size();
            int i, j;

            DateTime stepDate = model->getDate(step);
            DateTime setDate = inst->instSettle->settles(stepDate,inst->asset.get());
            double pvToSet    = inst->discount->pv(stepDate, setDate);

            // calculate step rebate for early exercise
            double reb = 0.0;
            if (!!inst->rebate){
                if(stepIsRebate[step])
                    reb = rebScale * inst->rebate->interpolate(stepDate);
            }

            if (stepIsKO[step]) {
                // Each slice prices the ladder option.
                // slice (i) is already achived i-1 ladder.
                // It must be caluclate from upper slice to consider state variable.
                for (i = numRungs-1; i>=0; i--){
                    for (j=bot; j<=top; j++) {
                        (p[i])[j] = priceLadder(s[j], i, p, j);
                    }                   
                }

                /* main tree calculation.  Using the above state variable trees.
                for (j=bot; j<=top; j++) {
                    (p[0])[j] = priceLadder(s[j], 0, price, j);
                } */                                  
            }
            if (stepCanExercise[step]) {
                // start fron numRungs, not numRungs-1, because allready breached tree
                // may early exercise.
                for (i = numRungs; i>=0; i--){
                    for (j=bot; j<=top; j++) {
                        (p[i])[j] = priceExeValue(s[j],i,p,j, pvToSet, reb);
                    }                   
                }
            }
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }
    
    // own methods
    double priceLadder(double spot, int iAchivedLevel, const vector< double * > & p, int j){
        int ilevel = iAchivedLevel;
        ilevel = iSearchLadderLevel(spot, iAchivedLevel);
        return (p[ilevel])[j];
    }

    double priceExeValue(double spot, int iAchivedLevel, const vector< double * > & p, int j, 
                        double pvToSettle, double reb){
        int ilevel = iAchivedLevel;
        double value = 0.0;
        if (inst->earlyExePayType == "DEFAULT"){
            value = GetIntrinsic(spot, koStrike[0], true, true);
            if (inst->isCapped) {
                value -= GetIntrinsic(spot,cap,true,true);
            }
        }
        ilevel = iSearchLadderLevel(spot, iAchivedLevel);
        if (ilevel>0){
            value = Maths::max( value, koRebate[ilevel-1]*pvToSettle );        // Locked the level as max so far.
        }
        value += reb*pvToSettle;
        value = Maths::max(value, (p[iAchivedLevel])[j]); // Compare intrinsinc v.s. option value.  (American exe).
        return value;
    }

    int iSearchLadderLevel(double spot, int iAchivedLevel){
        int ilevel = iAchivedLevel;
        double barrier = 0.0;
        for (int i = iAchivedLevel; i < koStrike.size(); i++){
            barrier = adjBarrier[i];
            if (spot > barrier *(1.0 - FP_MIN) || maxSoFar > koBarrier[i]*(1.0-FP_MIN))
                ilevel ++;
            else
                break;
        }
        return ilevel;
    }

private:
    const AmericanLadder* inst;
    vector<bool>  stepIsKO;
    vector<bool>  stepCanExercise;
    vector<bool>  stepIsRebate;

    double        maxSoFar;
    double        rebScale;
    double        spotScale;

    DoubleArray   koStrike;
    DoubleArray   koBarrier;
    DoubleArray   koRebate;

    double        vanStrike;
    double        cap;

    // barriers at a tree slice after discrete monitor adjustment
    vector< double > adjBarrier;
};

/** create a fd payoff tree - new for all fd/tree state variable interface */
FDProductSP AmericanLadder::createProduct(FDModel* model) const
{
    return FDProductSP( new AmericanLadderFDProd(this, model) );
}

class AmericanLadderHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("AmericanLadder");
        REGISTER(AmericanLadder, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        EMPTY_SHELL_METHOD(defaultAmericanLadder);
        FIELD(exerciseSchedule,        "Exercise Schedule");
        //FIELD(maturity, "maturity");
        //FIELD(strike, "strike");
        FIELD(rungs, "rungs");
        FIELD(ladderScales, "Scaling factors for Ladder payoff");
        FIELD(maxSoFar, "max so far");
        FIELD_MAKE_OPTIONAL(maxSoFar);
        FIELD(isCapped, "isCapped");
        FIELD(cap, "cap");
        FIELD_MAKE_OPTIONAL(cap);
        FIELD(spotAtMaturity, "spot @ maturity");
        FIELD_MAKE_OPTIONAL(spotAtMaturity);
        FIELD(canExerciseEarly, "Can option be exercised early");
        FIELD(earlyExePayType, "PayType for early Exercise.  DEFAULT (usual American), NO_INTRINSIC");
        FIELD_MAKE_OPTIONAL(earlyExePayType);
        FIELD(isExercised, "Indicates whether option has been exercised");
        FIELD_MAKE_OPTIONAL(isExercised);
        FIELD(dateExercised,        "Date on which option has been exercised");
        FIELD_MAKE_OPTIONAL(dateExercised);
        FIELD(rebate,        "rebate paid when it's early exercised.");
        FIELD_MAKE_OPTIONAL(rebate);
        FIELD(RebateNotScaled, "for fwd start option. true if rebate is not scaled, false(default)=it is scaled");
        FIELD_MAKE_OPTIONAL(RebateNotScaled);
        FIELD(IntraDayMonitor, "true for continuous monitoring, false if once a day");
        FIELD_MAKE_OPTIONAL(IntraDayMonitor);
        
    }

    static IObject* defaultAmericanLadder(){
        return new AmericanLadder();
    }
};

CClassConstSP const AmericanLadder::TYPE = CClass::registerClassLoadMethod(
    "AmericanLadder", typeid(AmericanLadder), AmericanLadderHelper::load);


/* for class loading */
bool AmericanLadderLoad() {
    return (AmericanLadder::TYPE != 0);
}


DRLIB_END_NAMESPACE

