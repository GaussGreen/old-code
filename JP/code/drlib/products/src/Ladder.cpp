//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Ladder.cpp
//
//   Description : it's the ladder (yawn)
// 
//   Author      : Andrew J Swain
//
//   Date        : 30 June 2003
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

class Ladder : public Generic1Factor,
               virtual public FDModel::IIntoProduct,
               virtual public LastSensDate
{
public:
    static CClassConstSP const TYPE;  

    virtual void Validate(){
        static const string method = "Ladder::Validate";
        
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

        if (Maths::isNegative(strike)) {
            throw ModelException(method, "strike ("+Format::toString(strike) + 
                                 ") is negative");
        }

        if (isCapped && cap < rungs[rungs.size()-1]) {
            throw ModelException(method, "cap (" + Format::toString(cap) + 
                                 ") is below the top rung (" + 
                                 Format::toString(rungs[rungs.size()-1])+")");
        }
            
        if (rungs[0] < strike) {
            throw ModelException(method, "first rung (" + 
                                 Format::toString(rungs[0]) + 
                                 ") is below the strike (" + 
                                 Format::toString(strike) + ")");
        }            

        for (int i = 1; i < rungs.size(); i++) {
            if (rungs[i] < rungs[i-1]) {
                throw ModelException(method, 
                                     "rungs  must be in ascending order - " 
                                     "rung " + Format::toString(i-1) + 
                                     "(" + Format::toString(rungs[i-1]) + 
                                     ") is higher than " + 
                                     "rung " + Format::toString(i) + 
                                     "(" + Format::toString(rungs[i]) +")");
            }
        }
    }

private:
    friend class LadderHelper;
    friend class LadderFDProd;

    Ladder():Generic1Factor(TYPE), strike(0.0), maxSoFar(0.0),
             isCapped(false), cap(0.0), spotAtMaturity(0.0) {}; 

    Ladder(const Ladder& rhs);
    Ladder& operator=(const Ladder& rhs);
   
    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    bool avoidVegaMatrix(const IModel*model){
        if (CTree1fLV::TYPE->isInstance(model)) {
            return true; // do pointwise instead
        }
        return false;
    }

    // make a LN vol request - not in real use as prefer LV tree
    CVolRequestLN* makeVolRequest() const {
        static const string method("Ladder::makeVolRequest");
        try {
            DateTime imntStartDate = fwdStarting ? startDate : valueDate;

            return new LinearStrikeVolRequest(strike,
                                              imntStartDate,
                                              maturity,
                                              fwdStarting);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Returns all strikes the Ladder is sensitive to  */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model){
        static const string method("Ladder::getSensitiveStrikes");
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
        DateTime instEnd  = instSettle->settles(maturity, asset.get());
        DateTime assetEnd = asset->settleDate(maturity);
        DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;    
    }

    bool sensShift(Theta* shift) {
        static const string method = "Ladder::sensShift";
        try  {
            DateTime newDate = shift->rollDate(valueDate);

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
                strike *= initialSpot;
                cap    *= initialSpot;
                for (int i = 0; i < rungs.size(); i++) {
                    rungs[i] *= initialSpot;
                }
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
    double finalPayoff() const {
        // find highest hit rung
        double highestRung = strike;

        for (int i = 0; i < rungs.size(); i++) {
            if (maxSoFar >= rungs[i]) {
                highestRung = rungs[i];
            }
        }

        // to make payoff equation simpler
        double capLevel = isCapped ? cap: spotAtMaturity;

        double payoff = Maths::max(highestRung-strike, 0.0) + 
            Maths::max(0.0, Maths::min(capLevel, spotAtMaturity) - highestRung);

        if (!oneContract) {
            payoff *= notional/initialSpot;
        }

        return payoff;
    }


    bool priceDeadInstrument(CControl* control, CResults* results) const{
        static const string method = "Ladder::priceDeadInstrument";
        try  {
            bool     deadInstrument = false;
            DateTime end = instSettle->settles(maturity,asset.get());

            if (valueDate >= end) {
                // all over, worth zero
                results->storePrice(0.0, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results);
                }
                              
                deadInstrument = true;  
            }
            else if (valueDate >= maturity ) {
                double value = finalPayoff()*discount->pv(end);
                                
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
            // FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           maturity,
                                           valueDate,
                                           asset.get());

            OutputRequest* request = NULL;
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArray paydates(1,instSettle->settles(maturity,asset.get()));
                OutputRequestUtil::recordPaymentDates(control,results,&paydates); 
            }
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                CashFlowArray cfl;
                if (valueDate >= maturity ) {
                    double value = finalPayoff();
                    DateTime pays = instSettle->settles(maturity, asset.get());
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
        }
        catch (exception&) {
            // don't die if any of these fail
        }
    }  

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

private:
    DateTime    maturity;
    double      strike;
    DoubleArray rungs;
    double      maxSoFar;
    bool        isCapped;
    double      cap;
    double      spotAtMaturity;
};

// product class
class LadderFDProd: public LatticeProdEDRIns
{
public:
    LadderFDProd( const Ladder * lad, FDModel * mdl ) :
        LatticeProdEDRIns( mdl, lad->rungs.size(), 1 + lad->rungs.size() + 1 + (lad->isCapped ? 1 : 0) ),
        inst( lad )
    {
        if( ! tree1f )
        {
            throw ModelException( "LadderFDProd::LadderFDProd", "Instrument of type "+
                                 inst->getClass()->getName() +
                                 " can be priced by CTree1f only" );
        }

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );

        // set up strike/barrier/rebate for each KO

        // how many price 'slices' are there ?
        // overall + 1 per rung, the high rung vanilla and maybe a cap

        koStrike.resize(numIns);
        koBarrier.resize(numIns);
        koRebate.resize(numIns);

        koStrike[0]  = inst->strike;
        koBarrier[0] = inst->rungs[0];
        koRebate[0]  = inst->rungs[0] - inst->strike;
  
        int i;
        for (i = 1; i < numIns; i++) {
            koStrike[i]  = inst->rungs[i-1];
            koBarrier[i] = inst->rungs[i];
            koRebate[i]  = inst->rungs[i] - inst->rungs[i-1];
        }
      
        // and the last vanilla
        vanStrike = inst->rungs[numIns-1];
        cap = inst->isCapped ? inst->cap : 0.;

        if (inst->fwdStarting) {
            double fwd = inst->asset->fwdValue(inst->startDate);
            for (i = 0; i < numIns; i++) {
                koStrike[i]  *= fwd;
                koBarrier[i] *= fwd;
                koRebate[i]  *= fwd;
            }            
            vanStrike *= fwd;
            if (inst->isCapped) {
                cap *= fwd;
            }
        }

        maxSoFar = inst->maxSoFar;

        settles = inst->instSettle->settles(inst->maturity, inst->asset.get());
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
        static const string method = "LadderFDProd::Init";
        try {
            // default to NODE_INSERTION smoothing
            if (tree1f->GetSmoothMethod() == CTree1f::DEFAULT) {
                tree1f->SetSmoothMethod(CTree1f::NODE_INSERTION);
            }

            tree1f->NumOfPrice = numPrices;
            tree1f->NumOfInsertNode =
                ( tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ? inst->rungs.size() : 0 );

            // define start and end points of tree
            DateTimeArray segDates(2);
            segDates[0] = getStartDate();
            segDates[1] = inst->maturity;

            if( inst->fwdStarting )
                tree1f->controlSameGridFwdStart(inst->ccyTreatment);
            
            // timeline configuration
            // 'density factor' for timeline
            IntArray density( 1, 1 );

            // don't bother with dollar divs
            tree1f->SetDivAmountTreatment(false);

            // kill off control var - unlikely to be very useful here
            tree1f->DEBUG_UseCtrlVar = false;

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
        static const string method = "LadderFDProd::InitProd";
        try {
            int i;
            // timeline exists at this point
            int lastStep = model->getLastStep();

            int numOfInsertNode = tree1f->NumOfInsertNode;

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
        static const string method("LadderFDProd::preCalc");
        try {
            int i;
            vector<double> vol;
            int idx = tree1f->getSliceIndex(step);

            for (i = 0; i < koBarrier.size(); i++) {
                adjBarrier[i] = koBarrier[i];

                // taking care of discrete monitoring, treat as daily monitor
                // get vol at barrier
                double level = adjBarrier[i];
                tree1f->GetStepVol(step, vol, &level, 0, 0);
                // please fix me .......
                Barrier::BarrierAdjustment(vol[0],true,adjBarrier[i]);

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

        // for iPrice = 0, all insert nodes is necessary.  
        // But shouldn't return "false", as it's going to no more move.
        if (iPrice == 0){
        }else{
            // activation is controlled by priority rather than flag.  Flag is currently working entirely.
            // by giving negative priority, we can switch off on the specific level of insStock.
            // need to operate on "1-idx" because RollTree looks 1-idx.
            for (int i = 0; i < (int)adjBarrier.size(); i++) {
                tree1f->SetInsertNode(1-idx, i, adjBarrier[i], (iPrice-1 == i ? 1 : -1));     
            }       
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
            prod_BWD_T(   s,
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
                // at top layer (i.e. vanilla option), the inserted node price is calculated by RollTree from 1-CurrIdx.
                // this could be not good and generate vial, because the higher moment is not consistent.
                // so, replaced by interpolated value from CurrIdx.
                const double *s_ins = insNodes->getValues();
                int idx = tree1f->getSliceIndex(step);
                for (int iPrice=(inst->isCapped ? pEnd-1 : pEnd); iPrice<=pEnd; iPrice++){
                    for (int j=0; j<tree1f->NumOfInsertNode; j++){
                        double interpPrice = tree1f->TreeInterp(s_ins[j], true, iPrice, 0);
                        //Use priority = 1, but it would be turned off at moveInsertNode.
                        tree1f->SetInsertNodeAndPrice(idx, j, s_ins[j], 1, iPrice, iPrice, interpPrice);    
                    }
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
        static const string method("LadderFDProd::prod_BWD_T");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int numRungs = koStrike.size();
            int i;
            int j;

            double settlePV = inst->instSettle->pvAdjust(inst->maturity,
                                                        inst->discount.get(), 
                                                        inst->asset.get());

            // we've 1 slice per rung + a vanilla + maybe a cap
            for (i = 1; i <= numRungs; i++) {
                for (j=bot; j<=top; j++) {
                    // price all the KO's
                    (p[i])[j] = priceKOAtMat(i-1, s[j])*settlePV;
                }
            }

            // then the top rung vanilla
            for (j=bot; j<=top; j++) {
                (p[numRungs+1])[j] = GetIntrinsic(s[j], vanStrike, true, true);
                (p[numRungs+1])[j] *= settlePV;
            }

            // and take off any cap
            if (inst->isCapped) {
                for (j=bot; j<=top; j++) {
                    (p[numRungs+2])[j] = -GetIntrinsic(s[j],cap,true,true);
                    (p[numRungs+2])[j] *= settlePV;
                }                
            }
            
            for (j=bot; j<=top; j++) {
                (p[0])[j] = 0.0;
                for (i = 1; i < numPrices; i++) {
                    (p[0])[j] += (p[i])[j];
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
        static const string method("LadderFDProd::prod_BWD");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int numRungs = koStrike.size();
            int i;
            if (!stepIsKO[step]) {
                return;
            }
            
            DateTime stepDate = model->getDate(step);
            double pvToEnd    = inst->discount->pv(stepDate, settles);

            for (i = 1; i <= numRungs; i++) {
                for (int j=bot; j<=top; j++) {
                    // price all the KO's
                    (p[i])[j] = priceKOEarly(i-1, (p[i])[j], s[j], pvToEnd);
                }
            }  

            for (int j=bot; j<=top; j++) {
                (p[0])[j] = 0.0;
                for (i = 1; i < numPrices; i++) {
                    (p[0])[j] += (p[i])[j];
                }
            }                   
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }
    
    // own methods
    double priceKOAtMat(int idx, double spot)
    {
        double value   = GetIntrinsic(spot, koStrike[idx], true, true);
        //double barrier = adjBarrier[idx];
        // no barrier adjustment to avoid double counting if
        // spot lies between adjusted and regular barrier
        double barrier = koBarrier[idx]; 
        if (spot > barrier-FP_MIN || maxSoFar > koBarrier[idx]-FP_MIN) {
            value = koRebate[idx];
        }
        return value;
    }

    double priceKOEarly(int idx, double live, double spot, double pv)
    {
        double value   = live;
        double barrier = adjBarrier[idx];
        if (spot > barrier-FP_MIN || maxSoFar > koBarrier[idx]-FP_MIN) {
            value = koRebate[idx]*pv;
        }
        return value;
    }

private:
    const Ladder* inst;
    vector<bool>  stepIsKO;

    DateTime      settles;
    double        maxSoFar;

    DoubleArray   koStrike;
    DoubleArray   koBarrier;
    DoubleArray   koRebate;

    double        vanStrike;
    double        cap;

    // barriers at a tree slice after discrete monitor adjustment
    vector< double > adjBarrier;
};

/** create a fd payoff tree - new for all fd/tree state variable interface */
FDProductSP Ladder::createProduct(FDModel* model) const
{
    return FDProductSP( new LadderFDProd(this, model) );
}

class LadderHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Ladder");
        REGISTER(Ladder, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        EMPTY_SHELL_METHOD(defaultLadder);
        FIELD(maturity, "maturity");
        FIELD(strike, "strike");
        FIELD(rungs, "rungs");
        FIELD(maxSoFar, "max so far");
        FIELD_MAKE_OPTIONAL(maxSoFar);
        FIELD(isCapped, "isCapped");
        FIELD(cap, "cap");
        FIELD_MAKE_OPTIONAL(cap);
        FIELD(spotAtMaturity, "spot @ maturity");
        FIELD_MAKE_OPTIONAL(spotAtMaturity);
    }

    static IObject* defaultLadder(){
        return new Ladder();
    }
};

CClassConstSP const Ladder::TYPE = CClass::registerClassLoadMethod(
    "Ladder", typeid(Ladder), LadderHelper::load);


/* for class loading */
bool LadderLoad() {
    return (Ladder::TYPE != 0);
}


DRLIB_END_NAMESPACE

