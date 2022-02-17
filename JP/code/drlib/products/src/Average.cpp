//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Average.hpp
//
//   Description : Average instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 9 March 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Average.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Black.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/LinearStrikeSpreadVolRequest.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Future.hpp"
#include "edginc/XCB.hpp"
#include "edginc/UnitXCBWithVol.hpp"
#include "edginc/PercXCBWithVol.hpp"

//FOR FD
#include "edginc/LastSensDate.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/VegaParallel.hpp"

#include "edginc/UtilFuncs.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/Schedule.hpp"

DRLIB_BEGIN_NAMESPACE

/* retrieve market data needed by Average - just valueDate, asset and
   discount yield curve */
void Average::GetMarket(const IModel*          model, 
                        const CMarketDataSP    market) {
    market->GetReferenceDate(valueDate);
    CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                               discount, asset);
    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
    if (premiumSettle.get()){
        premiumSettle->getMarket(model, market.get());
    }
}

/** what's today ? */
DateTime Average::getValueDate() const {
    return valueDate;
}

/** credit support */
CreditSupportSP Average::createCreditSupport(CMarketDataSP market){
    return CreditSupportSP(new AverageCreditSupport(this, market));
}

SampleListSP Average::getAveIn() const
{
    return SampleListSP(   );
}

/** when to stop tweaking */
DateTime Average::endDate(const Sensitivity* sensControl) const {
    DateTime end;

    if (RhoPointwise::TYPE->isInstance(sensControl)) {
        DateTime instEnd  = instSettle->settles(maturity, asset.get());
        DateTime assetEnd = asset->settleDate(maturity);
        end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    }
    else {
        end = maturity;
    }
    return end;
}

static double avgBlack(
    bool    isCall,
    double  avgForward,      
    double  effectiveStrike, 
    double  discFactor,      
    double  avgVariance) {
    static const string method = "avgBlack";
    try {
        double value;
        // Call Black with average variance and average forward
        // If option is a call and adjusted strike is negative then (S-K)
        // will always be positive, as K is negative. Therefore, premium
        // is simply PV(S-K).
        // If option is a put and strike is negative then premium is zero.

        if (isCall && Maths::isNegative(effectiveStrike)) {
            value = discFactor * (avgForward - effectiveStrike);
        }        
        else if (!isCall &&  effectiveStrike < -DBL_EPSILON) {
            value = 0.0;
        }
        else {
            value= Black::price(isCall, 
                                avgForward, 
                                effectiveStrike,
                                discFactor,
                                avgVariance);

        }  

        return value;
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

// might want to reconsider this as a static given the volume
// of class data coming across the interface
static void recordExtras(
    Control*              control,
    CResults*             results,
    const CAsset*         asset,
    const YieldCurve*     discount,
    const DateTime&       today,
    const DateTime&       maturity,
    double                strike,
    double                avgFwd,
    double                variance,
    double                premium,
    const InstrumentSettlement* premiumSettle,
    const InstrumentSettlement* instSettle,
    const SampleList*     avgOut) {
    static const string method = "recordExtras";
    try {
        OutputRequest* request = NULL;

        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       maturity,
                                       today,
                                       asset);
        
        if ((request = 
             control->requestsOutput(OutputRequest::EFFECTIVE_STRIKE))){
            results->storeRequestResult(request, strike);
        }
        if ((request = control->requestsOutput(OutputRequest::AVG_FWD))) {
            results->storeRequestResult(request, avgFwd);
        }
        if ((request = control->requestsOutput(OutputRequest::AVG_VARIANCE))) {
            results->storeRequestResult(request, variance);
        }
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTime paymentDate = instSettle->settles(maturity, asset);
            DateTimeArray date(1, paymentDate);
            OutputRequestUtil::recordPaymentDates(control,results,&date); 
        }
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            // if we're past the last sample then we know the payoff
            // given fv (fv = pv*c) => back out final payment as we have
            // all the data to hand
            DateTime end = avgOut->getLastDate();
            if (today.isGreater(end)) {
                DateTime paymentDate = instSettle->settles(maturity, asset);
                double pv = discount->pv(paymentDate);
                CashFlow cf(paymentDate, premium/pv);
                CashFlowArray cfl(1, cf);
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        &cfl);   
            }
        }


        // DELAY_PRICE
        InstrumentUtil::delayPriceHelper(control,
                                         results,
                                         premium,
                                         today,
                                         discount,
                                         asset,
                                         premiumSettle);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/*************************************   AverageSpot class  ****************************************/
class AverageSpot : public Average, 
                    public IAverageSpot,
                    virtual public CClosedFormLN::IIntoProduct,
                    virtual public IMCIntoProduct,
                    virtual public Theta::Shift,
                    public virtual FDModel::IIntoProduct{
private:
    // additional fields specific to average spot
    bool         fwdStarting;
    DateTime     startDate;
    bool         oneContract;
    double       notional;
    double       initialSpot;
    
    string       interpMethod; // FD: interpolation method
    int          dimAvgOutFD;  // FD: number of slices in average dimension for FD
    double       stateStdev;   // FD: number of standard deviations for average state variables, default=3.5
    bool         UseCtrlVar;   // FD: flag responsible for use of control variate technique 
public:
    static CClassConstSP const TYPE;
    
    AverageSpot();
    AverageSpot(bool                        isCall,
                const DateTime&             maturity,
                double                      strike,
                const SampleList*           avgOut,
                const InstrumentSettlement* instSettle,
                const InstrumentSettlement* premiumSettle,
                const CAsset*               asset,
                string                      ccyTreatment,
                const YieldCurve*           discount,
                const DateTime&             valueDate,
                bool                        fwdStarting,
                const DateTime&             startDate,
                bool                        oneContract,
                double                      notional,
                double                      initialSpot);

    // how we interpolate vol
    CVolRequestLN* volInterp(double strike) const;

    // indicative vol
    void indVol(CVolProcessedBS* vol, 
                OutputRequest*   request, 
                CResults*        results) const;

    virtual void Validate();
    
    /** Implementation of ClosedForm::IIntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    /** Implementation of IMCIntoProduct interface */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    virtual bool sensShift(Theta* shift);

    virtual CSensControl* AlterControl(
        const IModel*       modelParams,
        const CSensControl* sensControl) const;

    /** returns the strike which is relevant for the vol interpolation */
    double getAdjustedStrike() const;
     
    double calcAvgVariance( double strike ) const;
    void priceLN(
        Control* control,
        CResults* results,
        double strike,
        double avgVariance ) const;

private:
    friend class AverageSpotHelper;
    friend class AverageSpotClosedForm;
    friend class AverageSpotMC;

    friend class AverageSpotFD;      // FD: FD implementation of Average Spot 
    friend class AverageSpotState;   // FD: out-state variable

    void priceLN(Control* control, CResults* results) const;
    bool getFwdStartDate(DateTime &_startDate) const
    {
        if( fwdStarting ) _startDate = startDate; 
        return fwdStarting;
    }
};
/*************************************   AverageSpotClosedForm  ****************************************/
class AverageSpotClosedForm: public CClosedFormLN::IProduct{
private:
    const AverageSpot* avg; // a reference

public:
    AverageSpotClosedForm(const AverageSpot* avg): avg(avg){}

    void price(CClosedFormLN*  model,
               Control*        control, 
               CResults*       results) const{
        avg->priceLN(control, results);
    }
};
/************************ STATE VARIABLES *********************/
typedef TreeSliceLayer::StateSupport::InterpMethod AVG_INTERP;
class AverageSpotFD;

/* AverageSpotState class: averaging out state variable */
class AverageSpotState : public TreeSliceLayer::StateSupport{
public:
    static double calcAvg(const AverageSpot* inst){
        double avg = inst->avgOut->averageToDate(inst->startDate);
        return avg;
    }
    AverageSpotState(const AverageSpotFD* prod, const AverageSpot* inst, 
                     FDProductSP payoffIndex,   int dim, AVG_INTERP method) :
        StateSupport( "AverageSpotState", calcAvg(inst), method ),
        prod(prod),
        inst(inst),
        payoffIndex(payoffIndex),
        dim(dim),
        sampleCount(inst->avgOut->numDates(0)) {}

    /* state variable support using slice layers, TreeSliceLayer::StateSupport:: */
    // update state contents
    void update(int step, const TreeSlice& s);
    // set grid levels
    void setGrid( int step, bool init = false );
    // compute grid updates at a reset step. This is for each slice node level.
    virtual void transition(
        const vector< const TreeSlice * > & srcSlices,
        const vector< double > & gridIn,
        vector< double > & gridOut ) const;

    const AverageSpotFD*    prod;
    const AverageSpot*      inst;
    FDProductSP             payoffIndex;
    int                     dim;
    double                  avgSoFar;
    int                     firstResetStep;
    BoolArray               isResetStep;
    int                     sampleCount;
};
/************************ state operation methods *********************/
void AverageSpotState::setGrid( int step, bool init )
{
    if( init || step != firstResetStep )
    {
        // computing average range using stateStdev
        DateTimeArray sampleDates(inst->avgOut->getDates());
        DoubleArray values(inst->avgOut->getValues());
        DoubleArray weights(inst->avgOut->getWeights());
        sampleDates.resize(sampleCount);
        values.resize(sampleCount);
        weights.resize(sampleCount);

        SampleList samples(sampleDates, values, weights);
       
        // compute avg varince
        // computing average range using stateStdev
        ATMVolRequest volReq;
        CVolProcessedBSSP  volBS(inst->asset->getProcessedVol(&volReq));
        double avgVariance = samples.averageVariance(volBS.get(), inst->valueDate, true);

        double fwd = samples.expectedAverage(inst->asset.get(), inst->valueDate);
        double maxAvg = fwd*exp(-0.5*avgVariance + inst->stateStdev*sqrt(avgVariance));
        double minAvg = fwd*exp(-0.5*avgVariance - inst->stateStdev*sqrt(avgVariance));
      
        int n = dim; // *sqrt((double)sampleCount/inst->monitorDates.size());
        n = Maths::max(n, 2);
        // build a log array for average
        double increment = exp((log(maxAvg) - log(minAvg))/(n - 1.0));

        currGrid.resize(n);
        currGrid[0] = minAvg;
        for (int j = 1; j < n; ++j) {
            currGrid[j] = currGrid[j-1] * increment;
        }
    }
    else
    {
        // collapse grid to dim=1
        currGrid.resize(1); // !!! currGrid.resize(1, todayValue) does not work because resize down
        currGrid[0] = todayValue;
    }

    if( init )
        prevGrid = currGrid;

    populateGrid( step );
}
// compute grid transition through a reset/monitoring
void AverageSpotState::transition(
    const vector< const TreeSlice * > & srcSlices,
    const vector< double > & gridIn,
    vector< double > & gridOut ) const
{
    int nbGrid = gridOut.size();
    ASSERT( nbGrid == gridIn.size() );

    double s = srcSlices[ 0 ]->calc();

    for( int i = 0; i < nbGrid; ++i )
        gridOut[ i ] = ( gridIn[ i ] * ( sampleCount - 1 ) + s ) / sampleCount;
}
// interpolation etc.
void AverageSpotState::update(int step, const TreeSlice& s){
    if( isResetStep[step] ){
        setGrid( step );
        vector< const TreeSlice * > srcSlices( 1, &s );
        interpolate( step, srcSlices );
        --sampleCount;
    }
}
typedef refCountPtr<AverageSpotState > AverageSpotStateSP;
/************************ end state variables operation *********************/
/************************      FDProduct class          *********************/
class AverageSpotFD : public LatticeProdEDR, 
                      virtual public IFDProductLN {
private:
    friend class AverageSpotState;

    const AverageSpot*        inst;
    //double                    initialSpot;
    double                    strikeAvg; 
    //vector<bool>              stepExercise;
    BoolArray                 stepExercise;
    
    // price slices:
    // the price of the option, European, American or Bermudan
    TreeSliceSP value;
    // the price of the European option used for control variate technique
    TreeSliceSP valueCtrlVariate;
    
    // shorter for convenience
    AVG_INTERP interpMethod;

public: 
    AverageSpotStateSP avgOutState;

    /** constructor */
    AverageSpotFD(const AverageSpot* inst, FDModel* model) : 
        LatticeProdEDR(model), 
        inst(inst)
    {
        if (inst->fwdStarting) {
            strikeAvg = inst->strike * inst->asset->fwdValue(inst->startDate);
        }
        else
            strikeAvg = inst->strike; 

        double sumSoFar = inst->avgOut->sumToDate(inst->valueDate);

        strikeAvg = strikeAvg - sumSoFar;

        payoffIndex = model->createProduct(IProdCreatorSP(new
            IndexSpecEQ(inst->asset.getName(), inst->asset, inst->ccyTreatment)));

        if (inst->interpMethod == "QUADRATIC")
            interpMethod = TreeSliceLayer::StateSupport::INTERP_QUADRATIC;
        else if (inst->interpMethod == "LINEAR")
            interpMethod = TreeSliceLayer::StateSupport::INTERP_LINEAR;
        else
            throw ModelException("AverageSpotFD: ", "unknown interpMethod: " + inst->interpMethod);

        avgOutState = AverageSpotStateSP(new AverageSpotState(this, inst, payoffIndex, inst->dimAvgOutFD, interpMethod));  
    }

    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type);

     /** local method, product payoff method at maturity */
    void prod_BWD_T(const TreeSlice & s, int step);

    /** local method, this is payoff boundary condition, for KO, early exercise etc */
    void prod_BWD(const TreeSlice & s, int step);
    
    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(Control*  control) const;

    /** set flags for averaging time steps,
        return the first step of the monitorDates after valueDate */
    int setStepFlags(BoolArray& stepFlags, const DateTimeArray& stepDates, const DateTimeArray& monitorDates);
    
    /** initialising and setting product variables this is called per pricing call before each pricing */
    virtual void initProd();

    /** output prices and any additional results */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results);

    /** scale for notional */ 
    double scalePremium(const double& fairValue, YieldCurveConstSP disc);

    /** post price */
    double postPrice(double euroPrice, double exer_premium, YieldCurveConstSP disc);
    
    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const {return inst->ccyTreatment;}

    /** Returns start date */
    virtual DateTime getStartDate() const {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

    /** for the LogNormal model */
    CVolRequestLNSP getVolInterp(int iAsset) const {
        // get strike and maturity date from instrument
        DateTime matDate = inst->avgOut->getLastDate();
        
        double volStrike = (inst->fwdStarting) ? inst->strike : strikeAvg;

        DateTime imntStartDate = getStartDate();
        CVolRequestLNSP   reqarr;         
        reqarr = CVolRequestLNSP(new LinearStrikeVolRequest(volStrike, imntStartDate, matDate, inst->fwdStarting));        
        return reqarr;
    }
};

/** create a fd payoff product */
FDProductSP AverageSpot::createProduct(FDModel* model) const {
    return FDProductSP(new AverageSpotFD(this, model));
}
/* Update function*/
void AverageSpotFD::update(int& step, FDProduct::UpdateType type)
{
    const TreeSlice & s = payoffIndex->getValue(step);

    if( type == FDProduct::BWD_T )
    {
        avgOutState->setGrid( step, true );
        prod_BWD_T(s, step);
    }
    else if( type == FDProduct::BWD )
        prod_BWD(s, step);

    avgOutState->update(step, s);
}
/* local method, product payoff method at maturity */
void AverageSpotFD::prod_BWD_T(const TreeSlice & s, int step)
{
    static const string method = "AverageSpotFD::prod_BWD_T";
    try 
    {
        double settlementPV = inst->instSettle->pvAdjust(model->getDate(step), 
                                                                    inst->discount.get(), 
                                                                    inst->asset.get());
        if( ! inst->isCall )
            settlementPV = -settlementPV;

        const TreeSlice & outGrid = avgOutState->getGridSlice();
        *value = smax( 0., settlementPV * ( outGrid - strikeAvg ) );
        if ( inst->UseCtrlVar )
            *valueCtrlVariate = *value;
    } 
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}
/* local method, product payoff method before maturity */
void AverageSpotFD::prod_BWD(const TreeSlice & s, int step)
{
    static const string method = "AverageSpotFD::prod_BWD";
    try 
    {
        double settlementPV = inst->instSettle->pvAdjust(model->getDate(step), 
                                                                    inst->discount.get(), 
                                                                    inst->asset.get());
        if( ! inst->isCall )
            settlementPV = -settlementPV;

        // computes exercise price for american options
        if( !stepExercise[step] )
            return;

        const TreeSlice & outGrid = avgOutState->getGridSlice();
        *value = smax( *value, settlementPV * ( outGrid - strikeAvg ) );
    } 
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}
/** initialise tree1f - allow product customisation 
    must not init product variables here, use initProd() instead */
void AverageSpotFD::init(CControl* control) const
{
    static const string method = "AverageSpotFD::init()";
    try 
    {
        if (inst->dimAvgOutFD < 2) 
        {
            throw ModelException(method,"dimAvgOutFD must ba at least 2 when pricing with tree");
        }
        // customize tree parameters here and set up the tree
        DateTimeArray segDates;
        segDates.resize(2);
        segDates[0] = getStartDate();
        segDates[1] = inst->exerciseSchedule->lastDate();
        //segDates[1] = inst->maturity;
        IntArray density(1,1);

        // all averaging dates are copied to critical dates
        model->addCritDates(inst->avgOut->getDates());   
        // all exercise dates are copied to critical dates
        model->addCritDates(inst->exerciseSchedule->getDates());
        // prepare model set up
        model->initSegments(segDates, density);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/** set flags for averaging time steps,
    return the first step of the monitorDates after valueDate */
int AverageSpotFD::setStepFlags(BoolArray& stepFlags, const DateTimeArray& stepDates, const DateTimeArray& monitorDates) {

    static const string method = "AverageSpotFD::setStepAverage";
    try {
        int iMonitor = 0;
        while (monitorDates[iMonitor] < stepDates[0]) {
            iMonitor++;
        }

        int firstStep = -1;

        for (int iStep = 0; iStep < stepDates.size(); iStep++) {
            if (stepDates[iStep].equals(monitorDates[iMonitor], true)) {
                stepFlags[iStep] = true;
                if( firstStep < 0 ) firstStep = iStep;
                iMonitor++;
            }
            else {
                stepFlags[iStep] = false;
            }
        }
        return firstStep;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/** initialising and setting product variables
    this is called per pricing call before tree sweep call (after InitTree)  */
void AverageSpotFD::initProd()
{
    static const string method = "AverageSpotFD::initProd()";
    try 
    {
        int lastStep = model->getLastStep();
        stepExercise.resize(lastStep + 1);
       
        /*AssetUtil::setStepExercise(stepExercise,
                                   model->getDates(),
                                   inst->exerciseSchedule,
                                   inst->isAmerican,
                                   inst->asset.getSP());*/

        setStepFlags(stepExercise, model->getDates(), inst->exerciseSchedule->getDates());

        model->registerStateSupport( avgOutState.get() );

        startDEV( value = model->createSlice() ); 
        if ( inst->UseCtrlVar )
            startDEV( valueCtrlVariate = model->createSlice() ); 

        avgOutState->isResetStep.resize(lastStep + 1);
        avgOutState->firstResetStep =
            setStepFlags(avgOutState->isResetStep, model->getDates(), inst->avgOut->getDates());
        if(model->getDate(0) == inst->valueDate) // value date is already counted as history
            avgOutState->isResetStep[0] = false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/** output results */
void AverageSpotFD::recordOutput(Control* control, YieldCurveConstSP disc, Results* results) 
{
    double scaledPrice;
    double price = model->getPrice0( *value );
    if ( inst->UseCtrlVar ){
        double priceCtrlVariate = model->getPrice0( *valueCtrlVariate );
        scaledPrice = postPrice( price, price - priceCtrlVariate, disc);
    }
    else
        scaledPrice = scalePremium( price, disc );

    results->storePrice( scaledPrice, disc->getCcy() );

    // take care of additional outputs
    if (control && control->isPricing()) {

        DateTime matDate = inst->avgOut->getLastDate();

        OutputRequest* request = NULL;
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTime paymentDate = inst->instSettle->settles(matDate, inst->asset.get());
            DateTimeArray date(1, paymentDate);
            OutputRequestUtil::recordPaymentDates(control,results,&date); 
        }
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            if (inst->valueDate.isGreater(matDate)) {
                DateTime paymentDate = inst->instSettle->settles(matDate, inst->asset.get());
                double pv = disc->pv(paymentDate);
                CashFlow cf(paymentDate, scaledPrice/pv);
                CashFlowArray cfl(1, cf);
                OutputRequestUtil::recordKnownCashflows(control, results, disc->getCcy(), &cfl);   
            }
        }
        double        indVol;
        // calculate indicative vol
        if (matDate.isGreater(inst->valueDate)) {
            DateTime imntStartDate = getStartDate();

            // get vol request
            //double volStrike  = strike;
            double volStrike  = inst->fwdStarting ? inst->strike : strikeAvg;
            LinearStrikeVolRequest volRequest(volStrike, imntStartDate, matDate, inst->fwdStarting);

            try {
                // interpolate the vol
                CVolProcessedSP  vol(inst->asset->getProcessedVol(&volRequest));
                // cast to the BS vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                // this should never happen if our get market data has worked properly
                if (!volBS) {
                    throw ModelException("AverageSpotFD::recordOutput", "No Black Scholes Vol");
                }
                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
            }
            catch (exception& ) {
                indVol = 0.0;
            }
        }
        else {
            indVol = 0.0;
        }
    }
}
/* premium scaling */
double AverageSpotFD::scalePremium(const double& fairValue, YieldCurveConstSP disc) {

    double fwdAtStart = 0.0;
    double fwdStartDF = 1.0;
    if (inst->fwdStarting) {
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
// postPrice process
double AverageSpotFD::postPrice(double euroPrice, double exer_premium, YieldCurveConstSP disc)
{
    double price;

    AverageSpotClosedForm closedForm(inst);
    CClosedFormLN model;
    Control ctrl;
    CResults result;
    closedForm.price(&model, &ctrl, &result);
    double closedFormPrice = result.retrievePrice();

    price = scalePremium( exer_premium, disc ) + closedFormPrice;
    
    return price;
}
/************************ end of FDProduct class *********************/

static SimSeriesSP makeSimSeries(const SampleList* sample){
    SimSeriesSP simSeries(new SimSeries(1));
    sample->recordDates(simSeries.get());
    return simSeries;
}
// SampleList implements IPastValues
//IPastValuesConstSP((const IPastValues*)inst->avgOut.get()),

// Monte Carlo - use by all flavours of average
class Average::MC: public IMCProduct, virtual public IMCProductLN {
private:
    const Average*        inst;
    IRefLevel::IMCPathSP  avgOutSample;
    bool                  useRatio; // formula to use in payoff
    double                strike;
    double                mult;
    double                numDates; // cached for performance
public:
    
    // equivalent to InstIntoMCProduct 
    MC(const Average*            inst,
       const IRefLevelConstSP&   refLevel, // how to 'avg in'
       const DateTime&           simStartDate, // when simulation starts
       const IPastValuesConstSP& mcPastValues, // historic values
       bool                      useRatio, // which formula to use in payoff
       double                    strike,   // see payOff method
       double                    mult):
        IMCProduct(inst->asset.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  refLevel,
                  makeSimSeries(inst->avgOut.get()),
                  mcPastValues,
                  inst->instSettle.get(),
                  inst->maturity),
        inst(inst), useRatio(useRatio), strike(strike), mult(mult),
        numDates(inst->avgOut->numDates(0)){
        /* just need to create object which allows us to get sample to
           calculate level - we are treating the average out value as
           a 'reference level' */
        avgOutSample = 
            IRefLevel::IMCPathSP(
                inst->avgOut->createMCPath(inst->valueDate,
                                           DoubleArray(1), /* value is
                                                              irrelevant */
                                           inst->avgOut.get()));

    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /* express everything as either
       a) mult * Max(out/in - strike, 0)  or
       b) mult * Max(out - in * strike, 0) 
       with both mult and strike constants (for calls at least) */
    void payoff(const IPathGenerator*  pathGen,   // (I)
                IMCPrices&                prices) {  // (O)
        // only do future (or rather avoid doing partial past)
        if (pathGen->end(0) == numDates){
            // first start with 'in' value
            double inValue = pathGen->refLevel(0,0);
            // then do out value
            const double* path = pathGen->Path(0,0); // access path
            path += pathGen->begin(0); // skip over past
            double outValue = avgOutSample->refLevel(0, path);
            double price = useRatio? 
                (outValue/inValue - strike): (outValue - inValue * strike);
            if (!inst->isCall){
                price = -price;
            }
            prices.add(price < 0.0? 0.0: (price * mult));
        }
    }
    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {
        // choose how to interpolate the vol 
        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here
        /* we don't use the refLevel available from the pathGenerator since
           we do it ourselves (this methodology needs to be reviewed) */
        reqarr[0] = CVolRequestLNSP(inst->volInterp(inst->strike));
        return reqarr;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* AverageSpot::createProduct(const MonteCarlo* model) const {
    double mult, percStrike;
    bool   useRatio;
    double myInitialSpot = oneContract? 1.0: initialSpot;
    if (fwdStarting) {
        // forward starting => notional based, and % strike
        useRatio = true; // do payoff as (avg out/start value - % strike)
        mult = notional;
        percStrike = strike;
    } else {
        useRatio = false; // do payoff as scaling * (avg out - startVal * K)
        mult = oneContract? 1.0: notional/initialSpot;
        percStrike = strike/myInitialSpot;
    }
    // v simple RefLevel
    IRefLevelSP refLevel(IRefLevel::Util::makeTrivialAverage(startDate));
    // need to combine start date with avg out dates
    DateTimeArray allDates = avgOut->getDates(); // take copy
    allDates.insert(allDates.begin(), startDate); // add value
    // similarly for levels
    DoubleArray allVals = avgOut->getValues(); // take copy
    allVals.insert(allVals.begin(), initialSpot); // add value
    DoubleMatrix allValsAsMatrix(allVals);
    IPastValuesSP pastValues(IPastValues::Util::makeSimple(allDates, 
                                                           allValsAsMatrix));
    MC*  prod = new MC(this, refLevel, startDate, pastValues,
                       useRatio, percStrike, mult);
    return prod;
}
    
AverageSpot::AverageSpot() : Average(TYPE), notional(0.0), initialSpot(0.0),dimAvgOutFD(151), 
                             stateStdev(3.5), interpMethod("QUADRATIC"), UseCtrlVar(false) {}

AverageSpot::AverageSpot(bool                        isCall,
                         const DateTime&             maturity,
                         double                      strike,
                         const SampleList*           avgOut,
                         const InstrumentSettlement* instSettle,
                         const InstrumentSettlement* premiumSettle,
                         const CAsset*               asset,
                         string                      ccyTreatment,
                         const YieldCurve*           discount,
                         const DateTime&             valueDate,
                         bool                        fwdStarting,
                         const DateTime&             startDate,
                         bool                        oneContract,
                         double                      notional,            
                         double                      initialSpot) : 
    Average(TYPE), 
    fwdStarting(fwdStarting),
    startDate(startDate),
    oneContract(oneContract),
    notional(notional),
    initialSpot(initialSpot),
    interpMethod("QUADRATIC"),
    dimAvgOutFD(151),
    stateStdev(3.5),
    UseCtrlVar(false)
{
    this->isCall            = isCall;
    this->maturity          = maturity;
    this->strike            = strike;
    this->avgOut            = SampleListSP(copy(avgOut));
    this->exerciseSchedule  = exerciseSchedule;
    this->isAmerican        = isAmerican;
    this->UseCtrlVar        = UseCtrlVar;
    this->instSettle        = InstrumentSettlementSP(copy(instSettle));
    if (premiumSettle) {
        this->premiumSettle = InstrumentSettlementSP(copy(premiumSettle));
    }
    this->asset         = CAssetWrapper(copy(asset));
    this->ccyTreatment  = ccyTreatment;
    this->discount      = YieldCurveWrapper(copy(discount));
    this->valueDate     = valueDate;
}

bool AverageSpot::sensShift(Theta* shift) {
    static const string method = "AverageSpot::sensShift";
    try {
        // calculate the new date
        DateTime newDate = shift->rollDate(valueDate);
        bool useSpot = shift->useAssetFwds()? false : true;

        // If fwd start date falls between value date and theta date 
        // must multiply the percentage strike by spot (theta) or 
        // fwd at start date (theta forward spot) to convert to an 
        // absolute strike, as will no longer be forward starting.
        // Also have to set initial spot to same level.
        
        if ( fwdStarting                           && 
             startDate.isGreaterOrEqual(valueDate) && 
             newDate.isGreaterOrEqual(startDate)   ) {
            // returns the fwd price if shift is a theta fs
            initialSpot = asset->getThetaSpotOnDate(shift, startDate);
            strike *= initialSpot;

            // not fwd starting anymore 
            fwdStarting = false;
        }  

        // fill in any samples
        avgOut->roll(asset.get(), valueDate, newDate, useSpot);
                           
        // roll today 
        valueDate = newDate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
    return true; // our components have theta type sensitivity
}

/* determine how to interpolate on a vol surface for a non fwd starting
 * average price option. Use effective strike (strike adjusted by 
 * historic samples) scaled by 1/(total future weights) for interpolation.
 */
CVolRequestLN* AverageSpot::volInterp(double strike) const {
    static const string method = "AverageSpot::volInterp";
    try {
        CVolRequestLN* interp = 0;
        
        double sumSoFar = avgOut->sumToDate(valueDate);
        double adjStrike;

        if (avgOut->countFutureSamples(valueDate) > 0) {
            double futureWeight = avgOut->futureWeight(valueDate);
        
            if (Maths::isZero(futureWeight)) {
                throw ModelException(method, 
                                     "total of future avg out weights "
                                     "must be > 0.0");
            }

            adjStrike = (strike - sumSoFar)/futureWeight;
        }
        else {
            adjStrike = strike - sumSoFar;
        }
            
        // if effective strike is -ve, use ATM vol 
        if (Maths::isPositive(adjStrike)) {
            DateTime effectiveEnd;
            if (fwdStarting) {
                effectiveEnd = avgOut->expectedEndDate(valueDate);
            }

            interp = new LinearStrikeVolRequest(adjStrike, 
                                                fwdStarting ? startDate : 
                                                              valueDate, 
                                                fwdStarting ? effectiveEnd :
                                                              maturity,
                                                fwdStarting);
        }
        else {
            interp = new ATMVolRequest();
        }            
        return interp;
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

/* determine how to interpolate on a vol surface for a non fwd starting
 * average price option. Use effective strike (strike adjusted by 
 * historic samples) scaled by 1/(total future weights) for interpolation.
 */
CVolRequestLN* Average::volInterpSpot(const DateTime& valueDate,
                                      const DateTime& startDate,
                                      const DateTime& maturity,
                                      bool  fwdStarting,
                                      double strike,
                                      const SampleList* avgOut)
{
    static const string method = "AverageSpot::volInterp";
    try {
        CVolRequestLN* interp = 0;
        
        double sumSoFar = avgOut->sumToDate(valueDate);
        double adjStrike;

        if (avgOut->countFutureSamples(valueDate) > 0) {
            double futureWeight = avgOut->futureWeight(valueDate);
        
            if (Maths::isZero(futureWeight)) {
                throw ModelException(method, 
                                     "total of future avg out weights "
                                     "must be > 0.0");
            }

            adjStrike = (strike - sumSoFar)/futureWeight;
        }
        else {
            adjStrike = strike - sumSoFar;
        }
            
        // if effective strike is -ve, use ATM vol 
        if (Maths::isPositive(adjStrike)) {
            DateTime effectiveEnd;
            if (fwdStarting) {
                effectiveEnd = avgOut->expectedEndDate(valueDate);
            }

            interp = new LinearStrikeVolRequest(adjStrike, 
                                                fwdStarting ? startDate : 
                                                              valueDate, 
                                                fwdStarting ? effectiveEnd :
                                                              maturity,
                                                fwdStarting);
        }
        else {
            interp = new ATMVolRequest();
        }            
        return interp;
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

 
void AverageSpot::indVol(CVolProcessedBS* vol, 
                         OutputRequest*   request, 
                         CResults*        results) const {
    static const string method = "AverageSpot::indVol";
    // do 'indicative vol' - don't fail the whole pricing if this fails 
    try {
        double ivol;

        if (avgOut->countFutureSamples(valueDate) > 0) {
            DateTime endDate = avgOut->expectedEndDate(valueDate);

            ivol = vol->CalcVol(fwdStarting ? startDate:valueDate, endDate);
        }
        else {
            // all samples set
            ivol = 0.0;
        }
          
        results->storeRequestResult(request, ivol);
    }
    catch (exception& ) {
    }
}

void AverageSpot::Validate() {
    static const string method = "AverageSpot::Validate";
    try {
        if (!oneContract && !fwdStarting) {
            if (!Maths::isPositive(initialSpot)) {
                string m("initial spot (" + Format::toString(initialSpot) + 
                         ") <= 0.0");
                throw ModelException(method, m);
            }
        }

        if (fwdStarting && oneContract) {
            throw ModelException(method,
                                 "fwd starting options must be notional based");

        }

        AssetUtil::assetCrossValidate(asset.get(),
                                      fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);

        if (fwdStarting && valueDate.isGreater(startDate)) {
            string m("option is fwd starting, but start "
                     "date (" + startDate.toString() + ") is not after today ("
                     + valueDate.toString() + ")");
                throw ModelException(method, m);
        }

        DateTime lodate;
        DateTime hidate;

        avgOut->getBoundingDates(lodate, hidate);

        if (hidate.isGreater(maturity)) {
            string m("last average out date (" + hidate.toString() + 
                     ") must be <= maturity (" +  maturity.toString() + ")");
                throw ModelException(method, m); 
        }
         
        if (fwdStarting && startDate.isGreater(lodate)) {
            string m("first average out date (" + lodate.toString() + 
                     ") must be >= start date (" + startDate.toString() + 
                     ") for a forward starting option");
            throw ModelException(method, m);
        }

        if (!(!asset)) {
            // check that underlying is not a future
            if (Future::TYPE->isInstance(asset.get())) {
                throw ModelException(method,
                                     "Average options on Futures are not allowed");
            }
            // check that underlying is not an fx asset
            if (FXAsset::TYPE->isInstance(asset.get())) {
                throw ModelException(method,
                                     "Options on FX assets are not allowed yet");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

    
/** returns the strike which is relevant for the vol interpolation */
double AverageSpot::getAdjustedStrike() const
{
    static const string method = "AverageSpot::getAdjustedStrike";
    try {
        double sumSoFar = avgOut->sumToDate(valueDate);
        double adjStrike;

        if (avgOut->countFutureSamples(valueDate) > 0) {
            double futureWeight = avgOut->futureWeight(valueDate);
        
            if (Maths::isZero(futureWeight)) {
                throw ModelException(method, 
                                     "total of future avg out weights "
                                     "must be > 0.0");
            }
            adjStrike = (strike - sumSoFar)/futureWeight;
        }
        else {
            adjStrike = strike - sumSoFar;
        }
            
        // if effective strike is -ve, use current spot
        if (!Maths::isPositive(adjStrike)) {
            adjStrike = asset->getSpot();
        }

        return adjStrike;
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* AverageSpot::createProduct(
    CClosedFormLN* model) const{
    return new AverageSpotClosedForm(this);
}

/** "constructor" for average spot */
Average* Average::makeAvgSpot(
    bool                        isCall,
    const DateTime&             maturity,
    double                      strike,
    const SampleList*           avgOut,
    const InstrumentSettlement* instSettle,
    const InstrumentSettlement* premiumSettle,
    const CAsset*               asset,
    string                      ccyTreatment,
    const YieldCurve*           discount,
    const DateTime&             valueDate,
    bool                        fwdStarting,
    const DateTime&             startDate,
    bool                        oneContract,
    double                      notional,    
    double                      initialSpot) {

    return new AverageSpot(isCall,
                           maturity,
                           strike,
                           avgOut,
                           instSettle,
                           premiumSettle,
                           asset,
                           ccyTreatment,
                           discount,
                           valueDate,
                           fwdStarting,
                           startDate,
                           oneContract,
                           notional,
                           initialSpot);
}

/** Price a spread, lowStrike and highStrike are passed by reference so the
    calling function can have the strikes converted to absolute values
    for fwdStarting options without recalculating the asset fwd value */
double Average::priceSpotSpread(const DateTime& valueDate,
                                const DateTime& startDate,
                                const DateTime& maturityDate,
                                bool  fwdStarting,
                                bool  isCall,
                                bool  oneContract,
                                double notional,
                                double initialSpot,
                                double& lowStrike,
                                double& highStrike,
                                const InstrumentSettlement* instSettle,
                                const Asset* asset,
                                const YieldCurve* discount,
                                const SampleList* avgOut)
{
    
    static const string method = "Average::priceSpotSpread";
    try
    {
        // if we're past settlement, it's just worth nothing
        if (valueDate.isGreaterOrEqual(instSettle->settles(maturityDate, asset)))
        {
            return 0;
        }
        
        DateTime start = fwdStarting ? startDate : valueDate;
        double lowK = lowStrike;
        double highK = highStrike;

        // if it's forward starting, convert percentage strike 
        // to absolute value based on spot at start date
        double fwdAtStart = 0.0;
        if (fwdStarting) {
            fwdAtStart = asset->fwdValue(start);
            lowK *= fwdAtStart;
            highK *= fwdAtStart;
        }

        // estimate value of average forward 
        double avgFwd = avgOut->futureSampleSum(asset, valueDate);
        
        // effective strike 
        double sumSoFar = avgOut->sumToDate(valueDate);
        
        lowK -= sumSoFar;
        highK -= sumSoFar;

        // choose how to interpolate the vols 
        CVolRequestLNSP volRequestLow(Average::volInterpSpot(valueDate,
                                                             startDate,
                                                             maturityDate,
                                                             fwdStarting,
                                                             lowStrike,
                                                             avgOut));
        CVolRequestLNSP volRequestHigh(Average::volInterpSpot(valueDate,
                                                              startDate,
                                                              maturityDate,
                                                              fwdStarting,
                                                              highStrike,
                                                              avgOut));
        
        // interpolate the vol using our LN request
        CVolProcessedBSSP  volBSlow(asset->getProcessedVol(volRequestLow.get()));
        CVolProcessedBSSP  volBShigh(asset->getProcessedVol(volRequestHigh.get()));
        
        // variance 
        double avgVarianceLow = avgOut->averageVariance(volBSlow.get(),
                                                        fwdStarting ? startDate : 
                                                        valueDate,
                                                        false);

        double avgVarianceHigh = avgOut->averageVariance(volBShigh.get(),
                                                         fwdStarting ? startDate : 
                                                         valueDate,
                                                         false);
        
        
        // discounting 
        double pv = instSettle->pv(maturityDate, discount, asset);

        double highValue = avgBlack(isCall,
                                    avgFwd,
                                    highK,
                                    pv,
                                    avgVarianceHigh);

        double lowValue = avgBlack(isCall,
                                   avgFwd,
                                   lowK,
                                   pv,
                                   avgVarianceLow);
        
        if (!oneContract) 
        {
            // handle notional 
            if (fwdStarting) {
                if ( Maths::isZero(fwdAtStart) )
                {
                    throw ModelException("AverageSpot::priceSpread",
                                         "Forward at start is 0.0. Infinite premium.");
                }
                lowValue *= notional/fwdAtStart;
                highValue *= notional/fwdAtStart;
            }
            else  {
                // handle position 
                lowValue *= notional/initialSpot;
                highValue *= notional/initialSpot;
            }
        }
        
        // Now combine the premiums
        double spreadPrice = isCall?
            lowValue - highValue : /* call spread */
            highValue - lowValue;  /* put spread */

        if (fwdStarting) 
        {
            lowStrike *= fwdAtStart;
            highStrike *= fwdAtStart;
        }

        return spreadPrice;
        
    }
    catch (exception& e)
    {
        throw ModelException(&e, method);
    }
}

double AverageSpot::calcAvgVariance( double strike ) const
{
    // choose how to interpolate the vol
    CVolRequestLNSP volRequest(volInterp(strike));

    // interpolate the vol using our LN request
    CVolProcessedBSSP  volBS(asset->getProcessedVol(volRequest.get()));

    // variance
    return avgOut->averageVariance(volBS.get(),
                                   fwdStarting ? startDate :
                                                 valueDate,
                                   false);
}

void AverageSpot::priceLN(
    Control* control,
    CResults* results,
    double strike,
    double avgVariance ) const
{
    static const string method = "AverageSpot::priceLN";
    try {
        // if we're past settlement, it's just worth nothing
        if (valueDate.isGreaterOrEqual(instSettle->settles(maturity, asset.get()))) {
            results->storePrice(0.0, discount->getCcy());
            return;
        }

        DateTime start = fwdStarting ? startDate : valueDate;
        double   k = strike;

        // if it's forward starting, convert percentage strike
        // to absolute value based on spot at start date
        double fwdAtStart = 0.0;
        if (fwdStarting) {
            fwdAtStart = asset->fwdValue(start);
            k *= fwdAtStart;
        }

        // estimate value of average forward
        double avgFwd = avgOut->futureSampleSum(asset.get(), valueDate );

        // effective strike
        double sumSoFar = avgOut->sumToDate(valueDate);

        k -= sumSoFar;

        // choose how to interpolate the vol
        CVolRequestLNSP volRequest(volInterp(strike));

        // interpolate the vol using our LN request
        CVolProcessedBSSP  volBS(asset->getProcessedVol(volRequest.get()));

        // discounting
        double pv = instSettle->pv(maturity, discount.get(), asset.get());

        double value = avgBlack(isCall,
                                avgFwd,
                                k,
                                pv,
                                avgVariance);

        if (!oneContract) {
            // handle notional
            if (fwdStarting) {
              if ( Maths::isZero(fwdAtStart) )
              {
                  throw ModelException("AverageSpot::priceLN",
                                       "Forward at start is 0.0. Infinite premium.");
              }
              value*= notional/fwdAtStart;
            }
            else  {
                // handle position
                value *= notional/initialSpot;
            }
        }

        results->storePrice(value, discount->getCcy());

        if (control && control->isPricing() ) {
            OutputRequest* request =
                control->requestsOutput(OutputRequest::IND_VOL);
            if (request) {
                indVol(volBS.get(), request, results);
            }

            recordExtras(control,
                         results,
                         asset.get(),
                         discount.get(),
                         valueDate,
                         maturity,
                         k,
                         avgFwd,
                         avgVariance,
                         value,
                         premiumSettle.get(),
                         instSettle.get(),
                         avgOut.get());
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void AverageSpot::priceLN(Control* control, CResults* results) const
{
    // if we're past settlement, it's just worth nothing
    if (valueDate.isGreaterOrEqual(instSettle->settles(maturity, asset.get()))) {
        results->storePrice(0.0, discount->getCcy());
        return;
    }

    priceLN( control, results, strike, calcAvgVariance( strike ) );
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool AverageSpot::avoidVegaMatrix(const IModel* model)
{
    /* this should possibly be false for local vol models etc.
       We may require the model object here */
    return false;
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP AverageSpot::getSensitiveStrikes(OutputNameConstSP outputName,
                                               const IModel*      model)
{
    static const string method = "AverageSpot::getSensitiveStrikes";
    DoubleArraySP sensStrikes;
    try {
        sensStrikes = DoubleArraySP(new DoubleArray(0));
 
        if (avoidVegaMatrix(model)) {
            throw ModelException(method, 
                                 "VEGA_MATRIX is not valid for this instrument");
        }

        if (avgOut->countFutureSamples(valueDate) > 0) {
            CVolRequestSP volRequest(volInterp(strike));
       
            SensitiveStrikeDescriptor sensStrikeDesc;
            sensStrikeDesc.forwardOnly = false;

            asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                       sensStrikeDesc, sensStrikes);
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
    return sensStrikes;
}

// resize delta shift for fwd starting
CSensControl* AverageSpot::AlterControl(
    const IModel*       modelParams,
    const CSensControl* sensControl) const
{
    static const string method = "AverageSpot::AlterControl";
    try {
        SensControlPerName* alteredControl = 0;

        if (Delta::TYPE->isInstance(sensControl)) {
            const Delta* delta = 
                dynamic_cast<const Delta*>((IObject*)sensControl);
            double theStrike = fwdStarting? strike: getAdjustedStrike();
            ShiftSizeCollectorSP shiftSizeVisitor(new ShiftSizeCollector(
                delta,
                theStrike,
                fwdStarting?
                ShiftSizeCollector::FWD_START_ADJUSTMENT:
                ShiftSizeCollector::SPOT_START_ADJUSTMENT));
            asset->accept(shiftSizeVisitor.get());

            if (Maths::isPositive(shiftSizeVisitor->getShiftSize())) {
                alteredControl = new Delta(shiftSizeVisitor->getShiftSize());
                alteredControl->
                    setMarketDataName(sensControl->getMarketDataName());
            }
        }
        return alteredControl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/*************************************   AverageRatio class  ****************************************/
/*************************************   AverageRatio class  ****************************************/
/*************************************   AverageRatio class  ****************************************/
class AverageRatio : public Average, 
                     public virtual CClosedFormLN::IIntoProduct,
                     public virtual Theta::Shift,
                     public virtual FDModel::IIntoProduct {
private:
    // additional fields specific to average ratio
    double       notional;
    SampleListSP avgIn;

    string       interpMethod; // FD: interpolation method
    int          dimAvgOutFD;  // FD: number of slices in average out dimension for FD
    int          dimAvgInFD;   // FD: number of slices in average in dimension for FD
    double       stateStdev;   // FD: number of standard deviations for average state variables, default=3.5
    bool         UseCtrlVar;   // FD: flag responsible for use of control variate technique 
    
public:
    static CClassConstSP const TYPE;

    AverageRatio();
    AverageRatio(bool                       isCall,
                const DateTime&             maturity,
                double                      strike,
                const SampleList*           avgOut,
                const InstrumentSettlement* instSettle,
                const InstrumentSettlement* premiumSettle,
                const CAsset*               asset,
                string                      ccyTreatment,
                const YieldCurve*           discount,
                const DateTime&             valueDate,
                const SampleList*           avgIn,
                double                      notional);

    // how we interpolate vol
    CVolRequestLN* volInterp(double strike) const;

    // indicative vol
    void indVol(CVolProcessedBS* vol, 
                OutputRequest*   request, 
                CResults*        results) const;

    virtual void Validate();
    
    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    /** Implementation of IMCIntoProduct interface */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    virtual bool sensShift(Theta* shift);

    virtual CSensControl* AlterControl(const IModel*       modelParams,
                                       const CSensControl* sensControl) const;
      
    /** returns ave-in sample list. returns 0 if does not exist.
        default returns 0;*/
    virtual SampleListSP getAveIn() const{
        return avgIn;
    }

private:
    friend class AverageRatioHelper;
    friend class AverageRatioClosedForm;

    friend class AverageRatioFD;         // FD: FD implementation of Average Spot 
    friend class AverageRatioOutState;   // FD: out-state variable
    friend class AverageRatioInState;    // FD: in-state variable

    void priceLN(Control* control, CResults* results) const;
};
/*************************************   AverageSpotClosedForm  ****************************************/
class AverageRatioClosedForm: public CClosedFormLN::IProduct{
private:
    const AverageRatio* avg; // a reference

public:
    AverageRatioClosedForm(const AverageRatio* avg): avg(avg){}

    void price(CClosedFormLN*  model,
               Control*        control, 
               CResults*       results) const{
        avg->priceLN(control, results);
    }
};
/************************ STATE VARIABLES *********************/
typedef TreeSliceLayer::StateSupport::InterpMethod AVG_INTERP;
/************************ averaging out *********************/
class AverageRatioFD;
class AverageRatioOutState : public TreeSliceLayer::StateSupport{
public:
    static double calcAvg(const AverageRatio* inst){
        double avg = inst->avgOut->averageToDate(inst->valueDate);
        return avg;
    }
    AverageRatioOutState(const AverageRatioFD* prod, const AverageRatio* inst, 
                         FDProductSP payoffIndex,   int dim, AVG_INTERP method) :
        StateSupport( "AverageRatioOutState", calcAvg(inst), method ),
        prod(prod),
        inst(inst),
        payoffIndex(payoffIndex),
        dim(dim),
        sampleCount(inst->avgOut->numDates(0)) {}

    /* state variable support using slice layers, TreeSliceLayer::StateSupport:: */
    // update state contents
    void update(int step, const TreeSlice& s);
    // set grid levels
    void setGrid( int step, bool init = false );
    // compute grid updates at a reset step. This is for each slice node level.
    virtual void transition(
        const vector< const TreeSlice * > & srcSlices,
        const vector< double > & gridIn,
        vector< double > & gridOut ) const;

    const AverageRatioFD*    prod;
    const AverageRatio*      inst;
    FDProductSP              payoffIndex;
    int                      dim;
    double                   avgSoFar;
    int                      firstResetStep;
    BoolArray                isResetStep;
    int                      sampleCount;
};
/************************ state operation methods *********************/
void AverageRatioOutState::setGrid( int step, bool init )
{
    if( init || step != firstResetStep )
    {
        // computing average range using stateStdev
        DateTimeArray sampleDates(inst->avgOut->getDates());
        DoubleArray values(inst->avgOut->getValues());
        DoubleArray weights(inst->avgOut->getWeights());
        sampleDates.resize(sampleCount);
        values.resize(sampleCount);
        weights.resize(sampleCount);

        SampleList samples(sampleDates, values, weights);
       
        // compute avg varince
        // computing average range using stateStdev
        ATMVolRequest volReq;
        CVolProcessedBSSP  volBS(inst->asset->getProcessedVol(&volReq));
        double avgVariance = samples.averageVariance(volBS.get(), inst->valueDate, true);

        double fwd = samples.expectedAverage(inst->asset.get(), inst->valueDate);
        double maxAvg = fwd*exp(-0.5*avgVariance + inst->stateStdev*sqrt(avgVariance));
        double minAvg = fwd*exp(-0.5*avgVariance - inst->stateStdev*sqrt(avgVariance));
      
        int n = dim; // *sqrt((double)sampleCount/inst->monitorDates.size());
        n = Maths::max(n, 2);
        // build a log array for average
        double increment = exp((log(maxAvg) - log(minAvg))/(n - 1.0));

        currGrid.resize(n);
        currGrid[0] = minAvg;
        for (int j = 1; j < n; ++j) {
            currGrid[j] = currGrid[j-1] * increment;
        }
    }
    else {
        // collapse grid to dim=1
        currGrid.resize(1); // !!! currGrid.resize(1, todayValue) does not work because resize down
        currGrid[0] = todayValue;
    }

    if( init )
        prevGrid = currGrid;

    populateGrid( step );
}
// compute grid transition through a reset/monitoring
void AverageRatioOutState::transition(
    const vector< const TreeSlice * > & srcSlices,
    const vector< double > & gridIn,
    vector< double > & gridOut ) const
{
    int nbGrid = gridOut.size();
    ASSERT( nbGrid == gridIn.size() );

    double s = srcSlices[ 0 ]->calc();

    for( int i = 0; i < nbGrid; ++i )
        gridOut[ i ] = ( gridIn[ i ] * ( sampleCount - 1 ) + s ) / sampleCount;
}
// interpolation etc.
void AverageRatioOutState::update(int step, const TreeSlice& s){
    if( isResetStep[step] ){
        setGrid( step );
        vector< const TreeSlice * > srcSlices( 1, &s );
        interpolate( step, srcSlices );
        --sampleCount;
    }
}
typedef refCountPtr<AverageRatioOutState > AverageRatioOutStateSP;
/************************ averaging in *********************/
class AverageRatioInState : public TreeSliceLayer::StateSupport{
public:
    static double calcAvg(const AverageRatio* inst){
        double avg = inst->avgIn->averageToDate(inst->valueDate);
        return avg;
    }
    AverageRatioInState(const AverageRatioFD* prod, const AverageRatio* inst, 
                     FDProductSP payoffIndex,   int dim, AVG_INTERP method) :
        StateSupport( "AverageRatioInState", calcAvg(inst), method ),
        prod(prod),
        inst(inst),
        payoffIndex(payoffIndex),
        dim(dim),
        sampleCount(inst->avgIn->numDates(0)) {}

    /* state variable support using slice layers, TreeSliceLayer::StateSupport:: */
    // update state contents
    void update(int step, const TreeSlice& s);
    // set grid levels
    void setGrid( int step, bool init = false );
    // compute grid updates at a reset step. This is for each slice node level.
    virtual void transition(
        const vector< const TreeSlice * > & srcSlices,
        const vector< double > & gridIn,
        vector< double > & gridOut ) const;

    const AverageRatioFD*   prod;
    const AverageRatio*     inst;
    FDProductSP             payoffIndex;
    int                     dim;
    double                  avgSoFar;
    int                     firstResetStep;
    BoolArray               isResetStep;
    int                     sampleCount;
};
/************************ state operation methods *********************/
void AverageRatioInState::setGrid( int step, bool init )
{
    if( init || step != firstResetStep )
    {
        // computing average range using stateStdev
        DateTimeArray sampleDates(inst->avgIn->getDates());
        DoubleArray values(inst->avgIn->getValues());
        DoubleArray weights(inst->avgIn->getWeights());
        sampleDates.resize(sampleCount);
        values.resize(sampleCount);
        weights.resize(sampleCount);

        SampleList samples(sampleDates, values, weights);
       
        // compute avg varince
        // computing average range using stateStdev
        ATMVolRequest volReq;
        CVolProcessedBSSP  volBS(inst->asset->getProcessedVol(&volReq));
        double avgVariance = samples.averageVariance(volBS.get(), inst->valueDate, true);

        double fwd = samples.expectedAverage(inst->asset.get(), inst->valueDate);
        double maxAvg = fwd*exp(-0.5*avgVariance + inst->stateStdev*sqrt(avgVariance));
        double minAvg = fwd*exp(-0.5*avgVariance - inst->stateStdev*sqrt(avgVariance));
      
        int n = dim; // *sqrt((double)sampleCount/inst->monitorDates.size());
        n = Maths::max(n, 2);
        // build a log array for average
        double increment = exp((log(maxAvg) - log(minAvg))/(n - 1.0));

        currGrid.resize(n);
        currGrid[0] = minAvg;
        for (int j = 1; j < n; ++j) {
            currGrid[j] = currGrid[j-1] * increment;
        }
    }
    else {
        // collapse grid to dim=1
        currGrid.resize(1); // !!! currGrid.resize(1, todayValue) does not work because resize down
        currGrid[0] = todayValue;
    }

    if( init )
        prevGrid = currGrid;

    populateGrid( step );
}
// compute grid transition through a reset/monitoring
void AverageRatioInState::transition(
    const vector< const TreeSlice * > & srcSlices,
    const vector< double > & gridIn,
    vector< double > & gridOut ) const
{
    int nbGrid = gridOut.size();
    ASSERT( nbGrid == gridIn.size() );

    double s = srcSlices[ 0 ]->calc();

    for( int i = 0; i < nbGrid; ++i )
        gridOut[ i ] = ( gridIn[ i ] * ( sampleCount - 1 ) + s ) / sampleCount;
}
// interpolation etc.
void AverageRatioInState::update(int step, const TreeSlice& s){
    if( isResetStep[step] ){
        setGrid( step );
        vector< const TreeSlice * > srcSlices( 1, &s );
        interpolate( step, srcSlices );
        --sampleCount;
    }
}
typedef refCountPtr<AverageRatioInState > AverageRatioInStateSP;
/************************ end state variables operation *********************/
/************************      FDProduct class          *********************/
class AverageRatioFD : public LatticeProdEDR, 
                       virtual public IFDProductLN {
private:
    friend class AverageRatioOutState;
    friend class AverageRatioInState;

    const AverageRatio*        inst;
    //vector<bool>              stepExercise;
    BoolArray                 stepExercise;

    int                       avgInRemaining; 
    int                       avgOutRemaining;
    double                    avgOutSumSoFar;
    
    // price slices:
    // the price of the option, European, American or Bermudan
    TreeSliceSP value;
    // the price of the European option used for control variate technique
    TreeSliceSP valueCtrlVariate;
    
    // shorter for convenience
    AVG_INTERP interpMethod;
public: 
    AverageRatioOutStateSP avgOutState;
    AverageRatioInStateSP  avgInState;

    /** constructor */
    AverageRatioFD(const AverageRatio* inst, FDModel* model) : 
        LatticeProdEDR(model), 
        inst(inst)
    {
        payoffIndex = model->createProduct(IProdCreatorSP(new
            IndexSpecEQ(inst->asset.getName(), inst->asset, inst->ccyTreatment)));

        avgInRemaining  = inst->avgIn->countFutureSamples( inst->valueDate );
        avgOutRemaining = inst->avgOut->countFutureSamples( inst->valueDate );
        avgOutSumSoFar  = inst->avgOut->sumToDate( inst->valueDate );

        if (inst->interpMethod == "QUADRATIC")
            interpMethod = TreeSliceLayer::StateSupport::INTERP_QUADRATIC;
        else if (inst->interpMethod == "LINEAR")
            interpMethod = TreeSliceLayer::StateSupport::INTERP_LINEAR;
        else
            throw ModelException("AverageRatioFD: ", "unknown interpMethod: " + inst->interpMethod);

        avgOutState = AverageRatioOutStateSP(new AverageRatioOutState(this, inst, payoffIndex, inst->dimAvgOutFD, interpMethod));
        avgInState = AverageRatioInStateSP(new AverageRatioInState(this, inst, payoffIndex, inst->dimAvgInFD, interpMethod));
    }

    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type);

     /** local method, product payoff method at maturity */
    void prod_BWD_T(const TreeSlice & s, int step);

    /** local method, this is payoff boundary condition, for KO, early exercise etc */
    void prod_BWD(const TreeSlice & s, int step);
    
    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(Control*  control) const;

    /** set flags for averaging time steps,
        return the first step of the monitorDates after valueDate */
    int setStepFlags(BoolArray& stepFlags, const DateTimeArray& stepDates, const DateTimeArray& monitorDates);
    
    /** initialising and setting product variables this is called per pricing call before each pricing */
    virtual void initProd();

    /** output prices and any additional results */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results);

    /** scale for notional */ 
    double scalePremium(const double& fairValue);

    /** post price */
    double postPrice(double euroPrice, double exer_premium, YieldCurveConstSP disc);
    
    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const {return inst->ccyTreatment;}

    /** Returns start date */
    virtual DateTime getStartDate() const {
        return inst->valueDate;
    }

    /** for the LogNormal model */
    CVolRequestLNSP getVolInterp(int iAsset) const {
        // get strike and maturity date from instrument
        DateTime matDate = inst->avgOut->getLastDate();
        
        double volStrike = inst->strike;

        DateTime imntStartDate = getStartDate();
        CVolRequestLNSP   reqarr;         
        reqarr = CVolRequestLNSP(new LinearStrikeVolRequest(volStrike, imntStartDate, matDate, false));        
        return reqarr;
    }
};

/** create a fd payoff product */
FDProductSP AverageRatio::createProduct(FDModel* model) const {
    return FDProductSP(new AverageRatioFD(this, model));
}
/* Update function*/
void AverageRatioFD::update(int& step, FDProduct::UpdateType type)
{
    const TreeSlice & s = payoffIndex->getValue(step);

    if( type == FDProduct::BWD_T )
    {
        avgOutState->setGrid( step, true );
        avgInState->setGrid( step, true );
        
        prod_BWD_T(s, step);
    }
    else if( type == FDProduct::BWD )
        prod_BWD(s, step);

    avgOutState->update(step, s);
    avgInState->update(step, s);
}
/* local method, product payoff method at maturity */
void AverageRatioFD::prod_BWD_T(const TreeSlice & s, int step)
{
    static const string method = "AverageRatioFD::prod_BWD_T";
    try 
    {
        double settlementPV = inst->instSettle->pvAdjust(model->getDate(step), 
                                                         inst->discount.get(), 
                                                         inst->asset.get());
        if( ! inst->isCall )
            settlementPV = -settlementPV;

        const TreeSlice & outGrid = avgOutState->getGridSlice();
        const TreeSlice & inGrid = avgInState->getGridSlice();
        
        if ( avgInRemaining > 0 ) {
            *value = smax( 0., settlementPV * ( outGrid / inGrid - inst->strike ) );
        }
        else{
            if ( avgOutRemaining > 0 ) {
                //*value = smax( 0., settlementPV * ( outGrid / inGrid - ( inst->strike - avgOutSumSoFar ) ) );
                *value = inGrid;
            }
            else{
                *value = smax( 0., settlementPV * ( outGrid / inGrid -  inst->strike ) );
            }
        }
        if ( inst->UseCtrlVar )
                *valueCtrlVariate = *value;
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}
/* local method, product payoff method before maturity */
void AverageRatioFD::prod_BWD(const TreeSlice & s, int step)
{
    static const string method = "AverageRatioFD::prod_BWD";
    try 
    {
        double settlementPV = inst->instSettle->pvAdjust(model->getDate(step), 
                                                         inst->discount.get(), 
                                                         inst->asset.get());
        if( ! inst->isCall )
            settlementPV = -settlementPV;

        // computes exercise price for american options
        if( !stepExercise[step] )
            return;

        const TreeSlice & outGrid = avgOutState->getGridSlice();
        const TreeSlice & inGrid = avgInState->getGridSlice();
        
        if ( avgInRemaining > 0 ) {
            *value = smax( *value, settlementPV * ( outGrid / inGrid - inst->strike ) );
        }
        else{
            if ( avgOutRemaining > 0 ) {
                *value = smax( *value, settlementPV * ( outGrid / inGrid - ( inst->strike - avgOutSumSoFar ) ) );
            }
            else{
                *value = smax( *value, settlementPV * ( outGrid / inGrid -  inst->strike * inGrid) );
            }
        }
    } 
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}
/** initialise tree1f - allow product customisation 
    must not init product variables here, use initProd() instead */
void AverageRatioFD::init(CControl* control) const
{
    static const string method = "AverageRatioFD::init()";
    try 
    {
        if (inst->dimAvgOutFD < 2) 
        {
            throw ModelException(method,"dimAvgOutFD must ba at least 2 when pricing with tree");
        }
        // customize tree parameters here and set up the tree
        DateTimeArray segDates;
        segDates.resize(2);
        segDates[0] = getStartDate();
        segDates[1] = inst->exerciseSchedule->lastDate();
        IntArray density(1,1);

        // all averaging dates are copied to critical dates
        model->addCritDates(inst->avgOut->getDates()); 
        model->addCritDates(inst->avgIn->getDates()); 
        // all exercise dates are copied to critical dates
        model->addCritDates(inst->exerciseSchedule->getDates());
        // prepare model set up
        model->initSegments(segDates, density);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/** set flags for averaging time steps,
    return the first step of the monitorDates after valueDate */
int AverageRatioFD::setStepFlags(BoolArray& stepFlags, const DateTimeArray& stepDates, const DateTimeArray& monitorDates) {

    static const string method = "AverageRatioFD::setStepAverage";
    try {
        int iMonitor = 0;
        while (monitorDates[iMonitor] < stepDates[0]) {
            iMonitor++;
        }

        int firstStep = -1;

        for (int iStep = 0; iStep < stepDates.size(); iStep++) {
            if (stepDates[iStep].equals(monitorDates[iMonitor], true)) {
                stepFlags[iStep] = true;
                if( firstStep < 0 ) firstStep = iStep;
                iMonitor++;
            }
            else {
                stepFlags[iStep] = false;
            }
        }
        return firstStep;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/** initialising and setting product variables
    this is called per pricing call before tree sweep call (after InitTree)  */
void AverageRatioFD::initProd()
{
    static const string method = "AverageRatioFD::initProd()";
    try 
    {
        int lastStep = model->getLastStep();
        stepExercise.resize(lastStep + 1);
       
        /*AssetUtil::setStepExercise(stepExercise,
                                   model->getDates(),
                                   inst->exerciseSchedule,
                                   inst->isAmerican,
                                   inst->asset.getSP());*/

        setStepFlags(stepExercise, model->getDates(), inst->exerciseSchedule->getDates());

        model->registerStateSupport( avgOutState.get() );
        model->registerStateSupport( avgInState.get() );

        startDEV( value = model->createSlice() ); 
        if ( inst->UseCtrlVar )
            startDEV( valueCtrlVariate = model->createSlice() ); 

        avgOutState->isResetStep.resize(lastStep + 1);
        avgOutState->firstResetStep =
            setStepFlags(avgOutState->isResetStep, model->getDates(), inst->avgOut->getDates());
        if(model->getDate(0) == inst->valueDate) // value date is already counted as history
            avgOutState->isResetStep[0] = false;

        avgInState->isResetStep.resize(lastStep + 1);
        avgInState->firstResetStep =
            setStepFlags(avgInState->isResetStep, model->getDates(), inst->avgIn->getDates());
        if(model->getDate(0) == inst->valueDate) // value date is already counted as history
            avgInState->isResetStep[0] = false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/** output results */
void AverageRatioFD::recordOutput(Control* control, YieldCurveConstSP disc, Results* results) 
{
    double scaledPrice;
    double price = model->getPrice0( *value );
    if ( inst->UseCtrlVar ){
        double priceCtrlVariate = model->getPrice0( *valueCtrlVariate );
        scaledPrice = postPrice( price, price - priceCtrlVariate, disc);
    }
    else
        scaledPrice = scalePremium( price );

    results->storePrice( scaledPrice, disc->getCcy() );

    // take care of additional outputs
    if (control && control->isPricing()) {

        DateTime matDate = inst->avgOut->getLastDate();

        OutputRequest* request = NULL;
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTime paymentDate = inst->instSettle->settles(matDate, inst->asset.get());
            DateTimeArray date(1, paymentDate);
            OutputRequestUtil::recordPaymentDates(control,results,&date); 
        }
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            if (inst->valueDate.isGreater(matDate)) {
                DateTime paymentDate = inst->instSettle->settles(matDate, inst->asset.get());
                double pv = disc->pv(paymentDate);
                CashFlow cf(paymentDate, scaledPrice/pv);
                CashFlowArray cfl(1, cf);
                OutputRequestUtil::recordKnownCashflows(control, results, disc->getCcy(), &cfl);   
            }
        }
        double        indVol;
        // calculate indicative vol
        if (matDate.isGreater(inst->valueDate)) {
            DateTime imntStartDate = getStartDate();

            // get vol request
            //double volStrike  = strike;
            double volStrike  = inst->strike;
            LinearStrikeVolRequest volRequest(volStrike, imntStartDate, matDate, false);

            try {
                // interpolate the vol
                CVolProcessedSP  vol(inst->asset->getProcessedVol(&volRequest));
                // cast to the BS vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                // this should never happen if our get market data has worked properly
                if (!volBS) {
                    throw ModelException("AverageRatioFD::recordOutput", "No Black Scholes Vol");
                }
                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
            }
            catch (exception& ) {
                indVol = 0.0;
            }
        }
        else {
            indVol = 0.0;
        }
    }
}
/* premium scaling */
double AverageRatioFD::scalePremium(const double& fairValue ) {
    return (fairValue*inst->notional);
}
// postPrice process
double AverageRatioFD::postPrice(double euroPrice, double exer_premium, YieldCurveConstSP disc)
{
    double price;

    AverageRatioClosedForm closedForm(inst);
    CClosedFormLN model;
    Control ctrl;
    CResults result;
    closedForm.price(&model, &ctrl, &result);
    double closedFormPrice = result.retrievePrice();

    price = scalePremium( exer_premium ) + closedFormPrice;
    
    return price;
}
/************************ end of FDProduct class *********************/
/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* AverageRatio::createProduct(const MonteCarlo* model) const {
    // SampleList implements RefLevel and IPastValues. Since avgIn holds
    // the past values we don't need to populate the past values with them.
    IRefLevelConstSP refLevel(IRefLevelConstSP::attachToRef(avgIn.get())); 
    IPastValuesConstSP pastValues(
        IPastValuesConstSP::attachToRef(avgOut.get()));
    IMCProduct*  prod = new MC(this,
                              refLevel,
                              avgIn->getFirstDate(),
                              pastValues,
                              true, // use ratio
                              strike, notional);
    return prod;
}
    
AverageRatio::AverageRatio() : Average(TYPE),dimAvgOutFD(151),dimAvgInFD(151),
                               stateStdev(3.5),interpMethod("QUADRATIC"), UseCtrlVar(false) {}

void AverageRatio:: Validate() {
    static const string method = "AverageRatio::Validate";
    try {
        AssetUtil::assetCrossValidate(asset.get(),
                                      false,
                                      valueDate,
                                      valueDate,
                                      discount,
                                      this);

        DateTime loout;
        DateTime hiout;
        DateTime loin;
        DateTime hiin;

        avgOut->getBoundingDates(loout, hiout);
        avgIn->getBoundingDates(loin, hiin);

        if (hiout.isGreater(maturity)) {
                throw ModelException(method,
                                     "last average out date (" + 
                                     hiout.toString() + 
                                     ") must be <= maturity (" + 
                                     maturity.toString() + ")"); 
        }

        if (hiin.isGreater(loout)) {
                throw ModelException(method,
                                     "last average in date (" + 
                                     hiin.toString() + 
                                     ") must be < first average out date (" + 
                                     loout.toString() + ")"); 
        }
         
        // ensure no-one's trying to do anything clever with the sample weights
        if (!avgOut->weightsSumToOne()) {
                throw ModelException(method,
                                     "average out weights do not sum to one");
        }
        if (!avgIn->weightsSumToOne()) {
                throw ModelException(method,
                                     "average in weights do not sum to one");
        }

        if (!(!asset)) {
            // check that underlying is not a future
            if (Future::TYPE->isInstance(asset.get())) {
                throw ModelException(method,
                                     "Average options on Futures are not allowed");
            }
            // check that underlying is not an fx asset
            if (FXAsset::TYPE->isInstance(asset.get())) {
                throw ModelException(method,
                                     "Options on FX assets are not allowed yet");
            }

            // or any kind of basket
            if (XCB::TYPE->isInstance(asset.get()) ||
                UnitXCBWithVol::TYPE->isInstance(asset.get()) ||
                PercXCBWithVol::TYPE->isInstance(asset.get())) {
                throw ModelException(method,
                                     "avg ratio on basket not allowed");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
}

bool AverageRatio::sensShift(Theta* shift) {
    static const string method = "AverageRatio::sensShift";
    try {
        // calculate the new date
        DateTime newDate = shift->rollDate(valueDate);
        bool useSpot = shift->useAssetFwds() ? false : true;

        // fill in any samples 
        avgOut->roll(asset.get(), valueDate, newDate, useSpot);
        avgIn->roll(asset.get(), valueDate, newDate, useSpot);
                           
        // roll today 
        valueDate = newDate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
    return true; // our components have theta type sensitivity
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* AverageRatio::createProduct(
    CClosedFormLN* model) const{
    return new AverageRatioClosedForm(this);
}

CVolRequestLN* AverageRatio::volInterp(double strike) const {
    static const string method = "AverageRatio::volInterp";
    try {
        CVolRequestLN* interp = 0;

        DateTime avgInDate  = avgIn->expectedEndDate(valueDate);
        DateTime avgOutDate = avgOut->expectedEndDate(valueDate);

        // if all average-out samples are known, nothing to do 
        if (valueDate.isGreaterOrEqual(avgOutDate)) {
            interp = new ATMVolRequest();
        }
        else {
            double avgOutFutureWeight = avgOut->futureWeight(valueDate);
        
            if (Maths::isZero(avgOutFutureWeight)) {
                throw ModelException(method, 
                                     "total of future avg out weights "
                                     "must be > 0.0");
            }

            double avgOutSumSoFar = avgOut->sumToDate(valueDate);

            double avgInExpectedSum = avgIn->expectedAverage(asset.get(),
                                                             valueDate);

            double effectiveStrike = (strike*avgInExpectedSum-avgOutSumSoFar)/
                                      avgOutFutureWeight;

            // if all average-in samples are set, this is all we need to do 
            if (valueDate.isGreaterOrEqual(avgInDate)) {
                interp = new LinearStrikeVolRequest(effectiveStrike, 
                                                    valueDate, 
                                                    maturity,
                                                    false);

            }
            else {
                double fwdAtStrikeDate = asset->fwdValue(avgInDate);
                double spot = asset->getSpot();

                double interpLevel = spot * effectiveStrike/fwdAtStrikeDate;

                double volSpread;
                if (avgOutDate == avgInDate)
                {
                    volSpread = 0.0;
                }
                else
                {
                    DateTime effectiveDate = valueDate.add(avgOutDate.subtract(avgInDate));

                    ATMVolRequestSP atm(new ATMVolRequest());
                
                    LinearStrikeVolRequestSP k(new LinearStrikeVolRequest(interpLevel, 
                                                                          valueDate, 
                                                                          maturity,
                                                                          false));

                    CVolProcessedBSSP volATMBS(asset->getProcessedVol(atm.get()));
                    CVolProcessedBSSP volStrikeBS(asset->getProcessedVol(k.get()));
               
                    double atmVol    = volATMBS->CalcVol(valueDate, effectiveDate);
                    double strikeVol = volStrikeBS->CalcVol(valueDate, effectiveDate);

                    volSpread = strikeVol - atmVol;
                }

                interp = new LinearStrikeSpreadVolRequest(interpLevel/spot,
                                                          valueDate,
                                                          maturity,
                                                          volSpread);
            }
        }
        return interp;
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void AverageRatio::indVol(CVolProcessedBS* vol, 
                          OutputRequest*   request,
                          CResults*        results) const {
    static const string method = "AverageRatio::indVol";
    // do 'indicative vol' - don't fail the whole pricing if this fails 
    try {
        double ivol;

        if (avgOut->countFutureSamples(valueDate) > 0) {
            DateTime endDate = avgOut->expectedEndDate(valueDate);
            DateTime inDate  = avgIn->expectedEndDate(valueDate);

            DateTime start = inDate.isGreater(valueDate) ? inDate : valueDate;
            ivol = vol->CalcVol(start, endDate);
        }
        else {
            // all samples set
            ivol = 0.0;
        }
        results->storeRequestResult(request, ivol);
    }
    catch (exception& ) {
    }
}


void AverageRatio::priceLN(Control* control, CResults* results) const{
    static const string method = "AverageRatio::priceLN";
    try {
        // if we're past settlement, it's just worth nothing
        if (valueDate.isGreaterOrEqual(instSettle->settles(maturity, asset.get()))) {
            results->storePrice(0.0, discount->getCcy());
            return;
        }

        // estimate value of average forward
        double avgOutFutureSum = avgOut->futureSampleSum(asset.get(),
                                                         valueDate );

        double avgOutSumSoFar = avgOut->sumToDate(valueDate);

        double avgInFutureSum = avgIn->futureSampleSum(asset.get(), valueDate);
        double avgInSumSoFar  = avgIn->sumToDate(valueDate);

        // how to interpolate the vol 
        CVolRequestLNSP volRequest = CVolRequestLNSP(volInterp(strike));

        // interpolate the vol
        CVolProcessedBSSP  volBS(asset->getProcessedVol(volRequest.get()));

        // variance
        double avgOutVariance = avgOut->averageVariance(volBS.get(),
                                                        valueDate,
                                                        false);

        double avgInVariance  = avgIn->averageVariance(volBS.get(),
                                                       valueDate,
                                                       true);

        // covariance
        double covariance = avgOut->averageCovariance(avgIn.get(),
                                                      volBS.get(),
                                                      valueDate);
        // discounting 
        double pv = instSettle->pv(maturity, discount.get(), asset.get());

        int avgInRemaining  = avgIn->countFutureSamples(valueDate);
        int avgOutRemaining = avgOut->countFutureSamples(valueDate);

        double avgOutLevel = avgOutSumSoFar + avgOutFutureSum;
        double avgInLevel  = avgInSumSoFar  + avgInFutureSum;

        double value;
        double effectiveStrike;
        double avgFwd;
        double avgVariance;

        if (avgInRemaining > 0) {       
            effectiveStrike = strike * (avgInSumSoFar + avgInFutureSum);
            avgFwd = avgOutLevel;

            avgVariance = avgOutVariance + avgInVariance - 2 * covariance;

            double sqrtVar = sqrt(avgVariance);

            double fwd = avgOutLevel/avgInLevel * 
                         exp(0.5*avgVariance - 0.5*avgOutVariance + 
                             0.5*avgInVariance);

            double d1  = (log(fwd/strike) + 0.5*avgVariance)/sqrtVar;
            double d2  = d1 - sqrtVar;

            double nd1 = N1(d1 * (isCall ? 1.0 : -1.0));
            double nd2 = N1(d2 * (isCall ? 1.0 : -1.0));

            value = pv * (fwd*nd1 - strike*nd2) * (isCall ? 1.0 : -1.0);    
       }
        else {
            double k;

            if (avgOutRemaining > 0) {
                k = strike * avgInLevel - avgOutSumSoFar;
                avgFwd  = avgOutFutureSum;
            }
            else {
                k = strike * avgInLevel;
                avgFwd = avgOutLevel;
            }
        
            effectiveStrike = k; // for reporting

            avgVariance = avgOutVariance;

            value = avgBlack(isCall,
                             avgFwd,
                             k,
                             pv,
                             avgVariance);

            value /= avgInLevel;
        }
      
        // handle notional        
        value *= notional;
               
        results->storePrice(value, discount->getCcy());

        if (control && control->isPricing()) {
            OutputRequest* request;
            if ((request = control->requestsOutput(OutputRequest::IND_VOL))) {
                indVol(volBS.get(), request, results);
            }

            recordExtras(control,
                         results,
                         asset.get(),
                         discount.get(),
                         valueDate,
                         maturity,
                         effectiveStrike,
                         avgFwd,
                         avgVariance,
                         value,
                         premiumSettle.get(),
                         instSettle.get(),
                         avgOut.get());
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

bool AverageRatio::avoidVegaMatrix(const IModel* model)
{
    /* this should possibly be false for local vol models etc.
       We may require the model object here */
    return false;
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP AverageRatio::getSensitiveStrikes(OutputNameConstSP outputName,
                                                const IModel*      model)
{
    static const string method = "AverageRatio::getSensitiveStrikes";
    DoubleArraySP sensStrikes;
    try {
        sensStrikes = DoubleArraySP(new DoubleArray(0));
 
        if (avoidVegaMatrix(model)) {
            throw ModelException(method, 
                                 "VEGA_MATRIX is not valid for this instrument");
        }

        if (avgOut->countFutureSamples(valueDate) > 0 ||
            avgIn->countFutureSamples(valueDate)  > 0) {
            CVolRequestSP volRequest(volInterp(strike));
       
            SensitiveStrikeDescriptor sensStrikeDesc;
            sensStrikeDesc.forwardOnly = false;

            asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                       sensStrikeDesc, sensStrikes);
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
    return sensStrikes;
}

// resize delta shift for fwd starting
CSensControl* AverageRatio::AlterControl(
    const IModel*       modelParams,
    const CSensControl* sensControl) const
{
    static const string method = "AverageRatio::AlterControl";
    try {
        SensControlPerName* alteredControl = 0;
        // is it effectively forward starting ?
        if ((avgIn->countFutureSamples(valueDate) > 0) && 
            Delta::TYPE->isInstance(sensControl)) {
            const Delta* delta = 
                dynamic_cast<const Delta*>((IObject*)sensControl);
            ShiftSizeCollectorSP shiftSizeVisitor(new ShiftSizeCollector(
                delta,
                strike,
                ShiftSizeCollector::FWD_START_ADJUSTMENT));
            asset->accept(shiftSizeVisitor.get());

            if (Maths::isPositive(shiftSizeVisitor->getShiftSize())) {
                alteredControl = new Delta(shiftSizeVisitor->getShiftSize());
                alteredControl->
                    setMarketDataName(sensControl->getMarketDataName());
            }
        }

        return alteredControl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// AVERAGE HYBRID

class AverageHybrid : public Average, 
                      public virtual CClosedFormLN::IIntoProduct,
                      public virtual Theta::Shift {
private:
    // additional fields specific to average hybrid
    SampleListSP avgIn;

public:
    static CClassConstSP const TYPE;

    AverageHybrid();

    // how we interpolate vol
    CVolRequestLN* volInterp(double strike) const;

    // indicative vol
    void indVol(CVolProcessedBS* vol, 
                OutputRequest*   request, 
                CResults*        results) const;

    virtual void Validate();
    
    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** Implementation of IMCIntoProduct interface */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);
  
    virtual bool sensShift(Theta* shift);

    virtual CSensControl* AlterControl(
        const IModel*       modelParams,
        const CSensControl* sensControl) const;
     
    /** returns ave-in sample list. returns 0 if does not exist.
        default returns 0;*/
    virtual SampleListSP getAveIn() const
    {
        return avgIn;
    }

private:
    friend class AverageHybridHelper;
    friend class AverageHybridClosedForm;
    void priceLN(Control* control, CResults* results) const;
};

class AverageHybridClosedForm: public CClosedFormLN::IProduct{
private:
    const AverageHybrid* avg; // a reference

public:
    AverageHybridClosedForm(const AverageHybrid* avg): avg(avg){}

    void price(CClosedFormLN*  model,
               Control*        control, 
               CResults*       results) const{
        avg->priceLN(control, results);
    }
};
    
AverageHybrid::AverageHybrid() : Average(TYPE) {}

void AverageHybrid:: Validate() {
    static const string method = "AverageHybrid::Validate";
    try {
        AssetUtil::assetCrossValidate(asset.get(),
                                      false,
                                      valueDate,
                                      valueDate,
                                      discount,
                                      this);

        DateTime loout;
        DateTime hiout;
        DateTime loin;
        DateTime hiin;

        avgOut->getBoundingDates(loout, hiout);
        avgIn->getBoundingDates(loin, hiin);

        if (hiout.isGreater(maturity)) {
                throw ModelException(method,
                                     "last average out date (" + 
                                     hiout.toString() + 
                                     ") must be <= maturity (" + 
                                     maturity.toString() + ")"); 
        }

        if (hiin.isGreater(loout)) {
                throw ModelException(method,
                                     "last average in date (" + 
                                     hiin.toString() + 
                                     ") must be < first average out date (" + 
                                     loout.toString() + ")"); 
        }
        
        // ensure no-one's trying to do anything clever with the sample weights
        if (!avgOut->weightsSumToOne()) {
                throw ModelException(method,
                                     "average out weights do not sum to one");
        }
        if (!avgIn->weightsSumToOne()) {
                throw ModelException(method,
                                     "average in weights do not sum to one");
        }

        if (!(!asset)) {
            // check that underlying is not a future
            if (Future::TYPE->isInstance(asset.get())) {
                throw ModelException(method,
                                     "Average options on Futures are not allowed");
            }
            // check that underlying is not an fx asset
            if (FXAsset::TYPE->isInstance(asset.get())) {
                throw ModelException(method,
                                     "Options on FX assets are not allowed yet");
            }

            // or any kind of basket
            if (XCB::TYPE->isInstance(asset.get()) ||
                UnitXCBWithVol::TYPE->isInstance(asset.get()) ||
                PercXCBWithVol::TYPE->isInstance(asset.get())) {
                throw ModelException(method,
                                     "avg hybrid on basket not allowed");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
}

bool AverageHybrid::sensShift(Theta* shift) {
    static const string method = "AverageHybrid::sensShift";
    try {
        // calculate the new date
        DateTime newDate = shift->rollDate(valueDate);
        bool useSpot = shift->useAssetFwds() ? false : true;

        // fill in any samples 
        avgOut->roll(asset.get(), valueDate, newDate, useSpot);
        avgIn->roll(asset.get(), valueDate, newDate, useSpot);
                           
        // roll today 
        valueDate = newDate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
    return true; // our components have theta type sensitivity
}
/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* AverageHybrid::createProduct(
    CClosedFormLN* model) const{
    return new AverageHybridClosedForm(this);
}

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* AverageHybrid::createProduct(const MonteCarlo* model) const {
    // SampleList implements RefLevel and IPastValues. Since avgIn holds
    // the past values we don't need to populate the past values with them.
    IRefLevelConstSP refLevel(IRefLevelConstSP::attachToRef(avgIn.get())); 
    IPastValuesConstSP pastValues(
        IPastValuesConstSP::attachToRef(avgOut.get()));
    IMCProduct*  prod = new MC(this,
                              refLevel,
                              avgIn->getFirstDate(),
                              pastValues,
                              false, // don't use ratio
                              strike, 1.0);
    return prod;
}

CVolRequestLN* AverageHybrid::volInterp(double strike) const {
    static const string method = "AverageHybrid::volInterp";
    try {
        CVolRequestLN* interp = 0;

        DateTime avgInDate  = avgIn->expectedEndDate(valueDate);
        DateTime avgOutDate = avgOut->expectedEndDate(valueDate);

        // if all average-out samples are known, nothing to do 
        if (valueDate.isGreaterOrEqual(avgOutDate)) {
            interp = new ATMVolRequest();
        }
        else {
            double avgOutFutureWeight = avgOut->futureWeight(valueDate);
        
            if (Maths::isZero(avgOutFutureWeight)) {
                throw ModelException(method, 
                                     "total of future avg out weights "
                                     "must be > 0.0");
            }

            double avgOutSumSoFar = avgOut->sumToDate(valueDate);

            double avgInExpectedSum = avgIn->expectedAverage(asset.get(),
                                                             valueDate);

            double effectiveStrike = (strike*avgInExpectedSum-avgOutSumSoFar)/
                                      avgOutFutureWeight;

            // if all average-in samples are set, this is all we need to do 
            if (valueDate.isGreaterOrEqual(avgInDate)) {
                interp = new LinearStrikeVolRequest(effectiveStrike, 
                                                    valueDate, 
                                                    maturity,
                                                    false);

            }
            else {
                double fwdAtStrikeDate = asset->fwdValue(avgInDate);
                double spot = asset->getSpot();

                double interpLevel = spot * effectiveStrike/fwdAtStrikeDate;

                double volSpread;
                if (avgOutDate == avgInDate)
                {
                    volSpread = 0.0;
                }
                else
                {
                    DateTime effectiveDate = valueDate.add(avgOutDate.subtract(avgInDate));

                    ATMVolRequestSP atm(new ATMVolRequest());
                
                    LinearStrikeVolRequestSP k(new LinearStrikeVolRequest(interpLevel, 
                                                                          valueDate, 
                                                                          maturity,
                                                                          false));

                    CVolProcessedBSSP volATMBS(asset->getProcessedVol(atm.get()));
                    CVolProcessedBSSP volStrikeBS(asset->getProcessedVol(k.get()));

                    double atmVol    = volATMBS->CalcVol(valueDate, effectiveDate);
                    double strikeVol = volStrikeBS->CalcVol(valueDate, effectiveDate);

                    volSpread = strikeVol - atmVol;
                }

                interp = new LinearStrikeSpreadVolRequest(interpLevel/spot,
                                                          valueDate,
                                                          maturity,
                                                          volSpread);
            }
        }
        return interp;
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void AverageHybrid::indVol(CVolProcessedBS* vol, 
                           OutputRequest*   request,
                           CResults*        results) const {
    static const string method = "AverageHybrid::indVol";
    // do 'indicative vol' - don't fail the whole pricing if this fails 
    try {
        double ivol;

        if (avgOut->countFutureSamples(valueDate) > 0) {
            DateTime endDate = avgOut->expectedEndDate(valueDate);
            DateTime inDate  = avgIn->expectedEndDate(valueDate);

            DateTime start = inDate.isGreater(valueDate) ? inDate : valueDate;
            ivol = vol->CalcVol(start, endDate);
        }
        else {
            // all samples set
            ivol = 0.0;
        }
        results->storeRequestResult(request, ivol);
    }
    catch (exception& ) {
    }
}

void AverageHybrid::priceLN(Control* control, CResults* results) const{
    static const string method = "AverageHybrid::priceLN";
    try {
        // if we're past settlement, it's just worth nothing
        if (valueDate.isGreaterOrEqual(instSettle->settles(maturity, asset.get()))) {
            results->storePrice(0.0, discount->getCcy());
            return;
        }

        // estimate value of average forward
        double avgOutFutureSum = avgOut->futureSampleSum(asset.get(),
                                                         valueDate);

        // effective strike 
        double avgOutSumSoFar = avgOut->sumToDate(valueDate);

        double avgInFutureSum = avgIn->futureSampleSum(asset.get(), valueDate);
        double avgInSumSoFar  = avgIn->sumToDate(valueDate);

        // how to interpolate the vol 
        CVolRequestLNSP volRequest = CVolRequestLNSP(volInterp(strike));

        // interpolate the vol
        CVolProcessedBSSP  volBS(asset->getProcessedVol(volRequest.get()));

        // variance
        double avgOutVariance = avgOut->averageVariance(volBS.get(),
                                                        valueDate,
                                                        false);

        double avgInVariance  = avgIn->averageVariance(volBS.get(),
                                                       valueDate,
                                                       true);

        // covariance
        double covariance = avgOut->averageCovariance(avgIn.get(),
                                                      volBS.get(),
                                                      valueDate);
        // discounting 
        double pv = instSettle->pv(maturity, discount.get(), asset.get());

        int avgOutRemaining = avgOut->countFutureSamples(valueDate);

        double avgOutLevel = avgOutSumSoFar + avgOutFutureSum;
        double avgInLevel  = avgInSumSoFar  + avgInFutureSum;

        double value;

        double avgFwd;
        double avgVariance;
        double k;
        
        if (avgOutRemaining > 0) {
            k = strike * avgInLevel - avgOutSumSoFar;
            avgFwd  = avgOutFutureSum;
        }
        else {
            k = strike * avgInLevel;
            avgFwd = avgOutLevel;
        }

        avgVariance = avgOutVariance + avgInVariance - 2 * covariance;

        value = avgBlack(isCall,
                         avgFwd,
                         k,
                         pv,
                         avgVariance);

        results->storePrice(value, discount->getCcy());
 
        if (control && control->isPricing()) {
            OutputRequest* request;
            if ((request = control->requestsOutput(OutputRequest::IND_VOL))) {
                indVol(volBS.get(), request, results);
            }

            recordExtras(control,
                         results,
                         asset.get(),
                         discount.get(),
                         valueDate,
                         maturity,
                         k,
                         avgFwd,
                         avgVariance,
                         value,
                         premiumSettle.get(),
                         instSettle.get(),
                         avgOut.get());
        }
   }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

bool AverageHybrid::avoidVegaMatrix(const IModel* model)
{
    /* this should possibly be false for local vol models etc.
       We may require the model object here */
    return false;
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP AverageHybrid::getSensitiveStrikes(OutputNameConstSP outputName,
                                                 const IModel*      model)
{
    static const string method = "AverageHybrid::getSensitiveStrikes";
    DoubleArraySP sensStrikes;
    try {
        sensStrikes = DoubleArraySP(new DoubleArray(0));
 
        if (avoidVegaMatrix(model)) {
            throw ModelException(method, 
                                 "VEGA_MATRIX is not valid for this instrument");
        }

        if (avgOut->countFutureSamples(valueDate) > 0 ||
            avgIn->countFutureSamples(valueDate)  > 0) {
            CVolRequestSP volRequest(volInterp(strike));
       
            SensitiveStrikeDescriptor sensStrikeDesc;
            sensStrikeDesc.forwardOnly = false;

            asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                       sensStrikeDesc, sensStrikes);
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
    return sensStrikes;
}

// resize delta shift for fwd starting
CSensControl* AverageHybrid::AlterControl(
    const IModel*       modelParams,
    const CSensControl* sensControl) const
{
    static const string method = "AverageHybrid::AlterControl";
    try {
        SensControlPerName* alteredControl = 0;
        // is it effectively forward starting ?
        if ((avgIn->countFutureSamples(valueDate) > 0) && 
            Delta::TYPE->isInstance(sensControl)) {
            const Delta* delta = 
                dynamic_cast<const Delta*>((IObject*)sensControl);
            ShiftSizeCollectorSP shiftSizeVisitor(new ShiftSizeCollector(
                delta,
                strike,
                ShiftSizeCollector::FWD_START_ADJUSTMENT));
            asset->accept(shiftSizeVisitor.get());

            if (Maths::isPositive(shiftSizeVisitor->getShiftSize())) {
                alteredControl = new Delta(shiftSizeVisitor->getShiftSize());
                alteredControl->
                    setMarketDataName(sensControl->getMarketDataName());
            }
        }

        return alteredControl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

bool AverageLoad() {
    return (Average::TYPE && 
            AverageSpot::TYPE && 
            AverageRatio::TYPE && 
            AverageHybrid::TYPE);
}

/** expose below method for 3 hidden child classes */
DoubleArraySP Average::getSensitiveStrikes(Average*          inst,
                                           OutputNameConstSP outputName,
                                           const IModel*      model){
    ISensitiveStrikes* sensStrikes = inst;
    return sensStrikes->getSensitiveStrikes(outputName, model);
}

/** Returns the name of the instrument's discount currency. */
string Average::discountYieldCurveName() const {
    return discount.getName();
}


// for reflection
Average::Average(CClassConstSP clazz): CInstrument(clazz), isAmerican (false){}

/** Invoked when Class is 'loaded' */
void Average::load(CClassSP& clazz){
    REGISTER(Average, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(ISensitiveStrikes);
    FIELD(isCall, "is it a call option");
    FIELD(maturity, "option maturity");
    FIELD_MAKE_OPTIONAL(maturity);
    FIELD(strike, "strike");
    FIELD_MAKE_OPTIONAL(strike);
    FIELD(avgOut, "average-out samples");
    FIELD(exerciseSchedule, "exercise schedule");
    FIELD_MAKE_OPTIONAL(exerciseSchedule);
    FIELD(isAmerican, "isAmerican flag");
    FIELD_MAKE_OPTIONAL(isAmerican);
    FIELD(instSettle, "how instrument settles");
    FIELD(premiumSettle, "how the premium settles");
    FIELD_MAKE_OPTIONAL(premiumSettle);
    FIELD(asset, "identifies underlying");
    FIELD(ccyTreatment, "Currency Treatment");
    FIELD(discount, "identifies discount curve");
    FIELD(valueDate, "valuation Date");
    FIELD_MAKE_OPTIONAL(valueDate);
}

CClassConstSP const Average::TYPE = CClass::registerClassLoadMethod(
    "Average", typeid(Average), load);

class AverageSpotHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AverageSpot, clazz);
        SUPERCLASS(Average);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultAvgSpot);
        FIELD(fwdStarting, "is it fwd starting ?");
        FIELD(startDate, "start Date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(oneContract, "is it contract booking ?");
        FIELD(notional, "notional");
        FIELD_MAKE_OPTIONAL(notional);
        FIELD(initialSpot, "spot @ start date");
        FIELD_MAKE_OPTIONAL(initialSpot);
        FIELD(dimAvgOutFD, "dimension for averaging out period in FD");
        FIELD_MAKE_OPTIONAL(dimAvgOutFD);
        FIELD(stateStdev, "number of standard deviations for average state variables, default=3.5");
        FIELD_MAKE_OPTIONAL(stateStdev);
        FIELD(interpMethod, "interpolation method for state variable: QUADRATIC(default), LINEAR");
        FIELD_MAKE_OPTIONAL(interpMethod);
        FIELD(UseCtrlVar, "flag for using control variates");
        FIELD_MAKE_OPTIONAL(UseCtrlVar);
    }

    static IObject* defaultAvgSpot(){
        return new AverageSpot();
    }
};

CClassConstSP const AverageSpot::TYPE = CClass::registerClassLoadMethod(
    "AverageSpot", typeid(AverageSpot), AverageSpotHelper::load);
   

class AverageRatioHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AverageRatio, clazz);
        SUPERCLASS(Average);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultAvgRatio);
        FIELD(avgIn, "average-in samples");
        FIELD(notional, "notional");
        FIELD(dimAvgOutFD, "dimension for averaging out period in FD");
        FIELD_MAKE_OPTIONAL(dimAvgOutFD);
        FIELD(dimAvgInFD, "dimension for averaging in period in FD");
        FIELD_MAKE_OPTIONAL(dimAvgInFD);
        FIELD(stateStdev, "number of standard deviations for average state variables, default=3.5");
        FIELD_MAKE_OPTIONAL(stateStdev);
        FIELD(interpMethod, "interpolation method for state variable: QUADRATIC(default), LINEAR");
        FIELD_MAKE_OPTIONAL(interpMethod);
        FIELD(UseCtrlVar, "flag for using control variates");
        FIELD_MAKE_OPTIONAL(UseCtrlVar);
    }

    static IObject* defaultAvgRatio(){
        return new AverageRatio();
    }
};

CClassConstSP const AverageRatio::TYPE = CClass::registerClassLoadMethod(
    "AverageRatio", typeid(AverageRatio), AverageRatioHelper::load);

class AverageHybridHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AverageHybrid, clazz);
        SUPERCLASS(Average);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultAvgHybrid);
        FIELD(avgIn, "average-in samples");
    }

    static IObject* defaultAvgHybrid(){
        return new AverageHybrid();
    }
};

CClassConstSP const AverageHybrid::TYPE = CClass::registerClassLoadMethod(
    "AverageHybrid", typeid(AverageHybrid), AverageHybridHelper::load);

// --------------------------------------------------- //
// AvgOutPerf Class, extenal data Class implementation // 
// --------------------------------------------------- //
void AvgOutPerfMaker::validatePop2Object() {
    static const string method("CEqGainKONote::AvgOutPerfMaker::validatePop2Object");
    try {
        if (strike.empty()) {
            throw ModelException(method, "strike is empty");
        }     
            
        if (strike.size() != callPut.size()) {
            throw ModelException(method, "strike and call/put "
                                 "flags are different lengths");
        }     
        if (strike.size() != participation.size()) {
            throw ModelException(method, "strike and participation "
                                 "are different lengths");
        }     
        int i;
        for (i = 0; i < strike.size(); i++) {
            if (callPut[i] != 1 && callPut[i] != -1) {
                throw ModelException(method, "call/put " +
                                     Format::toString(i+1) + " (" + 
                                     Format::toString(callPut[i]) + 
                                     ") must be +1 (call) or -1 (put)");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

class AvgOutPerfMakerHelper {
public:
    // Invoked when Class is 'loaded' 
    // make all params optiona, so as to allow the empty in IMS.
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Average Out Peformance");
        REGISTER(AvgOutPerfMaker, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAvgOutPerfMaker);
        FIELD(callPut,       "1 for Call, -1 for Put.");
        FIELD_MAKE_OPTIONAL(callPut);
        FIELD(participation, "option participation array");
        FIELD_MAKE_OPTIONAL(participation);
        FIELD(strike,        "strike level array");
        FIELD_MAKE_OPTIONAL(strike);
        FIELD(avgSamples,    "Avarage Sampling Schedule and observed level (if in past)");
        FIELD_MAKE_OPTIONAL(avgSamples);
    }
    static IObject* defaultAvgOutPerfMaker(){
        return new AvgOutPerfMaker();
    }
};

CClassConstSP const AvgOutPerfMaker::TYPE = CClass::registerClassLoadMethod(
    "AvgOutPerfMaker", typeid(AvgOutPerfMaker), AvgOutPerfMakerHelper::load);

DateTime AvgOutPerfMaker::getFirstDate(){
    DateTimeArray dates = CashFlow::dates(avgSamples);
    return dates[0];
}

DateTime AvgOutPerfMaker::getLastDate(){
    DateTimeArray dates = CashFlow::dates(avgSamples);
    return dates[dates.size()-1];
}

DateTimeArray AvgOutPerfMaker::getAvgDates(){
    DateTimeArray dates = CashFlow::dates(avgSamples);
    return dates;
}

//class AvgOutPerfHelper {
//public:
//    static void load(CClassSP& clazz){
//        REGISTER(AvgOutPerf, clazz);
//        SUPERCLASS(CObject);
//        EMPTY_SHELL_METHOD(defaultAvgOutPerf);
//        FIELD(callPut,       "1 for Call, -1 for Put.");
//        FIELD(participation, "option participation array");
//        FIELD(strike,        "strike level array");
//        FIELD(avgOut,    "Avarage Sampling Schedule and observed level (if in past)");
//    }
//    static IObject* defaultAvgOutPerf(){
//        return new AvgOutPerf();
//    }
//};
//
//// AvgOutPerf Class, internal class implementation // 
//CClassConstSP const AvgOutPerf::TYPE = CClass::registerClassLoadMethod(
//    "AvgOutPerf", typeid(AvgOutPerf), AvgOutPerfHelper::load);

AvgOutPerf::AvgOutPerf(AvgOutPerfMakerSP avgOutPerf) {
    static const string method("AvgOutPerf::AvgOutPerf");
    try {

        // copy data from external data class
        callPut = avgOutPerf->callPut;
        participation = avgOutPerf->participation;
        strike = avgOutPerf->strike;
        
        int nSize = avgOutPerf->avgSamples.size();
        DateTimeArray avgDates(nSize);
        DoubleArray past(nSize);
        DoubleArray weights(nSize);
        for (int i = 0; i < nSize; i++){
            avgDates[i] = avgOutPerf->avgSamples[i].date;
            past[i] = avgOutPerf->avgSamples[i].amount;
            weights[i] = 1.0 / nSize;
        }
        avgOut = SampleListSP(new SampleList(avgDates, 
                                             past,
                                             weights));   
        
        isAlreadySetInst = false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// embedded the instrument information.
void AvgOutPerf::setInstInfo(const InstrumentSettlementSP instSettle,
                 const InstrumentSettlementSP premiumSettle,
                 const CAsset*               asset,
                 string                       ccyTreatment,
                 const YieldCurveConstSP      discount,
                 const DateTime               valueDate)
{
    this->instSettle = instSettle;
    this->premiumSettle = premiumSettle;
    this->asset = asset;
    this->ccyTreatment = ccyTreatment;
    this->discount = discount;
    this->valueDate = valueDate;
    isAlreadySetInst = true;
}

//AvgOutPerf::roll(){
//    static const string method = "AvgOutPerf::roll";
//    try{
//        if (isAlreadySetInst)
//            avgOut->roll(asset.get(), valueDate, newDate, useSpot);
//        else
//            throw ModelException(method, "setInstInfo is not called before this function." );
//    }
//    catch (exception& e) {
//        throw ModelException(e, method);
//    }
//}
//

double AvgOutPerf::AvgValue(double const strikeScale, 
                double const refLevel, 
                double const notl,
                const DateTime&             startDate) const {
    static const string method = "AvgOutPerf::AvgValue";
    try {     

        if(!isAlreadySetInst)
            throw ModelException(method, "setInstInfo is not called before this function." );

        double value = 0.0;

        CClosedFormLN cfln("VolPreferred");

        // use the pastAvgOut if there is.  To DO...
        // SampleListSP sample = !aop->pastAvgOut.get() ? aop->avgOut: aop->pastAvgOut;

        for (int i = 0; i < strike.size(); i++) {
            InstrumentSP avg(Average::makeAvgSpot(callPut[i]==1,
                                                  avgOut->getLastDate(),
                                                  strike[i]*strikeScale,
                                                  avgOut.get(),
                                                  instSettle.get(),
                                                  premiumSettle.get(),
                                                  asset,
                                                  ccyTreatment,
                                                  discount.get(),
                                                  valueDate,
                                                  startDate > valueDate,
                                                  startDate,
                                                  false, // always notional
                                                  notl,
                                                  refLevel));
                                   
            value += cfln.calcPrice(avg.get(), 0) * participation[i];
        }
        value /= discount->pv(startDate);        // return value at startDate, not valueDate.

        return value;                  
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// function to return private member.
IntArray AvgOutPerf::getCallPutArray()
{
    return callPut;
}

DoubleArray AvgOutPerf::getParticipation()
{
    return participation;
}

DoubleArray AvgOutPerf::getStrikeArray()
{
    return strike;
}

DateTimeArray AvgOutPerf::getAvgDates()
{
    return avgOut->getDates(); 
}

/*
// not fully tested //
double AvgOutPerf::getAdjFwd(double const strike, 
                            double const refLevel, 
                            double const notl,
                            const InstrumentSettlement* instSettle,
                            const InstrumentSettlement* premiumSettle,
                            const CAsset*               asset,
                            string                      ccyTreatment,
                            const YieldCurve*           discount,
                            const DateTime&             valueDate,
                            const SampleListSP          avgIn) const {
    static const string method = "AvgOutPerf::AvgValue";
    try {     
        DateTime maturity = avgOut->getLastDate();
        AverageRatio* avg = new AverageRatio(true,                   //isCall
                                           maturity,
                                           1.0,                     //strike,
                                           avgOut.get(),
                                           instSettle,
                                           premiumSettle,
                                           asset,
                                           ccyTreatment,
                                           discount,
                                           valueDate,
                                           avgIn.get(),
                                           notl);
                                   
        // estimate value of average forward
        double avgOutFutureSum = avgOut->futureSampleSum(asset,
                                                         valueDate );

        double avgOutSumSoFar = avgOut->sumToDate(valueDate);

        double avgInFutureSum = avgIn->futureSampleSum(asset, valueDate);
        double avgInSumSoFar  = avgIn->sumToDate(valueDate);

        // how to interpolate the vol 
        CVolRequestLNSP volRequest = CVolRequestLNSP(avg->volInterp(avg->strike));

        // interpolate the vol
        CVolProcessedBSSP  volBS(asset->getProcessedVol(volRequest.get()));

        // variance
        double avgOutVariance = avgOut->averageVariance(volBS.get(),
                                                        valueDate,
                                                        false);

        double avgInVariance  = avgIn->averageVariance(volBS.get(),
                                                       valueDate,
                                                       true);

        // covariance
        double covariance = avgOut->averageCovariance(avgIn.get(),
                                                      volBS.get(),
                                                      valueDate);
        // discounting 
        double pv = instSettle->pv(maturity, discount, asset);

        int avgInRemaining  = avgIn->countFutureSamples(valueDate);
        int avgOutRemaining = avgOut->countFutureSamples(valueDate);

        double avgOutLevel = avgOutSumSoFar + avgOutFutureSum;
        double avgInLevel  = avgInSumSoFar  + avgInFutureSum;

        double effectiveStrike;
        double avgFwd;
        double avgVariance;

        double fwd = 0.0;
        if (avgInRemaining > 0) {       
            effectiveStrike = strike * (avgInSumSoFar + avgInFutureSum);
            avgFwd = avgOutLevel;

            avgVariance = avgOutVariance + avgInVariance - 2 * covariance;

            double sqrtVar = sqrt(avgVariance);

            fwd = avgOutLevel/avgInLevel * 
                  exp(0.5*avgVariance - 0.5*avgOutVariance + 
                      0.5*avgInVariance);
        }
        else {
            double k;

            if (avgOutRemaining > 0) {
                k = strike * avgInLevel - avgOutSumSoFar;
                fwd  = avgOutFutureSum;
            }
            else {
                k = strike * avgInLevel;
                fwd = avgOutLevel;
            }
        }

        return fwd;                  
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
*/

/*
// AvgOutPerf Class, just having only one ATM avg out call // 
AvgOutPerf::AvgOutPerf(const DateTimeArray avgDates,
                       const DoubleArray past,
                       const DoubleArray weights){
    static const string method("AvgOutPerfMaker::setup");
    try {

        // copy data from external data class
        callPut = new IntArray(1,1);    // just one call
        participation = new DoubleArray(1, 1.0);
        strike = new DoubleArray(1, 1.0);
        
        avgOut = SampleListSP(new SampleList(avgDates, 
                                             past,
                                             weights));                                                    
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double AvgOutPerf::ExpectedAvgInTime(){
    static const string method("AvgOutPerf::ExpectedAvgInTime");
    try {


        // copy data from external data class
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}*/


// helper clas for AvgIn.  I think it's necessary to keep AvgIn object valid at tweaking.  No?
class AvgInHelper {
public:
    static void load(CClassSP& clazz){
        REGISTER(AvgIn, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAvgIn);
        FIELD(avgIn,       "average in sample list");
    }
    static IObject* defaultAvgIn(){
        return new AvgIn();
    }
};

// AvgIn Class, internal class implementation // 
CClassConstSP const AvgIn::TYPE = CClass::registerClassLoadMethod(
    "AvgIn", typeid(AvgIn), AvgInHelper::load);

// constructor
AvgIn::AvgIn(SampleListSP avgIn, const CAsset* asset):CObject(TYPE), avgIn(avgIn){
    this->asset = CAssetWrapper(copy(asset));
}

// return Fwd and ATM BS variance.
void AvgIn::getFwdAndBSVar(const DateTime valDate,     //(I) today
                           int iSample,                //(I) which sample date you want?
                           double &fwd,               //(O) fwd level
                           double &bsVar)              //(O) BS variance.   
{
    static const string method("AvgIn::getPeakSpot");
    try {
        DateTimeArray avgDates = avgIn->getDates();
        int nSize = avgDates.size();
        if (iSample >= nSize){
            throw ModelException(method, "Requested sample point, iSample[" + Format::toString(iSample+1) 
                                        + "] is larger than size of AvgIn["
                                        + Format::toString(nSize) +"]" );
        }
        CVolRequestLNSP     volRequest = CVolRequestLNSP(new ATMVolRequest());
        CVolProcessedBSSP   volBS(asset->getProcessedVol(volRequest.get()));                

        DateTimeArray tmpDt(2);
        tmpDt[0] = valDate;
        tmpDt[1] = avgDates[iSample]; //last avg in date
        DoubleArray vars(2);
        volBS->CalcVar(tmpDt, volBS->fromFirst, vars);

        bsVar = vars[0];       
        fwd = asset->fwdValue(tmpDt[1]);            
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// return expectation average in with known sStart and sEnd, using Brownian Bridge.
double AvgIn::getMeanCondAvgIn(const DateTime valDate,
                               const double   sStart, 
                               const double   sEnd, 
                               const double   varEnd){
    static const string method("AvgIn::getAvgIntMedDist");
    try {
        CVolRequestLNSP volRequest = CVolRequestLNSP(new ATMVolRequest());
        CVolProcessedBSSP  volBS(asset->getProcessedVol(volRequest.get()));        

        // --- calculate variance of AvgIn with Known sEnd ---
        DateTimeArray orgAvgDates   = avgIn->getDates();
        DoubleArray   orgPastLevels = avgIn->getValues();
        DoubleArray   orgWeights    = avgIn->getWeights();
        
        int nSize = orgAvgDates.size()-1;
        double dropedWeight = 1.0/(double)(nSize+1);

        DateTimeArray   avgIntMedDt(nSize);
        DoubleArray     pastLevels(nSize,0.0);
        DoubleArray     weights(nSize, 1.0/nSize);
        
        for (int iAvg=0;iAvg<nSize;iAvg++){
            avgIntMedDt[iAvg] = orgAvgDates[iAvg];
            pastLevels[iAvg]  = orgPastLevels[iAvg];
            orgWeights[iAvg]  = orgWeights[iAvg];            
        }
        SampleListSP avgIntMed = SampleListSP(new SampleList(avgIntMedDt, 
                                             pastLevels,
                                             weights));
        double varAvgIntMed = avgIntMed->averageVariance(volBS.get(), valDate, true);

        // --- set up Avg-In level ---
        // Brownian Bridge
        double kappa = varAvgIntMed/varEnd;
        double xEnd  = log(sEnd/sStart);
        double mean  = xEnd * kappa;
        
        //E[AvgIntMed | S_end]        
        double expectedAvgIntMed = sStart * exp(mean);
        
        // sum up AvgIntMed + sEnd
        double expectOfCondAvg = (1.0-dropedWeight)*expectedAvgIntMed + dropedWeight*sEnd;

        return expectOfCondAvg;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE
