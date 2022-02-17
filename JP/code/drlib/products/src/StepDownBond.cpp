//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StepDownBond.hpp
//
//   Description : Step down bond instrument. 
//                    if (LinearCoupon)
//                        coupon is scaled linearly between top and bottom rungs (must have 2 input rungs)
//                    else
//                        coupon is stepped according to Ladder rungs
//
//                    step coupon uses multiple payoffs of KI digitals on a single tree
//                    linear coupon uses 61 step coupons plus 2nd order Richardson extrapolation
//                    ladder must be in ascending order !!!
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/StepDownBond.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/VegaParallel.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

// helpers
class CStepDownBondHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(CStepDownBond, clazz);
        SUPERCLASS(Generic1Factor);
        EMPTY_SHELL_METHOD(defaultStepDownBond);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate);

        // Ladder parts
        FIELD(LadderLevelArray, "ladder levels, 1D array representing ladder level matrix");
        FIELD(RungNumber, "number of ladder rungs");
        FIELD(IntraDayMonitor, "true for hi/low continuous monitoring, false if once a day");
        FIELD(CurrentMinOrMax, "minimum for step down / max for step up, reached in current monitoring period - must be updated with every breach !!");
        FIELD_MAKE_OPTIONAL(CurrentMinOrMax);   //need to be remove
        FIELD(LadderCouponArray, "ladder coupons, 1D array representing ladder coupon matrix");

        // Fixed Leg
        FIELD(fixedLeg,     "Fixed Leg");
        
        // Libor parts
        FIELD(liborLeg,     "Libor Leg");

        // Optional parts.
        FIELD(isMonitorOnMax, "false(default)=step down type, true=step up");
        FIELD_MAKE_OPTIONAL(isMonitorOnMax);
        FIELD(LinearCoupon, "false(default)=step coupons, true=linearly interpolated coupons between ladder levels");
        FIELD_MAKE_OPTIONAL(LinearCoupon);
        FIELD(LinearGridSize, "grid size for linear coupon interp (default=61)");
        FIELD_MAKE_OPTIONAL(LinearGridSize);
        FIELD(VolStrike, "vol strike for log-normal only (default=spot). Volatility interpolated at this strike in log-normal model.");
        FIELD_MAKE_OPTIONAL(VolStrike);
        FIELD(SwitchPlainSwap, "Price with PlainSwap or Not.  Default = True (with PlainSwap).");
        FIELD_MAKE_OPTIONAL(SwitchPlainSwap);
        //  When SwitchPlainSwap = false, There is no KNOWN_CASH_FLOW info in this product.
        
        FIELD(pastValues,    "historical samples to know max/min values so as to determine Ladder Coupon.");
        FIELD_MAKE_OPTIONAL(pastValues);
        // 
        FIELD(DEBUG_num_segments, "For DR use only. Number of time line segments, default=1.");
        FIELD_MAKE_OPTIONAL(DEBUG_num_segments);
    }

    static IObject* defaultStepDownBond(){
        return new CStepDownBond();
    }
};


CClassConstSP const CStepDownBond::TYPE = CClass::registerClassLoadMethod(
    "StepDownBond", typeid(CStepDownBond), CStepDownBondHelper::load);
bool  CStepDownBondLoad() {
    return (CStepDownBond::TYPE != 0);
   }


// constructor
CStepDownBond::CStepDownBond(): Generic1Factor(TYPE)
{
    isMonitorOnMax = false;
    LinearCoupon = false;
    LinearGridSize = 61;
    VolStrike = 0.0;
    DEBUG_num_segments = 1;
    SwitchPlainSwap = true;
};

/** get historical max or min value arrays. If it's not sampled yet, return 0. */
CashFlowArray CStepDownBond::getHistMinOrMaxValues() const{
    
    int numPeriod = fixedLeg->getSize();
    DateTimeArray pastDates = CashFlow::dates(pastValues->getSamples());
    DoubleArraySP pastLevels = CashFlow::amounts(pastValues->getSamples());
    CashFlowArray histMinOrMax(numPeriod);
    
    DateTime theDate;
    double   minORmax;
    for (int j=0; j<numPeriod; j++){
        theDate = fixedLeg->AccrueStartDates[j];
        // looking for the max/min in the period only for the sample is started period.
        if (theDate <= valueDate){
            // Starting from the first date in the sampling period.
            // When the sampling doesn't have the sample of the accrue start date,
            // try to find out the first sample date before the value date.
            int i = Neighbour(theDate,pastDates,0,pastDates.size()-1,1);
            if (i<0 || pastDates[i] > valueDate){
                throw ModelException("CStepDownBondProd::getHistMinOrMaxValues", 
                                     "pastValues doesn't include any sample date for the period("
                                     + Format::toString(j) + ").  pastValues should have at least one sample between " 
                                     + theDate.toString() + " and " + valueDate.toString());
            }
            minORmax=(*pastLevels)[i];        
            theDate = pastDates[i];
            for (;i<pastDates.size();i++) {
                if (pastDates[i] <= fixedLeg->AccrueEndDates[j]
                    && pastDates[i] <  valueDate){
                    if (isMonitorOnMax && minORmax < (*pastLevels)[i]){
                        minORmax = (*pastLevels)[i];
                        theDate  = pastDates[i];
                    }
                    else if (!isMonitorOnMax && minORmax > (*pastLevels)[i]){
                        minORmax = (*pastLevels)[i];
                        theDate  = pastDates[i];
                    }                
                }
            }
            histMinOrMax[j] = CashFlow(theDate,minORmax);
        }
        else 
            histMinOrMax[j] = CashFlow(theDate,0.0);        
    }    
    return histMinOrMax;
}

// **** now a copy of Vanilla, to do : add more validations
void CStepDownBond::Validate()
{
    static const string method = "CStepDownBond::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object
    if (!asset){
        throw ModelException(method, "Asset is null");
    }
    if (!discount){
        throw ModelException(method, "Discount YC is null");
    }

    if (fixedLeg->AccrueEndDates.size() <1) {  // Need to check the size before get maturity.
        throw ModelException(method, "no ladder end schedule supplied.");
    }

    AssetUtil::assetCrossValidate(asset.get(),
                                  fwdStarting,
                                  startDate,
                                  valueDate,
                                  discount,
                                  this);
        
    if ( FXAsset::TYPE->isInstance(asset.get()) )
    {
        throw ModelException("CStepDownBond::Validate",
                             "Options on FX assets are not allowed yet");
    }

    // instrument specific validations
//    int firstLadderIndex = -1;
//    int reachedIndex = -1; // no ladder reached yet
//    double coup_gap = 0.0;

    if (fixedLeg->AccrueStartDates.size() <1)
        throw ModelException(method, "no ladder start schedule supplied.");

    if (RungNumber < 2)
        throw ModelException(method, "there at least two rungs.");

    if (LinearCoupon && RungNumber != 2)
        throw ModelException(method, "linear coupon payoff requires exactly 2 ladder rungs.");

    if (fwdStarting || startDate > valueDate)
        throw ModelException(method, "forward start not yet supported.");

    if (LinearCoupon && LinearGridSize <10 )
        throw ModelException(method, "Linear grid size too small for linear coupon calculation.");

    if (fixedLeg->AccrueEndDates.size() != fixedLeg->AccrueStartDates.size())
        throw ModelException(method, "fixedLeg->AccrueStartDates and LadderEndDate have different size array.");

//  Need to check this after remove all dummy data from input.
/*    DateTimeArray critDates; 
      for (i=0; i<fixedLeg->AccrueEndDates.size(); i++)
      {
      if (fixedLeg->AccrueStartDates[i] > fixedLeg->AccrueEndDates[i])
      throw ModelException(method, "fixedLeg->AccrueStartDates["+Format::toString(i)+"] is later than LadderEndDate["+Format::toString(i)+"] ");
      if (i>0 && fixedLeg->AccrueStartDates[i] < fixedLeg->AccrueEndDates[i-1])
      throw ModelException(method, "fixedLeg->AccrueStartDates["+Format::toString(i)+"] is earlier than previous LadderEndDate["+Format::toString(i-1)+"] ");
      if (fixedLeg->AccrueEndDates[i] > fixedLeg->PaymentDatesArray[i])
      throw ModelException(method, "fixedLeg->AccrueEndDates["+Format::toString(i)+"] is later than fixedLeg->PaymentDatesArray["+Format::toString(i)+"] ");
      }*/
    
    int coupNum = fixedLeg->PaymentDatesArray.size();
    int numsize = LadderLevelArray.size();
    if (numsize != coupNum*RungNumber)
        throw ModelException(method, "Ladder Level Should contains CoupNum * RungNumber size.    Your input are (size, CoupNum, RungNum) ="+Format::toString(numsize)+","+Format::toString(coupNum)+","+Format::toString(RungNumber));

    numsize = LadderCouponArray.size();
    if (numsize!= coupNum*(RungNumber+1))
        throw ModelException(method, "Ladder Level Should contains CoupNum * (RungNumber+1) size. Your input are (size, CoupNum, RungNum) ="+Format::toString(numsize)+","+Format::toString(coupNum)+","+Format::toString(RungNumber));

    if (DEBUG_num_segments <=0 || DEBUG_num_segments>10)
        throw ModelException(method, "DEBUG_num_segments must be from 1 to 10");

    if (valueDate > fixedLeg->AccrueStartDates[0]){
        if (!pastValues)
            throw ModelException(method, "Already Ladder sampling is started.  You need give pastValues");
    }

    DateTime matDate = fixedLeg->AccrueEndDates[fixedLeg->AccrueEndDates.size()-1];
    DateTime paymentDate = instSettle->settles(matDate, asset.get());    
    if (paymentDate < fixedLeg->PaymentDatesArray[fixedLeg->PaymentDatesArray.size()-1])
        throw ModelException(method, "Last PaymentDatesArray is later than paymentDate in settlement information");


}

// initiate GetMarket 
/** returns a vol request for log-normal vol */
CVolRequestConstSP CStepDownBond::makeVolRequest() const
{
    // get strike and maturity date from instrument
    DateTime        matDate = fixedLeg->AccrueEndDates[fixedLeg->AccrueEndDates.size()-1];

    double volStrike;
    if (VolStrike == 0)
        volStrike = asset->getSpot();
    else
        volStrike = VolStrike;


    DateTime imntStartDate = fwdStarting? 
                         startDate: valueDate;

    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(volStrike, imntStartDate, 
                                   matDate, fwdStarting));
    return volRequest;
}


/** returns the current value date */
DateTime CStepDownBond::getValueDate() const
{
   return valueDate;
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool CStepDownBond::priceDeadInstrument(CControl* control, CResults* results) const
{
    static string method = "CStepDownBond::priceDeadInstrument";
    DateTime matDate = fixedLeg->AccrueEndDates[fixedLeg->AccrueEndDates.size()-1];
    DateTime settledate = fixedLeg->PaymentDatesArray[fixedLeg->PaymentDatesArray.size()-1];
    double value =0.0;

    if (valueDate < matDate) 
        return false;
    else{
        if (valueDate < settledate) {
        
            value += liborLeg->getPV(valueDate, discount.get());
            value += fixedLeg->getPV(valueDate, discount.get());
            // for fixed leg, it's better to be changed to search Ladder Coupon from pastValue and LadderLevels, in future.
        }
        results->storePrice(value, discount->getCcy());        
        recordOutputRequests(control, results, value);
        return true;
    }
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool CStepDownBond::sensShift(Theta* shift)
{   
    // roll today 
    DateTime rollDate = shift->rollDate(valueDate);
    // Need to give fixing rate if valuation date is fixing dates.
    liborLeg->setFixingforThetaShift(valueDate,discount.get(),rollDate);
    
    valueDate = rollDate;

    return true;
};

/** when to stop tweaking */
CDateTime CStepDownBond::endDate(const Sensitivity* sensControl) const {
    DateTime lastCoupDate = fixedLeg->PaymentDatesArray[fixedLeg->PaymentDatesArray.size()-1];
    DateTime matDate = fixedLeg->AccrueEndDates[fixedLeg->AccrueEndDates.size()-1];
    DateTime instEnd  = instSettle->settles(lastCoupDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = (instEnd > assetEnd ? instEnd : assetEnd);
    return end;
}


void CStepDownBond::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime matDate = fixedLeg->AccrueEndDates[fixedLeg->AccrueEndDates.size()-1];
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

void CStepDownBond::recordOutputRequests(Control* control, Results* results, double fairValue) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime       matDate = fixedLeg->AccrueEndDates[fixedLeg->AccrueEndDates.size()-1];
        double         indVol;
        // calculate indicative vol
        try{
            if ( matDate.isGreater(valueDate) )
            {

                DateTime imntStartDate = fwdStarting? 
                                 startDate: valueDate;

                // get vol request
                CVolRequestConstSP lnVolRequest(makeVolRequest());
                
                // interpolate the vol
                CVolProcessedSP  vol(asset->getProcessedVol(lnVolRequest.get()));
                // cast to the type of vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                // this should never happen if our get market data has worked properly
                if (!vol){
                    throw ModelException("CStepDownBond::recordOutputRequests", 
                                         "No Black Scholes Vol");
                }

                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
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


        //calculate Swap related Value and Store Values
        //libor
        double legResult;
        int i,j;
        CashFlowArrayConstSP CflL = liborLeg->getCashFlowArray(valueDate, discount.get());
        if (CflL->size()>0)
        {
            CashFlowArraySP cfl (new CashFlowArray);
            for (i=0; i<liborLeg->getSize(); i++)
                    cfl->push_back((*CflL)[i]);
            results->storeGreek(cfl, Results::DEBUG_PACKET, OutputNameSP(new OutputName("liborLeg")));
        }
        legResult = liborLeg->getPV(valueDate, discount.get());
        results->storeScalarGreek(legResult, Results::DEBUG_PACKET, OutputNameSP(new OutputName("liborLegValue")));

        //fixed cash flow
        CashFlowArrayConstSP CflF = fixedLeg->getCashFlowArray();
        if (CflF->size()>0)
        {
            CashFlowArraySP cff (new CashFlowArray);
            for (j=0; j<fixedLeg->getSize(); j++)
                    cff->push_back((*CflF)[j]);
            results->storeGreek(cff, Results::DEBUG_PACKET, OutputNameSP(new OutputName("fixedLeg")));
        }
        legResult = fixedLeg->getPV(valueDate, discount.get());
        results->storeScalarGreek(legResult, Results::DEBUG_PACKET, OutputNameSP(new OutputName("fixedLegValue")));

        // Add PAYMENT_DATES
        OutputRequest* request = NULL;
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTime paymentDate = instSettle->settles(matDate, asset.get());
            DateTimeArray date(1, paymentDate);

            if (CflL->size()>0)
            {
                for (i=0; i<liborLeg->getSize(); i++) {
                    date.push_back((*CflL)[i].date);
                }
            }
            if (CflF->size()>0)
            {
                for (j=0; j<fixedLeg->getSize(); j++) {
                    date.push_back((*CflF)[j].date);
                }        
            }        

            OutputRequestUtil::recordPaymentDates(control,results,&date); 
        }
        // Add KNOW_CASHFLOWS
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            // these need to be in increasing date order to aggregate
            // hence the merge rather than push_back which handles this

            // Fixed Leg is not requires.  
            // Every known information should be calc from pastValues.

            // Pick up only known Libor Leg.
            CashFlowArraySP knownCflL (new CashFlowArray(0));
            if (CflL->size()>0 && SwitchPlainSwap)
            {
                for (i=0; i<liborLeg->getSize(); i++) {
                    if(valueDate > liborLeg->AccrualDates[i])
                        knownCflL->push_back((*CflL)[i]);
                }
            }

            // calculate known Ladder coupon
            CashFlowArraySP knownLdCpn(new CashFlowArray(0));
            if (valueDate > fixedLeg->AccrueStartDates[0]){
                CashFlowArray histMinOrMax = getHistMinOrMaxValues();
                if (histMinOrMax.getLength()>0){
                    for (i=0; i<histMinOrMax.getLength() ; i++){
                        if(valueDate > fixedLeg->AccrueEndDates[i])
                        {//find Ld Level achieved
                            for (j=0; j<RungNumber; j++){
                                //  LadderLevelArray[i] is always increasing order.  
                                //  LadderCouponArray[i] contains the Coupon if the extreame is between i and i+1.
                                if   (histMinOrMax[i].amount < LadderLevelArray[i*RungNumber+j])
                                    break;
                            }
                            double knownFixCpn = LadderCouponArray[i*(RungNumber+1)+j];
                            if (!SwitchPlainSwap)
                                knownFixCpn -= fixedLeg->CouponAmounts[i];
                            knownLdCpn->push_back(CashFlow(fixedLeg->PaymentDatesArray[i], knownFixCpn));
                        }
                    }
                }
            }

            // now glue them all together
            //CashFlowArraySP merge1(CashFlow::merge(knownCflL, knownCflF));
            //CashFlowArraySP merge2(CashFlow::merge(merge1, knownLdCpn));
            CashFlowArraySP knownCfl(CashFlow::merge(knownCflL, knownLdCpn));

            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    knownCfl.get());   
        }

        addOutputRequests(control,
                                        results,
                                        fairValue,
                                        indVol);
    }
}


/************************************************************************************************************************/
class CStepDownBondFDProd: public LatticeProdEDRIns 
{
public:

    CStepDownBondFDProd(const CStepDownBond* sdb, FDModel * model) :
        LatticeProdEDRIns( model, 1, sdb->LinearCoupon ? sdb->LinearGridSize : sdb->RungNumber + 1 ),
        inst(sdb)
    {
        if( ! tree1f )
        {
            throw ModelException( "CStepDownBondFDProd::CStepDownBondFDProd", "Instrument of type "+
                                 inst->getClass()->getName() +
                                 " can be priced by CTree1f only" );
        }

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );
    }

    virtual ~CStepDownBondFDProd(){}

    virtual CVolRequestConstSP GetLNRequest() const
    {
        CVolRequestConstSP volRequest(inst->makeVolRequest());
        return volRequest;
    }

    /** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

    
    virtual void init(CControl* control) const
    {
        static const string method = "CStepDownBondFDProd::Init";
        int i;
        
        try 
        { 
            if( tree1f )
            {
                // default to NODE_INSERTION smoothing
                if (tree1f->GetSmoothMethod() == CTree1f::DEFAULT) 
                    tree1f->SetSmoothMethod(CTree1f::NODE_INSERTION);
                
                tree1f->NumOfPrice = numPrices;
                tree1f->NumOfInsertNode =
                    ( tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ? numIns : 0 );    
                
                // default to DollarDivTreament = false !!! to be removed after EAS/EIS can handle it correctly !!!
                tree1f->SetDivAmountTreatment(false);    

                // this needs change if fwd start tree has to start today !!!
                if (inst->fwdStarting && inst->startDate>inst->valueDate)
                    tree1f->controlSameGridFwdStart(inst->ccyTreatment);                
            }

            /** customize tree parameters here and set up the tree */
            DateTimeArray segDates;
            segDates.resize(inst->DEBUG_num_segments+1); // to try using more segments

            segDates[0] = getStartDate();            

            DateTime mat =inst->fixedLeg->AccrueEndDates[inst->fixedLeg->AccrueEndDates.size()-1];
            IntArray density (inst->DEBUG_num_segments);
            
            // ********** testing simple density ratios, 3:2:1 ...need to change this after testing !!!
            double t = segDates[0].yearFrac(mat);

            for (i=0; i<inst->DEBUG_num_segments; i++)
            {
                density[i] = inst->DEBUG_num_segments-i;
                segDates[i+1] = segDates[i].rollDate((int)(365*t/inst->DEBUG_num_segments));
            }
            // make sure the last date is spot on
            segDates[inst->DEBUG_num_segments] = mat;                  

            // add critical dates
            model->addCritDates( inst->fixedLeg->AccrueStartDates );
            model->addCritDates( inst->fixedLeg->AccrueEndDates );
            model->addCritDates( inst->fixedLeg->PaymentDatesArray );

            // prepare timeline set up
            model->initSegments( segDates, density );
        }
        catch (exception& e) 
        {
            throw ModelException(e, method);
        }

    }

    /** initialise product specific data */
    virtual void initProd()
    {
        static const string method = "CStepDownBondFDPRod::InitProd";
        try
        {    
             int i, j;    
            int coupNum = inst->fixedLeg->PaymentDatesArray.size();    
            
            initSlices( numPrices );
            initInsertNode();
                                
            vector<vector<double> > tmpLdLevel(coupNum);     // num of ladder can be less, determined below
            vector<vector<double> > tmpLdLevelRaw(coupNum); // num of ladder can be less, determined below
            vector<vector<double> > tmpLdCoupon(coupNum);   // num of ladder can be less, determined below

            int ladder_size;
            int rung_num = inst->RungNumber;
            if (!inst->LinearCoupon)
                ladder_size = rung_num;
            else
            {// default = 61; Other case is multiple of 6(2*3)+1 for 2nd order Richardson extroplation
                if (inst->LinearGridSize != 61){
                    ladder_size = inst->LinearGridSize/6;
                    ladder_size = 6*(ladder_size+1)+1;
                }
                else
                    ladder_size = inst->LinearGridSize;
            }
            int real_ld_num =0;
            for (i=0; i<coupNum; i++)
            {
                tmpLdLevel[i].resize(ladder_size);
                tmpLdLevelRaw[i].resize(ladder_size);
                tmpLdCoupon[i].resize(ladder_size+1);
                for (j=0; j<ladder_size; j++)
                {
                    if ( j==0 || !inst->LinearCoupon)
                    {
                        tmpLdLevel[i][j] = inst->LadderLevelArray[j + i*rung_num];
                        tmpLdLevelRaw[i][j] = tmpLdLevel[i][j];
                        tmpLdCoupon[i][j] = inst->LadderCouponArray[j + i*(rung_num+1)];
                    }
                    else
                    {// Linear Coupon : There are only two ladder levels for linear case.  (it's already validated.)
                        tmpLdLevel[i][j] = tmpLdLevel[i][j-1]+(inst->LadderLevelArray[i*rung_num+1] - inst->LadderLevelArray[i*rung_num])*j/(rung_num-1);
                        tmpLdLevelRaw[i][j] = tmpLdLevel[i][j];
                        tmpLdCoupon[i][j] = tmpLdCoupon[i][j-1] + (inst->LadderCouponArray[i*(rung_num+1)+1]-inst->LadderCouponArray[i*(rung_num+1)])*j/rung_num; // coupon is for gap
                    }
                    if (j>0 && tmpLdLevel[i][j] < tmpLdLevel[i][j-1])
                        throw ModelException(method, "Ladder Level should be increasing order");
                }
                // one more for coupon
                if (inst->LinearCoupon)
                    tmpLdCoupon[i][j] = tmpLdCoupon[i][j] + (inst->LadderCouponArray[i*(rung_num+1)+1]-inst->LadderCouponArray[i*(rung_num+1)])*j/rung_num; // coupon is for gap
                else
                    tmpLdCoupon[i][j] = inst->LadderCouponArray[j + i * (rung_num+1)];
                // count real ladder sampling period.
                if (tmpLdLevel[i][0] != 0.0 || tmpLdLevel[i][ladder_size-1] != 0.0) // all zero ladder level means no coupon payment at this date.  So, it's dummy!!
                    real_ld_num ++;
            }
    
            if (inst->LinearCoupon) // create rungs (60 gaps default) for linear coupon   //Over wright RungNumber from now on.
            {
                LdGridSize = inst->LinearGridSize; // We must need to do this later.  Not Yet.
                rung_num = inst->RungNumber;
            }
            else
                LdGridSize = inst->RungNumber;


                                    // removing dummy data and re-arrange.//
            // Ladder Coupon PayDates and LiborPayDates are extract from fixedLeg->PaymentDatesArray, which contains all payment dates.
            // Ladder Coupon are rearranged to Base Coupon and Increment LdCoupon.
            DateTimeArray critDates; 
            LdLevel.resize(real_ld_num);    
            LdLevelRaw.resize(real_ld_num); 
            LdCoupon.resize(real_ld_num);   
            int n =0;   // Ladder Sample Number
            for (i=0; i<coupNum; i++)
            {   
                if (tmpLdLevel[i][0] != 0.0 || tmpLdLevel[i][rung_num-1] != 0.0) // all zero ladder level means no coupon payment at this date.  So, it's dummy!!
                {
                    // Ladder Start/End Pay dates are stored into Ld variable to remove meaningless dummy data.
                    LdCoupPayDates.push_back(inst->fixedLeg->PaymentDatesArray[i]);
                    LdStartDates.push_back(inst->fixedLeg->AccrueStartDates[i]);
                    LdEndDates.push_back(inst->fixedLeg->AccrueEndDates[i]);

                    //  Need to check this after remove all dummy data from input.
                    if (LdStartDates[n] > LdEndDates[n])
                        throw ModelException(method, "fixedLeg->AccrueStartDates["+Format::toString(i)+"] is later than LadderEndDate["+Format::toString(i)+"] Or set All Ladder Level = 0 to not use this date for ladder sampling");
                    if (n>0 && LdStartDates[n] < LdEndDates[n-1])
                        throw ModelException(method, "fixedLeg->AccrueStartDates["+Format::toString(i)+"] is earlier than previous LadderEndDate["+Format::toString(i-1)+"] Or set All Ladder Level = 0 to not use this date for ladder sampling");
                    if (LdEndDates[n] > inst->fixedLeg->PaymentDatesArray[i])
                        throw ModelException(method, "fixedLeg->AccrueEndDates["+Format::toString(i)+"] is later than fixedLeg->PaymentDatesArray["+Format::toString(i)+"] Or set All Ladder Level = 0 to not use this date for ladder sampling");
            
                    LdLevel[n].resize(LdGridSize);
                    LdLevelRaw[n].resize(LdGridSize);
                    LdCoupon[n].resize(LdGridSize+1);
                    for ( j= 0; j<LdGridSize; j++)
                    {// Store real Ladder data.
                        LdLevel[n][j] = tmpLdLevel[i][j];
                        LdLevelRaw[n][j] = tmpLdLevelRaw[i][j];
                        LdCoupon[n][j] = tmpLdCoupon[i][j];
                    }
                    //one more
                    LdCoupon[n][LdGridSize] = tmpLdCoupon[i][LdGridSize];
                    //LdCoupon[n][LdGridSize+1] = tmpLdCoupon[i][LdGridSize+1];

                    n ++;
                }

                // Libor is not always same as couopon payment dates.  
                // In Akasaka, it assume 0 date for Libor.  but I think it's only for VBA.
                // if from PYRAMID, it shouldn't work.
                /*
                LiborPayDates.push_back(inst->fixedLeg->PaymentDatesArray[i]);
                LiborAccrueDates.push_back(inst->liborLeg->PayDates[i]);
                if (i == coupNum-1)
                    LiborAccrueDates.push_back(inst->EndDates[i]);     // Accrue Date contains start date, end/start date and final end dates.
                */
                LiborAccrueDates = inst->liborLeg->AccrualDates;

                // Store Base coupons
                if (inst->isMonitorOnMax)
                    BaseCoupon.push_back(tmpLdCoupon[i][0]);
                else
                    BaseCoupon.push_back(tmpLdCoupon[i][rung_num]);
            }

            if (LdEndDates.size() != real_ld_num)
                throw ModelException(method,"Something wrong in Init Prod.");

            // mapping the pastValues to histMinOrMax
            histMinOrMax = inst->getHistMinOrMaxValues();

            // find out the starting idx. (no need past values for dummy Ladder period)
            int iHistStart = histMinOrMax.size()-real_ld_num;

            if (inst->valueDate > inst->fixedLeg->AccrueStartDates[0]){
                if (histMinOrMax.size()-iHistStart > real_ld_num)
                    throw ModelException(method,"pastValues don't contains the historical max/low");
            }

            // conditional coupons converted to KI coupon increment for pricing
            CouponAlreadyDetermined = 0.0;
            double discfact = 1.0;
            for (i=0; i<real_ld_num; i++)
            {// Final Price = - libor + BaseCoupon - Sum of LdCoupon with Barrier.
                for (j=0; j<rung_num; j++)
                {
                    if (inst->isMonitorOnMax)
                        LdCoupon[i][j] = -(LdCoupon[i][j+1] - LdCoupon[i][j]);      // adding coupon on base paymen for Step Up case 
                    else
                        LdCoupon[i][j] = LdCoupon[i][j+1] - LdCoupon[i][j];         // substracting coupon on base payment for Step Down case
                }

                LdCoupon[i].resize(rung_num); // last redundent coupon removed

                //Check the current MIN/MAX and add achived increment to Base Coupon.
                // <to Do>
                // Here, LdLevel should be LdLevelRaw if Ladder sampling was finished.
                // But need to be careful that this part treat the case that sampling is not finished case.
                if (inst->fixedLeg->AccrueStartDates[i] <= inst->valueDate && 
                    inst->valueDate < inst->fixedLeg->PaymentDatesArray[i] )
                {
                    //find out CurrentMinOrMax from pastValues
                    //exclude valDate = pastValues case.  Use current spot is safer as pastValues could have 0.
                    if (inst->valueDate > histMinOrMax[i+iHistStart].date)
                    {
                        CurrentMinOrMax = histMinOrMax[i+iHistStart].amount;
                        discfact = inst->discount->pv(LdCoupPayDates[i]);
                        for (j=0;j<rung_num;j++)
                        {
                            if (!inst->isMonitorOnMax && CurrentMinOrMax < LdLevel[i][j]*(1.0+FP_MIN))
                            {// step down 
                             //   BaseCoupon[i] -= LdCoupon[i][j];
                                CouponAlreadyDetermined -= discfact * LdCoupon[i][j];
                                LdCoupon[i][j] = 0.0;       // not a conditional coupon any more
                            }
                            else if (inst->isMonitorOnMax && CurrentMinOrMax > LdLevel[i][j]*(1.0-FP_MIN))
                            {// step up      
                             //   BaseCoupon[i] += LdCoupon[i][j];
                                CouponAlreadyDetermined += discfact * LdCoupon[i][j];
                                LdCoupon[i][j] = 0.0;       // not a conditional coupon any more
                            }
                        }
                    }
                }
            }

            // storage for inserted arrays since we reset every time in this model
            InStockArray.resize(2);
            InPriceArray.resize(2);
            // In***Array[][0] is actually not seems to be used in Akasaka.  So it would be good the size is just LdGridSize, no need to +1. (Keiji)
            InStockArray[0].resize(LdGridSize+1);
            InStockArray[1].resize(LdGridSize+1);
            InPriceArray[0].resize(LdGridSize+1);
            InPriceArray[1].resize(LdGridSize+1);

            CurrLadderPeriod = LdCoupPayDates.size()-1;     //period is start from 0 to size-1.
         }
         catch(exception& e)
         {// continue if indicative vol fail - non BS vol
             throw ModelException(e, method,"If the trade is already started, you may need give past value information.");
         }        
    }

    /** calculate at barriers for tree */
    virtual void preCalc(int step)
    {
        static const string method = "StepDownBondFDProd::preCalc";

        try 
        {
            int idx = tree1f->getSliceIndex(step);

          if (step != model->getLastStep())
          {
              // No need to calculate below if it is @ Maturity
             int iPrice;
             // determine if ladder change needed
             if (CurrLadderPeriod>0 && LdStartDates[CurrLadderPeriod]>model->getDate(step))
             {// shift ladder
                 CurrLadderPeriod --;
                 RefreshStatus = true; // waiting to refresh
             }

             // determine if ladder active
             if (   LdStartDates[CurrLadderPeriod] <= model->getDate(step) 
                    && LdEndDates[CurrLadderPeriod]  >= model->getDate(step))
             {
                LadderActive = true;
                if (RefreshStatus) // require coupon refreshing immediately when ladder becomes effective 
                    RefreshNow = true;
             }
             else
             {
                 // no ladder 
                  LadderActive = false;
                 // set inserted nodes to -1 if present so that they are not used at this step
                 if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION)
                 {
                     for (iPrice=0; iPrice<=LdGridSize; iPrice++)
                        InStockArray[idx][iPrice] = -1.0;
                 }
             }
          }
          
          DfToNextPayDate = inst->discount->pv(model->getDate(step), LdCoupPayDates[CurrLadderPeriod]);

          // to do: forward start adjustment
          // Prepare the adjusted ladder level taking care of discrete monitoring, treat for daily monitor
          if (!inst->IntraDayMonitor)
          {
              vector<double> vol;
              for (int iPrice =1; iPrice<=LdGridSize ; iPrice++)
              {
                  LdLevel[CurrLadderPeriod][iPrice-1] = LdLevelRaw[CurrLadderPeriod][iPrice-1];
                  tree1f->GetStepVol(step, vol, &LdLevel[CurrLadderPeriod][iPrice-1], 0, 0); // get vol at barrier
                  // adjust barrier if needed
                  if (inst->isMonitorOnMax)
                      Barrier::BarrierAdjustment(vol[0], true, LdLevel[CurrLadderPeriod][iPrice-1]);
                  else
                      Barrier::BarrierAdjustment(vol[0], false, LdLevel[CurrLadderPeriod][iPrice-1]);
              }
          }
        }
        catch (exception& e) 
        {
            throw ModelException(e, method);
        }
    }

    /** product payoff method at maturity */
    virtual void prod_BWD_T(const TreeSlice & spot,
                                  int step, 
                                  int bot, 
                                  int top, 
                                  int pStart, 
                                  int pEnd,
                                  const vector< TreeSliceSP > & price) 
    {
        static const string method = "CStepDownBondFDProd::prod_BWD_T";

        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        // For Insert Node case, nothing is done.  The judgement is not consistent and depend on numOfInsertNode!
        if (bot == 0 && top == 0)
            return;

        int iPrice, j;
        int currIdx = tree1f->GetSliceIdx();

        RefreshNow = false;
        RefreshStatus = false;
    
        // check Ladder
        if (   LdStartDates[CurrLadderPeriod] <= tree1f->getDate(step)
            && LdEndDates[CurrLadderPeriod]   >= tree1f->getDate(step) )
        {// now, ladder is active
            LadderActive = true;
        }
        else
        {
            throw ModelException(method, "Never reach here, Ladder Last Date should be maturity.");
            /*
            if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION)
            {
                for (iPrice=0; iPrice<=LdGridSize; iPrice++)
                    InStockArray[currIdx][iPrice] = -1.0;
            }*/
        }

        // calc payoff      
        for (iPrice=pStart; iPrice<=pEnd; iPrice++)
        {
            for (j = bot; j<= top; j++)      // Shoud I use BotClip/Top?
            {// floor and ceiling calculated
                if (iPrice == 0)
                    (p[iPrice])[j] = 0.0; // init Main Price [][0] array
                else if (  (!inst->isMonitorOnMax && s[j] < LdLevel[CurrLadderPeriod] [iPrice-1] *(1.0+ FP_MIN) )
                    ||( inst->isMonitorOnMax && s[j] > LdLevel[CurrLadderPeriod][iPrice-1] *(1.0- FP_MIN)) )
                {
                    (p[iPrice])[j] = DfToNextPayDate*LdCoupon[CurrLadderPeriod][iPrice-1];
                }
                else 
                    (p[iPrice])[j] = 0.0;
            }
            // price at barrier
            if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION && iPrice != 0)
            {
                InStockArray[currIdx][iPrice] = LdLevel[CurrLadderPeriod][iPrice-1];
                InPriceArray[currIdx][iPrice] = DfToNextPayDate*LdCoupon[CurrLadderPeriod][iPrice-1];
            }
        }
    }

    /** product payoff method at steps earlier than maturity */
    void prod_BWD(const TreeSlice & spot,
                         int step, 
                          int bot, 
                        int top, 
                         int pStart, 
                        int pEnd,
                        const vector< TreeSliceSP > & price) 
    {
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        int iPrice, j;

        // calc conditional payments
        for (iPrice=1; iPrice<=pEnd; iPrice++)
        {
            if (RefreshNow)
            {// subtract conditional payments
                for (j = bot; j<=top; j++)
                {
                    (p[0])[j] -= (p[iPrice])[j];            
                    (p[iPrice])[j] = 0.0;
                }
            }                    
            for (j = bot; j<=top; j++)
            {
                if (LadderActive && (  (!inst->isMonitorOnMax && s[j]<LdLevel[CurrLadderPeriod][iPrice-1]*(1.0+FP_MIN))
                                     ||( inst->isMonitorOnMax && s[j]>LdLevel[CurrLadderPeriod][iPrice-1]*(1.0-FP_MIN)) ))
                {//Checking Ladder is touched or not.
                    (p[iPrice])[j] = DfToNextPayDate*LdCoupon[CurrLadderPeriod][iPrice-1];
                }
            }
        }
    
        if(RefreshNow)
        {// refreshed
            RefreshNow = false;
            RefreshStatus = false;
        }

    }

    virtual string getCcyTreatment() const
    { 
        return inst->ccyTreatment;
    }

    /** premium scaling */
    virtual double scalePremium(vector<double> & P,
                                YieldCurveConstSP disc)
    {                
        int iPrice;
        double result = 0.0;

        // add already determined coupon
        result += CouponAlreadyDetermined;

        // combine result at step 0
        if (inst->SwitchPlainSwap)
        {
            result += inst->liborLeg->getPV(inst->valueDate, inst->discount.get());
            result += inst->fixedLeg->getPV(inst->valueDate, inst->discount.get());                  
        }

        if (inst->LinearCoupon)
        {// 2nd order Richardson extrapolation used
            for (iPrice=1; iPrice<=LdGridSize; iPrice+=6)
            {
                P[0] -= 3 * P[iPrice+1] - 3 * P[iPrice+2] + 6 * P[iPrice+3] - 3 * P[iPrice+4] + 3 * P[iPrice+5];                
            }
        }
        else
        {
            for (iPrice=1; iPrice<=LdGridSize; iPrice++)
            {
                P[0] -= P[iPrice];                
            }
        }

        P[0] += result;

        return P[0];
    }

    /** extra output requests */
    void recordOutput(Control* control, 
                      YieldCurveConstSP disc, 
                      Results* results)
    {
        // get prices at t=0
        int size = slices.size();
        vector< double > price0( size );
        for( int i = 0; i < size; ++i )
            price0[i] = model->getPrice0( *slices[i] );

        // save price
        double price = scalePremium(price0, disc);
        results->storePrice(price, disc->getCcy());

        // throw this back to the instrument itself
        inst->recordOutputRequests(control, results, price);
    }

    /** extra output requests */
    virtual bool positivePayoff() 
    {
        return Positive();
    }
    virtual bool Positive() 
    { 
        return false; 
    } 

    void SetInsertNodeAndPrice(int idx, int insPt, double insStock, int insPriority, int pStart, int pEnd, double insPrice)
    {
        if (insPt >= numIns) {
            throw ModelException("CStepDownBondFDProd::SetInsertNodeAndPrice", 
                                 "trying to set node beyond bound");
        }

        (*insNodes)[idx][insPt] = insStock;
        insOrders[idx][insPt] = insPriority;
        for (int i=pStart; i<=pEnd; i++) {
            (*insPrices)[idx][i][insPt] = insPrice;
        }
    }

    /** To make it possible to chnage InsertNode level, for Tree */
    virtual bool moveInsertNode(int currStep, int iPrice)
    {
        int currIdx = tree1f->getSliceIndex(currStep);

        if (iPrice == 0)
        {// for iPrice=0, there is no insert node.
            SetInsertNodeAndPrice(1-currIdx, 0, -100.0, 0, iPrice, iPrice, 9999999.9999);
        }
        else //if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION )//&& iPrice !=0 )
        {
            // copy previous insert node.
            SetInsertNodeAndPrice(1-currIdx, 0, InStockArray[1-currIdx][iPrice], 0, iPrice, iPrice, InPriceArray[1-currIdx][iPrice]);
            if (LadderActive)
            {// reset current array and keep it.
                InStockArray[currIdx][iPrice] = LdLevel[CurrLadderPeriod][iPrice-1]; // this is set for each ladder in effect
                InPriceArray[currIdx][iPrice] = DfToNextPayDate * LdCoupon[CurrLadderPeriod][iPrice-1];
                SetInsertNodeAndPrice(currIdx, 0, InStockArray[currIdx][iPrice], 0, iPrice, iPrice, InPriceArray[currIdx][iPrice]);
            }
        }

        return true;
    }

    void update(int & step, 
                FDProduct::UpdateType type)
    {
        // we assume just need one und level for spot here
        const TreeSlice & s = payoffIndex->getValue( step );
        int bot, top;
        s.getCalcRange( bot, top );

        const vector< TreeSliceSP > & price = slices;
        int pStart = 0, pEnd = price.size() - 1;

        if (type == FDProduct::BWD_T)
        {
            prod_BWD_T(s,
                       step,
                       bot,
                       top,
                       pStart, 
                       pEnd,
                       price);

            //insert nodes
            if (tree1f && tree1f->NumOfInsertNode>0)
            {
                prod_BWD_T(*insNodes,
                            step,
                            0,
                            tree1f->NumOfInsertNode-1,
                            pStart, 
                            pEnd,
                           *insPrices);
            }
        }
        else if(type == FDProduct::BWD)
        {
            prod_BWD(s,
                     step,
                     bot,
                     top,
                     pStart, 
                      pEnd,
                     price);

            //insert nodes
            if (tree1f && tree1f->NumOfInsertNode>0)
            {
                prod_BWD(*insNodes,
                          step,
                          0,
                           tree1f->NumOfInsertNode-1,
                          pStart, 
                          pEnd,
                         *insPrices);
            }
        }   
    }

protected:
 
    //All Ld variable has differnt array size of origianl input.  
    //Input Array contains dammy Date and Level if there is no ladder sampling.
    //At InitProd, those variable are resized for true lambda sampling size * ladder level.
    vector<vector<double> > LdLevel;        // Sorted Matrix by time by Ladder Level from LadderLevelArray
    vector<vector<double> > LdCoupon;       // Sorted Matrix by time by Ladder Level from LadderLevelArray  
    vector<vector<double> > LdLevelRaw;     // Sorted Matrix by time by Ladder Level from LadderLevelArray : To store original value
    int                     LdGridSize;     // Grid Size of LdLevel/Coupon.  This should be same num as NumOfPrice-1 in Tree

    CDateTimeArray          LdCoupPayDates; // Sorted along to Ladder Payment Dates.  Input Array include both of Libor and Ladder coupon dates.
    CDateTimeArray          LdStartDates;   // Ladder Start Time.  trimed all other dammy date from input array.
    CDateTimeArray          LdEndDates;   // Ladder  End  Time.  trimed all other dammy date from input array.
    
    double  CouponAlreadyDetermined;    // A coupon sum amount already touched ladder amount.
    vector<double> BaseCoupon;  // Unconditional Coupon to be payed.
    int  CurrLadderPeriod;   // Index of current Ladder Period.  0 to LdCoupPayDates-1
    bool LadderActive; // false(default) = off, true =on
    bool RefreshStatus, RefreshNow; // conditional coupon refresh

    CDateTimeArray           LiborAccrueDates;   //This array is extracted from From liborLeg->PayDates @ InitProd.
    CDateTimeArray           LiborPayDates;      //This array is extracted from From liborLeg->PayDates @ InitProd.

    vector<vector<double> > InStockArray;    // To change insert node depend on Price layer.
    vector<vector<double> > InPriceArray;       // To change insert node depend on Price layer.

private:
    const CStepDownBond* inst;
    // calc once in precalc to save time
    double DfToNextPayDate;
    double CurrentMinOrMax;                 //current min or max value to determine Ladder coupon

    CashFlowArray histMinOrMax;                // historical Min or Max array.

};

FDProductSP CStepDownBond::createProduct(FDModel* model) const
{
    return FDProductSP( new CStepDownBondFDProd(this, model) );
}

DRLIB_END_NAMESPACE
