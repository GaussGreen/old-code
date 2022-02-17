//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ExtendableNote.hpp
//
//   Description   Step down bond instrument. 
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/ExtendableNote.hpp"

DRLIB_BEGIN_NAMESPACE

// helpers
class CExtendableNoteHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(CExtendableNote, clazz);
        SUPERCLASS(CDblBarrier);
        EMPTY_SHELL_METHOD(defaultExtendableNote);
        
        FIELD(SpotFixing,         "Initial Spot Price to decide number of contract for last put/call option.");
        
        FIELD(fixedLeg, "cash flow of fixed payment");
        FIELD(liborLeg, "payment stream for Float Leg");

        // If the barrier are not on date case, Coupon Lock Date can be used
        FIELD(CouponLockDate, "All coupon before this date never KO even there is barrier.");
        FIELD_MAKE_OPTIONAL(CouponLockDate);    

        // Optionality to choice the output
        FIELD(SwitchPlainSwap, "Price with PlainSwap or Not.  Default = True (with PlainSwap).");
        FIELD_MAKE_OPTIONAL(SwitchPlainSwap);
        //  When SwitchPlainSwap = false, There is no KNOWN_CASH_FLOW info in this product.
        
        FIELD(RemoveSettlementOption, "default = false. If final put/fwd not needed, make it true");
        FIELD_MAKE_OPTIONAL(RemoveSettlementOption);    
        FIELD(RemoveCoupon, "default = false. If only need final put/fwd value, make it true");
        FIELD_MAKE_OPTIONAL(RemoveCoupon);    

        FIELD(isEquityPayer, "true : Pay equity payoff of Option/Fwd for. false(default) : Receive");
        FIELD_MAKE_OPTIONAL(isEquityPayer);    
        FIELD(isLiborPayer, "Recieve Fixed Leg for false.");
        FIELD_MAKE_OPTIONAL(isLiborPayer);    
    }

static IObject* defaultExtendableNote(){
        return new CExtendableNote();
    }
};

CClassConstSP const CExtendableNote::TYPE = CClass::registerClassLoadMethod(
    "ExtendableNote", typeid(CExtendableNote), CExtendableNoteHelper::load);

bool   CExtendableNoteLoad() {
    return (CExtendableNote::TYPE != 0);
}


// for reflection 
CExtendableNote::CExtendableNote(CClassConstSP clazz): CDblBarrier(clazz) {
    SwitchPlainSwap = true;
    RemoveSettlementOption = false;
    RemoveCoupon = false;
    CouponLockDate.rollDate(0);
    isEquityPayer = true;
    isLiborPayer = true;
}

// constructor
CExtendableNote::CExtendableNote(): CDblBarrier(TYPE)
{
    SwitchPlainSwap = true;
    RemoveSettlementOption = false;
    RemoveCoupon = false;
    CouponLockDate.rollDate(0);
    isEquityPayer = true;
    isLiborPayer = true;
};

/*--------------------------------------------



// **** now a copy of Vanilla, to do : add more validations */
void CExtendableNote::Validate()
{
    static const string method = "CExtendableNote::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object

    CDblBarrier::Validate();

    // Overwrite isCall
/*    if (PayoffMode == "PUT" || PayoffMode == "FORWARD")
        isCall = false;
    else
        isCall = true;
*/
	if (UpperBarrier->length() == 0 && LowerBarrier->length() == 0)     // Do I need it?
        throw ModelException(method, "no barrier schedule supplied.");

    const DateTime& matDate= exerciseSchedule->lastDate();
    int numItems = liborLeg->getSize();
    if (liborLeg->AccrualDates[numItems-1] > matDate)
        throw ModelException(method, "LiborAccrualDates["+Format::toString(numItems)+"] is later than maturity");
	numItems = fixedLeg->getSize();
    if (fixedLeg->AccrueStartDates[numItems-1] > matDate)
        throw ModelException(method, "FixedAccrueStartDates["+Format::toString(numItems)+"] is later than maturity");

    DateTime paymentDate = instSettle->settles(matDate, asset.get());    
    if (paymentDate < fixedLeg->PaymentDatesArray[numItems-1])
        throw ModelException(method, "Last PaymentDatesArray is later than paymentDate in settlement information");

    // to do:  Richardson iteration numbers.
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool CExtendableNote::sensShift(Theta* shift)
{    
    // roll today 
    DateTime rollDate = shift->rollDate(valueDate);
    // Need to give fixing rate if valuation date is fixing dates.
    liborLeg->setFixingforThetaShift(valueDate,discount.get(),rollDate);
    
    valueDate = rollDate;
    return true;
};


void CExtendableNote::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    
    int i, j;
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime matDate= exerciseSchedule->lastDate();
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

        //calculate Swap related Value and Store Values
        //libor
        double legResult;
        CashFlowArrayConstSP CflL = liborLeg->getCashFlowArray(valueDate, discount.get());
        CashFlowArraySP cfl (new CashFlowArray(0));
        if (CflL->size()>0)
       {
            for (i=0; i<liborLeg->getSize(); i++)
                    cfl->push_back((*CflL)[i]);
		    results->storeGreek(cfl, Results::DEBUG_PACKET, OutputNameSP(new OutputName("liborLeg")));
        }
        legResult = liborLeg->getPV(valueDate, discount.get());
        results->storeScalarGreek(legResult, Results::DEBUG_PACKET, OutputNameSP(new OutputName("liborLegValue")));

        //fixed cash flow
        CashFlowArrayConstSP CflF = fixedLeg->getCashFlowArray();
        CashFlowArraySP cff (new CashFlowArray(0));
        if (CflF->size()>0)
        {
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

            CashFlowArraySP finalCFL(new CashFlowArray(0));;                              
            // if there is known equity payoff, need to handle it.
            if (valueDate>=matDate) {
                double strike = exerciseSchedule->lastValue();
                double payout = notional/SpotFixing * (strike - spotAtMaturity);
                CashFlow final(instSettle->settles(matDate,asset.get()),payout);
                finalCFL->push_back(final);
            }


            // Pick up only known Fixed Leg.
            CashFlowArraySP knownCflF (new CashFlowArray(0));
            if (CflF->size()>0 && SwitchPlainSwap)
            {
                for (j=0; j<fixedLeg->getSize(); j++)
                {
                    if(valueDate > fixedLeg->AccrueStartDates[j])
                        knownCflF->push_back((*CflF)[j]);
                }
            }
            // Pick up only known Libor Leg.
            CashFlowArraySP knownCflL (new CashFlowArray(0));
            if (CflL->size()>0 && SwitchPlainSwap)
            {
                for (i=0; i<liborLeg->getSize(); i++) {
                    if(valueDate > liborLeg->AccrualDates[i])
                        knownCflL->push_back((*CflL)[i]);
                }
            }

            // now glue them all together
            CashFlowArraySP merge1(CashFlow::merge(knownCflL, knownCflF));
            CashFlowArraySP knownCfl(CashFlow::merge(merge1, finalCFL));

            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    knownCfl.get());   
        }

        // Add BARRIER_LEVEL
        request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (request && !fwdStarting) {
            // report barrier levels over a date range
            DateTime upperDate = BarrierLevel::barrierWindow(valueDate);

            BarrierLevelArraySP levels(new BarrierLevelArray(0));

            if (LowerBarType != "NA") {
                // use economic barrier (if it exists)
                Schedule* s = lowerEcoBarrier.get() ? lowerEcoBarrier.get(): 
                    LowerBarrier.get();
                CashFlowArraySP subset(s->subset(valueDate, upperDate));
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(false,(*subset)[i].date,(*subset)[i].amount);
                    levels->push_back(bl);
                }
            }

            if (UpperBarType != "NA") {
                // use economic barrier (if it exists)
                Schedule* s = upperEcoBarrier.get() ? upperEcoBarrier.get(): 
                    UpperBarrier.get();
                CashFlowArraySP subset(s->subset(valueDate, upperDate));
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(true,(*subset)[i].date,(*subset)[i].amount);
                    levels->push_back(bl);
                }
            }

            if (!levels->empty()) {
                OutputRequestUtil::recordBarrierLevels(control,
                                                       results,
                                                       asset->getTrueName(),
                                                       levels.get());
            }
        }            

    }
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool CExtendableNote::priceDeadInstrument(CControl* control, CResults* results) const
{
    double    value         = 0.0;
    DateTime  exerDate;

    static string method = "CExtendableNote::priceDeadInstrument";
    if (IntraDayMonitor && !fwdStarting)
    {
        double s = asset->getSpot();

        // check breach of barriers on value date for KO
		if ((UpperBarType == "KO") && (UpperBarBreached == false)
            && UpperBarrier->getDates()[0] <= valueDate
            && UpperBarrier->lastDate() >= valueDate) // barrier started and not finished.
		{
			if(UpperBarrier->getInterp() == "N")
			{
               for(int i = 0; i < UpperBarrier->length(); i++)
			   {
				   if (UpperBarrier->getDates()[i].equals(valueDate,false))
				   {
					   if (s > UpperBarrier->interpolate(valueDate)*(1.0-FP_MIN))
					   {	  
					       UpperBarBreached = true;
			               UpperBarBreachDate = valueDate;
					   }
					   break;
				   }
			   }
			}
			else
			{
               if (s > UpperBarrier->interpolate(valueDate)*(1.0-FP_MIN))
			   {
                       UpperBarBreached = true;
			           UpperBarBreachDate = valueDate;
			   }
			}
		}
		if ((LowerBarType == "KO") && (LowerBarBreached == false)
            && LowerBarrier->getDates()[0] <= valueDate
            && LowerBarrier->lastDate() >= valueDate) // barrier started and not finished.
		{
			if(LowerBarrier->getInterp() == "N")
			{
               for(int i = 0; i < LowerBarrier->length(); i++)
			   {
				   if (LowerBarrier->getDates()[i].equals(valueDate,false))
				   {
					   if (s < LowerBarrier->interpolate(valueDate)*(1.0+FP_MIN))
					   {
					       LowerBarBreached = true;
			               LowerBarBreachDate = valueDate;
					   }
					   break;
				   }
			   }
			}
			else
			{
               if (s < LowerBarrier->interpolate(valueDate)*(1.0+FP_MIN))
			   {
                   LowerBarBreached = true;
			       LowerBarBreachDate = valueDate;
			   }
			}
		}
    }  

    if (isExercised && ((LowerBarBreached  && LowerBarType == "KO") ||
                        (UpperBarBreached  && UpperBarType == "KO")))
        throw ModelException(method, 
                             "Can not breach KO barrier and exercise.");


    DateTime matDate = exerciseSchedule->lastDate();

    bool lowerOut = LowerBarBreached  && LowerBarType == "KO";
    bool upperOut = UpperBarBreached  && UpperBarType == "KO";

    bool isOut = ((lowerOut || upperOut) && BarrierDependence != "TWO_TOUCH")
                 || (lowerOut && upperOut);

    DateTimeArray payDates = fixedLeg->PaymentDatesArray;
    //DoubleArray part= ParticipationSchedule->getValues();
    //DateTimeArray exeDates = exerciseSchedule->getDates();
    //DoubleArray strikes = exerciseSchedule->getValues();
        
    //bool isDead = valueDate >= matDate || (isExercised && canExerciseEarly) || isOut;
    bool isDead = valueDate >= matDate || isOut;
    if (!isDead)
        return false; // not dead yet

    DateTime settlementDate = instSettle->settles(matDate, asset.get());

    int k = -9; // initialize with Neighbour failer case.

    // check already settled before final settlement of Strucutore
    if (isOut)
    {//find out next payment dates
        if((k = Neighbour(lowerOut ? LowerBarBreachDate : UpperBarBreachDate,
                            payDates, 0, payDates.size()-1, 1)) >= 0){
            settlementDate = payDates[k];
        }
        //if cannot find, use settlemenDate of whole structure.
    }

    if (valueDate >= settlementDate)
    {// settled already
        results->storePrice(0.0, discount->getCcy());
        addOutputRequests(control, results, 0.0, 0.0);
        return true;
    }

    // sort out ko case first
    if (isOut) // two touch case has problem choosing rebate - use lower for now
    {        
        if (settlementDate >= valueDate) // should really be > but then can't record cash flow if no delay
        {
            if (lowerOut && !!LowerRebate)
                value = LowerRebate->interpolate(LowerBarBreachDate);
            else if (upperOut && !!UpperRebate)
                value = UpperRebate->interpolate(UpperBarBreachDate);
        }
        else{// already terminated.  Just return 0.0.
            results->storePrice(0.0, discount->getCcy());
            return true;
        }
    }
    else if (settlementDate < valueDate){
        results->storePrice(0.0, discount->getCcy());
        return true;
    }
    else if (valueDate>=matDate) 
    {// if there is known equity payoff, need to handle it.
        double strike = exerciseSchedule->lastValue();
        if (isEquityPayer)
            value = -notional/SpotFixing * (strike - spotAtMaturity);
        else
            value =  notional/SpotFixing * (strike - spotAtMaturity);
    }

    // pv from settlement to today
    value *= discount->pv(valueDate, settlementDate);

    // add Swap Leg (Pay FixLeg / Recieve Libor Leg)
    if (k>=0){
        if (valueDate == payDates[k]){          //if valueDate = payDate, it's ignored.
            CashFlowArrayConstSP libor = liborLeg->getCashFlowArray(valueDate, discount.get());
            value += (*libor)[k].amount;  // Just Cash Flow of Today.            
            value += fixedLeg->CouponAmounts[k];
        }
        value += fixedLeg->getPV(valueDate ,valueDate ,payDates[k], discount.get());
        value += liborLeg->getPV(valueDate ,valueDate ,payDates[k] , discount.get());
    }
    else{// if k<0, then should be valueDate != settlementDate.
        value += fixedLeg->getPV(valueDate ,valueDate ,settlementDate , discount.get());
        value += liborLeg->getPV(valueDate ,valueDate ,settlementDate , discount.get());
    }

     if (!SwitchPlainSwap)
        value = 0.0;    //  If substruct PlainSwap, then there is no KNOWN_CASH_FLOW.

    // store results
    results->storePrice(value, discount->getCcy());
    addOutputRequests(control, results, value, 0.0);

    return true;
}

/*******************************************************************************************************************************/

/** create a fd payoff product **/
FDProductSP CExtendableNote::createProduct(FDModel* model) const
{
    return FDProductSP( new CExtendableNoteFDProd(this, model) );
}

void CExtendableNoteFDProd::recordOutput(Control*          control, 
										 YieldCurveConstSP disc, 
										 Results*          results)
{
    // get prices at t=0
    // save price
    double price = scalePremium(model->getPrice0( *slices[0] ), model->getPrice0( *slices[1] ), disc);   
    results->storePrice(price, disc->getCcy());

    if (control && control->isPricing()) 
	{
        DateTime       matDate = inst->exerciseSchedule->lastDate();
        double         indVol;
        // calculate indicative vol
        try 
		{
            if ( matDate.isGreater(inst->valueDate) )
            {

                DateTime imntStartDate = inst->fwdStarting? 
                    inst->startDate: inst->valueDate;

                // get vol request
                CVolRequestConstSP lnVolRequest = GetLNRequest();

                // interpolate the vol
                CVolProcessedSP  vol(inst->asset->getProcessedVol(lnVolRequest.get()));
                // cast to the type of vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                // this should never happen if our get market data has worked properly
                if (!vol)
				{
                    throw ModelException("CExtendableNoteFDProd::recordOutput", 
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

		// Store Settlement Option Price
        results->storeScalarGreek(SettleOptionPrice, Results::DEBUG_PACKET, 
                                    OutputNameSP(new OutputName("SettleOption")));    
        // Store Conditional Leg
        results->storeScalarGreek(KOSwapPrice, Results::DEBUG_PACKET, 
                                    OutputNameSP(new OutputName("KOSwap")));             

        inst->addOutputRequests(control,
                                results,
                                price,
                                indVol);
	}
}

/** premium scaling */
double CExtendableNoteFDProd::scalePremium(double P0, 
										   double P1,
										   YieldCurveConstSP disc)										  
{
	static const string method = "CExtendableNoteFDProd::ScalePremium";

    SettleOptionPrice = P1;     //Store the settlement option price
    KOSwapPrice = P0;   

	// combine results if needed
	if (!inst->RemoveSettlementOption) // combine with option
	{
        if (inst->isEquityPayer)
        {
			P0 -= P1;		
            SettleOptionPrice *= -1.0;
        }
        else
        {
		    P0 += P1;
        }
	}

	// substract PlainSwap
    if (!inst->SwitchPlainSwap)
    {
        P0 -= inst->liborLeg->getPV(inst->valueDate, inst->discount.get());
        P0 -= inst->fixedLeg->getPV(inst->valueDate, inst->discount.get());        
    }

	return P0;
}

void CExtendableNoteFDProd::initProd()
{
    static const string method = "CExtendableNoteFDProd::InitProd";
    try
    {    
        DblBarrierFDProd::initProd();    
        CalcStepCoupon();
     }
     catch(exception& e)
     {// continue if indicative vol fail - non BS vol
         throw ModelException(e, method,"Hey Keiji.  Your Init Prod is rubish.  Check it Again!!");
     }
}

void CExtendableNoteFDProd::CalcStepCoupon()
{
	int j, k;
    int lastStep = model->getLastStep();
    
    CashFlowArrayConstSP libor = inst->liborLeg->getCashFlowArray(inst->valueDate, inst->discount.get());    
    CashFlowArrayConstSP fixedflow = inst->fixedLeg->getCashFlowArray();

	Coupon.resize(lastStep+1);
	double ExCoupon = 0.0;

	for (j=0; j<=lastStep; j++)
		Coupon[j] = 0.0; // init
    if (!inst->RemoveCoupon)
    {
        int numFixedCoup = (*fixedflow).size();
        // fixed leg
        for(j=0; j< numFixedCoup; j++) // assuming all arrays are the same size
	    {
		    if ( (*fixedflow)[j].date > inst->valueDate)
		    { // only coupons that pay date end in the future are included
			    if ((k = Neighbour(inst->fixedLeg->AccrueStartDates[j], model->getDates(), 0, lastStep, -1)) > 0)
			    {
				    if (inst->fixedLeg->AccrueEndDates[j] > inst->CouponLockDate)
					    Coupon[k] += (*fixedflow)[j].amount*inst->discount->pv(model->getDate(k), (*fixedflow)[j].date);
				    else if(model->getDate(k) <= inst->CouponLockDate)		//coupon before Rock Date are included in ExCoupon   keiji
					    ExCoupon += (*fixedflow)[j].amount*inst->discount->pv(model->getDate(0), (*fixedflow)[j].date);                    
			    }
			    else // PV'ed to step 1 (not 0) for theta calc
				    ExCoupon += (*fixedflow)[j].amount*inst->discount->pv(model->getDate(0), (*fixedflow)[j].date);                
		    }
        }
		// floating leg
        int numFloatCoup = (*libor).size();
        for(j=0; j< numFloatCoup; j++)
        {
            if ( (*libor)[j].date > inst->valueDate)
		    { // only coupons that pay date end in the future are included
			    if ((k = Neighbour(inst->liborLeg->AccrualDates[j], model->getDates(), 0, lastStep, -1)) > 0)
			    {
				    if (inst->liborLeg->AccrualDates[j] > inst->CouponLockDate) 	
					    Coupon[k] += (*libor)[j].amount * inst->discount->pv(model->getDate(k),(*libor)[j].date);
				    else if(model->getDate(k) <= inst->CouponLockDate)		//coupon before Rock Date are included in ExCoupon   keiji
					    ExCoupon += (*libor)[j].amount * inst->discount->pv((*libor)[j].date);
                }
			    else // PV'ed to step 1 (not 0) for theta calc
				    ExCoupon += (*libor)[j].amount * inst->discount->pv((*libor)[j].date);
		    }
        }
    }

    // ex coupon is added to coupon at first step
    Coupon[0] += ExCoupon;
}

/** product payoff method at maturity */
void CExtendableNoteFDProd::prod_BWD_T(const TreeSlice & spot,
									         int step, 
											 int bot, 
											 int top, 
											 int pStart, 
											 int pEnd,
											 const vector< TreeSliceSP > & price)
{
    static const string method = "CExtendableNoteFDProd::prod_BWD_T";
    
    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );

    int j;
    double claim = (inst->isCall? 1.0:-1.0);
    
    // calc terminal price at T first
	
	int lastStep = model->getLastStep();
	
    double strikeT = stepStrike[lastStep]; //  Suppose this is Strike.Exercise.Interpolate(TradeTime[NumOfStep]);
	for (j = bot; j<=top; j++)
	{// floor and ceiling calculated
		if ((UType == KO && s[j] >= UpperBar * (1.0- FP_MIN)) 
				|| (LType == KO && s[j] <= LowerBar * (1.0+ FP_MIN))) // dead
		{
			(p[1])[j] = (p[0])[j] = 0.0;
			(p[2])[j] = (p[3])[j] = 0.0;		//keiji
		}
		else
		{			
            // geared settlement option  (EuroPrice)
            (p[3])[j] = claim*(s[j]-strikeT);
			if (inst->PayoffMode == "CALL" || inst->PayoffMode == "PUT")
                (p[3])[j] = Maths::max((p[3])[j],0.0); // option settlement
            (p[3])[j] *= inst->notional/inst->SpotFixing; // Leveleged.
            
            (p[2])[j] = Coupon[lastStep];
			if ((UType == KI && s[j] < UpperBar * (1.0+FP_MIN))		//<-actulally this is not needed
					|| (LType == 1 && s[j] > LowerBar * (1.0-FP_MIN)))
				(p[1])[j] = (p[0])[j] = 0.0;		//failed to KI
			else
			{
				(p[1])[j] = (p[3])[j];		
				(p[0])[j] = (p[2])[j];
			}
		}
    }      
}

/** product payoff method at steps earlier than maturity */
void CExtendableNoteFDProd::prod_BWD(const TreeSlice & spot,
											int step,
										    int bot,
											int top,
											int pStart,
										    int pEnd,
										    const vector< TreeSliceSP > & price)
{
	static const string method = "CExtendableNoteFDProd::prod_BWD";

    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );

	for (int j = bot; j<=top; j++)
	{
		// barrier pay off
		if((UType == KO && s[j] > UpperBar * (1.0- FP_MIN)) || (LType == KO && s[j] < LowerBar * (1.0+ FP_MIN))) // dead
		{
			(p[3])[j] = (p[2])[j] = 0.0;	//keiji
		}
		else
		{
			(p[2])[j] += Coupon[step];
		}

		if (s[j] < UpperBar * (1.0- FP_MIN) && s[j] > LowerBar * (1.0+ FP_MIN)) // still within barriers
		{	
			if (LType != KI)	//with DI, there is no coupon	keiji
				(p[0])[j] += Coupon[step];
		}
		else 
		{
			(p[0])[j] = (p[2])[j];
			(p[1])[j] = (p[3])[j];
		}
	}
}

/** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
void CExtendableNoteFDProd::update(int& step, FDProduct::UpdateType type)
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
		    prod_BWD_T(   *insNodes,
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

DRLIB_END_NAMESPACE
