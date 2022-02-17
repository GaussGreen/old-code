//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CMCDSFloatingLeg.cpp
//
//   Description : Methods to price a CMCDS floating leg
//
//   Author      : Mehdi Chaabouni
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/CMCDSFloatingLeg.hpp"
#include "edginc/IDiscountCurve.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/Q3MQQuasiPricer.hpp"
#include "edginc/CDSVolRequestSimpleEuropean.hpp"
#include "edginc/CDSVolProcessedSimpleEuropean.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/ICDS.hpp"
#include "edginc/IFixedRateCreditFeeLeg.hpp"
#include "edginc/Settlement.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/DayCountConventionFactory.hpp"

#ifndef EDG_CMCDSFloatingLeg_CPP
#define EDG_CMCDSFloatingLeg_CPP

DRLIB_BEGIN_NAMESPACE

/*Make sure the class links*/
bool CMCDSFloatingLegLoad(){
	return (CMCDSFloatingLeg::TYPE!=0);
}

/*Default Constructor, used by reflection*/
CMCDSFloatingLeg::CMCDSFloatingLeg() : Generic1FactorCredit(TYPE),
				       RichardsonPolyOrder(9),
				       isCashSettle(false),
				       floatAccrualOnDefault(false)
				       
{
  rollDates = DateTimeArray();
  beforeRollSpreads = DoubleArray();
  afterRollSpreads = DoubleArray();
  cap = DoubleArray();
  floor = DoubleArray();
  shift = DoubleArray();
};

string  CMCDSFloatingLeg::discountYieldCurveName() const{
  return discount.getName();
        }
        
/** Check the inputs to validate product*/
void CMCDSFloatingLeg::Validate() {
  
  static const string method = "CMCDSFloatingLeg::validate";
     
  //Check that arrays have same size
  int size = ulResets.size();
  if (size ==0)
	  throw ModelException(method,"Reset Dates must have a positive size");
  if (ulStartDate.size()!=size ||
      floatPayDate.size()!=size ||
      ulMaturityDate.size()!=size ||
      floatAccStart.size()!=size ||
      floatAccEnd.size()!=size)
    {
    throw ModelException(method,"Dates sizes must be the same");
     }

  // Check beforeRoll and afterRoll have the same size
  if (beforeRollSpreads.size() != afterRollSpreads.size()
      ||beforeRollSpreads.size() != rollDates.size() )
    throw ModelException(method,"rollDates, beforeRollSpreads and afterRollSpreads  must have the same size");

  if (beforeRollSpreads.size() !=0)
  {
	  for (int j = 0; j < beforeRollSpreads.size() ; j++)
  {
	  if (rollDates[j] > valueDate)
			throw ModelException(method,"rollDates must be in the past");
	  if (beforeRollSpreads[j] < 0.0 || afterRollSpreads[j] < 0.0)
			throw ModelException(method,"beforeRollSpreads and afterRollSpreads  must be positive");
  }
  }
  if (cap.size() != floor.size())
    throw ModelException(method,"cap and floor arrays  must have the same size");

  if (cap.size() != 0 && cap.size() != size)
    throw ModelException(method,"cap array must have the same size as the reset dates or be omitted");

  if (shift.size() != 0 && shift.size() != size)
    throw ModelException(method,"shift array must have the same size as the reset dates or be omitted");

  if (shift.size() == 0)
    shift = DoubleArray(size,0.0);
  
  //Check order of dates and cap >= floor
  int i = 0;
  for (i=0; i< size;i++)
    {
    if (ulMaturityDate[i] <= ulStartDate[i])
      throw ModelException(method,"Maturities must be after start dates");
    if (ulResets[i].date > ulStartDate[i])
      throw ModelException(method,"Reset Dates must be before start dates");
    if (floatPayDate[i] < ulResets[i].date)
      throw ModelException(method,"Float Pay Dates must be after reset dates");
	if ( i > 0 && ulStartDate[i] < ulStartDate[i-1])
		throw ModelException(method,"Start Dates must be in increasing order");

    if (cap.size() ==size)
      {
    if (cap[i] < floor[i])
      throw ModelException(method,"Cap levels must be greater than (or equal to) the floor ones");
      }
    }
 
  //Build underlying if not input
  if (!ulPayAccOnDefault) ulPayAccOnDefault.reset(CBool::create(true));
  if(!ulCouponFreq)   ulCouponFreq.reset(CInt::create(cdsParSpreads->getSwapFrequency()));
  DayCountConventionSP ulCouponDccObj;
  if (ulCouponDcc.empty()) ulCouponDccObj =DayCountConventionSP(cdsParSpreads->dayCountConv());
  else ulCouponDccObj = DayCountConventionSP(DayCountConventionFactory::make(ulCouponDcc));
  ulAccDcc = DayCountConventionSP(ulCouponDccObj.clone());
  ulCouponBdc = BadDayConventionSP(BadDayConventionFactory::make("F"));
  ulAccBdc = BadDayConventionSP(BadDayConventionFactory::make("F"));
  ulCouponHol=HolidayWrapper(cdsParSpreads->getHolidays().clone());
  ulAccHol=HolidayWrapper(cdsParSpreads->getHolidays().clone());
  if (ulStub.empty())  ulStub = "sf";
  ulRecovery = cdsParSpreads->getRecovery();
      
      if (!underlying)
	{
	  underlying = VanillaCDSSP( new VanillaCDS(valueDate,
						ulStartDate[0],              //protection start date 
						ulMaturityDate[0],        //maturity date i.e. protection end date 
						true,          //knock out cds if the underlying credit defaults before issue date. 
						discount,      //discount curve 
						cdsParSpreads,//risky curve 
						    ulPayAccOnDefault->boolValue(),              //pay accrued fee on default or not 
						1.0,                  //fee 
						    MaturityPeriodSP(new MaturityPeriod(ulCouponFreq->intValue())),    //payment frequency 
						ulCouponDccObj,     //payment day count convention 
						ulAccDcc, //accrual day count convention 
						ulCouponBdc,       //payment bad day convention 
						ulAccBdc,   //accrual bad day convention 
						ulCouponBdc, //valuation bad day convention 
						ulCouponHol,           //payment holiday(s) 
						ulAccHol,       //accrual holiday(s) 
						ulStub,            //sf (short front stub), lf(long front stub) 
						ulRecovery,             //recovery rate 
						0.0,                    //delay between default date and settlement upon default 
						SettlementSP(new RollingSettlement())));  
  
    }
  else
    {
      // if the underlying is input, we need to validate it.
      underlying->Validate();
    }
}        
//GetMarket
void CMCDSFloatingLeg::GetMarket(
     const IModel*        model, 
     const CMarketDataSP  market)
{
  //Get market data for underlying

    market->GetReferenceDate(valueDate);
 
   
    /*=============================================================
     *Get discount Curve
     *==============================================================*/
    discount.getData(model, market);

    /*==================================================================
     *Get Credit Spread Curve
     *==================================================================*/
      cdsParSpreads.getData(model,market);
    
    //underlying->GetMarket(model,market);
    // call instrument specific getMarket routine if applicable
    if (IGetMarket::TYPE->isInstance(this))
      {
        IGetMarket* imnt = dynamic_cast<IGetMarket*>(this);
        imnt->getMarket(model,market.get());
      }

    
}



//GetValueDate
DateTime CMCDSFloatingLeg::getValueDate() const 
{
    return valueDate;
}

//computePastJumps
double CMCDSFloatingLeg::computePastJumps(DateTime d) const
{
  static const string method = "CMCDSFloatingLeg::computePastJumps";

  if (beforeRollSpreads.size() == 0) return 0.0;
  
  double pastJumps = 0.0;
  int i =0;
  for (i=0; i < beforeRollSpreads.size() ; i++)
    {
      if ( rollDates[i] < d)
      pastJumps += afterRollSpreads[i] - beforeRollSpreads[i];
    }
  return pastJumps;
  
}
              
//Price
void CMCDSFloatingLeg::quasiPrice(
			     CResults* results,
			     Control* control,
			     const CModel* model) const
   
{
     static const string method = "CMCDSFloatingLeg::quasiPrice";
     try {

        // retrieve the model to be used in the calculation
        // of fee leg forward rates
        IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>
            (const_cast<Model*>(model));
        if (!ihfrp)
        {
            throw ModelException(method,
                "Model must implement IHasForwardRateModel");
        }
        IForwardRatePricerSP frModel = ihfrp->getForwardRatePricer();

       OutputRequest* request = 0;
       OutputRequest* requestKcf = 0;
       const string& ccy = cdsParSpreads->getCcy();
       double recovery = cdsParSpreads->getRecovery();
    IDecretionCurveConstSP prepay = cdsParSpreads->getPrepayCurve();
       
     //Price and outputs to store here
     double totalValue = 0.0;
     DoubleArray fwdRates(ulStartDate.size(),0.0);
     DoubleArray adjRates(ulStartDate.size());
     DoubleArray floorPVs(ulStartDate.size());
     DoubleArray capPVs(ulStartDate.size());
     DoubleArray survProbs(ulStartDate.size(),1.0);
     DoubleArray discFactors(ulStartDate.size(),1.0);
     DoubleArray matAdjs(ulStartDate.size(),0.0);
     DoubleArray fwAdjs(ulStartDate.size(),0.0);
     DoubleArray accValues(ulStartDate.size(),0.0);
     DoubleArray atmVols(ulStartDate.size(),0.0);


     
     requestKcf = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
               
     //This is the incremental maturity adjustment due to the composite fixing.
     double pastJumps = computePastJumps(valueDate);
     double adjustmentMaturity = pastJumps;
     double forwardSpread = 0.0;
     double forwardSpread1 = 0.0;
     double forwardSpread2 = 0.0;
     
     //
     DayCountConventionSP floatDccObj(DayCountConventionFactory::make(floatDCC));
     bool isDefaulted = cdsParSpreads->defaulted();
      
     // Base CDS we will use to build a CDS with new dates in the loop.
     ICDSConvention * changeableCDS =
       dynamic_cast<ICDSConvention *>(underlying.get());

     if (isDefaulted)
       {
	 /* For now, return 0 is defaulted
	    Need better tackle this. Generate the accrual on default.
	    Not sure what to do if  flag is defaulted for an index*/
	 results->storePrice(0.0,ccy);
	 return;
       }
     
    for(int i=0 ; i < ulStartDate.size() ; i++)
    {
      double value =0.0;
      double adjRate =.0;
      double floorPV = .0;
      double capPV = .0;
      
      // year fraction applying to the rate
      double yearFraction = floatDccObj->years(floatAccStart[i],floatAccEnd[i]);
      
      if (ulResets[i].date < valueDate)
	{
	  //We observe a fixing which we move by a shift and the past jumps then we apply the floor and the cap
	  adjRate =  ulResets[i].amount;
	  if (cap.size() != 0)
	    capPV = Maths::max( adjRate + shift[i]- computePastJumps(ulResets[i].date)-cap[i],0.0);
	  if (floor.size() != 0)
	    floorPV = Maths::max(floor[i] - ( adjRate + shift[i]- computePastJumps(ulResets[i].date)),0.0);
	  value = adjRate + shift[i]+ floorPV - capPV - computePastJumps(ulResets[i].date) ;

	  //The cash flow is known if it is in the past or if it is in the future and not risky.
	  if (requestKcf && ( (valueDate >= floatAccEnd[i]) || !isFloatRisky ))
	     {
	       CashFlow cf(floatPayDate[i],
			   value * yearFraction * floatNotional[i]);
	       if (!kcf) kcf=CashFlowArraySP(new CashFlowArray(1,cf));
	       else kcf->push_back(cf);
	     }
	} else
	{
      double adjustmentForward=0.0;
      /*=========================================================================     
       * Build CDS to be observed with the right dates and compute forward
       *==========================================================================*/
            
      ICDSSP extendedCDS = changeableCDS->generateCDS(ulStartDate[i],
						      ulMaturityDate[i],
						      1);

      IBadDayAdjuster* bdAdj = (IBadDayAdjuster*)dynamic_cast<IBadDayAdjuster*>(extendedCDS.get());
      if (!bdAdj)
      {
        throw ModelException(method, "extendedCDS must implement IBadDayAdjuster");
      }
     IBadDayAdjusterConstSP bda(IBadDayAdjusterConstSP::attachToRef(bdAdj));
     ICreditContingentLegSP ctgLeg = extendedCDS->getContingentLeg();
     double fwdProtVal = ctgLeg->getContingentLegPV(ulStartDate[i],*(cdsParSpreads.get()),bda)
       / extendedCDS->getNotional();

    //fee leg must be of a fixed rate variety for this to work
    ICreditFeeLegSP feeLeg = extendedCDS->getFeeLeg();
    IFixedRateCreditFeeLeg* frcfl = dynamic_cast<IFixedRateCreditFeeLeg*>(feeLeg.get());
    if (!frcfl)
    {
        throw ModelException(method, "Fee leg is not of a fixed fee variety "
                                     "(must implement IFixedRateCreditFeeLeg)");
    }

    DateTime earliestRiskyDate = ctgLeg->firstObservationStartDate();
    DateTime latestRiskyDate = ctgLeg->lastObservationEndDate();
    double fwdAnnuityVal = feeLeg->getFeeLegPV(ulStartDate[i],earliestRiskyDate,latestRiskyDate,*(cdsParSpreads.get()),*(cdsParSpreads.get()),prepay,true,extendedCDS->getAccrualDcc(),frModel)
       / (extendedCDS->getNotional() * frcfl->getRate());
     
     if (fwdAnnuityVal <=0)
       {
	 throw ModelException(method,
			      "CDS must have a positive annuity");
       }
     else 
       {
	 forwardSpread = fwdProtVal / fwdAnnuityVal;
	 fwdRates[i] = forwardSpread;
       }
     
    /*=========================================================================     
     *Compute Maturity adjustment increment - Deterministic approximation
     *==========================================================================*/
     if (isMaturityAdjusted && i>0)
       {
     // calculate the forward spread for a single name for the reduced maturity
	 ICDSSP CDSForMatAdj1 = changeableCDS->generateCDS(ulStartDate[i-1],
							 ulMaturityDate[i-1],
							1);
     IBadDayAdjuster* bdAdj1 = (IBadDayAdjuster*)dynamic_cast<IBadDayAdjuster*>(CDSForMatAdj1.get());
      if (!bdAdj1)
      {
        throw ModelException(method, "CDSForMatAdj1 must implement IBadDayAdjuster");
      }
      IBadDayAdjusterConstSP bda1(IBadDayAdjusterConstSP::attachToRef(bdAdj1));

    ICreditContingentLegSP adj1CtgLeg = CDSForMatAdj1->getContingentLeg();
	 double fwdProtVal1 = adj1CtgLeg->getContingentLegPV(ulStartDate[i-1],*(cdsParSpreads.get()),bda1)
	   /CDSForMatAdj1->getNotional();

    //fee leg must be of a fixed rate variety for this to work
    ICreditFeeLegSP adj1FeeLeg = CDSForMatAdj1->getFeeLeg();
    IFixedRateCreditFeeLeg* adj1frcfl = dynamic_cast<IFixedRateCreditFeeLeg*>(adj1FeeLeg.get());
    if (!adj1frcfl)
    {
        throw ModelException(method, "Fee leg is not of a fixed fee variety "
                                     "(must implement IFixedRateCreditFeeLeg)");
    }
    earliestRiskyDate = adj1CtgLeg->firstObservationStartDate();
    latestRiskyDate = adj1CtgLeg->lastObservationEndDate();
    double fwdAnnuityVal1 = adj1FeeLeg->getFeeLegPV(ulStartDate[i-1],earliestRiskyDate,latestRiskyDate,*(cdsParSpreads.get()),*(cdsParSpreads.get()),prepay,true,CDSForMatAdj1->getAccrualDcc(),frModel)
	   /(CDSForMatAdj1->getNotional() * adj1frcfl->getRate());
     
     if (fwdAnnuityVal1 <=0)
       {
	 throw ModelException(method,
			      "CDS must have a positive annuity");
       }
     else 
       {
	 forwardSpread1 = fwdProtVal1 / fwdAnnuityVal1;                             
       }
           // calculate the forward spread for a single name for the reduced maturity
	 ICDSSP CDSForMatAdj2 = changeableCDS->generateCDS(ulStartDate[i-1],
							 ulMaturityDate[i],
							1);
     IBadDayAdjuster* bdAdj2 = (IBadDayAdjuster*)dynamic_cast<IBadDayAdjuster*>(CDSForMatAdj2.get());
      if (!bdAdj2)
      {
        throw ModelException(method, "CDSForMatAdj2 must implement IBadDayAdjuster");
      }
     IBadDayAdjusterConstSP bda2(IBadDayAdjusterConstSP::attachToRef(bdAdj2));

     
     ICreditContingentLegSP adj2CtgLeg = CDSForMatAdj2->getContingentLeg();
	 double fwdProtVal2 = adj2CtgLeg->getContingentLegPV(ulStartDate[i-1],*(cdsParSpreads.get()),bda2)
	   /CDSForMatAdj2->getNotional();

    //fee leg must be of a fixed rate variety for this to work
    ICreditFeeLegSP adj2FeeLeg = CDSForMatAdj2->getFeeLeg();
    IFixedRateCreditFeeLeg* adj2frcfl = dynamic_cast<IFixedRateCreditFeeLeg*>(adj2FeeLeg.get());
    if (!adj2frcfl)
    {
        throw ModelException(method, "Fee leg is not of a fixed fee variety "
                                     "(must implement IFixedRateCreditFeeLeg)");
    }
    earliestRiskyDate = adj2CtgLeg->firstObservationStartDate();
    latestRiskyDate = adj2CtgLeg->lastObservationEndDate();
	 double fwdAnnuityVal2 = adj2FeeLeg->getFeeLegPV(ulStartDate[i-1],earliestRiskyDate,latestRiskyDate,*(cdsParSpreads.get()),*(cdsParSpreads.get()),prepay,true,CDSForMatAdj2->getAccrualDcc(),frModel)
	   /(CDSForMatAdj2->getNotional() * adj2frcfl->getRate());
     
     if (fwdAnnuityVal2 <=0)
       {
	 throw ModelException(method,
			      "CDS must have a positive annuity");
       }
     else 
       {
	 forwardSpread2 = fwdProtVal2 / fwdAnnuityVal2;                             
       }                   
     //Maturity adjustment
     adjustmentMaturity += (forwardSpread2 - forwardSpread1);
      
       }
     
    /*==========================================================================
     * No-Ko adjustment for indices - Deterministic approximation
     *===========================================================================*/
     if (isForwardAdjusted && ulResets[i].date > valueDate)
       {
     double survProbToReset = cdsParSpreads->survivalProb(valueDate, ulResets[i].date);
	 double survProbToStart  = cdsParSpreads->survivalProb(valueDate, ulStartDate[i]);
	 double dfToReset = discount->pv(valueDate, ulResets[i].date);
	 double dfToStart =discount->pv(valueDate, ulStartDate[i]);
	 
     adjustmentForward = (1-recovery) * (1-survProbToReset) * dfToReset
		 / (fwdAnnuityVal * survProbToStart * dfToStart);
	 //adjustmentForward = (1-recovery) * (1-survProbToReset)
		// / (fwdAnnuityVal);
       }
     
    /*===========================================================================
     * Fetch MultiQ and mean reversion from the market
     *===========================================================================*/
     CDSVolRequestSimpleEuropean volReq(
					forwardSpread, forwardSpread,true, 
					true, ulStartDate[i],ulMaturityDate[i], 
					true );
     CDSVolProcessedSimpleEuropeanSP 
       procVolSP(dynamic_cast<CDSVolProcessedSimpleEuropean*>(cdsParSpreads->getProcessedVol(&volReq))); 
     
     MultiQDistributionSP multiq = procVolSP->multiQ();
     atmVols[i] = procVolSP->atmVolatility();
     double beta = procVolSP->meanReversion();
       
    /*=========================================================================     
     * Build the Q3 Pricer object and compute expectations under new measure.
     *==========================================================================*/
     Q3MQQuasiPricer   quasiPricer(valueDate,
				   ulResets[i].date,
				   ulStartDate[i],
				   ulMaturityDate[i],
				   floatPayDate[i],
				   VanillaCDSSP(dynamic_cast<VanillaCDS *>(extendedCDS.get())),
				   isFloatRisky,
				   ulResets[i].date < valueDate,
				   ulResets[i].amount,
				   beta,
				   isCashSettle,
				   multiq,
				   cdsParSpreads.getSP(),
				   discount.getSP(),
				   RichardsonPolyOrder);
     
     // Expectation of the rate under FA measure
      adjRate = quasiPricer.adjRate();
      
      
      //Expectation of the floor option under FA measure
      if (floor.size() != 0)
	{
      double floorStrike = floor[i] - shift[i]  - 
                         (adjustmentForward - adjustmentMaturity);
      
      floorPV =  quasiPricer.optionPrice(false, floorStrike);
	}
      else floorPV =0.0;
      
      // Expectation of the cap option under FA measure
      if (cap.size() != 0)
	{
      double capStrike = cap[i] - shift[i] - 
                         (adjustmentForward - adjustmentMaturity);
                                    
      capPV = quasiPricer.optionPrice(true, capStrike);
	}
      else capPV = 0.0;
                                    
      //Bring together different components                                    
      value += adjRate + shift[i] + adjustmentForward 
	-adjustmentMaturity + floorPV - capPV;

      fwAdjs[i]=adjustmentForward;
      matAdjs[i]=adjustmentMaturity - pastJumps;
	}
      //Store expectations to be output
      adjRates[i]=adjRate;
      capPVs[i] = capPV;
      floorPVs[i]=floorPV;
      
      
    /*============================================================================
     * Calculate discounted value
     *===========================================================================*/
      
     
	
      double discountFactor =(floatPayDate[i] > valueDate ?
			      discount->pv(valueDate, floatPayDate[i]):
			      0.0);
     
      
     double survivalProbability = (isFloatRisky ?
				   cdsParSpreads->survivalProb(valueDate, floatPayDate[i]):
				   1.);
   
     
     value *= yearFraction * discountFactor * survivalProbability * floatNotional[i];

    /*=========================================================================     
     * Compute Accruals on default - Deterministic rough approximation
     * For the moment suppose no accrual on default as legal docs not clear
     * about this point and Trader does not know.
     *==========================================================================*/
//      double accDefProb = (floatAccrualOnDefault && isFloatRisky ?
// 			  - cdsParSpreads->survivalProb(valueDate,floatAccEnd[i]) +
// 			  (floatAccStart[i]>valueDate?cdsParSpreads->survivalProb(valueDate,floatAccStart[i]):1.0):
// 			  0.0);
//      double accYearFrac =  floatDccObj->years(floatAccStart[i]> valueDate ?floatAccStart[i]:valueDate,
// 						 floatAccEnd[i]);
//      double accValue = accDefProb * accYearFrac * .5 * value;
     double accValue = 0.0;
     value += accValue ;

     discFactors[i]=discountFactor;
     survProbs[i] = survivalProbability;
     accValues[i] = accValue;

     // Final incremental value
     totalValue += value;
     }
     
    /*==========================================================================
     * Outputting the results
     *============================================================================*/
    results->storePrice(totalValue,ccy);
    
    request = control->requestsOutput(OutputRequest::FORWARD_CDS_SPREAD);
    DoubleArraySP fr(new DoubleArray(fwdRates));
    if (request) {
      results->storeRequestResult(request, fr);
    }
    
    request = control->requestsOutput(OutputRequest::CMCDS_ADJUSTED_RATES);
    DoubleArraySP ar(new DoubleArray(adjRates));
    if (request) {
      results->storeRequestResult(request, ar);
    }

    request = control->requestsOutput(OutputRequest::CMCDS_CAP_PV);
    DoubleArraySP cpv(new DoubleArray(capPVs));
    if (request) {
      results->storeRequestResult(request, cpv);
    }

    request = control->requestsOutput(OutputRequest::CMCDS_FLOOR_PV);
    DoubleArraySP fpv(new DoubleArray(floorPVs));
    if (request) {
      results->storeRequestResult(request, fpv);
    }
 request = control->requestsOutput(OutputRequest::CMCDS_DISCOUNT_FACTORS);
    DoubleArraySP cfact(new DoubleArray(discFactors));
    if (request) {
      results->storeRequestResult(request, cfact);
    }
     request = control->requestsOutput(OutputRequest::CMCDS_SURV_PROBS);
    DoubleArraySP sprob(new DoubleArray(survProbs));
    if (request) {
      results->storeRequestResult(request, sprob);
    }
     request = control->requestsOutput(OutputRequest::CMCDS_FORWARD_ADJS);
    DoubleArraySP fadj(new DoubleArray(fwAdjs));
    if (request) {
      results->storeRequestResult(request, fadj);
    }
     request = control->requestsOutput(OutputRequest::CMCDS_MAT_ADJS);
    DoubleArraySP madj(new DoubleArray(matAdjs));
    if (request) {
      results->storeRequestResult(request, madj);
    }
     request = control->requestsOutput(OutputRequest::CMCDS_ACCS_ON_DEFAULT);
    DoubleArraySP aod(new DoubleArray(accValues));
    if (request) {
      results->storeRequestResult(request, aod);
    }
     request = control->requestsOutput(OutputRequest::IND_VOL);
    DoubleArraySP atv(new DoubleArray(atmVols));
    if (request) {
      results->storeRequestResult(request, atv);
    }
    //kKNOW_CASHFLOWS
    if (requestKcf) 
      OutputRequestUtil::recordKnownCashflows(control,
					      results,
					      ccy,
					      kcf.get());
	// PAYMENT_DATES
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
	  DateTimeArraySP dates(new DateTimeArray(floatPayDate));
	  OutputRequestUtil::recordPaymentDates(control,results,dates.get());
        }
	
     }catch (exception &e) {
       throw ModelException(&e,method);
     }
}


/*=============================================================================
 * Reflection, loading, etc.
 *===========================================================================*/
class CMCDSFloatingLegHelper {
public:
    static IObject* defaultCMCDSFloatingLeg();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CMCDSFloatingLeg, clazz);
        SUPERCLASS(Generic1FactorCredit);
        IMPLEMENTS(MultiQQuasiSmile::IIntoProduct);
	IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultCMCDSFloatingLeg);

        // Fields
        
	FIELD(  ulResets,"Cashflow array of when the spread will reset and the fixing rate if known");
	FIELD(  ulStartDate,"Array of start dates for the CDS's whose spread is going to be observed");
	FIELD(  ulMaturityDate,"Maturity dates of the CDS's whose spread is going to be observed");
	FIELD(  floatPayDate,"Array of dates specifying when the observed spread is paid");
	FIELD(  floatAccStart,"Array of dates specifying when the spread starts accruing");
	FIELD(  floatAccEnd,"Array of dates specifying when the spread accruing ends.");
	FIELD(      floatNotional,"Notional for each payment");
	FIELD(  floatDCC,"DCC for each payment. This together with accrual start and end will transform the rate into a payment.");
	FIELD(    isFloatRisky,"whether payment is contingent on default");
	FIELD(isForwardAdjusted,"whether to adjust the forward. Relevant for indices");
	FIELD(isMaturityAdjusted,"For composite fixing, whether to take maturity adjustments into account");
	
	
	//FIELD(    cdsParSpreads,"name of the spread curve");
	//FIELD(       discount,"name of the yield curve");

	//Transient Fields
	FIELD(kcf,"Known cashflows");
	FIELD_MAKE_TRANSIENT(kcf);

	//Optional Fields
	FIELD(rollDates,"Dates on which roll has happened. Only past dates will be taken into account.");
	FIELD(beforeRollSpreads,"for each roll date, the spread before the roll. This is used in computing the past jumps");
	FIELD(afterRollSpreads,"for each roll date, the spread after the roll. This is used in computing the past jumps");
	FIELD(   cap,"Array of caps on the observed spread");
	FIELD(   floor,"Array of floors on the observed spread ");
	FIELD(   shift,"Array of shifts to the observed Spread");
	FIELD(  floatAccrualOnDefault,"whether accrual on the floater is paid upon default. NOT USED as term sheet unclear.");
	FIELD(isCashSettle,"whether compute cvxity and delay adjustments using a flat or a curvey curve");
	FIELD(ulPayAccOnDefault,"Observed underlying CDS pays accrued on default on not. Default = true");
	FIELD(ulCouponFreq,"Observed underlying CDS coupon frequency. Default from market");
	FIELD(ulCouponDcc,"Observed underlying CDS coupon DCC. Default from market");
	FIELD(ulStub,"Observed underlying CDS stub. Default = sf");
	FIELD(     underlying,"Example CDS whose spread is going to be observed");
	FIELD(RichardsonPolyOrder,"Parameter for numerical integration. Defaulted to 9");

	FIELD_MAKE_OPTIONAL(underlying);
	FIELD_MAKE_OPTIONAL(ulPayAccOnDefault);
	FIELD_MAKE_OPTIONAL(ulCouponFreq);
	FIELD_MAKE_OPTIONAL(ulCouponDcc);
	FIELD_MAKE_OPTIONAL(ulStub);
	FIELD_MAKE_OPTIONAL(rollDates);
	FIELD_MAKE_OPTIONAL(beforeRollSpreads);
	FIELD_MAKE_OPTIONAL(afterRollSpreads);
	FIELD_MAKE_OPTIONAL(cap);
	FIELD_MAKE_OPTIONAL(floor);
	FIELD_MAKE_OPTIONAL(shift);
	FIELD_MAKE_OPTIONAL(isCashSettle);
	FIELD_MAKE_OPTIONAL(RichardsonPolyOrder);
	FIELD_MAKE_OPTIONAL(floatAccrualOnDefault);
    }

};


IObject* CMCDSFloatingLegHelper::defaultCMCDSFloatingLeg() {
    return new CMCDSFloatingLeg();
}


//TYPE
CClassConstSP const CMCDSFloatingLeg::TYPE = 
    CClass::registerClassLoadMethod("CMCDSFloatingLeg", typeid(CMCDSFloatingLeg),CMCDSFloatingLegHelper::load);


/*=============================================================================
 * Pricing Models
 *===========================================================================*/
 
 // MultiQ Quasi smile pricer
class CMCDSFloatingLegMultiQQuasi : public MultiQQuasiSmile::IProduct{
private:
    const CMCDSFloatingLeg* cmcds; // a reference
    const MultiQQuasiSmile* model;

public:
    CMCDSFloatingLegMultiQQuasi(const CMCDSFloatingLeg* cmcds, const MultiQQuasiSmile* model): 
                                      cmcds(cmcds), model(model) {}
    void price(MultiQQuasiSmile* model,
               Control*         control, 
               CResults*        results) const {
        cmcds->quasiPrice(results, control, model);
    }
};
    
MultiQQuasiSmile::IProduct* CMCDSFloatingLeg::createProduct(MultiQQuasiSmile* model) const{
    return new CMCDSFloatingLegMultiQQuasi(this, model);
}

/*=============================================================================
 * Theta::Shift Interface
 *===========================================================================*/
bool CMCDSFloatingLeg::sensShift(Theta* theta) {
    try {
        valueDate = theta->rollDate(valueDate);
    } catch (exception& e) {
        throw ModelException(e, "CMCDSFloatingLeg::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}


DRLIB_END_NAMESPACE

#endif

