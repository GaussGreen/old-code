//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : Q3MQQuasiPricer.cpp
//
//   Description : Pricer using numerical multiQ under adjusted forward measure.
//             
//   Author      : Mehdi Chaabouni
//
//   Date        : january 2006
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/Q3MQQuasiPricer.hpp"
#include "edginc/IDiscountCurve.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/CDSVolRequestSimpleEuropean.hpp"
#include "edginc/CDSVolProcessedSimpleEuropean.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/ICDS.hpp"
#include "edginc/IFixedRateCreditFeeLeg.hpp"
#include "edginc/Addin.hpp"  
#include "edginc/Settlement.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Range.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"

#include "edginc/Integrator.hpp"


DRLIB_BEGIN_NAMESPACE

#define Q3_BDRY 1E-2
#define Q3_TINY 1E-6
#define Q3BOND_PRICE(c,f,m,y) \
         (fabs(y)<Q3_TINY) ? (1 + (m) * (c)) : ((c) / (y) + \
         (1- (c) / (y)) * pow(1 + (y)/(f), -(m) * (f)))

//const string Q3MQQuasiPricer::CALL_OPTION = "CALL_OPTION";
//const string Q3MQQuasiPricer::PUT_OPTION = "PUT_OPTION";
  




// Constructors:
Q3MQQuasiPricer::Q3MQQuasiPricer() : CInstrument(TYPE) {
};

Q3MQQuasiPricer::Q3MQQuasiPricer(DateTime valueDate,
				 DateTime resetDate,
				 DateTime startDate,
				 DateTime maturityDate,
				 DateTime payDate,
				 VanillaCDSSP underlying,
				 bool isFloatRisky,
				 bool isFixed,
				 double fixing,
				 double meanReversion,
				 bool isCashSettle,
				 MultiQDistributionSP multiq,
				 ICDSParSpreadsConstSP        crv,
				 YieldCurveConstSP            rfCrv,
				 int RichardsonPolyOrder) :
  CInstrument(TYPE),
  valueDate(valueDate),
  resetDate(resetDate),
  startDate(startDate),
  maturityDate(maturityDate),
  payDate(payDate),
  cds(underlying),
  isFloatRisky(isFloatRisky),
  isFixed(isFixed),
  fixing(fixing),
  meanReversion(meanReversion),
  isCashSettle(isCashSettle),
  mq(multiq),
  RichardsonPolyOrder(RichardsonPolyOrder),
  calibrated(false),
  crv(crv),
  rfCrv(rfCrv)
{
  Validate();
}


double Q3MQQuasiPricer::adjRate()
{
  if (isFixed) return fixing;

  Validate();
  return faMq.forward();
}

double Q3MQQuasiPricer::optionPrice(bool isCall, double strike)
{
  if (isFixed) return (isCall ? (max(fixing - strike,0.0)) : (max(strike - fixing,0.0)));

  Validate();
  return faMq.vanillaOptionPrice(isCall,strike);
}




 
/************************************************************************
******** CInstrument ***************************************************
*************************************************************/
//GetMarket
void Q3MQQuasiPricer::GetMarket(const IModel* model, const CMarketDataSP market)
{
//  market->GetReferenceDate(valueDate);

//     /*=========================================================================
//      * GET THE DISCOUNT CURVE 
//      *=======================================================================*/
//     discount = cds->getYieldCurveWrapper();
//     discount.getData(model, market);
//     rfCrv = YieldCurveSP(discount.get());

//     /*=========================================================================
//      * GET THE CREDIT SPREAD CURVE
//      *=======================================================================*/
//     cdsParSpreads = cds->getParSpreadsWrapper();
//     ICDSParSpreads::getMarketData(
//         model,
//         market.get(),
//         discount.getName(),
//         cdsParSpreads);
//     crv = ICDSParSpreadsSP(cdsParSpreads.get());

}

// discountYieldCurveName
string Q3MQQuasiPricer::discountYieldCurveName() const
{
  return "not Implemented" ; //discount.getName();
}
//getValueDate
DateTime Q3MQQuasiPricer::getValueDate() const
{
  return valueDate;
}

//Validate
void Q3MQQuasiPricer::Validate() 
{
  static const char* method = "Q3MQQuasiPricer::Validate";

  // if already calibrated do nothing.
  if (calibrated) return;

    // Require a model for calculating fees
    // so just provide the default closed form
    ClosedFormForwardRatePricerSP cfPricer =
        ClosedFormForwardRatePricerSP(
            new ClosedFormForwardRatePricer());
    // retrieve the pre-payment curve
    IDecretionCurveConstSP prepay = crv->getPrepayCurve();

  //otherwise calibrate
  // Compute mutable internal parameters:
   volStart = 0;
   expiry = valueDate.yearFrac(resetDate);
   rateStart = valueDate.yearFrac(startDate);
   delayZeroMaturity = startDate.yearFrac(payDate);
   convexityZeroMaturity = startDate.yearFrac(maturityDate) ;
   swapMaturity = startDate.yearFrac(maturityDate);
   
   delayZeroRate =(delayZeroMaturity !=0.0 ?
		   - log(crv->survivalProb(startDate,payDate))/delayZeroMaturity
		   : 0.0);
   convexityZeroRate = (convexityZeroMaturity !=0.0 ?
			-log(crv->survivalProb(startDate,maturityDate))/convexityZeroMaturity
			: 0.0);

   
   
   double fwdProtVal = cds->getContingentLegPV(startDate,*crv)
       / cds->getNotional();

    //the fee leg must be a fixed fee type for this to work
    ICreditFeeLegSP feeLeg = cds->getFeeLeg();
    IFixedRateCreditFeeLeg* frcfl = dynamic_cast<IFixedRateCreditFeeLeg*>(feeLeg.get());
    if (!frcfl)
    {
	    throw ModelException(method,
			                 "Fee leg must be a fixed rate type");
    }
    DateTime earliestRiskyDate = startDate;
    DateTime latestRiskyDate = feeLeg->getLastPayDate();
   double fwdAnnuityVal = feeLeg->getFeeLegPV(startDate,earliestRiskyDate,latestRiskyDate,*crv,*crv,prepay,true,cds->getAccrualDcc(),cfPricer)
       / (cds->getNotional() * frcfl->getRate());
     
     if (fwdAnnuityVal <=0)
       {
	 throw ModelException(method,
			      "CDS must have a positive annuity");
       }
     else 
       {
	 swapFwdSpread = fwdProtVal / fwdAnnuityVal;                             
       }
     convexityRecovery = crv->getRecovery();
     fwdAnnuity = fwdAnnuityVal;
     swapFreq = double(crv->getSwapFrequency()); // double(cds->paymentFreq->annualFrequency());
     MaturityPeriod maturityPeriod((int)swapFreq);
     swapRate = swapFwdSpread / (1. - convexityRecovery) +
       rfCrv->couponRate(startDate,maturityDate, maturityPeriod,false,DayCountConventionFactory::make("Act/360"));
     vnfmFwdAnnuity = (1.0 -  crv->pv(startDate,maturityDate))/swapRate;
     CDSVolRequestSimpleEuropean volReq(
					swapFwdSpread,swapFwdSpread,true, 
					true, startDate,maturityDate, 
					true);
     CDSVolProcessedSimpleEuropeanSP 
       procVolSP(dynamic_cast<CDSVolProcessedSimpleEuropean*>(crv->getProcessedVol(&volReq))); 
     sigATM = procVolSP->atmVolatility();
   
     
   convexityRiskFreeRate = (convexityZeroMaturity !=0.0 ?
			    - log(rfCrv->pv(startDate,maturityDate))/convexityZeroMaturity
			    : 0.0);

VNFMSP vnfm = VNFMSP(new VNFMCMCDS(volStart,
								expiry,
								rateStart,
								swapMaturity,
								swapFreq,
								meanReversion,
								isCashSettle,
								swapRate,
								sigATM,
								vnfmFwdAnnuity,
								swapFwdSpread,
								convexityRecovery));

vnfm->Q3VNFMZero2Swap(convexityZeroMaturity,
					  convexityZeroRate,
					  convexAlpha,
					  convexPower,
					  convexT);


if (isFloatRisky)
{
vnfm->Q3VNFMZero2Swap(delayZeroMaturity,
					  delayZeroRate,
					  delayAlpha,
					  delayPower,
					  delayT);
} else
{
	delayAlpha = 0.0;
	delayPower = 1.0;
}
 
 //Use convexity and delay adjustments to compute the new adjusted distribution			   
 faMq =  FAMultiQDistribution(mq,
			      convexAlpha,
			      convexPower,
			      convexityRecovery,
			      convexityRiskFreeRate,
			      convexityZeroMaturity,
			      delayAlpha,
			      delayPower,
			      delayZeroMaturity,
			      RichardsonPolyOrder);

  
  calibrated = true;
  return;

}
//Not used for the moment.
double Q3MQQuasiPricer::convexityAndDelayAdjustmentFunction(double y) const
{
 double yCutoff = Q3_BDRY * mq->forward();
    double yCut = (y > Maths::max(Q3_TINY,yCutoff) ? y : yCutoff);
    
    //Calculate delay adjustment
    double delayLambda;
    double delayAdjustment;
    if (y < Q3_TINY)
      {
	delayAdjustment = 1.0;
      }else{
     delayLambda = delayAlpha * pow(yCut,delayPower);
     delayAdjustment = pow(1+delayLambda,-delayZeroMaturity);
     double  p = Maths::min(y/yCutoff,1.0);
     delayAdjustment = p * delayAdjustment + (1-p);
    }

    //Calculate convexity adjustment
  
    double convexityLambda = convexAlpha * pow(yCut,convexPower);
    double convexityAdjustment = (1-convexityRecovery) / yCut
      * convexityLambda /(convexityLambda + convexityRiskFreeRate)
      * (1- pow((1+convexityLambda)*(1+convexityRiskFreeRate),-convexityZeroMaturity));

    //Return adjustement function
    return  delayAdjustment / convexityAdjustment;
}

			       
		

// Loading
class Q3MQQuasiPricerHelper{
public:
static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Q3MQQuasiPricer, clazz);
	SUPERCLASS(CInstrument);
	EMPTY_SHELL_METHOD(defaultQ3MQQuasiPricer);
	//Fields
	FIELD(valueDate,"");
	FIELD(resetDate,"");
	FIELD(startDate,"");
	FIELD(maturityDate,"");
	FIELD(payDate,"");
	FIELD(cds,"");
	FIELD_MAKE_OPTIONAL(cds);
	FIELD(isFloatRisky,"");
	FIELD(isFixed,"");
	FIELD(fixing,"");
	FIELD(meanReversion,"");
	FIELD(isCashSettle,"");
	FIELD(mq,"");
	FIELD(crv,"");
	FIELD(rfCrv,"");
	FIELD(RichardsonPolyOrder,"");

    }

	    static IObject* defaultQ3MQQuasiPricer() {
        return new Q3MQQuasiPricer();
    }
};

//Loading
CClassConstSP const Q3MQQuasiPricer::TYPE = CClass::registerClassLoadMethod("Q3MQQuasiPricer", 
									    typeid(Q3MQQuasiPricer),
									    Q3MQQuasiPricerHelper::load);

// for linker
bool   Q3MQQuasiPricerLoad() {
    return (Q3MQQuasiPricer::TYPE != 0);
   }


/*=============================================================================
 * Class to provide add-in functions - values and option prices
 *===========================================================================*/
class Q3MQQuasiPricerAddin2 : public CObject {
public:
    static CClassConstSP const TYPE;

    // addin parameters
    Q3MQQuasiPricerSP q3;
    double value;
    string callOrPut;

    double forward() {
    return q3->adjRate();
  }

    double vanillaPrice() {
        bool isCall = false;
        if (callOrPut.size()<1) 
            throw ModelException("Q3MQQuasiPricerAddin2::vanillaPrice",
            "callOrPut may not be blank and must begin with 'C' for call or 'P' for put.");
        char cp = toupper(callOrPut[0]);
        if (cp!='C' && cp!='P')
            throw ModelException("Q3MQQuasiPricerAddin2::vanillaPrice",
            "callOrPut must begin with 'C' for call or 'P' for put.");

        if (cp=='C') isCall = true;

        return q3->optionPrice(isCall, value);
    };

    
    Q3MQQuasiPricerAddin2() : CObject(TYPE), value(0), callOrPut("") {};

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(Q3MQQuasiPricerAddin2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultQ3MQQuasiPricerAddin2);
        FIELD(q3, "Q3MQQuasiPricer Object");
        FIELD(value,"Input value or strike for addin functions");
        FIELD(callOrPut,"Should begin with 'C' for call or 'P' for put.");
        FIELD_MAKE_OPTIONAL(callOrPut);

       

        Addin::registerDoubleMethod("Q3MQQUASI_VANILLA_PRICE",
            Addin::UTILITIES,
            "Prices a vanilla call or put with given strike under the forward adjusted measure.",
            &Q3MQQuasiPricerAddin2::vanillaPrice);

	Addin::registerDoubleMethod("Q3MQQUASI_FORWARD",
				    Addin::UTILITIES,
				    "Returns the expectation of the rate",
				    &Q3MQQuasiPricerAddin2::forward);
    }
  
  static IObject* defaultQ3MQQuasiPricerAddin2() {
        return new Q3MQQuasiPricerAddin2();
  }
};

CClassConstSP const Q3MQQuasiPricerAddin2::TYPE =
  CClass::registerClassLoadMethod("Q3MQQuasiPricerAddin2", typeid(Q3MQQuasiPricerAddin2), Q3MQQuasiPricerAddin2::load);


DRLIB_END_NAMESPACE


