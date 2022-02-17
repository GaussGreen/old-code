//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : Q3MQQuasiPricer.hpp
//
//   Description : Pricer using numerical multiQ under adjusted forward measure.
//             
//   Author      : Mehdi Chaabouni
//
//   Date        : january 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDG_Q3MQQuasiPricer_HPP
#define EDG_Q3MQQuasiPricer_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/MultiQDistribution.hpp"
#include "edginc/FAMultiQDistribution.hpp"
#include "edginc/VanillaCDS.hpp"
#include "edginc/VNFMCMCDS.hpp"

DRLIB_BEGIN_NAMESPACE


FORWARD_DECLARE(Settlement)
FORWARD_DECLARE(Schedule)
FORWARD_DECLARE(ICDS)
FORWARD_DECLARE(MultiQDistribution)
FORWARD_DECLARE(FAMultiQDistribution)
FORWARD_DECLARE(Expiry)
FORWARD_DECLARE(Control)
FORWARD_DECLARE(Q3MQQuasiPricer)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)
FORWARD_DECLARE_WRAPPER(YieldCurve)
  
struct FAData
{
  //Fields
  double zeroVol;
  double swapVol;
  double zeroSwapCorr;
  double alpha;
  double power;
  //should contain other distribution data (multiQ) ??? or just
  // how to change measures ?
  //Default Constructor
  FAData():
    zeroVol(0),
    swapVol(0),
    zeroSwapCorr(0),
    alpha(0),
    power(0) {}
  
  //Constructor
  FAData(double zeroVol,
	 double swapVol,
	 double zeroSwapCorr,
	 double alpha,
	 double power):
    zeroVol(zeroVol),
    swapVol(swapVol),
    zeroSwapCorr(zeroSwapCorr),
    alpha(alpha),
    power(power) {}

};


class PRODUCTS_DLL Q3MQQuasiPricer : public CInstrument
{
public:
  static CClassConstSP const TYPE;
  
  //Constructors
  Q3MQQuasiPricer();
  Q3MQQuasiPricer(DateTime valueDate,
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
		  MultiQDistributionSP mq,
		  ICDSParSpreadsConstSP        crv,
		  YieldCurveConstSP            rfCrv,
		  int RichardsonPolyOrder);
       
      
    /*======================================================================
     * CInstrument Implementations
     *======================================================================*/
    virtual void GetMarket(const IModel*, const CMarketDataSP);
    virtual  DateTime getValueDate() const;
    virtual void Validate();
    virtual string discountYieldCurveName() const;
  
    /*======================================================================
     * Pricing
     *=======================================================================*/
    // Get the FA adjusted rate
    double adjRate();
  
  //Get the expectation of payoff under FA measure
  // static const string CALL_OPTION;
  //static const string PUT_OPTION;
   double optionPrice(bool isCall, double strike);
  
protected:
  // internally calibrated members
  mutable double volStart;
  mutable double expiry;
  mutable double rateStart;
  mutable double delayZeroMaturity;
  mutable double convexityZeroMaturity;
  mutable double swapMaturity;
  mutable double swapRate;
  mutable double swapFwdSpread;
  mutable double delayZeroRate;
  mutable double convexityZeroRate;
  mutable double fwdAnnuity;
  mutable double vnfmFwdAnnuity;
  mutable double sigATM;
  mutable double swapFreq;
  mutable bool Q3VNFMZero2SwapInputsComputed;
  mutable FAMultiQDistribution faMq;
  mutable double convexityRecovery;
  mutable double convexityRiskFreeRate;

  double mutable convexAlpha;
  double mutable convexPower;
  double mutable convexT;
  double mutable delayAlpha; 
  double mutable delayPower;
  double mutable delayT;

  void computeQ3VNFMZero2SwapCRInputs();

private:
  
    /*=====================================================================
     * Fields
     *====================================================================*/   
  DateTime valueDate;  
  DateTime resetDate;
  /** start Date of the CDS*/
  DateTime startDate;
  /** maturity Date of the CDS*/
  DateTime maturityDate;
	/** PayDate is the date when the flow is paid*/
  DateTime payDate;
  /** This object englobes all the conventions related to the CDS apart from its start and end date*/
  VanillaCDSSP cds;
     
  /** Determines if we do the delay adjustment*/
  bool isFloatRisky;
  /** True if the the payment to be made is known and has already been fixed*/
  bool isFixed;
  /** only relevant if isFixed is TRUE*/
  double fixing;
  
  double meanReversion;
 bool isCashSettle;
  //Those will be taken from the market in the future. But need to change
  //  CDSVolCubeMultiQSmile
  MultiQDistributionSP mq;
  
  //Richardson Poly order passed to FAMultiQDistribution to do the numerical integration
  int RichardsonPolyOrder;
   
  /***************************************/
  //internally calibrated parameters
 
  //Identity function
  double identity (double x)
  { return x;}
  double convexityAndDelayAdjustmentFunction(double y) const;
 
  // replace by Validate()
  // void calibrate() ;
  mutable bool calibrated;
 

// these all come from the underlying discount and credit curves and
    // are always overwritten by those in the underlying.
   ICDSParSpreadsConstSP        crv;
    YieldCurveConstSP            rfCrv;
   //  ICDSParSpreadsWrapper           cdsParSpreads;    
//     YieldCurveWrapper               discount;    



  friend class Q3MQQuasiPricerHelper;
  friend class Q3MQQuasiPricerFunction;
 
};



typedef smartPtr<Q3MQQuasiPricer>    Q3MQQuasiPricerSP;

DRLIB_END_NAMESPACE
#endif

