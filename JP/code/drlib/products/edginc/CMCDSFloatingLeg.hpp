//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CMCDSFloatingLeg.hpp
//
//   Description : Constant Maturity CDS Floating Payment
//             
//   Author      : Mehdi Chaabouni
//
//   Date        : january 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDG_CMCDSFloatingLeg_HPP
#define EDG_CMCDSFloatingLeg_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/MultiQQuasiSmile.hpp"
#include "edginc/Model.hpp"
#include "edginc/Generic1FactorCredit.hpp"

DRLIB_BEGIN_NAMESPACE


FORWARD_DECLARE(Settlement)
  FORWARD_DECLARE(Schedule)
  FORWARD_DECLARE(DayCountConvention)
  FORWARD_DECLARE(ICDS)
  FORWARD_DECLARE(Expiry)
  FORWARD_DECLARE(Control)
  FORWARD_DECLARE_WRAPPER(ICDSParSpreads)
  FORWARD_DECLARE_WRAPPER(YieldCurve)
  
/**Class defining a Constant Maturity CDS. */
  class PRODUCTS_DLL CMCDSFloatingLeg :   public Generic1FactorCredit,
					  public virtual MultiQQuasiSmile::IIntoProduct,
					  public virtual Theta::Shift
{    
public:
    static CClassConstSP const TYPE;
    
     /*=========================================================================
     * I/CInstrument Interface
     *=======================================================================*/
    // Underlying Market
    virtual void GetMarket(const IModel* model, const CMarketDataSP);
    virtual void Validate();
    virtual DateTime getValueDate() const;
    virtual string discountYieldCurveName() const;
    
    /*=========================================================================
     * Pricing Model Interface Implementation
     *=======================================================================*/
    virtual MultiQQuasiSmile::IProduct* createProduct(MultiQQuasiSmile * model) const;
    
    void quasiPrice(
         CResults* results,
         Control* control,
         const CModel* model) const;
         
  /*=========================================================================
   * Tweak interface methods
   *=======================================================================*/
  bool sensShift(Theta* theta);



protected:
  CMCDSFloatingLeg();   

private:
  
  double computePastJumps(DateTime) const;
  
  /*=======================================================================
    * Underlying - Conventions for which spread to observe
    *=====================================================================*/
  ICDSSP underlying;
  CBoolSP ulPayAccOnDefault;
  CIntSP ulCouponFreq;
  string ulCouponDcc;
  mutable DayCountConventionSP ulAccDcc;
  mutable BadDayConventionSP ulCouponBdc;
  mutable BadDayConventionSP ulAccBdc;
  mutable HolidayWrapper ulCouponHol;
  mutable HolidayWrapper ulAccHol;
  string ulStub;
  mutable double ulRecovery;
        /*=========================================================================
         * Data Fields
          *=======================================================================*/
        //Dates
  //  DateTime valueDate; // Comes with Generic1FactorCredit
  /** Fixing dates and fixing values (if any)*/
  CashFlowArray ulResets;
  //  DateTimeArray ulResetDate;
  /** start Dates*/
  DateTimeArray ulStartDate;
  /** When the floating spread is paid*/
  DateTimeArray floatPayDate;
  /** Start date for the accrual*/
  DateTimeArray floatAccStart;
  /** End date for the accrual*/
  DateTimeArray floatAccEnd;
  /** End Dates of the observed CDSs*/
  DateTimeArray ulMaturityDate;
  
  
  //Rate calculation definition
  /** Payoff = Max(Min(observed Spread at reset + shift,  cap), floor) */
  DoubleArray cap;
  DoubleArray floor;
  DoubleArray shift;
  
  //Fixings
        /** an array of booleans to say if the floating has been fixed.
	    This usually happens if the observation date is in the past */
  //BoolArray isFixed;
  /** used only when the corresponding boolean isFixed is true */
  //DoubleArray fixing;
  
  //Past Events introducing parallel shifts on the fixing
  /** For indices, we need to keep track of names replacement
      We store the sum of changes at roll dates. This will be added to the payoff */
  //double pastIndexMaturityAdj;
  DateTimeArray rollDates;
  DoubleArray beforeRollSpreads;
  DoubleArray afterRollSpreads;
 
    
  //Floating amount conventions
  /** Once you know the expected rate under the forward adjusted measure,
      you need to multiply it by the DCC and the notional*/
  DoubleArray floatNotional;
  //  DayCountConventionArray floatingDCC;
  string floatDCC;
  bool floatAccrualOnDefault;
  
  
 
  //Flags
  bool isFloatRisky;
  bool isForwardAdjusted;
  bool isMaturityAdjusted;
  
   
 
  // defines the way to compute the delay and convexity adjustments
  bool isCashSettle;
  
  /** These will be used to determine the forward, the adjustments, the discounts...*/
  //These now come from the Generic1FactorCredit parent class.
  //ICDSParSpreadsWrapper cdsParSpreads;
  //YieldCurveWrapper discount;
  
  //Cash Flow array for KNOWN_CASHFLOWS
  mutable CashFlowArraySP kcf;
  
  // Parameter for accuracy of numerical integration. Defaulted to 9.
  int RichardsonPolyOrder;
     /*=========================================================================
      * Friends
      *=======================================================================*/
  friend class CMCDSFloatingLegHelper;
  friend class CMCDSFLoatingLegMultiQQuasi;
  
  /*=========================================================================
   * For Reflection
     *=======================================================================*/
  static void load(CClassSP& clazz);
  static IObject* defaultCMCDSFLoatingLeg();
};

//typedef smartConstPtr<CMCDSFloatingLeg>  CMCDSFloatingLegConstSP;
//typedef smartPtr<CMCDSFloatingLeg>       CMCDSFloatingLegSP;


DRLIB_END_NAMESPACE

#endif
