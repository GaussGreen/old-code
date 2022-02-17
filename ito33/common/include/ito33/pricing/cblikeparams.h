/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cblikeparams.h
// Purpose:     params class for CB-like
// Created:     2004/08/19
// RCS-ID:      $Id: cblikeparams.h,v 1.66 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cblikeparams.h
    @brief params class for CB-like

    Implementation of the params class for CB-like.
 */

#ifndef _ITO33_PRICING_CBLIKEPARAMS_H_
#define _ITO33_PRICING_CBLIKEPARAMS_H_

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/list.h"
#include "ito33/date.h"

#include "ito33/pricing/cbeventmanager.h"
#include "ito33/pricing/cbevent.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/params.h"
#include "ito33/pricing/cashflows.h"

#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/bondlike/triggeraspercentageof.h"

namespace ito33 
{
 
namespace pricing
{

  class PathDepStructure;
  class Model;


/**
   Params class for CB-like.

   Pricing parameters for CB-like
 */
class CBLikeParams : public Params
{
public: 

  /**
     Ctor by pricing CB-Like contract. 

     Other members must be initialized by the SetXXX() functions

     @param cbLike reference to CB-like contract objet
   */
  CBLikeParams(CBLike& cbLike) : Params(cbLike), m_cbLike(cbLike),
                                 m_dFXRate(1.), 
                                 m_bIsStartOfYear(false) { }
  
  /**
     Creates Params object by CB-like contract and common financial and 
     numerical datas. 

     @param cbLike reference to CB-like contract objet
     @param sessionData reference to financial SessionData
     @param pNumParams numeric parameters
     @param pMeshParams parameters for mesh builder
   */
  CBLikeParams(CBLike& cbLike,
               const finance::SessionData& sessionData,
               const shared_ptr<numeric::NumParams>& pNumParams,
               const shared_ptr<numeric::MeshParams>& pMeshParams)
             : Params(cbLike, sessionData, pNumParams, pMeshParams),
               m_cbLike(cbLike),
               m_dFXRate(1.), m_bIsStartOfYear(false) { }

  /// virtual dtor for base class
  virtual ~CBLikeParams() { }

  /// @name implmentation of virtual functions
  //@{
  virtual void Init();
  
  // Path dependent pricing needs multiple copies of the pricing
  // objects.  Since this (cblikeparams) base class is passed
  // to cbricer, need a function for the derived classes to 
  // imlemement for copying/cloning the real object.
  virtual CBLikeParams* Clone() const = 0;

  // The default implementation take the list of special times as 
  // two parts: one from contract and one from session. but this 
  // becomes not true with New Share feature of Convertible Bond
  void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const;
  
  /**
     Updates local current time and the eventmanager.

     @param dTime is a double standing for a time
   */
  virtual void Update(double dTime);
  
  /**
     Updates time at maturity for different members (time, eventmanager etc).

     @param dTime should be the time at maturity
   */
  virtual void SetInitialState(double dTime);

  /**
     Updates the index of the mono date event.

     @param pCBEvent the event to  treat
     @param bConstraintsNeedUpdate true if update needed    
   */
  virtual void MonoDateEventIndexTreatment(const CBEvent* pCBEvent, 
                                           bool& bConstraintsNeedUpdate);
  
  //@}
  
  /**
     Gets the array of recovery values at current time.
     
     @param pdS the spots at which the recovery values will be computed.
     @param nNbS array size
     @param pdRecoveryValues array of recovery values
   */
  void 
  GetRecoveryValues(const double* pdS, size_t nNbS, double *pdRecoveryValues);
    
  /**
     Updates the index of the (mono date) cb events at current time.
     (MonoDatePut, MonoDateCall, monoDateConversion, coupon)
    
     @return true if the constraints need update because of some index changes
             false otherwise
   */
  bool UpdateMonoDateEventIndex();

  ///@name member functions for getting some important pointers
  //@{

  /// Gets a reference to the CB-like contract
  CBLike& GetCBLike() const { return m_cbLike; }
  
  /**
     Gets the Cashflows in the CB-like contract.

     @return the Cashflows in the CB-like contract
   */
  CashFlows* GetCashFlows() const { return m_cbLike.GetCashFlows(); }
  
  /**
     Gets the CallProvision in the CB-like contract.

     @return the callProvision in the CB-like contract
   */
  CallProvisions* GetCalls() const 
  { 
    return m_cbLike.GetCalls(); 
  }

  /**
     Gets the ConversionProvision in the CB-like contract.

     @return the ConversionProvisions in the CB-like contract
   */
  ConversionProvisions* GetConversions() const 
  { 
    return m_cbLike.GetConversions(); 
  }

  /**
     Gets the Puts in the CB-like contract.

     @return the Puts in the CB-like contract
   */
  CBPuts* GetPuts() const { return m_cbLike.GetPuts(); }
 
  //@} // member functions for getting some important pointers

  /**
     Gets the financial payoff as initial values for the price.
     
     This is probably a temporay solution. This function can be used to 
     pass around model independant informations, but for passing the whole
     greek structure(price, vega...), it needs to be done at pricer level
     (or other place?)

     @param pdS the spots at which the values will be computed
     @param nNbS the size of the array
     @param pdValues the computed values
   */
  void GetInitialValues(const double* pdS, size_t nNbS, double *pdValues);

  /// @name new share related functions
  //@{

  /** 
     Sets to false the boolean m_bIsStartOfYear
   */
  void DisableStartOfYear()
  {
    m_bIsStartOfYear = false;
  }
  
  /**
     Indicates if we are at a fiscal year start time.
   
     @return true if we are at a fiscal year start time, false otherwise.
   */
  bool IsStartOfYear() const { return m_bIsStartOfYear; }
   
  /**
     Indicates if new share method must be used.
     
     @return true if new share method must be used, false otherwise
   */
  bool HasNewShare() const;
  
  /**
     Gets a special cb for the new share case.

     This special cb is used to compute new share prices at maturity of
     our contract.

     @return the special cb for the new share case
   */
  AutoPtr<CB> GetNewShareCB();

  /**
     Gets the cb params for new share.

     @param pdSpots array of spots
     @param nNbSpots number of spots

     @return a cb like params
   */
  shared_ptr<CBLikeParams> 
  GetNewShareParams(const double *pdSpots, const size_t nNbSpots);

  //@}    //new share related functions

  /// @name event handle related
  //@{

  /**
     Gets the next CB event occuring at the current time. 
  
     @return the next cbevent occuring at the current time (if any), or NULL if 
             no more cbevents occur at the current time.
   */
  const CBEvent* GetCBEventAtCurrentTime() const 
  {
    return m_CBEventManager.GetEvent();
  }
 
  /**
     Gets all the events stored in the CB event manager.
     
     @return the container of all the CB events stored n the cb event manager
   */
  const std::list< shared_ptr<CBEvent> >& GetAllCBEvents() const 
  { 
    return m_CBEventManager.GetAll(); 
  } 
  
  bool HasEventsNow() const 
  { 
    return m_eventManager.HasEventsNow() || m_CBEventManager.HasEventsNow(); 
  }

  //@} // event handle related


  /// @name functions for getting current state
  //@{
  
  /**
     Gets the index of call

     Note that when this function is called, the state "index of call" must
     have been well defined. Even though this becomes a source of code error,
     we don't like always use the functions with "dTime" and "bPlus" arguments:
     they are time consuming. That is why we have a debug boolean
     "m_bStateFunctionShouldBeCalledWithTimeArgument" to help developer to find
     errors at compiling time.

     However "m_bStateFunctionShouldBeCalledWithTimeArgument" can't solve all
     problems related to state get functions. Developer still need to pay
     many attention when using them.

     The boolean is not check in GetIndexCall and GetIndexConversion, as they
     are needed in Conversions::ComputeRoots() which is called only in
     MeshManager::SetupMe(). In this circonstance, computeRoots() is called
     only when the grid period have both conversion and call, thus we have
     to calculate their index in advance.

     @return index of call or INVALIDINDEX
   */
  size_t GetIndexCall() const { return m_nIdxCall; }

  /**
     Gets the index of conversion.

     @sa GetIndexCall

     @return index of conversion or INVALIDINDEX
   */
  size_t GetIndexConversion() const { return m_nIdxConversion; }

  /**
     Gets the index of Put.

     @sa GetIndexCall

     @return index of Put or INVALIDINDEX
   */
  size_t GetIndexPut() const 
  {    
    ASSERT(!m_bStateFunctionShouldBeCalledWithTimeArgument);
    return m_nIdxPut; 
  }

  /**
     Gets the accrued interest value at the current time.

     @sa GetIndexCall

     @return the accrued interest value at current time 
   */
  double GetAccruedInterest() const;

  /**
     Gets the coupon amount at current time.

     @sa GetIndexCall

     @return coupon amount if there is a coupon at current time, 0 otherwise
   */
  double GetCouponAmount() const
  {
    ASSERT(!m_bStateFunctionShouldBeCalledWithTimeArgument);
    return m_dCouponAmount;
  }
  
  /**
     Gets current make-whole premium

     @sa GetIndexCall
 
     @return current make-whole premium
   */
  double GetCurrentMakeWholePositivePremium() const
  {
    return   m_dMakeWholeCurrentPremium > 0
           ? m_dMakeWholeCurrentPremium
           : 0;
  }
  
  /**
     Gets the claim value at current time.

     @sa GetIndexCall

     @return the claim value at current time
   */
  double GetClaim() const
  {
    ASSERT(!m_bStateFunctionShouldBeCalledWithTimeArgument);
    return m_dClaim;
  }  
  

  /**
     Gets the claim for call or puts.

     @param dTime current time
     @param dYield guaranteed yield
     @param bPlus true if right side, false if not.

     @return the claim at the given time dTime.
   */
  double GetClaimFromYield(double dTime, double dYield, bool bPlus=true) const;


  /**
     Gets the claim for call or puts at the current time.

     @param dYield guaranteed yield
     @param bPlus true if right side, false if not.

     @return the claim value at the current time.
   */
  double GetClaimFromYield(double dYield, bool bPlus=true) const;

  /**
     Gets the conversion price at current time.

     @sa GetIndexCall
     
     @return the conversion price at current time
   */
  double GetConversionPrice() const
  {
    ASSERT(!m_bStateFunctionShouldBeCalledWithTimeArgument);
    return m_dConversionPrice;
  }

  /**
     Gets the conversion price for call trigger at current time.
     
     @return the conversion price for call trigger at current time
   */
  double GetConversionPriceForCall() const
  {
    ASSERT(!m_bStateFunctionShouldBeCalledWithTimeArgument);
    return m_dConversionPriceForCall;
  }

  //@}

  /// @name constraints related
  //@{

  /**
     Gets the constraint values imposed by the call (if any).

     @param pdS the spots at which the constraint values will be computed
     @param nNbS the size of the arrays
     @param pdValues the computed constraint values
     @param nIdxStartConversion index of the point issuer call
                                to force conversion                                
     @param pdNewSharePrices the new share prices
     @param bHasNotice if getting values for call notice
   */
  bool GetCallConstraintValues
       (const double* pdS, size_t nNbS, 
        double* pdValues, size_t& nIdxStartConversion,
        const double* pdNewSharePrices, 
        bool bHasNotice);

  /**
     Gets the constraint values imposed by the conversion (if any).

     @param pdS the spots at which the constraint values will be computed
     @param nNbS the size of the arrays
     @param pdNewSharePrices the new share prices
     @param pdValues the computed constraint values
   */
  bool GetConversionConstraintValues
       (const double* pdS, size_t nNbS, 
        const double* pdNewSharePrices,
        double* pdValues);

  /**
     Gets the constraint values imposed by the put (if any).

     @param pdS the spots at which the constraint values will be computed
     @param nNbS the size of the arrays
     @param pdValues the computed constraint values
   */
  bool GetPutConstraintValues
       (const double* pdS, size_t nNbS, double* pdValues);
  
  //@} // constraint related

  /**
     Checks if the call has notice period.

     @return true if the its call has notice period
   */
  bool HasNoticePeriod() const
  {
    return GetCalls()->HasNoticePeriod();
  }

  /**
     Gets the notice period as a fraction of year.

     @return Notice period in days;
   */
  double GetNoticePeriod() const
  {
    return GetCalls()->GetNoticePeriod();
  }

  /**
     This function checks whether we have a non monodate call at a given time.
     If yes, we set the m_nIdxCall to its index and the function returns true.
     otherwise, m_nIdxCall is set to be invalid and the function returns false.
        
     @param dTime the time at which we want to know if we have a non monodate 
                   call
     @param bPlus at which side we are looking for existence of call
     @return true if we have a non monodate call at the given time
   */
  bool HasContinousCallAt(double dTime, bool bPlus)
  {
    m_nIdxCall = GetCalls()->GetIndexContinousCallAt(dTime, bPlus);

    return m_nIdxCall != INVALIDINDEX;
  }

  /**
     This function checks whether we have a non monodate Conversion at
     a given time.
     If yes, we set the m_nIdxConversion to its index and
             the function returns true.
     otherwise, m_nIdxConversion is set to be invalid and
             the function returns false.
        
     @param dTime the time at which we want to know if we have a non monodate
                   Conversion
     @param at which side we are looking for existence of conversion
     @return true if we have a non monodate Conversion at the given time
   */
  bool HasContinousConversionAt(double dTime, bool bPlus)
  {
    m_nIdxConversion = 
      GetConversions()->GetIndexContinousConversionAt(dTime, bPlus);

    return m_nIdxConversion != INVALIDINDEX;
  }

  /**
     Gets the call notice cb params for the current time.
   */
  virtual CBLikeParams* GetCallNoticeParams();
   
  /**
     Sets the foreign exchange (FX) rate at the given value dCurrentFXRate

     @param dCurrentFXRate the current FX rate.
   */
  void SetFXRate(double dCurrentFXRate){ m_dFXRate = dCurrentFXRate; }
   
  /**
     Gets the current foreign exchange (FX) rate.

     In the case of a fixed quanto, the current fx rate is the fixed FX rate.

     @return the current FX rate
   */
  double GetFXRate() const { return m_dFXRate; }
   
  /**
     Gets the foreign exchange (FX) rate at the time dTime.

     In the case of a fixed quanto, the current fx rate is the fixed FX rate.

     @param dTime the time at which we need the FX rate
     @return the current FX rate
   */
  inline double GetFXRate(double dTime) const;

  /// @name path dependent related
  //@{

  /**
     Checks if there is coco clause that requires a path dependent pricing.
   */
  bool HasPathDepCoCo() const;

  /**
     Checks if there is call clause that requires a path dependent pricing.
   */
  bool HasPathDepCall() const;

  /**
     Checks if the pricing of the convertible like is path dependent.
   */
  virtual bool IsPathDependent() const
  {
    return HasPathDepCoCo() || HasPathDepCall(); 
  }
  
  virtual std::vector<double> 
  ConstructPathDepGrid(Model& model) const;

  virtual std::vector< AutoPtr<CBLikeParams> >
  ConstructPaths(const std::vector<double>& grid) const;

  virtual std::list< shared_ptr<PathDepEvent> >
  ConstructPathDepEvents(const std::vector< AutoPtr<CBLikeParams> >& paths) const;

  virtual size_t GetPathToSave(const std::vector<double>& grid) const;

  virtual void InitPaths(PathDepStructure& pathDepStruct);

  //@}

# ifndef NDEBUG
  
  bool m_bStateFunctionShouldBeCalledWithTimeArgument;

# endif


protected:

  ///
  virtual void StateValuesMustBeReset()
  {
    #ifndef NDEBUG
    {
      m_bAccruedInterestUpdated = false;
      m_bClaimUpdated = false;
      m_bConversionPriceUpdated = false;
      m_bConversionPriceForCallUpdated = false;
      m_bFXRateUpdated = false;
    }
    #endif
  }

  #ifndef NDEBUG
    bool
      m_bAccruedInterestUpdated,
      m_bClaimUpdated,
      m_bConversionPriceUpdated,
      m_bConversionPriceForCallUpdated,
      m_bFXRateUpdated;
  #endif

protected:
 
  /**
     Updates the fx rate at current time.
   */
  void UpdateFXRate()
  {
    m_dFXRate = GetFXRate(m_dCurrentTime);

    #ifndef NDEBUG
      m_bFXRateUpdated = true;
    #endif
  }

  /**
     Updates Claim at given side of current time 

     @param bPlus true if right side, false if not.
   */
  void UpdateClaim(bool bPlus = true);

  /**
     Updates the conversion price at current time.

     Must be called after the function UpdateClaim().
     Must be called after the function UpdateAccruedInterest().
   */
  void UpdateConversionPrice()
  {
    ASSERT_MSG
      (
        m_bAccruedInterestUpdated,
        "Accrued interest is not initialized."
      );
    ASSERT_MSG
      (
        m_bClaimUpdated,
        "claim is not initialized."
      );
    ASSERT_MSG
      (
        !m_cbLike.IsCrossCurrency() || m_bFXRateUpdated,
        "FxRate is not initialized."
      );

    if ( GetConversions()->GetNbConversions() != 0 )
    {
      finance::TriggerAsPercentageOf 
        of = GetConversions()->GetTriggerAsPercentageOf();

      m_dConversionPrice = GetConversions()->GetConversionPrice(of);
      
      if ( GetCalls()->GetNbCalls() != 0  )
      {
        finance::TriggerAsPercentageOf 
          of = GetCalls()->GetTriggerAsPercentageOf();
         
        m_dConversionPriceForCall = GetConversions()->GetConversionPrice(of);
      }
    }

    #ifndef NDEBUG
      m_bConversionPriceUpdated = true;
      m_bConversionPriceForCallUpdated       = true;
    #endif
  }

  /**
     Updates the accrued interest at current time.
   */
  void UpdateAccruedInterest()
  {
    if ( m_dSlopeCoupon != 0 )
      m_dAccruedInterest = m_dSlopeCoupon  
          * (m_dCurrentTime - GetCashFlows()->GetTime(m_nIdxNextCoupon - 1));
    else
      m_dAccruedInterest = 0;

    #ifndef NDEBUG
      m_bAccruedInterestUpdated = true;
    #endif
  }

  /**
     Sets up the initial call notice params. The objet will be reused by 
     GetCallNoticePeriod instead of returning each time a new one.
   */
  void SetupCallNoticeParams();

  /**
     Updates purely discrete Cb Events related data(coupon amount, put index).

     Normally, the data only makes sense for DoEvents in InstData. But we still
     need to set them to invalid/non signifcant values outside of DoEvents 
     with the actual code. 
   */
  void UpdateDiscreteCBEventData();

  CBEventManager m_CBEventManager;

  /// helper variable indicating whether we are at end of grid
  double m_dFormerTime;

  /**
     The forward discount factor of yield curve associated to the currency of
     instrument.  

     We may also need this paramater for other instruments. But for now, 
     it is calculated thanks to m_dFormerTime, so stays (temporarily?) in
     CBLikeParams.
   */
  double m_dDerivativeYieldCuveFDF;

  /// @name identities at current time
  //@{

  bool m_bPlus;
  double m_dClaim;
  double m_dConversionPrice;
  double m_dConversionPriceForCall;
  double m_dAccruedInterest;
  double m_dCouponAmount;
  double m_dSlopeCoupon;
  double m_dMakeWholeCurrentPremium;
  double m_dFXRate;

  size_t m_nIdxCall;
  size_t m_nIdxConversion;
  size_t m_nIdxPut;
  size_t m_nIdxNextCoupon;

  //@}

  /// @name Make-Whole concerned
  //@{
  
  /// initialization of make-whole concerned values if possible
  void InitMakeWholeData();

  /// artificial make-whole final premium.
  double m_dMakeWholeFinalPremium;

  /// Index of coupon whose date is after the soft call
  size_t m_nIndexAfterLastCouponMakeWhole;
  //@}

  // Boolean use for new share case: true if start of year, false otherwise
  bool m_bIsStartOfYear;

  /// CBLikeParams for solving the small cb problem (notice)
  AutoPtr<CBLikeParams> m_pCallNoticeParams;

  CBLike& m_cbLike;


private:

  NO_COPY_CLASS(CBLikeParams);

}; // class CBLikeParams

inline double CBLikeParams::GetFXRate(double dTime) const
{
  if ( GetCBLike().IsFixedQuanto() )
    return GetConversions()->GetFixedFXRate();

  double dFXRate;   
  
  double pdTimes[2];
  double pdRates[2];
  double pdDerivativeRates[2];

  pdTimes[0] = m_dValuationTime;
  pdTimes[1] = dTime;

  GetYieldCurve()->GetCompoundFactor(pdTimes, pdRates, 2);   
  
  m_cbLike.GetDerivativeCurve()->GetCompoundFactor
                                 (pdTimes, pdDerivativeRates, 2);
  
  dFXRate = GetCBLike().GetSpotFXRate() * pdRates[0] / pdDerivativeRates[0];
  
  dFXRate *= pdDerivativeRates[1] / pdRates[1];

  return dFXRate;
}

inline double CBLikeParams::GetAccruedInterest() const
{ 
  ASSERT(!m_bStateFunctionShouldBeCalledWithTimeArgument);

  return m_dAccruedInterest; 
}

inline bool CBLikeParams::HasNewShare() const
{  
  ConversionProvisions& conv = *(GetConversions());
  size_t nNbConversions = conv.GetNbConversions();
  if ( nNbConversions == 0 )
    return false;
  else
    return m_cbLike.HasNewShare();
}


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBLIKEPARAMS_H_
