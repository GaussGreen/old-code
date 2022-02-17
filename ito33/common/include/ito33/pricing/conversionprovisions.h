/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/conversionprovisions.h
// Purpose:     base class for conversion
// Author:      Laurence
// Created:     2004/03/12
// RCS-ID:      $Id: conversionprovisions.h,v 1.54 2006/05/03 10:14:25 nabil Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_CONVERSIONPROVISIONS_H_
#define _ITO33_PRICING_CONVERSIONPROVISIONS_H_

#include "ito33/vector.h"
#include "ito33/list.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/triggeraspercentageof.h"
#include "ito33/finance/bondlike/triggerincurrencyof.h"
#include "ito33/finance/bondlike/cocotype.h"


namespace ito33
{

  class ITO33_DLLDECL Date;

namespace numeric
{
  namespace mesh
  {
    class Root;
  }
}

namespace pricing
{

class CBLikeParams;

class ConversionProvisions
{ 
public:

  /// default ctor
  ConversionProvisions()
    : m_bKeepAccrued(true), m_bForfeitCoupon(false), 
      m_nNbConversions(0),
      m_triggerInCurrencyOf(finance::TriggerInCurrencyOf_Max),
      m_triggerAsPercentageOf
          (finance::TriggerAsPercentageOf_Principal),
      m_dFixedFXRate(0.)
  {
  }
  
  ConversionProvisions(bool bKeepAccrued, bool bForfeitCoupon)
    : m_bKeepAccrued(bKeepAccrued), m_bForfeitCoupon(bForfeitCoupon), 
      m_nNbConversions(0),
      m_triggerInCurrencyOf(finance::TriggerInCurrencyOf_Max),
      m_triggerAsPercentageOf
          (finance::TriggerAsPercentageOf_Principal),
      m_dFixedFXRate(0.)
  {
  }

  /// virtual dtor
  virtual ~ConversionProvisions() { }

  /**
     Sets the trigger as percentage of.
    
     @param the trigger as percentage of
   */
  void SetTriggerAsPercentageOf(finance::TriggerAsPercentageOf of)
  {
    m_triggerAsPercentageOf = of; 
  }

  /**
     Sets the trigger in currency of.
    
     @param the trigger as percentage of
   */
  void SetTriggerInCurrencyOf(finance::TriggerInCurrencyOf of)
  {
    m_triggerInCurrencyOf = of; 
  }

  /**
     Sets the fixed FX rate (used in cross currency case).
    
     @param the fixed FX rate of the security
   */
  void SetFixedFXRate(double dFixedFXRate) { m_dFixedFXRate = dFixedFXRate; }

  /**
     Gets the fixed FX rate (in cross currency case).
    
     @return the fixed FX rate of the security
   */
  double GetFixedFXRate() { return m_dFixedFXRate; }

  /// @name member functions for get conversion values in different cases
  //@{
  
  /**
     Evaluates the gross parity(r * S). Accured is not included for now.
     
     The resulted values are expressed in the currency of security for a cross
     currency since they are used as constraints so to be compared against 
     prices.

     @param pdS given spot array 
     @param nNbS @sa GetConversionConstraintValues
     @param pdNewSharePrices new share prices 
     @param pdValues gross parity for given spots
     @return true when we have conversion clause here. otherwise false,
             and does nothing for pdValues
   */
  virtual bool GetGrossParities
               (const double* pdS, size_t nNbS,
                const double* pdNewSharePrices,
                double* pdValues) const = 0;
  
  /**
     Evaluates conversion constraint values, if any, for given spots.
 
     @TODO : should the following functions return bool or void? That is,
             should we check index inside or outside the function

     @param pdS given spot array
     @param nNbS size of the given spots array
     @param pdNewSharePrices new share prices
     @param pdValues (output) conversion constraint values
     @return true when we have conversion clause here. otherwise false,
             and does nothing for pdValues
   */
  bool GetConversionConstraintValues
       ( const double* pdS, size_t nNbS, 
         const double* pdNewSharePrices, 
         double* pdValues ) const;

  /**
     Evaluates the conversion values forced by early redemption (call).
     therefore make whole is taken in account and trigger is ignored.

     @param pdS given spot array
     @param nNbS size of the given spots array
     @param pdNewSharePrices new share prices 
     @param pdValues the forced conversion values(trigger needs not be met)
     @return true when we have conversion clause here. otherwise false,
             and does nothing for pdValues
   */
  bool GetForcedConversionValues
       (const double* pdS, size_t nNbS, 
        const double* pdNewSharePrices, 
        double* pdValues) const;

  //@}

  virtual void ComputeRoots(double dTime, 
    size_t &nNbRoots, numeric::mesh::Root *pRoots, bool bPlus = true) = 0;

  /*
     The conversion price at given time, is always expressed in the currency of 
     stock even if the security is a cross currency since the conversion price 
     is used to work with the triggers.
  */
  virtual double 
  GetConversionPrice
  (finance::TriggerAsPercentageOf triggerAsPercentageOf ) const = 0;

  virtual double 
  GetConversionPrice
  (double dTime, finance::TriggerAsPercentageOf triggerAsPercentageOf, 
   bool bPlus = true) const = 0;
  
  bool IsKeepAccrued() const { return m_bKeepAccrued; }

  bool IsForfeitCoupon() const { return m_bForfeitCoupon; }
  
  size_t GetNbConversions() const { return m_nNbConversions; }

  double GetStartTime(size_t nIdx) const
  {
    return m_pdStartTimes[nIdx];
  }
  
  double GetEndTime(size_t nIdx) const
  {
    return m_pdEndTimes[nIdx];
  }  

  /**
     All conversion base ratios are resetted
     to dRatio.

     @param dRatio base conversion ratio for all conversion period
   */
  virtual void SetRatios(double /* dRatio */) = 0;


  /**
     Helper function to check if the conversion is continous, so not monodate.
   */
  bool IsContinuous(size_t nIdx) const 
  {
    return m_pdStartTimes[nIdx] != m_pdEndTimes[nIdx];
  }

  /**
     Sets pointer to related params object

     @param pParams pointer to params
   */
  void SetParams(const CBLikeParams* pParams)
  {
    m_pParams = pParams;
  }

  /**
     Gets the Index of continous conversion at given time and at the given side,
     returns INVALIDINDEX if there is no continuous conversion.

     @param dTime given time
     @param bPlus at which side we are looking for index
     @return index of continous Conversion at time dTime
   */
  size_t GetIndexContinousConversionAt(double dTime, bool bPlus);
  
  double GetTriggerLevel() const;

  /*
     Some virtual functions for path dependent, note that if the trigger 
     information is stored in the base class, it's possible to implement 
     these functions in the base, so not virtual.
  */

  /// @name pathdependent functions
  //@{

  void DeactivateDiscretePeriods();
  
  finance::CoCoType GetPeriodType(size_t nIdx) const { return m_pCoCoTypes[nIdx]; }

  double GetTriggerRate(size_t nIdx) const { return m_pdTriggerRates[nIdx]; }
    
  void RemoveTriggers();

  /// Is there an active period with path dependent CoCo feature
  bool HasPathDepCoCo(double dValuationTime, double dMaturityTime);

  finance::TriggerAsPercentageOf GetTriggerAsPercentageOf() const
  {
    return m_triggerAsPercentageOf;
  }

  //@}

protected:
  
  /// keep accrued flag
  bool m_bKeepAccrued;

  /// forfeit coupon
  bool m_bForfeitCoupon;

  /// number of conversion windows    
  size_t m_nNbConversions;
  
  /// start time array
  std::vector<double> m_pdStartTimes;
  
  /// end time array
  std::vector<double> m_pdEndTimes;

  /// Trigger is a percentage of fixed or accreted conversion price
  finance::TriggerAsPercentageOf m_triggerAsPercentageOf;

  /// Trigger is in currency of underlying or security
  finance::TriggerInCurrencyOf m_triggerInCurrencyOf;

  /// The fixed fx rate (for trigger of a cross currency or fixed quanto)
  double m_dFixedFXRate;

  /// The conversion trigger rates (to be multiplied by the base amount)
  std::vector<double> m_pdTriggerRates;

  /// The type of window (quarterly, anytime, etc)
  std::vector<finance::CoCoType> m_pCoCoTypes;

  /// The cash amount paid if conversion occurs
  std::vector<double> m_pdCashs;

  // Discrete periods are not used, except when forced by a call
  /// Flag indicating if the current period is active 
  std::vector<bool> m_pbIsActive;

  // we save the pointer here, otherwise almost all
  // member function should take params as argument
  /// pointer to params 
  const CBLikeParams* m_pParams;

}; // class ConversionProvisions


} // namespace pricing

} // namespace ito33

#endif  //  #ifndef _ITO33_PRICING_CONVERSIONPROVISIONS_H_
