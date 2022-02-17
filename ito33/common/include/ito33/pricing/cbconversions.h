/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cbconversions.h
// Purpose:     convertible bond conversions class
// Author:      Laurence
// Created:     2004/03/12
// RCS-ID:      $Id: cbconversions.h,v 1.42 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_CBCONVERSIONS_H_
#define _ITO33_PRICING_CBCONVERSIONS_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/numeric/mesh/roots.h"

#include "ito33/pricing/conversionprovisions.h"
#include "ito33/finance/bondlike/cocotype.h"
#include "ito33/finance/bondlike/conversionschedule.h"

namespace ito33
{
  class Date;

namespace finance
{
  class ITO33_DLLDECL ConversionPeriod;
  class ITO33_DLLDECL ConversionSchedule;
}

namespace pricing
{

/// Declaration of the CBConversions class
class CBConversions : public ConversionProvisions
{
public:

  CBConversions(): ConversionProvisions() { };

  CBConversions(const shared_ptr<finance::ConversionSchedule>& pConversions,
                Date valuationDate);

  // default dtor is ok

  virtual ~CBConversions() {}

  virtual void 
  ComputeRoots(double dTime, size_t& nNbRoots, numeric::mesh::Root* pRoots,
               bool bPlus = true);

  /// @name implement virtual functions
  //@{

  virtual bool GetGrossParities(const double* pdS, size_t nNbS, 
                        const double* pdNewSharePrices,
                        double* pdValues) const;
  
  virtual double 
  GetConversionPrice
  (finance::TriggerAsPercentageOf triggerAsPercentageOf) const;

  virtual double 
  GetConversionPrice
  (double dTime, finance::TriggerAsPercentageOf triggerAsPercentageOf, 
   bool bPlus = true) const;

  virtual double GetRatio(size_t nIdx) const
  {
    return m_pdRatios[nIdx];
  }
  
  void SetRatios(double dRatio)
  {
    size_t nNbRatios = m_pdRatios.size();
    for (size_t nIdx = 0; nIdx < nNbRatios; nIdx++)
      m_pdRatios[nIdx] = dRatio;
  }
   
  /**
     Trigger rate is always a function of the base conversion ratio.
     When SetRatio is called the trigger rates need to
     still be based on the initial base conversion ratio.

     @param dRatio current conversion ratio
   */
  void ChangeTriggerRates(double dRatio)
  {
    size_t nNbTriggerRate = m_pdTriggerRates.size();
    
    for (size_t nIdx = 0; nIdx < nNbTriggerRate; nIdx++)
      m_pdTriggerRates[nIdx] *= dRatio / (m_pdRatios[nIdx] + 1.e-100);
  }

  //@}


protected:
  
  /**
    Check if the current window must be split.  For example, the
    window is longer than one quarter, and the trigger rate changes.

    Pass in data specifically (instead of contained in a conversion period)
    so that derived classes can also use the function.  For example,
    sharedependentconverion uses this function, but does not use
    conversion periods.

    @param startDate window start date
    @param endDate window end date
    @param cocoType the type of CoCo window
    @param dTriggerRate the trigger rate
    @param dExtremeTriggerRate the extreme trigger rate
    @param dRateChange the (yearly) change in the trigger rate
    @param bIsLastTriggerMet was the last trigger condition met
  */
  bool IsWindowTooLong(Date startDate,
                       Date endDate,
                       finance::CoCoType cocoType,
                       double dTriggerRate,
                       double dExtremeRate,
                       double dRateChange,
                       bool bIsLastTriggerMet);

  /**
    Break up the current window into quarterly chunks.  Adjust the trigger
    rates if needed.

    Pass in data specifically (instead of contained in a conversion period)
    so that derived classes can also use the function.  For example,
    sharedependentconverion uses this function, but does not use
    conversion periods.

    @param windoStartDate window start date
    @param finalEndDate window end date
    @param dRatio the conversion ratio 
    @param dCash the cash paid upon conversion
    @param cocoType the type of CoCo window
    @param dTriggerRate the trigger rate
    @param dExtremeTriggerRate the extreme trigger rate
    @param dRateChange the (yearly) change in the trigger rate
    @param bIsLastTriggerMet was the last trigger condition met
    @param valuationDate the valuaiton date
  */
  void BreakupWindow(Date windowStartDate,
                     Date finalEndDate,
                     double dRatio,
                     double dCash,
                     finance::CoCoType cocoType,
                     double dTriggerRate,
                     double dExtremeRate,
                     double dRateChange,                     
                     bool bIsLastTriggerMet,
                     Date valuationDate);



  /// The ratios
  std::vector<double> m_pdRatios;  

}; // class CBConversions


} // namespace pricing

} // namespace ito33

#endif  //  #ifndef _ITO33_PRICING_CBCONVERSIONS_H_
