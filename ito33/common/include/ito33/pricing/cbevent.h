/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cbevent.h
// Purpose:     convertible bond event class
// Author:      Nabil
// Created:     2003/11/20
// RCS-ID:      $Id: cbevent.h,v 1.24 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_CBEVENT_H_
#define _ITO33_PRICING_CBEVENT_H_

#include "ito33/beforestd.h"
#include <functional>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/numeric/predicatetime.h"

namespace ito33
{

namespace pricing
{

enum CBEventType 
{
  CBET_Void,
  // For cb option contract
  // these events must declared before all the other and in this order
  // to guarantee the good work of the update for the data of the cb option
  // at the CBET_CBOptionMaturity event.
  CBET_CBOptionMaturity,
  CBET_CBOptionFloatingCoupon,

  CBET_EndConversion,
  CBET_EndCall,
  CBET_PayDividend,
  CBET_StartOfYear,
  CBET_PayCoupon,
  CBET_MonoDatePut,
  CBET_MonoDateCall,
  CBET_MonoDateConversion,
  CBET_Maturity,
  CBET_StartCall,
  CBET_StartConversion,
  CBET_ChangeOfRoot,
  CBET_StartAccretion,

  CBET_Max
};


class CBEvent
{
public:

  size_t
    m_nIndex;
  
  double
    m_dTime;
  
  CBEventType
    m_cbeventType;

  CBEvent() {}
  
  CBEvent(double dTime, CBEventType cbeventType) 
  { 
    SetEvent(dTime, cbeventType); 
  }
  
  CBEvent(double dTime, CBEventType cbeventType, size_t nIndex) 
  {
    SetEvent(dTime, cbeventType, nIndex);
  }

  ~CBEvent() {}

  void SetEvent(double dTime, CBEventType cbeventType) 
  {
    m_dTime = dTime; 
    m_cbeventType = cbeventType;
  }

  void SetEvent(double dTime, CBEventType cbeventType, size_t nIndex) 
  {
    m_dTime = dTime; 
    m_cbeventType = cbeventType; 
    m_nIndex = nIndex;
  }

  double GetTime() const { return m_dTime; }
    
  /**
      less comparison operator.
    
      @param cbevent1 the first event
      @param cbevent2 the second event to compare
   */
  friend bool operator < (const CBEvent& cbevent1, const CBEvent& cbevent2)
  {
    // If events occur at the same time, check event type for the order
    if ( numeric::AreTimesEqual(cbevent1.m_dTime, cbevent2.m_dTime) )
      return  cbevent1.m_cbeventType < cbevent2.m_cbeventType;
    else
      return cbevent1.m_dTime < cbevent2.m_dTime;
  }

  /**
      great comparison operator.
    
      @param cbevent1 the first event
      @param cbevent2 the second event to compare
   */
  friend bool operator > (const CBEvent& cbevent1, const CBEvent& cbevent2)
  {
    // If events occur at the same time, check event type for the order
    if ( numeric::AreTimesEqual(cbevent1.m_dTime, cbevent2.m_dTime) )
      return  cbevent1.m_cbeventType > cbevent2.m_cbeventType;
    else
      return cbevent1.m_dTime > cbevent2.m_dTime;
  }

};

} // namespace pricing

} // namespace ito33

namespace std
{

/**
    Used for sorting a list of pointers on events of type CBEvent.
    Pay attention to the name greater, it is used for ascending sort here
 */
template<>
inline bool std::greater< ito33::shared_ptr<ito33::pricing::CBEvent> >::operator()
     (const ito33::shared_ptr<ito33::pricing::CBEvent>& pcbEvent1, 
      const ito33::shared_ptr<ito33::pricing::CBEvent>& pcbEvent2) const
{
  return *pcbEvent1 < *pcbEvent2;
}

} // namespace std

#endif // #ifndef _ITO33_PRICING_CBEVENT_H_
