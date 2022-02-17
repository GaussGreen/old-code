/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/onetouches.cpp
// Purpose:     Implementation of ForwardOption class
// Created:     2006/02/23
// RCS-ID:      $Id: onetouches.cpp,v 1.3 2006/08/19 22:40:30 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/onetouches.h"

// operator helps to sort the one touches in descending order of barrier
template<>
inline bool 
std::greater< ito33::shared_ptr<ito33::finance::OneTouch> >::operator()
     (const ito33::shared_ptr<ito33::finance::OneTouch>& pOneTouch1, 
      const ito33::shared_ptr<ito33::finance::OneTouch>& pOneTouch2) const
{
  return pOneTouch1->GetBarrier() > pOneTouch2->GetBarrier();
}

namespace ito33
{

namespace finance
{

OneTouches::OneTouches(const Elements& oneTouchList)
                     : m_oneTouchList(oneTouchList)
{
  ASSERT( !m_oneTouchList.empty() );
  BarrierType barrierType = m_oneTouchList.front()->GetBarrierType();
  
  m_oneTouchList.sort( std::greater< shared_ptr<OneTouch> >() );

  // Choose the one touch so that prices at transformed spots won't need
  // extrapolation
  if ( barrierType == Barrier_UpAndOut )
    m_oneTouch = m_oneTouchList.back();
  else if ( barrierType == Barrier_DownAndOut )
    m_oneTouch = m_oneTouchList.back();

  Date maturityDate = m_oneTouch->GetMaturityDate();
  double dS0 = m_oneTouch->GetSessionData()->GetSpotSharePrice();
  double dBarrier = m_oneTouch->GetBarrier();

  Elements::const_iterator iter = m_oneTouchList.begin();
  for (; iter != m_oneTouchList.end(); ++iter)
  {
    ASSERT( (*iter)->GetMaturityDate() == maturityDate );

    ASSERT( (*iter)->GetBarrierType() == barrierType );  
    
    m_pdSpots.push_back( dS0 * dBarrier / (*iter)->GetBarrier() );

    if ( (*iter)->HasMarketPrice() )
      m_pdMarketPrices.push_back( (*iter)->GetMarketPrice() );
  }  
}

} // namespace finance

} // namespace ito33
