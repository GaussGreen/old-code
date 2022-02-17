/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/onetouches.h
// Purpose:     Class holding multiple one touches having the same maturity
// Created:     2006/02/23
// RCS-ID:      $Id: onetouches.h,v 1.2 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/onetouches.h
    @brief declaration of class holding multiple one touches having the same
           maturity and the same barrier type
 */

#ifndef _ITO33_FINANCE_EXOTICOPTION_ONETOUCHES_H_
#define _ITO33_FINANCE_EXOTICOPTION_ONETOUCHES_H_

#include "ito33/sharedptr.h"
#include "ito33/list.h"
#include "ito33/vector.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL OneTouch;

/**
    OneTouches can be used internally to price at the same time a list of 
    one touches having the same maturity and barrier type.
 */
class OneTouches
{
public:

  typedef std::list< shared_ptr<OneTouch> > Elements;

  /**
      Creates a OneTouches object.

      We may not respect the original ordering of the list since the one 
      touches are internally sorted according to the barrier level.

      @param optionList list of options
   */
  OneTouches(const Elements& oneTouchList);

  // Default dtor is ok

  /**
      @name Methods for accessing the OneTouches.
   */
  //@{

  /**
      Gets the spots corresponding to the list of one touches.

      @return spots(ascending order) corresponding to the list of one touches
   */
  const std::vector<double>& GetSpots() const
  {
    return m_pdSpots;
  }
  
  /** 
      Gets the market prices corresponding to the list of one touches.

      @return market prices corresponding to the list of one touches
   */
  const std::vector<double>& GetMarketPrices() const
  {
    return m_pdMarketPrices;
  }
  
  /**
      Gets the reference one touch.

      @return the reference one touch that will be used as a normal one touch
              pricing
   */
  shared_ptr<OneTouch> GetOneTouch() const
  { 
    return m_oneTouch;
  }

  //@}


protected:

  Elements m_oneTouchList;
  
  shared_ptr<OneTouch> m_oneTouch;

  std::vector<double> m_pdSpots;

  std::vector<double> m_pdMarketPrices;

}; // class OneTouches


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_EXOTICOPTION_ONETOUCHES_H_
