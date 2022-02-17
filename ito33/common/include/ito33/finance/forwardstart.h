/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/forwardstart.h
// Purpose:     financial forward start class
// Author:      wang
// Created:     18/12/2003
// RCS-ID:      $Id: forwardstart.h,v 1.3 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/forwardstart.h
    @brief declaration of the financial forward start class    
 */

#ifndef _ITO33_FINANCE_FORWARDSTART_H_
#define _ITO33_FINANCE_FORWARDSTART_H_

#include "ito33/useexception.h"
#include "ito33/date.h"

#include "ito33/finance/derivative.h"

namespace ito33
{

namespace finance
{


/**
    ForwardStart represents the financial aspects of a forward start instrument.

    It doesn't derive from option because a forward start instrument has a percent
    strike instead of a strike for an option.
 */
class ForwardStart : public Derivative
{
public:
  
  /**
     Create a empty forward start objet
     and initialize the members to invalid values

     Use the setting methods below to initialize the forward start
   */
  ForwardStart() : Derivative(),
                   m_dPercentStrike(-1.)
  { }

  /// Default dtor is ok since no class will derive from ForwardStart ATM

  /**
    Set the percent strike of the forward start

    @param dPercentStrike the percent strike of the option

    @todo Does the percent strike need only to be positif or we should give
          a upper bound for it.
   */
  void SetPercentStrike(double dPercentStrike)
  {
    if (dPercentStrike < 0.)
      throw EXCEPTION_MSG
           (
            ITO33_BAD_DATA,
            TRANS("ForwardStart definition: Setting negatif percent strike")
           );

    if (m_dPercentStrike != dPercentStrike)
    {
      m_dPercentStrike = dPercentStrike; 

      Invalidate();
    }
  }

  /**
     Set the start date of the forward start

     @param StartDate the start date of the forward start
   */
  void SetStartDate(Date StartDate)
  {
    if ( !StartDate.IsValid() )
      throw EXCEPTION_MSG
           (
            ITO33_BAD_DATA,
            TRANS("ForwardStart definition: Setting invalid start date")
           );

    if (m_StartDate != StartDate)
    {
      StartDate = StartDate;

      Invalidate();
    }
  }

  /**
     Get the percent strike of the forward start

     @return The percent strike of the forward start
   */
  double GetPercentStrike() const 
  { 
    if (m_dPercentStrike < 0.)
      throw EXCEPTION_MSG
           (
            ITO33_UNDEF_DATA,
            TRANS("ForwardStart definition: The percent strike is not set.")
           );

    return m_dPercentStrike; 
  }

  /**
     Get the start date of the forward start

     @return the start date of the forward start
   */
  Date GetStartDate() const 
  {
    if ( !m_StartDate.IsValid() )
      throw EXCEPTION_MSG
           (
            ITO33_UNDEF_DATA,
            TRANS("ForwardStart definition: The start date is not set.")
           );

    return m_StartDate; 
  }


protected:
  
  virtual void DoValidate() const
  {
    // forward start is special in the sense that the pricer doesn't care 
    // about the issue date
    if ( !m_MaturityDate.IsValid() )
      throw EXCEPTION_MSG
           (
            ITO33_UNDEF_DATA,
            TRANS("Derivative definition: Setting invalid maturity.")
           );

    // If the Issue date is set, it needs to be smaller than the maturity date.
    if ( m_IssueDate.IsValid() && !(m_IssueDate < m_MaturityDate) )
      throw EXCEPTION_MSG
           (
            ITO33_BAD_DATA,
            TRANS("Derivative definition: The issue date is set but not"
                  "before the maturity.")
           );    

    // forward start specific data validation
    if (m_dPercentStrike < 0.)
      throw EXCEPTION_MSG
           (
            ITO33_UNDEF_DATA,
            TRANS("ForwardStart definition: The percent strike is not set.")
           );
 
    if ( !m_StartDate.IsValid() )
      throw EXCEPTION_MSG
           (
            ITO33_UNDEF_DATA,
            TRANS("ForwardStart definition: The start date is not set.")
           );

    // Start date needs to be smaller than the maturity
    if ( !(m_StartDate < m_MaturityDate) )
      throw EXCEPTION_MSG
           (
            ITO33_BAD_DATA,
            TRANS("ForwardStart definition: The start date is not"
                  "before the maturity date.")
           );    
  }

  double m_dPercentStrike;

  Date m_StartDate;

}; // class ForwardStart

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_FORWARDSTART_H_
