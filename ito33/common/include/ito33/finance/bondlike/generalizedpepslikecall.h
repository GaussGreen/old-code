 /////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/generalizedpepslike.h
// Purpose:     PEPS-like financial class
// Author:      Wang
// Created:     2004/08/17
// RCS-ID:      $Id: generalizedpepslikecall.h,v 1.4 2006/03/23 09:27:14 yann Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/generalizedpepslikecall.h
    @brief declaration of the financial PEPS-like class   
 */

#ifndef _ITO33_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_H_
#define _ITO33_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_H_

#include "ito33/date.h"
#include "ito33/dlldecl.h"

#include "ito33/finance/bondlike/call.h"
#include "ito33/finance/bondlike/generalizedpepslikecalltype.h"

namespace ito33
{

namespace XML
{
  class Tag;
}
namespace finance
{

/**
   This class describes the conversion at issuer's option
   for a GeneralizedPEPSLike insturment.
 */
class ITO33_DLLDECL GeneralizedPEPSLikeCall : public Call
{
public:
  /**
     Constructs the call using all the required data

     @param startDate The start date of the call.
     @param endDate The end date of the call.
     @param dTriggerRate the trigger expressed as percentage of higher strike
            which should be greater than or equal to 100%.
     @param type variable share or fixed share?
   */
  GeneralizedPEPSLikeCall(Date startDate, Date endDate,
            double dTriggerRate, GeneralizedPEPSLikeCallType type);

  // copy constructor is ok

  /// @name Getter functions
  // @{

  /**
     Gets the call period start date.

     @return The call period start date.
   */
  Date GetStartDate() const { return m_startDate; }

  /**
     Gets the call period end date.

     @return The call period end date.
   */
  Date GetEndDate() const { return m_endDate; }

  /**
     The trigger expressed as a percentage of the higher strike of the peps.

     @return The trigger
   */
  double GetTrigger() const { return m_dTriggerRate; }

  /**
     The type (Fixed share or variable share).

     @return the type of the call
   */
  GeneralizedPEPSLikeCallType GetType() const
  {
    return m_type;
  }

  //@}

  // implement virtual function  
  virtual bool HasMakeWhole() const;

  /**
     XML Dump
     @internal

     @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;

private:

  Date m_startDate;
  
  Date m_endDate;
  
  double m_dTriggerRate;

  GeneralizedPEPSLikeCallType m_type;
};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_H_
