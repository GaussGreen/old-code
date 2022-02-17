/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/callfixedshare.h
// Purpose:     fixed share call for PEPS-like
// Author:      Wang
// Created:     2004/08/17
// RCS-ID:      $Id: callfixedshare.h,v 1.18 2005/03/31 10:24:34 pedro Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/callfixedshare.h
    @brief declaration of the fixed share call class, used for PEPS-like.
    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CALLFIXEDSHARE_H_
#define _ITO33_FINANCE_BONDLIKE_CALLFIXEDSHARE_H_

#include "ito33/date.h"

#include "ito33/finance/bondlike/call.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
   When the issuer calls back the instrument, it will be converted into a 
   fixed number of shares.
 */
class ITO33_DLLDECL CallFixedShare : public Call
{
public:

  /**
     Constructs the call using all the required data

     @param startDate The start date of a call period.
     @param endDate The end date of a call period.
     @param dRatio The applicable conversion ratio during this call period.
   */
  CallFixedShare(Date startDate, Date endDate, double dRatio);

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
     Gets the applicable conversion ratio during the call period

     @return The conversion ratio.
   */
  double GetRatio() const { return m_dRatio; }

  /**
     The call trigger expressed as a percentage of the conversion price

     @param dTriggerRate The trigger (expressed as a percentage of the 
            conversion price)
   */
  void SetTrigger(double dTriggerRate);

  /**
     The call trigger expressed as a percentage of the conversion price

     @return The trigger (expressed as a percentage of the conversion price)
   */
  double GetTrigger() const { return m_dTriggerRate; }

  /**
     Indicates if it is a soft call

     @return true if it is a soft call, false otherwise
   */
  bool IsSoft() const { return m_dTriggerRate > 0.; }

  //@}

  // implement virtual function  
  virtual bool HasMakeWhole() const;

  /**
    @internal

    @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;

private:

  Date m_startDate;
  
  Date m_endDate;
  
  double m_dRatio;
  
  double m_dTriggerRate;

}; // class CallFixedShare


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CALLFIXEDSHARE_H_
