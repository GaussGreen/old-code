/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/impliedcdsspreads.h
// Purpose:     Class for computing implied CDS spreads
// Author:      ITO33
// Created:     2004/11/19
// RCS-ID:      $Id: impliedcdsspreads.h,v 1.9 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/impliedcdsspreads.h
    @brief Class for computing implied CDS spreads
*/

#ifndef _ITO33_IHG_IMPLIEDCDSSPREADS_H_
#define _ITO33_IHG_IMPLIEDCDSSPREADS_H_

#include "ito33/vector.h"
#include "ito33/useexception.h"

#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/finance/frequency.h"

#include "ito33/ihg/common.h"

extern const ito33::Error ITO33_BAD_DATA;

namespace ito33 
{
  
namespace finance
{
  class ITO33_DLLDECL SessionData;
}

namespace ihg 
{

  class ITO33_IHG_DLLDECL TheoreticalModel;

/**
   Class for computing implied CDS spreads for a term structure
   of maturities.
 */
class ImpliedCDSSpreads
{

public:

  /**
     Ctor initializes the spread calculator with all the spread stream
     details except the spread values, which will be computed to give prices
     of zero.

     @param contractingDate contracting date
     @param firstPaymentDate first payment date
     @param dcc day count convention
     @param freq payment frequency (monthly, semi-annual, etc)
     @param dRecoveryRate recovery rate for the CDS contracts
   */
  ImpliedCDSSpreads(Date contractingDate,
                    Date firstPaymentDate,
                    Date::DayCountConvention dcc,
                    finance::Frequency freq,
                    double dRecoveryRate)                         
           : m_contractingDate(contractingDate),
             m_firstPaymentDate(firstPaymentDate),
             m_dcc(dcc),
             m_freq(freq),
             m_dRecoveryRate(dRecoveryRate)
  { 
    // Verify the data
    if ( dRecoveryRate < 0. || dRecoveryRate > 1. )
      throw EXCEPTION_MSG
            (
              ITO33_BAD_DATA,
              TRANS("Implied CDS spread calculator: Setting invalid recovery rate.")
            );

    if ( m_firstPaymentDate <= contractingDate )
	    throw EXCEPTION_MSG
            ( 
              ITO33_BAD_DATA,
              TRANS("Implied CDS spread calculator: First payment date equal or before issue date.")
            );

    if ( !m_freq ) 
      throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Implied CDS spread calculator: Invalid frequency for payment.")
          );
  }

  // Default dtor is ok

  /**
      Compute the implied spreads.

      @param model the theroretical used or pricing
      @param sessionData the session data used for pricing
      @param pMaturityDates the maturity dates of the CDS term structure
      @return vector of spread values corresponding to the maturity dates
  */
  std::vector<double> Compute( const ihg::TheoreticalModel& model, 
                               const shared_ptr<finance::SessionData>& sessionData, 
                               const std::vector<Date>& pMaturityDates);


protected:

  /// contracting date for spread stream
  Date m_contractingDate;

  /// first payment date for spread stream
  Date m_firstPaymentDate;

  /// day count convention for spread stream
  Date::DayCountConvention m_dcc;

  /// payment frequency for spread stream
  finance::Frequency m_freq;

  /// The recovery rate
  double m_dRecoveryRate;

  NO_COPY_CLASS(ImpliedCDSSpreads);

}; // class ImpliedCDSSpreads


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_IMPLIEDCDSSPREADS_H_
