/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/impliedcdsspreads.h
// Purpose:     Class for computing implied CDS spreads
// Created:     2005/02/08
// RCS-ID:      $Id: impliedcdsspreads.h,v 1.7 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/impliedcdsspreads.h
    @brief Class for computing implied CDS spreads
 */

#ifndef _ITO33_FINANCE_IMPLIEDCDSSPREADS_H_
#define _ITO33_FINANCE_IMPLIEDCDSSPREADS_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/finance/frequency.h"

#ifdef __CPP2ANY__
#include "ito33/finance/theoreticalmodel.h"
#endif

#include "ito33/dlldecl.h"

namespace ito33 
{
  
namespace finance
{
  class ITO33_DLLDECL SessionData;
  class ITO33_DLLDECL TheoreticalModel;


/**
   Class for computing implied CDS spreads for a term structure
   of maturities.
 */
class ITO33_DLLDECL ImpliedCDSSpreads
{

public:

  /**
     Ctor initializes the spread calculator with all the spread stream
     details except the spreads, which will be computed to give prices
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
                    double dRecoveryRate) ;

  // Default dtor is ok

  /**
      Computes the implied CDS spreads.

      @param model the theroretical used or pricing
      @param sessionData the session data used for pricing
      @param pMaturityDates the maturity dates of the CDS term structure
      @return vector of spread values corresponding to the maturity dates
   */
  std::vector<double> 
  Compute( const shared_ptr<TheoreticalModel>& model, 
           const shared_ptr<SessionData>& sessionData, 
           const std::vector<Date>& pMaturityDates );


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

}; // class ImpliedCDSSpreads


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_IMPLIEDCDSSPREADS_H_
