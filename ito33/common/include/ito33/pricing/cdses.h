/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cdses.h
// Purpose:     contracts class for a list of cds contracts
// Author:      David
// Created:     2004/03/31
// RCS-ID:      $Id: cdses.h,v 1.17 2006/08/21 16:09:46 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cdses.h
    @brief The declaration of the cdses class.

    CDSes contain a list of cds-like contracts 
 */

#ifndef _ITO33_PRICING_CDSES_H_
#define _ITO33_PRICING_CDSES_H_

#include "ito33/sharedptr.h"
#include "ito33/list.h"
#include "ito33/dateutils.h"

#include "ito33/pricing/contracts.h"

namespace ito33
{

class ITO33_DLLDECL Date;

namespace finance
{
  class ITO33_DLLDECL CDSLike;
  class ITO33_DLLDECL CashFlowStreamUniform;
}


namespace pricing
{

/// The declaration of the CDSes class.
class CDSes : public Contracts 
{
public:

  /**
      ctor accepting list of shared pointers.

      @param pCDSes a sorted list of cds-like contracts
   */
  CDSes(const std::list< shared_ptr<finance::CDSLike> >& pCDSes);

  /**
      ctor accepting list of raw pointers.

      @param pCDSes a sorted list of cds-like contracts
   */
  CDSes(const std::list< finance::CDSLike* >& pCDSes);

  // default dtor is ok

  void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const;

  double GetMaturityTime() const 
  { 
    return GetDoubleFrom( m_pMaturityDates.back() );
  }

  /**
      Returns the recovery rate, which is enforced to be the same for all 
      underlying CDS contracts
   */
  double GetRecoveryRate() const { return m_dRecoveryRate; }

  /** 
      Returns the spread stream. For now, only one payment stream is
      assumed, meaning all underlying CDS contracts must use the same
      spread dates (but can have different maturities).
   */
  const finance::CashFlowStreamUniform& 
  GetSpreadStream() const { return *m_pSpreadStream; }


protected:

  /// The recovery rate, enforced to be the same for all CDS contracts
  double m_dRecoveryRate;

  /// The maturity dates for the different CDS contracts
  std::list<Date> m_pMaturityDates;

  // The spread stream
  shared_ptr<finance::CashFlowStreamUniform> m_pSpreadStream;

  NO_COPY_CLASS(CDSes);

}; // class CDSes;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CDSES_H_
