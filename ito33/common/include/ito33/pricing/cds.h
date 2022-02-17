/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cds.h
// Purpose:     contracts class for cds (backward)
// Created:     2004/03/02
// RCS-ID:      $Id: cds.h,v 1.18 2006/08/21 14:43:39 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cds.h
    @brief The declaration of the cds contracts class.
 */

#ifndef _ITO33_PRICING_CDS_H_
#define _ITO33_PRICING_CDS_H_

#include "ito33/sharedptr.h"

#include "ito33/pricing/contract.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL CDSLike;
  class ITO33_DLLDECL CashFlowStreamUniform;
}

namespace pricing
{

/// The declaration of the (backward) cds contract class.
class CDS : public Contract 
{
public:
  
  /**
      The ctor takes a financial cds liek object.
    
      @param cds a reference to an object of type finance::CDSLike
      @todo{Treat the cross currency case}
   */
  CDS(const finance::CDSLike &cds);

  // Default dtor is ok 

  /// Thev recovery rate when default occurs
  double GetRecoveryRate() const { return m_dRecoveryRate; }

  /// Gets the spread stream
  const finance::CashFlowStreamUniform& GetSpreadStream() const 
  {
    return *m_pSpreadStream; 
  }


private:

  double m_dRecoveryRate;

  shared_ptr<finance::CashFlowStreamUniform> m_pSpreadStream;

  NO_COPY_CLASS(CDS);

}; // class CDS;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CDS_H_
