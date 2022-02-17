/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/eds.h
// Purpose:     Contracts class for EDS
// Created:     2005/01/26
// RCS-ID:      $Id: eds.h,v 1.3 2005/02/08 13:37:11 zhang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/eds.h
   @brief The declaration of the EDS contracts class.
 */

#ifndef _ITO33_PRICING_EDS_H_
#define _ITO33_PRICING_EDS_H_

#include "ito33/dateutils.h"

#include "ito33/finance/eds.h"

#include "ito33/pricing/contract.h"

namespace ito33
{

namespace pricing
{


/// The declaration of the EDS contract class.
class EDS : public Contract 
{
public:
  
  EDS(const finance::EDS& eds);

  // Default dtor is ok 

  double GetRecoveryRate() const { return m_dRecoveryRate; }

  const finance::CashFlowStreamUniform& GetSpreadStream() const 
  {
    return m_spreadStream; 
  }

  double GetBarrier() const { return m_dBarrier; }


private:

  double m_dRecoveryRate;

  const finance::CashFlowStreamUniform& m_spreadStream;

  double m_dBarrier;

  NO_COPY_CLASS(EDS);

}; // class EDS;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_EDS_H_

