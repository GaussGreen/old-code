/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/parbond.h
// Purpose:     contracts class for parbond (backward)
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbond.h,v 1.2 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/parbond.h
   @brief The declaration of the parbond contracts class.
 */

#ifndef _ITO33_PRICING_PARBOND_H_
#define _ITO33_PRICING_PARBOND_H_

#include "ito33/list.h"
#include "ito33/dateutils.h"

#include "ito33/finance/parbond.h"

#include "ito33/pricing/contract.h"

namespace ito33
{

namespace pricing
{


/// The declaration of the (backward) parbond contract class.
class ParBond : public Contract 
{
public:
  
  /**
    The ctor of the class
    
    @param parbond a reference to an object of type finance::ParBond
    @todo{Treat the cross currency case}
   */
  ParBond(const finance::ParBond &parbond);

  // Default dtor is ok 

  double GetRecoveryRate() const { return m_dRecoveryRate; }

  const shared_ptr<finance::CashFlowStream>& GetSpreadStream() const 
  {
    return m_pSpreadStream; 
  }


private:

  double m_dRecoveryRate;

  const shared_ptr<finance::CashFlowStream> m_pSpreadStream;

  NO_COPY_CLASS(ParBond);

}; // class ParBond;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_PARBOND_H_

