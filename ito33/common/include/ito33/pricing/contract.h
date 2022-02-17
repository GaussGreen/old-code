/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/contract.h
// Purpose:     a single contract class
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: contract.h,v 1.9 2006/01/10 16:53:16 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/contract.h
    @brief The declaration of the contract class.

    The base class for a single contract.  
 */

#ifndef _ITO33_PRICING_CONTRACT_H_
#define _ITO33_PRICING_CONTRACT_H_

#include "ito33/pricing/contracts.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivative;
}

namespace pricing
{


/** 
   The declaration of the contract class.

   The base class for a single contract
 */
class Contract : public Contracts
{

public:

  Contract() : Contracts() { }

  Contract(const finance::Derivative& derivative);

  virtual ~Contract() { }

  /**
     Sets the maturity (as a fraction of year) of the contracts
    
     @param dMaturityTime the maturity (as a fraction of year) of the contract
   */
  void SetMaturityTime(double dMaturityTime) { m_dMaturityTime = dMaturityTime; }

  double GetMaturityTime() const { return m_dMaturityTime; }


protected:

  /// The maturity date of the contract
  double m_dMaturityTime;

}; // class Contract;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CONTRACT_H_

