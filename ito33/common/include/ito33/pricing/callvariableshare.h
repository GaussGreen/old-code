/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/callvariableshare.h
// Purpose:     variable share call class for PEPS-like
// Author:      Wang
// Created:     2004/08/23
// RCS-ID:      $Id: callvariableshare.h,v 1.4 2006/03/22 13:10:23 yann Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_CALLVARIABLESHARE_H_
#define _ITO33_PRICING_CALLVARIABLESHARE_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/date.h"

#include "ito33/pricing/callprovisions.h"
#include "ito33/pricing/mandatory_payoff_structure.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL GeneralizedPEPSLike;
}

namespace pricing
{
/// Declaration of the call variable share
class CallVariableShare : public CallProvisions
{

public:

  CallVariableShare(): CallProvisions() { }
  
  CallVariableShare(const finance::GeneralizedPEPSLike& peps);

  // Default dtor is ok

  //  ---------  implement virtual functions  --------------
  virtual void GetCallStrikes
               ( const double* pdS, size_t nNbS, 
                 const double* pdNewSharePrices, double* pdStrikes ) const;

  virtual void GetCallStrikesWithoutCoupon
               ( double dTime, const double* pdS, size_t nNbS, 
                 const double* pdNewSharePrices, double* pdValues,
                 bool bPlus = true ) const;

  virtual void Clear();
  //  ---------  implement virtual functions  -------------- end

private:

  std::vector<MandatoryPayoffStructure> m_pData;


}; // class CallVariableShare



} // namespace pricing

} // namespace ito33

#endif  //  #ifndef _ITO33_PRICING_CALLVARIABLESHARE_H_
