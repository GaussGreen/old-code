/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/callfixedshare.h
// Purpose:     Fixed share call class for PEPS-like
// Author:      Wang
// Created:     2004/08/23
// RCS-ID:      $Id: callfixedshare.h,v 1.12 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_CALLFIXEDSHARE_H_
#define _ITO33_PRICING_CALLFIXEDSHARE_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/date.h"

#include "ito33/pricing/callprovisions.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL CallFixedShare;
  class ITO33_DLLDECL GeneralizedPEPSLike;
}

namespace pricing
{


class CallFixedShare : public CallProvisions
{

public:

  CallFixedShare(): CallProvisions() { }
  
  CallFixedShare(const shared_ptr<finance::CallFixedShare>& pCall,
                 Date maturityDate);

  CallFixedShare(const finance::GeneralizedPEPSLike& peps);

  // Default dtor is ok

  //  ---------  implement virtual functions  --------------
  virtual void GetCallStrikes
               (const double* pdS, size_t nNbS, 
                const double* pdNewSharePrices, double* pdStrikes) const;

  virtual void GetCallStrikesWithoutCoupon
  (double dTime, const double* pdS, size_t nNbS, 
   const double* pdNewSharePrices, double* pdValues, 
   bool bPlus = true) const;

  virtual void Clear();
  //  ---------  implement virtual functions  -------------- end

private:

  std::vector<double> m_pdRatios;

}; // class CallFixedShare



} // namespace pricing

} // namespace ito33

#endif  //  #ifndef _ITO33_PRICING_CALLFIXEDSHARE_H_
