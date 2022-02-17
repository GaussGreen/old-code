/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cbcalls.h
// Purpose:     convertible bond calls class
// Author:      Nabil
// Created:     2004/03/11
// RCS-ID:      $Id: cbcalls.h,v 1.43 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_CBCALLS_H_
#define _ITO33_PRICING_CBCALLS_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/pricing/callprovisions.h"

namespace ito33
{

class ITO33_DLLDECL Date;

namespace finance
{
  class ITO33_DLLDECL CallSchedule;
}

namespace pricing
{

/// Declaration of the CBCalls class
class CBCalls : public CallProvisions
{

public:

  CBCalls(): CallProvisions() { }
  
  CBCalls(const shared_ptr<finance::CallSchedule>& pCall);

  // default dtor is ok

  virtual void GetCallStrikes
               ( const double* pdS, size_t nNbS, 
                 const double* pdNewSharePrices, double* pdValues) const;

  virtual void GetCallStrikesWithoutCoupon
               ( double dTime, const double* pdS, size_t nNbS, 
                 const double* pdNewSharePrices, double* pdValues, 
                 bool bPlus = true ) const;

  double GetStrikeWithoutCoupon(double dTime, bool bMonis = true) const;

  void Clear();

  /**
    Get the strike rate at the specified index.  Returning 0.0 means
    the strike is not defined (must be a yield)

    @param nIdx the call index
    @return the strike rate corresponding to the given index
   */
  double GetStrikeRate(size_t nIdx) const { return m_pdStrikeRates[nIdx]; }

protected:

  /**
     Calculates strike amount at current time.

     @return strike amount
   */
  double GetStrike() const;  

  std::vector<double> m_pdStrikeRates;  


}; // class CBCalls


} // namespace pricing

} // namespace ito33

#endif  //  #ifndef _ITO33_PRICING_CBCALLS_H_
