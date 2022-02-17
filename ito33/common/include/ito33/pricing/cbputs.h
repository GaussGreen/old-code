/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cbputs.h
// Purpose:     convertible bond puts class
// Author:      Laurence
// Created:     2004/03/11
// RCS-ID:      $Id: cbputs.h,v 1.30 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_CBPUTS_H_
#define _ITO33_PRICING_CBPUTS_H_

#include "ito33/common.h"
#include "ito33/vector.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL PutSchedule;
}

namespace pricing
{

class CBLikeParams;

/**
  Class storing the put data.  Currently:
  - puts occur at discrete times
  - the put is either a strike rate applied to the current principal
    or a yield to maturity

  For simplicity, store vectors for the strikes and yields.  One
  of them must have a non-zero entry corresponding to each put
  date. The other must have a zero value.
*/

class CBPuts
{

public:
  
  CBPuts() { DefaultInit(); }
  
  CBPuts(const shared_ptr<finance::PutSchedule>& pPuts);

  ~CBPuts(){ }

  ///@name methods for accessing puts
  //@{

  /**
    Get the time of the specified put index

    @param nIdx The put index
    @return the put time corresponding to the given index
   */
  double GetTime(size_t nIdx) const { return m_pdTimes[nIdx]; }

  /**
    Get the number of put dates

    @return the number of put dates
   */
  size_t GetNbPuts() const { return m_nNbPuts; } 
  
  /**
    Get the strike rate at the specified index.  Returning 0.0 means
    the strike is not defined (must be a yield)

    @param nIdx the put index
    @return the strike rate corresponding to the given index
   */
  double GetStrikeRate(size_t nIdx) const { return m_pdStrikeRates[nIdx]; }

  /**
    Get the yield to put at the specified index. relevant only when
    strike rate at given index is 0.

    @param nIdx the put index
    @return the yield corresponding to the given index
   */
  double GetYieldToPut(size_t nIdx) const { return m_pdYieldToPut[nIdx]; }

  /**
    Get keep accrued flag

    @return keep accrued flag
    */
  bool GetKeepAccrued() const { return m_bKeepAccrued; }
  
  /**
    Get forfeit coupon flag

    @return forfeit coupon flag
    */
  bool GetForfeitCoupon() const { return m_bForfeitCoupon; }
  
  //@}  // name methods for accessing puts

  /**
     Evaluates Put constraint values, if any, for given spots

     @param pdS given spot array
     @param nNbS @sa number given spots 
     @param pdValues (output) Put constraint values
     @return true when we have Put clause here. otherwise false,
            and does nothing for pdValues
   */
  bool GetPutConstraintValues(const double* pdS, size_t nNbS, 
                              double* pdValues) const;

  /**
     Helper function to remove all the puts.
   */
  void Clear() 
  {
    CBPuts puts;

    *this = puts;
  }

  /**
     Sets pointer to related params object

     @param pParams pointer to params
   */
  void SetParams(const CBLikeParams* pParams)
  {
    m_pParams = pParams;
  }


private:
  
  void DefaultInit()
  {
    m_nNbPuts        = 0;
    m_bKeepAccrued   = true;
    m_bForfeitCoupon = false;
  }

  /// The number of put dates
  size_t m_nNbPuts;

  /// The put dates
  std::vector<double> m_pdTimes;

  /// The strike rates (a zero value means undefined)
  std::vector<double> m_pdStrikeRates;

  /// The yield to put values (a zero value means undefined)
  std::vector<double> m_pdYieldToPut;

  /// Keep accrued flag
  bool m_bKeepAccrued;

  /// Forfeit coupon flag
  bool m_bForfeitCoupon;

  // We save the pointer here, otherwise almost all
  // member function should take params as argument
  /// pointer to params 
  const CBLikeParams* m_pParams;

}; // class CBPuts


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBPUTS_H_
