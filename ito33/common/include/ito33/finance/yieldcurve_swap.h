/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/yieldcurve_swap.h
// Purpose:     YieldCurveSwap class
// Author:      ZHANG Yunzhi
// Created:     2004/04/11
// RCS-ID:      $Id: yieldcurve_swap.h,v 1.13 2006/08/21 14:27:26 zhang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/yieldcurve_swap.h
    @brief The declaration of YieldCurveSwap class.

    Class for swap curve.
 */

#ifndef _ITO33_FINANCE_YIELDCURVE_SWAP_H_
#define _ITO33_FINANCE_YIELDCURVE_SWAP_H_

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/cashrates.h"

#include "ito33/finance/swaprates.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"

namespace ito33
{

namespace finance
{

/**
    The YieldCurveSwap class allows to construct a yield curve composed of
    interbank deposit rates and par swap rates.
 */
class ITO33_DLLDECL YieldCurveSwap : public YieldCurve
{
public:
  
  /**
      Default constructor.
      By default the cash basis is 30/360 and swap basis is ACT/365.

      NOTE: the Reference date must be specified later.

      @noexport
   */
  YieldCurveSwap()
  {
  }

  /**
      Initializes the new swap curve with its reference date.
      By default the cash basis is 30/360 and swap basis is ACT/365.

      The legs have to be added to the yield curve later.

      @param ReferenceDate the Reference date which must be valid
   */
  YieldCurveSwap(const Date& ReferenceDate)
                    : YieldCurve(ReferenceDate)
  {
  }

  /// @name implement the base class pure virtuals
  //@{
  virtual void GetZeroRates(const double *pdDays,
                            double *pdValues,
                            size_t nPoints) const;

  virtual shared_ptr<YieldCurve> Perturb(double dShift) const;
  //@}

  /// cash rates terms
  void SetCashRates(const CashRates& cashRates)
  { 
    m_cashRates = cashRates;

    // we will have to validate the yield curve data.
    Invalidate();
  }

  /// Swap rate terms
  void SetSwapRates(const SwapRates& swapRates)
  { 
    m_swapRates = swapRates;

    // we will have to validate the yield curve data.
    Invalidate();
  }

  /// cash rates terms
  const CashRates& GetCashRates() const
  { 
    return m_cashRates;
  }

  /// Swap rate terms
  const SwapRates& GetSwapRates() const
  { 
    return m_swapRates;
  }

  /**
      Infers the zero curve from the specified swap curve
      using the bootstrapping methodology.
   */
  shared_ptr<YieldCurveAnnuallyCompounded> BootStrap() const;

  /**
      Removes all the legs from this yield curve.
   */
  void Clear()
  {
    m_cashRates = CashRates();
    m_swapRates = SwapRates();

    Invalidate();
  };

  virtual XML::Tag Dump(const char *name, XML::Tag& tagParent) const;

  // implement base class pure virtuals
  virtual void Visit(YieldCurveVisitor& visitor) const;

private:

  /// DoValidate() helper
  bool CheckData();

  /// the Validate() core
  virtual void DoValidate();

  /// cash rates
  CashRates m_cashRates;

  /// swap rates
  SwapRates m_swapRates;

  shared_ptr<YieldCurveAnnuallyCompounded> m_pYieldCurve;

  NO_COPY_CLASS(YieldCurveSwap);

};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_YIELDCURVE_SWAP_H_
