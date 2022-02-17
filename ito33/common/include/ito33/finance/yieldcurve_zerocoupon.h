/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/yieldcurve_zerocoupon.h
// Purpose:     YieldCurveZeroCoupon class
// Author:      ZHANG Yunzhi
// Created:     2006/08/17
// RCS-ID:      $Id: yieldcurve_zerocoupon.h,v 1.19 2006/08/21 16:41:06 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/yieldcurve_zerocoupon.h
    @brief The declaration of YieldCurveZeroCoupon class.

    Class for annualized annually compounded yield curve.
 */

#ifndef _ITO33_FINANCE_YIELDCURVE_ZEROCOUPON_H_
#define _ITO33_FINANCE_YIELDCURVE_ZEROCOUPON_H_

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/yieldcurveleg.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL YieldCurveAnnuallyCompounded;

/**
    A zero coupon rate curve.

    Note that this class is very similar to YieldCurveAnnuallyCompounded.
    In fact, YieldCurveAnnuallyCompounded is a short version of this class.

    @iterator GetAll
 */
class ITO33_DLLDECL YieldCurveZeroCoupon : public YieldCurve
{
public:
  
  // implement the base class pure virtuals
  virtual void GetZeroRates(const double *pdDays,
                            double *pdValues,
                            size_t nPoints) const;

  virtual shared_ptr<YieldCurve> Perturb(double dShift) const;

  /// the type of the container we use for storing the legs
  typedef std::vector<YieldCurveLeg> Legs;

  /**
      Default constructor.

      NOTE: the Reference date must be specified later.
      
      @noexport
   */
  YieldCurveZeroCoupon() : m_dcc(Date::DayCountConvention_Act365)
  {
    m_legs.reserve(20);
  }

  /**
      Constructor initializes the new yield curve with its reference date.

      The legs have to be added to the yield curve later.

      @param ReferenceDate the Reference date which must be valid
   */
  YieldCurveZeroCoupon(const Date& ReferenceDate)
                     : YieldCurve(ReferenceDate),
                       m_dcc(Date::DayCountConvention_Act365)
  {
    m_legs.reserve(20);
  }

  /**
      Constructor initializes the new yield curve with its reference date.

      The legs have to be added to the yield curve later.

      @param ReferenceDate the Reference date which must be valid
      @param nNumberLegs the approximate number of legs in the curve,
                         used for storage allocation optimization purposes only
   */
  YieldCurveZeroCoupon(const Date& ReferenceDate, size_t nNumberLegs)
                     : YieldCurve(ReferenceDate),
                       m_dcc(Date::DayCountConvention_Act365)
  {
    m_legs.reserve(nNumberLegs);
  }

  // default dtor are ok


  /// @name Operations
  //@{

  /**
      Removes all the legs from this yield curve.

      This method also resets the validated status and so even the curve had
      been Validate()'d before, new legs can be added to it again now.
   */
  void Clear()
  {
    m_legs.clear();

    Invalidate();
  };

  /**
      Adds a leg to the yield curve.

      @param nMaturityDuration the actual number of TimeUnit from reference date
      @param The maturityUnit unit in which the maturity is expressed.
      @param dRate the annually compounded rate of the leg
   */
  void AddLeg(size_t nMaturityDuration, TimeUnit maturityUnit, double dRate);

  /**
       Day-count basis on which the rates are expressed.

       @param dcc Day-count basis on which the rates are expressed.
   */
  void SetDayCountConvention(Date::DayCountConvention dcc);

  //@}

  /**
      Gets direct access to the container of the legs.

      This function calls Validate() first if it hadn't been done yet and so
      may throw if the curve data is invalid. If Validate() had already been
      called before, it doesn't throw.
   */
  const Legs& GetAll() const
  {
    Validate();

    return m_legs;
  }

  /**
       Day-count basis on which the rates are expressed.

       @return Day-count basis on which the rates are expressed.
   */
  Date::DayCountConvention GetDayCountConvention() const
  {
    return m_dcc;
  }

  virtual XML::Tag Dump(const char *name, XML::Tag& tagParent) const;

  // implement base class pure virtuals
  virtual void Visit(YieldCurveVisitor& visitor) const;

private:

  /// the container with all the legs (unsorted) which were added
  Legs m_legs;

  /// the day basis on which the rates are expressed.
  Date::DayCountConvention m_dcc;

  /// Implementation of the class, available after Validate() call.
  shared_ptr<YieldCurveAnnuallyCompounded> m_pImpl;

private:

  /// the Validate() core
  virtual void DoValidate();

  NO_COPY_CLASS(YieldCurveZeroCoupon);
};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_YIELDCURVE_ZEROCOUPON_H_
