/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/yieldcurve_annuallycompounded.h
// Purpose:     YieldCurveAnnuallyCompounded class
// Author:      Wang
// Created:     2004/06/17
// RCS-ID:      $Id: yieldcurve_annuallycompounded.h,v 1.24 2006/08/21 14:27:26 zhang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/yieldcurve_annuallycompounded.h
    @brief The declaration of YieldCurveAnnuallyCompounded class.

    Class for annualized annually compounded yield curve.
 */

#ifndef _ITO33_FINANCE_YIELDCURVE_ANNUALLYCOMPOUNDED_H_
#define _ITO33_FINANCE_YIELDCURVE_ANNUALLYCOMPOUNDED_H_

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/zerorate.h"
#include "ito33/finance/yieldcurve.h"

namespace ito33
{

namespace finance
{


/**
    A yield curve defined as a collection of legs (day/annually compounded
    rate pairs) with linear interpolation.

    @iterator GetAll
 */
class ITO33_DLLDECL YieldCurveAnnuallyCompounded : public YieldCurve
{
public:
  
  // implement the base class pure virtuals
  virtual void GetZeroRates(const double *pdDays,
                            double *pdValues,
                            size_t nPoints) const;

  virtual shared_ptr<YieldCurve> Perturb(double dShift) const;

  /// the type of the container we use for storing the legs
  typedef std::vector<ZeroRate> Legs;

  /**
      Default constructor.

      NOTE: the Reference date must be specified later.

      @noexport
   */
  YieldCurveAnnuallyCompounded()
  {
    m_legs.reserve(20);
  }

  /**
      Constructor initializes the new yield curve with its reference date.

      The legs have to be added to the yield curve later.

      @param ReferenceDate the Reference date which must be valid
   */
  YieldCurveAnnuallyCompounded(const Date& ReferenceDate)
                             : YieldCurve(ReferenceDate)
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
  YieldCurveAnnuallyCompounded(const Date& ReferenceDate, size_t nNumberLegs)
                             : YieldCurve(ReferenceDate)
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

      @param nDay the actual number of days from reference date
      @param dRate the value of the leg for this day
   */
  void AddLeg(size_t nDay, double dRate)
  {
    if(dRate < 0.)
      ThrowInvalidRate();

    m_legs.push_back(ZeroRate(dRate, nDay));

    // data changed, modify the validation state.
    Invalidate();
  }

  /**
      Sets all of the yield curve legs at once.

      This method automatically calls Validate() and so after calling it no
      more legs can be added to the curve. On the other hand, it can be
      called even if Validate() had already been called because it calls
      Clear() first anyhow.

      It may throw if Validate() throws.

      @param pnDays array of durations of each leg
      @param pdRates array of rates of each leg
      @noexport Java
   */
  void SetLegs(const std::vector<size_t>& pnDays,
               const std::vector<double>& pdRates);

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

  virtual XML::Tag Dump(const char *name, XML::Tag& tagParent) const;

  // implement base class pure virtuals
  virtual void Visit(YieldCurveVisitor& visitor) const;

private:


  /// the container with all the legs (unsorted) which were added
  Legs m_legs;

  // the fields below are valid only after (successful) call to Validate()

  // the arrays containing sorted days
  std::vector<double> m_pdDays;

  /// the array containing rate values corresponding to the days.
  std::vector<double> m_pdRates;

private:

  /// DoValidate() helper
  bool CheckData();

  /// the Validate() core
  virtual void DoValidate();

  /// throws invalid rate exception 
  static void ThrowInvalidRate();

  /// throws exception for invalid data of SetLegs
  static void ThrowInvalidSetLegsData();

  NO_COPY_CLASS(YieldCurveAnnuallyCompounded);
};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_YIELDCURVE_ANNUALLYCOMPOUNDED_H_
