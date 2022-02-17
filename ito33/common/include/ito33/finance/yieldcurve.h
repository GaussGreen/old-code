/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/yieldcurve.h
// Purpose:     YieldCurve and related classes
// Author:      Vadim Zeitlin
// Created:     12.02.04
// RCS-ID:      $Id: yieldcurve.h,v 1.50 2006/08/21 14:27:26 zhang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/yieldcurve.h
    @brief The declaration of YieldCurve class and some derived classes.

    YieldCurve defines an abstract interface to a yield curve, the concrete
    classes deriving from it provide several possible implementations of this
    interface.
 */

#ifndef _ITO33_FINANCE_YIELDCURVE_H_
#define _ITO33_FINANCE_YIELDCURVE_H_

#include "ito33/date.h"
#include "ito33/sharedptr.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

class ITO33_DLLDECL YieldCurveVisitor;

/**
    Base class for all yield curve implementations.

    Note that the user can define his own yield curves easily simply by
    deriving from this class and implementing the pure virtual functions it
    defines. Doing so is slightly more difficult than using one of the
    predefined concrete classes below, but has the advantage of completely
    avoiding any sort of interpolation for the yield curve points as
    GetZeroRates() could then directly return the exact values.

    @nocreate
 */
class ITO33_DLLDECL YieldCurve
{
public:

  // default copy ctor and assignment operator are ok

  /// Virtual dtor for any base class
  virtual ~YieldCurve() {}


  /**
      The yield curve reference date.

      @return the yield curve reference date
   */
  const Date& GetReferenceDate() const;

  /**
      The yield curve reference date.

      @referenceDate a valid reference date.
   */
  void SetReferenceDate(const Date& referenceDate);

  /**
      Call this function to validate the yield curve data. In general this
      signals the end of adding legs to the yield curve.

      Note that it is not necessary to call Validate(), it will be called by
      the functions accessing the curve data anyhow internally if necessary.
      You may also call Validate() multiple times, it is harmless.
   */
  void Validate() const;

  /**
      @internal
      @brief Gets the yield curve values for the given days.

      For efficiency, we return many values at once instead of requiring a
      separate function call for each date -- this is important for numerical
      code. For the same reason we don't use std::vector<> here: returning a
      new such object each time would be more expensive than simply writing in
      a C arrays like we do now.

      This function is pure virtual and must be implemented by the derived
      classes in the most efficient way possible.

      @param nPoints number of elements in the following arrays
      @param pdDays the number of days to retrieve the values for, as offsets, in
                     fractions of year, from the Reference date
      @param pdValues the array filled with the rate values on return
      
      @noexport
   */
  virtual void GetZeroRates(const double *pdDays,
                            double *pdValues,
                            size_t nPoints) const = 0;

  /**
      @internal
      @brief Convenient wrapper around GetZeroRates() for returning a single value.

      Return a single rate value for the given date.

      @param dDay number of days (as offset to m_ReferenceDate) to retrieve the value for
      @return the rate value for this date
      
      @noexport
   */
  double GetZeroRate(double dDay) const
  {
    double dValue;
    GetZeroRates(&dDay, &dValue, 1);
    return dValue;
  }

  /**
      @internal
      @brief Returns a perturbed yield curve.

      To calculate the derivatives with respect to the yield curve, we must be
      able to shift it by some (small) amount. This method must be implemented
      by the derived classes to do this.

      Note that it doesn't modify this object at all but rather returns a new,
      perturbed, copy of it.

      @param dShift the perturbation parameter, its meaning varies for
                    different derived classes
      @return a new perturbed yield curve
      
      @noexport
   */
  virtual shared_ptr<YieldCurve> Perturb(double dShift) const = 0;
 
  /**
      @internal
      @brief Calculates annualized continously compounded Interest Rate at
             maturity.

      @param dMaturity the maturity, represented as fraction of year

      @return annualized continously compounded Interest Rate. We have
              \f$  e ^ { R_c ( T - t_{ref} ) } 
                 = e^{ \int_{t_{ref}}^T r(s) ds } 
              \f$  or \f$ R_c = log(1 + R_0) \f$   

      @noexport
   */
  double GetContinuousRate(double dMaturity) const;
  
  /**
      @internal
      @brief Calculates annualized continously compounded Interest Rate at
             maturity.

      @param maturity the maturity date 

      @noexport
   */
  double GetContinuousRate(Date maturity) const;

  /**
      @internal
      @brief Calculates Discount Factor at an ordered serie of maturity.
     
      We have \f$ pdDF[i] = e^{ - \int_{t_{ref}}^{T_i} r(s) ds } \f$ 
      with \f$ i = 0 ... N - 1 \f$   

      @param pdMaturities the maturities, represented as fraction of year
      @param pdDF the calculated discount factor, should be allocated by caller
      @param nNb the number of the maturities

      @noexport
   */
  void GetDiscountFactor(const double *pdMaturities, double *pdDF, 
                         size_t nNb) const;

  /**
      @internal
      @brief Calculates compound factor (Inverse of discount factor) at an 
             ordered serie of maturity.

      We have \f$ pdCF[i] = e^{ \int_{t_{ref}}^{T_i} r(s) ds } \f$ 
      with \f$ i = 0 ... N - 1 \f$
    
      @param pdMaturities the maturities, represented as fraction of year
      @param pdCF the calculated compound factor, should be allocated by caller
      @param nNb the number of the maturities

      @noexport
   */
  void GetCompoundFactor(const double *pdMaturities, double *pdCF, 
                         size_t nNb) const;
  
  /**
      @internal 
      @brief Calculates forward discount factor between two maturities.
    
      @param dMaturity1 the first maturity
      @param dMaturity2 the second maturity, should not be less than the first

      @return \f$ e^{ - \int_{T_1}^{T_2} r(s) ds } \f$

      @noexport
   */
  double GetForwardDiscountFactor(double dMaturity1, double dMaturity2) const;
 
  /**
      @internal 
      @brief Calculates forward discount factor between two maturities.
    
      @param date1 the first maturity
      @param date2 the second maturity, should not be less than the first

      @return \f$ e^{ - \int_{T_1}^{T_2} r(s) ds } \f$

      @noexport
   */
  double GetForwardDiscountFactor(Date date1, Date date2) const;

  /**
      @internal 
      @brief Calculates forward discount factor between an ordered serie of
             maturities.

      We have \f$ pdDFF[i] = e^{ - \int_{T_i}^{T_{i + 1}} r(s) ds } \f$ 
      with \f$ i = 0 ... N - 2 \f$

      @param pdMaturities the maturities, represented as fraction of year
      @param pdDFF the calculated forward discount factor, 
                   should be allocated by caller
      @param nNb the number of the maturities

      @noexport
   */
  void GetForwardDiscountFactor(const double *pdMaturities, double *pdDFF, 
                                size_t nNb) const;

  /**
      @internal
      @brief Dumps all data of this yield curve in XML format.

      Unlike Derivative::Dump() and Session::Dump(), this method takes a name
      as we typically have several yield curves and so have to use different
      names for them instead of always the same as in the other cases.

      @param name the name of the tag to create
      @param tagParent the parent tag under which our tag(s) should be created
      @return the tag we created (so that caller could add more stuff to it
              before closing it)
      
      @noexport
   */
  virtual XML::Tag Dump(const char *name, XML::Tag& tagParent) const = 0;

  /**
      @internal
      @brief This method is part of the implementation of the visitor pattern.

      Visitor pattern allows to do different things for different kinds of
      derivatives without adding virtual functions to do all of them in the
      base Derivative class itself nor doing the dreaded switchs on typeid of
      the concrete type. Instead, define a new class deriving from
      DerivativeVisitor and do whatever is required in its methods.

      @noexport
   */
  virtual void Visit(YieldCurveVisitor& visitor) const = 0;

protected:

  /**
      Constructor initializes the new yield curve with a reference date.

      All days in other methods are counted from the Reference.

      @param ReferenceDate a valid Reference date
   */
  YieldCurve(const Date& ReferenceDate);

  /**
      Constructor initializes the new yield curve without specifying the
      reference date. The reference date must be set later.
   */
  YieldCurve();

  /**
      Call this function to be sure that Reference date has been defined.
   */
  void CheckReferenceDate() const;

  /**
      Invalidates the yield curve
   */
  void Invalidate() { m_bValidated = false; }

  /**
      Whether the yield curve is valid.

      @return whether the yield curve is valid
   */
  bool IsValid() const { return m_bValidated; }

  /**
      Really validates the yield curve data. This function implements Validate().

      By default it does nothing.
   */
  virtual void DoValidate() {}

private:

  // the Reference date, all days are counted from it
  Date m_ReferenceDate;

  // the Reference date, represented as a fraction of year
  double m_dReferenceDate;

  /// true if Validate() has been called
  bool m_bValidated;

}; // class YieldCurve;

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_YIELDCURVE_H_
