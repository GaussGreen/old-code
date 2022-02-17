/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/issuer.h
// Purpose:     Issuer class declaration
// Created:     18.12.2003
// RCS-ID:      $Id: issuer.h,v 1.29 2006/06/21 16:09:00 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/issuer.h
    @brief Declaration of the issuer class.
 */

#ifndef _ITO33_FINANCE_ISSUER_H_
#define _ITO33_FINANCE_ISSUER_H_

#include "ito33/vector.h"
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

/**
    Class defining the issuer of an instrument.
 */
class ITO33_DLLDECL Issuer 
{

public:
  /**
      Constructs an empty Issuer object.

      Sets a default value(Jan 1) for the fiscal year start date.
   */
  Issuer() 
  {
    m_fiscalYearStartDate = Date(2000, Date::Jan, 1);
  }

  /**
      The fiscal year start date.

      @param fiscalYearStartDate the fiscal year start date
   */
  void SetFiscalYearStartDate(Date fiscalYearStartDate);

  /**
      Sets the default intensity of this issuer.

      The default intensity will only be used for an exchangeable bond, and
      is a piecewise constant function of time.

      @param pDates a vector of increasing dates
      @param pdDefaultIntensities a vector of default intensities
   */
  void SetDefaultIntensity(const std::vector<Date>& pDates,
                           const std::vector<double>& pdDefaultIntensities);

  /**
      The fiscal year start date.

      @return The fiscal year start date
   */
  const Date& GetFiscalYearStartDate() const 
  { 
    return m_fiscalYearStartDate; 
  }
  
  /**
      Checks if the default intensity of the issuer is defined.

      @return true if the default instensity is defined, false otherwise
   */
  bool HasDefaultIntensity() const;

  /**
      Gets the time array of the default intensity (a piecewise constant
      function of time).

      @return the time array of the default intensity function
   */
  const std::vector<Date>& GetDatesOfDefaultIntensity() const;

  /**
      Gets the value array of the default intensity (a piecewise constant
      function of time).

      @return the value array of the default intensity function
   */
  const std::vector<double>& GetValuesOfDefaultIntensity() const;

  /**
      @internal

      @noexport
   */
  void Dump(XML::Tag& tagParent) const;

private:

  /// The fiscal year start date
  Date m_fiscalYearStartDate;

  /// Time array of the default intensity
  std::vector<Date> m_pDates;

  /// Value array of the default intensity
  std::vector<double> m_pdDefaultIntensities;

}; // class Issuer


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_ISSUER_H_
