/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/domain.h
// Purpose:     Class defines the grid of (time, spot) points
// Author:      Vadim Zeitlin
// Created:     2004-05-03
// RCS-ID:      $Id: domain.h,v 1.17 2006/01/09 09:47:55 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/domain.h
    @brief Domain class defines the grid of (time, spot) points.
 */
#ifndef _ITO33_FINANCE_DOMAIN_H_
#define _ITO33_FINANCE_DOMAIN_H_

#include "ito33/date.h"
#include "ito33/vector.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    Domain class defines the grid of points for which surface values are
    defined.

    @nocreate
 */
class ITO33_DLLDECL Domain
{
public:
  
  /// The container used for the dates
  typedef std::vector<Date> Dates;

  /// The container used for the spots
  typedef std::vector<double> Spots;

  // default and copy ctors and operator=() are ok

  /// Dummy virtual dtor for base class
  virtual ~Domain() { }

  /**
      Gets all points of the time axis for which we have values.

      @return reference to internal array of dates
   */
  const Dates& GetDates() const { return m_dates; }

  /**
      The underlying share prices for which we currently have the values.

      If spots are not strictly increasing, an exception is thrown.

      @param spots the spots we're interested in in strictly increasing order
   */
  void SetUnderlyingSharePrices(const Spots& spots);

  /**
      The underlying share prices for which we currently have the values.

      If SetSpots() had been called, the array returned will contain the same
      elements as were passed to SetSpots(), so it isn't very useful to call
      this function in this case. OTOH, if SetSpots() had not been called, then
      the class determines the spots automatically and you need to use this
      method to retrieve them.

      Also notice that the reference returned by this method will be changed if
      SetSpots() is called [again], so you need to copy the values if you want
      to keep them.

      @return reference to internal array of spots
   */
  const Spots& GetUnderlyingSharePrices() const; 


protected:

  /// the (fixed) array of dates
  Dates m_dates;

  /// the array of spots
  Spots m_spots;

}; // class Domain


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_DOMAIN_H_
