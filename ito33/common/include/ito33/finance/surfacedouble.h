/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/surfacedouble.h
// Purpose:     Classes representing a double value depending on time and spot
// Created:     2004/05/04
// RCS-ID:      $Id: surfacedouble.h,v 1.17 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/surfacedouble.h
    @brief Financial double value Surface class.

    A financial SurfaceDouble class represents any double value depending on 
    time and spot.
    
    Such Surface has an associated numerical surface which implements 
    its methods. The domain associated to the surface is obtained from
    the numeric surface's domain.
 */
#ifndef _ITO33_FINANCE_SURFACEDOUBLE_H_
#define _ITO33_FINANCE_SURFACEDOUBLE_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace numeric
{
  class SurfaceDouble;
}

namespace finance
{

  class ITO33_DLLDECL Domain;

/**
    Surface class associating a double value to each point of the associated
    domain.

    @nocreate
 */
class ITO33_DLLDECL SurfaceDouble
{
public:

  /// Container of values returned by GetValuesAt()
  typedef std::vector<double> Doubles;

  /**
      Constructor associates the Surface with a numeric double value surface.

      @param pSurface shared pointer to numeric surface
   */
  SurfaceDouble(const shared_ptr<numeric::SurfaceDouble>& pSurface)
    : m_pSurface(pSurface) { }

  /**
      Retrieves all values for the given date index.

      This is a less efficient but easier to export to other languages
      alternative to the method above. Unless RVO is used, you should prefer
      the above method in C++ code.

      @param nDate the index of the date we want to retrieve values for
      @return the array containing the values for each of the spots
   */
  Doubles GetValuesAt(size_t nDate) const;

  /**
      @internal
      @brief Retrieves all values for the given date index.

      The time index refers to the array returned by Domain::GetDates(), an
      exception is thrown if it is out of bounds. The array elements are the
      values for each domain spot (ie each element of Domain::GetSpots()) for
      the given time.

      This function is more efficient than the overload below unless you
      construct a new vector object from the result returned because it avoids
      copying data.

      @param nDate the index of the date we want to retrieve values for
      @param values the array that we will retrieve values to
      
      @noexport
   */
  void GetValuesAt(size_t nDate, Doubles& values) const;

  /**
      Retrieves a single value for the given date index and spot.

      The date is identified by index in Domain::GetDates() array while spot is
      arbitrary.

      This function is very inefficient compared to GetValuesAt() but may be
      sometimes more convenient to use.

      @param nDate the index of the date we want to retrieve value for
      @param dSpot the value of spot to retrieve value for
      
      @return the value at the given date for the given spot
   */
  double GetValueAt(size_t nDate, double dSpot) const;

  /**
      @internal
      @brief Gets the domain associated with this surface.

      @return a shared pointer to a finance domain

      @noexport
   */
  shared_ptr<Domain> GetDomain() const;

  /**
      @internal
      @brief Gets the numerical surface associated with this surface.

      @return a shared pointer to a numerical surface

      @noexport
   */
  shared_ptr<numeric::SurfaceDouble> GetImpl() const;


private:

  shared_ptr<numeric::SurfaceDouble> m_pSurface;

}; // class SurfaceDouble

inline
SurfaceDouble::Doubles SurfaceDouble::GetValuesAt(size_t nDate) const
{
  Doubles values;
  GetValuesAt(nDate, values);
  return values;
}

} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_SURFACEDOUBLE_H_
