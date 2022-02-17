/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/surfaceflag.h
// Purpose:     Classes representing a flag value depending on time and spot
// Author:      Wang
// Created:     2004/05/05
// RCS-ID:      $Id: surfaceflag.h,v 1.17 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/surfaceflag.h
    @brief Financial flag value Surface class.

    A financial SurfaceFlag class represents any flag value depending on 
    time and spot.

    A flag indicates the constraint type at a (time, spot) point.
    
    Such Surface has an associated numerical surface which implements 
    its methods. The domain associated to the surface is obtained from
    the numeric surface's domain.
 */
#ifndef _ITO33_FINANCE_SURFACEFLAG_H_
#define _ITO33_FINANCE_SURFACEFLAG_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace numeric
{
  class SurfaceFlag;
}

namespace finance
{

  class ITO33_DLLDECL Domain;

/**
    Surface class associating a flag value to each point of the associated
    domain.

    @nocreate
 */
class ITO33_DLLDECL SurfaceFlag
{
public:

  /// the type used to represent a flag
  typedef int Flag;
  
  /// Container of values returned by GetValuesAt()
  typedef std::vector<int /* and not Flag as c2a chokes on this */ > Flags;


  /**
      Constructor associates the Surface with a numeric flag value surface.

      @param pSurface shared pointer to numeric surface

      @noexport
   */
  SurfaceFlag(const shared_ptr<numeric::SurfaceFlag>& pSurface)
    : m_pSurface(pSurface) { }
 
  /**
      Retrieves all values for the given date index.

      This is a less efficient but easier to export to other languages
      alternative to the method above. Unless RVO is used, you should prefer
      the above method in C++ code.

      @param nDate the index of the date we want to retrieve values for
      @return the array containing the values for each of the spots
   */
  Flags GetValuesAt(size_t nDate) const;

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
  void GetValuesAt(size_t nDate, Flags& values) const;

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
  Flag GetValueAt(size_t nDate, double dSpot) const;

  /**
      @internal
      @brief Gets the domain associated with this surface.

      @return a shared pointer to a finance domain

      @noexport
   */
  shared_ptr<Domain> GetDomain() const;


private:

  shared_ptr<numeric::SurfaceFlag> m_pSurface;

}; // class SurfaceFlag

inline
SurfaceFlag::Flags SurfaceFlag::GetValuesAt(size_t nDate) const
{
  Flags values;
  GetValuesAt(nDate, values);
  return values;
}


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_SURFACEFLAG_H_
