/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/surfaceflag.h
// Purpose:     Classes representing a flag depending on time and spot
// Created:     2004/05/05
// RCS-ID:      $Id: surfaceflag.h,v 1.10 2006/08/19 22:26:42 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/surfaceflag.h
    @brief Numerical flag value Surface for price.
 */
#ifndef _ITO33_NUMERIC_SURFACEFLAG_H_
#define _ITO33_NUMERIC_SURFACEFLAG_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/numeric/surfacedata.h"

namespace ito33
{

namespace numeric
{

  class Domain;


/**
    Surface class associates a flag to each point of the associated Domain.
    
    The surface can be discoutinuous on time. At a given time, there may have
    three values that we call respectively FirstValue(s), Value(s), 
    LastValue(s) 
    
    This class is really for backward problem having constraints since there 
    hasn't yet forward problems that have constraints. 

    We won't/can't do interpolation on time, which is different 
    from SurfaceDouble.
 */
class SurfaceFlag
{
public:

  /// type used to represent a flag
  typedef int Flag;

  /// Container of values returned by GetValuesAt()
  typedef std::vector<Flag> Flags;

  /// container of values representing an array of spots
  typedef std::vector<double> Spots;

  /**
      Constructor associates the Surface with a numerical domain.

      @param pDomain a shared pointer to the numerical domain of the surface.
   */
  SurfaceFlag(const shared_ptr<Domain>& pDomain) 
    : m_pDomain(pDomain), m_surface(pDomain->GetTimes().size())
  {
  }

  /**
      Retrieves output values we have for the given time(identified by an
      index) and spots.

      The time index refers to the time array, an exception is thrown if it is 
      out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for values for 
      @param values the retrieved values 
   */
  void GetValuesAt(size_t nIdxT, const Spots& pdS, Flags& values) const;

  /*
    Unlike for the class SurfaceDouble, it doesn't make much sense, neither 
    is required to get the FirstValues and LastValues for any spot at a 
    given time index.
  */
  /**
      @name Helper functions to help the interpolation
   */
  //@{  

  /**
      Retrieves the FirstValues we have for the given time (identified by an 
      index). 

      The time index refers to the time array, an exception is thrown if it is 
      out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
   */
  const Flags& GetFirstValuesAt(size_t nIdxT) const
  {
    return m_surface.GetFirstValuesAt(nIdxT);
  }

  /**
      Retrieves the LastValues we have for the given time (identified by an
      index).

      The time index refers to the time array, an exception is thrown if it is 
      out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
   */
  const Flags& GetLastValuesAt(size_t nIdxT) const
  {
    return m_surface.GetLastValuesAt(nIdxT);
  }

  /**
      Retrieves output values we have for the given time(identified by an
      index) and spots.

      The time index refers to the time array, an exception is thrown if it is 
      out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
      
      @return the flag values at index nIdxT
   */
  const Flags& GetValuesAt(size_t nIdxT) const
  {
    return m_surface.GetValuesAt(nIdxT);
  }

  //@}  



  /**
      Adds the values at the given time.

      @param pFlags values to be added to the surface
      @param bInitValuesForNextStep whether it is at end(begin) of grid
   */
  void Add(const Flags& pFlags, bool bInitValuesForNextStep)
  { 
    if(bInitValuesForNextStep)
      m_surface.AddInitValuesForNextStep(pFlags);
    else
      m_surface.Add(pFlags, m_pDomain->GetCurrentIndex());
  }  

  // we don't need the Insert() functions as in surfaceDouble, as we don't need to
  // calculate Delta, Gamma, theta etc. for flag surface.

  /**
      Gets the domain associated with this surface.

      @return const reference to the internal shared domain pointer. 
   */
  const shared_ptr<Domain>& GetDomain() const
  {
    return m_pDomain;
  }
  
  /**
      Helper function to dumps the numerical result to separate files for
      gnuplot.
   */
  void DumpToFiles() const;


private:

  shared_ptr<Domain> m_pDomain;

  SurfaceDataInt m_surface;

}; // class SurfaceFlag


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_SURFACEFLAG_H_

