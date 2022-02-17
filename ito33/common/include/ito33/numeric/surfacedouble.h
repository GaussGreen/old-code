/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/surfacedouble.h
// Purpose:     Classes representing a double depending on time and spot
// Created:     2004/05/04
// RCS-ID:      $Id: surfacedouble.h,v 1.13 2006/08/19 22:26:42 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/surfacedouble.h
    @brief Numerical double value Surface for price and greeks.
 */
#ifndef _ITO33_NUMERIC_SURFACEDOUBLE_H_
#define _ITO33_NUMERIC_SURFACEDOUBLE_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

namespace ito33
{

namespace numeric
{

  class Domain;


/**
    Surface class associates a double to each point of the associated Domain.
    
    The surface can be discontinuous on time. At a given time, there may have
    three values that we call respectively FirstValue(s), Value(s), 
    LastValue(s) 
    
    We are not using Left and Right because difference between
    backward and forward. FirstValue for backward is RightValue, but is 
    LeftValue for forward.
 */
class SurfaceDouble
{
public:

  /// The container used for the spots
  typedef std::vector<double> Spots;

  /// Container of values returned by GetValuesAt()
  typedef std::vector<double> Doubles;

  /**
      Constructor associates the Surface with a numerical domain.

      @param pDomain a shared pointer to the numerical domain of the surface.
   */
  SurfaceDouble(const shared_ptr<Domain>& pDomain) : m_pDomain(pDomain) { }

  /// Dummy virtual dtor for base class
  virtual ~SurfaceDouble() { }
  
  /**
      Retrieves output values we have for the given time(identified by an
      index) and spots.

      The time index refers to the time array. an exception is thrown if it is 
      out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for 
      @param values the retrieved values 
   */
  virtual void GetValuesAt
               (size_t nIdxT, const Spots& pdS, Doubles& values) const = 0;

  /**
      Retrieves the FirstValues we have for the given time (identified by an 
      index) and spots.

      The time index refers to the time array, an exception is thrown if it is 
      out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for
      @param values the retrieved values 
   */
  virtual void GetFirstValuesAt
               (size_t nIdxT, const Spots& pdS, Doubles& values) const = 0;

  /**
      Retrieves the LastValues we have for the given time (identified by an
      index) and spots.

      The time index refers to the time array, an exception is thrown if it is 
      out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for  
      @param values the retrieved values 
   */
  virtual void GetLastValuesAt
               (size_t nIdxT, const Spots& pdS, Doubles& values) const = 0;

  /*
    Unlike the class SurfaceFlag, we don't have function as
    virtual const Doubles & GetFirstValuesAt(size_t nIdxT) const = 0;
    since there may have surface that the values depend only on time, 
    so that the values need not to be stored as Doubles.

    Anyway, What are really needed is the following general functions 
    from numerical point of view.
  */

  /**
      Retrieves FirstValues we have for the given time and spots.

      Values will be zero if the time is out of the time mesh. 

      @param dTime the time we want to retrieve values for
      @param pdS the spots we want to retrive values for 
      @param values the retrieved values 
   */
  virtual void GetFirstValuesAt
               (double dTime, const Spots& pdS, Doubles& values) const = 0;

  /**
      Retrieves LastValues we have for the given time and spots.

      Values will be zero if the time is out of the time mesh. 

      @param dTime the time we want to retrieve values for
      @param pdS the spots we want to retrive values for 
      @param values the retrieved values 
   */
  virtual void GetLastValuesAt
               (double dTime, const Spots& pdS, Doubles& values) const = 0;

  /**
      Gets the domain associated with this surface.

      @return const reference to the internal shared domain pointer. 
   */
  const shared_ptr<Domain>& GetDomain() const
  {
    return m_pDomain;
  }


  /**
      Constructs delta and gamma surfaces from the existing surface.

      Both delta and gamma are computed in one call for efficiency.

      @param pDeltaSurface the delta surface (initially empty)
      @param pGammaSurface the gamma surface (initially empty)
   */
  virtual void 
  GetDeltaAndGamma(shared_ptr<SurfaceDouble>& pDeltaSurface,
                   shared_ptr<SurfaceDouble>& pGammaSurface) const = 0;

  /**
      Constructs theta from the existing surface.

      Similar to delta and gamma, except the derivative is taken
      in the "time" direction.

      @param pThetaSurface the theta surface (initially empty)
   */
  virtual void 
  GetThetaBackwardOnly(shared_ptr<SurfaceDouble>& pThetaSurface) const = 0;
  
  /**
      Performs a finite difference calculation.

      Both surfaces are assumed to have the exact same structure.

      @param shiftedSurface The values of the shifted surface
      @param dInverseShift The inverse of the finite difference shift value 
      @return surface containing the derivative values computed via finite diffs 
   */
  virtual shared_ptr<SurfaceDouble> 
  ComputeFiniteDifference
  (const SurfaceDouble& shiftedSurface, double dInverseShift) const = 0;

  /**
      Helper function to dumps all the numerical result to separate files for
      gnuplot.
   */
  virtual void DumpToFiles() const { }


protected:

  shared_ptr<Domain> m_pDomain;

}; // class SurfaceDouble


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_SURFACEDOUBLE_H_
