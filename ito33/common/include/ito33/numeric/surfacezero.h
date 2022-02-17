/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/surfacezero.h
// Purpose:     Classes representing a zero value on time and spot
// Created:     2004/05/05
// RCS-ID:      $Id: surfacezero.h,v 1.11 2006/08/19 22:26:42 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/surfacezero.h
    @brief Numerical double value Surface for price and greeks.

    Hope that the name won't create confusion. The name SurfaceDoubleZero
    just seems so ridiculous and there won't have a corresponding version
    for flag value surface.
 */
#ifndef _ITO33_NUMERIC_SURFACEZERO_H_
#define _ITO33_NUMERIC_SURFACEZERO_H_

#include "ito33/vector.h"

#include "ito33/numeric/surfacedouble.h"

namespace ito33
{

namespace numeric
{

/**
    Surface with only zero as value for every (time, spot) point.

    This class doesn't make much sense from financial point of view, 
    but we can use it to create easily a surface for delta, gamma of time only 
    problem. 
    
    A general flat surface class is not implemented, because it doesn't make
    much sense, and it will complicate the implementation for the value outside
    of the time mesh(should be zero normally). 
 */
class SurfaceZero : public SurfaceDouble
{
public:

  /**
      Constructor calls the base class constructor.

      @param pDomain a shared pointer to the numerical domain of the surface.
   */
  SurfaceZero(const shared_ptr<Domain>& pDomain) : SurfaceDouble(pDomain) { }

  /**
      Retrieves output values we have for the given time(identified by an
      index) and spots.

      We don't verify here if the index is out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for values for 
      @param values the retrieved values 
   */
  virtual void
  GetValuesAt(size_t nIdxT, const Spots& pdS, Doubles& values) const;

  /**
      Retrieves the FirstValues we have for the given time (identified by an 
      index) and spots.

      We don't verify here if the index is out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for values for 
      @param values the retrieved values 
   */
  void GetFirstValuesAt
       (size_t nIdxT, const Spots& pdS, Doubles& values) const;

  /**
      Retrieves the LastValues we have for the given time (identified by an
      index) and spots.

      We don't verify here if the index is out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for values for 
      @param values the retrieved values 
   */
  void GetLastValuesAt
       (size_t nIdxT, const Spots& pdS, Doubles& values) const;

  /**
      Retrieves FirstValues we have for the given time and spots.

      Values will be zero if the time is out of the time mesh. 

      @param dTime the time we want to retrieve values for
      @param pdS the spots we want to retrive values for values for 
      @param values the retrieved values 
   */
  void GetFirstValuesAt
       (double dTime, const Spots& pdS, Doubles& values) const;

  /**
      Retrieves LastValues we have for the given time and spots.

      Values will be zero if the time is out of the time mesh. 

      @param dTime the time we want to retrieve values for
      @param pdS the spots we want to retrive values for values for 
      @param values the retrieved values 
   */
  void GetLastValuesAt
       (double dTime, const Spots& pdS, Doubles& values) const;

  /**
      Constructs delta and gamma surfaces from the existing surface

      Both delta and gamma are computed in one call for efficiency

      @param pDeltaSurface the delta surface (initially empty)
      @param pGammaSurface the gamma surface (initially empty)
   */
  virtual void 
  GetDeltaAndGamma(shared_ptr<SurfaceDouble>& pDeltaSurface,
                   shared_ptr<SurfaceDouble>& pGammaSurface) const;


  /**
      Constructs theta from the existing surface.

      Dummy implementation to allow the code to compile.

      Similar to delta and gamma, except the derivative is taken
      in the "time" direction.

      @param pThetaSurface the theta surface (initially empty)
   */
  void GetThetaBackwardOnly(shared_ptr<SurfaceDouble>& pThetaSurface) const
  {
    pThetaSurface = make_ptr( new SurfaceZero(m_pDomain) );
  }  
  
  /// See base class
  shared_ptr<SurfaceDouble> 
  ComputeFiniteDifference(const SurfaceDouble& shiftedSurface,
                          double dInverseShift) const;

}; // class SurfaceZero


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_SURFACEZERO_H_
