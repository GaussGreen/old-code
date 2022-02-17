/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/surfacegeneral.h
// Purpose:     Classes representing a general surface for price etc
// Created:     2004/05/05
// RCS-ID:      $Id: surfacegeneral.h,v 1.20 2006/08/19 22:26:42 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/surfacegeneral.h
    @brief General Numerical Surface for price and greeks.

    @todo the storage for the values at a given time is not ideal. 
 */
#ifndef _ITO33_NUMERIC_SURFACEGENERAL_H_
#define _ITO33_NUMERIC_SURFACEGENERAL_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/surfacedata.h"
#include "ito33/numeric/surfacedouble.h"
#include "ito33/numeric/domain.h"

namespace ito33
{

namespace numeric
{


/**
    General Numerical Surface for price and greeks.

    We store at most three values at a given time.
 */
class SurfaceGeneral : public SurfaceDouble
{
public:

  /**
      Constructor calls the base class constructor.

      @param pDomain a shared pointer to the numerical domain of the surface.
   */
  SurfaceGeneral(const shared_ptr<Domain>& pDomain) 
    : SurfaceDouble(pDomain), m_surface(pDomain->GetTimes().size())
  {
  }

  /**
      @name implementation of the interface defined in the base class 
            SurfaceDouble
   */
  //@{
  
  /**
      Retrieves output values we have for the given time(identified by an
      index) and spots.

      The actual convention is the middle or the last values.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for
      @param values the retrieved values 
   */
  void GetValuesAt(size_t nIdxT, const Spots& pdS, Doubles& values) const;

  /**
      Retrieves the FirstValues we have for the given time (identified by an 
      index) and spots.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for
      @param values the retrieved values 
   */
  void GetFirstValuesAt
       (size_t nIdxT, const Spots& pdS, Doubles& values) const;

  /**
      Retrieves the LastValues we have for the given time (identified by an
      index) and spots.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for
      @param values the retrieved values 
   */
  void GetLastValuesAt
       (size_t nIdxT, const Spots& pdS, Doubles& values) const;

  /**
      Retrieves FirstValues we have for the given time and spots.

      Values will be zero if the time is out of the time mesh. 

      @param dTime the time we want to retrieve values for
      @param pdS the spots we want to retrive values for
      @param values the retrieved values 
   */
  void GetFirstValuesAt
       (double dTime, const Spots& pdS, Doubles& values) const;

  /**
      Retrieves LastValues we have for the given time and spots.

      Values will be zero if the time is out of the time mesh. 

      @param dTime the time we want to retrieve values for
      @param pdS the spots we want to retrive values for
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
  void GetDeltaAndGamma(shared_ptr<SurfaceDouble>& pDeltaSurface,
                        shared_ptr<SurfaceDouble>& pGammaSurface) const;


  /**
      Constructs theta from the existing surface.

      Similar to delta and gamma, except the derivative is taken
      in the "time" direction.

      @param pThetaSurface the theta surface (initially empty)
   */
  void GetThetaBackwardOnly(shared_ptr<SurfaceDouble>& pThetaSurface) const;

  // See base class
  shared_ptr<SurfaceDouble> 
  ComputeFiniteDifference
  (const SurfaceDouble& shiftedSurface, double dInverseShift) const;

  //@}


  /** 
      @name functions help to construct the surface
   */
  //@{

  /**
      Adds the values at the given time

      @param pdValues values to be added to the surface
      @param bIsInitValuesForNextStep is true if the given values is 
             special and is "initial values for next step"
   */
  void Add(const Doubles& pdValues, bool bIsInitValuesForNextStep = false)
  { 
    if (bIsInitValuesForNextStep)
      m_surface.AddInitValuesForNextStep(pdValues);
    else
      m_surface.Add( pdValues, m_pDomain->GetCurrentIndex());
  }

  //@}


  /**
      @name Helper functions to help the interpolation
   */
  //@{

  const Doubles& GetValuesAt(size_t nIdxT) const
  {
    return m_surface.GetValuesAt(nIdxT);
  }

  const Doubles& GetFirstValuesAt(size_t nIdxT) const
  {
    return m_surface.GetFirstValuesAt(nIdxT);
  }

  const Doubles& GetLastValuesAt(size_t nIdxT) const
  {
    return m_surface.GetLastValuesAt(nIdxT);
  }

  //@}

  /**
      @name Miscellaneous functions 
   */
  //@{

  /**
      Applies a function to each data value.

      F must have operator () defined as double F(double) const.

      @param functor The function object to apply
   */
  template <class F>
  void ApplyFunctor( const F& functor )
  {
    m_surface.ApplyFunctor( functor );
  }

  //@}

  void DumpToFiles() const;


protected:
  
  /** 
      @name functions help to construct the surface
   */
  //@{
  
  /**
      Adds the values at the given time.

      @param pdValues values to be added to the surface
      @param nWhere position where values should be added
   */
  void Insert(const Doubles& pdValues, size_t nWhere)
  { 
    m_surface.Insert( pdValues, nWhere );
  }


  /**
      Adds the values at the given time.

      @param pdValues1 FirstValues to be added to the surface
      @param pdValues2 LastValues to be added to the surface
      @param nWhere position where values should be added
   */
  void Insert(const Doubles& pdValues1, const Doubles& pdValues2, size_t nWhere)
  {
    m_surface.Insert(pdValues1, pdValues2, nWhere);
  }

  /**
      Adds the initial condition values at the given position.

      @param pdValues FirstValues to be added to the surface
      @param nWhere position where values should be added
   */
  void InsertInitConditionForNextStep
            (
              const Doubles& pdValues, 
              size_t nWhere
            )
  {
    m_surface.InsertInitConditionForNextStep(pdValues, nWhere);
  }

  //@}

private:

  SurfaceDataDouble m_surface;

}; // class SurfaceGeneral


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_SURFACEGENERAL_H_
