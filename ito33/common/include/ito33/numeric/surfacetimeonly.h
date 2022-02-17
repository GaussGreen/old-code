/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/surfacetimeonly.h
// Purpose:     Classes representing a surface depending only on time
// Created:     2004/05/07
// RCS-ID:      $Id: surfacetimeonly.h,v 1.13 2006/08/19 22:26:42 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/surfacetimeonly.h
    @brief Numerical double value Surface for price and greeks.

    When the value depend only on time, there will have at most two values.
 */
#ifndef _ITO33_NUMERIC_SURFACETIMEONLY_H_
#define _ITO33_NUMERIC_SURFACETIMEONLY_H_

#include "ito33/vector.h"

#include "ito33/numeric/surfacedouble.h"

namespace ito33
{

namespace numeric
{

/**
    Surface with value depending only on time.
 */
class SurfaceTimeOnly : public SurfaceDouble
{
public:

  /**
      Constructor calls the base class constructor.

      @param pDomain a shared pointer to the numerical domain of the surface.
   */
  SurfaceTimeOnly(const shared_ptr<Domain>& pDomain) 
    : SurfaceDouble(pDomain) { }

  /**
      @name implementation of the interface defined in the base class 
            SurfaceDouble
   */
  //@{

  /**
      Retrieves output values we have for the given time(identified by an
      index) and spots.

      We don't verify here if the index is out of bounds.

      @param nIdxT the index of the time we want to retrieve values for
      @param pdS the spots we want to retrive values for values for 
      @param values the retrieved values 
   */
  void GetValuesAt(size_t nIdxT, const Spots& pdS, Doubles& values) const;

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
      Constructs delta and gamma surfaces from the existing surface.

      Both delta and gamma are computed in one call for efficiency

      @param pDeltaSurface the delta surface (initially empty)
      @param pGammaSurface the gamma surface (initially empty)
   */
  void GetDeltaAndGamma(shared_ptr<SurfaceDouble>& pDeltaSurface,
                        shared_ptr<SurfaceDouble>& pGammaSurface) const;

  void GetThetaBackwardOnly(shared_ptr<SurfaceDouble>& pThetaSurface) const;
  
  // See base class
  shared_ptr<SurfaceDouble> 
  ComputeFiniteDifference(const SurfaceDouble& shiftedSurface,
                          double dInverseShift) const;

  //@}

  /**
      @name SurfaceTimeOnly specific interface
   */
  //@{

  /**
      Retrieves output values we have for the given time(identified by an
      index). 

      The actual convention is lastvalue

      @param nIdxT the index of the time we want to retrieve values for

      @return the retrieved value 
   */
  double GetValueAt(size_t nIdxT) const { return GetLastValueAt(nIdxT); }

  /**
      Retrieves the FirstValue we have for the given time (identified by an 
      index).

      @param nIdxT the index of the time we want to retrieve values for
      @return the retrieved value
   */
  double GetFirstValueAt(size_t nIdxT) const { return m_ppdValues[nIdxT][0]; }

  /**
      Retrieves the LastValue we have for the given time (identified by an 
      index).

      @param nIdxT the index of the time we want to retrieve values for
      @return the retrieved value
   */
  double GetLastValueAt(size_t nIdxT) const
  {
    return m_ppdValues[nIdxT].back(); 
  }

  /**
      Retrieves FirstValue we have for the given time.

      Values will be zero if the time is out of the time mesh. 

      @param dTime the time we want to retrieve values for
      @return the retrieved value 
   */
  double GetFirstValueAt(double dTime) const;

  /**
      Retrieves LastValues we have for the given time and spots.

      Values will be zero if the time is out of the time mesh. 

      @param dTime the time we want to retrieve values for
      @return the retrieved value
   */
  double GetLastValueAt(double dTime) const;

  //@}


  /** 
      @name functions help to construct the surface
   */
  //@{

  /**
      Adds a single value at a given time, so the surface is continuous 
      at current time if no more value added.

      @param dValue the value at current time
   */
  void Add(double dValue) 
  {
    // If the underlying domain has been advanced, then add a new entry.
    // Otherwise, add to the end of the exisiting list, assuming this means
    // multiple data is occuring at the same time. This assumes that the
    // domain is ALWAYS updated before the surface
    size_t nNbTimes = m_ppdValues.size();
    
    if ( m_pDomain->IsValidTimeIndex(nNbTimes) )
    {
      // new entry
      std::vector<double> pdValues;
      
      pdValues.push_back(dValue); 
      
      m_ppdValues.push_back(pdValues);    
    }
    else
    {
      // add to old entries, since domain has not been advanced
      ASSERT_MSG(m_ppdValues[nNbTimes-1].size() < 3, 
                 "More than 3 surface entries at current time");
      
      m_ppdValues.back().push_back(dValue);
    }
  }

  /**
      Adds two value at a given time, so the surface is discontinuous 
      at current time if the two values are not equals.

      @param dValue1 the first value at current time
      @param dValue2 the last value at current time
   */
  void Add(double dValue1, double dValue2)
  { 
    std::vector<double> pdValues(2);
    
    pdValues[0] = dValue1;
    pdValues[0] = dValue2;

    m_ppdValues.push_back(pdValues);
  }

  //@}


private:

  /// the storage of the surface values
  std::vector< std::vector<double> > m_ppdValues;

}; // class SurfaceTimeOnly


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_SURFACETIMEONLY_H_
