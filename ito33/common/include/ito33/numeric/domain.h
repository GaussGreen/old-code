/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/domain.h
// Purpose:     Domain class defines the times and spots 
// Author:      WANG Xuewen
// Created:     2004/05/03
// RCS-ID:      $Id: domain.h,v 1.13 2005/09/14 14:39:32 nabil Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/domain.h
    @brief Numeric Domain class defines the grid of (time, spot) points 
 */

#ifndef _ITO33_NUMERIC_DOMAIN_H_
#define _ITO33_NUMERIC_DOMAIN_H_

#include <cmath>

#include "ito33/vector.h"
#include "ito33/constants.h"

#include "ito33/finance/domain.h"

namespace ito33
{

namespace numeric
{


/**
   Domain class defines the times and spots.

   Note that there isn't a class for TimeOnly domain since it's 
   rarely(or never ?) used directly either by end user or by numerical code. 
   
   More, a domain without space component sounds a bit strange.

   A time only domain can be treated as a DomeFixedSpaceMesh with 
   a trivially constructed fixed space mesh, so that an uniformed
   interface can be obtained for all kinds of domain.
 */
class Domain : public finance::Domain
{
  
public:

  /// Default constructor
  Domain() : finance::Domain() 
  {
    #ifndef NDEBUG
    {
      m_direction = Direction_Max;
    }
    #endif
  }
  /// Dummy virtual dtor for base class
  virtual ~Domain() { }

  /**
     Generates output dates for end user output. 

     The numerical time mesh can have more than one points for a given date,
     which we would like to hide for the output.

     This function needs to be called before the construction of the output
     surface. 
   */
  void GenerateOutputDates();

  /**
     Gets the index on the time points for the date with index nIdxDate

     @param nIdxDate the index of the date for whom the index at the time array 
                     is needed

     @return the index of the date at the time array
   */
  size_t GetTimeIndexFromDateIndex(size_t nIdxDate) const
  {
    if ( nIdxDate < m_pnIdxDates.size() )
      return m_pnIdxDates[nIdxDate];

    ThrowInvalidDateIndex();

    return 0;
  }

  /**
     Gets the first space mesh at the given time index

     @param nIdxT the index of the time at which the space mesh is required

     @return the space mesh at current time
   */
  virtual const Spots& GetFirstSpaceMeshAt(size_t nIdxT) const = 0;
 
  /**
     Gets the last space mesh at the given time index

     @param nIdxT the index of the time at which the space mesh is required

     @return the space mesh at current time
   */
  virtual const Spots& GetLastSpaceMeshAt(size_t nIdxT) const = 0;

  /**
     Gets the output space mesh at the given time index

     @param nIdxT the index of the time at which the space mesh is required

     @return the space mesh at current time
   */
  virtual const Spots& GetOutputSpaceMeshAt(size_t nIdxT) const = 0;


  /**
    Gets the current time index, ie the largest index

    @return current index
    */
  size_t GetCurrentIndex() const
  {
    return m_pdTimes.size() - 1;
  }
  /** 
     Check if the time index is valid

     @param nIdxT the index to be checked

     @return true if the index is valid, false otherwise
   */
  bool  IsValidTimeIndex(size_t nIdxT) const
  {
    return nIdxT < m_pdTimes.size();
  }

  /** 
     Gets the time at the given time index

     @param nIdxT the index of the time at which the time is required

     @return the time at the specified index
   */
  double GetTimeAt(size_t nIdxT) const
  {
    return m_pdTimes[nIdxT];
  }

  /**
    Gets the time array

    @return time array
    */
  const std::vector<double>& GetTimes() const { return m_pdTimes;  }


protected:

  /**
     Add a time point to the numerical surface

     @param dTime A time point to be added to the surface
   */
  void AddTime(double dTime)
  {
    if ( m_pdTimes.empty() == true)
      m_pdTimes.push_back(dTime);
    else
    {
      // Only add the time if it is new, and larger than the last time
      double dLastTime = m_pdTimes.back();

      // check if the times are in order. Need to know order in which times
      // are added
      if ( fabs( dTime - dLastTime ) > TIMETOLERANCE )
      {
        m_pdTimes.push_back(dTime);
  
        ASSERT_MSG
          (
            m_direction == Direction_Max ? // return always true
                  ( m_direction = dTime > dLastTime + TIMETOLERANCE
                                ? Direction_Forward
                                : Direction_Backward) != Direction_Max
                  :
                  m_direction == Direction_Forward ? 
                    ( dTime > dLastTime + TIMETOLERANCE )
                    :
                    ( dTime < dLastTime - TIMETOLERANCE ),
            "Times out of order in domain_fixedspacemesh"
            
          );
      }
    }
  }

  /// The time points of the domain, can be different than m_dates 
  std::vector<double> m_pdTimes;

  /// Helper vector for the index of m_dates on m_pdTimes
  std::vector<size_t> m_pnIdxDates;

#ifndef NDEBUG
  enum Direction
  {
    Direction_Forward,
    Direction_Backward,
    Direction_Max
  } m_direction;
#endif

private:

  static void ThrowInvalidDateIndex();

}; // class Domain


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_DOMAIN_H_
