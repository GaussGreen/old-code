/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/timestruct.h
// Purpose:     time array classes
// Author:      (z)
// Created:     03/11/04
// RCS-ID:      $Id: timestruct.h,v 1.12 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NUMERIC_TIMESTRUCT_H_
#define _ITO33_NUMERIC_TIMESTRUCT_H_

#include "ito33/useexception.h"
#include "ito33/array.h"

namespace ito33
{

namespace numeric
{

/**
     time array class
  */
class TimeStruct
{
public:

  /**
     ctor
     
     @param pdTimes time array
     @param nNumber number of elements
  */
  TimeStruct(const double *pdTimes, size_t nNumber);

  /// dtor
  virtual ~TimeStruct()
  {
  }

  /**
     Set Initial time and direction

     @param bIsBackward search direction
     @param dTime initial time value
    */
  virtual void SetInitialTime(bool bIsBackward, double dTime) = 0;
  
  /**
     Get the index of the interval where find given time

     @param dTime given time value
   */
  virtual size_t GetIndex(double dTime) = 0;

  /** 
      Gets the size of time array

      @return the size
    */
  size_t GetSize() const { return m_nNumber; }
  
  /**
      Gets the time values

      @param pBegin iterator pointing to the beginning of times
      @param pEnd iterator pointing one past the end of times
   */
  template <typename T>
  void GetTimes(T pBegin, T pEnd) const
  {
    size_t n = 0;
    for ( T pTimes = pBegin; pTimes != pEnd; ++pTimes, ++n)
    {
      *pTimes = m_pdTimes[n];
    }
  }
protected:
  /// time array
  Array<double> m_pdTimes;
  /// number of elements in the time array
  size_t m_nNumber;

  /// direction
  bool m_bIsBackward;

  /// index
  size_t m_nIndex;
};

/**
  time array where the function is piece-wise constant

  the function f(s) defined by a time array t[n] and a value array v[n], n<N is

  f(s) = v[0]      if s < t[0]
         v[n]      if t[n-1] <= s < t[n] for 0 <= n < N
         v[N-1]    if t[N - 1] <= s
  */
class TimeStructConstant : public TimeStruct
{
public:
  
  /**
     ctor
     
     @param pdTimes time array
     @param nNumber number of elements
  */
  TimeStructConstant(const double *pdTimes, size_t nNumber) :
      TimeStruct(pdTimes, nNumber) {}
  
  /**
     Set Initial time and direction

     @param bIsBackward search direction
     @param dTime initial time value
    */
  void SetInitialTime(bool bIsBackward, double dTime)
  {
    m_bIsBackward = bIsBackward;

    if(m_bIsBackward)
    {
      for(m_nIndex= m_nNumber - 1; 
        m_nIndex >= 1 && dTime < m_pdTimes[m_nIndex - 1]; m_nIndex--)
          ;
    }
    else
    {
      for(m_nIndex= 0; 
        m_nIndex < m_nNumber && dTime >= m_pdTimes[m_nIndex]; m_nIndex++)
          ;
    }
  }

  /**
     Get the index of the interval where find given time

     @param dTime given time value
   */
  size_t GetIndex(double dTime)
  {
    if(m_bIsBackward)
    {
      for(; m_nIndex >= 1 && dTime < m_pdTimes[m_nIndex - 1]; m_nIndex--)
          ;
    }
    else
    {
      for(; m_nIndex < m_nNumber && dTime >= m_pdTimes[m_nIndex]; m_nIndex++)
          ;
    }

    return m_nIndex;
  }
};


} // namespace numeric

} // namespace numeric

#endif // _ITO33_NUMERIC_TIMESTRUCT_H_
