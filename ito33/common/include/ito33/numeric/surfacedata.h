/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/surfacedata.h
// Purpose:     Helper functions to handle a SurfaceDataAtSameTime type
// Author:      ZHANHG Yunzhi
// Created:     2004/10/08
// RCS-ID:      $Id: surfacedata.h,v 1.8 2006/03/22 13:10:22 yann Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/surfacedata.h
    @brief Helper functions to handle a SurfaceDataAtSameTime type

 */
#ifndef _ITO33_NUMERIC_SURFACEDATA_H_
#define _ITO33_NUMERIC_SURFACEDATA_H_

#include "ito33/vector.h"
#include "ito33/debug.h"

namespace ito33
{

namespace numeric
{

/**
  At same time, we have two types of values. 

  standard ones: result after having resolved the PDE at this time step or
                 result after application of events. Note that user wanted
                 result is always the last standard values at this time point.

  special ones: this is not common. If the constraints before and after some
                time point are different, this is the case when a call period
                ends here and another one starts here, before we start a new
                time step, we must apply the new constraint to the prices.
                Thus we call them special values. We also call it "initial
                values for next step".
  */
template <class T>
class SurfaceDataAtSameTime
{
public:

  typedef T DataType;

  typedef std::vector<DataType> DataLine;

  /**
     Adds the values at the given time

     @param values values to be added to the surface
   */
  void Add(const DataLine& values)
  {
    ASSERT_MSG(m_ppdStandardValues.size() < 2, 
           "More than 2 surface entries at current time");

    ASSERT_MSG(!HasInitValuesForNextStep(),
           "Add() must be called after AddInitValuesForNextStep().");

    m_ppdStandardValues.push_back(values);
  }

  /**
    Adds special values at given time

    @param values special values
   */
  void AddInitValuesForNextStep(const DataLine& values)
  {
    ASSERT_MSG(!HasInitValuesForNextStep(),
               "More than 1 Initial values entries at current time");

    ASSERT_MSG
      (
        m_ppdStandardValues.size() > 0,
        "InitValuesForNextStep can't stand alone without standard values."
      );

    m_pdInitValuesForNextStep = values;
  }

  /**
     @name Helper functions to get data out
   */
  //@{

  /**
    Gets the values. This is the last standard values

    @return user required values
    */
  const DataLine& GetValues() const
  {
    SelfCheck();

    return m_ppdStandardValues.back();
  }

  /**
    Gets the first commit values.

    @return the first added values
    */
  const DataLine& GetFirstValues() const
  {
    SelfCheck();

    return m_ppdStandardValues[0];
  }

  /**
    Gets the last commit values.

    @return the last add values
    */
  const DataLine& GetLastValues() const
  {
    SelfCheck();

    return   m_pdInitValuesForNextStep.size() 
           ? m_pdInitValuesForNextStep
           : m_ppdStandardValues.back();
  }

  //@}

  /**
    Gets the size

    @return the size
    */
  size_t SizeOfStandardValues() const
  {
    return m_ppdStandardValues.size();
  }

  /**
    Check for special values.

    @return true, if there are special values at current time, false otherwise
  */
  bool HasInitValuesForNextStep() const
  {
    return m_pdInitValuesForNextStep.size() > 0;
  }

  /**
    Compute finite differences of the current data.

    @param shifted The shifted data values
    @param dInverseShift The inverse of the shift amount
    @param result The computed finite differences
  */
  void ComputeFiniteDifference
          (
            const SurfaceDataAtSameTime<DataType>& shifted,
            double dInverseShift,
            SurfaceDataAtSameTime<DataType>& result
          ) const
  {
    ASSERT_MSG(SizeOfStandardValues() == shifted.SizeOfStandardValues (),
      "Surfaces do not have the same structure.");
    
    ASSERT_MSG
        (
          HasInitValuesForNextStep() == shifted.HasInitValuesForNextStep(),
          "Surfaces do not have the same structure."
        );

    result.m_ppdStandardValues.resize(m_ppdStandardValues.size());
    for(size_t nIdx = 0; nIdx < m_ppdStandardValues.size(); nIdx++)
    {
      ASSERT_MSG
        (
          m_ppdStandardValues[nIdx].size() 
              == shifted.m_ppdStandardValues[nIdx].size(),
          "Surfaces do not have the same structure."
        );  

      result.m_ppdStandardValues[nIdx].resize(
                                          m_ppdStandardValues[nIdx].size());

      for (size_t nXIdx = 0; nXIdx < m_ppdStandardValues[nIdx].size(); nXIdx++)
      {
        result.m_ppdStandardValues[nIdx][nXIdx] 
          = dInverseShift * 
            ( shifted.m_ppdStandardValues[nIdx][nXIdx]
                  - m_ppdStandardValues[nIdx][nXIdx] );
      } // loop over the actual data at a given time (eg. over spot)
    }

    if(HasInitValuesForNextStep())
    {
      result.m_pdInitValuesForNextStep.resize
                      (m_pdInitValuesForNextStep.size());
      for (size_t nXIdx = 0; nXIdx < m_pdInitValuesForNextStep.size(); nXIdx++)
      {
        result.m_pdInitValuesForNextStep[nXIdx]
              = dInverseShift * 
                ( shifted.m_pdInitValuesForNextStep[nXIdx]
                  - m_pdInitValuesForNextStep[nXIdx] );
      } 
    }
  }

  /**
    Apply a function to each data value.

    F must have operator () defined as DataType F(DataType) const.

    @param functor The function object to apply
  */
  template <class F>
  void ApplyFunctor( const F& functor )
  {
    // modify the special data
    if ( HasInitValuesForNextStep() )
    {
      for (size_t nIdx = 0; nIdx < m_pdInitValuesForNextStep.size(); nIdx++)
      {
        DataType currentVal = m_pdInitValuesForNextStep[nIdx];
        m_pdInitValuesForNextStep[nIdx] = functor(currentVal);
      }
    } // if special data

    // modify the normal data
    for (size_t nIdx = 0; nIdx < m_ppdStandardValues.size(); nIdx++)
    {
      for (size_t nIdy = 0; nIdy < m_ppdStandardValues[nIdx].size(); nIdy++)
      {
        DataType currentVal = m_ppdStandardValues[nIdx][nIdy];
        m_ppdStandardValues[nIdx][nIdy] = functor(currentVal);
      }
    } // loop over normal data

  }


private:

  void SelfCheck() const
  {
    ASSERT_MSG(m_ppdStandardValues.size() <= 2,
               "StandardValues can't hold more than 2 results");
  }

  std::vector< DataLine > m_ppdStandardValues;

  DataLine m_pdInitValuesForNextStep;
};

/// Declaration of the SurfaceData class
template <class T>
class SurfaceData
{
public:

  typedef T DataType;

  typedef SurfaceDataAtSameTime<DataType> DataAtSameTime;

  typedef std::vector<DataType> DataLine;

  /**
    default constructor

    @param nSize the size of the surface (time points)
    */
  SurfaceData(size_t nSize = 0) : m_pData(nSize) {}

  /** 
     @name functions to add data
   */
  //@{

  /**
     Inserts the values at the given position

     @param values values to be added to the surface
     @param nWhere position where values should be added
   */
  void Add(const DataLine& values, size_t nWhere)
  {
    ASSERT_MSG
      (
           (m_pData.size() == 0 && nWhere <= m_pData.size() )
        || nWhere == m_pData.size() || nWhere == m_pData.size() - 1,
        "Wrong position to Add surface data. Please check if "
        "each domain.AddTime() call is followed by surface.Add() call."
      );

    if(nWhere == m_pData.size())
      m_pData.push_back( SurfaceDataAtSameTime<DataType>() );

    m_pData[nWhere].Add(values);
  }

  void AddInitValuesForNextStep(const DataLine& values)
  {
    m_pData.back().AddInitValuesForNextStep(values);
  }

  /**
     Inserts the values at the given position

     @param values values to be added to the surface
     @param nWhere position where values should be added
   */
  void Insert(const DataLine& values, size_t nWhere)
  {
#ifndef NDEBUG
    CheckPositionForInsert(nWhere);
#endif
    m_pData[nWhere].Add(values);
  }

  /**
     Inserts the values at the given position

     @param values1 FirstValues to be added to the surface
     @param values2 LastValues to be added to the surface
     @param nWhere position where values should be added
   */
  void Insert(const DataLine& values1,
              const DataLine& values2,
              size_t nWhere)
  {
#ifndef NDEBUG
    CheckPositionForInsert(nWhere);
#endif

    m_pData[nWhere].Add(values1);
    m_pData[nWhere].Add(values2);
  }

  /**
     Add the initial condition values at the given position

     @param pdValues initial condition values to be added to the surface
     @param nWhere position where values should be added
   */
  void InsertInitConditionForNextStep
          (
            const DataLine& values,
            size_t nWhere
          )
  {
#ifndef NDEBUG
    CheckPositionForInsert(nWhere);
#endif

    m_pData[nWhere].AddInitValuesForNextStep(values);
  }

  //@}

  /**
     @name Helper functions to get data out
   */
  //@{
  const DataLine& GetValuesAt(size_t nIdxT) const
  {
    return m_pData[nIdxT].GetValues();
  }

  const DataLine& GetFirstValuesAt(size_t nIdxT) const
  {
    return m_pData[nIdxT].GetFirstValues();
  }

  const DataLine& GetLastValuesAt(size_t nIdxT) const
  {
    return m_pData[nIdxT].GetLastValues();
  }

  //@}

  /** 
     @name functions to modify data
   */
  //@{

  /**
    Apply a function to each data value.

    F must have operator () defined as DataType F(DataType) const.

    @param functor The function object to apply
  */
  template <class F>
  void ApplyFunctor( const F& functor )
  {
    for (size_t nIdxT = 0; nIdxT < m_pData.size(); nIdxT++)
      m_pData[nIdxT].ApplyFunctor( functor );
  }

  //@}

  typedef typename std::vector< DataAtSameTime >::const_iterator
    const_iterator;

  const_iterator begin() const { return m_pData.begin(); }

  const_iterator end() const { return m_pData.end(); }

  
  typedef typename std::vector< DataAtSameTime >::iterator iterator;

  iterator begin() { return m_pData.begin(); }

  iterator end() { return m_pData.end(); }

private:
  
#ifndef NDEBUG
  void CheckPositionForInsert(size_t nWhere)
  {
    ASSERT_MSG
      (
        nWhere < m_pData.size(),
        "Wrong position to insert surface data. Please check if "
        "each domain.AddTime() call is followed by surface.Add() call."
      );
  }
#endif
  
  std::vector< DataAtSameTime > m_pData;

}; // class SurfaceDataAtSameTime

typedef SurfaceData<double> SurfaceDataDouble;


typedef SurfaceData<int> SurfaceDataInt;

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_SURFACEDATA_H_

