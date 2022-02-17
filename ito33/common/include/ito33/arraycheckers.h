/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/arraycheckers.h
// Purpose:     Some useful checkers for array
// Author:      Wang
// Created:     2004/06/04
// RCS-ID:      $Id: arraycheckers.h,v 1.9 2005/02/24 14:33:59 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_ARRAYCHECKERS_H_
#define _ITO33_ARRAYCHECKERS_H_

#include "ito33/vector.h"
#include "ito33/useexception.h"

extern const ito33::Error ITO33_ARRAY_ZEROSIZE;
extern const ito33::Error ITO33_ARRAY_NOTINCREASING;

namespace ito33
{

/**
  Validate non nul array size. Otherwise, throw an exception.

  @param nSize the size to be checked.
  @param pcArrayName name of the array to be appeared in the error message
  @return validated size
  */
inline size_t CheckArraySize(size_t nSize, const std::string& pcArrayName = "")
{
  if(nSize == 0)
  {
    if(pcArrayName == "")
      throw EXCEPTION(ITO33_ARRAY_ZEROSIZE);
    else
    {
      std::string msg(TRANS(pcArrayName.c_str()));
      msg += ": ";
      msg += ITO33_ARRAY_ZEROSIZE.GetMessage();
      throw EXCEPTION_MSG
            (
              ITO33_ARRAY_ZEROSIZE,
              msg.c_str()
            );
    }
  }

  return nSize;
}

/**
   Check a standard container size, throw an exception if empty.

   @param container The container to be checked
   @param pcArrayName name of the array to be appeared in the error message
   
   @return validated container
 */
template <class T>
inline const T& CheckSize(const T& container, const std::string& pcArrayName = "")
{
  if ( container.empty() )
  {
    if ( pcArrayName == "" )
      throw EXCEPTION(ITO33_ARRAY_ZEROSIZE);
    else
    {
      std::string msg(TRANS(pcArrayName.c_str()));
      msg += ": ";
      msg += ITO33_ARRAY_ZEROSIZE.GetMessage();
      throw EXCEPTION_MSG
            (
              ITO33_ARRAY_ZEROSIZE,
              msg.c_str()
            );
    }
  }

  return container;
}

/**
   Validate an array of non negative values

   Throw an exception when the array is empty or has negative element.
   
   @param pdValues the array to be validated
   @param nNumber the size of the array
  @param pcArrayName name of the array to be appeared in the error message

   @return the validated array
 */
const double* 
CheckNonNegativity(const double* pdValues, 
                   size_t nNumber,
                   const std::string& pcArrayName = "");

/**
   Validate an vector of non negative values

   Throw an exception when the vector is empty or has negative element.
   
   @param pdValues the vector to be validated
   @param pcArrayName name of the array to be appeared in the error message

   @return a const reference to the validated vector
 */
const std::vector<double>& 
CheckNonNegativity(const std::vector<double>& pdValues,
                   const std::string& pcArrayName = "");

/**
   Validate an array of values which should be strictly positive

   Throw an exception when the array is empty or has non positive element.
   
   @param pdValues the array to be validated
   @param nNumber the size of the array
   @param pcArrayName name of the array to be appeared in the error message

   @return the validated array
 */
const double* 
CheckPositivity(const double* pdValues,
                size_t nNumber,
                const std::string& pcArrayName = "");

/**
   Validate a vector of values which should be strictly positive

   Throw an exception when the vector is empty or has non positive element.
   
   @param pdValues the vector to be validated
   @param pcArrayName name of the array to be appeared in the error message

   @return a const reference to the validated vector
 */
const std::vector<double>& 
CheckPositivity(const std::vector<double>& pdValues,
                const std::string& pcArrayName = "");

/**
   Validate a sorted (<) array. 
   
   Throw an exception when the array is zero size or has not been sorted.

   @param pValues pointer to the array
   @param nNumber the size of the array
   @param pcArrayName name of the array to be appeared in the error message
   @return the validated pointer
 */
template <class T>
const T* CheckIncreasingOrder(const T* pValues,
                              size_t nNumber,
                              const std::string& pcArrayName = "")
{
  CheckArraySize(nNumber);

  for(size_t nIdx = 1; nIdx < nNumber; nIdx++)
    if(pValues[nIdx - 1] >= pValues[nIdx])
    {
      if(pcArrayName == "")
        throw EXCEPTION_MSG
                (
                ITO33_ARRAY_NOTINCREASING,
                ITO33_ARRAY_NOTINCREASING.GetMessage()
                );
      else
      {
        std::string msg(TRANS(pcArrayName.c_str()));
        msg += ": ";
        msg += ITO33_ARRAY_NOTINCREASING.GetMessage();
        throw EXCEPTION_MSG
                (
                ITO33_ARRAY_NOTINCREASING,
                msg.c_str()
                );
      }
    }

  return pValues;
}

/**
   Validate a sorted vector 
   Throw an exception when the vector is empty or has not been sorted.
   
   @param pValues the vector to be validated
   @param pcArrayName name of the array to be appeared in the error message

   @return a const reference to the validated vector
 */
template <class T>
const std::vector<T>& CheckIncreasingOrder(const std::vector<T>& pValues,
                                           const std::string& pcArrayName = "")
{
  CheckIncreasingOrder(&pValues[0], pValues.size(), pcArrayName);
  return pValues;
}


} // namespace ito33

#endif // #ifndef _ITO33_ARRAYCHECKERS_H_
