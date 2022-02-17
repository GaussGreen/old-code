/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/removeduplicates.h
// Purpose:     Some useful functions for the containers to remove duplicates
// Author:      Nabil
// Created:     2004/02/23
// RCS-ID:      $Id: removeduplicates.h,v 1.5 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/removeduplicates.h
    @brief Some useful functions for the containers to remove duplicates

 */

#ifndef _ITO33_REMOVEDUPLICATES_H_
#define _ITO33_REMOVEDUPLICATES_H_

#include "ito33/beforestd.h"
#include <list>
#include <cmath>
#include <algorithm>
#include "ito33/afterstd.h"

namespace ito33
{

/**
    This template function removes duplicates for a sorted list of elements. 

    The type of the elements is T and the predicat used must be the
    predicat defined for this same type.
    
    @param pElements List of elements
    @param predicat the predicat (the comparison criterion)

    @remarks{ The container must be sorted. }
 */
template<class T, class Predicat>
void RemoveDuplicates(std::list<T>& pElements, Predicat predicat)
{
  std::list<T>::iterator
    Result;
  
  Result = std::unique(pElements.begin(), pElements.end(), predicat);

  //Use erase to remove no needed elements.  
  pElements.erase(Result, pElements.end());  
}

} //namespace ito33

#endif //_ITO33_REMOVEDUPLICATES_H_
