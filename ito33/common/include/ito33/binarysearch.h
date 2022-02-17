
/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/binarysearch.h
// Purpose:     binary search
// Author:      Yann d'Halluin
// Created:     22/08/2004
// RCS-ID:      $Id: binarysearch.h,v 1.5 2004/10/28 11:36:39 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_BINARYSEARCH_H_
#define _ITO33_BINARYSEARCH_H_

#include <cstddef> // for size_t

namespace ito33
{

  /**
     Binary search 

     @return 
     if x value lies between A[l] <= x <= A[u] returns index = u
     if x < A[0], index = 1
     if x > A[n-1], index = n-1
   */
  size_t BinSearch(const double* pdA, size_t nSizeA, double dx);

}

#endif // #ifndef _ITO33_BINARYSEARCH_H_
