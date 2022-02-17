/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/binarysearch.cpp
// Purpose:     binary search
// Author:      Yann d'Halluin
// Created:     22/08/2004
// RCS-ID:      $Id: binarysearch.cpp,v 1.6 2006/07/21 21:55:18 dave Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/binarysearch.h"
#include "ito33/debug.h"

namespace ito33 
{
size_t BinSearch(const double *pdA, const size_t nSizeA, const double dx)
{
   ASSERT_MSG( pdA !=  0,"Binary Search error, input array null");
   ASSERT_MSG( nSizeA  >= 2, "Binary Search error, array should be greater than 2"); 

   size_t l, u, m;
     // x value lies between
     // A[l] <= x <= A[u]
     // returns index = u
     //  if x < A[0], index = 0
     //  if x > A[n-1], index = n-1
     //  this means that if linear interpolation
     //  is used in the caller, this will use the nearest
     //   A[] values to extrapolate off the end of the grid
    size_t nSize = nSizeA;
    double dVal  = dx;

    dVal = std::max(dVal,pdA[0]);
    dVal = std::min(dVal,pdA[nSize-1]);
    
   l = 0;
   u = nSize-1;
   while ( l+1 < u) {
      m = (l+u)/2;
      if ( pdA[m] <= dVal) {
         l = m;
      } else {
         u = m;
      }
   }// end while

   return u;
}

} //namespace ito33
