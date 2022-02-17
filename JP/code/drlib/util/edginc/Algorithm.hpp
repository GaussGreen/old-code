//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Algorithm.hpp
//
//   Description : Various useful algorithm functions etc
//
//   Author      : Mark A Robson
//
//   Date        : 10 Oct 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ALGORITHM_HPP
#define EDG_ALGORITHM_HPP
#include "edginc/AtomicArray.hpp"
#include "edginc/IDoubleArray.hpp"

DRLIB_BEGIN_NAMESPACE

/** Various useful algorithm functions etc */

class UTIL_DLL Algorithm {

public:
    /* sorts array into ascending order */
    static void shellSort(DoubleArray& doublesToSort);

    /* sorts given portion of array into ascending order */
    static void shellSort(
        DoubleArray& doublesToSort,   /* array to sort */
        int          begin,   /* where in the array to start */
        int          end);    /* where in the array to stop  */
    
    static void shellSort(
        IDoubleArray& doublesToSort,   /* array to sort */
        int           begin,   /* where in the array to start */
        int           end);    /* where in the array to stop  */
};

DRLIB_END_NAMESPACE

#endif

