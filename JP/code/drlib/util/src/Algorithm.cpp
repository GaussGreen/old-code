//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Algorithm.cpp
//
//   Description : Various useful algorithm type functions etc
//
//   Author      : Mark A Robson
//
//   Date        : 10 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Algorithm.hpp"

DRLIB_BEGIN_NAMESPACE
/* sorts array into ascending order */
void Algorithm::shellSort(DoubleArray& doublesToSort){
    shellSort(doublesToSort, 0, doublesToSort.size()-1);
}

/* sorts given portion of List into ascending order */
void Algorithm::shellSort(
    DoubleArray& List,    /* array to sort */
    int          begin,   /* where in the array to start */
    int          end){    /* where in the array to stop  */
    /// may want to replace this with stl version (or have another method)
    double temp;

    int gap;            /* gap between values being compared   */
    int This;           /* this value to be compared....       */
    int next;           /*   ....with this one!                */
    int start;          /* just a counter                      */
    int values;         /* number of values to sort            */

    /* initialise pointers  */
    values = end - begin + 1;
    gap = values;

        /* sort */
    do {
        gap /= 2;
        if (gap > 0) {
            for (start = 0; start < values - gap; start++) {
                This = start;
                while (This >= 0) {
                    next = This + gap;
                    /* primary comparison between 'This' and 'next'   */
                    if (List[This] < List[next]) {
                                /*** swap them  ***/
                        temp = List[This];
                        List[This] = List[next];
                        List[next] = temp;
                    } else {
                        This = 0;
                    }

                    // secondary comparison if primary comparison swapped 
                    This -= gap;
                }
            }
        }
    } while (gap != 0);
}
 
/* sorts given portion of List into ascending order */
void Algorithm::shellSort(
    IDoubleArray& List,    /* array to sort */
    int           begin,   /* where in the array to start */
    int           end){    /* where in the array to stop  */
    /// may want to replace this with stl version (or have another method)
    double temp;

    int gap;            /* gap between values being compared   */
    int This;           /* this value to be compared....       */
    int next;           /*   ....with this one!                */
    int start;          /* just a counter                      */
    int values;         /* number of values to sort            */

    /* initialise pointers  */
    values = end - begin + 1;
    gap = values;

        /* sort */
    do {
        gap /= 2;
        if (gap > 0) {
            for (start = 0; start < values - gap; start++) {
                This = start;
                while (This >= 0) {
                    next = This + gap;
                    /* primary comparison between 'This' and 'next'   */
                    if (List[This] < List[next]) {
                                /*** swap them  ***/
                        temp = List[This];
                        List[This] = List[next];
                        List[next] = temp;
                    } else {
                        This = 0;
                    }

                    // secondary comparison if primary comparison swapped 
                    This -= gap;
                }
            }
        }
    } while (gap != 0);
}

DRLIB_END_NAMESPACE
