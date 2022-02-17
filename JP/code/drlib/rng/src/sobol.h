
#ifndef _SC_SOBOL_H
#define _SC_SOBOL_H

#include "edginc/coreConfig.hpp"
CORE_BEGIN_NAMESPACE


//#include <stdio.h>
//#include <stdlib.h>

enum { DIMENSION = 396};

/******************************************************
 This generator can generate at most 2^30 quasi
 random numbers in each dimension.
 To get more, just simply modify the definition of
 MAXBIT
******************************************************/
enum { MAXBIT=30};

/******************************************************
 The MAXORDPOLY specifies the maximum order of the
 primitive polynomial.
 To add more polynomials, MAXORDPOLY, the order of
 each polynomial deg[], and initdat need to be changed
*******************************************************/
enum { MAXORDPOLY=12};

int sobinit( int dimension,
             int **Sobol_gen,
             int *Sobol_xn,
             int *Sobol_n,
             int *Sobol_nstream /*pass by reference */
           );
int sobvect( int dimension,
             double *x,
             int **Sobol_gen,
             int *Sobol_xn,
             int *Sobol_n,
             int *Sobol_nstream /*pass by reference */
           );
int sobstream( int **Sobol_gen,
               int *Sobol_xn,
               int *Sobol_n,
               int *Sobol_nstream /*pass by reference */
             );
double sobolseq(int streamid,
                int **Sobol_gen,
                int *Sobol_xn,
                int *Sobol_n,
                int *Sobol_nstream /*pass by reference */
               );
CORE_END_NAMESPACE
#endif // _SC_SOBOL_H
