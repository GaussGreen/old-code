/* ======================================================
   FILENAME:  num_h_gaussj.h

   PURPOSE:   inversion de matrice par methodes de GAUSS
   ====================================================== */

#ifndef NUM_H_GAUSSJ_H
#define NUM_H_GAUSSJ_H

#include "uterror.h"

/*
        Solves by Gauss-Jordan elimination for Ax1 = B1      , Ax2 = B2 ,... On
   output      , **a is replaced by its inverse      , and **b (all the vectors)
   is replaced by the solution vectors x1      , x2      ,...
*/

Err gaussj(double **a, int n, double **b, int m);

#endif