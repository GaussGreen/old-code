/*
 * $Log: matlib.h,v $
 * Revision 1.3  2003/08/05 07:57:07  jmprie
 * ajout de GaussLegendre()
 *
 * Revision 1.2  2002/11/25 16:40:49  mab
 * Formatting
 *
 */

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
#ifndef _MATLIB_H
#define _MATLIB_H



#include <stdio.h>
#include <math.h>
#include <string.h>
#include "linalg.h" 
#include "armglob.h"
#ifdef unix
#include <ieeefp.h>
#else
#include <float.h>
#endif



#define EPS 1.0e-9
#define INF exp(3000)


typedef double (*T_FUNC) (ARM_Matrix* , 
                          ARM_Vector*,
                          void **);

extern double* default_op(double *opt=NULL);

extern void searchq(int *pcnt,
                    double *fnew,
                    double oldx,
                    ARM_Vector *matl,
                    ARM_Vector *matx,
                    double gdold,
                    double *stepsize);   
            
extern double cubici1(double fnew,double fold,double graddnew,
                      double graddold,double stepsize);

extern double cubici2(double graddold,
                      ARM_Vector *matl,
                      ARM_Vector *matx);

extern double cubici3(double fnew,double fold,double graddnew,
                      double graddold,double *stepsize);

extern int fminu(ARM_Vector *p,int *CST,
                 void **parameters,
                 ARM_Matrix *data,
                 T_FUNC func,
                 double *OP);

extern void updhess(ARM_Vector xnew,ARM_Vector xold,ARM_Vector gradxnew,
                    ARM_Vector gradxold,ARM_Matrix *invhess,double *para,
                    ARM_Vector *SD);

extern void Hessien(ARM_Vector* x,
                    ARM_Matrix* hessien,
                    void ** parameters,
                    ARM_Matrix* data,
                    T_FUNC func);

extern double leastsq(ARM_Vector *p,int *CST,
                     void **parameters,
                     ARM_Matrix *data,
                     T_FUNC func,
                     double *OP);

extern void GaussLegendre(int n, double* x, double* w);

#endif
/*-------------------------------------------------------------------*/
/*---- End Of File ----*/
