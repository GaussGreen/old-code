/*
 * $Log: xdfpmin.h,v $
 * Revision 1.2  2002/06/03 11:45:21  mab
 * ADDED : compilation directive
 *
 */


#ifndef _XDFPMIN
#define _XDFPMIN


/* Driver for routine dfpmin */

#include <stdio.h>
#include <math.h>
#include <string.h>


#include "linalg.h" 


#define EPSM 1.0e-9
#define EPSD 1.0e-8
#define GTOL 1.0e-5
#define ITMAX 500
#define ITMAXDB 50
#define STPMX 200.0
#define ALF 1.0e-9
#define TOLX 1.0e-4
#define TOL 2.0e-4
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define CGOLD 0.3819660
#define ZEPS 1.0e-10


#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);




#define FMAX(a,b) (((a) < (b)) ? (b) : (a))

 
typedef double (*T_FUNC) (ARM_Matrix* , ARM_Vector*, void **);
 
 


#define SQR2(a) ((a)*(a))



#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);



typedef double (*T_FUNC) (ARM_Matrix* , 
                          ARM_Vector*,
                          void **);


extern long dfpmin(ARM_Vector* ,
                   ARM_Matrix* ,
                   int *,
                   double *,
                   void ** ,
                   ARM_Matrix* ,
                   T_FUNC func  
                   );



extern void lnsrch(ARM_Vector* ,
            double ,
            ARM_Vector* ,
            ARM_Vector* ,
            ARM_Vector* ,
            double* ,
            double ,
            int *,
            void ** ,
            ARM_Matrix* ,
            T_FUNC func
            );

double dbrent(double ax, double bx, double cx, double tol, double *xmin,
                ARM_Vector *p,
                ARM_Vector *descente,
                void ** parameters,
                ARM_Matrix* data,
                T_FUNC func);


extern void NumJac(ARM_Vector* ,
            double ,
            ARM_Vector* ,
            void ** ,
            ARM_Matrix* ,
            T_FUNC func );

void NumHess(ARM_Vector* x,
            ARM_Matrix* hessien,
            void ** parameters,
            ARM_Matrix* data,
            T_FUNC func );            

void powell(ARM_Vector* p,
            ARM_Matrix* xi,
            double ftol,
            int *iter,
            double *fret,
            void ** parameters,
            ARM_Matrix* data,
            T_FUNC func );

void linmin(ARM_Vector* p,
            ARM_Vector* xi,
            double *fret, 
            void ** parameters,
            ARM_Matrix* data,
            T_FUNC func );

double f1dim(ARM_Vector* pcom,
                 ARM_Vector* xicom,
                 double x,
                 void** parameters,
                 ARM_Matrix* data,
                 T_FUNC func);

void mnbrak(double* ax, double* bx, double* cx, double* fa,
            double* fb, double* fc,
            void** parameters,
            ARM_Matrix* data,
            T_FUNC func,
            double (*func2)(ARM_Vector*,ARM_Vector*,double,
                            void** ,ARM_Matrix*,T_FUNC func),
                            ARM_Vector* pcom,
                            ARM_Vector* xicom);

double brent(double ax, double bx, double cx,
             void ** parameters, 
             ARM_Matrix* data,   
             T_FUNC func,        
             double (*func2)(ARM_Vector*, ARM_Vector*, double,
                     void** ,ARM_Matrix*, T_FUNC func),
             ARM_Vector* pcom,
             ARM_Vector* xicom,
             double tol, double* xmin);


#endif


/*--------------------------------------------------------------*/
/*---- End Of File ----*/
