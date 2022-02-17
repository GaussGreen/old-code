/*********************************************************************************
 * Transition matrix utils
 *
 ********************************************************************************/


#ifndef _UTIL_H_
#define _UTIL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdarg.h>
#include "transition.h"
#include "srm123.h"

/*********************************************************************************
 * Date routines : convert from YYYYMMDD long format.
 *                              
 * Converts a DR date type (encoded as YYYYMMDD in a long)
 * to a TDate. Returns -1 if failed.
 *
 ********************************************************************************/
TDate DrDate2TDate(long drDate);

/*********************************************************************************
 * Date routines : convert to YYYYMMDD long format.
 *                              
 * Converts a date in TDate type to DR date (encoded as YYYYMMDD in a long)
 *
 ********************************************************************************/
long TDate2DrDate(TDate tDate);

/*********************************************************************************
 * Convert DR DCC to Alib DCC
 *
 ********************************************************************************/
int AlibDCC2DrDCC(
    char *dcc,                            /* (O) Dr DCC                         */
    long AlibDCC);                        /* (I) Alib DCC                       */
    
/*********************************************************************************
 * Convert Dr T_CURVE to Alib TCurve
 *
 ********************************************************************************/
int TCurve2DrCurve(
    T_CURVE  *DrCurve,                    /* (O) Dr curve                       */
    TDate     today,                      /* (I) today                          */
    TCurve   *AlibCurve);                 /* (I) Alib curve                     */

    
/*********************************************************************************
 * n-th eigenvector for Upper-Triangular matrix, given eigenval
 * which is (n,n) element of the matrix
 *
 ********************************************************************************/
int MatNthEigVecUT(
    double        *evec,                  /* (O) ptr to eigenvec                */
    int            idx,                   /* (I) index of the evec to calc      */
    double         eval,                  /* (I) eigenval specified             */
    Mat           *B,                     /* (I) matrix                         */
    double        *pEvec);                /* (I) array to eigenvecs             */

/*********************************************************************************
 * set size for mat
 *
 ********************************************************************************/
int MatSetSize(
    Mat    *A,                            /* (I/O) output matrix                */
    int    size);                         /* (I) size                           */

/*********************************************************************************
 * A = B * C
 *
 ********************************************************************************/
int MatMult(
    Mat    *A,                            /* (O) output matrix                  */
    Mat    *B,                            /* (I) input matrix                   */
    Mat    *C);                           /* (I) input matrix                   */

/*********************************************************************************
 * A = B * C
 * Upper Triangle matrix multiplication only 
 *
 ********************************************************************************/
int MatMultUT(
    Mat    *A,                            /* (O) product                        */
    Mat    *B,                            /* (I) matrix B                       */
    Mat    *C);                           /* (I) matrix C                       */

/*********************************************************************************
 * X = B * Y
 * X: Nx1, B: NxN, Y: Nx1
 *
 ********************************************************************************/
int MatMultVec(
    double     *X,                        /* (O) vector Y                       */
    Mat        *B,                        /* (I) mat                            */
    double     *Y);                       /* (I) vector                         */

/*********************************************************************************
 * X = B * Y
 * X: Nx1, B: NxN, Y: Nx1
 *
 ********************************************************************************/
int MatMultVecUT(
    double    *X,                         /* (O) vector Y                       */
    Mat       *B,                         /* (I) mat                            */
    double    *Y);                        /* (I) vector                         */
    
/*********************************************************************************
 *     Computes exp(t*H), the matrix exponential of a general matrix in
 *     full, using the irreducible rational Pade approximation to the 
 *     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
 *     combined with scaling-and-squaring.
 *     adapted from expokit.
 *
 ********************************************************************************/
int MatExp(
    double       *expH,                   /* (O) output exp(tH)                 */
    int           ideg,                   /* (I) degree of diagonal Pade.       */
    double        t,                      /* (I) time-scale                     */
    double       *H,                      /* (I) matrix                         */
    int           m);                     /* (I) order of H                     */

/*********************************************************************************
 * scale an array
 *
 ********************************************************************************/
int PtrScale(
    double    *ptr,                       /* (I/O) array of doubles             */
    int        size,                      /* (I) size of ptr                    */
    double     scale);                    /* (I) scale factor                   */

/*********************************************************************************
 * copy one ptr to another
 *
 ********************************************************************************/
int PtrCopy(
    double    *to,                        /* (O) destination                    */
    double    *from,                      /* (I) source                         */
    int        size);                     /* (I) size of array                  */

/*********************************************************************************
 * set an arry to a certain value
 *
 ********************************************************************************/
int PtrSet(
    double       *ptr,                    /* (O) array                          */
    int           size,                   /* (I) size of array                  */
    double        val);                   /* (I) value to be set to             */

    
typedef int (*DRDiagnosticMsgCB)(const char*, va_list);
    
/*********************************************************************************
 *       Retrieve the error message in a string buffer
 *
 ********************************************************************************/
const char* CRXErrorRetrieve();

/*********************************************************************************
 *   Print an error message.
 *
 ********************************************************************************/
void CRXError(const char* format, ...);

/*********************************************************************************
 *       Store the error message in a string buffer
 *
 ********************************************************************************/
int CRXErrorCallback(const char* format, va_list ap);

/*********************************************************************************
 * set the error callback function as porposed by Jerry.
 *
 ********************************************************************************/
void CRXErrorSetCallBack(DRDiagnosticMsgCB errorCallbackFnp);

/*f----------------------------------------------------------------------------
 * Double precision cumulative normal function
 */
double NormCum (double x);

/*f----------------------------------------------------------------------------
 * Inverse cumulative normal distribution
 *
 * Based on Risk Magazine
 */
double NormCumInv (double prob);

/*********************************************************************************
 * Returns the value of the beta function B(z, w)
 *
 ********************************************************************************/
double Beta(double z, double w);
    
/*********************************************************************************
 * Returns the value of the beta dist function B(z, w) pdf
 *
 ********************************************************************************/
double BetaPDF(double x, double z, double w);

/*********************************************************************************
 * Returns the incom lete beta function I x (a, b).
 *
 ********************************************************************************/
double Betai(double x, double a, double b);
    
    
#ifdef __cplusplus
}
#endif

#endif
