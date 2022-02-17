#ifndef ESL_UTIL_DOT_H
#define ESL_UTIL_DOT_H

/** NOTE: This file should be only included through 'esl_util.c'
 */

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
#include "esl_date.h"
#include "esl_alloc.h"


#ifdef  __cplusplus
extern "C" {
#endif

/*  convert SwapDCC(360, 365, ACT) to CurveDCC(3,5,0) */
int Curve_DCC(const char* swapdcc_str, char* crvdcc);

/*  convert CurveDCC(3,5,0) to SwapDCC(360, 365, ACT) */
int Swap_DCC(char crvdcc, char* swapdcc_str);


/*****  Smooth_Step  ********************************************************/
/**
*       Julia's smooth step function (cf. corresponding memo).
*/
int     Smooth_Step (   double	*SmoothValue,  /**< (O) Output smoothed value */
                        double	UpValue,       /**< (I) Up value              */
                        double	DownValue,     /**< (I) Down value            */
                        double	Index,         /**< (I) Index value           */
                        double	Barrier,       /**< (I) Barrier level         */
                        double	IndexStep);    /**< (I) Index step            */

/*****  DrSmoothStep  ********************************************************/
/**
 *      Smooth version of the step function.
 *      Returns DownValue           if Index <= barrier - step
 *              'smooth curve'      if barrier - step < Index < barrier + step
 *              UpValue             if Index >= barrier + step
 *
 */
double DrSmoothStep(double UpValue,  /**< (I) used if Index > Barrier + step  */
                    double DownValue,/**< (I) used if Index < Barrier - step  */
                    double Index,    /**< (I) index value                     */
                    double Barrier,  /**< (I) barrier value                   */
                    double Step);    /**< (I) smooth between Barrier +/- step */

/*****  DrSmoothMax  ********************************************************/
/**
 *      Smooth version of the max(a,b)
 *      This is the based on the integral of the DrSmoothStep function
 *      Returns b                   if a <= b - step
 *              'smooth curve'      if b - step < a < b + step
 *              a                   if a >= b + step
 *
 */
double  DrSmoothMax(double  a,     /**< (I) Argument 1                 */
                    double  b,     /**< (I) Argument 2                 */
                    double  Step); /**< (I) smooth in [-step, step]    */

/*****  Conv_Freq  **********************************************************/
/**
*       Convert the character Freq ('A'nnual, 'S'emi-Annual, 'Q'uarterly,
*       'M'onthly to an integer
*/
int     Conv_Freq (char   Freq);

/*****  Conv_Index  *********************************************************/
/**
*       Convert the index name into maturity, frequency and base.
*/
int     Conv_Index (int*                    IndexMat,  /**< (O) Index maturity in months */
                    char*                   IndexF,    /**< (O) Frequency as a character */
                    char*                   IndexBase, /**< (O) Index base               */
                    char const*             Index,     /**< (I) Index name               */
                    T_CURVE const*    t_curve);  /**< (I) t_curve data             */


int     Conv_Index_s ( int        *IndexMat,  /* (O) Index maturity in months */
                       char       *IndexF,    /* (O) Frequency as a character */
                       char       *IndexBase, /* (O) Index base               */
                       char const *Index,     /* (I) Index name               */
		       char       SwapF,      /* (I) swap freq convention     */
                       char       SwapDCC,    /* (I) swap dcc convention      */
                       char       LiborConv); /* (I) Money market convention  */

/*****  Coupon_Accrued  ******************************************************/
/**
*       Calculate the accrued interests using the convention and the length
*       of the coupon period. If full coupon period, calculate the coupon.
*/
double  Coupon_Accrued (
		double CouponRate,  /**< Coupon rate as input */
		int    F,           /**< Frequency of payment as an integer */
		char   CouponBase,  /**< Base convention for coupon */
		long   StartDate,   /**< Start of the coupon period */
		long   EndDate,     /**< End of the coupon period */
		long   CurrentDate);/**< Current date */


/*****  Conv_Index_Flex  *****************************************************/
/**
*       Convert the index name into maturity, frequency and base.
*       This version differs from Conv_Index in that any tenor is
*       broken down into <number><period> so 7m rates, etc can be
*		handled.  The function also indicates whether the rate is
*		cash or swap with 'C' or 'S'.
*/
int Conv_Index_Flex(
	 int*                   IndexMat,   /**< (O) Index maturity in months */
     char*                  IndexF,     /**< (O) Frequency as a character */
     char*                  IndexBase,  /**< (O) Index base               */
     char*                  IndexType,  /**< (O) Index type		*/
     char const*            Index,      /**< (I) Index name               */
     T_CURVE const*   t_curve);   /**< (I) t_curve data             */


/*****  NR_Poly  ************************************************************/
/**
*       Newton-Raphson for polynomials.
*/
double  NR_Poly (   double  Guess,  /**< First guess for the root       */
                    double  *a,     /**< Coefficients of the polynomial */
                    int     n);     /**< Degree of the polynomial       */

/*****  Quadratic_Solve  **********************************************/
/**
*       Find analytic roots of quadratic eqn: a[2]x^2+a[1]x+a[0]=0.
*       Uses Numerical Recipes method, p184.
*       Returns failure unless there is at least one real root.
*       If a[2]= 0 then second root is set to be same as first root. 
*/
int  Quadratic_Solve (double  *Root,  /**< (O) Roots of quadratic       */
                      double  *a);    /**< (I) Quadratic coefficients   */


/*****  gaussj  *************************************************************/
/**
*       Solve 1x1 or 2x2 system of linear equations.
*/
int     gaussj (double  Matrix[3][3], /**< (I)    Martix of linear system  */
                int     dim,          /**< (I)    Dimension of the problem */
                double  *vec);        /**< (I/O)  RHS of system / solution */


/*****  polint  *************************************************************/
/**
*       Polynomial interpolation of arbitrary degree n. 
*/
int     polint (double 	xa[],    /**< (I) Set of x values */
                double 	ya[],    /**< (I) Set of y values */
                int     n,       /**< (I) Interpolation order */
                double 	x,       /**< (I) x coordinate of required point */
                double 	*y,      /**< (O) y value */
                double 	*dy);    /**< (O) error estimate */


/*****  dlinterp  ***********************************************************/
/**
*       Linear interpolation between 2 Dates.
*/
int     dlinterp (  long d, 
                    double *y, 
                    long d0, 
                    long d1, 
                    double y0, 
                    double y1);

/*****  linterp  ************************************************************/
/**
*       Linear interpolation between 2 points.
*/
int     linterp (double x, 
                 double *y, 
                 double x0, 
                 double x1, 
                 double	y0, 
                 double y1);

/*****  qinterp  ************************************************************/
/**
*       Quadratic interpolation based on numerical recipes.
*/
void qinterp (double xa[],         /**< (I) Set of x values */
              double ya[],         /**< (I) Set of y values */
              double x,            /**< (I) x coordinate of required point */
              double *y,           /**< (O) y value */
              double DefaultValue);/**< (I) Default value if error is too big */

/*****  d4interp  ***********************************************************/
/**
*       4th degree interpolation between 5 points.
*/
void    d4interp (	double 	x,                  
                    double 	*y, 
                    double 	*xa, 
                    double 	*ya);

/*****  sqinterp  ***********************************************************/
/**
*       Safe quadratic interpolation. The output value is "capped/floored"
*       to within a range proportional to  the difference between the max
*       input ordinate and the min input ordinate.
*/ 
int    sqinterp(double  x0,         /**< (I) Set of three x values          */
                double  x1,
                double  x2,
                double  y0,         /**< (I) Set of three y values          */
                double  y1,
                double  y2,
                double  d0,         /**< (I) Set of three quadratic coefs   */
                double  d1,
                double  d2,
                double  x,          /**< (I) x coordinate of required point */
                double  *y);        /**< (O) required value                 */


/*************** Triangulation ***********************************************/
/**Produces the "usual" Lower triangular matrix for orthogonalising corrolated 
   factors.                                                                    
   NOTE: If Nbfac < 3 then unused matrix elements are set to zero.             
                                                                               
 *****************************************************************************/

int Triangulation (
        double TriangMtx[3][3],/**<(O) Lower triangular matrix transformation*/
        int    Nbfac,          /**<(I) Number of factors                     */
        double *Rho);          /**<(I) correlation coefficients              */

/*****  ExtendSpotVol  *****************************************************/
/**
        Extends a spot vol curve using the flat interp method and modifies
        the date indexing of the extended curve. It copes with cases where
        a period in the extended curve covers multiple periods in the source
        curve.
  
        On input:
        ---------
        SrcSpotVol[t] is the spot vol applicable between ExpDate[t-1] to
        ExpDate[t]. (when t=0, it will be between Today and ExpDate[0])
  
        Note: ExpDates must be in strictly ascending order
  
        On output:
        ----------
        ExtSpotVol[t] is the spot vol applicable between TimePt[t] and 
        TimePt[t+1].
  
        Note: 1  We assume TimePt[0] to be Today
              2  Size of ExtSpotVol must be the same as that of TimePts
              3  We assume ExtSpotVol[NbTimePts-1] is never used because
                 it is meaningless.
              4  TimePts must be in strictly ascending order
  
 */
int  ExtendSpotVol (
        double *ExtSpotVol,    /**< (O) Mem must be alloc'ed on entry */
        long    NbSrcSpotVols, /**< (I) Size of the source curve      */
        long   *ExpDate,       /**< (I) Vol expiry date               */
        double *SrcSpotVol,    /**< (I) Source spot vols              */
        long    NbTimePts,     /**< (I) Size of extended timeline     */
        long   *TimePt);       /**< (I) Extended timeline             */


/*****  tableinterp  ********************************************************/
/**
*       Interpolation based on a array: linear inside, flat outside.
*/
void    tableinterp (double x, 
                     double *y, 
                     double *xVal, 
                     double *yVal, 
                     int    nbVal);
	
/*****  AnnPV  *********************************************************/
/**
*       Evalute value of annuity
*/
double  AnnPV (double    AnnRate, /**< (I) Rate of annuity       */
               double    IRR,     /**< (I) Discounting yield     */
               int       Tenor,   /**< (I) Nb of periods remaing */
               int       Freq);   /**< (I) Frequency             */

int Mapping_2Q(double   *Value,     /* (O)*/
               double   QLeft,      /* (I) */
               double   QRight,     /* (I) */   
               double   FwdShift,
               double   InputX);    /* (I)*/

/**f---------------------------------------------------------------------
 * Vector unary operations.
 *
 * <br><br>
 * Performs a unary vector operation on a vector 
 * <br> (v<sub>0</sub>, ..., v<sub>n-1</sub>).<br>
 * with <i> vectArg</i> as argument (assumed to have same length as <i> v</i>).
 * The argument <i> operation</i> can be one of the follwing: <br><br>
 * <i> "="</i>: copies <i> vectarg</i> to <i> v</i>. <br>
 * <i> "+"</i>: adds <i> vectarg</i> to <i> v</i>. <br>
 * <i> "-"</i>: subtracts <i> vectarg</i> form <i> v</i>. <br>
 * <i> "*"</i>: mutiplies <i> v</i> by <i> vectArg</i> (element by element). <br>
 * <i> "/"</i>: divides <i> v</i> by <i> vectArg</i> (element by element). <br>
 * Returns SUCCESS/FAILURE.
 */

int
VectUnaryOper(
    double *v,        /**< (B) vector [0..n-1]            */
    int n,            /**< (I) number of elements         */
    char *operation,  /**< (I) operation type (=,+,-,*,/) */
    double *vectArg); /**< (I) vector argument [0..n-1]   */


/**f---------------------------------------------------------------------
 * Interpolation 1D spline on doubles: initialization.
 * 
 * <br><br>
 * Initialize a spline interpolation of an array of points
 * $(x^a_i,y^a_i)_{i=0,\dots,m-1}$.
 * The values of ${dy/ dx}$ at $x_0$ and $X_{n-1}$
 * should be given (pass 1e12 for natural splines, i.e. ${d^2y/ dx^2}=0$).
 * On exit, <i> y2</i> contains the coefficients
 * to be used in the function <i> SplineInterp</i>
 * to interpolate.
 * If $x&lt;= x^a_0$ computes $y = y^a_0$
 * (and if $x&gt;= x^a_{m-1}$, $y = y^a_{n-1}$).
 * Returns 0 iff successful.
 */

int
SplineInterp1dInit(
    double *x,      /**< (I) array of X(i) (i=0,..,n-1) */
    double *y,      /**< (I) array of Y(i) (i=0,..,n-1) */
    int n,          /**< (I) # of points */
    double yp1,     /**< (I) dx/dy at X(0) */
    double ypn,     /**< (I) dx/dy at X(n-1) */
    double *y2);    /**< (O) should be allocated on entry */

/*****  AnnIRR  ****************************************************/
/**
*      Find IRR for annuity
*/
int     AnnIRR (double    Outs,            /**< (I) Outstanding            */
                double    AnnRate,         /**< (I) Rate of annuity        */
                int       Tenor,           /**< (I) Nb of periods remaing  */
                int       Freq,            /**< (I) Frequency              */
                double   *IRR);            /**< (O) IRR to be found        */


/*****  DR_ATOF  ****************************************************/
/*
*      Convert a string into double if the string starts with a number
*      otherwise, return error
*/
int DR_ATOF(char    *InputString,    /* (I) */
            double  *OutputDouble);   /* (O) */

#ifdef  __cplusplus
}
#endif

#endif

