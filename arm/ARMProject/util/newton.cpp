/* 
 * $Log: newton.cpp,v $
 * Revision 1.4  2003/07/18 06:33:46  ebenhamou
 * corrected some code... but can tell this is not used
 *
 * Revision 1.3  2001/12/05 14:49:39  arm
 * Rajout de : comments $Log: newton.cpp,v $
 * Rajout de : comments Revision 1.4  2003/07/18 06:33:46  ebenhamou
 * Rajout de : comments corrected some code... but can tell this is not used
 * Rajout de : comments
 *
 */


/*----------------------------------------------------------------------------*
     newton.cpp
 
    Newton-Raphson routines to find roots of functions

*----------------------------------------------------------------------------*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "expt.h"
#include "armglob.h"

#include "newton.h"
#include "linalg.h"






/*----------------------------------------------------------------------------*
    SYNOPSIS
        double newtonRoot(
            double (*function)(double &, double &, double, void **), 
            double x1, double x2, 
            void     **fixedParams, 
            double     rootTol, 
            int        iterMax)


    Uses the Newton-Raphson method to find a root of a function.

    This is a relatively failsafe Newton-Raphson routine.  It
    requires that the user bracket the desired root in a interval
    between 'x1' and 'x2' thus insuring that at least something
    is known about the function.  After first checking to make
    sure that a root is bracketed (i.e. f(x1) * f(x2) < 0), it
    starts a loop that takes a Newton-Raphson step each iteration
    unless the step would take it outside the bracketing boundaries.
    In this case a bisection step is taken.  After each step the
    new value becomes a new boundary.  Thus, in the worst case
    this routine becomes a simple bisection search.

    Input:  
        double (*function)() - a pointer to a routine that calculates the
              value of the function and its derivative. Synopsis must be 
        double function(double &f, double &dfdx, double x, void **fixedParams)
         with x the variable, f the function value, dfdx the function 
       derivative and fixedParams any fixed parameters required for function
                                evaluation. 
      double    x1, x2          -  x1 and x2 must bracket a root
      void    **fixedParams      - an array that carries any other 
           fixed parameters that the function and its derivative 
              may depend on. 
      double    rootTol            - the error tolerance for the root
      int        iterMax            - the maximum number of iterations allowed.

    Output:
        The routine returns K_HUGE_DOUBLE in any of the following situations:

            1)   f(x1) or f(x2) >= K_HUGE_DOUBLE.

            2)   'x1' and 'x2' do not bracket a root, i.e. f(x1) * f(x2) > 0

            3)   any intermediate value of f(x) or f'(x) >= K_HUGE_DOUBLE.

            4)   the maximum number of iterations in the loop has been 
                 exceeded.  

         Otherwise it returns the value of the root having satisfied
         the condition;

                1)   abs((x(n) - x(n-1))) < 'rootTol' + 4*K_DOUBLE_TOL*x(n)

         If the interval contains an odd number of roots, only one will
         be found.  If it contains an even number then of course it will
         not satisfy the bracketing condition.

*----------------------------------------------------------------------------*/

double newtonRoot(double (*function)(double &, double &, double, void **), 
                  double     x1,
                  double     x2, 
                  void*      fixedParams[], 
                  double     rootTol, 
                  int        iterMax)
{
    double swap, xl, xh, fl, fh, approxfl, approxfh,
           xnew, dxold, dx, f, dfdx;

    int    count=0;


    // Check to make sure that the root is bracketed.  If one
    //  of the end points is close enough, return that value.
    //  Also make sure that the function has not returned values
    //  too large to consider.

    fl = (*function)(fl, dfdx, x1, fixedParams);
    fh = (*function)(fh, dfdx, x2, fixedParams);

    if (( fabs(fl) >= K_HUGE_DOUBLE ) 
        ||
        ( fabs(fh) >= K_HUGE_DOUBLE )
       )
    {
       throw MathException(__LINE__, __FILE__, ERR_FUNCTION_EVAL_OVFL,
            "Function value out of range");
    }
    else if ( fabs(fl) <= 2.0 * K_DOUBLE_TOL )
       return(x1);
    else if ( fabs(fh) <= 2.0 * K_DOUBLE_TOL )
       return(x2);
    else if ((fl/fabs(fl))*(fh/fabs(fh)) > 0.0 )
    {
       throw MathException(__LINE__, __FILE__, ERR_INITIAL_VALUE_PB,
            "Root is not bracketed");
    }

     //    Orient the search so that fl < 0 < fh.  This allows the
     //    bisection step to work if the newton-raphson fails.

    if ( fl < 0.0 )
    {
       xl = x1;
       xh = x2;
    }
    else
    {
       xh = x1;
       xl = x2;
       swap = fl;
       fl = fh;
       fh = swap;
    }

    //    Initialize the parameters.

    xnew = 0.5 * (xh + xl);
    dxold = fabs(xh - xl);
    dx = dxold;

    f = (*function)(f, dfdx, xnew, fixedParams);

    if (( f >= K_HUGE_DOUBLE ) || ( dfdx >= K_HUGE_DOUBLE ))
    {
       throw MathException(__LINE__, __FILE__, ERR_FUNCTION_EVAL_OVFL,
            "Function value out of range");
    }

    count = 0;

    while((fabs(dx) > (rootTol + 4.0*K_DOUBLE_TOL*xnew))
          && (f != 0.0)  && (count++ < iterMax)) 
    {

        approxfl = f - (xnew-xl)*dfdx;
        approxfh = f - (xnew-xh)*dfdx;

        //  If    1) the sum of the logs of the approximate function values
        //     'approxfl' and 'approxfh' is greater than K_LOG_HUGE_DOUBLE
        //
        //  or    2) the newton-raphson would try to step outside the
        //           bracketing boundaries,
        //
        //  or    3) the step it would take is smaller than that the
        //           bisection step would take,
        //
        //  then take a bisection step.

        if ((log(fabs(approxfl)) + log(fabs(approxfh))) > K_LOG_HUGE_DOUBLE ||
               approxfl * approxfh >= 0.0 || fabs(2.0*f) > fabs(dxold*dfdx)) 
        {
            dxold = dx;
            dx = 0.5 * (xh-xl);
            xnew = xl + dx;
        }
        else 
        {    // newton-raphson is ok
            dxold = dx;
            dx = f/dfdx;
            xnew  -= dx;
        }

        f = (*function)(f, dfdx, xnew, fixedParams);

        if (f >= K_HUGE_DOUBLE || dfdx >= K_HUGE_DOUBLE)
        {
            throw MathException(__LINE__, __FILE__, ERR_FUNCTION_EVAL_OVFL,
                "Function value out of range");
        }


        //  Reset bracketing boundaries for bisection step test.

        if ( f < 0.0 )
        {
           xl = xnew;
           fl = f;
        }
        else
        {
           xh = xnew;
           fh = f;
        }
    }

    //  Check to see if number of iterations has been exceeded.
    //  Note that this is almost unecessary since at a minimum
    //  iterMax bisection steps would be taken if the error tolerance
    //  has not been met.

    if ( count >= iterMax )
    {
       throw MathException(__LINE__, __FILE__, ERR_MAX_ITER_NUM_EXD,
            "Maximum number of iterations reached");
    }
    
    return(xnew);

}




/*----------------------------------------------------------------------------*
    SYNOPSIS
        int  multiDimNewtonRoot(
        int    (*function)(ARM_Vector *, ARM_Matrix *, ARM_Vector *, void **),
            ARM_Vector     &root,
            void        **fixedParams, 
            double        tol, 
            int            itermax)
        
    multi dim Newton-Raphson algorithm to find root of 
            function of several variables. 
    
    Input:
        ARM_Vector    *(*function)(ARM_Vector *, ARM_Matrix *, 
                         ARM_Vector *, void **),
                    - user input function. Prototype should be
                    ARM_Vector    *f(ARM_Vector *val,
                    ARM_Matrix *grad, ARM_Vector *x, void **fixedParams)
                            with x the variables, 
                    val the function value, grad the function gradient
                and fixedParams any fixed parameters necessary to
                                    function evaluation.
        ARM_Vector    *root            - first guess for root
        void    **fixedParams        - fixed parameters for function evaluation
        double tol                    - tolerance to use
        int    itermax                    - maximum number of iterations
    Output:
        ARM_Vector    *root            - found root
        The routine returns a pointer to root if ok, NULL otherwise

*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
	EB: 
	
	WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

	NOT EVER TESTED AND USED 
	since I had to add the initialisation code for the Gradient
	fGrad= new ARM_Matrix(size,size);

	SHOULD PROBABLY BE CLEANED OUT

*----------------------------------------------------------------------------*/


int multiDimNewtonRoot(
    int    (*function)(ARM_Vector *, ARM_Matrix *, ARM_Vector *, void **),
    ARM_Vector     *root,
    void        **fixedParams, 
    double        tol, 
    int            itermax)
{
    int k, i, size;
    double errx, errf, det;
    ARM_Vector    *fVal;
    ARM_Matrix    *fGrad;



    //    allocate memory
    size = root->GetSize();
    fVal = new ARM_Vector(size);
	fGrad= new ARM_Matrix(size,size);
   

    //    Iterate to solve for zeros of function    
    
    for (k=1; k<=itermax; k++) 
    {
        //    evaluate function and gradient    
        
        if (!function(fVal, fGrad, root, fixedParams))
        {
            throw MathException(__LINE__, __FILE__, ERR_FUNCTION_EVAL_OVFL,
            "Function evaluation failed");
        }
        
        
        for (i=0; i<size; i++) 
        {
            (*fVal)[i] =  -(*fVal)[i];
        }

        //    check for convergence    
        errf=0.0;
        for (i=0; i<size; i++) 
        {
            errf += fabs((*fVal)[i]);
        }

        if (errf <= tol) break;

        //    solve linear system    
        fGrad->LinSolve(fVal, det);

        //    update root    
        errx=0.0;

        for (i=0; i<size; i++) 
        {
            errx += fabs((*fVal)[i]);
            (*root)[i] += (*fVal)[i];
        }
    
        //    check for convergence    
        if (errx <= tol) break;
    }

    
    if (k==itermax)
    {
        throw MathException(__LINE__, __FILE__, ERR_MAX_ITER_NUM_EXD,
            "Maximum number of iterations reached");
    }
    
    
    return(1);
}

 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/    
