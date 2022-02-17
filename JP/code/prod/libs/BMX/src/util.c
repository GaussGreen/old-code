/****************************************************************************/
/*      Utility routines.                                                   */
/****************************************************************************/
/*      UTIL.c                                                              */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>  
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "bmx123head.h"



/*****  GetIndexStep  *******************************************************/
/*                                                                           
 *      Obtains the difference between values of an index at adjacent nodes.
 *      This value is then used in the smoothing algorithm (Smooth_Step).
 *      This function hides the dimension generality from  the  caller, 
 *      but it has no means to check that the adequate amount of space has 
 *      been allocated under the void * being passed.
 */
double   GetIndexStep ( double     *Index,      /* (I) Index pointer       */
                        int         Dim,        /* (I) Index dimension     */
                        int         i,          /* (I) Node indices        */
                        int         j,
                        int         k,
                        int         t,          /* (I) Current time point  */
                        TREE_DATA   *tree_data) /* (I) Tree data structure */
{

    double  *IndexL;            /* Local pointer */

    double  IndexStep;          /* Output index step */
    double  IndexVal;           /* Index value at mid node */

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */



    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    IndexStep = ERROR;  /* To avoid division by 0 */
                
    switch (Dim)
    {
        case 1:
        {
            IndexL = Index + Node_Offset(1, 0, 0, t, tree_data);

            if (i > Bottom1)
                IndexStep = MAX (IndexStep, fabs (IndexL[i-1] - IndexL[i]));
            if (i < Top1)                                                           
                IndexStep = MAX (IndexStep, fabs (IndexL[i+1] - IndexL[i]));

            break;
        }                    
        case 2:
        {
            IndexL = Index + Node_Offset(2, i, 0, t, tree_data);

            IndexVal = IndexL[j];

            if (j > Bottom2[i])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[j-1] - IndexVal));
            }
            if (j < Top2[i])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[j+1] - IndexVal));
            }

            if (i > Bottom1)
            {
                IndexL = Index + Node_Offset(2, i-1, 0, t, tree_data);

                if ((j>=Bottom2[i-1])&&(j<=Top2[i-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }
            if (i < Top1)
            {         
                IndexL = Index + Node_Offset(2, i+1, 0, t, tree_data);

                if ((j>=Bottom2[i+1])&&(j<=Top2[i+1]))
                {                                                  
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }

            break;
        }                    
        case 3:
        {
            IndexL = Index + Node_Offset(3, i, j, t, tree_data);

            IndexVal = IndexL[k];

            if (k > Bottom3[i][j])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[k-1] - IndexVal));
            }
            if (k < Top3[i][j])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[k+1] - IndexVal));
            }

            if (i > Bottom1)
            {
                if ((j>=Bottom2[i-1])&&(j<=Top2[i-1]))
                {
                    IndexL = Index + Node_Offset(3, i-1, j, t, tree_data);

                    if ((k>=Bottom3[i-1][j])&&(k<=Top3[i-1][j]))
                    {
                        IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                    }
                }
            }
            if (i < Top1)    
            {
                if ((j>=Bottom2[i+1])&&(j<=Top2[i+1]))
                {
                    IndexL = Index + Node_Offset(3, i+1, j, t, tree_data);

                    if ((k>=Bottom3[i+1][j])&&(k<=Top3[i+1][j]))
                    {                                                       
                        IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                    }
                }
            }
            if (j > Bottom2[i])
            {
                IndexL = Index + Node_Offset(3, i, j-1, t, tree_data);

                if ((k>=Bottom3[i][j-1])&&(k<=Top3[i][j-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                }
            }
            if (j < Top2[i])
            {
                IndexL = Index + Node_Offset(3, i, j+1, t, tree_data);

                if ((k>=Bottom3[i][j+1])&&(k<=Top3[i][j+1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                }
            }
        }                    
        default:
        {
            break;
        }                    
    }  /* switch */


    return (IndexStep);

}  /* GetIndexStep */


 
/*****  Smooth_Step  ********************************************************/
/*
*       Julia's smooth step function (cf. corresponding memo).
*/
int     Smooth_Step (   double  *SmoothValue,  /* (O) Output smoothed value */
                        double  UpValue,       /* (I) Up value              */
                        double  DownValue,     /* (I) Down value            */
                        double  Index,         /* (I) Index value           */
                        double  Barrier,       /* (I) Barrier level         */
                        double  IndexStep)     /* (I) Index step            */
{

    double  x, y;
                     
                     
    x = (Index - Barrier) / IndexStep;
    
    if (x < -1.)                    /* Step function = 0 below step */
    { 
        *SmoothValue = DownValue;
    }
    else if (x > 1.)                /* Step function = 1 above step */
    { 
        *SmoothValue = UpValue;
    }
    else                            /* Step function: polynomial in between */
    {   
        y = 0.5 + x / 16. * (15. + x * x * (-10. + 3. * x * x));
        
        *SmoothValue = (UpValue - DownValue) * y + DownValue;
        
    }  /* if then else */   


    return (SUCCESS);

}  /* Smooth_Step */


/*****  DrSmoothStep  ********************************************************/
/*
 *      Smooth version of the step function.
 *      Returns DownValue           if Index <= barrier - step
 *              'smooth curve'      if barrier - step < Index < barrier + step
 *              UpValue             if Index >= barrier + step
 *
 */
double DrSmoothStep(double  UpValue,  /* (I) used if Index > Barrier + step  */
                    double  DownValue,/* (I) used if Index < Barrier - step  */
                    double  Index,    /* (I) index value                     */
                    double  Barrier,  /* (I) barrier value                   */
                    double  Step)     /* (I) smooth between Barrier +/- step */
{
    double smoothValue, xSquare, x;

    Step = ABS(Step);

    /* do the limiting case when step = 0.0 */
    if (Step < TINY)
    {
        smoothValue = (Index < Barrier) ? DownValue : UpValue;
        return (smoothValue);
    }

    /* do the smooth case */
    x = (Index - Barrier)/Step;

    if (x <= -1.0)
    {
        smoothValue = DownValue;
    }
    else if (x >= 1.0)
    {
        smoothValue = UpValue;
    }
    else
    {
        xSquare = x * x;
        smoothValue = ((3.*xSquare-10.) * xSquare + 15.) * x *.0625 +.5;
        smoothValue = DownValue + (UpValue - DownValue) * smoothValue;
    } 
    return (smoothValue);

} /* DrSmoothStep */


/*****  DrSmoothMax  ********************************************************/
/*
 *      Smooth version of the max(a,b)
 *      This is the based on the integral of the DrSmoothStep function
 *      Returns b                   if a <= b - step
 *              'smooth curve'      if b - step < a < b + step
 *              a                   if a >= b + step
 *
 */
double  DrSmoothMax(double  a,     /* (I) Argument 1                 */
                    double  b,     /* (I) Argument 2                 */
                    double  Step)  /* (I) smooth in [-step, step]    */
{
    double smoothValue, x, xSquare;

    Step = ABS(Step);

    /* do the limiting case when step = 0.0 */
    if (Step < TINY)
    {
        smoothValue = (a < b) ? b : a;
        return (smoothValue);
    }

    /* do the smooth case */
    x = (a - b)/Step;

    if (x <= -1.0)
    {
        smoothValue = b;
    }
    else if (x >= 1.0)
    {
        smoothValue = a;
    }
    else
    {
        xSquare = x * x;
        smoothValue = ((xSquare-5.) * xSquare + 15.) * xSquare * 0.03125 + 
                      0.5 * x + 0.15625;
        smoothValue = Step * smoothValue + b;
    } 
    return (smoothValue);
}


/*****  Conv_Freq  **********************************************************/
/*
*       Convert the character Freq ('A'nnual, 'S'emi-Annual, 'Q'uarterly,
*       'M'onthly to an integer
*/
int     Conv_Freq (char   Freq)
{
    switch (Freq)
    {
        case 'A':
            return (1);
        case 'S':
            return (2);
        case 'Q':
            return (4);
        case 'M':
            return (12);
    }

    return (0);

}  /* Conv_Freq */



/*****  Conv_Index  *********************************************************/
/*
*       Convert the index name into maturity, frequency and base.
*/
int     Conv_Index ( int         *IndexMat,  /* (O) Index maturity in months */
                     char        *IndexF,    /* (O) Frequency as a character */
                     char        *IndexBase, /* (O) Index base               */
                     char        *Index,     /* (I) Index name               */
                     T_CURVE     *t_curve)   /* (I) t_curve data             */
{

    char    SwapF;              /* Underlying benchmark swap frequency */
    char    LiborDCC;          /* Base of yield curve Libor           */
    char    SwapDCC;           /* Base of yield curve benchmark swap  */

    int     status = FAILURE;   /* Error status = FAILURE initially    */

        
    if (t_curve->MMB == 360)
        LiborDCC = '0';
    else        
        LiborDCC = '5';
            
    if (!strcmp (t_curve->SwapDCC, "ACT"))
        SwapDCC = '3';
    else if (!strcmp (t_curve->SwapDCC, "360"))
        SwapDCC = '0';
    else 
        SwapDCC = '5';
        
    SwapF = t_curve->SwapFreq;
        
    if (     !strcmp (Index, "1m"))
    {
        *IndexMat  = 1;
        *IndexF    = 'M';
        *IndexBase = LiborDCC;         
    }   
    else if (!strcmp (Index, "3m"))
    {
        *IndexMat  = 3;
        *IndexF    = 'Q';
        *IndexBase = LiborDCC;         
    }   
    else if (!strcmp (Index, "6m"))
    {
        *IndexMat  = 6;
        *IndexF    = 'S';
        *IndexBase = LiborDCC;         
    }   
    else if (!strcmp (Index, "12m"))
    {
        *IndexMat  = 12;
        *IndexF    = 'A';
        *IndexBase = LiborDCC;         
    }   
    else if (!strcmp (Index, "1y"))
    {
        *IndexMat  = 12;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "2y"))
    {
        *IndexMat  = 24;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "3y"))
    {
        *IndexMat  = 36;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "4y"))
    {
        *IndexMat  = 48;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "5y"))
    {
        *IndexMat  = 60;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "6y"))
    {
        *IndexMat  = 72;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "7y"))
    {
        *IndexMat  = 84;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "8y"))
    {
        *IndexMat  = 96;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "9y"))
    {
        *IndexMat  = 108;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "10y"))
    {
        *IndexMat  = 120;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "12y"))
    {
        *IndexMat  = 144;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "15y"))
    {
        *IndexMat  = 180;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "20y"))
    {
        *IndexMat  = 240;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "25y"))
    {
        *IndexMat  = 300;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else if (!strcmp (Index, "30y"))
    {
        *IndexMat  = 360;
        *IndexF    = SwapF;
        *IndexBase = SwapDCC;         
    }   
    else    
    {
        DR_Error ("Conv_Index3: incorrect index name!");
        goto RETURN;
        
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    return (status);


}  /* Conv_Index */



/*****  Conv_Index_Flex  *********************************************************/
/*
*       Convert the index name into maturity, frequency and base.
*       This version differs from Conv_Index in that any tenor is
*       broken down into <number><period> so 7m rates, etc can be
*       handled.  The function also indicates whether the rate is
*       cash or swap with 'C' or 'S'.
*/
int Conv_Index_Flex(
    int                 *IndexMat,          /* (O) Index maturity in months */
    char                *IndexF,            /* (O) Frequency as a character */
    char                *IndexBase,         /* (O) Index base               */
    char                *IndexType,         /* (O) Index type               */
    char                *Index,             /* (I) Index name               */
    T_CURVE             *t_curve)           /* (I) t_curve data             */
{
    char                SwapF;              /* Underlying benchmark swap frequency */
    char                LiborDCC;           /* Base of yield curve Libor           */
    char                SwapDCC;            /* Base of yield curve benchmark swap  */
    int                 status = FAILURE;   /* Error status = FAILURE initially    */

    /* Libor daycount convention */
    if (t_curve->MMB == 360)
        LiborDCC = '0';
    else        
        LiborDCC = '5';
            
    /* Swap daycount convention */
    if (!strcmp (t_curve->SwapDCC, "ACT"))
        SwapDCC = '3';
    else if (!strcmp (t_curve->SwapDCC, "360"))
        SwapDCC = '0';
    else 
        SwapDCC = '5';

    /* Swap payment frequency */
    SwapF = t_curve->SwapFreq;

    /* Type code */
    switch (tolower(Index[strlen(Index)-1]))
    {
    case 'm':   /* Months (cash) */
                *IndexType = CASH_RATE;
                *IndexBase = LiborDCC;         
                *IndexMat  = atoi(Index);
                *IndexF    = 'A';
                break;

    case 'y':   /* Years (swap) */
                *IndexType = SWAP_RATE;
                *IndexF    = SwapF;
                *IndexBase = SwapDCC;         
                *IndexMat  = atoi(Index) * 12;
                break;

    default:    /* Unknown index */
                DR_Error ((char*)"Conv_Index_Flex: incorrect index name!");
                goto RETURN;
    }

    /* OK if we get here */
    status = SUCCESS;

RETURN:

    /* Exit */
    return status;
} 


/*****  NR_Poly  ************************************************************/
/*
*       Newton-Raphson for polynomials.
*/
double  NR_Poly (   double  Guess,  /* First guess for the root       */
                    double  *a,     /* Coefficients of the polynomial */
                    int     n)      /* Degree of the polynomial       */
{

    double  x;    /* Working value of the root */
    double  dx;   /* Variation of x from one iteration to the other */
    double  P;    /* Value of the polynomial */
    double  dP;   /* Value of the derivative of the polynomial */

    int     i;
    double  j;    /* Iteration index */


    x  = Guess;

    for (j = 0; j < 20; j++)   /* Stop after 20 iterations */
    {
        P  = a[n];
        dP = 0.;

        /* Evaluate P and dP at the current value x */
        for (i = n - 1; i >= 0; i--)
        {
            dP = P + dP * x;
            P  = a[i] + P * x;
        }

        if (dP == 0.)
            DR_Error ("Derivative should not vanish. (NR_Poly)");

        dx = P/dP;
        x -= dx;    /* Next guess in the Newton-Raphson */

        if (fabs (dx) < TINY)
            return (x);

    }  /* for j */

    DR_Error ("NR_Poly: maximum number of iterations exceeded");

    return (0.);

}  /* NR_Poly */



/*****  Quadratic_Solve  **********************************************/
/*
*       Find analytic roots of quadratic eqn: a[2]x^2+a[1]x+a[0]=0.
*       Uses Numerical Recipes method, p184.
*       Returns failure unless there is at least one real root.
*       If a[2]= 0 then second root is set to be same as first root. 
*/
int  Quadratic_Solve (double  *Root,  /* (O) Roots of quadratic       */
                      double  *a)     /* (I) Quadratic coefficients   */
{
    double  a0;
    double  a1;
    double  a2;
    int     status = FAILURE;   /* Error status */

    if (Root == NULL || a == NULL) goto RETURN;

    a0 = a[0];
    a1 = a[1];
    a2 = a[2];

    if (fabs(a2) < TINY)
    {
        if (fabs(a1) < TINY) goto RETURN; /* ill-defined quadratic */
        else                              /* equation is linear    */
        {
            Root[0] = Root[1] = -a0/a1;
        }
    }
    else
    {
        double Discriminant;
        Discriminant = a1*a1 - 4*a0*a2; 
        
        if ( (SIGN(Discriminant) < 0) && (sqrt(fabs(Discriminant)) >= TINY) )
        {
            goto RETURN;                          /* no real roots    */
        }
        else if (sqrt(fabs(Discriminant)) < TINY) /* single real root */
        {
            Root[0] = Root[1] = - 0.5 * a1 / a2;
        }
        else                                      /* two real roots   */
        {
            double  q;

            q = -0.5*(sqrt(Discriminant) * SIGN(a1) + a1);
            
            Root[0] = q / a2;
            Root[1] = a0 / q;
        }
    }

    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
    {
        DR_Error ("Quadratic_Solve: at least one real root does not exist");
    }
    
    return (status);

}  /* Quadratic_Solve */


/*****  gaussj  *************************************************************/
/*
*       Solve 1x1 or 2x2 system of linear equations.
*/
int     gaussj (double  Matrix[3][3], /* (I)    Martix of linear system  */
                int     dim,          /* (I)    Dimension of the problem */
                double  *vector)      /* (I/O)  RHS of system / solution */
{

    double  det;                /* Determinant               */
    double  maxelt;             /* Maximal element in matrix */
    double  v0;

    int     status = FAILURE;   /* Error status  */
    char    ErrorMsg[MAXBUFF];  /* Error message */


    if (dim == 1)        /* 1x1 system */
    {
        if (fabs(Matrix[0][0]) < ERROR)
        {
            sprintf (ErrorMsg, "gaussj: linear system is unstable!");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        vector[0] /= Matrix[0][0];
    }
    else if (dim == 2)
    {
        det = Matrix[1][1] * Matrix[0][0] - Matrix[0][1] * Matrix[1][0];
        maxelt = MAX(fabs(Matrix[0][0]),fabs(Matrix[0][1]));
        maxelt = MAX(maxelt,fabs(Matrix[1][0]));
        maxelt = MAX(maxelt,fabs(Matrix[1][1]));

        if (fabs(det)/maxelt < ERROR)
        {
            sprintf (ErrorMsg, "gaussj: linear system is unstable!");
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        v0 = vector[0];
        vector[0] = (vector[0] * Matrix[1][1] - vector[1] * Matrix[0][1])/det;
        vector[1] = (vector[1] * Matrix[0][0] - v0        * Matrix[1][0])/det;
    }
    else 
    {
        sprintf (ErrorMsg, "gaussj: unsupported dimension of linear system!");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* gaussj */



/*****  polint  *************************************************************/
/*
*       Polynomial interpolation of arbitrary degree n. 
*/
int     polint (double  xa[],    /* (I) Set of x values */
                double  ya[],    /* (I) Set of y values */
                int     n,       /* (I) Interpolation order */
                double  x,       /* (I) x coordinate of required point */
                double  *y,      /* (O) y value */
                double  *dy)     /* (O) error estimate */
{

    double  *c,*d;

    int     i,m,ns=1;
    double  den,dif,dift,ho,hp,w;

    int     status = FAILURE;   /* Error status  */



    c = (double *) DR_Array (DOUBLE, 1, n);
    d = (double *) DR_Array (DOUBLE, 1, n);

    if ((c == NULL) || (d == NULL))
    {
        DR_Error ("polint: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;        
    }

    dif=fabs(x-xa[1]);
    for (i=1;i<=n;i++) {
        if ( (dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=1;i<=n-m;i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ((den=ho-hp) == 0.0)
            {
                DR_Error ("Error in routine polint");
                goto FREE_MEM_AND_RETURN;
            }
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Free_DR_Array (c, DOUBLE, 1, n);
    Free_DR_Array (d, DOUBLE, 1, n);

    return (status);

}  /* polint */



/*****  dlinterp  ***********************************************************/
/*
*       Linear interpolation between 2 Dates.
*/
int     dlinterp (  long d, 
                    double *y, 
                    long d0, 
                    long d1, 
                    double y0, 
                    double y1)
{
    int    status = FAILURE;
    double a,b,l0,l1;

    /* basic check */
    if (d0 == d1)
    {
        DR_Error("dlinterp: points supplied for "
                 "interpolation fall on the same date.\n");
        goto RETURN;
    }

    if ((d == d0) || IS_EQUAL(y0,y1))
    {
        *y = y0;
    }
    else if (d == d1)
    {
        *y = y1;
    }
    else
    {
        a=(double) Daysact(d ,d1);
        b=(double) Daysact(d0,d1);
        l0=a/b;

        a=(double) Daysact(d0,d );
        l1=a/b;

        *y=l0*y0+l1*y1;
    }

    status = SUCCESS;

RETURN:

    return(status);

}  /* dlinterp */



/*****  linterp  ************************************************************/
/*
*       Linear interpolation between 2 points.
*/
int     linterp (double x, 
                 double *y, 
                 double x0, 
                 double x1, 
                 double y0, 
                 double y1)
{
    int    status = FAILURE;
    double a,b,l0,l1;

    /* basic check */    
    if ((x == x0) || (y0 == y1))
    {
        *y = y0;
    }
    else if (x == x1)
    {
        *y = y1;
    }
    else
    {
        a=x-x1;
        b=x0-x1;
        l0=a/b;

        a=x-x0;
        b=x1-x0;
        l1=a/b;

        *y=l0*y0+l1*y1;
    }

    status = SUCCESS;

RETURN:

    return(status);

}  /* linterp */



/*****  linterp2d  ************************************************************/
/*
*       Linear interpolation between 4 points in 2d.
*/
int     linterp2d (double x, 
                   double y,
                   double *z,
                   double x0, 
                   double x1, 
                   double y0, 
                   double y1,
                   double z00,
                   double z01,
                   double z10,
                   double z11)
{
    int    status = FAILURE;
    double l00,l01,l10,l11;

    /* basic check */
    if (IS_EQUAL(x0,x1) || IS_EQUAL(y0,y1))
    {
        DR_Error("linterp: x,y inputs supplied for "
                 "interpolation are the same.\n");
        goto RETURN;
    }

    l00 = (x-x1)*(y-y1)/((x0-x1)*(y0-y1));
    l01 = (x-x1)*(y-y0)/((x0-x1)*(y1-y0));
    l10 = (x-x0)*(y-y1)/((x1-x0)*(y0-y1));
    l11 = (x-x0)*(y-y0)/((x1-x0)*(y1-y0));

    *z=l00*z00+l01*z01+l10*z10+l11*z11;

    status = SUCCESS;

RETURN:

    return(status);

} /* linterp2d */



/*****  qinterp  ************************************************************/
/*
*       Quadratic interpolation based on numerical recipes.
*/
void qinterp (double  xa[],         /* (I) Set of x values */
              double  ya[],         /* (I) Set of y values */
              double  x,            /* (I) x coordinate of required point */
              double  *y,           /* (O) y value */
              double  DefaultValue) /* (I) Default value if error is too big */
{

    double  dy;              /* Error estimate */
    double  d1, d2, d3;      /* Distance on x axis */
    double  R12, R13, R23;   /* Intermediate ratios for efficiency */
    double  P12, P23;        /* Linear interp of 2 of the input points */


    d1 = xa[1] - x;
    d2 = xa[2] - x;
    d3 = xa[3] - x;
        
    R12 = 1. / (xa[1] - xa[2]);
    R13 = 1. / (xa[1] - xa[3]);
    R23 = 1. / (xa[2] - xa[3]);
                       
    P12 = R12 * (d1 * ya[2] - d2 * ya[1]);
    P23 = R23 * (d2 * ya[3] - d3 * ya[2]);
        
    *y = R13 * (d1 * P23 - d3 * P12);
        
    if (fabs (d1) < fabs (d3))
        dy = *y - P12;
    else 
        dy = *y - P23;

    if (fabs (dy) > .1 * fabs (*y))
        *y = DefaultValue;

    return;

}  /* qinterp */



/*****  sqinterp  ***********************************************************/
/*
*       Safe quadratic interpolation. The output value is "capped/floored"
*       to within a range proportional to  the difference between the max
*       input ordinate and the min input ordinate.
*/ 
int    sqinterp(double  x0,         /* (I) Set of three x values          */
                double  x1,
                double  x2,
                double  y0,         /* (I) Set of three y values          */
                double  y1,
                double  y2,
                double  d0,         /* (I) Set of three quadratic coefs   */
                double  d1,
                double  d2,
                double  x,          /* (I) x coordinate of required point */
                double  *y)         /* (O) required value                 */
{

    #define   K      0.5

    double    ymin,  ymax;
    double    h;
    double    aux;

    aux = (x - x1) * (x - x2) * d0 * y0    
        + (x - x0) * (x - x2) * d1 * y1
        + (x - x0) * (x - x1) * d2 * y2;

    /* Check if quadratic interp output inside [ymin,ymax] */
    ymin = MIN(y0,MIN(y1,y2));
    ymax = MAX(y0,MAX(y1,y2));
    h    = ymax - ymin;
      
    aux =  MAX(ymin-K*h, MIN(ymax + K*h, aux));

    *y = aux;
    return (SUCCESS);
    
}



/*****  d4interp  ***********************************************************/
/*
*       4th degree interpolation between 5 points.
*/
void    d4interp (  double  x,                  
                    double  *y, 
                    double  *xa, 
                    double  *ya)
{
    double 
            a,b,l0,l1,l2,l3,l4;

    a=(x-xa[1])*(x-xa[2])*(x-xa[3])*(x-xa[4]);
    b=(xa[0]-xa[1])*(xa[0]-xa[2])*(xa[0]-xa[3])*(xa[0]-xa[4]);
    l0=a/b;

    a=(x-xa[0])*(x-xa[2])*(x-xa[3])*(x-xa[4]);
    b=(xa[1]-xa[0])*(xa[1]-xa[2])*(xa[1]-xa[3])*(xa[1]-xa[4]);
    l1=a/b;

    a=(x-xa[0])*(x-xa[1])*(x-xa[3])*(x-xa[4]);
    b=(xa[2]-xa[0])*(xa[2]-xa[1])*(xa[2]-xa[3])*(xa[2]-xa[4]);
    l2=a/b;

    a=(x-xa[0])*(x-xa[1])*(x-xa[2])*(x-xa[4]);
    b=(xa[3]-xa[0])*(xa[3]-xa[1])*(xa[3]-xa[2])*(xa[3]-xa[4]);
    l3=a/b;

    a=(x-xa[0])*(x-xa[1])*(x-xa[2])*(x-xa[3]);
    b=(xa[4]-xa[0])*(xa[4]-xa[1])*(xa[4]-xa[2])*(xa[4]-xa[3]);
    l4=a/b;

    /*  4th degree interpolation.  */
    *y=ya[0]*l0+ya[1]*l1+ya[2]*l2+ya[3]*l3+ya[4]*l4;

    return;

}  /* d4interp */



/*****  tableinterp  ********************************************************/
/*
*       Interpolation based on a array: linear inside, flat outside.
*/
void    tableinterp (double x, 
                     double *y, 
                     double *xVal, 
                     double *yVal, 
                     int    nbVal)
{
    int k=0;

    if (x - xVal[0] < TINY) *y = yVal[0];
    else
    if (x - xVal[nbVal-1] > -TINY) *y = yVal[nbVal-1];
    else
    {
        while ((x > xVal[k]) && (k < nbVal - 1)) k++;

        if (k <=0 || k >= nbVal) 
        {
            printf ("OUCH\n");
        }

        linterp (x,y,xVal[k-1],xVal[k],yVal[k-1],yVal[k]);
    }

    return;

}  /* tableinterp */



/*****  matrixinterp  *******************************************************/
/*
*       Interpolation based on a matrix: linear inside, flat outside.
*/
void    matrixinterp (double x, 
                      double y,
                      double *z,
                      double *xVal, 
                      double *yVal, 
                      double **zVal,
                      int    xNbVal,
                      int    yNbVal)
{
    int    k=0,l=0;
    double x0, x1, y0, y1, z00, z01, z10, z11;

    if (x <= xVal[0])  
    {
        x0 = xVal[0];
        k  = 0;
    }
    else if (x >= xVal[xNbVal-1]) 
    {
        x0 = xVal[xNbVal-1];
        k  = xNbVal-2;
    }
    else
    {
        x0 = x;
        while ((k < xNbVal - 1) && (x > xVal[k+1])) k++;
    }

    if (y <= yVal[0])
    {
        y0 = yVal[0];
        l  = 0;
    }
    else if (y >= yVal[yNbVal-1])
    {
        y0 = yVal[yNbVal-1];
        l  = yNbVal-2;
    }
    else
    {
        y0 = y;
        while ((l < yNbVal - 1) && (y > yVal[l+1])) l++;
    }

    linterp2d (x0,
               y0,
               z,
               xVal[k],
               xVal[k+1],
               yVal[l],
               yVal[l+1],
               zVal[k][l],
               zVal[k][l+1],
               zVal[k+1][l],
               zVal[k+1][l+1]);
    
    return;

}  /* tableinterp */


/*****  AnnPV  *********************************************************/
/*
*       Evalute value of annuity
*/
double  AnnPV (double    AnnRate,         /* (I) Rate of annuity        */
               double    IRR,             /* (I) Discounting yield      */
               int       Tenor,           /* (I) Nb of periods remaing  */
               int       Freq)            /* (I) Frequency              */
{

    if (fabs(IRR) < TINY) 
    {
        return (AnnRate*Tenor)/Freq;
    }
    else
    {
        return AnnRate*(1.-pow(1.+IRR/Freq,-Tenor))/IRR;
    }
}


/*****  AnnIRR  ****************************************************/
/*
*      Find IRR for annuity
*/
#define IRRDELTA     0.0001
#define ANNPVERR     0.00000001
#define IRRMAX       9.99

int     AnnIRR (double    Outs,            /* (I) Outstanding            */
                double    AnnRate,         /* (I) Rate of annuity        */
                int       Tenor,           /* (I) Nb of periods remaing  */
                int       Freq,            /* (I) Frequency              */
                double   *IRR)             /* (O) IRR to be found        */
{
    int
            count = 0,
            found,
            status = FAILURE;      /* Error status = FAILURE initially */

    double 
            Yield,
            aPV, aPV2,
            y0, y1, err0, err1;


    if (Outs > (AnnRate * Tenor / Freq)) goto RETURN;
        
    if (Outs < AnnPV(AnnRate,IRRMAX,Tenor,Freq))
    {
        *IRR   = IRRMAX;
        status = SUCCESS;
        goto RETURN;
    }

    y0    = 2*Freq*(pow((Outs*Freq)/(AnnRate*Tenor),-1/(Tenor+1))-1);
    y1    = AnnRate/Outs;
    err0  = Outs-AnnPV(AnnRate,y0,Tenor,Freq);
    err1  = Outs-AnnPV(AnnRate,y1,Tenor,Freq);  
    Yield = (fabs(err0)<fabs(err1)) ? y0 : y1;
    
    found = FALSE;
    do 
    {
        aPV  = AnnPV (AnnRate,Yield,Tenor,Freq);
        aPV2 = AnnPV (AnnRate,Yield+IRRDELTA,Tenor,Freq);
        
        if (fabs(aPV-Outs) > ANNPVERR)
        {
            Yield  += (Outs - aPV) * IRRDELTA / (aPV2 - aPV);
        }
        else
        {
            found = TRUE;
        }
        count++ ;
    }
    while (!found && (count <= 10));

    if (count <= 10)
    {
        *IRR   = Yield;
        status = SUCCESS;
    }

    RETURN:                            

    return (status);

}  /* AnnIRR */




/*****  DrlSplineInterp1d  ****************************************************/
/*
*      Spline interpolation routine
*/

#define NMAX    1000

int DrlSplineInterp1d(double *xa,       /* (I) array of X(i) (i=0,..,n-1) */
                      double *ya,       /* (I) array of X(i) (i=0,..,n-1) */
                      int n,            /* (I) # of points                */
                      double yp1,       /* (I) dx/dy at X(0)              */
                      double ypn,       /* (I) dx/dy at X(n-1)            */
                      double x,         /* (I) point to intepolate        */
                      double *y)        /* (O) interpolated value         */
{

    double     y2[NMAX];
    double     xa2[NMAX];
    double     ya2[NMAX];
    int        i, counter = 1;
    int        status = FAILURE;
    
    if (n >= NMAX) {
        DR_Error(" Array too large !");
        goto RETURN;
    }

    xa2[0] = xa[0];
    ya2[0] = ya[0];

    /* Make sure the array xa is stricly increasing */
    for (i=1; i < n; i++)
    {
        if (xa[i] > xa[i-1])
        {
            xa2[counter] = xa[i];
            ya2[counter] = ya[i];
            counter ++;
        }
    }


    if (DrlSplineInterp1dInit(xa2, ya2, counter, yp1, ypn, y2) == FAILURE)
        goto RETURN;

    if (DrlSplineInterp1dInterp(xa2, ya2, y2, counter, x, y) == FAILURE)
        goto RETURN;

    status = SUCCESS;

RETURN:

    return(status);

}



/**********************************************************************
 * Interpolation 1D spline on doubles: initialization.
 */

 #define    HUGE    (0.99e11)

int DrlSplineInterp1dInit(double *x,        /* (I) array of X(i) (i=0,..,n-1)   */
                          double *y,        /* (I) array of Y(i) (i=0,..,n-1)   */
                          int n,            /* (I) # of points                  */
                          double yp1,       /* (I) dx/dy at X(0)                */
                          double ypn,       /* (I) dx/dy at X(n-1)              */
                          double *y2)       /* (O) should be allocated on entry */
{
    int i,k;
    double  p,qn,sig,un;
    double  *u;
    int status = FAILURE;


    u = (double *) DR_Array(DOUBLE, 0, n);
    if (u == NULL) return(status);

    /* nrc conventions */
    x--; y--; y2--; u--;


    if (yp1 > HUGE)
        y2[1]=u[1]=0.0;
    else {
        y2[1] = -0.5;
        u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {

        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > HUGE)
        qn=un=0.0;
    else {
        qn=0.5;
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
    for (k=n-1;k>=1;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];

    
    status = SUCCESS;
    ++u;
    Free_DR_Array(u, DOUBLE, 0, n);

    return(status);
    
}


/***********************************************************************
 * Interpolation 1D spline on doubles: interpolation.
 */

int DrlSplineInterp1dInterp(double *xa,          /* (I) array of X(i) (i=0,..,n-1) */
                            double *ya,          /* (I) array of X(i) (i=0,..,n-1) */
                            double *y2a,         /* (I) from SplineInit            */
                            int n,               /* (I) # of points                */
                            double x,            /* (I) point to intepolate        */
                            double *y)           /* (O) interpolated value         */
{
    int    klo,khi,k;
    double h,b,a;

    /* nrc2 convetions */
    xa--; ya--; y2a--;


    if (x <= xa[1]) {
        *y = ya[1];
        return(SUCCESS);
    } 
    else if (x >= xa[n])
    {
        *y = ya[n];


        return(SUCCESS);
    }



    klo=1;
    khi=n;
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xa[k] > x) khi=k;
        else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) {
        DR_Error("bad xa input.\n!");
        return(FAILURE);
    }
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+
        (b*b*b-b)*y2a[khi])*(h*h)/6.0;

    if ((*y < ya[klo]) || (*y > ya[khi]) ) 
    {
        linterp (x,y,xa[klo],xa[khi],ya[klo],ya[khi]);
    }


    return(SUCCESS);
}







