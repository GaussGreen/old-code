#include <math.h>
#include <stdlib.h>
#include <stdio.h>
/*#include <config.h>*/
#include "error2.h"
#include "proba_utils.h"
#include "integration_gsl.h"

/**************************************************************/
/*    UTILITIES                                               */
/**************************************************************/
static /*inline*/ void qpsrt (gsl_integration_workspace * workspace);

#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_FN_EVAL(F,x) (*((F)->PTR_FUNCTION))(x,(F)->params)

#undef HAVE_EXTENDED_PRECISION_REGISTERS

#ifdef HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

static double rescale_error (double err, const double result_abs, const double result_asc)
{
    err = fabs(err) ;

    if (result_asc != 0 && err != 0)
    {
	    double scale = pow((200 * err / result_asc), 1.5) ;
	
	    if (scale < 1)
	    {
	        err = result_asc * scale ;
	    }
	    else 
	    {
	        err = result_asc ;
	    }
    }
    if (result_abs > GSL_DBL_MIN / (50 * GSL_DBL_EPSILON))
    {
        double min_err = 50 * GSL_DBL_EPSILON * result_abs ;

        if (min_err > err) 
	    {
	        err = min_err ;
	    }
    }
  
    return err ;
}


/* Definition of an integration rule */
typedef void gsl_integration_rule ( const PTR_FUNCTION f,
				                    double a,
                                    double b,
				                    double *result,
                                    double *abserr,
				                    double *defabs,
                                    double *resabs,
                                    void *param);

/****************************************************/
/* INTEGRATION METHODS                              */
/****************************************************/
void gsl_integration_qk15 ( const PTR_FUNCTION f,double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc, void *param);

void gsl_integration_qk21 (const PTR_FUNCTION f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc, void *param);

void gsl_integration_qk31 (const PTR_FUNCTION f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc, void *param);

void gsl_integration_qk41 (const PTR_FUNCTION f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc, void *param);

void gsl_integration_qk51 (const PTR_FUNCTION f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc, void *param);

void gsl_integration_qk61 (const PTR_FUNCTION f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc, void *param);

void gsl_integration_qcheb (PTR_FUNCTION f, double a, double b, 
                            double *cheb12, double *cheb24, void *param);

/* The low-level integration rules in QUADPACK are identified by small
   integers (1-6). We'll use symbolic constants to refer to them.  */

enum
{
    GSL_INTEG_GAUSS15 = 1,	/* 15 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS21 = 2,	/* 21 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS31 = 3,	/* 31 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS41 = 4,	/* 41 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS51 = 5,	/* 51 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS61 = 6	/* 61 point Gauss-Kronrod rule */
};

/****************************************************/
/* WORKSPACE PTR_FUNCTIONS                          */
/****************************************************/
/**********************************/
/* gsl_integration_workspace_alloc*/
/**********************************/
gsl_integration_workspace * gsl_integration_workspace_alloc (const size_t n) 
{
    static char routine[] = "gsl_integration_workspace_alloc";
    int status = FAILURE;
    gsl_integration_workspace * w = NULL;
  
    if (n == 0)
    {
        DR_Error ("workspace length n must be positive integer");
        goto RETURN;
    }

    w = (gsl_integration_workspace *) 
    malloc (sizeof (gsl_integration_workspace));
    if (w==NULL)
    {
        DR_Error ("failed to allocate space for workspace struct");
        goto RETURN;
    }

    w->alist = (double *) malloc (n * sizeof (double));
    if(w->alist==NULL)
    {
        DR_Error("failed to allocate space for alist ranges");
        goto RETURN;
    }

    w->blist = (double *) malloc (n * sizeof (double));
    if(w->blist==NULL)
    {
        DR_Error("failed to allocate space for blist ranges");
        goto RETURN;
    }

    w->rlist = (double *) malloc (n * sizeof (double));
    if(w->rlist==NULL)
    {
        DR_Error("failed to allocate space for rlist ranges");
        goto RETURN;
    }

    w->elist = (double *) malloc (n * sizeof (double));
    if(w->elist==NULL)
    {
        DR_Error("failed to allocate space for elist ranges");
        goto RETURN;
    }

    w->order = (size_t *) malloc (n * sizeof (size_t));
    if(w->order==NULL)
    {
        DR_Error("failed to allocate space for order ranges");
        goto RETURN;
    }

    w->level = (size_t *) malloc (n * sizeof (size_t));
    if(w->level==NULL)
    {
        DR_Error("failed to allocate space for level ranges");
        goto RETURN;
    }

    w->size = 0 ;
    w->limit = n ;
    w->maximum_level = 0 ;
    status = SUCCESS;

RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        if(w)
        {
            if(w->level) free (w->level);
            if(w->order) free (w->order);
            if(w->elist) free (w->elist);
            if(w->rlist) free (w->rlist);
            if(w->blist) free (w->blist);
            if(w->alist) free (w->alist);
            free (w);
            w = NULL;
        }

    }
    return w ;
}

void gsl_integration_workspace_free (gsl_integration_workspace * w)
{
    if(w)
    {
        if(w->level) free (w->level);
        if(w->order) free (w->order);
        if(w->elist) free (w->elist);
        if(w->rlist) free (w->rlist);
        if(w->blist) free (w->blist);
        if(w->alist) free (w->alist);
        free (w);
    }
}



/****************************/
/* update                   */
/****************************/
static /*inline*/ void update ( gsl_integration_workspace * workspace,
	                            double a1,
                                double b1,
                                double area1,
                                double error1,
	                            double a2,
                                double b2,
                                double area2,
                                double error2)
{
    double * alist = workspace->alist ;
    double * blist = workspace->blist ;
    double * rlist = workspace->rlist ;
    double * elist = workspace->elist ;
    size_t * level = workspace->level ;
    
    const size_t i_max = workspace->i ;
    const size_t i_new = workspace->size ;
    
    const size_t new_level = workspace->level[i_max] + 1;

    /* append the newly-created intervals to the list */
    if (error2 > error1)
    {
        alist[i_max] = a2;	/* blist[maxerr] is already == b2 */
        rlist[i_max] = area2;
        elist[i_max] = error2;
        level[i_max] = new_level;
        
        alist[i_new] = a1;
        blist[i_new] = b1;
        rlist[i_new] = area1;
        elist[i_new] = error1;
        level[i_new] = new_level;
    }
    else
    {
        blist[i_max] = b1;	/* alist[maxerr] is already == a1 */
        rlist[i_max] = area1;
        elist[i_max] = error1;
        level[i_max] = new_level;
        
        alist[i_new] = a2;
        blist[i_new] = b2;
        rlist[i_new] = area2;
        elist[i_new] = error2;
        level[i_new] = new_level;
    }
  
    workspace->size++;

    if (new_level > workspace->maximum_level)
    {
        workspace->maximum_level = new_level;
    }

    qpsrt (workspace) ;
}

/****************************/
/* retrieve                 */
/****************************/
static /*inline*/ void retrieve (   const gsl_integration_workspace * workspace, 
	                            double *a,
                                double *b,
                                double *r,
                                double *e)
{
    const size_t i = workspace->i;
    double * alist = workspace->alist;
    double * blist = workspace->blist;
    double * rlist = workspace->rlist;
    double * elist = workspace->elist;
    
    *a = alist[i] ;
    *b = blist[i] ;
    *r = rlist[i] ;
    *e = elist[i] ;
}

/****************************/
/* sum_results              */
/****************************/
static /*inline*/ double sum_results ( const gsl_integration_workspace * workspace)
{
    const double * const rlist = workspace->rlist ;
    const size_t n = workspace->size;
    
    size_t k;
    double result_sum = 0;
    
    for (k = 0; k < n; k++)
    {
        result_sum += rlist[k];
    }
  
    return result_sum;
}

/****************************/
/* subinterval_too_small    */
/****************************/
static /*inline*/ int subinterval_too_small (double a1, double a2, double b2)
{
    const double e = GSL_DBL_EPSILON;
    const double u = GSL_DBL_MIN;
    
    double tmp = (1 + 100 * e) * (fabs (a2) + 1000 * u);
    
    int status = fabs (a1) <= tmp && fabs (b2) <= tmp;
    
    return status;
}

/****************************/
/* initialise               */
/****************************/
static /*inline*/ void initialise ( gsl_integration_workspace * workspace,
                                    double a,
                                    double b)
{
    workspace->size = 0;
    workspace->nrmax = 0;
    workspace->i = 0;
    workspace->alist[0] = a;
    workspace->blist[0] = b;
    workspace->rlist[0] = 0.0;
    workspace->elist[0] = 0.0;
    workspace->order[0] = 0;
    workspace->level[0] = 0;
    
    workspace->maximum_level = 0;
}

/****************************/
/* set_initial_result       */
/****************************/
static /*inline*/ void set_initial_result ( gsl_integration_workspace * workspace, 
			                                double result,
                                            double error)
{
    workspace->size = 1;
    workspace->rlist[0] = result;
    workspace->elist[0] = error;
}

/****************************/
/* increase_nrmax           */
/****************************/
static int increase_nrmax (gsl_integration_workspace * workspace)
{
    int k;
    int id = workspace->nrmax;
    int jupbnd;
    
    const size_t * level = workspace->level;
    const size_t * order = workspace->order;
    
    size_t limit = workspace->limit ;
    size_t last = workspace->size - 1 ;
    
    if (last > (1 + limit / 2))
    {
        jupbnd = limit + 1 - last;
    }
    else
    {
        jupbnd = last;
    }
  
    for (k = id; k <= jupbnd; k++)
    {
        size_t i_max = order[workspace->nrmax];
      
        workspace->i = i_max ;

        if (level[i_max] < workspace->maximum_level)
	    {
	        return 1;
	    }

        workspace->nrmax++;
    }
  return 0;
}

/****************************/
/* reset_nrmax              */
/****************************/
static /*inline*/ void reset_nrmax (gsl_integration_workspace * workspace)
{
    workspace->nrmax = 0;
    workspace->i = workspace->order[0] ;
}

/****************************/
/* large_interval           */
/****************************/
static int large_interval (gsl_integration_workspace * workspace)
{
    size_t i = workspace->i ;
    const size_t * level = workspace->level;
  
    if (level[i] < workspace->maximum_level)
    {
        return 1 ;
    }
    else
    {
        return 0 ;
    }
}

/****************************/
/* test_positivity          */
/****************************/
/* Compare the integral of f(x) with the integral of |f(x)|
   to determine if f(x) covers both positive and negative values */
static /*inline*/ int test_positivity ( double result,
                                        double resabs)
{
    int status = (fabs (result) >= (1 - 50 * GSL_DBL_EPSILON) * resabs);

    return status;
}


/************************************************************/
/* extrapolation table utilities                            */
/************************************************************/
struct extrapolation_table
{
    size_t n;
    double rlist2[52];
    size_t nres;
    double res3la[3];
};

/****************************/
/* initialise_table         */
/****************************/
static void initialise_table (struct extrapolation_table *table)
{
    table->n = 0;
    table->nres = 0;
}
#ifdef JUNK
static void initialise_table (  struct extrapolation_table *table,
                                double y)
{
    table->n = 0;
    table->rlist2[0] = y;
    table->nres = 0;
}
#endif
/****************************/
/* append_table             */
/****************************/
static void append_table (  struct extrapolation_table *table,
                            double y)
{
    size_t n;
    n = table->n;
    table->rlist2[n] = y;
    table->n++;
}



/****************************/
/* qelg                     */
/****************************/

static /*inline*/ void qelg (struct extrapolation_table *table,
                             double *result,
                             double *abserr)
{
    double *epstab = table->rlist2;
    double *res3la = table->res3la;
    const size_t n = table->n - 1;
    
    const double current = epstab[n];
    
    double absolute = GSL_DBL_MAX;
    double relative = 5 * GSL_DBL_EPSILON * fabs (current);
    
    const size_t newelm = n / 2;
    const size_t n_orig = n;
    size_t n_final = n;
    size_t i;
    
    const size_t nres_orig = table->nres;
    
    *result = current;
    *abserr = GSL_DBL_MAX;
    
    if (n < 2)
    {
        *result = current;
        *abserr = GSL_MAX (absolute, relative);
        return;
    }

    epstab[n + 2] = epstab[n];
    epstab[n] = GSL_DBL_MAX;
    
    for (i = 0; i < newelm; i++)
    {
        double res = epstab[n - 2 * i + 2];
        double e0 = epstab[n - 2 * i - 2];
        double e1 = epstab[n - 2 * i - 1];
        double e2 = res;
        
        double e1abs = fabs (e1);
        double delta2 = e2 - e1;
        double err2 = fabs (delta2);
        double tol2 = GSL_MAX (fabs (e2), e1abs) * GSL_DBL_EPSILON;
        double delta3 = e1 - e0;
        double err3 = fabs (delta3);
        double tol3 = GSL_MAX (e1abs, fabs (e0)) * GSL_DBL_EPSILON;
        
        double e3, delta1, err1, tol1, ss;
        
        if (err2 <= tol2 && err3 <= tol3)
	    {
	    /* If e0, e1 and e2 are equal to within machine accuracy,
	        convergence is assumed.  */

	        *result = res;
	        absolute = err2 + err3;
	        relative = 5 * GSL_DBL_EPSILON * fabs (res);
	        *abserr = GSL_MAX (absolute, relative);
	        return;
	    }

        e3 = epstab[n - 2 * i];
        epstab[n - 2 * i] = e1;
        delta1 = e1 - e3;
        err1 = fabs (delta1);
        tol1 = GSL_MAX (e1abs, fabs (e3)) * GSL_DBL_EPSILON;

        /* If two elements are very close to each other, omit a part of
            the table by adjusting the value of n */

        if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
	    {
	        n_final = 2 * i;
	        break;
	    }

        ss = (1 / delta1 + 1 / delta2) - 1 / delta3;

        /* Test to detect irregular behaviour in the table, and
            eventually omit a part of the table by adjusting the value of
            n. */

        if (fabs (ss * e1) <= 0.0001)
	    {
	        n_final = 2 * i;
	        break;
	    }

        /* Compute a new element and eventually adjust the value of
            result. */

        res = e1 + 1 / ss;
        epstab[n - 2 * i] = res;

        {
	        const double error = err2 + fabs (res - e2) + err3;

	        if (error <= *abserr)
	        {
	            *abserr = error;
	            *result = res;
	        }
        }
    }

    /* Shift the table */

    {
        const size_t limexp = 50 - 1;

        if (n_final == limexp)
        {
	        n_final = 2 * (limexp / 2);
        }
    }

    if (n_orig % 2 == 1)
    {
        for (i = 0; i <= newelm; i++)
	    {
	        epstab[1 + i * 2] = epstab[i * 2 + 3];
	    }
    }
    else
    {
        for (i = 0; i <= newelm; i++)
	    {
	        epstab[i * 2] = epstab[i * 2 + 2];
	    }
    }

    if (n_orig != n_final)
    {
        for (i = 0; i <= n_final; i++)
	    {
	        epstab[i] = epstab[n_orig - n_final + i];
	    }
    }

    table->n = n_final + 1;

    if (nres_orig < 3)
    {
        res3la[nres_orig] = *result;
        *abserr = GSL_DBL_MAX;
    }
    else
    {				/* Compute error estimate */
        *abserr = (fabs (*result - res3la[2]) + fabs (*result - res3la[1])
		 + fabs (*result - res3la[0]));

        res3la[0] = res3la[1];
        res3la[1] = res3la[2];
        res3la[2] = *result;
    }

  /* In QUADPACK the variable table->nres is incremented at the top of
     qelg, so it increases on every call. This leads to the array
     res3la being accessed when its elements are still undefined, so I
     have moved the update to this point so that its value more
     useful. */

    table->nres = nres_orig + 1;  

    *abserr = GSL_MAX (*abserr, 5 * GSL_DBL_EPSILON * fabs (*result));

    return;
}


/********************************************************************/
/* INTEGRATION ON A FINITE INTERVAL                                 */
/********************************************************************/

static int qag (    const PTR_FUNCTION f,
                    const double a,
                    const double b,
                    const double epsabs,
                    const double epsrel,
                    const size_t limit,
                    gsl_integration_workspace * workspace,
                    double * result,
                    double * abserr,
                    gsl_integration_rule * q,
                    void *param) ;
/****************************/
/* gsl_integration_qag      */
/****************************/
int gsl_integration_qag (   const PTR_FUNCTION f,
	            	        double a,
                            double b,
		                    double epsabs, 
                            double epsrel,
                            size_t limit,
		                    int key,
		                    gsl_integration_workspace * workspace,
		                    double * result,
                            double * abserr,
                            void *param)
{
    static char routine[] = "gsl_integration_qag";
    int status = FAILURE;
    gsl_integration_rule * integration_rule = gsl_integration_qk15 ;

    if (key < GSL_INTEG_GAUSS15)
    {
        key = GSL_INTEG_GAUSS15 ;
    } 
    else if (key > GSL_INTEG_GAUSS61) 
    {
      key = GSL_INTEG_GAUSS61 ;
    }

    switch (key) 
    {
    case GSL_INTEG_GAUSS15:
        integration_rule = gsl_integration_qk15 ;
        break ;
    case GSL_INTEG_GAUSS21:
        integration_rule = gsl_integration_qk21 ;
        break ;
    case GSL_INTEG_GAUSS31:
        integration_rule = gsl_integration_qk31 ; 
        break ;
    case GSL_INTEG_GAUSS41:
        integration_rule = gsl_integration_qk41 ;
        break ;      
    case GSL_INTEG_GAUSS51:
        integration_rule = gsl_integration_qk51 ;
        break ;      
    case GSL_INTEG_GAUSS61:
        integration_rule = gsl_integration_qk61 ;
        break ;      
    default:
        DR_Error("value of key does specify a known integration rule") ;
        goto RETURN;
    }

    status = qag (f, a, b, epsabs, epsrel, limit,
                workspace, 
                result, abserr, 
                integration_rule,
                param) ;
  
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
  return status ;
}

/****************************/
/* qag                      */
/****************************/
static int qag (const PTR_FUNCTION f,
                const double a,
                const double b,
                const double epsabs,
                const double epsrel,
                const size_t limit,
                gsl_integration_workspace * workspace,
                double *result, double *abserr,
                gsl_integration_rule * q,
                void *param)
{
    static char routine[] = "qag";
    double area, errsum;
    double result0, abserr0, resabs0, resasc0;
    double tolerance;
    size_t iteration = 0;
    int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;
    int status = FAILURE;
    double round_off;	
    
    /* Initialize results */
    
    initialise (workspace, a, b);
    
    *result = 0;
    *abserr = 0;
    
    if (limit > workspace->limit)
    {
        DR_Error ("iteration limit exceeds available workspace") ;
        goto RETURN;
    }

    if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28))
    {
        DR_Error ("tolerance cannot be acheived with given epsabs and epsrel");
        goto RETURN;
    }

    /* perform the first integration */

    q (f, a, b, &result0, &abserr0, &resabs0, &resasc0, param);

    set_initial_result (workspace, result0, abserr0);

    /* Test on accuracy */

    tolerance = GSL_MAX (epsabs, epsrel * fabs (result0));

    /* need IEEE rounding here to match original quadpack behavior */

    round_off = GSL_COERCE_DBL (50 * GSL_DBL_EPSILON * resabs0);

    if (abserr0 <= round_off && abserr0 > tolerance)
    {
        *result = result0;
        *abserr = abserr0;

        DR_Error ("cannot reach tolerance because of roundoff error "
		 "on first attempt");
        goto RETURN;
    }
    else if ((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0.0)
    {
        *result = result0;
        *abserr = abserr0;

        status =  SUCCESS;
        goto RETURN;
    }
    else if (limit == 1)
    {
        *result = result0;
        *abserr = abserr0;

        DR_Error ("a maximum of one iteration was insufficient");
        goto RETURN;
    }

    area = result0;
    errsum = abserr0;

    iteration = 1;

    do
    {
        double a1, b1, a2, b2;
        double a_i, b_i, r_i, e_i;
        double area1 = 0, area2 = 0, area12 = 0;
        double error1 = 0, error2 = 0, error12 = 0;
        double resasc1, resasc2;
        double resabs1, resabs2;

        /* Bisect the subinterval with the largest error estimate */

        retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

        a1 = a_i; 
        b1 = 0.5 * (a_i + b_i);
        a2 = b1;
        b2 = b_i;
        
        q (f, a1, b1, &area1, &error1, &resabs1, &resasc1, param);
        q (f, a2, b2, &area2, &error2, &resabs2, &resasc2, param);
        
        area12 = area1 + area2;
        error12 = error1 + error2;
        
        errsum += (error12 - e_i);
        area += area12 - r_i;
        
        if (resasc1 != error1 && resasc2 != error2)
	    {
	        double delta = r_i - area12;
            
	        if (fabs (delta) <= 1.0e-5 * fabs (area12) && error12 >= 0.99 * e_i)
	        {
	            roundoff_type1++;
	        }
	        if (iteration >= 10 && error12 > e_i)
	        {
	            roundoff_type2++;
	        }
	    }

        tolerance = GSL_MAX (epsabs, epsrel * fabs (area));

        if (errsum > tolerance)
	    {
	        if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
	        {
	            error_type = 2;	/* round off error */
	        }

	        /* set error flag in the case of bad integrand behaviour at
	        a point of the integration range */

	        if (subinterval_too_small (a1, a2, b2))
	        {
	            error_type = 3;
	        }
	    }

        update (workspace, a1, b1, area1, error1, a2, b2, area2, error2);

        retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

        iteration++;

    }
    while (iteration < limit && !error_type && errsum > tolerance);

    *result = sum_results (workspace);
    *abserr = errsum;

    if (errsum <= tolerance)
    {
        status = SUCCESS;
        goto RETURN;
    }
    else if (error_type == 2)
    {
        DR_Error ("roundoff error prevents tolerance from being achieved");
        status = FAILURE;
        goto RETURN;
    }
    else if (error_type == 3)
    {
        DR_Error ("bad integrand behavior found in the integration interval");
        status = FAILURE;
        goto RETURN;
    }
    else if (iteration == limit)
    {
        DR_Error ("maximum number of subdivisions reached");
        status = FAILURE;
        goto RETURN;
    }
    else
    {
        DR_Error ("could not integrate PTR_FUNCTION");
        status = FAILURE;
        goto RETURN;
    }

RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return status;
}


/****************************/
/* qpsrt                    */
/****************************/
static /*inline*/ void qpsrt (gsl_integration_workspace * workspace)
{
    const size_t last = workspace->size - 1;
    const size_t limit = workspace->limit;
    double * elist = workspace->elist;
    size_t * order = workspace->order;
    double errmax ;
    double errmin ;
    int i, k, top;
    
    size_t i_nrmax = workspace->nrmax;
    size_t i_maxerr = order[i_nrmax] ;
    
    /* Check whether the list contains more than two error estimates */
    
    if (last < 2) 
    {
        order[0] = 0 ;
        order[1] = 1 ;
        workspace->i = i_maxerr ;
        return ;
    }

    errmax = elist[i_maxerr] ;
    
    /* This part of the routine is only executed if, due to a difficult
       integrand, subdivision increased the error estimate. In the normal
       case the insert procedure should start after the nrmax-th largest
       error estimate. */
    
    while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) 
    {
        order[i_nrmax] = order[i_nrmax - 1] ;
        i_nrmax-- ;
    } 

    /* Compute the number of elements in the list to be maintained in
       descending order. This number depends on the number of
       subdivisions still allowed. */
  
    if(last < (limit/2 + 2)) 
    {
        top = last ;
    }
    else
    {
        top = limit - last + 1;
    }
  
    /* Insert errmax by traversing the list top-down, starting
       comparison from the element elist(order(i_nrmax+1)). */
  
    i = i_nrmax + 1 ;
  
    /* The order of the tests in the following line is important to
       prevent a segmentation fault */

    while (i < top && errmax < elist[order[i]])
    {
        order[i-1] = order[i] ;
        i++ ;
    }
  
    order[i-1] = i_maxerr ;
  
    /* Insert errmin by traversing the list bottom-up */
  
    errmin = elist[last] ;
  
    k = top - 1 ;
  
    while (k > i - 2 && errmin >= elist[order[k]])
    {
        order[k+1] = order[k] ;
        k-- ;
    }
  
    order[k+1] = last ;

    /* Set i_max and e_max */

    i_maxerr = order[i_nrmax] ;
  
    workspace->i = i_maxerr ;
    workspace->nrmax = i_nrmax ;
}

/********************************************************************/
/* INTEGRATION ON AN INFINITE INTERVAL                              */
/********************************************************************/
/****************************/
/* gsl_integration_qags     */
/****************************/
static int qags (   const PTR_FUNCTION f,
                    const double a,
                    const double b,
                    const double epsabs,
                    const double epsrel,
                    const size_t limit,
                    gsl_integration_workspace * workspace,
                    double *result,
                    double *abserr,
                    gsl_integration_rule * q,
                    void *param);

int gsl_integration_qags (  const PTR_FUNCTION f,
		                    double a,
                            double b,
		                    double epsabs,
                            double epsrel,
                            size_t limit,
		                    gsl_integration_workspace * workspace,
		                    double * result,
                            double * abserr,
                            void *param)
{
    static char routine[] = "gsl_integration_qags";
    int status = qags (f, a, b, epsabs, epsrel, limit,
                     workspace, 
                     result, abserr, 
                     &gsl_integration_qk21, param) ;
    
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status ;
}

/* QAGI: evaluate an integral over an infinite range using the
   transformation

   integrate(f(x),-Inf,Inf) = integrate((f((1-t)/t) + f(-(1-t)/t))/t^2,0,1)

   */
typedef struct
{
    void *param;
    PTR_FUNCTION f;
} Transform_Param_i;

int i_transform (double t, double *y, void *param)
{
    Transform_Param_i *f_transform_param = (Transform_Param_i *) param;
    double y1,y2;
    void *f_param = f_transform_param->param;
    PTR_FUNCTION f = (PTR_FUNCTION) f_transform_param->f;
    double x = (1 - t) / t;
    (*f)(x,&y1,f_param);
    (*f)(-x,&y2,f_param);
    *y = y1 + y2;
    *y = ((*y) / t) / t;
    return SUCCESS;
}



/****************************/
/* gsl_integration_qagi     */
/****************************/
int gsl_integration_qagi (  PTR_FUNCTION f,
		                    double epsabs,
                            double epsrel,
                            size_t limit,
		                    gsl_integration_workspace * workspace,
		                    double *result,
                            double *abserr,
                            void *param)
{
    int status;
    
    PTR_FUNCTION f_transform = &i_transform;
    Transform_Param_i f_transform_param;
    f_transform_param.f = f;
    f_transform_param.param = param;
    
    status = qags (f_transform, 0.0, 1.0, 
                   epsabs, epsrel, limit,
                   workspace,
                   result, abserr,
                   &gsl_integration_qk15,
                   &f_transform_param);
    
    return status;
}



/* QAGIL: Evaluate an integral over an infinite range using the
   transformation,
   
   integrate(f(x),-Inf,b) = integrate(f(b-(1-t)/t)/t^2,0,1)

   */

typedef struct
{
    void *param;
    double b;
    PTR_FUNCTION f;
} Transform_Param_il;

int il_transform (double t,
                  double *y,
                  void *param)
{
    Transform_Param_il *transform_param = (Transform_Param_il *) param;
    double b = transform_param->b;
    PTR_FUNCTION f = transform_param->f;
    void *f_param = transform_param->param;
    double x = b - (1 - t) / t;
    (*f)(x,y,f_param);
    *y =  (*y / t) / t;
    return SUCCESS;

}

/****************************/
/* gsl_integration_qagil    */
/****************************/
int gsl_integration_qagil ( PTR_FUNCTION f,
		                    double b,
            		        double epsabs,
                            double epsrel,
                            size_t limit,
		                    gsl_integration_workspace * workspace,
		                    double *result,
                            double *abserr,
                            void *param)
{
    static char routine[] = "gsl_integration_qagil";
    int status = FAILURE;
    PTR_FUNCTION f_transform = &il_transform;
    Transform_Param_il f_transform_param;
    f_transform_param.f = f;
    f_transform_param.param = param;
    f_transform_param.b = b;
    
    status = qags (f_transform, 0.0, 1.0, 
                   epsabs, epsrel, limit,
                   workspace,
                   result, abserr,
                   &gsl_integration_qk15,
                   &f_transform_param);
    
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}



/* QAGIU: Evaluate an integral over an infinite range using the
   transformation

   integrate(f(x),a,Inf) = integrate(f(a+(1-t)/t)/t^2,0,1)

   */

typedef struct
{
    void *param;
    double a;
    PTR_FUNCTION f;
} Transform_Param_iu;

int iu_transform (    double t,
                      double *y,
                      void *param)
{
    Transform_Param_iu *transform_param = (Transform_Param_iu *) param;
    double a = transform_param->a;
    PTR_FUNCTION f = transform_param->f;
    void *f_param = transform_param->param;
    double x = a + (1 - t) / t;
    (*f)(x,y,f_param);
    *y =  (*y / t) / t;
    return SUCCESS;

}

/****************************/
/* gsl_integration_qagiu    */
/****************************/
int gsl_integration_qagiu ( PTR_FUNCTION f,
		                    double a,
		                    double epsabs,
                            double epsrel,
                            size_t limit,
		                    gsl_integration_workspace * workspace,
		                    double *result,
                            double *abserr,
                            void *param)
{
    static char routine[] = "gsl_integration_qagiu";
    int status;
    
    PTR_FUNCTION f_transform = &iu_transform;
    Transform_Param_iu f_transform_param;
    f_transform_param.f = f;
    f_transform_param.param = param;
    f_transform_param.a = a;
    
    status = qags (f_transform, 0.0, 1.0, 
                   epsabs, epsrel, limit,
                   workspace,
                   result, abserr,
                   &gsl_integration_qk15,
                   &f_transform_param);

    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}



/* Main integration function */
/****************************/
/* qags                     */
/****************************/
static int qags (const PTR_FUNCTION f,
      const double a, const double b,
      const double epsabs, const double epsrel,
      const size_t limit,
      gsl_integration_workspace * workspace,
      double *result, double *abserr,
      gsl_integration_rule * q,
      void * param)
{
    static char routine[] = "qags";
    int status = FAILURE;
    double area, errsum;
    double res_ext, err_ext;
    double result0, abserr0, resabs0, resasc0;
    double tolerance;
    double ertest = 0;
    double error_over_large_intervals = 0;
    double reseps = 0, abseps = 0, correc = 0;
    size_t ktmin = 0;
    int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 = 0;
    int error_type = 0, error_type2 = 0;
    
    size_t iteration = 0;
    
    int positive_integrand = 0;
    int extrapolate = 0;
    int disallow_extrapolation = 0;
    
    struct extrapolation_table table;
    
    /* Initialize results */
    
    initialise (workspace, a, b);
    
    *result = 0;
    *abserr = 0;
    
    if (limit > workspace->limit)
    {
        DR_Error ("iteration limit exceeds available workspace") ;
        goto RETURN;
    }
    
    /* Test on accuracy */
    
    if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28))
    {
        DR_Error ("tolerance cannot be acheived with given epsabs and epsrel");
        goto RETURN;
    }
    
    /* Perform the first integration */
    
    q (f, a, b, &result0, &abserr0, &resabs0, &resasc0, param);
    
    set_initial_result (workspace, result0, abserr0);
    
    tolerance = GSL_MAX (epsabs, epsrel * fabs (result0));
    
    if (abserr0 <= 100 * GSL_DBL_EPSILON * resabs0 && abserr0 > tolerance)
    {
        *result = result0;
        *abserr = abserr0;
    
        DR_Error ("cannot reach tolerance because of roundoff error"
    		 "on first attempt");
        goto RETURN;
    }
    else if ((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0.0)
    {
        *result = result0;
        *abserr = abserr0;
    
        return SUCCESS;
    }
    else if (limit == 1)
    {
        *result = result0;
        *abserr = abserr0;
    
        DR_Error ("a maximum of one iteration was insufficient");
        goto RETURN;
    }
    
    /* Initialization */
    
    initialise_table (&table);
    append_table (&table, result0);
    
    area = result0;
    errsum = abserr0;
    
    res_ext = result0;
    err_ext = GSL_DBL_MAX;
    
    positive_integrand = test_positivity (result0, resabs0);
    
    iteration = 1;
    
    do
    {
        size_t current_level;
        double a1, b1, a2, b2;
        double a_i, b_i, r_i, e_i;
        double area1 = 0, area2 = 0, area12 = 0;
        double error1 = 0, error2 = 0, error12 = 0;
        double resasc1, resasc2;
        double resabs1, resabs2;
        double last_e_i;
    
        /* Bisect the subinterval with the largest error estimate */
    
        retrieve (workspace, &a_i, &b_i, &r_i, &e_i);
    
        current_level = workspace->level[workspace->i] + 1;
    
        a1 = a_i;
        b1 = 0.5 * (a_i + b_i);
        a2 = b1;
        b2 = b_i;
    
        iteration++;
    
        q (f, a1, b1, &area1, &error1, &resabs1, &resasc1, param);
        q (f, a2, b2, &area2, &error2, &resabs2, &resasc2, param);
    
        area12 = area1 + area2;
        error12 = error1 + error2;
        last_e_i = e_i;
    
        /* Improve previous approximations to the integral and test for
           accuracy.
    
           We write these expressions in the same way as the original
           QUADPACK code so that the rounding errors are the same, which
           makes testing easier. */
    
        errsum = errsum + error12 - e_i;
        area = area + area12 - r_i;
    
        tolerance = GSL_MAX (epsabs, epsrel * fabs (area));
    
        if (resasc1 != error1 && resasc2 != error2)
    	{
    	    double delta = r_i - area12;
    
    	    if (fabs (delta) <= 1.0e-5 * fabs (area12) && error12 >= 0.99 * e_i)
    	    {
    	        if (!extrapolate)
    		    {
    		        roundoff_type1++;
    		    }
    	        else
    		    {
    		        roundoff_type2++;
    		    }
    	    }
    	    if (iteration > 10 && error12 > e_i)
    	    {
    	        roundoff_type3++;
    	    }
    	}
    
        /* Test for roundoff and eventually set error flag */
    
        if (roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20)
    	{
    	    error_type = 2;	/* round off error */
    	}
    
        if (roundoff_type2 >= 5)
    	{
    	    error_type2 = 1;
    	}
    
        /* set error flag in the case of bad integrand behaviour at
           a point of the integration range */
    
        if (subinterval_too_small (a1, a2, b2))
    	{
    	    error_type = 4;
    	}
    
        /* append the newly-created intervals to the list */
    
        update (workspace, a1, b1, area1, error1, a2, b2, area2, error2);
    
        if (errsum <= tolerance)
    	{
    	    goto compute_result;
    	}
    
        if (error_type)
    	{
    	    break;
    	}
    
        if (iteration >= limit - 1)
    	{
    	    error_type = 1;
    	    break;
    	}
    
        if (iteration == 2)	/* set up variables on first iteration */
    	{
    	    error_over_large_intervals = errsum;
    	    ertest = tolerance;
    	    append_table (&table, area);
    	    continue;
    	}
    
        if (disallow_extrapolation)
    	{
    	    continue;
    	}
    
        error_over_large_intervals += -last_e_i;
    
        if (current_level < workspace->maximum_level)
    	{
    	    error_over_large_intervals += error12;
    	}
    
        if (!extrapolate)
    	{
    	    /* test whether the interval to be bisected next is the
    	       smallest interval. */
    
    	    if (large_interval (workspace))
    	        continue;
    
    	    extrapolate = 1;
    	    workspace->nrmax = 1;
    	}
    
        if (!error_type2 && error_over_large_intervals > ertest)
    	{
    	    if (increase_nrmax (workspace))
    	        continue;
    	}
    
        /* Perform extrapolation */
    
        append_table (&table, area);
    
        qelg (&table, &reseps, &abseps);
    
        ktmin++;
    
        if (ktmin > 5 && err_ext < 0.001 * errsum)
    	{
    	    error_type = 5;
    	}
    
        if (abseps < err_ext)
    	{
    	    ktmin = 0;
    	    err_ext = abseps;
    	    res_ext = reseps;
    	    correc = error_over_large_intervals;
    	    ertest = GSL_MAX (epsabs, epsrel * fabs (reseps));
    	    if (err_ext <= ertest)
    	    break;
    	}
    
        /* Prepare bisection of the smallest interval. */
    
        if (table.n == 1)
    	{
    	    disallow_extrapolation = 1;
    	}
    
        if (error_type == 5)
    	{
    	    break;
    	}
    
        /* work on interval with largest error */
    
        reset_nrmax (workspace);
        extrapolate = 0;
        error_over_large_intervals = errsum;
    
    }
    while (iteration < limit);
    
    *result = res_ext;
    *abserr = err_ext;
    
    if (err_ext == GSL_DBL_MAX)
        goto compute_result;
    
    if (error_type || error_type2)
    {
        if (error_type2)
    	{
    	    err_ext += correc;
    	}
    
        if (error_type == 0)
    	    error_type = 3;
    
        if (res_ext != 0.0 && area != 0.0)
    	{
    	    if (err_ext / fabs (res_ext) > errsum / fabs (area))
    	    goto compute_result;
    	}
        else if (err_ext > errsum)
    	{
    	    goto compute_result;
    	}
        else if (area == 0.0)
    	{
    	    goto RETURN;
    	}
    }
    
    /*  Test on divergence. */
    
    {
        double max_area = GSL_MAX (fabs (res_ext), fabs (area));
    
        if (!positive_integrand && max_area < 0.01 * resabs0)
            goto RETURN;
    }
    
    {
        double ratio = res_ext / area;
    
        if (ratio < 0.01 || ratio > 100.0 || errsum > fabs (area))
            error_type = 6;
    }
    
    goto RETURN;

compute_result:

  *result = sum_results (workspace);
  *abserr = errsum;

RETURN:

    if (error_type > 2)
    error_type--;

    if (error_type == 0) 
    {
      status = SUCCESS;
    }
    else if (error_type == 1)
    {
      DR_Error ("number of iterations was insufficient");
    }
    else if (error_type == 2)
    {
      DR_Error ("cannot reach tolerance because of roundoff error");

    }
    else if (error_type == 3)
    {
      DR_Error ("bad integrand behavior found in the integration interval");

    }
    else if (error_type == 4)
    {
      DR_Error ("roundoff error detected in the extrapolation table");
    }
    else if (error_type == 5)
    {
      DR_Error ("integral is divergent, or slowly convergent");
    }
    else
    {
      DR_Error ("could not integrate function");
    }

    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}


/****************************************************************/
/* GENERAL INTEGRATION METHODS                                  */
/****************************************************************/
void gsl_integration_qk (   const int n, 
                            const double xgk[],
                            const double wg[],
                            const double wgk[],
                            double fv1[],
                            double fv2[],
                            PTR_FUNCTION f,
                            double a,
                            double b,
                            double *result,
                            double *abserr,
                            double *resabs,
                            double *resasc,
                            void *param)
{

    const double center = 0.5 * (a + b);
    const double half_length = 0.5 * (b - a);
    const double abs_half_length = fabs (half_length);
    double f_center;
    double result_gauss = 0;
    double result_kronrod;
    double result_abs;
    double result_asc = 0;
    double mean = 0, err = 0;
    int j;
    
    (*f)( center, &f_center, param);
    result_kronrod = f_center * wgk[n - 1];
    result_abs = fabs (result_kronrod);
    
    
    if (n % 2 == 0)
    {
        result_gauss = f_center * wg[n / 2 - 1];
    }

    for (j = 0; j < (n - 1) / 2; j++)
    {
        const int jtw = j * 2 + 1;	/* j=1,2,3 jtw=2,4,6 */
        const double abscissa = half_length * xgk[jtw];
        double fval1;
        double fval2;
        double fsum;
        (*f)( center - abscissa, &fval1, param);
        (*f)( center + abscissa, &fval2, param);
        fsum = fval1 + fval2;
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        result_gauss += wg[j] * fsum;
        result_kronrod += wgk[jtw] * fsum;
        result_abs += wgk[jtw] * (fabs (fval1) + fabs (fval2));
    }

    for (j = 0; j < n / 2; j++)
    {
        int jtwm1 = j * 2;
        const double abscissa = half_length * xgk[jtwm1];
        double fval1;
        double fval2;
        (*f)( center - abscissa, &fval1, param);
        (*f)( center + abscissa, &fval2, param);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        result_kronrod += wgk[jtwm1] * (fval1 + fval2);
        result_abs += wgk[jtwm1] * (fabs (fval1) + fabs (fval2));
    };

    mean = result_kronrod * 0.5;

    result_asc = wgk[n - 1] * fabs (f_center - mean);

    for (j = 0; j < n - 1; j++)
    {
      result_asc += wgk[j] * (fabs (fv1[j] - mean) + fabs (fv2[j] - mean));
    }

    /* scale by the width of the integration region */

    err = (result_kronrod - result_gauss) * half_length;

    result_kronrod *= half_length;
    result_abs *= abs_half_length;
    result_asc *= abs_half_length;

    *result = result_kronrod;
    *resabs = result_abs;
    *resasc = result_asc;
    *abserr = rescale_error (err, result_abs, result_asc);
}

/************************************/
/*  gsl_integration_qk15            */
/************************************/ 
/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */



void gsl_integration_qk15 ( const PTR_FUNCTION f,
                            double a,
                            double b,
                            double *result,
                            double *abserr,
                            double *resabs,
                            double *resasc,
                            void *param)
{
    const double xgk[8] =	/* abscissae of the 15-point kronrod rule */
    {
      0.991455371120812639206854697526329,
      0.949107912342758524526189684047851,
      0.864864423359769072789712788640926,
      0.741531185599394439863864773280788,
      0.586087235467691130294144838258730,
      0.405845151377397166906606412076961,
      0.207784955007898467600689403773245,
      0.000000000000000000000000000000000
    };

    /* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
       xgk[0], xgk[2], ... abscissae to optimally extend the 7-point gauss rule */

    const double wg[4] =	/* weights of the 7-point gauss rule */
    {
      0.129484966168869693270611432679082,
      0.279705391489276667901467771423780,
      0.381830050505118944950369775488975,
      0.417959183673469387755102040816327
    };

    const double wgk[8] =	/* weights of the 15-point kronrod rule */
    {
      0.022935322010529224963732008058970,
      0.063092092629978553290700663189204,
      0.104790010322250183839876322541518,
      0.140653259715525918745189590510238,
      0.169004726639267902826583426598550,
      0.190350578064785409913256402421014,
      0.204432940075298892414161999234649,
      0.209482141084727828012999174891714
    };

    double fv1[8], fv2[8];
    gsl_integration_qk (8, xgk, wg, wgk, fv1, fv2, f, a, b, result, abserr, resabs, resasc, param);
}

 
/************************************/
/*  gsl_integration_qk21            */
/************************************/ 
/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */

void gsl_integration_qk21 (const PTR_FUNCTION f, double a, double b,
		      double *result, double *abserr,
		      double *resabs, double *resasc,
              void *param)
{
    const double xgk[11] =	/* abscissae of the 21-point kronrod rule */
    {
      0.995657163025808080735527280689003,
      0.973906528517171720077964012084452,
      0.930157491355708226001207180059508,
      0.865063366688984510732096688423493,
      0.780817726586416897063717578345042,
      0.679409568299024406234327365114874,
      0.562757134668604683339000099272694,
      0.433395394129247190799265943165784,
      0.294392862701460198131126603103866,
      0.148874338981631210884826001129720,
      0.000000000000000000000000000000000
    };

    /* xgk[1], xgk[3], ... abscissae of the 10-point gauss rule. 
       xgk[0], xgk[2], ... abscissae to optimally extend the 10-point gauss rule */

    const double wg[5] =	/* weights of the 10-point gauss rule */
    {
      0.066671344308688137593568809893332,
      0.149451349150580593145776339657697,
      0.219086362515982043995534934228163,
      0.269266719309996355091226921569469,
      0.295524224714752870173892994651338
    };

    const double wgk[11] =	/* weights of the 21-point kronrod rule */
    {
      0.011694638867371874278064396062192,
      0.032558162307964727478818972459390,
      0.054755896574351996031381300244580,
      0.075039674810919952767043140916190,
      0.093125454583697605535065465083366,
      0.109387158802297641899210590325805,
      0.123491976262065851077958109831074,
      0.134709217311473325928054001771707,
      0.142775938577060080797094273138717,
      0.147739104901338491374841515972068,
      0.149445554002916905664936468389821
    };
    
    double fv1[11], fv2[11];
    gsl_integration_qk (11, xgk, wg, wgk, fv1, fv2, f, a, b, result, abserr, resabs, resasc, param);
}


/************************************/
/*  gsl_integration_qk31            */
/************************************/   
/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */


void gsl_integration_qk31 (const PTR_FUNCTION f, double a, double b,
      double *result, double *abserr,
      double *resabs, double *resasc,
      void *param)
{

    const double xgk[16] =	/* abscissae of the 31-point kronrod rule */
    {
      0.998002298693397060285172840152271,
      0.987992518020485428489565718586613,
      0.967739075679139134257347978784337,
      0.937273392400705904307758947710209,
      0.897264532344081900882509656454496,
      0.848206583410427216200648320774217,
      0.790418501442465932967649294817947,
      0.724417731360170047416186054613938,
      0.650996741297416970533735895313275,
      0.570972172608538847537226737253911,
      0.485081863640239680693655740232351,
      0.394151347077563369897207370981045,
      0.299180007153168812166780024266389,
      0.201194093997434522300628303394596,
      0.101142066918717499027074231447392,
      0.000000000000000000000000000000000
    };

    /* xgk[1], xgk[3], ... abscissae of the 15-point gauss rule. 
       xgk[0], xgk[2], ... abscissae to optimally extend the 15-point gauss rule */

    const double wg[8] =	/* weights of the 15-point gauss rule */
    {
      0.030753241996117268354628393577204,
      0.070366047488108124709267416450667,
      0.107159220467171935011869546685869,
      0.139570677926154314447804794511028,
      0.166269205816993933553200860481209,
      0.186161000015562211026800561866423,
      0.198431485327111576456118326443839,
      0.202578241925561272880620199967519
    };

    const double wgk[16] =	/* weights of the 31-point kronrod rule */
    {
      0.005377479872923348987792051430128,
      0.015007947329316122538374763075807,
      0.025460847326715320186874001019653,
      0.035346360791375846222037948478360,
      0.044589751324764876608227299373280,
      0.053481524690928087265343147239430,
      0.062009567800670640285139230960803,
      0.069854121318728258709520077099147,
      0.076849680757720378894432777482659,
      0.083080502823133021038289247286104,
      0.088564443056211770647275443693774,
      0.093126598170825321225486872747346,
      0.096642726983623678505179907627589,
      0.099173598721791959332393173484603,
      0.100769845523875595044946662617570,
      0.101330007014791549017374792767493
    };

    double fv1[16], fv2[16];
    gsl_integration_qk (16, xgk, wg, wgk, fv1, fv2, f, a, b, result, abserr, resabs, resasc, param);
}

/************************************/
/*  gsl_integration_qk41            */
/************************************/ 
/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */



void gsl_integration_qk41 (const PTR_FUNCTION f, double a, double b,
		      double *result, double *abserr,
		      double *resabs, double *resasc,
              void *param)
{

    const double xgk[21] =	/* abscissae of the 41-point kronrod rule */
    {
      0.998859031588277663838315576545863,
      0.993128599185094924786122388471320,
      0.981507877450250259193342994720217,
      0.963971927277913791267666131197277,
      0.940822633831754753519982722212443,
      0.912234428251325905867752441203298,
      0.878276811252281976077442995113078,
      0.839116971822218823394529061701521,
      0.795041428837551198350638833272788,
      0.746331906460150792614305070355642,
      0.693237656334751384805490711845932,
      0.636053680726515025452836696226286,
      0.575140446819710315342946036586425,
      0.510867001950827098004364050955251,
      0.443593175238725103199992213492640,
      0.373706088715419560672548177024927,
      0.301627868114913004320555356858592,
      0.227785851141645078080496195368575,
      0.152605465240922675505220241022678,
      0.076526521133497333754640409398838,
      0.000000000000000000000000000000000
    };

    /* xgk[1], xgk[3], ... abscissae of the 20-point gauss rule. 
       xgk[0], xgk[2], ... abscissae to optimally extend the 20-point gauss rule */

    const double wg[11] =	/* weights of the 20-point gauss rule */
    {
      0.017614007139152118311861962351853,
      0.040601429800386941331039952274932,
      0.062672048334109063569506535187042,
      0.083276741576704748724758143222046,
      0.101930119817240435036750135480350,
      0.118194531961518417312377377711382,
      0.131688638449176626898494499748163,
      0.142096109318382051329298325067165,
      0.149172986472603746787828737001969,
      0.152753387130725850698084331955098
    };

    const double wgk[21] =	/* weights of the 41-point kronrod rule */
    {
      0.003073583718520531501218293246031,
      0.008600269855642942198661787950102,
      0.014626169256971252983787960308868,
      0.020388373461266523598010231432755,
      0.025882133604951158834505067096153,
      0.031287306777032798958543119323801,
      0.036600169758200798030557240707211,
      0.041668873327973686263788305936895,
      0.046434821867497674720231880926108,
      0.050944573923728691932707670050345,
      0.055195105348285994744832372419777,
      0.059111400880639572374967220648594,
      0.062653237554781168025870122174255,
      0.065834597133618422111563556969398,
      0.068648672928521619345623411885368,
      0.071054423553444068305790361723210,
      0.073030690332786667495189417658913,
      0.074582875400499188986581418362488,
      0.075704497684556674659542775376617,
      0.076377867672080736705502835038061,
      0.076600711917999656445049901530102
    };

    double fv1[21], fv2[21];
    gsl_integration_qk (21, xgk, wg, wgk, fv1, fv2, f, a, b, result, abserr, resabs, resasc, param);
}


/************************************/
/*  gsl_integration_qk51            */
/************************************/ 
/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */



/* wgk[25] was calculated from the values of wgk[0..24] */

void gsl_integration_qk51 (const PTR_FUNCTION f, double a, double b,
		      double *result, double *abserr,
		      double *resabs, double *resasc,
              void *param)
{

    const double xgk[26] =	/* abscissae of the 51-point kronrod rule */
    {
      0.999262104992609834193457486540341,
      0.995556969790498097908784946893902,
      0.988035794534077247637331014577406,
      0.976663921459517511498315386479594,
      0.961614986425842512418130033660167,
      0.942974571228974339414011169658471,
      0.920747115281701561746346084546331,
      0.894991997878275368851042006782805,
      0.865847065293275595448996969588340,
      0.833442628760834001421021108693570,
      0.797873797998500059410410904994307,
      0.759259263037357630577282865204361,
      0.717766406813084388186654079773298,
      0.673566368473468364485120633247622,
      0.626810099010317412788122681624518,
      0.577662930241222967723689841612654,
      0.526325284334719182599623778158010,
      0.473002731445714960522182115009192,
      0.417885382193037748851814394594572,
      0.361172305809387837735821730127641,
      0.303089538931107830167478909980339,
      0.243866883720988432045190362797452,
      0.183718939421048892015969888759528,
      0.122864692610710396387359818808037,
      0.061544483005685078886546392366797,
      0.000000000000000000000000000000000
    };

    /* xgk[1], xgk[3], ... abscissae of the 25-point gauss rule. 
       xgk[0], xgk[2], ... abscissae to optimally extend the 25-point gauss rule */

    const double wg[13] =	/* weights of the 25-point gauss rule */
    {
      0.011393798501026287947902964113235,
      0.026354986615032137261901815295299,
      0.040939156701306312655623487711646,
      0.054904695975835191925936891540473,
      0.068038333812356917207187185656708,
      0.080140700335001018013234959669111,
      0.091028261982963649811497220702892,
      0.100535949067050644202206890392686,
      0.108519624474263653116093957050117,
      0.114858259145711648339325545869556,
      0.119455763535784772228178126512901,
      0.122242442990310041688959518945852,
      0.123176053726715451203902873079050
    };

    const double wgk[26] =	/* weights of the 51-point kronrod rule */
    {
      0.001987383892330315926507851882843,
      0.005561932135356713758040236901066,
      0.009473973386174151607207710523655,
      0.013236229195571674813656405846976,
      0.016847817709128298231516667536336,
      0.020435371145882835456568292235939,
      0.024009945606953216220092489164881,
      0.027475317587851737802948455517811,
      0.030792300167387488891109020215229,
      0.034002130274329337836748795229551,
      0.037116271483415543560330625367620,
      0.040083825504032382074839284467076,
      0.042872845020170049476895792439495,
      0.045502913049921788909870584752660,
      0.047982537138836713906392255756915,
      0.050277679080715671963325259433440,
      0.052362885806407475864366712137873,
      0.054251129888545490144543370459876,
      0.055950811220412317308240686382747,
      0.057437116361567832853582693939506,
      0.058689680022394207961974175856788,
      0.059720340324174059979099291932562,
      0.060539455376045862945360267517565,
      0.061128509717053048305859030416293,
      0.061471189871425316661544131965264,
      0.061580818067832935078759824240066
    };
    double fv1[26], fv2[26];
    gsl_integration_qk (26, xgk, wg, wgk, fv1, fv2, f, a, b, result, abserr, resabs, resasc, param);
}



/************************************/
/*  gsl_integration_qk61            */
/************************************/ 
/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */


void gsl_integration_qk61 (const PTR_FUNCTION f, double a, double b,
		      double *result, double *abserr,
		      double *resabs, double *resasc,
              void *param)
{

    const double xgk[31] =	/* abscissae of the 61-point kronrod rule */
    {
      0.999484410050490637571325895705811,
      0.996893484074649540271630050918695,
      0.991630996870404594858628366109486,
      0.983668123279747209970032581605663,
      0.973116322501126268374693868423707,
      0.960021864968307512216871025581798,
      0.944374444748559979415831324037439,
      0.926200047429274325879324277080474,
      0.905573307699907798546522558925958,
      0.882560535792052681543116462530226,
      0.857205233546061098958658510658944,
      0.829565762382768397442898119732502,
      0.799727835821839083013668942322683,
      0.767777432104826194917977340974503,
      0.733790062453226804726171131369528,
      0.697850494793315796932292388026640,
      0.660061064126626961370053668149271,
      0.620526182989242861140477556431189,
      0.579345235826361691756024932172540,
      0.536624148142019899264169793311073,
      0.492480467861778574993693061207709,
      0.447033769538089176780609900322854,
      0.400401254830394392535476211542661,
      0.352704725530878113471037207089374,
      0.304073202273625077372677107199257,
      0.254636926167889846439805129817805,
      0.204525116682309891438957671002025,
      0.153869913608583546963794672743256,
      0.102806937966737030147096751318001,
      0.051471842555317695833025213166723,
      0.000000000000000000000000000000000
    };

    /* xgk[1], xgk[3], ... abscissae of the 30-point gauss rule. 
       xgk[0], xgk[2], ... abscissae to optimally extend the 30-point gauss rule */

    const double wg[15] =	/* weights of the 30-point gauss rule */
    {
      0.007968192496166605615465883474674,
      0.018466468311090959142302131912047,
      0.028784707883323369349719179611292,
      0.038799192569627049596801936446348,
      0.048402672830594052902938140422808,
      0.057493156217619066481721689402056,
      0.065974229882180495128128515115962,
      0.073755974737705206268243850022191,
      0.080755895229420215354694938460530,
      0.086899787201082979802387530715126,
      0.092122522237786128717632707087619,
      0.096368737174644259639468626351810,
      0.099593420586795267062780282103569,
      0.101762389748405504596428952168554,
      0.102852652893558840341285636705415
    };

    const double wgk[31] =	/* weights of the 61-point kronrod rule */
    {
      0.001389013698677007624551591226760,
      0.003890461127099884051267201844516,
      0.006630703915931292173319826369750,
      0.009273279659517763428441146892024,
      0.011823015253496341742232898853251,
      0.014369729507045804812451432443580,
      0.016920889189053272627572289420322,
      0.019414141193942381173408951050128,
      0.021828035821609192297167485738339,
      0.024191162078080601365686370725232,
      0.026509954882333101610601709335075,
      0.028754048765041292843978785354334,
      0.030907257562387762472884252943092,
      0.032981447057483726031814191016854,
      0.034979338028060024137499670731468,
      0.036882364651821229223911065617136,
      0.038678945624727592950348651532281,
      0.040374538951535959111995279752468,
      0.041969810215164246147147541285970,
      0.043452539701356069316831728117073,
      0.044814800133162663192355551616723,
      0.046059238271006988116271735559374,
      0.047185546569299153945261478181099,
      0.048185861757087129140779492298305,
      0.049055434555029778887528165367238,
      0.049795683427074206357811569379942,
      0.050405921402782346840893085653585,
      0.050881795898749606492297473049805,
      0.051221547849258772170656282604944,
      0.051426128537459025933862879215781,
      0.051494729429451567558340433647099
    };

    double fv1[31], fv2[31];
    gsl_integration_qk (31, xgk, wg, wgk, fv1, fv2, f, a, b, result, abserr, resabs, resasc, param);
}

