#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_user_fcn_least_squares (Mfloat (*fcn) (Mint, Mfloat), Mint, Mint, Mfloat[], Mfloat[], va_list);
#else
static VA_LIST_HACK l_user_fcn_least_squares ();
#endif

static Mfloat *lv_coefficients = NULL;
#ifdef ANSI
Mfloat     *imsl_f_user_fcn_least_squares (Mfloat (*fcn) (Mint, Mfloat), Mint nbasis, Mint ndata, Mfloat xdata[], Mfloat ydata[],...)
#else
Mfloat     *imsl_f_user_fcn_least_squares (fcn, nbasis, ndata, xdata, ydata, va_alist)
    Mfloat      (*fcn) ();
    Mint        nbasis;
    Mint        ndata;
    Mfloat      xdata[], ydata[];
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, ydata);
#ifdef DOUBLE
    imsl_e1psh ("imsl_d_user_fcn_least_squares");
#else
    imsl_e1psh ("imsl_f_user_fcn_least_squares");
#endif
    lv_coefficients = NULL;
    IMSL_CALL (l_user_fcn_least_squares (fcn, nbasis, ndata, xdata, ydata, argptr));
    va_end (argptr);
#ifdef DOUBLE
    imsl_e1pop ("imsl_d_user_fcn_least_squares");
#else
    imsl_e1pop ("imsl_f_user_fcn_least_squares");
#endif
    return lv_coefficients;
}




#ifdef ANSI
static VA_LIST_HACK l_user_fcn_least_squares (Mfloat (*fcn) (Mint, Mfloat), Mint nbasis, Mint ndata, Mfloat xdata[], Mfloat ydata[], va_list argptr)
#else
static VA_LIST_HACK l_user_fcn_least_squares (fcn, nbasis, ndata, xdata, ydata, argptr)
    Mfloat      (*fcn) ();
    Mint        nbasis;
    Mint        ndata;
    Mfloat      xdata[], ydata[];
    va_list     argptr;
#endif
{

    Mint        code = 1;
    Mint        arg_number = 5;
    Mfloat     *intercept = NULL;
    Mfloat     *weights = NULL;
    Mfloat     *ssq = NULL;
    Mfloat      local_ssq = F_ZERO;
    Mint        free_space = 0;
    Mint        users_space = 0;
    Mint        iwt = 0;
    Mint        use_intercept = 0;
    Mint        intcep = 0;
    Mint        wants_ssq = 0;
    Mfloat     *coefficients = NULL;
    Mfloat     *a = NULL;
    Mfloat     *wk = NULL;
    Mint        temp_int;
    Mint	i;
    while (code > 0) {
	code = va_arg (argptr, int);
	++arg_number;
	switch (code) {
	case IMSL_RETURN_USER:
	    coefficients = va_arg (argptr, Mfloat *);
	    ++arg_number;
	    users_space = 1;
	    free_space = 1;
	    break;
	case IMSL_WEIGHTS:
	    weights = va_arg (argptr, Mfloat *);
	    ++arg_number;
	    iwt = 1;
	    break;
	case IMSL_INTERCEPT:
	    intercept = va_arg (argptr, Mfloat *);
	    ++arg_number;
	    use_intercept = 1;
	    intcep = 1;
	    break;
	case IMSL_SSE:
	    ssq = va_arg (argptr, Mfloat *);
	    ++arg_number;
	    wants_ssq = 1;
	    break;
	case 0:
	    break;
	default:
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    goto RETURN;
	}
    }
    /* Check NBASIS */
    if (nbasis < 1) {
	imsl_e1sti (1, nbasis);
	imsl_ermes (IMSL_TERMINAL, IMSL_NBASIS_FCNS_TOO_SMALL);
    }
    if (imsl_n1rty (0) > 0)
	goto RETURN;
    /* GET SPACE FOR THE COEFFICIENTS AND WORKSPACE */
    a = (Mfloat *) imsl_malloc ((nbasis + 1) * sizeof (*a));
    temp_int = (intcep + nbasis) * (intcep + nbasis) + 4 * (intcep + nbasis) + iwt + 1;
    wk = (Mfloat *) imsl_malloc ((temp_int) * sizeof (*wk));
    if ((a == NULL) || (wk == NULL)) {
	free_space = 1;
	imsl_e1stl (1, "nbasis");
	imsl_e1sti (1, nbasis);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }
    /* CALL THE ROUTINE */
    if (wants_ssq)
	imsl_f2lsq (fcn, &intcep, &nbasis, &ndata, xdata, ydata, &iwt,
	    weights, a, ssq, wk);
    else
	imsl_f2lsq (fcn, &intcep, &nbasis, &ndata, xdata, ydata, &iwt,
	    weights, a, &local_ssq, wk);
    if (imsl_n1rty (1) > 3) {
	free_space = 1;
	goto FREE_SPACE;
    }

    if (users_space == 1) {
	if (intcep == 1) {
	    *intercept = a[0];
	    scopy (nbasis, &a[1], 1, coefficients, 1);
	    lv_coefficients = coefficients;
	}
	else
	    scopy (nbasis, &a[0], 1, coefficients, 1);
	lv_coefficients = coefficients;
    }
    else {
	if (intcep == 1) {
	    *intercept = a[0];
/* This assignment was replaced with the loop and assignment because
   lv_coefficients = &(a[1]) means that we can not legally free
   lv_coefficients (it does not point the the start of an allocated
   block).  We overwrite *a, but it was not returned anyway.

	    lv_coefficients = &(a[1]);
*/
	    for (i=0; i<nbasis; i++)
		a[i] = a[i+1];
	    lv_coefficients = a;
	}
	else
	    lv_coefficients = &(a[0]);
    }
FREE_SPACE:
    if ((free_space == 1) && (a != NULL))
	imsl_free (a);
    if (wk != NULL)
	imsl_free (wk);
RETURN:
    return (argptr);
}
