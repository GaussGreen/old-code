#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif


static VA_LIST_HACK PROTO (l_gauss_quad_rule, (Mint n, Mfloat weights[],
	            Mfloat points[], va_list argptr));
#ifdef ANSI
    void        imsl_f_gauss_quad_rule (Mint n, Mfloat weights[],
                Mfloat points[],...)
#else
    void        imsl_f_gauss_quad_rule (n, weights, points, va_alist)
    Mint        n;
    Mfloat      weights[];
    Mfloat      points[];
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, points);
    E1PSH ("imsl_f_gauss_quad_rule", "imsl_d_gauss_quad_rule");
    IMSL_CALL (l_gauss_quad_rule (n, weights, points, argptr));
    va_end (argptr);
    E1POP ("imsl_f_gauss_quad_rule", "imsl_d_gauss_quad_rule");
    return;
}


#ifdef ANSI
static VA_LIST_HACK l_gauss_quad_rule (Mint n, Mfloat weights[],
                Mfloat points[], va_list argptr)
#else
static VA_LIST_HACK l_gauss_quad_rule (n, weights, points, argptr)
    Mint        n;
    Mfloat      weights[];
    Mfloat      points[];
    va_list     argptr;
#endif
{
    Mint        arg_number = 3;
    Mint        code;
    Mint        iweigh = 1;	/* DEFAULT VALUE */
    Mint        nfix = 0;	/* DEFAULT VALUE */
    Mfloat      fixed_points[2];
    Mfloat      alpha = F_ZERO;
    Mfloat      beta = F_ZERO;
    Mfloat     *wk = NULL;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, int);
	arg_number++;
	switch (code) {
	case IMSL_LEGENDRE:
	    iweigh = 1;
	    break;
	case IMSL_CHEBYSHEV_FIRST:
	    iweigh = 2;
	    break;
	case IMSL_CHEBYSHEV_SECOND:
	    iweigh = 3;
	    break;
	case IMSL_HERMITE:
	    iweigh = 4;
	    break;
	case IMSL_COSH:
	    iweigh = 7;
	    break;

	case IMSL_JACOBI:
	    iweigh = 5;
	    alpha = (Mfloat) va_arg (argptr, double);
	    beta = (Mfloat) va_arg (argptr, double);
	    break;
	case IMSL_JACOBI_ADR:
	    iweigh = 5;
	    alpha = *( va_arg (argptr, Mfloat *));
	    beta =  *( va_arg (argptr, Mfloat *));
	    break;


	case IMSL_GEN_LAGUERRE:
	    iweigh = 6;
	    alpha = (Mfloat) va_arg (argptr, double);
	    break;
	case IMSL_GEN_LAGUERRE_ADR:
	    iweigh = 6;
	    alpha = *( va_arg (argptr, Mfloat *));
	    break;


	case IMSL_FIXED_POINT:
	    nfix = 1;
	    fixed_points[0] = (Mfloat) va_arg (argptr, double);
	    break;
	case IMSL_FIXED_POINT_ADR:
	    nfix = 1;
	    fixed_points[0] = *( va_arg (argptr, Mfloat *));
	    break;


	case IMSL_TWO_FIXED_POINTS:
	    nfix = 2;
	    fixed_points[0] = (Mfloat) va_arg (argptr, double);
	    fixed_points[1] = (Mfloat) va_arg (argptr, double);
	    break;
	case IMSL_TWO_FIXED_POINTS_ADR:
	    nfix = 2;
	    fixed_points[0] = *( va_arg (argptr, Mfloat *));
	    fixed_points[1] = *( va_arg (argptr, Mfloat *));
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

    if (n < 1) {
	imsl_e1sti (1, n);
	/*
	 * (5, 1, "The number of quadrature points n = %(i1).  N must be at
	 * least 1.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NUM_QUADRATURE_POINTS);
	goto RETURN;
    }

    if (weights == NULL) {
	imsl_e1stl (1, "weights");
	/* (5, 1, "The required argument %(L1) is NULL."); */
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
    }
    if (points == NULL) {
	imsl_e1stl (1, "points");
	/* (5, 1, "The required argument %(L1) is NULL."); */
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
    }
    if (imsl_n1rty (0))
	goto RETURN;


/* GET THE WORKSPACE  */
    wk = (Mfloat *) imsl_malloc (n * sizeof (*wk));

/* CALL TO ROUTINE    */
    imsl_g2rul (&n, &iweigh, &alpha, &beta, &nfix, fixed_points, points, weights, wk);

/* FREE THE SPACE     */
    if (wk != NULL)
	imsl_free (wk);
RETURN:
    return argptr;
}
