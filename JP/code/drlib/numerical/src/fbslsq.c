#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static Mf_spline *pp = NULL;
static VA_LIST_HACK PROTO (l_spline_least_squares, (Mint ndata, Mfloat xdata[],
	            Mfloat fdata[], Mint spline_space_dim,
	            va_list argptr));
#ifdef ANSI
    Mf_spline  *imsl_f_spline_least_squares (Mint ndata, Mfloat xdata[], Mfloat fdata[],
                Mint spline_space_dim,...)
#else
    Mf_spline  *imsl_f_spline_least_squares (ndata, xdata, fdata, spline_space_dim, va_alist)
    Mint        ndata;
    Mfloat      xdata[];
    Mfloat      fdata[];
    Mint        spline_space_dim;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, spline_space_dim);
#ifdef DOUBLE
    imsl_e1psh ("imsl_d_spline_least_squares");
#else
    imsl_e1psh ("imsl_f_spline_least_squares");
#endif
    IMSL_CALL (l_spline_least_squares (ndata, xdata, fdata, spline_space_dim, argptr));
    va_end (argptr);
#ifdef DOUBLE
    imsl_e1pop ("imsl_d_spline_least_squares");
#else
    imsl_e1pop ("imsl_f_spline_least_squares");
#endif
    return pp;
}
#ifdef ANSI
static VA_LIST_HACK l_spline_least_squares (Mint ndata, Mfloat xdata[], Mfloat fdata[],
                Mint spline_space_dim, va_list argptr)
#else
static VA_LIST_HACK l_spline_least_squares (ndata, xdata, fdata, spline_space_dim, argptr)
    Mint        ndata;
    Mfloat      xdata[];
    Mfloat      fdata[];
    Mint        spline_space_dim;
    va_list     argptr;
#endif
{
    Mint        code;
    Mint        users_xorder;
    Mint        i;
    Mint        domain_dim;
    Mint        target_dim;
    Mint        free_the_structure = 0;
    Mint        arg_number = 4;
    Mint        four = 4;
    Mfloat     *users_knots = NULL;
    Mfloat     *weights = NULL;
    Mint       *num_coefs = NULL;
    Mint       *order = NULL;
    Mint        order_given = 0;
    Mint        knots_given = 0;
    Mint        weights_given = 0;
    Mint        sse_wanted = 0;
    Mint        use_bsvls = 0;	/* THE DEFAULT IS TO USE BSLSQ */
    Mint        derivative_x = 0;
    Mfloat     *sse_return = NULL;
    Mfloat     *wk1_b2lsq = NULL;
    Mfloat     *wk2_b2lsq = NULL;
    Mfloat     *wk3_b2lsq = NULL;
    Mfloat     *wk4_b2lsq = NULL;
    Mint       *iwk1_b2lsq = NULL;
    Mfloat     *wk1_b2vls = NULL;
    Mint       *iwk1_b2vls = NULL;
    Mfloat     *xguess = NULL;
    Mfloat      sse;
    Mfloat      temp_float;
    Mint        temp_int;
    Mfloat      x_small;
    Mfloat      x_big;
    Mfloat      interval_size;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_ORDER:
	    order_given = 1;
	    users_xorder = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_KNOTS:
	    knots_given = 1;
	    users_knots = va_arg (argptr, Mfloat *);
	    arg_number++;
	    break;
	case IMSL_WEIGHTS:
	    weights_given = 1;
	    weights = va_arg (argptr, Mfloat *);
	    arg_number++;
	    break;
	case IMSL_SSE:
	    sse_wanted = 1;
	    sse_return = va_arg (argptr, Mfloat *);
	    arg_number++;
	    break;
	case IMSL_OPTIMIZE:
	    use_bsvls = 1;
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
    /* SET PARAMETERS */
    /* DOMAIN & RANGE DIMENSIONS */
    domain_dim = 1;
    target_dim = 1;
    /* ORDERS */
    if (order_given == 0) {
	order = &four;
    }
    else {
	order = &users_xorder;
    }
    /* NUM_COEFS */
    num_coefs = &spline_space_dim;
    /* CHECK KORDER */
    if (*order < 1) {
	imsl_e1sti (1, *order);
	imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
	goto RETURN;
    }
    /* CHECK NCOEF */
    if (*num_coefs < *order) {
	imsl_e1sti (1, *num_coefs);
	imsl_e1sti (2, *order);
	imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS);
	goto RETURN;
    }
    /* CHECK ARGUMENT NCOEF */
    if (*num_coefs > ndata) {
	imsl_e1sti (1, *num_coefs);
	imsl_e1sti (2, ndata);
	imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS_2);
	goto RETURN;
    }

    /* WEIGHTS */
    if (weights_given == 0) {
	weights = (Mfloat *) imsl_malloc (ndata * sizeof (Mfloat));
	if (weights == NULL) {
	    imsl_e1stl (1, "ndata");
	    imsl_e1sti (1, ndata);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto RETURN;
	}
	for (i = 0; i < ndata; i++) {
	    (weights[i]) = F_ONE;
	}
    }
    /* CREATE THE STRUCTURE */
    if (knots_given == 0)
	pp = imsl_f_spline_create (domain_dim, target_dim, order, num_coefs, 0);
    else
	pp = imsl_f_spline_create (domain_dim, target_dim, order, num_coefs, IMSL_KNOTS, &users_knots, 0);
    if (imsl_n1rty (1) == 4) {
	imsl_e1mes (0, 0, " ");
	free_the_structure = 1;
	imsl_e1stl (1, "ndata");
	imsl_e1sti (1, ndata);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto RETURN;
    }
    /*
     * COMPUTE THE KNOTS IF THE USER DID NOT SUPPLY THEM.  THESE KNOTS ARE
     * EQUALLY SPACED IN THE INTERVAL AND STACKED APPROPRIATELY AT THE
     * ENDPOINTS.
     */

    if (knots_given == 0) {
	x_big = xdata[ndata-1];
	x_small = xdata[0];
	for (i = 1; i < ndata; i++) {
	    if (xdata[i] < x_small) {
		x_small = xdata[i];
	    }
	    else if (xdata[i] > x_big) {
		x_big = xdata[i];
	    }
	}
	interval_size = fabs (x_big - x_small);

	for (i = 1; i <= (*num_coefs - *order + 2); i++) {
	    pp->knots[0][i + *order - 2] = x_small + interval_size * ((Mfloat) (i - 1) / (Mfloat) (*num_coefs -
		    *order + 1));
	}
	pp->knots[0][*num_coefs] += 0.001;
	/* Stack knots. */
	for (i = 1; i <= (*order - 1); i++) {
	    pp->knots[0][i - 1] = pp->knots[0][*order - 1];
	    pp->knots[0][i + *num_coefs] = pp->knots[0][*num_coefs];
	}
    }
    /*
     * IF THE USER SUPPLIED THE OPTIONAL ARGUMENT IMSL_OPTIMIZE, WE USE
     * BSVLS.  OTHERWISE, WE USE BSLSQ
     */

    if (use_bsvls == 1) {
	xguess = (Mfloat *) imsl_malloc ((spline_space_dim + *order) * sizeof (*xguess));
	if ((xguess == NULL)) {
	    free_the_structure = 1;
	    imsl_e1stl (1, "order");
	    imsl_e1sti (1, *order);
	    imsl_e1stl (2, "spline_space_dim");
	    imsl_e1sti (2, spline_space_dim);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	    goto B2VLS_FREE_SPACE;
	}
	for (i = 0; i < *order + spline_space_dim; i++) {
	    xguess[i] = pp->knots[0][i];
	}
	temp_int = spline_space_dim * (6 + 2 ** order) + (*order) * (7 - (*order)) + 3 * ndata + 3;
	wk1_b2vls = (Mfloat *) imsl_malloc (temp_int * sizeof (*wk1_b2vls));
	iwk1_b2vls = (Mint *) imsl_malloc (ndata * sizeof (*iwk1_b2vls));
	if ((wk1_b2vls == NULL) || (iwk1_b2vls == NULL)) {
	    free_the_structure = 1;
	    imsl_e1stl (1, "ndata");
	    imsl_e1sti (1, ndata);
	    imsl_e1stl (2, "order");
	    imsl_e1sti (2, *order);
	    imsl_e1stl (3, "spline_space_dim");
	    imsl_e1sti (3, spline_space_dim);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_3);
	    goto B2VLS_FREE_SPACE;
	}

	imsl_b2vls (&ndata, xdata, fdata, weights, pp->order, pp->num_coef, xguess, pp->knots[0],
	    pp->coef[0], &sse, iwk1_b2vls, wk1_b2vls);
	if (imsl_n1rty (1) > 3) {
	    free_the_structure = 1;
	    goto B2VLS_FREE_SPACE;
	}
	if (sse_wanted == 1) {
	    *sse_return = (sse) * (sse);
	}
	goto B2VLS_FREE_SPACE;
    }
    /*
     * THIS ELSE CLAUSE WILL BE EXECUTED IN ALL CASES EXCEPT WHEN
     * IMSL_OPTIMIZE IS SPECIFIED.
     */
    else {
	/* GET THE WORKSPACE NEEDED IN B2LSQ    */
	wk1_b2lsq = (Mfloat *) imsl_malloc (((3 + (*num_coefs)) * (*order)) * sizeof (*wk1_b2lsq));
	wk2_b2lsq = (Mfloat *) imsl_malloc (ndata * sizeof (*wk2_b2lsq));
	wk3_b2lsq = (Mfloat *) imsl_malloc (ndata * sizeof (*wk3_b2lsq));
	wk4_b2lsq = (Mfloat *) imsl_malloc (ndata * sizeof (*wk4_b2lsq));
	iwk1_b2lsq = (Mint *) imsl_malloc (ndata * sizeof (*iwk1_b2lsq));
	if ((wk1_b2lsq == NULL) || (wk2_b2lsq == NULL) || (wk3_b2lsq == NULL) || (wk4_b2lsq == NULL) || (iwk1_b2lsq == NULL)) {
	    free_the_structure = 1;
	    imsl_e1stl (1, "ndata");
	    imsl_e1sti (1, ndata);
	    imsl_e1stl (3, "order");
	    imsl_e1sti (3, *order);
	    imsl_e1stl (5, "spline_space_dim");
	    imsl_e1sti (5, spline_space_dim);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_3);
	    goto B2LSQ_FREE_SPACE;
	}
	/* CALL THE SPLINE ROUTINE */
	imsl_b2lsq (&ndata, xdata, fdata, weights, (pp->order),
	    pp->knots[0], (pp->num_coef), pp->coef[0],
	    wk1_b2lsq, wk2_b2lsq, wk3_b2lsq, wk4_b2lsq, iwk1_b2lsq);
	/*
	 * IF THE SPLINE COULD NOT BE COMPUTED, THEN FREE THE STRUCTURE AND
	 * RETURN NULL
	 */
	if (imsl_n1rty (1) > 3) {
	    free_the_structure = 1;
	    goto B2LSQ_FREE_SPACE;
	}
	/*
	 * IF THE USER WANTS THE SUM OF THE SQUARES OF THE ERRORS AT THE
	 * ORIGINAL DATA POINTS, COMPUTE IT
	 */
	if (sse_wanted == 1) {
	    sse = F_ZERO;
	    for (i = 0; i < ndata; i++) {
		temp_float = imsl_b2der (&derivative_x, &(xdata[i]), (pp->order), pp->knots[0],
		    (pp->num_coef), (pp->coef[0]), wk1_b2lsq, wk2_b2lsq, wk3_b2lsq);
		sse += (temp_float - fdata[i]) * (temp_float - fdata[i]);
	    }
	    *sse_return = (sse);
	}
    }
    /* FREE THE WORKSPACE USED */
B2LSQ_FREE_SPACE:
    if (wk1_b2lsq != NULL)
	imsl_free (wk1_b2lsq);
    if (wk2_b2lsq != NULL)
	imsl_free (wk2_b2lsq);
    if (wk3_b2lsq != NULL)
	imsl_free (wk3_b2lsq);
    if (wk4_b2lsq != NULL)
	imsl_free (wk4_b2lsq);
    if (iwk1_b2lsq != NULL)
	imsl_free (iwk1_b2lsq);
    if (free_the_structure == 1) {
	if (pp != NULL)
	    imsl_free (pp);
	pp = NULL;
    }
B2VLS_FREE_SPACE:
    if ((weights_given == 0)&&(weights != NULL))
	imsl_free (weights);
    if ((use_bsvls == 1)&&(xguess != NULL))
	imsl_free (xguess);
    if (wk1_b2vls != NULL)
	imsl_free (wk1_b2vls);
    if (iwk1_b2vls != NULL)
	imsl_free (iwk1_b2vls);
    if (free_the_structure == 1) {
	if (pp != NULL)
	    imsl_free (pp);
	pp = NULL;
    }

RETURN:
    return argptr;
}
