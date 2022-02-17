#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static VA_LIST_HACK PROTO (l_spline_knots, (Mint ndata, Mfloat xdata[], va_list argptr));

static Mfloat *knots = NULL;
#ifdef ANSI
Mfloat     *imsl_f_spline_knots (Mint ndata, Mfloat xdata[],...)
#else
Mfloat     *imsl_f_spline_knots (ndata, xdata, va_alist)
    Mint        ndata;
    Mfloat      xdata[];
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, xdata);
#ifdef DOUBLE
    imsl_e1psh ("imsl_d_spline_knots");
#else
    imsl_e1psh ("imsl_f_spline_knots");
#endif
    IMSL_CALL (l_spline_knots (ndata, xdata, argptr));
    va_end (argptr);
#ifdef DOUBLE
    imsl_e1pop ("imsl_d_spline_knots");
#else
    imsl_e1pop ("imsl_f_spline_knots");
#endif
    return (knots);

}
#ifdef ANSI
static VA_LIST_HACK l_spline_knots (Mint ndata, Mfloat xdata[], va_list argptr)
#else
static VA_LIST_HACK l_spline_knots (ndata, xdata, argptr)
    Mint        ndata;
    Mfloat      xdata[];
    va_list     argptr;
#endif
{
    Mint        arg_number = 2;
    Mint        code;
    Mint        temp_int;
    Mint        order = 4;
    Mint        wants_b2opk = 0;
    Mint        use_users_knot_space = 0;
    Mint        free_knots = 0;
    Mfloat     *users_knot_space;
    Mint       *iwk_b2nak = NULL;
    Mint        maxit = 10;
    Mfloat     *wk_b2opk = NULL;
    Mint       *iwk_b2opk = NULL;
    Mfloat     *wk_b2nak = NULL;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, int);
	arg_number++;
	switch (code) {
	case IMSL_ORDER:
	    order = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_OPT:
	    wants_b2opk = 1;
	    break;
	case IMSL_MIN_PROJECTION:
	    /* NOT YET IMPLEMENTED */
	    break;
	case IMSL_RETURN_USER:
	    use_users_knot_space = 1;
	    users_knot_space = va_arg (argptr, Mfloat *);
	    arg_number++;
	    break;
	case IMSL_OPT_ITMAX:
	    maxit = va_arg (argptr, Mint);
	    arg_number++;
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
    if (wants_b2opk == 1) {
	/* Check KORDER */
	if (order < 3) {
	    imsl_e1sti (1, order);

/*              imsl_ermes(5, 1, "The order of the spline must be at least 3 while KORDER = %(i1) is given.");
*/
	    imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER_2);
	}
	/* Check NDATA */
	if (ndata < order) {
	    imsl_e1sti (1, ndata);
	    imsl_e1sti (2, order);

/*              imsl_ermes(5, 2, "The number of data points must be at least as large as the order of the spline while NDATA = %(i1) and KORDER = %(i2) are given.");
*/
	    imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_NEED_DATA_PTS);
	}
	if (imsl_n1rty (0) != 0)
	    goto RETURN;

	temp_int = (ndata - order) * (3 * order - 2) + 6 * ndata + 2 * order + 5;
	wk_b2opk = (Mfloat *) imsl_malloc (temp_int * sizeof (*wk_b2opk));
	iwk_b2opk = (Mint *) imsl_malloc (ndata * sizeof (*iwk_b2opk));
	if ((wk_b2opk == NULL) || (iwk_b2opk == NULL)) {
	    imsl_e1stl (1, "ndata");
	    imsl_e1sti (1, ndata);
	    imsl_e1stl (2, "order");
	    imsl_e1sti (2, order);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	    goto FREE_SPACE_B2OPK;
	}

	if (use_users_knot_space == 1) {
	    imsl_b2opk (&ndata, xdata, &order, users_knot_space, &maxit, wk_b2opk, iwk_b2opk);
	    if (imsl_n1rty (1) < 4)
		knots = users_knot_space;
	}
	else {
	    knots = (Mfloat *) imsl_malloc ((ndata + order) * sizeof (*knots));
	    if (knots == NULL) {
		imsl_e1stl (1, "ndata");
		imsl_e1sti (1, ndata);
		imsl_e1stl (2, "order");
		imsl_e1sti (2, order);
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
		goto FREE_SPACE_B2OPK;
	    }

	    imsl_b2opk (&ndata, xdata, &order, knots, &maxit, wk_b2opk, iwk_b2opk);
	    if (imsl_n1rty (1) > 3) {
		free_knots = 1;
	    }

	}
	goto FREE_SPACE_B2OPK;
    }


    else {

	if (order <= 1) {
	    imsl_e1sti (1, order);

/*              imsl_ermes(5, 1, "The order of the spline must be at least 2 while KORDER = %(i1) is given.");
*/
	    imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER_1);
	}
	/* Check NDATA */
	if (ndata < order) {
	    imsl_e1sti (1, ndata);
	    imsl_e1sti (2, order);

/*              imsl_ermes(5, 2, "The number of data points must be at least as large as the order of the spline while NDATA = %(i1) and KORDER = %(i2) are given.");
*/
	    imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_NEED_DATA_PTS);
	}
	if (imsl_n1rty (0) != 0)
	    goto RETURN;


	wk_b2nak = (Mfloat *) imsl_malloc (ndata * sizeof (*wk_b2nak));
	iwk_b2nak = (Mint *) imsl_malloc (ndata * sizeof (*iwk_b2nak));

	if ((wk_b2nak == NULL) || (iwk_b2nak == NULL)) {
	    imsl_e1stl (1, "ndata");
	    imsl_e1sti (1, ndata);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE_B2NAK;
	}

	if (use_users_knot_space == 1) {
	    imsl_b2nak (&ndata, xdata, &order, users_knot_space, wk_b2nak, iwk_b2nak);
	    if (imsl_n1rty (1) < 4)
		knots = users_knot_space;
	}
	else {
	    knots = (Mfloat *) imsl_malloc ((ndata + order) * sizeof (*knots));
	    if (knots == NULL) {
		imsl_e1stl (1, "ndata");
		imsl_e1sti (1, ndata);
		imsl_e1stl (2, "order");
		imsl_e1sti (2, order);
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
		goto FREE_SPACE_B2NAK;
	    }
	    imsl_b2nak (&ndata, xdata, &order, knots, wk_b2nak, iwk_b2nak);
	    if (imsl_n1rty (1) > 3) {
		free_knots = 1;
	    }
	}
	goto FREE_SPACE_B2NAK;
    }
FREE_SPACE_B2NAK:
    if ((knots != NULL) && (free_knots == 1)) {
	imsl_free (knots);
	knots = NULL;
    }
    if (iwk_b2nak != NULL)
	imsl_free (iwk_b2nak);
    if (wk_b2nak != NULL)
	imsl_free (wk_b2nak);
    goto RETURN;

FREE_SPACE_B2OPK:
    if ((knots != NULL) && (free_knots == 1)) {
	imsl_free (knots);
	knots = NULL;
    }
    if (iwk_b2opk != NULL)
	imsl_free (iwk_b2opk);
    if (wk_b2opk != NULL)
	imsl_free (wk_b2opk);
    goto RETURN;

RETURN:
    return (argptr);
}
