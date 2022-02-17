#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static Mf_spline *pp = NULL;

static VA_LIST_HACK PROTO (l_spline_interp, (Mint ndata, Mfloat xdata[], Mfloat fdata[], va_list argptr));
#ifdef ANSI
Mf_spline  *imsl_f_spline_interp (Mint ndata, Mfloat xdata[], Mfloat fdata[],...)
#else
Mf_spline  *imsl_f_spline_interp (ndata, xdata, fdata, va_alist)
    Mint        ndata;
    Mfloat      xdata[];
    Mfloat      fdata[];
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, fdata);
#ifdef DOUBLE
    imsl_e1psh ("imsl_d_spline_interp");
#else
    imsl_e1psh ("imsl_f_spline_interp");
#endif
    IMSL_CALL (l_spline_interp (ndata, xdata, fdata, argptr));
    va_end (argptr);
#ifdef DOUBLE
    imsl_e1pop ("imsl_d_spline_interp");
#else
    imsl_e1pop ("imsl_f_spline_interp");
#endif
    return pp;
}
#ifdef ANSI
static VA_LIST_HACK l_spline_interp (Mint ndata, Mfloat xdata[], Mfloat fdata[], va_list argptr)
#else
static VA_LIST_HACK l_spline_interp (ndata, xdata, fdata, argptr)
    Mint        ndata;
    Mfloat      xdata[];
    Mfloat      fdata[];
    va_list     argptr;
#endif
{
    Mint        arg_number = 3;
    Mint        four = 4;
    Mint        code;
    Mint        domain_dim;
    Mint        target_dim;
    Mint        order_given = 0;
    Mint        knots_given = 0;
    Mfloat     *users_knots;
    Mfloat     *xsrt;
    Mfloat     *wk1;
    Mfloat     *wk2;
    Mfloat     *wk3;
    Mint       *iwk_b2nak;
    Mint       *iwk_b2int;
    Mint        users_order;
    Mint       *order;
    Mint       *num_coefs;
    Mint        tempint;
    Mint        free_the_structure = 0;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, int);
	arg_number++;
	switch (code) {
	case IMSL_ORDER:
	    order_given = 1;
	    users_order = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_KNOTS:
	    knots_given = 1;
	    users_knots = va_arg (argptr, Mfloat *);
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

    domain_dim = 1;
    target_dim = 1;

    if (order_given == 0)
	order = &four;
    else
	order = &users_order;
    tempint = ndata;
    num_coefs = &tempint;
    if (*order < 1) {
	imsl_e1sti (1, *order);

	imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
	goto RETURN;
    }
    /* CHECK NDATA */
    if (ndata < *order) {
	imsl_e1sti (1, ndata);
	imsl_e1sti (2, *order);

	imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_NEED_DATA_PTS);
	goto RETURN;
    }
    /* CREATE THE STRUCTURE */
    if (knots_given == 0)
	pp = imsl_f_spline_create (domain_dim, target_dim, order, num_coefs, 0);
    else
	pp = imsl_f_spline_create (domain_dim, target_dim, order, num_coefs, IMSL_KNOTS, &users_knots, 0);
    if (imsl_n1rty (1) == 4) {
	imsl_e1mes (0, 0, " ");
	imsl_e1stl (1, "ndata");
	imsl_e1sti (1, ndata);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto RETURN;
    }
    /* GET THE WORKSPACE NEEDED IN B2NAK    */
    xsrt = (Mfloat *) imsl_malloc (ndata * sizeof (*xsrt));
    iwk_b2nak = (Mint *) imsl_malloc (ndata * sizeof (*iwk_b2nak));

    if ((iwk_b2nak == NULL) || (xsrt == NULL)) {
	free_the_structure = 1;
	imsl_e1stl (1, "ndata");
	imsl_e1sti (1, ndata);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto B2NAK_FREE_SPACE;
    }

    if (imsl_n1rty (1))
	goto RETURN;
    /*
     * COMPUTE THE KNOTS IF THE USER DID NOT SUPPLY THEM
     */
    if (knots_given == 0)
	imsl_b2nak (&ndata, xdata, pp->order, pp->knots[0], xsrt, iwk_b2nak);

    /* GET THE WORKSPACE NEEDED IN B2INT    */
    wk1 = (Mfloat *) imsl_malloc (((5 ** order) - 2) * ndata * sizeof (*wk1));
    wk2 = (Mfloat *) imsl_malloc (ndata * sizeof (*wk2));
    wk3 = (Mfloat *) imsl_malloc (ndata * sizeof (*wk3));
    iwk_b2int = (Mint *) imsl_malloc (ndata * sizeof (*iwk_b2int));

    if ((iwk_b2int == NULL) || (wk1 == NULL) || (wk2 == NULL) || (wk3 == NULL)) {
	free_the_structure = 1;
	imsl_e1stl (1, "ndata");
	imsl_e1sti (1, ndata);
	imsl_e1stl (1, "order");
	imsl_e1sti (1, *order);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	goto B2INT_FREE_SPACE;
    }
    if (imsl_n1rty (1))
	goto RETURN;
    /* CALL THE SPLINE ROUTINE */
    imsl_b2int (&ndata, xdata, fdata, pp->order, pp->knots[0], pp->coef[0], wk1, wk2, wk3, iwk_b2int);
   if (imsl_n1rty(1)>3) {
             free_the_structure = 1;
             goto B2INT_FREE_SPACE;
   }

    /* FREE THE WORKSPACE USED */
B2INT_FREE_SPACE:
    if (iwk_b2int != NULL)
	imsl_free (iwk_b2int);
    if (wk3 != NULL)
	imsl_free (wk3);
    if (wk2 != NULL)
	imsl_free (wk2);
    if (wk1 != NULL)
	imsl_free (wk1);
B2NAK_FREE_SPACE:
    if ((free_the_structure == 1) && (pp != NULL))
	imsl_free (pp);
    if (iwk_b2nak != NULL)
	imsl_free (iwk_b2nak);
    if (xsrt != NULL)
	imsl_free (xsrt);
    goto RETURN;
RETURN:
    return argptr;
}
