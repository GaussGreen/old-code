#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static Mfloat *lv_value = NULL;

 /* work_basis is a global static variable that is used to point to some
    global (within this file) workspace  */

static Mfloat *work_basis = NULL;
 
 /* Prototypes for local codes */

static VA_LIST_HACK PROTO (l_radial_evaluate, (Mint n, Mfloat *x, 
	Mf_radial_basis_fit * radial_struct, va_list argptr));

 /* TOP LEVEL OF SUPER-CODE */

#ifdef ANSI
    Mfloat     *imsl_f_radial_evaluate (Mint n, Mfloat *x,
                	Mf_radial_basis_fit *radial_struct,...)
#else
    Mfloat     *imsl_f_radial_evaluate (n, x, radial_struct, va_alist)
    Mint        n;
    Mfloat     *x;
    Mf_radial_basis_fit *radial_struct;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, radial_struct);
    E1PSH ("imsl_f_radial_evaluate", "imsl_d_radial_evaluate");
    IMSL_CALL (l_radial_evaluate (n, x, radial_struct, argptr));
    va_end (argptr);
    E1POP ("imsl_f_radial_evaluate", "imsl_d_radial_evaluate");
    return lv_value;
}

 /* 2ND OF  LEVEL SUPER-CODE */

#ifdef ANSI
static VA_LIST_HACK l_radial_evaluate (Mint n, Mfloat *x, 
		Mf_radial_basis_fit *radial_struct, va_list argptr)
#else
static VA_LIST_HACK l_radial_evaluate (n, x, radial_struct, argptr)
    Mint        n;
    Mfloat     *x;
    Mf_radial_basis_fit *radial_struct;
    va_list     argptr;
#endif
{
    Mint        code;
    Mint        dimension;
    Mint        num_centers=0;
    Mint        arg_number = 3;
    Mint        return_user = 0;
    Mint        i;
    Mint        j;
    Mint        k;
    Mint        which_point;
    static Mfloat      l_work[50];
    Mfloat      distance;
    Mfloat      xnew = 0.0;

    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_RETURN_USER:
	    lv_value = va_arg (argptr, Mfloat *);
	    if (!lv_value) {
		imsl_e1stl (1, "value");
		imsl_e1stl (2, "IMSL_RETURN_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		goto RETURN;
	    }
	    ++arg_number;
	    return_user = 1;
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


    /* ERROR CHECKING OF INPUT ARGUMENTS */

    /* Check n < 1 */

    if (n < 1) {
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NUM_POINTS_LESS_THAN_ONE);
	goto RETURN;
      /*printf ("PROBLEM with n < 1 in l_radial_evaluate.\n");*/
    }

    /* Check pointer to the structure. */

    if (radial_struct == NULL) {
	imsl_ermes (IMSL_TERMINAL, IMSL_INVALID_RADIAL_STRUCT);
	goto RETURN;
      /*printf ("PROBLEM with radial_struct == NULL in l_radial_evaluate.\n");*/
    }

    /* Check dimension */

    if (radial_struct->dimension < 1) {
	imsl_e1sti (1, radial_struct->dimension);
	imsl_ermes (IMSL_TERMINAL, IMSL_DIM_LESS_THAN_ONE);
	goto RETURN;
      /*printf ("PROBLEM with dimension in l_radial_evaluate \n");*/
    }

    /* Check num_centers */

    if (radial_struct->num_centers < 1) {
	imsl_e1sti (1, radial_struct->num_centers);
	imsl_ermes (IMSL_TERMINAL, IMSL_NUM_CENTERS_LESS_THAN_ONE);
	goto RETURN;
      /*printf ("PROBLEM with num_centers in l_radial_evaluate \n");*/
    }

     /* Check additional_terms. The following three cases are valid: no
        additional terms:radial_struct->additional_terms = 0 A constant term
        in the model: radial_struct->additional_terms = 1 A linear   term in
        the model: radial_struct->additional_terms =
        (radial_struct->dimension)+1) */

    if ((radial_struct->additional_terms < 0) &&
	(radial_struct->additional_terms != 1) &&
	(radial_struct->additional_terms != ((radial_struct->dimension) + 1))) {
      /*printf ("PROBLEM with additionl_terms in l_radial_evaluate \n");*/
	goto RETURN;
    }

    /* Check pointer to centers */

    if (radial_struct->centers == NULL) {
	imsl_ermes (IMSL_TERMINAL, IMSL_INVALID_CENTERS_POINTER);
	goto RETURN;
      /*printf ("PROBLEM with centers in l_radial_evaluate \n");*/
    }

    /* Check pointer to coefs */

    if (radial_struct->coefficients == NULL) {
	imsl_ermes (IMSL_TERMINAL, IMSL_INVALID_COEFF_POINTER);
	goto RETURN;
      /*printf ("PROBLEM with centers in l_radial_evaluate \n");*/
    }

    /* Check the pointer to radial function. */

    if (radial_struct->radial_function == NULL) {
	imsl_ermes (IMSL_TERMINAL, IMSL_INVALID_RAD_FUNC_POINTER);
	goto RETURN;
      /*printf ("PROBLEM with radial function in l_radial_evaluate \n");*/
    }

    /* GET SPACE FOR RETURNED VALUES IF NEEDED  */

     /* if the user did not supply space for the returned values, malloc the
        space */

    if (!return_user) {
	lv_value = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));
	if (lv_value == NULL) {
	    imsl_e1stl (1, "n");
	    imsl_e1sti (1, n);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}
    }
    
     /* work_basis is a global(within this file) static variable that is used
        to point to some global workspace.  To save some time, there is an
        array l_work, that can be used for this space if dimension is less
        than 50.  This elliminates the need to malloc every time this routine
        is called */
   
    if (radial_struct->dimension <= 50) {
	work_basis = l_work;
    }
    else {
	work_basis = (Mfloat *) imsl_malloc ((radial_struct->dimension) *
	    sizeof (*work_basis));
	if ((work_basis == NULL)) {
	    imsl_e1stl (1, "radial_struct->dimension");
	    imsl_e1sti (1, radial_struct->dimension);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}
    }
    
     /* Start computing the value of the radial basis fit at the input point */
    
    for (which_point = 0; which_point < n; ++which_point) {
	xnew = 0.0;
	for (j = 0; j < radial_struct->num_centers; ++j) {
	    for (k = 0; k < radial_struct->dimension; ++k)
		work_basis[k] =
		    x[(radial_struct->dimension) * which_point + k] -
		    (radial_struct->centers)[j * radial_struct->dimension + k];
	    distance = imsl_snrm2 (radial_struct->dimension, work_basis, 1);
	    xnew += (radial_struct->radial_function)
		(distance) * ((radial_struct->coefficients)[j]);

	}

	
	 /* If a constant or linear term was used, then compute their part in
	    the answer */
	 

	if (radial_struct->additional_terms > 0)
	    xnew += radial_struct->coefficients[num_centers];
	for (i = 1; i < radial_struct->additional_terms; ++i) {
	    xnew += radial_struct->coefficients[num_centers + i] * x[i - 1];
	}

	/* On errors return NaN. */
	if (imsl_n1rty (0) > 3)
	    lv_value[which_point] = imsl_amach (6);
	else
	    lv_value[which_point] = xnew;
    }
FREE_SPACE:
    if ((work_basis != NULL) && (radial_struct->dimension > 50)) {
	imsl_free (work_basis);
	work_basis = NULL;
    }
    if ((imsl_n1rty (0) > 3) && (lv_value != NULL)) {
	imsl_free (lv_value);
	lv_value = NULL;
    }

RETURN:
    return argptr;
}
