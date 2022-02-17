#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static Mf_ppoly *pp = NULL;

static VA_LIST_HACK PROTO (l_ppoly_create, (Mint domain_dim, Mint target_dim, Mint *orders,
	            Mint *num_breakpoints, va_list argptr));



#ifdef ANSI
    Mf_ppoly   *imsl_f_ppoly_create (Mint domain_dim, Mint target_dim, Mint *orders,
                Mint *num_breakpoints,...)
#else
    Mf_ppoly   *imsl_f_ppoly_create (domain_dim, target_dim, orders, num_breakpoints, va_alist)
    Mint        domain_dim;
    Mint        target_dim;
    Mint       *orders;
    Mint       *num_breakpoints;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, num_breakpoints);
#ifdef DOUBLE
    imsl_e1psh ("imsl_d_ppoly_create");
#else
    imsl_e1psh ("imsl_f_ppoly_create");
#endif
    IMSL_CALL (l_ppoly_create (domain_dim, target_dim, orders, num_breakpoints, argptr));
    va_end (argptr);
#ifdef DOUBLE
    imsl_e1pop ("imsl_d_ppoly_create");
#else
    imsl_e1pop ("imsl_f_ppoly_create");
#endif
    return pp;
}
#ifdef ANSI
static VA_LIST_HACK l_ppoly_create (Mint domain_dim, Mint target_dim, Mint *orders,
                Mint *num_breakpoints, va_list argptr)
#else
static VA_LIST_HACK l_ppoly_create (domain_dim, target_dim, orders, num_breakpoints, argptr)
    Mint        domain_dim;
    Mint        target_dim;
    Mint       *orders;
    Mint       *num_breakpoints;
    va_list     argptr;
#endif
{
    Mint        arg_number = 4;
    Mvoid      *tmp_ptr = NULL;
    Mvoid      *tmp_ptr_orig = NULL;
    Mint        i;
    Mint        j;
    Mint        code;
    Mint        num_ints;
    Mint        num_floats;
    Mint        num_pointers;
    Mint        temp_product = 1;
    Mint        supplied_breakpoints = 0;
    Mint        not_enough_memory = 0;
    Mint        supplied_coefs = 0;
    Mfloat    **callers_breakpoints = NULL;
    Mfloat    **callers_coefs = NULL;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, int);
	arg_number++;
	switch (code) {
	case IMSL_BREAKPOINTS:
	    supplied_breakpoints = 1;
	    callers_breakpoints = va_arg (argptr, Mfloat **);
	    arg_number++;
	    break;
	case IMSL_COEFS:
	    supplied_coefs = 1;
	    callers_coefs = va_arg (argptr, Mfloat **);
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

    if (imsl_n1rty (0))
	goto RETURN;
/* COMPUTE THE TOTOAL AMOUTNT OF MEMORY NEEDED FOR THE SPLINE STRUCTURE */
    num_ints = 3 * (domain_dim) + 2;
#if (defined(COMPUTER_DECOSF) && defined(CMATH_MINT_TO_LONG))
    /* Since decosf wants longs aligned on 8-byte boundries. */
    num_ints +=1;
#endif
/* THE FOLLOWING IF BLOCK WAS INSERTED BECAUSE THE HP98C
   REQUIRES DOUBLE TO BE ALIGNED ON 8-BYTE BOUNDARIES. */
#if (defined(COMPUTER_HP98C) || defined(COMPUTER_HP97C) || defined(COMPUTER_MIPNT)) || defined(COMPUTER_SUN5) || defined(COMPUTER_GSUN4) && defined(DOUBLE)
    if (mod( (domain_dim + target_dim), 2) == 0)
        num_pointers = 5 + domain_dim + target_dim + 1;
    else
        num_pointers = 5 + domain_dim + target_dim;
#else
    num_pointers = 5 + domain_dim + target_dim;
#endif
    num_floats = 0;
    /* FLOATS FOR THE BREAKPOINTS */
    for (i = 0; i < domain_dim; i++)
	num_floats += (num_breakpoints[i]);
    /* FLOATS FOR THE COEFFICIENTS */
    for (i = 0; i < domain_dim; i++)
	temp_product *= orders[i] * ((num_breakpoints[i]) - 1);
    for (i = 0; i < target_dim; i++)
	num_floats += (target_dim * (temp_product));
    /* GET THE TOTAL SPACE */
    tmp_ptr_orig = (Mvoid *) imsl_malloc (num_ints * sizeof (Mint) +
	num_floats * sizeof (Mfloat) +
	num_pointers * sizeof (tmp_ptr));
    tmp_ptr = tmp_ptr_orig;
    if (tmp_ptr == NULL) {
	not_enough_memory = 1;
	imsl_e1stl (1, "domain_dim");
	imsl_e1sti (1, domain_dim);
	imsl_e1stl (2, "target_dim");
	imsl_e1sti (2, target_dim);
	imsl_ermes (IMSL_FATAL, IMSL_OUT_OF_MEMORY_2);
	goto FREE_SPACE;
    }

    /*
     * Set PP to point to the beginning of the malloc'ed space.
     */
    pp = (Mf_ppoly *) tmp_ptr;
    /* ADJUST THE TEMPORARY POINTER */
    tmp_ptr = (Mvoid *) (pp + 1);
    /* SET POINTER TO POINTER FOR BREAKPOINTS */
    pp->breakpoints = (Mfloat **) tmp_ptr;
    /* ADJUST THE TEMPORARY POINTER */
    tmp_ptr = (Mvoid *) (pp->breakpoints + domain_dim);
    /* SET POINTER TO POINTER FOR COEF */
    pp->coef = (Mfloat **) tmp_ptr;
    /* ADJUST THE TEMPORARY POINTER */
/* THE FOLLOWING IF BLOCK WAS INSERTED BECAUSE THE HP98C
   REQUIRES DOUBLE TO BE ALIGNED ON 8-BYTE BOUNDARIES. */
#if (defined(COMPUTER_HP98C) || defined(COMPUTER_HP97C) || defined(COMPUTER_MIPNT)) || defined(COMPUTER_SUN5) || defined(COMPUTER_GSUN4) && defined(DOUBLE)
    if (mod( (domain_dim + target_dim), 2) == 0)
	tmp_ptr = (Mvoid *) (pp->coef + target_dim + 1);
    else
	tmp_ptr = (Mvoid *) (pp->coef + target_dim);
#else
    tmp_ptr = (Mvoid *) (pp->coef + target_dim);
#endif
    /* SET THE POINTER TO THE SETS OF BREAKPOINTS */
    for (i = 0; i < domain_dim; i++) {
	pp->breakpoints[i] = (Mfloat *) tmp_ptr;
	/* ADJUST THE TEMPORARY POINTER */
	tmp_ptr = (Mvoid *) ((pp->breakpoints[i]) + (num_breakpoints[i]));
    }
    /* SET THE POINTER TO THE SETS OF COEFS */
    for (i = 0; i < target_dim; i++) {
	pp->coef[i] = (Mfloat *) tmp_ptr;
	/* ADJUST THE TEMPORARY POINTER */
	tmp_ptr = (Mvoid *) (pp->coef[i] + temp_product);
    }
    /* SET THE POINTER TO THE SET OF ORDERS */
#if (defined(COMPUTER_DECOSF) && defined(CMATH_MINT_TO_LONG))
    /* Since decosf wants longs aligned on 8-byte boundries. */
    if (mod( ((long)tmp_ptr), 8) == 0)
      pp->order = (Mint *) ((int *)tmp_ptr);
    else
      pp->order = (Mint *) ((int *)tmp_ptr + 1);
#else
    pp->order = (Mint *) ((int *)tmp_ptr);
#endif

    /* Adjust the temporary pointers */
    tmp_ptr = (Mvoid *) (pp->order + domain_dim);
    /* SET THE POINTER TO THE SET OF NUM_COEFS */
    pp->num_coef = (Mint *) tmp_ptr;
    tmp_ptr = (Mvoid *) (pp->num_coef + domain_dim);
    /* SET THE POINTER TO THE SET OF NUM_BREAKPOINTS */
    pp->num_breakpoints = (Mint *) tmp_ptr;
    /* COPY INPUT INFORMATION INTO THE STRUCTURE */
    /* SET DOMAIN_DIM AND RANGE_DIM */
    pp->domain_dim = domain_dim;
    pp->target_dim = target_dim;
    temp_product = 1;
    for (i = 0; i < domain_dim; i++)
	temp_product *= orders[i] * ((num_breakpoints[i]) - 1);
    for (i = 0; i < domain_dim; i++) {
	pp->order[i] = orders[i];
	pp->num_breakpoints[i] = num_breakpoints[i];
	pp->num_coef[i] = temp_product;
    }
    /* COPY THE INPUT BREAKPOINTS IF THEY WERE SUPPLIED */
    if (supplied_breakpoints == 1) {
	for (i = 0; i < domain_dim; i++)
	    for (j = 0; j < pp->num_breakpoints[i]; j++)
		pp->breakpoints[i][j] = callers_breakpoints[i][j];
    }
    /* COPY THE INPUT COEFFICIENTS IF THEY WERE SUPPLIED */
    if (supplied_coefs == 1) {
	for (i = 0; i < target_dim; i++)
	    for (j = 0; j < pp->num_coef[i]; j++)
		pp->coef[i][j] = callers_coefs[i][j];
    }
FREE_SPACE:
    if (not_enough_memory) {
	if (tmp_ptr_orig != NULL)
	    imsl_free (tmp_ptr_orig);
    }
RETURN:
    return (argptr);
}
