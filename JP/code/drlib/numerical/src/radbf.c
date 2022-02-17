#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

 /* Pointer to space returned */
static Mfloat *coefficients = NULL;
static Mf_radial_basis_fit *returned_structure = NULL;
 /*
  * Pointer to some workspace that is in various routines in this file.
  * Rather than malloc it many times, it is going to be available to all
  * routines in this file.
  */
static Mfloat *work_basis = NULL;

static Mint which_solver = 0; /* svd */
static Mfloat tol;

 /* --------------------------------------------------- */
 /*
  * In order to pass a function as an optional argument, the following
  * typedef is introduced.
  */
 /* --------------------------------------------------- */

 /* --------------------------------------------------- */
 /* Prototypes for local codes                          */
 /* --------------------------------------------------- */

static VA_LIST_HACK PROTO (l_radial_scattered_fit, (Mint dimension,
	            Mint num_points,
	            Mfloat *abscissae,
	            Mfloat *fdata,
	            Mint num_centers,
	            va_list argptr));
    static Mvoid PROTO (l_r3dbf,
                (Mint *, Mint *, Mint *, Mfloat *, Mfloat *, Mint *, Mfloat *));
    Mfloat      PROTO (imsl_f_radial_function,
                (Mfloat));

    static void PROTO (l_r6dbf,
                (Mint *, Mint *, Mint *, Mfloat *, Mfloat *, Mfloat *));
    static Mfloat *PROTO (l_r4dbf,
                (Mfloat (*basis_function) (Mfloat),
	            Mint *, Mint *, Mint *, Mfloat *, Mfloat *, Mfloat *,
	            Mfloat[], Mint *));


 /* --------------------------------------------------- */
 /* TOP LEVEL OF SUPER-CODE              */
 /* --------------------------------------------------- */
#ifdef ANSI
    Mf_radial_basis_fit *imsl_f_radial_scattered_fit (Mint dimension,
                Mint num_points,
                Mfloat *abscissae,
                Mfloat *fdata,
                Mint num_centers,...)
#else
    Mf_radial_basis_fit *imsl_f_radial_scattered_fit (dimension, num_points,
                abscissae, fdata,
                num_centers,
                va_alist)
    Mint        dimension;
    Mint        num_points;
    Mfloat     *abscissae;
    Mfloat     *fdata;
    Mint        num_centers;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, num_centers);
    E1PSH ("imsl_f_radial_scattered_fit", "imsl_d_radial_scattered_fit");
    IMSL_CALL (l_radial_scattered_fit (dimension, num_points,
	    abscissae, fdata,
	    num_centers,
	    argptr));
    va_end (argptr);
    E1POP ("imsl_f_radial_scattered_fit", "imsl_d_radial_scattered_fit");
    return returned_structure;
}



 /* --------------------------------------------------- */
 /* 2ND OF  LEVEL SUPER-CODE             */
 /* --------------------------------------------------- */


#ifdef ANSI
static VA_LIST_HACK l_radial_scattered_fit (Mint dimension, Mint num_points,
                Mfloat *abscissae, Mfloat *fdata,
                Mint num_centers,
                va_list argptr)
#else
static VA_LIST_HACK l_radial_scattered_fit (dimension, num_points, abscissae,
                fdata, num_centers, argptr)
    Mint        dimension;
    Mint        num_points;
    Mfloat     *abscissae;
    Mfloat     *fdata;
    Mint        num_centers;
    va_list     argptr;
#endif
{
    Mint        arg_number = 5;
    Mint        i;
    Mint        weights_given = 0;
    Mint        basis_given = 0;
    Mint        delta_given = 0;
    Mint        centers_given = 0;
    Mint        random_seed_set = 0;
    Mint        coefficients_return_user = 0;
    Mint        linear_term = 0;
    Mint        constant_term = 0;
    Mint        additional_columns = 0;
    Mint        code = 0;
    Mint        num_weights_zero = 0;
    Mint        random_seed = 234579;
    Mfloat     *weights = NULL;
    Mfloat     *users_centers = NULL;
    Mfloat     *least_squares_matrix = NULL;
    Mfloat     *work_rhs = NULL;
    Mfloat      delta = 1.0;
    Mfloat      ratio = 0.5;
    Mf_radial_fcn basis = NULL;
    Mf_radial_fcn users_basis = NULL;


    Mint        num_ints = 0;
    Mint        num_floats = 0;
    Mint        num_pointers = 0;
    Mvoid      *tmp_ptr = NULL;
    Mvoid      *tmp_ptr_orig = NULL;

    tol = sqrt (imsl_amach (4));

    code = 1;
    while (code > 0) {
	code = va_arg (argptr, int);
	arg_number++;
	switch (code) {
	case IMSL_RETURN_USER:
	    coefficients_return_user = 1;
	    coefficients = va_arg (argptr, Mfloat *);
	    arg_number++;
	    break;
	case IMSL_CENTERS:
	    centers_given = 1;
	    users_centers = va_arg (argptr, Mfloat *);
	    arg_number++;
	    break;
	case IMSL_SUPPLY_BASIS:
	    basis_given = 1;
	    users_basis = va_arg (argptr, Mf_radial_fcn);
	    arg_number++;
	    break;
	case IMSL_SUPPLY_DELTA:
	    delta_given = 1;
	    delta = va_arg (argptr, Mdouble);
	    arg_number++;
	    break;
	case IMSL_SUPPLY_DELTA_ADR:
	    delta_given = 1;
	    delta = *va_arg (argptr, Mfloat *);
	    arg_number++;
	    break;
	case IMSL_WEIGHTS:
	    weights_given = 1;
	    weights = va_arg (argptr, Mfloat *);
	    arg_number++;
	    break;
	case IMSL_CENTERS_RATIO_ADR:
	    ratio = *va_arg (argptr, Mfloat *);
	    arg_number++;
	    break;
	case IMSL_CENTERS_RATIO:
	    ratio = va_arg (argptr, Mdouble);
	    arg_number++;
	    break;
	case IMSL_LINEAR_TERM:
	    linear_term = 1;
	    arg_number++;
	    break;
	case IMSL_CONSTANT_TERM:
	    constant_term = 1;
	    arg_number++;
	    break;
	case IMSL_RANDOM_SEED:
	    random_seed_set = 1;
	    random_seed = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_NO_SVD:
	    which_solver = 1;
	    break;
	case IMSL_TOL:
/*
 * GCC 4.01 doesn't like this and warns that the program will abort
 * if control ever gets here.  I think this is because of implicit
 * conversion of float to double in ANSI C function calls.
 */

/*	    tol = va_arg (argptr, Mfloat); */
	    tol = va_arg (argptr, double); 
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

    /* ----------ERROR CHECKING OF INPUT ARGUMENTS ----------------------- */
    /* THESE CHECKS WILL BE REPLACED WITH ACTUAL ERROR-HANDLER STUFF LATER. */
    /* Check dimension */
    if (dimension < 1) {
	imsl_e1sti (1, dimension);
	imsl_ermes (IMSL_TERMINAL, IMSL_DIM_LESS_THAN_ONE);
	goto RETURN;
    }
    /* Check num_points */
    if (num_points < 1) {
	imsl_e1sti (1, num_points);
	imsl_ermes (IMSL_TERMINAL, IMSL_NUM_POINTS_LESS_THAN_ONE);
	goto RETURN;
    }
    /* Check num_centers */
    if (num_centers < 1) {
	imsl_e1sti (1, num_centers);
	imsl_ermes (IMSL_TERMINAL, IMSL_NUM_CENTERS_LESS_THAN_ONE);
	goto RETURN;
    }
    if (num_centers > num_points) {
	imsl_e1sti (1, num_centers);
	imsl_e1sti (2, num_points);
	imsl_ermes (IMSL_TERMINAL, IMSL_CENTERS_GT_POINTS);
	goto RETURN;
    }
    /* Check pointer to abscissae */
    if (abscissae == NULL) {
	imsl_e1stl (1, "abscissae");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
	goto RETURN;
    }
    /* Check pointer to fdata */
    if (fdata == NULL) {
	imsl_e1stl (1, "fdata");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
	goto RETURN;
    }
    /* Check the random seed. */
    if (random_seed < 0) {
	imsl_e1sti (1, random_seed);
	imsl_ermes (IMSL_TERMINAL, IMSL_SEED_NEGATIVE);
	goto RETURN;
    }
    /* Check that user did not specify both LINEAR_TERM and CONSTANT_TERM */
    if ((linear_term) && (constant_term)) {
	goto RETURN;
    }
    /* Check the pointer to users centers if they were specified. */
    if ((centers_given) && (users_centers == NULL)) {
	imsl_e1stl (1, "centers");
	imsl_e1stl (2, "IMSL_CENTERS");
	imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	goto RETURN;
    }
    /* Check the pointer to users weights if they were specified. */
    if ((weights_given) && (weights == NULL)) {
	imsl_e1stl (1, "weights");
	imsl_e1stl (2, "IMSL_WEIGHTS");
	imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	goto RETURN;
    }
    /* Check the pointer to users return space if it was specified. */
    if ((coefficients_return_user) && (coefficients == NULL)) {
	/*printf ("PROBLEM with users return space in l_radial_scattered_fit \n");*/
	goto RETURN;
    }
    /* Check the pointer to users radial function if it was specified. */
    if ((basis_given) && (users_basis == NULL)) {
	imsl_e1stl (1, "radial_function");
	imsl_e1stl (2, "IMSL_SUPPLY_BASIS");
	imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	goto RETURN;
    }
    /* ------------ CHECKS COMPLETED -------------------------------- */
    /*
     * Compute the number of addional columns that will be needed in the
     * coefficient matrix h because of either CONSTANT_TERM  or LINEAR_TERM
     * was specified.
     */
    if (constant_term)
	additional_columns = 1;
    if (linear_term)
	additional_columns = dimension + 1;

    /*
     * If the user supplied a basis function, then use it, otherwise use our
     * own, which is named imsl_f_radial_function, and is located within this
     * file..
     */
    if (!basis_given) {
	basis = imsl_f_radial_function;
    }
    else {
	basis = users_basis;
    }


    /*
     * Generate the space for the returned structure and partition it into
     * its appropriate parts.  By getting all the space with one malloc, it
     * will only take one free to release the space.
     */

    num_ints = 3;
    num_floats = 1 + (dimension + 1) * num_centers + additional_columns;
    num_pointers = 3;

    tmp_ptr_orig = (Mvoid *) imsl_malloc (num_ints * sizeof (Mint) +
	num_floats * sizeof (Mfloat) +
	num_pointers * sizeof (tmp_ptr));
    tmp_ptr = tmp_ptr_orig;
    if (tmp_ptr == NULL) {
	imsl_e1stl (1, "dimension");
	imsl_e1sti (1, dimension);
	imsl_e1stl (2, "num_center");
	imsl_e1sti (2, num_centers);
	imsl_ermes (IMSL_FATAL, IMSL_OUT_OF_MEMORY_2);
	goto FREE_SPACE;
    }


    /*
     * Assign the address of the starting point of the space to the return
     * value of this function.
     */
    returned_structure = (Mf_radial_basis_fit *) tmp_ptr_orig;
    /*
     * Move the temporary pointer to the address immediately following the
     * space for the definition of the strucure.
     */
    tmp_ptr = (Mvoid *) (returned_structure + 1);

    /* Set the pointer to the centers */
    returned_structure->centers = (Mfloat *) tmp_ptr;

    /*
     * Move the temporary pointer to the address immediately following the
     * space for the centers.
     */
    tmp_ptr = (Mvoid *) (returned_structure->centers + (dimension * num_centers));

    /* Set the pointer to the coefficients */
    returned_structure->coefficients = (Mfloat *) tmp_ptr;

    /*
     * Fill out the rest of the stucture with the information we currently
     * have available.
     */

    returned_structure->dimension = dimension;
    returned_structure->num_centers = num_centers;
    returned_structure->additional_terms = additional_columns;
    returned_structure->delta = delta;
    returned_structure->radial_function = basis;





    /* ------------ CHECK OR COMPUTE WEIGHTS ------------------------ */

    /*
     * SET UP THE WEIGHTS IF THEY WERE NOT SUPPLIED BY THE USER.  IF THEY
     * WERE SUPPLIED, CHECK THAT THEY ARE ALL NONNEGATIVE AND AT LEAST ON IS
     * POSITIVE.
     */
    if (weights_given) {
	for (i = 0; i < num_points; i++) {
	    if (weights[i] == F_ZERO)
		num_weights_zero++;
	    if (weights[i] < F_ZERO) {
		imsl_e1sti (1, i);
		imsl_e1str (1, weights[i]);
		imsl_e1stl (1, "X");
		imsl_ermes (IMSL_FATAL, IMSL_NEGATIVE_WEIGHTS);
		goto RETURN;
	    }
	}
	if (num_weights_zero == num_points) {
	    imsl_ermes (IMSL_TERMINAL,
		IMSL_SPLINE_NO_POS_ELMNT);
	    goto RETURN;
	}
    }
    /*
     * If this else is reached, we need to get space for the weights and fill
     * it with ones.
     */
    else {
	weights = (Mfloat *) imsl_malloc ((num_points) * sizeof (*weights));
	if (weights != NULL)
	    sset (num_points, F_ONE, weights, 1);
	if (weights == NULL) {
	    imsl_e1stl (1, "num_points");
	    imsl_e1sti (1, num_points);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}
    }
    /* ------------ FINISHED WITH WEIGHTS --------------------------- */

    /* ------------ COMPUTE THE CENTERS ----------------------------- */

    /*
     * If the user did not supply centers to use, then compute our own based
     * on the original data, and the ratio of gridded to randomly placed
     * centers.  This ratio is stored in the variable ratio.
     */


    if (!centers_given) {
	l_r3dbf (&num_centers, &dimension, &num_points, abscissae,
	    &ratio, &random_seed,
	    (Mfloat *) returned_structure->centers);
	/* If an error occured in l_r3bdf, free space and exit.  */
	if (imsl_n1rty (0) > 3)
	    goto FREE_SPACE;
    }
    else {
	for (i = 0; i < dimension * num_centers; ++i) {
	    returned_structure->centers[i] = users_centers[i];
	}

    }
    /* ------------ CENTERS COMPUTED -------------------------------- */

    /* ------------ COMPUTE THE LEAST SQUARES  MATRIX --------------- */
    /* Get some workspace. */
    work_basis = (Mfloat *) imsl_malloc ((dimension) * sizeof
	(*work_basis));
    if ((work_basis == NULL)) {
	imsl_e1stl (1, "dimension");
	imsl_e1sti (1, dimension);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }
    /*
     * If the user supplied a basis function, then use id, otherwise use our
     * own.
     */
    if (!basis_given) {
	basis = imsl_f_radial_function;
    }
    else {
	basis = users_basis;
    }

    /* l_r4dbf returns a pointer to the coefficien matrix h.  */
    least_squares_matrix = l_r4dbf (basis, &dimension, &num_points,
	&num_centers, &delta, abscissae,
	returned_structure->centers,
	weights, &additional_columns);
    /* If an error occured in l_r4bdf, free space and exit.  */
    if (imsl_n1rty (0) > 3)
	goto FREE_SPACE;
    /* ------------ LEAST_SQUARES MATRIX COMPUTED ------------------- */
    /* ------------ COMPUTE THE RHS VECTOR -------------------------- */
    work_rhs = (Mfloat *) imsl_malloc ((num_points) * sizeof (*work_rhs));
    if ((work_rhs == NULL)) {
	imsl_e1stl (1, "num_points");
	imsl_e1sti (1, num_points);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }
    /* Scale the right hand side to allow for the weights */
    for (i = 0; i < num_points; ++i) {
	work_rhs[i] = sqrt (weights[i]) * fdata[i];
    }
    /* ------------ COMPUTE THE COEFFICIENTS ------------------------ */
    l_r6dbf (&num_points, &num_centers, &additional_columns, least_squares_matrix, work_rhs,
	returned_structure->coefficients);

    goto FREE_SPACE;
    /* ------------ COEFFICIENTS COMPUTED --------------------------- */

FREE_SPACE:
    if ((imsl_n1rty (0) > 3) && (returned_structure != NULL)) {
	imsl_free (returned_structure);
	returned_structure = NULL;
    }

    if ((!weights_given) && (weights != NULL)) {
	imsl_free (weights);
	weights = NULL;
    }
    if (least_squares_matrix != NULL) {
	imsl_free (least_squares_matrix);
	least_squares_matrix = NULL;
    }
    if (work_rhs != NULL) {
	imsl_free (work_rhs);
	work_rhs = NULL;
    }
    if (work_basis != NULL) {
	imsl_free (work_basis);
	work_basis = NULL;
    }


RETURN:
    return argptr;
}



/*  -----------------------------------------------------------------------
    IMSL Name:  l_r3dbf

    Revised:    July, 1, 1991

    Purpose:    Comute the values of the default centers used by the
                radial basis function.

    Usage:      l_radial_centers (Mint *num_centers,
                                  Mint *dimension,
                                  Mint *num_points,
                                  Mfloat *abscissae
                                  Mfloat *ratio)

    Arguments:
  *num_centers- Pointer to the number of centers to be computed.  (Input)
   *dimension - Number of domain dimensions.  (Input)
  *num_points - Pointer to the number of points the centers are to be computed
                with respect to.  (Input)
 *abscissae   - Pointer to abscissa values the centers are to be
                computed with respect to.  (Input)
      *ratio  - Pointer to a  floating point value between 0.0 and 1.0
                defined by: ratio = (number_gridded)/num_centers).  (Input)
 *random_seed - Pointer to the random seed to use.  (Input)
   *centers   - The space into which the centers will be stored.  (Output)

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mvoid l_r3dbf (Mint *num_centers, Mint *dimension, Mint *num_points,
                Mfloat *abscissae, Mfloat *ratio, Mint *random_seed,
                Mfloat *centers)
#else
static Mvoid l_r3dbf (num_centers, dimension, num_points, abscissae,
                ratio, random_seed, centers)
    Mint       *num_centers;
    Mint       *dimension;
    Mint       *num_points;
    Mfloat     *abscissae;
    Mfloat     *ratio;
    Mint       *random_seed;
    Mfloat     *centers;
#endif
{
    Mint        i;
    Mint        j;
    Mint        kk;
    Mfloat     *abscissae_small = NULL;
    Mfloat     *abscissae_big = NULL;
    Mfloat     *interval_size = NULL;
    Mint       *m = NULL;
    Mint       *p = NULL;
    Mint        r;
    Mfloat     *random_points = NULL;
    Mfloat     *c_t = NULL;
    Mfloat      temp_abscissae;
    Mint        old_random_seed;
    Mint        ng0;
    Mint        ng1;
    Mint        nn;
    Mint        num_gridded;
    Mint        num_random;
    imsl_e1psh ("l_r3dbf");
    /*
     * GET THE VALUE OF THE CURRENT RANDOM NUMBER SEED SO WE CAN RESET IT
     * WHEN WE EXIT HIS ROUTINE.
     */
    old_random_seed = imsl_random_seed_get ();
    /* SET THE RANDOM NUMBER SEED WITH THE INPUT SEED. */
    imsl_random_seed_set (*random_seed);
    /*
     * GET SOME SPACE TO HOLD VECTORS OF VALUES RELATING TO THE INTERVALS FOR
     * EACH DIMENSION.
     */
    abscissae_small = (Mfloat *) imsl_malloc (((*dimension)) *
	sizeof (*abscissae_small));
    abscissae_big = (Mfloat *) imsl_malloc (((*dimension)) *
	sizeof (*abscissae_big));
    m = (Mint *) imsl_malloc (((*dimension)) *
	sizeof (*m));
    interval_size = (Mfloat *) imsl_malloc (((*dimension)) *
	sizeof (*interval_size));
    if ((abscissae_small == NULL) || (abscissae_big == NULL)
	|| (m == NULL) || (interval_size == NULL)) {
	imsl_e1stl (1, "num_centers");
	imsl_e1sti (1, (*num_centers));
	imsl_e1stl (2, "dimension");
	imsl_e1sti (2, (*dimension));
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	goto FREE_SPACE;
    }
    /*
     * FIGURE OUT HOW MANY OF EACH TYPE OF CENTER IS NEEDED NUMBER OF
     * GRIDDED.  THE METHOD USED IS AS FOLLOWS: 1.  ng0 = the desired ratio
     * of gridded centers versus random centers. 2.  ng1 = the number of
     * centers values in each dimension.  We impose the restriction that
     * there will be the same number of centers values in each dimension. So,
     * for example, if (*dimension) = 2, (*ratio) = .75 and num centers = 35,
     * then ng0 =26, so there will be AT most 26 gridded centers.  ng1 then
     * is computed as the (dimension)th root of 26, ng1 = 5. 3.  So
     * num_gridded is 25, and num_random is 1.
     */

    ng0 = (*ratio) * (*num_centers);
    ng1 = pow ((Mdouble) ng0, 1.0 / (Mdouble) (*dimension));
    num_gridded = pow ((Mdouble) ng1, (Mdouble) (*dimension));
    /* NUMBER OF RANDOM */
    num_random = (*num_centers) - num_gridded;
    /* COMPUTE SOME INFORMATION ABOUT THE INTERVALS IN EACH DIMENSION. */


    for (i = 0; i < (*dimension); ++i) {
	m[i] = 1;
	abscissae_small[i] = abscissae[i];
	abscissae_big[i] = abscissae[i];

	for (j = i; j < (*num_points); j++) {
	  abscissae_small[i] = (abscissae_small[i] < abscissae[j * (*dimension) + i]) ? abscissae_small[i] : abscissae[j * (*dimension) + i];
	  abscissae_big[i] = (abscissae_big[i] > abscissae[j * (*dimension) + i]) ? abscissae_big[i] : abscissae[j * (*dimension) + i];

	}
	interval_size[i] = fabs (abscissae_big[i] - abscissae_small[i]);
    }

    if (num_gridded == 1) {
	for (i = 0; i < (*dimension); ++i)
	    centers[i] = abscissae_small[i] + .5 * interval_size[i];
	goto COMPUTE_RANDOM;
    }
    /*
     * THE FOLLOWING LOOP(S) ARE A NON_RECURSIVE IMPLEMENTATION OF dimension
     * NESTED LOOPS.  THIS IS WHERE THE num_gridded CENTERS ARE COMPUTED.
     */

    kk = -1;
    c_t = centers;
    nn = ng1;
    m[0] = 0;
L_30:m[0] = m[0] + 1;
    if ((*dimension) == 1)
	goto L_50;
    for (i = 1; i < (*dimension); ++i) {
	if (m[i - 1] <= nn) {
	    kk++;
	    goto L_50;
	}
	kk++;
	m[i - 1] = 1;
	m[i] = m[i] + 1;
    }
L_50:if (m[(*dimension) - 1] > nn)
	goto COMPUTE_RANDOM;
    for (i = 0; i < (*dimension); ++i) {
	j = m[i];
	c_t[i] = abscissae_small[i] + interval_size[i] * ((Mfloat) (mod (j, ng1)) / (Mfloat) (ng1 -
		1));
    }
    c_t = &(c_t[(*dimension)]);
    goto L_30;
    /*
     * AT THIS POINT, ALL THE GRIDDED CENTERS HAVE BEEN COMPUTED, IF ANY. WE
     * NOW START TO COMPUTE THE RANDOMLY DISTRIBUTED CENTERS.  THESE CENTERS
     * ARE RANDOMLY DISTRIBUTED THROUGHOUT THE DOATA POINTS.
     */
COMPUTE_RANDOM:

    if (num_random > 0) {


/* GENERATE A RANDOM SUBSET OF SIZE NUM_RANDOM OF THE
   ORIGINAL DATA POINTS.  SEE BOOK COMBINATORIAL ALGORITHMS BY
   REINGOLD, PAGE 189. */

	p = (Mint *) imsl_malloc (((*num_points)) * sizeof (*p));
	random_points = (Mfloat *) imsl_malloc ((num_random) * (*dimension) * sizeof (*random_points));

	if ((p == NULL) || (random_points == NULL)) {
	    imsl_e1stl (1, "num_points");
	    imsl_e1sti (1, (*num_points));
	    imsl_e1stl (2, "num_random");
	    imsl_e1sti (2, num_random);
	    imsl_e1stl (3, "dimension");
	    imsl_e1sti (3, (*dimension));
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_3);
	    goto FREE_SPACE;
	}

	for (j = 0; j < (*num_points); ++j) {
	    p[j] = j;
	}
	for (j = 0; j < num_random; ++j) {
	    imsl_f_random_uniform (1, IMSL_RETURN_USER, &temp_abscissae, 0);
	    r = j + (Mint) ((*num_points) - j) * temp_abscissae;
	    for (i = 0; i < (*dimension); ++i) {
		random_points[j * (*dimension) + i] = abscissae[(*dimension) * p[r] + i];
	    }
	    p[r] = p[j];

	}


    }

/* COPY THE RANDOM SET OF POINTS TO CENTERS. */
    for (i = 0; i < num_random; i++) {
	for (j = 0; j < (*dimension); j++) {
	    centers[(num_gridded + i) * (*dimension) + j] = random_points[(*dimension) * i + j];
	}
    }


FREE_SPACE:

    if (abscissae_small != NULL) {
	imsl_free (abscissae_small);
	abscissae_small = NULL;
    }
    if (abscissae_big != NULL) {
	imsl_free (abscissae_big);
	abscissae_big = NULL;
    }
    if (m != NULL) {
	imsl_free (m);
	m = NULL;
    }
    if (interval_size != NULL) {
	imsl_free (interval_size);
	interval_size = NULL;
    }
    if (p != NULL) {
	imsl_free (p);
	p = NULL;
    }
    if (random_points != NULL) {
	imsl_free (random_points);
	random_points = NULL;
    }
    /* RESET THE RANDOM NUMBER SEED TO THE VALUE BEFOR ENTERING THIS ROUTINE. */
    imsl_random_seed_set (old_random_seed);
    imsl_e1pop ("l_r3dbf");
    return;
}



/*  -----------------------------------------------------------------------
    IMSL Name:  l_r6dbf

    Revised:    July, 1, 1991

    Purpose:    Compute the radial basis function coefficients by solving
                a least squares problem.

    Usage:      l_r6dbf        (Mint *num_points, Mint *num_centers,
                                Mint *additional_columns,
                                Mfloat *least_squares_matrix,
                                Mfloat *rhs)

    Arguments:
   *num_points - The number of points in the least squares fitting problem.
                (Input)
   *num_centers-pointer to the number of columns in the radial basis fit.
                (Input)
   *additional_columns - pointer the the number of additional columns of the
                least squares matrix as a result of either a constant or
                linear term being added to the model.  (Input)
   *least_squares_matrix  - The num_points by (num_centers+addtional_columns)
                matrix containing the
                coefficient matrix of the least squares system to be
                solved.  (Input)
        *rhs  - array of length num_points containing the right-hand
                side of the least squares system.

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r6dbf (Mint *num_points,
                Mint *num_centers,
                Mint *additional_columns,
                Mfloat *least_squares_matrix,
                Mfloat *rhs, Mfloat *solution)
#else
static void l_r6dbf (num_points, num_centers, additional_columns,
                least_squares_matrix, rhs, solution)
    Mint       *num_points;
    Mint       *num_centers;
    Mint       *additional_columns;
    Mfloat     *least_squares_matrix;
    Mfloat     *rhs;
    Mfloat     *solution;
#endif
{
    Mfloat      tol = sqrt (imsl_amach (4));
    Mfloat     *s;
    Mfloat     *u;
    Mfloat     *v;
    Mfloat     *accum;
    Mint        i;
    Mint        min;
    imsl_e1psh ("l_r6dbf");


    if (which_solver) {
        solution = imsl_f_lin_least_squares_gen (*num_points,
	    *num_centers + (*additional_columns),
	    least_squares_matrix,
	    rhs,
/*          IMSL_RESIDUAL, &residual,*/
	    IMSL_RETURN_USER, solution, 0);
/*      imsl_f_write_matrix("residual", *num_points,1,residual,0);*/
    }
    else {

    s = imsl_f_lin_svd_gen (*num_points, *num_centers + *additional_columns,
	least_squares_matrix,
	IMSL_U, &u,
	IMSL_V, &v,
	0);
    if (imsl_error_code () >= 4L)
	goto RETURN;

    min = imsl_i_min (*num_points, *num_centers + *additional_columns);
    for (i = 0; i < min; i++) {
	if (s[i] > tol)
	    s[i] = 1.0 / s[i];
	else
	    s[i] = 0.0;
    }

    accum = imsl_f_mat_mul_rect ("trans(A)*x", IMSL_A_MATRIX,
	*num_points, min, u,
	IMSL_X_VECTOR, *num_points, rhs,
	0);

    for (i = 0; i < min; i++)
	accum[i] *= s[i];

    imsl_f_mat_mul_rect ("A*x", IMSL_A_MATRIX,
	*num_centers + *additional_columns, min, v,
	IMSL_X_VECTOR, min, accum,
	IMSL_RETURN_USER, solution,
	0);

RETURN:

    if (s != NULL)
	free (s);
    if (u != NULL)
	free (u);
    if (v != NULL)
	free (v);
    if (accum != NULL)
	free (accum);
    }

    imsl_e1pop ("l_r6dbf");
    return;
}
/*  -----------------------------------------------------------------------
    IMSL Name: imsl_f_radial_function

    Revised:    July, 1, 1991

    Purpose:    Compute the value of the default radial basis function.

    Usage:      imsl_f_basis_function (Mint dim, Mfloat *x, Mfloat *center,
                                Mfloat delta)

    Arguments:
      distance-
      delta   - Scalar value used in the default radial basis funtion.
              (Input)


    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */


#ifdef ANSI
Mfloat      imsl_f_radial_function (Mfloat l_distance)
#else
Mfloat      imsl_f_radial_function (l_distance)
    Mfloat      l_distance;
#endif
{
    Mfloat      distance;
    Mfloat      delta;
    distance = (Mfloat) l_distance;
    delta = (Mfloat) (returned_structure->delta);
    return ((Mfloat) sqrt (delta * delta + distance * distance));
}


/*  -----------------------------------------------------------------------
    IMSL Name:  l_radial_build_h

    Revised:    July, 1, 1991

    Purpose:    Compute the entries in the coefficient matrix used in the
                least squares problem to compute the radial basis function
                coefficients.


    Usage:      l_radial_build_h (Mfloat (*basis_function)(Mfloat, Mfloat),
                                  Mint *dim, Mint *num_points,
                                  Mint *num_centers, Mfloat *delta,
                                  Mfloat *x, Mfloat *centers)

    Arguments:
   basis_function - Radial basis function used to compute the entries in
                    the returned matrix.
   *dimension - Number of domain dimensions.  (Input)
  *num_points - The number of points in the least squares fitting problem.
  *num_centers- The number of centers to use. (Input)
     *delta   - Scalar value used in the default radial basis funtion.
                (Input)
   *abscissae  - Pointer to abscissa values the radial basis function is to
                be evaluated at.  (Input)
   *centers   - Pointer to the centers to be used to compute the fit.  (Input)
   *weights   - Weights of the problem. (Input)
  additional_columns - the number of additional columns in the model
                       as a result of either a constant or linear term.
                      (Input)


    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
*/
#ifdef ANSI
static Mfloat *l_r4dbf (Mfloat (*basis_function) (Mfloat),
                Mint *dimension, Mint *num_points, Mint *num_centers,
                Mfloat *delta, Mfloat *abscissae, Mfloat *centers,
                Mfloat weights[], Mint *additional_columns)
#else
static Mfloat *l_r4dbf (basis_function, dimension, num_points,
                num_centers, delta, abscissae,
                centers, weights, additional_columns)
/*    Mf_radial_fcn     (*basis_function;*/
    Mf_radial_fcn basis_function;
    Mint       *dimension;
    Mint       *num_points;
    Mint       *num_centers;
    Mfloat     *delta;
    Mfloat     *abscissae;
    Mfloat     *centers;
    Mfloat      weights[];
    Mint       *additional_columns;
#endif
{
    Mint        i;
    Mint        k;
    Mint        j;
    Mint        num_rows;
    Mint        num_columns;
    Mfloat     *least_squares_matrix = NULL;
    Mfloat      distance;
    Mfloat      weight_factor;
    Mfloat     *abscissa_to_use = NULL;
    imsl_e1psh ("r4dbf");

    /* Get the space for the retrned least_squares_matrix. */
    least_squares_matrix = (Mfloat *) imsl_malloc (
	(*num_centers + *additional_columns) ** num_points * sizeof (*least_squares_matrix));
    if ((least_squares_matrix == NULL)) {
	imsl_e1stl (1, "num_centers");
	imsl_e1sti (1, *num_centers);
	imsl_e1stl (2, "num_points");
	imsl_e1sti (2, *num_points);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	goto FREE_SPACE;
    }
    num_rows = *num_points;
    num_columns = *num_centers + *additional_columns;



    for (j = 0; j < *num_points; ++j) {
	/* Each row uses the same abscissa */
	abscissa_to_use = &(abscissae[*dimension * j]);
	for (i = 0; i < *num_centers; ++i) {
	    weight_factor = sqrt (weights[j]);
	    /* Loop through all centers for each row. */
	    distance = F_ZERO;
	    /* Compute the distance between the current abscissa and center */
	    for (k = 0; k < *dimension; ++k)
		*(work_basis + k) =
		    *(abscissa_to_use + k) -
		    *(centers + i ** dimension + k);
	    distance = imsl_snrm2 (*dimension, work_basis, 1);
	    /* Fill least_squares_matrix[i][j] */
	    least_squares_matrix[j * num_columns + i] = weight_factor *
		basis_function (
		distance);

	}
    }

    /* Compute the additional columns of H if needed.   */
    /*
     * If linear or constant term is desired, then the first extra column
     * will be filled with 1.0*sqrt(weight).
     */
    if (*additional_columns >= 1) {
	for (j = 0; j < *num_points; ++j) {
	    least_squares_matrix[*num_centers * (j + 1)] = sqrt (weights[j]);
	}
    }

    /*
     * If a linear term is desired, we have to fill out ndim more columns.
     */
    if (*additional_columns > 1) {
	for (i = 0; i < *dimension; ++i) {
	    for (j = 0; j < *num_points; ++j) {
		least_squares_matrix[*num_centers * (j + 1) + i + 1] = sqrt (weights[j]) *
		    abscissae[*dimension * j + i];
	    }

	}
    }
    goto RETURN;
FREE_SPACE:
    if (least_squares_matrix != NULL) {
	imsl_free (least_squares_matrix);
	least_squares_matrix = NULL;
    }
RETURN:
    imsl_e1pop ("r4dbf");
    return (least_squares_matrix);

}
