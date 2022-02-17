#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#undef FALSE

#define FALSE 0

static Mint *orig_itype = NULL;



static VA_LIST_HACK PROTO(l_spline_lsq_constrained, (Mint ndata, Mfloat xdata[],

	Mfloat fdata[], Mint ncoef, Mint num_con_pts, Mf_constraint_struct

	constraints[], va_list argptr));

static void PROTO (l_c2nft, (Mint *ndata, Mfloat xdata[], Mfloat fdata[],

	            Mfloat weight[], Mint *nxval, Mfloat xval[],

	            Mint *nhard, Mint *ider, Mint *itype, Mfloat bl[],

	            Mfloat bu[], Mint *korder, Mfloat xknot[],

	            Mint *ncoef, Mfloat bscoef[], Mfloat h[],

	            Mfloat g[], Mfloat a[], Mfloat rhs[], Mfloat wk[],

	            Mint *iperm, Mint *iwk));

    static void PROTO (l_c3nft, (Mint *nxval, Mfloat xval[], Mint itype[],

	            Mint ider[], Mfloat bl[], Mfloat bu[], Mfloat xknot[],

	            Mint *ncoef, Mint *korder, Mint iperm[]));

    static void PROTO (l_c4nft, (Mint *ndata, Mfloat xdata[], Mfloat fdata[],

	            Mfloat weight[], Mint *korder, Mfloat xknot[],

	            Mint *ncoef, Mfloat *h, Mfloat g[], Mfloat biatx[],

	            Mfloat wk1[], Mfloat wk2[]));

    static void PROTO (l_c5nft, (Mint *ndata, Mfloat xdata[], Mint *nxval,

	            Mfloat xval[], Mint *ider, Mint *itype, Mfloat bl[],

	            Mfloat bu[], Mint *korder, Mfloat xknot[], Mint *ncoef,

	            Mint iperm[], Mfloat *a, Mint *lda, Mfloat rhs[],

	            Mint *irow, Mint *neq, Mfloat wk1[],

	            Mfloat wk2[], Mfloat wk3[], Mfloat wk4[]));

    static void PROTO (l_c6nft, (Mfloat *x, Mint *nxval, Mint *korder,

	            Mint *ider, Mfloat xknot[], Mint *ncoef, Mint *iflag,

	            Mfloat *a, Mint *lda, Mfloat rhs[], Mint *irow, Mfloat *bound,

	            Mfloat *wka, Mfloat *dbiatx, Mint *iconst));

    static void PROTO (l_c7nft, (Mfloat *s, Mfloat *t, Mint *nxval, Mint *ider,

	            Mint *korder, Mfloat xknot[], Mint *ncoef, Mint *iflag,

	            Mfloat *a, Mint *lda, Mfloat rhs[], Mint *irow, Mfloat *bound,

	            Mint *iconst, Mfloat wk1[], Mfloat wk2[], Mfloat wk3[],

	            Mfloat wk4[]));

    static void PROTO (l_c8nft, (Mint *nxval, Mint *nhard, Mint ider[],

	            Mint itype[], Mint *korder, Mint *ncoef,

	            Mfloat bscoef[], Mfloat *h, Mfloat g[], Mfloat *a, Mint *lda,

	            Mfloat rhs[], Mint iperm[], Mint *neq, Mint *irow, Mfloat wk1[],

	            Mfloat wk2[], Mfloat wk3[], Mint iwk1[], Mint iwk2[]));

    static void PROTO (l_c9nft, (Mint *ndata, Mfloat xdata[], Mint *korder,

	            Mfloat xknot[], Mfloat weight[], Mint *ncoef,

	            Mint *scalar));

    static void PROTO (l_c10ft, (Mint *m, Mint *nvar, Mfloat *a, Mint *lda,

	            Mfloat bl[], Mfloat bu[], Mint itype[], Mfloat xlb[],

	            Mfloat xub[], Mint *iflag, Mfloat wk[], Mint iwk[]));

    static void PROTO (l_permu, (Mint *n, Mfloat *x, Mint *ipermu, Mint *ipath,

	            Mfloat *xpermu));

static void PROTO( l_b32gd, (Mfloat t[], Mint *k, Mfloat *x, Mint *left,

                                Mfloat *a, Mfloat *dbiatx, Mint *nderiv));

static void PROTO( l_b42gd, (Mfloat t[], Mint *jhigh, Mint *index, Mfloat *x,

                   Mint * left, Mfloat biatx[]));



static Mf_spline *pp = NULL;



#ifdef ANSI

Mf_spline *imsl_f_spline_lsq_constrained(Mint ndata, Mfloat xdata[],

	Mfloat fdata[], Mint ncoef, Mint num_con_pts,

	Mf_constraint_struct constraints[], ...)

#else

Mf_spline *imsl_f_spline_lsq_constrained(ndata, xdata, fdata, ncoef,

	num_con_pts, constraints, va_alist)

	Mint ndata;

	Mfloat xdata[];

	Mfloat fdata[];

	Mint ncoef;

	Mint num_con_pts;

	Mf_constraint_struct constraints[];

	va_dcl

#endif



{

	va_list	argptr;

	VA_START (argptr, constraints);

	E1PSH ("imsl_f_spline_lsq_constrained", 

		"imsl_d_spline_lsq_constrained");

	pp = NULL;

	IMSL_CALL (l_spline_lsq_constrained (ndata, xdata, fdata, ncoef,

		num_con_pts, constraints, argptr));

	va_end (argptr);

	E1POP ("imsl_f_spline_lsq_constrained", 

		"imsl_d_spline_lsq_constrained");

	return pp;

}





#ifdef ANSI

static VA_LIST_HACK l_spline_lsq_constrained (Mint ndata, Mfloat xdata[],

	Mfloat fdata[], Mint ncoef, Mint num_con_pts, Mf_constraint_struct

	constraints[], va_list argptr)

#else

static VA_LIST_HACK l_spline_lsq_constrained (ndata, xdata, fdata, ncoef,

	num_con_pts, constraints, argptr)

	Mint 	ndata;

	Mfloat	xdata[];

	Mfloat	fdata[];

	Mint	ncoef;

	Mint	num_con_pts;

	Mf_constraint_struct	constraints[];

	va_list	argptr;

#endif

{

	Mint	arg_number = 5;

	Mint	four = 4;

	Mint	code;

	Mint 	domain_dim = 1;

	Mint	target_dim = 1;

	Mint	order_given = 0;

	Mint	knots_given = 0;

	Mint	weight_given = 0;

	Mfloat	*weight = NULL;

	Mfloat	*users_knots;

	Mfloat	x_small;

	Mfloat	x_big;

	Mfloat	interval_size;

	Mint	users_order;

	Mint	*order;

	Mint	korder;

	Mint	nxval = num_con_pts;

	Mint	i;

	Mint	*num_coefs;

	Mint	nhard = 0;

	Mint	free_the_structure = 0;



	Mfloat	*xval	= NULL;

	Mint	*ider	= NULL;

	Mint	*itype	= NULL;

	Mfloat	*bl	= NULL;

	Mfloat	*bu	= NULL;



	Mfloat	*h	= NULL;

	Mfloat	*g	= NULL;

	Mfloat 	*a	= NULL;

	Mfloat	*rhs	= NULL;

	Mfloat	*wk	= NULL;

	Mint	*iperm	= NULL;

	Mint	*iwk	= NULL;

	Mint	length_of_orig_itype;

	Mint	length_of_h;

	Mint	length_of_g;

	Mint	length_of_a;

	Mint	length_of_rhs;

	Mint	length_of_wk;

	Mint	length_of_iperm;

	Mint	length_of_iwk;











	code = 1;

	while (code > 0) {

		code = va_arg(argptr, Mint);

		arg_number++;

		switch (code) {

			case IMSL_ORDER:

				order_given = 1;

				users_order = va_arg (argptr, Mint);

				arg_number++;

				break;

			case IMSL_KNOTS:

				knots_given = 1;

				users_knots = va_arg (argptr, Mfloat*);

				arg_number++;

				break;

			case IMSL_WEIGHTS:

				weight_given = 1;

				weight = va_arg (argptr, Mfloat*);

				arg_number++;

				break;

			case IMSL_NHARD:

				nhard = va_arg (argptr, Mint);

				arg_number++;

				break;

			case 0:

				break;

			default:

				imsl_e1sti (1, code);

				imsl_e1sti (2, arg_number);

				imsl_ermes (IMSL_TERMINAL,IMSL_UNKNOWN_OPTION);

				goto RETURN;

			

		}

	}



	if (!order_given)

		order = &four;

	else

		order = &users_order;

	if (*order < 1) {

		imsl_e1sti (1, *order);

		imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);

		goto RETURN;

	}

    if (num_con_pts < 0) {

        imsl_e1sti (1, num_con_pts);

        imsl_ermes (IMSL_TERMINAL, IMSL_NXVAL_POSITIVE);

        goto RETURN;

    }



/* check ndata */

	if (ndata < *order) {

		imsl_e1sti (1, ndata);

		imsl_e1sti (2, *order);

		imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_NEED_DATA_PTS);

		goto RETURN;

	}



	num_coefs = &ncoef;

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



/* create the structure */

	if (!knots_given)

		pp = imsl_f_spline_create (domain_dim, target_dim, order,

			num_coefs, 0);

	else

		pp = imsl_f_spline_create (domain_dim, target_dim, order,

			num_coefs, IMSL_KNOTS, &users_knots, 0);

	

	if (imsl_n1rty(1) == 4) {

		imsl_e1mes (0, 0, " ");

		imsl_e1stl (1, "ndata");

		imsl_e1sti (1, ndata);

		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

		goto RETURN;

	}



/*

   compute the knots if user did not supply them. these knots are

   equally spaced in the interval and stacked appropriately at the

   endpoints 

 */



	if (!knots_given) {

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

			pp->knots[0][i + *order - 2] = x_small + 

				interval_size * ((Mfloat) (i - 1) / 

				(Mfloat) (*num_coefs - *order + 1));

		}

		pp->knots[0][*num_coefs] += 0.001;



/* stack knots */



		for (i = 1; i <= (*order - 1); i++) {

			pp->knots[0][i - 1] = pp->knots[0][*order - 1];

			pp->knots[0][i + *num_coefs] = pp->knots[0][*num_coefs];

		}

	}



/* get necessary space for c2nft */



	korder = *order;





	length_of_orig_itype = nxval;

	length_of_h = ncoef*ncoef;

	length_of_g = ncoef;

	length_of_a = (2*num_con_pts+korder) * (ncoef+1);

	length_of_rhs = 2*num_con_pts + korder;

	length_of_wk = (korder+1)*(2*korder+1) + (3*ncoef*ncoef+13*ncoef)/2

			+ 2*num_con_pts + 3*ndata + korder 

			+ (2*num_con_pts+korder)*(2*num_con_pts + korder

			+ 29) + 4*ncoef + 1;

	length_of_iperm = num_con_pts;

	length_of_iwk = ndata + 30*(2*num_con_pts+korder) + 4*ncoef;



	if (length_of_orig_itype == 0) length_of_orig_itype = 1;

	orig_itype = (Mint*) imsl_malloc (length_of_orig_itype*sizeof(*orig_itype));

	h = (Mfloat*) imsl_malloc (length_of_h*sizeof(*h));

	g = (Mfloat*) imsl_malloc (length_of_g*sizeof(*g));

	a = (Mfloat*) imsl_malloc (length_of_a*sizeof(*a));

	rhs = (Mfloat*) imsl_malloc (length_of_rhs*sizeof(*rhs));

	wk = (Mfloat*) imsl_malloc (length_of_wk*sizeof(*wk));

	if (length_of_iperm == 0) length_of_iperm = 1;

	iperm = (Mint*) imsl_malloc (length_of_iperm*sizeof(*iperm));

	iwk = (Mint*) imsl_malloc (length_of_iwk*sizeof(*iwk));

	if (!weight_given)

		weight = (Mfloat*) imsl_malloc (ndata*sizeof(*weight));



	if ( (h == NULL) || (g == NULL) || (a == NULL) || (rhs == NULL) ||

		(wk == NULL) || (iperm == NULL) || (iwk == NULL) ||

		(weight == NULL) ) {

		free_the_structure = 1;

		imsl_e1stl (1, "ndata");

		imsl_e1sti (1, ndata);

		imsl_e1stl (2, "num_con_pts");

		imsl_e1sti (2, num_con_pts);

		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);

		goto RETURN;

	}



	nxval = num_con_pts;

	if (nxval == 0) nxval = 1;

	xval = (Mfloat*) imsl_malloc (nxval*sizeof(*xval));

	ider = (Mint*) imsl_malloc (nxval*sizeof(*ider));

	itype = (Mint*) imsl_malloc (nxval*sizeof(*itype));

	bl = (Mfloat*) imsl_malloc (nxval*sizeof(*bl));

	bu = (Mfloat*) imsl_malloc (nxval*sizeof(*bu));

	nxval = num_con_pts;



	if ( (xval == NULL) || (ider == NULL) || (itype == NULL) ||

		(bl == NULL) || (bu == NULL) ) {

		free_the_structure = 1;

		imsl_e1stl (1, "nxval");

		imsl_e1sti (1, nxval);

		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

		goto RETURN;

	}

	if (imsl_n1rty(1))

		goto RETURN;



/* fill c2nft input vectors */



	if (!weight_given)

		for (i=0; i<ndata; i++) 

			*(weight+i) = 1.0;



	for (i=0; i<nxval; i++) {

		*(xval+i) = (constraints+i)->xval;

		*(ider+i) = (constraints+i)->der;

		*(itype+i) = (constraints+i)->type;

                *(orig_itype+i) = (constraints+i)->type;

		*(bl+i) = (constraints+i)->bl;

		*(bu+i) = (constraints+i)->bu;

	}

	for (i = 1; i <= nxval; ++i)

	  {

	/* Check ITYPE */

	if ((((itype[i - 1]) < 1) || (abs (itype[i - 1]) > 12)) &&

	    ((itype[i - 1] != 20) && (itype[i - 1] != 99))) {

	    imsl_e1sti (1, orig_itype[i - 1]);

	    imsl_e1sti (2, i - 1);

            imsl_ermes (IMSL_TERMINAL, IMSL_BAD_CNSTR_TYPE);

	    goto RETURN;

	}

	    

	  }



	for (i=0; i<nxval; i++) {

		if ( (constraints+i)->type == 5) {

			*(itype+i) = 1;

			*(ider+i) = -1;

		}

		if ( (constraints+i)->type == 6) {

			*(itype+i) = 2;

 			*(ider+i) = -1;

		}

		if ( (constraints+i)->type == 7) {

			*(itype+i) = 3;

			*(ider+i) = -1;

		}

		if ( (constraints+i)->type == 8) {

			*(itype+i) = 4;

			*(ider+i) = -1;

		}



		if ( (constraints+i)->type == 9)

			*(itype+i) = -1;



		if ( (constraints+i)->type == 10)

			*(itype+i) = -2;



		if ( (constraints+i)->type == 11)

			*(itype+i) = -3;



		if ( (constraints+i)->type == 12) 

			*(itype+i) = -4;



		if ( (constraints+i)->type == 20)

			*(itype+i) = 10;

	}	



/* do the work */





	l_c2nft (&ndata, xdata, fdata, weight, &nxval, xval,

		&nhard, ider, itype, bl, bu, pp->order, pp->knots[0],

		&ncoef, pp->coef[0], h, g, a, rhs, wk, iperm, iwk);



RETURN:



	if (free_the_structure && pp != NULL) {

		imsl_free(pp);

		pp = NULL;

	}

	if (orig_itype != NULL) { imsl_free(orig_itype); orig_itype = NULL;}

	if (h != NULL) imsl_free(h);

	if (g != NULL) imsl_free(g);

	if (a != NULL) imsl_free(a);

	if (rhs != NULL) imsl_free(rhs);

	if (wk != NULL) imsl_free(wk);

	if (iperm != NULL) imsl_free(iperm);

	if (iwk != NULL) imsl_free(iwk);

	if (!weight_given && weight != NULL) imsl_free (weight);

	if (xval != NULL) imsl_free(xval);

	if (ider != NULL) imsl_free(ider);

	if (itype != NULL) imsl_free(itype);

	if (bl != NULL) imsl_free(bl);

	if (bu != NULL) imsl_free(bu);



	return argptr;

}





/*  -----------------------------------------------------------------------

    IMSL Name:  C2NFT/DC2NFT (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    March 1, 1990



    Purpose:    Compute the least squares constrained spline

                approximation, returning the B-spline coefficients.



    Usage:      CALL C2NFT (NDATA, XDATA, FDATA, WEIGHT, NXVAL,

                            XVAL, NHARD, IDER, ITYPE, BL, BU,

                            KORDER, XKNOT, NCOEF, BSCOEF, H, G, A,

                            RHS, WK, IPERM, IWK)



    Arguments:  (See CONFT)



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_c2nft (Mint *ndata, Mfloat xdata[], Mfloat fdata[],

                Mfloat weight[], Mint *nxval, Mfloat xval[],

                Mint *nhard, Mint *ider, Mint *itype, Mfloat bl[],

                Mfloat bu[], Mint *korder, Mfloat xknot[],

                Mint *ncoef, Mfloat bscoef[], Mfloat h[],

                Mfloat g[], Mfloat a[], Mfloat rhs[], Mfloat wk[],

                Mint *iperm, Mint *iwk)

#else

static void l_c2nft (ndata, xdata, fdata, weight, nxval, xval,

                nhard, ider, itype, bl, bu, korder, xknot,

                ncoef, bscoef, h, g, a, rhs, wk, iperm, iwk)

    Mint       *ndata;

    Mfloat      xdata[], fdata[], weight[];

    Mint       *nxval;

    Mfloat      xval[];

    Mint       *nhard, ider[], itype[];

    Mfloat      bl[], bu[];

    Mint       *korder;

    Mfloat      xknot[];

    Mint       *ncoef;

    Mfloat      bscoef[], h[], g[], a[], rhs[], wk[];

    Mint        iperm[], iwk[];

#endif

{

    Mint        _l0, i, irow, kwk1, kwk2, kwk3, kwk4, kwk5, kwk6, kwk7,

                kwk8, kwk9, lda, neq;





    imsl_e1psh ("C2NFT");



    /* Check NXVAL */

    if (*nxval < 0) {

	imsl_e1sti (1, *nxval);

        imsl_ermes (IMSL_TERMINAL, IMSL_NXVAL_POSITIVE);

	goto L_9000;

    }

    /* Check NXVAL .ge. NHARD */

    if (*nxval < *nhard) {

	imsl_e1sti (1, *nxval);

	imsl_e1sti (2, *nhard);

        imsl_ermes (IMSL_TERMINAL, IMSL_NHARD_GT_TOTAL);

    }

    /* Check arguments. */

    _l0 = FALSE;

    l_c9nft (ndata, xdata, korder, xknot, weight, ncoef, &_l0);

    /* Check for errors */

    if (imsl_n1rty (0) != 0)

	goto L_9000;

    /*

     * Partition workspace

     */

    kwk1 = 1;

    kwk2 = kwk1 + *korder + 1;

    kwk3 = kwk2 + (*korder + 1) ** korder;

    kwk4 = kwk3 + (*korder + 1) ** korder;

    kwk5 = kwk4 + *ncoef;

    kwk6 = kwk5 + (3 ** ncoef ** ncoef + 11 ** ncoef) / 2 + 2 ** nxval + *korder;

    kwk7 = kwk6 + *ndata;

    kwk8 = kwk7 + *ndata;

    kwk9 = kwk8 + *ndata;

    /*

     * If NXVAL is equal to zero then call B2LSQ for the unconstrained

     * least-squares fit.

     */

    if (*nxval == 0) {

	imsl_b2lsq (ndata, xdata, fdata, weight, korder, xknot, ncoef,

	    bscoef, &wk[kwk5 - 1], &wk[kwk6 - 1], &wk[kwk7 - 1], &wk[kwk8 - 1],

	    iwk);

    }

    else {

	/*

	 * Sort the XDATA array and arrange the FDATA array correspondingly.

	 */

	for (i = 1; i <= *ndata; i++) {

	    iwk[i - 1] = i;

	}

	imsl_svrgp (*ndata, xdata, &wk[kwk6 - 1], iwk);

	_l0 = 1;

	l_permu (ndata, fdata, iwk, &_l0, &wk[kwk7 - 1]);



	/*

	 * Permute input constraints so that the equality constraints come

	 * first. the permutation is stored in iperm. C3NFT also does some

	 * checks on the input constraints.

	 */

	l_c3nft (nxval, xval, itype, ider, bl, bu, xknot, ncoef, korder,

	    iperm);

	if (imsl_n1rty (1) != 0)

	    goto L_9000;

	/*

	 * Compute the matrix H and vector G for the objective function to be

	 * used in the call to Q2ROG.

	 */

	l_c4nft (ndata, &wk[kwk6 - 1], &wk[kwk7 - 1], weight, korder,

	    xknot, ncoef, h, g, &wk[kwk1 - 1], &wk[kwk2 - 1], &wk[kwk3 - 1]);

	/* Initialize A and RHS to zero. */

	sset ((2 ** nxval + *korder) ** ncoef, 0.0, a, 1);

	sset (2 ** nxval + *korder, 0.0, rhs, 1);

	/*

	 * Compute the constraint matrix a and the right hand side, RHS, to

	 * be used in the call to Q2ROG.

	 */

	lda = 2 ** nxval + *korder;

	l_c5nft (ndata, &wk[kwk6 - 1], nxval, xval, ider, itype, bl,

	    bu, korder, xknot, ncoef, iperm, a, &lda, rhs, &irow, &neq,

	    &wk[kwk1 - 1], &wk[kwk5 - 1], &wk[kwk3 - 1], &wk[kwk2 - 1]);

	if (imsl_n1rty (1) != 0)

	    goto L_9000;

	/*

	 * Solve the quadratic programming problem made up of the matrix H

	 * and the vector G for the objective function, the matrix A as the

	 * constraint matrix, and RHS as the right hand side of the system of

	 * constraints.

	 */

	l_c8nft (nxval, nhard, ider, itype, korder, ncoef, bscoef, h,

	    g, a, &lda, rhs, iperm, &neq, &irow, &wk[kwk4 - 1], &wk[kwk5 - 1],

	    &wk[kwk9 - 1], &iwk[0], &iwk[*ncoef]);

    }



L_9000:

    imsl_e1pop ("C2NFT");

    return;

}				/* end of function */

#undef FALSE



/*  -----------------------------------------------------------------------

    IMSL Name:  C10FT/DC10FT (Single/Double precision version)



    Computer:   SUN4/SINGLE



    Revised:    March 1, 1990



    Purpose:    Check if a set of constraints are consistent by calling

                a linear programming routine which uses the revised

                simplex algorithm.



    Usage:      CALL C10FT (M, NVAR, A, LDA, BL, BU, ITYPE, XLB, XUB,

                            IFLAG, WK, IWK)



    Arguments:

       M      - Number of constraints.  (Input)

       NVAR   - Number of variables.  (Input)

       A      - M by NVAR matrix containing the coefficients of the

                M constraints.  (Input)

       LDA    - Leading dimension of A exactly as specified in the

                dimension statement of the calling program.  (Input)

       BL     - Vector of length M containing the lower limit of the

                general constraints; if there is no lower limit on the

                I-th constraint, then BL(I) is not referenced.

                (Input)

       BU     - Vector of length M containing the upper limit of the

                general constraints; if there is no upper limit on

                the I-th constraint, then BU(I) is not referenced;

                if there is no range constraint, BL and BU can share

                the same storage locations.  (Input)

       ITYPE  - Vector of length M indicating the types of general

                constraints in the matrix A.  (Input)

                Let

                   R(I) = A(I,1)*XSOL(1) + ... + A(I,NVAR)*XSOL(NVAR).

                Then the value of IRTYPE(I) signifies the following.

                   IRTYPE(I)           I-th CONSTRAINT

                  -----------    ----------------------------

                       0          BL(I) .EQ. R(I) .EQ. BU(I)

                       1          R(I) .LE. BU(I)

                       2          R(I) .GE. BL(I)

                       3          BL(I) .LE. R(I) .LE. BU(I)

       XLB    - Vector of length NVAR containing the lower bound on the

                variables; if there is no lower bound on a variable,

                then 1.0E30 should be set as the lower bound.  (Input)

       XUB    - Vector of length NVAR containing the upper bound on the

                variables; if there is no upper bound on a variable,

                then -1.0E30 should be set as the upper bound.  (Input)

       IFLAG  - Iflag is set to 99 if the constaints were found to be

                inconsistent, and one otherwise.  (Output)

       WK     - Real ork array of size M*(M+29)+2*NVAR+1.

       IWK    - Integer work array of size 29*M+3*NVAR.



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_c10ft (Mint *m, Mint *nvar, Mfloat *a, Mint *lda,

                Mfloat bl[], Mfloat bu[], Mint itype[], Mfloat xlb[],

                Mfloat xub[], Mint *iflag, Mfloat wk[], Mint iwk[])

#else

static void l_c10ft (m, nvar, a, lda, bl, bu, itype, xlb, xub,

                iflag, wk, iwk)

    Mint       *m, *nvar;

    Mfloat     *a;

    Mint       *lda;

    Mfloat      bl[], bu[];

    Mint        itype[];

    Mfloat      xlb[], xub[];

    Mint       *iflag;

    Mfloat      wk[];

    Mint        iwk[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

    Mint        kr1, kr2, kr3, kr4, kr5;

    Mfloat      obj;

    Mint	*ibasis = NULL;

    Mint	*ibb = NULL;







    imsl_e1psh ("l_c10ft");

    /* Compute workspace offsets. */

    kr1 = 1;

    kr2 = kr1 + *nvar;

    kr3 = kr2 + *m;

    kr4 = kr3 + 1;

    kr5 = kr4 + *m * (*m + 28);

    /*

     * Set the coefficients of the objective function sent to D2PRS to zero.

     * This will cause D2PRS to stop after phase one.

     */

    sset (*nvar, 0.0e0, &wk[kr5 - 1], 1);



    /* call D2PRS. */



#if 0

    imsl_d2prs (*m, *nvar, a, *lda, bl, bu, &wk[kr5 - 1], itype, xlb, xub,

	&obj, &wk[kr1 - 1], &wk[kr2 - 1], &wk[kr4 - 1],

	iwk, 10000);

#endif



        ibasis = (Mint *) imsl_malloc((*m + *nvar)*sizeof(*ibasis));

        ibb = (Mint *) imsl_malloc((*m + *nvar)*sizeof(*ibb));

	if (ibasis == NULL || ibb == NULL) {

                imsl_e1stl (1, "m");

                imsl_e1sti (1, *m);

                imsl_e1stl (2, "nvar");

                imsl_e1sti (2, *nvar);

                imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);

        }



    imsl_d2prs (*m, *nvar, a, *lda, bl, bu, &wk[kr5 - 1], itype, xlb, xub,

	0, ibasis, ibb,

	&obj, &wk[kr1 - 1], &wk[kr2 - 1], &wk[kr4 - 1],

	iwk, 10000);



	if (ibasis != NULL) imsl_free (ibasis);

	if (ibb != NULL) imsl_free (ibb);



    /*

     * See if an error state has been set in D2PRS.  If the constraints were

     * found to be inconsistent, or made the problem infeasible,  set IFLAG

     * to 99, and remove the error state. Otherwise, set IFLAG to 1 and do

     * not remove the error state.

     */

/*    if (((imsl_n1rty (1) == 4) && (imsl_n1rcd (1) == 5)) || ((imsl_n1rty (1) ==

		3) && (imsl_n1rcd (1) == 3))) {*/

    if ((imsl_error_code () == IMSL_BOUNDS_INCONSISTENT) || (imsl_error_code () == IMSL_PROB_INFEASIBLE)) {

	*iflag = 99;

	imsl_e1mes (0, 0, " ");

    }

    else {

	*iflag = 1;

    }

    imsl_e1pop ("l_c10ft");

    return;

}				/* end of function */



#undef A



/*  -----------------------------------------------------------------------

    IMSL Name:  C3NFT/DC3NFT (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    March 1, 1990



    Purpose:    Compute a permutation vector that arranges constraints

                for CONFT so that the equality constraints are first and

                make some checks on the constraints.



    Usage:      CALL C3NFT (NXVAL, XVAL, ITYPE, IDER, BL, BU,

                            XKNOT, NCOEF, KORDER, IPERM)



    Arguments:

       NXVAL  - Number of elements to be permuted.  (Input)

       XVAL   - Array of length NXVAL containing the abscissas at

                which the fit is to be constrained.  (Input)

       ITYPE  - Array of length NXVAL containing the types of

                constraints to be permuted.  (Input)

       IDER   - Array of length NXVAL containing the derivative

                value of the spline  which is to be constrained.  (Input)

       BL     - Array of length NXVAL containing the lower limit of the

                general constraints; if there is no lower limit on the

                I-th constraint, then BL(I) is not referenced.

                (Input)

       BU     - Array of length NXVAL containing the upper limit of the

                general constraints; if there is no upper limit on

                the I-th constraint, then BU(I) is not referenced;

                if there is no range constraint, BL and BU can share

                the same storage locations.  (Input)

       XKNOT  - Array of length NCOEF+KORDER containing the knot

                sequence.  (Input)

                XKNOT must be nondecreasing.

       NCOEF  - Number of B-spline coefficients.  (Input)

       KORDER - Order of the spline.  (Input)

       IPERM  - Array of length NXVAL containing the permutation

                vector.  (Output)



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_c3nft (Mint *nxval, Mfloat xval[], Mint itype[],

                Mint ider[], Mfloat bl[], Mfloat bu[], Mfloat xknot[],

                Mint *ncoef, Mint *korder, Mint iperm[])

#else

static void l_c3nft (nxval, xval, itype, ider, bl, bu, xknot, ncoef,

                korder, iperm)

    Mint       *nxval;

    Mfloat      xval[];

    Mint        itype[], ider[];

    Mfloat      bl[], bu[], xknot[];

    Mint       *ncoef, *korder, iperm[];

#endif

{

    Mint        i, j, k, l, tmp1, tmp2;





    imsl_e1psh ("C3NFT");

    /*

     * J is used as a counter of the number of equality constraints.

     */

    j = 1;

    /*

     * This next line protects against IPERM(1) being undefined when it is

     * accessed later in this routine. It is set to -1 sothare is no way to

     * confuse it with the index of a constraint.

     */

    if (*nxval >= 1)

	iperm[0] = -1;

    /*

     * The following DO loop makes some checks on all of the constraints. It

     * also puts the indices of all the equality constraints in the first

     * entries of IPERM.

     */

    for (i = 1; i <= *nxval; i++) {

	/*

	 * Make sure XVAL(I) is within the given knot sequence.

	 */

	if ((xval[i - 1] < xknot[*korder - 1]) || (xval[i - 1] > xknot[*ncoef])) {

	    imsl_e1sti (1, i-1);

	    imsl_e1sti (2, *ncoef);

	    imsl_e1str (1, xknot[*korder - 1]);

	    imsl_e1str (2, xval[i - 1]);

	    imsl_e1str (3, xknot[*ncoef]);

            imsl_ermes (IMSL_FATAL, IMSL_XVAL_WITHIN_KNOTS);

	    goto L_9000;

	}

	/* Check IDER */

	if (ider[i - 1] < -1) {

	    imsl_e1sti (1, i - 1);

	    imsl_e1sti (2, ider[i - 1]);

            imsl_ermes (IMSL_TERMINAL, IMSL_DER_TOO_SMALL);

	    goto L_9000;

	}

	/* Check ITYPE */

	if (((abs (itype[i - 1]) < 1) || (abs (itype[i - 1]) > 4)) &&

	    ((itype[i - 1] != 10) && (itype[i - 1] != 99))) {

	    imsl_e1sti (1, orig_itype[i - 1]);

	    imsl_e1sti (2, i - 1);

            imsl_ermes (IMSL_TERMINAL, IMSL_BAD_CNSTR_TYPE);

	    goto L_9000;

	}

	else {

	    /*

	     * Place the index, w.r.t. itype, of the constraints which

	     * involve equality constaints. ITYPE equal to 1 or 10 implies an

	     * equality constraint.

	     */

	    if ((abs (itype[i - 1]) == 1) || (abs (itype[i - 1]) ==

		    10)) {

		iperm[j - 1] = i;

		j += 1;

	    }

	}

    }

    /*

     * Make more checks on costraints. The checks are arranged as follows: 1.

     * (BL .lt. BU) for range constraints 2. (IDER(i) .ge. 0) if ITYPE(i)

     * .eq. 10 3. Five checks done in case IDER(i).eq.-1. 4. Three checks

     * done if ITYPE(i) .eq. -1, -2, -3, or -4

     */

    i = 1;

L_20:

    if (i <= *nxval) {

	/*

	 * Check BL .lt. BU for range constraints.

	 */

	if (abs (itype[i - 1]) == 4) {

	    if (bl[i - 1] > bu[i - 1]) {

		imsl_e1str (1, bl[i - 1]);

		imsl_e1str (2, bu[i - 1]);

		imsl_e1sti (3, i - 1);

                imsl_ermes(IMSL_FATAL, IMSL_BAD_RANGE);

		goto L_9000;

	    }

	}

	/* Check IDER for ITYPE(I) = 10 */

	if ((itype[i - 1] == 10) && (ider[i - 1] < 0)) {

	    imsl_e1sti (1, i - 1);

	    imsl_e1sti (2, 20);

	    imsl_e1sti (3, ider[i - 1]);

            imsl_ermes (IMSL_TERMINAL, IMSL_BAD_PERIODIC_DER);

	    goto L_9000;

	}

	else if (itype[i - 1] == 10) {

	    i += 1;

	    goto L_20;

	}

	/* Check constraints on integrals. */

	if (ider[i - 1] == -1) {

	    if (i == *nxval) {

		imsl_e1sti (1, i - 1);

                imsl_ermes(IMSL_TERMINAL, IMSL_HALF_CONSTRAINT);

		goto L_9000;

	    }

	    /*

	     * Give error message in the case ITYPE(i) .ne. ITYPE(i+1) when

	     * it is required that they be equal.

	     */

	    if (itype[i - 1] != itype[i]) {

                tmp1 = orig_itype[i-1];

                tmp2 = orig_itype[i];

		imsl_e1sti (1, tmp1);

		imsl_e1sti (2, tmp2);

		imsl_e1sti (3, i - 1);

		imsl_e1sti (4, i);

                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_TYPE_PAIR);

		goto L_9000;

	    }

	    /* ITYPE must 1,2,3, or 4. */

	    if ((itype[i - 1] < 0) || (itype[i - 1] == 10)) {

		imsl_e1sti (1, i - 1);

		imsl_e1sti (2, itype[i - 1]);

                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_INTGR_CNSTR);

		goto L_9000;

	    }

	    /*

	     * Check to see that the left endpoint is less than or equal to

	     * the right end point.

	     */

	    if (xval[i - 1] >= xval[i]) {

		imsl_e1sti (1, i - 1);

		imsl_e1sti (2, i);

		imsl_e1str (1, xval[i - 1]);

		imsl_e1str (2, xval[i]);

                imsl_ermes(IMSL_FATAL, IMSL_BAD_INTGR_RANGE);

		goto L_9000;

	    }

	    /* Check IDER(I)=IDER(I+1)=-1 */

	    if (ider[i] != -1) {

		imsl_e1sti (1, i - 1);

		imsl_e1sti (2, ider[i - 1]);

		imsl_e1sti (3, i);

		imsl_e1sti (4, ider[i]);

                imsl_ermes(IMSL_FATAL, IMSL_BAD_INTGR_DER);

		goto L_9000;

	    }

	    /*

	     * If it makes it this far then the integral constraint appears

	     * to be O.K.

	     */

	    i += 2;

	    goto L_20;

	    /* Check constraints for ITYPE .lt. 0 */

	}

	else if (itype[i - 1] < 0) {

	    if (i == *nxval) {

		imsl_e1sti (1, i - 1);

                imsl_ermes ( IMSL_TERMINAL, IMSL_ONLY_HALF_CONSTR);

		goto L_9000;

	    }

	    /*

	     * Give error message in the case ITYPE(i) .ne. ITYPE(i+1) when

	     * it is required that they be equal.

	     */

	    if (itype[i - 1] != itype[i]) {

                tmp1 = orig_itype[i-1];

                tmp2 = orig_itype[i];

		imsl_e1sti (1, tmp1);

		imsl_e1sti (2, tmp2);

		imsl_e1sti (3, i - 1);

		imsl_e1sti (4, i);

                imsl_ermes (IMSL_TERMINAL, IMSL_BAD_TYPE_PAIR);

		goto L_9000;

	    }

	    /* Check IDER(I+1) .ge.0 */

	    if (ider[i] < 0) {

                tmp1 = orig_itype[i - 1];

		imsl_e1sti (1, i - 1);

		imsl_e1sti (4, ider[i]);

		imsl_e1sti (3, i);

		imsl_e1sti (2, tmp1);

                imsl_ermes (IMSL_TERMINAL, IMSL_DER_GE_ZERO);

		goto L_9000;

	    }

	    /*

	     * If it makes it this far then the two point constraint appears

	     * to be O.K.

	     */

	    i += 2;

	    goto L_20;

	}

	else {

	    /*

	     * If it makes to here, the constraint is either an ignored

	     * constraint, or is a one point constraint.  If it is a one

	     * point constraint, the necessary checks were done in the above

	     * DO loop.

	     */

	    i += 1;

	    goto L_20;

	}

    }

    k = iperm[0];

    l = 1;

    /*

     * At this point all of the equality constraints are represented in the

     * first (j-1) elements of IPERM. all that are left now are inequality

     * constaints.

     */

    for (i = 1; i <= *nxval; i++) {

	/*

	 * Check if the ith constraint was an equality constraint.  If it

	 * was, then skip to the next one.

	 */

	if (i == k) {

	    k = iperm[l];

	    l += 1;

	    /*

	     * If the ith constraint was not an equality constraint place the

	     * index of the constraint in IPERM.

	     */

	}

	else {

	    iperm[j - 1] = i;

	    j += 1;

	}

    }

L_9000:

    imsl_e1pop ("C3NFT");

    return;

}				/* end of function */



/*  -----------------------------------------------------------------------

    IMSL Name:  C4NFT/DC4NFT (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    March 1, 1990



    Purpose:    Compute the Hessian matrix H and vector G to be used

                by CONFT in a call to QPROG.



    Usage:      CALL C4NFT (NDATA, XDATA, FDATA, WEIGHT,

                            KORDER, XKNOT, NCOEF, H, G,

                            BIATX, WK1, WK2)



    Arguments:

       NDATA  - Number of data points.  (Input)

       XDATA  - Array of length NDATA containing the data point

                abscissas.  (Input)

       FDATA  - Array of length NDATA containing the values to be

                approximated.  (Input)

       WEIGHT - Array of length NDATA containing the weights.  (Input)

       KORDER - Order of the spline.  (Input)

       XKNOT  - Array of length NCOEF+KORDER containg the knot

                sequence.  (Input)

       NCOEF  - Number of B-spline coefficients.  (Input)

       H      - Array of size NCOEF by NCOEF containing the Hessian

                matrix of the objective function to be used in a

                call to QPROG.  (Output)

       G      - Array of length NCOEF containing the coefficients of

                the linear term of the objective function used in a

                call to QPROG.  (Output)

       BIATX  - Array of length KORDER used as workspace.

       WK1    - Array of length KORDER used as workspace.

       WK2    - Array of length KORDER used as workspace.



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_c4nft (Mint *ndata, Mfloat xdata[], Mfloat fdata[],

                Mfloat weight[], Mint *korder, Mfloat xknot[],

                Mint *ncoef, Mfloat *h, Mfloat g[], Mfloat biatx[],

                Mfloat wk1[], Mfloat wk2[])

#else

static void l_c4nft (ndata, xdata, fdata, weight, korder, xknot,

                ncoef, h, g, biatx, wk1, wk2)

    Mint       *ndata;

    Mfloat      xdata[], fdata[], weight[];

    Mint       *korder;

    Mfloat      xknot[];

    Mint       *ncoef;

    Mfloat     *h, g[], biatx[], wk1[], wk2[];

#endif

{

#define H(I_,J_)	(h+(I_)*(*ncoef)+(J_))

    Mint        _l0, i, j, k, left, mflag;

    Mfloat      cj;

    /*

     * Zero out the vector G and the matrix H.

     */

    sset (*ncoef, 0.0, g, 1);

    sset (*ncoef ** ncoef, 0.0, h, 1);

    /*

     * Run the outer loop for each data point.

     */

    for (i = 1; i <= *ndata; i++) {

	/*

	 * Find the values of all possibly nonzero B-splines.

	 */

	_l0 = *korder + *ncoef;

	imsl_b4der (xknot, &_l0, &xdata[i - 1], &left,

	    &mflag);

	imsl_b4int (xknot, korder, &xdata[i - 1], &left, biatx, wk1, wk2);

	/*

	 * Do 20 loop imsl_runs only for the KORDER possibly nonzero

	 * B-splines.

	 */

	for (j = left - *korder + 1; j <= left; j++) {

	    cj = weight[i - 1] * biatx[j - left + *korder - 1];

	    g[j - 1] += -fdata[i - 1] * cj;

	    /*

	     * Do 10 loop imsl_runs only for the KORDER possibly nonzero

	     * B-splines.

	     */

	    for (k = left - *korder + 1; k <= left; k++) {

		*H (k - 1, j - 1) += cj * biatx[k - left + *korder - 1];

	    }

	}

    }

    return;

}				/* end of function */



#undef H



/*  -----------------------------------------------------------------------

    IMSL Name:  C5NFT/DC5NFT (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    March 1, 1990



    Purpose:    Generate the rows of the constraint matrix A.



    Usage:      CALL C5NFT (NDATA, XDATA, NXVAL, XVAL,

                            IDER, ITYPE, BL, BU, KORDER, XKNOT,

                            NCOEF, IPERM, A, LDA, RHS, IROW,

                            NEQ, WK1, WK2, WK3, WK4)



    Arguments:

       NDATA  - Number of data points.  (Input)

       XDATA  - Arrray of length NDATA containing the data point

                abscissas.  (Input)

       NXVAL  - Number of points in the vector XVAL.  (Input)

       XVAL   - Array of length NXVAL containing the abscissas at

                which the fit is to constrained.  (Input)

       IDER   - Array of length NXVAL containing the derivative

                value of the spline which is to be fit.  (Input)

       ITYPE  - Array of length NXVAL indicating the types of general

                constraints in the matrix A.  (Input)

       BL     - Vector of length NXVAL containing the lower limit of the

                general constraints; if there is no lower limit on the

                I-th constraint, then BL(I) is not referenced.

                (Input)

       BU     - Vector of length NXVAL containing the upper limit of the

                general constraints; if there is no upper limit on

                the I-th constraint, then BU(I) is not referenced;

                if there is no range constraint, BL and BU can share

                the same storage locations.  (Input)

       KORDER - Order of the spline.  (Input)

       XKNOT  - Array of length NCOEF+KORDER containing the knot

                sequence.  (Input)

                XKNOT must be nondecreasing.

       NCOEF  - Number of B-spline coefficients.  (Input)

       IPERM  - Array of length NXVAL that holds the permutation of

                the constraints.  (Input)

       A      - The constraint matrix.  (Output)

       LDA    - Leading dimension of A exactly as specified in

                the dimension statement of the calling program.  (Input)

                LDA must be at least 2*NXVAL+KORDER.

       RHS    - Array of size at least 2*NXVAL+KORDER that is the

                right hand side of the system of constraints in the

                matrix A.  (Output)

       IROW   - The number of rows in the matrix A upon exit.  (Output)

                IROW will be less than or equal to 2*NXVAL+KORDER.

       NEQ    - The number of equality constraints in the matrix A.

                (Output)

       WK1    - Work array of size KORDER+1.  (Input)

       WK2    - Work array of size (KORDER+1)*KORDER+NCOEF.  (Input)

       WK3    - Work array of size (KORDER+1)*KORDER.  (Input)

       WK4    - Work array of size KORDER+1.  (Input)



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_c5nft (Mint *ndata, Mfloat xdata[], Mint *nxval,

                Mfloat xval[], Mint *ider, Mint *itype, Mfloat bl[],

                Mfloat bu[], Mint *korder, Mfloat xknot[], Mint *ncoef,

                Mint iperm[], Mfloat *a, Mint *lda, Mfloat rhs[],

                Mint *irow, Mint *neq, Mfloat wk1[],

                Mfloat wk2[], Mfloat wk3[], Mfloat wk4[])

#else

static void l_c5nft (ndata, xdata, nxval, xval, ider, itype, bl,

                bu, korder, xknot, ncoef, iperm, a, lda, rhs, irow, neq, wk1,

                wk2, wk3, wk4)

    Mint       *ndata;

    Mfloat      xdata[];

    Mint       *nxval;

    Mfloat      xval[];

    Mint        ider[], itype[];

    Mfloat      bl[], bu[];

    Mint       *korder;

    Mfloat      xknot[];

    Mint       *ncoef, iperm[];

    Mfloat     *a;

    Mint       *lda;

    Mfloat      rhs[];

    Mint       *irow, *neq;

    Mfloat      wk1[], wk2[], wk3[], wk4[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

    Mint        icount, iflag, k;

    Mfloat      _f0, bound;





    imsl_e1psh ("C5NFT");

    /* INITIALIZE RELEVANT COUNTERS */

    *neq = 0;

    *irow = 1;

    icount = 1;

    /*

     * Compute a row(s) of the matrix A. Note that if the constraint is a

     * range constraint, (ITYPE = 4 or ITYPE = -4) then two rows of the

     * matrix A are computed.

     */

L_10:

    if (icount <= *nxval) {

	/* ITYPE equal to 10 */

	if (itype[iperm[icount - 1] - 1] == 10) {

	    for (k = 0; k <= imsl_i_min (ider[iperm[icount - 1] - 1], *korder); k++) {

		iflag = 0;

		_f0 = 0.0;

		l_c6nft (&xdata[0], nxval, korder, &k, xknot, ncoef,

		    &iflag, a, lda, rhs, irow, &_f0, wk2, wk3,

		    &iperm[icount - 1]);

		iflag = 1;

		l_c6nft (&xdata[*ndata - 1], nxval, korder, &k, xknot,

		    ncoef, &iflag, a, lda, rhs, irow, &_f0, wk2,

		    wk3, &iperm[icount - 1]);

		*irow += 1;

		*neq += 1;

	    }

	    icount += 1;

	    goto L_10;

	    /* ITYPE equal to 99 */

	}

	else if (itype[iperm[icount - 1] - 1] == 99) {

	    *A (*ncoef, *irow - 1) = iperm[icount - 1];

	    *irow += 1;

	    icount += 1;

	    goto L_10;

	}

	/*

	 * IDER .eq. -1 implies the constraint is on an integral of the

	 * curve.

	 */

	if (ider[iperm[icount - 1] - 1] == -1) {

	    /* IDER = (-1) and ITYPE = 4 */

	    if (itype[iperm[icount - 1] - 1] == 4) {

		iflag = 0;

		l_c7nft (&xval[iperm[icount - 1] - 1], &xval[iperm[icount] -

			1], nxval, &ider[iperm[icount - 1] - 1], korder,

		    xknot, ncoef, &iflag, a, lda, rhs, irow, &bl[iperm[icount - 1] -

			1], &iperm[icount - 1], wk1, wk4, wk3, wk2);

		*irow += 1;

		iflag = 1;

		l_c7nft (&xval[iperm[icount - 1] - 1], &xval[iperm[icount] -

			1], nxval, &ider[iperm[icount - 1] - 1], korder,

		    xknot, ncoef, &iflag, a, lda, rhs, irow, &bu[iperm[icount - 1] -

			1], &iperm[icount - 1], wk1, wk4, wk3, wk2);

		*irow += 1;

		icount += 2;

		goto L_10;

	    }

	    else {

		/* IDER = (-1) and ITYPE = 1, 2, or 3 */

		if (itype[iperm[icount - 1] - 1] == 2) {

		    iflag = 1;

		    bound = bu[iperm[icount - 1] - 1];

		}

		else {

		    iflag = 0;

		    bound = bl[iperm[icount - 1] - 1];

		}

		l_c7nft (&xval[iperm[icount - 1] - 1], &xval[iperm[icount] -

			1], nxval, &ider[iperm[icount - 1] - 1], korder,

		    xknot, ncoef, &iflag, a, lda, rhs, irow, &bound,

		    &iperm[icount - 1], wk1, wk4, wk3, wk2);

		if (itype[iperm[icount - 1] - 1] == 1)

		    *neq += 1;

		*irow += 1;

		icount += 2;

		goto L_10;

	    }

	}

	else {

	    /*

	     * Start of constraints on derivatives. IDER .ge. (0) and ITYPE =

	     * 1, 2, or 3

	     */

	    if ((itype[iperm[icount - 1] - 1] >= 1) && (itype[iperm[icount - 1] -

			1] <= 3)) {

		if (itype[iperm[icount - 1] - 1] == 2) {

		    iflag = 1;

		    bound = bu[iperm[icount - 1] - 1];

		}

		else {

		    iflag = 0;

		    bound = bl[iperm[icount - 1] - 1];

		}

		l_c6nft (&xval[iperm[icount - 1] - 1], nxval, korder,

		    &ider[iperm[icount - 1] - 1], xknot, ncoef, &iflag,

		    a, lda, rhs, irow, &bound, wk2, wk3, &iperm[icount - 1]);

		if (itype[iperm[icount - 1] - 1] == 1)

		    *neq += 1;

		*irow += 1;

		icount += 1;

		goto L_10;

		/* IDER .ge. (0) and ITYPE = 4 */

	    }

	    else if (itype[iperm[icount - 1] - 1] == 4) {

		iflag = 1;

		l_c6nft (&xval[iperm[icount - 1] - 1], nxval, korder,

		    &ider[iperm[icount - 1] - 1], xknot, ncoef, &iflag,

		    a, lda, rhs, irow, &bu[iperm[icount - 1] - 1], wk2,

		    wk3, &iperm[icount - 1]);

		*irow += 1;

		iflag = 0;

		l_c6nft (&xval[iperm[icount - 1] - 1], nxval, korder,

		    &ider[iperm[icount - 1] - 1], xknot, ncoef, &iflag,

		    a, lda, rhs, irow, &bl[iperm[icount - 1] - 1], wk2,

		    wk3, &iperm[icount - 1]);

		*irow += 1;

		icount += 1;

		goto L_10;

		/* IDER .ge. (0) and ITYPE = (-4) */

	    }

	    else if (itype[iperm[icount - 1] - 1] == -4) {

		iflag = 0;

		l_c6nft (&xval[iperm[icount - 1] - 1], nxval, korder,

		    &ider[iperm[icount - 1] - 1], xknot, ncoef, &iflag,

		    a, lda, rhs, irow, &bl[iperm[icount - 1] - 1], wk2,

		    wk3, &iperm[icount - 1]);

		iflag = 1;

		_f0 = 0.0;

		l_c6nft (&xval[iperm[icount] - 1], nxval, korder, &ider[iperm[icount] -

			1], xknot, ncoef, &iflag, a, lda, rhs, irow, &_f0,

		    wk2, wk3, &iperm[icount - 1]);

		*irow += 1;

		iflag = 1;

		l_c6nft (&xval[iperm[icount - 1] - 1], nxval, korder,

		    &ider[iperm[icount - 1] - 1], xknot, ncoef, &iflag,

		    a, lda, rhs, irow, &bu[iperm[icount - 1] - 1], wk2,

		    wk3, &iperm[icount - 1]);

		icount += 1;

		iflag = 0;

		_f0 = 0.0;

		l_c6nft (&xval[iperm[icount - 1] - 1], nxval, korder,

		    &ider[iperm[icount - 1] - 1], xknot, ncoef, &iflag,

		    a, lda, rhs, irow, &_f0, wk2, wk3, &iperm[icount - 2]);

		icount += 1;

		*irow += 1;

		goto L_10;

	    }

	    else {

		/*

		 * IDER .ge. (0) and ITYPE = (-1), (-2), or (-3).

		 */

		if (itype[iperm[icount - 1] - 1] == -2) {

		    iflag = 0;

		    bound = bu[iperm[icount - 1] - 1];

		}

		else {

		    iflag = 1;

		    bound = bl[iperm[icount - 1] - 1];

		}

		l_c6nft (&xval[iperm[icount - 1] - 1], nxval, korder,

		    &ider[iperm[icount - 1] - 1], xknot, ncoef, &iflag,

		    a, lda, rhs, irow, &bound, wk2, wk3, &iperm[icount - 2]);

		icount += 1;

		iflag = mod (iflag, 1);

		bound = 0.0;

		l_c6nft (&xval[iperm[icount - 1] - 1], nxval, korder,

		    &ider[iperm[icount - 1] - 1], xknot, ncoef, &iflag,

		    a, lda, rhs, irow, &bound, wk2, wk3, &iperm[icount - 2]);

		if (itype[iperm[icount - 1] - 1] == -1)

		    *neq += 1;

		icount += 1;

		*irow += 1;

		goto L_10;

	    }

	}

    }

    else {

	/*

	 * If all constraints have been implemented exit C5NFT.

	 */

	goto L_9000;

    }



L_9000:

    *irow -= 1;

    imsl_e1pop ("C5NFT");

    return;

}				/* end of function */



#undef A



/*  -----------------------------------------------------------------------

    IMSL Name:  C6NFT/DC6NFT (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    March 1, 1990



    Purpose:    Compute one row of the constraint matrix A for a

                derivative constraint.



    Usage:      CALL C6NFT (X, NXVAL, KORDER, IDER, XKNOT, NCOEF,

                            IFLAG, A, LDA, RHS, IROW, BOUND,

                            WKA, DBIATX, ICONST )



    Arguments:

       X      - Value at which the constraint is made.  (Input)

       NXVAL  - Number of points in the vector XVAL.  (Input)

       KORDER - Order of the spline.  (Input)

       IDER   - Array of length NXVAL containing the derivative

                value of the spline which is to be fit.  (Input)

       XKNOT  - Array of length NCOEF+KORDER containing the knot

                sequence.  (Input)

                XKNOT must be nondecreasing.

       NCOEF  - Number of B-spline coefficients.  (Input)

       IFLAG  - Integer flag indicatingthe type of constraint.  (Input)

                IFLAG should be set to one for .LE. and set to

                zero otherwise.

       A      - The constraint matrix.  (Output)

       LDA    - Leading dimension of A exactly as specified in the

                dimension statement of the calling program.  (Input)

       RHS    - Array of size at least 2*NXVAL+KORDER that is the

                right hand side of the system of constraints in the

                matrix A.  (Output)

       IROW   - The row of A the constraint should be entered.  (Input)

       BOUND  - The bound for the constraint.  (Input)

       WKA    - Work array of size KORDER*KORDER.  (Input)

       DBIATX - Work array of size (KORDER+1)*KORDER.  (Input)

       ICONST - Integer indicating which constraint is being inserted

                in the matrix A.  (Input)

                This value is stored in the (NCOEF+1)st column of the

                matrix A.



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_c6nft (Mfloat *x, Mint *nxval, Mint *korder,

                Mint *ider, Mfloat xknot[], Mint *ncoef, Mint *iflag,

                Mfloat *a, Mint *lda, Mfloat rhs[], Mint *irow, Mfloat *bound,

                Mfloat *wka, Mfloat *dbiatx, Mint *iconst)

#else

static void l_c6nft (x, nxval, korder, ider, xknot, ncoef, iflag,

                a, lda, rhs, irow, bound, wka, dbiatx, iconst)

    Mfloat     *x;

    Mint       *nxval, *korder, *ider;

    Mfloat      xknot[];

    Mint       *ncoef, *iflag;

    Mfloat     *a;

    Mint       *lda;

    Mfloat      rhs[];

    Mint       *irow;

    Mfloat     *bound, *wka, *dbiatx;

    Mint       *iconst;

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

#define WKA(I_,J_)	(wka+(I_)*(*korder)+(J_))

#define DBIATX(I_,J_)	(dbiatx+(I_)*(*korder)+(J_))

    Mint        _l0, i, left, mflag;

    Mfloat      const_;

    /*

     * Find constant to transform inequality to .ge. for Q2ROG

     */

    if (*iflag == 1) {

	const_ = -1.0;

    }

    else {

	const_ = 1.0;

    }

    sset (*korder * (*ider + 1), 0.0e0, dbiatx, 1);

    /* Find the first knot to the left of X */

    _l0 = *ncoef + *korder;

    imsl_b4der (xknot, &_l0, x, &left, &mflag);

    /*

     * Compute the values of the possibly nonzero derivatives.

     */

    _l0 = *ider + 1;

    l_b32gd (xknot, korder, x, &left, wka, dbiatx, &_l0);

    /*

     * Compute the row of the A matrix. loop over all possibly nonzero

     * derivatives

     */

    for (i = left - *korder + 1; i <= left; i++) {

	*A (i - 1, *irow - 1) += const_ ** DBIATX (*ider, i - left + *korder - 1);

    }

    /*

     * Compute the right hand side of the equality/inequality

     */

    rhs[*irow - 1] += *bound * const_;

    /*

     * Insert the consraint number in the last column of the matrix A.

     */

    *A (*ncoef, *irow - 1) = (Mfloat) (*iconst);



    return;

}				/* end of function */



#undef A

#undef WKA

#undef DBIATX



/*  -----------------------------------------------------------------------

    IMSL Name:  C7NFT/DC7NFT (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    March 1, 1990



    Purpose:    Compute a one row of the constraint matrix A for an

                integral constraint.



    Usage:      CALL C7NFT (S, T, NXVAL, IDER, KORDER, XKNOT, NCOEF,

                            IFLAG, A, LDA, RHS, IROW, BOUND,

                            ICONST, WK1, WK2, WK3, WK4 )



    Arguments:

       S      - Lower limit of the constrained integral.  (Input)

       T      - Upper limit of the constrained integral.  (Input)

       NXVAL  - Number of points in the vector XVAL.  (Input)

       IDER   - Array of length NXVAL containing the derivative

                value of the spline which is to be fit.  (Input)

       KORDER - Order of the spline.  (Input)

       XKNOT  - Array of length NCOEF+KORDER containing the knot

                sequence.  (Input)

                XKNOT must be nondecreasing.

       NCOEF  - Number of B-spline coefficients.  (Input)

       IFLAG  - Integer flag indicatingthe type of constraint.  (Input)

                IFLAG should be set to one for .LE. and set to

                zero otherwise.

       A      - The constraint matrix.  (Output)

       LDA    - Leading dimension of A exactly as specified in the

                dimension statement of the calling program.  (Input)

       RHS    - Array of size at least 2*NXVAL+KORDER that is the

                right hand side of the system of constraints in the

                matrix A.  (Output)

       IROW   - The row of A the constraint should be entered.  (Input)

       BOUND  - The bound for the constraint.  (Input)

       ICONST - Integer indicating which constraint is being inserted

                in the matrix A.  (Input)

                This value is stored in the (NCOEF+1)st column of the

                matrix A.

       WK1    - Work array of size KORDER+1.  (Input)

       WK2    - Work array of size KORDER+1.  (Input)

       WK3    - Work array of size KORDER+1.  (Input)

       WK4    - Work array of size KORDER+NCOEF+1.  (Input)



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_c7nft (Mfloat *s, Mfloat *t, Mint *nxval, Mint *ider,

                Mint *korder, Mfloat xknot[], Mint *ncoef, Mint *iflag,

                Mfloat *a, Mint *lda, Mfloat rhs[], Mint *irow, Mfloat *bound,

                Mint *iconst, Mfloat wk1[], Mfloat wk2[], Mfloat wk3[],

                Mfloat wk4[])

#else

static void l_c7nft (s, t, nxval, ider, korder, xknot, ncoef, iflag,

                a, lda, rhs, irow, bound, iconst, wk1, wk2, wk3, wk4)

    Mfloat     *s, *t;

    Mint       *nxval, *ider, *korder;

    Mfloat      xknot[];

    Mint       *ncoef, *iflag;

    Mfloat     *a;

    Mint       *lda;

    Mfloat      rhs[];

    Mint       *irow;

    Mfloat     *bound;

    Mint       *iconst;

    Mfloat      wk1[], wk2[], wk3[], wk4[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

    Mint        _l0, i, lefts, leftt, mflag;

    Mfloat      const_;

    /*

     * Find constant to transform inequality to .ge. for imsl_q2rog.

     */

    if (*iflag == 1) {

	const_ = -1.0;

    }

    else {

	const_ = 1.0;

    }

    /* Find LEFTS */

    _l0 = *ncoef + *korder;

    imsl_b4der (xknot, &_l0, s, &lefts, &mflag);

    /* Find LEFTT */

    _l0 = *ncoef + *korder;

    imsl_b4der (xknot, &_l0, t, &leftt, &mflag);

    sset (*ncoef, 0.0, &wk4[*korder + 1], 1);

    /*

     * Loop over all possibly nonzero B-spline values

     */

    for (i = lefts - *korder + 1; i <= leftt; i++) {

	wk4[*korder + i] = 1.0;

	*A (i - 1, *irow - 1) += const_ * imsl_b2itg (s, t, korder, xknot,

	    ncoef, &wk4[*korder + 1], wk1, wk2, wk3, wk4);

	wk4[*korder + i] = 0.0;

    }

    /* Compute right hand side */

    rhs[*irow - 1] += *bound * const_;

    /*

     * Insert the constraint number in the last column of the matrix A.

     */

    *A (*ncoef, *irow - 1) = (Mfloat) (*iconst);

    return;

}				/* end of function */



#undef A



/*  -----------------------------------------------------------------------

    IMSL Name:  C8NFT/DC8NFT (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    March 1, 1990



    Purpose:    Solve the quadratic programming problem in CONFT.



    Usage:      CALL C8NFT (NXVAL, NHARD, IDER, ITYPE, KORDER,

                            NCOEF, BSCOEF, H, G, A, LDA, RHS, IPERM,

                            NEQ, IROW, WK1, WK2, WK3, IWK1, IWK2 )



    Arguments:

       NXVAL  - Number of points in the vector XVAL.  (Input)

       NHARD  - The number of points involved in the 'hard' constraints.

                (Input)

       IDER   - Array of length NXVAL containing the derivative

                value of the spline which is to be fit.  (Input)

       ITYPE  - Array of length NXVAL indicating the types of general

                constraints in the matrix A.  (Input)

       KORDER - Order of the spline.  (Input)

       NCOEF  - Number of B-spline coefficients.  (Input)

       BSCOEF - Array of length NCOEF containing the B-spline

                coefficients.  (Output)

       H      - Array of size NCOEF by NCOEF containing the Hessian

                matrix of the objective function to be used in a

                call to QPROG.  (Input)

       G      - Array of length NCOEF containing the coefficients of

                the linear term of the objective function used in a

                call to QPROG.  (Input)

       A      - The constraint matrix.  (Input)

                The index of the constraints are stored in the

                (NCOEF+1)st column of A.

       LDA    - Leading dimension of A exactly as specified in the

                dimension statement of the calling program.  (Input)

       RHS    - Array of size at least 2*NXVAL+KORDER that is the

                right hand side of the system of constraints in the

                matrix A.  (Input)

       IPERM  - Array of length NXVAL that holds the permutation of

                the constraints.  (Input)

       NEQ    - The number of equality constraints in A.  (Input)

                The equality constraints are in the first NEQ rows the

                matrix A.

       IROW   - The number of rows of constraints in the matrix A.

                (Input)

       WK1    - Work array of size NCOEF.  (Input)

       WK2    - Work array of size (3*NCOEF*NCOEF+11*NCOEF)/2 +

                2*NXVAL+KORDER.  (Input)

       WK3    - Work array of size 4*NCOEF+IROW*(IROW+30)+1.  (Input)

       IWK1   - Integer work array of length NCOEF.  (Input)

       IWK2   - Integer work array of length 30*IROW+3*NCOEF.  (Input)



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_c8nft (Mint *nxval, Mint *nhard, Mint ider[],

                Mint itype[], Mint *korder, Mint *ncoef,

                Mfloat bscoef[], Mfloat *h, Mfloat g[], Mfloat *a, Mint *lda,

                Mfloat rhs[], Mint iperm[], Mint *neq, Mint *irow, Mfloat wk1[],

                Mfloat wk2[], Mfloat wk3[], Mint iwk1[], Mint iwk2[])

#else

static void l_c8nft (nxval, nhard, ider, itype, korder, ncoef,

                bscoef, h, g, a, lda, rhs, iperm, neq, irow, wk1, wk2, wk3, iwk1,

                iwk2)

    Mint       *nxval, *nhard, ider[], itype[], *korder, *ncoef;

    Mfloat      bscoef[], *h, g[], *a;

    Mint       *lda;

    Mfloat      rhs[];

    Mint        iperm[], *neq, *irow;

    Mfloat      wk1[], wk2[], wk3[];

    Mint        iwk1[], iwk2[];

#endif

{

#define H(I_,J_)	(h+(I_)*(*ncoef)+(J_))

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

    Mint        iflag, k, kfirst, ki1, ki2, kr1, kr2, kr3, krow, ldh, nact;

    Mfloat      diag, ermax, xmtch1, xmtch2;





    imsl_e1psh ("l_c8nft");

    ermax = imsl_amach (4) * 1.0e02;

    /*

     * Put data in correct form to send to C10FT. Compute workspace offsets.

     */

    ki1 = 1;

    ki2 = ki1 + *irow;

    kr1 = 1;

    kr2 = kr1 + *ncoef;

    kr3 = kr2 + *ncoef;

    /*

     * Copy the entries of RHS to the workspace that is going to be used as

     * BL in the call to C10FT. Set the workspace that is sent to C10FT for

     * ITYPE to have NEQ equality constraints and the rest greater than or

     * equal.

     */

    iset (*irow, 2, &iwk2[ki1 - 1], 1);

    iset (*neq, 0, &iwk2[ki1 - 1], 1);

    /*

     * Set the workspace sent to R10FT for XLB and XUB to values that will

     * not put bounds on the variables.

     */

    sset (*ncoef, 1.0e30, &wk3[kr1 - 1], 1);

    sset (*ncoef, -1.0e30, &wk3[kr2 - 1], 1);

    k = *nxval;

L_10:

    if (k >= *nhard) {

	/* See if the problem is feasible. */

	l_c10ft (irow, ncoef, a, lda, rhs, rhs, &iwk2[ki1 - 1], &wk3[kr1 - 1],

	    &wk3[kr2 - 1], &iflag, &wk3[kr3 - 1], &iwk2[ki2 - 1]);

	/*

	 * If IFLAG is not equal to one, then remove a constraint from A by

	 * zeroing out the row(s) of A, entries of RHS, and entries of the

	 * workspace corresponding to the BL int he call to C10FT which

	 * correspond to the lowest priority constraint left.

	 */

L_20:

	if (iflag != 1) {

	    if (k <= *nhard)

		goto L_50;

	    /*

	     * Compute which values to look for in the last column of A.

	     */

	    if ((ider[k - 1] == -1) || (itype[k - 1] < 0)) {

		xmtch1 = (Mfloat) (k);

		xmtch2 = (Mfloat) (k - 1);

		k -= 1;

	    }

	    else {

		xmtch1 = (Mfloat) (k);

		xmtch2 = (Mfloat) (k);

	    }

	    kfirst = 1;

    L_30:

	    if ((*A (*ncoef, kfirst - 1) != xmtch1) && (*A (*ncoef, kfirst - 1) !=

		    xmtch2)) {

		kfirst += 1;

		goto L_30;

	    }

	    /*

	     * KFIRST is the first row of a that contains at least a part of

	     * the constraint. We now zero out all the rows that make up the

	     * constaint.

	     */

	    krow = kfirst;

	    /* Remove the constraint. */

    L_40:

	    sset (*ncoef, 0.0, A (0, krow - 1), *lda);

	    rhs[krow - 1] = 0.0;

	    if (krow < *irow) {

		if ((fabs (*A (*ncoef, krow) - xmtch1) < ermax) ||

		    (fabs (*A (*ncoef, krow) - xmtch2) < ermax)) {

		    krow += 1;

		    goto L_40;

		}

	    }

	    k -= 1;

	    goto L_10;

	}

	else {

	    /* Try to solve the QP problem. */

	    ldh = *ncoef;

	    imsl_q2rog (*ncoef, *irow, *neq, a, *lda, rhs, g, h, ldh, &diag,

		bscoef, &nact, iwk1, wk1, wk2);

	    /*

	     * Tell how many constraints were removed.

	     */

	    if ((*nxval - k) > 0) {

		imsl_e1sti (1, *nxval - k);

                imsl_ermes ( IMSL_WARNING, IMSL_CONSTR_REMOVED);

	    }

	    goto L_9000;

	}

    }

    /*

     * Exit with error saying the hard constraints could not be met.

     */

L_50:

    imsl_e1sti (1, *nhard);

    imsl_e1sti (2, *nxval - k);

    imsl_ermes(IMSL_FATAL, IMSL_NO_FIT_OBTAINED);

L_9000:

    imsl_e1pop ("l_c8nft");

    return;

}				/* end of function */



#undef H

#undef A



/*  -----------------------------------------------------------------------

    IMSL Name:  C9NFT/DC9NFT (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    March 1, 1990



    Purpose:    Check B-spline parameters.



    Usage:      CALL C9NFT (NDATA, XDATA, KORDER,

                            XKNOT, WEIGHT, NCOEF, SCALAR)



    Arguments:

       NDATA  - Number of data points.  (Input)

       XDATA  - Array of length NDATA containing the data

                point abscissas.  (Input)

       KORDER - Order of the spline.  (Input)

       XKNOT  - Array of length NCOEF+KORDER containing the knot

                sequence.  (Input)

                XKNOT must be nondecreasing.

       WEIGHT - Array of length NDATA containing the weights.  (Input)

       NCOEF  - Length of BSCOEF.  (Input)

       SCALAR - Logical variable telling if only scalars are to be

                checked.  (Input)

                SCALAR should be set to .TRUE. if only scalars

                are to be checked.



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_c9nft (Mint *ndata, Mfloat xdata[], Mint *korder,

                Mfloat xknot[], Mfloat weight[], Mint *ncoef,

                Mint *scalar)

#else

static void l_c9nft (ndata, xdata, korder, xknot, weight, ncoef,

                scalar)

    Mint       *ndata;

    Mfloat      xdata[];

    Mint       *korder;

    Mfloat      xknot[], weight[];

    Mint       *ncoef;

    Mint       *scalar;

#endif

{

    Mint        i, mult;

    /* Check KORDER */

    if (*korder < 1) {

	imsl_e1sti (1, *korder);

	imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);

	goto L_9000;

    }

    /* Check NCOEF */

    if (*ncoef < *korder) {

	imsl_e1sti (1, *ncoef);

	imsl_e1sti (2, *korder);

	imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_COEFF_X);

	goto L_9000;

    }

    /* Check NCOEF */

    if (*ncoef > *ndata) {

	imsl_e1sti (1, *ncoef);

	imsl_e1sti (2, *ndata);

        imsl_ermes (IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS_2);

	goto L_9000;

    }

    /* Check knot sequence */

    if (!*scalar) {

	mult = 1;

	for (i = 2; i <= (*ncoef + *korder); i++) {

	    if (xknot[i - 1] == xknot[i - 2]) {

		mult += 1;

		if (mult > *korder) {

                    imsl_e1stl (1,"X");

		    imsl_e1sti (1, i - mult);

		    imsl_e1sti (2, i - 1);

		    imsl_e1str (1, xknot[i - 1]);

		    imsl_e1sti (3, *korder);

	            imsl_ermes (IMSL_FATAL, IMSL_KNOT_MULTIPLICITY);

		    goto L_9000;

		}

	    }

	    else if (xknot[i - 1] < xknot[i - 2]) {

                imsl_e1stl (1,"X");

		imsl_e1sti (1, i - 2);

		imsl_e1sti (2, i - 1);

		imsl_e1str (1, xknot[i - 2]);

		imsl_e1str (2, xknot[i - 1]);

	         imsl_ermes (IMSL_FATAL, IMSL_KNOT_NOT_INCREASING);

		goto L_9000;

	    }

	    else {

		mult = 1;

	    }

	}

	/* Test XDATA(1) .GE. XKNOT(KORDER) */

	if (xdata[0] < xknot[*korder - 1]) {

            imsl_e1stl (1,"X");

	    imsl_e1str (1, xdata[0]);

	    imsl_e1str (2, xknot[*korder - 1]);

	    imsl_ermes (IMSL_FATAL, IMSL_XDATA_TOO_SMALL);

	    goto L_9000;

	}

	/* Test XDATA(NDATA) .LE. XKNOT(NDATA+1) */

	if (xdata[*ndata - 1] > xknot[*ncoef]) {

            imsl_e1stl (1,"X");

	    imsl_e1str (1, xdata[*ndata - 1]);

	    imsl_e1str (2, xknot[*ncoef]);

	    imsl_ermes (IMSL_FATAL, IMSL_XDATA_TOO_LARGE);

	    goto L_9000;

	}

	/* Check weights */

	for (i = 1; i <= *ndata; i++) {

	    if (weight[i - 1] <= 0.0e0) {

		imsl_e1sti (1, i - 1);

		imsl_e1str (1, weight[i - 1]);

                imsl_ermes ( IMSL_FATAL, IMSL_WEIGHT_LE_ZERO); 

		goto L_9000;

	    }

	}

    }



L_9000:

    ;

    return;

}				/* end of function */

/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 16:08:50

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  PERMU/DPERMU (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    July 30, 1985



    Purpose:    Rearrange the elements of an array as specified by a

                permutation.



    Usage:      CALL PERMU (N, X, IPERMU, IPATH, XPERMU)



    Arguments:

       N      - Length of the arrays X and XPERMU.  (Input)

       X      - Real vector of length N containing the array to be

                permuted.  (Input)

       IPERMU - Integer vector of length N containing a permutation

                IPERMU(1), ..., IPERMU(N) of the integers 1, ..., N.

                (Input)

       IPATH  - Integer flag.  (Input)

                IPATH = 1 means IPERMU represents a forward permutation,

                          i.e., X(IPERMU(I)) is moved to XPERMU(I).

                IPATH = 2 means IPERMU represents a backward permutation,

                          i.e., X(I) is moved to XPERMU(IPERMU(I)).

       XPERMU - Real vector of length N containing the array X permuted.

                (Output)

                If X is not needed, X and XPERMU can share the same

                storage locations.



    Keywords:   Utilities; Forward permutation; Backward permutation



    GAMS:       N8



    Chapters:   MATH/LIBRARY Utilities

                STAT/LIBRARY Utilities



    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_permu (Mint *n, Mfloat *x, Mint *ipermu, Mint *ipath, Mfloat *xpermu)

#else

static void l_permu (n, x, ipermu, ipath, xpermu)

    Mint       *n;

    Mfloat      x[];

    Mint        ipermu[], *ipath;

    Mfloat      xpermu[];

#endif

{

    Mint        i, j, k;

    Mfloat      temp;





    imsl_e1psh ("l_permu");



    if (*n <= 0) {

	imsl_e1sti (1, *n);



	/*

	 * (5, 1, "The length of the arrays X and XPERMU must be positive

	 * while N = %(i1) is given.");

	 */

	imsl_ermes (IMSL_TERMINAL, IMSL_X_AND_XPERMU_LENGTH_LE_0);

    }

    if (*ipath != 1 && *ipath != 2) {

	imsl_e1sti (1, *ipath);



	/*

	 * (5, 2, "IPATH must be equal to 1 or 2 while a value of %(i1) is

	 * given.");

	 */

	imsl_ermes (IMSL_TERMINAL, IMSL_IPATH_RANGE_3);

    }

    if (imsl_n1rcd (0) != 0)

	goto L_9000;



    /*

     * MAKE A COPY OF X IN XPERMU AND WORK WITH XPERMU

     */

    scopy (*n, x, 1, xpermu, 1);

    if (*n == 1)

	goto L_9000;

    for (i = 1; i <= *n; i++) {

	if (ipermu[i - 1] < 1 || ipermu[i - 1] > *n) {

	    imsl_e1sti (1, i);

	    imsl_e1sti (2, *n);

	    imsl_e1sti (3, ipermu[i - 1]);



	    /*

	     * (5, 3, "IPERMU(%(i1)) = %(i3) is not allowed.  It must be

	     * between 1 and N = %(i2).");

	     */

	    imsl_ermes (IMSL_TERMINAL, IMSL_IPERMU_RANGE);

	}

	else {

	    ipermu[i - 1] = -ipermu[i - 1];

	}

    }

    if (imsl_n1rcd (0) != 0)

	goto L_9000;



    if (*ipath == 1) {



	/*

	 * FORWARD PERMUTATION

	 */

	for (i = 1; i <= *n; i++) {

	    if (ipermu[i - 1] <= 0) {

		j = i;

		ipermu[j - 1] = -ipermu[j - 1];

		k = ipermu[j - 1];

	L_20:

		;

		if (ipermu[k - 1] <= 0) {

		    temp = xpermu[j - 1];

		    xpermu[j - 1] = xpermu[k - 1];

		    xpermu[k - 1] = temp;

		    ipermu[k - 1] = -ipermu[k - 1];

		    j = k;

		    k = ipermu[k - 1];

		    goto L_20;

		}

	    }

	}

    }

    else {



	/*

	 * BACKWARD PERMUTATION

	 */

	for (i = 1; i <= *n; i++) {

	    if (ipermu[i - 1] <= 0) {

		ipermu[i - 1] = -ipermu[i - 1];

		j = ipermu[i - 1];

	L_40:

		;

		if (j != i) {

		    temp = xpermu[i - 1];

		    xpermu[i - 1] = xpermu[j - 1];

		    xpermu[j - 1] = temp;

		    ipermu[j - 1] = -ipermu[j - 1];

		    j = ipermu[j - 1];

		    goto L_40;

		}

	    }

	}

    }

L_9000:

    imsl_e1pop ("l_permu");

    return;

}

/*  -----------------------------------------------------------------------

    IMSL Name:  B32GD/DB32GD (Single/Double precision version)



    Computer:   SUN/SINGLE



    Revised:    January 31, 1990



    Purpose:    Calculate the value and derivatives of all B-splines

                which do not vanish at X.



    Usage:      CALL B32GD (T, K, X, LEFT, A, DBIATX, NDERIV)



    Arguments:

       T      - Array of length LEFT+K containing the knot sequence.

                (Input)

       K      - Order of the spline.  (Input)

       X      - Point at which derivative is to be evaluated.  (Input)

       LEFT   - The left endpoint of the interval of interest.  (Input)

                The K B-splines whose support contains the interval

                (T(LEFT,T(LEFT+1)) are to be considered.

       A      - Work array of order (K,K).  (Input)

       DBIATX - Array of order (K,NDERIV).  (Output)

                Entry (I,M) contains the value of (M-1)st derivative

                of (LEFT-K+I)th B-spline of order K for knot sequence

                T, I = 1, ... , NDERIV.

       NDERIV - Integer indicating that values of B-splines and their

                derivatives up to but not including the NDERIV-th are

                requested.  (Input)



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_b32gd (Mfloat t[], Mint *k, Mfloat *x, Mint *left,

                                Mfloat *a, Mfloat *dbiatx, Mint *nderiv)

#else

static void l_b32gd (t, k, x, left, a, dbiatx, nderiv)

    Mfloat       t[];

    Mint        *k;

    Mfloat      *x;

    Mint        *left;

    Mfloat      *a, *dbiatx;

    Mint        *nderiv;

#endif

{

#define A(I_,J_)	(a+(I_)*(*k)+(J_))

#define DBIATX(I_,J_)	(dbiatx+(I_)*(*k)+(J_))

    Mint         _l0, _l1, i, ideriv, il, j, jlow, jp1mid, kp1, kp1mm, ldummy, m,

                mhigh;

#ifdef COMPUTER_MIPNT

    Mint	star_k = *k;

    Mint	jm1, im1, im2;

#endif

    Mfloat       factor, fkp1mm, sum;





    mhigh = imsl_i_max (imsl_i_min (*nderiv, *k), 1);



    kp1 = *k + 1;

    _l0 = kp1-mhigh;

    _l1 = 1;

    l_b42gd (t, &_l0, &_l1, x, left, dbiatx);

    if (mhigh == 1)

	goto L_90;



    ideriv = mhigh;

    for (m = 2; m <= mhigh; m++) {

	jp1mid = 1;

	for (j = ideriv; j <= *k; j++) {

	    *DBIATX (ideriv - 1, j - 1) = *DBIATX (0, jp1mid - 1);

	    jp1mid += 1;

	}

	ideriv -= 1;

    _l0 = kp1-ideriv;

    _l1 = 2;

	l_b42gd (t, &_l0, &_l1, x, left, dbiatx);

    }



    jlow = 1;

    for (i = 1; i <= *k; i++) {

	for (j = jlow; j <= *k; j++) {

	    *A (i - 1, j - 1) = 0.0;

	}

	jlow = i;

	*A (i - 1, i - 1) = 1.0;

    }



    for (m = 2; m <= mhigh; m++) {

	kp1mm = kp1 - m;

	fkp1mm = (float) (kp1mm);

	il = *left;

	i = *k;

	for (ldummy = 1; ldummy <= kp1mm; ldummy++) {

	    factor = fkp1mm / (t[il + kp1mm - 1] - t[il - 1]);

	    for (j = 1; j <= i; j++) {

#ifdef COMPUTER_MIPNT

		jm1 = j - 1;

		im1 = i - 1;

		im2 = i - 2;

		*(a+jm1*star_k+im1) = (*(a+jm1*star_k+im1) - *(a+jm1*star_k+im2)) * factor;

#else

		*A (j - 1, i - 1) = (*A (j - 1, i - 1) - *A (j - 1, i - 2)) *

		    factor;

#endif

	    }

	    il -= 1;

	    i -= 1;

	}



	for (i = 1; i <= *k; i++) {

	    sum = 0.0;

	    jlow = imsl_i_max (i, m);

	    for (j = jlow; j <= *k; j++) {

		sum += *A (i - 1, j - 1) ** DBIATX (m - 1, j - 1);

	    }

	    *DBIATX (m - 1, i - 1) = sum;

	}

    }

L_90:

    return;

}				/* end of function */



#undef A

#undef DBIATX



/*----------------------------------------------------------------------- */



/*  IMSL Name:  B42GD/DB42GD (Single/Double precision version)



    Computer:   SUN/SINGLE



    Revised:    January 31, 1990



    Purpose:    Calculate the value of all possibly nonzero B-splines

                at X of order JOUT = MAX(JHIGH,(J+1*(INDEX-1)) with

                knot sequence T.



    Usage:      CALL B42GD (T, JHIGH, INDEX, X, LEFT, BIATX)



    Arguments:

       T      - Array of length JOUT+LEFT containing the knot sequence.

                (Input)

       JHIGH  - Integer used to determine the order.  (Input)

       INDEX  - Integer used to determine the order.  (Input)

       X      - Point at which the B-splines are to be evaluated.

                (Input)

       LEFT   - Integer chosen so that T(LEFT).LE.X.LE.T(LEFT+1).

                (Input)

       BIATX  - Array of length JOUT, with BIATX(I) containing the

                value at X of the polynomial of order JOUT which agrees

                with the B-spline B(LEFT-JOUT+I,JOUT,T) on the interval

                (T(LEFT),T(LEFT+1)).  (Output)



    Chapter:    MATH/LIBRARY Interpolation and Approximation



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#define	JMAX	20

#ifdef ANSI

static void l_b42gd (Mfloat t[], Mint *jhigh, Mint *index, Mfloat *x,

                   Mint * left, Mfloat biatx[])

#else

static void l_b42gd (t, jhigh, index, x, left, biatx)

    Mfloat       t[];

    Mint        *jhigh, *index;

    Mfloat      *x;

    Mint        *left;

    Mfloat       biatx[];

#endif

{

    Mint         i, jp1;

    static Mfloat       deltal[JMAX], deltar[JMAX], saved, term;

    static Mint  j = 1;







    if (*index == 1) {

	j = 1;

	biatx[0] = 1.0;

	if (j >= *jhigh)

	    goto L_30;

    }

L_10:

    jp1 = j + 1;

    deltar[j - 1] = t[*left + j - 1] - *x;

    deltal[j - 1] = *x - t[*left - j];

    saved = 0.0;

    for (i = 1; i <= j; i++) {

	term = biatx[i - 1] / (deltar[i - 1] + deltal[jp1 - i - 1]);

	biatx[i - 1] = saved + deltar[i - 1] * term;

	saved = deltal[jp1 - i - 1] * term;

    }

    biatx[jp1 - 1] = saved;

    j = jp1;

    if (j < *jhigh)

	goto L_10;

L_30:

    return;

}				/* end of function */

#undef JMAX



