#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif



#define max(a,b) ( ((a) > (b)) ? (a) : (b))

#ifdef TRUE

#undef TRUE

#define TRUE 1

#else

#define TRUE 1

#endif



static void     PROTO (l_permu, (Mint, Mfloat*, Mint*, Mint, Mfloat*));

static void     PROTO (l_qmdmrg, (Mint[], Mint[], Mint[],

         Mint[], Mint[], Mint[], Mint*,

         Mint*, Mint[], Mint[], Mint[]));

static void     PROTO (l_qmdqt, (Mint*, Mint[], Mint[],

         Mint[], Mint*, Mint[], Mint[]));

static void     PROTO (l_qmdrch, (Mint*, Mint[], Mint[],

         Mint[], Mint[], Mint*, Mint[],

         Mint*, Mint[]));

static void     PROTO (l_qmdupd, (Mint[], Mint[], Mint*,

         Mint[], Mint[], Mint[], Mint[],

         Mint[], Mint[], Mint[]));

static void     PROTO (l_standard_to_adjacency, (Mint, Mf_sparse_elem*,

         Mint, Mint*, Mint*, Mint*));

static void     PROTO (l_CSC_to_adjacency, (Mint, Mint*, Mint*, Mint*, Mint*,

         Mint*));

static void     PROTO (l_smbfct, (Mint*, Mint[], Mint[],

         Mint[], Mint[], Mint[], Mint*,

         Mint[], Mint[], Mint*, Mint*));

static void     PROTO (l_genqmd, (Mint*, Mint[], Mint[], Mint[], Mint[], Mint*));

static void     PROTO (l_sym_fact_with_ordering, (Mint, Mint, Mf_sparse_elem*,

         Mint*, Mint**, Mint*, Mint*, Mint*, Mint*, Mint*, Mint, Mint, Mint*, Mint*));

static void     PROTO (l_coordinate_to_compressed, (Mint, Mf_sparse_elem*,

         Mint*, Mint*, Mint*, Mint*, Mfloat*, Mfloat*, Mint*));

static void     PROTO (l_CSC_to_compressed, (Mint, Mint*, Mint*, Mfloat*,

         Mint*, Mint*, Mint*, Mint*, Mfloat*, Mfloat*, Mint*));

static void     PROTO (l_gsfct, (Mint*, Mint[], Mfloat[],

         Mint[], Mint[], Mfloat[], Mint*));

static void     PROTO (l_gsslv, (Mint*, Mint[], Mfloat[],

         Mint[], Mint[], Mint[], Mfloat[], Mfloat[]));

static void     PROTO (l_post_order, (Mint[], Mint[],

         Mint[], Mint[], Mint*));

static void     PROTO (l_match_vect, (Mint*, Mint[], Mint[],

         Mint[]));

static void     PROTO (l_compute_frontal_storage, (Mint*, Mint[], Mint[],

         Mint[], Mint*));

static void     PROTO (l_mffct, (Mint*, Mint[], Mfloat[],

         Mint[], Mint[], Mfloat[],

         Mint*, Mint*));

static void     PROTO (l_l6fxd, (Mint*, Mfloat[], Mint*,

         Mfloat[], Mint[]));

static void     PROTO (l_l8fxd, (Mint*, Mfloat[], Mint*,

         Mint*));

static void     PROTO (l_l9fxd, (Mint*, Mfloat[], Mint*,

         Mint[]));

static VA_LIST_HACK  PROTO (l_lin_sol_posdef_coordinate, (Mint, Mint,

         Mf_sparse_elem*, Mfloat*, va_list argptr));



static Mfloat *lv_x;

static Mint  multifrontal_space;



#ifdef ANSI

Mfloat *imsl_f_lin_sol_posdef_coordinate (Mint n, Mint nz, Mf_sparse_elem *a,

	Mfloat *rhs, ...)

#else

Mfloat     *imsl_f_lin_sol_posdef_coordinate (n, nz, a, rhs, va_alist)

    Mint         n;

    Mint         nz;

    Mf_sparse_elem *a;

    Mfloat     *rhs;

va_dcl

#endif

{

    va_list     argptr;

    VA_START (argptr, rhs);

    lv_x = NULL;



    E1PSH ("imsl_f_lin_sol_posdef_coordinate", "imsl_d_lin_sol_posdef_coordinate");
#ifndef LINUX64
    argptr = l_lin_sol_posdef_coordinate (n, nz, a, rhs, argptr);
#else
    l_lin_sol_posdef_coordinate (n, nz, a, rhs, argptr);
#endif

    E1POP ("imsl_f_lin_sol_posdef_coordinate", "imsl_d_lin_sol_posdef_coordinate");



    va_end (argptr);

    return lv_x;

}



#ifdef ANSI

static VA_LIST_HACK l_lin_sol_posdef_coordinate (Mint n, Mint nz, Mf_sparse_elem *a,

	Mfloat *rhs, va_list argptr)

#else

static VA_LIST_HACK l_lin_sol_posdef_coordinate (n, nz, a, rhs, argptr)

    Mint         n;

    Mint         nz;

    Mf_sparse_elem *a;

    Mfloat     *rhs;

    va_list     argptr;

#endif

{

    Mint         maxsub;



    Mint        *tmp_vec = NULL;

    static Mint *perm = NULL;

    static Mint *invp = NULL;



    static Mint *xlnz = NULL;

    Mint         maxlnz;

    static Mint *xnzsub = NULL;

    static Mint *nzsub = NULL;



    static Mfloat *diag = NULL;

    static Mfloat *alnz = NULL;

    Mint         flag = 0;



    Mint         code = 1;

    Mint         arg_number = 6;



    Imsl_symbolic_factor *p_symbolic_factor_struct;

    Mf_numeric_factor *p_numeric_factor_struct;



    Mint         return_symbolic_factor = 0;

    Mint         supply_symbolic_factor = 0;

    Mint         symbolic_factor_only = 0;

    Mint         compute_symbolic_factor = 1;



    Mint         return_numeric_factor = 0;

    Mint         supply_numeric_factor = 0;

    Mint         numeric_factor_only = 0;

    Mint         compute_numeric_factor = 1;



    Mint         return_smallest_diag = 0;

    Mint         return_largest_diag = 0;

    Mint        return_num_factor_nonzeros = 0;

    Mfloat      smallest_diag;

    Mfloat      largest_diag;

    Mint       *p_num_factor_nonzeros;

    Mfloat     *p_smallest_diag;

    Mfloat     *p_largest_diag;



    Mint         multifrontal = 0;

    Mint         error = 0;

    Mint         i;

    Mint         ir;

    Mint         jc;



    Mint        CSC_format = 0;

    Mint       *CSC_col_ptr = NULL;

    Mint       *CSC_row_ind = NULL;

    Mfloat     *CSC_values = NULL;

 

    Mint		data_incremented = 0;



    while (code > 0) {

	code = va_arg (argptr, int);

	arg_number++;

	switch (code) {

	case IMSL_RETURN_SYMBOLIC_FACTOR:

	    return_symbolic_factor = 1;

	    p_symbolic_factor_struct = va_arg (argptr,

		Imsl_symbolic_factor *);

	    arg_number++;

	    break;



	case IMSL_SUPPLY_SYMBOLIC_FACTOR:

	    supply_symbolic_factor = 1;

	    compute_symbolic_factor = 0;

	    p_symbolic_factor_struct = va_arg (argptr,

		Imsl_symbolic_factor *);

	    arg_number++;

	    break;



	case IMSL_SYMBOLIC_FACTOR_ONLY:

	    symbolic_factor_only = 1;

	    compute_numeric_factor = 0;

	    break;



	case IMSL_RETURN_NUMERIC_FACTOR:

	    return_numeric_factor = 1;

	    p_numeric_factor_struct = va_arg (argptr,

		Mf_numeric_factor *);

	    arg_number++;

	    break;



	case IMSL_SUPPLY_NUMERIC_FACTOR:

	    supply_numeric_factor = 1;

	    compute_numeric_factor = 0;

	    p_numeric_factor_struct = va_arg (argptr,

		Mf_numeric_factor *);

	    arg_number++;

	    break;



	case IMSL_NUMERIC_FACTOR_ONLY:

	    numeric_factor_only = 1;

	    break;



	case IMSL_SOLVE_ONLY:

	    compute_symbolic_factor = 0;

	    compute_numeric_factor = 0;



	case IMSL_MULTIFRONTAL_FACTORIZATION:

	    multifrontal = 1;

	    break;



	case IMSL_RETURN_USER:

	    lv_x = va_arg (argptr, Mfloat *);

	    arg_number++;

	    break;



	case IMSL_SMALLEST_DIAGONAL_ELEMENT:

	    p_smallest_diag = va_arg (argptr, Mfloat *);

	    return_smallest_diag = 1;

	    arg_number++;

	    break;



	case IMSL_LARGEST_DIAGONAL_ELEMENT:

	    p_largest_diag = va_arg (argptr, Mfloat *);

	    return_largest_diag = 1;

	    arg_number++;

	    break;



	case IMSL_NUM_NONZEROS_IN_FACTOR:

	    p_num_factor_nonzeros = va_arg (argptr, Mint *);

	    return_num_factor_nonzeros = 1;

	    arg_number++;

	    break;



        case IMSL_CSC_FORMAT:

            CSC_format = 1;

            CSC_col_ptr = va_arg (argptr, Mint *);

            CSC_row_ind = va_arg (argptr, Mint *);

            CSC_values = va_arg (argptr, Mfloat *);

            arg_number += 3;

            break;



	case 0:

	    break;



	default:

	    error++;

	    imsl_e1sti (1, code);

	    imsl_e1sti (2, arg_number);

	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);

	    break;

	}

    }



    if (!(supply_symbolic_factor || supply_numeric_factor)) {

    if (n < 1) {

	error++;

	imsl_e1sti (1, n);

	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);

    }



    if (nz < n) {

	error++;

	imsl_e1sti (1, n);

	imsl_e1sti (2, nz);

	imsl_ermes (IMSL_TERMINAL, IMSL_NZ_LESS_THAN_N);

    }



    if (error)

	goto RETURN;



    if (!CSC_format) {



    for (i = 0; i < nz; i++) {

	a[i].row++;

	a[i].col++;

    }

    data_incremented = 1;



/*  check for entries in the upper triangle  */



    for (i = 0; i < nz; i++)

	if (a[i].row < a[i].col) {

	    error++;

	    imsl_e1sti (1, i);

	    imsl_e1sti (2, a[i].row);

	    imsl_e1sti (3, a[i].col);

	    imsl_ermes (IMSL_TERMINAL, IMSL_NONZERO_IN_UP_TRIANGLE);

	}



/*  check for duplicate entries in A  */



    tmp_vec = (Mint *) malloc (nz * sizeof (*tmp_vec));

    if (!tmp_vec) {

	error++;

	imsl_e1stl (1, "nz");

	imsl_e1sti (1, nz);

	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

    }

    else {

	for (i = 0; i < nz; i++)

	    tmp_vec[i] = n * (a[i].col - 1) + a[i].row;

	imsl_i_sort (nz, tmp_vec, IMSL_RETURN_USER, tmp_vec, 0);

	for (i = 0; i < nz - 1; i++) {

	    if (tmp_vec[i] == tmp_vec[i + 1]) {

		error++;

		ir = tmp_vec[i] % n;

		jc = (tmp_vec[i] - ir) / n + 1;

		imsl_e1sti (1, ir);

		imsl_e1sti (2, jc);

		imsl_ermes (IMSL_TERMINAL, IMSL_BAD_NONZERO);

	    }

	}

    }

    if (tmp_vec != NULL) {

	imsl_free (tmp_vec);

	tmp_vec = NULL;

    }

    }

    }



    if (return_symbolic_factor && supply_symbolic_factor) {

	error++;

	imsl_ermes (IMSL_TERMINAL, IMSL_EXCLUSIVE_SYMFAC_REQUEST);

    }



    if (return_numeric_factor && supply_numeric_factor) {

	error++;

	imsl_ermes (IMSL_TERMINAL, IMSL_EXCLUSIVE_NUMFAC_REQUEST);

    }



    if (symbolic_factor_only && (return_smallest_diag ||

	    return_largest_diag)) {

	error++;

	imsl_ermes (IMSL_TERMINAL, IMSL_BAD_DIAG_SIZE_REQUEST);

    }



    if (!compute_symbolic_factor && return_num_factor_nonzeros) {

	error++;

	imsl_ermes (IMSL_TERMINAL, IMSL_BAD_RETURN_NUMNZ_REQUEST);

    }



    if (error)

	goto RETURN;



    if (compute_symbolic_factor) {

/*  allocate space for the permutation and its inverse, to be returned by

    genqmd */



	perm = (Mint *) malloc (n * sizeof (*perm));

	invp = (Mint *) malloc (n * sizeof (*invp));





/*  allocate space to hold the return from sbmfct -- the off-diagonal

     non-zeros of the Cholesky factor in compressed subscript format */



	xlnz = (Mint *) malloc ((n + 1) * sizeof (*xlnz));

	xnzsub = (Mint *) malloc ((n + 1) * sizeof (*xnzsub));

	nzsub = (Mint *) malloc ((maxsub = nz) * sizeof (*nzsub));

	if (!perm || !invp || !xlnz || !xnzsub || !nzsub) {

	    if (nzsub != NULL) {

		free (nzsub);

		nzsub = NULL;

	    }

	    if (xnzsub != NULL) {

		free (xnzsub);

		xnzsub = NULL;

	    }

	    if (xlnz != NULL) {

		free (xlnz);

		xlnz = NULL;

	    }

	    if (invp != NULL) {

		free (invp);

		invp = NULL;

	    }

	    if (perm != NULL) {

		free (perm);

		perm = NULL;

	    }

	    error = 1;

	    imsl_e1stl (1, "n");

	    imsl_e1stl (2, "nz");

	    imsl_e1sti (1, n);

	    imsl_e1sti (2, nz);

	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);

	}



	if (error)

	    goto RETURN;



/*  now symbolic factorization */



	l_sym_fact_with_ordering (n, nz, a, &maxsub, &nzsub,

	    xnzsub, &maxlnz, xlnz, perm, invp, multifrontal, 

	    CSC_format, CSC_col_ptr, CSC_row_ind);



	if (return_num_factor_nonzeros)

	    *p_num_factor_nonzeros = maxlnz + n;



    }			     /* if compute_symbolic_factor  */



/*  setup to return the symbolic factorization, if necessary  */

    if (return_symbolic_factor) {

	p_symbolic_factor_struct->nzsub = &nzsub;

	p_symbolic_factor_struct->xnzsub = &xnzsub;

	p_symbolic_factor_struct->maxsub = maxsub;

	p_symbolic_factor_struct->xlnz = &xlnz;

	p_symbolic_factor_struct->maxlnz = maxlnz;

	p_symbolic_factor_struct->perm = &perm;

	p_symbolic_factor_struct->invp = &invp;

	p_symbolic_factor_struct->multifrontal_space = multifrontal_space;

    }



/*  get data from symbolic factor structure  */



    if (supply_symbolic_factor) {

	nzsub = *p_symbolic_factor_struct->nzsub;

	xnzsub = *p_symbolic_factor_struct->xnzsub;

	nzsub = *p_symbolic_factor_struct->nzsub;

	maxsub = p_symbolic_factor_struct->maxsub;

	xlnz = *p_symbolic_factor_struct->xlnz;

	maxlnz = p_symbolic_factor_struct->maxlnz;

	perm = *p_symbolic_factor_struct->perm;

	invp = *p_symbolic_factor_struct->invp;

	multifrontal_space = p_symbolic_factor_struct->multifrontal_space;

    }



    if (!symbolic_factor_only) {



	if (compute_numeric_factor) {



/*  change data from coordinate to compressed storage format; allocate

    storage to hold data structure */



	    diag = (Mfloat *) malloc (n * sizeof (*diag));

	    if (maxlnz > 0)

	    	alnz = (Mfloat *) calloc (maxlnz, sizeof (*alnz));

	    else

	    	alnz = (Mfloat *) calloc (1, sizeof (*alnz));

	    if (!alnz || !diag) {

		error++;

		imsl_e1stl (1, "n");

		imsl_e1sti (1, n);

		imsl_ermes (IMSL_TERMINAL,

		    IMSL_OUT_OF_MEMORY_1);

	    }



	    if (error) {

		if (alnz != NULL)

		    free (alnz);

		alnz = NULL;

		goto RETURN;

	    }



	    if (CSC_format)

		l_CSC_to_compressed (n, CSC_col_ptr, CSC_row_ind, CSC_values,

                    nzsub, xnzsub, xlnz, invp, alnz, diag, &flag);

	    else

	    	l_coordinate_to_compressed (nz, a, nzsub,

		    xnzsub, xlnz, invp, alnz, diag, &flag);



	    if (!(return_symbolic_factor))

		if (invp != NULL) {

		    free (invp);

		    invp = NULL;

		}



	    flag = 0;

	    if (multifrontal)

/*  perform multifrontal  numeric factorization, over the input data in alnz */



		l_mffct (&n, xlnz, alnz, xnzsub, nzsub, diag,

		    &multifrontal_space, &flag);

	    else

/*  perform standard numeric factorization, over the input data in alnz */



		l_gsfct (&n, xlnz, alnz, xnzsub, nzsub, diag,

		    &flag);

	}

	if (flag != 0) {

	    error++;

	    imsl_ermes (IMSL_TERMINAL, IMSL_BAD_SQUARE_ROOT);

	}

	if (error) {

	    if (diag != NULL) {

		free (diag);

		diag = NULL;

	    }

	    if (alnz != NULL) {

		free (alnz);

		alnz = NULL;

	    }

            if (nzsub != NULL) {

                free (nzsub);

                nzsub = NULL;

            }

            if (xnzsub != NULL) {

                free (xnzsub);

                xnzsub = NULL;

            }

            if (xlnz != NULL) {

                free (xlnz);

                xlnz = NULL;

            }

            if (invp != NULL) {

                free (invp);

                invp = NULL;

            }

            if (perm != NULL) {

                free (perm);

                perm = NULL;

            }

	    goto RETURN;

	}



/*  setup to return the numeric factorization, if necessary  */



	if (return_numeric_factor) {

	    p_numeric_factor_struct->nzsub = &nzsub;

	    p_numeric_factor_struct->xnzsub = &xnzsub;

	    p_numeric_factor_struct->xlnz = &xlnz;

	    p_numeric_factor_struct->alnz = &alnz;

	    p_numeric_factor_struct->perm = &perm;

	    p_numeric_factor_struct->diag = &diag;

	}



/*  get data from numeric factor_structure, if provided  */



	if (supply_numeric_factor) {

	    nzsub = *p_numeric_factor_struct->nzsub;

	    xnzsub = *p_numeric_factor_struct->xnzsub;

	    xlnz = *p_numeric_factor_struct->xlnz;

	    alnz = *p_numeric_factor_struct->alnz;

	    perm = *p_numeric_factor_struct->perm;

	    diag = *p_numeric_factor_struct->diag;

	}



	if (!numeric_factor_only) {



/*  gsslv overwrite the right hand side, so use the return value here  */



	    if (!lv_x)

		lv_x = (Mfloat *) malloc (n * sizeof (*lv_x));

	    memcpy ((char*)lv_x, (char*)rhs, n * sizeof (*lv_x));



/*  do forward and back solves */



	    l_gsslv (&n, xlnz, alnz, xnzsub, nzsub, perm, diag, lv_x);

	}



	if (return_smallest_diag) {

	    smallest_diag = imsl_f_machine (2);

	    for (i = 0; i < n; i++)

		smallest_diag = imsl_f_min (smallest_diag,

		    diag[i]);

	    *p_smallest_diag = smallest_diag;

	}



	if (return_largest_diag) {

	    largest_diag = 0.0;

	    for (i = 0; i < n; i++)

		largest_diag = imsl_f_max (largest_diag,

		    diag[i]);

	    *p_largest_diag = largest_diag;

	}



	if (!(return_numeric_factor || supply_numeric_factor)) {

	    if (alnz != NULL) {

		free (alnz);

		alnz = NULL;

	    }

	    if (diag != NULL) {

		free (diag);

		diag = NULL;

	    }

	}



    }			     /* if !symbolic_factor_only */



    if (!(return_symbolic_factor || supply_symbolic_factor)

	&& !(return_numeric_factor || supply_numeric_factor)) {

	if (nzsub != NULL) {

	    free (nzsub);

	    nzsub = NULL;

	}

	if (xnzsub != NULL) {

	    free (xnzsub);

	    xnzsub = NULL;

	}

	if (xlnz != NULL) {

	    free (xlnz);

	    xlnz = NULL;

	}

	if (perm != NULL) {

	    free (perm);

	    perm = NULL;

	}

    }



RETURN:;



    if (data_incremented) 

    for (i = 0; i < nz; i++) {

        a[i].row--;

        a[i].col--;

    }

    return argptr;

}



#ifdef ANSI

static void l_sym_fact_with_ordering (Mint n, Mint nz, Mf_sparse_elem *a,

	Mint *maxsub, Mint **nzsub, Mint *xnzsub, Mint *maxlnz, Mint *xlnz, 

	Mint *perm, Mint *invp, Mint multifrontal, Mint CSC_format,

	Mint *CSC_col_ptr, Mint *CSC_row_ind)

#else

static void l_sym_fact_with_ordering (n, nz, a, maxsub, nzsub, xnzsub,

                maxlnz, xlnz, perm, invp, multifrontal, CSC_format,

		CSC_col_ptr, CSC_row_ind)

    Mint         n;

    Mint         nz;

    Mf_sparse_elem *a;

    Mint        *maxsub;

    Mint       **nzsub;

    Mint        *xnzsub;

    Mint        *maxlnz;

    Mint        *xlnz;

    Mint        *perm;

    Mint        *invp;

    Mint         multifrontal;

    Mint        CSC_format;

    Mint       *CSC_col_ptr;

    Mint       *CSC_row_ind;

#endif

{

    Mint         iflag;

    Mint        *xadj = NULL;

    Mint        *adjncy = NULL;

    Mint        *xadj_copy = NULL;

    Mint        *adjncy_copy = NULL;



    Mint         nofsub;

    Mint         maxtmp;



    imsl_e1psh ("l_sym_fact_with_ordering");



/*  get space for the adjancecny structure; standard_to_adjacency converts

    our input format to adjacency format */



    xadj = (Mint *) malloc ((n + 1) * sizeof (*xadj));

    if (nz == n)

    	adjncy = (Mint *) malloc (sizeof (*adjncy));

    else

    	adjncy = (Mint *) malloc (2 * (nz - n) * sizeof (*adjncy));



    if (xadj == NULL || adjncy == NULL) {

	imsl_e1stl (1, "n");

        imsl_e1sti (1, n);

	imsl_e1stl (2, "nz");

        imsl_e1sti (2, nz);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);

	goto RETURN;

    }

	

    if (CSC_format)

        l_CSC_to_adjacency (n, CSC_col_ptr, CSC_row_ind, xadj,

            adjncy, &iflag);

    else

        l_standard_to_adjacency (nz, a, n, xadj, adjncy,

	    &iflag);



/* the minimum degree function destroys the adjacency structure,

   so make a copy for use in the symbolic factorization */



    xadj_copy = (Mint *) malloc ((n + 1) * sizeof (*xadj_copy));

    if (xadj_copy == NULL) {

	imsl_e1stl (1, "n");

        imsl_e1sti (1, n);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

	goto RETURN;

    }

    memcpy ((char*)xadj_copy, (char*)xadj, (n + 1) * sizeof (*xadj));



    if (nz == n) 

    	adjncy_copy = (Mint *) malloc (sizeof (*adjncy_copy));

    else {

    	adjncy_copy = (Mint *) malloc (2 * (nz - n) * sizeof (*adjncy_copy));

	if (adjncy_copy == NULL) {

            imsl_e1stl (1, "n");

            imsl_e1sti (1, n); 

	    imsl_e1stl (2, "nz");

            imsl_e1sti (2, nz);

            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);

            goto RETURN;

	}

    	memcpy ((char*)adjncy_copy, (char*)adjncy, 2 * (nz - n)

	    * sizeof (*adjncy));

    }



/*  do minimum degree */



    l_genqmd (&n, xadj_copy, adjncy_copy, perm, invp, &nofsub);



/*  free the extra copy of the adjancency structure */



    free (xadj_copy);

    free (adjncy_copy);



/*  now symbolic factorization */



    iflag = 1;

    while (iflag) {

	l_smbfct (&n, xadj, adjncy, perm, invp, xlnz, maxlnz,

	    xnzsub, *nzsub, maxsub, &iflag);

	if (iflag) {

	    *maxsub += n;

	    *nzsub = (Mint *) realloc (*nzsub, *maxsub * sizeof (**nzsub));

    	    if (*nzsub == NULL) {

        	imsl_e1stl (1, "n");

        	imsl_e1sti (1, n);

        	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        	goto RETURN;

    	    }

	}

    }



    if (multifrontal) {

    	maxtmp = *maxsub;

	if (nz > n) {

	    l_post_order (*nzsub, xnzsub, invp, perm, &n);

	    iflag = 1;

	    while (iflag) {

		l_smbfct (&n, xadj, adjncy, perm, invp, xlnz,

		    maxlnz, xnzsub, *nzsub, &maxtmp, &iflag);

		if (iflag) {

		    *maxsub = maxtmp += n;

		    *nzsub = (Mint *) realloc (*nzsub, maxtmp * sizeof (**nzsub));

    	    	    if (*nzsub == NULL) {

        		imsl_e1stl (1, "n");

        		imsl_e1sti (1, n);

        		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        		goto RETURN;

    	    	    }

		}

	    }

	}

	l_compute_frontal_storage (&n, xlnz, xnzsub, *nzsub,

	    &multifrontal_space);

    }



RETURN:



/*  we don't need the adjancency structure anymore */



    if (xadj != NULL) free (xadj);

    if (adjncy != NULL) free (adjncy);

    imsl_e1pop ("l_sym_fact_with_ordering");

}



/*Translated by FOR_C, v3.4 (P), on 02/04/93 at 10:41:13 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/04/93 at 10:35:44 */



/*  l_qmdmrg - quotient minimum degree merge

    

    this function merges indistinguisable nodes in the minimum degree

    ordering algorithm.  ia also computes the new degrees of these new

    supernodes.



    input - (xadj, adjncy) - the adjacency structure

            deg - the number of nodes in the given set

            (nhdsze, nbrhd) - the set of eliminated supernodes adjacent 

                              to some nodes in the set



    updated parameters - 

    deg - the degree vector

    qsize - size of indistinguishable nodes

    qlink - linked list for indistinguishable nodes

    marker - the given set is given by those nodes with marker value

             set to 1. those nodes with degree updated will have value 

             set to 2.

    

    working parameters -

    rchset - the reachable set

    ovrlp - temp vector to store the intersection of two reachable sets.

*/



#ifdef ANSI

static void l_qmdmrg (Mint xadj[], Mint adjncy[], Mint deg[],

                Mint qsize[], Mint qlink[], Mint marker[], Mint *deg0,

                Mint *nhdsze, Mint nbrhd[], Mint rchset[], Mint ovrlp[])

#else

static void l_qmdmrg (xadj, adjncy, deg, qsize, qlink, marker, deg0, nhdsze,

                nbrhd, rchset, ovrlp)

    Mint         xadj[];

    Mint         adjncy[];

    Mint         deg[];

    Mint         qsize[];

    Mint         qlink[];

    Mint         marker[];

    Mint        *deg0;

    Mint        *nhdsze;

    Mint         nbrhd[];

    Mint         rchset[];

    Mint         ovrlp[];

#endif

{

    Mint         deg1, head, inhd, inhd_, iov, iov_, irch, irch_, j, j_, jstop,

                jstrt, link, lnode, mark, mrgsze, nabor, node, novrlp, rchsze,

                root;

/*  initialization  */



    if (*nhdsze > 0) {

	for (inhd = 1; inhd <= *nhdsze; inhd++) {

	    inhd_ = inhd - 1;

	    root = nbrhd[inhd_];

	    marker[root - 1] = 0;

	}



/*  loop through each eliminated supernode in the set (nhdsze,nbrhd)  */



	for (inhd = 1; inhd <= *nhdsze; inhd++) {

	    inhd_ = inhd - 1;

	    root = nbrhd[inhd_];

	    marker[root - 1] = -1;

	    rchsze = 0;

	    novrlp = 0;

	    deg1 = 0;

    L_200:

	    jstrt = xadj[root - 1];

	    jstop = xadj[root] - 1;



/*  determine the reachable set and its intersection with the   */

/*  input reachable set                                         */



	    for (j = jstrt; j <= jstop; j++) {

		j_ = j - 1;

		nabor = adjncy[j_];

		root = -nabor;

		if (nabor < 0)

		    goto L_200;

		if (!nabor)

		    goto L_700;



		mark = marker[nabor - 1];

		if (mark >= 0) {

		    if (mark <= 0) {

			rchsze += 1;

			rchset[rchsze - 1] = nabor;

			deg1 += qsize[nabor - 1];

			marker[nabor - 1] = 1;

		    }

		    else if (mark <= 1) {

			novrlp += 1;

			ovrlp[novrlp - 1] = nabor;

			marker[nabor - 1] = 2;

		    }

		}

	    }



/*  from the overlapped set, determine the nodes that can be merged together */



    L_700:

	    head = 0;

	    mrgsze = 0;

	    for (iov = 1; iov <= novrlp; iov++) {

		iov_ = iov - 1;

		node = ovrlp[iov_];

		jstrt = xadj[node - 1];

		jstop = xadj[node] - 1;

		for (j = jstrt; j <= jstop; j++) {

		    j_ = j - 1;

		    nabor = adjncy[j_];

		    if (!marker[nabor - 1])

			goto L_1401;

		}



/*  node belongs to the new merged supernode.  update the vectors  */

/*  qlink and qsize                                                */



		mrgsze += qsize[node - 1];

		marker[node - 1] = -1;

		lnode = node;

		while (TRUE) {

		    link = qlink[lnode - 1];

		    if (link <= 0)

			goto L_1000;

		    lnode = link;

		}

	L_1000:

		qlink[lnode - 1] = head;

		head = node;

		goto L_1402;

	L_1401:

		marker[node - 1] = 1;

	L_1402:

		;

	    }

	    if (head > 0) {

		qsize[head - 1] = mrgsze;

		deg[head - 1] = *deg0 + deg1 - 1;

		marker[head - 1] = 2;

	    }



/*  reset marker values  */



	    root = nbrhd[inhd_];

	    marker[root - 1] = 0;

	    if (rchsze > 0) {

		for (irch = 1; irch <= rchsze; irch++) {

		    irch_ = irch - 1;

		    node = rchset[irch_];

		    marker[node - 1] = 0;

		}

	    }

	}

    }

    return;

}



/*Translated by FOR_C, v3.4 (P), on 02/04/93 at 10:41:15 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/04/93 at 10:35:49

   Options SET: none



     QMDQT  ..... QUOT MIN DEG QUOT TRANSFORM



     PURPOSE - THIS SUBROUTINE PERFORMS THE QUOTIENT GRAPH

        TRANSFORMATION AFTER A NODE HAS BEEN ELIMINATED.



     INPUT PARAMETERS -

        ROOT - THE NODE JUST ELIMINATED. IT BECOMES THE

               REPRESENTATIVE OF THE NEW SUPERNODE.

        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.

        (RCHSZE, RCHSET) - THE REACHABLE SET OF ROOT IN THE

               OLD QUOTIENT GRAPH.

        NBRHD - THE NEIGHBORHOOD SET WHICH WILL BE MERGED

               WITH ROOT TO FORM THE NEW SUPERNODE.

        MARKER - THE MARKER VECTOR.



     UPDATED PARAMETER -

        ADJNCY - BECOMES THE ADJNCY OF THE QUOTIENT GRAPH.

*/



#ifdef ANSI

static void l_qmdqt (Mint *root, Mint xadj[], Mint adjncy[],

                Mint marker[], Mint *rchsze, Mint rchset[], Mint nbrhd[])

#else

static void l_qmdqt (root, xadj, adjncy, marker, rchsze, rchset, nbrhd)

    Mint        *root;

    Mint         xadj[];

    Mint         adjncy[];

    Mint         marker[];

    Mint        *rchsze;

    Mint         rchset[];

    Mint         nbrhd[];

#endif

{

    Mint         inhd, irch, irch_, j, j_, jstop, jstrt, link, nabor, node;

    irch = 0;

    inhd = 0;

    node = *root;

    while (TRUE) {

	jstrt = xadj[node - 1];

	jstop = xadj[node] - 2;

	if (jstop >= jstrt) {



/*  place reach nodes into the adjacent list of node  */



	    for (j = jstrt; j <= jstop; j++) {

		j_ = j - 1;

		irch += 1;

		adjncy[j_] = rchset[irch - 1];

		if (irch >= *rchsze)

		    goto L_400;

	    }

	}



/*  link to other space provided by the nbrhd set  */



	link = adjncy[jstop];

	node = -link;

	if (link >= 0) {

	    inhd += 1;

	    node = nbrhd[inhd - 1];

	    adjncy[jstop] = -node;

	}

    }



/*  all reachable nodes have been saved. end the adj list   */

/*  add root to the nbr list of each node in the reach set  */



L_400:

    adjncy[j] = 0;

    for (irch = 1; irch <= *rchsze; irch++) {

	irch_ = irch - 1;

	node = rchset[irch_];

	if (marker[node - 1] >= 0) {

	    jstrt = xadj[node - 1];

	    jstop = xadj[node] - 1;

	    for (j = jstrt; j <= jstop; j++) {

		j_ = j - 1;

		nabor = adjncy[j_];

		if (marker[nabor - 1] < 0)

		    goto L_602;

	    }

	    goto L_601;

    L_602:

	    adjncy[j - 1] = *root;

	}

L_601:

	;

    }

    return;

}



/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/04/93 at 10:35:51

 *  Options SET: none



     QMDRCH ..... QUOT MIN DEG REACH SET



     PURPOSE - THIS SUBROUTINE DETERMINES THE REACHABLE SET OF

        A NODE THROUGH A GIVEN SUBSET.  THE ADJACENCY STRUCTURE

        IS ASSUMED TO BE STORED IN A QUOTIENT GRAPH FORMAT.



     INPUT PARAMETERS -

        ROOT - THE GIVEN NODE NOT IN THE SUBSET.

        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE PAIR.

        DEG - THE DEGREE VECTOR.  DEG(I) LT 0 MEANS THE NODE

               BELONGS TO THE GIVEN SUBSET.



     OUTPUT PARAMETERS -

        (RCHSZE, RCHSET) - THE REACHABLE SET.

        (NHDSZE, NBRHD) - THE NEIGHBORHOOD SET.



     UPDATED PARAMETERS -

        MARKER - THE MARKER VECTOR FOR REACH AND NBRHD SETS.

               GT 0 MEANS THE NODE IS IN REACH SET.

               LT 0 MEANS THE NODE HAS BEEN MERGED WITH

               OTHERS IN THE QUOTIENT OR IT IS IN NBRHD SET.

*/



#ifdef ANSI

static void l_qmdrch (Mint *root, Mint xadj[], Mint adjncy[],

                Mint deg[], Mint marker[], Mint *rchsze, Mint rchset[],

                Mint *nhdsze, Mint nbrhd[])

#else

static void l_qmdrch (root, xadj, adjncy, deg, marker, rchsze, rchset, nhdsze,

                nbrhd)

    Mint        *root;

    Mint         xadj[];

    Mint         adjncy[];

    Mint         deg[];

    Mint         marker[];

    Mint        *rchsze;

    Mint         rchset[];

    Mint        *nhdsze;

    Mint         nbrhd[];

#endif

{

    Mint         i, i_, istop, istrt, j, j_, jstop, jstrt, nabor, node;



/*  loop through the neighbors of root in the quotient graph  */



    *nhdsze = 0;

    *rchsze = 0;

    istrt = xadj[*root - 1];

    istop = xadj[*root] - 1;

    if (istop >= istrt) {

	for (i = istrt; i <= istop; i++) {

	    i_ = i - 1;

	    nabor = adjncy[i_];

	    if (!nabor)

		goto L_601;

	    if (!marker[nabor - 1]) {

		if (deg[nabor - 1] < 0) {



/*  nabor has been eliminated, find nodes reachable from it  */



		    marker[nabor - 1] = -1;

		    *nhdsze += 1;

		    nbrhd[*nhdsze - 1] = nabor;

	    L_300:

		    jstrt = xadj[nabor - 1];

		    jstop = xadj[nabor] - 1;

		    for (j = jstrt; j <= jstop; j++) {

			j_ = j - 1;

			node = adjncy[j_];

			nabor = -node;

			if (node < 0)

			    goto L_300;

			if (!node)

			    goto L_602;

			if (!marker[node - 1]) {

			    *rchsze += 1;

			    rchset[*rchsze - 1] = node;

			    marker[node - 1] = 1;

			}

		    }

		}

		else {



/*  include nabor into the reachable set  */



		    *rchsze += 1;

		    rchset[*rchsze - 1] = nabor;

		    marker[nabor - 1] = 1;

		}

	    }

    L_602:

	    ;

	}

    }

L_601:

    return;

}



/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/04/93 at 10:35:53

 *  Options SET: none



     QMDUPD ..... QUOT MIN DEG UPDATE

     PURPOSE - THIS ROUTINE PERFORMS DEGREE UPDATE FOR A SET

        OF NODES IN THE MINIMUM DEGREE ALGORITHM.



     INPUT PARAMETERS -

        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.

        (NLIST, LIST) - THE LIST OF NODES WHOSE DEGREE HAS TO

               BE UPDATED.



     UPDATED PARAMETERS -

        DEG - THE DEGREE VECTOR.

        QSIZE - SIZE OF INDISTINGUISHABLE SUPERNODES.

        QLINK - LINKED LIST FOR INDISTINGUISHABLE NODE.

        MARKER - USED TO MARK THOSE NODES IN REACH/NBRHD SETS.



     WORKING PARAMETERS -

        RCHSET - THE REACHABLE SET.

        NBRHD -  THE NEIGHBORHOOD SET. 



     PROGRAM SUBROUTINES - 

        QMDMRG.

*/



#ifdef ANSI

static void l_qmdupd (Mint xadj[], Mint adjncy[], Mint *nlist,

                Mint list[], Mint deg[], Mint qsize[], Mint qlink[],

                Mint marker[], Mint rchset[], Mint nbrhd[])

#else

static void l_qmdupd (xadj, adjncy, nlist, list, deg, qsize, qlink, marker,

                rchset, nbrhd)

    Mint         xadj[];

    Mint         adjncy[];

    Mint        *nlist;

    Mint         list[];

    Mint         deg[];

    Mint         qsize[];

    Mint         qlink[];

    Mint         marker[];

    Mint         rchset[];

    Mint         nbrhd[];

#endif

{

    Mint         deg0, deg1, il, il_, inhd, inhd_, inode, irch, irch_, j,

                j_, jstop, jstrt, mark, nabor, nhdsze, node, rchsze;



/*  find all eliminated supernodes that are adjacent to some  */

/*  nodes in the given list.  put them into (nhdsze, nbrhd).  */

/*  dego contains the number of nodes in the list             */



    if (*nlist > 0) {

	deg0 = 0;

	nhdsze = 0;

	for (il = 1; il <= *nlist; il++) {

	    il_ = il - 1;

	    node = list[il_];

	    deg0 += qsize[node - 1];

	    jstrt = xadj[node - 1];

	    jstop = xadj[node] - 1;

	    for (j = jstrt; j <= jstop; j++) {

		j_ = j - 1;

		nabor = adjncy[j_];

		if (!(marker[nabor - 1] || deg[nabor - 1] >= 0)) {

		    marker[nabor - 1] = -1;

		    nhdsze += 1;

		    nbrhd[nhdsze - 1] = nabor;

		}

	    }

	}



/*  merge indistinguishable nodes in the list by calling qmdmrg  */



	if (nhdsze > 0)

	    l_qmdmrg (xadj, adjncy, deg, qsize, qlink, marker, &deg0,

		&nhdsze, nbrhd, rchset, &nbrhd[nhdsze]);



/*  find the new degrees of the nodes that havenot been merged  */



	for (il = 1; il <= *nlist; il++) {

	    il_ = il - 1;

	    node = list[il_];

	    mark = marker[node - 1];

	    if (!(mark > 1 || mark < 0)) {

		marker[node - 1] = 2;

		l_qmdrch (&node, xadj, adjncy, deg, marker, &rchsze,

		    rchset, &nhdsze, nbrhd);

		deg1 = deg0;

		if (rchsze > 0) {

		    for (irch = 1; irch <= rchsze; irch++) {

			irch_ = irch - 1;

			inode = rchset[irch_];

			deg1 += qsize[inode - 1];

			marker[inode - 1] = 0;

		    }

		}

		deg[node - 1] = deg1 - 1;

		if (nhdsze > 0) {

		    for (inhd = 1; inhd <= nhdsze; inhd++) {

			inhd_ = inhd - 1;

			inode = nbrhd[inhd_];

			marker[inode - 1] = 0;

		    }

		}

	    }

	}

    }

    return;

}



#ifdef ANSI

static void l_standard_to_adjacency (Mint nz, Mf_sparse_elem *a, Mint n,

	Mint *xadj, Mint *adjncy, Mint *iflag)

#else

static void l_standard_to_adjacency (nz, a, n, xadj, adjncy,

                iflag)

    Mint         nz;

    Mf_sparse_elem *a;

    Mint         n;

    Mint        *xadj;

    Mint        *adjncy;

    Mint        *iflag;

#endif

{

    Mint         i;

    Mint         ir;

    Mint         ic;

    Mint        *ixwork = NULL;



    imsl_e1psh ("l_standard_to_adjacency");

    ixwork = (Mint *) malloc (n * sizeof (*ixwork));

    if (ixwork == NULL) {

        imsl_e1stl (1, "n");

        imsl_e1sti (1, n);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        goto RETURN;

    }



    *iflag = 0;

    if (nz <= 0) {

    	imsl_e1pop ("l_standard_to_adjacency");

	return;

    }

    for (i = 0; i < n + 1; i++)

	(xadj)[i] = -2;

    for (i = 0; i < n; i++)

	ixwork[i] = 0;

    for (i = 0; i < nz; i++) {

	(xadj)[a[i].row] += 1;

	(xadj)[a[i].col] += 1;

    }

    (xadj)[0] = 1;

    for (i = 1; i < n + 1; i++)

	(xadj)[i] += (xadj)[i - 1];



    for (i = 0; i < nz; i++) {

	if (a[i].row != a[i].col) {

	    ir = a[i].row - 1;

	    ic = a[i].col - 1;



	    (adjncy)[(xadj)[ir] + ixwork[ir] - 1] = ic + 1;

	    ixwork[ir] += 1;



	    (adjncy)[(xadj)[ic] + ixwork[ic] - 1] = ir + 1;

	    ixwork[ic] += 1;

	}

    }



RETURN:

    if (ixwork != NULL) free (ixwork);

    imsl_e1pop ("l_standard_to_adjacency");

}



#ifdef ANSI

static void l_CSC_to_adjacency (Mint n, Mint *CSC_col_ptr, Mint *CSC_row_ind,

	Mint *xadj, Mint *adjncy, Mint *iflag)

#else

static void l_CSC_to_adjacency (n, CSC_col_ptr, CSC_row_ind, xadj,

		adjncy, iflag)

    Mint         n;

    Mint	       *CSC_col_ptr;

    Mint	       *CSC_row_ind;

    Mint        *xadj;

    Mint        *adjncy;

    Mint        *iflag;

#endif

{

    Mint		start;

    Mint		stop;

    Mint         i;

    Mint		j;

    Mint         ir;

    Mint         ic;

    Mint        *ixwork = NULL;



    imsl_e1psh ("l_CSC_to_adjacency");

    ixwork = (Mint *) malloc (n * sizeof (*ixwork));

    if (ixwork == NULL) {

        imsl_e1stl (1, "n");

        imsl_e1sti (1, n);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        goto RETURN;

    }

 

    *iflag = 0;

    for (i = 0; i < n + 1; i++)

        (xadj)[i] = -2;

    for (i = 0; i < n; i++)

        ixwork[i] = 0;



    for (i = 0; i < n; i++) {

	start = CSC_col_ptr[i];

	stop = CSC_col_ptr[i + 1];

	for (j = start; j < stop; j++) {

	    ir = CSC_row_ind[j];

	    ic = i;

            (xadj)[ir+1] += 1;

            (xadj)[ic+1] += 1;

	}

    }



    (xadj)[0] = 1;

    for (i = 1; i < n + 1; i++)

        (xadj)[i] += (xadj)[i - 1];

 

    for (i = 0; i < n; i++) {

	start = CSC_col_ptr[i];

	stop = CSC_col_ptr[i + 1];

	for (j = start; j < stop; j++) {

	    ir = CSC_row_ind[j];

	    ic = i;

            if (ir != ic) {

               

            	(adjncy)[(xadj)[ir] + ixwork[ir] - 1] = ic + 1;

            	ixwork[ir] += 1;

 

            	(adjncy)[(xadj)[ic] + ixwork[ic] - 1] = ir + 1;

            	ixwork[ic] += 1;

	    }

        }

    }



RETURN:

    if (ixwork != NULL) free (ixwork);

    imsl_e1pop ("l_CSC_to_adjacency");

}



/*Translated by FOR_C, v3.4 (P), on 02/04/93 at 10:41:29 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/04/93 at 10:36:01

 *  Options SET: none



      SMBFCT ..... SYMBOLIC FACTORIZATION



      PURPOSE - THIS ROUTINE PERFORMS SYMBOLIC FACTORIZATION 

         ON A PERMUTED LINEAR SYSTEM AND IT ALSO SETS UP THE

         COMPRESSED DATA STRUCTURE FOR THE SYSTEM.



      INPUT PARAMETERS -

         NEQNS - NUMBER OF EQUATIONS.

         (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.

         (PERM, INVP) - THE PERMUTATION VECTOR AND ITS INVERSE.



      UPDATED PARAMETERS -

         MAXSUB - SIZE OF THE SUBSCRIPT ARRAY NZSUB.  ON RETURN,

                IT CONTAINS THE NUMBER OF SUBSCRIPTS USED



      OUTPUT PARAMETERS -

         XLNZ - INDEX INTO THE NONZERO STORAGE VECTOR LNZ.

         (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT VECTORS.

         MAXLNZ - THE NUMBER OF NONZEROS FOUND.

         FLAG - ERROR FLAG.  POSITIVE VALUE INDICATES THAT.

                NZSUB ARRAY IS TOO SMALL.



      WORKING PARAMETERS -

         MRGLNK - A VECTOR OF SIZE NEQNS.  AT THE KTH STEP,

                MRGLNK(K), MRGLNK(MRGLNK(K)) , .........

                IS A LIST CONTAINING ALL THOSE COLUMNS L(*,J)

                WITH J LESS THAN K, SUCH THAT ITS FIRST OFF-

                DIAGONAL NONZERO IS L(K,J).  THUS, THE

                NONZERO STRUCTURE OF COLUMN L(*,K) CAN BE FOUND

                BY MERGING THAT OF SUCH COLUMNS L(*,J) WITH

                THE STRUCTURE OF A(*,K).

         RCHLNK - A VECTOR OF SIZE NEQNS.  IT IS USED TO ACCUMULATE

                THE STRUCTURE OF EACH COLUMN L(*,K).  AT THE

                END OF THE KTH STEP,

                    RCHLNK(K), RCHLNK(RCHLNK(K)), ........

                IS THE LIST OF POSITIONS OF NONZEROS IN COLUMN K

                OF THE FACTOR L.

         MARKER  - AN INTEGER VECTOR OF LENGTH NEQNS. IT IS USED 

               TO TEST IF MASS SYMBOLIC ELIMINATION CAN BE

                PERFORMED.  THAT IS, IT IS USED TO CHECK WHETHER

                THE STRUCTURE OF THE CURRENT COLUMN K BEING

                PROCESSED IS COMPLETELY DETERMINED BY THE SINGLE

                COLUMN MRGLNK(K).



*/



#ifdef ANSI

static void l_smbfct (Mint *neqns, Mint xadj[], Mint adjncy[],

                Mint perm[], Mint invp[], Mint xlnz[], Mint *maxlnz,

                Mint xnzsub[], Mint nzsub[], Mint *maxsub, Mint *flag)

#else

static void l_smbfct (neqns, xadj, adjncy, perm, invp, xlnz, maxlnz, xnzsub,

                nzsub, maxsub, flag)

    Mint        *neqns;

    Mint         xadj[];

    Mint         adjncy[];

    Mint         perm[];

    Mint         invp[];

    Mint         xlnz[];

    Mint        *maxlnz;

    Mint         xnzsub[];

    Mint         nzsub[];

    Mint        *maxsub;

    Mint        *flag;

#endif

{

    Mint         i, inz, j, j_, jstop, jstrt, jstrt_, k, k_, knz, kxsub, lmax,

                m, mrgk, mrkflg, nabor, node, np1, nzbeg, nzend, rchm;



    Mint        *rchlnk = NULL, *mrglnk = NULL, *marker = NULL;



    imsl_e1psh ("l_smbfct");

    rchlnk = (Mint *) malloc (*neqns * sizeof (*rchlnk));

    mrglnk = (Mint *) malloc (*neqns * sizeof (*mrglnk));

    marker = (Mint *) malloc (*neqns * sizeof (*marker));

    if (rchlnk == NULL || mrglnk == NULL || marker == NULL) {

        imsl_e1stl (1, "n");

        imsl_e1sti (1, *neqns);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        goto RETURN;

    }



/*  initialization  */



    nzbeg = 1;

    nzend = 0;

    xlnz[0] = 1;

    for (k = 1; k <= *neqns; k++) {

	k_ = k - 1;

	mrglnk[k_] = 0;

	marker[k_] = 0;

    }

/*  for each column, knz counts the number of nonzeros in column   */

/*  k accumulated in rchlnk                                        */



    np1 = *neqns + 1;

    for (k = 1; k <= *neqns; k++) {

	k_ = k - 1;

	knz = 0;

	mrgk = mrglnk[k_];

	mrkflg = 0;

	marker[k_] = k;

	if (mrgk)

	    marker[k_] = marker[mrgk - 1];

	xnzsub[k_] = nzend;

	node = perm[k_];

	jstrt = xadj[node - 1];

	jstop = xadj[node] - 1;

	if (jstrt <= jstop) {



/*  use rchlnk to link through the structure of A(*, k) below diagonal  */



	    rchlnk[k_] = np1;

	    for (j = jstrt; j <= jstop; j++) {

		j_ = j - 1;

		nabor = adjncy[j_];

		nabor = invp[nabor - 1];

		if (nabor > k) {

		    rchm = k;

		    while (TRUE) {

			m = rchm;

			rchm = rchlnk[m - 1];

			if (rchm > nabor)

			    goto L_1601;

		    }

	    L_1601:

		    knz += 1;

		    rchlnk[m - 1] = nabor;

		    rchlnk[nabor - 1] = rchm;

		    if (marker[nabor - 1] != marker[k_])

			mrkflg = 1;

		}

	    }



/*  test for mass symbolic elimination  */



	    lmax = 0;

	    if (!(mrkflg || !mrgk)) {

		if (!mrglnk[mrgk - 1]) {

		    xnzsub[k_] = xnzsub[mrgk - 1] + 1;

		    knz = xlnz[mrgk] - (xlnz[mrgk - 1] + 1);

		    goto L_1400;

		}

	    }



/*  link through each column i that affects L(*,k) */



	    i = k;

	    while (TRUE) {

		i = mrglnk[i - 1];

		if (!i)

		    goto L_800;

		inz = xlnz[i] - (xlnz[i - 1] + 1);

		jstrt = xnzsub[i - 1] + 1;

		jstop = xnzsub[i - 1] + inz;

		if (inz > lmax) {

		    lmax = inz;

		    xnzsub[k_] = jstrt;

		}



/*  merge structure of (L(*,i) in nzsub into rchlnk  */



		rchm = k;

		for (j = jstrt; j <= jstop; j++) {

		    j_ = j - 1;

		    nabor = nzsub[j_];

		    while (TRUE) {

			m = rchm;

			rchm = rchlnk[m - 1];

			if (rchm >= nabor)

			    goto L_1602;

		    }

	    L_1602:

		    if (rchm != nabor) {

			knz += 1;

			rchlnk[m - 1] = nabor;

			rchlnk[nabor - 1] = rchm;

			rchm = nabor;

		    }

		}

	    }



/*  check if subscripts duplicate those of another column  */



    L_800:

	    if (knz != lmax) {



/*  or if tail of k-1st column matches head of k-th  */



		if (nzbeg <= nzend) {

		    i = rchlnk[k_];

		    for (jstrt = nzbeg; jstrt <= nzend; jstrt++) {

			jstrt_ = jstrt - 1;

			if (nzsub[jstrt_] == i)

			    goto L_1000;

			if (nzsub[jstrt_] > i)

			    goto L_1200;

		    }

		    goto L_1200;

	    L_1000:

		    xnzsub[k_] = jstrt;

		    for (j = jstrt; j <= nzend; j++) {

			j_ = j - 1;

			if (nzsub[j_] != i)

			    goto L_1200;

			i = rchlnk[i - 1];

			if (i > *neqns)

			    goto L_1400;

		    }

		    nzend = jstrt - 1;

		}



/*  copy the structure of L(*,k) from rchlnk to the data structure  */

/*  (xnzsub, nzsub)                                                 */



	L_1200:

		nzbeg = nzend + 1;

		nzend += knz;

		if (nzend > *maxsub)

		    goto L_1600;

		i = k;

		for (j = nzbeg; j <= nzend; j++) {

		    j_ = j - 1;

		    i = rchlnk[i - 1];

		    nzsub[j_] = i;

		    marker[i - 1] = k;

		}

		xnzsub[k_] = nzbeg;

		marker[k_] = k;

	    }



/*  update the vector mrglnk.  note column L(*,k) just found is required to  */

/*  determine column L(*,j), where L(j,k) is the first nonzero in L(*,k)     */

/*  below the diagonal                                                       */



    L_1400:

	    if (knz > 1) {

		kxsub = xnzsub[k_];

		i = nzsub[kxsub - 1];

		mrglnk[k_] = mrglnk[i - 1];

		mrglnk[i - 1] = k;

	    }

	}

	xlnz[k_ + 1] = xlnz[k_] + knz;

    }

    *maxlnz = xlnz[*neqns - 1] - 1;

    *maxsub = xnzsub[*neqns - 1];

    xnzsub[*neqns] = xnzsub[*neqns - 1];

    *flag = 0;

    goto L_1603;



/*  error - insufficient storage for nonzero subscripts  */



L_1600:

    *flag = 1;

L_1603:



RETURN:

    if (rchlnk != NULL) free (rchlnk);

    if (mrglnk != NULL) free (mrglnk);

    if (marker != NULL) free (marker);

    imsl_e1pop ("l_smbfct");

}



/*Translated by FOR_C, v3.4 (P), on 02/04/93 at 10:41:06 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/04/93 at 10:35:38

 *  Options SET: none



     GENQMD ..... QUOT MIN DEGREE ORDERING



      PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE

         ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENT-

         ATION OF THE ELIMINATION GRAPHS BY QUOTIENT GRAPHS,

         AND THE NOTION OF INDISTINGUISHABLE NODES.

         CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE

         DESTROYED.



      INPUT PARAMETERS -

         NEQNS - NUMBER OF EQUATIONS.

         (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.



      OUTPUT PARAMETERS -

         PERM - THE MINIMUM DEGREE ORDERING.

         INVP - THE INVERSE OF PERM.



      WORKING PARAMETERS -

         DEG - THE DEGREE VECTOR. DEG(I) IS NEGATIVE MEANS

                NODE I HAS BEEN NUMBERED.

         MARKER - A MARKER VECTOR, WHERE MARKER(I) IS 

                NEGATIVE MEANS NODE I HAS BEEN MERGED WITH

                ANOTHER NODE AND THUS CAN BE IGNORED.

         RCHSET - VECTOR USED FOR THE REACHABLE SET.

         NBRHD - VECTOR USED FOR THE NEIGHBORHOOD SET.

         QSIZE - VECTOR USED TO STORE THE SIZE OF

                INDISTINGUISHABLE SUPERNODES.

         QLINK - VECTOR TO STORE INDISTINGUISHABLE NODES,

                I, QLINK(I), QLINK(QLINK(I)) ... ARE THE

                MEMBERS OF THE SUPERNODE REPRESENTED BY I.



      PROGRAM SUBROUTINES -

         QMDRCH, QMDQT, QMDUPD.



 */



#ifdef ANSI

static void l_genqmd (Mint *neqns, Mint xadj[], Mint adjncy[],

                Mint perm[], Mint invp[],

                Mint *nofsub)

#else

static void l_genqmd (neqns, xadj, adjncy, perm, invp, nofsub)

    Mint        *neqns;

    Mint         xadj[];

    Mint         adjncy[];

    Mint        *perm;

    Mint        *invp;

    Mint        *nofsub;

#endif

{

    Mint         inode, ip, irch, irch_, j, j_, mindeg, ndeg, nhdsze, node,

                node_, np, num, nump1, nxnode, rchsze, search, thresh;



    Mint        *deg = NULL;

    Mint        *marker = NULL;

    Mint        *rchset = NULL;

    Mint        *nbrhd = NULL;

    Mint        *qsize = NULL;

    Mint        *qlink = NULL;



/*  initialize degree vector and other working variables */



/* working storage */



    imsl_e1psh ("l_genqmd");

    deg = (Mint *) malloc (*neqns * sizeof (*deg));

    marker = (Mint *) malloc (*neqns * sizeof (*marker));

    rchset = (Mint *) malloc (*neqns * sizeof (*rchset));

    nbrhd = (Mint *) malloc (*neqns * sizeof (*nbrhd));

    qsize = (Mint *) malloc (*neqns * sizeof (*qsize));

    qlink = (Mint *) malloc ((*neqns + 1) * sizeof (*qlink));

    if (deg == NULL || marker == NULL || rchset == NULL || nbrhd == NULL

	|| qsize == NULL || qlink == NULL) {

        imsl_e1stl (1, "n");

        imsl_e1sti (1, *neqns);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        goto RETURN;

    }



    mindeg = *neqns;

    *nofsub = 0;

    for (node = 1; node <= *neqns; node++) {

	node_ = node - 1;

	(perm)[node_] = node;

	(invp)[node_] = node;

	marker[node_] = 0;

	qsize[node_] = 1;

	qlink[node_] = 0;

	ndeg = xadj[node_ + 1] - xadj[node_];

	deg[node_] = ndeg;

	if (ndeg < mindeg)

	    mindeg = ndeg;

    }

    num = 0;



/*  perform threshold search to get a node of min degree.

    variable search points to where search should start */



    while (TRUE) {

	search = 1;

	thresh = mindeg;

	mindeg = *neqns;

L_300:

	nump1 = num + 1;

	if (nump1 > search)

	    search = nump1;

	for (j = search; j <= *neqns; j++) {

	    j_ = j - 1;

	    node = (perm)[j_];

	    if (marker[node - 1] >= 0) {

		ndeg = deg[node - 1];

		if (ndeg <= thresh)

		    goto L_500;

		if (ndeg < mindeg)

		    mindeg = ndeg;

	    }

	}

	goto L_801;



/*  node has minimum degree -- find its reachable sets by calling

    qmdrch */



L_500:

	search = j;

	*nofsub += deg[node - 1];

	marker[node - 1] = 1;

	l_qmdrch (&node, xadj, adjncy, deg, marker, &rchsze, rchset,

	    &nhdsze, nbrhd);



/*  eliminate all nodes indistinguishable from node -- they are

    given by node, qlink(node), ...  */



	nxnode = node;

	while (TRUE) {

	    num += 1;

	    np = (invp)[nxnode - 1];

	    ip = (perm)[num - 1];

	    (perm)[np - 1] = ip;

	    (invp)[ip - 1] = np;

	    (perm)[num - 1] = nxnode;

	    (invp)[nxnode - 1] = num;

	    deg[nxnode - 1] = -1;

	    nxnode = qlink[nxnode - 1];

	    if (nxnode <= 0)

		goto L_802;

	}

	/* 113. */

L_802:

	if (rchsze > 0) {



/*  update the degrees of the nodes in the reachable set

    and identify the indistinguishable nodes */



	    l_qmdupd (xadj, adjncy, &rchsze, rchset, deg, qsize, qlink,

		marker, &rchset[rchsze], &nbrhd[nhdsze]);



/*  reset marker value of nodes in reach set.

    update threshold value for cyclic search.

    also call qmdqt to form new quotient graph  */



	    marker[node - 1] = 0;

	    for (irch = 1; irch <= rchsze; irch++) {

		irch_ = irch - 1;

		inode = rchset[irch_];

		if (marker[inode - 1] >= 0) {

		    marker[inode - 1] = 0;

		    ndeg = deg[inode - 1];

		    if (ndeg < mindeg)

			mindeg = ndeg;

		    if (ndeg <= thresh) {

			mindeg = thresh;

			thresh = ndeg;

			search = (invp)[inode - 1];

		    }

		}

	    }

	    if (nhdsze > 0)

		l_qmdqt (&node, xadj, adjncy, marker, &rchsze, rchset,

		    nbrhd);

	}

	if (num < *neqns)

	    goto L_300;

	goto L_803;

L_801:

	;

    }

L_803:



RETURN:



/* free the work vectors */



    if (deg != NULL) free (deg);

    if (marker != NULL) free (marker);

    if (rchset != NULL) free (rchset);

    if (nbrhd != NULL) free (nbrhd);

    if (qsize != NULL) free (qsize);

    if (qlink != NULL) free (qlink);

    imsl_e1pop ("l_genqmd");

    return;

}





#ifdef ANSI

static void l_coordinate_to_compressed (Mint nz, Mf_sparse_elem *a,

	Mint *nzsub, Mint *xnzsub, Mint *xlnz, Mint *invp, Mfloat *alnz,

	Mfloat *diag, Mint *iflag)

#else

static void l_coordinate_to_compressed (nz, a, nzsub, xnzsub,

                xlnz, invp, alnz, diag, iflag)

    Mint         nz;

    Mf_sparse_elem *a;

    Mint        *nzsub;

    Mint        *xnzsub;

    Mint        *xlnz;

    Mint        *invp;

    Mfloat     *alnz;

    Mfloat     *diag;

    Mint        *iflag;

#endif

{

    Mint         i;

    Mint         ir, ic;

    Mint         irp, icp;

    Mint         itemp, jtemp;

    Mint         ii;

    Mint         k;

    *iflag = 0;

    for (i = 0; i < nz; i++) {

	ir = a[i].row;

	ic = a[i].col;



/* check for diagonal */



	if (ir == ic) {

	    jtemp = invp[ic - 1];

	    diag[jtemp - 1] = a[i].val;

	}

	else {

	    irp = invp[ir - 1];

	    icp = invp[ic - 1];

	    if (irp <= icp) {

		jtemp = irp;

		irp = icp;

		icp = jtemp;

	    }

	    itemp = xnzsub[icp - 1];

	    k = xlnz[icp] - xlnz[icp - 1];

	    for (ii = itemp; ii <= itemp + k - 1; ii++) {

		if (nzsub[ii - 1] == irp) {

		    alnz[xlnz[icp - 1] + ii - itemp - 1] = a[i].val;

		    *iflag = 1;

		}

	    }

	    if (*iflag != 1) {

		*iflag = 2;

		break;

	    }

	}

    }

}



#ifdef ANSI

static void l_CSC_to_compressed (Mint n, Mint *CSC_col_ptr, Mint *CSC_row_ind,

	Mfloat *CSC_values, Mint *nzsub, Mint *xnzsub, Mint *xlnz, Mint *invp,

	Mfloat *alnz, Mfloat *diag, Mint *iflag)

#else

static void l_CSC_to_compressed (n, CSC_col_ptr, CSC_row_ind, CSC_values, 

		nzsub, xnzsub, xlnz, invp, alnz, diag, iflag)

    Mint         n;

    Mint	       *CSC_col_ptr;

    Mint	       *CSC_row_ind;

    Mfloat     *CSC_values;

    Mint        *nzsub;

    Mint        *xnzsub;

    Mint        *xlnz;

    Mint        *invp;

    Mfloat     *alnz;

    Mfloat     *diag;

    Mint        *iflag;

#endif

{

    Mint         i;

    Mint         j;

    Mint         ir, ic;

    Mint         irp, icp;

    Mint         itemp, jtemp;

    Mint         ii;

    Mint         k;

    Mint		start;

    Mint		stop;



    *iflag = 0;

    for (i = 0; i < n; i++) {

	start = CSC_col_ptr[i];

	stop = CSC_col_ptr[i + 1];

	for (j = start; j < stop; j++) {



            ir = CSC_row_ind[j] + 1;

            ic = i + 1;

 

/* check for diagonal */

 

            if (ir == ic) {

            	jtemp = invp[ic - 1];

            	diag[jtemp - 1] = CSC_values[j];

            }

            else {

            	irp = invp[ir - 1];

            	icp = invp[ic - 1];

            	if (irp <= icp) {

                    jtemp = irp;

                    irp = icp;

                    icp = jtemp;

            	}

            	itemp = xnzsub[icp - 1];

            	k = xlnz[icp] - xlnz[icp - 1];

            	for (ii = itemp; ii <= itemp + k - 1; ii++) {

                    if (nzsub[ii - 1] == irp) {

                    	alnz[xlnz[icp - 1] + ii - itemp - 1] = CSC_values[j];

                    	*iflag = 1;

                    }

            	}

            	if (*iflag != 1) {

                    *iflag = 2;

                    break;

            	}

            }

        }

    }

}



/*Translated by FOR_C, v3.4 (P), on 02/11/93 at 09:41:42 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/11/93 at 09:38:48

 *  Options SET: fmt=t s=n

 *-----------------------------------------------------------------------

 *  IMSL Name:  L6CXD (Single precision version)

 *

 *  Computer:   SUN4/SINGLE

 *

 *  Revised:    May 11, 1990

 *

 *  Purpose:    Post order the topologically ordered tree generated by

 *              L5CXD.

 *

 *  Usage:      CALL L6CXD (MAXSUB, NZSUB, NZSUBX, INVP, IPERM, N,

 *                          ICHILD, ISKCNT, ISKNOD, ISKOST, IXCHLD)

 *

 *  Arguments:  See L5CXD.

 *

 *  Chapter:    MATH/LIBRARY Linear Systems

 *

 *  Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

 *

 *  Warranty:   IMSL warrants only that IMSL testing has been applied

 *              to this code.  No other warranty, expressed or implied,

 *              is applicable.

 *

 *-----------------------------------------------------------------------

 * */

#ifdef ANSI

static void l_post_order (Mint nzsub[], Mint nzsubx[],

                Mint invp[], Mint iperm[], Mint *n)

#else

static void l_post_order (nzsub, nzsubx, invp, iperm, n)

    Mint         nzsub[];

    Mint         nzsubx[];

    Mint         invp[];

    Mint         iperm[];

    Mint        *n;

#endif

{

    Mint         i, i_, ioldnd, ip, iroot, iskptr, istart, itemp1, itemp2,

                iter, j, j_, label, newnde;

    Mint        *ichild = NULL;

    Mint        *iskcnt = NULL;

    Mint        *isknod = NULL;

    Mint        *iskost = NULL;

    Mint        *ixchld = NULL;



    imsl_e1psh ("l_post_order");



/*  allocate work vectors  */



    ichild = (Mint *) malloc ((*n + 1) * sizeof (*ichild));

    iskcnt = (Mint *) malloc ((*n + 1) * sizeof (*iskcnt));

    isknod = (Mint *) malloc ((*n + 1) * sizeof (*isknod));

    iskost = (Mint *) malloc ((*n + 1) * sizeof (*iskost));

    ixchld = (Mint *) malloc ((*n + 1) * sizeof (*ixchld));

    if (ichild == NULL || iskcnt == NULL || isknod == NULL ||

	iskost == NULL || ixchld == NULL) {

        imsl_e1stl (1, "n");

        imsl_e1sti (1, *n);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        goto RETURN;

    }



/*  set up child information, node i will have ixchild[i+1]-ixchild[i]

    children starting at ichild[ixchild[i]]  */



    memset ((char*)ixchld, 0, (*n + 1) * sizeof (*ixchld));

    for (j = 1; j <= *n; j++) {

	j_ = j - 1;

	ip = nzsub[nzsubx[j_] - 1];

	if (ip && ip != j)

	    ixchld[ip - 1] += 1;

    }



    itemp1 = ixchld[0];

    itemp2 = ixchld[1];

    ixchld[0] = 1;

    for (j = 1; j <= (*n - 1); j++) {

	j_ = j - 1;

	ixchld[j_ + 1] = ixchld[j_] + itemp1;

	itemp1 = itemp2;

	itemp2 = ixchld[j_ + 2];

    }

    ixchld[*n] = ixchld[*n - 1] + itemp1;



/*  use iskopt as temporary storage vector  */



    memcpy ((char*)iskost, (char*)ixchld, *n * sizeof (*iskost));

    for (j = 1; j <= *n; j++) {

	j_ = j - 1;

	ip = nzsub[nzsubx[j_] - 1];

	if (ip && ip != j) {

	    ichild[iskost[ip - 1] - 1] = j;

	    iskost[ip - 1] += 1;

	}

    }



    iter = 1;

    for (j = 1; j <= (*n - 1); j++) {

	j_ = j - 1;

	ip = nzsub[nzsubx[j_] - 1];

	if (!ip || ip == j)

	    iter += 1;

    }



    label = 1;

    istart = 1;

    for (j = 1; j <= iter; j++) {

	j_ = j - 1;



/*  find root of tree  */



	if (j == 1) {

	    iroot = *n;

	}

	else {

	    for (i = istart; i <= *n; i++) {

		i_ = i - 1;

		ip = nzsub[nzsubx[i_] - 1];

		if (!ip || ip == i) {

		    iroot = i;

		    goto L_60;

		}

	    }

    L_60:

	    ;

	    istart = iroot + 1;

	}



/*  ispkpt - current std position

    isknod - contains the nodes

    iskcnt - contains the number of childre that have not been pushed yet

    iskost - index offset into ichild indicating which child is on the stack */



	iskptr = 0;



/*  push root onto stack  */



	iskptr = 1;

	isknod[iskptr - 1] = iroot;

	iskcnt[iskptr - 1] = ixchld[iroot] - ixchld[iroot - 1];

	iskost[iskptr - 1] = -1;



L_70:

	;



	if (iskcnt[iskptr - 1]) {



/*  if a node has children that have not been processed, decrement iskcnt,

    increment offset.  next push the current child (as indicated by offset)

    along with iskcnt and iskost  */



	    iskcnt[iskptr - 1] -= 1;

	    iskost[iskptr - 1] += 1;

	    iskptr += 1;

	    isknod[iskptr - 1] = ichild[ixchld[isknod[iskptr - 2] -

		    1] + iskost[iskptr - 2] - 1];

	    iskcnt[iskptr - 1] = ixchld[isknod[iskptr - 1]] - ixchld[isknod[iskptr - 1] -

		1];

	    iskost[iskptr - 1] = -1;

	    if (!ixchld[isknod[iskptr - 1] - 1])

		iskcnt[iskptr - 1] = 0;

	}

	else {





/*  if a node has count zero, all children have been processed.  process that

    node.  save the old node and rename the node  */



	    if (iskptr > 1) {

		ioldnd = ichild[ixchld[isknod[iskptr - 2] - 1] + iskost[iskptr - 2] -

		    1];

		if (ioldnd == label)

		    goto L_80;

		ichild[ixchld[isknod[iskptr - 2] - 1] + iskost[iskptr - 2] -

		    1] = label;

	    }

	    else {

		ioldnd = iroot;

	    }

	    newnde = label;



/*  make the necessary changes in invp  */



	    invp[iperm[ioldnd - 1] - 1] = newnde;



    L_80:

	    label += 1;



/*  if at the bottom of stack, get out; else pop and continue  */



	    if (iskptr == 1)

		goto L_90;

	    iskptr -= 1;

	}

	goto L_70;

L_90:

	;

    }



    for (i = 1; i <= *n; i++) {

	i_ = i - 1;

	iperm[invp[i_] - 1] = i;

    }



RETURN:



    if (ichild != NULL) free (ichild);

    if (iskcnt != NULL) free (iskcnt);

    if (isknod != NULL) free (isknod);

    if (iskost != NULL) free (iskost);

    if (ixchld != NULL) free (ixchld);



    imsl_e1pop ("l_post_order");

    return;

}



/*Translated by FOR_C, v3.4 (P), on 02/11/93 at 09:41:28 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/11/93 at 09:38:39

 *  Options SET: fmt=t s=n

 *-----------------------------------------------------------------------

 *  IMSL Name:  L11XD (Single precision version)

 *

 *  Computer:   SUN4/SINGLE

 *

 *  Revised:    February 28, 1991

 *

 *  Purpose:    Determine the required size of working storage

 *              for stack of frontal matrices.

 *

 *  Usage:      CALL L11XD (NEQNS, IXLNZ, NZSUBX, NZSUB, ISTACK, MAXWS)

 *

 *  Arguments:

 *     NEQNS  - Number of equation.  (Input)

 *     IXLNZ  - Index vector for LNZ.  (Input)

 *     NZSUBX - Part of the compressed subscript data structure for L.

 *              (Input)

 *     NZSUB  - Part of the compressed subscript data structure for L.

 *              (Input)

 *     ISTACK - Stack of nodes, each corresponds to a frontal matrix

 *              on the stack storage WSTORE.  (Input)

 *     MAXWS  - Required size of working storage for stack of

 *              frontal matrices.  (Output)

 *

 *  Chapter:    MATH/LIBRARY Linear Systems

 *

 *  Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

 *

 *  Warranty:   IMSL warrants only that IMSL testing has been applied

 *              to this code.  No other warranty, expressed or implied,

 *              is applicable.

 *

 *-----------------------------------------------------------------------

 * */

#ifdef ANSI

static void l_compute_frontal_storage (Mint *neqns, Mint ixlnz[], Mint nzsubx[],

                Mint nzsub[], Mint *maxws)

#else

static void l_compute_frontal_storage (neqns, ixlnz, nzsubx, nzsub, maxws)

    Mint        *neqns;

    Mint         ixlnz[];

    Mint         nzsubx[];

    Mint         nzsub[];

    Mint        *maxws;

#endif

{

    Mint         i, i_, ifront, im1frn, im1par, ipop1, isknod, iskpar, isktop,

                ixfree, j, lastj, loc, ndloc, nelim, nextj, nmat, nmatnz,

                nsub, nsubnz;



    Mint        *istack = NULL;



    imsl_e1psh ("l_compute_frontal_storage");

    istack = (Mint *) malloc (*neqns * sizeof (*istack));

    if (istack == NULL) {

        imsl_e1stl (1, "n");

        imsl_e1sti (1, *neqns);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        goto RETURN;

    }

    if (*neqns <= 0) {

    	imsl_e1pop ("l_compute_frontal_storage");

	return;

    }



/*  initialize working variables: isktop points to top of stack, isknod is

    the top node at stack whose parent is given by iskpar, ixfree points

    to start of free space */



    *maxws = 0;

    istack[0] = 0;

    isktop = 2;

    isknod = 0;

    iskpar = 0;

    ixfree = 1;

    j = 1;



L_10:

    ;

    nextj = j + 1;

    nmat = ixlnz[nextj - 1] - ixlnz[j - 1] + 1;

    nmatnz = nmat * (nmat + 1) / 2;

    if (iskpar == j)

	goto L_20;

    *maxws = max (*maxws, ixfree + nmatnz - 1);

    goto L_50;

L_20:

    ;

    ndloc = nzsubx[j - 1] - 1;

    ipop1 = 1;

L_30:

    ;

    nsub = ixlnz[isknod] - ixlnz[isknod - 1];

    nsubnz = nsub * (nsub + 1) / 2;

    loc = nzsubx[isknod - 1];

    l_match_vect (&nsub, &nzsub[loc - 1], &nzsub[ndloc - 1], &istack[isktop - 1]);

    ixfree -= nsubnz;

    if (!ipop1)

	goto L_40;

    *maxws = max (*maxws, ixfree + nmatnz - 1);

    ipop1 = 0;

L_40:

    ;

    isktop -= 1;

    isknod = istack[isktop - 2];

    iskpar = 0;

    if (isknod <= 0)

	goto L_50;

    loc = nzsubx[isknod - 1];

    iskpar = nzsub[loc - 1];

    if (iskpar == j)

	goto L_30;

L_50:

    ;

    nelim = 1;

    im1frn = ixlnz[nextj - 1] - ixlnz[j - 1];

    if (im1frn <= 0)

	goto L_70;

    loc = nzsubx[j - 1];

    im1par = nzsub[loc - 1];

    for (i = nextj; i <= *neqns; i++) {

	i_ = i - 1;

	if (im1par != i || iskpar == i)

	    goto L_70;

	ifront = ixlnz[i_ + 1] - ixlnz[i_];

	if (ifront != im1frn - 1)

	    goto L_70;

	nelim += 1;

	if (ifront <= 0)

	    goto L_70;

	loc = nzsubx[i_];

	im1par = nzsub[loc - 1];

	im1frn = ifront;

    }



L_70:

    ;

    lastj = j + nelim - 1;

    nextj = lastj + 1;

    nsub = nmat - nelim;

    if (nsub <= 0)

	goto L_80;

    nsubnz = nsub * (nsub + 1) / 2;

    ixfree += nsubnz;

    istack[isktop - 1] = lastj;

    isktop += 1;

    isknod = lastj;

    loc = nzsubx[lastj - 1];

    iskpar = nzsub[loc - 1];

L_80:

    ;

    j = nextj;

    if (j <= *neqns)

	goto L_10;



RETURN:



    if (istack != NULL) free (istack);

    imsl_e1pop ("l_compute_frontal_storage");

    return;

}



/*Translated by FOR_C, v3.4 (P), on 02/11/93 at 09:41:46 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/11/93 at 09:38:50

 *  Options SET: fmt=t s=n

 *-----------------------------------------------------------------------

 *  IMSL Name:  L7FXD (Single precision version)

 *

 *  Computer:   SUN4/SINGLE

 *

 *  Revised:    March 7, 1990

 *

 *  Purpose:    Match a subvector with another vector.  The locations

 *              of its values in the super vector will be returned.  The

 *              first entries of subvector and vector are assumed to

 *              be the same.

 *

 *  Usage:      CALL L7FXD (NSUB, ISUBVC, MSTVEC, IXTEND)

*

 *  Arguments:

 *     NSUB   - Size of subvector.  (Input)

 *     ISUBVC - Subvector, values are in ascending sequence.  (Input)

 *     MSTVEC - The master vector, values are in ascending sequence.

 *              (Input)

 *     IXTEND - Subvector location in super vector, length is NSUB.

 *              (Output)

 *

 *  Chapter:    MATH/LIBRARY Linear Systems

 *

 *  Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

 *

 *  Warranty:   IMSL warrants only that IMSL testing has been applied

 *              to this code.  No other warranty, expressed or implied,

 *              is applicable.

 *

 *-----------------------------------------------------------------------

 * */

#ifdef ANSI

static void l_match_vect (Mint *nsub, Mint isubvc[], Mint mstvec[],

                Mint ixtend[])

#else

static void l_match_vect (nsub, isubvc, mstvec, ixtend)

    Mint        *nsub;

    Mint         isubvc[];

    Mint         mstvec[];

    Mint         ixtend[];

#endif

{

    Mint         isub, isub_, isubvl, ivec;

    ixtend[0] = 1;

    if (*nsub <= 1)

	return;



    ivec = 1;

    for (isub = 2; isub <= *nsub; isub++) {

	isub_ = isub - 1;

	isubvl = isubvc[isub_];

L_10:

	ivec += 1;

	if (mstvec[ivec - 1] < isubvl)

	    goto L_10;

	ixtend[isub_] = ivec;

    }

    return;

}



/*Translated by FOR_C, v3.4 (P), on 02/04/93 at 10:41:10 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/04/93 at 10:35:41

 *  Options SET: none



      GSFCT ..... GENERAL SPARSE SYMMETRIC FACT



      PURPOSE - THIS SUBROUTINE PERFORMS THE SYMMETRIC

         FACTORIZATION FOR A GENERAL SPARSE SYSTEM, STORED IN 

         THE COMPRESSED SUBSCRIPT DATA FORMAT.



      INPUT PARAMETERS -

         NEQNS - NUMBER OF EQUATIONS.

         XLNZ - INDEX VECTOR FOR LNZ.  XLNZ(I) POINTS TO THE

                START OF NONZEROS IN COLUMN I OF FACTOR L.

         (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT DATA

                STRUCTURE FOR FACTOR L.



      UPDATED PARAMETERS -

         LNZ - ON INPUT, CONTAINS NONZEROS OF A, AND ON

                RETURN, THE NONZEROS OF L.

         DIAG - THE DIAGONAL OF L OVERWRITES THAT OF A.

         IFLAG - THE ERROR FLAG.  IT IS SET TO 1 IF A ZERO OR

                NEGATIVE SQUARE ROOT OCCURS DURING THE

                FACTORIZATION.

         OPS   - A DOUBLE PRECISION COMMON PARAMETER THAT IS

                 INCREMENTED BY THE NUMBER OF OPERATIONS

                 PERFORMED BY THE SUBROUTINE.



      WORKING PARAMETERS -

         LINK - AT STEP J, THE LIST IN

                   LINK(J), LINK(LINK(J)), ...........

                CONSISTS OF THOSE COLUMNS THAT WILL MODIFY

                THE COLUMN L(*,J).

         FIRST - TEMPORARY VECTOR TO POINT TO THE FIRST

                NONZERO IN EACH COLUMN THAT WILL BE USED

                NEXT FOR MODIFICATION.

         TEMP - A TEMPORARY VECTOR TO ACCUMULATE MODIFICATIONS.



*/



#ifdef ANSI

static void l_gsfct (Mint *neqns, Mint xlnz[], Mfloat lnz[],

                Mint xnzsub[], Mint nzsub[], Mfloat diag[],

                Mint *iflag)

#else

static void l_gsfct (neqns, xlnz, lnz, xnzsub, nzsub, diag, iflag)

    Mint        *neqns;

    Mint         xlnz[];

    Mfloat      lnz[];

    Mint         xnzsub[];

    Mint         nzsub[];

    Mfloat      diag[];

    Mint        *iflag;

#endif

{

    Mint         i, i_, ii, ii_, istop, istrt, isub, j, j_, k, kfirst, newk;

    Mfloat      diagj, ljk;

    Mint        *link = NULL;

    Mint        *first = NULL;

    Mfloat     *temp = NULL;



    imsl_e1psh ("l_gsfct");



    link = (Mint *) malloc (*neqns * sizeof (*link));

    first = (Mint *) malloc (*neqns * sizeof (*first));

    temp = (Mfloat *) malloc (*neqns * sizeof (*temp));

    if (link == NULL || first == NULL || temp == NULL) {

        imsl_e1stl (1, "n");

        imsl_e1sti (1, *neqns);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        goto RETURN;

    }



/*  initialize work vectors */



    for (i = 1; i <= *neqns; i++) {

	i_ = i - 1;

	link[i_] = 0;

	temp[i_] = 0.0e0;

    }



/*  compute column L(*, j) for j = 1, ..., neqns  */



    for (j = 1; j <= *neqns; j++) {

	j_ = j - 1;



/*  for each column L(*,k) tjat affects L(*,j) */



	diagj = 0.0e0;

	newk = link[j_];

	while (TRUE) {

	    k = newk;

	    if (!k)

		goto L_400;

	    newk = link[k - 1];



/*  outer product modification of L(*,j) by */

/*  L(*,k) starting ar first[k] of L(*,k)   */



	    kfirst = first[k - 1];

	    ljk = lnz[kfirst - 1];

	    diagj += ljk * ljk;

	    istrt = kfirst + 1;

	    istop = xlnz[k] - 1;

	    if (istop >= istrt) {



/*  before modification, update vectors first,  */

/* and link for future modification steps       */



		first[k - 1] = istrt;

		i = xnzsub[k - 1] + (kfirst - xlnz[k - 1]) + 1;

		isub = nzsub[i - 1];

		link[k - 1] = link[isub - 1];

		link[isub - 1] = k;



/*  the actual mod is saved in vector temp  */



		for (ii = istrt; ii <= istop; ii++) {

		    ii_ = ii - 1;

		    isub = nzsub[i - 1];

		    temp[isub - 1] += lnz[ii_] * ljk;

		    i += 1;

		}

	    }

	}



/*  apply the modifications accumulated in temp to column L(*,j)  */



L_400:

	diagj = diag[j_] - diagj;

	if (diagj <= 0.0e0)

	    goto L_700;

	diagj = sqrt (diagj);

	diag[j_] = diagj;

	istrt = xlnz[j_];

	istop = xlnz[j_ + 1] - 1;

	if (istop >= istrt) {

	    first[j_] = istrt;

	    i = xnzsub[j_];

	    isub = nzsub[i - 1];

	    link[j_] = link[isub - 1];

	    link[isub - 1] = j;

	    for (ii = istrt; ii <= istop; ii++) {

		ii_ = ii - 1;

		isub = nzsub[i - 1];

		lnz[ii_] = (lnz[ii_] - temp[isub - 1]) / diagj;

		temp[isub - 1] = 0.0e0;

		i += 1;

	    }

	}

    }

    goto L_701;



/* error - zero or negative square root in factorization */



L_700:

    *iflag = 1;

L_701:



RETURN:



    if (link != NULL) free (link);

    if (first != NULL) free (first);

    if (temp != NULL) free (temp);

    imsl_e1pop ("l_gsfct");

    return;

}





/*Translated by FOR_C, v3.4 (P), on 02/04/93 at 10:41:12 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/04/93 at 10:35:42

 *  Options SET: none



      GSSLV ..... GENERAL SPARSE SYMMETRIC SOLVE



      PURPOSE - TO PERFORM SOLUTION OF A FACTORED SYSTEM, WHERE

         THE MATRIX IS STORED IN THE COMPRESSED SUBSCRIPT

         SPARSE FORMAT.



      INPUT PARAMETERS -

         NEQNS - NUMBER OF EQUATIONS.

         (XLNZ, LNZ) - STRUCTURE OF NONZEROS IN L.

         (XNZSUB, NZSUB) - COMPRESSED SUBSCRIPT STRUCTURE.

         DIAG - DIAGONAL COMPONENTS OF L.



      UPDATED PARAMETER -

         RHS - ON INPUT, IT CONTAINS THE RHS VECTOR, AND ON

                OUTPUT, THE SOLUTION VECTOR.



*/



#ifdef ANSI

static void l_gsslv (Mint *neqns, Mint xlnz[], Mfloat lnz[],

                Mint xnzsub[], Mint nzsub[], Mint perm[], Mfloat diag[], Mfloat rhs[])

#else

static void l_gsslv (neqns, xlnz, lnz, xnzsub, nzsub, perm, diag, rhs)

    Mint        *neqns;

    Mint         xlnz[];

    Mfloat      lnz[];

    Mint         xnzsub[];

    Mint         nzsub[];

    Mint         perm[];

    Mfloat      diag[];

    Mfloat      rhs[];

#endif

{

    Mint         i, ii, ii_, istop, istrt, isub, j, j_, jj;

    Mfloat      rhsj, s;

/*  forward substitution  */



/* apply permutation */



    l_permu (*neqns, rhs, perm, 1, rhs);



    for (j = 1; j <= *neqns; j++) {

	j_ = j - 1;

	rhsj = rhs[j_] / diag[j_];

	rhs[j_] = rhsj;

	istrt = xlnz[j_];

	istop = xlnz[j_ + 1] - 1;

	if (istop >= istrt) {

	    i = xnzsub[j_];

	    for (ii = istrt; ii <= istop; ii++) {

		ii_ = ii - 1;

		isub = nzsub[i - 1];

		rhs[isub - 1] += -lnz[ii_] * rhsj;

		i += 1;

	    }

	}

    }



/*  backward substitution  */



    j = *neqns;

    for (jj = 1; jj <= *neqns; jj++) {

	s = rhs[j - 1];

	istrt = xlnz[j - 1];

	istop = xlnz[j] - 1;

	if (istop >= istrt) {

	    i = xnzsub[j - 1];

	    for (ii = istrt; ii <= istop; ii++) {

		ii_ = ii - 1;

		isub = nzsub[i - 1];

		s += -lnz[ii_] * rhs[isub - 1];

		i += 1;

	    }

	}

	rhs[j - 1] = s / diag[j - 1];

	j -= 1;

    }



    l_permu (*neqns, rhs, perm, 2, rhs);

    return;

}



#ifdef ANSI

static void l_permu (Mint n, Mfloat *x, Mint *ipermu, Mint ipath, Mfloat *xpermu)

#else

static void l_permu (n, x, ipermu, ipath, xpermu)

    Mint         n;

    Mfloat      x[];

    Mint         ipermu[], ipath;

    Mfloat      xpermu[];

#endif

{

    Mint         i, j, k;

    Mfloat      temp;

/*  make a copy of x in xpermu and work with xpermu  */



    memcpy ((char*)xpermu, (char*)x, n * sizeof (*x));



    if (n == 1)

	goto L_9000;



    for (i = 0; i < n; i++)

	ipermu[i] = -ipermu[i];



    if (ipath == 1) {



/* forward permutation */



	for (i = 1; i <= n; i++) {

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



/* backward permutation */



	for (i = 1; i <= n; i++) {

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

    return;

}



/*Translated by FOR_C, v3.4 (P), on 02/11/93 at 09:41:41 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/11/93 at 09:38:47

 *  Options SET: fmt=t s=n

 *-----------------------------------------------------------------------

 *  IMSL Name:  L5FXD/DL5FXD (Single/Double precision version)

 *

 *  Computer:   SUN4/SINGLE

 *

 *  Revised:    March 7, 1990

 *

 *  Purpose:    Perform the multifrontal factorization for a sparse

 *              symmetric system, stored in the compress subscript

 *              data format.

 *

 *  Usage:      CALL L5FXD (NEQNS, IXLNZ, ALNZ, NZSUBX, NZSUB, DIAG,

 *                          ISTACK, WSTORE, MAXWS, IFLAG)

 *

 *  Arguments:

 *     NEQNS  - Number of equation.  (Input)

 *     IXLNZ  - Index vector for ALNZ.  (Input)

 *     ALNZ   - On input contains nonzeros of A, on output contains

 *              the nonzeros of L.  (Input/Output)

 *     NZSUBX - Part of the compressed subscript data structure for L.

 *              (Input)

 *     NZSUB  - Part of the compressed subscript data structure for L.

 *              (Input)

 *     DIAG   - On input contains the diagonal of A, on output contains

 *              the diagonal of L.  (Input/Output)

 *     ISTACK - Stack of nodes, each corresponds to a frontal matrix

 *              on the stack storage WSTORE.  (Input/Output)

 *              Tail portion is also used in the routine for extension

 *              vector.

 *     WSTORE - Working storage for stack of frontal matrices.

 *              (Input/Output)

 *     MAXWS  - Required size of wstore vector.  (Input)

 *     IFLAG  - Error flag set to 1 is a zero or negative square root

 *              occurs during the factorization.  (Output)

 *              Set to 2 if stack storage is not large enough.

 *

 *  Chapter:    MATH/LIBRARY Linear Systems

 *

 *  Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

 *

 *  Warranty:   IMSL warrants only that IMSL testing has been applied

 *              to this code.  No other warranty, expressed or implied,

 *              is applicable.

 *

 *-----------------------------------------------------------------------

 * */

#ifdef ANSI

static void l_mffct (Mint *neqns, Mint ixlnz[], Mfloat alnz[],

                Mint nzsubx[], Mint nzsub[], Mfloat diag[],

                Mint *maxws, Mint *iflag)

#else

static void l_mffct (neqns, ixlnz, alnz, nzsubx, nzsub, diag,

                maxws, iflag)

    Mint        *neqns;

    Mint         ixlnz[];

    Mfloat      alnz[];

    Mint         nzsubx[];

    Mint         nzsub[];

    Mfloat      diag[];

    Mint        *maxws;

    Mint        *iflag;

#endif

{

    Mint         i, i_, ifront, im1frn, im1par, ipop1, isknod, iskpar, isktop,

                ixfree, ixmat, j, k, k_, kstop, kstrt, lastj, loc, matloc,

                ndloc, nelim, nextj, nmat, nmatnz, nsub, nsubnz;

    Mint         loop;



    Mint        *istack = NULL;

    Mfloat     *wstore = NULL;



    imsl_e1psh ("l_mffct");



/*  initialize working variables:  isktop points to top of stack,

    sknode is the top nade at stack whose parent is given by iskpar.

    ixfree points to start of free space  */



    istack = (Mint *) malloc (*neqns * sizeof (*neqns));

    wstore = (Mfloat *) malloc (*maxws * sizeof (*wstore));

    if (istack == NULL || wstore == NULL) {

        imsl_e1stl (1, "n");

        imsl_e1sti (1, *neqns);

        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

        goto RETURN;

    }



    istack[0] = 0;

    isktop = 2;

    isknod = 0;

    iskpar = 0;

    ixfree = 1;

    j = 1;



/*  for each node j, do ...  */



L_10:

    ;

    nextj = j + 1;

    nmat = ixlnz[nextj - 1] - ixlnz[j - 1] + 1;

    nmatnz = nmat * (nmat + 1) / 2;

    if (iskpar == j)

	goto L_20;





/*  if j is a leaf, allocate space for its frontal matrix  */



    matloc = ixfree;

    if (ixfree + nmatnz - 1 > *maxws)

	goto L_150;

/*

	sset( &nmatnz, ADR(_f0,0.0e0), &wstore[ixfree - 1], ADR(_l0,1) );

*/

    for (loop = 0; loop < nmatnz; loop++)

	wstore[ixfree - 1 + loop] = 0.0;



    goto L_60;

L_20:

    ;



/*  otherwise, extend the top frontal matrix from stack and add subsequent

    frontal matrix from stack if the parent node of stack matrix is 'j'  */



    ndloc = nzsubx[j - 1] - 1;

    ipop1 = 1;

L_30:

    ;

    nsub = ixlnz[isknod] - ixlnz[isknod - 1];

    nsubnz = nsub * (nsub + 1) / 2;

    loc = nzsubx[isknod - 1];

    l_match_vect (&nsub, &nzsub[loc - 1], &nzsub[ndloc - 1], &istack[isktop - 1]);

    ixfree -= nsubnz;

    if (!ipop1)

	goto L_40;

    matloc = ixfree;

    if (ixfree + nmatnz - 1 > *maxws)

	goto L_150;

    l_l9fxd (&nmat, &wstore[ixfree - 1], &nsub, &istack[isktop - 1]);

    ipop1 = 0;

    goto L_50;

L_40:

    ;

    l_l6fxd (&nmat, &wstore[matloc - 1], &nsub, &wstore[ixfree - 1],

	&istack[isktop - 1]);

L_50:

    ;

    isktop -= 1;

    isknod = istack[isktop - 2];

    iskpar = 0;

    if (isknod <= 0)

	goto L_60;

    loc = nzsubx[isknod - 1];

    iskpar = nzsub[loc - 1];

    if (iskpar == j)

	goto L_30;

L_60:

    ;



/*  determine the size of mass elim associated with col j  */



    nelim = 1;

    im1frn = ixlnz[nextj - 1] - ixlnz[j - 1];

    if (im1frn <= 0)

	goto L_80;

    loc = nzsubx[j - 1];

    im1par = nzsub[loc - 1];

    for (i = nextj; i <= *neqns; i++) {

	i_ = i - 1;

	if (im1par != i || iskpar == i)

	    goto L_80;

	ifront = ixlnz[i_ + 1] - ixlnz[i_];

	if (ifront != im1frn - 1)

	    goto L_80;

	nelim += 1;

	if (ifront <= 0)

	    goto L_80;

	loc = nzsubx[i_];

	im1par = nzsub[loc - 1];

	im1frn = ifront;

    }



/*  add contributions from original matrix columns to the current frontal

    matrix (at wstore[ixfree]).  then perform 'nelim' steps of eliminations  */



L_80:

    ;

    lastj = j + nelim - 1;

    ixmat = matloc;

    for (i = j; i <= lastj; i++) {

	i_ = i - 1;

	wstore[ixmat - 1] += diag[i_];

	ixmat += 1;

	kstrt = ixlnz[i_];

	kstop = ixlnz[i_ + 1] - 1;

	if (kstop < kstrt)

	    goto L_100;

	for (k = kstrt; k <= kstop; k++) {

	    k_ = k - 1;

	    wstore[ixmat - 1] += alnz[k_];

	    ixmat += 1;

	}

L_100:

	;

    }

    l_l8fxd (&nmat, &wstore[matloc - 1], &nelim, iflag);

    if (*iflag == 1)

	goto L_140;



/*  save the computed factor coulmns into storage for L  */



    ixmat = matloc;

    for (i = j; i <= lastj; i++) {

	i_ = i - 1;

	diag[i_] = wstore[ixmat - 1];

	ixmat += 1;

	kstrt = ixlnz[i_];

	kstop = ixlnz[i_ + 1] - 1;

	if (kstop < kstrt)

	    goto L_120;

	for (k = kstrt; k <= kstop; k++) {

	    k_ = k - 1;

	    alnz[k_] = wstore[ixmat - 1];

	    ixmat += 1;

	}

L_120:

	;

    }



/*  push remaining frontal matrix into stack  */



    nextj = lastj + 1;

    nsub = nmat - nelim;

    if (nsub <= 0)

	goto L_130;

    nsubnz = nsub * (nsub + 1) / 2;

/*

	scopy( &nsubnz, &wstore[ixmat - 1], ADR(_l0,1), &wstore[ixfree - 1],

	 ADR(_l1,1) );

*/

    memcpy ((char*)&wstore[ixfree - 1], (char*)&wstore[ixmat - 1],

	nsubnz * sizeof (*wstore));

    ixfree += nsubnz;

    istack[isktop - 1] = lastj;

    isktop += 1;

    isknod = lastj;

    loc = nzsubx[lastj - 1];

    iskpar = nzsub[loc - 1];

L_130:

    ;

    j = nextj;

    if (j <= *neqns)

	goto L_10;

L_140:

    if (istack != NULL) free (istack);

    if (wstore != NULL) free (wstore);

    imsl_e1pop ("l_mffct");

    return;



/*  error -- insufficient stack storage  */



L_150:

    *iflag = 2;



RETURN:

    if (istack != NULL) free (istack);

    if (wstore != NULL) free (wstore);

    imsl_e1pop ("l_mffct");

    return;

}



/*Translated by FOR_C, v3.4 (P), on 02/11/93 at 09:41:44 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/11/93 at 09:38:49

 *  Options SET: fmt=t s=n

 *-----------------------------------------------------------------------

 *  IMSL Name:  L6FXD/DL6FXD (Single/Double precision version)

 *

 *  Computer:   SUN4/SINGLE

 *

 *  Revised:    March 7, 1990

 *

 *  Purpose:    Add a frontal submatrix into another frontal matrix.

 *              The lower triangular part of a frontal matrix is stored

 *              column by column in contiguous storage.

 *

 *  Usage:      CALL L6FXD (NMAT, AMAT, NSUB, SUBMAT, IXTEND)

 *

 *  Arguments:

 *     NMAT   - Size of frontal matrix.  (Input)

 *     AMAT   - Frontal matrix.  (Input/Output)

 *     NSUB   - Size of frontal matrix.  (Input)

 *     SUBMAT - Contents of frontal matrix.  (Input)

 *     IXTEND - Extension vector of submatrix to frontal matrix.

 *              (Input)

 *

 *  Chapter:    MATH/LIBRARY Linear Systems

 *

 *  Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

 *

 *  Warranty:   IMSL warrants only that IMSL testing has been applied

 *              to this code.  No other warranty, expressed or implied,

 *              is applicable.

 *

 *-----------------------------------------------------------------------

 * */

#ifdef ANSI

static void l_l6fxd (Mint *nmat, Mfloat amat[], Mint *nsub,

                Mfloat submat[], Mint ixtend[])

#else

static void l_l6fxd (nmat, amat, nsub, submat, ixtend)

    Mint        *nmat;

    Mfloat      amat[];

    Mint        *nsub;

    Mfloat      submat[];

    Mint         ixtend[];

#endif

{

    Mint         isub, isub_, ixsub, jmat, jorign, jsub, jsub_, matloc;



/*  for an 'n'-by-'n' lower frontal matrix, function 'idispl' gives the

    displacement so that the (i,j)-th entry in this lower triangular 

    matrix is stored at location 'isispl(n,j) + i' in storage vector */



#ifdef COMPUTER_DECOSF

#define IDISPL(n,j)	((Mint)(((j) - 1)*(2*(n) - (j))/2))

#else

#define IDISPL(n,j)	((long)(((j) - 1)*(2*(n) - (j))/2))

#endif

    if (*nsub <= 0)

	return;

    ixsub = 1;

    for (jsub = 1; jsub <= *nsub; jsub++) {

	jsub_ = jsub - 1;

	jmat = ixtend[jsub_];

	jorign = IDISPL (*nmat, jmat);

	for (isub = jsub; isub <= *nsub; isub++) {

	    isub_ = isub - 1;

	    matloc = jorign + ixtend[isub_];

	    amat[matloc - 1] += submat[ixsub - 1];

	    ixsub += 1;

	}

    }

    return;

#undef	IDISPL

}



/*Translated by FOR_C, v3.4 (P), on 02/11/93 at 09:41:49 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/11/93 at 09:38:51

 *  Options SET: fmt=t s=n

 *-----------------------------------------------------------------------

 *  IMSL Name:  L8FXD/DL8FXD (Single/Double precision version)

 *

 *  Computer:   SUN4/SINGLE

 *

 *  Revised:    March 7, 1990

 *

 *  Purpose:    Perform numeric elimination of frontal matrix.

 *

 *  Usage:      CALL L8FXD (NMAT, AMAT, NELIM, IERROR)

 *

 *  Arguments:

 *     NMAT   - Number of equations.  (Input)

 *     AMAT   - Triangular matrix.  (Input/Output)

 *     NELIM  - Number of eliminations to be performed.  (Input)

 *     IERROR - Error flag set to 1 is negative or zero pivots.

 *              (Output)

 *

 *  Chapter:    MATH/LIBRARY Linear Systems

 *

 *  Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

 *

 *  Warranty:   IMSL warrants only that IMSL testing has been applied

 *              to this code.  No other warranty, expressed or implied,

 *              is applicable.

 *

 *-----------------------------------------------------------------------

 * */

#ifdef ANSI

static void l_l8fxd (Mint *nmat, Mfloat amat[], Mint *nelim,

                Mint *ierror)

#else

static void l_l8fxd (nmat, amat, nelim, ierror)

    Mint        *nmat;

    Mfloat      amat[];

    Mint        *nelim;

    Mint        *ierror;

#endif

{

    Mint         i, i_, ixdiag, ixmat, j, k, k_, kstop, kstrt, nrows;

    Mfloat      aljk, diag;



    *ierror = 1;

    nrows = *nmat;

    ixdiag = 1;

    for (j = 1; j <= *nelim; j++) {

	diag = amat[ixdiag - 1];

	if (diag <= 0.0)

	    goto L_50;

	diag = sqrt (diag);

	amat[ixdiag - 1] = diag;

	nrows -= 1;

	if (nrows <= 0)

	    goto L_40;



/*  divide column j by sqrt of its diagonal element */



	kstrt = ixdiag + 1;

	kstop = ixdiag + nrows;

	for (k = kstrt; k <= kstop; k++) {

	    k_ = k - 1;

	    amat[k_] /= diag;

	}



/*  outer-product modification of remaining matrix by col j  */



	ixdiag = kstop + 1;

	ixmat = ixdiag;

	for (k = kstrt; k <= kstop; k++) {

	    k_ = k - 1;

	    aljk = amat[k_];

	    for (i = k; i <= kstop; i++) {

		i_ = i - 1;

		amat[ixmat - 1] += -amat[i_] * aljk;

		ixmat += 1;

	    }

	}

L_40:

	;

    }

    *ierror = 0;

L_50:

    return;

}



/*Translated by FOR_C, v3.4 (P), on 02/11/93 at 09:41:51 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn s=dvowa str=ln */

/* Structured by FOR_STRUCT, v2.0, on 02/11/93 at 09:38:52

 *  Options SET: fmt=t s=n

 *-----------------------------------------------------------------------

 *  IMSL Name:  L9FXD/DL9FXD (Single/Double precision version)

 *

 *  Computer:   SUN4/SINGLE

 *

 *  Revised:    March 7, 1990

 *

 *  Purpose:    Perform an in-place extension of a frontal matrix.  The

 *              starting entries of the input frontal matrix and the

 *              extended matrix are assumed to be in location 1 of

 *              the storage vector.

 *

 *  Usage:      CALL L9FXD (NMAT, AMAT, NSUB, IXTEND)

 *

 *  Arguments:

 *     NMAT   - Size of the frontal matrix.  (Input)

 *     AMAT   - Frontal matrix.  (Input/Output)

 *     NSUB   - Size of submatrix.  (Input)

 *     IXTEND - Extension vector of submatrix to frontal matrix.

 *              (Input)

 *

 *  Chapter:    MATH/LIBRARY Linear Systems

 *

 *  Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

 *

 *  Warranty:   IMSL warrants only that IMSL testing has been applied

 *              to this code.  No other warranty, expressed or implied,

 *              is applicable.

 *

 *-----------------------------------------------------------------------

 * */

#ifdef ANSI

static void l_l9fxd (Mint *nmat, Mfloat amat[], Mint *nsub,

                Mint ixtend[])

#else

static void l_l9fxd (nmat, amat, nsub, ixtend)

    Mint        *nmat;

    Mfloat      amat[];

    Mint        *nsub;

    Mint         ixtend[];

#endif

{

    Mint         i, isub, ixsub, j, jmat, jorign, jsub, matloc, matsze;

    Mint         loop;



#ifdef COMPUTER_DECOSF

#define IDISPL(n,j)	((Mint)(((j) - 1)*(2*(n) - (j))/2))

#else

#define IDISPL(n,j)	((long)(((j) - 1)*(2*(n) - (j))/2))

#endif

    if (*nsub <= 0 || *nsub == *nmat)

	return;



/*  initialization -- set tail part of storage vector to zero */



    ixsub = *nsub * (*nsub + 1) / 2;

    matsze = *nmat * (*nmat + 1) / 2;

    matloc = ixsub + 1;

/*

	sset( ADR(_l0,matsze - ixsub), ADR(_f0,0.0e0), &amat[matloc - 1],

	 ADR(_l1,1) );

*/

    for (loop = 0; loop < matsze - ixsub; loop++)

	amat[matloc - 1 + loop] = 0.0;



/*  perform backward extension on columns 'nsub' downto 1 */



    jsub = *nsub;

    for (j = 1; j <= *nsub; j++) {

	jmat = ixtend[jsub - 1];

	jorign = IDISPL (*nmat, jmat);

	isub = *nsub;

	for (i = jsub; i <= *nsub; i++) {

	    matloc = jorign + ixtend[isub - 1];

	    if (ixsub == matloc)

		goto L_10;

	    amat[matloc - 1] = amat[ixsub - 1];

	    amat[ixsub - 1] = 0.0;

    L_10:

	    ixsub -= 1;

	    isub -= 1;

	}

	jsub -= 1;

    }

    return;

#undef	IDISPL

}

