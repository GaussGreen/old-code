#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#define MAX_TOKENS  3
#define ARRAY(P)         lv_array[7&lv_parse[P]]
#define IS_MATRIX(P)     ((7&lv_parse[P])==1 || (7&lv_parse[P])==2)
#define IS_VECTOR(P)     ((7&lv_parse[P])==3 || (7&lv_parse[P])==4)
#define IS_TRANSPOSE(P)  (lv_parse[P]&8)
#define IS_CTRANSPOSE(P) (lv_parse[P]&16)
#define IS_TRANS(P)      (lv_parse[P]&24)
#define IS_SET(P)        ARRAY(P).set

static VA_LIST_HACK PROTO (l_mat_mul_rect, (Mchar *string, va_list argptr));
static Mint PROTO (l_parse_string, (Mchar *string));
static void PROTO (l_compute, (Mint ntoken, Mchar *string));
static void PROTO (l_ccopy_conj, (Mint n, Mf_complex cx[], Mint incx, Mf_complex cy[],
	            Mint incy));
    static void PROTO (l_conj, (Mint n, Mf_complex cx[], Mint incx));

    static Mint lv_parse[MAX_TOKENS+1];
    static struct {
	Mint        nrow;
	Mint        ncol;
	Mint        maxcol;
	Mint        set;
	Mf_complex *value;
    }           lv_array[5];

    static Mf_complex *lv_result;
    static Mint lv_result_maxcol;
    static Mchar *lv_name[] = {NULL, "A", "B", "x", "y"};
#ifdef ANSI
    Mf_complex *imsl_c_mat_mul_rect (Mchar *string,...)
#else
    Mf_complex *imsl_c_mat_mul_rect (string, va_alist)
    Mchar      *string;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, string);

    E1PSH ("imsl_c_mat_mul_rect", "imsl_z_mat_mul_rect");
    lv_result = NULL;
    IMSL_CALL (l_mat_mul_rect (string, argptr));
    E1POP ("imsl_c_mat_mul_rect", "imsl_z_mat_mul_rect");
    return lv_result;
}


#ifdef ANSI
static VA_LIST_HACK l_mat_mul_rect (Mchar *string, va_list argptr)
#else
static VA_LIST_HACK l_mat_mul_rect (string, argptr)
    Mchar      *string;
    va_list     argptr;
#endif
{
    Mint        code = 1;
    Mint        arg_number = 1;
    Mint        ntoken;
    Mint        k;
    lv_array[1].set = lv_array[2].set = lv_array[3].set = lv_array[4].set = 0;
    lv_array[1].maxcol = lv_array[2].maxcol = -1;
    lv_array[3].ncol = lv_array[4].ncol = 1;
    lv_array[3].maxcol = lv_array[4].maxcol = 1;
    lv_result_maxcol = -1;

    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_A_MATRIX:
	    arg_number += 3;
	    lv_array[1].nrow = va_arg (argptr, Mint);
	    lv_array[1].ncol = va_arg (argptr, Mint);
	    lv_array[1].value = va_arg (argptr, Mf_complex *);
	    lv_array[1].set = 1;
	    if (!lv_array[1].value) {
		imsl_e1stl (1, "A");
		imsl_e1stl (2, "IMSL_A_MATRIX");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	    }
	    break;
	case IMSL_A_COL_DIM:
	    lv_array[1].maxcol = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_B_MATRIX:
	    arg_number += 3;
	    lv_array[2].nrow = va_arg (argptr, Mint);
	    lv_array[2].ncol = va_arg (argptr, Mint);
	    lv_array[2].value = va_arg (argptr, Mf_complex *);
	    lv_array[2].set = 1;
	    if (!lv_array[2].value) {
		imsl_e1stl (1, "B");
		imsl_e1stl (2, "IMSL_B_MATRIX");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	    }
	    break;
	case IMSL_B_COL_DIM:
	    lv_array[2].maxcol = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_X_VECTOR:
	    arg_number += 2;
	    lv_array[3].nrow = va_arg (argptr, Mint);
	    lv_array[3].value = va_arg (argptr, Mf_complex *);
	    lv_array[3].set = 1;
	    if (!lv_array[3].value) {
		imsl_e1stl (1, "X");
		imsl_e1stl (2, "IMSL_X_VECTOR");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	    }
	    break;
	case IMSL_Y_VECTOR:
	    arg_number += 2;
	    lv_array[4].nrow = va_arg (argptr, Mint);
	    lv_array[4].value = va_arg (argptr, Mf_complex *);
	    lv_array[4].set = 1;
	    if (!lv_array[4].value) {
		imsl_e1stl (1, "Y");
		imsl_e1stl (2, "IMSL_Y_VECTOR");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	    }
	    break;
	case IMSL_RETURN_USER:
	    arg_number++;
	    lv_result = va_arg (argptr, Mf_complex *);
	    if (!lv_result) {
		imsl_e1stl (1, "ans");
		imsl_e1stl (2, "IMSL_RETURN_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	    }
	    break;
	case IMSL_RETURN_COL_DIM:
	    arg_number++;
	    lv_result_maxcol = va_arg (argptr, Mint);
	    break;
	case 0:
	    break;
	default:
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    break;
	}
    }
    if (imsl_n1rty (0))
	goto RETURN;

    for (k = 1; k < 3; k++) {
	if (lv_array[k].set)
	    if (lv_array[k].maxcol == -1) {
		lv_array[k].maxcol = lv_array[k].ncol;
	    }
	    else if (lv_array[k].maxcol < lv_array[k].ncol) {
		/*
		 * The number of columns in the matrix must be less than or
		 * equal to its column dimension while for %(L1) the number
		 * of columns is %(I1) and the column dimension is %(I2).
		 */
		imsl_e1stl (1, (k == 1) ? "A" : "B");
		imsl_e1sti (2, lv_array[k].maxcol);
		imsl_e1sti (1, lv_array[k].ncol);
		imsl_ermes (IMSL_TERMINAL, IMSL_COL_DIM_LESS_COL);
	    }
    }
    if (imsl_n1rty (0))
	goto RETURN;

    ntoken = l_parse_string (string);
    if (ntoken) {
	Mint        col_dim = IS_TRANSPOSE (ntoken - 1) ? ARRAY (ntoken - 1).nrow : ARRAY (ntoken - 1).ncol;
	if (lv_result_maxcol > -1 && (lv_result_maxcol < col_dim)) {
	    imsl_e1sti (1, lv_result_maxcol);
	    imsl_e1sti (2, col_dim);
	    imsl_e1stl (1, string);
	    imsl_ermes (IMSL_TERMINAL, IMSL_BAD_RESULT_COL_DIM);
	}
	else
	    l_compute (ntoken, string);
    }
RETURN:
    return (argptr);
}
#ifdef ANSI
static Mint l_parse_string (Mchar *string)
#else
static Mint l_parse_string (string)
    Mchar      *string;
#endif
{
    Mchar      *p;
    Mchar      *q;
    Mint        ntoken = 0;
    Mint        len;
    Mint        error = 0;
    p = q = string;
    while (q != NULL && !error) {
	if (ntoken > MAX_TOKENS) {
	    error = 1;
	}
	else {
	    q = strchr (p, '*');
	    if (q == NULL)
		len = strlen (p);
	    else
		len = q - p;
	    if (len == 1) {
		switch (*p) {
		case 'A':
		    lv_parse[ntoken++] = 1;
		    break;
		case 'B':
		    lv_parse[ntoken++] = 2;
		    break;
		case 'x':
		    lv_parse[ntoken++] = 3;
		    break;
		case 'y':
		    lv_parse[ntoken++] = 4;
		    break;
		default:
		    error = 1;
		}
	    }
	    else if (len == 8 && strncmp (p, "trans(", 6) == 0) {
		switch (*(p + 6)) {
		case 'A':
		    lv_parse[ntoken++] = 9;
		    break;
		case 'B':
		    lv_parse[ntoken++] = 10;
		    break;
		case 'x':
		    lv_parse[ntoken++] = 11;
		    break;
		case 'y':
		    lv_parse[ntoken++] = 12;
		    break;
		default:
		    error = 1;
		}
	    }
	    else if (len == 9 && strncmp (p, "ctrans(", 7) == 0) {
		switch (*(p + 7)) {
		case 'A':
		    lv_parse[ntoken++] = 17;
		    break;
		case 'B':
		    lv_parse[ntoken++] = 18;
		    break;
		case 'x':
		    lv_parse[ntoken++] = 19;
		    break;
		case 'y':
		    lv_parse[ntoken++] = 20;
		    break;
		default:
		    error = 1;
		}
	    }
	    else {
		error = 1;
	    }
	    p = q + 1;
	}
    }
    if (error) {
	imsl_e1stl (1, string);
	imsl_ermes (IMSL_TERMINAL, IMSL_INVALID_MULT_STRING);
	ntoken = 0;
	goto RETURN;
    }

RETURN:
    return ntoken;
}
#ifdef ANSI
static void l_compute (Mint ntoken, Mchar *string)
#else
static void l_compute (ntoken, string)
    Mint        ntoken;
    Mchar      *string;
#endif
{
    Mint        narow;
    Mint        nacol;
    Mint        nbrow;
    Mint        nbcol;
    Mint        ncrow;
    Mint        nccol;
    Mint        ia;
    Mint        ix;
    Mint        k;
    Mchar      *trans_a;
    static Mf_complex czero = {0.0, 0.0};
    static Mf_complex cone = {1.0, 0.0};
    static Mint ione = 1;
    for (ia = 0; ia < ntoken; ia++) {
	if (!IS_SET (ia)) {
	    /*
	     * The array %(L1) is used but has not been set via an optional
	     * argument.
	     */
	    imsl_e1stl (1, lv_name[7 & lv_parse[ia]]);
	    imsl_ermes (IMSL_TERMINAL, IMSL_MAT_MUL_UNDEFINED);
	    return;
	}
    }
    switch (ntoken) {
    case 1:

	narow = (IS_TRANS (0)) ? ARRAY (0).ncol : ARRAY (0).nrow;
	nacol = (IS_TRANS (0)) ? ARRAY (0).nrow : ARRAY (0).ncol;
	if (lv_result_maxcol == -1)
	    lv_result_maxcol = nacol;

	/* allocate workspace */
	if (lv_result == NULL) {
	    lv_result = (Mf_complex *) imsl_malloc (lv_result_maxcol * narow * sizeof (Mf_complex));
	    if (lv_result == NULL) {
		imsl_ermes(IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
		return;
	    }
	}
	if (IS_VECTOR (0)) {
	    if (IS_TRANSPOSE (0))
		imsl_ccopy (&ARRAY (0).nrow, ARRAY (0).value, &ione, lv_result, &ione);
	    else
		l_ccopy_conj (ARRAY (0).nrow, ARRAY (0).value, ione, lv_result, ione);
	}
	else if (!IS_TRANS (0)) {
	    Mf_complex *a = ARRAY (0).value;
	    Mf_complex *r = lv_result;
	    Mint        k;
	    for (k = 0; k < ARRAY (0).nrow; k++) {
		imsl_ccopy (&ARRAY (0).ncol, a, &ione, r, &ione);
		a += ARRAY (0).maxcol;
		r += lv_result_maxcol;
	    }
	}
	else if (IS_TRANSPOSE (0)) {	/* return transpose */
	    Mf_complex *a = ARRAY (0).value;
	    Mf_complex *r = lv_result;
	    for (k = 0; k < ARRAY (0).nrow; k++) {
		imsl_ccopy (&ARRAY (0).ncol, a, &ione, r, &lv_result_maxcol);
		a += ARRAY (0).maxcol;
		r++;
	    }
	}
	else {			/* return conjugate transpose */
	    Mf_complex *a = ARRAY (0).value;
	    Mf_complex *r = lv_result;
	    for (k = 0; k < ARRAY (0).nrow; k++) {
		l_ccopy_conj (ARRAY (0).ncol, a, ione, r, lv_result_maxcol);
		a += ARRAY (0).maxcol;
		r++;
	    }
	}
	for (k = nacol; k < lv_result_maxcol; k++) {
	    imsl_cset (&narow, &czero, (lv_result + k), &lv_result_maxcol);
	}
	break;
    case 2:
	nacol = IS_TRANS (0) ? ARRAY (0).nrow : ARRAY (0).ncol;
	nbrow = IS_TRANS (1) ? ARRAY (1).ncol : ARRAY (1).nrow;
	narow = IS_TRANS (0) ? ARRAY (0).ncol : ARRAY (0).nrow;
	nbcol = IS_TRANS (1) ? ARRAY (1).nrow : ARRAY (1).ncol;
	/* check that the dimensions match */
	if (nacol != nbrow) {
	    /*
	     * Cannot multiply a %(I1) by %(I2) matrix and a %(I3) by %(I4)
	     * matrix.
	     */
	    imsl_e1sti (1, narow);
	    imsl_e1sti (2, nacol);
	    imsl_e1sti (3, nbrow);
	    imsl_e1sti (4, nbcol);
	    imsl_ermes (IMSL_TERMINAL, IMSL_MATMUL_DIM_MISMATCH_2);
	    goto RETURN;
	}
	if (lv_result_maxcol == -1)
	    lv_result_maxcol = nbcol;
	/* allocate workspace */
	if (lv_result == NULL) {
	    lv_result = (Mf_complex *)
		imsl_malloc (narow * lv_result_maxcol * sizeof (Mf_complex));
	    if (lv_result == NULL) {
		imsl_ermes(IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
		return;
	    }
	}
	for (k = nbcol; k < lv_result_maxcol; k++) {
	    imsl_cset (&narow, &czero, (lv_result + k), &lv_result_maxcol);
	}
	/* matrix * vector   or  vector * matrix  */
	if (IS_MATRIX (0) && IS_VECTOR (1) || IS_MATRIX (1) && IS_VECTOR (0)) {
	    Mint        x_conj;
	    Mint        y_conj = 0;
	    Mint        incr = 1; 
	    if (IS_MATRIX (0) && IS_VECTOR (1)) {
		ia = 0;
		ix = 1;
		x_conj = IS_CTRANSPOSE (ix);
		trans_a = (IS_TRANS (ia)) ? "N" : "T";
		if (IS_CTRANSPOSE (ia)) {
		    x_conj = !x_conj;
		    y_conj = 1;
		}
	    }
	    else {
		ia = 1;
		ix = 0;
		x_conj = IS_CTRANSPOSE (ix);
		trans_a = (IS_TRANS (ia)) ? "T" : "N";
		if (IS_CTRANSPOSE (ia))
		    trans_a = "C";
	    }
	    if (narow > 1){
	      incr = lv_result_maxcol;
	    }
	    if (x_conj)
		l_conj (ARRAY (ix).nrow, ARRAY (ix).value, 1);
	    imsl_cgemv (trans_a, ione, &ARRAY (ia).ncol, &ARRAY (ia).nrow, &cone,
		ARRAY (ia).value, &ARRAY (ia).maxcol, ARRAY (ix).value, &ione,
		&czero, lv_result, &incr);
	    if (x_conj)
		l_conj (ARRAY (ix).nrow, ARRAY (ix).value, 1);
	    if (y_conj)
		l_conj (narow, lv_result, lv_result_maxcol);
	    break;
	}
	/* trans(vector) * vector */
	if (IS_VECTOR (0) && IS_VECTOR (1) &&
	    IS_TRANSPOSE (0) && !IS_TRANS (1)) {
	    *lv_result = imsl_cdotu (&ARRAY (0).nrow, ARRAY (0).value, &ione,
		ARRAY (1).value, &ione);
	    break;
	}
	if (IS_VECTOR (0) && IS_VECTOR (1) &&
	    IS_CTRANSPOSE (0) && !IS_TRANS (1)) {
	    *lv_result = imsl_cdotc (&ARRAY (0).nrow, ARRAY (0).value, &ione,
		ARRAY (1).value, &ione);
	    break;
	}
	/* general case */
	{
	    Mint        cdb = (IS_TRANS (1)) ? 1 : ARRAY (1).maxcol;
	    Mint        incb = (IS_TRANS (1)) ? ARRAY (1).maxcol : 1;
	    Mint        cda;
	    Mf_complex *a_hold;
	    Mf_complex *b = ARRAY (1).value;
	    Mf_complex *b1;
	    Mf_complex *r = lv_result;
	    Mint        m, n, k, inc_b1, nb;
	    if (IS_CTRANSPOSE (0)) {
		Mf_complex *a, *r;

		trans_a = "C";
		m = ARRAY (0).nrow;
		n = ARRAY (0).ncol;
		cda = m;
		nb = m;
				/* Store ctrans(ARRAY(0)) in copy */
		a = a_hold = ARRAY(0).value;
		r = ARRAY(0).value =
		    (Mf_complex *)imsl_malloc (m * n * sizeof (Mf_complex));
		if (a_hold == NULL) {
		    imsl_ermes(IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
		    return;
		}
		for (k = 0; k < m; k++) {
		    imsl_ccopy(&n, a, &ione, r, &m);
		    a += ARRAY(0).maxcol;
		    r++;
		}
	    }
	    else {
		trans_a = (IS_TRANSPOSE (0)) ? "N" : "T";
		m = ARRAY (0).ncol;
		n = ARRAY (0).nrow;
		cda = ARRAY (0).maxcol;
		nb = n;
	    }
	    if (IS_CTRANSPOSE (1)) {
		b1 = (Mf_complex *)imsl_malloc (nb * sizeof (Mf_complex));
		if (b1 == NULL) {
		    imsl_ermes(IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
		    return;
		}
		inc_b1 = 1;
	    }
	    else {
		inc_b1 = cdb;
	    }
	    for (k = 0; k < nbcol; k++) {
		if (IS_CTRANSPOSE (1)) {
		    l_ccopy_conj (nb, b, cdb, b1, inc_b1);
		}
		else {
		    b1 = b;
		}
		imsl_cgemv (trans_a, ione, &m, &n,
		    &cone, ARRAY (0).value, &cda,
		    b1, &inc_b1, &czero, r, &lv_result_maxcol);
		b += incb;
		r++;
	    }
	    if (strcmp (trans_a, "C") == 0) {
		imsl_free(ARRAY(0).value);
		ARRAY(0).value = a_hold;
	    }
	    if (IS_CTRANSPOSE (1))
		imsl_free (b1);
	    break;
	}
    case 3:
	narow = IS_TRANS (0) ? ARRAY (0).ncol : ARRAY (0).nrow;
	nacol = IS_TRANS (0) ? ARRAY (0).nrow : ARRAY (0).ncol;
	nbrow = IS_TRANS (1) ? ARRAY (1).ncol : ARRAY (1).nrow;
	nbcol = IS_TRANS (1) ? ARRAY (1).nrow : ARRAY (1).ncol;
	ncrow = IS_TRANS (2) ? ARRAY (2).ncol : ARRAY (2).nrow;
	nccol = IS_TRANS (2) ? ARRAY (2).nrow : ARRAY (2).ncol;
	if (nacol != nbrow || nbcol != ncrow) {
	    /*
	     * Cannot multiply a %(I1) by %(I2) matrix times a %(I3) by %(I4)
	     * matrix times a %(I5) by %(I6) matrix.
	     */
	    imsl_e1sti (1, narow);
	    imsl_e1sti (2, nacol);
	    imsl_e1sti (3, nbrow);
	    imsl_e1sti (4, nbcol);
	    imsl_e1sti (5, ncrow);
	    imsl_e1sti (6, nccol);
	    imsl_ermes (IMSL_TERMINAL, IMSL_MATMUL_DIM_MISMATCH_3);
	    goto RETURN;
	}
	/* trans(vector) * matrix * vector   (bilinear form) */
	if (IS_VECTOR (0) && IS_MATRIX (1) && IS_VECTOR (2) &&
	    IS_TRANSPOSE (0) && !IS_TRANSPOSE (1) && !IS_TRANSPOSE (2)) {
	    Mf_complex *row;
	    Mf_complex  sum;
	    Mf_complex  term;
	    sum.re = sum.im = F_ZERO;
	    for (row = ARRAY (1).value, k = 0; k < ARRAY (0).nrow; k++) {
		term = imsl_cdotu (&ARRAY (1).ncol, row, &ione,
		    ARRAY (2).value, &ione);
		term = imsl_c_mul (term, ARRAY (0).value[k]);
		sum.re += term.re;
		sum.im += term.im;
		row += ARRAY (1).maxcol;
	    }
	    if (lv_result_maxcol == -1) {
		lv_result_maxcol = 1;
	    }
	    if (lv_result == NULL) {
		lv_result = (Mf_complex *) imsl_malloc (lv_result_maxcol * sizeof (Mf_complex));
		if (lv_result == NULL) {
		    imsl_ermes(IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
		    return;
		}
	    }
	    *lv_result = sum;
	    if (lv_result_maxcol > 1) {
		ia = lv_result_maxcol - 1;
		imsl_cset (&ia, &czero, (lv_result + 1), &ione);
	    }
	    break;
	} else {
				/* compute A*B*C as temp=B*C, A*temp */
	    Mf_complex	*temp;
	    Mint	temp_parse = (lv_parse[0]&7 != 2) ? 2 : 1;
	    Mint	parse0 = lv_parse[0];
	    Mint	result_maxcol = lv_result_maxcol;

	    lv_parse[0] = lv_parse[1];
	    lv_parse[1] = lv_parse[2];
	    lv_result_maxcol = nccol;
	    l_compute(2, string);
	    lv_parse[0] = parse0;
	    lv_parse[1] = temp_parse;
	    lv_result_maxcol = result_maxcol;
				/* set A or B to temp */
	    lv_array[temp_parse].nrow   = nbrow;
	    lv_array[temp_parse].ncol   = nccol;
	    lv_array[temp_parse].maxcol = nccol;
	    lv_array[temp_parse].value = temp = lv_result;
	    lv_array[temp_parse].set = 1;
	    l_compute(2, string);
	    imsl_free(temp);
	    break;
	}
    }
RETURN:
    return;
}


#ifdef ANSI
static void l_ccopy_conj (Mint n, Mf_complex cx[], Mint incx, Mf_complex cy[],
                Mint incy)
#else
static void l_ccopy_conj (n, cx, incx, cy, incy)
    Mint        n;
    Mf_complex  cx[];
    Mint        incx;
    Mf_complex  cy[];
    Mint        incy;
#endif
{
    Mint        i, ix, iy;
    if (n <= 0)
	return;
    if (incx != 1 || incy != 1) {
	ix = 1;
	iy = 1;
	if (incx < 0)
	    ix = (-n + 1) * incx + 1;
	if (incy < 0)
	    iy = (-n + 1) * incy + 1;
	for (i = 0; i < n; i++) {
	    cy[iy - 1].re = cx[ix - 1].re;
	    cy[iy - 1].im = -cx[ix - 1].im;
	    ix += incx;
	    iy += incy;
	}
    }
    else {
	for (i = 0; i < n; i++) {
	    cy[i].re = cx[i].re;
	    cy[i].im = -cx[i].im;
	}
    }
}



#ifdef ANSI
static void l_conj (Mint n, Mf_complex cx[], Mint incx)
#else
static void l_conj (n, cx, incx)
    Mint        n;
    Mf_complex  cx[];
    Mint        incx;
#endif
{
    Mint        k;
    Mint        j;
    for (j = 0, k = 0; k < n; k++) {
	cx[j].im = -cx[j].im;
	j += incx;
    }
}
