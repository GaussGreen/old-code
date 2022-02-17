#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_sort (Mint, Mfloat*, va_list);
static void l_svrgp (Mint, Mfloat[], Mfloat[], Mint[]);
static void l_svrbp (Mint*, Mfloat[], Mfloat[], Mint[]);
static void l_svrgn (Mint, Mfloat[], Mfloat[]);
static void l_svrbn (Mint*, Mfloat[], Mfloat[]);
#else
static VA_LIST_HACK l_sort ();
static void l_svrgp ();
static void l_svrbp ();
static void l_svrgn ();
static void l_svrbn ();
#endif

static Mfloat *lv_sort = NULL;
#ifdef ANSI
Mfloat     *imsl_f_sort (Mint n, Mfloat *x, ...)
#else
Mfloat     *imsl_f_sort (n, x, va_alist)
    Mint        n;
    Mfloat     *x;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, x);
    E1PSH ("imsl_f_sort", "imsl_d_sort");
    lv_sort = NULL;
    IMSL_CALL (l_sort (n, x, argptr));
    va_end (argptr);
    E1POP ("imsl_f_sort", "imsl_d_sort");
    return lv_sort;
}


#ifdef ANSI
static VA_LIST_HACK l_sort (Mint n, Mfloat *x, va_list argptr)
#else
static VA_LIST_HACK l_sort (n, x, argptr)
    Mint        n;
    Mfloat     *x;
    va_list     argptr;
#endif
{
    Mint        code;
    Mint        arg_number = 2;
    Mint	algebraic 	= 1;
    Mint	absolute	= 0;
    Mint	permutation	= 0;
    Mint	user_supplied_perm_space 	= 0;
    Mint	user_supplied_sort_space	= 0;
    Mint	**iperm		= NULL;
    Mint	*user_iperm	= NULL;
    Mint	i;

    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_ABSOLUTE:
	    algebraic = 0;
	    absolute = 1;
	    break;
	case IMSL_PERMUTATION:
	    iperm = va_arg(argptr, Mint **);
	    permutation = 1;
	    arg_number++;
	    break;
	case IMSL_PERMUTATION_USER:
	    user_iperm = va_arg(argptr, Mint *);
	    user_supplied_perm_space = 1;
	    permutation = 1;
	    arg_number++;
	    break;
	case IMSL_RETURN_USER:
	    lv_sort = va_arg(argptr, Mfloat *);
	    user_supplied_sort_space = 1;
	    arg_number++;
	case 0:
	    break;
	default:
	    /* Argument number %(I2) is an unknown */
	    /* optional argument %(I1). */
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    break;
	}
    }
 
    if (n < 1) {
        imsl_e1sti (1, n);
        imsl_ermes (IMSL_TERMINAL, IMSL_LARGER_N_VALUE_NEEDED);
    }


    if (imsl_n1rty (0))
	goto FREE_SPACE;

    if (!user_supplied_sort_space) 
	lv_sort = (Mfloat*) imsl_malloc(n*sizeof(*lv_sort));

    if (permutation && !user_supplied_perm_space) 
	*iperm = (Mint*) imsl_malloc(n*sizeof(**iperm));

    if (lv_sort == NULL || (permutation && !user_supplied_perm_space &&
		*iperm == NULL) ) {
	imsl_e1stl(1,"n");
	imsl_e1sti(1, n);
	imsl_ermes(IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }
    if (permutation) {
	if (!user_supplied_perm_space) {
		for (i=0; i<n; i++)
			*(*iperm+i) = i;
    		if (algebraic)
			l_svrgp(n, x, lv_sort, *iperm);
		if (absolute)
			l_svrbp(&n, x, lv_sort, *iperm);
	} else {
		for (i=0; i<n; i++)
			*(user_iperm+i) = i;
    		if (algebraic)
			l_svrgp(n, x, lv_sort, user_iperm);
		if (absolute)
			l_svrbp(&n, x, lv_sort, user_iperm);
	}
    }
		
    if (!permutation) {
	if (algebraic)
		l_svrgn(n, x, lv_sort);
	if (absolute)
		l_svrbn(&n, x, lv_sort);
    }


	
FREE_SPACE:
    if (imsl_n1rty(0) > 3) {
	if (lv_sort != NULL) imsl_free(lv_sort);
	lv_sort = NULL;
    }
RETURN:
    return (argptr);
}




/*Translated by FOR_C++, v0.1, on 08/08/91 at 15:40:06 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/08/91 at 15:40:03
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SVRGP/DSVRGP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 30, 1991

    Purpose:    Sort a real array by algebraically increasing value and
                return the permutation that rearranges the array.

    Usage:      CALL SVRGP (N, RA, RB, IPERM)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       RA     - Vector of length N containing the array to be sorted.
                (Input)
       RB     - Vector of length N containing the sorted array.
                (Output)
                If RA is not needed, RA and RB can share the same
                storage locations.
       IPERM  - Vector of length N.  (Input/Output)
                On input, IPERM should be initialized to the values 1,
                2, ..., N.  On output, IPERM contains a record of
                permutations made on the vector RA.

    Remark:
       For wider applicability, integers (1,2,...,N) that are to be
       associated with RA(I) for I = 1, 2, ..., N may be entered into
       IPERM(I) in any order.  Note that these integers must be unique.

    Keyword:    Utilities

    GAMS:       N6a1b; N6a2b

    Chapters:   MATH/LIBRARY Utilities
                STAT/LIBRARY Utilities

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_svrgp (Mint n, Mfloat ra[], Mfloat rb[], Mint iperm[])
#else
static void l_svrgp (n, ra, rb, iperm)
    Mint         n;
    Mfloat       ra[], rb[];
    Mint         iperm[];
#endif
{
    Mint         i, ij, il[21], it, itt, iu[21], j, k, l, m;
    Mfloat       r, t, tt;


    imsl_e1psh ("SVRGP ");
    scopy (n, ra, 1, rb, 1);

    m = 1;
    i = 1;
    j = n;
    r = .375;
L_10:
    if (i == j)
	goto L_70;
    if (r <= .5898437) {
	r += 3.90625e-2;
    }
    else {
	r -= .21875;
    }
L_20:
    k = i;
    /*
     * SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION T
     */
    ij = i + (j - i) * r;
    t = rb[ij - 1];
    it = iperm[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
     */
    if (rb[i - 1] > t) {
	rb[ij - 1] = rb[i - 1];
	rb[i - 1] = t;
	t = rb[ij - 1];
	iperm[ij - 1] = iperm[i - 1];
	iperm[i - 1] = it;
	it = iperm[ij - 1];
    }
    l = j;
    /*
     * IF LAST ELEMENT OF ARRAY IS LESS THAN T, INTERCHANGE WITH T
     */
    if (rb[j - 1] >= t)
	goto L_40;
    rb[ij - 1] = rb[j - 1];
    rb[j - 1] = t;
    t = rb[ij - 1];
    iperm[ij - 1] = iperm[j - 1];
    iperm[j - 1] = it;
    it = iperm[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
     */
    if (rb[i - 1] <= t)
	goto L_40;
    rb[ij - 1] = rb[i - 1];
    rb[i - 1] = t;
    t = rb[ij - 1];
    iperm[ij - 1] = iperm[i - 1];
    iperm[i - 1] = it;
    it = iperm[ij - 1];
    goto L_40;
L_30:
    if (rb[l - 1] == rb[k - 1])
	goto L_40;
    tt = rb[l - 1];
    rb[l - 1] = rb[k - 1];
    rb[k - 1] = tt;
    itt = iperm[l - 1];
    iperm[l - 1] = iperm[k - 1];
    iperm[k - 1] = itt;
    /*
     * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER THAN
     * T
     */
L_40:
    l -= 1;
    if (rb[l - 1] > t)
	goto L_40;
    /*
     * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER THAN T
     */
L_50:
    k += 1;
    if (rb[k - 1] < t)
	goto L_50;
    /* INTERCHANGE THESE ELEMENTS */
    if (k <= l)
	goto L_30;
    /*
     * SAVE UPPER AND LOWER SUBSCRIPTS OF THE ARRAY YET TO BE SORTED
     */
    if (l - i <= j - k)
	goto L_60;
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    m += 1;
    goto L_80;
L_60:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    m += 1;
    goto L_80;
    /*
     * BEGIN AGAIN ON ANOTHER PORTION OF THE UNSORTED ARRAY
     */
L_70:
    m -= 1;
    if (m == 0)
	goto L_9000;
    i = il[m - 1];
    j = iu[m - 1];
L_80:
    if (j - i >= 11)
	goto L_20;
    if (i == 1)
	goto L_10;
    i -= 1;
L_90:
    i += 1;
    if (i == j)
	goto L_70;
    t = rb[i];
    it = iperm[i];
    if (rb[i - 1] <= t)
	goto L_90;
    k = i;
L_100:
    rb[k] = rb[k - 1];
    iperm[k] = iperm[k - 1];
    k -= 1;
    if (t < rb[k - 1])
	goto L_100;
    rb[k] = t;
    iperm[k] = it;
    goto L_90;

L_9000:
    imsl_e1pop ("SVRGP ");
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  SVRBP/DSVRBP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 31, 1991

    Purpose:    Sort a real array by nondecreasing absolute value and
                return the permutation that rearranges the array.

    Usage:      CALL SVRBP (N, RA, RB, IPERM)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       RA     - Vector of length N containing the array to be sorted.
                (Input)
       RB     - Vector of length N containing the sorted array.
                (Output)
                If RA is not needed, RA and RB can share the same
                storage locations.
       IPERM  - Vector of length N.  (Input/Output)
                On input, IPERM should be initialized to the values 1,
                2, ..., N.  On output, IPERM contains a record of
                permutations made on the vector IA.

    Remark:
       For wider applicability, integers (1,2,...,N) that are to be
       associated with RA(I) for I = 1, 2, ..., N may be entered into
       IPERM(I) in any order.  Note that these integers must be unique.

    Keyword:    Utilities

    GAMS:       N6a1b; N6a2b

    Chapter:    MATH/LIBRARY Utilities

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_svrbp (Mint *n, Mfloat ra[], Mfloat rb[], Mint iperm[])
#else
static void l_svrbp (n, ra, rb, iperm)
    Mint        *n;
    Mfloat       ra[], rb[];
    Mint         iperm[];
#endif
{
    Mint         i, ij, ik, il[21], it, itt, iu[21], j, k, l, m;
    Mfloat       r, t, tt;


    imsl_e1psh ("SVRBP ");
    scopy (*n, ra, 1, rb, 1);
    /* FIND ABSOLUTE VALUES OF ARRAY RA */
    for (i = 1; i <= *n; i++) {
	if (rb[i - 1] < 0.0) {
	    iperm[i - 1] = -iperm[i - 1];
	    rb[i - 1] = -rb[i - 1];
	}
    }

    m = 1;
    i = 1;
    j = *n;
    r = .375;
L_20:
    if (i == j)
	goto L_90;

    if (r <= .5898437) {
	r += 3.90625e-2;
    }
    else {
	r -= .21875;
    }

L_40:
    k = i;
    /*
     * SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION T
     */
    ij = i + (j - i) * r;
    t = rb[ij - 1];
    it = iperm[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
     */
    if (rb[i - 1] > t) {
	rb[ij - 1] = rb[i - 1];
	rb[i - 1] = t;
	t = rb[ij - 1];
	iperm[ij - 1] = iperm[i - 1];
	iperm[i - 1] = it;
	it = iperm[ij - 1];
    }
    l = j;
    /*
     * IF LAST ELEMENT OF ARRAY IS LESS THAN T, INTERCHANGE WITH T
     */
    if (rb[j - 1] >= t)
	goto L_60;
    rb[ij - 1] = rb[j - 1];
    rb[j - 1] = t;
    t = rb[ij - 1];
    iperm[ij - 1] = iperm[j - 1];
    iperm[j - 1] = it;
    it = iperm[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
     */
    if (rb[i - 1] <= t)
	goto L_60;
    rb[ij - 1] = rb[i - 1];
    rb[i - 1] = t;
    t = rb[ij - 1];
    iperm[ij - 1] = iperm[i - 1];
    iperm[i - 1] = it;
    it = iperm[ij - 1];
    goto L_60;
L_50:
    if (rb[l - 1] == rb[k - 1])
	goto L_60;
    tt = rb[l - 1];
    rb[l - 1] = rb[k - 1];
    rb[k - 1] = tt;
    itt = iperm[l - 1];
    iperm[l - 1] = iperm[k - 1];
    iperm[k - 1] = itt;
    /*
     * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER THAN
     * T
     */
L_60:
    l -= 1;
    if (rb[l - 1] > t)
	goto L_60;
    /*
     * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER THAN T
     */
L_70:
    k += 1;
    if (rb[k - 1] < t)
	goto L_70;
    /* INTERCHANGE THESE ELEMENTS */
    if (k <= l)
	goto L_50;
    /*
     * SAVE UPPER AND LOWER SUBSCRIPTS OF THE ARRAY YET TO BE SORTED
     */
    if (l - i <= j - k)
	goto L_80;
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    m += 1;
    goto L_100;
L_80:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    m += 1;
    goto L_100;
    /*
     * BEGIN AGAIN ON ANOTHER PORTION OF THE UNSORTED ARRAY
     */
L_90:
    m -= 1;
    if (m == 0)
	goto L_130;
    i = il[m - 1];
    j = iu[m - 1];
L_100:
    if (j - i >= 11)
	goto L_40;
    if (i == 1)
	goto L_20;
    i -= 1;
L_110:
    i += 1;
    if (i == j)
	goto L_90;
    t = rb[i];
    it = iperm[i];
    if (rb[i - 1] <= t)
	goto L_110;
    k = i;
L_120:
    rb[k] = rb[k - 1];
    iperm[k] = iperm[k - 1];
    k -= 1;
    if (t < rb[k - 1])
	goto L_120;
    rb[k] = t;
    iperm[k] = it;
    goto L_110;

L_130:
    for (ik = 1; ik <= *n; ik++) {
	if (iperm[ik - 1] < 0) {
	    iperm[ik - 1] = -iperm[ik - 1];
	    rb[ik - 1] = -rb[ik - 1];
	}
    }

L_9000:
    imsl_e1pop ("SVRBP ");
    return;
}				/* end of function */




/*Translated by FOR_C++, v0.1, on 08/09/91 at 11:26:30 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/91 at 11:26:26
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SVRGN/DSVRGN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 19, 1991

    Purpose:    Sort a real array by algebraically increasing value.

    Usage:      CALL SVRGN (N, RA, RB)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       RA     - Vector of length N containing the array to be sorted.
                (Input)
       RB     - Vector of length N containing the sorted array.
                (Output)
                If RA is not needed, RA and RB can share the same
                storage locations.

    Keyword:    Utilities

    GAMS:       N6a2b

    Chapters:   MATH/LIBRARY Utilities
                STAT/LIBRARY Utilities

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_svrgn (Mint n, Mfloat ra[], Mfloat rb[])
#else
static void l_svrgn (n, ra, rb)
    Mint         n;
    Mfloat       ra[], rb[];
#endif
{
    Mint         i, ij, il[21], iu[21], j, k, l, m;
    Mfloat       r, t, tt;


    imsl_e1psh ("SVRGN ");

    scopy (n, ra, 1, rb, 1);
    /* Check input argument */

    m = 1;
    i = 1;
    j = n;
    r = .375;
L_10:
    if (i == j)
	goto L_80;
L_20:
    if (r <= .5898437) {
	r += 3.90625e-2;
    }
    else {
	r -= .21875;
    }
L_30:
    k = i;
    /*
     * SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION T
     */
    ij = i + (j - i) * r;
    t = rb[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
     */
    if (rb[i - 1] > t) {
	rb[ij - 1] = rb[i - 1];
	rb[i - 1] = t;
	t = rb[ij - 1];
    }
    l = j;
    /*
     * IF LAST ELEMENT OF ARRAY IS LESS THAN T, INTERCHANGE WITH T
     */
    if (rb[j - 1] >= t)
	goto L_50;
    rb[ij - 1] = rb[j - 1];
    rb[j - 1] = t;
    t = rb[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
     */
    if (rb[i - 1] <= t)
	goto L_50;
    rb[ij - 1] = rb[i - 1];
    rb[i - 1] = t;
    t = rb[ij - 1];
    goto L_50;
L_40:
    if (rb[l - 1] == rb[k - 1])
	goto L_50;
    tt = rb[l - 1];
    rb[l - 1] = rb[k - 1];
    rb[k - 1] = tt;
    /*
     * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER THAN
     * T
     */
L_50:
    l -= 1;
    if (rb[l - 1] > t)
	goto L_50;
    /*
     * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER THAN T
     */
L_60:
    k += 1;
    if (rb[k - 1] < t)
	goto L_60;
    /* INTERCHANGE THESE ELEMENTS */
    if (k <= l)
	goto L_40;
    /*
     * SAVE UPPER AND LOWER SUBSCRIPTS OF THE ARRAY YET TO BE SORTED
     */
    if (l - i <= j - k)
	goto L_70;
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    m += 1;
    goto L_90;
L_70:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    m += 1;
    goto L_90;
    /*
     * BEGIN AGAIN ON ANOTHER PORTION OF THE UNSORTED ARRAY
     */
L_80:
    m -= 1;
    if (m == 0)
	goto L_9000;
    i = il[m - 1];
    j = iu[m - 1];
L_90:
    if (j - i >= 11)
	goto L_30;
    if (i == 1)
	goto L_10;
    i -= 1;
L_100:
    i += 1;
    if (i == j)
	goto L_80;
    t = rb[i];
    if (rb[i - 1] <= t)
	goto L_100;
    k = i;
L_110:
    rb[k] = rb[k - 1];
    k -= 1;
    if (t < rb[k - 1])
	goto L_110;
    rb[k] = t;
    goto L_100;

L_9000:
    imsl_e1pop ("SVRGN ");

    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  SVRBN/DSVRBN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 19, 1991

    Purpose:    Sort a real array by nondecreasing absolute value.

    Usage:      CALL SVRBN (N, RA, RB)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       RA     - Vector of length N containing the array to be sorted.
                (Input)
       RB     - Vector of length N containing the sorted array.
                (Output)
                If RA is not needed, RA and RB can share the same
                storage locations.

    Keyword:    Utilities

    GAMS:       N6a2b

    Chapters:   MATH/LIBRARY Utilities

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_svrbn (Mint *n, Mfloat ra[], Mfloat rb[])
#else
static void l_svrbn (n, ra, rb)
    Mint        *n;
    Mfloat       ra[], rb[];
#endif
{
    Mint         i, ij, il[21], iu[21], j, k, l, m;
    Mfloat       r, t, tt;


    imsl_e1psh ("SVRBN ");
    scopy (*n, ra, 1, rb, 1);

    m = 1;
    i = 1;
    j = *n;
    r = .375;
L_10:
    if (i == j)
	goto L_90;
L_20:
    if (r <= .5898437) {
	r += 3.90625e-2;
    }
    else {
	r -= .21875;
    }
L_30:
    k = i;
    /*
     * SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION T
     */
    ij = i + (j - i) * r;
    t = rb[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
     */
    if (fabs (rb[i - 1]) <= fabs (t))
	goto L_40;
    rb[ij - 1] = rb[i - 1];
    rb[i - 1] = t;
    t = rb[ij - 1];
L_40:
    l = j;
    /*
     * IF LAST ELEMENT OF ARRAY IS LESS THAN T, INTERCHANGE WITH T
     */
    if (fabs (rb[j - 1]) >= fabs (t))
	goto L_60;
    rb[ij - 1] = rb[j - 1];
    rb[j - 1] = t;
    t = rb[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
     */
    if (fabs (rb[i - 1]) <= fabs (t))
	goto L_60;
    rb[ij - 1] = rb[i - 1];
    rb[i - 1] = t;
    t = rb[ij - 1];
    goto L_60;
L_50:
    if (fabs (rb[l - 1]) == fabs (rb[k - 1]))
	goto L_60;
    tt = rb[l - 1];
    rb[l - 1] = rb[k - 1];
    rb[k - 1] = tt;
    /*
     * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER THAN
     * T
     */
L_60:
    l -= 1;
    if (fabs (rb[l - 1]) > fabs (t))
	goto L_60;
    /*
     * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER THAN T
     */
L_70:
    k += 1;
    if (fabs (rb[k - 1]) < fabs (t))
	goto L_70;
    /* INTERCHANGE THESE ELEMENTS */
    if (k <= l)
	goto L_50;
    /*
     * SAVE UPPER AND LOWER SUBSCRIPTS OF THE ARRAY YET TO BE SORTED
     */
    if (l - i <= j - k)
	goto L_80;
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    m += 1;
    goto L_100;
L_80:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    m += 1;
    goto L_100;
    /*
     * BEGIN AGAIN ON ANOTHER PORTION OF THE UNSORTED ARRAY
     */
L_90:
    m -= 1;
    if (m == 0)
	goto L_9000;
    i = il[m - 1];
    j = iu[m - 1];
L_100:
    if (j - i >= 11)
	goto L_30;
    if (i == 1)
	goto L_10;
    i -= 1;
L_110:
    i += 1;
    if (i == j)
	goto L_90;
    t = rb[i];
    if (fabs (rb[i - 1]) <= fabs (t))
	goto L_110;
    k = i;
L_120:
    rb[k] = rb[k - 1];
    k -= 1;
    if (fabs (t) < fabs (rb[k - 1]))
	goto L_120;
    rb[k] = t;
    goto L_110;

L_9000:
    imsl_e1pop ("SVRBN ");
    return;
}				/* end of function */
