#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_sort (Mint, Mint*, va_list);
static void l_svigp (Mint*, Mint[], Mint[], Mint[]);
static void l_svibp (Mint*, Mint[], Mint[], Mint[]);
static void l_svign (Mint, Mint[], Mint[]);
static void l_svibn (Mint*, Mint[], Mint[]);
#else
static VA_LIST_HACK l_sort ();
static void l_svigp ();
static void l_svibp ();
static void l_svign ();
static void l_svibn ();
#endif

static Mint *lv_sort = NULL;
#ifdef ANSI
Mint     *imsl_i_sort (Mint n, Mint *x, ...)
#else
Mint     *imsl_i_sort (n, x, va_alist)
    Mint        n;
    Mint     *x;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, x);
    E1PSH ("imsl_i_sort", "imsl_i_sort");
    lv_sort = NULL;
    IMSL_CALL (l_sort (n, x, argptr));
    va_end (argptr);
    E1POP ("imsl_i_sort", "imsl_i_sort");
    return lv_sort;
}


#ifdef ANSI
static VA_LIST_HACK l_sort (Mint n, Mint *x, va_list argptr)
#else
static VA_LIST_HACK l_sort (n, x, argptr)
    Mint        n;
    Mint     *x;
    va_list     argptr;
#endif
{
    Mint        code;
    Mint        arg_number = 2;
    Mint        algebraic       = 1;
    Mint        absolute        = 0;
    Mint        permutation     = 0;
    Mint        user_supplied_perm_space        = 0;
    Mint	user_supplied_sort_space	= 0;
    Mint        **iperm         = NULL;
    Mint        *user_iperm     = NULL;
    Mint        i;
 
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
            lv_sort = va_arg(argptr, Mint *);
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
        lv_sort = (Mint*) imsl_malloc(n*sizeof(*lv_sort));

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
                        l_svigp(&n, x, lv_sort, *iperm);
                if (absolute)
                        l_svibp(&n, x, lv_sort, *iperm);
        } else {
                for (i=0; i<n; i++)
                        *(user_iperm+i) = i;
                if (algebraic)
                        l_svigp(&n, x, lv_sort, user_iperm);
                if (absolute)
                        l_svibp(&n, x, lv_sort, user_iperm);
        }
    }
                
    if (!permutation) {
        if (algebraic)
                l_svign(n, x, lv_sort);
        if (absolute)
                l_svibn(&n, x, lv_sort);
    }
    

FREE_SPACE:
    if (imsl_n1rty(0) > 3) {
        if (lv_sort != NULL) imsl_free(lv_sort);
        lv_sort = NULL;
    }
RETURN:
    return (argptr);
}




/*Translated by FOR_C++, v0.1, on 08/08/91 at 15:41:06 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/08/91 at 15:41:04
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SVIGP

    Computer:   FORC/SINGLE

    Revised:    May 30, 1991

    Purpose:    Sort an integer array by algebraically increasing value
                and return the permutation that rearranges the array.

    Usage:      CALL SVIGP (N, IA, IB, IPERM)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       IA     - Integer vector of length N containing the array to be
                sorted.  (Input)
       IB     - Integer vector of length N containing the sorted array.
                (Output)
                If IA is not needed, IA and IB can share the same
                storage locations.
       IPERM  - Vector of length N.  (Input/Output)
                On input, IPERM should be initialized to the values 1,
                2, ..., N.  On output, IPERM contains a record of
                permutations made on the vector IA.

    Remark:
       For wider applicability, integers (1,2,...,N) that are to be
       associated with IA(I) for I = 1, 2, ..., N may be entered into
       IPERM(I) in any order.  Note that these integers must be unique.

    Keyword:    Utilities

    GAMS:       N6a1a; N6a2a

    Chapters:   MATH/LIBRARY Utilities
                STAT/LIBRARY Utilities

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_svigp (Mint *n, Mint ia[], Mint ib[], Mint iperm[])
#else
static void l_svigp (n, ia, ib, iperm)
    Mint        *n, ia[], ib[], iperm[];
#endif
{
    Mint         i, ij, il[21], it, itt, iu[21], j, k, l, m, mt, mtt;
    Mfloat       r;


    imsl_e1psh ("SVIGP ");
    icopy (*n, ia, 1, ib, 1);

    m = 1;
    i = 1;
    j = *n;
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
     * SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION MT
     */
    ij = i + (j - i) * r;
    mt = ib[ij - 1];
    it = iperm[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN MT, INTERCHANGE WITH MT
     */
    if (ib[i - 1] > mt) {
	ib[ij - 1] = ib[i - 1];
	ib[i - 1] = mt;
	mt = ib[ij - 1];
	iperm[ij - 1] = iperm[i - 1];
	iperm[i - 1] = it;
	it = iperm[ij - 1];
    }
    l = j;
    /*
     * IF LAST ELEMENT OF ARRAY IS LESS THAN MT, INTERCHANGE WITH MT
     */
    if (ib[j - 1] >= mt)
	goto L_40;
    ib[ij - 1] = ib[j - 1];
    ib[j - 1] = mt;
    mt = ib[ij - 1];
    iperm[ij - 1] = iperm[j - 1];
    iperm[j - 1] = it;
    it = iperm[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN MT, INTERCHANGE WITH MT
     */
    if (ib[i - 1] <= mt)
	goto L_40;
    ib[ij - 1] = ib[i - 1];
    ib[i - 1] = mt;
    mt = ib[ij - 1];
    iperm[ij - 1] = iperm[i - 1];
    iperm[i - 1] = it;
    it = iperm[ij - 1];
    goto L_40;
L_30:
    if (ib[l - 1] == ib[k - 1])
	goto L_40;
    mtt = ib[l - 1];
    ib[l - 1] = ib[k - 1];
    ib[k - 1] = mtt;
    itt = iperm[l - 1];
    iperm[l - 1] = iperm[k - 1];
    iperm[k - 1] = itt;
    /*
     * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER THAN
     * MT
     */
L_40:
    l -= 1;
    if (ib[l - 1] > mt)
	goto L_40;
    /*
     * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER THAN
     * MT
     */
L_50:
    k += 1;
    if (ib[k - 1] < mt)
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
    mt = ib[i];
    it = iperm[i];
    if (ib[i - 1] <= mt)
	goto L_90;
    k = i;
L_100:
    ib[k] = ib[k - 1];
    iperm[k] = iperm[k - 1];
    k -= 1;
    if (mt < ib[k - 1])
	goto L_100;
    ib[k] = mt;
    iperm[k] = it;
    goto L_90;

L_9000:
    imsl_e1pop ("SVIGP ");
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  SVIBP

    Computer:   FORC/SINGLE

    Revised:    May 31, 1991

    Purpose:    Sort an integer array by nondecreasing absolute value
                and return the permutation that rearranges the array.

    Usage:      CALL SVIBP (N, IA, IB, IPERM)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       IA     - Integer vector of length N containing the array to be
                sorted.  (Input)
       IB     - Integer vector of length N containing the sorted array.
                (Output)
                If IA is not needed, IA and IB can share the same
                storage locations.
       IPERM  - Vector of length N.  (Input/Output)
                On input, IPERM should be initialized to the values 1,
                2, ..., N.  On output, IPERM contains a record of
                permutations made on the vector IA.

    Remark:
       For wider applicability, integers (1,2,...,N) that are to be
       associated with IA(I) for I = 1, 2, ..., N may be entered into
       IPERM(I) in any order.  Note that these integers must be unique.

    Keyword:    Utilities

    GAMS:       N6a1a; N6a2a

    Chapter:    MATH/LIBRARY Utilities

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_svibp (Mint *n, Mint ia[], Mint ib[], Mint iperm[])
#else
static void l_svibp (n, ia, ib, iperm)
    Mint        *n, ia[], ib[], iperm[];
#endif
{
    Mint         i, ij, ik, il[21], it, itt, iu[21], j, k, l, m, mt, mtt;
    Mfloat       r;


    imsl_e1psh ("SVIBP ");
    icopy (*n, ia, 1, ib, 1);
    /* FIND ABSOLUTE VALUES OF ARRAY IA */
    for (i = 1; i <= *n; i++) {
	if (ib[i - 1] < 0) {
	    iperm[i - 1] = -iperm[i - 1];
	    ib[i - 1] = -ib[i - 1];
	}
    }

    m = 1;
    i = 1;
    j = *n;
    r = .375;
L_20:
    if (i == j)
	goto L_90;

L_30:
    if (r <= .5898437) {
	r += 3.90625e-2;
    }
    else {
	r -= .21875;
    }

L_40:
    k = i;
    /*
     * SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION MT
     */
    ij = i + (j - i) * r;
    mt = ib[ij - 1];
    it = iperm[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN MT, INTERCHANGE WITH MT
     */
    if (ib[i - 1] > mt) {
	ib[ij - 1] = ib[i - 1];
	ib[i - 1] = mt;
	mt = ib[ij - 1];
	iperm[ij - 1] = iperm[i - 1];
	iperm[i - 1] = it;
	it = iperm[ij - 1];
    }
    l = j;
    /*
     * IF LAST ELEMENT OF ARRAY IS LESS THAN MT, INTERCHANGE WITH MT
     */
    if (ib[j - 1] >= mt)
	goto L_60;
    ib[ij - 1] = ib[j - 1];
    ib[j - 1] = mt;
    mt = ib[ij - 1];
    iperm[ij - 1] = iperm[j - 1];
    iperm[j - 1] = it;
    it = iperm[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN MT, INTERCHANGE WITH MT
     */
    if (ib[i - 1] <= mt)
	goto L_60;
    ib[ij - 1] = ib[i - 1];
    ib[i - 1] = mt;
    mt = ib[ij - 1];
    iperm[ij - 1] = iperm[i - 1];
    iperm[i - 1] = it;
    it = iperm[ij - 1];
    goto L_60;
L_50:
    if (ib[l - 1] == ib[k - 1])
	goto L_60;
    mtt = ib[l - 1];
    ib[l - 1] = ib[k - 1];
    ib[k - 1] = mtt;
    itt = iperm[l - 1];
    iperm[l - 1] = iperm[k - 1];
    iperm[k - 1] = itt;
    /*
     * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER THAN
     * MT
     */
L_60:
    l -= 1;
    if (ib[l - 1] > mt)
	goto L_60;
    /*
     * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER THAN
     * MT
     */
L_70:
    k += 1;
    if (ib[k - 1] < mt)
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
    mt = ib[i];
    it = iperm[i];
    if (ib[i - 1] <= mt)
	goto L_110;
    k = i;
L_120:
    ib[k] = ib[k - 1];
    iperm[k] = iperm[k - 1];
    k -= 1;
    if (mt < ib[k - 1])
	goto L_120;
    ib[k] = mt;
    iperm[k] = it;
    goto L_110;

L_130:
    for (ik = 1; ik <= *n; ik++) {
	if (iperm[ik - 1] < 0) {
	    iperm[ik - 1] = -iperm[ik - 1];
	    ib[ik - 1] = -ib[ik - 1];
	}
    }

L_9000:
    imsl_e1pop ("SVIBP ");
    return;
}				/* end of function */



/* Structured by FOR_STRUCT, v0.2, on 08/09/91 at 14:53:28
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SVIGN

    Computer:   FORC/SINGLE

    Revised:    March 19, 1991

    Purpose:    Sort an integer array by algebraically increasing
                value.

    Usage:      CALL SVIGN (N, IA, IB)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       IA     - Integer vector of length N containing the array to be
                sorted.  (Input)
       IB     - Integer vector of length N containing the sorted array.
                (Output)
                If IA is not needed, IA and IB can share the same
                storage locations.

    Keyword:    Utilities

    GAMS:       N6a2a

    Chapters:   MATH/LIBRARY Utilities
                STAT/LIBRARY Utilities

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_svign (Mint n, Mint ia[], Mint ib[])
#else
static void l_svign (n, ia, ib)
    Mint         n, ia[], ib[];
#endif
{
    Mint         i, ij, il[21], iu[21], j, k, l, m, mt, mtt;
    Mfloat       r;


    imsl_e1psh ("SVIGN ");
    icopy (n, ia, 1, ib, 1);

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
     * SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION MT
     */
    ij = i + (j - i) * r;
    mt = ib[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN MT, INTERCHANGE WITH MT
     */
    if (ib[i - 1] > mt) {
	ib[ij - 1] = ib[i - 1];
	ib[i - 1] = mt;
	mt = ib[ij - 1];
    }
    l = j;
    /*
     * IF LAST ELEMENT OF ARRAY IS LESS THAN MT, INTERCHANGE WITH MT
     */
    if (ib[j - 1] >= mt)
	goto L_50;
    ib[ij - 1] = ib[j - 1];
    ib[j - 1] = mt;
    mt = ib[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN MT, INTERCHANGE WITH MT
     */
    if (ib[i - 1] <= mt)
	goto L_50;
    ib[ij - 1] = ib[i - 1];
    ib[i - 1] = mt;
    mt = ib[ij - 1];
    goto L_50;
L_40:
    if (ib[l - 1] == ib[k - 1])
	goto L_50;
    mtt = ib[l - 1];
    ib[l - 1] = ib[k - 1];
    ib[k - 1] = mtt;
    /*
     * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER THAN
     * MT
     */
L_50:
    l -= 1;
    if (ib[l - 1] > mt)
	goto L_50;
    /*
     * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER THAN
     * MT
     */
L_60:
    k += 1;
    if (ib[k - 1] < mt)
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
    mt = ib[i];
    if (ib[i - 1] <= mt)
	goto L_100;
    k = i;
L_110:
    ib[k] = ib[k - 1];
    k -= 1;
    if (mt < ib[k - 1])
	goto L_110;
    ib[k] = mt;
    goto L_100;

L_9000:
    imsl_e1pop ("SVIGN ");
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  SVIBN

    Computer:   FORC/SINGLE

    Revised:    March 19, 1991

    Purpose:    Sort an integer array by nondecreasing absolute value.

    Usage:      CALL SVIBN (N, IA, IB)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       IA     - Integer vector of length N containing the array to be
                sorted.  (Input)
       IB     - Integer vector of length N containing the sorted array.
                (Output)
                If IA is not needed, IA and IB can share the same
                storage locations.

    Keyword:    Utilities

    GAMS:       N6a2a

    Chapters:   MATH/LIBRARY Utilities

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_svibn (Mint *n, Mint ia[], Mint ib[])
#else
static void l_svibn (n, ia, ib)
    Mint        *n, ia[], ib[];
#endif
{
    Mint         i, ij, il[21], iu[21], j, k, l, m, mt, mtt;
    Mfloat       r;


    imsl_e1psh ("SVIBN ");
    icopy (*n, ia, 1, ib, 1);

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
     * SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION MT
     */
    ij = i + (j - i) * r;
    mt = ib[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN MT, INTERCHANGE WITH MT
     */
    if (abs (ib[i - 1]) <= abs (mt))
	goto L_40;
    ib[ij - 1] = ib[i - 1];
    ib[i - 1] = mt;
    mt = ib[ij - 1];
L_40:
    l = j;
    /*
     * IF LAST ELEMENT OF ARRAY IS LESS THAN MT, INTERCHANGE WITH MT
     */
    if (abs (ib[j - 1]) >= abs (mt))
	goto L_60;
    ib[ij - 1] = ib[j - 1];
    ib[j - 1] = mt;
    mt = ib[ij - 1];
    /*
     * IF FIRST ELEMENT OF ARRAY IS GREATER THAN MT, INTERCHANGE WITH MT
     */
    if (abs (ib[i - 1]) <= abs (mt))
	goto L_60;
    ib[ij - 1] = ib[i - 1];
    ib[i - 1] = mt;
    mt = ib[ij - 1];
    goto L_60;
L_50:
    if (abs (ib[l - 1]) == abs (ib[k - 1]))
	goto L_60;
    mtt = ib[l - 1];
    ib[l - 1] = ib[k - 1];
    ib[k - 1] = mtt;
    /*
     * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER THAN
     * MT
     */
L_60:
    l -= 1;
    if (abs (ib[l - 1]) > abs (mt))
	goto L_60;
    /*
     * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER THAN
     * MT
     */
L_70:
    k += 1;
    if (abs (ib[k - 1]) < abs (mt))
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
    mt = ib[i];
    if (abs (ib[i - 1]) <= abs (mt))
	goto L_110;
    k = i;
L_120:
    ib[k] = ib[k - 1];
    k -= 1;
    if (abs (mt) < abs (ib[k - 1]))
	goto L_120;
    ib[k] = mt;
    goto L_110;

L_9000:
    imsl_e1pop ("SVIBN ");
    return;
}				/* end of function */
