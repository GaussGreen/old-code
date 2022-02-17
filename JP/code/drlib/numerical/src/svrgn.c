#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 04/30/90 at 17:11:43
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SVRGN/DSVRGN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 18, 1988

    Purpose:    Sort a real array by algebraically increasing value.

    Usage:      CALL SVRGN (N, RA, RB)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       RA     - Vector of length N containing the array to be sorted.
                (Input)
       RB     - Vector of length N containing the sorted array.  (Output)
                If RA is not needed, RA and RB can share the same
                storage locations.

    Keyword:    Utilities

    GAMS:       N6a2b

    Chapters:   MATH/LIBRARY Utilities
                STAT/LIBRARY Utilities

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_svrgn(Mint n, Mfloat *ra, Mfloat *rb)
#else
void imsl_svrgn(n, ra, rb)
	Mint             n;
	Mfloat           ra[], rb[];
#endif
{
	Mint             i, ij, il[21], iu[21], j, k, l, m;
	Mfloat           r, t, tt;


	imsl_e1psh("SVRGN");

	scopy(n, ra, 1, rb, 1);
	/* Check input argument */
	if (n <= 0) {
		imsl_e1sti(1, n);

/*		(5, 1, "N = %(i1).  The length of RA and RB, N, must be greater than or equal to one.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LENGTH_OF_VECTORS_2);
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;

	m = 1;
	i = 1;
	j = n;
	r = .375;
L_10:
	if (i == j)
		goto L_80;
	if (r <= .5898437) {
		r += 3.90625e-2;
	} else {
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
	 * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER
	 * THAN T
	 */
L_50:
	l -= 1;
	if (rb[l - 1] > t)
		goto L_50;
	/*
	 * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER
	 * THAN T
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
	imsl_e1pop("SVRGN");

	return;
}				/* end of function */
