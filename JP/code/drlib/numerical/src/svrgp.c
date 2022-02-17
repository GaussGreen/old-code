#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/*
  -----------------------------------------------------------------------
    IMSL Name:  SVRGP/DSVRGP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 19, 1985

    Purpose:    Sort a real array by algebraically increasing value
                and return the permutation that rearranges the array.

    Usage:      CALL SVRGP (N, RA, RB, IPERM)

    Arguments:
       N      - Number of elements in the array to be sorted.  (Input)
       RA     - Vector of length N containing the array to be sorted.
                (Input)
       RB     - Vector of length N containing the sorted array.  (Output)
                If RA is not needed, RA and RB can share the same
                storage locations.
       IPERM  - Vector of length N.  (Input/Output)
                On input IPERM should be initialized to the values
                1, 2, ..., N.  On output IPERM contains a record of
                permutations made on the vector RA.

    Remark:
       For wider applicability, integers (1, 2, ..., N) that are to be
       associated with RA(I) for I = 1, 2, ..., N may be entered into
       IPERM(I) in any order.  Note that these integers must be unique.

    Keyword:    Utilities

    GAMS:       N6a1b; N6a2b

    Chapters:   MATH/LIBRARY Utilities
                STAT/LIBRARY Utilities

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_svrgp(Mint n, Mfloat ra[], Mfloat rb[], Mint iperm[])
#else
void imsl_svrgp(n, ra, rb, iperm)
	Mint            n;
	Mfloat          ra[], rb[];
	Mint            iperm[];
#endif
{
	Mint            i, ij, il[21], it, itt, iu[21], j, k, l, m;
	Mfloat          r, t, tt;


	imsl_e1psh("imsl_svrgp");
	scopy(n, ra, 1, rb, 1);
	/* CHECK FOR INPUT ERRORS */
	if (n <= 0) {
		imsl_e1sti(1, n);

/*		(5, 1, "The length of vectors RA, RB and IPERM must be greater than or equal to one.  On input the length, N, is given as %(i1).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LENGTH_OF_VECTORS_3);
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;

	m = 1;
	i = 1;
	j = n;
	r = .375;
L_10:
	if (i == j)
		goto L_70;
	if (r <= .5898437) {
		r += 3.90625e-2;
	} else {
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
	 * FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER
	 * THAN T
	 */
L_40:
	l -= 1;
	if (rb[l - 1] > t)
		goto L_40;
	/*
	 * FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER
	 * THAN T
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
	imsl_e1pop("imsl_svrgp");
	return;
}				/* end of function */

