#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#include <limits.h>

#ifdef ANSI
static void l_f_mod_m1ran(Mint, Mint, Mfloat *, Mfloat *);
static Mint l_prod_mod(Mint, Mint, Mint);
#else
static void l_f_mod_m1ran();
static Mint l_prod_mod();
#endif

/* -----------------------------------------------------------------------
    Purpose:    Transpose of a rectangular matrix.

    Usage:      CALL M1RAN (NRA, NCA, A, ATRAN)

    Arguments:
       NRA    - Number of rows of A.  (Input)
       NCA    - Number of columns of A.  (Input)
       A      - Real vector of length NRA*NCA containing the NRA by NCA
                matrix to be transposed.  (Input)
       ATRAN  - Real vector of length NRA*NCA containing the NCA by NRA
                transposed matrix.  (Input)
                If A is not needed, ATRAN can share the same storage
                locations as A.
  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_f_m1ran(Mint nra, Mint nca, Mfloat *a, Mfloat *atran)
#else
void imsl_f_m1ran(nra, nca, a, atran)
	Mint             nra, nca;
	Mfloat           a[], atran[];
#endif
{
	Mfloat 		*tmp_tran, temp;
	Mint		i, j, nram, k, nrap1, nram1, l;

                                        /* CHECK INPUT */
        imsl_e1psh("M1RAN_F");
        if (nra <= 0) {
                imsl_e1sti(1, nra);
/*              (5, 1, "The number of rows must be greater than zero.  NRA = %(i1)");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NOT_ENOUGH_ROWS);
        }
        if (nca <= 0) {
                imsl_e1sti(1, nca);
/*              (5, 2, "The number of columns must be greater than zero.  NCA = %(i1)");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NOT_ENOUGH_COLUMNS);
        }
        if (imsl_n1rcd(0) != 0)  goto L_9000;

                                /* CHECK TO SEE IF A IS SQUARE */
        nram = nra * nca - 1;
        if (nra != nca)
                goto L_30;
                                /* SQUARE MATRICES ARE HANDLED SEPARATELY */

        if (atran != a) {
#if defined(COMPUTER_PMXUX)
                memcpy(atran, a, sizeof(Mfloat)*nra*nca);
#elif defined(ANSI)
                memcpy((void*)atran, (const void*)a, sizeof(Mfloat)*nra*nca);
#else
                memcpy(atran, a, sizeof(Mfloat)*nra*nca);
#endif
        }

        k = 2;
        nrap1 = nra + 1;
        nram1 = nra - 1;
                                /* PERFORM TRANSPOSE OPERATION */
        for (i = nra; i<=nram; i += nra) {
                l = k + nram1;
                for (j = k; j <= i; j++) {
                        temp = atran[j-1];
                        atran[j-1] = atran[l-1];
                        atran[l-1] = temp;
                        l += nra;
                }
                k += nrap1;
        }
        goto L_9000;

	/* A is not square */
L_30:
	if (a == atran) {
		tmp_tran = (Mfloat*) imsl_malloc(nca*nra*sizeof(Mfloat));
		if (!tmp_tran) {

			/* Not enough space for temporary copy of transpose.
			   We will switch to modulo arithmetic method of
		 	   transposing a */

			l_f_mod_m1ran(nra, nca, a, atran);
			goto L_9000;
		}
	} else {
		tmp_tran = atran;
	}

        for (i = 0; i < nra; i++) {
                for (j = 0; j < nca; j++) 
			tmp_tran[i+nra*j] = a[j+nca*i];
        }
	if (a == atran) {
#if defined(COMPUTER_PMXUX)
                memcpy(atran, tmp_tran, sizeof(*a)*nra*nca);
#elif defined(ANSI)
                memcpy((void*)atran, (const void*)tmp_tran, sizeof(*a)*nra*nca);
#else
                memcpy(atran, tmp_tran, sizeof(*a)*nra*nca);
#endif
		imsl_free(tmp_tran);
	} 

L_9000:
	imsl_e1pop("M1RAN_F");
	return;
}

/* -----------------------------------------------------------------------
    Purpose:    Transpose of a rectangular matrix, using modulo
		arithmetic.

    Usage:      CALL M1RAN (NRA, NCA, A, ATRAN)

    Arguments:
       NRA    - Number of rows of A.  (Input)
       NCA    - Number of columns of A.  (Input)
       A      - Real vector of length NRA*NCA containing the NRA by NCA
                matrix to be transposed.  (Input)
       ATRAN  - Real vector of length NRA*NCA containing the NCA by NRA
                transposed matrix.  (Input)
                If A is not needed, ATRAN can share the same storage
                locations as A.
  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f_mod_m1ran(Mint nca, Mint nra, Mfloat *a, Mfloat *atran)
#else
static void l_f_mod_m1ran(nca, nra, a, atran)
	Mint             nra, nca;
	Mfloat           a[], atran[];
#endif
{
	Mint            i, id, iexp[13], ipf[13],
                        ipwr[13], isd, istart, it, itemp[13], j, k, kk, l,
                        ll, npf, nram, nram1, nramis, nrap1, too_big = 0;
	Mfloat           temp, temp1;

					/* CHECK INPUT */
	imsl_e1psh("M1RAN_F");
	if (nra <= 0) {
		imsl_e1sti(1, nra);
/*		(5, 1, "The number of rows must be greater than zero.  NRA = %(i1)");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NOT_ENOUGH_ROWS);
	}
	if (nca <= 0) {
		imsl_e1sti(1, nca);
/*		(5, 2, "The number of columns must be greater than zero.  NCA = %(i1)");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NOT_ENOUGH_COLUMNS);
	}
	if (imsl_n1rcd(0) != 0)  goto L_9000;
	if (nra*nca > imsl_i_machine(5)/nra || 
		nca*nra > imsl_i_machine(5)/nca) too_big = 1;
					/* COPY A TO ATRAN */

	if ( &atran[0] != &a[0]) {
/*		for (i=0;i<nra*nca;i++)  atran[i] = a[i];  */
#if defined(COMPUTER_PMXUX)
		memcpy(atran, a, sizeof(*a)*nra*nca);
#elif defined(ANSI) 
		memcpy((void*)atran, (const void*)a, sizeof(*a)*nra*nca);
#else
		memcpy(atran, a, sizeof(*a)*nra*nca);
#endif 
	}
				/* IF A IS ONE COLUMN OR ONE ROW RETURN */
	if (nra <= 1 || nca <= 1) goto L_9000;
				/* CHECK TO SEE IF A IS SQUARE */
	nram = nra * nca - 1;
	if (nra != nca)
		goto L_30;
				/* SQUARE MATRICES ARE HANDLED SEPARATELY */
	k = 2;
	nrap1 = nra + 1;
	nram1 = nra - 1;
				/* PERFORM TRANSPOSE OPERATION */
	for (i = nra; i<=nram; i += nra) {
		l = k + nram1;
		for (j = k; j <= i; j++) {
			temp = atran[j-1];
			atran[j-1] = atran[l-1];
			atran[l-1] = temp;
			l += nra;
		}
		k += nrap1;
	}
	goto L_9000;
	/*
	 * A IS NOT SQUARE. TRANSPOSITION OF AN NRA BY NCA MATRIX AMOUNTS TO
	 * REPLACING THE ELEMENT AT VECTOR POSITION I WITH THE ELEMENT AT
	 * POSITION NRA*I (MOD NRA*NCA-1). EACH SUBCYCLE OF THIS PERMUTATION
	 * IS COMPLETED IN ORDER. DECOMPOSE NRAM INTO ITS PRIME FACTORS
	 */
L_30:
	imsl_prime(nram, &npf, ipf, iexp, ipwr);
	for (i = 1; i <= npf; i++) {
		itemp[i - 1] = 0;
	}
				/* GENERATE DIVISORS OF NRAM LESS THAN NRAM/2 */
	id = 1;
L_50:
	if (id >= nram / 2)
		goto L_9000;
	j = nram / id;
	for (i = 1; i <= npf; i++) {
		if (itemp[i - 1] != iexp[i - 1])
			j = (j / ipf[i - 1]) * (ipf[i - 1] - 1);
	}
	/*
	 * THE STARTING POINT OF A SUBCYCLE IS DIVISIBLE ONLY BY ID AND MUST
	 * NOT APPEAR IN ANY OTHER SUBCYCLE
	 */
	istart = id;
L_70:
	nramis = nram - istart;
	if (istart == id)
		goto L_100;
	isd = istart / id;
	for (i = 1; i <= npf; i++) {
		if (itemp[i - 1] == iexp[i - 1])
			goto L_80;
		if (mod(isd, ipf[i - 1]) == 0)
			goto L_140;
L_80:
		;
	}
	it = istart;
L_90:
	if (too_big)
		it = l_prod_mod(nra, it, nram);
	else
		it = mod(nra * it, nram); 
	if (it < istart || it > nramis)
		goto L_140;
	if (it > istart && it < nramis)
		goto L_90;
L_100:
	temp = atran[istart];
	temp1 = atran[nramis];
	k = istart;
L_110:
	if (too_big)
		kk = l_prod_mod(nra, k, nram);
	else
		kk = mod(nra * k, nram);
	l = nram - k;
	ll = nram - kk;
	j -= 2;
	/*
	 * M1VE TWO ELEMENTS. THE SECOND FROM THE NEGATIVE SUBCYCLE. CHECK TO
	 * SEE IF THE SUBCYCLE IS COMPLETE.
	 */
	if (kk == istart)
		goto L_120;
	if (ll == istart)
		goto L_130;
	/*
	 * SUBCYCLE NOT COMPLETE. M1VE ELEMENTS AND UPDATE POINTERS
	 */
	atran[k] = atran[kk];
	atran[l] = atran[ll];
	k = kk;
	goto L_110;
	/* SUBCYCLE COMPLETE. RECOMPUTE ID */
L_120:
	atran[k] = temp;
	atran[l] = temp1;
	goto L_140;
L_130:
	atran[k] = temp1;
	atran[l] = temp;
L_140:
	istart += id;
	if (j > 0)
		goto L_70;
	for (i = 1; i <= npf; i++) {
		if (itemp[i - 1] == iexp[i - 1])
			goto L_150;
		itemp[i - 1] += 1;
		id *= ipf[i - 1];
		goto L_50;
L_150:
		itemp[i - 1] = 0;
		id /= ipwr[i - 1];
	}

L_9000:
	imsl_e1pop("M1RAN_F");
	return;
}


#ifdef ANSI
void imsl_c_m1ran(Mint nca, Mint nra, Mf_complex *a, Mf_complex *atran)
#else
void imsl_c_m1ran(nca, nra, a, atran)
	Mint             nra, nca;
	Mf_complex        a[], atran[];
#endif
{
	Mint            i, id, iexp[13], ipf[13],
                        ipwr[13], isd, istart, it, itemp[13], j, k, kk, l,
                        ll, npf, nram, nram1, nramis, nrap1;
	Mint too_big = 0;
	Mf_complex        temp, temp1;

					/* CHECK INPUT */
	imsl_e1psh("M1RAN_C");
	if (nra <= 0) {
		imsl_e1sti(1, nra);
/*		(5, 1, "The number of rows must be greater than zero.  NRA = %(i1)");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NOT_ENOUGH_ROWS);
	}
	if (nca <= 0) {
		imsl_e1sti(1, nca);
/*		(5, 2, "The number of columns must be greater than zero.  NCA = %(i1)");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NOT_ENOUGH_COLUMNS);
	}
	if (imsl_n1rcd(0) != 0)  goto L_9000;
	if (nra*nca > imsl_i_machine(5)/nra || 
		nca*nra > imsl_i_machine(5)/nca) too_big = 1;

					/* COPY A TO ATRAN */
	if (&atran[0] != &a[0]) {
/*		for (i=0;i<nra*nca;i++)  atran[i] = a[i]; */
#if defined(COMPUTER_PMXUX)
                memcpy(atran, a, sizeof(*a)*nra*nca);
#elif defined(ANSI)
                memcpy((void*)atran, (const void*)a, sizeof(*a)*nra*nca);
#else
                memcpy(atran, a, sizeof(*a)*nra*nca);
#endif
	}
				/* IF A IS ONE COLUMN OR ONE ROW RETURN */
	if (nra <= 1 || nca <= 1) goto L_9000;
				/* CHECK TO SEE IF A IS SQUARE */
	nram = nra * nca - 1;
	if (nra != nca) goto L_30;
				/* SQUARE MATRICES ARE HANDLED SEPARATELY */
	k = 2;
	nrap1 = nra + 1;
	nram1 = nra - 1;
				/* PERFORM TRANSPOSE OPERATION */
	for (i = nra; i<=nram; i += nra) {
		l = k + nram1;
		for (j = k; j <= i; j++) {
			temp = atran[j-1];
			atran[j-1] = atran[l-1];
			atran[l-1] = temp;
			l += nra;
		}
		k += nrap1;
	}
	goto L_9000;
	/*
	 * A IS NOT SQUARE. TRANSPOSITION OF AN NRA BY NCA MATRIX AMOUNTS TO
	 * REPLACING THE ELEMENT AT VECTOR POSITION I WITH THE ELEMENT AT
	 * POSITION NRA*I (MOD NRA*NCA-1). EACH SUBCYCLE OF THIS PERMUTATION
	 * IS COMPLETED IN ORDER. DECOMPOSE NRAM INTO ITS PRIME FACTORS
	 */
L_30:
	imsl_prime(nram, &npf, ipf, iexp, ipwr);
	for (i = 1; i <= npf; i++) {
		itemp[i - 1] = 0;
	}
				/* GENERATE DIVISORS OF NRAM LESS THAN NRAM/2 */
	id = 1;
L_50:
	if (id >= nram / 2) goto L_9000;
	j = nram / id;
	for (i = 1; i <= npf; i++) {
		if (itemp[i-1] != iexp[i-1])
			j = (j / ipf[i - 1]) * (ipf[i - 1] - 1);
	}
	/*
	 * THE STARTING POINT OF A SUBCYCLE IS DIVISIBLE ONLY BY ID AND MUST
	 * NOT APPEAR IN ANY OTHER SUBCYCLE
	 */
	istart = id;
L_70:
	nramis = nram - istart;
	if (istart == id) goto L_100;
	isd = istart / id;
	for (i = 1; i <= npf; i++) {
		if (itemp[i - 1] == iexp[i - 1]) goto L_80;
		if (mod(isd, ipf[i - 1]) == 0) goto L_140;
L_80:   ;
	}
	it = istart;
L_90:
	if (too_big)
		it = l_prod_mod(nra, it, nram);
	else
		it = mod(nra * it, nram); 
	if (it < istart || it > nramis) goto L_140;
	if (it > istart && it < nramis) goto L_90;
L_100:
	temp = atran[istart];
	temp1 = atran[nramis];
	k = istart;
L_110:
	if (too_big)
		kk = l_prod_mod(nra, k, nram);
	else
		kk = mod(nra*k, nram); 
	l = nram - k;
	ll = nram - kk;
	j -= 2;
	/*
	 * M1VE TWO ELEMENTS. THE SECOND FROM THE NEGATIVE SUBCYCLE. CHECK TO
	 * SEE IF THE SUBCYCLE IS COMPLETE.
	 */
	if (kk == istart) goto L_120;
	if (ll == istart) goto L_130;
	/*
	 * SUBCYCLE NOT COMPLETE. M1VE ELEMENTS AND UPDATE POINTERS
	 */
	atran[k] = atran[kk];
	atran[l] = atran[ll];
	k = kk;
	goto L_110;
	                            /* SUBCYCLE COMPLETE. RECOMPUTE ID */
L_120:  atran[k] = temp;
	atran[l] = temp1;
	goto L_140;
L_130:  atran[k] = temp1;
	atran[l] = temp;
L_140:  istart += id;
	if (j > 0) goto L_70;
	for (i = 1; i <= npf; i++) {
		if (itemp[i - 1] == iexp[i - 1]) goto L_150;
		itemp[i - 1] += 1;
		id *= ipf[i - 1];
		goto L_50;
L_150:          itemp[i - 1] = 0;
		id /= ipwr[i - 1];
	}

L_9000:
	imsl_e1pop("M1RAN_C");
	return;
}

#ifndef DOUBLE

#ifdef ANSI
void imsl_i_m1ran(Mint nca, Mint nra, Mint *a, Mint *atran)
#else
void imsl_i_m1ran(nca, nra, a, atran)
	Mint             nra, nca;
	Mint             a[], atran[];
#endif
{
	Mint            i, id, iexp[13], ipf[13],
                        ipwr[13], isd, istart, it, itemp[13], j, k, kk, l,
                        ll, npf, nram, nram1, nramis, nrap1;
	Mint too_big = 0;
	Mint             temp, temp1;

					/* CHECK INPUT */
	imsl_e1psh("M1RAN_I");
	if (nra <= 0) {
		imsl_e1sti(1, nra);
/*		(5, 1, "The number of rows must be greater than zero.  NRA = %(i1)");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NOT_ENOUGH_ROWS);
	}
	if (nca <= 0) {
		imsl_e1sti(1, nca);
/*		(5, 2, "The number of columns must be greater than zero.  NCA = %(i1)");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NOT_ENOUGH_COLUMNS);
	}
	if (imsl_n1rcd(0) != 0)  goto L_9000;
	if (nra*nca > imsl_i_machine(5)/nra || 
		nca*nra > imsl_i_machine(5)/nca) too_big = 1;
					/* COPY A TO ATRAN */
        if (&atran[0] != &a[0]) {
/*		for (i=0;i<nra*nca;i++)  atran[i] = a[i]; */
#if defined(COMPUTER_PMXUX)
		memcpy(atran, a, sizeof(*a)*nra*nca);
#elif defined(ANSI)
		memcpy((void*)atran, (const void*)a, sizeof(*a)*nra*nca);
#else
		memcpy(atran, a, sizeof(*a)*nra*nca);
#endif
	}

				/* IF A IS ONE COLUMN OR ONE ROW RETURN */
	if (nra <= 1 || nca <= 1) goto L_9000;
				/* CHECK TO SEE IF A IS SQUARE */
	nram = nra * nca - 1;
	if (nra != nca) goto L_30;
				/* SQUARE MATRICES ARE HANDLED SEPARATELY */
	k = 2;
	nrap1 = nra + 1;
	nram1 = nra - 1;
				/* PERFORM TRANSPOSE OPERATION */
	for (i = nra; i<=nram; i += nra) {
		l = k + nram1;
		for (j = k; j <= i; j++) {
			temp = atran[j-1];
			atran[j-1] = atran[l-1];
			atran[l-1] = temp;
			l += nra;
		}
		k += nrap1;
	}
	goto L_9000;
	/*
	 * A IS NOT SQUARE. TRANSPOSITION OF AN NRA BY NCA MATRIX AMOUNTS TO
	 * REPLACING THE ELEMENT AT VECTOR POSITION I WITH THE ELEMENT AT
	 * POSITION NRA*I (MOD NRA*NCA-1). EACH SUBCYCLE OF THIS PERMUTATION
	 * IS COMPLETED IN ORDER. DECOMPOSE NRAM INTO ITS PRIME FACTORS
	 */
L_30:
	imsl_prime(nram, &npf, ipf, iexp, ipwr);
	for (i = 1; i <= npf; i++) {
		itemp[i - 1] = 0;
	}
				/* GENERATE DIVISORS OF NRAM LESS THAN NRAM/2 */
	id = 1;
L_50:
	if (id >= nram / 2) goto L_9000;
	j = nram / id;
	for (i = 1; i <= npf; i++) {
		if (itemp[i-1] != iexp[i-1])
			j = (j / ipf[i - 1]) * (ipf[i - 1] - 1);
	}
	/*
	 * THE STARTING POINT OF A SUBCYCLE IS DIVISIBLE ONLY BY ID AND MUST
	 * NOT APPEAR IN ANY OTHER SUBCYCLE
	 */
	istart = id;
L_70:
	nramis = nram - istart;
	if (istart == id) goto L_100;
	isd = istart / id;
	for (i = 1; i <= npf; i++) {
		if (itemp[i - 1] == iexp[i - 1]) goto L_80;
		if (mod(isd, ipf[i - 1]) == 0) goto L_140;
L_80:   ;
	}
	it = istart;
L_90:
	if (too_big)
		it = l_prod_mod(nra, it, nram);
	else
		it = mod(nra * it, nram); 
	if (it < istart || it > nramis) goto L_140;
	if (it > istart && it < nramis) goto L_90;
L_100:
	temp = atran[istart];
	temp1 = atran[nramis];
	k = istart;
L_110:
	if (too_big)
		kk = l_prod_mod(nra, k, nram);
	else
		kk = mod(nra*k, nram); 
	l = nram - k;
	ll = nram - kk;
	j -= 2;
	/*
	 * M1VE TWO ELEMENTS. THE SECOND FROM THE NEGATIVE SUBCYCLE. CHECK TO
	 * SEE IF THE SUBCYCLE IS COMPLETE.
	 */
	if (kk == istart) goto L_120;
	if (ll == istart) goto L_130;
	/*
	 * SUBCYCLE NOT COMPLETE. M1VE ELEMENTS AND UPDATE POINTERS
	 */
	atran[k] = atran[kk];
	atran[l] = atran[ll];
	k = kk;
	goto L_110;
	                            /* SUBCYCLE COMPLETE. RECOMPUTE ID */
L_120:  atran[k] = temp;
	atran[l] = temp1;
	goto L_140;
L_130:  atran[k] = temp1;
	atran[l] = temp;
L_140:  istart += id;
	if (j > 0) goto L_70;
	for (i = 1; i <= npf; i++) {
		if (itemp[i - 1] == iexp[i - 1]) goto L_150;
		itemp[i - 1] += 1;
		id *= ipf[i - 1];
		goto L_50;
L_150:          itemp[i - 1] = 0;
		id /= ipwr[i - 1];
	}

L_9000:
	imsl_e1pop("M1RAN_I");
	return;
}


/* Structured by FOR_STRUCT, v0.2, on 06/12/90 at 10:21:40
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  PRIME (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 30, 1985

    Purpose:    Decompose an integer into its prime factors.

    Usage:      CALL PRIME (N, NPF, IPF, IEXP, IPW)

    Arguments:
       N      - Integer to be decomposed.  (Input)
       NPF    - Number of different prime factors of ABS(N).  (Output)
                If N is equal to -1, 0, or 1, NPF is set to 0.
       IPF    - Integer vector of length 13.  (Output)
                IPF(I) contains the prime factors of the absolute value
                of N, for I=1,...,NPF.  The remaining 13-NPF locations
                are not used.
       IEXP   - Integer vector of length 13.  (Output)
                IEXP(I) is the exponent of IPF(I), for I=1,...,NPF.  The
                remaining 13-NPF locations are not used.
       IPW    - Integer vector of length 13.  (Output)
                IPW(I) contains the quantity IPF(I)**IEXP(I),
                for I=1,...,NPF.  The remaining 13-NPF locations are not
                used.

    Remark:
       The output from PRIME should be interpreted in the following way:
          ABS(N) = IPF(1)**IEXP(1) * .... * IPF(NPF)**IEXP(NPF).

    Keyword:    Utilities

    GAMS:       A6c

    Chapter:    MATH/LIBRARY Utilities

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_prime(Mint n, Mint *npf, Mint *ipf, Mint *iexp, Mint *ipw)
#else
void imsl_prime(n, npf, ipf, iexp, ipw)
	Mint             n, *npf, ipf[], iexp[], ipw[];
#endif
{
	Mint             i, id, ifc, iq, k, nf;


	i = 0;
	nf = abs(n);
	for (k = 1; k <= 13; k++) {
		ipf[k - 1] = 0;
		iexp[k - 1] = 0;
		ipw[k - 1] = 0;
	}
	/* IF N IS -1, 0, OR 1, SET NPF TO 0 */
	if (nf <= 1)
		goto L_50;
	ifc = 0;
	id = 2;
L_20:
	iq = nf / id;
	if (nf != id * iq)
		goto L_30;
	/* ID IS A FACTOR OF NF */
	nf = iq;
	if (id > ifc) {
		/* ID IS A NEW DIVISOR */
		i += 1;
		/*
		 * UPDATE THE PRIME FACTOR VECTOR, THE POWER VECTOR, AND SET
		 * THE EXPONENT TO 1
		 */
		ipf[i - 1] = id;
		ipw[i - 1] = id;
		iexp[i - 1] = 1;
		ifc = id;
		/*
		 * ID FACTORS N MORE THAN ONCE. UPDATE THE POWER VECTOR AND
		 * THE EXPONENT VECTOR.
		 */
	} else {
		ipw[i - 1] *= id;
		iexp[i - 1] += 1;
	}
	goto L_20;
	/*
	 * ID IS NOT A FACTOR OF NF. IF NF IS NOT PRIME, CONTINUE TO FACTOR.
	 */
L_30:
	if (iq <= id)
		goto L_40;
	/* UPDATE THE DIVISOR */
	if (id <= 2) {
		id = 3;
	} else {
		id += 2;
	}
	goto L_20;
	/* NF IS PRIME OR 1 */
L_40:
	if (nf > 1) {
		/*
		 * NF IS PRIME. UPDATE THE PRIME FACTOR VECTOR, THE POWER
		 * VECTOR, AND THE EXPONENT VECTOR.
		 */
		i += 1;
		ipf[i - 1] = nf;
		ipw[i - 1] = nf;
		iexp[i - 1] = 1;
	}
L_50:
	*npf = i;

	return;
}
#endif


				/* compute (ij) mod m */
#ifdef ANSI
static Mint l_prod_mod(Mint i, Mint j, Mint m)
#else
static Mint l_prod_mod(i, j, m)
    Mint		i;
    Mint		j;
    Mint		m;
#endif
{
				/* this code assumes that a double */
				/* has at least a few more bits of */
				/* precision than an int. */
				/* should check that eps_double*max_int < 0.001 */
    Mdouble	di, dj;
    Mdouble	i0, i1, j0, j1, m0, m1, q0, q1;
    Mdouble	dr, r0, r1, r2;
#ifdef COMPUTER_DECOSF
    unsigned int 	   q;
    static unsigned int    h, mask;
#else
    unsigned long 	   q;
    static unsigned long   h, mask;
#endif
    static Mint 		   shift=0;


				/* initialize */
    if (shift == 0) {
				/* compute h = sqrt(biggest int) */
#ifdef COMPUTER_DECOSF
	shift = CHAR_BIT*sizeof(int)/2;
#else
	shift = CHAR_BIT*sizeof(long int)/2;
#endif
	h     = 1 << shift;
	mask  = h - 1;
    }
				/* reduce i and j */
    di = i = i%m;
    dj = j = j%m;
				/* approximate q = floor[(ij)/m] */
    q = (int) (di*dj/m);
				/* split q, m, i and j */
				/* q = q1*h + q0, etc. */
    q0 = q & mask;	q1 = q >> shift;
    m0 = m & mask;	m1 = m >> shift;
    i0 = i & mask;	i1 = i >> shift;
    j0 = j & mask;	j1 = j >> shift;
				/* compute split form of mod */
				/* r = i*j - q*m */
    /* r2*h*h + r1*h + r0 = (i1*h+i0)*(j1*h+j0) - (q1*h+q0)*(m1*h+m0) */
    r0 = i0*j0 - q0*m0;
    r1 = (i1*j0+i0*j1) - (q1*m0+q0*m1);
    r2 = i1*j1 - q1*m1;
				/* sum parts */
    dr = (r2*h + r1) * h + r0;
				/* fix in case q was off by a little */
    while (dr > m) dr -= m;
    while (dr < 0) dr += m;
				/* round to int and return */
    return (int) (dr+0.5);
}
