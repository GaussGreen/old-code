#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK	l_quadratic_prog(Mint, Mint, Mint, Mfloat*, Mfloat*,
                                 Mfloat*, Mfloat*, va_list argptr);
static Mfloat   l_a1ot(Mint, Mfloat*, Mint, Mfloat*, Mint);
#else
static VA_LIST_HACK  l_quadratic_prog();
static Mfloat   l_a1ot();
#endif

static Mfloat	*lv_x;

#ifdef ANSI
Mfloat *imsl_f_quadratic_prog(Mint m, Mint n, Mint meq, Mfloat *a,
                              Mfloat *b, Mfloat *g, Mfloat *h, ...)
#else
Mfloat *imsl_f_quadratic_prog(m, n, meq, a, b, g, h, va_alist)
    Mint	m;
    Mint	n;
    Mint	meq;
    Mfloat	*a;
    Mfloat	*b;
    Mfloat	*g;
    Mfloat	*h;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr,h);

    E1PSH("imsl_f_quadratic_prog","imsl_d_quadratic_prog");
    lv_x = NULL;
    IMSL_CALL(l_quadratic_prog(m, n, meq, a, b, g, h, argptr));
    va_end(argptr);
    E1POP("imsl_f_quadratic_prog","imsl_d_quadratic_prog");
    return lv_x;
}


#ifdef ANSI
static VA_LIST_HACK l_quadratic_prog(Mint m, Mint n, Mint meq, Mfloat *a,
                                Mfloat *b, Mfloat *g, Mfloat *h,
                                va_list argptr)
#else
static VA_LIST_HACK l_quadratic_prog(m, n, meq, a, b, g, h, argptr)
    Mint	m;
    Mint	n;
    Mint	meq;
    Mfloat	*a;
    Mfloat	*b;
    Mfloat	*g;
    Mfloat	*h;
    va_list	argptr;
#endif
{
    Mint	    code;
    Mint	    arg_number  = 7;
    Mint	    a_col_dim   = n;
    Mint	    h_col_dim   = n;
    Mfloat	    *d          = NULL;
    Mfloat	    **pd        = NULL;
    Mint	    user_d      = 0;
    Mint	    user_pd     = 0;
    Mfloat	    *diag       = NULL;
    Mfloat	    *obj        = NULL;
    Mint	    user_obj    = 0;
    Mfloat	    *at  	= NULL;
    Mfloat	    *work	= NULL;
    Mint	    *iwork      = NULL;
    Mint	    nact;
    Mint	    i, k;
    Mint	    return_user = 0;
    Mfloat	    tmp;

    diag = &tmp;
    code = 1;
    while (code > 0) {
           code = va_arg(argptr, Mint);
           arg_number++;
           switch (code) {
              case IMSL_RETURN_USER:
                  lv_x = va_arg(argptr, Mfloat*);
                  arg_number++;
                  return_user = 1;
                  break;
              case IMSL_DUAL_USER:
                  user_d = 1;
                  d = va_arg(argptr, Mfloat*);
                  arg_number++;
                  break;
              case IMSL_DUAL:
                  user_pd = 1;
                  pd = va_arg(argptr, Mfloat**);
                  arg_number++;
                  break;
              case IMSL_A_COL_DIM:
                  a_col_dim = va_arg(argptr, Mint);
                  arg_number++;
                  break;
              case IMSL_H_COL_DIM:
                  h_col_dim = va_arg(argptr, Mint);
                  arg_number++;
                  break;
              case IMSL_ADD_TO_DIAG_H:
                  diag = va_arg(argptr, Mfloat*);
                  arg_number++;
                  break;
              case IMSL_OBJ:
                  user_obj = 1;
                  obj = va_arg(argptr, Mfloat*);
                  arg_number++;
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

    if (imsl_n1rty(0)) goto RETURN;

    if (n <= 0) {
        imsl_e1sti(1, n);
        imsl_ermes(IMSL_TERMINAL, IMSL_N_MUST_BE_POSITIVE);

    } else if (n > a_col_dim) {
        imsl_e1sti(1, n);
        imsl_e1sti(2, a_col_dim);
        imsl_e1stl(1, "a");
        /* ermes(5, "The order of the matrix must be <= its column dim.") */
        imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
    }
    if (m <= 0) {
        imsl_e1sti(1, m);
        /* ermes(5, "The number of constraints, must be at least zero") */
        imsl_ermes(IMSL_TERMINAL, IMSL_LINEAR_CONSTRAINT_VALUE);
    }
    if (meq > m) {
        imsl_e1sti(1, meq);
        imsl_e1sti(2, m);
        imsl_ermes(IMSL_TERMINAL, IMSL_NEQ_CANNOT_BE_GT_NCON);
    }
    if (imsl_n1rty(0)) goto RETURN;

    if (!user_d)   d = (Mfloat *) imsl_malloc (n*sizeof(*d));

    at      = (Mfloat *) imsl_malloc (m*n*sizeof(*at));
    work    = (Mfloat *) imsl_malloc ((n*(n*3+11)/2+m)*sizeof(*work));
    iwork   = (Mint *) imsl_malloc (n*sizeof(*iwork));

    if (at==NULL || work==NULL || iwork==NULL || d==NULL) {
        imsl_e1stl(1, "n");
        imsl_e1sti(1, n);
        imsl_e1stl(2, "m");
        imsl_e1sti(2, m);
        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
        goto FREE_SPACE;
    }
    if (lv_x == NULL) {
	lv_x = (Mfloat *)imsl_malloc (n*sizeof(*lv_x));
	if (lv_x == NULL) {
	    imsl_e1stl(1, "n");
	    imsl_e1sti(1, n);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}
    }
    for (i = 1; i <= n; i++) {
         k = (i - 1) * m;
         scopy (m, &a[i - 1], a_col_dim, &at[k], 1);
    }

    imsl_q2rog(n, m, meq, at, m, b, g, h, h_col_dim, diag, lv_x, &nact,
            iwork, d, work);

    if (user_obj) {
        tmp = F_ZERO;
        for (i = 1; i <= n; i++) {
             k = (i - 1) * n;
             tmp += lv_x[i - 1] * imsl_sdot(n, &h[k], 1, lv_x, 1);
        }
        *obj = imsl_sdot(n, g, 1, lv_x, 1) + F_HALF * tmp;
    }

    if (user_pd) *pd = d;

FREE_SPACE:
    if (!user_d && !user_pd && d != NULL)
        imsl_free(d);
    else {
        sset(m,F_ZERO,at,1);
        for (i = 1; i <= nact; i++) {
             k = iwork[i - 1];
             at[k - 1] = d[i - 1];
        }
        scopy (m, at, 1, d, 1);
    }

    if (at != NULL)     imsl_free(at);
    if (work != NULL)   imsl_free(work);
    if (iwork != NULL)  imsl_free(iwork);

RETURN:
    if (imsl_n1rty(0) > 3) {
        if (!return_user && lv_x != NULL)
            imsl_free(lv_x);
        lv_x = NULL;
    }
    return (argptr);
}

/* -----------------------------------------------------------------------
    IMSL Name:  Q2ROG/DQ2ROG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 15, 1985

    Purpose:    Solve a quadratic programming problem subject to linear
                equality/inequality constraints.

    Usage:      CALL Q2ROG (NVAR, NCON, NEQ, A, LDA, B, GRAD, H, LDH,
                            DIAG, SOL, NACT, IACT, ALAMDA, WK)

    Arguments:  (See QPROG)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_q2rog(Mint nvar, Mint ncon, Mint neq, Mfloat *a,
                    Mint lda, Mfloat *b, Mfloat *grad, Mfloat *h,
                    Mint ldh, Mfloat *diag, Mfloat *sol, Mint *nact,
                    Mint *iact, Mfloat *alamda, Mfloat *wk)
#else
void imsl_q2rog(nvar, ncon, neq, a, lda, b, grad, h, ldh,
                    diag, sol, nact, iact, alamda, wk)
	Mint           nvar, ncon, neq;
	Mfloat         *a;
	Mint           lda;
	Mfloat          b[], grad[], *h;
	Mint           ldh;
	Mfloat         *diag, sol[];
	Mint           *nact, iact[];
	Mfloat          alamda[], wk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
#define H(I_,J_)	(h+(I_)*(ldh)+(J_))
	Mint            _d_l, _d_m, _do0, _do1, i, ia, id, ifinc,
	                iflag, ii, il, info, ir, ira, irb, is,
	                iterc, itref, iu, iw, iwa, iwd, iwr, iws, iww,
	                iwx, iwy, iwz, ix, iz, iza, j, jfinc, jflag, jl,
	                ju, k, kdrop, kfinc, kflag, kk, knext, lflag, m,
	                meq, mflag, n, nflag, nm, nu;
	Mfloat          cvmax, diagr, fdiff, fdiffa, ga, gb,
	                parinc, parnew, ratio, res, step, sum, suma, sumb,
	                sumc, sumx, sumy, temp, tempa, vsmall, xmag, xmagr,
	                zero;


	imsl_e1psh("Q2ROG ");
	zero = F_ZERO;
	n = nvar;
	m = ncon;
	meq = neq;
	iwz = n;
	iwr = iwz + n * n;
	iww = iwr + (n * (n + 3)) / 2;
	iwd = iww + n;
	iwx = iwd + n;
	iwa = iwx + n;
	/*
	 * Set some parameters. Numbers less than VSMALL are assumed to be
	 * negligible. The multiple of I that is addid to H is at most DIAGR
	 * times the least multiple of I that gives positive definiteness. X
	 * is re-initialized if its magnitude is reduced by the factor XMAGR.
	 * A check is made for an increase in F every IFINC iterations, after
	 * KFINC iterations are completed.
	 */
	vsmall = imsl_amach(4);
	diagr = F_TWO;
	xmagr = 0.01;
	ifinc = 3;
	kfinc = imsl_i_max(10, n);
	/*
	 * Find the reciprocals of the lengths of the constraint normals.
	 * Return id a constraint is infeasible due to a zero normal.
	 */
	*nact = 0;
	if (m <= 0)
		goto L_40;
	for (k = 1; k <= m; k++) {
		sum = F_ZERO;
		sum += imsl_sdot(n, A(0, k - 1), lda, A(0, k - 1), lda);
		if (sum > F_ZERO)
			goto L_10;
		if (b[k - 1] == F_ZERO)
			goto L_20;
		info = -k;
		if (k <= meq)
			goto L_9000;
		if (b[k - 1] > 0)
			goto L_9000;
		goto L_20;
L_10:
		sum = F_ONE / sqrt(sum);
L_20:
		ia = iwa + k;
		wk[ia - 1] = sum;
	}
	/*
	 * If necessary increase the DIAGonal elements of H.
	 */
L_40:
	*diag = F_ZERO;
	for (i = 1; i <= n; i++) {
		id = iwd + i;
		wk[id - 1] = *H(i - 1, i - 1);
		*diag = imsl_f_max(*diag, vsmall - wk[id - 1]);
		if (i == n)
			goto L_60;
		ii = i + 1;
		for (j = ii; j <= n; j++) {
			ga = -imsl_f_min(wk[id - 1], *H(j - 1, j - 1));
			gb = fabs(wk[id - 1] - *H(j - 1, j - 1)) + fabs(*H(j - 1, i - 1));
			if (gb > F_ZERO)
				ga += imsl_fi_power(*H(j - 1, i - 1), 2) / gb;
			*diag = imsl_f_max(*diag, ga);
		}
L_60:
		;
	}
	if (*diag <= F_ZERO)
		goto L_90;
L_70:
	*diag *= diagr;
	for (i = 1; i <= n; i++) {
		id = iwd + i;
		*H(i - 1, i - 1) = *diag + wk[id - 1];
	}
	/*
	 * Form the cholesky factorization of H. The transpose of the factor
	 * will be placed in the R-partition of W.
	 */
L_90:
	ir = iwr;
	for (j = 1; j <= n; j++) {
		ira = iwr;
		irb = ir + 1;
		for (i = 1; i <= j; i++) {
			temp = *H(j - 1, i - 1);
			if (i != 1) {
				temp -= imsl_sdot(ir - irb + 1, &wk[irb - 1], 1, &wk[ira],
						  1);
				ira += ir - irb + 1;
			}
			ir += 1;
			ira += 1;
			if (i < j)
				wk[ir - 1] = temp / wk[ira - 1];
		}
		if (temp < vsmall)
			goto L_130;
		wk[ir - 1] = sqrt(temp);
	}
	goto L_160;
	/*
	 * Increase further the DIAGonal elements of H.
	 */
L_130:
	wk[j - 1] = F_ONE;
	sumx = F_ONE;
	k = j;
L_140:
	sum = F_ZERO;
	ira = ir - 1;
	for (i = k; i <= j; i++) {
		sum += -wk[ira - 1] * wk[i - 1];
		ira += i;
	}
	ir -= k;
	k -= 1;
	wk[k - 1] = sum / wk[ir - 1];
	sumx += imsl_fi_power(wk[k - 1], 2);
	if (k >= 2)
		goto L_140;
	*diag += vsmall - temp / sumx;
	goto L_70;
	/*
	 * Set Z to the inverse of the matrix in R.
	 */
L_160:
	nm = n - 1;
	for (i = 1; i <= n; i++) {
		iz = iwz + i;
		sset(i - 1, F_ZERO, &wk[iz - 1], n);
		iz += (i - 1) * n;
		ir = iwr + (i + i * i) / 2;
		wk[iz - 1] = F_ONE / wk[ir - 1];
		if (i == n)
			goto L_190;
		iza = iz;
		for (j = i; j <= nm; j++) {
			ir += i;
			sum = F_ZERO;
			for (k = iza, _do0 = DOCNT(iza, iz, _do1 = n); _do0 > 0; k += _do1, _do0--) {
				sum += wk[k - 1] * wk[ir - 1];
				ir += 1;
			}
			iz += n;
			wk[iz - 1] = -sum / wk[ir - 1];
		}
L_190:
		;
	}
	/*
	 * Set the initial values of some variables. ITERC counts the number
	 * of iterations. ITREF is set to one when iterative refinement is
	 * required. JFINC indicates when to test for an increase in F.
	 */
	iterc = 1;
	itref = 0;
	jfinc = -kfinc;
	/*
	 * Set X to zero and set the corresponding residuals of the
	 * KUHN-TUCKER conditions.
	 */
L_200:
	iflag = 1;
	iws = iww - n;
	for (i = 1; i <= n; i++) {
		sol[i - 1] = F_ZERO;
		iw = iww + i;
		wk[iw - 1] = grad[i - 1];
		if (i > *nact)
			goto L_210;
		wk[i - 1] = F_ZERO;
		is = iws + i;
		k = iact[i - 1];
		wk[is - 1] = b[k - 1];
L_210:
		;
	}
	xmag = F_ZERO;
	if (*nact > 0)
		goto L_260;
	goto L_290;
	/*
	 * Set the residuals of the KUHN-TUCKER conditions for general X.
	 */
L_220:
	iflag = 2;
	iws = iww - n;
	for (i = 1; i <= n; i++) {
		iw = iww + i;
		wk[iw - 1] = grad[i - 1];
		wk[iw - 1] += imsl_sdot(n, H(0, i - 1), ldh, &sol[0], 1);
	}
	if (*nact == 0)
		goto L_290;
	for (k = 1; k <= *nact; k++) {
		kk = iact[k - 1];
		is = iws + k;
		wk[is - 1] = b[kk - 1];
		for (i = 1; i <= n; i++) {
			iw = iww + i;
			wk[iw - 1] += -wk[k - 1] ** A(i - 1, kk - 1);
			wk[is - 1] += -sol[i - 1] ** A(i - 1, kk - 1);
		}
	}
	/*
	 * Pre-multiply the vector in the S-partition of W by the inverse of
	 * R transpose.
	 */
L_260:
	ir = iwr;
	il = iws + 1;
	iu = iws + *nact;
	for (i = il; i <= iu; i++) {
		sum = F_ZERO;
		ju = i - 1;
		sum += imsl_sdot(ju - il + 1, &wk[ir], 1, &wk[il - 1], 1);
		ir += (ju - il + 1) + 1;
		wk[i - 1] = (wk[i - 1] - sum) / wk[ir - 1];
	}
	/*
	 * Shift X to satisfy the active constraints and make the
	 * corresponding change to the gradient residuals.
	 */
	for (i = 1; i <= n; i++) {
		iz = iwz + i;
		sum = F_ZERO;
		sum += imsl_sdot(iu - il + 1, &wk[il - 1], 1, &wk[iz - 1], n);
		sol[i - 1] += sum;
		saxpy(n, sum, H(0, i - 1), ldh, &wk[iww], 1);
	}
	/*
	 * Form the scalar product of the current gradient residuals with
	 * each column of Z.
	 */
L_290:
	kflag = 1;
	goto L_850;
L_300:
	if (*nact == n)
		goto L_320;
	/*
	 * Shift X so that it satisfies the remaining KUHN-TUCKER conditions.
	 */
	il = iws + *nact + 1;
	iza = iwz + *nact * n;
	for (i = 1; i <= n; i++) {
		sum = F_ZERO;
		iz = iza + i;
		sum += imsl_sdot(iww - il + 1, &wk[iz - 1], n, &wk[il - 1], 1);
		sol[i - 1] -= sum;
	}
	info = iterc;
	if (m == 0)
		goto L_640;
	if (*nact == 0)
		goto L_350;
	/* Update the LAGRANGE multipliers. */
L_320:
	lflag = 3;
	goto L_650;
L_330:
	for (k = 1; k <= *nact; k++) {
		iw = iww + k;
		wk[k - 1] += wk[iw - 1];
	}
	/*
	 * Revise the value of XMAG. Branch if iterative refinement is
	 * required.
	 */
L_350:
	jflag = 1;
	goto L_830;
L_360:
	if (iflag == itref)
		goto L_220;
	/*
	 * Delete a constraint if a LAGRANGE multiplier of an inequality
	 * constraint is negative.
	 */
	kdrop = 0;
	goto L_380;
L_370:
	kdrop += 1;
	if (wk[kdrop - 1] >= F_ZERO)
		goto L_380;
	if (iact[kdrop - 1] <= meq)
		goto L_380;
	nu = *nact;
	mflag = 1;
	goto L_720;
L_380:
	if (kdrop < *nact)
		goto L_370;
	/*
	 * Seek the greatest normalized constraint violation, disregarding
	 * any that may be due to computer rounding errors.
	 */
L_390:
	cvmax = F_ZERO;
	for (k = 1; k <= m; k++) {
		ia = iwa + k;
		if (wk[ia - 1] <= F_ZERO)
			goto L_400;
		sum = -b[k - 1] + imsl_sdot(n, &sol[0], 1, A(0, k - 1), lda);
		sumx = -sum * wk[ia - 1];
		if (k <= meq)
			sumx = fabs(sumx);
		if (sumx <= cvmax)
			goto L_400;
		temp = fabs(b[k - 1]) + l_a1ot(n, &sol[0], 1, A(0, k - 1),
						  lda);
		tempa = temp + fabs(sum);
		if (tempa <= temp)
			goto L_400;
		temp += 1.5 * fabs(sum);
		if (temp <= tempa)
			goto L_400;
		cvmax = sumx;
		res = sum;
		knext = k;
L_400:
		;
	}
	/* Test for convergence. */
	info = iterc;
	if (cvmax <= vsmall)
		goto L_630;
	/*
	 * Return if, due to rounding errors, the actual change in X may not
	 * increase the objective function.
	 */
	jfinc += 1;
	if (jfinc == 0)
		goto L_430;
	if (jfinc != ifinc)
		goto L_440;
	fdiff = F_ZERO;
	fdiffa = F_ZERO;
	for (i = 1; i <= n; i++) {
		sum = F_TWO * grad[i - 1];
		sumx = fabs(sum);
		for (j = 1; j <= n; j++) {
			ix = iwx + j;
			temp = *H(j - 1, i - 1) * (wk[ix - 1] + sol[j - 1]);
			sum += temp;
			sumx += fabs(temp);
		}
		ix = iwx + i;
		fdiff += sum * (sol[i - 1] - wk[ix - 1]);
		fdiffa += sumx * fabs(sol[i - 1] - wk[ix - 1]);
	}
	info = 0;
	sum = fdiffa + fdiff;
	if (sum <= fdiffa)
		goto L_630;
	temp = fdiffa + 1.5 * fdiff;
	if (temp <= sum)
		goto L_630;
	jfinc = 0;
L_430:
	scopy(n, &sol[0], 1, &wk[iwx], 1);
	/*
	 * Form the scalar product of the new constraint normal with each
	 * column of Z. PARNEW will become the LAGRANGE multiplier of the new
	 * constraint.
	 */
L_440:
	iterc += 1;
	scopy(n, A(0, knext - 1), lda, &wk[iww], 1);
	iws = iwr + (*nact + *nact ** nact) / 2;
	kflag = 2;
	goto L_850;
L_450:
	parnew = F_ZERO;
	/*
	 * Apply GIVENS rotations to make the last (N-NACT-2) scalar products
	 * equal to zero.
	 */
	if (*nact == n)
		goto L_480;
	nu = n;
	nflag = 1;
	goto L_780;
	/*
	 * Branch if there is no need to delete a constraint.
	 */
L_460:
	is = iws + *nact;
	if (*nact == 0)
		goto L_580;
	iz = iwz + *nact * n;
	suma = F_ZERO;
	sumb = F_ZERO;
	sumc = F_ZERO;
	suma += imsl_sdot(n, &wk[iww], 1, &wk[iz], 1);
	sumb += l_a1ot(n, &wk[iww], 1, &wk[iz], 1);
	sumc += imsl_sdot(n, &wk[iz], 1, &wk[iz], 1);
	temp = sumb + 0.1 * fabs(suma);
	tempa = sumb + 0.2 * fabs(suma);
	if (temp <= sumb)
		goto L_480;
	if (tempa <= temp)
		goto L_480;
	ia = iwa + knext;
	sumc = sqrt(sumc) / wk[ia - 1];
	temp = sumc + 0.1 * fabs(suma);
	tempa = sumc + 0.2 * fabs(suma);
	if (temp <= sumc)
		goto L_470;
	if (tempa <= temp)
		goto L_470;
	goto L_580;
	/*
	 * Calculate the multipliers for the new constraint normal expressed
	 * in terms of the active constraint normals. Then work out which
	 * constraint to drop.
	 */
L_470:
	lflag = 4;
	goto L_650;
L_480:
	lflag = 1;
	goto L_650;
	/*
	 * Complete the test for linearly dependent constraints
	 */
L_490:
	for (i = 1; i <= n; i++) {
		suma = *A(i - 1, knext - 1);
		sumb = fabs(suma);
		for (k = 1; k <= *nact; k++) {
			kk = iact[k - 1];
			if (kk <= ncon)
				goto L_498;
			kk -= ncon;
			temp = zero;
			if (kk == i)
				temp = wk[iww + kk - 1];
			kk -= nvar;
			if (kk == i)
				temp = -wk[iww + kk - 1];
			goto L_499;
	L_498:
			;
			iw = iww + k;
			temp = wk[iw - 1] ** A(i - 1, kk - 1);
	L_499:
			;
			suma -= temp;
			sumb += fabs(temp);
		}
		temp = sumb + 0.1 * fabs(suma);
		tempa = sumb + 0.2 * fabs(suma);
		if (temp <= sumb)
			goto L_510;
		if (tempa <= temp)
			goto L_510;
		goto L_570;
L_510:
		;
	}
	lflag = 1;
	goto L_690;
	/*
	 * Branch if the constraints are inconsistent.
	 */
L_520:
	info = -knext;
	if (kdrop == 0)
		goto L_630;
	parinc = ratio;
	parnew = parinc;
	/*
	 * Revise the LAGRANGE multipliers of the active constraints.
	 */
L_530:
	for (k = 1; k <= *nact; k++) {
		iw = iww + k;
		wk[k - 1] += -parinc * wk[iw - 1];
		if (iact[k - 1] > meq)
			wk[k - 1] = imsl_f_max(F_ZERO, wk[k - 1]);
	}
	if (kdrop == 0)
		goto L_610;
	/*
	 * Delete the constraint to be dropped. Shift the vector of scalar
	 * products. Then, if appropriate, make one more scalar product zero.
	 */
	nu = *nact + 1;
	mflag = 2;
	goto L_720;
L_550:
	iws += -*nact - 1;
	nu = imsl_i_min(n, nu);
	for (i = 1; i <= nu; i++) {
		is = iws + i;
		j = is + *nact;
		wk[is - 1] = wk[j];
	}
	nflag = 2;
	goto L_780;
	/*
	 * Calculate the step to the violated constraint.
	 */
L_570:
	is = iws + *nact;
L_580:
	sumy = wk[is];
	step = -res / sumy;
	parinc = step / sumy;
	if (*nact == 0)
		goto L_600;
	/*
	 * Calculate the changes to the LAGRANGE multipliers, and reduce the
	 * step along the new search direction if necessary.
	 */
	lflag = 2;
	goto L_650;
L_590:
	if (kdrop == 0)
		goto L_600;
	temp = F_ONE - ratio / parinc;
	if (temp <= F_ZERO)
		kdrop = 0;
	if (kdrop == 0)
		goto L_600;
	step = ratio * sumy;
	parinc = ratio;
	res *= temp;
	/*
	 * Update X and the LAGRANGE multipliers. Drop a constraint if the
	 * full step is not taken.
	 */
L_600:
	iwy = iwz + *nact * n;
	/*
	 * DO 760  I=1, N IY     = IWY + I SOL(I) = SOL(I) + STEP*WK(IY) 760
	 * CONTINUE
	 */
	saxpy(n, step, &wk[iwy], 1, &sol[0], 1);
	parnew += parinc;
	if (*nact >= 1)
		goto L_530;
	/*
	 * Add the new constraint to the active set.
	 */
L_610:
	*nact += 1;
	wk[*nact - 1] = parnew;
	iact[*nact - 1] = knext;
	ia = iwa + knext;
	wk[ia - 1] = -wk[ia - 1];
	/*
	 * Estimate the magnitude of X. Then begin a new iteration,
	 * re-initializing X if this magnitude is small.
	 */
	jflag = 2;
	goto L_830;
L_620:
	if (sum < xmagr * xmag)
		goto L_200;
	if (itref > 0)
		goto L_220;
	goto L_390;
	/*
	 * Initiate iterative refinement if it has not yet been used, or
	 * return after restoring the DIAGonal elements of H.
	 */
L_630:
	itref += 1;
	jfinc = -1;
	if (itref == 1)
		goto L_220;
L_640:
	scopy(n, &wk[iwd], 1, H(0, 0), ldh + 1);
	/*
	 * 800 DO 810  I=1, N ID     = IWD + I H(I,I) = WK(ID) 810 CONTINUE
	 * Due to rounding errors, the actual change in the variable may not
	 * increase the objective function.
	 */
	if (info == 0) {
               /* imsl_ermes(3, "Due to the effect of computer rounding       */
               /*                error, a change in the variables fail to     */
               /*                improve the objective function value.        */
               /*                Usually the solution is close to optimum."); */
                imsl_ermes(IMSL_WARNING, IMSL_NO_MORE_PROGRESS);

	} else if (info < 0) {
		/* The constraint with index ABS(INFO) and the constraints  */
		/* whose indices are IACT(K),K=1,.,NACT, are inconsistent.  */

                /* imsl_ermes(4, "The system of equations is inconsistent. */
                /*                There is no solution.");                 */
                imsl_ermes(IMSL_FATAL, IMSL_SYSTEM_INCONSISTENT);

	} else {
                /* The calculation is successful */
		scopy(*nact, wk, 1, alamda, 1);
		sset(nvar - *nact, F_ZERO, &alamda[*nact], 1);
	}
L_9000:
	imsl_e1pop("Q2ROG ");
	return;
	/*
	 * The remaining instructions are used as subroutines. Calculate the
	 * LAGRANGE multipliers by pre-multiplying the vector in the
	 * S-partition of W by the inverse of R.
	 */
L_650:
	ir = iwr + (*nact + *nact ** nact) / 2;
	i = *nact;
	sum = F_ZERO;
	goto L_680;
L_660:
	ira = ir - 1;
	sum = F_ZERO;
	for (j = i; j <= *nact; j++) {
		iw = iww + j;
		sum += wk[ira - 1] * wk[iw - 1];
		ira += j;
	}
	ir -= i;
	i -= 1;
L_680:
	iw = iww + i;
	is = iws + i;
	wk[iw - 1] = (wk[is - 1] - sum) / wk[ir - 1];
	if (i > 1)
		goto L_660;
	if (lflag == 3)
		goto L_330;
	if (lflag == 4)
		goto L_490;
	/*
	 * Calculate the next constraint to drop.
	 */
L_690:
	kdrop = 0;
	for (k = 1; k <= *nact; k++) {
		if (iact[k - 1] <= meq)
			goto L_710;
		iw = iww + k;
		if (res * wk[iw - 1] >= F_ZERO)
			goto L_710;
		temp = wk[k - 1] / wk[iw - 1];
		if (kdrop == 0)
			goto L_700;
		if (fabs(temp) >= fabs(ratio))
			goto L_710;
L_700:
		kdrop = k;
		ratio = temp;
L_710:
		;
	}
	if (lflag == 1)
		goto L_520;
	if (lflag == 2)
		goto L_590;
	/*
	 * Drop the constraint in position KDROP in the active set.
	 */
L_720:
	ia = iwa + iact[kdrop - 1];
	wk[ia - 1] = -wk[ia - 1];
	if (kdrop == *nact)
		goto L_770;
	/*
	 * Set some indices and calculate the elements of the next GIVENS
	 * rotation.
	 */
	iz = iwz + kdrop * n;
	ir = iwr + (kdrop + kdrop * kdrop) / 2;
L_730:
	ira = ir;
	ir += kdrop + 1;
	temp = imsl_f_max(fabs(wk[ir - 2]), fabs(wk[ir - 1]));
	sum = temp * sqrt(imsl_fi_power(wk[ir - 2] / temp, 2) + imsl_fi_power(wk[ir - 1] / temp, 2));
	ga = wk[ir - 2] / sum;
	gb = wk[ir - 1] / sum;
	/* Exchange the columns or R. */
	for (i = 1; i <= kdrop; i++) {
		ira += 1;
		j = ira - kdrop;
		temp = wk[ira - 1];
		wk[ira - 1] = wk[j - 1];
		wk[j - 1] = temp;
	}
	wk[ir - 1] = F_ZERO;
	/* Apply the rotation to the rows of R. */
	wk[j - 1] = sum;
	kdrop += 1;
	for (i = kdrop; i <= nu; i++) {
		temp = ga * wk[ira - 1] + gb * wk[ira];
		wk[ira] = ga * wk[ira] - gb * wk[ira - 1];
		wk[ira - 1] = temp;
		ira += i;
	}
	/*
	 * Apply the rotation to the columns of Z.
	 */
	for (i = 1; i <= n; i++) {
		iz += 1;
		j = iz - n;
		temp = ga * wk[j - 1] + gb * wk[iz - 1];
		wk[iz - 1] = ga * wk[iz - 1] - gb * wk[j - 1];
		wk[j - 1] = temp;
	}
	/*
	 * Revise IACT and the LAGRANGE multipliers.
	 */
	iact[kdrop - 2] = iact[kdrop - 1];
	wk[kdrop - 2] = wk[kdrop - 1];
	if (kdrop < *nact)
		goto L_730;
L_770:
	*nact -= 1;
	if (mflag == 1)
		goto L_220;
	if (mflag == 2)
		goto L_550;
	/*
	 * Apply GIVENS rotations to reduce some of the scalar products in
	 * the S-partition of W to zero.
	 */
L_780:
	iz = iwz + nu * n;
L_790:
	iz -= n;
L_800:
	is = iws + nu;
	nu -= 1;
	if (nu == *nact)
		goto L_820;
	if (wk[is - 1] == F_ZERO)
		goto L_790;
	temp = imsl_f_max(fabs(wk[is - 2]), fabs(wk[is - 1]));
	sum = temp * sqrt(imsl_fi_power(wk[is - 2] / temp, 2) + imsl_fi_power(wk[is - 1] / temp, 2));
	ga = wk[is - 2] / sum;
	gb = wk[is - 1] / sum;
	wk[is - 2] = sum;
	for (i = 1; i <= n; i++) {
		k = iz + n;
		temp = ga * wk[iz - 1] + gb * wk[k - 1];
		wk[k - 1] = ga * wk[k - 1] - gb * wk[iz - 1];
		wk[iz - 1] = temp;
		iz -= 1;
	}
	goto L_800;
L_820:
	if (nflag == 1)
		goto L_460;
	if (nflag == 2)
		goto L_570;
	/*
	 * Calculate the magnitude of X and revise XMAG.
	 */
L_830:
	sum = F_ZERO;
	for (i = 1; i <= n; i++) {
		sum += fabs(sol[i - 1]) * (fabs(grad[i - 1]) + fabs(*H(i - 1, i - 1) *
							       sol[i - 1]));
	}
	xmag = imsl_f_max(xmag, sum);
	if (jflag == 1)
		goto L_360;
	if (jflag == 2)
		goto L_620;
	/*
	 * Pre-multiply the vector in the W-partition of W by Z transpose.
	 */
L_850:
	jl = iww + 1;
	iz = iwz;
	for (i = 1; i <= n; i++) {
		is = iws + i;
		wk[is - 1] = F_ZERO;
		wk[is - 1] += imsl_sdot(iwd - jl + 1, &wk[iz], 1, &wk[jl - 1],
					1);
		iz += iwd - jl + 1;
	}
	if (kflag == 1)
		goto L_300;
	if (kflag == 2)
		goto L_450;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/23/90 at 08:25:32 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/23/90 at 08:25:31
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  A1OT/DA1OT (Single/Double precision version)

    Computer:   $COMPUTER/$PRECISION

    Revised:    August 9, 1986

    Purpose:    Compute sum of absolute values of products.

    Usage:      A1OT(N, SX, INCX, SY, INCY)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be.. SX(1+(I-1)*INCX) if INCX .GE. 0
                or SX(1+(I-N)*INCX) if INCX .LT. 0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).  (Input)
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be.. SY(1+(I-1)*INCY) if INCY .GE. 0
                or SY(1+(I-N)*INCY) if INCY .LT. 0.
       A1OT   - Sum from I=1 to N of ABS(X(I)*Y(I)).  (Output)
                X(I) and Y(I) refer to specific elements of SX and SY,
                respectively.  See INCX and INCY argument descriptions.

    GAMS:       D1a4

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_a1ot(Mint n, Mfloat *sx, Mint incx, Mfloat *sy,
                    Mint incy)
#else
static Mfloat l_a1ot(n, sx, incx, sy, incy)
	Mint            n;
	Mfloat          sx[];
	Mint            incx;
	Mfloat          sy[];
	Mint            incy;
#endif
{
	Mint            i, ix, iy;
	Mfloat          a1ot_v;


	a1ot_v = F_ZERO;
	if (n > 0) {
		if (incx != 1 || incy != 1) {
			/* CODE FOR UNEQUAL INCREMENTS */
			ix = 1;
			iy = 1;
			if (incx < 0)
				ix = (-n + 1) * incx + 1;
			if (incy < 0)
				iy = (-n + 1) * incy + 1;
			for (i = 1; i <= n; i++) {
				a1ot_v += fabs(sx[ix - 1] * sy[iy - 1]);
				ix += incx;
				iy += incy;
			}
		} else {
			for (i = 1; i <= n; i++) {
				a1ot_v += fabs(sx[i - 1] * sy[i - 1]);
			}
		}
	}
	return (a1ot_v);
}				/* end of function */
