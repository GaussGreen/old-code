#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_covariances (Mint, Mint, Mfloat*, va_list argptr);
static void l_c2rvc(Mint, Mint, Mint, Mfloat*, Mint, Mint, Mint, Mint, Mint,
                    Mfloat*, Mfloat*, Mint, Mint*, Mint, Mint*, Mint*, Mfloat*,
                    Mfloat[]);
static Mint l_i1nan(Mint , Mfloat *, Mint);
#else
static VA_LIST_HACK l_covariances();
static void l_c2rvc();
static Mint l_i1nan();
#endif

static Mfloat *lv_covariances = NULL;

#ifdef ANSI
Mfloat *imsl_f_covariances (Mint n_observations, Mint n_variables, Mfloat *x, ...)
#else
Mfloat *imsl_f_covariances (n_observations, n_variables, x, va_alist)
    Mint        n_observations, n_variables;
    Mfloat      *x;
    va_dcl
#endif
{
    va_list     argptr;
    VA_START(argptr, x);
    E1PSH("imsl_f_covariances", "imsl_d_covariances");
    lv_covariances = NULL;
    IMSL_CALL(l_covariances (n_observations, n_variables, x, argptr));
    va_end(argptr);
    E1POP("imsl_f_covariances", "imsl_d_covariances");
    return lv_covariances;
}

#ifdef ANSI
static VA_LIST_HACK l_covariances (Mint n_observations, Mint n_variables, Mfloat *x, va_list argptr)
#else
static VA_LIST_HACK l_covariances (n_observations, n_variables, x, argptr)
    Mint        n_observations, n_variables;
    Mfloat      *x;
    va_list     argptr;
#endif
{
    Mint        col_dim_x               = n_variables;
    Mint        col_dim_covariances     = n_variables;
    Mint	compute_option		= 0;
    Mfloat    **means = NULL;
    Mfloat     *means_ptr = NULL;

    Mint         code = 1, icmpt = 0;
    Mint   return_means = 0, user_means = 0;
    Mint   n_arguements = 3, user_covariances = 0;
    Mfloat      *wk, sumwt;
    Mint	incd[1], nobs, nmiss;
 
    while (code>0){
        code = va_arg(argptr, Mint);
        ++n_arguements;
        switch (code){
            case IMSL_RETURN_USER:
                lv_covariances = va_arg(argptr, Mfloat*);
                user_covariances = 1;
                ++n_arguements;
                break;
            case IMSL_X_COL_DIM:
                col_dim_x = va_arg(argptr, Mint);
                ++n_arguements;
                break;
            case IMSL_COVARIANCE_COL_DIM:
                col_dim_covariances = va_arg(argptr, Mint);
                ++n_arguements;
                break;
	    case IMSL_VARIANCE_COVARIANCE_MATRIX:
                if (!icmpt) {
                   compute_option = 0;
                   icmpt = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_COV_OPTIONS);
                }
	    case IMSL_CORRECTED_SSCP_MATRIX:
                if (!icmpt) {
                   compute_option = 1;
                   icmpt = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_COV_OPTIONS);
                   goto RETURN;
                }
	    case IMSL_CORRELATION_MATRIX:
                if (!icmpt) {
                   compute_option = 2;
                   icmpt = 2;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_COV_OPTIONS);
                   goto RETURN;
                }
	    case IMSL_STDEV_CORRELATION_MATRIX:
                if (!icmpt) {
                   compute_option = 3;
                   icmpt = 3;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_COV_OPTIONS);
                   goto RETURN;
                }
	    case IMSL_MEANS:
		means        = va_arg(argptr, Mfloat**);
		user_means   = 0;
		return_means = 1;
		break;
	    case IMSL_MEANS_USER:
		means_ptr    = va_arg(argptr, Mfloat*);
		user_means   = 1;
		return_means = 1;
		break;
            case 0:
                break;
            default:
                imsl_e1sti (1, code);
                imsl_e1sti (2, n_arguements);
                imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
                break;
        
	}
    }
    if (imsl_n1rty(0) > 0) goto RETURN;
    
    if (x==NULL) {
	imsl_e1stl (1, "x");
        imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
    }

    if (col_dim_covariances < n_variables){
        imsl_e1stl(1, "covariance");
        imsl_e1sti(1, n_variables);
        imsl_e1sti(2, col_dim_covariances);
        imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_LESS_COL);
    }

    if (col_dim_x < n_variables){
        imsl_e1stl(1, "x");
        imsl_e1sti(1, n_variables);
        imsl_e1sti(2, col_dim_x);
        imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_LESS_COL);
    }

    if (n_observations < 1){
        imsl_e1stl(1, "n_observations");
        imsl_e1sti(1, n_observations);
        imsl_e1stl(2, "1");
        imsl_ermes(IMSL_TERMINAL, IMSL_CHOOSE_S1_GREATER_S2);
    }

    if (n_variables < 1){
        imsl_e1stl(1, "n_variables");
        imsl_e1sti(1, n_variables);
        imsl_e1stl(2, "0");
        imsl_ermes(IMSL_TERMINAL, IMSL_CHOOSE_S1_GREATER_S2);
    }
    if (imsl_n1rty(0) > 0) goto RETURN;
    switch (compute_option){
	case 1:
	    wk = (Mfloat *) imsl_malloc (n_variables*(n_variables+2)*sizeof(Mfloat));
	    break;
	case 2:
	    wk = (Mfloat *) imsl_malloc (n_variables*(n_variables+2)*sizeof(Mfloat));
	    break;
	case 3:
	    wk = (Mfloat *) imsl_malloc (2*n_variables*(n_variables+1)*sizeof(Mfloat));
	    break;
	default: /* IMSL_VARIANCE_COVARIANCE_MATRIX */
	    wk = (Mfloat *) imsl_malloc (3*n_variables*sizeof(Mfloat));
	    break;
    }
    if (!wk){
        imsl_e1sti (1,n_variables);
        imsl_e1stl (1,"n_variables");
        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
        goto RETURN;
    }

    if (!user_covariances) {
        lv_covariances = (Mfloat *)imsl_malloc (n_variables*col_dim_covariances*sizeof(Mfloat));
        if (!lv_covariances) {
           imsl_e1sti (1,n_variables);
           imsl_e1sti (2,col_dim_covariances);
           imsl_e1stl (1,"n_variables");
           imsl_e1stl (2,"covariance_col_dim");
           imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
            goto FREE_SPACE;
        }
    }
    if (!user_means){
	means_ptr = (Mfloat *) imsl_malloc (n_variables*sizeof(Mfloat));
	if (!means_ptr){
            imsl_e1sti (1,n_variables);
            imsl_e1stl (1,"n_variables");
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    if (!user_covariances){
		if (lv_covariances) imsl_free (lv_covariances);
		lv_covariances = NULL;
	    }
            goto FREE_SPACE;	    
	}
    }

    imsl_f_m1ran(n_observations, col_dim_x, x, x);

    l_c2rvc(0, n_observations, n_variables, x, n_observations, 0, 0, 0, compute_option, means_ptr, lv_covariances, col_dim_covariances, incd, 1, &nobs, &nmiss, &sumwt, wk);

    imsl_f_m1ran(col_dim_x, n_observations, x, x);

    if ((code=imsl_n1rty(1))>3 && code!=6 &&!user_covariances){
        if (lv_covariances) imsl_free (lv_covariances);
	return_means = 0;
        lv_covariances = NULL;
    }
    else if (return_means && !user_means) *means = means_ptr;

FREE_SPACE:
    if (!return_means && !user_means && means_ptr) imsl_free(means_ptr);
    if (wk) imsl_free(wk);

RETURN:

    return argptr;
}

/*
  -----------------------------------------------------------------------
    IMSL Name:  C2RVC/DC2RVC (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 22, 1988

    Purpose:    Compute the means and the variance-covariance or
                correlation matrix.

    Usage:      CALL C2RVC (IDO, NROW, NVAR, X, LDX, IFRQ, IWT, MOPT,
                            ICOPT, XMEAN, COV, LDCOV, INCD, LDINCD, NOBS,
                            NMISS, SUMWT, WK)

    Arguments:  See CORVC

    Chapter:    STAT/LIBRARY Correlation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
static void l_c2rvc(ido, nrow, nvar, x, ldx, ifrq, iwt, mopt, icopt, xmean, cov, ldcov, incd, ldincd, nobs, nmiss, sumwt, wk)
	Mint            ido, nrow, nvar;
	Mfloat         *x;
	Mint            ldx, ifrq, iwt, mopt, icopt;
	Mfloat          xmean[], *cov;
	Mint            ldcov, *incd, ldincd, *nobs, *nmiss;
	Mfloat         *sumwt, wk[];
{
#define X(I_,J_)	(x+(I_)*(ldx)+(J_))
#define COV(I_,J_)	(cov+(I_)*(ldcov)+(J_))
#define INCD(I_,J_)	(incd+(I_)*(ldincd)+(J_))
	Mint            i, igo, imean, imiss, irow, ivar, j, jj, jmean, jvar, k, kk, kmean, kvar, l, l1, l2, ner, nrowp, num, nunits;
	static Mint     icall = 1, ipairm, ipairv, ipairw, ncol, nvar2, nvarsq;
	Mfloat          cnt, cntm1, covl, denom, rat, temp1, temp2, tj, tjk, tkj, wi, wts, xmiss, xnum, xwt;


	imsl_e1psh("l_c2rvc");
	ner = 1;
	xmiss = imsl_amach(6);
	imsl_c1iarg(ido, "ido", 0, 3, &ner);
	imsl_c1iarg(nvar, "nvar", 0, -1, &ner);
	imsl_c1iarg(ifrq, "ifrq", 0, -1, &ner);
	imsl_c1iarg(iwt, "iwt", 0, -1, &ner);
	imsl_c1iarg(mopt, "mopt", 0, 3, &ner);
	imsl_c1iarg(icopt, "icopt", 0, 3, &ner);
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	nrowp = (nrow > 0 ? nrow : -nrow);

	if (ido == 3 && nrow == 0)
		goto L_90;

	imsl_c1dim(0, nrowp, "abs(nrow)", ldx, "ldx", &ner);
	if (imsl_n1rcd(0) >= 4)
		goto L_9000;
	/* Get some constants */
	ncol = nvar;
	if (ifrq > 0)
		ncol += 1;
	if (iwt > 0)
		ncol += 1;
	nvar2 = 2 * nvar;
	nvarsq = nvar * nvar;
	ipairm = nvar2;
	ipairw = nvar2 + nvarsq;
	ipairv = nvar2 + nvarsq;
	if (mopt == 3 && iwt > F_ZERO)
		ipairv += (nvarsq + nvar) / 2;

	if (ido <= 1) {
		/*
		 * Initialize XMEAN, COV, INCD, WK and ICALL
		 */
		icall = 1;
		*INCD(0, 0) = 0;
		sset(nvar, F_ZERO, xmean, 1);
		for (j = 1; j <= nvar; j++) {
			sset(nvar, F_ZERO, COV(j - 1, 0), 1);
			if (mopt > 0)
				iset(nvar, 0, INCD(j - 1, 0), 1);
		}
		/* Determine length of workspace. */
		if (mopt < 1) {
			if (iwt <= 0) {
				nunits = 3 * nvar;
			} else {
				nunits = 4 * nvar;
			}
		} else if (mopt == 1 || mopt == 2) {
			if (iwt <= 0) {
				nunits = nvar * (nvar + 2);
			} else {
				nunits = (nvar * (3 * nvar + 5)) / 2;
			}
		} else if (mopt == 3) {
			if (iwt <= 0) {
				nunits = 2 * nvar * (nvar + 1);
			} else {
				nunits = (5 * nvar * (nvar + 1)) / 2;
			}
		}
		sset(nunits, F_ZERO, wk, 1);
		jj = 1;
		/* Initialize constant indicator */
		sset(nvar, imsl_amach(6), &wk[nvar], 1);
		*nobs = 0;
		*nmiss = 0;
		*sumwt = F_ZERO;
	} else {
		icall += 1;
	}

	if (nrow < 0) {
		irow = -1;
	} else {
		irow = 1;
	}
	/*
	 * Accumulate means and adjusted sum of squares and cross products.
	 * Test for constant observation set.
	 */
	xnum = F_ONE;
	for (i = 1; i <= nrowp; i++) {
		imsl_c1wfr(ido, icall, x, ldx, i, irow, ifrq, iwt, xmiss, nmiss,
		      &xnum, &wi, &igo);
		if (igo == 3)
			goto L_9000;
		if (igo == 1 || igo == 2)
			goto L_80;
		num = nint(xnum);
		xwt = wi * num;
		imiss = l_i1nan(ncol, X(0, i - 1), ldx);
		if (imiss != 0)
			*nmiss += 1;
		if (mopt < 1) {
			/*
			 * No missing values or else listwise deletion option
			 */
			if (imiss != 0)
				goto L_80;
			*sumwt += xwt;
			*nobs += num;
			if (*nobs < 0) {

/*				(4, 10, "More observations have been deleted than were originally entered.");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_TOO_MANY_OBS_DELETED);
				goto L_9000;
			}
			if (*sumwt == F_ZERO)
				goto L_80;
			rat = xwt / *sumwt;
			jj = 0;
			temp1 = xwt * (F_ONE - rat);
			for (j = 1; j <= nvar; j++) {
		L_20:
				jj += 1;
				if (jj == ifrq || jj == iwt)
					goto L_20;
				/* Check for constant */
				if (imsl_ifnan(wk[nvar + j - 1])) {
					if ((wi > F_ZERO) && (!imsl_ifnan(*X(jj - 1, i - 1))))
						wk[nvar + j - 1] = *X(jj - 1, i - 1);
				} else {
					if ((wi > F_ZERO) && (*X(jj - 1, i - 1) != wk[nvar + j - 1]))
						wk[j - 1] = F_ONE;
				}
				tj = *X(jj - 1, i - 1) - xmean[j - 1];
				wk[nvar2 + j - 1] = tj;
				temp2 = temp1 * tj;
				saxpy(j, temp2, &wk[nvar2], 1, COV(0, j - 1), ldcov);
				xmean[j - 1] += tj * rat;
			}
			*INCD(0, 0) = *nobs;
		} else {
			/* Pairwise deletion option */
			*sumwt += xwt;
			*nobs += num;
			if (*nobs < 0) {

/*				(4, 10, "More observations have been deleted than were originally entered.");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_TOO_MANY_OBS_DELETED);
				goto L_9000;
			}
			jmean = ipairm - nvar;
			jvar = ipairv - nvar;
			jj = 0;
			for (j = 1; j <= nvar; j++) {
		L_40:
				jj += 1;
				if (jj == ifrq || jj == iwt)
					goto L_40;
				jmean += nvar;
				jvar += nvar;
				if (imsl_ifnan(*X(jj - 1, i - 1)))
					goto L_70;
				if (imsl_ifnan(wk[nvar + j - 1])) {
					if ((wi > F_ZERO) && (!imsl_ifnan(*X(jj - 1, i - 1))))
						wk[nvar + j - 1] = *X(jj - 1, i - 1);
				} else {
					if ((wi > F_ZERO) && (*X(jj - 1, i - 1) != wk[nvar + j - 1]))
						wk[j - 1] = F_ONE;
				}
				kmean = ipairm - nvar;
				kvar = ipairv - nvar;
				kk = 0;
				l = j * (j - 1) / 2;
				for (k = 1; k <= j; k++) {
			L_50:
					kk += 1;
					if (kk == ifrq || kk == iwt)
						goto L_50;
					l += 1;
					kmean += nvar;
					kvar += nvar;
					if (!imsl_ifnan(*X(kk - 1, i - 1))) {
						*INCD(k - 1, j - 1) += num;
						if (*INCD(k - 1, j - 1) < 0) {
							imsl_e1sti(1, j);
							imsl_e1sti(2, k);

/*							(4, 11, "More observations are being deleted from COV(%(i1),%(I2)) than were originally entered.  (INCD(%(i1),%(i2)) is less than zero.)");
*/

imsl_ermes(IMSL_FATAL, IMSL_MORE_OBS_DEL_THAN_ENTERED);
							goto L_9000;
						}
						if (iwt > F_ZERO) {
							wk[ipairw + l - 1] += xwt;
							if (wk[ipairw + l - 1] == F_ZERO)
								goto L_60;
							rat = xwt / wk[ipairw + l - 1];
						} else {
							if (*INCD(k - 1, j - 1) == 0)
								goto L_60;
							rat = xwt / *INCD(k - 1, j - 1);
						}
						tjk = *X(kk - 1, i - 1) - wk[jmean + k - 1];
						tkj = *X(jj - 1, i - 1) - wk[kmean + j - 1];
						*COV(k - 1, j - 1) += xwt * tkj * tjk * (F_ONE - rat);
						wk[jmean + k - 1] += tjk * rat;
						if (mopt == 3)
							wk[jvar + k - 1] += xwt * tjk * tjk * (F_ONE -
								       rat);
						if (j != k) {
							wk[kmean + j - 1] += tkj * rat;
							if (mopt == 3)
								wk[kvar + j - 1] += xwt * tkj * tkj * (F_ONE -
								       rat);
						}
					}
			L_60:
					;
				}
				if (iwt > 0) {
					wts = wk[ipairw + j * (j + 1) / 2 - 1];
					if (wts == F_ZERO) {
						xmean[j - 1] = F_ZERO;
					} else {
						xmean[j - 1] += xwt * (*X(jj - 1, i - 1) - xmean[j - 1]) /
							wts;
					}
				} else {
					if (*INCD(j - 1, j - 1) == F_ZERO) {
						xmean[j - 1] = F_ZERO;
					} else {
						xmean[j - 1] += xwt * (*X(jj - 1, i - 1) - xmean[j - 1]) /
							*INCD(j - 1, j - 1);
					}
				}
		L_70:
				;
			}
		}
L_80:
		;
	}
L_90:
	if (ido == 0 || ido == 3) {
		/* Postprocessing begins here */
		for (i = 1; i <= nvar; i++) {
			if (*COV(i - 1, i - 1) < -100.0 * imsl_amach(4)) {
				imsl_e1sti(1, i);

/*				(4, 18, "Different observations are being deleted from COV(%(i1),%(i1)) then were originally entered.  (COV(%(i1),%(i1)) is less than zero.)");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_DIFFERENT_OBS_DELETED);
				goto L_9000;
			}
		}
		cntm1 = *nobs - 1;
		/* Sum of the weights is zero */
		if (*sumwt < F_TEN * imsl_amach(4)) {
			sset(nvar, xmiss, xmean, 1);
			if (icopt <= 1) {
				for (i = 1; i <= nvar; i++) {
					sset(nvar, xmiss, COV(i - 1, 0), 1);
				}

/*				(3, 12, "The sum of the weights is zero.  The means, variances and covariances are set to NaN (not a number).");
*/
                                imsl_ermes(IMSL_WARNING,
				IMSL_ZERO_SUM_OF_WEIGHTS_2);
			} else {
				for (i = 1; i <= nvar; i++) {
					sset(nvar, xmiss, COV(i - 1, 0), 1);
				}

/*				(3, 13, "The sum of the weights is zero.  The means and correlations are set to NaN (not a number).");
*/
                                imsl_ermes(IMSL_WARNING,
				IMSL_ZERO_SUM_OF_WEIGHTS_3);
			}
			goto L_9000;
		}
		/* Check for constant values */
		for (i = 1; i <= nvar; i++) {
			if (wk[i - 1] == F_ZERO) {
				*COV(i - 1, i - 1) = F_ZERO;
				if (icopt >= 2) {
					imsl_e1sti(1, i);

/*					(3, 14, "Correlations are requested but the observations on variable %(i1) are constant.  The pertinent correlation coefficients are set to NaN (not a number).");
*/
                                        imsl_ermes(IMSL_WARNING,
					IMSL_CONSTANT_VARIABLE);
				}
			}
		}
		if (mopt == 1) {
			/*
			 * Adjust sum of squares and cross-products matrix.
			 */
			l = 0;
			imean = ipairm - nvar;
			for (i = 1; i <= nvar; i++) {
				jmean = ipairm - nvar;
				imean += nvar;
				for (j = 1; j <= i; j++) {
					l += 1;
					jmean += nvar;
					if (iwt > 0) {
						cnt = wk[ipairw + l - 1];
					} else {
						cnt = *INCD(j - 1, i - 1);
					}
					*COV(j - 1, i - 1) += cnt * (xmean[i - 1] * xmean[j - 1] +
								     wk[jmean + i - 1] * wk[imean + j - 1] - xmean[i - 1] *
								     wk[imean + j - 1] - xmean[j - 1] * wk[jmean + i - 1]);
				}
				if (*COV(i - 1, i - 1) < F_ZERO)
					*COV(i - 1, i - 1) = F_ZERO;
			}
		}
		if (icopt == 0) {
			if (mopt == 0) {
				/*
				 * Compute the variance-covariance matrix
				 */
				if (*INCD(0, 0) >= 2) {
					temp1 = F_ONE / cntm1;
					for (i = 1; i <= nvar; i++) {
						sscal(i, temp1, COV(0, i - 1), ldcov);
					}
					goto L_270;
				}
/*			        (3, 15, "Variances and covariances are requested but fewer than two valid observations are present for a variable.  The pertinent statistics are set to NaN (not a number).");
*/
                                imsl_ermes(IMSL_WARNING,
				IMSL_INSUFFICIENT_DATA);
				for (i = 1; i <= nvar; i++) {
					sset(i, xmiss, COV(0, i - 1), ldcov);
				}
				goto L_270;
			}
			l = 0;
			for (i = 1; i <= nvar; i++) {
				for (j = 1; j <= i; j++) {
					if (*INCD(j - 1, i - 1) >= 2) {
						*COV(j - 1, i - 1) /= *INCD(j - 1, i - 1) -
							F_ONE;
					} else {
						*COV(j - 1, i - 1) = xmiss;
						l = 1;
					}
				}
			}
			if (l == 1) {

/*				(3, 15, "Variances and covariances are requested but fewer than two valid observations are present for a variable.  The pertinent statistics are set to NaN (not a number).");
*/
                                imsl_ermes(IMSL_WARNING,
				IMSL_INSUFFICIENT_DATA);
			}
			if (mopt == 3) {
				/*
				 * Change the pairwise sum of products to
				 * covariances
				 */
				l2 = -nvar;
				for (i = 1; i <= nvar; i++) {
					k = i - 1;
					l1 = i - nvar;
					l2 += -l1 + 1;
					for (j = 1; j <= k; j++) {
						l1 += nvar;
						l2 += 1;
						if (*INCD(j - 1, i - 1) >= 2) {
							wk[ipairv + l1 - 1] /= *INCD(j - 1, i - 1) -
								1;
							wk[ipairv + l2 - 1] /= *INCD(j - 1, i - 1) -
								1;
						} else {
							wk[ipairv + l1 - 1] = F_ZERO;
							wk[ipairv + l2 - 1] = F_ZERO;
						}
					}
					l2 += 1;
					if (*INCD(i - 1, i - 1) >= 2) {
						wk[ipairv + l2 - 1] /= *INCD(i - 1, i - 1) -
							1;
					} else {
						wk[ipairv + l2 - 1] = F_ZERO;
					}
				}
			}
		} else if (icopt >= 2) {
			/* Compute the correlation matrix */
			if (nvar > 1) {
				jvar = ipairv;
				for (j = 2; j <= nvar; j++) {
					jvar += nvar;
					k = j - 1;
					ivar = ipairv - nvar;
					for (i = 1; i <= k; i++) {
						ivar += nvar;
						if (*COV(i - 1, i - 1) == F_ZERO || *COV(j - 1, j - 1) ==
						    F_ZERO) {
							*COV(i - 1, j - 1) = xmiss;
							if (wk[i - 1] != F_ZERO && wk[j - 1] != F_ZERO) {

/*								(3, 16, "Correlations are requested but the observations on a variable are constant.  The pertinent correlation coefficients are set to NaN (not a number).");
*/

imsl_ermes(IMSL_WARNING, IMSL_CONSTANT_VARIABLE);
							}
							goto L_230;
						}
						if (mopt <= 0) {
							if (cntm1 < 1) {
								*COV(i - 1, j - 1) = xmiss;
								if (wk[i - 1] != F_ZERO && wk[j - 1] !=
								    F_ZERO) {

/*							(3, 17, "Correlations are requested but fewer than two valid observations are present for a variable.  The pertinent correlation coefficients are set to NaN (not a number).");
*/
imsl_ermes(IMSL_WARNING, IMSL_TOO_FEW_VALID_OBS_CORREL);
								}
								goto L_230;
							}
							*COV(i - 1, j - 1) /= sqrt(*COV(j - 1, j - 1)) *
								sqrt(*COV(i - 1, i - 1));
							goto L_220;
						}
						if (*INCD(i - 1, j - 1) < 2) {
							*COV(i - 1, j - 1) = xmiss;
							if (wk[i - 1] != F_ZERO && wk[j - 1] != F_ZERO) {

/*							 (3, 17, "Correlations are requested but fewer than two valid observations are present for a variable.  The pertinent correlation coefficients are set to NaN (not a number).");
*/
imsl_ermes(IMSL_WARNING, IMSL_TOO_FEW_VALID_OBS_CORREL);
							}
							goto L_230;
						}
						if (mopt != 3)
							*COV(i - 1, j - 1) = (*COV(i - 1, j - 1) *
									      sqrt((*INCD(j - 1, j - 1) - F_ONE) * (*INCD(i - 1, i - 1) -
														  F_ONE))) / ((*INCD(i - 1, j - 1) - F_ONE) * sqrt(*COV(i - 1, i - 1)) *
															    sqrt(*COV(j - 1, j - 1)));
						if (mopt == 3) {
							if (wk[jvar + i - 1] <= F_ZERO || wk[ivar + j - 1] <=
							    F_ZERO) {
								*COV(i - 1, j - 1) = xmiss;
								if (wk[i - 1] != F_ZERO && wk[j - 1] !=
								    F_ZERO) {

/*								(3, 16, "Correlations are requested but the observations on a variable are constant.  The pertinent correlation coefficients are set to NaN (not a number).");
*/
imsl_ermes(IMSL_WARNING, IMSL_TOO_FEW_VALID_OBS_CORREL);
								}
								goto L_230;
							}
							denom = sqrt(wk[jvar + i - 1] * wk[ivar + j - 1]);
							imsl_c1div(*COV(i - 1, j - 1), denom, COV(i - 1, j - 1));
							if (imsl_ifnan(*COV(i - 1, j - 1))) {
								if (wk[i - 1] != F_ZERO && wk[j - 1] !=
								    F_ZERO) {

/*								(3, 16, "Correlations are requested but the observations on a variable are constant.  The pertinent correlation coefficients are set to NaN (not a number).");
*/
imsl_ermes(IMSL_WARNING, IMSL_TOO_FEW_VALID_OBS_CORREL);
								}
								goto L_230;
							}
						}
				L_220:
						if (!imsl_ifnan(*COV(i - 1, j - 1))) {
							if (*COV(i - 1, j - 1) > F_ONE)
								*COV(i - 1, j - 1) = F_ONE;
							if (*COV(i - 1, j - 1) < -F_ONE)
								*COV(i - 1, j - 1) = -F_ONE;
						}
				L_230:
						;
					}
				}
				if (icopt == 3) {
					/*
					 * Compute standard deviations and
					 * store in diagonals
					 */
					for (i = 1; i <= nvar; i++) {
						covl = *COV(i - 1, i - 1);
						*COV(i - 1, i - 1) = xmiss;
						if (mopt <= 0 && cntm1 >= F_ONE)
							*COV(i - 1, i - 1) = sqrt(covl / cntm1);
						if (mopt > 0) {
							if (*INCD(i - 1, i - 1) > 1)
								*COV(i - 1, i - 1) = sqrt(covl / (*INCD(i - 1, i - 1) -
								      F_ONE));
						}
					}
				} else {
					/* Set diagonals to 1.0 */
					for (i = 1; i <= nvar; i++) {
						*COV(i - 1, i - 1) = F_ONE;
					}
				}
			}
		}
	}
	/*
	 * Fill out the upper triangular portion of COV and INCD
	 */
L_270:
	for (i = 1; i <= (nvar - 1); i++) {
		scopy(nvar - i, COV(i - 1, i), 1, COV(i, i - 1), ldcov);
		if (mopt > 0)
			icopy(nvar - i, INCD(i - 1, i), 1, INCD(i, i - 1), ldincd);
	}

L_9000:
	imsl_e1pop("l_c2rvc");
	return;
}				/* end of function */
#ifdef ANSI
static Mint l_i1nan(Mint n, Mfloat *sx, Mint incx)
#else
static Mint l_i1nan(n, sx, incx)
	Mint            n;
	Mfloat          *sx;
	Mint            incx;
#endif
{

	Mint             i, isnan_v, k;
	isnan_v = 0;
	for ( k=1,i=1; k > n; k++){
	  if (imsl_ifnan(*(sx+i-1))) {
	    isnan_v = k;
	    break;
	  }
	  i += incx;
	}
	return(isnan_v);
      }
