#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static void l_u2sta(Mint, Mint, Mint, Mfloat*, Mint, Mint, Mint, Mint, Mfloat,
                    Mfloat, Mint, Mfloat*, Mint, Mint*, Mfloat[]);
static void l_c1div(Mfloat, Mfloat, Mfloat*);
static VA_LIST_HACK l_simple_statistics (Mint, Mint, Mfloat*, va_list);
static Mint l_ismax(Mint, Mfloat[], Mint);
static void l_median(Mint, Mint, Mfloat*, Mint, Mfloat*, Mint*, Mint*, Mint);
static void l_get_mad(Mint, Mint, Mfloat*, Mint, Mfloat*, Mint*, Mint*, Mint);
#else
static void l_u2sta();
static void l_c1div();
static VA_LIST_HACK l_simple_statistics();
static Mint l_ismax();
static void l_median();
static void l_get_mad();
#endif

static Mfloat *lv_simple_statistics = NULL;

#ifdef ANSI
Mfloat *imsl_f_simple_statistics (Mint n_observations, Mint n_variables, Mfloat *x, ...)
#else
Mfloat *imsl_f_simple_statistics (n_observations, n_variables, x, va_alist)
    Mint	n_observations, n_variables;
    Mfloat	*x;
    va_dcl
#endif
{
    va_list	argptr;
    VA_START(argptr, x);
    E1PSH("imsl_f_simple_statistics", "imsl_d_simple_statistics");

    lv_simple_statistics = NULL;
    IMSL_CALL(l_simple_statistics (n_observations, n_variables, x, argptr));
    va_end(argptr);

    E1POP("imsl_f_simple_statistics", "imsl_d_simple_statistics");
    return lv_simple_statistics;
}

#ifdef ANSI
static VA_LIST_HACK l_simple_statistics (Mint n_observations, Mint n_variables, Mfloat *x, va_list argptr)
#else
static VA_LIST_HACK l_simple_statistics (n_observations, n_variables, x, argptr)
    Mint	n_observations, n_variables;
    Mfloat	*x;
    va_list	argptr;
#endif
{
    Mfloat      *road_wk                = NULL;
    Mint        *iperm                  = NULL;
    Mint        *counts                 = NULL;
    Mint	col_dim_x		= n_variables;
    Mint	col_dim_stat		= n_variables;
    Mfloat	confidence_means	= 95.0e0;
    Mfloat	confidence_variances	= 95.0e0;

    Mint	code = 1, ner=0, i, ldstat;
    Mint	n_arguements = 3, user_statistics = 0;
    Mfloat	*wk;
    Mint	n_missing_values = 0;
    Mint        b_median = 0, mutual=0;

    while (code>0){
	code = va_arg(argptr, Mint);
        ++n_arguements;
	switch (code){
            case IMSL_MEDIAN:
                if (b_median) mutual = 1;
                b_median = 1;
                break;
            case IMSL_MEDIAN_AND_SCALE:
                if (b_median) mutual = 1;
                b_median = 2;
                break;
	    case IMSL_RETURN_USER:
		lv_simple_statistics = va_arg(argptr, Mfloat*);
		user_statistics = 1;
		++n_arguements;
		break;
	    case IMSL_X_COL_DIM:
		col_dim_x = va_arg(argptr, Mint);
		++n_arguements;
		break;
	    case IMSL_STAT_COL_DIM:
		col_dim_stat = va_arg(argptr, Mint);
		++n_arguements;
		break;
	    case IMSL_CONFIDENCE_MEANS:
		confidence_means = (Mfloat) va_arg(argptr, Mdouble);
		++n_arguements;
		break;
	    case IMSL_CONFIDENCE_VARIANCES:
		confidence_variances = (Mfloat) va_arg(argptr, Mdouble);
		++n_arguements;
		break;
	    case IMSL_CONFIDENCE_MEANS_ADR:
		confidence_means = *(va_arg(argptr, Mfloat *));
		++n_arguements;
		break;
	    case IMSL_CONFIDENCE_VARIANCES_ADR:
		confidence_variances = *(va_arg(argptr, Mfloat *));
		++n_arguements;
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

    imsl_c1iarg(n_variables, "n_variables", 1, -1, &ner);
    imsl_c1iarg(n_observations, "n_observations", 1, -1, &ner);
    if (imsl_n1rty(0) > 0) goto RETURN;

    if (mutual) {
        imsl_e1stl(1, "IMSL_MEDIAN and IMSL_MEDIAN_AND_SCALE");
	imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXCLUSIVE);
    }

    wk = (Mfloat *) imsl_malloc (2*n_variables*sizeof(Mfloat));
    if (!wk){
        imsl_e1sti (1,n_variables);
        imsl_e1stl (1,"n_variables");
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
        goto RETURN;
    }
    if (b_median>0) {
        road_wk = (Mfloat*)imsl_malloc(n_observations*sizeof(Mfloat));
        iperm = (Mint*) imsl_malloc(n_observations*sizeof(Mint));
        counts = (Mint*) imsl_malloc(n_variables*sizeof(Mint));
    }
    if (b_median==0) {
        ldstat = 14;
    } else if (b_median==1) {
        ldstat = 15;
    } else {
        ldstat = 17;
    }
    if (!user_statistics) {
	    lv_simple_statistics = (Mfloat*)imsl_malloc(ldstat*col_dim_stat*sizeof(Mfloat));
	if (!lv_simple_statistics || (b_median && (!road_wk || !iperm || !counts))) {
            imsl_e1sti (1,n_observations);
            imsl_e1sti (2,n_variables);
            imsl_e1stl (1,"n_variables");
            imsl_e1stl (2,"n_observations");
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	    goto FREE_SPACE;
	}
    }

    imsl_f_m1ran(n_observations, col_dim_x, x, x);

    l_u2sta (0, n_observations, n_variables, x, n_observations, 0, 0, 0, confidence_means, confidence_variances, 0, lv_simple_statistics, ldstat, &n_missing_values, wk);

    if (n_missing_values){
/*        (6,1,"At least one row of X contained NaN (a missing	value).  These rows were excluded from calculations.");
*/
          imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_ROW_OF_X_CONTAINED_NAN);
    }


    if ((code=imsl_n1rty(1))>3 && code!=6 &&!user_statistics){
	imsl_free (lv_simple_statistics);
        lv_simple_statistics = NULL;
        imsl_f_m1ran(col_dim_x, n_observations, x, x);
    }
    else{
	imsl_f_m1ran(col_dim_stat,ldstat,lv_simple_statistics,lv_simple_statistics);
        if (b_median) {
           for (i=0; i<n_variables; i++){
	       *(counts+i) =
	       (Mint)*(lv_simple_statistics+col_dim_stat*9+i);
	   }
           l_median (n_observations, n_variables, x, col_dim_x, road_wk,
	             iperm, counts, col_dim_stat);
        }

        if (b_median==2) {
           l_get_mad (n_observations, n_variables, x, col_dim_x, road_wk,
	             iperm, counts, col_dim_stat);
	}

        imsl_f_m1ran(col_dim_x, n_observations, x, x);
    }

FREE_SPACE:

    if (iperm!=NULL)imsl_free(iperm);
    if (road_wk!=NULL)imsl_free(road_wk);
    if (wk!=NULL)imsl_free(wk);
    if (counts!=NULL)imsl_free(counts);

RETURN:
    if (imsl_n1rty(0) > 3) {
	if(!user_statistics) {
	    if (lv_simple_statistics != NULL) imsl_free(lv_simple_statistics);
	}
        lv_simple_statistics = NULL;
    }
    return argptr;
}

/*
  -----------------------------------------------------------------------
    IMSL Name:  U2STA/DU2STA (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 5, 1990

    Purpose:    Compute basic univariate statistics.

    Usage:      CALL U2STA (IDO, NROW, NVAR, X, LDX, IFRQ, IWT, MOPT,
                            CONPRM, CONPRV, IPRINT, STAT, LDSTAT, NRMISS,
                            WK)

    Arguments:  See UVSTA.

    Chapter:    STAT/LIBRARY Basic Statistics

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u2sta(Mint ido, Mint nrow, Mint nvar, Mfloat *x, Mint ldx, Mint ifrq, Mint iwt, Mint mopt, Mfloat conprm, Mfloat conprv, Mint iprint, Mfloat *stat, Mint ldstat, Mint *nrmiss, Mfloat wk[])
#else
static void l_u2sta(ido, nrow, nvar, x, ldx, ifrq, iwt, mopt, conprm, conprv, iprint, stat, ldstat, nrmiss, wk)
	Mint            ido, nrow, nvar;
	Mfloat         *x;
	Mint            ldx, ifrq, iwt, mopt;
	Mfloat          conprm, conprv;
	Mint            iprint;
	Mfloat         *stat;
	Mint            ldstat, *nrmiss;
	Mfloat          wk[];
#endif
{
#define X(I_,J_)	(x+(I_)*(ldx)+(J_))
#define STAT(I_,J_)	(stat+(I_)*(ldstat)+(J_))
	Mint             i, igo, imaxc, imaxr, imaxw, imiss, irow, j, jj, ner, nrowp;
	static Mint     icall;
	Mfloat          cnt, df, e, etj, frq, fw, fwome, halfl, ome, q, q2, stat2j, stat3j, sumwt, tj, tj2, tj3, tj4, wi, xij, xmiss;

	imsl_e1psh("l_u2sta");
	xmiss = imsl_amach(6);
	ner = 1;
	imsl_c1iarg(ido, "ido", 0, 3, &ner);
	imsl_c1iarg(nvar, "nvar", 1, -1, &ner);
	imsl_c1iarg(ifrq, "ifrq", 0, -1, &ner);
	imsl_c1iarg(iwt, "iwt", 0, -1, &ner);
	imsl_c1iarg(iprint, "iprint", 0, 2, &ner);
	if (imsl_n1rty(0) >= 4)
		goto L_9000;
	if (ido == 3 && nrow == 0)
		goto L_80;
	nrowp = nrow;
	if (nrowp < 0) nrowp = -nrowp;
	if (ido <= 1) {
		if (nrow <= 0) {
			imsl_e1sti(1, nrow);
			imsl_e1sti(2, ido);

/*			imsl_ermes(5, 7, "NROW = %(i1) and IDO = %(i2).  If NROW is nonpositive, IDO must be equal to 2 or 3.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_NONPOSITIVE_NROW_VALUE);
		}
		/* Initialize statistics */
		icall = 1;
		*nrmiss = 0;
		sset(nvar, xmiss, wk, 1);
		sset(nvar, F_ZERO, &wk[nvar], 1);
		sset(nvar, F_ZERO, STAT(0, 0), ldstat);
		sset(nvar, F_ZERO, STAT(0, 1), ldstat);
		sset(nvar, F_ZERO, STAT(0, 2), ldstat);
		sset(nvar, F_ZERO, STAT(0, 3), ldstat);
		sset(nvar, F_ZERO, STAT(0, 4), ldstat);
		sset(nvar, xmiss, STAT(0, 5), ldstat);
		sset(nvar, xmiss, STAT(0, 6), ldstat);
		sset(nvar, F_ZERO, STAT(0, 8), ldstat);
		sset(nvar, F_ZERO, STAT(0, 9), ldstat);
		if (iwt > 0)
			sset(nvar, F_ZERO, STAT(0, 14), ldstat);
	} else {
		icall += 1;
	}
	if (ido <= 1 || ido == 3) {
		if (conprm >= 100.0) {
			imsl_e1str(1, conprm);

/*			imsl_ermes(5, 8, "CONPRM = %(r1).  The confidence level for the mean, CONPRM, must be less than 100.0.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_CONPRM_VALUE_TOO_BIG);
		}
		if (conprv >= 100.0) {
			imsl_e1str(1, conprv);

/*			imsl_ermes(5, 9, "CONPRV = %(r1).  The confidence level for the variance, CONPRV, must be less than 100.0.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_CONPRV_VALUE_TOO_BIG);
		}
		if (imsl_n1rty(0) > 0)
			goto L_9000;
	}
	if (nrow < 0) {
		irow = -1;
	} else {
		irow = 1;
	}
	/*
	 * Accumulate means and adjusted sums of squares and test for
	 * constant observation set.
	 */
	if (mopt <= 0) {
		/* Listwise deletion option */
		cnt = *STAT(0, 9);
		if (iwt > 0) {
			sumwt = *STAT(0, 14);
		} else {
			sumwt = cnt;
		}
		for (i = 1; i <= nrowp; i++) {
			imsl_c1wfr(ido, icall, x, ldx, i, irow, ifrq, iwt, xmiss,
			      nrmiss, &frq, &wi, &igo);
			if (igo == 3)
				goto L_9000;
			if (igo == 1 || igo == 2)
				goto L_30;
			jj = 1;
			if (iwt <= 0)
				wi = F_ONE;
			for (j = 1; j <= nvar; j++) {
				if (jj == ifrq || jj == iwt)
					jj += 1;
				if (imsl_ifnan(*X(jj - 1, i - 1))) {
					if (nrow > 0) {
						*nrmiss += 1;
					} else {
						*nrmiss -= 1;
					}
					goto L_30;
				}
				jj += 1;
			}
			fw = frq * wi;
			sumwt += fw;
			cnt += frq;
			if (cnt < F_ZERO) {

/*				imsl_ermes(5, 10, "The number of nonmissing observations is less than zero.");
*/
                                imsl_ermes(IMSL_TERMINAL,
				IMSL_NUM_NONMISS_OBS_LT_ZERO);
				goto L_9000;
			} else if (cnt == F_ZERO) {
				sset(nvar, F_ZERO, STAT(0, 0), ldstat);
				sset(nvar, F_ZERO, STAT(0, 1), ldstat);
				sset(nvar, F_ZERO, STAT(0, 3), ldstat);
				sset(nvar, F_ZERO, STAT(0, 4), ldstat);
				if (iwt > 0)
					sset(nvar, F_ZERO, STAT(0, 14), ldstat);
				goto L_30;
			}
			if (wi > F_ZERO) {
				e = fw / sumwt;
				ome = F_ONE - e;
				fwome = fw * ome;
			}
			jj = 1;
			for (j = 1; j <= nvar; j++) {
				if (jj == ifrq || jj == iwt)
					jj += 1;
				xij = *X(jj - 1, i - 1);
				jj += 1;
				if (nrow > 0) {
					/* UPDATE MINIMA, MAXIMA */
					if (imsl_ifnan(wk[j - 1])) {
						*STAT(j - 1, 5) = xij;
						*STAT(j - 1, 6) = xij;
						wk[j - 1] = xij;
					} else {
						if (xij < *STAT(j - 1, 5))
							*STAT(j - 1, 5) = xij;
						if (xij > *STAT(j - 1, 6))
							*STAT(j - 1, 6) = xij;
						if (xij != wk[j - 1])
							wk[nvar + j - 1] = F_ONE;
					}
				}
				if (wi > F_ZERO) {
					/* Update/downdate moments */
					tj = xij - *STAT(j - 1, 0);
					tj2 = tj * tj;
					tj3 = tj * tj2;
					tj4 = tj * tj3;
					etj = e * tj;
					*STAT(j - 1, 4) += -F_FOUR * etj ** STAT(j - 1, 3) + F_SIX *
						etj * etj ** STAT(j - 1, 1) + fwome * tj4 * (F_ONE - F_THREE *
								   e * ome);
					*STAT(j - 1, 3) += -F_THREE * etj ** STAT(j - 1, 1) + fwome *
						tj3 * (F_ONE - F_TWO * e);
					*STAT(j - 1, 1) += fwome * tj2;
					*STAT(j - 1, 0) += etj;
				}
			}
	L_30:
			;
		}
		for (j = 1; j <= nvar; j++) {
			*STAT(j - 1, 9) = cnt;
			if (iwt > 0)
				*STAT(j - 1, 14) = sumwt;
		}
	} else {
		/* Elementwise deletion option */
		for (i = 1; i <= nrowp; i++) {
			imiss = 0;
			imsl_c1wfr(ido, icall, x, ldx, i, irow, ifrq, iwt, xmiss,
			      nrmiss, &frq, &wi, &igo);
			if (igo == 3)
				goto L_9000;
			if (igo == 1 || igo == 2)
				goto L_60;
			if (iwt <= 0)
				wi = F_ONE;
			jj = 1;
			for (j = 1; j <= nvar; j++) {
				if (jj == ifrq || jj == iwt)
					jj += 1;
				xij = *X(jj - 1, i - 1);
				jj += 1;
				if (imsl_ifnan(xij)) {
					imiss += 1;
					goto L_50;
				}
				if (nrow > 0) {
					/* UPDATE MINIMA, MAXIMA */
					if (imsl_ifnan(wk[j - 1])) {
						*STAT(j - 1, 5) = xij;
						*STAT(j - 1, 6) = xij;
						wk[j - 1] = xij;
					} else {
						if (xij < *STAT(j - 1, 5))
							*STAT(j - 1, 5) = xij;
						if (xij > *STAT(j - 1, 6))
							*STAT(j - 1, 6) = xij;
						if (xij != wk[j - 1])
							wk[nvar + j - 1] = F_ONE;
					}
				}
				*STAT(j - 1, 9) += frq;
				if (*STAT(j - 1, 9) < F_ZERO) {

/*					imsl_ermes(5, 10, "The number of nonmissing observations is less than zero.");
*/
                                        imsl_ermes(IMSL_TERMINAL,
					IMSL_NUM_NONMISS_OBS_LT_ZERO);
					goto L_9000;
				} else if (*STAT(j - 1, 9) == F_ZERO) {
					*STAT(j - 1, 0) = F_ZERO;
					*STAT(j - 1, 1) = F_ZERO;
					*STAT(j - 1, 3) = F_ZERO;
					*STAT(j - 1, 4) = F_ZERO;
					if (iwt > F_ZERO)
						*STAT(j - 1, 14) = F_ZERO;
					goto L_60;
				}
				if (wi > F_ZERO) {
					/* Update/downdate moments */
					if (iwt > 0) {
						sumwt = *STAT(j - 1, 14);
					} else {
						sumwt = *STAT(j - 1, 9) - frq;
					}
					fw = frq * wi;
					sumwt += fw;
					e = fw / sumwt;
					ome = F_ONE - e;
					fwome = fw * ome;
					tj = xij - *STAT(j - 1, 0);
					tj2 = tj * tj;
					tj3 = tj * tj2;
					tj4 = tj * tj3;
					etj = e * tj;
					*STAT(j - 1, 4) += -F_FOUR * etj ** STAT(j - 1, 3) + F_SIX *
						etj * etj ** STAT(j - 1, 1) + fwome * tj4 * (F_ONE - F_THREE *
								   e * ome);
					*STAT(j - 1, 3) += -F_THREE * etj ** STAT(j - 1, 1) + fwome *
						tj3 * (F_ONE - F_TWO * e);
					*STAT(j - 1, 1) += fwome * tj2;
					*STAT(j - 1, 0) += etj;
					if (iwt > 0)
						*STAT(j - 1, 14) = sumwt;
				}
		L_50:
				;
			}
			if (imiss > 0) {
				if (nrow > 0) {
					*nrmiss += 1;
				} else {
					*nrmiss -= 1;
				}
			}
	L_60:
			;
		}
	}
	if (ido == 1 || ido == 2) goto L_9000;
	/*
	 * Postprocessing begins here
	 */
L_80:
	;
	imaxc = l_ismax(nvar, STAT(0, 9), ldstat);
	if (iwt > 0)
		imaxw = l_ismax(nvar, STAT(0, 14), ldstat);
	if (*STAT(imaxc - 1, 9) <= F_ONE) {
		/*
		 * Fewer than two valid observations for all variables
		 */
		for (j = 1; j <= nvar; j++) {
			sset(4, xmiss, STAT(j - 1, 1), 1);
			*STAT(j - 1, 8) = xmiss;
			if (conprm > F_ZERO) {
				*STAT(j - 1, 10) = xmiss;
				*STAT(j - 1, 11) = xmiss;
			}
			if (conprv > F_ZERO) {
				*STAT(j - 1, 12) = xmiss;
				*STAT(j - 1, 13) = xmiss;
			}
		}

/*		imsl_ermes(6, 1, "Fewer than two valid observations are present.  The corresponding statistics are set to NaN (not a number), (except for the mean, which is not correct if no valid observations ");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE,
		IMSL_LESS_THAN_TWO_VALID_OBS);
		/*
		 * &               'are present, or is correct if one
		 * observation '// &               'is present).')
		 */
		goto L_9000;
	}
	if (iwt > 0) {
		if (*STAT(imaxw - 1, 14) < F_TEN * imsl_amach(4)) {
			/*
			 * Sum of the weights is zero for all variables
			 */
			for (j = 1; j <= nvar; j++) {
				sset(5, xmiss, STAT(j - 1, 0), 1);
				*STAT(j - 1, 8) = xmiss;
				if (conprm > F_ZERO) {
					*STAT(j - 1, 10) = xmiss;
					*STAT(j - 1, 11) = xmiss;
				}
				if (conprv > F_ZERO) {
					*STAT(j - 1, 12) = xmiss;
					*STAT(j - 1, 13) = xmiss;
				}
			}

/*			imsl_ermes(6, 15, "The sum of the weights is zero. The statistics, except for the minima, maxima, ranges and counts, are set to NaN (not a number).");
*/
                        imsl_ermes(IMSL_WARNING_IMMEDIATE,
			IMSL_ZERO_SUM_OF_WEIGHTS);
			goto L_9000;
		}
	}
	for (j = 1; j <= nvar; j++) {
		if (*STAT(j - 1, 9) > F_ONE) {
			if (iwt > 0 && *STAT(j - 1, 14) < F_TEN * imsl_amach(4)) {
				sset(5, xmiss, STAT(j - 1, 0), 1);
				*STAT(j - 1, 8) = xmiss;
				if (conprm > F_ZERO) {
					*STAT(j - 1, 10) = xmiss;
					*STAT(j - 1, 11) = xmiss;
				}
				if (conprv > F_ZERO) {
					*STAT(j - 1, 12) = xmiss;
					*STAT(j - 1, 13) = xmiss;
				}
				imsl_e1sti(1, j);

/*				imsl_ermes(6, 15, "The sum of the weights for variable %(i1) is zero. The statistics, except for the minima, maxima, ranges and counts, are set to NaN (not a number).");
*/
                                imsl_ermes(IMSL_WARNING_IMMEDIATE,
				IMSL_SUM_OF_WEIGHTS_ZERO);
			}
		}
	}
	/* Compute the range */
	for (j = 1; j <= nvar; j++) {
		if (imsl_ifnan(*STAT(j - 1, 6)) || imsl_ifnan(*STAT(j - 1, 5))) {
			*STAT(j - 1, 7) = F_ZERO;
		} else {
			*STAT(j - 1, 7) = *STAT(j - 1, 6) - *STAT(j - 1, 5);
		}
	}
	imaxr = l_ismax(nvar, STAT(0, 7), ldstat);
	if (*STAT(imaxr - 1, 7) < F_ZERO) {
		for (j = 1; j <= nvar; j++) {
			sset(4, xmiss, STAT(j - 1, 1), 1);
			*STAT(j - 1, 7) = xmiss;
			*STAT(j - 1, 8) = xmiss;
			if (conprm > F_ZERO) {
				*STAT(j - 1, 10) = xmiss;
				*STAT(j - 1, 11) = xmiss;
			}
			if (conprv > F_ZERO) {
				*STAT(j - 1, 12) = xmiss;
				*STAT(j - 1, 13) = xmiss;
			}
		}

/*		imsl_ermes(6, 4, "The maximum value is less than the minimum value.  The corresponding statistics are  set to NaN (not a number).");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_MAX_LESS_THAN_MIN);
		goto L_9000;
	} else {
		for (j = 1; j <= nvar; j++) {
			if (*STAT(j - 1, 7) < F_ZERO) {
				sset(4, xmiss, STAT(j - 1, 1), 1);
				*STAT(j - 1, 7) = xmiss;
				*STAT(j - 1, 8) = xmiss;
				if (conprm > F_ZERO) {
					*STAT(j - 1, 10) = xmiss;
					*STAT(j - 1, 11) = xmiss;
				}
				if (conprv > F_ZERO) {
					*STAT(j - 1, 12) = xmiss;
					*STAT(j - 1, 13) = xmiss;
				}
				imsl_e1sti(1, j);

/*				imsl_ermes(6, 4, "The maximum value is less than the minimum value for variable %(i1).  The corresponding statistics are  set to NaN (not a number).");
*/
                                imsl_ermes(IMSL_WARNING_IMMEDIATE,
				IMSL_MIN_GREATER_THAN_MAX);
			}
		}
	}
	/*
	 * Compute the variance, standard deviation, skewness and kurtosis
	 */
	for (j = 1; j <= nvar; j++) {
		if (!imsl_ifnan(*STAT(j - 1, 7))) {
			if (*STAT(j - 1, 9) <= F_ONE) {
				sset(4, xmiss, STAT(j - 1, 1), 1);
				*STAT(j - 1, 8) = xmiss;
				if (conprm > F_ZERO) {
					*STAT(j - 1, 10) = xmiss;
					*STAT(j - 1, 11) = xmiss;
				}
				if (conprv > F_ZERO) {
					*STAT(j - 1, 12) = xmiss;
					*STAT(j - 1, 13) = xmiss;
				}
				imsl_e1sti(1, j);

/*				imsl_ermes(6, 1, "Fewer than two valid observations are
present for variable %(i1).  The corresponding statistics are set to NaN (not a number)
, (except for the mean, which is not correct if no valid observations are present, or
is correct if one observation is present)."); */
                                imsl_ermes(IMSL_WARNING_IMMEDIATE,
				IMSL_NOT_ENOUGH_OBSERVATIONS);
				goto L_150;
			} else if (wk[nvar + j - 1] != F_ZERO && *STAT(j - 1, 1) > F_ZERO) {
				stat2j = *STAT(j - 1, 1) / *STAT(j - 1, 9);
				stat3j = sqrt(stat2j);
				if (*STAT(j - 1, 7) < sqrt(imsl_amach(1))) {
					*STAT(j - 1, 1) = F_ZERO;
					*STAT(j - 1, 2) = F_ZERO;
					*STAT(j - 1, 3) = xmiss;
					*STAT(j - 1, 4) = xmiss;
					imsl_e1sti(1, j);

/*					imsl_ermes(6, 12, "Since the range of variable
%(i1) is very small, the variance for this variable underflows.  Therefore, the variance
and standard deviation are set to 0, and the skewness and kurtosis are set to NaN (not
a number).");
*/
                                        imsl_ermes(IMSL_WARNING_IMMEDIATE,
					IMSL_VARIANCE_UNDERFLOW);
				} else {
					*STAT(j - 1, 1) /= *STAT(j - 1, 9) - F_ONE;
					*STAT(j - 1, 2) = sqrt(*STAT(j - 1, 1));
					if (*STAT(j - 1, 7) < pow(imsl_amach(1), F_ONE / F_THREE)) {
						*STAT(j - 1, 3) = xmiss;
						*STAT(j - 1, 4) = xmiss;
						imsl_e1sti(1, j);

/*						imsl_ermes(6, 13, "Since the range of variable %(i1) is very small, the higher order moments for this variable underflow.  Therefore, the skewness and kurtosis are set to NaN (not a number).");
*/

imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_HIGH_ORDER_UNDERFLOW);
					} else if (F_THREE * log(stat3j) < log(imsl_amach(1))) {
						*STAT(j - 1, 3) /= *STAT(j - 1, 9) * stat3j;
						*STAT(j - 1, 3) /= stat3j;
						*STAT(j - 1, 3) /= stat3j;
					} else {
						*STAT(j - 1, 3) /= *STAT(j - 1, 9) * pow(stat3j, F_THREE);
					}
					if (*STAT(j - 1, 7) < pow(imsl_amach(1), .25)) {
						*STAT(j - 1, 4) = xmiss;
						imsl_e1sti(1, j);

/*						imsl_ermes(6, 14, "Since the range of variable %(i1) is very small, the fourth order moment for this variable underflows.  Therefore, the kurtosis is set to NaN (not a number).");
*/

imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_FOURTH_ORDER_UNDERFLOW);
					} else {
						if (F_TWO * log(stat2j) < log(imsl_amach(1))) {
							*STAT(j - 1, 4) /= *STAT(j - 1, 9) * stat2j;
							*STAT(j - 1, 4) = *STAT(j - 1, 4) / stat2j -
								F_THREE;
						} else {
							*STAT(j - 1, 4) = *STAT(j - 1, 4) / (*STAT(j - 1, 9) *
											     pow(stat2j, F_TWO)) - F_THREE;
						}
					}
				}
			} else if (wk[nvar + j - 1] == 0) {
				imsl_e1sti(1, j);

/*				imsl_ermes(6, 2, "The observations on variable %(i1) are constant.");
*/
                                imsl_ermes(IMSL_WARNING_IMMEDIATE,
				IMSL_CONSTANT_OBSERVATIONS);
				*STAT(j - 1, 1) = F_ZERO;
				*STAT(j - 1, 2) = F_ZERO;
				*STAT(j - 1, 3) = F_ZERO;
				*STAT(j - 1, 4) = F_ZERO;
			} else if (!imsl_ifnan(*STAT(j - 1, 1))) {
				*STAT(j - 1, 2) = *STAT(j - 1, 1);
				imsl_e1sti(1, j);

/*				imsl_ermes(6, 3, "The variance and standard deviation are negative for variable %(i1).");
*/
                                imsl_ermes(IMSL_WARNING_IMMEDIATE,
				IMSL_VAR_AND_STD_ARE_NEGATIVE);
			}
			/* Compute coefficients of variation */
			if (!imsl_ifnan(*STAT(j - 1, 8))) {
				if (imsl_ifnan(*STAT(j - 1, 0)) || fabs(*STAT(j - 1, 0)) <
				    imsl_amach(4)) {
					*STAT(j - 1, 8) = xmiss;
				} else {
					l_c1div(*STAT(j - 1, 2), *STAT(j - 1, 0), STAT(j - 1, 8));
				}
				if (imsl_ifnan(*STAT(j - 1, 8))) {
					imsl_e1sti(1, j);

/*					imsl_ermes(6, 5, "The coefficient of variation is set to NaN (not a number) for variable %(i1).  This is due to the value of the mean, standard deviation, or both.");
*/
                                        imsl_ermes(IMSL_WARNING_IMMEDIATE,
					IMSL_COEFF_OF_VARIATION_NAN);
				}
			}
			/* Compute confidence intervals for mean */
			if (conprm > F_ZERO) {
				q = F_HALF + conprm / 200.0;
				df = *STAT(j - 1, 9) - F_ONE;
				if (*STAT(j - 1, 2) < F_ZERO) {
					imsl_e1sti(1, j);

/*					imsl_ermes(6, 7, "The standard deviation is negative for variable %(i1).  The corresponding confidence limits for the mean are set to NaN (not a number).");
*/
                                        imsl_ermes(IMSL_WARNING_IMMEDIATE,
					IMSL_NEGATIVE_STD_VALUE);
					*STAT(j - 1, 10) = xmiss;
					*STAT(j - 1, 11) = xmiss;
				} else {
					halfl = imsl_f_t_inverse_cdf(q, df);
					if (imsl_n1rty(1) >= 5) {

/*						imsl_ermes(6, 8, "An error occured in determining the t statistic.  The confidence limits are set to NaN (not a number).");
*/

imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_ERROR_IN_T_STATISTIC);
						*STAT(j - 1, 10) = xmiss;
						*STAT(j - 1, 11) = xmiss;
					} else {
						if (imsl_ifnan(*STAT(j - 1, 2))) {
							*STAT(j - 1, 10) = xmiss;
							*STAT(j - 1, 11) = xmiss;
						} else {
							halfl = halfl ** STAT(j - 1, 2) / sqrt(*STAT(j - 1, 9));
							*STAT(j - 1, 10) = *STAT(j - 1, 0) - halfl;
							*STAT(j - 1, 11) = *STAT(j - 1, 0) + halfl;
						}
					}
				}
			}
			/*
			 * Compute confidence intervals for variances
			 */
			if (conprv > F_ZERO) {
				q = F_HALF * (F_ONE - conprv / 100.0);
				q2 = F_ONE - q;
				df = *STAT(j - 1, 9) - F_ONE;
				if (*STAT(j - 1, 1) < F_ZERO) {
					imsl_e1sti(1, j);

/*					imsl_ermes(6, 9, "The variance is negative for variable %(i1).  The corresponding confidence limits for the variance are set to NaN (not a number).");
*/
                                        imsl_ermes(IMSL_WARNING_IMMEDIATE,
					IMSL_NEGATIVE_VARIANCE);
					*STAT(j - 1, 12) = xmiss;
					*STAT(j - 1, 13) = xmiss;
				} else if (imsl_ifnan(*STAT(j - 1, 1))) {
					*STAT(j - 1, 12) = xmiss;
					*STAT(j - 1, 13) = xmiss;
				} else {
					halfl = imsl_f_chi_squared_inverse_cdf(q2,df);
					if (imsl_n1rty(1) >= 4) {

/*						imsl_ermes(6, 10, "An error occured in determining the chi-squared statistic.  The lower confidence limit for the variance is set to NaN (not a number).");
*/

imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_CHI_SQUARED_STAT_ERROR);
						*STAT(j - 1, 12) = xmiss;
					} else {
						*STAT(j - 1, 12) = *STAT(j - 1, 1) * df / halfl;
					}
					halfl = imsl_f_chi_squared_inverse_cdf(q,df);
					if (imsl_n1rty(1) >= 4) {

/*						imsl_ermes(6, 11, "An error occured in determining the chi-squared statistic.  The upper confidence limit for the variance is set to NaN (not a number).");
*/

imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_CHI_SQUARED_STAT_ERROR);
						*STAT(j - 1, 13) = xmiss;
					} else {
						*STAT(j - 1, 13) = *STAT(j - 1, 1) * df / halfl;
					}
				}
			}
		}

L_150:
	;
	}
L_9000:
	;
	imsl_e1pop("l_u2sta");
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  C1DIV/DC1DIV (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 21, 1985

    Purpose:    Perform a divide check, and output the quotient.

    Usage:      CALL C1DIV (TOP, BOTTOM, QUOT)

    Arguments:
       TOP    - Numerator.  (Input)
       BOTTOM - Denominator.  (Input)
       QUOT   - Quotient.  (Output)

    Chapter:    STAT/LIBRARY Utilities (not documented)

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c1div(Mfloat top, Mfloat bottom, Mfloat *quot)
#else
static void l_c1div(top, bottom, quot)
        Mfloat      top, bottom, *quot;
#endif
{
        Mfloat      absbot, anan, big, small;


        anan = imsl_amach(6);
        if (imsl_ifnan(top) || imsl_ifnan(bottom)) {
                *quot = anan;
                return;
        }
        absbot = fabs(bottom);
        if (absbot <= F_ONE) {
                /*
                 * ABS(QUOT) is greater or equal to AMACH(1) or equal to 0.0.
                 */
                big = imsl_amach(2);
                if (fabs(top) < big * absbot)
                        goto L_10;
                /*
                 * TOP and BOTTOM are 0.0 or ABS(QUOT) is greater than or
                 * equal to BIG.
                 */
                if (top == F_ZERO) {
                        *quot = anan;
                } else if (bottom >= F_ZERO) {
                        if (top >= F_ZERO) {
                                *quot = imsl_amach(7);
                        } else {
                                *quot = imsl_amach(8);
                        }
                } else {
                        if (top >= F_ZERO) {
                                *quot = imsl_amach(8);
                        } else {
                                *quot = imsl_amach(7);
                        }
                }
                return;
        } else {
                /*
                 * ABS(QUOT) is less than AMACH(2) and greater than or equal
                 * to 0.0.
                 */
                small = imsl_amach(1);
                if (fabs(top) >= small * absbot)
                        goto L_10;
                /*
                 * ABS(QUOT) is greater than or equal to 0.0 and less than
                 * AMACH(1).
                 */
                *quot = F_ZERO;
                return;
        }
L_10:
        *quot = top / bottom;
        return;
}                               /* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  ISMAX (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Find the smallest index of the component of a
                single-precision vector having maximum value.

    Usage:      ISMAX(N, SX, INCX)

    Arguments:
       N      - Length of vector X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than zero.
       ISMAX  - The smallest index I such that X(I)  is the maximum of
                X(J) for J=1 to N.  (Output)
                X(I) refers to a specific element of SX. See INCX
                argument description.

    Keyword:    Level 1 BLAS

    GAMS:       D1a2

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mint l_ismax(Mint n, Mfloat sx[], Mint incx)
#else
static Mint l_ismax(n, sx, incx)
	Mint           n;
	Mfloat           sx[];
	Mint           incx;
#endif
{
	Mint            i, ismax_v, ix;
	Mfloat           smax;


	ismax_v = 0;
	if (n >= 1) {
		ismax_v = 1;
		if (n != 1) {
			if (incx != 1) {
				/* CODE FOR INCREMENT NOT EQUAL TO 1 */
				ix = 1;
				smax = sx[0];
				ix += incx;
				for (i = 2; i <= n; i++) {
					if (sx[ix - 1] > smax) {
						ismax_v = i;
						smax = sx[ix - 1];
					}
					ix += incx;
				}
			} else {
				/* CODE FOR INCREMENT EQUAL TO 1 */
				smax = sx[0];
				for (i = 2; i <= n; i++) {
					if (sx[i - 1] > smax) {
						ismax_v = i;
						smax = sx[i - 1];
					}
				}
			}
		}
	}
	return (ismax_v);
}				/* end of function */
/* X_COL_DIM is not used , but keep the calling sequence intact. */
#ifdef ANSI
static void l_median(Mint nobs, Mint nvar, Mfloat *x, Mint x_col_dim, Mfloat *road_wk,  Mint *iperm, Mint *counts, Mint col_dim_stat)
#else
static void l_median(nobs, nvar, x, x_col_dim, road_wk, iperm, counts, col_dim_stat)
  Mint    nobs, nvar, x_col_dim, *iperm, *counts, col_dim_stat;
  Mfloat  *x, *road_wk;
#endif
{
  Mint    var, quot, rem;

  for (var=0; var<nvar; var++){
    if (counts[var] == nobs) {
      scopy(nobs, &x[var*nobs], 1, &road_wk[0], 1);
      imsl_svrgp(nobs, road_wk, road_wk, iperm);
      quot =  nobs/2;
      rem  =  nobs%2;
      if (rem != 0) {
         lv_simple_statistics[var+14*col_dim_stat] = road_wk[quot];
      } else {
         lv_simple_statistics[var+14*col_dim_stat] = (road_wk[quot-1] + road_wk[quot])/F_TWO;
      }
    } else {
	lv_simple_statistics[var+14*col_dim_stat] = imsl_amach(6);
    }
  }
  return;
}




#ifdef ANSI
static void l_get_mad(Mint nobs, Mint nvar, Mfloat *x, Mint x_col_dim, Mfloat *road_wk,  Mint *iperm, Mint *counts, Mint col_dim_stat)
#else
static void l_get_mad(nobs, nvar, x, x_col_dim, road_wk, iperm, counts, col_dim_stat)
  Mint    nobs, nvar, x_col_dim, *iperm, *counts, col_dim_stat;
  Mfloat  *x, *road_wk;
#endif
{
  Mint    i, var, quot, rem, index, start;
  Mfloat  med, standardize;
  standardize = imsl_f_normal_inverse_cdf(0.75);
  for (var=0; var<nvar; var++){
     med = *(lv_simple_statistics+col_dim_stat*14+var);
     start = var*x_col_dim;
     for (i=0; i<nobs; i++) {
       index = start + i;
       road_wk[i] = fabs(*(x+index) - med);
     }

     if (counts[var]==nobs) {
       imsl_svrgp(nobs, road_wk, road_wk, iperm);
       quot =  nobs/2;
       rem  =  nobs%2;

       if (rem != 0) {
         lv_simple_statistics[var+15*col_dim_stat] = road_wk[quot];
       } else {
         lv_simple_statistics[var+15*col_dim_stat] = (road_wk[quot-1] + road_wk[quot])/F_TWO;
       }
       lv_simple_statistics[var+16*col_dim_stat] =
       lv_simple_statistics[var+15*col_dim_stat]/standardize;

     } else {
       lv_simple_statistics[var+15*col_dim_stat] = imsl_amach(6);
       lv_simple_statistics[var+16*col_dim_stat] = imsl_amach(6);
     }
  }
}

/* IFNAN/DIFNAN */

#ifdef COMPUTER_APLC

#ifdef imsl_ifnan
#undef imsl_ifnan
#endif

#ifndef DOUBLE

#ifdef ANSI
Mint imsl_ifnan(Mfloat x)
#else
Mint imsl_ifnan(x)
Mfloat x;
#endif /*ANSI*/
{
	Mfloat tmpaplc;
        Mfloat tmpaplc2;
	tmpaplc = x;
        tmpaplc2 = x;
        if (tmpaplc != tmpaplc2) return 1;
        return 0;
}

#else /*if DOUBLE*/

#ifdef ANSI
Mint imsl_difnan(Mdouble x)
#else
Mint imsl_difnan(x)
Mdouble x;
#endif /*ANSI*/
{
        Mdouble tmpaplc;
        Mdouble tmpaplc2;
        tmpaplc = x;
        tmpaplc2 = x;
        if (tmpaplc != tmpaplc2) return 1;
        return 0;
}

#endif /* float/double */

#endif /* COMPUTER_APLC */

#if defined(IMSL_MACHINE_80X86) || defined(IMSL_MACHINE_NT)

/*  NaNs have the maximum biased exponent -- 255 in float,
    2047 in double -- and a non-zero fraction.  The sign
    bit is ignored.  If the high-order bit of the fractional
    part is 1, it is a quiet NaN (QNaN); otherwise it is a
    signaling NaN.  */

#define F_MASK_EXP	0x7f800000	/* float exponent */
#define F_MASK_FRAC	0x007fffff	/* float fraction */
#ifndef COMPUTER_LINUX
typedef
#endif
union _if {
    Mlong   i;
    Mfloat  f;
};

#define D_MASK_EXP	0x7ff00000	/* double exponent */
#define D_MASK_FRAC	0x000fffff	/* double fraction */
#ifndef COMPUTER_LINUX
typedef
#endif
union _id {
    Mlong   i[2];
    Mdouble d;
};

#ifndef DOUBLE

#ifdef ANSI
Mint imsl_ifnan(Mfloat x)
#else
Mint imsl_ifnan(x)
Mfloat x;
#endif /*ANSI*/
{
    union _if	_lx;

    _lx.f = x;

    if ((_lx.i & F_MASK_EXP) == F_MASK_EXP &&
	(_lx.i & F_MASK_FRAC) != 0) return 1;

    return 0;
}

#else /*if DOUBLE*/

#ifdef ANSI
Mint imsl_difnan(Mdouble x)
#else
Mint imsl_difnan(x)
Mdouble x;
#endif /*ANSI*/
{
    union _id	_lx;

    _lx.d = x;

    if ((_lx.i[1] & D_MASK_EXP) == D_MASK_EXP &&
	((_lx.i[1] & D_MASK_FRAC) != 0 || _lx.i[0] != 0)) return 1;

    return 0;
}

#endif /* float/double */

#endif /* IMSL_MACHINE_80X86 */
