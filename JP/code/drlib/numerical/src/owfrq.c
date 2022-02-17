#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static void PROTO(l_owfrq,(Mint *nobs, Mfloat x[], Mint *k, Mint *iopt, 
        Mfloat *xlo, Mfloat *xhi, Mfloat *clhw, Mfloat div[],Mfloat
	table[]));
static Mint PROTO(l_isanan,(Mint n, Mfloat sx[], Mint incx));
static Mint PROTO(l_ismax,(Mint, Mfloat[], Mint));
static Mint PROTO(l_ismin,(Mint n, Mfloat sx[], Mint incx));
static VA_LIST_HACK PROTO(l_table_oneway, (Mint nobs, Mfloat y[], Mint n_intervals, va_list argptr));

/****************************** USER INTERFACE ************************/

static Mfloat *lv_table;

#ifdef ANSI
Mfloat *imsl_f_table_oneway(Mint nobs, Mfloat y[], Mint n_intervals, ...)
#else
Mfloat *imsl_f_table_oneway(nobs, y, n_intervals, va_alist)
    Mint        nobs, n_intervals;
    Mfloat      y[];
    va_dcl
#endif
{
    va_list     argptr;
    VA_START(argptr,n_intervals);
    E1PSH("imsl_f_table_oneway","imsl_d_table_oneway");
    lv_table = NULL;
    IMSL_CALL(l_table_oneway(nobs, y, n_intervals, argptr));
    va_end(argptr);
    E1POP("imsl_f_table_oneway","imsl_d_table_oneway");
    return lv_table;
}

/******************************* SUPER CODE ****************************/
#ifdef ANSI
static VA_LIST_HACK l_table_oneway(Mint nobs, Mfloat y[], Mint n_intervals, va_list argptr)
#else
static VA_LIST_HACK l_table_oneway(nobs, y, n_intervals, argptr)
    Mint        nobs, n_intervals;
    Mfloat      y[];
    va_list     argptr;
#endif
{
    Mint    code       = 1;
    Mint    arg_number = 3;
    Mint    user_table = 0;
    Mint    return_bounds = 0;
    Mfloat  *lv_min = NULL;
    Mfloat  *lv_max = NULL;
    Mint    choice_made=0;

    Mint    iopt       = 0; /*default*/
    Mfloat  xlo,xhi, clhw;
    Mfloat  *div = NULL;  /*length depends on iopt and n_intervals*/
    Mint	user_div = 0;

    while (code > 0) {
        code = va_arg(argptr, int);
        ++arg_number;
        switch (code) {
            case IMSL_RETURN_USER:
                lv_table = va_arg(argptr, Mfloat*);
                user_table  = 1;
                arg_number++;
                break;
            case IMSL_DATA_BOUNDS:
                lv_min = va_arg(argptr, Mfloat*);
                lv_max = va_arg(argptr, Mfloat*);
                iopt = 0; return_bounds = 1;
                arg_number += 2;
                choice_made++;
                break;
            case IMSL_KNOWN_BOUNDS:
                xlo  = (Mfloat) va_arg(argptr, Mdouble);
                xhi  = (Mfloat) va_arg(argptr, Mdouble);
                iopt = 1;
                arg_number += 2;
                choice_made++;
                break;
            case IMSL_KNOWN_BOUNDS_ADR:
                xlo  = *(va_arg(argptr, Mfloat *));
                xhi  = *(va_arg(argptr, Mfloat *));
                iopt = 1;
                arg_number += 2;
                choice_made++;
                break;
            case IMSL_CUTPOINTS:
                div = va_arg(argptr, Mfloat*);
		user_div = 1;
                iopt = 2;
                arg_number++;
                choice_made++;
                break;
            case IMSL_CLASS_MARKS:
                div = va_arg(argptr, Mfloat*);
		user_div = 1;
                iopt = 3;
                arg_number++;
                choice_made++;
                break;
            case 0:
                break;
            default:
                imsl_e1sti (1, code);
                imsl_e1sti (2, arg_number);
                imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
                return argptr;
        }
    }
    if (choice_made > 1) {
        imsl_ermes (IMSL_TERMINAL,IMSL_MUT_EXC_TALLY_OPT);
        goto RETURN;
    }

    if (iopt == 1 && n_intervals < 3) {
        imsl_e1sti(1, n_intervals);
        imsl_ermes(IMSL_TERMINAL, IMSL_IOPT_NEED_3_INTERVALS);
	goto RETURN;
    } else if (iopt == 2 && n_intervals < 2) {
        imsl_e1sti(1, n_intervals);
        imsl_e1sti(2, iopt);
        imsl_ermes(IMSL_TERMINAL, IMSL_IOPT_NEED_2_INTERVALS);
	goto RETURN;
    } else if (n_intervals < 1) {
        imsl_e1sti(1, n_intervals);
        imsl_ermes(IMSL_TERMINAL, IMSL_N_INTERVALS_LT_1);
	goto RETURN;
    }

    if (!user_table) {
        lv_table = (Mfloat *) imsl_malloc (n_intervals*sizeof(Mfloat));
        if (!lv_table){
            imsl_e1sti (1,n_intervals);
            imsl_e1stl (1,"n_intervals");
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto RETURN;
        }
    }

    if (iopt == 2 || iopt == 3) {
        if (n_intervals < 2) {
                imsl_e1sti (1, n_intervals);
                imsl_e1sti (2, iopt);
                imsl_ermes (IMSL_TERMINAL, IMSL_IOPT_NEED_2_INTERVALS);
                goto RETURN;
	}
        clhw = (div[1]-div[0])/2.0;
    }

    if (iopt == 0 || iopt == 1) {
	div = (Mfloat *) imsl_malloc (n_intervals*sizeof(Mfloat));
        if (!div){
            imsl_e1sti (1,n_intervals);
            imsl_e1stl (1,"n_intervals");
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto RETURN;
        }
    }

    l_owfrq(&nobs, y, &n_intervals, &iopt, &xlo, &xhi, &clhw, div, lv_table);
    if (imsl_n1rty(0) > 3) goto RETURN;

    if (iopt==0 && return_bounds==1) {
	*lv_min = y[l_ismin(nobs,y,1)-1];
	*lv_max = y[l_ismax(nobs,y,1)-1];
    }

RETURN:
    if (imsl_n1rty(0) > 3) {
        if(!user_table) {
            if (lv_table != NULL) imsl_free (lv_table);
        }
        lv_table = NULL;
        if (div !=NULL && !user_div) imsl_free (div);  div = NULL;
    }
    if (iopt == 0 || iopt == 1) {
        /*since div is not returned*/
        if (div !=NULL) imsl_free (div);  
        div =  NULL;
    }

    return (argptr);
}


/*Translated by FOR_C++, v0.1, on 09/09/91 at 08:58:13 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/*  -----------------------------------------------------------------------
    IMSL Name:  OWFRQ/DOWFRQ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 25, 1991

    Purpose:    Tally observations into a one-way frequency table.

    Usage:      CALL OWFRQ (NOBS, X, K, IOPT, XLO, XHI, CLHW, DIV, TABLE)

    Arguments:
       NOBS   - Number of observations.  (Input)
       X      - Vector of length NOBS containing the data.  (Input)
       K      - Number of intervals.  (Input)
       IOPT   - Tallying option.  (Input)
                 IOPT   Action
                   0    Intervals of equal length, determined from the
                        data, are used.  Let XMIN and XMAX be the
                        minimum and maximum values in X, respectively.
                        Then, TABLE(1) is the tally of
                        observations less than or equal to XMIN+(XMAX-
                        XMIN)/K, TABLE(2) is the tally of observations
                        greater than XMIN+(XMAX-XMIN)/K and less than
                        or equal to XMIN+2*(XMAX-XMIN)/K, and so on.
                        TABLE(K) is the tally of observations greater
                        than XMAX-(XMAX-XMIN)/K.
                   1    Intervals of equal length are used just as
                        in the case of IOPT  =  0, except the upper and
                        lower bounds are taken as the user-supplied
                        variables XLO and XHI, instead of the actual
                        minimum and maximum in the data.  Therefore,
                        the first and the last intervals are semi-
                        infinite in length.  K must be greater than 2.
                   2    K - 1 cutpoints are input in DIV.  The tally in
                        TABLE(1) is the number of observations less
                        than or equal to DIV(1).  For I greater than 1
                        and less than K, the tally in TABLE(I) is
                        the number of observations greater than DIV(I-
                        1) and less than or equal to DIV(I).  The
                        tally in TABLE(K) is the
                        number of observations greater than DIV(K - 1).
                        K must be greater than 1.
                   3    Class marks are input in DIV and a constant
                        class half-width is input in CLHW.  The total
                        of the elements in TABLE may be less than
                        NOBS.  The tally in TABLE(I) is the
                        number of observations between DIV(I)-CLHW and
                        DIV(I)+CLHW.
       XLO    - If IOPT = 1, XLO is the lower bound at which to begin
                forming the class intervals.  (Input)
                XLO is used only if IOPT = 1.
       XHI    - If IOPT = 1, XHI is the upper bound to use in forming
                the class intervals.  (Input)
                XHI is used only if IOPT = 1.
       CLHW   - If IOPT = 3, CLHW is the half-width of the class
                intervals.  (Input)
                CLHW is not used if IOPT is not equal to 3.
       DIV    - Vector of varying length and contents depending on
                IOPT.  (Input if IOPT = 2 or 3; output if IOPT = 0 or
                1.)
                The contents of DIV are in ascending order.
                 IOPT   Contents
                   0    DIV is of length K containing interval
                        midpoints.  (DIV is output.)
                   1    DIV is of length K containing interval
                        midpoints.  Since the first and last intervals
                        are semi-infinite in length, DIV(1) contains
                        XLO minus half the interval length, and DIV(K)
                        contains XHI plus half the interval length.
                        (DIV is output.)
                   2    DIV is a vector of length K-1 containing
                        cutpoints.  (DIV is input.)
                   3    DIV is of length K containing classmarks.
                        (DIV is input.)
       TABLE  - Vector of length K containing the counts.  (Output)

    Keywords:   Categorical data; Univariate; Class; Interval; Class
                mark; Class boundary; Cutpoint; Midpoint; Half-width;
                Count; Group; Descriptive statistics; Summary statistics;
                Continuous data

    GAMS:       L2b

    Chapter:    STAT/LIBRARY Basic Statistics

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_owfrq(Mint *nobs, Mfloat x[], Mint *k, Mint *iopt, 
        Mfloat *xlo, Mfloat *xhi, Mfloat *clhw, Mfloat div[],Mfloat table[])
#else
static void l_owfrq(nobs, x, k, iopt, xlo, xhi, clhw, div, table)
	Mint            *nobs;
	Mfloat           x[];
	Mint            *k, *iopt;
	Mfloat          *xlo, *xhi, *clhw, div[], table[];
#endif
{
	Mint             l0, l1, i, imax, imin, j, ner;
	Mfloat           ep2, xlen, xmax, xmin, xx;


	imsl_e1psh("OWFRQ ");
	if (*nobs < 1) {
		imsl_e1sti(1, *nobs);
                imsl_ermes(IMSL_TERMINAL, IMSL_NOBS_LESS_THAN_ONE);
	}
	ner = 2;
        l0 = 0; l1 = 3;
	imsl_c1iarg(*iopt, "IOPT", l0, l1, &ner);
	if (*iopt == 1) {
		if (*xlo >= *xhi) {
			imsl_e1str(1, *xlo);
			imsl_e1str(2, *xhi);
                        imsl_ermes(IMSL_TERMINAL, IMSL_XHI_LT_XLO);
		}
	}
	if (*iopt == 2 || *iopt == 3) {
		for (i = 2; i <= *k; i++) {
			if (*iopt == 2 && i == *k)
				goto L_10;
			if (div[i - 1] <= div[i - 2]) {
                                if (*iopt==2) {
                                imsl_e1stl(1, "cutpoints");
				} else {
                                imsl_e1stl(1, "class_marks");
				}
				imsl_e1sti(1, i);
				imsl_e1sti(2, i - 1);
				imsl_e1str(1, div[i - 1]);
				imsl_e1str(2, div[i - 2]);
                                imsl_ermes(IMSL_TERMINAL, IMSL_DIV_NOT_MONOTONIC);
			}
	L_10:
			;
		}
	}
        if (*iopt == 1 && *k < 3) {
                imsl_e1sti(1, *k);
                imsl_ermes(IMSL_TERMINAL, IMSL_IOPT_NEED_3_INTERVALS);
        } else if (*iopt == 2 && *k < 2) {
                imsl_e1sti(1, *k);
                imsl_e1sti(2, *iopt);
                imsl_ermes(IMSL_TERMINAL, IMSL_IOPT_NEED_2_INTERVALS);
        } else if (*k < 1) {
                imsl_e1sti(1, *k);
                imsl_ermes(IMSL_TERMINAL, IMSL_N_INTERVALS_LT_1);
        }
	if (*nobs > 0) {
                l0 = 1;
		if (l_isanan(*nobs, x, l0) > 0) {
                        imsl_ermes(IMSL_TERMINAL, IMSL_NAN_NOT_ALLOWED);
		}
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	/* Initialize TABLE */
	sset(*k, 0.0, &table[0], 1);
	if (*iopt == 0) {
		/* Determine interval bounds and length */
                l0 = 1;
		imin = l_ismin(*nobs, &x[0], l0);
		imax = l_ismax(*nobs, &x[0], l0);
		xmin = x[imin - 1];
		xmax = x[imax - 1];
		xlen = (xmax - xmin) / (Mfloat) (*k);
		/*
		 * Determine which interval contains X(I)
		 */
		ep2 = 10.0 * imsl_amach(4) * xlen;
		for (i = 1; i <= *nobs; i++) {
			xx = x[i - 1] - xmin;
			for (j = 1; j <= (*k - 1); j++) {
				if (j == 1) {
					if (xx - xlen <= ep2) {
						table[0] += 1.0;
						goto L_30;
					}
				} else {
					if (xx > (Mfloat) (j - 1) * xlen && xx - (Mfloat) (j) *
					    xlen <= ep2) {
						table[j - 1] += 1.0;
						goto L_30;
					}
				}
			}
			table[*k - 1] += 1.0;
	L_30:
			;
		}
		for (i = 1; i <= *k; i++) {
			div[i - 1] = xmin + xlen * ((Mfloat) (i - 1) + .5);
		}
	} else if (*iopt == 1) {
		xlen = (*xhi - *xlo) / (Mfloat) (*k - 2);
		for (i = 1; i <= *nobs; i++) {
			if (x[i - 1] <= *xlo) {
				table[0] += 1.0;
			} else if (x[i - 1] > *xhi) {
				table[*k - 1] += 1.0;
			} else {
				xx = x[i - 1] - *xlo;
				for (j = 2; j <= (*k - 1); j++) {
					if (xx > (Mfloat) (j - 2) * xlen && xx <= (Mfloat) (j -
								  1) * xlen)
						table[j - 1] += 1.0;
				}
			}
		}
		div[0] = *xlo - xlen / 2.0;
		div[*k - 1] = *xhi + xlen / 2.0;
		for (i = 1; i <= (*k - 2); i++) {
			div[i] = *xlo + xlen * ((Mfloat) (i - 1) + .5);
		}
		/* Cutpoints contained in DIV are used */
	} else if (*iopt == 2) {
		for (i = 1; i <= (*k - 1); i++) {
			for (j = 1; j <= *nobs; j++) {
				if (i == 1) {
					if (x[j - 1] <= div[0])
						table[0] += 1.0;
					if (x[j - 1] > div[*k - 2])
						table[*k - 1] += 1.0;
				} else {
					if (x[j - 1] > div[i - 2] && x[j - 1] <= div[i - 1])
						table[i - 1] += 1.0;
				}
			}
		}
	} else {
		/* Classmarks contained in DIV are used */
		for (i = 1; i <= *nobs; i++) {
			for (j = 1; j <= *k; j++) {
				if (fabs(x[i - 1] - div[j - 1]) <= *clhw) {
					table[j - 1] += 1.0;
					goto L_110;
				}
			}
	L_110:
			;
		}
	}
L_9000:
	imsl_e1pop("OWFRQ ");
	return;
}				/* end of function */
/* 
  -----------------------------------------------------------------------
    IMSL Name:  ISANAN (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    May 9, 1991

    Purpose:    Find the smallest index of a single-precision vector
                with corresponding element equal to NaN (not a number).

    Usage:      ISANAN(N, SX, INCX)

    Arguments:
       N      - Length of vector X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX).  INCX must be
                greater than zero.
       ISANAN - The smallest index I such that X(I) equals NaN
                (not a number).  (Output)
                X(I) refers to a specific element of SX.
                See INCX argument description.  ISANAN = 0 means
                no element of X equals NaN.

    GAMS:       N

    Chapter:    STAT/LIBRARY Mathematical Support

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mint l_isanan(Mint n, Mfloat sx[], Mint incx)
#else
static Mint l_isanan(n, sx, incx)
	Mint            n;
	Mfloat          sx[];
	Mint            incx;
#endif
{
	Mint            i, isanan_v, k, ner;

	ner = 1;
	isanan_v = 0;
	if (n >= 1 && incx >= 1) {
		i = 1;
		k = 1;
		/*
		 * Repeat until NaN is found or the vector has been checked
		 */
L_10:
		if (imsl_ifnan(sx[i-1]))
			isanan_v = k;
		i += incx;
		k += 1;
		if (isanan_v == 0 && k <= n)
			goto L_10;
	}
L_9000:
	return (isanan_v);
}				/* end of function */


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
}                               /* end of function */

#ifdef ANSI
static Mint l_ismin(Mint n, Mfloat sx[], Mint incx)
#else
static Mint l_ismin(n, sx, incx)
        Mint             n;
        Mfloat           sx[];
        Mint             incx;
#endif
{
        Mint             i, ismin_v, ix;
        Mfloat           smin;


        ismin_v = 0;
        if (n >= 1) {
                ismin_v = 1;
                if (n != 1) {
                        if (incx != 1) {
                                /* CODE FOR INCREMENT NOT EQUAL TO 1 */
                                ix = 1;
                                smin = sx[0];
                                ix += incx;
                                for (i = 2; i <= n; i++) {
                                        if (sx[ix - 1] < smin) {
                                                ismin_v = i;
                                                smin = sx[ix - 1];
                                        }
                                        ix += incx;
                                }
                        } else {
                                /* CODE FOR INCREMENT EQUAL TO 1 */
                                smin = sx[0];
                                for (i = 2; i <= n; i++) {
                                        if (sx[i - 1] < smin) {
                                                ismin_v = i;
                                                smin = sx[i - 1];
                                        }
                                }
                        }
                }
        }
        return (ismin_v);
}                               /* end of function */
