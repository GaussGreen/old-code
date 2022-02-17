#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK  l_chi_squared_test (Mfloat (*fcn) (Mfloat), Mint, Mint,
                    Mfloat *, va_list);
static void     l_chigf (Mfloat (*fcn) (Mfloat), Mint, Mfloat[], Mfloat[], Mint, Mfloat[], Mint, Mfloat[], Mfloat[], Mfloat[], Mfloat[], Mfloat *, Mfloat *);
#if defined(COMPUTER_HP97C)
static Mfloat   l_gfnin (Mfloat (*fcn) (Mfloat), Mint, Mfloat, Mfloat, Mfloat);
#else
static Mfloat   l_gfnin (Mfloat (*fcn) (Mfloat), Mfloat, Mfloat, Mfloat);
#endif
static Mfloat   l_ssum (Mint, Mfloat[], Mint);
#else
static VA_LIST_HACK  l_chi_squared_test ();
static void     l_chigf ();
static Mfloat   l_gfnin ();
static Mfloat   l_ssum ();
#endif

static Mfloat   lv_chi_squared_statistics;


#ifdef ANSI
Mfloat          imsl_f_chi_squared_test (Mfloat (*cdf) (Mfloat), Mint n_observations,
                    Mint n_catagories, Mfloat *x,...)
#else
Mfloat          imsl_f_chi_squared_test (cdf, n_observations, n_catagories, x, va_alist)
Mfloat          (*cdf) ();
Mint            n_observations, n_catagories;
Mfloat         *x;
va_dcl
#endif
{
    va_list         argptr;
    VA_START (argptr, x);

    E1PSH ("imsl_f_chi_squared_test", "imsl_d_chi_squared_test");
    lv_chi_squared_statistics = imsl_amach (6);
    IMSL_CALL(l_chi_squared_test (cdf, n_observations, n_catagories, x, argptr));
    va_end (argptr);

    E1POP ("imsl_f_chi_squared_test", "imsl_d_chi_squared_test");

    return lv_chi_squared_statistics;
}
#ifdef ANSI
static VA_LIST_HACK  l_chi_squared_test (Mfloat (*cdf) (Mfloat), Mint n_observations, Mint n_categories, Mfloat *x, va_list argptr)
#else
static VA_LIST_HACK  l_chi_squared_test (cdf, n_observations, n_categories, x, argptr)
Mfloat          (*cdf) ();
Mint            n_observations, n_categories;
Mfloat         *x;
va_list         argptr;
#endif
{
    Mfloat         *frequencies = NULL;
    Mfloat          bounds[2];
    Mint            n_parameters = 0;
    Mfloat        **cutpoints;
    Mfloat         *cutpoints_ptr = NULL;
    Mfloat        **cell_counts;
    Mfloat         *cell_counts_ptr = NULL;
    Mfloat        **cell_chi_squared;
    Mfloat         *cell_chi_squared_ptr = NULL;
    Mfloat        **cell_expected;
    Mfloat         *cell_expected_ptr = NULL;
    Mfloat         *degrees_of_freedom = NULL;
    Mfloat         *chi_squared = NULL;
    Mint            n_cutpoints = n_categories;

    Mint             code = 1, arg_number = 4;
    Mint       return_cutpoints = 0;
    Mint       return_cell_expected = 0;
    Mint       return_cell_counts = 0;
    Mint       return_cell_chi_squared = 0;
    Mint       return_degrees_of_freedom = 0;
    Mint       return_chi_squared = 0;
    Mint       user_cutpoints = 0;
    Mint       user_cell_expected = 0;
    Mint       user_cell_counts = 0;
    Mint       user_cell_chi_squared = 0;
    Mint       equal_cutpoints = 1;
    Mint       error = 0;
    Mfloat          frequency_array[1];
    Mfloat         *chi_squared_ptr = NULL;
    Mfloat          p, df;


    bounds[0] = F_ZERO;
    bounds[1] = F_ZERO;
    frequency_array[0] = -F_ONE;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	++arg_number;
	switch (code) {
	    case IMSL_FREQUENCIES:
		frequencies = va_arg (argptr, Mfloat *);
		if (!frequencies) {
		    imsl_e1stl (1, "frequencies");
		    imsl_e1stl (2, "IMSL_FREQUENCIES_USER");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		    error = 1;
		}
		++arg_number;
		break;
	    case IMSL_BOUNDS:
		bounds[0] = (Mfloat) va_arg (argptr, Mdouble);
		++arg_number;
		bounds[1] = (Mfloat) va_arg (argptr, Mdouble);
		++arg_number;
		break;
	    case IMSL_BOUNDS_ADR:
		bounds[0] = *(va_arg (argptr, Mfloat *));
		++arg_number;
		bounds[1] = *(va_arg (argptr, Mfloat *));
		++arg_number;
		break;
	    case IMSL_N_PARAMETERS_ESTIMATED:
		n_parameters = va_arg (argptr, Mint);
		++arg_number;
		break;

	    case IMSL_CUTPOINTS:
		cutpoints = va_arg (argptr, Mfloat **);
		cutpoints_ptr = *cutpoints;
		return_cutpoints = 1;
		user_cutpoints = 0;
		++arg_number;
		break;
	    case IMSL_CUTPOINTS_USER:
		cutpoints_ptr = va_arg (argptr, Mfloat *);
		if (!cutpoints_ptr) {
		    imsl_e1stl (1, "cutpoints");
		    imsl_e1stl (2, "IMSL_CUTPOINTS_USER");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		    error = 1;
		}
		return_cutpoints = 1;
		user_cutpoints = 1;
		equal_cutpoints = 0;
		n_cutpoints = n_categories;
		++arg_number;
		break;
	    case IMSL_CUTPOINTS_EQUAL:
		equal_cutpoints = 1;
		break;

	    case IMSL_CELL_COUNTS:
		cell_counts = va_arg (argptr, Mfloat **);
		cell_counts_ptr = *cell_counts;
		user_cell_counts = 0;
		return_cell_counts = 1;
		++arg_number;
		break;
	    case IMSL_CELL_COUNTS_USER:
		cell_counts_ptr = va_arg (argptr, Mfloat *);
		if (!cell_counts_ptr) {
		    imsl_e1stl (1, "cell_counts");
		    imsl_e1stl (2, "IMSL_CELL_COUNTS_USER");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		    error = 1;
		}
		return_cell_counts = 1;
		user_cell_counts = 1;
		++arg_number;
		break;

	    case IMSL_CELL_EXPECTED:
		cell_expected = va_arg (argptr, Mfloat **);
		cell_expected_ptr = *cell_expected;
		return_cell_expected = 1;
		user_cell_expected = 0;
		++arg_number;
		break;
	    case IMSL_CELL_EXPECTED_USER:
		cell_expected_ptr = va_arg (argptr, Mfloat *);
		if (!cell_expected_ptr) {
		    imsl_e1stl (1, "cell_expected");
		    imsl_e1stl (2, "IMSL_EXPECTED_USER");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		    error = 1;
		}
		return_cell_expected = 1;
		user_cell_expected = 1;
		++arg_number;
		break;

	    case IMSL_CELL_CHI_SQUARED:
		cell_chi_squared = va_arg (argptr, Mfloat **);
		cell_chi_squared_ptr = *cell_chi_squared;
		return_cell_chi_squared = 1;
		user_cell_chi_squared = 0;
		++arg_number;
		break;
	    case IMSL_CELL_CHI_SQUARED_USER:
		cell_chi_squared_ptr = va_arg (argptr, Mfloat *);
		if (!cell_chi_squared_ptr) {
		    imsl_e1stl (1, "cell_chi_squared");
		    imsl_e1stl (2, "IMSL_CELL_CHI_SQUARED_USER");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		    error = 1;
		}
		return_cell_chi_squared = 1;
		user_cell_chi_squared = 1;
		++arg_number;
		break;

	    case IMSL_CHI_SQUARED:
		chi_squared = va_arg (argptr, Mfloat *);
		if (!chi_squared) {
		    imsl_e1stl (1, "cell_chi_squared");
		    imsl_e1stl (2, "IMSL_CHI_SQUARED");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		    error = 1;
		}
		return_chi_squared = 1;
		++arg_number;
		break;
	    case IMSL_DEGREES_OF_FREEDOM:
		degrees_of_freedom = va_arg (argptr, Mfloat *);
		if (!degrees_of_freedom) {
		    imsl_e1stl (1, "df");
		    imsl_e1stl (2, "IMSL_DEGREES_OF_FREEDOM");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		    error = 1;
		}
		return_degrees_of_freedom = 1;
		++arg_number;
		break;

	    case 0:
		break;
	    default:
		imsl_e1sti (1, code);
		imsl_e1sti (2, arg_number);
		imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
		goto RETURN;
	}
    }

    if (n_parameters < 0) {
	imsl_e1sti (1, n_parameters);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_N_PARAMETERS_GE_ZERO);
    }
    if (n_categories < 2) {
	imsl_e1sti (1, n_categories);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_N_CATEGORIES_GE_TWO);
    }
    if (imsl_n1rty (0))
	goto RETURN;

    if (equal_cutpoints) {
	n_cutpoints = -n_cutpoints;
    }
    if (!frequencies)
	frequencies = frequency_array;

    if (!user_cutpoints)
	cutpoints_ptr = (Mfloat *) imsl_malloc ((n_categories - 1) * sizeof (Mfloat));

    if (!user_cell_expected)
	cell_expected_ptr = (Mfloat *) imsl_malloc (n_categories * sizeof (Mfloat));

    if (!user_cell_counts)
	cell_counts_ptr = (Mfloat *) imsl_malloc (n_categories * sizeof (Mfloat));

    chi_squared_ptr = (Mfloat *) imsl_malloc ((n_categories + 1) * sizeof (Mfloat));

    if (!cell_counts_ptr || !cell_expected_ptr || !cutpoints_ptr || !chi_squared_ptr) {
        imsl_e1sti (1,n_categories);
        imsl_e1stl (1,"n_categories");
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	return_cell_chi_squared = 0;
	return_cell_counts = 0;
	return_cell_expected = 0;
	return_cutpoints = 0;
	goto FREE_SPACE;
    }

    l_chigf (cdf, n_observations, x, frequencies, n_cutpoints, bounds, n_parameters, cutpoints_ptr, cell_counts_ptr, cell_expected_ptr, chi_squared_ptr, &p, &df);

    if ((code = imsl_n1rty (1)) > 3 && code != 6) {
	return_cell_chi_squared = 0;
	return_cell_counts = 0;
	return_cell_expected = 0;
	return_cutpoints = 0;
    }
    else {
	lv_chi_squared_statistics = p;
	if (return_chi_squared)
	    *chi_squared = *(chi_squared_ptr + n_categories);
	if (return_degrees_of_freedom)
	    *degrees_of_freedom = df;
	if (!user_cell_chi_squared && return_cell_chi_squared) {
/*	    realloc (chi_squared_ptr, n_categories); */
	    cell_chi_squared_ptr = chi_squared_ptr;
	}
	else if (user_cell_chi_squared) {
	    scopy (n_categories, chi_squared_ptr, 1, cell_chi_squared_ptr, 1);
	    if (chi_squared_ptr) imsl_free (chi_squared_ptr);
	}
    }


FREE_SPACE:
    if (!return_cell_expected && !user_cell_expected) {
	if (cell_expected_ptr) imsl_free (cell_expected_ptr);
	cell_expected_ptr = NULL;
    }
    else if (!user_cell_expected)
	*cell_expected = cell_expected_ptr;

    if (!return_cutpoints && !user_cutpoints) {
	if (cutpoints_ptr) imsl_free (cutpoints_ptr);
	cutpoints_ptr = NULL;
    }
    else if (!user_cutpoints)
	*cutpoints = cutpoints_ptr;

    if (!return_cell_chi_squared && !user_cell_chi_squared) {
	if (chi_squared_ptr) imsl_free (chi_squared_ptr);
	cell_chi_squared_ptr = NULL;
	chi_squared_ptr = NULL;
    }
    else if (!user_cell_chi_squared)
	*cell_chi_squared = cell_chi_squared_ptr;

    if (!return_cell_counts && !user_cell_counts) {
	if (cell_counts_ptr) imsl_free (cell_counts_ptr);
	cell_counts_ptr = NULL;
    }
    else if (!user_cell_counts)
	*cell_counts = cell_counts_ptr;

RETURN:
    return argptr;
}

/*
  -----------------------------------------------------------------------
    IMSL Name:  CHIGF/DCHIGF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 19, 1986

    Purpose:    Perform a chi-squared goodness-of-fit test.

    Usage:      CALL CHIGF (IDO, CDF, NELM, X, FREQ, NCAT, RNGE, NDFEST,
                            CUTP, COUNTS, EXPECT, CHISQ, P, DF)

    Arguments:
       IDO    - Processing option.  (Input)
                IDO  Action
                 0   This is the only call to CHIGF, and all of the data
                     are input on this call.
                 1   This is the first call to CHIGF, and additional
                     calls to CHIGF will be made.  Initialization and
                     updating for the data in X are performed.
                 2   This is an intermediate call to CHIGF.  Updating for
                     the data in X is performed.
                 3   This is the final call to CHIGF.  Updating for the
                     data in X and wrap-up computations are performed.
                Calls to CHIGF with IDO = 2 or 3 may be intermixed.
                It is permissable for a call with IDO = 2 to follow a
                call with IDO = 3.
       CDF    - User-supplied FUNCTION to compute the cumulative
                distribution function (CDF) at a given value.
                The form is
                CDF(Y), where
                Y      - Value at which the CDF is to be evaluated.
                         (Input)
                CDF    - Value of the CDF at Y.  (Output)
                CDF must be declared EXTERNAL in the calling program.
       NELM   - The absolute value of NELM is the number of data elements
                currently input in X.  (Input)
                NELM may be positive, zero, or negative.  Negative NELM
                means delete the -NELM data elements from the analysis.
       X      - Vector of length IABS(NELM) containing the data elements
                for this call.  (Input)
                If the data element is missing (NaN, not a number) then
                the observation is ignored.
       FREQ   - Vector containing the frequencies.  (Input)
                If the first element of FREQ is -1.0 then all frequencies
                are taken to be 1 and FREQ is of length 1.  Otherwise
                FREQ is of length IABS(NELM) and the elements in FREQ
                contain the frequency of the corresponding observation
                in X.  If the frequency is missing (NaN, not a number)
                (and FREQ(1) is not -1.0) the observation is ignored.
       NCAT   - The absolute value of NCAT is the number of cells into
                which the observations are to be tallied.  (Input)
                If NCAT is negative, then CHIGF chooses the cutpoints in
                CUTP so that the cells are equiprobable in continuous
                distributions.  NCAT should not be negative in discrete
                distributions.  The user must be careful to define
                cutpoints in discrete distributions, since no error
                message can be generated in this situation if NCAT is
                negative.
       RNGE   - Vector of length 2 containing the lower and upper
                endpoints of the range of the distribution, respectively.
                (Input)
                If the lower and upper endpoints are equal, a range on
                the whole real line is used.  If the lower and upper
                endpoints are different, points outside of the range
                are ignored so that distributions conditional on the
                range can be used.  In this case, the point RNGE(1) is
                excluded from the first interval, but the point
                RNGE(2) is included in the last interval.
       NDFEST - Number of parameters estimated in computing the CDF.
                (Input)
       CUTP   - Vector of length IABS(NCAT)-1 containing the cutpoints
                defining the cells.  (Input if NCAT is positive; output
                otherwise)
                IABS(NCAT)-1 cutpoints define the cells to be used.  If
                NCAT is positive then the cutpoints are input by the
                user.  The intervals defined by the cutpoints are such
                that the lower endpoint is not included while the upper
                endpoint is included in the interval.
       COUNTS - Vector of length IABS(NCAT) containing the counts in each
                of the cells.  (Output if IDO = 0 or 1; input/output if
                IDO > 1)
       EXPECT - Vector of length IABS(NCAT) containing the expected
                count in each cell.  (Output if IDO = 0 or 3; not
                referenced otherwise)
       CHISQ  - Vector of length IABS(NCAT)+1 containing the
                contributions to chi-squared.  (Output if IDO=0 or 3;
                not referenced otherwise)
                Elements 1 through IABS(NCAT) contain the contributions
                to chi-squared for the corresponding cell.  Element
                IABS(NCAT)+1 contains the total chi-squared statistic.
       P      - p-value for the chi-squared statistic in
                CHISQ(IABS(NCAT)+1).  (Output)
                This chi-squared statistic has DF degrees of freedom.
       DF     - Degrees of freedom in chi-squared.  (Output)

    Remark:
       Informational errors
       Type Code
         4   4  There are more observations deleted from a cell
                than added.
         4   5  All observations are missing.
         3   6  An expected value is less than 1.
         3   7  An expected value is less than 5.
         4   8  The function CDF is not a cumulative distribution
                function.
         4   9  The probability of the range of the distribution
                is not positive.
         4   10 An error has occurred when inverting the cumulative
                distribution function.  This function must be continuous
                and defined over the whole real line.  If all else fails,
                you must specify the cutpoints (i.e., NCAT must be
                positive).

    Keywords:   Chi-squared test for normality; Hypothesis test;
                Cumulative distribution function; CDF; P-value;
                Frequency; Tally; Cell; Count; Cutpoint; Nonparametric
                statistics; Goodness of fit

    GAMS:       L4a1c

    Chapter:    STAT/LIBRARY Tests of Goodness of Fit and Randomness

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_chigf (Mfloat (*cdf) (Mfloat), Mint nelm, Mfloat x[],
                    Mfloat freq[], Mint ncat, Mfloat rnge[], Mint ndfest,
                    Mfloat cutp[], Mfloat counts[], Mfloat expect[],
                    Mfloat chisq[], Mfloat *p, Mfloat *df)
#else
static void     l_chigf (cdf, nelm, x, freq, ncat, rnge, ndfest, cutp, counts, expect, chisq, p, df)
Mfloat          (*cdf) ();
Mint            nelm;
Mfloat          x[], freq[];
Mint            ncat;
Mfloat          rnge[];
Mint            ndfest;
Mfloat          cutp[], counts[], expect[], chisq[], *p, *df;
#endif
{
    Mint            i, imi, irn, j, m, irnge, abs_ncat, abs_nelm;
    Mfloat          cell1, e, en, frq = F_ONE, pdif, plower, prob, pupper,
                    t1, t2, xx = F_ZERO, xxo, df2, eps;
    imsl_e1psh ("l_chigf");

    eps = F_TEN * imsl_amach (4);
    abs_ncat = (ncat > 0) ? ncat : -ncat;
    sset (abs_ncat, F_ZERO, counts, 1);
 /* Check CUTP and RNGE */
    irn = 0;
    imi = 0;
    t1 = fabs (rnge[1] - rnge[0]);
    t2 = imsl_f_max (fabs (rnge[0]), fabs (rnge[1]));
    if (t2 == F_ZERO)
	t2 = F_ONE;
    if (t1 <= eps * t2) {
	irnge = 0;
	pupper = F_ONE;
	plower = F_ZERO;
	pdif = F_ONE;
    }
    else {
	irnge = 1;
	imsl_e1usr ("ON");
	pupper = (*cdf) (rnge[1]);
	imsl_e1usr ("OFF");
	imsl_e1usr ("ON");
	plower = (*cdf) (rnge[0]);
	imsl_e1usr ("OFF");
	if (pupper > F_ONE + eps || plower < F_ZERO - eps) {
	    imsl_e1str (1, plower);
	    imsl_e1str (2, pupper);

/*			imsl_ermes(4, 8, "CDF(LOWER_BOUND) = %(r1) and CDF(UPPER_BOUND) = %(r2).  The function CDF is not a cumulative distribution function.");
*/
	    imsl_ermes (IMSL_FATAL, IMSL_INCORRECT_CDF_1);
	    goto L_9000;
	}
	pdif = pupper - plower;
	if (pdif <= 0) {
	    imsl_e1str (1, pdif);

/*			imsl_ermes(4, 9, "The probability of the range is %(r1).  This value is not positive.");
*/
	    imsl_ermes (IMSL_FATAL, IMSL_INCORRECT_CDF_2);
	    goto L_9000;
	}
    }
 /* Generate counts */
    abs_nelm = (nelm > 0 ? nelm : -nelm);
    if (freq[0] == -F_ONE) {
	if (nelm > 0) {
	    frq = F_ONE;
	}
	else {
	    frq = -F_ONE;
	}
    }
    for (i = 1; i <= abs_nelm; i++) {
	if (freq[0] != -F_ONE) {
	    frq = freq[i - 1];
	}
	xx = x[i - 1];
	if (irnge) {
	    if (xx > rnge[1] || xx <= rnge[0]) {
		if (irn == 0) {
		    imsl_e1sti (1, i);
		    imsl_e1str (2, xx);

/*					imsl_ermes(6, 11, "Row %(i1) of X contains a value which is out of range.");
*/
		    imsl_ermes (IMSL_WARNING_IMMEDIATE,
			IMSL_X_VALUE_OUT_OF_RANGE);
		    irn = 1;
		}
		goto L_20;
	    }
	}
	if (imsl_ifnan (frq) || imsl_ifnan (xx)) {
	    if (imi == 0) {

/*				imsl_ermes(6, 3, "At least one data element is missing.");
*/
		imsl_ermes (IMSL_WARNING_IMMEDIATE,
		    IMSL_MISSING_DATA_ELEMENT);
		imi = 1;
	    }
	    goto L_20;
	}
	if (ncat > 0) {
	    m = 1;
	    for (j = 1; j <= (ncat - 1); j++) {
		if (xx > cutp[j - 1]) {
		    m += 1;
		}
	    }
	}
	else {
	    imsl_e1usr ("ON");
	    xx = (*cdf) (xx);
	    imsl_e1usr ("OFF");
	    if (xx < plower - eps || xx > pupper + eps) {
		imsl_e1str (1, plower);
		imsl_e1str (2, pupper);
		imsl_e1str (3, xx);
		imsl_e1sti (1, i-1);

/*				imsl_ermes(4, 8, "CDF(LOWER_BOUND) = %(r1) and CDF(UPPER_BOUND) = %(r2), but CDF(X(%(i1))) = %(r3).  The function CDF is not a cumulative distribution function.");
*/
		imsl_ermes (IMSL_FATAL, IMSL_INCORRECT_CDF_3);
		goto L_9000;
	    }
	    prob = (xx - plower) / pdif;
	    m = abs_ncat * prob + F_ONE - eps;
	    if (m == 0)
		m = 1;
	}
	counts[m - 1] += frq;
	if (counts[m - 1] < -100.0 * imsl_amach (4)) {

/*			imsl_ermes(4, 4, "There are more observations deleted from the cell than added.");
*/
	    imsl_ermes (IMSL_FATAL,
		IMSL_TOO_MANY_CELL_DELETIONS);

	    goto L_9000;
	}
L_20:
	;
    }
 /*
  * Determine the cutpoints if NCAT is negative
  */
    if (ncat < 0) {
	xx = F_ZERO;
	for (i = 1; i <= (abs_ncat - 1); i++) {
	    *p = i;
	    *p = (*p / abs_ncat) * pdif + plower;
#if defined(COMPUTER_HP97C)
	    cutp[i - 1] = l_gfnin (cdf, 1, *p, 0.0001, xx);
#else
	    cutp[i - 1] = l_gfnin (cdf, *p, 0.0001, xx);
#endif
	    if (imsl_n1rty (1) == 4) {
/*				    imsl_ermes(4, 10, "An error occurred when inverting the cumulative distribution function.  This function must be continuous and defined on the whole real line.  If all else fails, specify the cutpoints.");
*/
		imsl_ermes (IMSL_FATAL,
		    IMSL_INCORRECT_CDF_5);
		goto L_9000;
	    }
	    xx = cutp[i - 1];
	    if (i != 1) {
		t1 = fabs (cutp[i - 1] - cutp[i - 2]);
		if (t1 < eps)
		    cutp[i - 1] = cutp[i - 2];
	    }
	}
    }
 /* Determine expected values of COUNTS */
    en = l_ssum (abs_ncat, counts, 1);
    if (en <= 0) {

/*			imsl_ermes(4, 5, "All observations are missing (NaN, not a number) values.");
*/
	imsl_ermes (IMSL_FATAL,
	    IMSL_ALL_OBSERVATIONS_MISSING);
	goto L_9000;
    }
    *df = abs_ncat;
    xxo = plower;
    for (i = 1; i <= abs_ncat; i++) {
	if (i < abs_ncat) {
	    imsl_e1usr ("ON");
	    xx = (*cdf) (cutp[i - 1]);
	    imsl_e1usr ("OFF");
	    if (xx < plower - eps || xx > pupper + eps) {
		imsl_e1str (1, plower);
		imsl_e1str (2, pupper);
		imsl_e1str (3, xx);
		imsl_e1sti (4, i-1);
/*				        imsl_ermes(4, 8, "CDF(LOWER_BOUND) = %(r1) and CDF(UPPER_BOUND) = %(r2), but CDF(CUT_POINT(%(i4))) = %(r3).  The function CDF is not a cumulative distribution fuction.");
*/
		imsl_ermes (IMSL_FATAL,
		    IMSL_INCORRECT_CDF_4);
		goto L_9000;
	    }
	    e = (xx - xxo) / pdif;
	}
	else {
	    e = (pupper - xxo) / pdif;
	}
	xxo = xx;
	e *= en;
	expect[i - 1] = e;
	if (e < F_FIVE) {

/*				imsl_ermes(3, 7, "An expected value is less than 5.");
*/
	    imsl_ermes (IMSL_WARNING,
		IMSL_EXPECTED_VAL_LESS_THAN_5);
	}
	if (e < F_ONE) {

/*				imsl_ermes(3, 6, "An expected value is less than 1.");
*/
	    imsl_ermes (IMSL_WARNING,
		IMSL_EXPECTED_VAL_LESS_THAN_1);
	}
	cell1 = counts[i - 1] - e;
    /* Find CHISQ */
	if (e > F_ZERO) {
	    chisq[i - 1] = cell1 * cell1 / e;
	}
	else {
	    *df -= F_ONE;
	    chisq[i - 1] = F_ZERO;
	}
    }
    chisq[abs_ncat] = l_ssum (abs_ncat, chisq, 1);
    *df += -1 - ndfest;
 /* Find P */
    if (*df > F_HALF) {
	df2 = *df;
	*p = F_ONE - imsl_f_chi_squared_cdf (chisq[abs_ncat], df2);
    }
    else {
	*p = imsl_amach (6);
    }
L_9000:
    imsl_e1pop ("l_chigf");
    return;
}	       /* end of function */

/*
  -----------------------------------------------------------------------
    IMSL Name:  GFNIN/DGFNIN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 11, 1988

    Purpose:    Evaluate the inverse of a general continuous cumulative
                distribution function given in a subprogram.

    Usage:      GFNIN(F, P, EPS, GUESS)

    Arguments:
       F      - User-supplied FUNCTION to be inverted.  F must be
                continuous and strictly monotone.
                The form is
                F(X), where
                X      - The argument to the function.  (Input)
                F      - The value of the function at X.  (Output)
                F must be declared EXTERNAL in the calling program.
       P      - The point at which the inverse of F is desired.  (Input)
       EPS    - Convergence criterion.  (Input)
                When the relative change in GFNIN from one iteration to
                the next is less than EPS, convergence is assumed.
                A common value for EPS is 0.0001.  Another common value
                is 100 times the machine epsilon.
       GUESS  - An initial estimate of the inverse of F at P.  (Input)
       GFNIN  - The inverse of the function F at the point P.  (Output)
                F(GFNIN) is "close" to P.

    Remarks:
    1. Informational errors
       Type Code
         4   1  After 100 attempts, a bound for the inverse can not be
                determined.  Try again with a different initial estimate.
         4   2  No unique inverse exists.
         4   3  Over 100 iterations have occurred without convergence.
                Convergence is assumed.

    2. The function to be inverted need not be a distribution function;
       it can be any continuous, monotonic function.

    Keywords:   Percentage point; Probability distribution; Continuous
                random variable; Percentile; Fractile

    GAMS:       L5a2; F1a2

    Chapter:    STAT/LIBRARY Probability Distribution Functions and
                             Inverses

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static Mfloat   l_gfnin (Mfloat (*f) (Mfloat), Mint idum, Mfloat p, Mfloat eps, Mfloat guess)
#else
static Mfloat   l_gfnin (Mfloat (*f) (Mfloat), Mfloat p, Mfloat eps, Mfloat guess)
#endif /* COMPUTER_HP97C */
#else
static Mfloat   l_gfnin (f, p, eps, guess)
Mfloat          (*f) (), p, eps, guess;
#endif
{
    short int       ibisec;
    Mint            iter;
    Mfloat          delta, f1, f2, f3, fd, imsl_gfnin_v, slope, x1, x2, x3,
                    xd, xm;
    imsl_gfnin_v = imsl_amach (6);
    x1 = guess;
    imsl_e1usr ("ON");
    f1 = (*f) (x1) - p;
    imsl_e1usr ("OFF");
    if (f1 == F_ZERO) {
	x2 = x1;
	goto L_30;
    }
    if (fabs (guess) >= F_ONE) {
	x2 = guess * 1.05;
	xd = x2 - x1;
    }
    else {
	x2 = guess + 0.05;
	xd = 0.05;
    }
    imsl_e1usr ("ON");
    f2 = (*f) (x2) - p;
    imsl_e1usr ("OFF");
 /*
  * Need to bracket estimate and make sure that it is feasible (i.e., that a
  * solution exists.)
  */
    slope = imsl_f_max (.01, (f2 - f1) / xd);
    delta = -f1 / slope;
    iter = 0;
L_10:
    delta *= F_TWO;
    iter += 1;
    if (iter > 100) {

/*		imsl_ermes(4, 1, "After 100 attempts, a bound for the inverse cannot be determined.  Try again with a different initial estimate.");
*/
	imsl_ermes (IMSL_FATAL, IMSL_NO_BOUND_AFTER_100_TRYS);
	imsl_gfnin_v = imsl_amach (6);
	goto L_9000;
    }
    x2 = x1 + delta;
    imsl_e1usr ("ON");
    f2 = (*f) (x2) - p;
    imsl_e1usr ("OFF");
    if (f1 * f2 >= F_ZERO) {
	x1 = x2;
	goto L_10;
    }
 /* Use regula falsi to get the estimate */
    ibisec = 0;
    for (iter = 1; iter <= 100; iter++) {
	xm = (x1 + x2) / F_TWO;
	fd = f2 - f1;
	xd = x2 - x1;
	if (xd != F_ZERO && fd == F_ZERO) {
	    imsl_e1str (1, x1);
	    imsl_e1str (2, x2);
	    imsl_e1str (3, f1);

/*			imsl_ermes(4, 2, "F(%(r1)) = F(%(r2)) = %(r3).  No unique inverse exists.");
*/
	    imsl_ermes (IMSL_FATAL, IMSL_NO_UNIQUE_INVERSE_EXISTS);
	    goto L_9000;
	}
	if (fabs (xm) != F_ZERO) {
	    if (fabs (xd / xm) < eps)
		goto L_30;
	}
	else {
	    if (fabs (xd) < eps)
		goto L_30;
	}
    /*
     * X3 is the intersection of the secant line through(X1,F1), (X2,F2) and
     * the X axis
     */
	if (ibisec) {
	    x3 = xm;
	}
	else {
	    x3 = x2 - f2 * xd / fd;
	}
	ibisec = 0;
	imsl_e1usr ("ON");
	f3 = (*f) (x3) - p;
	imsl_e1usr ("OFF");
	if (f3 * f2 <= F_ZERO) {
	/* Root was trapped, so use regula falsi */
	    x1 = x2;
	    f1 = f2;
	    x2 = x3;
	    f2 = f3;
	}
	else {
	/*
	 * Root was not trapped use Illinois modification
	 */
	    x2 = x3;
	    f2 = f3;
	    f1 /= F_TWO;
	    if (fabs (f2) > fabs (f1)) {
	    /* Use bisection */
		f1 *= F_TWO;
		ibisec = 1;
	    }
	}
    }

/*	imsl_ermes(4, 3, "Over 100 iterations have occurred without convergence.  Convergence is assumed.");
*/
    imsl_ermes (IMSL_FATAL, IMSL_CONVERGENCE_ASSUMED);
L_30:
    imsl_gfnin_v = (x1 + x2) / F_TWO;

L_9000:
    return imsl_gfnin_v;
}	       /* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  SSUM (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Sum the values of a single precision vector.

    Usage:      SSUM(N, SX, INCX)

    Arguments:
       N      - Length of vectors X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than 0.
       SSUM   - Single precision sum from I=1 to N of X(I).  (Output)
                X(I) refers to a specific element of SX.

    Keyword:    Level 1 BLAS

    GAMS:       D1a

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat   l_ssum (Mint n, Mfloat sx[], Mint incx)
#else
static Mfloat   l_ssum (n, sx, incx)
Mint            n;
Mfloat          sx[];
Mint            incx;
#endif
{
    Mint            _d_l, _d_m, _do0, _do1, i, nincx;
    Mfloat          ssum_v;


    ssum_v = F_ZERO;
    if (n > 0) {
	if (incx != 1) {
	/* CODE FOR INCREMENT NOT EQUAL TO 1 */
	    nincx = n * incx;
	    for (i = 1, _do0 = DOCNT (1, nincx, _do1 = incx); _do0 > 0; i += _do1, _do0--) {
		ssum_v += sx[i - 1];
	    }
	}
	else {
	    for (i = 1; i <= n; i++) {
		ssum_v += sx[i - 1];
	    }
	}
    }
    return (ssum_v);
}	       /* end of function */
