#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#define     SIGN(A,B) ((B) < (0) ? (-fabs(A)) : (fabs(A)))
#define	    IABS(A)	(int)fabs((double)A)

static VA_LIST_HACK	PROTO(l_regression,(Mint, Mint, Mfloat*, Mfloat[],
		    va_list));
static void	PROTO(l_r2se,(Mint, Mfloat[], Mint, Mfloat*, Mint, Mint, 
		    Mfloat[], Mfloat*, Mfloat*, Mfloat[], Mint, Mfloat*,
		    Mint*, Mfloat[]));
static void	PROTO(l_r3se,(Mint, Mfloat[], Mint, Mfloat*, Mint, Mint,
		    Mfloat[], Mfloat*, Mfloat*, Mfloat[], Mint, Mfloat*,
		    Mint*, Mfloat[], Mfloat[], Mfloat[], Mfloat[], Mfloat[]));
static void	PROTO(l_r3ivn,(Mint, Mint, Mint, Mint, Mint, Mint, Mint*,
		    Mint, Mint*, Mint, Mint, Mint, Mfloat, Mint, Mint, Mint));
static void	PROTO(l_girts,(Mint, Mfloat[], Mint, Mint, Mfloat[], Mint, 
		    Mint, Mint*, Mfloat[], Mint, Mfloat[], Mint));
static Mint     PROTO(l_i1nan,(Mint , Mfloat *, Mint));

static Mfloat   *lv_coef_vector;
static Mfloat   lv_tolerance;

#ifdef ANSI
Mfloat *imsl_f_regression(Mint n_observations, Mint n_independent,
        Mfloat *x, Mfloat y[], ...)
#else
Mfloat *imsl_f_regression(n_observations, n_independent, x, y, va_alist)
    Mint        n_observations, n_independent;
    Mfloat      *x, y[];
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr,y);

    E1PSH("imsl_f_regression","imsl_d_regression");

    lv_coef_vector = NULL;
    IMSL_CALL(l_regression(n_observations, n_independent, x, y, argptr));
    va_end(argptr);

    E1POP("imsl_f_regression","imsl_d_regression");

    return lv_coef_vector;
}

#ifdef ANSI
static VA_LIST_HACK l_regression(Mint n_observations, Mint n_independent,
               Mfloat *x, Mfloat y[], va_list argptr)
#else
static VA_LIST_HACK l_regression(n_observations, n_independent, x, y, argptr)
    Mint        n_observations, n_independent;
    Mfloat      *x, y[];
    va_list     argptr;
#endif
{
    Mint            x_col_dim     = n_independent;
    Mint            col_dim_covb  = 0;
    Mfloat	  **covb_ptr      = NULL;
    Mfloat	   *covb	  = NULL;
    Mint	   *rank          = NULL;
    Mint	    intercept     = 1;
    Mfloat	  **anova_ptr     = NULL;
    Mfloat	   *anova	  = NULL;
    Mfloat	   *x_mean	  = NULL;
    Mfloat	  **x_mean_ptr	  = NULL;
    Mfloat	   *residual	  = NULL;
    Mfloat	  **residual_ptr  = NULL;

    Mint		code		    = 1;
    Mint	arg_number	    = 4;
    Mint	return_covb	    = 0;
    Mint	return_anova	    = 0;
    Mint	user_coef_vector    = 0;
    Mint	user_covb	    = 0;
    Mint	user_anova	    = 0;
    Mint   user_x_mean	    = 0;
    Mint   return_x_mean	    = 0;
    Mint	user_residual	    = 0;
    Mint	return_residual	    = 0;
    Mint   return_rank         = 0;

    Mint        n_coef, n_missing, n_work, i;
    Mfloat      dfe, sst, sse, *wk = NULL, *r = NULL; 
    Mint         error = 0, ner=1;

    code = 1;
    lv_tolerance = 100.0*imsl_amach(4);
    while (code > 0) {
        code = va_arg(argptr, int);
        ++arg_number;
        switch (code) {
            case IMSL_RETURN_USER:
                lv_coef_vector = va_arg(argptr, Mfloat*);
		user_coef_vector = 1;
		++arg_number;
                break;
            case IMSL_X_COL_DIM:
                x_col_dim = va_arg(argptr, Mint);
		++arg_number;
                break;
            case IMSL_COV_COL_DIM:
		col_dim_covb = va_arg(argptr, Mint);
		++arg_number;
		break;
            case IMSL_COEF_COVARIANCES_USER:
                covb    = va_arg(argptr, Mfloat*);
                return_covb = 1;
		user_covb   = 1;
		++arg_number;
                break;
            case IMSL_COEF_COVARIANCES:
                covb_ptr    = va_arg(argptr, Mfloat**);
                return_covb = 1;
		user_covb   = 0;
		++arg_number;
                break;
	    case IMSL_RANK:
		rank = va_arg(argptr, Mint*);
                return_rank = 1;
		++arg_number;
		break;
	    case IMSL_NO_INTERCEPT:
		intercept = 0;
		break;
            case IMSL_ANOVA_TABLE_USER:
                anova    = va_arg(argptr, Mfloat*);
                return_anova = 1;
		user_anova   = 1;
		++arg_number;
                break;
            case IMSL_ANOVA_TABLE:
                anova_ptr  = va_arg(argptr, Mfloat**);
                return_anova = 1;
		user_anova   = 0;
		++arg_number;
                break;
	    case IMSL_TOLERANCE:
		lv_tolerance = (Mfloat) va_arg(argptr, Mdouble);
		break;
	    case IMSL_TOLERANCE_ADR:
		lv_tolerance = *(va_arg(argptr, Mfloat *));
		break;
	    case IMSL_X_MEAN:
		x_mean_ptr = va_arg(argptr, Mfloat**);
		return_x_mean = 1;
		user_x_mean   = 0;
		++arg_number;
		break;
	    case IMSL_X_MEAN_USER:
		x_mean = va_arg(argptr, Mfloat*);
		return_x_mean = 1;
		user_x_mean   = 1;
		++arg_number;
		break;
	    case IMSL_RESIDUAL:
		residual_ptr = va_arg(argptr, Mfloat**);
		return_residual = 1;
		user_residual   = 0;
		++arg_number;
		break;
	    case IMSL_RESIDUAL_USER:
		residual = va_arg(argptr, Mfloat*);
		return_residual = 1;
		user_residual   = 1;
		++arg_number;
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

    if (x == NULL) {
        imsl_e1stl (1, "x");
        imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
        return argptr;
    }

    if (y == NULL) {
        imsl_e1stl (1, "y");
        imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
        return argptr;
    }

    imsl_c1iarg(n_observations, "n_observations", 1, -2, &ner);
    imsl_c1iarg(n_independent, "n_independent", 1, -2, &ner);
    if (imsl_n1rty(0) != 0) return argptr;
    
    if (x_col_dim < n_independent) {
        imsl_e1stl(1, "x");
	imsl_e1sti(1, n_independent);
	imsl_e1sti(2, x_col_dim);
        imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_LESS_COL);
        return argptr;
    }

    n_coef = n_independent + intercept;
    
    if (col_dim_covb && !return_covb) {
/*	(5, 3, "IMSL_COV_COL_DIM must be used with IMSL_COEF_COVARIANCES or IMSL_COEF_COVARIANCES_USER.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_COVARIANCE_SPECIFIERS);
	return argptr;
    }
    else if (!col_dim_covb) {
	col_dim_covb = n_coef;
    }

    r  = (Mfloat *) imsl_malloc (n_coef*n_coef*sizeof(Mfloat));
    n_work = 5*n_independent + 4*intercept + 2;
    if (!r || !(wk = (Mfloat *) imsl_malloc (n_work*sizeof(Mfloat)))) {
        imsl_e1sti (1,n_independent);
        imsl_e1stl (1,"n_independent");
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	error = 1;
    } else if (!user_coef_vector){
	lv_coef_vector = (Mfloat *) imsl_malloc (n_coef*sizeof(Mfloat));
	if (!lv_coef_vector){
            imsl_e1sti (1,n_independent);
            imsl_e1stl (1,"n_independent");
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    error = 1;
	}
    }

    if (!error){
	imsl_f_m1ran(n_observations,x_col_dim,x,x);
	if (!(error = imsl_n1rty(1) > 3 ? 1: 0)) {
	
	    l_r2se(n_observations, y, n_independent, x, n_observations, intercept, lv_coef_vector, &sst, &sse, r, n_coef, &dfe, &n_missing, wk);
           if (return_rank) *rank = n_observations - n_missing-dfe;

	    if ((error = (imsl_n1rty(1)>3 && code!=6) ? 1 : 0)) {
		if (!user_coef_vector) {
		    imsl_free(lv_coef_vector);
		    lv_coef_vector = NULL;
		} 
	    }
	    if (wk) imsl_free (wk); 
            wk = NULL;

	    if (!error){
		Mint    count, j;

		if (return_covb) {
		    Mfloat  s2;
		    if (!user_covb) {
			covb = (Mfloat *) imsl_malloc(n_coef*col_dim_covb*sizeof(Mfloat));
			if (!covb){
/*			    imsl_ermes (5, 2, "Not enough memory to compute the estimated variances and covariances of the estimated coefficients.");
*/
                            imsl_ermes(IMSL_TERMINAL,
			    IMSL_COVARIANCE_MEMORY_REQ);
			    error = 1;
			}
		    }
		    if (!error){
			imsl_c1div(sse,dfe,&s2);
			imsl_rcovb(n_coef, r, n_coef, s2, covb, col_dim_covb);
			error = (imsl_n1rty(1)>3 && code!=6) ? 2: 0;
			if (error && !user_covb) {
			    if (covb) imsl_free(covb);
			    covb = NULL;
			} else if (!user_covb){
			    *covb_ptr = covb;
			    covb      = NULL;
			}
		    }
		}

		if (r) imsl_free(r); r = NULL;
		error = 0;

		if (return_anova){
		    if (!user_anova){
			anova = (Mfloat *) imsl_malloc(15*sizeof(Mfloat));
			if (!anova){
/*			    imsl_ermes (5, 2, "Not enough memory to compute the analysis of variance statistics.");
*/
                            imsl_ermes(IMSL_TERMINAL,
			    IMSL_NO_MEMORY_VARIANCE_STATS);
			    error = 1;
			}
		    }
		    if (!error) {
			Mfloat y_bar = 0.0;
			Mfloat ssr, dfr;

			ssr = imsl_f_max(sst - sse, 0.0);
			dfr = (Mfloat) (n_observations - intercept) - dfe;
			count = 0.0;

			for (i = 0; i < n_observations; i++) {
			    if (!imsl_ifnan(y[i])) {
				y_bar += y[i];
				++count;
			    }
			}
			y_bar /= (Mfloat) count;

			imsl_g1aov(dfr, ssr, dfe, sse, y_bar, anova);
			error = imsl_n1rty(1) > 3 ?  1: 0;
			if (error && !user_anova) {
			    if (anova) imsl_free(anova);
			    anova = NULL;
			} else if (!user_anova) {
			    *anova_ptr = anova;
			    anova      = NULL;
			}
		    }
		}

		error = 0;

		if (return_x_mean){
		    if (!user_x_mean){
			x_mean = (Mfloat *)imsl_malloc(n_independent*sizeof(Mfloat));
			if (!x_mean){
/*			    imsl_ermes (5, 2, "Not enough memory to compute the vector of means for X.");
*/
                            imsl_ermes(IMSL_TERMINAL,
			    IMSL_VECTOR_OF_MEANS_X_MEMORY);
			    error = 1;
			}
		    }
		    if (!error){
			Mint irow = 0;
			Mint iobs;

			for (j = 0; j < n_independent; j++){
			    x_mean[j] = 0.0;
			    count     = 0;
			    for (i = 0, iobs = irow; i < n_observations; i++, iobs++){
				if (!imsl_ifnan(x[iobs])) {
				    x_mean[j] += x[iobs];
				    count++;
				}
			    }
			    if (count > 0) x_mean[j] /= (Mfloat) count;
			    irow = iobs;
			}
			if (!user_x_mean) {
			    *x_mean_ptr = x_mean;
			    x_mean = NULL;
			}
		    }
		}

		error = 0;

		if (return_residual) {
		    if (!user_residual){
			residual = (Mfloat *) imsl_malloc (n_observations*sizeof(Mfloat));
			if (!residual) {
/*			    imsl_ermes (5, 2, "Not enough memory to compute the vector of residuals.");
*/
                            imsl_ermes(IMSL_TERMINAL,
			    IMSL_RESIDUAL_VECTOR_MEMORY);
			    error = 1;
			}
		    }
		    if (!error){
			Mfloat	*coef_ptr;

			coef_ptr = lv_coef_vector + intercept;
			for (i = 0; i < n_observations; i++){
			    if (!imsl_ifnan(y[i])){
				residual[i] = y[i] - imsl_sdot(n_independent, coef_ptr, 1, (x+i), n_observations);
				if (intercept) residual[i] -= *lv_coef_vector;
			    } else {
				residual[i] = imsl_amach(6);
			    }
			}
			if (!user_residual) {
			    *residual_ptr = residual;
			    residual      = NULL;
			}
		    }
		}
	    }
	    imsl_f_m1ran(x_col_dim,n_observations,x,x);
	}
    }
    if (wk) imsl_free(wk);
    if (r)  imsl_free(r);

    if ((imsl_n1rty(0)>3) && (imsl_n1rty(0)!=6)) {
	if (!user_coef_vector) {
	    if (lv_coef_vector != NULL) imsl_free(lv_coef_vector);
	}
        lv_coef_vector = NULL;
    }
    return (argptr);
}
/*
  -----------------------------------------------------------------------
    IMSL Name:  R2SE/DR2SE (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 23, 1989

    Purpose:    Fit a multiple linear regression model using least
                squares.

    Usage:      CALL R2SE (NOBS, Y, NIND, X, LDX, INTCEP, B, SST, SSE,
                           R, LDR, DFE, NRMISS, WK)

    Arguments:  See RLSE.

    Chapters:   STAT/LIBRARY Regression
                MATH/LIBRARY Interpolation and Approximation

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r2se(Mint nobs, Mfloat y[], Mint nind, Mfloat *x, Mint ldx,
            Mint intcep, Mfloat b[], Mfloat *sst, Mfloat *sse, Mfloat r[],
            Mint ldr, Mfloat *dfe, Mint *nrmiss, Mfloat wk[])
#else
static void l_r2se(nobs, y, nind, x, ldx, intcep, b, sst, sse, r, ldr, dfe,
            nrmiss, wk)
	Mint            nobs;
	Mfloat          y[];
	Mint            nind;
	Mfloat          *x;
	Mint            ldx, intcep;
	Mfloat          b[], *sst, *sse, *r;
	Mint            ldr;
	Mfloat         *dfe;
	Mint           *nrmiss;
	Mfloat          wk[];
#endif
{
	Mint            idwk, iwk, ixmax, ixmin, nreg;


	imsl_e1psh("l_r2se");

	if (intcep != 0 && intcep != 1) {
		imsl_e1sti(1, intcep);
/*		(5, 1, "INTCEP = %(i1).  The intercept option, INTCEP, must be equal to 0 or 1.");
NOTE: This message can be triggered only in the debugger*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_INTCEP_VALUE);
	}
	if (intcep + nind <= 0) {
		imsl_e1sti(1, intcep);
		imsl_e1sti(2, nind);

/*		(5, 2, "INTCEP = %(i1) and NIND = %(i2).  There must be at least one regressor in the model.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MODEL_REGRESSOR_REQUIRED);
	} else {
		if (ldr < intcep + nind) {
			imsl_e1sti(1, ldr);
			imsl_e1sti(2, intcep);
			imsl_e1sti(3, nind);

/*			imsl_ermes(5, 6, "LDR = %(i1), INTCEP = %(i2) and NIND = %(i3).  The leading dimension  of R, LDR, must be greater than or equal  to the number of regressors, NIND+INTCEP.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_NEED_LARGER_LDR_VALUE);
		}
	}
	if (nobs <= 0) {
		imsl_e1sti(1, nobs);

/*		imsl_ermes(5, 4, "NOBS = %(i1).  The number of observations, NOBS, must be greater than zero.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NOBS_GE_ZERO);
	} else {
		if (ldx < nobs) {
			imsl_e1sti(1, nobs);
			imsl_e1sti(2, ldx);

/*			imsl_ermes(5, 5, "NOBS = %(i1) and LDX = %(i2).  The leading dimension of X, LDX, must be greater than or equal to the number of observations, NOBS.");
*/
                        imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDX_GE_NOBS);
		}
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;

	nreg = nind + intcep;
	idwk = nind + 2;
	ixmin = idwk + nreg;
	ixmax = ixmin + nreg;
	iwk = ixmax + nreg;
	l_r3se(nobs, y, nind, x, ldx, intcep, b, sst, sse, r, ldr, dfe,
	     nrmiss, wk, &wk[idwk - 1], &wk[ixmin - 1], &wk[ixmax - 1], &wk[iwk - 1]);

L_9000:
	imsl_e1pop("l_r2se");
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  R3SE/DR3SE (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 23, 1989

    Purpose:    Fit a multiple linear regression model using least
                squares.

    Usage:      CALL R3SE (NOBS, Y, NIND, X, LDX, INTCEP, B, SST, SSE,
                           R, LDR, DFE, NRMISS, XWK, D, XMIN, XMAX,
                           WK)

    Arguments:
       NOBS   - Number of observations.  (Input)
       Y      - Vector of length NOBS containing the dependent (response)
                variable.  (Input)
       NIND   - Number of independent (explanatory) variables.  (Input)
       X      - NOBS by NIND matrix containing the independent
                (explanatory) variables.  (Input)
       LDX    - Leading dimension of X exactly as specified in the
                dimension statement in the calling program.  (Input)
       INTCEP - Intercept option.  (Input)
                INTCEP  Action
                   0    An intercept is not in the model.
                   1    An intercept is in the model.
       B      - Vector of length INTCEP + NIND containing a least
                squares solution for the regression coefficients.
                (Output)
                For INTCEP = 0, the fitted value for observation I is
                   B(1)*X(I,1) + B(2)*X(I,2) + ... + B(NIND)*X(I,NIND)
                For INTCEP = 1, the fitted value for observation I is
                   B(1) + B(2)*X(I,1) + ... + B(NIND+1)*X(I,NIND)
       SST    - Total sum of squares.  (Output)
                If INTCEP = 1, the total sum of squares is corrected
                for the mena.
       SSE    - Sum of squares for error.  (Output)
       R      - INTCEP + NIND by INTCEP + NIND upper triangular matrix
                containing the R matrix from a QR decomposition
                of the matrix of regressors.  (Output)
                All of the diagonal element of R are taken to be
                nonnegative.  The rank of the matrix of regressors is
                the number of positive diagonal elements.
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement in the calling program.  (Input)
       DFE    - Degrees of freedom for error.  (Output)
       NRMISS - Number of rows of data containing NaN (not a number)
                for the dependent or independent variables.  (Output)
                If a row of data contains NaN for any of these variables,
                that row is excluded from the computations.
       XWK    - Work vector of length NIND + 1.  (Output)
       D      - Vector of length INTCEP + NIND containing scale
                factors for fast Givens transformations.  (Output)
       XMIN   - Vector of length INTCEP + NIND containing the
                minimum values for each of the regressors.  (Output)
       XMAX   - Vector of length INTCEP + NIND containing the
                maximum values for each of the regressors.  (Output)
       WK     - Work vector of length NIND + INTCEP + 1.  (Output)

    Chapters:   STAT/LIBRARY Regression
                MATH/LIBRARY Interpolation and Approximation

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r3se(Mint nobs, Mfloat y[], Mint nind, Mfloat *x, Mint ldx,
                    Mint intcep, Mfloat b[], Mfloat *sst, Mfloat *sse,
                    Mfloat r[], Mint ldr, Mfloat *dfe, Mint *nrmiss,
                    Mfloat xwk[], Mfloat d[], Mfloat xmin[], Mfloat xmax[],
                    Mfloat wk[])
#else
static void l_r3se(nobs, y, nind, x, ldx, intcep, b, sst, sse, r, ldr, 
                   dfe, nrmiss, xwk, d, xmin, xmax, wk)
	Mint            nobs;
	Mfloat          y[];
	Mint            nind;
	Mfloat         *x;
	Mint            ldx, intcep;
	Mfloat          b[], *sst, *sse, *r;
	Mint            ldr;
	Mfloat         *dfe;
	Mint           *nrmiss;
	Mfloat          xwk[], d[], xmin[], xmax[], wk[];
#endif
{
#define X(I_,J_)	(x+(I_)*(ldx)+(J_))
#define R(I_,J_)	(r+(I_)*(ldr)+(J_))
	Mint            idum[1], iobs, irank, isub, k, ldb;
	Mfloat          ssr, tol;


	imsl_e1psh("l_r3se");
	/* Initialization */
	tol = lv_tolerance;
	ldb = nind + intcep;
	if (intcep == 0) {
		isub = 0;
	} else {
		isub = 1;
	}
	imsl_r2ivn(1, 0, nind + 1, xwk, 1, intcep, -nind, idum, -1, idum,
	      0, 0, isub, tol, b, ldb, r, ldr, d, &irank, dfe, sse,
	      1, nrmiss, xmin, xmax, wk);
	/*
	 * Intermediate steps Adding one row of X at a time
	 */
	for (iobs = 1; iobs <= nobs; iobs++) {
		scopy(nind, X(0, iobs - 1), ldx, xwk, 1);
		xwk[nind] = y[iobs - 1];
		imsl_r2ivn(2, 1, nind + 1, xwk, 1, intcep, -nind, idum, -1, idum,
		      0, 0, isub, tol, b, ldb, r, ldr, d, &irank, dfe, sse,
		      1, nrmiss, xmin, xmax, wk);
	}
	/* Compute SSR before final computations */
	if (intcep == 0) {
		k = 1;
	} else {
		k = 2;
	}
	ssr = imsl_sxyz(nind, &b[k - 1], 1, &d[k - 1], 1, &b[k - 1], 1);
	/* Compute SST */
	*sst = ssr + *sse;
	/* Wrap up computations */
	imsl_r2ivn(3, 0, nind + 1, xwk, 1, intcep, -nind, idum, -1, idum,
	      0, 0, isub, tol, b, ldb, r, ldr, d, &irank, dfe, sse,
	      1, nrmiss, xmin, xmax, wk);

	/* Check for a unique solution */
	if (irank < intcep + nind) {
		imsl_e1sti(1, irank);

/*		imsl_ermes(3, 1, "The model is not full rank.  There is not a unique least squares solution.  The rank of the matrix of regressors is %(i1).");
*/
                imsl_ermes(IMSL_WARNING, IMSL_RANK_DEFICIENT);
	}
	imsl_e1pop("l_r3se");
	return;

#undef X
#undef R
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  R2IVN/DR2IVN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 1, 1985

    Purpose:    Fit a multivariate linear regression model via fast
                Givens transformations.

    Usage:      CALL R2IVN (IDO, NROW, NVAR, X, LDX, INTCEP, IIND,
                            INDIND, IDEP, INDDEP, IFRQ, IWT, ISUB,
                            TOL, B, LDB, R, LDR, D, IRANK, DFE, SCPE,
                            LDSCPE, NRMISS, XMIN, XMAX, WK)

    Arguments:  See RGIVN.

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_r2ivn(Mint ido, Mint nrow, Mint nvar, Mfloat *x, Mint ldx,
                Mint intcep, Mint iind, Mint indind[], Mint idep,
                Mint inddep[], Mint ifrq, Mint iwt, Mint isub, Mfloat tol,
                Mfloat *b, Mint ldb, Mfloat *r, Mint ldr, Mfloat d[],
                Mint *irank, Mfloat *dfe, Mfloat *scpe, Mint ldscpe,
                Mint *nrmiss, Mfloat xmin[], Mfloat xmax[], Mfloat wk[])
#else
void imsl_r2ivn(ido, nrow, nvar, x, ldx, intcep, iind, indind, idep,
                inddep, ifrq, iwt, isub, tol, b, ldb, r, ldr, d, irank,
                dfe, scpe, ldscpe, nrmiss, xmin, xmax, wk)
    Mint            ido, nrow, nvar;
    Mfloat          x[];
    Mint            ldx, intcep, iind, indind[], idep, inddep[], ifrq;
    Mint            iwt, isub;
    Mfloat          tol, *b;
    Mint            ldb;
    Mfloat         *r;
    Mint            ldr;
    Mfloat          d[];
    Mint           *irank;
    Mfloat         *dfe, *scpe;
    Mint            ldscpe, *nrmiss;
    Mfloat          xmin[], xmax[], wk[];
#endif
{
#define B(I_,J_)	(b+(I_)*(ldb)+(J_))
#define R(I_,J_)	(r+(I_)*(ldr)+(J_))
#define SCPE(I_,J_)	(scpe+(I_)*(ldscpe)+(J_))
	Mint           i, i1, idepx, igo, intp1, iobs, irow, j, jdepix, jdepjx, jdepx, k, ldep, ncoef, nconst, ndep, nind, nobs;
	static Mint     icall;
	Mfloat          frq, sd2, sparam[5], sumwt, temp, tolsq, wt;


	imsl_e1psh("imsl_r2ivn");
	/*
	 * Check for terminal errors.
	 */
	l_r3ivn(ido, nrow, nvar, ldx, intcep, iind, indind, idep, inddep,
	      ifrq, iwt, isub, tol, ldb, ldr, ldscpe);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	ndep = IABS(idep);
	nind = IABS(iind);
	ncoef = intcep + nind;
	intp1 = intcep + 1;
	idepx = ncoef + 1;
	tolsq = pow(tol, 2.0);
	if (ido <= 1) {
		/* Initialize ICALL, NRMISS, and DFE. */
		icall = 1;
		*nrmiss = 0;
		*dfe = 0.0;
		/*
		 * Initialize elements of R, B, D, and SCPE.
		 */
		for (i = 1; i <= ncoef; i++) {
			sset(ncoef, 0.0, R(i - 1, 0), 1);
		}
		for (i = 1; i <= ndep; i++) {
			sset(ncoef, 0.0, B(i - 1, 0), 1);
		}
		sset(ncoef, 1.0, d, 1);
		for (i = 1; i <= ndep; i++) {
			sset(i, 0.0, SCPE(i - 1, 0), 1);
		}
		/*
		 * Initialize XMIN and XMAX to not-a-number.
		 */
		sset(ncoef, imsl_amach(6), xmin, 1);
		sset(ncoef, imsl_amach(6), xmax, 1);
	} else {
		icall += 1;
	}
	if (nrow < 0) {
		/*
		 * Rows of data are to be deleted from from analysis.
		 */
		nobs = -nrow;
		irow = -1;
	} else {
		/*
		 * Rows of data are to be added to analysis.
		 */
		nobs = nrow;
		irow = 1;
	}
	if (isub == 0) {
		i1 = 1;
	} else {
		i1 = 2;
	}
	for (iobs = 1; iobs <= nobs; iobs++) {
		/*
		 * Check frequency and weight.
		 */
		imsl_c1wfr(ido, icall, x, ldx, iobs, irow, ifrq, iwt, imsl_amach(6),
		      nrmiss, &frq, &wt, &igo);
		if (igo == 3)
			goto L_9000;
		if (igo == 2)
			goto L_130;
		if (igo == 1)
			goto L_130;
		/*
		 * Gather independent and dependent variables in row IOBS of
		 * X into WK. Check for missing value code not-a-number.
		 */
		if (intcep == 1) {
			wk[0] = 1.0;
		}
		for (i = 1; i <= iind; i++) {
			wk[intcep + i - 1] = x[(indind[i - 1] - 1) * ldx + iobs -
					       1];
		}
		for (i = 1; i <= (-iind); i++) {
			wk[intcep + i - 1] = x[(i - 1) * ldx + iobs - 1];
		}
		if (l_i1nan(nind, &wk[intp1 - 1], 1) > 0) {
			*nrmiss += irow;
			goto L_130;
		}
		jdepx = idepx;
		for (i = 1; i <= idep; i++) {
			wk[jdepx - 1] = x[(inddep[i - 1] - 1) * ldx + iobs - 1];
			jdepx += 1;
		}
		for (i = idep + 1; i <= 0; i++) {
			wk[jdepx - 1] = x[(nvar + i - 1) * ldx + iobs - 1];
			jdepx += 1;
		}
		if (l_i1nan(ndep, &wk[idepx - 1], 1) > 0) {
			*nrmiss += irow;
			goto L_130;
		}
		/* Update degrees of freedom. */
		*dfe += frq;
		if (irow == 1) {
			if (ncoef > 0) {
				if (imsl_ifnan(xmin[0])) {
					/*
					 * Initialize XMIN and XMAX to first
					 * nonmissing observation.
					 */
					scopy(ncoef, wk, 1, xmin, 1);
					scopy(ncoef, wk, 1, xmax, 1);
				}
			}
			for (i = 1; i <= ncoef; i++) {
				/* Update XMIN and XMAX. */
				temp = wk[i - 1];
				if (temp < xmin[i - 1])
					xmin[i - 1] = temp;
				if (temp > xmax[i - 1])
					xmax[i - 1] = temp;
			}
		} else {
			for (i = intp1; i <= ncoef; i++) {
				temp = wk[i - 1];
				if (temp == xmin[i - 1]) {
					imsl_e1sti(1, i);
					imsl_e1str(1, temp);
					imsl_e1stl(1, "xmin");
/*					imsl_ermes(6, 5, "Downdating is requested, but XMIN(%(i1)) equals the current value (= %(r1)) of the associated regressor variable.  Downdating of XMIN cannot occur.");
*/
                                        imsl_ermes(IMSL_WARNING_IMMEDIATE,
					IMSL_DOWNDATING_REQUESTED);
                                        
				}
				if (temp == xmax[i - 1]) {
					imsl_e1sti(1, i);
					imsl_e1str(1, temp);
					imsl_e1stl(1, "xmax");
/*					imsl_ermes(6, 6, "Downdating is requested, but XMAX(%(i1)) equals the current value (= %(r1)) of the associated regressor variable.  Downdating of XMAX cannot occur.");
*/
                                        imsl_ermes(IMSL_WARNING_IMMEDIATE,
					IMSL_DOWNDATING_REQUESTED);
				}
			}
		}
		if (wt == 0.0)
			goto L_130;
		sd2 = wt * frq;
		if (isub == 1) {
			/* Update sums. */
			sumwt = *R(0, 0);
			*R(0, 0) = sumwt + sd2;
			if (nind > 0) {
				saxpy(nind, sd2, &wk[1], 1, R(1, 0), ldr);
			}
			saxpy(ndep, sd2, &wk[idepx - 1], 1, b, ldb);

			/* Center variables. */
			if (*R(0, 0) == 0.0) {
				goto L_130;
			}
			d[0] = 1.0 / *R(0, 0);
			if (nind > 0) {
				saxpy(nind, -d[0], R(1, 0), ldr, &wk[1], 1);
			}
			saxpy(ndep, -d[0], B(0, 0), ldb, &wk[idepx - 1], 1);
			if (sumwt == 0.0) {
				goto L_130;
			}
			sd2 = sd2 ** R(0, 0) / sumwt;
		}
		/*
		 * Update D, R, and B via fast Givens transformations.
		 */
		for (i = i1; i <= ncoef; i++) {
			imsl_srotmg(&d[i - 1], &sd2, R(i - 1, i - 1), &wk[i - 1], sparam);
			if (ndep > 0)
				imsl_srotm(ndep, B(0, i - 1), ldb, &wk[idepx - 1], 1,
				      sparam);
			if (i == ncoef)
				goto L_100;
			imsl_srotm(ncoef - i, R(i, i - 1), ldr, &wk[i], 1, sparam);
	L_100:
			;
		}
		/* Update upper triangle of SCPE. */
		jdepjx = idepx;
		for (j = 1; j <= ndep; j++) {
			jdepix = idepx;
			for (i = 1; i <= j; i++) {
				*SCPE(j - 1, i - 1) += wk[jdepix - 1] * sd2 * wk[jdepjx - 1];
				jdepix += 1;
			}
			jdepjx += 1;
		}
L_130:
		;
	}

	if (ido == 0 || ido == 3) {
		/*
		 * Remove regressors which are linearly dependent.
		 */
		nconst = 0;
		for (i = 1; i <= ncoef; i++) {
			ldep = 0;
			if (xmin[i - 1] == xmax[i - 1]) {
				if (xmin[i - 1] == 0.0) {
					ldep = 1;
				} else {
					nconst += 1;
					if (nconst > 1)
						ldep = 1;
				}
			} else if (i > intp1) {
				/*
				 * Scale factor for each independent variable
				 * is SUM(X-XBAR)**2. If INTCEP = 0, XBAR is
				 * taken to be zero in the formula above.
				 */
				temp = imsl_sxyz(i - intcep, R(i - 1, intp1 - 1), 1, &d[intp1 - 1], 1, R(i - 1, intp1 - 1), 1);
				if (d[i - 1] * pow(*R(i - 1, i - 1), 2.0) <= tolsq * temp)
					ldep = 1;
			}
			if (ldep == 1) {
				for (j = i + 1; j <= ncoef; j++) {
					imsl_srotmg(&d[j - 1], &d[i - 1], R(j - 1, j - 1),
					       R(j - 1, i - 1), sparam);
					if (ndep > 0)
						imsl_srotm(ndep, B(0, j - 1), ldb, B(0, i - 1), ldb, sparam);
					if (j == ncoef)
						goto L_140;
					imsl_srotm(ncoef - j, R(j, j - 1), ldr, R(j, i - 1),ldr, sparam);
			L_140:
					;
				}
				for (j = 1; j <= ndep; j++) {
					for (k = 1; k <= j; k++) {
						*SCPE(j - 1, k - 1) += *B(k - 1, i - 1) * d[i - 1] *
							*B(j - 1, i - 1);
					}
				}
				sset(ndep, 0.0, B(0, i - 1), ldb);
				sset(ncoef - i + 1, 0.0, R(i - 1, i - 1), ldr);
			}
		}
		/*
		 * Compute B by back substitution.
		 */
		l_girts(ncoef, r, ldr, ndep, b, ldb, 1, irank, b, ldb, r, ldr);
		*dfe -= *irank;
		if (*dfe <= 0.0) {
			imsl_e1str(1, *dfe);

/*			imsl_ermes(6, 4, "DFE = %(r1).  Statistical inference is not possible.  More observations are needed.");
*/
                        imsl_ermes(IMSL_WARNING_IMMEDIATE,
			IMSL_NO_STAT_INFERENCE);
		}
		/*
		 * Set R = SQRT(DIAG(D))*R and then put diag(D) = 1.
		 */
		for (i = 1; i <= ncoef; i++) {
			d[i - 1] = SIGN(sqrt(d[i - 1]), *R(i - 1, i - 1));
		}
		for (i = 1; i <= ncoef; i++) {
			sscal(ncoef - i + 1, d[i - 1], R(i - 1, i - 1), ldr);
		}
		sset(ncoef, 1.0, d, 1);
	}
	/*
	 * Copy upper triangle of SCPE into lower triangle.
	 */
	if (ndep > 1) {
	    for (i = 1; i <= ndep; i++) {
		scopy(ndep - i, SCPE(i, i - 1), ldscpe, SCPE(i - 1, i), 1);
	     }
	}
	/* Exit section */
L_9000:
	imsl_e1pop("imsl_r2ivn");
	return;

#undef B
#undef R
#undef SCPE
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  R3IVN/DR3IVN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 4, 1985

    Purpose:    Check for user errors in calls to RGIVN and R2IVN.

    Usage:      CALL R3IVN (IDO, NROW, NCOL, LDX, INTCEP, IIND, INDIND,
                            IDEP, INDDEP, IFRQ, IWT, ISUB, TOL, LDB,
                            LDR, LDSCPE)

    Arguments:  See RGIVN.

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r3ivn(Mint ido, Mint nrow, Mint ncol, Mint ldx, Mint intcep,
                    Mint iind, Mint *indind, Mint idep, 
                    Mint *inddep, Mint ifrq, Mint iwt, Mint isub, Mfloat tol,
                    Mint ldb, Mint ldr, Mint ldscpe)
#else
static void l_r3ivn(ido, nrow, ncol, ldx, intcep, iind, indind, idep, 
                    inddep, ifrq, iwt, isub, tol, ldb, ldr, ldscpe)
	Mint            ido, nrow, ncol, ldx, intcep, iind, indind[], idep, inddep[], ifrq, iwt, isub;
	Mfloat          tol;
	Mint            ldb, ldr, ldscpe;
#endif
{
	Mint            i, ndep, ner, nind;


	ner = 1;
	imsl_c1iarg(ido, "ido", 0, 3, &ner);
	imsl_c1iarg(ncol, "ncol", 0, -1, &ner);

	imsl_c1dim(0, IABS(nrow), "iabs(nrow)", ldx, "ldx", &ner);
	imsl_c1iarg(intcep, "intcep", 0, 1, &ner);
	imsl_c1ind(0, ifrq, "ifrq", ncol, "ncol", &ner);
	imsl_c1ind(0, iwt, "iwt", ncol, "ncol", &ner);
	if (tol < 0.0 || tol > 1.0) {
		imsl_e1str(1, tol);

/*		imsl_ermes(5, ner, "TOL = %(r1).  TOL must be between 0.0 and 1.0, inclusive.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_VALUE_OF_TOL);
	}
	ner += 1;
	if (intcep == 0) {
		if (isub != 0) {
			imsl_e1sti(1, isub);

/*			imsl_ermes(5, ner, "INTCEP = 0 and ISUB = %(i1).  When INTCEP = 0, ISUB must equal 0.");
*/
                        imsl_ermes(IMSL_TERMINAL, IMSL_ISUB_SHOULD_BE_ZERO);
		}
		ner += 1;
	} else {
		imsl_c1iarg(isub, "isub", 0, 1, &ner);
	}
	nind = IABS(iind);
	ndep = IABS(idep);
	if (ndep > 0) {

		imsl_c1dim(0, intcep + nind, "intcep+iabs(iind)", ldb, "tdb", &ner);
		imsl_c1dim(0, intcep + nind, "*intcep+iabs(iind)", ldr, "tdr", &ner);
		imsl_c1dim(0, ndep, "iabs(idep)", ldscpe, "tdscpe", &ner);
	} else {
		ner += 3;
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	for (i = 1; i <= idep; i++) {
		if (inddep[i - 1] <= 0 || inddep[i - 1] > ncol) {
			imsl_e1sti(1, i);
			imsl_e1sti(2, inddep[i - 1]);
			imsl_e1sti(3, ncol);

/*			imsl_ermes(5, ner, "INDDEP(%(i1)) = %(i2) and NCOL = %(i3).  INDDEP(%(i1)) must be greater than or equal 1 and less than or equal to NCOL.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_NEED_ONE_LE_INDDEP_LE_NCOL);
		}
	}
	ner += 1;
	for (i = 1; i <= iind; i++) {
		if (indind[i - 1] <= 0 || indind[i - 1] > ncol) {
			imsl_e1sti(1, i);
			imsl_e1sti(2, indind[i - 1]);
			imsl_e1sti(3, ncol);

/*			imsl_ermes(5, ner, "INDIND(%(i1)) = %(i2) and NCOL = %(i3).  INDIND(%(i1)) must be greater than or equal 1 and less than or equal to NCOL.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_NEED_ONE_LE_INDIND_LE_NCOL);
		}
	}
L_9000:
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  C1WFR/DC1WFR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 2, 1985

    Purpose:    Check weights and frequencies in statistics routines.

    Usage:      CALL C1WFR (IDO, ICALL, X, LDX, IOBS, IROW, IFRQ, IWT,
                            XMISS, NMISS, FRQ, WT, IGO)

    Arguments:
       IDO    - Processing option.  (Input)
                This is the standard processing option parameter.  It is
                used in this subroutine only to determine the form of an
                error message.
       ICALL  - The number of the current invocation of the routine that
                calls C1WFR.  (Input)
       X      - Matrix containing the data.  (Input)
       LDX    - Row dimension of X exactly as specified in the
                declarations statement in the calling program.  (Input)
       IOBS   - Index of the current row in X.  (Input)
       IROW   - Indicator of whether current observation is being added
                or deleted.  (Input)
                IROW = 1 means addition; IROW = -1 means deletion.
       IFRQ   - Frequency option.  (Input)
                IFRQ=0 means that all frequencies are 1.0.  For positive
                IFRQ, column number IFRQ of X contains the frequencies.
       IWT    - Weighting option.  (Input)
                IWT=0 means that all weights are 1.0.  For positive IWT,
                column IWT of X contains the weights.
       XMISS  - Missing value indicator (not-a-number).  (Input).
       NMISS  - Number of rows of data encountered that contain any
                missing values.  (Input/output)
       FRQ    - The frequency of the current row.  (Output)
                FRQ is negative for rows being deleted and, in that case,
                -FRQ is the frequency.
       WT     - The weight of the current row.  (Output)
       IGO    - Status indicator.  (Output)
                IGO = 0 means the frequency is positive and nonmissing
                and the weight is nonnegative and nonmissing.
                IGO = 1 means the frequency is zero.
                IGO = 2 means the frequency is missing and the weight is
                not negative or else the weight is missing and the
                frequency is positive.  NMISS is updated.
                IGO = 3 means the frequency is negative or else the
                weight is negative and the frequency is either positive
                or missing.  In this case, a type 4 error message is
                issued.

    Remark:
       The manner in which IGO is set in C1WFR is consistent with the
       treatment of missing, negative, or zero weights and/or frequencies
       in STAT/LIBRARY.  In general, a zero frequency means that the
       row is to be eliminated from the analysis; no further processing,
       counting of missing values, or error checking is done on the row.
       Although it is not required that frequencies be integral, the
       logic of their treatment implicitly assumes that they are.
       Weights, on the other hand, are allowed to be continuous; they
       are not thought of as replication factors.  A weight of zero
       results in the row being counted, and updates are made of order
       statistics and of the number of missing values.  A missing value
       for the frequency or a missing value for the weight when the
       frequency is nonzero results in the row being deleted from the
       analysis; but even in that case, if one is nonmissing, it is an
       error for that nonmissing weight or frequency to be negative.
       The following table gives the values of IGO for all possibilities:
                                 FRQ
                      Miss   Neg.    Zero   Pos.
                    ------------------------------
              Miss !   2      3       1      2
          W   Neg. !   3      3       1      3
          T   Zero !   2      3       1      0
              Pos. !   2      3       1      0

        Generally, the routine calling C1WFR would, following the call,
             IF (IGO .EQ. 3) branch to end of routine
             IF (IGO .EQ. 2) treat the row as missing and so branch
                             to the end of the loop handling the row
             IF (IGO .EQ. 1) branch to the end of the loop handling the
                             row
             ELSE            check for other missing values, update
                             counts, update order statistics, etc.,
                             then IF (WT .EQ. 0.0) skip over any
                             processing for which zero weights result
                             in no change.

    Copyright:  1985 by IMSL, Inc.  All rights reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  ----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c1wfr(Mint ido, Mint icall, Mfloat *x, Mint ldx, Mint iobs,
                    Mint irow, Mint ifrq, Mint iwt, Mfloat xmiss, Mint  *nmiss,
                    Mfloat *frq, Mfloat *wt, Mint *igo)
#else
void imsl_c1wfr(ido, icall, x, ldx, iobs, irow, ifrq, iwt, xmiss, nmiss,
                    frq, wt, igo)
    Mint            ido, icall;
    Mfloat          *x;
    Mint            ldx, iobs, irow, ifrq, iwt;
    Mfloat          xmiss;
    Mint           *nmiss;
    Mfloat         *frq, *wt;
    Mint           *igo;
#endif
{

	*igo = 0;
	if (ifrq > 0) {
                *frq = x[(ifrq-1)*ldx+iobs-1];
		if (imsl_ifnan(*frq)) {
			*nmiss += irow;
			*igo = 2;
		} else if (*frq == 0.0) {
			*igo = 1;
			goto L_9000;
		}
	}
	if (iwt > 0) {
                *wt = x[ldx*(iwt-1)+iobs-1];
		if (imsl_ifnan(*wt)) {
			if (*igo != 2) {
				*nmiss += irow;
				*igo = 2;
			}
		}
	}
	if (ifrq > 0) {
		if (!imsl_ifnan(*frq)) {
			if (*frq < 0.0) {
				imsl_e1sti(1, iobs);
				imsl_e1str(1, *frq);
				if (ido > 0) {
					imsl_e1sti(2, icall);

/*					imsl_ermes(4, 2, "The frequency for row %(i1) of X on invocation number %(i2) of this routine is %(r1).  Frequencies must be nonnegative.");
*/
                                        imsl_ermes(IMSL_FATAL,
					IMSL_NONNEG_FREQ_REQUEST_1);
				} else {

/*					imsl_ermes(4, 2, "The frequency for row %(i1) of X is %(r1).  Frequencies must be nonnegative.");
*/
                                        imsl_ermes(IMSL_FATAL,
					IMSL_NONNEG_FREQ_REQUEST_2);
				}
				*igo = 3;
				goto L_9000;
			}
		}
	} else {
		*frq = 1.0;
	}
	if (irow == -1)
		*frq = -*frq;
	if (iwt > 0) {
		if (!imsl_ifnan(*wt)) {
			if (*wt < 0.0) {
				imsl_e1sti(1, iobs);
				imsl_e1str(1, *wt);
				if (ido > 0) {
					imsl_e1sti(2, icall);

/*					imsl_ermes(4, 1, "The weight for row %(i1) of X on invocation number %(i2) of this routine was %(r1).  Weights must be nonnegative.");
*/
                                        imsl_ermes(IMSL_FATAL,
					IMSL_NONNEG_WEIGHT_REQUEST_1);
				} else {

/*					imsl_ermes(4, 1, "The weight for row %(i1) of X was %(r1).  Weights must be nonnegative.");
*/
                                        imsl_ermes(IMSL_FATAL,
					IMSL_NONNEG_WEIGHT_REQUEST_2);
				}
				*igo = 3;
				goto L_9000;
			}
		}
	} else {
		*wt = 1.0;
	}
L_9000:
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  SROTMG (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Construct a modified Givens plane rotation in single
                precision.

    Usage:      CALL SROTMG (SD1, SD2, SX1, SY1, SPARAM)

    Arguments:
       SD1    - Scale factor.  (Input/Output)
                On input, SD1 contains the first scale factor.  On
                output, SD1 contains the updated scale factor.
       SD2    - Scale factor.  (Input/Output)
                On input, SD2 contains the second scale factor.  On
                output, SD2 contains the updated scale factor.
       SX1    - On input, SX1 contains the first component of the vector
                to be rotated.  On output SX1 contains the rotated value
                of the first component .  (Input/Output)
       SY1    - Second component of the vector to be rotated.  (Input)
                Since this component is zeroed by the rotation, it is
                left unchanged in storage.
       SPARAM - Real vector of length 5 which defines the rotation matrix
                H.  (Input/Output)
                See remark.

    Remark:
       SROTMG constructs a modified Givens rotation H and updates the
       scale factors SD1 and SD2 which zero SY1. The transformed value of
       SD1 replaces SD1, i.e.
       On input,  SW1  =  SQRT(SD1)*SX1
                  SZ1  =  SQRT(SD2)*SY1.
       On output, ( C  S ) (SW1) = (C*SW1 + S*SZ1) = (SQRT(SD1)*SX1)
                  (-S  C ) (SZ1) = (      0      ) = (      0       )
       where C and S define a Givens rotation. H takes the form
       SPARAM(1) = -2.0
         SPARAM(2) = UNCHANGED   SPARAM(4) = UNCHANGED
         SPARAM(3) = UNCHANGED   SPARAM(5) = UNCHANGED
       SPARAM(1) = -1.0
         SPARAM(2) = H11         SPARAM(4) = H12
         SPARAM(3) = H21         SPARAM(5) = H22
       SPARAM(1) = 0.0
         SPARAM(2) = UNCHANGED   SPARAM(4) = H12
         SPARAM(3) = H21         SPARAM(5) = UNCHANGED
       SPARAM(1) = 1.0
         SPARAM(2) = H11         SPARAM(4) = UNCHANGED
         SPARAM(3) = UNCHANGED   SPARAM(5) = H22

    Keywords:   Level 1 BLAS; SROTMG

    GAMS:       D1a8

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_srotmg(Mfloat *sd1, Mfloat *sd2, Mfloat *sx1, Mfloat *sy1, Mfloat sparam[])
#else
void imsl_srotmg(sd1, sd2, sx1, sy1, sparam)
	Mfloat        *sd1, *sd2, *sx1, *sy1, sparam[];
#endif
{
	Mint           igo;
	Mfloat         sflag, sh11, sh12, sh21, sh22, sp1, sp2, sq1, sq2, stemp, su;
	static Mfloat  zero = 0.0;
	static Mfloat  one = 1.0;
	static Mfloat  two = 2.0;
	static Mfloat  gam = 4096.0;
	static Mfloat  gamsq = 1.678e7;
	static Mfloat  rgamsq = 5.96e-8;



	if (*sd1 >= zero)
		goto L_10;
	/* GO ZERO-H-D-AND-SX1.. */
	goto L_60;
L_10:
	;
	/* CASE-SD1-NONNEGATIVE */
	sp2 = *sd2 ** sy1;
	if (sp2 != zero)
		goto L_20;
	sflag = -two;
	goto L_270;
	/* REGULAR-CASE.. */
L_20:
	;
	sp1 = *sd1 ** sx1;
	sq2 = sp2 ** sy1;
	sq1 = sp1 ** sx1;

	if (fabs(sq1) <= fabs(sq2))
		goto L_40;
	sh21 = -*sy1 / *sx1;
	sh12 = sp2 / sp1;

	su = one - sh12 * sh21;

	if (su > zero)
		goto L_30;
	/* GO ZERO-H-D-AND-SX1.. */
	goto L_60;
L_30:
	;
	sflag = zero;
	*sd1 /= su;
	*sd2 /= su;
	*sx1 *= su;
	/* GO SCALE-CHECK.. */
	goto L_100;
L_40:
	;
	if (sq2 >= zero)
		goto L_50;
	/* GO ZERO-H-D-AND-SX1.. */
	goto L_60;
L_50:
	;
	sflag = one;
	sh11 = sp1 / sp2;
	sh22 = *sx1 / *sy1;
	su = one + sh11 * sh22;
	stemp = *sd2 / su;
	*sd2 = *sd1 / su;
	*sd1 = stemp;
	*sx1 = *sy1 * su;
	/* GO SCALE-CHECK */
	goto L_100;
	/* PROCEDURE..ZERO-H-D-AND-SX1.. */
L_60:
	;
	sflag = -one;
	sh11 = zero;
	sh12 = zero;
	sh21 = zero;
	sh22 = zero;

	*sd1 = zero;
	*sd2 = zero;
	*sx1 = zero;
	/* RETURN.. */
	goto L_230;
	/* PROCEDURE..FIX-H.. */
L_70:
	;
	if (sflag < zero)
		goto L_90;

	if (sflag != zero)
		goto L_80;
	sh11 = one;
	sh22 = one;
	sflag = -one;
	goto L_90;
L_80:
	;
	sh21 = -one;
	sh12 = one;
	sflag = -one;
L_90:
	;
	if (igo == 130) {
		goto L_130;
	} else if (igo == 160) {
		goto L_160;
	} else if (igo == 190) {
		goto L_190;
	} else if (igo == 220) {
		goto L_220;
	}
	/* PROCEDURE..SCALE-CHECK */
L_100:
	;
L_110:
	;
	if (*sd1 > rgamsq)
		goto L_140;
	if (*sd1 == zero)
		goto L_170;
	igo = 130;
	/* FIX-H.. */
	goto L_70;
L_130:
	;
	*sd1 *= pow(gam, 2.0);
	*sx1 /= gam;
	sh11 /= gam;
	sh12 /= gam;
	goto L_110;
L_140:
	;
L_150:
	;
	if (*sd1 < gamsq)
		goto L_170;
	igo = 160;
	/* FIX-H.. */
	goto L_70;
L_160:
	;
	*sd1 = *sd1 / pow(gam, 2.0);
	*sx1 *= gam;
	sh11 *= gam;
	sh12 *= gam;
	goto L_150;
L_170:
	;
L_180:
	;
	if (fabs(*sd2) > rgamsq)
		goto L_200;
	if (*sd2 == zero)
		goto L_230;
	igo = 190;
	/* FIX-H.. */
	goto L_70;
L_190:
	;
	*sd2 *= pow(gam, 2.0);
	sh21 /= gam;
	sh22 /= gam;
	goto L_180;
L_200:
	;
L_210:
	;
	if (fabs(*sd2) < gamsq)
		goto L_230;
	igo = 220;
	/* FIX-H.. */
	goto L_70;
L_220:
	;
	*sd2 = *sd2 / pow(gam, 2.0);
	sh21 *= gam;
	sh22 *= gam;
	goto L_210;
L_230:
	;
	if (sflag < zero) {
		goto L_260;
	} else if (sflag == zero) {
		goto L_240;
	} else {
		goto L_250;
	}
L_240:
	;
	sparam[2] = sh21;
	sparam[3] = sh12;
	goto L_270;
L_250:
	;
	sparam[1] = sh11;
	sparam[4] = sh22;
	goto L_270;
L_260:
	;
	sparam[1] = sh11;
	sparam[2] = sh21;
	sparam[3] = sh12;
	sparam[4] = sh22;
L_270:
	;
	sparam[0] = sflag;
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  SROTM (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Apply a modified Givens plane rotation in single
                precision.

    Usage:      CALL SROTM (N, SX, INCX, SY, INCY, SPARAM)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).
                (Input/Output)
                SROTM replaces X(I) with H11*X(I)+H12*Y(I) for I=1,...,N.
                X(I) and Y(I) refer to specific elements of SX and SY.
                The H components refer to the rotation defined by SPARAM.
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be
                   SX(1+(I-1)*INCX) if INCX.GE.0  or
                   SX(1+(I-N)*INCX) if INCX.LT.0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).
                (Input/Output)
                SROTM replaces Y(I) with H21*X(I)+H22*Y(I) for I=1,...,N.
                X(I) and Y(I) refer to specific elements of SX and SY.
                The H components refer to the rotation defined by SPARAM.
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be
                   SY(1+(I-1)*INCY) if INCY.GE.0  or
                   SY(1+(I-N)*INCY) if INCY.LT.0.
       SPARAM - Real vector of length 5 which defines the rotation matrix
                H.  (Input)
                See remark.

    Remark:
       SROTM applies the modified Givens rotation H to the 2 by N matrix
                           ( X(1) ... X(N) )
                           ( Y(1) ... Y(N) )
       H takes one of the following forms,
       SPARAM(1) = -2.0
         H11 = 1.0        H12 = 0.0
         H21 = 0.0        H22 = 1.0
       SPARAM(1) = -1.0
         H11 = SPARAM(2)  H12 = SPARAM(4)
         H21 = SPARAM(3)  H22 = SPARAM(5)
       SPARAM(1) = 0.0
         H11 = 1.0        H12 = SPARAM(4)
         H21 = SPARAM(3)  H22 = 1.0
       SPARAM(1) = 1.0
         H11 = SPARAM(2)  H12 = 1.0
         H21 = -1.0       H22 = SPARAM(5)

    Keywords:   Level 1 BLAS; SROTM

    GAMS:       D1a8

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_srotm(Mint n, Mfloat sx[], Mint incx, Mfloat sy[], Mint incy, Mfloat sparam[])
#else
void imsl_srotm(n, sx, incx, sy, incy, sparam)
    Mint           n;
    Mfloat         sx[];
    Mint           incx;
    Mfloat         sy[];
    Mint           incy;
    Mfloat         sparam[];
#endif
{
	Mint           i, kx, ky;
	Mfloat         sflag, sh11, sh12, sh21, sh22, w, z;


	sflag = sparam[0];
	if (n > 0 && (sflag != -2.0)) {
		/* CODE FOR BOTH INCREMENTS EQUAL TO 1 */
		if (incx == 1 && incy == 1) {
			if (sflag == 0.0) {
				sh12 = sparam[3];
				sh21 = sparam[2];
				for (i = 1; i <= n; i++) {
					w = sx[i - 1];
					z = sy[i - 1];
					sx[i - 1] = w + z * sh12;
					sy[i - 1] = w * sh21 + z;
				}
			} else if (sflag > 0.0) {
				sh11 = sparam[1];
				sh22 = sparam[4];
				for (i = 1; i <= n; i++) {
					w = sx[i - 1];
					z = sy[i - 1];
					sx[i - 1] = w * sh11 + z;
					sy[i - 1] = -w + sh22 * z;
				}
			} else if (sflag < 0.0) {
				sh11 = sparam[1];
				sh12 = sparam[3];
				sh21 = sparam[2];
				sh22 = sparam[4];
				for (i = 1; i <= n; i++) {
					w = sx[i - 1];
					z = sy[i - 1];
					sx[i - 1] = w * sh11 + z * sh12;
					sy[i - 1] = w * sh21 + z * sh22;
				}
			}
			/*
			 * CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
			 * NOT EQUAL TO 1
			 */
		} else {
			kx = 1;
			ky = 1;
			if (incx < 0)
				kx = 1 + (1 - n) * incx;
			if (incy < 0)
				ky = 1 + (1 - n) * incy;
			if (sflag == 0.0) {
				sh12 = sparam[3];
				sh21 = sparam[2];
				for (i = 1; i <= n; i++) {
					w = sx[kx - 1];
					z = sy[ky - 1];
					sx[kx - 1] = w + z * sh12;
					sy[ky - 1] = w * sh21 + z;
					kx += incx;
					ky += incy;
				}
			} else if (sflag > 0.0) {
				sh11 = sparam[1];
				sh22 = sparam[4];
				for (i = 1; i <= n; i++) {
					w = sx[kx - 1];
					z = sy[ky - 1];
					sx[kx - 1] = w * sh11 + z;
					sy[ky - 1] = -w + sh22 * z;
					kx += incx;
					ky += incy;
				}
			} else if (sflag < 0.0) {
				sh11 = sparam[1];
				sh12 = sparam[3];
				sh21 = sparam[2];
				sh22 = sparam[4];
				for (i = 1; i <= n; i++) {
					w = sx[kx - 1];
					z = sy[ky - 1];
					sx[kx - 1] = w * sh11 + z * sh12;
					sy[ky - 1] = w * sh21 + z * sh22;
					kx += incx;
					ky += incy;
				}
			}
		}
	}
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  GIRTS/DGIRTS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 10, 1986

    Purpose:    Solve a triangular (possibly singular) set of linear
                systems and/or compute a generalized inverse of an upper
                triangular matrix.

    Usage:      CALL GIRTS (N, R, LDR, NB, B, LDB, IPATH, IRANK, X, LDX,
                            RINV, LDRINV)

    Arguments:
       N      - Order of the upper triangular matrix R.  (Input)
       R      - N by N upper triangular matrix.  (Input)
                If R contains a zero along the diagonal, the remaining
                elements of the row must also be zero.  Only
                the upper triangle of R is referenced.
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement of the calling program.  (Input)
       NB     - Number of columns in B.  (Input)
                NB must be nonnegative.  If NB is zero, no linear systems
                are solved.
       B      - N by NB matrix containing the right hand sides of the
                linear system.  (Input, if NB .GT. 0)
                If NB = 0, B is not referenced and can be a vector of
                length one.
       LDB    - Leading dimension of B exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPATH  - Path option.  (Input)
                IPATH  Action
                  1    Solve R*X = B.
                  2    Solve trans(R)*X = B.
                  3    Solve R*X = B and compute RINV.
                  4    Solve trans(R)*X = B and compute RINV.
       IRANK  - Rank of R.  (Output)
       X      - N by NB matrix containing the solution matrix
                corresponding to the right hand side B.  (Output, if
                NB .GT.0)
                If B is not needed, then X and B can share the same
                storage locations.  If NB = 0, X is not referenced
                and can be a vector of length one.
       LDX    - Leading dimension of X exactly as specified in the
                dimension statement of the calling program.  (Input)
       RINV   - N by N upper triangular matrix that is the inverse of R
                when R is nonsingular.  (Output, if IPATH equals 3 or 4)
                (When R is singular, RINV is a g3 inverse.  See Remark 3
                for an explanation of g3 inverses.)   If IPATH = 1 or 2,
                RINV is not referenced and can be a vector of length one.
                If IPATH = 3 or 4 and R is not needed, then R and RINV
                can share the same storage locations.
       LDRINV - Leading dimension of RINV exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remarks:
    1. Informational error
       Type Code
         3   1  The linear system of equations is inconsistent.

    2. GIRTS assumes that a singular R is represented by zero rows
       in R.  No other forms of singularity in R are allowed.

    3. RINV is a g3 inverse means that it satisfies conditions 1, 2, and
       3 for the Moore-Penrose inverse but generally fails condition 4.
       The four conditions for AINV to be a Moore-Penrose inverse of A
       are as follows:
          1.  A*AINV*A = A
          2.  AINV*A*AINV = AINV
          3.  A*AINV is symmetric
          4.  AINV*A is symmetric.

    GAMS:       D9

    Chapter:    STAT/LIBRARY Mathematical Support

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_girts(Mint n, Mfloat r[], Mint ldr, Mint nb, Mfloat b[], Mint ldb,
        Mint ipath, Mint *irank, Mfloat x[], Mint ldx, Mfloat rinv[], Mint ldrinv)
#else
static void l_girts(n, r, ldr, nb, b, ldb, ipath, irank, x, ldx, rinv, ldrinv)
    Mint           n;
    Mfloat         r[];
    Mint           ldr, nb;
    Mfloat         b[];
    Mint           ldb, ipath, *irank;
    Mfloat         x[];
    Mint           ldx;
    Mfloat         rinv[];
    Mint           ldrinv;
#endif
{
	Mint           i, j, k, ner;
	Mfloat         temp1, temp2;


	imsl_e1psh("l_girts");
	/* Check for terminal errors */
	ner = 1;
	imsl_c1dim(1, n, "n", ldr, "ldr", &ner);
	imsl_c1iarg(nb, "nb", 0, -1, &ner);
	imsl_c1iarg(ipath, "ipath", 1, 4, &ner);
	if (nb > 0) {
		imsl_c1dim(1, nb, "*nb", ldb, "tdb", &ner);
		imsl_c1dim(1, nb, "*nb", ldx, "tdx", &ner);
	} else {
		ner += 6;
	}
	if (ipath == 3 || ipath == 4) {

		imsl_c1dim(1, n, "*n", ldrinv, "tdrinv", &ner);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * Linear dependent rows of R must be represented only by rows whose
	 * elements are zero.
	 */
	imsl_c1r(n, r, ldr, &ner);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Get rank */
	*irank = 0;
	for (i = 1; i <= n; i++) {
		if (r[i + ldr * (i - 1) - 1] != 0.0)
			*irank += 1;
	}
	/*
	 * Make A copy of B in X and work with X
	 */
	for (j = 1; j <= nb; j++) {
		scopy(n, &b[ldb * (j - 1)], 1, &x[ldx * (j - 1)], 1);
	}

	if (ipath == 1 || ipath == 3) {
		/* Solve R*X = B */
		if (*irank < n) {
			for (i = 1; i <= nb; i++) {
				for (j = n; j >= 1; j--) {
					if (r[j + ldr * (j - 1) - 1] == 0.0) {
						if (x[j + ldx * (i - 1) - 1] != 0.0) {
							imsl_e1sti(1, j);
							imsl_e1sti(2, i);
							imsl_e1str(1, x[j + ldx * (i - 1) - 1]);

/*							imsl_ermes(3, 1, "The linear system of equations is inconsistent within the computed tolerance.  Elements of row %(i1) are zero, but B(%(i1),%(i2)) = %(r1).  X(%(i1),%(i2)) is set to zero.");
*/

                                                        imsl_ermes(IMSL_WARNING, IMSL_TOLERALNCE_INCONSISTENT);
						}
						x[j + ldx * (i - 1) - 1] = 0.0;
					} else {
						x[j + ldx * (i - 1) - 1] /= r[j + ldr * (j - 1) - 1];
						temp1 = -x[j + ldx * (i - 1) - 1];
						saxpy(j - 1, temp1, &r[ldr * (j - 1)], 1, &x[ldx * (i - 1)],
						      1);
					}
				}
			}
		} else {
			for (i = 1; i <= nb; i++) {

				imsl_strsv("UPPER", "NOT-TRANS", "NOT-DIAG", n, r, ldr, &x[ldx * (i - 1)], 1);
			}
		}

	} else if (ipath == 2 || ipath == 4) {
		/* Solve R'*X=B */
		if (*irank < n) {
			/* Case of singular R */
			for (i = 1; i <= nb; i++) {
				for (j = 1; j <= n; j++) {
					temp1 = x[j + ldx * (i - 1) - 1] - imsl_sdot(j - 1,
										  &r[ldr * (j - 1)], 1, &x[ldx * (i - 1)], 1);
					if (r[j + ldr * (j - 1) - 1] == 0.0) {
						temp2 = fabs(x[j + ldx * (i - 1) - 1]) + imsl_a1ot(j -
											      1, &r[ldr * (j - 1)], 1, &x[ldx * (i - 1)], 1);
						temp2 *= 200.0 * imsl_amach(4);
						if (fabs(temp1) > temp2) {
							imsl_e1sti(1, j);
							imsl_e1sti(2, i);

/*							imsl_ermes(3, 1, "The linear system of equations is inconsistent within the computed tolerance.  X(%(i1),%(i2)) is set to zero.");
*/

imsl_ermes(IMSL_WARNING, IMSL_TOLERALNCE_INCONSISTENT_2);
						}
						x[j + ldx * (i - 1) - 1] = 0.0;
					} else {
						x[j + ldx * (i - 1) - 1] = temp1 / r[j + ldr * (j - 1) - 1];
					}
				}
			}
		} else {
			/* Case of nonsingular R */
			for (i = 1; i <= nb; i++) {

				imsl_strsv("UPPER", "TRANSPOSE", "NOT-UNIT", n, r, ldr, &x[ldx * (i - 1)], 1);
			}
		}
	}
	if (ipath == 3 || ipath == 4) {
		/* Invert R */
		for (j = 1; j <= n; j++) {
			scopy(j, &r[ldr * (j - 1)], 1, &rinv[ldrinv * (j - 1)], 1);
		}
		for (k = 1; k <= n; k++) {
			if (rinv[k + ldrinv * (k - 1) - 1] == 0.0) {
				sset(k, 0.0, &rinv[ldrinv * (k - 1)], 1);
				sset(n - k, 0.0, &rinv[k + ldrinv * k - 1], ldrinv);
			} else {
				rinv[k + ldrinv * (k - 1) - 1] = 1.0 / rinv[k + ldrinv * (k - 1) - 1];
				temp1 = -rinv[k + ldrinv * (k - 1) - 1];
				sscal(k - 1, temp1, &rinv[ldrinv * (k - 1)], 1);
				if (k < n) {
					imsl_sger(k - 1, n - k, 1.0, &rinv[ldrinv * (k - 1)],
						1, &rinv[k + ldrinv * k - 1], ldrinv, &rinv[ldrinv * k],
						ldrinv);
					sscal(n - k, rinv[k + ldrinv * (k - 1) - 1], &rinv[k + ldrinv * k - 1],
					      ldrinv);
				}
			}
		}
		/*
		 * Fill lower triangle of RINV with zeros
		 */
		for (i = 1; i <= (n - 1); i++) {
		    sset(n - i, 0.0, &rinv[(i + 1) + ldr * (i - 1) - 1], 1);
		}
	}
	/* Exit section */
L_9000:
	imsl_e1pop("l_girts");
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  C1R/DC1R (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 17, 1985

    Purpose:    Call E1MES if a diagonal element of the Cholesky factor R
                is zero but a remaining element in that row of R is not.

    Usage:      CALL C1R (N, R, LDR, NER)

    Arguments:
       N      - Order of R.  (Input)
       R      - Upper triagular matrix, N by N, with the Cholesky factor.
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement of the calling program.  (Input)
       NER    - Error code used in call to E1MES.  (Input/Output)
                On output NER is incremented by one.

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c1r(Mint n, Mfloat r[], Mint ldr, Mint *ner)
#else
void imsl_c1r(n, r, ldr, ner)
    Mint            n;
    Mfloat          r[];
    Mint            ldr, *ner;
#endif
{
	Mint            i, j;


	for (i = 1; i <= n; i++) {
		if (r[i + ldr * (i - 1) - 1] == 0.0) {
			for (j = i + 1; j <= n; j++) {
				if (r[i + ldr * (j - 1) - 1] != 0.0) {
					imsl_e1sti(1, i);
					imsl_e1sti(2, j);
					imsl_e1str(1, r[i + ldr * (i - 1) - 1]);
					imsl_e1str(2, r[i + ldr * (j - 1) - 1]);

/*					imsl_ermes(5, *ner, "R(%(i1),%(i1)) = %(r1).  Remaining elements for the row must also be zero, but R(%(i1),%(i2)) = %(r2) is not.");
*/
                                        imsl_ermes(IMSL_TERMINAL,
					IMSL_REMAINING_ELMNTS_NOT_ZERO);
					goto L_9000;
				}
			}
		}
	}
	*ner += 1;
L_9000:
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  RCOVB/DRCOVB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 1, 1985

    Purpose:    Compute the estimated variance-covariance matrix of the
                estimated regression coefficients given the R matrix.

    Usage:      CALL RCOVB (NCOEF, R, LDR, S2, COVB, LDCOVB)

    Arguments:
       NCOEF  - Number of regression coefficients in the model.
                (Input)
       R      - NCOEF by NCOEF upper triangular matrix containing the R
                matrix.  (Input)
                The R matrix can come from a regression fit based on a QR
                decomposition of the matrix of regressors or based on a
                Cholesky factorization trans(R)*R of the matrix of sums
                of squares and crossproducts of the regressors.  Elements
                to the right of a diagonal element of R that is zero must
                also be zero.  A zero row indicates a nonfull rank model.
                For an R matrix that comes from a regression fit with
                linear equality restrictions on the parameters, each row
                of R corresponding to a restriction must have a
                corresponding diagonal element that is negative. The
                remaining rows of R must have positive diagonal elements.
                Only the upper triangle of R is referenced.
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement in the calling program.  (Input)
       S2     - s-squared, the estimated variance of the error in the
                regression model.  (Input)
                S2 is the error mean square from the regression fit.
       COVB   - NCOEF by NCOEF matrix that is the estimated
                variance-covariance matrix of the estimated regression
                coefficients when R is nonsingular and is from an
                unrestricted regression fit.  (Output)
                See Remark for an explanation of COVB when R is singular
                or R is from a restricted regression fit.  If R is not
                needed, COVB and R can share the same storage locations.
       LDCOVB - Leading dimension of COVB exactly as specified in the
                dimension statement in the calling program.  (Input)

    Remark:
       When R is nonsingular and comes from an unrestricted regression
       fit, COVB is the estimated variance-covariance matrix of the
       estimated regression coefficients, and COVB =  S2*inv(trans(R)*R).
       Otherwise, variances and covariances of estimable functions of the
       regression coefficients can be obtained using COVB, and
       COVB = S2*G*D*trans(G).  Here, D is the diagonal matrix with
       diagonal elements equal to 0 if the corresponding rows of R
       are restrictions, and with diagonal elements equal to one
       otherwise.  Also, G is a particular generalized inverse of R.
       See the Algorithm section.

    Keywords:   Statistical inference; General linear model analysis

    GAMS:       L8a10; L8h

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_rcovb(Mint ncoef, Mfloat *r, Mint ldr, Mfloat s2, Mfloat *covb, Mint ldcovb)
#else
void imsl_rcovb(ncoef, r, ldr, s2, covb, ldcovb)
	Mint            ncoef;
	Mfloat         *r;
	Mint            ldr;
	Mfloat          s2, *covb;
	Mint            ldcovb;
#endif
{
#define R(I_,J_)	(r+(I_)*(ldr)+(J_))
#define COVB(I_,J_)	(covb+(I_)*(ldcovb)+(J_))
	Mint            irank, j, k, ner;
	Mfloat          t, wk[1];


	imsl_e1psh("imsl_rcovb");
	/* Check for terminal errors */
	ner = 1;

	imsl_c1dim(1, ncoef, "m", ldr, "ldr", &ner);

	imsl_c1dim(1, ncoef, "*m", ldcovb, "cov_col_dim", &ner);
	imsl_c1ge0(s2, "s2", &ner);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	imsl_c1r(ncoef, r, ldr, &ner);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * Compute a g3 inverse of R
	 */
	l_girts(ncoef, r, ldr, 0, wk, 1, 3, &irank, wk, 1, covb, ldcovb);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Form G3INV(R)*G3INV(R)' */
	for (j = 1; j <= ncoef; j++) {
		if (*COVB(j - 1, j - 1) > 0.0) {
			for (k = 1; k <= (j - 1); k++) {
				t = *COVB(j - 1, k - 1);
				saxpy(k, t, COVB(j - 1, 0), 1, COVB(k - 1, 0), 1);
			}
			t = *COVB(j - 1, j - 1);
			sscal(j, t, COVB(j - 1, 0), 1);
		} else {
			sset(j, 0.0, COVB(j - 1, 0), 1);
		}
	}
	/* Multiply by S2 */
	for (j = 1; j <= ncoef; j++) {
		sscal(j, s2, COVB(j - 1, 0), 1);
	}
	/* Fill lower portion of COVB */
	for (j = 1; j < ncoef; j++){
	    scopy (ncoef-j, COVB(j, j - 1), ldcovb, COVB(j - 1, j), 1);
	}
	/* Exit section */
L_9000:
	imsl_e1pop("imsl_rcovb");
	return;

#undef R
#undef COVB
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  A1OT/DA1OT (Single/Double precision version)

    Computer:   FORC/SINGLE

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
Mfloat imsl_a1ot(Mint n, Mfloat sx[], Mint incx, Mfloat sy[], Mint incy)
#else
Mfloat imsl_a1ot(n, sx, incx, sy, incy)
    Mint            n;
    Mfloat          sx[];
    Mint            incx;
    Mfloat          sy[];
    Mint            incy;
#endif
{
	Mint            i, ix, iy;
	Mfloat          a1ot_v;


	a1ot_v = 0.0;
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
/*
  -----------------------------------------------------------------------
    IMSL Name:  C1GE0/DC1GE0 (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 25, 1985

    Purpose:    Call E1MES if a real argument is negative.

    Usage:      CALL C1GE0 (RARG, NMARG, NER)

    Arguments:
       RARG   - Argument to be checked.  (Input)
       NMARG  - Character string that contains the name of the
                argument to be checked.
       NER    - Error code.  (Input/Output)
                The input NER is the error code used by E1MES.  NER is
                incremented by 1 on output.

    Copyright:  1985 by IMSL, Inc.  All rights reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c1ge0(Mfloat rarg, Mchar *nmarg, Mint *ner)
#else
void imsl_c1ge0(rarg, nmarg, ner)
	Mfloat          rarg;
	Mchar           *nmarg;
	Mint            *ner;
#endif
{

	if (rarg < 0.0) {
		imsl_e1str(1, rarg);
		imsl_e1stl(1, nmarg);

/*		imsl_ermes(5, *ner, "%(l1) = %(r1).  %(l1) must be greater than or equal to zero.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_RARG_GE_ZERO);
	}
	*ner += 1;
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
void imsl_c1div(Mfloat top, Mfloat bottom, Mfloat *quot)
#else
void imsl_c1div(top, bottom, quot)
	Mfloat      top, bottom, *quot;
#endif
{
	Mfloat	    absbot, anan, big, small;


	anan = imsl_amach(6);
	if (imsl_ifnan(top) || imsl_ifnan(bottom)) {
		*quot = anan;
		return;
	}
	absbot = fabs(bottom);
	if (absbot <= 1.0) {
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
		if (top == 0.0) {
			*quot = anan;
		} else if (bottom >= 0.0) {
			if (top >= 0.0) {
				*quot = imsl_amach(7);
			} else {
				*quot = imsl_amach(8);
			}
		} else {
			if (top >= 0.0) {
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
		*quot = 0.0;
		return;
	}
L_10:
	*quot = top / bottom;
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  G1AOV/DG1AOV (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 20, 1985

    Purpose:    Get an analysis of variance table and related statistics.

    Usage:      CALL G1AOV (DFR, SSR, DFE, SSE, GMEAN, AOV)

    Arguments:
       DFR    - Degrees of freedom for regression.  (Input)
       SSR    - Sum of squares for regression.  (Input)
       DFE    - Degrees of freedom for error.  (Input)
       SSE    - Sum of squares for error.  (Input)
       GMEAN  - Grand mean.  (Input)
                If the grand mean is not known it may be set to
                not-a-number.
       AOV    - Vector of length 15 that contains statistics relating
                to the analysis of variance.  (Output)
                  I            AOV(I)
                  1      Degrees of freedom for regression
                  2      Degrees of freedom for error
                  3      Total degrees of freedom
                  4      Sum of squares for regression
                  5      Sum of squares for error
                  6      Total sum of squares
                  7      Regression mean square
                  8      Error mean square
                  9      F statistic
                 10      p-value
                 11      R-squared (in percent)
                 12      Adjusted R-squared (in percent)
                 13      Estimated standard deviation of the model error
                 14      Mean of the response (dependent variable)
                 15      Coefficient of variation (in percent)

    Remark:
       For writing the analysis of variance table see WRAOV/DWRAOV.

    Keyword:    F test

    GAMS:       L8a10

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_g1aov(Mfloat dfr, Mfloat ssr, Mfloat dfe, Mfloat sse,
                    Mfloat gmean, Mfloat aov[])
#else
void imsl_g1aov(dfr, ssr, dfe, sse, gmean, aov)
	Mfloat          dfr, ssr, dfe, sse, gmean, aov[];
#endif
{
	Mint            ner;
	Mfloat          anan, cv, r2, r2adj, rms, tms;


	imsl_e1psh("imsl_g1aov");
	/* Check for terminal errors */
	ner = 1;
	imsl_c1ge0(dfr, "DFR", &ner);
	imsl_c1ge0(ssr, "SSR", &ner);
	imsl_c1ge0(dfe, "DFE", &ner);
	imsl_c1ge0(sse, "SSE", &ner);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	anan = imsl_amach(6);

	/* Get degrees of freedom */
	aov[0] = dfr;
	aov[1] = dfe;
	aov[2] = dfr + dfe;

	/* Get sums of squares */
	aov[3] = ssr;
	aov[4] = sse;
	aov[5] = ssr + sse;

	/* Get mean squares */
	imsl_c1div(ssr, dfr, &aov[6]);
	imsl_c1div(sse, dfe, &aov[7]);

	/* Get F value and p-value */
	imsl_c1div(aov[6], aov[7], &aov[8]);
	if (imsl_ifnan(aov[8])) {
	    aov[9] = imsl_amach(6);
	}
	else if (aov[8] == imsl_amach(7)) {
	    aov[9] = 1.0;
	}
	else{
	    aov[9] = 1.0 - imsl_f_F_cdf(aov[8],dfr,dfe);
	}

	/* Compute R-squared */
	imsl_c1div(100.0 * ssr, aov[5], &r2);
	aov[10] = r2;
	/*
	 * Compute adjusted R-squared = MAX(0.0,1.0-(ERROR MS)/(TOTAL MS))
	 */
	imsl_c1div(aov[5], aov[2], &tms);
	imsl_c1div(aov[7], tms, &r2adj);
	if (!imsl_ifnan(r2adj)) {
		if (r2adj > 1.0) {
			r2adj = 0.0;
		} else {
			r2adj = 100.0 * (1.0 - r2adj);
		}
	}
	aov[11] = r2adj;
	/*
	 * Compute std. dev. of error and coef. of variation.
	 */
	if (imsl_ifnan(aov[7])) {
		aov[12] = anan;
		aov[14] = anan;
	} else {
		rms = sqrt(aov[7]);
		aov[12] = rms;
		imsl_c1div(100.0 * rms, gmean, &cv);
		aov[14] = cv;
	}
	/* Store grand mean. */
	aov[13] = gmean;
	/* Exit section */
L_9000:
	imsl_e1pop("imsl_g1aov");
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  SXYZ (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Compute a single precision XYZ product.

    Usage:      SXYZ(N, SX, INCX, SY, INCY, SZ, INCZ)

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
       SZ     - Real vector of length MAX(N*IABS(INCZ),1).  (Input)
       INCZ   - Displacement between elements of SZ.  (Input)
                Z(I) is defined to be.. SZ(1+(I-1)*INCZ) if INCZ .GE. 0
                or SZ(1+(I-N)*INCZ) if INCZ .LT. 0.
       SXYZ   - Sum from I=1 to N of X(I)*Y(I)*Z(I).  (Output)
                X(I), Y(I) and Z(I) refer to specific elements of SX,
                SY and SZ, respectively.  See INCX, INCY and INCZ
                argument descriptions.

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
Mfloat imsl_sxyz(Mint n, Mfloat sx[], Mint incx, Mfloat sy[], Mint incy, Mfloat sz[], Mint incz)
#else
Mfloat imsl_sxyz(n, sx, incx, sy, incy, sz, incz)
    Mint             n;
    Mfloat           sx[];
    Mint             incx;
    Mfloat           sy[];
    Mint             incy;
    Mfloat           sz[];
    Mint             incz;
#endif
{
	Mint             i, ix, iy, iz;
	Mfloat           sxyz_v;


	sxyz_v = 0.0;
	if (n <= 0)
		goto L_9000;
	if ((incx != 1 || incy != 1) || incz != 1) {

		/*
		 * CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
		 * TO 1
		 */
		ix = 1;
		iy = 1;
		iz = 1;
		if (incx < 0)
			ix = (-n + 1) * incx + 1;
		if (incy < 0)
			iy = (-n + 1) * incy + 1;
		if (incz < 0)
			iz = (-n + 1) * incz + 1;
		for (i = 1; i <= n; i++) {
			sxyz_v += sx[ix - 1] * sy[iy - 1] * sz[iz - 1];
			ix += incx;
			iy += incy;
			iz += incz;
		}
	} else {
		for (i = 1; i <= n; i++) {
			sxyz_v += sx[i - 1] * sy[i - 1] * sz[i - 1];
		}
	}
L_9000:
	return (sxyz_v);
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

