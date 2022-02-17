#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static VA_LIST_HACK	PROTO(l_poly_regression,(Mint, Mfloat[], Mfloat[], Mint,
			va_list));
static void	PROTO(l_r2orp,(Mint, Mint, Mfloat*, Mint, Mint, Mint, Mint,
			Mint, Mint, Mint, Mfloat, Mint, Mint*, Mfloat*,
			Mfloat*, Mfloat[], Mfloat[], Mfloat[], Mfloat[],
			Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mint*, Mfloat[],
			Mint[]));
static void	PROTO(l_r3orp,(Mint, Mfloat[], Mfloat[], Mfloat[], Mfloat[],
			Mint, Mint, Mfloat, Mint, Mint*, Mfloat*, Mfloat*, 
			Mfloat[], Mfloat[], Mfloat[], Mfloat[], Mfloat*, 
			Mfloat*, Mfloat*, Mfloat*, Mint[], Mfloat*));
static void	PROTO(l_r4orp,(Mint, Mfloat*, Mint, Mint, Mint, Mint, Mint, 
			Mint, Mint, Mfloat, Mint, Mint*, Mfloat*, Mfloat*,
			Mfloat[], Mfloat[], Mfloat[], Mfloat[], Mfloat*,
			Mfloat*, Mfloat*, Mfloat*, Mint*, Mfloat[], Mfloat[], 
			Mfloat[], Mfloat[], Mint[], Mfloat[]));
static void	PROTO(l_r2tap,(Mint, Mfloat[], Mfloat[], Mfloat, Mfloat, 
			Mfloat[], Mfloat[], Mfloat, Mfloat, Mint, Mfloat, 
			Mfloat, Mint, Mfloat[], Mfloat*, Mint, Mfloat*, Mint, 
			Mfloat*, Mint, Mfloat[]));
static void	PROTO(l_r3tap,(Mint, Mfloat[], Mfloat[], Mfloat[], Mfloat, 
			Mfloat, Mfloat*, Mint));
static void     PROTO(l_r4tap,(Mint, Mfloat[], Mfloat[], Mfloat[], Mfloat, 
			Mfloat, Mfloat*, Mint));
static void	PROTO(l_r5tap,(Mint, Mfloat[], Mfloat[], Mfloat, Mfloat, 
			Mfloat[], Mfloat[], Mfloat, Mfloat, Mint, Mfloat, 
			Mfloat, Mint, Mfloat[], Mfloat*, Mint, Mfloat*, Mint, 
			Mfloat*, Mint, Mfloat[], Mfloat[]));
static void	PROTO(l_c1f,( Mfloat, Mfloat, Mfloat, Mfloat, Mfloat*, 
			Mfloat*));
static void	PROTO(l_c1t,(Mfloat, Mfloat, Mfloat, Mfloat*, Mfloat*));
static Mint	PROTO(l_ismax,(Mint, Mfloat[], Mint));
static Mint	PROTO(l_ismin,(Mint, Mfloat[], Mint));
static Mfloat	PROTO(l_ssum,(Mint, Mfloat[], Mint));
static void	PROTO(l_rcoef,(Mint, Mfloat[], Mfloat[], Mfloat, Mfloat*, Mint));


static Mfloat   *lv_coefficients;
 
#ifdef ANSI
Mfloat *imsl_f_poly_regression(Mint nobs, Mfloat x_data[], Mfloat y_data[], Mint degree, ...)
#else
Mfloat *imsl_f_poly_regression(nobs, x_data, y_data, degree, va_alist)
    Mint        nobs;
    Mfloat      x_data[], y_data[];
    Mint	degree;
    va_dcl
#endif
{
    va_list     argptr;
    VA_START(argptr,degree);

    E1PSH("imsl_f_poly_regression","imsl_d_poly_regression");

    lv_coefficients = NULL;
    IMSL_CALL(l_poly_regression(nobs, x_data, y_data, degree, argptr));
    va_end(argptr);

    E1POP("imsl_f_poly_regression","imsl_d_poly_regression");

    return lv_coefficients;
}

#ifdef ANSI
static VA_LIST_HACK l_poly_regression(Mint nobs, Mfloat x_data[], Mfloat y_data[], Mint degree, va_list argptr)
#else
static VA_LIST_HACK l_poly_regression(nobs, x_data, y_data, degree, argptr)
    Mint        nobs;
    Mfloat      x_data[], y_data[];
    Mint	degree;
    va_list     argptr;
#endif
{

    Mfloat	*weights	= NULL;
    Mfloat	*sspoly		= NULL;
    Mfloat	**sspoly_ptr	= NULL;
    Mint	sspoly_col_dim  = 4;
    Mfloat	*sslof		= NULL;
    Mfloat	**sslof_ptr     = NULL;
    Mint	sslof_col_dim	= 4;
    Mfloat	*anova		= NULL;
    Mfloat	**anova_ptr	= NULL;
    Mfloat	*residual	= NULL;
    Mfloat	**residual_ptr  = NULL;
    Mfloat	*x_mean		= NULL;
    Mfloat	*x_variance	= NULL;
    Mint	*df_pure	= NULL;
    Mfloat	*ssq_pure	= NULL;

    Mint	code		= 1;
    Mint	arg_number	= 4;
    Mint        user_sspoly	= 0;
    Mint	return_sspoly	= 0;
    Mint        user_sslof	= 0;	
    Mint        return_sslof	= 0;
    Mint        user_anova	= 0;
    Mint        return_anova	= 0;
    Mint        user_residual	= 0;
    Mint        return_residual = 0;
    Mint        user_coef       = 0;
    Mint        lof		= 0;

    Mint	*iwk    = NULL;
    Mfloat	*x	= NULL;
    Mfloat	*a	= NULL;
    Mfloat	*b	= NULL; 
    Mfloat	*wk	= NULL;
    Mfloat	*scoef	= NULL;
    Mfloat	*d	= NULL;
    Mint	error   = 0;
    Mint        air     = 0;
    Mint	ncol, iwt, itemp;

    while (code > 0) {    
	code = va_arg(argptr, int);
        ++arg_number;
	switch (code) {
	    case IMSL_RETURN_USER:
		lv_coefficients = va_arg(argptr, Mfloat*);
		user_coef  = 1;
		++arg_number;
		break;
	    case IMSL_WEIGHTS:
		weights = va_arg(argptr, Mfloat*);
		++arg_number;
		break;
	    case IMSL_SSQ_POLY:
		sspoly_ptr = va_arg(argptr, Mfloat**);
		user_sspoly   = 0;
		return_sspoly = 1;
		++arg_number;
		break;
	    case IMSL_SSQ_POLY_USER:
		sspoly = va_arg(argptr, Mfloat*);
		user_sspoly   = 1;
		return_sspoly = 1;
		++arg_number;
		break;
	    case IMSL_SSQ_POLY_COL_DIM:
		sspoly_col_dim = va_arg(argptr, Mint);
		++arg_number;
		break;
	    case IMSL_SSQ_LOF:
		sslof_ptr = va_arg(argptr, Mfloat**);
		user_sslof   = 0;
		return_sslof = 1;
		lof          = 1;
		++arg_number;
		break;
	    case IMSL_SSQ_LOF_USER:
		sslof = va_arg(argptr, Mfloat*);
		user_sslof   = 1;
		return_sslof = 1;
		lof          = 1;
		++arg_number;
		break;
	    case IMSL_SSQ_LOF_COL_DIM:
		sslof_col_dim = va_arg(argptr, Mint);
		++arg_number;
		break;
	    case IMSL_ANOVA_TABLE:
		anova_ptr = va_arg(argptr, Mfloat**);
		return_anova = 1;
		user_anova   = 0;
		++arg_number;
		break;
	    case IMSL_ANOVA_TABLE_USER:
		anova = va_arg(argptr, Mfloat*);
		return_anova = 1;
		user_anova   = 1;
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
	    case IMSL_X_MEAN:
		x_mean = va_arg(argptr, Mfloat*);
		++arg_number;
		break;
	    case IMSL_X_VARIANCE:
		x_variance = va_arg(argptr, Mfloat*);
	        ++arg_number;
		break;
	    case IMSL_SSQ_PURE_ERROR:
		ssq_pure = va_arg(argptr, Mfloat*);
		++arg_number;
		break;
	    case IMSL_DF_PURE_ERROR:
		df_pure = va_arg(argptr, Mint*);
		++arg_number;
		break;
	    case 0:
		break;
	    default:
                imsl_e1sti (1, code);
                imsl_e1sti (2, arg_number);
                imsl_ermes(IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
		return argptr;
	}
    }
    if (sslof_col_dim < 4) {
        imsl_e1sti(1, sslof_col_dim);
        imsl_ermes(IMSL_TERMINAL, IMSL_LOF_COL_DIM_4); 
	air = 1;
    }
    if (sspoly_col_dim < 4) {
        imsl_e1sti(1, sspoly_col_dim);
        imsl_ermes(IMSL_TERMINAL, IMSL_POLY_COL_DIM_4); 
	air = 1;
    }

    if (x_data == NULL) {
        imsl_e1stl (1, "x");
        imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
        air = 1;
    }

    if (y_data == NULL) {
        imsl_e1stl (1, "y");
        imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
        air = 1;
    }

    if (nobs < 1) {
	imsl_e1sti(1, nobs);
        imsl_ermes(IMSL_TERMINAL, IMSL_NOBS_LESS_THAN_ONE);
	air = 1;
    }
	/* Check NDEG */
    if (degree < 0) {
	imsl_e1sti(1, degree);
/*	imsl_ermes(5, 2, "DEGREE = %(i1).  The degree of the polynomial to be fit must be greater than or equal to 0.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_POLY_DEGREE_LT_ZERO);
	air = 1;
    }
    if (air==1) return argptr;

    if (!weights){
	ncol = 2;
	iwt  = 0;
    }
    else{
	ncol = 3;
	iwt  = 3;
    }
    x  = (Mfloat *) imsl_malloc(ncol*nobs*sizeof(Mfloat));
    /*The extra +1 is needed for the allocations of a and b */
    /*because it is illegal to allocate 0 workspace*/
    a     = (Mfloat *) imsl_malloc ((degree+1)*sizeof(Mfloat));
    b     = (Mfloat *) imsl_malloc ((degree+1)*sizeof(Mfloat));
    scoef = (Mfloat *) imsl_malloc ((degree+1)*sizeof(Mfloat));
    d     = (Mfloat *) imsl_malloc ((degree+1)*sizeof(Mfloat));
    wk    = (Mfloat *) imsl_malloc (8*nobs*sizeof(Mfloat));
    iwk = (Mint *) imsl_malloc (nobs*sizeof(Mint));

    if (!x || !a || !b || !scoef || !d || !wk || !iwk) {
                imsl_e1sti (1,degree);
                imsl_e1stl (1,"degree");
                imsl_e1sti (2,nobs);
                imsl_e1stl (2,"n_observations");
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	error = 1;
    } else {
	if (!user_coef) {
	    lv_coefficients = (Mfloat *) imsl_malloc((degree+1)*sizeof(Mfloat));
	    if (!lv_coefficients){
                imsl_e1sti (1,degree);
                imsl_e1stl (1,"degree");
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		error = 1;
	    }
	}

	if (!error){
	    Mint    ldx	    = nobs;
	    Mint    ifrq    = 0;
	    Mint    icrit   = 0;
	    Mint    ind	    = 1;
	    Mint    irsp    = 2;
	    Mint    nrmiss;
	    Mint    deg_out;
	    Mfloat  crit=0.0, smultc, saddc, dfe, sse, dfpe=0.0, sspe=0.0,  dft;

	    scopy (nobs, x_data, 1, x, 1);
	    scopy (nobs, y_data, 1, x+nobs, 1);
	    if (weights) scopy (nobs, weights, 1, x+2*nobs, 1);

	    l_r2orp(nobs, ncol, x, ldx, irsp, ind, ifrq, iwt, degree, icrit, 
		    crit, lof, &deg_out, &smultc, &saddc, a, b, scoef, d, 
		    &dfe, &sse, &dfpe, &sspe, &nrmiss, wk, iwk);
            if ((imsl_n1rty(1) > 3) && (imsl_n1rty(1)!=6)) {
		error = 1;
	    } else {
		error = 0;
	    }

	    if (!error){
		Mint	ldsslof  = degree;
		Mint	ldsspoly = degree;
		Mint    ldcoef = deg_out + 1;
		Mint	iprint = 0;
		Mint    i, j, dplus1 = degree+1;
		Mfloat  *coef, yhat, *res;

		if (x) imsl_free(x);   x = NULL;
		if (wk) imsl_free(wk);  wk = NULL;
		if (iwk) imsl_free(iwk); iwk = NULL;

                if (return_residual) {
                  /***Computations for residuals***/
                  wk = (Mfloat *) imsl_malloc(dplus1*sizeof(Mfloat));
                  if (!user_residual) {
                      *residual_ptr = (Mfloat*) imsl_malloc 
                                      (nobs*sizeof(Mfloat));
                      res = *residual_ptr;
		  } else {
		      res = residual;
		  }
                  if (!wk || !res) {
                      imsl_e1sti (1,degree);
                      imsl_e1stl (1,"degree");
                      imsl_e1sti (2,nobs);
                      imsl_e1stl (2,"n_observations");
		      imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
                    if (wk) imsl_free(wk);
                    if (res) imsl_free(res);
                    goto RETURN;
		  }
	          for (i=0;i<nobs;i++) {
                        if (!imsl_ifnan(x_data[i]) && !imsl_ifnan(y_data[i])) {
                              wk[0] = 1;
                              for (j=2;j<=dplus1;j++){
				 if (j==2) {
                                     wk[j-1] = (smultc*x_data[i] +saddc) - a[j-2];
                                 } else {
				     wk[j-1] = ((smultc*x_data[i] +saddc) - a[j-2])*wk[j-2]-b[j-2]*wk[j-3]; 
				 }
			      }
                              yhat = imsl_sdot(dplus1, scoef, 1, wk, 1);
                              *(res+i) = y_data[i] - yhat; 
			} else {
			    *(res+i) = imsl_amach(6);
			}
		  }
		  if (wk) imsl_free(wk); wk = NULL;
                 /********************************/
		}


                if (degree !=0) {
		if (!user_sspoly){
/*The following line is removed because we are now allocating degree amount of*/
/*space no matter what the deg_out degree of the polynomial is.               */
         	    /*if (!return_sspoly) ldsspoly = deg_out;*/ 
		    sspoly = (Mfloat *) imsl_malloc (ldsspoly*sspoly_col_dim*sizeof(Mfloat));
                    if (!sspoly) {
                      imsl_e1sti (1,degree);
                      imsl_e1stl (1,"degree");
                      imsl_e1sti (2,sspoly_col_dim);
                      imsl_e1stl (2,"ssq_poly_col_dim");
		      imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
                      goto RETURN;
		    }
		  }
                 sset(degree, 0.0, &sspoly[0], 1);
                 sset(degree, imsl_amach(6), &sspoly[2*degree], 1);
	      } /*end if (degree != 0)*/
 

		if (!user_sslof && return_sslof) {
                    if (degree != 0) {
		      sslof = (Mfloat *) imsl_malloc (ldsslof*sslof_col_dim*sizeof(Mfloat));
                      if (!sslof) {
                        imsl_e1sti (1,degree);
                        imsl_e1stl (1,"degree");
                        imsl_e1sti (2,sslof_col_dim);
                        imsl_e1stl (2,"ssq_lof_col_dim");
		        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
                        goto RETURN;
		      }
		    }
		} else if (!return_sslof){
		    lof = 0;
		    ldsslof = 1;
		}

                if (degree != 0) {
                  if (return_sslof) {
                    sset(degree, 0.0, &sslof[0], 1);
                    sset(degree, imsl_amach(6), &sslof[degree*2], 1);
		  }
		}

		if (!user_anova){
		    anova = (Mfloat *) imsl_malloc (15*sizeof(Mfloat));
                    if (!anova) {
	              imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
                      goto RETURN;
		  }
		}
		coef = (Mfloat *) imsl_malloc (ldcoef*4*sizeof(Mfloat));
		wk = (Mfloat *) imsl_malloc ((deg_out+1)*(deg_out+7)*sizeof(Mfloat));
		if (!coef || !wk) error = 1;

		if (!error){
		    l_r2tap(deg_out, a, b, smultc, saddc, scoef, d, dfe, sse, lof, dfpe, sspe, iprint, anova, sspoly, ldsspoly, coef, ldcoef, sslof, ldsslof, wk);
 
                    if ((imsl_n1rty(1) > 3) && (imsl_n1rty(1)!=6)) {
		        error = 1;
	            } else {
		        error = 0;
	            }

		    if (!error) {



			if (ssq_pure) *ssq_pure = sspe;
			if (df_pure) *df_pure = (Mint) dfpe;

			scopy (deg_out+1, coef, 1, lv_coefficients, 1);
			for (i=deg_out+1; i<=degree; i++){
			    *(lv_coefficients+i) = F_ZERO;
			}

			dft = (Mfloat) (nobs - nrmiss - 1);

			if (x_mean || x_variance) {
			    Mfloat  mean;
			    imsl_c1div(*a - saddc, smultc, &mean);

			    if (x_mean) *x_mean = mean;

			    if (x_variance) {
				*x_variance = F_ZERO;
				for (i = 0; i < nobs; i++) {
				    if (!imsl_ifnan(x_data[i]) && !imsl_ifnan(y_data[i])) {
					*x_variance += imsl_fi_power(x_data[i]-mean, 2);
				    }
				}
			    }
			    imsl_c1div(*x_variance, dft, x_variance);
			}
		    }
		}
	        if (coef) imsl_free (coef);

		if (!user_sspoly){
		    if (!return_sspoly || error){
			if (sspoly) imsl_free(sspoly);
		    } else {
                        if (degree != 0) {
		          imsl_f_m1ran(sspoly_col_dim,degree,sspoly,sspoly);
                          itemp = sslof_col_dim * deg_out;
                          sset(degree-deg_out, 0.0, &sspoly[itemp], sspoly_col_dim);
                          sset(degree-deg_out, 0.0, &sspoly[itemp+1], sspoly_col_dim);
                          sset(degree-deg_out, imsl_amach(6), &sspoly[itemp+2], sspoly_col_dim);
                          sset(degree-deg_out, imsl_amach(6), &sspoly[itemp+3], sspoly_col_dim);
			}
			*sspoly_ptr = sspoly;
		    }
		    sspoly = NULL;
		} else {
                        if (degree != 0) {
		          imsl_f_m1ran(sspoly_col_dim,degree,sspoly,sspoly);
                          itemp = sslof_col_dim * deg_out;
                          sset(degree-deg_out, 0.0, &sspoly[itemp], sspoly_col_dim);
                          sset(degree-deg_out, 0.0, &sspoly[itemp+1], sspoly_col_dim);
                          sset(degree-deg_out, imsl_amach(6), &sspoly[itemp+2], sspoly_col_dim);
                          sset(degree-deg_out, imsl_amach(6), &sspoly[itemp+3], sspoly_col_dim);
			}
		}
		
		if (!user_sslof){
		    if (!return_sslof || error){
			if (sslof) imsl_free(sslof);
		    } else {
                        if (degree != 0) {                       
			  imsl_f_m1ran(sslof_col_dim,degree,sslof,sslof);
                          itemp = sslof_col_dim * deg_out;
                          sset(degree-deg_out, 0.0, &sslof[itemp], sslof_col_dim);
                          sset(degree-deg_out, 0.0, &sslof[itemp+1], sslof_col_dim);
                          sset(degree-deg_out, imsl_amach(6), &sslof[itemp+2], sslof_col_dim);
                          sset(degree-deg_out, imsl_amach(6), &sslof[itemp+3], sslof_col_dim);
			}
			*sslof_ptr = sslof;
		    }
		    sslof = NULL;
		} else{
                        if (degree != 0) {
		          imsl_f_m1ran(sslof_col_dim,degree,sslof,sslof);
                          itemp = sslof_col_dim * deg_out;
                          sset(degree-deg_out, 0.0, &sslof[itemp], sslof_col_dim);
                          sset(degree-deg_out, 0.0, &sslof[itemp+1], sslof_col_dim);
                          sset(degree-deg_out, imsl_amach(6), &sslof[itemp+2], sslof_col_dim);
                          sset(degree-deg_out, imsl_amach(6), &sslof[itemp+3], sslof_col_dim);
			}
		}

		if (!user_anova){
		    if (!return_anova || error){
			if (anova) imsl_free(anova);
		    } else {
			*anova_ptr = anova;
		    }
		    anova = NULL;
		}
	    }
	}
    }

  RETURN:
    if (wk) imsl_free (wk);
    if (iwk) imsl_free (iwk);
    if (scoef) imsl_free (scoef);
    if (d) imsl_free (d);
    if (b) imsl_free (b);
    if (a) imsl_free (a);
    if (x) imsl_free (x);
    if (error) {
	if (lv_coefficients && !user_coef) {
		imsl_free(lv_coefficients);
		lv_coefficients = NULL;
	}
    }

    return(argptr);
}
/* 
  -----------------------------------------------------------------------
    IMSL Name:  R2ORP/DR2ORP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 24, 1986

    Purpose:    Fit an orthogonal polynomial regression model.

    Usage:      CALL R2ORP (NOBS, NCOL, X, LDX, IRSP, IND, IFRQ, IWT,
                            MAXDEG, ICRIT, CRIT, LOF, NDEG, SMULTC,
                            SADDC, A, B, SCOEF, SSORP, DFE, SSE, DFPE,
                            SSPE, NRMISS, WK, IWK)

    Arguments:  See RFORP

    Remarks:    See RFORP

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r2orp(Mint nobs, Mint ncol, Mfloat *x, Mint ldx, Mint irsp, 
       Mint ind, Mint ifrq, Mint iwt, Mint maxdeg, Mint icrit, Mfloat crit, 
       Mint lof, Mint *ndeg, Mfloat *smultc, Mfloat *saddc, Mfloat a[], 
       Mfloat b[], Mfloat scoef[], Mfloat ssorp[], Mfloat *dfe, Mfloat *sse, 
       Mfloat *dfpe, Mfloat *sspe, Mint *nrmiss, Mfloat wk[], Mint iwk[])
#else
static void l_r2orp(nobs, ncol, x, ldx, irsp, ind, ifrq, iwt, maxdeg, icrit, crit, lof, ndeg, smultc, saddc, a, b, scoef, ssorp, dfe, sse, dfpe, sspe, nrmiss, wk, iwk)
	Mint            nobs, ncol;
	Mfloat         *x;
	Mint            ldx, irsp, ind, ifrq, iwt, maxdeg, icrit;
	Mfloat          crit;
	Mint            lof, *ndeg;
	Mfloat         *smultc, *saddc, a[], b[], scoef[], ssorp[], *dfe,
	               *sse, *dfpe, *sspe;
	Mint           *nrmiss;
	Mfloat          wk[];
	Mint            iwk[];
#endif
{
#define X(I_,J_)	(x+(I_)*(ldx)+(J_))
	Mint            ner;


	imsl_e1psh("l_r2orp");

	if (nobs <= 1) {
		imsl_e1sti(1, nobs);

/*		imsl_ermes(5, 1, "NOBS = %(i1).  The number of observations, NOBS, must be greater than 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LARGER_NOBS_REQUIRED); 
	}
	if (ncol < 1) {
		imsl_e1sti(1, ncol);

/*		imsl_ermes(5, 2, "NCOL = %(i1).  The number of columns in X, NCOL, must be greater than or equal to 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NCOL_MUST_BE_GE_ONE); 
	}
	ner = 4;
	imsl_c1ind(1, irsp, "IRSP", ncol, "NCOL", &ner);
	imsl_c1ind(0, ind, "IND", ncol, "NCOL", &ner);
	imsl_c1ind(0, ifrq, "IFRQ", ncol, "NCOL", &ner);
	imsl_c1ind(0, iwt, "IWT", ncol, "NCOL", &ner);
	if (maxdeg < 0) {
		imsl_e1sti(1, maxdeg);

/*		imsl_ermes(5, 4, "MAXDEG = %(i1).  The maximum degree of the polynomial to be fit, MAXDEG, must be greater than or equal to 0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_MAXDEG_GE_ZERO);
	}
	if (icrit < 0 || icrit > 2) {
		imsl_e1sti(1, icrit);

/*		imsl_ermes(5, 5, "ICRIT = %(i1).  The criterion option, ICRIT, must be either 0, 1, or 2.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_ICRIT_VALUE);
	}
	if (lof < 0 || lof > 1) {
		imsl_e1sti(1, lof);

/*		imsl_ermes(5, 6, "LOF = %(i1).  The lack of fit option, LOF, must be either 0 or
1.");*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_LOF_VALUE);
	}
	if (icrit == 2) {
		if (lof != 1) {
			imsl_e1sti(1, icrit);
			imsl_e1sti(2, lof);

/*			imsl_ermes(5, 7, "ICRIT = %(i1) and LOF = %(i2).  When the criterion option, ICRIT, equals 2, the lack of fit option, LOF, must be equal to 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_ICRIT_OR_LOF_VALU); 
		}
		if (crit < F_ZERO || crit > 100.0) {
			imsl_e1sti(1, icrit);
			imsl_e1str(1, crit);

/*			imsl_ermes(5, 8, "CRIT = %(r1) when ICRIT = %(i1).  The significance level for the lack of fit test, CRIT, must be greater than or equal to 0 and less than or equal to 100 when the criterion option, ICRIT, is equal  to 2.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_CRIT_OR_ICRIT_VALU); 
		}
	}
	if (icrit == 1) {
		if (crit <= F_ZERO || crit > 100.0) {
			imsl_e1sti(1, icrit);
			imsl_e1str(1, crit);

/*			imsl_ermes(5, 9, "CRIT = %(r1) when ICRIT = %(i1).  The R-squared that the fitted polynomial must achieve, CRIT, must be greater than 0 and less than or equal to 100 when the criterion option, ICRIT, is equal to 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_CRIT_IF_ICRIT_EQ_1); 
		}
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	l_r4orp(nobs, x, ldx, irsp, ind, ifrq, iwt, maxdeg, icrit, crit,
	       lof, ndeg, smultc, saddc, a, b, scoef, ssorp, dfe, sse, dfpe,
		   sspe, nrmiss, wk, &wk[nobs], &wk[nobs * 2], &wk[nobs * 3], iwk, &wk[nobs * 4]);

	/* Exit section */
L_9000:
	imsl_e1pop("l_r2orp");
	return;
#undef X
}				/* end of function */

/*
  -----------------------------------------------------------------------
    IMSL Name:  R3ORP/DR3ORP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 10, 1986

    Purpose:    Nuclei of RFORP using orthogonal polynomials for a
                regression fit.

    Usage:      CALL R3ORP (NOBS, X, Y, WT, FRQ, MAXDEG, ICRIT, CRIT,
                            LOF, NDEG, SMULTC, SADDC, A, B, SCOEF,
                            D, DFE, SSE, DFPE, SSPE, IWK, PWK)

    Arguments:  See RFORP

    Remarks:    See RFORP

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r3orp(Mint nobs, Mfloat x[], Mfloat y[], Mfloat wt[], Mfloat frq[],
                    Mint maxdeg, Mint icrit, Mfloat crit, Mint lof, Mint *ndeg,
                    Mfloat *smultc, Mfloat *saddc, Mfloat a[], Mfloat b[],
                    Mfloat scoef[], Mfloat d[], Mfloat *dfe, Mfloat *sse,
                    Mfloat *dfpe, Mfloat *sspe, Mint iwk[], Mfloat *pwk)
#else
static void l_r3orp(nobs, x, y, wt, frq, maxdeg, icrit, crit, lof, ndeg, smultc, saddc, 
                    a, b, scoef, d, dfe, sse, dfpe, sspe, iwk, pwk)
	Mint            nobs;
	Mfloat          x[], y[], wt[], frq[];
	Mint            maxdeg, icrit;
	Mfloat          crit;
	Mint            lof, *ndeg;
	Mfloat         *smultc, *saddc, a[], b[], scoef[], d[], *dfe, *sse, *dfpe, *sspe;
	Mint            iwk[];
	Mfloat         *pwk;
#endif
{
#define PWK(I_,J_)	(pwk+(I_)*(nobs)+(J_))
	Mint            i, i1, i2, i3, i4, icnt, ii, ixmax, ixmin, iycons,
	                j, jj, k, kk, lof1, mxdeg1, ngrp;
	Mfloat          amslf, amspe, crit1, f, prob, r, sn, sslf, ssorpi,
	                sspei, sspf, ssr, sst, sum, sumfrq, swtfrq, sy,
	                sy2, x1, xn, ybar, yn, yy;


	imsl_e1psh("l_r3orp");

	/* Find maximum and minimum X */
	ixmax = l_ismax(nobs, x, 1);
	ixmin = l_ismin(nobs, x, 1);
	if (ixmin == ixmax) {

/*		imsl_ermes(4, 7, "Each value of the independent variable is the same.  The independent variable, X(*,IND), cannot be constant.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_CONSTANT_XVALUES); 
		goto L_9000;
	}
	iycons = 1;
	for (i = 2; i <= nobs; i++) {
		if (y[i - 2] == y[i - 1])
			goto L_10;
		iycons = 0;
		goto L_20;
L_10:
		;
	}
	/* Find distinct groups of X */
L_20:
	iset(nobs, 0, iwk, 1);
	ngrp = 1;
	icnt = 1;
	for (i = 1; i <= (nobs - 1); i++) {
		if (x[i - 1] == x[i]) {
			icnt += 1;
		} else {
			iwk[ngrp - 1] = icnt;
			icnt = 1;
			ngrp += 1;
		}
	}
	iwk[ngrp - 1] = icnt;
	sumfrq = l_ssum(nobs, frq, 1);
	lof1 = lof;
	if (icrit == 0) {
		crit1 = 100.0;
	} else {
		crit1 = crit;
	}
	r = crit1 * .01;
	if (lof == 1) {
		/* Compute pure error. */
		*sspe = F_ZERO;
		k = 1;
		kk = 1;
		for (ii = 1; ii <= ngrp; ii++) {
			sum = F_ZERO;
			swtfrq = F_ZERO;
			sspei = F_ZERO;
			for (jj = 1; jj <= iwk[ii - 1]; jj++) {
				sum += y[k - 1] * wt[k - 1];
				swtfrq += wt[k - 1];
				k += 1;
			}
			ybar = sum / swtfrq;
			for (jj = 1; jj <= iwk[ii - 1]; jj++) {
				sspei += wt[kk - 1] * (y[kk - 1] - ybar) * (y[kk - 1] -
								      ybar);
				kk += 1;
			}
			*sspe += sspei;
		}
		*dfpe = sumfrq - ngrp;
		if (*dfpe == F_ZERO) {
			/* Lack of fit test does not exist. */
			lof1 = 0;
		}
	}
	/* Scale X to (-2,2) */
	yy = F_HALF * (x[ixmin - 1] + x[ixmax - 1]);
	yn = F_FOUR / (x[ixmax - 1] - x[ixmin - 1]);
	i1 = 1;
	i2 = 2;
	i3 = 3;
	i4 = 4;
	sy = imsl_sdot(nobs, y, 1, wt, 1);
	sy2 = imsl_sxyz(nobs, y, 1, y, 1, wt, 1);
	xn = l_ssum(nobs, wt, 1);
	for (i = 1; i <= nobs; i++) {
		*PWK(i4 - 1, i - 1) = sqrt(wt[i - 1]);
	}
	sset(nobs, F_ZERO, PWK(i1 - 1, 0), 1);
	sset(nobs, F_ZERO, PWK(i3 - 1, 0), 1);
	sset(nobs, F_ONE, PWK(i2 - 1, 0), 1);
	for (i = 1; i <= nobs; i++) {
		x[i - 1] = yn * (x[i - 1] - yy);
	}
	x1 = F_ZERO;
	sset(maxdeg + 1, F_ZERO, scoef, 1);
	ssr = F_ZERO;
	/* Find first coefficient */
	d[0] = xn;
	scoef[0] = sy / xn;
	*ndeg = 0;
	mxdeg1 = imsl_i_min(maxdeg + 1, ngrp);
	for (i = 1; i <= mxdeg1; i++) {
		j = i + 1;
		ssorpi = sy * sy / xn;
		sy2 -= ssorpi;
		if (sy2 <= F_ZERO)
			sy2 = F_ZERO;
		if (j <= 2) {
			sst = sy2;
			if (iycons == 1)
				sst = F_ZERO;
			sspf = (F_ONE - F_TEN * imsl_amach(4)) * sst;
		} else {
			ssr += ssorpi;
			if (ssr > sst)
				ssr = sst;
		}
		if (lof1 == 1 && j > 2) {
			/* Compute lack of fit test */
			if (i >= ngrp) {
				prob = F_ONE;
				goto L_90;
			}
			sslf = (sst - ssr) - *sspe;
			if (sslf < F_ZERO)
				sslf = F_ZERO;
			amslf = sslf / (ngrp - i);
			imsl_c1div(*sspe, *dfpe, &amspe);
			l_c1f(amslf, amspe, (float) (ngrp - i), *dfpe, &f, &prob);
		}
		/* Check for stopping criterion */
L_90:
		if (j <= 2) {
			if (sst <= F_ZERO)
				goto L_130;
		} else {
			if (lof1 == 1 && icrit == 2) {
				if (prob > r)
					goto L_130;
			} else {
				if (ssr >= sst * r)
					goto L_130;
				if (ssr >= sspf)
					goto L_130;
			}
		}
		if (i == mxdeg1)
			goto L_130;
		*ndeg = i;
		/*
		 * Compute constants
		 */
		sn = imsl_sxyz(nobs, x, 1, PWK(i4 - 1, 0), 1, PWK(i4 - 1, 0), 1);
		sn /= xn;
		if (i != 1)
			x1 = xn / d[i - 2];
		/* Find next orthogonal polynomial */
		for (k = 1; k <= nobs; k++) {
			*PWK(i1 - 1, k - 1) = (x[k - 1] - sn) ** PWK(i2 - 1, k - 1) -
				x1 ** PWK(i1 - 1, k - 1);
			*PWK(i3 - 1, k - 1) = *PWK(i1 - 1, k - 1) * sqrt(wt[k - 1]);
		}
		a[i - 1] = sn;
		b[i - 1] = x1;
		sy = F_ZERO;
		for (k = 1; k <= nobs; k++) {
			sy += *PWK(i3 - 1, k - 1) * (y[k - 1] - *PWK(i2 - 1, k - 1)) *
				sqrt(wt[k - 1]);
		}
		xn = imsl_sdot(nobs, PWK(i3 - 1, 0), 1, PWK(i3 - 1, 0), 1);
		/* Compute new coefficients */
		d[j - 1] = xn;
		scoef[j - 1] = sy / xn;
		ii = i1;
		i1 = i2;
		i2 = ii;
		ii = i3;
		i3 = i4;
		i4 = ii;
	}
	/*
	 * Compute multiplicative and additive constants, SSE and DFE
	 */
L_130:
	*sse = sy2;
	*smultc = yn;
	*saddc = -yn * yy;
	*dfe = sumfrq - (*ndeg + 1);
	if (*dfe == F_ZERO)
		*sse = F_ZERO;
	if (*sse == F_ZERO || ssr >= sspf) {
		imsl_e1sti(1, *ndeg);

/*		imsl_ermes(6, 11, "The degree %(i1) polynomial fit is a perfect fit within machine precision.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_PERFECT_FIT_POLY);
	}
	if (sst <= F_ZERO) {

/*		imsl_ermes(3, 4, "The response variable, X(*,IRSP), is constant.  A zero degree polynomial is fit.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_CONSTANT_YVALUES); 
		*sse = F_ZERO;
		if (lof == 1)
			*sspe = F_ZERO;
	} else if (*ndeg == ngrp - 1 && *ndeg < maxdeg) {
		imsl_e1sti(1, *ndeg);
		imsl_e1sti(2, ngrp);

/*		imsl_ermes(3, 5, "A degree %(i1) polynomial is fitted.  There are only %(i2) distinct abscissas, so a higher degree polynomial fit cannot be computed.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_FEW_DISTINCT_XVALUES);
	} else if ((icrit == 0 && ssr >= sspf) && *ndeg < maxdeg) {
		/* Perfect fit */
		imsl_e1sti(1, maxdeg);
		imsl_e1sti(2, *ndeg);

/*		imsl_ermes(3, 6, "Although a degree %(i1) fit was requested, a degree %(i2) polynomial is fitted because the degree %(i2) polynomial fit is a perfect fit within machine precision.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_PERFECT_FIT);
	}
	/* Exit section */
L_9000:
	imsl_e1pop("l_r3orp");
	return;
#undef PWK
}				/* end of function */

/*
  -----------------------------------------------------------------------
    IMSL Name:  R4ORP/DR4ORP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 24, 1986

    Purpose:    Fit an orthogonal polynomial regression model.

    Usage:      CALL R4ORP (NOBS, X, LDX, IRSP, IND, IFRQ, IWT, MAXDEG,
                            ICRIT, CRIT, LOF, NDEG, SMULTC, SADDC, A, B,
                            SCOEF, SSORP, DFE, SSE, DFPE, SSPE, NRMISS,
                            XWK, YWK, WWK, FWK, IWK, PWK)

    Arguments:  See RFORP

    Remarks:    See RFORP

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r4orp(Mint nobs, Mfloat *x, Mint ldx, Mint irsp, Mint ind,
                    Mint ifrq, Mint iwt, Mint maxdeg, Mint icrit, Mfloat crit,
                    Mint lof, Mint *ndeg, Mfloat *smultc, Mfloat *saddc,
                    Mfloat a[], Mfloat b[], Mfloat scoef[], Mfloat ssorp[],
                    Mfloat *dfe, Mfloat *sse, Mfloat *dfpe, Mfloat *sspe,
                    Mint *nrmiss, Mfloat xwk[], Mfloat ywk[], Mfloat wwk[],
                    Mfloat fwk[], Mint iwk[], Mfloat pwk[])
#else
static void l_r4orp(nobs, x, ldx, irsp, ind, ifrq, iwt, maxdeg, icrit, crit,
                    lof, ndeg, smultc, saddc, a, b, scoef, ssorp, dfe, sse,
                    dfpe, sspe, nrmiss, xwk, ywk, wwk, fwk, iwk, pwk)
	Mint            nobs;
	Mfloat         *x;
	Mint            ldx, irsp, ind, ifrq, iwt, maxdeg, icrit;
	Mfloat          crit;
	Mint            lof, *ndeg;
	Mfloat         *smultc, *saddc, a[], b[], scoef[], ssorp[], *dfe,
	               *sse, *dfpe, *sspe;
	Mint           *nrmiss;
	Mfloat          xwk[], ywk[], wwk[], fwk[];
	Mint            iwk[];
	Mfloat         *pwk;
#endif
{
#define X(I_,J_)	(x+(I_)*(ldx)+(J_))
#define PWK(I_,J_)	(pwk+(I_)*(nobs)+(J_))
	Mint             i, icall, igo, nmiss, nob1, _l0=0, _l1=1;
	Mfloat           frq, wt;


	imsl_e1psh("l_r4orp");
	/*
	 * Eliminate all missing values and bad values for weights and freqs
	 */
	icall = 1;
	*nrmiss = 0;
	nob1 = 0;
	for (i = 1; i <= nobs; i++) {
		if (!imsl_ifnan(*X(ind - 1, i - 1))) {
			if (!imsl_ifnan(*X(irsp - 1, i - 1))) {
				imsl_c1wfr(_l0, icall, x, ldx, i, _l1, ifrq, iwt, imsl_amach(6), &nmiss, &frq, &wt, &igo);
				if (igo == 3) {
					goto L_9000;
				} else if (igo != 2) {
					nob1 += 1;
					xwk[nob1 - 1] = *X(ind - 1, i - 1);
					*PWK(0, nob1 - 1) = *X(irsp - 1, i - 1);
					*PWK(1, nob1 - 1) = wt;
					*PWK(2, nob1 - 1) = frq;
				} else {
					*nrmiss += 1;
				}
			} else {
				*nrmiss += 1;
			}
		} else {
			*nrmiss += 1;
		}
		icall += 1;
	}
	if (nob1 == 0) {

/*		imsl_ermes(4, 3, "Each row of X contains a missing value in either the IRSP, IFRQ, IND, or IWT column.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_ALL_OBSERVATIONS_MISSING);
		goto L_9000;
	}
	/* Sort X */
	for (i = 1; i <= nob1; i++) {
		iwk[i - 1] = i;
	}
	imsl_svrgp(nob1, xwk, xwk, iwk);
	/*
	 * Put weights, freqs and response variables in proper order along
	 * with X
	 */
	for (i = 1; i <= nob1; i++) {
		ywk[i - 1] = *PWK(0, iwk[i - 1] - 1);
		wwk[i - 1] = *PWK(1, iwk[i - 1] - 1);
		fwk[i - 1] = *PWK(2, iwk[i - 1] - 1);
	}
	for (i = 1; i <= nob1; i++) {
		wwk[i - 1] *= fwk[i - 1];
	}
	/*
	 * Call computational routine
	 */
	l_r3orp(nob1, xwk, ywk, wwk, fwk, maxdeg, icrit, crit, lof, ndeg,
	       smultc, saddc, a, b, scoef, ssorp, dfe, sse, dfpe, sspe, iwk,
		   pwk);

	/* Exit section */
L_9000:
	imsl_e1pop("l_r4orp");
	return;

#undef X
#undef PWK
}				/* end of function */


/*
  -----------------------------------------------------------------------
    IMSL Name:  R2TAP/DR2TAP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 24, 1986

    Purpose:    Compute summary statistics for a polynomial regression
                model given the fit based on orthogonal polynomials.

    Usage:      CALL R2TAP (NDEG, A, B, SMULTC, SADDC, SCOEF, D, DFE,
                            SSE, LOF, DFPE, SSPE, IPRINT, AOV, SQSS,
                            LDSQSS, COEF, LDCOEF, TLOF, LDTLOF, WK)

    Arguments:  See RSTAP.

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r2tap(Mint ndeg, Mfloat a[], Mfloat b[], Mfloat smultc, Mfloat saddc, 
		    Mfloat scoef[], Mfloat d[], Mfloat dfe, Mfloat sse, Mint lof, 
		    Mfloat dfpe, Mfloat sspe, Mint iprint, Mfloat aov[], Mfloat *sqss,
		    Mint ldsqss, Mfloat *coef, Mint ldcoef, Mfloat *tlof, Mint ldtlof,
		    Mfloat wk[])
#else
static void l_r2tap(ndeg, a, b, smultc, saddc, scoef, d, dfe, sse, lof, dfpe, sspe, iprint,
		    aov, sqss, ldsqss, coef, ldcoef, tlof, ldtlof, wk)
	Mint            ndeg;
	Mfloat          a[], b[], smultc, saddc, scoef[], d[], dfe, sse;
	Mint            lof;
	Mfloat          dfpe, sspe;
	Mint            iprint;
	Mfloat          aov[], *sqss;
	Mint            ldsqss;
	Mfloat         *coef;
	Mint            ldcoef;
	Mfloat         *tlof;
	Mint            ldtlof;
	Mfloat          wk[];
#endif
{
#define SQSS(I_,J_)	(sqss+(I_)*(ldsqss)+(J_))
#define COEF(I_,J_)	(coef+(I_)*(ldcoef)+(J_))
#define TLOF(I_,J_)	(tlof+(I_)*(ldtlof)+(J_))
	Mint            i;


	imsl_e1psh("l_r2tap");

	if (ndeg < 0) {
		imsl_e1sti(1, ndeg);

/*		imsl_ermes(5, 1, "NDEG = %(i1).  The degree of the polynomial regression, NDEG, must be greater than or equal to 0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_NDEG_VALUE);
	} else {
		for (i = 1; i <= (ndeg + 1); i++) {
			if (d[i - 1] <= F_ZERO) {
				imsl_e1sti(1, i);
				imsl_e1str(1, d[i - 1]);

/*				imsl_ermes(5, 2, "D(%(i1)) = %(r1).  Each element of D, the vector containing the diagonal elements of the sum of squares and cross-products matrix, must be positive.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_D_ELMNTS);
			}
		}
	}
	if (dfe < F_ZERO) {
		imsl_e1str(1, dfe);

/*		imsl_ermes(5, 3, "DFE = %(r1).  The degrees of freedom for error, DFE, must be greater than or equal to 0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_DFE_VALUE);
	}
	if (sse < F_ZERO) {
		imsl_e1str(1, sse);

/*		imsl_ermes(5, 4, "SSE = %(r1).  The sum of squares for error, SSE, must be greater than or equal to 0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_SSE_VALUE);
	}
	if (lof != 0 && lof != 1) {
		imsl_e1sti(1, lof);

/*		imsl_ermes(5, 5, "LOF = %(i1).  The lack of fit test option, LOF, must be either 0 or 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_LOF_VALUE); 
	}
	if (lof == 1) {
		if (dfpe < F_ZERO) {
			imsl_e1str(1, dfpe);

/*			imsl_ermes(5, 6, "DFPE = %(r1).  The degrees of freedom for pure error, DFPE, must be greater than or equal to 0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_DFPE_VALUE);
		}
		if (sspe < F_ZERO) {
			imsl_e1str(1, sspe);

/*			imsl_ermes(5, 7, "SSPE = %(r1).  The sums of squares for pure error, SSPE, must be greater than or equal to 0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_SSPE_VALUE);
		}
	}
	if (iprint != 0 && iprint != 1) {
		imsl_e1sti(1, iprint);

/*		imsl_ermes(5, 8, "IPRINT = %(i1).  The printing option, IPRINT, must be equal to either 0 or 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_IPRINT_VALUE);
	}
	if (ndeg > 0) {
		if (ldsqss < ndeg) {
			imsl_e1sti(1, ndeg);
			imsl_e1sti(2, ldsqss);

/*			imsl_ermes(5, 9, "NDEG = %(i1) and LDSQSS = %(i2).  The leading dimension of SQSS, LDSQSS, must be greater than or equal to the degree of the polynomial regression, NDEG.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDSQSS_GE_NDEG);
		}
		if (lof == 1) {
			if (ldtlof < ndeg) {
				imsl_e1sti(1, ndeg);
				imsl_e1sti(2, ldtlof);

/*				imsl_ermes(5, 11, "NDEG = %(i1) and LDTLOF = %(i2).  The leading dimension of TLOF, LDTLOF, must be greater than or equal to the degree of the polynomial regression, NDEG.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDTLOF_GE_NDEG);
			}
		}
	}
	if (ndeg >= 0) {
		if (ldcoef < ndeg + 1) {
			imsl_e1sti(1, ndeg);
			imsl_e1sti(2, ldcoef);

/*			imsl_ermes(5, 10, "NDEG = %(i1) and LDCOEF = %(i2).  The leading dimension of COEF, LDCOEF, must be greater than or equal to 1 plus the degree of the polynomial regression, NDEG.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LARGER_LDCOEF);
		}
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;

	l_r5tap(ndeg, a, b, smultc, saddc, scoef, d, dfe, sse, lof, dfpe,
		sspe, iprint, aov, sqss, ldsqss, coef, ldcoef, tlof, ldtlof,
		   wk, &wk[(ndeg + 1) * 4]);

L_9000:
	imsl_e1pop("l_r2tap");
	return;
#undef SQSS
#undef COEF
#undef TLOF
}				/* end of function */


/*
  -----------------------------------------------------------------------
    IMSL Name:  R3TAP/DR3TAP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 10, 1986

    Purpose:    Nuclei of RSTAP used to uncode the coefficients for
                the orthogonal polynomials for the original data.

    Usage:      CALL R3TAP (NDEG, COEF, A, B, SMULTC, SADDC, PWK, LDPWK)

    Arguments:
       NDEG   - Degree of the final model.  (Input)
       COEF   - Vector of length NDEG+1.  (Input/Output)
                On input COEF contains the coefficients for the
                scaled model.  On output COEF contains the coefficients
                for the original data.
       A      - See RSTAP.
       B      - See RSTAP.
       SMULTC - Multiplicative constant used to compute a scaled version
                of the data.  (Input)
       SADDC  - Additive constant used to compute a scaled version of the
                data.  (Input)
       PWK    - NDEG+1 by 4 work matrix.
       LDPWK  - Leading dimension of PWK exactly as specified in the
                dimension statement of the calling program.  (Input)

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r3tap(Mint ndeg, Mfloat coef[], Mfloat a[], Mfloat b[], Mfloat smultc, Mfloat saddc, Mfloat *pwk, Mint ldpwk)
#else
static void l_r3tap(ndeg, coef, a, b, smultc, saddc, pwk, ldpwk)
	Mint             ndeg;
	Mfloat           coef[], a[], b[], smultc, saddc, *pwk;
	Mint             ldpwk;
#endif
{
#define PWK(I_,J_)	(pwk+(I_)*(ldpwk)+(J_))
	Mint             i, ii, ii1, j, k;
	Mfloat           s, skl, sl, sss, tt;

	/*
	 * Perform partial decoding - compute coefficients for the model
	 * expressed as a polynomial function of the scaled independent
	 * variable
	 */
	scopy(ndeg + 1, coef, 1, PWK(3, 0), 1);
	for (i = 1; i <= 3; i++) {
		sset(ndeg + 1, F_ZERO, PWK(i - 1, 0), 1);
	}
	for (i = 2; i <= (ndeg + 1); i++) {
		*PWK(1, i - 1) = F_ONE;
		k = i - 1;
		for (ii = 2; ii <= i; ii++) {
			ii1 = ii - 1;
			*PWK(2, ii - 1) = *PWK(1, ii1 - 1) - *PWK(1, ii - 1) * a[k - 1] -
				b[k - 1] ** PWK(0, ii - 1);
			*PWK(3, ii1 - 1) += coef[i - 1] ** PWK(2, ii - 1);
		}
		if (i == ndeg + 1)
			goto L_50;
		for (ii = 1; ii <= i; ii++) {
			*PWK(0, ii - 1) = *PWK(1, ii - 1);
			*PWK(1, ii - 1) = *PWK(2, ii - 1);
		}
	}
	/*
	 * Finish decoding - compute coefficients for the model expressed as
	 * a polynomial function of the original independent variable
	 */
L_50:
	sss = F_ONE;
	s = F_ONE;
	sl = F_ONE;
	skl = F_ONE;
	for (j = 1; j <= (ndeg + 1); j++) {
		tt = sss ** PWK(3, j - 1);
		k = j + 1;
		skl = j;
		sl = F_ONE;
		for (i = k; i <= (ndeg + 1); i++) {
			s = saddc * s * skl / sl;
			tt += *PWK(3, i - 1) * s;
			if (j == 1)
				goto L_60;
			skl += F_ONE;
			sl += F_ONE;
	L_60:
			;
		}
		coef[j - 1] = tt;
		sss *= smultc;
		s = sss;
	}
	return;

#undef PWK
}				/* end of function */

/*
  -----------------------------------------------------------------------
    IMSL Name:  R4TAP/DR4TAP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 7, 1986

    Purpose:    Nuclei of RSTAP used to compute variance estimates
                for the uncoded orthogonal polynomial regression
                coefficients.

    Usage:      CALL R4TAP (NDEG, COEF, A, B, SMULTC, SADDC, TWK, LDTWK)

    Arguments:
       NDEG   - Degree of the final model.  (Input)
       COEF   - Vector of length NDEG+1.  (Input/Output)
                On input COEF contains the variances for the
                coded model.  On output COEF contains the variances
                for the uncoded model.
       A      - See RSTAP.
       B      - See RSTAP.
       SMULTC - Multiplicative constant used to compute a scaled version
                of the data.  (Input)
       SADDC  - Additive constant used to compute a scaled version of the
                data.  (Input)
       TWK    - NDEG+1 by NDEG+3 work matrix.
       LDTWK  - Leading dimension of TWK exactly as specified in the
                dimension statement of the calling program.  (Input)

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r4tap(Mint ndeg, Mfloat coef[], Mfloat a[], Mfloat b[], Mfloat smultc, Mfloat saddc, Mfloat *twk, Mint ldtwk)
#else
static void l_r4tap(ndeg, coef, a, b, smultc, saddc, twk, ldtwk)
	Mint            ndeg;
	Mfloat          coef[], a[], b[], smultc, saddc, *twk;
	Mint            ldtwk;
#endif
{
#define TWK(I_,J_)	(twk+(I_)*(ldtwk)+(J_))
	Mint            i, j, k=0, k1, l;
	Mfloat          imsl_beta, s, skl, sl, ss, sss, sum;


	if (ndeg > 0) {
		*TWK(1, 0) = -a[0];
	}
	for (i = 1; i <= (ndeg + 1); i++) {
		/*
		 * Calculate coefficients of the coded regression coefficient
		 * estimates, used to obtain partially decoded estimates
		 */
		*TWK(i - 1, i - 1) = F_ONE;
		if (i > 2) {
			*TWK(i - 1, 0) = -*TWK(k - 1, 0) * a[k - 1] - *TWK(i - 3, 0) *
				b[k - 1];
		}
		k = i;
	}
	for (j = 2; j <= ndeg; j++) {
		k = j + 1;
		*TWK(k - 1, j - 1) = *TWK(j - 1, j - 2) - a[j - 1];
		for (i = j + 2; i <= (ndeg + 1); i++) {
			*TWK(i - 1, j - 1) = *TWK(k - 1, j - 2) - a[k - 1] ** TWK(k - 1, j - 1) -
				b[k - 1] ** TWK(k - 2, j - 1);
			k = i;
		}
	}
	/*
	 * Calculate variances and covariances of partially decoded estimates
	 */
	for (i = 1; i <= (ndeg + 1); i++) {
		*TWK(ndeg + 1, i - 1) = coef[i - 1] + imsl_sxyz(ndeg - i + 1, &coef[i],
			     1, TWK(i, i - 1), ldtwk, TWK(i, i - 1), ldtwk);
	}
	for (i = 1; i <= ndeg; i++) {
		for (j = i + 1; j <= (ndeg + 1); j++) {
			*TWK(i - 1, j - 1) = imsl_sxyz(ndeg - j + 2, &coef[j - 1],
						       1, TWK(j - 1, i - 1), ldtwk, TWK(j - 1, j - 1), ldtwk);
		}
	}
	/*
	 * Calculate coefficients of the partially decoded regression
	 * coefficient estimates, used to obtain completely decoded estimates
	 */
	imsl_beta = smultc;
	sss = F_ONE;
	ss = saddc;
	s = F_ONE;
	sl = F_ONE;
	skl = F_ONE;
	for (j = 1; j <= (ndeg + 1); j++) {
		k1 = j;
		k = j + 1;
		*TWK(k1 - 1, j - 1) = sss;
		if (j != 1) {
			skl = j;
			sl = F_ONE;
		}
		for (i = k; i <= (ndeg + 1); i++) {
			s = ss * s * skl / sl;
			k1 += 1;
			*TWK(k1 - 1, j - 1) = s;
			if (j == 1)
				goto L_70;
			skl += F_ONE;
			sl += F_ONE;
	L_70:
			;
		}
		sss *= imsl_beta;
		s = sss;
	}
	/*
	 * CALCULATE VARIANCES OF COMPLETELY DECODED REGRESSION COEFFICIENTS
	 */
	for (i = 1; i <= (ndeg + 1); i++) {
		*TWK(ndeg + 2, i - 1) = imsl_sxyz(ndeg - i + 2, TWK(ndeg + 1, i - 1),
		     1, TWK(i - 1, i - 1), ldtwk, TWK(i - 1, i - 1), ldtwk);
	}
	for (i = 1; i <= ndeg; i++) {
		sum = F_ZERO;
		for (l = i; l <= ndeg; l++) {
			for (j = l + 1; j <= (ndeg + 1); j++) {
				sum += *TWK(l - 1, i - 1) ** TWK(j - 1, i - 1) ** TWK(l - 1, j - 1);
			}
		}
		coef[i - 1] = *TWK(ndeg + 2, i - 1) + sum + sum;
	}
	coef[ndeg] = *TWK(ndeg + 2, ndeg);

	return;
#undef TWK
}				/* end of function */


/*
  -----------------------------------------------------------------------
    IMSL Name:  R5TAP/DR5TAP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 24, 1986

    Purpose:    Compute summary statistics for a polynomial regression
                model given the fit based on orthogonal polynomials.

    Usage:      CALL R5TAP (NDEG, A, B, SMULTC, SADDC, SCOEF, D, DFE,
                            SSE, LOF, DFPE, SSPE, IPRINT, AOV, SQSS,
                            LDSQSS, COEF, LDCOEF, TLOF, LDTLOF, PWK, TWK)

    Arguments:  See RSTAP.

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* IPRINT was not used, but leave the calling sequence intact. */
#ifdef ANSI
static void l_r5tap(Mint ndeg, Mfloat a[], Mfloat b[], Mfloat smultc, 
		    Mfloat saddc, Mfloat scoef[], Mfloat d[], Mfloat dfe, 
		    Mfloat sse, Mint lof, Mfloat dfpe, Mfloat sspe,
		    Mint iprint, Mfloat aov[], Mfloat *sqss, Mint ldsqss,
		    Mfloat *coef, Mint ldcoef, Mfloat *tlof, Mint ldtlof,
		    Mfloat pwk[], Mfloat twk[])
#else
static void l_r5tap(ndeg, a, b, smultc, saddc, scoef, d, dfe, sse, lof, 
		    dfpe, sspe, iprint, aov, sqss, ldsqss, coef, ldcoef, 
		    tlof, ldtlof, pwk, twk)
	Mint            ndeg;
	Mfloat          a[], b[], smultc, saddc, scoef[], d[], dfe, sse;
	Mint            lof;
	Mfloat          dfpe, sspe;
	Mint            iprint;
	Mfloat          aov[], *sqss;
	Mint            ldsqss;
	Mfloat         *coef;
	Mint            ldcoef;
	Mfloat         *tlof;
	Mint            ldtlof;
	Mfloat          pwk[], twk[];
#endif
{
#define SQSS(I_,J_)	(sqss+(I_)*(ldsqss)+(J_))
#define COEF(I_,J_)	(coef+(I_)*(ldcoef)+(J_))
#define TLOF(I_,J_)	(tlof+(I_)*(ldtlof)+(J_))

	Mint            nan, i, ldpwk, ldtwk;
	Mfloat          amse, amslf, amspe, dfr, gmean, s2, sqms, ssr;

	/* Compute sequential SS and AOV vector */
	for (i = 1; i <= ndeg; i++) {
		*SQSS(1, i - 1) = d[i] * imsl_fi_power(scoef[i], 2);
	}
	ssr = l_ssum(ndeg, SQSS(1, 0), 1);
	dfr = ndeg;
	gmean = scoef[0];
	imsl_g1aov(dfr, ssr, dfe, sse, gmean, aov);

	/* Compute uncoded coefficients */
	ldpwk = ndeg + 1;
	scopy(ndeg + 1, scoef, 1, coef, 1);
	l_r3tap(ndeg, coef, a, b, smultc, saddc, pwk, ldpwk);

	/*
	 * Compute variances of the coded coefficients
	 */
	imsl_c1div(sse, dfe, &s2);
	nan = imsl_ifnan(s2);
	if (!nan) {
		for (i = 1; i <= (ndeg + 1); i++) {
			*COEF(1, i - 1) = (F_ONE / d[i - 1]) * s2;
		}
	}
	/*
	 * Compute variances of the uncoded coefficients
	 */
	if (nan) {
		sset(ndeg + 1, imsl_amach(6), COEF(1, 0), 1);
	} else {
		ldtwk = ndeg + 1;
		l_r4tap(ndeg, COEF(1, 0), a, b, smultc, saddc, twk, ldtwk);
		for (i = 1; i <= (ndeg + 1); i++) {
			if (*COEF(1, i - 1) > F_ZERO) {
				*COEF(1, i - 1) = sqrt(*COEF(1, i - 1));
			} else {
				*COEF(1, i - 1) = F_ZERO;
			}
		}
	}
	/* Compute COEF matrix */
	if (!nan) {
		l_rcoef(ndeg + 1, coef, COEF(1, 0), dfe, coef, ldcoef);
	} else {
		sset(ndeg + 1, imsl_amach(6), COEF(2, 0), 1);
		sset(ndeg + 1, imsl_amach(6), COEF(3, 0), 1);
	}

	if (ndeg > 0) {
		/* Compute SQSS matrix */
		sset(ndeg, F_ONE, sqss, 1);
		imsl_c1div(sse, dfe, &amse);
		for (i = 1; i <= ndeg; i++) {
			imsl_c1div(*SQSS(1, i - 1), *SQSS(0, i - 1), &sqms);
			l_c1f(sqms, amse, *SQSS(0, i - 1), dfe, SQSS(2, i - 1), SQSS(3, i - 1));
		}
	}
	if (lof == 1 && ndeg > 0) {
		/* COMPUTE LOF TEST */
		*TLOF(0, ndeg - 1) = dfe - dfpe;
		*TLOF(1, ndeg - 1) = sse - sspe;
		for (i = ndeg - 1; i >= 1; i--) {
			*TLOF(0, i - 1) = *TLOF(0, i) + F_ONE;
			*TLOF(1, i - 1) = *TLOF(1, i) + *SQSS(1, i);
		}
		imsl_c1div(sspe, dfpe, &amspe);
		for (i = 1; i <= ndeg; i++) {
			imsl_c1div(*TLOF(1, i - 1), *TLOF(0, i - 1), &amslf);
			l_c1f(amslf, amspe, *TLOF(0, i - 1), dfpe, TLOF(2, i - 1),
				 TLOF(3, i - 1));
		}
	}
	return;
#undef SQSS
#undef COEF
#undef TLOF
}				/* end of function */

 
/*
  -----------------------------------------------------------------------
    IMSL Name:  C1F/DC1F (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 2, 1985

    Purpose:    Perform a divide check, and output an F statistic and
                p-value.

    Usage:      CALL C1F (FNUM, FDEN, DF1, DF2, F, PVALUE)

    Arguments:
       FNUM   - Numerator of the F statistic.  (Input)
       FDEN   - Denominator of the F statistic.  (Input)
       DF1    - Degrees of freedom for numerator.  (Input)
       DF2    - Degrees of freedom for denominator.  (Input)
       F      - F statistic.  (Output)
       PVALUE - P-value.  (Output)

    Chapter:    STAT/LIBRARY Utilities (not documented)

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c1f(Mfloat fnum, Mfloat fden, Mfloat df1, Mfloat df2, Mfloat *f, Mfloat *pvalue)
#else
static void l_c1f(fnum, fden, df1, df2, f, pvalue)
	Mfloat           fnum, fden, df1, df2, *f, *pvalue;
#endif
{
	Mfloat           anan, big;


	anan = imsl_amach(6);
	imsl_c1div(fnum, fden, f);
	if ((imsl_ifnan(*f) || df1 <= F_ZERO) || df2 <= F_ZERO) {
		*pvalue = anan;
	} else{
	    big = imsl_amach(7);
	    if (*f == big) {
		*pvalue = F_ZERO;
	    } else {
		*pvalue = F_ONE - imsl_f_F_cdf(*f, df1, df2);
	    }
	}

	return;
}				/* end of function */


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
static Mfloat l_ssum(Mint n, Mfloat sx[], Mint incx)
#else
static Mfloat l_ssum(n, sx, incx)
	Mint             n;
	Mfloat           sx[];
	Mint             incx;
#endif
{
	Mint             _d_l, _d_m, _do0, _do1, i, nincx;
	Mfloat           ssum_v;


	ssum_v = F_ZERO;
	if (n > 0) {
		if (incx != 1) {
			/* CODE FOR INCREMENT NOT EQUAL TO 1 */
			nincx = n * incx;
			for (i = 1, _do0 = DOCNT(1, nincx, _do1 = incx); _do0 > 0; i += _do1, _do0--) {
				ssum_v += sx[i - 1];
			}
		} else {
			for (i = 1; i <= n; i++) {
				ssum_v += sx[i - 1];
			}
		}
	}
	return (ssum_v);
}				/* end of function */

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
	Mint             n;
	Mfloat           sx[];
	Mint             incx;
#endif
{
	Mint             i, ismax_v, ix;
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

/*
  -----------------------------------------------------------------------
    IMSL Name:  ISMIN (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Find the smallest index of the component of a
                single-precision vector having minimum value.

    Usage:      ISMIN(N, SX, INCX)

    Arguments:
       N      - Length of vector X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than zero.
       ISMIN  - The smallest index I such that X(I) is the minimum of
                of X(J) for J=1 to N.  (Output)
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
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  RCOEF/DRCOEF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the t-values and p-values for the regression
                coefficients given their estimates and their estimated
                standard errors.

    Usage:      CALL RCOEF (NCOEF, B, STDB, DFE, COEF, LDCOEF)

    Arguments:
       NCOEF  - Number of regression coefficients in the model.
                (Input)
       B      - Vector of length NCOEF containing a least squares
                solution for the regression coefficients.  (Input)
       STDB   - Vector of length NCOEF containing the estimated
                standard errors of the estimated regression coefficients.
                (Input)
       DFE    - Degrees of freedom for error.  (Input)
       COEF   - Matrix, NCOEF by 4, containing B in column 1, STDB in
                column 2, corresponding t-values in column 3, and
                corresponding p-values in column 4.  (Output)
                The t-values computed are for the null hypothesis that
                the corresponding coefficient is zero.  The p-values
                correspond to two-sided alternative hypotheses, i.e.,
                the p-value is the probability of a t random variable
                exceeding the computed t-value in column 3 in absolute
                value.  B and column 1 of COEF can share the same
                storage locations.  STDB and column 2 of COEF can
                share the same storage locations.
       LDCOEF - Leading dimension of COEF exactly as specified in the
                dimension statement of the calling program.

    GAMS:       L8h

    Chapter:    STAT/LIBRARY Regression

    Copyright:  1985 by IMSL, Inc.  All rights reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_rcoef(Mint ncoef, Mfloat b[], Mfloat stdb[], Mfloat dfe, Mfloat *coef, Mint ldcoef)
#else
static void l_rcoef(ncoef, b, stdb, dfe, coef, ldcoef)
	Mint            ncoef;
	Mfloat          b[], stdb[], dfe, *coef;
	Mint            ldcoef;
#endif
{
#define COEF(I_,J_)	(coef+(I_)*(ldcoef)+(J_))
	Mint            i, ner;


	imsl_e1psh("l_rcoef");
	ner = 1;

	imsl_c1dim(1, ncoef, "NCOEF", ldcoef, "LDCOEF", &ner);
	if (dfe < F_ZERO) {
		imsl_e1str(1, dfe);

/*       		imsl_ermes(5, ner, "DFE = %(r1).  It must be nonnegative.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_DFE_VALUE);
	}
	ner += 1;
	for (i = 1; i <= ncoef; i++) {
		if (!imsl_ifnan(stdb[i - 1])) {
			if (stdb[i - 1] < F_ZERO) {
				imsl_e1sti(1, i);
				imsl_e1str(1, stdb[i - 1]);

/*				imsl_ermes(5, ner, "STDB(%(i1)) = %(r1).  It must be nonnegative.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_STDB_VALUE);
			}
		}
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	scopy(ncoef, b, 1, COEF(0, 0), 1);
	scopy(ncoef, stdb, 1, COEF(1, 0), 1);
	for (i = 1; i <= ncoef; i++) {
		l_c1t(*COEF(0, i - 1), *COEF(1, i - 1), dfe, COEF(2, i - 1), COEF(3, i - 1));
	}
L_9000:
	imsl_e1pop("l_rcoef");

	return;

#undef COEF
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  C1T/DC1T (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 2, 1985

    Purpose:    Perform a divide check, and output a t statistic and
                p-value.

    Usage:      CALL C1T (TNUM, TDEN, DF, T, PVALUE)

    Arguments:
       TNUM   - Numerator of the t statistic.  (Input)
       TDEN   - Denominator of the t statistic.  (Input)
                TDEN must be greater than or equal to 0.
       DF     - Degrees of freedom.  (Input)
                Fractional degrees of freedom, as necessary in the
                Behrens-Fisher problem, are permitted.
       T      - T statistic.  (Output)
       PVALUE - P-value.  (Output)

    Chapter:    STAT/LIBRARY Utilities (not documented)

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c1t(Mfloat tnum, Mfloat tden, Mfloat df, Mfloat *t, Mfloat *pvalue)
#else
static void l_c1t(tnum, tden, df, t, pvalue)
	Mfloat           tnum, tden, df, *t, *pvalue;
#endif
{
	Mfloat           abst, temp;


	imsl_c1div(tnum, tden, t);
	if (imsl_ifnan(*t) || df == F_ZERO) {
		*pvalue = imsl_amach(6);
	} else {
	    temp = imsl_amach(8);
	    if ((fabs(*t) == imsl_amach(7) && *t > F_ZERO) || (fabs(*t) == -temp &&
							*t < F_ZERO)) {
		*pvalue = F_ZERO;
	    } else {
		abst = fabs(*t);
		*pvalue = F_TWO * (F_ONE - imsl_f_t_cdf(abst, df));
	    }
	}

	return;
}				/* end of function */
