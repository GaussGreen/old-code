/* ======================================================================= 

	FILENAME:   srt_f_twofaccorrel.c

	PURPOSE:    relevant routines for a two factor correlation calibration 

	AUTHOR:     O. Van Eyseren and A. Savine

	FUNCTIONS:  Err     stdev_correl(...)
	            double  model_correl(...)
	            Err     interp_corr_tenor(...)
   
   ======================================================================= */
#include "math.h"
#include "srt_h_all.h"
#include "grf_h_all.h"
#include "srtaccess.h"

/* -----------------------------------------------------------------------
	Given model parameters (alpha, beta, rho)
	returns the correlation implied by the model
   ----------------------------------------------------------------------- */
double model_correl (
		double mat1, 
		double mat2, 
		double alpha, 
		double beta, 
		double rho)
{

	double normvar1;
	double normvar2;
	double normcovar;
	double correl;

	normvar1 = 	1 
		+ 2*rho*alpha*exp(-beta*mat1) 
		+ alpha*alpha*exp(-2*beta*mat1) ;
	normvar2 = 1 
		+ 2*rho*alpha*exp(-beta*mat2) 
		+ alpha*alpha*exp(-2*beta*mat2) ;
	normcovar = 1 
		+ rho*alpha*(exp(-beta*mat1)+exp(-beta*mat2))
			+ alpha*alpha*exp(-beta*(mat1+mat2)); 

	correl = normcovar / sqrt( normvar1 * normvar2);

	return correl;
}

/* ========================================================================= */

/* -----------------------------------------------------------------------
	Given model parameters (alpha, beta, rho)
	Fills the full correlation matrix (model_corr) implied by the model, 
	for a given set of ifr tenor matrurities (in years)
   ----------------------------------------------------------------------- */

Err model_corr_matrix(
	double **model_corr,
	double *tenor_mat,
	long num_tenor,
	double alpha,
	double beta,
	double rho)
{		
	long i, j;
	
	if (!model_corr)
		return serror("Corr matrix not initialised in model_corr_matrix");
	for (i=0;i<num_tenor;i++)
	{
		for (j=0; j<i;j++)
		{
			model_corr[i][j] = model_correl(
					tenor_mat[i],
					tenor_mat[j],
					alpha,
					beta,
					rho);
			model_corr[j][i] = model_corr[i][j];
		}
	}
	
	return NULL;
}

/* -----------------------------------------------------------------------
	Given model parameters (alpha, beta, rho) and a set of dates (in Y)
	return the confidence level explain by the two eigen vectors 
   ----------------------------------------------------------------------- */

Err TwoFactorEigenConf(double	*pdSetOfDates,
					   int		iNumberOfDates,
					   double	dAlpha,
					   double	dGamma,
					   double	dRho,
					   /* return a vector of two elements */
					   double	pdConfidenceLevel[2])
{
double		**CorrelMatrix= NULL,
			*EigenValues = NULL,
			**EigenVectors = NULL;
int			i;
Err			err;

	CorrelMatrix = dmatrix(0,
						   iNumberOfDates-1,
						   0,
						   iNumberOfDates-1);   
	EigenVectors = dmatrix(0,
						   iNumberOfDates-1,
						   0,
						   iNumberOfDates-1);   

	EigenValues  = (double *)malloc(sizeof(double) * (iNumberOfDates-1));

	if (err = model_corr_matrix(CorrelMatrix,
								pdSetOfDates,
								iNumberOfDates-1,
								dAlpha,
								dGamma,
								dRho))
		return err;

	/*  force the diagonal to be equal to 1 */
	for(i=0; i<iNumberOfDates; i++)
		CorrelMatrix[i][i] = 1.0;
	
	if ( err=jacobi_diagonalisation(CorrelMatrix , 
									iNumberOfDates, 
									EigenValues, 
									EigenVectors, 
									&i))
		return err;

	pdConfidenceLevel[0] = EigenValues[0] / (EigenValues[0]+EigenValues[1]);
	pdConfidenceLevel[1] = EigenValues[1] / (EigenValues[0]+EigenValues[1]);
	
	return err;
}





/* ======================================================================== */

/* ------------------------------------------------------------------------
	Returns in chisq(sum of errors squared ), the difference 
	between model and historical correlation matrix, for given
	model parameters (alpha, beta, rho) and ifr tenors 
   ------------------------------------------------------------------------ */

Err chisq_correl (
		double **hist_corr,
		double *tenor_mat,
		long num_tenor,		
		double alpha,
		double beta,
		double rho, 
		double *chisq)
{
	long i, j;
	Err err = NULL;
	double diff;

	double **model_corr;

	model_corr = dmatrix(0, num_tenor-1, 0, num_tenor-1);
	err = model_corr_matrix(
			model_corr,
			tenor_mat,
			num_tenor,
			alpha,
			beta,
			rho);
	if (err)
	{
		free_dmatrix(model_corr, 0, num_tenor-1, 0, num_tenor-1);
		return err;
	}
	
	*chisq = 0.0;
	for (i=0;i<num_tenor-1;i++)
	{
		for (j=0;j<num_tenor-1;j++)
		{
			diff = model_corr[i][j]-hist_corr[i][j];
			*chisq += diff * diff;
		}
	}		

	free_dmatrix(model_corr, 0, num_tenor-1, 0, num_tenor-1);
	
	return NULL;
}


/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------
	Returns in stdev (sqrt of average of errors squared ), the difference 
	between model and historical correlation matrix, for given
	model parameters (alpha, beta, rho) and ifr tenors 
   ------------------------------------------------------------------------ */

Err stdev_correl (
		double **hist_corr,
		double *tenor_mat,
		long num_tenor,		
		double alpha,
		double beta,
		double rho,
		double *stdev)
{
	double chisq;
	Err err = NULL;

	err = chisq_correl ( hist_corr, tenor_mat, num_tenor, 
			alpha, beta, rho, &chisq);
	if (err)
		return err;

	*stdev = sqrt (chisq / (num_tenor * num_tenor));

	return NULL;
}

/* ------------------------------------------------------------------------ */

/* -----------------------------------------------------------------------
	Given a tenor string (1m, 1y, 1y6m)
	returns the corresponding maturity for the model
   ----------------------------------------------------------------------- */

Err interp_corr_tenor(String fra_tenor, double *fra_mat)
{
	int num1 = 0;
	char ten1;
	int num2 = 0;
	char ten2;
	String endstr;

/* Make them capital letters without blanks */
	strupper(fra_tenor);
	strip_white_space(fra_tenor);

/* First number and first letter (Y or M) */
	num1 = (int)strtol(fra_tenor, &endstr, 10);
	if ( !(*endstr) || (num1<0))
		return serror ("Bad CORR tenor %s",fra_tenor);
	ten1 = *endstr;
	if ( ten1 != 'M' && ten1 != 'Y')
		return serror ("Bad CORR tenor %s",fra_tenor);
	endstr++;

/* Second number and second letter (M) */
	if (*endstr)
	{
		num2 = (int)strtol(endstr, &endstr, 10);
		if ( !(*endstr) || (num2<0))
			return serror ("Bad CORR tenor %s",fra_tenor);
		ten2 = *endstr;
		if ( ten2 != 'M')
			return serror ("Bad CORR tenor %s",fra_tenor);
  	endstr++;
	}

/* Check there is nothing after */
	if (*endstr)
			return serror ("Bad CORR tenor %s: %s not allowed",fra_tenor, endstr);

/* Sets the mat as being number of years plus number of month / 12.0 */
	if (ten1 == 'M')
	{
		*fra_mat = (double)num1 / 12.0;
	}
	else
	{
		*fra_mat = (double)num1 + (double)num2 / 12.0;
	}
	 
/* Return a success message */
	return NULL;
	
} /* END interp_corr_tenor() */ 

/* ======================================================================= */

Err srt_f_twofac_ifr_correl (
		String    tenor1, 
		String    tenor2, 
		String    und_name,
		double    *correl)

{
	Err                  err;
	SrtUndPtr            und;
	
	double               mat1;
	double               mat2;
	double               alpha;
	double               beta;
	double               rho;
	TermStruct           *ts;
	SrtUnderlyingType    und_type;
	SrtMdlType	         mdl_type;
	SrtMdlDim	         mdl_dim;
	
	
	if (!(und = lookup_und (und_name)))
	{
		return serror(" Couldn't find underlying %s ", und_name);
	}
	
	und_type = get_underlying_type(und);

	if (und_type != INTEREST_RATE_UND)
	{
		return serror(" Cannot display correlation of a non IR underlying");
	}
		
	if( err = get_underlying_mdldim(und,&mdl_dim))
	{
		return err;
	}

	if (mdl_dim == ONE_FAC)
	{
		*correl = 1.0;
		return NULL;
	}
	
	if( err = get_underlying_mdltype(und,&mdl_type))
	{
		return err;
	}

	
	if (err = get_underlying_ts (und, &ts))
	{
		return err;
	}

	alpha = ((TwoFacIrmTermStructVal*) (ts->head->element->val.pval))->exp_ts[1].sig
			/ ((TwoFacIrmTermStructVal*) (ts->head->element->val.pval))->exp_ts[0].sig;

	beta = (1.0 / ((TwoFacIrmTermStructVal*) (ts->head->element->val.pval))->exp_ts[1].tau)
			- (1.0 / ((TwoFacIrmTermStructVal*) (ts->head->element->val.pval))->exp_ts[0].tau);

	rho = ((TwoFacIrmTermStructVal*) (ts->head->element->val.pval))->rho;

	if (err = interp_corr_tenor (tenor1, &mat1))
	{
		return err;
	}
	
	if (err = interp_corr_tenor (tenor2, &mat2))
	{
		return err;
	}

	*correl = model_correl (mat1, mat2, alpha, beta, rho);

	return NULL;
}


/* ======================================================================= */

	
/*

  Given any model, returns the normal correl between swap rate 1 and swap rate 2
  senn at some observation date through a grfn tableau.

*/

Err swap_correl (Date obs_date, 
				 Date start1, Date end1, String cmp1, String basis1, 
				 Date start2, Date end2, String cmp2, String basis2,  
				 String und, 
				 int mdl_rows, char** paramStrings, char **valueStrings,
				 String n_ln,
				 double *correl)
{

	char	***tableau = smatrix_size (0, 0, 0, 1, GRFN_DEF_ARGBUFSZ),
			swap1[GRFN_DEF_ARGBUFSZ],
			swap2[GRFN_DEF_ARGBUFSZ];
	double	**ss = dmatrix (0, 0, 0, 1),
			temp_price,
			temp_stdev;
	long	i;
	
	double	e1 = 0.0,
			e2 = 0.0,
			var1 = 0.0,
			var2 = 0.0,
			covar = 0.0;

	int		**mask = imatrix (0, 0, 0, 1); 
	
	SrtHistData		*S1 = NULL , *S2 = NULL;

	Err		err;

	mask[0][0] = mask[0][1] = GRFNSCELL; 

	strupper(cmp1);
	strip_white_space(cmp1);
	strupper(cmp2);
	strip_white_space(cmp2);
	strupper(basis1);
	strip_white_space(basis1);
	strupper(basis2);
	strip_white_space(basis2);
	strupper(n_ln);
	strip_white_space(n_ln);

	if (*n_ln == 'N')
	{
		sprintf (swap1, "SWAP(%d,%d,\"%s\",\"%s\")", 
			(int) start1, (int) end1, cmp1, basis1);
		sprintf (swap2, "SWAP(%d,%d,\"%s\",\"%s\")", 
			(int) start2, (int) end2, cmp2, basis2);
	}
	else if (*n_ln == 'L')
	{
		sprintf (swap1, "LOG(SWAP(%d,%d,\"%s\",\"%s\"))", 
			(int) start1, (int) end1, cmp1, basis1);
		sprintf (swap2, "LOG(SWAP(%d,%d,\"%s\",\"%s\"))", 
			(int) start2, (int) end2, cmp2, basis2);
	}
	else return serror ("parameter %s should be normal or lognormal", n_ln);

	strcpy (tableau[0][0], "(");
	strcat (tableau[0][0], swap1);
	strcat (tableau[0][0], ")|HIST(\"X\")");
	
	strcpy (tableau[0][1], "(");
	strcat (tableau[0][1], swap2);
	strcat (tableau[0][1], ")|HIST(\"Y\")");
	
	if (err = SrtGrfnMain (und, mdl_rows, paramStrings, valueStrings, 1, &obs_date,
		1, 2, tableau, mask, 0, 0, NULL, &temp_price, &temp_stdev, ss,NULL ))
	{	
		return err;
	}
	
	S1 = srt_f_gethistdata ("X");
	S2 = srt_f_gethistdata ("Y");
	
	for (i=1; i<=S1->num_path; i++)
	{
		e1 += S1->data[i];
		e2 += S2->data[i];
		var1 += (S1->data[i])*(S1->data[i]);
		var2 += (S2->data[i])*(S2->data[i]);
		covar += (S1->data[i])*(S2->data[i]);
	}

	e1 /= S1->num_path;
	e2 /= S1->num_path;
	var1 /= S1->num_path;
	var2 /= S1->num_path;
	covar /= S1->num_path;

	var1 -= e1 * e1;
	var2 -= e2 * e2;
	covar -= e1 * e2;
	
	*correl = covar / sqrt (var1 * var2);

	free_smatrix_size (tableau, 0, 0, 0, 1, GRFN_DEF_ARGBUFSZ);
	free_dmatrix (ss, 0, 0, 0, 1);
	
	return NULL;
}

/*

  Given any model, returns the normal correl between instr 1 and instr 2
  senn at some observation date through a grfn tableau.

*/

Err generic_correl (Date obs_date, 
					String instr1, 
					String instr2,  
					String und, 
					int mdl_rows, char** paramStrings, char **valueStrings,
					double *correl)
{

	char	***tableau = smatrix_size (0, 0, 0, 1, GRFN_DEF_ARGBUFSZ);
	double	**ss = dmatrix (0, 0, 0, 1),
			temp_price,
			temp_stdev;
	long	i;
	
	double	e1 = 0.0,
			e2 = 0.0,
			var1 = 0.0,
			var2 = 0.0,
			covar = 0.0;

	int		**mask = imatrix (0, 0, 0, 1); 
	
	SrtHistData		*S1 = NULL , *S2 = NULL;

	Err		err;

	mask[0][0] = mask[0][1] = GRFNSCELL; 

	strcpy (tableau[0][0], "(");
	strcat (tableau[0][0], instr1);
	strcat (tableau[0][0], ")|HIST(\"X\")");
	
	strcpy (tableau[0][1], "(");
	strcat (tableau[0][1], instr2);
	strcat (tableau[0][1], ")|HIST(\"Y\")");
	
	if (err = SrtGrfnMain (und, mdl_rows, paramStrings, valueStrings, 1, &obs_date,
		1, 2, tableau, mask, 0, 0, NULL, &temp_price, &temp_stdev,ss, NULL))
	{	
		return err;
	}
	
	S1 = srt_f_gethistdata ("X");
	S2 = srt_f_gethistdata ("Y");
	
	for (i=0; i<=S1->num_path-1; i++)
	{
		e1 += S1->data[i];
		e2 += S2->data[i];
		var1 += (S1->data[i])*(S1->data[i]);
		var2 += (S2->data[i])*(S2->data[i]);
		covar += (S1->data[i])*(S2->data[i]);
	}

	e1 /= S1->num_path;
	e2 /= S1->num_path;
	var1 /= S1->num_path;
	var2 /= S1->num_path;
	covar /= S1->num_path;

	var1 -= e1 * e1;
	var2 -= e2 * e2;
	covar -= e1 * e2;
	
	*correl = covar / sqrt (var1 * var2);

	free_smatrix_size (tableau, 0, 0, 0, 1, GRFN_DEF_ARGBUFSZ);
	free_dmatrix (ss, 0, 0, 0, 1);
	
	return NULL;
}
