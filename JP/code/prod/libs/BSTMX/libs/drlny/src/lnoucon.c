/************************************************************************
 * Module:	DRL
 * Submodule:	LNO - Linear and Nonlinear Optimization
 * File:	
 * Function:	Optimisation Routines
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>

#include "drlio.h"
#include "drlmem.h"

#include "drllno.h"		/* prototype consistency */


#if defined(USEIMSL) 
#include <imsl.h>

static	double	_DrlLnoUconMin_imsl_f(
	int n,			/* (I) number of variables */
	double *x);		/* (I) point where funct evaluated x[0..n-1] */

static	const	char* _DrlLnoUconMin_imsl_code(int errCode);
#endif


#define	MAXVAR	256

/* Static variable */
static	TDrlLnoObjFunc		objfunS;
static	void			*userDataS;
static	int			isErrorS;
static	int			nxS;
static	double			xS[MAXVAR];
static	int			k2iS[MAXVAR];


/*f---------------------------------------------------------------------
 * A generic nonlinear unconstrained optimization routine.
 * 
 * The argument optimMethod can have teh following values:\\
 * "IMSL" for imsl.\\
 * WARNING: SINCE THE ROUTINE USES STATIC VARIABLES,
 * IT CANNOT BE CALLED RECURSIVELY.
 *
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlLnoUnconMin(
	int nx,			/* (I) number of dimensions */
	double *x,		/* (B) starting/optimal point [0..nx-1] */
	int *optx,		/* (I) optimize (true/false) [0..nx-1] */
	double *f,		/* (O) final value of objective */
	int havedf,		/* (I) gradient provided (TRUE, FALSE) */
	TDrlLnoObjFunc objfun,	/* (I) objective function (see type definition) */
	void *userData,		/* (I) user data (for callback) */
	char *optimMethod,	/* (I) optimization method */
	TLNOParams *params) 	/* (I) optimization parameters (or NULL=def) */
{
static	char	routine[] = "DrlLnoUconMin";
	int	status = FAILURE;
	int	ny,		/* number of optimized variables */
		i,		/* index in input variable [0..nx-1] */
		k;		/* index in optimization [0..ny-1] */
	int	itermax = 100;

#if defined(USEIMSL) 
	int		errCodeImsl;
	double		*y0 = NULL,	/* initial guess point for the optim */
			*ySol = NULL;	/* final point */
#endif

	TLNOParams	defParams;	/* default parameters */


	/* Set default parameters */
	IF_FAILED_DONE(TLNOParamsSetDefault(&defParams));
	if (params == NULL) params = &defParams;

	/*
	 * Map indices if only subset of variable optimized:
	 * The k-th  optmization variable corresponds
	 * to the the i-th input variable
	 */
	if (optx) {
		ny = 0;
		k = 0;
		for (i=0; i<nx; i++) {
			if (optx[i] != 0) {
				k2iS[k] = i;
				/* one more opt variable */
				k++;
			}
		}
		ny = k;
	} else {
		/* Optimize everything */
		ny = nx;
		for (k=0; k<nx; k++) {
			k2iS[k] = k;
		}
	}

	if (ny == 0) {
		*f = 0e0;
		status = SUCCESS;
		goto done;
	}


#ifdef	DEBUG
	for (i=0; i<nx; i++)
		printf("%2d/%2d  %lf  %2d\n", i, nx, x[i], optx[i]);

	for (k=0; k<ny; k++)
		printf("%2d %2d  k:%2d  k2i:%2d \n", nx, ny, k, k2iS[k]);
#endif

	/*
	 * Set static variables used in callbacks
	 */
	objfunS = objfun;
	userDataS = userData;
	isErrorS = FALSE;
	nxS = nx;
	for (i=0; i<nx; i++) xS[i] = x[i];

	/* Initial solution */
	if ((y0 = DrlDoubleVectAlloc(0, ny-1)) == NULL) 
		goto done;
	for (k=0; k<ny; k++) y0[k] = x[k2iS[k]];

#ifdef	DEBUG
	printf("%s: initial\n", routine);
	for (k=0; k<ny; k++)
		printf("Y0[%2d] %lf\n", k, y0[k]);
#endif



#if defined(USEIMSL) 
	if (optimMethod[0] == 'I') {
		/*==============================================
		 * USING IMSL 
		 *==============================================*/

		/* IMSL error setting */
		imsl_error_options(
			IMSL_SET_STOP, IMSL_FATAL, 0,
			IMSL_SET_STOP, IMSL_TERMINAL, 0,
			IMSL_SET_ERROR_FILE, stdout,
			0);


		/* Call to IMSL routine */
		ySol = imsl_d_min_uncon_multivar(
			_DrlLnoUconMin_imsl_f,
			ny,

			IMSL_XGUESS, y0,
			IMSL_MAX_FCN, (int) itermax,
			IMSL_MAX_ITN, (int) itermax,
			IMSL_FVALUE, (double*) f,
			0);
		errCodeImsl = imsl_error_code();
	
		if ((ySol == NULL) || (errCodeImsl != 0L)) {
			GtoErrMsg("%s: optimization failed (code %d, %s)\n",
				routine, errCodeImsl,
				_DrlLnoUconMin_imsl_code(errCodeImsl));
			if (ySol != NULL) {
				for (k=0; k<ny; k++) x[k2iS[k]] = ySol[k];
			}
			goto done;
		} else {
			for (k=0; k<ny; k++) x[k2iS[k]] = ySol[k];
		}
#endif	/* USEIMSL */
	} else {
	    GtoErrMsg("%s: unavailable method `%s'.\n", routine, optimMethod);
	}


	/* check for error in user defined functions */
	if (isErrorS != 0) {
	    GtoErrMsg("%s: error occured in callback.\n", routine);
	    goto done;
	}

	/*
	 * End
	 */
	status = SUCCESS;
done:
	if (ySol != NULL) free((char*) ySol);
	DrlDoubleVectFree(y0,  0, ny-1);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}





#if defined(USEIMSL) 

/*----------------------------------------------------------------------
 * Callback for imsl libraries
 */

static	double
_DrlLnoUconMin_imsl_f(
	int ny,			/* (I) number of variables */
	double *y)		/* (I) point where funct eval y[0..ny-1] */
{
static	char	routine[] = "_DrlLnoUconMin_imsl_f";
	int	errCode;
static	int	tmpNeeddf[MAXVAR];
	int	i, k;
	double	fvalue;

	/**
	 ** Objective function
	 **/
	/* Set optimizaed variable */
	for (k=0; k<ny; k++) xS[k2iS[k]] = y[k];

	for (i=0; i<nxS; i++) tmpNeeddf[i] = 0;

	errCode = objfunS(
		nxS,
		xS,
		1, 	/* 0=unused, 1=needed, 2=provided */
		&fvalue,
		tmpNeeddf,
		NULL,
		userDataS);/* (I) user supplied data */


#ifdef	DEBUG
	printf("%s: ", routine);
	for (i=0; i<nxS; i++) printf(" %lf ", xS[i]);
	printf("\t  FVALUE %lf\n", fvalue);
#endif


	if (errCode != 0) {
		isErrorS = TRUE;
		GtoErrMsg("%s: objective function failed (warning).\n",
				routine);
	}

	return (fvalue);
}


static	const	char*
_DrlLnoUconMin_imsl_code(int errCode)
{
#undef	CASE
#define	CASE(code)	case code: return(#code);
	switch (errCode) {
	CASE(IMSL_STEP_TOLERANCE)

	CASE(IMSL_REL_FCN_TOLERANCE)
	CASE(IMSL_TOO_MANY_ITN)
	CASE(IMSL_TOO_MANY_FCN_EVAL)
	CASE(IMSL_TOO_MANY_GRAD_EVAL)
	CASE(IMSL_UNBOUNDED)
	CASE(IMSL_NO_FURTHER_PROGRESS)

	CASE(IMSL_FALSE_CONVERGENCE)
	}
	return ("(unknown)");
}

#endif	/* USEIMSL */





