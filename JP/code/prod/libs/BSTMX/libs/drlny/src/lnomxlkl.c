/************************************************************************
 * Module:	DRL
 * Submodule:	LNO - Linear and Nonlinear Optimization
 * File:	
 * Function:	
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "drlmem.h"
#include "drlstr.h"
#include "drlio.h"
#include "drlgpars.h"		/* GParEval() */
#include "drlvtype.h"

#include "drllno.h"

#undef	__DEBUG__
#define	__DEBUG__

#define	MAXVAR	10


static	int	ndata,
		nvar,
 		ncnln;
static	double	xvec[MAXVAR+1],
		*Xi,
		*Yi,
		*blsc,
		*busc,
		*blcnln,
		*bucnln;
static	char	*fopt;

static	int	confun(int *NCNLN, int *N, int *MODE,
			int *NEEDC, double *X,
			void *userData, double *C, double *cjac1);
static	int	objfun(int *MODE, int *N, double *X,
			void *userData, double *OBJF, double *OBJGRD);


/*f---------------------------------------------------------------------
 * Performs a maximum likehood estimation of the form
 * <br>
 *   \min_{z_1,\dots,z_n}
 *          \Bigl(y_i - \phi(x_i,z_1,\dots,z_n)\Bigr)
 * <br>
 */

DLL_EXPORT(int)
DrlMaxLikehoodL(
	long *ndataL,		/* 'L'  1 (I) # data */
	double *XiL,		/* 'F'  2 (I) array of X data [0..nData-1] */
	double *YiL,		/* 'F'  3 (I) array of Y data [0..nData-1] */
	long *nvarL,		/* 'L'  4 (I) # variables */
#ifdef	_SKIP
	long *optFlagL,		/* 'L'  5 (I) starting point [0..nVar-1] */
#endif
	double *zstartL,	/* 'F'  6 (I) starting point [0..nVar-1] */
	double *blscL,		/* 'F'  7 (I) low  bound [0..nVar-1] */
	double *buscL,		/* 'F'  8 (I) high bound [0..nVar-1] */
 	long *ncnlnL,		/* 'L'  9 (I) # non lin constrains */
	double *blcnlnL,	/* 'F'  0 (I) low  bound nonlin constr */
	double *bucnlnL,	/* 'F' 11 (I) high bound nonlin constr */
	char *fcnlnL,		/* 'C' 12 (I) nonlin constraints */
	long *itermaxL,		/* 'L' 13 (I) max # of iterations */
	char *foptL,		/* 'C' 14 (I) functional form */
	double *valoptL,	/* 'F' 15 (O) optimal max lik [0] */
	double *zoptL)		/* 'F' 16 (O) optimal point [0..nVar-1] */
{
static	char	routine[] = "DrlMaxLikehoodL";
	int	status = FAILURE, errCode;
	int	i;
	double	Z[MAXVAR];

#undef	CHECK
#define	CHECK(cond)	{if (!(cond)) {GtoErrMsg("%s: assertion failed (%s)\n",\
			routine, #cond); goto done;}}


        /* Set log file */
#ifdef	__DEBUG__
        if (DrlStdoutFileSet("stdout.log", "w") < 0) return(FAILURE);
	fprintf(stdout, "%s:\n", routine);
#endif


	/* check arguments */
	CHECK(ARGSIZE(ndataL) == 1);
	CHECK(ARGSIZE(XiL) >= ARGSIZE(ndataL));
	CHECK(ARGSIZE(YiL) >= ARGSIZE(ndataL));

	ndata = (int) ndataL[1];
	Xi = XiL+1;
	Yi = YiL+1;

	CHECK(ARGSIZE(nvarL) == 1);
	CHECK(ARGSIZE(nvarL) <= MAXVAR);
	/*CHECK(ARGSIZE(optFlagL) >= ARGSIZE(nvarL));*/
	CHECK(ARGSIZE(zstartL) >= ARGSIZE(nvarL));
	CHECK(ARGSIZE(zoptL) >= ARGSIZE(nvarL));
	CHECK(ARGSIZE(blscL) >= ARGSIZE(nvarL));
	CHECK(ARGSIZE(buscL) >= ARGSIZE(nvarL));

	nvar = (int) nvarL[1];
	blsc = blscL+1;
	busc = buscL+1;

	CHECK(ARGSIZE(valoptL) == 1);
	CHECK(ARGSIZE(itermaxL) == 1);

	CHECK(ARGSIZE(foptL) == 1);
	fopt = foptL+1;

	for(i=0; i<=nvar-1;i++)
		Z[i] = zstartL[i+1];


	/*
	 * Check syntax
	 */
	xvec[0] = 1.234567;
	for(i=0; i<=nvar-1;i++) xvec[i+1] = Z[i];
	if (DrlGParEval(
		nvar+1,
		xvec,
		fopt,
		valoptL+1) != 0) {
		    GtoErrMsg("%s: parsing failed\n", routine);
		    goto done;
	}



	/* do the optimization */
	errCode = DrlNLinProg(
		nvar, blsc, busc,
		0, NULL, NULL, NULL, 	/* nclin */
 		/*ncnln, blcnln, bucnln,*/
 		0, NULL, NULL,
		confun,
		objfun,
		(int) itermaxL[1],
		NULL,	/* initial value of cnln */
		valoptL+1,
		Z,
		NULL,
		DRL_LNO_IMSL_NLN,
		0);

	if (errCode != 0) {
		GtoErrMsg("%s: optimization failed.\n", routine);
		goto done;
	} else {
		for(i=0; i<=nvar-1;i++)
			zoptL[i+1] = Z[i];
	}


	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}



/*ARGSUSED*/
static	int
objfun(
	int *MODE,		/* (I) NOT USED */
	int *N,			/* (I) number of variables */
	double *X,		/* (I) point where obective valued */
	void *userData,		/* (I) user data */
	double *OBJF,		/* (O) value of the objective */
	double *OBJGRD)		/* (O) NOT USED */
{
	double	retVal;
	int	i, k;

#ifdef	__DEBUG__
		fprintf(stdout, "objfun:\n");
		for(i=0; i<=nvar-1;i++)
			fprintf(stdout, "%lf\t", X[i]);
		fprintf(stdout, "\n");
#endif
	*OBJF = 0e0;
	for (k=0; k<=ndata-1; k++) {
		xvec[0] = Xi[k];
		for(i=0; i<=nvar-1;i++) xvec[i+1] = X[i];

		if (DrlGParEval(
			nvar+1,
			xvec,
			fopt,
			&retVal) != 0) {
		    GtoErrMsg("objfun: parsing failed\n");
		    *OBJF = 0e0;
		    return(FAILURE);
		}
#ifdef	__DEBUG__
		fprintf(stdout, "data[%3d] val=%lf\tYi=%lf\n",
			k, retVal, Yi[k]);
#endif
		*OBJF += (retVal - Yi[k])*(retVal-Yi[k]);
	}
#ifdef	__DEBUG__
	fprintf(stdout, "OBJF=%lf\n", *OBJF);
#endif
	return(SUCCESS);
}


/*ARGSUSED*/
static	int
confun(
	int *NCNLN,		/* (I) number of nonlinear constraints */
	int *N,			/* (I) number of variables */
	int *MODE,		/* (I) NOT USED */
	int *NEEDC,		/* (I) array indicate if constraint needed */
	double *X,		/* (I) point wher constraints valued */
	void *userData,		/* (I) user data */
	double *C,		/* (O) array of values of the constraints */
	double *cjac1)		/* (O) NOT USED */
{
	return(SUCCESS);
}


