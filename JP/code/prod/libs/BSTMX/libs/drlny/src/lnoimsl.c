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
#ifdef	__DEBUG__
	int	_imslPrintLevel = 1;
#else
	int	_imslPrintLevel = 0;
#endif

#define	DRL_NLINPROG_MAXVAR	256		/* max # of variables */
#define	DRL_NLINPROG_MAXCONST	256		/* max # of constraints */

static	int	_nIter;
static	int	_ncnln,
		_nclin;
static	double	*_blclin,
		*_buclin,
		*_blcnln,
		*_bucnln,
		**_matclin;
static	void	*_userData;

static	void
_imsl_patch_objfun(
	int	m, 	/* (I) total number of constraints */
	int	meq,	/* (I) number of equality constraints */
	int	n,	/* (I) number of variables */
	double	*x,	/* (I) point where funct evaluated x[0..n-1] */
	int *active,	/* (I) array [0..max(1,m)-1] of active constraint */
	double	*f,	/* (O) value computed at point x */
	double	*g);	/* (I) array [0..max(1,m)-1] of values of constraint */

#endif

static	TNLinProgConFunc	confun;
static	TNLinProgObjFunc	objfun;

static	int		nlinprog_errno;



/*f---------------------------------------------------------------------
 * A wrapper for the call to the "imsl_d_min_con_nonlin" routine
 * of the IMSL library.
 */

DLL_EXPORT(int)
DrlNLinProgImslMinConNonlin(
	int	NVAR,		/* (I) number of dimensions */
	double	*blsc,		/* (I) lbounds state variables */
	double	*busc,		/* (I) ubounds state variables */
	int	NCLIN,		/* (I) number of lin constr */
	double	*blclin,	/* (I) lbounds on lin constr */
	double	*buclin,	/* (I) ubounds on lin constr */
	double	**a,		/* (I) coeff of lin constr */
	int	NCNLN,		/* (I) number of nlin constr */
	double	*blcnln,	/* (I) lbounds for nlin constr */
	double	*bucnln,	/* (I) ubounds for nlin constr */
	TNLinProgConFunc nlnc,	/* (I) C-routine for nlin constr */
	TNLinProgObjFunc obj,	/* (I) C-routine for obj */
	int	ITER,		/* (I) maximum number of iterations */
	double	*C,		/* (I) initial values of nlin constr
				 * (O) final values of nlin constr */
	double	*OBJF,		/* (I) initial value of obj
				 * (O) final value of obj */
	double	*X,		/* (I) starting point
				 * (O) optimal point */
	void	*userData)	/* (I) user data */
{
static	char	routine[] = "DrlNLinProgImslMinConNonlin";
	int	status = FAILURE;

#if defined(USEIMSL) 
	/****************************************************************
	 * Interface for C-Math Imsl Library
	 ****************************************************************/

	int	i,
		m, 		/* total number of constraints */
		meq,		/* number of equality constraints */
		n,		/* number of variables */
		ibtype;		/* scalar indicating type of bounds on var */
	long	imslErrCode;
	double	*x0 = NULL,	/* initial guess point for the optim */
		*xSol = NULL,	/* final point */
		*xlb = NULL,	/* array of low bounds on variable */
		*xub = NULL;	/* array of high bounds on variable */ 

	/*
	 *
	 */
	nlinprog_errno = 0;

	n = NVAR;		/* number of variables */
	m = 2*(NCLIN + NCNLN);	/* total number of >= constraints */
	meq = 0;		/* no equality constraint */
	ibtype = 0;		/* user supplies all the bounds */

	_nclin = NCLIN;
	_ncnln = NCNLN;
	_nIter=0;

	_blclin = blclin;
	_buclin = buclin;
	_matclin = a;

	_blcnln = blcnln;
	_bucnln = bucnln;

	objfun = obj ;
	confun = nlnc ;

	_userData = userData;

	/*
	 * Parameter checking
	 */
	if ((NVAR <= 0) || (NVAR > DRL_NLINPROG_MAXVAR)) {
		GtoErrMsg("%s: [imsl] too large or too low NVAR (%d).\n",
			routine, NVAR);
		goto done;
	}

	if ((ITER <= 0) || (ITER > 40000)) {
		GtoErrMsg("%s: [imsl] too large or too low ITER (%d).\n",
			routine, ITER);
		goto done;
	}

	if ((NCLIN < 0) || (NCLIN > DRL_NLINPROG_MAXCONST)) {
		GtoErrMsg("%s: [imsl] too large or too low NCLIN (%d).\n",
			routine, NCLIN);
		goto done;
	}

	if ((NCNLN < 0) || (NCNLN > DRL_NLINPROG_MAXCONST)) {
		GtoErrMsg("%s: [imsl] too large or too low NCNLN (%d).\n",
			routine, NCNLN);
		goto done;
	}


	/*
	 *
	 */

	if ((x0 = DrlDoubleVectAlloc(0, n-1)) == NULL) 
		goto done;
	for (i=0; i<=n-1; i++) x0[i] = X[i];

	if ((xlb = DrlDoubleVectAlloc(0, n-1)) == NULL) 
		goto done;

	if ((xub = DrlDoubleVectAlloc(0, n)) == NULL) 
		goto done;

	/* Variable bounds constraints */
	for (i=0; i<=n-1; i++) {
		xlb[i] = blsc[i];
		xub[i] = busc[i];
	}

	/*
	 *
	 */
	if (_imslPrintLevel > 0) {
		int	i1, j1;
		fprintf(stdout, "%s:\n", routine);
		fprintf(stdout, "NVAR = %d\n", NVAR);
		for (i1=0; i1<=n-1; i1++) {
			fprintf(stdout,
				"\t[%2d]\txLow=%lf\txHigh=%lf\tx0=%lf\n", 
				i1, xlb[i1], xub[i1], x0[i1]);
		}
		fprintf(stdout, "Linear constraints (NCLIN=%d)\n", NCLIN);
		fprintf(stdout, "Low bound:\n");
		DrlFPrintDoubleVect(stdout, NULL, _blclin, NCLIN);
		fprintf(stdout, "High bound:\n");
		DrlFPrintDoubleVect(stdout, NULL, _buclin, NCLIN);
		fprintf(stdout, "Matrix:\n");
		DrlFPrintDoubleMatr(stdout, NULL, _matclin, NCLIN, NVAR);

		fprintf(stdout, "Nonlinear constraints (NCNLN=%d)\n", NCNLN);
		fprintf(stdout, "Low bound:\n");
		DrlFPrintDoubleVect(stdout, NULL, _blcnln, NCNLN);
		fprintf(stdout, "High bound:\n");
		DrlFPrintDoubleVect(stdout, NULL, _bucnln, NCNLN);

		fprintf(stdout, "Total positive constraints (m=%d)\n", m);
		fprintf(stdout, "Total equality constraints (meq=%d)\n", meq);
	}

	/*
	 * IMSL error setting
	 */
	imsl_error_options(
		IMSL_SET_STOP, IMSL_FATAL, 0,
		IMSL_SET_STOP, IMSL_TERMINAL, 0,
		IMSL_SET_ERROR_FILE, stdout,
		0);


	/*
	 * Call to IMSL routine
	 */
	xSol = imsl_d_min_con_nonlin(
			_imsl_patch_objfun,
			m, meq, n, ibtype, xlb, xub,
			IMSL_XGUESS, x0,
			IMSL_PRINT, _imslPrintLevel,
			IMSL_ITMAX, ITER,
			/*IMSL_RETURN_USER, X,*/
			IMSL_OBJ, OBJF,
			0);
	imslErrCode = imsl_error_code();
	if (_imslPrintLevel > 0) {
		fprintf(stdout, "%s: imsl err code = %d\n",
			routine, imslErrCode);
	}

	if ((xSol == NULL) || (imslErrCode != 0L)) {
		if (_imslPrintLevel > 0) {
		    fprintf(stdout, "IMSL ERROR code = %d\n", imslErrCode);
		}
		GtoErrMsg("%s: optimization failed (code %d)\n",
			routine, imslErrCode);
		GtoErrMsg("%s: %s\n", routine,
			DrlNLinProgImslMinConNonlinErrMsg((int) imslErrCode));
		if (xSol != NULL) {
			for (i=0; i<=n-1; i++) X[i] = xSol[i];
		} else {
			for (i=0; i<=n-1; i++) X[i] = 0.0;
		}
		goto done;
	} else {
		for (i=0; i<=n-1; i++) X[i] = xSol[i];
	}

	/* check for error in user defined functions */
	if (nlinprog_errno != 0) {
	    GtoErrMsg("%s: error (code %d) occured in callback.\n",
		routine, nlinprog_errno);
	    goto done;
	}

	/*
	 * End
	 */
	status = SUCCESS;
done:
	if (xSol != NULL) free((char*) xSol);
	DrlDoubleVectFree(x0,  0, n-1);
	DrlDoubleVectFree(xlb, 0, n-1);
	DrlDoubleVectFree(xub, 0, n);

	fflush(stdout);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);


#else
	GtoErrMsg("%s: routine not implemented\n", routine);
	return(-1);
#endif
}



#if defined(USEIMSL) 
/*----------------------------------------------------------------------
 * Objective and Constraint Patch Functions for IMSL
 * There is a total of m = 2*(NCLIN+NCNLN) positive constraints:
 *	m=0		1st lin LOW
 *	m=1		1st lin HIGH
 *	m=2		2nd lin LOW
 *	etc..
 *	m=2*NCNLN	1st nlin LOW
 *	m=2*NCNLN+1	1st nlin HIGH
 *	rtc...
 */

static	void
_imsl_patch_objfun(
	int	m, 		/* (I) total number of constraints */
	int	meq,		/* (I) number of equality constraints */
	int	n,		/* (I) number of variables */
	double	*x,		/* (I) point where funct evaluated x[0..n-1] */
	int	*active,	/* (I) array [0..max(1,m)-1] active constr */
	double	*f,		/* (O) value computed at point x */
	double	*g)		/* (I) array [0..max(1,m)-1] constr  values */
{
static	char	routine[] = "_imsl_patch_objfun";
	double	val;
	int	mode = 0,
		i, j;
	int	errCode;

static	double	tmpConstr[DRL_NLINPROG_MAXCONST+1];
static	int	tmpNeedc[DRL_NLINPROG_MAXCONST+1];

	/**
	 ** Objective function
	 **/
	errCode = objfun(&mode, &n, x, _userData, f, NULL);
	nlinprog_errno |= errCode;
	if (errCode != 0) {
	    GtoErrMsg("%s: objective function failed (warning).\n", routine);
	}


	/**
	 ** Linear constraints
	 **/
	for (i=0; i<=_nclin-1; i++) {
	   /* check if constraint active */
	   if ((active[2*i]) || (active[2*i+1])) {
		val = 0.;
		for (j=0; j<=n-1; j++)
			val += _matclin[i][j] * x[j];

	  	if (active[2*i  ]) g[2*i  ] =   val - _blclin[i];
	  	if (active[2*i+1]) g[2*i+1] = - val + _buclin[i];
	   }
	}

	/**
	 ** Nonlinear constraints
	 **/
	for (i=0; i<=_ncnln-1; i++) {
		tmpNeedc[i] = ((active[2*(_nclin+i)]) ||
			       (active[2*(_nclin+i)+1]));
	}

	/* call nonlinear constraints function */
	errCode = confun(&_ncnln, &n, &mode, tmpNeedc,
				x, _userData, tmpConstr, (double*) NULL);
	nlinprog_errno |= errCode;
	if (errCode != 0) {
	    GtoErrMsg("%s: constraint function failed (warning).\n", routine);
	}




	if (_imslPrintLevel >= 2) {
		fprintf(stdout, "NCNLN: ");
		for (i=0; i<=_ncnln-1; i++)
			fprintf(stdout, "\t%.3g", tmpConstr[i]);
		fprintf(stdout, "\n");

		fprintf(stdout, "BLCNLN:");
		for (i=0; i<=_ncnln-1; i++)
			fprintf(stdout, "\t%.3g", _blcnln[i]);
		fprintf(stdout, "\n");

		fprintf(stdout, "BUCNLN:");
		for (i=0; i<=_ncnln-1; i++)
			fprintf(stdout, "\t%.3g", _bucnln[i]);
		fprintf(stdout, "\n");

	}

	for (i=0; i<=_ncnln-1; i++) {
	  if ((active[2*(_nclin+i)]) || (active[2*(_nclin+i)+1])) {
	  	if (active[2*(_nclin+i)  ] != 0)
			g[2*(_nclin+i)  ] =   tmpConstr[i] - _blcnln[i];
	  	if (active[2*(_nclin+i)+1] != 0)
			g[2*(_nclin+i)+1] = - tmpConstr[i] + _bucnln[i];
	  }
	}

	if (_imslPrintLevel >= 2) {
		fprintf(stdout, "<%3d>", ++_nIter);
		for (i=0; i<=m-1; i++) fprintf(stdout, "\t%d", i);
		fprintf(stdout, "\n");

		fprintf(stdout, "Active:");
		for (i=0; i<=m-1; i++) fprintf(stdout, "\t%d", active[i]);
		fprintf(stdout, "\n");

		fprintf(stdout, "Value:");
		for (i=0; i<=m-1; i++) fprintf(stdout, "\t%.3g", g[i]);
		fprintf(stdout, "\n");
	}


	/*
	 *
	 */
	return;
}


#endif


/*----------------------------------------------------------------------
 * Error message printout
 */

DLL_EXPORT(char*)
DrlNLinProgImslMinConNonlinErrMsg(int errCode)
{

#if defined(USEIMSL) 

static	char	buf[256];

	switch (errCode) {
	case 0:
		return(" - Optimal solution found.");
	case IMSL_TOO_MANY_ITN:
		return(" - [imsl] too many iterations");
	case IMSL_UPHILL_DIRECTION:
		return(" - [imsl] uphill search direction");
	case IMSL_TOO_MANY_LINESEARCH:
		return(" - [imsl] line search took more that "
		"5 function calls (probably unconsistent constraints)");
	case IMSL_NO_PROGRESS_MADE:
		return(" - [imsl] search direction close to zero");
	case IMSL_QP_INCONSISTENT:
		return(" - [imsl] the constraints for the QP "
		"problem are inconsistent");
	default:
		sprintf(buf, " - [imsl] unknown error code %d", errCode);
		return(buf);
	}
#else
	return(" - ERROR: unknown platform ");
#endif

	return(NULL);	/* never reached */
}









