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

#if defined(USENAG)
/*
 * Fortran wrapper routines
 */

extern	void	c_04ucf(int *N, double *blsc, double* busc,
			int *NCLIN, double *blclin, double *buclin,
				double *a1,
			int *NCNLN, double *blcnln, double *bucnln,
			int *ITER, double *C,
			double *OBJF, double *X, int *IFAIL);

	void	c_e04ucf_objfun(int *MODE, int *N, double *X,
			double *OBJF, double *OBJGRD);
	void	c_e04ucf_confun(int *NCNLN, int *N, int *MODE,
			int *NEEDC, double *X, double *C, double *cjac1);

#endif

static	TNLinProgConFunc	nagE04UCF_confun;
static	TNLinProgObjFunc	nagE04UCF_objfun;


/*f---------------------------------------------------------------------
 * C wrapper for the call to the nonlinear optmization routine
 * E04UCF of the NAG Fortran Library.
 */

DLL_EXPORT(int)
DrlNLinProgNagE04UCF(
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
static	char	routine[] = "DrlNLinProg";
	int	status = FAILURE;

#if defined(USENAG)
	/***********************************************************
	 * Interface for NAG routine E04UCF 
	 **********************************************************/
	double	*a1 ;
	int	iFail , 
		i, j ;

	nagE04UCF_objfun = obj ;
	nagE04UCF_confun = nlnc ;

	if ((a1 = DrlDoubleVectorAlloc((NVAR+1)*(NCLIN+1))) == NULL) {
		GtoErrMsg("%s: memory allocation failed (NVAR=%d, NCLIN=%d)\n",
			routine, NVAR, NCLIN);
		return(-4);
	}


	for (i=0; i<=NCLIN-1; i++)
	for (j=0; j<=NVAR-1; j++)
		a1[i*NVAR+j] = a[i][j] ;



	c_04ucf(&NVAR, blsc, busc,
			&NCLIN, blclin, buclin, a1,
			&NCNLN, blcnln, bucnln,
			&ITER, C, OBJF, X, &iFail) ;

	DrlDoubleVectorFree(a1) ;


	if (iFail != 0) {
		GtoErrMsg("DrlNLinProg: %s\n", DrlNLinProgErrMsg(iFail));
	}

	return(iFail) ;

#else
	GtoErrMsg("%s: routine not implemented\n", routine);
	return(-1);
#endif
}



#if defined(USENAG)

/*---------------------------------------------------------------------
 * Objective and Constraint Patch Functions for Nag
 */

void
c_e04ucf_objfun( int *MODE, int *N, double *X, double *OBJF, double *OBJGRD )
{
	nagE04UCF_objfun(MODE, N, X, OBJF, OBJGRD) ;
}



/*
 *	Constraint Function:
 *		cjac(i,j), 1<=i<=ncnln, 1<=j<=n
 *		si c_i(x_1,..,x_n) est la contrainte, cjac(i,j) = dc_i/dx_j
 *C	DO 400 I=1, NCNLN
 *C	DO 410 I=1, N
 *C		CJAC(I,J) = cjac1(N*(I-1)+J-1)
 */

void
c_e04ucf_confun(int *NCNLN, int *N, int *MODE,
		int *NEEDC, double *X, double *C, double *cjac1 )
{
	nagE04UCF_confun(NCNLN, N, MODE, NEEDC, X, C, cjac1) ;
}

#endif


/*----------------------------------------------------------------------
 * Error message printout
 */

DLL_EXPORT(char*)
DrlNLinProgNagE04UCFErrMsg(int errCode)
{
	switch (errCode) {
	case 0:
	return(" - Optimal solution found.");
	case 1:
	return(" - Optimal solution found, "
		" but requested accuracy not achieved.");
	case 2:
	return(" - No feasible point for the linear constraints.");
	case 3:
	return(" - No feasible point for the nonlinear constraints.");
	case 4:
	return(" - Too many major iterations.             ");
	case 5:
	return(" - Problem is unbounded (or badly scaled).");
	case 6:
	return(" - Current point cannot be improved upon. ");
	case 7:
	return(" - Large errors found in the derivatives. ");
	case 8:
	return(" - Error type 8.");
	case 9:
	return(" - Errors found in the input parameters. Problem abandoned.");
	default:
		if (errCode < 0) 
			return(" - user requested termination.");
		else
			return(" - Unknown error.");
		break;
	}

	return(NULL);	/* never reached */
}


