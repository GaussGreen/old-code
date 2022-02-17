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



typedef	struct	{
	TNLinProgObjFunc	objFunc;
	TNLinProgConFunc	conFunc;
	void			*userData;
} TLnoObjConFuncData;

static	int	drlLnoObjConFunc(
	int NVAR,		/* (I) number of variables */
	double *X,		/* (I) point where obective valued [NVAR] */
	int NEEDO,		/* (I) indicates if objective  needed */
	double *OBJF,		/* (O) value of objective */
	int NCNLN,		/* (I) number of nlin constr */
	int *NEEDC,		/* (I) array indicate if const needed [NCNLN] */
	double *CNLN,		/* (O) array of values of constr [NCNLN] */
	void *userData);	/* (I) user data */


/*f---------------------------------------------------------------------
 * A generic nonlinear optimization routine with
 * linear and non-linear constraints and no derivatives provided. \\
 * It performs a minimization of a function
 * 	objective(x0,...,xNVAR-1)
 * of NVAR variables under \\
 * 1. NVAR state variable constraints \\
 *      blsc(i) <=  xi <= busc(i) for i=0,..,NVAR-1 \\
 * 2. NCLIN linear constraints on the state variables of the form \\
 *      blclin(i) <=  sum a(i,j) xj <= buclin(i) for i=0,..,NCLIN-1 \\
 * 3. NCNLN nonlinear constraints of the form
 *      blcnln(i) <= constr(i,x0,...,xNVAR-1) xj <= bucnln(i)
 *                                              for i=0,..,NCNLN-1 \\
 * The optimization method currently available are: \\
 * DRL_LNO_IMSL_NLN: uses the IMSL CMATH routine min_con_nonlin
 *     (when the library is complied with the USEIMSL flag
 *     and linked with the IMSL CMATH library).\\
 * DRL_LNO_NAG: uses the NAG routine EC04UCF (when the library
 *     is complied with the USENAG flag and linked with the NAG Fortran
 *     library).\\
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlNLinProg(
	int NVAR,		/* (I) number of dimensions */
	double *blsc,		/* (I) lbounds state variables */
	double *busc,		/* (I) ubounds state variables */
	int NCLIN,		/* (I) number of lin constr */
	double *blclin,		/* (I) lower bounds on lin constr */
	double *buclin,		/* (I) upper bounds on lin constr */
	double **a,		/* (I) matrix coeff of lin constr */
	int NCNLN,		/* (I) number of nlin constr */
	double *blcnln,		/* (I) lower bounds for nlin constr */
	double *bucnln,		/* (I) upper bounds for nlin constr */
	TNLinProgConFunc nlnc,	/* (I) C-routine for nlin constr */
	TNLinProgObjFunc obj,	/* (I) C-routine for obj */
	int ITER,		/* (I) maximum number of iterations */
	double *C,		/* (I) initial values of nlin constr
				 * (O) final values of nlin constr */
	double *OBJF,		/* (I) initial value of objective
				 * (O) final value of objective */
	double *X,		/* (I) starting point
				 * (O) optimal point */
	void *userData,		/* (I) user data */
	long optimMethod,	/* (I) optimization method */
	TLNOParams *params) 	/* (I) optimization parameters (or NULL=def) */
{
static	char	routine[] = "DrlNLinProg";
	int	status = FAILURE;
	TLnoObjConFuncData	data;
	TLNOParams	defParams;	/* default parameters */


	/* Set default parameters */
	IF_FAILED_DONE(TLNOParamsSetDefault(&defParams));
	if (params == NULL) params = &defParams;

	/*
	 *
	 */
	switch (optimMethod) {
	case DRL_LNO_IMSL_NLN:
		/* Imsl min_non_conlin routine
		 */

		IF_FAILED_DONE( DrlNLinProgImslMinConNonlin(
			NVAR,
			blsc,
			busc,
			NCLIN,
			blclin,
			buclin,
			a,
			NCNLN,
			blcnln,
			bucnln,
			nlnc,
			obj,
			ITER,
			C,
			OBJF,
			X,
			userData));
		break;

#ifdef	_SKIP
		IF_FAILED_DONE( DrlNLinOptimQMC(
			NVAR,
			blsc,
			busc,
			NCLIN,
			blclin,
			buclin,
			a,
			NCNLN,
			blcnln,
			bucnln,
			nlnc,
			obj,
			1000, /*ITERMAX,*/
			C,
			OBJF,
			X,
			userData,
			NULL));

	IF_FAILED_DONE( DrlNLinProg(
			NVAR,
			blsc,
			busc,
			NCLIN,
			blclin,
			buclin,
			a,
			NCNLN,
			blcnln,
			bucnln,
			nlnc,
			obj,
			ITER,
			C,
			OBJF,
			X,
			userData,
			NULL));

#endif


	case DRL_LNO_COMPBOX:
		if (NCLIN != 0) {
			GtoErrMsg("%s: COMPBOX does not support linear"
				" constraints.\n", routine);
			goto done;
		}
		data.objFunc = obj;
		data.conFunc = nlnc;
		data.userData = userData;


		IF_FAILED_DONE(	DrlNLinProgComplexBox(
			NVAR,
			blsc,
			busc,
			NCNLN,
			blcnln,
			bucnln,
			drlLnoObjConFunc,
			ITER,
			C,
			OBJF,
			X,
			(void*) &data,
			params));


		break;
	default:
		GtoErrMsg("%s: unknown optimization method %ld.\n",
			routine, optimMethod);
		goto done;
	}


	/*GtoErrMsg("%s: OPTIMAL VALUE %e\n", routine, *OBJF);
	DrlFPrintDoubleVect(NULL, "\t%g", X, NVAR);*/


	/*
	 * End
	 */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}




/*f---------------------------------------------------------------------
 * Patch routine for complex box method.
 */


static	int
drlLnoObjConFunc(
	int NVAR,		/* (I) number of variables */
	double *X,		/* (I) point where obective valued [NVAR] */
	int NEEDO,		/* (I) indicates if objective  needed */
	double *OBJF,		/* (O) value of objective */
	int NCNLN,		/* (I) number of nlin constr */
	int *NEEDC,		/* (I) array indicate if const needed [NCNLN] */
	double *CNLN,		/* (O) array of values of constr [NCNLN] */
	void *userData)		/* (I) user data */
{
	TLnoObjConFuncData	*data =  (TLnoObjConFuncData*) userData;
	int			retCode;

	retCode = data->objFunc(
		NULL,
		&NVAR,
		X,
		data->userData,
		OBJF,
		NULL);
	if (retCode != SUCCESS) return(retCode);

	retCode = data->conFunc(
		&NCNLN,
		&NVAR,
		NULL,
		NEEDC,
		X,
		data->userData,
		CNLN,
		NULL);
	if (retCode != SUCCESS) return(retCode);


	return(SUCCESS);
}





/*f---------------------------------------------------------------------
 * Sets the default values for the optimization parameters.
 */

int
TLNOParamsSetDefault(TLNOParams *params)
{
	params->optimMethod = DRL_LNO_IMSL_NLN;

	params->xtol = DBL_EPSILON;
	params->ftol = DBL_EPSILON;
	params->btol = DBL_EPSILON;
	params->restart = 0;

	return(SUCCESS);
}


