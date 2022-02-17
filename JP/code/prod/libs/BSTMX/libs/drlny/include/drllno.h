/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	LNO - Linear and Nonlinear Optimization
 * File:	drllno.h
 * Function:	Optimization routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drllno_H
#define	_drllno_H

#include "drlstd.h"

/*==============================================================*
 * Prototype and optimization parameters definitions.
 *==============================================================*/

/*
 * Optimization types
 */

#define	DRL_LNO_NONE		(0L)
#define	DRL_LNO_IMSL		(14L)
#define	DRL_LNO_IMSL_NLN	(10L)
#define	DRL_LNO_QMC		(12L)
#define	DRL_LNO_COMPBOX		(20L)

/*t-@CDOC(idxn=TLNOParams,catn=structdef)
 * Structure for optimization parameters.
 */
typedef	struct {
	long	optimMethod;	/* optimization method */
	double	xtol;		/* tolerance on variable */
	double	ftol;		/* tolerance on objective function */
	double	btol;		/* tolerance on constraints */
	int	restart;	/* number of restarts (default 0) */
} TLNOParams;
/*e*/

extern	int	TLNOParamsSetDefault(TLNOParams *params);


/*==============================================================*
 *
 * UNCONSTRAINED NON-LINEAR OPTIMIZATION
 * 
 *==============================================================*/

/*t-@CDOC(idxn=TDrlLnoObjFunc,catn=structdef)
 * Prototype for a user-suplied  objective function for an unconstrained
 * nonlinear optimization problem.
 */
typedef	int	(*TDrlLnoObjFunc)(
	int nx,		/* (I) number of variables */
	double *x,	/* (I) point where obective valued [nx] */
	int needf,	/* (I) need objective: 0=no, 1=yes, 2=given */
	double *f,	/* (B) value of objective */
	int *needdf,	/* (I) need grad: 0=no, 1=yes, 2=given [0..nx-1] */
	double *df,	/* (B) gradient of objective [0..nx-1] */
	void *userData);/* (I) user supplied data */


/*
 * General unconstrained optimization routine
 */
DLL_EXPORT(int)	DrlLnoUnconMin(
	int nx,			/* (I) number of dimensions */
	double *x,		/* (B) starting/optimal point */
	int *optx,		/* (I) optimize (true/false) [0..nx-1] */
	double *f,		/* (O) final value of objective */
	int havedf,		/* (I) gradient provided (TRUE, FALSE) */
	TDrlLnoObjFunc objfun,	/* (I) objective function */
	void *userData,		/* (I) user data (for callback) */
	char *optimMethod,	/* (I) optimization method */
	TLNOParams *params); 	/* (I) optimization parameters (or NULL=def) */


/*==============================================================*
 *
 * CONSTRAINED NON-LINEAR OPTIMIZATION
 * 
 *==============================================================*/

/*t-@CDOC(idxn=TLNOObjConFunc,catn=structdef)
 * Prototype for user supplied routine for a 
 * nonlinear optimization with nonlinear constraints.
 */

typedef	int	(*TLNOObjConFunc)(
	int NVAR,		/* (I) number of variables */
	double *X,		/* (I) point where obective valued [NVAR] */
	int NEEDO,		/* (I) indicates if objective  needed */
	double *OBJF,		/* (O) value of objective */
	int NCNLN,		/* (I) number of nlin constr */
	int *NEEDC,		/* (I) array indicate if const needed [NCNLN] */
	double *CNLN,		/* (O) array of values of constr [NCNLN] */
	void *userData);	/* (I) user data */
/*e*/




/*t-@CDOC(idxn=TNLinProgObjFunc,catn=structdef)
 * These types define generic objective and constraint functions
 * for a nolinear problem with nonlinear constrains where
 * derivatives are not available.\\
 * \vspace{2mm}
 * The type {\tt NLinProgObjFunc} \index{NLinProgObjFunc}
 * defines a routine that computes an objective function
 * $$F(x_1,\dots,x_N)$$
 * of $N$ variables.
 * \vspace{2mm}
 * {\tt NLinProgConFunc} \index{NLinProgConFunc} is a routine
 * that compute a set of $\mbox{NCNLN}$ constraints
 * $$G_1(x_1,\dots,x_N),\dots,G_{\mbox{NCNLN}}(x_1,\dots,x_N).$$
 * On entry, the argument {\tt NEEDC} is ana array of length
 * $\mbox{NCNLN}$ of booleans containing TRUE or FALSE whether the
 * corresponding constraints need to be computed.
 */
typedef	int	(*TNLinProgObjFunc)(
	int	*MODE,		/* (I) not used */
	int	*N,		/* (I) number of variables */
	double	*X,		/* (I) point where obective valued */
	void	*userData,	/* (I) user data */
	double	*OBJF,		/* (O) value of the objective */
	double	*OBJGRD);	/* (O) not used */

typedef	int	(*TNLinProgConFunc)(
	int	*NCNLN,		/* (I) number of nonlinear constraints */
	int	*N,		/* (I) number of variables */
	int	*MODE,		/* (I) not used */
	int	*NEEDC,		/* (I) array indicate if constraint needed */
	double	*X,		/* (I) point wher constraints valued */
	void	*userData,	/* (I) user data */
	double	*C,		/* (O) array of values of the constraints */
	double	*cjac1);	/* (O) not used */

/*e*/




/*==============================================================*
 * Optimization routines
 *==============================================================*/


extern	DLL_EXPORT(int)	DrlNLinProg(
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
	TLNOParams *params);	/* (I) optimization parameters */



/*
 * Complex Box method.
 */

extern	DLL_EXPORT(int)	DrlNLinProgComplexBox(
	int NVAR,		/* (I) number of dimensions */
	double *blsc,		/* (I) lbounds state variables */
	double *busc,		/* (I) ubounds state variables */
	int NCNLN,		/* (I) number of nlin constr */
	double *blcnln,		/* (I) lbounds for nlin constr */
	double *bucnln,		/* (I) ubounds for nlin constr */
	TLNOObjConFunc objCon,	/* (I) objective and constraints */
	int ITER,		/* (I) maximum number of iterations */
	double *C,		/* (B) init/final  values of nlin constr */
	double *OBJF,		/* (B) initial value of obj */
	double *X,		/* (B) starting/optimal point */
	void *userData,		/* (I) user data */
	TLNOParams *params);	/* (I) optimization parameters */


/*
 * Wrappers to libraries: NAG, IMSL, etc.
 */

extern	DLL_EXPORT(int)	DrlNLinProgNagE04UCF(
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
	void	*userData);	/* (I) user data */

extern	DLL_EXPORT(char*)	DrlNLinProgNagE04UCFErrMsg(int errCode);


extern	DLL_EXPORT(int)	DrlNLinProgImslMinConNonlin(
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
	void	*userData);	/* (I) user data */

extern	DLL_EXPORT(char*)	DrlNLinProgImslMinConNonlinErrMsg(int errCode);










/*--------------------------------------------------------------
 *
 */

extern	DLL_EXPORT(int)	DrlMaxLikehoodL(
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
	double *zoptL);		/* 'F' 16 (O) optimal point [0..nVar-1] */




#endif	/* _drllno_H */
