/************************************************************************
 * Module:	DCU
 * Submodule:	LNO - Linear and Nonlinear Optimization
 * File:	
 * Function:	Optimisation Routines
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include "cerror.h"
#include "cgeneral.h"
#include "macros.h"

#include "drlio.h"			/* I/O rutines */
#include "drlmem.h"

#include "drllno.h"			/* prototype consistency */

static	int	drlSobolSequence(int *n, double *x);

typedef	struct	{
	int			fNcnln,
				fNclin,
				fNIter,
				fPrintLevel;
	double			*fBlclin,
				*fBuclin,
				*fBlcnln,
				*fBucnln,
				**fMatclin;
	void			*fUserData;

	TNLinProgObjFunc	fObjfun;
	TNLinProgConFunc	fConfun;

} TDrlNLinOptimQMCData;


static	int
drlNLinOptimQMCConstFunc(
	int n,		/* (I) number of variables */
	int m, 		/* (I) total number of constraints */
	double *x,	/* (I) point where funct evaluated x[0..n-1] */
	int *active,	/* (I) array[0..max(1,m)-1] of active constraint */
	double *f,	/* (O) value computed at point x */
	double *g,	/* (I) array [0..max(1,m)-1]of constraint value */
	TDrlNLinOptimQMCData *oData);



/*f---------------------------------------------------------------------
 * Does a QMC search in a constrained optimization problem
 * <br>
 *   \min f(x_1,\dots,x_n)
 * <br>
 * subject to the constraints
 * <br>
 *   x^{l}_i &&lt;=& x_i &lt;= x^{u}_i,
 *        i=0,\dots,n \\
 *   b^{ll}_k &&lt;=& \sum_j a_{kj}x_j &lt;= b^{ul}_k,
 *        k=0,\dots,\mbox{NCLIN} \\
 *   b^{ln}_k &&lt;=& h_k(x_1,\dots,x_n) &lt;= b^{un}_k,
 *        k=0,\dots,\mbox{NCLIN} 
 * <br>
 */

DLL_EXPORT(int)
DrlNLinOptimQMC(
	int NVAR,		/* (I) number of dimensions */
	double *blsc,		/* (I) lbounds state variables */
	double *busc,		/* (I) ubounds state variables */
	int NCLIN,		/* (I) number of lin constr */
	double *blclin,		/* (I) lbounds on lin constr */
	double *buclin,		/* (I) ubounds on lin constr */
	double **a,		/* (I) coeff of lin constr */
	int NCNLN,		/* (I) number of nlin constr */
	double *blcnln,		/* (I) lbounds for nlin constr */
	double *bucnln,		/* (I) ubounds for nlin constr */
	TNLinProgConFunc nlnc,	/* (I) C-routine for nlin constr */
	TNLinProgObjFunc obj,	/* (I) C-routine for obj */
	int ITERMAX,		/* (I) maximum number of iterations */
	double *C,		/* (I) initial values of nlin constr */
				/* (O) final values of nlin constr */
	double *OBJF,		/* (B) initial/final value of obj */
	double *X,		/* (B) initial/final value */
	void *userData,		/* (I) user data */
	...)			/* (I)  other options (last arg should be 0) */
{
static	char	routine[] = "DrlNLinOptimQMC";
	int	status = FAILURE;
	int	i, k,
		NCTOT, 		/* total number of constraints */
		meq,		/* number of equality constraints */
		n;		/* number of variables */
	double	*x0 = NULL,	/* initial guess point for the optim */
		fx0,
		*x2 = NULL,
		fx2,
		*xSol = NULL,	/* final point */
		fxSol,
		*g = NULL,	/* constraints */
		*xlb = NULL,	/* array of low bounds on variable */
		*xub = NULL;	/* array of high bounds on variable */ 
	int	lFlag,
		*active = NULL;	/* vector of active flags constr */

	int	nsbdim,
		nsbiter;
	double	sbsq[30];
	TDrlNLinOptimQMCData	*oData = NULL;

	/*
	 * Parameter checking
	 */
	if ((NVAR <= 0) || (NVAR > 100)) {
		GtoErrMsg("%s: bad NVAR=%d\n", routine, NVAR);
		goto done;
	}

	if (ITERMAX <= 0) {
		GtoErrMsg("%s: bad ITERMAX=%d\n", routine, ITERMAX);
		goto done;
	}

	if ((NCLIN < 0) || (NCLIN > 1000)) {
		GtoErrMsg("%s: bad NCLIN=%d\n", routine, NCLIN);
		goto done;
	}

	if ((NCNLN < 0) || (NCNLN > 1000)) {
		GtoErrMsg("%s: bad NCNLN=%d\n", routine, NCNLN);
		goto done;
	}

	n = NVAR;		/* number of variables */
	NCTOT = 2*(NCLIN + NCNLN);	/* total number of >= constraints */
	meq = 0;		/* no equality constraint */

	if ((oData = NEW(TDrlNLinOptimQMCData)) == NULL)
		goto done;

	oData->fNclin = NCLIN;		/* static var */
	oData->fNcnln = NCNLN;
	oData->fNIter=0;

	oData->fBlclin = blclin;
	oData->fBuclin = buclin;
	oData->fMatclin = a;

	oData->fBlcnln = blcnln;
	oData->fBucnln = bucnln;

	oData->fObjfun = obj;
	oData->fConfun = nlnc;
	oData->fUserData = userData;

	oData->fPrintLevel = 0;

	fxSol = fx2 = 1e12;

	/*
	 *
	 */
	if ((x0     = DrlDoubleVectAlloc(0, n-1)) == NULL) goto done;
	if ((x2     = DrlDoubleVectAlloc(0, n-1)) == NULL) goto done;
	if ((xSol   = DrlDoubleVectAlloc(0, n-1)) == NULL) goto done;
	if ((xlb    = DrlDoubleVectAlloc(0, n-1)) == NULL) goto done;
	if ((xub    = DrlDoubleVectAlloc(0, n-1)) == NULL) goto done;
	if (NCTOT > 0) {
	    if ((g      = DrlDoubleVectAlloc(0, NCTOT-1)) == NULL) goto done;
	    if ((active = DrlIntVectAlloc(0, NCTOT-1)) == NULL) goto done;
	}



	for (i=0; i<=n-1; i++) x0[i] = X[i];


	/* Variable bounds constraints */
	for (i=0; i<=n-1; i++) {
		xlb[i] = blsc[i];
		xub[i] = busc[i];
	}

	/*
	 *
	 */
#ifndef	_WINDLL
	if (oData->fPrintLevel > 0) {
		int	i1;
		fprintf(stdout, "%s:\n", routine);
		fprintf(stdout, "NVAR = %d\n", NVAR);
		for (i1=0; i1<=n-1; i1++) {
			fprintf(stdout,
				"\t[%2d]\txLow=%lf\txHigh=%lf\tx0=%lf\n", 
				i1, xlb[i1], xub[i1], x0[i1]);
		}
		fprintf(stdout, "Linear constraints (NCLIN=%d)\n", NCLIN);
		fprintf(stdout, "Low bound:\n");
		DrlFPrintDoubleVect(stdout, NULL, oData->fBlclin, NCLIN);
		fprintf(stdout, "High bound:\n");
		DrlFPrintDoubleVect(stdout, NULL, oData->fBuclin, NCLIN);
		fprintf(stdout, "Matrix:\n");
		DrlFPrintDoubleMatr(stdout, NULL, oData->fMatclin, NCLIN, NVAR);

		fprintf(stdout, "Nonlinear constraints (NCNLN=%d)\n", NCNLN);
		fprintf(stdout, "Low bound:\n");
		DrlFPrintDoubleVect(stdout, NULL, oData->fBlcnln, NCNLN);
		fprintf(stdout, "High bound:\n");
		DrlFPrintDoubleVect(stdout, NULL, oData->fBucnln, NCNLN);

		fprintf(stdout, "Total positive constraints (NCTOT=%d)\n",
			NCTOT);
		fprintf(stdout, "Total equality constraints (meq=%d)\n", meq);
	}
#endif


	/*
	 *
	 */

	/* intialize sequence */
	nsbdim = -1;
	if (drlSobolSequence(&nsbdim, NULL) != SUCCESS)
		goto done;

	nsbdim = n;
	for (nsbiter = 0; nsbiter <= ITERMAX-1; nsbiter++) {

		/* get Sobol deviates */
		if (drlSobolSequence(&nsbdim, sbsq) != SUCCESS)
			goto done;

		/* build try point */
		for (i=0; i<=n-1; i++) {
			x0[i] = xlb[i] + (xub[i] - xlb[i]) * sbsq[i];
		}


		/* check linear constraints */
		for (k=0; k<=2*NCLIN-1; k++)
	   		active[k] = 1;
		for (k=2*NCLIN; k<=NCTOT-1; k++)
	   		active[k] = 0;

		if (drlNLinOptimQMCConstFunc(n, NCTOT, x0, active,
			&fx0, g, oData) != SUCCESS)
				goto done;

		lFlag = 0;
		for (k=0; k<=2*NCLIN-1; k++)
			lFlag += (g[k] < 0e0);
		if (lFlag > 0) {
			continue;
		}

		/* check linear constraints */
		for (k=0; k<=2*NCLIN-1; k++)
	   		active[k] = 0;
		for (k=2*NCLIN; k<=NCTOT-1; k++)
	   		active[k] = 1;

		if (drlNLinOptimQMCConstFunc(n, NCTOT, x0, active,
			&fx0, g, oData) != SUCCESS)
				goto done;

		lFlag = 0;
		for (k=2*NCLIN; k<=NCTOT-1; k++)
			lFlag += (g[k]  < 0e0);
		if (lFlag > 0) {
			continue;
		}


		/* point is feasible: check value  */
		for (k=0; k<=NCTOT-1; k++)
	   		active[k] = 0;
		if (drlNLinOptimQMCConstFunc(n, NCTOT, x0, active,
			&fx0, g, oData) != SUCCESS)
				goto done;

		if (fx0 < fxSol) {
			/* save 2nd best */
			for (i=0; i<=n-1; i++) x2[i] = x0[i];
			fx2 = fx0;
			for (i=0; i<=n-1; i++) xSol[i] = x0[i];
			fxSol = fx0;
		} else if (fx0 < fx2) {
			for (i=0; i<=n-1; i++) x2[i] = x0[i];
			fx2 = fx0;
		}
	}

	/*
	 *
	 */
#ifndef	_WINDLL
	if (oData->fPrintLevel > 0) {
		fprintf(stdout, "\tOptim\t2ndBest\tDiff\n");
		for (i=0; i<=n-1; i++)  {
			fprintf(stdout, "\t%lf\t%lf\t%g\n",
				xSol[i], x2[i], xSol[i] - x2[i]);
			fflush(stdout);
		}
	}
#endif
	if (fxSol > 1e10) {
		GtoErrMsg("%s: no feasible point found.\n", routine);
		goto done;
	} else {
		for (i=0; i<=n-1; i++) 
			X[i] = xSol[i];
	}


	/* made it through OK */
	status = SUCCESS;
done:
	if (oData) FREE(oData);
	DrlDoubleVectFree(x0, 0, n-1);
	DrlDoubleVectFree(x2, 0, n-1);
	DrlDoubleVectFree(xSol, 0, n-1);
	DrlDoubleVectFree(xlb, 0, n-1);
	DrlDoubleVectFree(xub, 0, n-1);
	DrlDoubleVectFree(g, 0, NCTOT-1);
	DrlIntVectFree(active, 0, NCTOT-1);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*----------------------------------------------------------------------
 * Constraint patch function:
 * there is a total of m = 2*(NCLIN+NCNLN) positive constraints:
 *	m=0		1st lin LOW
 *	m=1		1st lin HIGH
 *	m=2		2nd lin LOW
 *	etc..
 *	m=2*NCNLN	1st nlin LOW
 *	m=2*NCNLN+1	1st nlin HIGH
 *	rtc...
 */

static	int
drlNLinOptimQMCConstFunc(
	int n,		/* (I) number of variables */
	int NCTOT,	/* (I) total number of constraints */
	double *x,	/* (I) point where funct evaluated x[0..n-1] */
	int *active,	/* (I) array[0..max(1,m)-1] of active constraint */
	double *f,	/* (O) value computed at point x */
	double *g,	/* (I) array [0..max(1,m)-1]of constraint value */
	TDrlNLinOptimQMCData *oData)	/* (I) user data */
{
static	char	routine[] = "drlNLinOptimQMCConstFunc";
	double	val;
	int	MODE = 0,
		errCode = 0,
		i, j;

static	double	tmpConstr[100];
static	int	tmpNeedc[100];

	/*
	 * Computation of the objective
	 */
	errCode = oData->fObjfun(
		&MODE,
		&n,
		x,
		oData->fUserData,
		f,
		NULL);
	if (errCode != 0) {
	     GtoErrMsg("%s: objective function returned "
		    "code %d (warning).\n", routine, errCode);
	}


	/*
	 * Linear constraints
	 */
	for (i=0; i<=oData->fNclin-1; i++) {
	   /* check if constraint active */
	   if ((active[2*i]) || (active[2*i+1])) {
		val = 0e0;
		for (j=0; j<=n-1; j++)
			val += oData->fMatclin[i][j] * x[j];

	  	if (active[2*i  ]) g[2*i  ] =   val - oData->fBlclin[i];
	  	if (active[2*i+1]) g[2*i+1] = - val + oData->fBuclin[i];
	   }
	}

	/*
	 * Nonlinear constraints
	 */
	for (i=0; i<=oData->fNcnln-1; i++) {
		tmpNeedc[i] = ((active[2*(oData->fNclin+i)]) ||
					(active[2*(oData->fNclin+i)+1]));
	}

	errCode = oData->fConfun(
		&oData->fNcnln,
		&n,
		&MODE,
		tmpNeedc,
		x,
		oData->fUserData,
		tmpConstr,
		(double*) NULL);

	if (errCode != 0) {
	     GtoErrMsg("%s: objective function returned "
		    "code %d (warning).\n", routine, errCode);
	}

#ifndef	_WINDLL
	if (oData->fPrintLevel >= 2) {
		fprintf(stdout, "NCNLN: ");
		for (i=0; i<=oData->fNcnln-1; i++)
			fprintf(stdout, "\t%.3g", tmpConstr[i]);
		fprintf(stdout, "\n");

		fprintf(stdout, "BLCNLN:");
		for (i=0; i<=oData->fNcnln-1; i++)
			fprintf(stdout, "\t%.3g", oData->fBlcnln[i]);
		fprintf(stdout, "\n");

		fprintf(stdout, "BUCNLN:");
		for (i=0; i<=oData->fNcnln-1; i++)
			fprintf(stdout, "\t%.3g", oData->fBucnln[i]);
		fprintf(stdout, "\n");

	}
#endif
	for (i=0; i<=oData->fNcnln-1; i++) {
	  if ((active[2*(oData->fNclin+i)]) ||
	      (active[2*(oData->fNclin+i)+1])) {
	  	if (active[2*(oData->fNclin+i)  ] != 0)
			g[2*(oData->fNclin+i)  ] =
				tmpConstr[i] - oData->fBlcnln[i];
	  	if (active[2*(oData->fNclin+i)+1] != 0)
			g[2*(oData->fNclin+i)+1] =
				- tmpConstr[i] + oData->fBucnln[i];
	  }
	}

#ifndef	_WINDLL
	if (oData->fPrintLevel >= 2) {
		fprintf(stdout, "<%3d>", ++oData->fNIter);
		for (i=0; i<=n-1; i++) fprintf(stdout, "\t%.3g", x[i]);
		fprintf(stdout, "\n");

		fprintf(stdout, "CSTNUM");
		for (i=0; i<=NCTOT-1; i++) fprintf(stdout, "\t%d", i);
		fprintf(stdout, "\n");

		fprintf(stdout, "CSTACT:");
		for (i=0; i<=NCTOT-1; i++) fprintf(stdout, "\t%d", active[i]);
		fprintf(stdout, "\n");

		fprintf(stdout, "CSTVAL");
		for (i=0; i<=NCTOT-1; i++) fprintf(stdout, "\t%.3g", g[i]);
		fprintf(stdout, "\n");
	}
#endif


	/*
	 *
	 */
	return(0);
}





/************************************************************************
 * Module:	DRL - RAND
 * Function:	Random Number Generation and Simulation
 *		Sobol Sequence
 * Author:	NRC2 - modified by C. Daher
 * Revision:	$Header$
 ************************************************************************/

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define MAXBIT 30
#define MAXDIM 6
static unsigned long	inS,
			ixS[MAXDIM+1],
			*iuS[MAXBIT+1];
static unsigned long	mdegS[MAXDIM+1]={0,1,2,3,3,4,4};
static unsigned long	ipS[MAXDIM+1]={0,0,1,1,2,1,4};
static unsigned long	ivS[MAXDIM*MAXBIT+1]={
				0,1,1,1,1,1,1,3,1,3,3,
				1,1,5,7,7,3,3,5,15,11,
				5,15,13,9};


/*----------------------------------------------------------------------
 * Generates pseudo random Sobol sequence numbers.
 * When $n < 0$, internally initialize a set of direction numbers.
 * When $0<n<\mbox{MAXDIM}=30$, returns as vector <i> x[0..n-1]</i> the
 * next values form $n$ of these sequences
 * ($n$ must not be changed between initializations).
 * Returns 0 iff OK.
 */

static	int
drlSobolSequence(int *n, double *x)
{
static	char	routine[] = "drlSobolSequence";
	int		j,k,l;
	unsigned long	i,im,ipp;
static double		fac;
static unsigned long	in,ix[MAXDIM+1],*iu[MAXBIT+1];
static unsigned long	mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
static unsigned long	ip[MAXDIM+1]={0,0,1,1,2,1,4};
static unsigned long	iv[MAXDIM*MAXBIT+1]={ 0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

	/* Because NRC vectors are [1..n] */
	x--;

#ifdef	__DEBUG__
	printf("Sobol: n=%d");
#endif

	if (*n < 0) {
		/*
		 * Reset all static vectors
		 */
#undef	COPYVECT
#define	COPYVECT(x,y,n)	{for(i=0; i<=(n)-1; i++) *((x)+i) = *((y)+i);}
		in = inS;
		COPYVECT(ix, ixS, MAXDIM+1);
		COPYVECT(iu, iuS, MAXBIT+1);
		COPYVECT(mdeg, mdegS, MAXDIM+1);
		COPYVECT(ip, ipS, MAXDIM+1);
		COPYVECT(iv, ivS, MAXDIM*MAXBIT+1);
#undef	COPYVECT


		/*
		 * Initialize
		 */
		for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) {
			iu[j] = &iv[k];
		}
		for (k=1;k<=MAXDIM;k++) {
			for (j=1;j<=(int)mdeg[k];j++) {
				iu[j][k] <<= (MAXBIT-j);
			}
			for (j=mdeg[k]+1;j<=MAXBIT;j++) {
				ipp=ip[k];
				i=iu[j-mdeg[k]][k];
				i ^= (i >> mdeg[k]);
				for (l=mdeg[k]-1;l>=1;l--) {
					if (ipp & 1) i ^= iu[j-l][k];
					ipp >>= 1;
				}
				iu[j][k]=i;
			}
		}
		fac=1.0/(1L << MAXBIT);
		in=0;
	} else {
		im=in;
		for (j=1;j<=MAXBIT;j++) {
			if (!(im & 1)) break;
			im >>= 1;
		}
		if (j > MAXBIT) {
			GtoErrMsg("%s: MAXBIT too small in sobseq\n",
				routine);
			return(1);
		}
		im=(j-1)*MAXDIM;
		for (k=1;k<=IMIN(*n,MAXDIM);k++) {
			ix[k] ^= iv[im+k];
			x[k]=ix[k]*fac;
		}
		in++;
	}
	return(0);
}
#undef MAXBIT
#undef MAXDIM
#undef NRANSI



