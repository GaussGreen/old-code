/****************************************************************
 * Purpose:	Test Driver
 * Author:	C. Daher
 *
 ****************************************************************/
#include "drlstd.h"		/* otherwise DBL_EPSLION redefined by ALIB ! */

#include "drllno.h"
#include "drlmem.h"
#include "drlvtype.h"		/* DrlLilVectLogging() */

#include "drlmedis.h"		/* Prototype consistency */


/*----------------------------------------------------------------------
 *
 */

typedef	struct {
	int np;			/* number of probabilit points */
	int nc;			/* number of constraints */
	double **a;		/* linear constraints [0..nc-1][0..np-1] */
	double *c;		/* linear constraints [0..nc-1] */
	double *lm;		/* Lagrange multipliers [0..nc-1] */
	double *q;		/* prior distribution [0..np-1] */
} TMedisData;



/*----------------------------------------------------------------------
 *
 */


static	int
medisImslCallback(
        int nx,         /* (I) number of variables */
        double *x,      /* (I) point where obective valued [nx] */
        int needf,      /* (I) 0=unused, 1=needed, 2=provided */
        double *fvalue,	/* (B) value of objective */
        int *needdf,    /* (I) 0=unused, 1=needed, 2=provided [0..nx-1] */
        double *df,     /* (B) gradient of objective [0..nx-1] */
        void *userData) /* (I) user supplied data */
{
	TMedisData *medisData = (TMedisData*) userData;
	int	i, k,
		np = medisData->np;
	double	zl, e1;			/* tmp variables */
	double	*q = medisData->q,
		*c = medisData->c,
		**a = medisData->a;

	
	zl = 0e0;
	for (i=0; i<np; i++) {
		e1 = 0e0;
		for (k=0; k<nx; k++) {
			e1 += x[k] * a[k][i];
		}
		zl += q[i] * exp(e1);
	}
	*fvalue = log(zl);

	for (k=0; k<nx; k++) {
		(*fvalue) -= x[k] * c[k];
	}

	if (GtoLoggingGet() >= 1) {
		GtoErrMsg("Callback:\n");
		for (k=0; k<nx; k++) {
			GtoErrMsg("%lf ", x[k]);
		}
		GtoErrMsg(": %lf\n", (*fvalue));
	}

	return (SUCCESS);
}



/*f---------------------------------------------------------------------
 * Minimum entropy discrete distribution with linear side conditions.
 *
 * Determines the probability distribution p[i] such that the relative
 * entropy to a prior distribution q[i] (i.e. Kullback-Leiber or 
 * relative entropy) \\
 *	sum_i p[i] log(p[i] / q[i]) \\
 * is minimized under a series of nc linear side conditions of the form \\
 * 	sum_i a[k][i] * p[i] = c[k], for k=0,..,nc-1.\\
 * On exit, the optimal entropy, Lagrange multipliers and distribution
 * are returned. \\
 */

int
DrlMedisDiscrete(
	int np,			/* (I) number proba */
	double *q,		/* (I) prior distribution [0..np-1] */
	int nc,			/* (I) number of constraints */
	double **a,		/* (I) linear constraints coef [0..nc-1][0..np-1] */
	double *c,		/* (I) linear constraints value [0..nc-1] */
	int *opt,		/* (I) enable constraint flags (TRUE/FALSE) [0..nc-1] */
	double *en,		/* (O) entropy */
	double *lm,		/* (B) start/output Lagrange multipliers */
	double *p)		/* (O) output distribution [0..np-1] */
{
static	char	routine[] = "DrlMedisDiscrete";
	int	status = FAILURE;

	TMedisData	medisData;
	int		i, k;
	double		e1, zl, fvalue;

	medisData.np = np;
	medisData.nc = nc;
	medisData.a = a;
	medisData.c = c;
	medisData.lm = lm;
	medisData.q = q;

	/*
	 * Log input
	 */
	if (GtoLoggingGet() >= 1) {
		GtoErrMsg("A:\n");
		for (i=0; i<np; i++) {
			for (k=0; k<nc; k++) {
				GtoErrMsg("%lf ", a[k][i]);
			}
			GtoErrMsg("\n");
		}
	}

	/*
	 * Check input
	 */
	e1 = 0e0;
	for (i=0; i<np; i++) {
		if (q[i] < 0e0) {
			GtoErrMsg("%s: input prior dist q[%d] (%lf) != 0.\n",
				routine, i, q[i]);
			goto done;
		}
		e1 += q[i];
	}
#ifdef	SKIP
	if (fabs(e1 - 1e0) > 1e-6) {
		GtoErrMsg("%s: input prior dist does not sum to 1 "
				"(diff %e).\n", routine, e1-1e0);
		goto done;
	}
#endif

	/*
	 * Solve for the lagrange multipliers
	 */
	if (DrlLnoUnconMin(
		nc,
		lm,
		opt,
		&fvalue,
		FALSE,		/* no gradient */
		medisImslCallback, /* callback */
		(void*) &medisData,
		"IMSL",
		NULL) != SUCCESS) {
			GtoErrMsg("%s: lagrange mult optimization failed.\n",
				routine);
			goto done;
	}
	(*en) = fvalue;

	/* Logging */
	if (GtoLoggingGet() >= 1) {
		GtoErrMsg("Lagrange multpliers:\n");
		for (k=0; k<nc; k++) {
			GtoErrMsg("%lf ", lm[k]);
		}
		GtoErrMsg("\n");
	}

	/*
	 * Compute optimal distribution
	 */

	zl = 0e0;
	for (i=0; i<np; i++) {
		e1 = 0e0;
		for (k=0; k<nc; k++) {
			e1 += lm[k] * a[k][i];
		}
		zl += q[i] * exp(e1);
	}

	for (i=0; i<np; i++) {
		e1 = 0e0;
		for (k=0; k<nc; k++) {
			e1 += lm[k] * a[k][i];
		}
		p[i] = q[i] * exp(e1) / zl;
	}

	/*
	 * Verify
	 */
	for (k=0; k<nc; k++) {
	    if (opt[k]) {
		e1 = 0e0;
		for (i=0; i<np; i++) {
			e1 += a[k][i] * p[i];
		}
		e1 -= c[k];
		if (fabs(e1) >= 1e-4) {
			GtoErrMsg("%s: constraint # %d not satisfied "
				"(diff %lf).\n", routine,
				k, e1);
			goto done;
		}
	    }
	}


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}




/*f---------------------------------------------------------------------
 * LIL interface of DrlMedisDiscrete.
 */


int
DrlMedisDiscreteL(
	double *qL,	/* (I) 00 'F' prior distribution [0..np-1] */
	double *aL,	/* (I) 01 'F' linear constraints [0..nc-1][0..np-1] */
	double *cL,	/* (I) 02 'F' linear constraints [0..nc-1] */
	double *lmL,	/* (I) 03 'F' input Lagrange multipliers */
	long *optL,	/* (I) 04 'L' optimize flags */
	double *enOL,	/* (O) 05 'F' entropy */
	double *lmOL,	/* (O) 06 'F' output Lagrange multipliers */
	double *pOL)	/* (O) 07 'F' output distribution [0..np-1] */
{
static	char	routine[] = "DrlMedisDiscreteL";
	int	status = FAILURE;

	int	np;
	int	nc;
	double	*q;
	double	**a = NULL;
	int	*opt = NULL;
	double	*c;
	double	*lm;
	double	*p;
	int	i, k;

	/*
	 * Check inputs
	 */
	WRAP_CHECK_VECTOR( qL);
	WRAP_CHECK_VECTOR( aL);
	WRAP_CHECK_VECTOR( cL);
	WRAP_CHECK_VECTOR( lmL);
	WRAP_CHECK_VECTOR( optL);
	WRAP_CHECK_VECTOR( enOL);
	WRAP_CHECK_VECTOR( lmOL);
	WRAP_CHECK_VECTOR( pOL);

	/* Log wrapper inputs  */
	if (GtoLoggingGet() >= 1) {
	    DrlLilVectLoggingFp(NULL, "MEDIS_D",
		DRL_DOUBLE_L, qL,	"DIST_PRIOR",
		DRL_DOUBLE_L, aL,	"CONST_COEFF",
		DRL_DOUBLE_L, cL,	"CONST_VALUE",
		DRL_DOUBLE_L, lmL,	"LAGR_IN",
		DRL_LONG_L,   optL,	"OPTFLAG_IN",
		DRL_DOUBLE_L, enOL,	"ENT_OUT",
		DRL_DOUBLE_L, lmOL,	"LAGR_OUT",
		DRL_DOUBLE_L, pOL,	"DIST_OUT",
		(TVType)0L);
	}


	/* Unwrap arguments */
	np = ARGSIZE(qL);
	nc = ARGSIZE(cL);
	ASSERT_OR_DONE(ARGSIZE(aL) == nc*np);
	ASSERT_OR_DONE(ARGSIZE(lmL) == nc);
	ASSERT_OR_DONE(ARGSIZE(optL) == nc);

	ASSERT_OR_DONE(ARGSIZE(enOL) >= 1);
	ASSERT_OR_DONE(ARGSIZE(lmOL) >= nc);
	ASSERT_OR_DONE(ARGSIZE(pOL) >= np);

	q = &(qL[1]);
	ASSERT_OR_DONE((a = DrlDoubleMatrAlloc(0, nc-1, 0, np-1)) != NULL);
	ASSERT_OR_DONE((opt = DrlIntVectAlloc(0, nc-1)) != NULL);
	c = &(cL[1]);
	lm = &(lmOL[1]);
	for (k=0; k<nc; k++) {
		lm[k] = lmL[k+1];
	}
	for (k=0; k<nc; k++) {
		opt[k] = (int) optL[k+1];
	}
	p = &(pOL[1]);

        for (i=np-1; i>=0 ;i--)
        for (k=nc-1; k>=0; k--) {
                a[k][i] = aL[ WRAP_MATR_IDX(np, nc, i, k) ];

        }




	/* Make call to routine */
	IF_FAILED_DONE( DrlMedisDiscrete(
		np,
		q,
		nc,
		a,
		c,
		opt,
		&(enOL[1]),
		lm,
		p));


	/* output */


	/* OK */
	status = SUCCESS;
done:
	DrlDoubleMatrFree(a, 0, nc-1, 0, np-1);
	DrlIntVectFree(opt, 0, nc-1);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


