/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	MEDIS - Minimum entropy distribution
 * File:	drlmedis.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlmedis_H
#define	_drlmedis_H
#include "drlstd.h"

/*
 * Minimum entropy discrete distribution with linear side conditions.
 */
int	DrlMedisDiscrete(
	int np,			/* (I) number proba */
	double *q,		/* (I) prior distribution [0..np-1] */
	int nc,			/* (I) number of constraints */
	double **a,		/* (I) linear constraints [0..nc-1][0..np-1] */
	double *c,		/* (I) linear constraints [0..nc-1] */
	int *opt,		/* (I) optimize flags */
	double *en,		/* (O) entropy */
	double *lm,		/* (B) output Lagrange multipliers */
	double *p);		/* (O) output distribution [0..np-1] */


/*
 * LIL interface of DrlMedisDiscrete.
 */
int	DrlMedisDiscreteL(
	double *qL,	/* (I) 00 'F' prior distribution [0..np-1] */
	double *aL,	/* (I) 01 'F' linear constraints [0..nc-1][0..np-1] */
	double *cL,	/* (I) 02 'F' linear constraints [0..nc-1] */
	double *lmL,	/* (I) 03 'F' input Lagrange multipliers */
	long *optL,	/* (I) 04 'L' optimize flags */
	double *enOL,	/* (O) 05 'F' entropy */
	double *lmOL,	/* (O) 06 'F' output Lagrange multipliers */
	double *pOL);	/* (O) 07 'F' output distribution [0..np-1] */



#endif	/* _drlmedis_H */


