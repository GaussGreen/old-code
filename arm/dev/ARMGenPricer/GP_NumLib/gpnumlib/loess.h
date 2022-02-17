/*!
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *	\file loess.h
 *
 *  \brief 
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 
 */

#ifndef LOESS_H
#define LOESS_H

#include "gpbase/port.h"
#include "gpbase/env.h"

CC_BEGIN_NAMESPACE( ARM )

// Structure for Loess

struct loess_struct {
	struct {
		long    n;
	        long    p;
	        double  *y;
	        double  *x;
		double	*weights;
	} in;
	struct {
	        double  span;
	        long    degree;
	        long    normalize;
	        long    parametric[8];
	        long    drop_square[8];
	        char    *family;
	} model;
	struct {
	        char    *surface;
	        char    *statistics;
	        double  cell;
	        char    *trace_hat;
	        long    iterations;
	} control;
	struct {
		long	*parameter;
		long	*a;
		double	*xi;
		double	*vert;
		double	*vval;
	} kd_tree;
	struct {
		double	*fitted_values;
	        double  *fitted_residuals;
		double  enp;
		double	s;
		double  one_delta;
		double	two_delta;
		double	*pseudovalues;
		double	trace_hat;
		double	*diagonal;
		double	*robust;
		double  *divisor;
	} out;
};

struct pred_struct {
	double	*fit;
	double	*se_fit;
	double  residual_scale;
	double  df;
};

struct anova_struct {
	double	dfn;
	double	dfd;
	double  F_value;
	double  Pr_F;
};

struct ci_struct {
	double	*fit;
	double	*upper;
	double  *lower;
};

int comp(double *d1, double *d2);

// Useful functions

void loess_setup(double *x, double *y, long n, long p, struct  loess_struct* lo);
void loess(loess_struct *lo);
void loess_summary(struct	loess_struct	*lo);
void loess_free_mem(struct	loess_struct *lo);
void loess_free();

void predict(double *eval, long m, struct loess_struct *lo, struct pred_struct *pre, long se);
void pred_free_mem(struct pred_struct	*pre);
void anova(struct loess_struct *one, struct loess_struct *two, struct anova_struct *out);
void pointwise(struct  pred_struct *pre, long m, double coverage, struct  ci_struct *ci);
void pw_free_mem(struct ci_struct *ci);

CC_END_NAMESPACE()

#endif