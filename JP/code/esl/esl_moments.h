#ifndef ESL_MOMENTS_DOT_H
#define ESL_MOMENTS_DOT_H

/** NOTE: This file should be only included through 'esl_moments.c'
 */

#ifdef  __cplusplus
extern "C" {
#endif

int TreeMean_IrExt (double *, int, double *, double, double *, double *, double *);
int TreeMean_IrNum (double *, int);
int TreeMean_FxExt (double *, int, double *, double *, double, double *, double *, double, double *, double *, double *, double *, double *);
int TreeMean_FxNum (double *, int, double *, double *, double, double *, double *, double, double *, double *, double *);
int TreeMean_EqExt (double *, int, double *, double *, double, double *, double *, double *);
int TreeMean_EqNum (double *, int, double *, double *, double, double *);

int TreeCovariance_IrIr (double *, int, double *, double, double *, double, double *, double *);
int TreeCovariance_IrFx (double *, int, double *, double, double *, double *, double, double *, double *, double, double *, double *, double *, double *, double *);
int TreeCovariance_IrEq (double *, int, double *, double, double *, double *, double, double *, double *, double *, double *);
int TreeCovariance_FxFx (double *, int, double *, double *, double, double *, double *, double, double *, double *, double, double *, double *, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int TreeCovariance_FxEq (double *, int, double *, double *, double, double *, double *, double, double *, double *, double, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int TreeCovariance_EqEq (double *, int, double *, double *, double, double *, double *, double, double *, double *, double *, double *, double *, double *, double *);

#ifdef  __cplusplus
}
#endif

#endif
