// AffineModelGen.h : generic affine model option pricing

#ifndef __AFFINEMODELGEN_H__
#define __AFFINEMODELGEN_H__

typedef struct _SAffineCoef {
  SFxSVUndFull *o;
  double T_pay;
  int idx_d, idx_f, idx_x;

  double **B, ***D;
  int n;
} SAffinCoef;

Err AffineCalcPhi(SFxSVComm_InvFT *q, double u_re, double u_im, double *h_re,
                  double *h_im);
Err AffineCalcMoments(SFxSVComm_InvFT *q, double *mean, double *std);

#endif // #ifndef __AFFINEMODELGEN_H__