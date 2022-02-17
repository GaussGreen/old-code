/****************************************************************************/
/*      HYB4_LIB.h                                                          */
/****************************************************************************/

#include "eslhead.h"

#ifndef _HYB4_LIB_H
#define _HYB4_LIB_H

#include "cupsmodl.h"
#include "hyb4_model.h"

#ifdef __cplusplus
extern "C" {
#endif

/* hyb4_alloc.c */

void    Hyb4_Tree_Init         (HYB4_TREE_DATA *);
int     Hyb4_Tree_Alloc        (HYB4_TREE_DATA *, HYB3_TREE_DATA *);
int     Hyb4_Tree_Free         (HYB4_TREE_DATA *);
int     Hyb4_Tree_Limits_Alloc (HYB4_TREE_DATA *, HYB3_TREE_DATA *);
int     Hyb4_Tree_Limits_Free  (HYB4_TREE_DATA *);
int     Hyb4_Tree_Limits_Assign(HYB4_TREE_DATA *, HYB3_TREE_DATA *);
int     Hyb4_Tree_Limits_Unassign(HYB4_TREE_DATA *);
int     Hyb4_Offset_Alloc      (HYB4_TREE_DATA *);
int     Hyb4_Offset_Free       (HYB4_TREE_DATA *);
void    Hyb4_Dev_Init          (HYB4_DEV_DATA *);
int     Hyb4_Dev_Alloc         (HYB4_DEV_DATA *, HYB4_TREE_DATA *);
int     Hyb4_Dev_Free          (HYB4_DEV_DATA *, HYB4_TREE_DATA *);
void    Hyb4_Asset_IR_Init     (ASSET_IR *);
int     Hyb4_Asset_IR_Alloc    (ASSET_IR *, HYB4_TREE_DATA *);
int     Hyb4_Asset_IR_Free     (ASSET_IR *, HYB4_TREE_DATA *);
void    Hyb4_Asset_EQ_Init     (ASSET_EQ *);
int     Hyb4_Asset_EQ_Alloc    (ASSET_EQ *, HYB4_TREE_DATA *);
int     Hyb4_Asset_EQ_Free     (ASSET_EQ *, HYB4_TREE_DATA *);
void    Hyb4_Asset_FX_Init     (ASSET_FX *);
int     Hyb4_Asset_FX_Alloc    (ASSET_FX *, HYB4_TREE_DATA *);
int     Hyb4_Asset_FX_Free     (ASSET_FX *, HYB4_TREE_DATA *);
TSLICE  Hyb4_Alloc_Slice       (HYB4_TREE_DATA *,int);      
void    Hyb4_Free_Slice        (TSLICE,HYB4_TREE_DATA *,int);

/* hyb4_bond.c */

int     Hyb4_Bond_Price(TSLICE,double,long,double,long,int,int,int,int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);

/* hyb4_claimbank.c */

void    Hyb4_CbkInit(CLAIM_BANK *);
int     Hyb4_CbkAlloc(CLAIM_BANK *,int,HYB4_TREE_DATA *, int);
int     Hyb4_CbkFree(CLAIM_BANK *,HYB4_TREE_DATA *, int);
int     Hyb4_CbkSizeFromDL(long *,int,long *,int);
int     Hyb4_CbkCalcSize(int,CRIT_DATE *,int,int);
int     Hyb4_CbkDev(CLAIM_BANK *,long,int,int,int, int, HYB4_DEV_DATA *,HYB4_TREE_DATA *);
TSLICE  Hyb4_CbkPopSlice(CLAIM_BANK *);
int     Hyb4_CbkPushSlice(CLAIM_BANK *, TSLICE,long,long);
TSLICE  Hyb4_CbkReadSlice(CLAIM_BANK *,long);
int     Hyb4_CbkGetOffset(CLAIM_BANK *,long,int);

/* hyb4_copy.c */

int     Hyb4_CopyTreeCups (HYB4_TREE_DATA *, HYB3_TREE_DATA *);
int     Hyb4_CopyTreeEq   (HYB4_TREE_DATA *, HYB3_TREE_DATA *, int, int);
int     Hyb4_CopyTreeFx   (HYB4_TREE_DATA *, HYB3_TREE_DATA *, int);
int     Hyb4_CopyTreeFx2FxEq(HYB4_TREE_DATA *, HYB3_TREE_DATA *, int, double);
int     Hyb4_CopyTreeEq3d (HYB4_TREE_DATA *, HYB3_TREE_DATA *, int, int);
int     Hyb4_CopyTreeFxEq (HYB4_TREE_DATA *, HYB3_TREE_DATA *, int);
int     Hyb4_CopyAssetIR  (HYB4_TREE_DATA *, HYB3_TREE_DATA *, ASSET_IR *, MKTVOL_DATA *, int);
int     Hyb4_CopyAssetEQ  (HYB4_TREE_DATA *, HYB3_TREE_DATA *, ASSET_EQ *);
int     Hyb4_CopyAssetFX  (HYB4_TREE_DATA *, HYB3_TREE_DATA *, ASSET_FX *);

/* hyb4_dev.c */

int     Hyb4_Dev(TSLICE,int,int,int,int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);

/* hyb4_floater.c */

int     Hyb4_Floater_Price (TSLICE,int,TSLICE,TSLICE,long,double,double,char,char,double,double,long,int,int,int,int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);
int     Hyb4_FwdFloater_t(TSLICE,TSLICE,int,TSLICE,int,TSLICE,long,double,double,long,char,double,char,long,char,char,long,int,int,int,int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);
int     Hyb4_Collaret_t(TSLICE,TSLICE,TSLICE,int,double,double,double,char,double,char,double,int,int,int,int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);
int     Hyb4_Floater_Pmt (TSLICE,int,TSLICE,TSLICE,double,double,double,double,int,int,HYB4_TREE_DATA *);

/* hyb4_koopt.c */

int     Hyb4_KoOption_t (TSLICE,int,TSLICE,long,double,double,int,TSLICE,char,char,int,int,HYB4_TREE_DATA *);

/* hyb4_lattice.c */

int     Hyb4_ConfigureLatticeCode(HYB4_TREE_DATA *, LATTICE_PROG *);
int     Hyb4_Lattice(HYB4_DEV_DATA *, HYB4_TREE_DATA *, LATTICE_PROG *, int, int);

/* hyb4_moments.c */

int     Hyb4_Moments_Analytic  (HYB4_TREE_DATA *);
int     Hyb4_Moments_Empirical (HYB4_TREE_DATA *, LATTICE_PROG *);
int     Hyb4_Moments_Print     (HYB4_TREE_DATA *, char *, char *);
int     Hyb4_StateProbs_Print  (HYB4_TREE_DATA *, char *);
int     Hyb4_Limits_Print      (HYB4_TREE_DATA *, char *);
int     Hyb4_AssetPayoff_Empirical (HYB4_TREE_DATA *, LATTICE_PROG *, int, int, double *);

/* hyb4_opt.c */

int     Hyb4_Option_t(TSLICE,TSLICE,double,double,long,int,int,int,int,int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);

/* hyb4_paryield.c */

int     Hyb4_Par_Yield_t(TSLICE,int,int,TSLICE *,long const*,long,long,long,int,char,char,double,int,HYB4_TREE_DATA const*);

/* hyb4_slice.c */

int     Hyb4_Init_Slice (double *,int,double,HYB4_TREE_DATA *);
int     Hyb4_SetSlice(TSLICE,int,double,int,HYB4_TREE_DATA *);
int     Hyb4_CopySlice(TSLICE,TSLICE,int,int,HYB4_TREE_DATA *);
int     Hyb4_ExpandSlice(TSLICE,int,TSLICE,int,int,HYB4_TREE_DATA *);
int     Hyb4_AddTwoSlices(TSLICE,int,TSLICE,TSLICE,int,HYB4_TREE_DATA *);
int     Hyb4_MultiplyTwoSlices(TSLICE,int,TSLICE,TSLICE,int,HYB4_TREE_DATA *);
int     Hyb4_LCombTwoSlices(TSLICE,int,TSLICE,double,TSLICE,double,int,HYB4_TREE_DATA *);
int     Hyb4_AddScalar(TSLICE,int,double,int,HYB4_TREE_DATA *);
int     Hyb4_MultiplyScalar(TSLICE,int,double,int,HYB4_TREE_DATA *);
int     Hyb4_InvertSlice(TSLICE,int,int,HYB4_TREE_DATA *);
int     Hyb4_MinOnSlice(TSLICE,int,double,int,HYB4_TREE_DATA *);
int     Hyb4_MaxMinOnSlice(TSLICE,int,double,double,int,HYB4_TREE_DATA *);
double  Hyb4_GetValueAtNode(int,TSLICE,int,int,int,int,int,HYB4_TREE_DATA *);
int     Hyb4_SmoothStepUp(double*,double*,double,double,int,int,int,HYB4_TREE_DATA *);
int     Hyb4_Offset_Calc (int, HYB4_TREE_DATA *);

/* hyb4_stateprices.c */

int     Hyb4_UpdateStatePrices4D(int, int, HYB4_TREE_DATA *, HYB4_DEV_DATA *, LATTICE_PROG *, TSLICE, TSLICE);

/* hyb4_turbo.c */

int     Hyb4_Turbo_Price(TSLICE,TSLICE,double,double,char,TSLICE,double,double,char,double,double,double,double,TSLICE,int,char,int,double,int,double,TSLICE,int,int,int,HYB4_DEV_DATA   *,HYB4_TREE_DATA  *);
int     Hyb4_FwdTurbo_t(TSLICE,TSLICE,TSLICE,double,double,char,TSLICE,double,double,char,double,double,double,double,TSLICE,int,long,char,char,char,int,double,char,char,long,int,int,int,HYB4_DEV_DATA  *,HYB4_TREE_DATA *) ;

/* hyb4_util.c */

double  Hyb4_GetIndexStep (TSLICE,int,int,int,int,int,int,HYB4_TREE_DATA *);

/* hyb4_zerobank.c */

int     Hyb4_ZbkUpdate(CLAIM_BANK *,int,long,long,int,int, int, int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);
int     Hyb4_FXZbkUpdate(CLAIM_BANK *,int,long,long,int,int, int, int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);
double *Hyb4_ZbkReadZero(CLAIM_BANK *,long,int,long, int, int,HYB4_TREE_DATA *);
int     Hyb4_ZbkParYield_t(double *,double *, int,CLAIM_BANK *,long,long,int,char,char,double,int,HYB4_TREE_DATA *);
int     Hyb4_ZbkAnnuity_t(double *,int,CLAIM_BANK *,long,long,int,char,char,int,HYB4_TREE_DATA *);

/* hyb4_zeros.c */

int     Hyb4_Zero_t(TSLICE,long,int,int,int,int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);
int     Hyb4_FXZero_t (TSLICE,long,int,int,int,int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);
int     Hyb4_Zero_Bank(TSLICE *,long *,int *,int,long,long,int,int,int,int,HYB4_DEV_DATA *,HYB4_TREE_DATA *);

/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif /* HYB4_LIB_H */

