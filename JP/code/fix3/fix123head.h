/****************************************************************************/
/*      Function templates for library.                                     */
/****************************************************************************/
/*      FIX123HEAD.h                                                        */
/****************************************************************************/

/*
$Header$
*/


/* Use safe muliple inclusion for Jerry Cohen's stream libraries */

#ifndef _fix123head_h
#define _fix123head_h

#include "eslhead.h"
#include "fix123.h"

#ifdef __cplusplus
extern "C" {
#endif

/* alloc.c */
void    Fix3_Tree_Init (FIX3_TREE_DATA *);
int     Fix3_Tree_Alloc(FIX3_TREE_DATA *);
int     Fix3_Tree_Free (FIX3_TREE_DATA *);
void    Fix3_Dev_Init  (FIX3_DEV_DATA    *dev_data);
int     Fix3_Dev_Alloc (FIX3_DEV_DATA *,FIX3_TREE_DATA const*);
int     Fix3_Dev_Free (FIX3_DEV_DATA *,FIX3_TREE_DATA const*);
double  *Fix3_Alloc_Slice (FIX3_TREE_DATA const*);
int     Fix3_Free_Slice (double *,FIX3_TREE_DATA const*);
int     *Fix3_Alloc_Slice_Int (FIX3_TREE_DATA const*);
int     Fix3_Free_Slice_Int (int *,FIX3_TREE_DATA const*);

/* amortloan.c */
int     Fix3_AnnAmortSwap_t (double **,double *,double *,double *,char,char,double,double,double,long,long,long,int,double,double,double **,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* bond.c */
int     Fix3_Bond_t   (double *,double,long,double,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_FwdBond_t(double *,double *,double *,long,double,double,char,long,double,double,long,double,char,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* callzero.c */
int     Fix3_Callzero_t (double *,double *,double *,double *,double *,long,long,long,double,double,double,char,char,char,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* cap.c */
int     Fix3_Cap_t (double *,double *,double *,int,long,double,double,char,double,char,double,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_Caplet_t (double *,double *,double *,int,long,double,double,char,double,char,double,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_Collaret_t (double *,double *,double *,long,double,double,double,char,double,char,double,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* cet.c */
int     Fix3_Cet_Classic(int,T_CURVE const*,MKTVOL_DATA *,FIX3_TREE_DATA const*);
int     Fix3_Cet_Manager(MKTVOL_DATA const*,FIX3_TREE_DATA const*,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Cet_Model_Manager_Classic(MKTVOL_DATA const*,FIX3_TREE_DATA const*,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Cet_Model_Manager_TimeDep(MKTVOL_DATA const*,FIX3_TREE_DATA const*,MKTVOL_DATA *,FIX3_TREE_DATA *); 
int     Fix3_Cet_Model_Manager_Smd(MKTVOL_DATA const*,FIX3_TREE_DATA const*,MKTVOL_DATA *,FIX3_TREE_DATA *); 
int     Fix3_Cet_Schedule(long, CET_OUT_DATA *,MKTVOL_DATA *,FIX3_TREE_DATA const*,FIX3_TREE_DATA *); 
int     Fix3_Cet_Calc(MKTVOL_DATA *,FIX3_TREE_DATA *,CET_OUT_DATA *);
void    Fix3_Cet_Tree_Reset(FIX3_TREE_DATA *);
int     Fix3_Cet_Print_Vol_File(CET_OUT_DATA *,MKTVOL_DATA *,int,int);
int     Fix3_Cet_Model_Manager_E2Q(MKTVOL_DATA const* , FIX3_TREE_DATA const*, MKTVOL_DATA* ,FIX3_TREE_DATA* );     

/* cet_tmx.c */
int     Fix3_Cet_Tmx(int,T_CURVE const*,MKTVOL_DATA *,FIX3_TREE_DATA const*);
int     Cet_Calc (T_CURVE const*,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Cet_Manager (FIX3_TREE_DATA const*,FIX3_TREE_DATA *); 
int     Cet_NmrVanl (long,T_CURVE const*,MKTVOL_DATA *,int);
int     Cet_Print (int,double *,MKTVOL_DATA *);
int     Cet_Schedule (long,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Cet_SwapVols (int,int , double *, double *, double *, double *, MKTVOL_DATA *, FIX3_TREE_DATA *, FIX3_DEV_DATA *, int); 
int     Cet_UpdSwapSml(long,MKTVOL_DATA *);

/* chooser.c */
int     Fix3_Chooser_t (double **,double *,long,int,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_Survivor_t (double **,double *,double *,long,int,int,double,char,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* claimbank.c */
void    Fix3_CbkInit(CLAIM_BANK *);
int     Fix3_CbkAlloc(CLAIM_BANK *,int,FIX3_TREE_DATA *);
int     Fix3_CbkFree(CLAIM_BANK *,FIX3_TREE_DATA *);
int     Fix3_CbkSizeFromDL(long *,int,long *,int);
int     Fix3_CbkCalcSize(int,CRIT_DATE *,int,int);
int     Fix3_CbkDev(CLAIM_BANK *,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
double *Fix3_CbkPopSlice(CLAIM_BANK *);
int     Fix3_CbkPushSlice(CLAIM_BANK *,double *,long,long);
double *Fix3_CbkReadSlice(CLAIM_BANK *,long);
int     Fix3_CbkGetOffset(CLAIM_BANK const*,long,int);

/* dev.c */int     Fix3_Dev_Classic (double *,int,int,int,FIX3_DEV_DATA const*,FIX3_TREE_DATA const*);
int     Fix3_Ev  (double *,int,int,    FIX3_DEV_DATA const*,FIX3_TREE_DATA const*);

/* dev_tmx.c */
int     Fix3_Dev_Tmx (double *,int,int,int,FIX3_DEV_DATA const*,FIX3_TREE_DATA const*);

/* drift.c */
int     Fix3_Drift_1D_Classic (MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Drift_2D_Classic (MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Drift_3D_Classic (MKTVOL_DATA *,FIX3_TREE_DATA *);

/* drift_timedep.c */
int     Fix3_Drift_1D_TimeDep (MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Drift_2D_TimeDep (MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Drift_3D_TimeDep (MKTVOL_DATA *,FIX3_TREE_DATA *);

/* drift_smd.c */
int     Fix3_Drift_2D_Smd (MKTVOL_DATA *,FIX3_TREE_DATA *);

/* drift_tmx.c */
int     Fix3_Drift_1D_Tmx (MKTVOL_DATA *,FIX3_TREE_DATA *);

/* drift_e2q.c */
void	computeMandS ( double, double, double, double, const char*, double, double, double, double*, double*, double*);
double	getU ( double, double, double, double, double);
void	ExpansionParameters (double, double, double, double, double, double*, double*, double*, double*, double*, double*);
void	computePolynomialCoeff (double, double, double, double, double, double, double, int, int, int, double*	const, double*	const, double	P[4]);					
int     Fix3_Drift_1D_E2Q (MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Drift_2D_E2Q (MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Drift_3D_E2Q (MKTVOL_DATA *,FIX3_TREE_DATA *);

/* warning.c */
void    DR_Warning (int, char *, ...);

/* fix3model.c -- vol bootstrapping interface */
extern int (*Fix3_BFactor) (double*, long ,long, char, char, MKTVOL_DATA* , T_CURVE const*);
extern int (*Fix3_SpotVol) (MKTVOL_DATA*, T_CURVE const*);
extern int (*Fix3_Interp_SpotVol) (FIX3_TREE_DATA*, MKTVOL_DATA*);
extern int (*Fix3_IndexVol) (double *, long, long, char, char,MKTVOL_DATA*,T_CURVE const*);
extern int (*Fix3_Filtered_SpotVol) (MKTVOL_DATA *, T_CURVE const*);
extern int (*Fix3_GenIndexVol) (double *, long, long, long, long, char, char, MKTVOL_DATA *, T_CURVE const* );   
extern int (*Fix3_Cet_Main) (int,T_CURVE const*,MKTVOL_DATA *,FIX3_TREE_DATA const*);
extern int (*Fix3_IndexLimits) (double* ,double*, long, long, long, int, char, char,  MKTVOL_DATA*, T_CURVE const*);
/* fix3model.c -- tree building interface */
extern int (*Fix3_Tree_Limits) (int *,int *,int *,int *,int **,int **,int ***,int ***,int, FIX3_TREE_DATA *tree_data, MKTVOL_DATA *mktvol_data);
extern int (*Fix3_Drift_1D) (MKTVOL_DATA *,FIX3_TREE_DATA *);
extern int (*Fix3_Drift_2D) (MKTVOL_DATA *,FIX3_TREE_DATA *);
extern int (*Fix3_Drift_3D) (MKTVOL_DATA *,FIX3_TREE_DATA *);
extern int (*Fix3_Lattice) (FIX3_DEV_DATA *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *);
extern int (*Fix3_Dev) (double *,int,int,int,FIX3_DEV_DATA const*,FIX3_TREE_DATA const*);
/* fix3model.c -- model I/O interface */
extern int (*Fix3_Param_Check) (int,MKTVOL_DATA *,FIX3_TREE_DATA *);
extern int (*Fix3_Cet_Model_Manager) (MKTVOL_DATA const*,FIX3_TREE_DATA const*,MKTVOL_DATA *,FIX3_TREE_DATA *);  
extern int (*Fix3_PackModelData) (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
extern int (*Fix3_UnPackModelData) (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int         Fix3_Model_Interface_Init (int);

/* flexswap.c */
int     Fix3_Flexswap_t (double **,double **,long *,char *, long,double,int,double *,double,double,double,double,long,double *,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* floater.c */
int     Fix3_Floater_t (double *,double *,double *,long,double,double,char,char,double,double,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_FwdFloater_t (double *,double *,double *,double *,long,double,double,long,char,double,char,long,char,char,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_CmpFloater_t (double *,double *,double *,long,long,double,double,char,char,double,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* idxam.c */
int     Fix3_IdxamSwap_t (double **,long,long,long,double *,double *,double,double,double,double,double,double,double *,double *,int,char,char,char,int,double,double,double **,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);


/* irrloan.c */
int     Fix3_IRRLoanSwap_t (double **,long,long,double *,double *,double *,double,double,double,double,double,double,double,char,char,int,int,double,double,double **,int,FIX3_TREE_DATA *);


/* iou.c */
int     Fix3_Ioucap_t (double **,double *,double *,long,int,double,double,char,char,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_Iouswap_t(double **,double *,double *,long,int,double,double,double,double **,char,char,char,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);


/* koopt.c */
int     Fix3_KoOption_t (double *,double *,long,double,double,double,char,char,int,FIX3_TREE_DATA *);
int     Fix3_KoOptionEps_t (double *, double *, long, double, double, double, double, double, char, char, int, FIX3_TREE_DATA   *tree_data);
int     Fix3_KiOption_t (double *,double *,double *,long,double,double,char,char,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_Trigger_t (double *,double *,double *,double,double,long,int,double,double,char,char,int,int,int,FIX3_DEV_DATA *, FIX3_TREE_DATA *);
int     Fix3_Top_t(double *,double *,double *,double,double,long,int,double,double,char,char,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* kostream.c */
int     Fix3_KoCap_t (double *,double *,double *,double *,double *,long,double,double,double,char,char,char,long,long,double,long,char,char,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_KoSwap_t (double *,double *,double *,double *,double *,double *,long,double,double,double,char,char,char,long,long,double,long,char,long,double,double,long,char,char,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_KoOptionVarRebate_t (double *,double *,double *,long,double,double,char,char,int,int,FIX3_TREE_DATA *);
int     Fix3_KoFwdSwap_t(double *,double *,double *,double *,double *,double *,double *,long,double,double,double,char,char,char,char,long,long,double,double,long,char,long,double,double,long,char,char,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* ladder.c */
int  Fix3_LadderInterpolate_t(double **, double const*, double const*, double, double, double, char, double, int, double const*, double const* PrevStates, int, FIX3_TREE_DATA const* tree_data);

int     Fix3_LadderSwap_t (double **,long,long,double *, double *, double *,double,double,double,double,char,char,int,double *,double *,int,FIX3_TREE_DATA *);
int     Fix3_LadderStep_t (double *,double *,double,double,double,double,double,char,int,FIX3_TREE_DATA *);

/* lattice.c */
int     Fix3_Lattice_Classic (FIX3_DEV_DATA *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *);

/* lattice_timedep.c */
int     Fix3_Lattice_TimeDep (FIX3_DEV_DATA *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *);

/* lattice_smd.c */
int     Fix3_Lattice_Smd     (FIX3_DEV_DATA *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *);

/* lattice_tmx.c */
int     Fix3_Lattice_Tmx (FIX3_DEV_DATA *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_NmrToCcy    (FIX3_DEV_DATA *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_CcyToNmr    (FIX3_DEV_DATA *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *);

/* lattice_e2q.c */
int     Fix3_Lattice_E2Q     (FIX3_DEV_DATA *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *);

/* multiserial.cpp */
int     Fix3_Multiserial_t(double **, int, double *, double, double, char, int, double **, double, double, int, FIX3_TREE_DATA *);

/* numer.c */
int     Nmr_Alloc(MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Nmr_Schedule(long,long,MKTVOL_DATA *);
int     Nmr_Calc(T_CURVE const*,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Nmr_Interp(double *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     TMX_StateVar(double *, double *, long, double, long, long, char, char, MKTVOL_DATA *, T_CURVE *);

/* opt.c */
int     Fix3_Option_t (double *,double *,double,double,long,int,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_OptionPlus_t (double *,double *,double *,double *,double *,double,double,long,int,int,double *,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_OptionStat_t (double *,double *,double *,double *,double *,double,double,long,int,int,double *,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* paryield.c */
int     Fix3_Par_Yield_t (double *,int,double **,long *,long,long,long,int,char,char,double,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* printdebug.c */
void    Fix3_Print_FIX3_TREE_DATA(FILE *, FIX3_TREE_DATA *);

/* slice.c */
int     Fix3_Dump_Slice (double const*,int,FIX3_TREE_DATA const*);
int     Fix3_Set_Slice (double *,double,int,FIX3_TREE_DATA const*);
int     Fix3_Copy_Slice (double *,double *,int,FIX3_TREE_DATA const*);
int     Fix3_LCombTwoSlices(double *,double *,double,double *,double,int,FIX3_TREE_DATA const*);
int     Fix3_AddTwoSlices(double *,double *,double *,int,FIX3_TREE_DATA const*);
int     Fix3_MultiplyTwoSlices(double *,double *,double *,int,FIX3_TREE_DATA const*);
int     Fix3_MultTwoSlicesAddAll(double *,double *,double *,int,FIX3_TREE_DATA const*);
int     Fix3_AddScalar(double *,double,int,FIX3_TREE_DATA const*);
int     Fix3_MultiplyScalar(double *,double,int,FIX3_TREE_DATA const*);
int     Fix3_MaxMinOnSlice(double *,double,double,int,FIX3_TREE_DATA const*);
int     Fix3_MaxOnSlice(double *,double,int,FIX3_TREE_DATA const*);
int     Fix3_MinOnSlice(double *,double,int,FIX3_TREE_DATA const*);
int     Fix3_SmoothMinOnSlice(double *,double,int,int,FIX3_TREE_DATA const*);
int     Fix3_Node_Offset (int,int,int,int, FIX3_TREE_DATA const*);
int     Fix3_Init_Slice (double *,double,FIX3_TREE_DATA const*);
int     Fix3_MaxOfTwoSlices (double *,double *,double *,int,FIX3_TREE_DATA const*);
int     Fix3_lSmoothMaxOfTwoSlices (double *,double *,double *,int,int,FIX3_TREE_DATA const*);
int     Fix3_MinOfTwoSlices (double *,double *,double *,int,FIX3_TREE_DATA const*);
int     Fix3_SmoothMinOfTwoSlices (double *,double *,double *,int,int,FIX3_TREE_DATA const*);
int     Fix3_AccrueSlice(double *, double, int, int, FIX3_TREE_DATA const*);
int     Fix3_DivideTwoSlices (double *,double *,double *,int,FIX3_TREE_DATA const*);
int     Fix3_StepUp(double *Slice, double const* Index, double Up, double Barrier, int Smooth, int t, FIX3_TREE_DATA const* tree);
int     Fix3_KOIndicator(double *indicator, double const* index, double low, double high, char io, int smooth, int t, FIX3_TREE_DATA const* tree);
int     Fix3_ExIndicator(double *indicator, double const* index, double barrier, char ba, int smooth, int t, FIX3_TREE_DATA const* tree);
int     Fix3_Expectation(double* result, double const* x, double y, double const* p, int t, FIX3_TREE_DATA const* tree_data);
int     Fix3_LogicalOr(double* result, double const* x, double const* y, int t, FIX3_TREE_DATA const* tree_data);
int     Fix3_AddScalar2(double *result, double const* arg, double scalar, int t, FIX3_TREE_DATA const* tree);
int     Fix3_SliceTimesScalar (double *,double *,double,int,FIX3_TREE_DATA const*);

/* statevar.c */
double  Fix3_StateVar_Interp(int, int, double *, double **, double[MAXNBSTATES][3], double);
int     Fix3_StateVar_Generate(double **, int, double, double, char);

/* stdinput.c */
int     Fix3_Param_Input (MKTVOL_DATA *,FIX3_TREE_DATA *,int,char[6][MAXBUFF],char const*);
int     Fix3_Param_Input_Classic (MKTVOL_DATA *, FIX3_TREE_DATA *, char const*, char const* );            
int     Fix3_Param_Check_Original(int,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Param_Check_Classic(int, MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_MktVolAndModel_Input (MKTVOL_DATA *,FIX3_TREE_DATA *,T_CURVE *,char[2][8],char[6][MAXBUFF]);
int     Fix3_MktAndModel_Input (MKTVOL_DATA *,FIX3_TREE_DATA *,T_CURVE *,char[2][8],char[6][MAXBUFF]);
int     Fix3_Model_Input_Classic (MKTVOL_DATA *, FIX3_TREE_DATA*, char const* );    
int     Fix3_Num_Input (MKTVOL_DATA *,FIX3_TREE_DATA *, char [6][MAXBUFF], char const* ); 
int     Fix3_Num_Input_New (MKTVOL_DATA *,FIX3_TREE_DATA *, char const*); 
int     Fix3_Env_Logging (char[2][8],MKTVOL_DATA *,FIX3_TREE_DATA *,T_CURVE *);               
int     Fix3_ZeroCurve_Input(T_CURVE *, MKTVOL_DATA *, FIX3_TREE_DATA *); 
                           
/* stdinput_timedep.c */
int     Fix3_Param_Input_TimeDep (MKTVOL_DATA *,FIX3_TREE_DATA *, char const*, char const*);
int     Fix3_Param_Check_TimeDep (int, MKTVOL_DATA *, FIX3_TREE_DATA *);
int     Fix3_Model_Input_TimeDep (MKTVOL_DATA *, FIX3_TREE_DATA*, char const* );                           

/* stdinput_smd.c */
int     Fix3_Param_Check_Smd (int,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Param_Input_Smd (MKTVOL_DATA *,FIX3_TREE_DATA *, char const*, char const*);
int     Fix3_Model_Input_Smd (MKTVOL_DATA *, FIX3_TREE_DATA *, char const* );

/* stdinput_tmx.c */
int     Fix3_Param_Input_Tmx_Old (MKTVOL_DATA *,FIX3_TREE_DATA *,int,char[6][MAXBUFF],char const*,char const *);
int     Fix3_Param_Check_Tmx (int,MKTVOL_DATA *,FIX3_TREE_DATA *); 
int     Fix3_MktVol_Input_W_Tmx (MKTVOL_DATA *,char *,T_CURVE const*,char[MAXBUFF],char *,char *,char *,char *,char *);
int     Convert_VoV(MKTVOL_DATA *,T_CURVE const*);
int     Fix3_Param_Input_Tmx (   MKTVOL_DATA *,FIX3_TREE_DATA *, char const*, char const*); 
int     Fix3_Model_Input_Tmx (MKTVOL_DATA *, FIX3_TREE_DATA *,  char const*);
                  
/* stdinput_e2q.c */
int     Fix3_Param_Check_E2Q(int,MKTVOL_DATA *,FIX3_TREE_DATA *);
int     Fix3_Param_Input_E2Q (MKTVOL_DATA *,FIX3_TREE_DATA *,char const*, char const*);
int     Fix3_Model_Input_E2Q (MKTVOL_DATA *,FIX3_TREE_DATA *, char const* );


/* sticky.c */
int     Fix3_Sticky_t(double**,long,double*,double*,double,double,double,double,double,double,double,double,double,double,double,double,int,int,int,int,int,double,double,double**,int,int,int,FIX3_DEV_DATA*,FIX3_TREE_DATA*); 
int     Fix3_StickySwap_t(double **,long,long,double *,double *,double *,double *,double,double,double,double,double,double,double,double,double,double,double,double,double,int,int,int,int,int,double,double,double **,int,FIX3_TREE_DATA *);
int     Fix3_VolBond_t(double **,long,double *,double *,double,double,double,double,double,double,int,int,double,double,double **,int,FIX3_TREE_DATA *);
int     Fix3_VolSwap_t(double **,long,long,double *,double *,double *,double *,double,double,double,double,double,double,double,double,double,double,double,int,int,double,double,double **,int,FIX3_TREE_DATA *);

/* swap.c */
int     Fix3_ParSwap_t (double *,double *,double *,double,long,double,char,long,double,long,double,char,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_AmortSwap_t (double *,double *,double *,double *,double *,double *,long,double,double,double,long,char,long,double,double,long,char,double,char,long,double,char,char,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_Swaplet_t (double *,double *,double *,long,double,double,char,double,char,double,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA*);
int     Fix3_FwdSwap_t (double *,double *,double *,double *,double *,double *,double *,long,double,double,double,long,char,long,double,double,long,char,char,long,double,char,char,char,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_FwdSwapFlows_t(double *,CLAIM_BANK *,CLAIM_BANK *,double *,double,double,double,long,long,long,char,char,double *,double,double *,double,double *,double,long,long,long,long,char,char,char,double,long,int,FIX3_TREE_DATA *);

/* swapvol.c */
double	Fix3_ExpDecay (double,double);
int     Fix3_BFactor_Classic (double*, long ,long, char, char, MKTVOL_DATA* , T_CURVE const*);
int     Fix3_SpotVol_Classic (MKTVOL_DATA*, T_CURVE const*);
int     Fix3_Interp_SpotVol_Classic (FIX3_TREE_DATA*, MKTVOL_DATA*);
int     Fix3_IndexVol_Classic (double *, long, long, char, char,MKTVOL_DATA*,T_CURVE const*);
int     Fix3_Filtered_SpotVol_Classic (MKTVOL_DATA *, T_CURVE const*);
int     Fix3_GenIndexVol_Classic (double *, long, long, long, long, char, char, MKTVOL_DATA *, T_CURVE const* );     

/* swapvol_timedep.c */
int     Fix3_BFactor_TimeDep (double*, long ,long, char, char, MKTVOL_DATA* , T_CURVE const*);
int     Fix3_SpotVol_TimeDep (MKTVOL_DATA*, T_CURVE const*);
int     Fix3_Interp_SpotVol_TimeDep (FIX3_TREE_DATA*, MKTVOL_DATA*);
int     Fix3_IndexVol_TimeDep (double *, long, long, char, char,MKTVOL_DATA*,T_CURVE const*);
int     Fix3_Interp_SpotVol_TimeDepOld (FIX3_TREE_DATA*, MKTVOL_DATA*);
int     Fix3_Filtered_SpotVol_TimeDep(MKTVOL_DATA *, T_CURVE const*);
int     Fix3_GenIndexVol_TimeDep (double *, long, long, long, long, char, char, MKTVOL_DATA *, T_CURVE const* );     
int     QInterp ( double ,double *, long  ,long*, double *, int  );

/* swapvol_smd.c */
double  Smd_MeanDrift (double, double, double, double);
double  Fix3_Tanh (double);
double  Fix3_ExpInt (double,double,double,double);
int     Fix3_SpotVol_Smd (MKTVOL_DATA*, T_CURVE const*);
int     Fix3_Interp_SpotVol_Smd (FIX3_TREE_DATA*, MKTVOL_DATA*);
int     Fix3_IndexVol_Smd (double *, long, long, char, char,MKTVOL_DATA*,T_CURVE const*);
int     Smd_SpotVol_2F_TimeDep (double [MAXNBDATE],long ,int ,long *,double *,int *,char ,char ,long *,long *,double ,double ,double ,double ,double ,double ,int ,double *,double [MAXNBDATE],double [MAXNBDATE],double [MAXNBDATE],int ,int ,MKTVOL_DATA *,T_CURVE const*);          
int     Smd_IndexVol_2F_TimeDep (double *,long ,long ,char ,char ,int ,long ,long *,double [MAXNBDATE],double ,double ,double ,double ,double ,double ,int ,double *,double *,double *,double *,MKTVOL_DATA *,T_CURVE const*);       

/* swapvol_tmx.c */
int     Fix3_BFactor_Tmx (double*, long ,long, char, char, MKTVOL_DATA* , T_CURVE const*);
int     Fix3_SpotVol_Tmx (MKTVOL_DATA*, T_CURVE const*);
int     Fix3_Interp_SpotVol_Tmx (FIX3_TREE_DATA*, MKTVOL_DATA*);
int     Fix3_IndexVol_Tmx (double *, long, long, char, char,MKTVOL_DATA*,T_CURVE const*);
int     Fix3_Filtered_SpotVol_Tmx (MKTVOL_DATA *, T_CURVE const*);
int     Fix3_GenIndexVol_Tmx (double *, long, long, long, long, char, char, MKTVOL_DATA *, T_CURVE const* );     
int     Fix3_IndexLimits_Tmx (double* ,double*, long, long, long, int, char, char,  MKTVOL_DATA*, T_CURVE const*);

/* swapvol_e2q.c */
int     Fix3_BFactor_E2Q (double*, long ,long, char, char, MKTVOL_DATA* , T_CURVE const*);
int     Fix3_SpotVol_E2Q (MKTVOL_DATA*, T_CURVE const*);
int     Fix3_Interp_SpotVol_E2Q (FIX3_TREE_DATA*, MKTVOL_DATA*);
int     Fix3_IndexVol_E2Q (double *, long, long, char, char,MKTVOL_DATA*,T_CURVE const*);
int     Fix3_Filtered_SpotVol_E2Q(MKTVOL_DATA *, T_CURVE const*);
int     Fix3_GenIndexVol_E2Q (double*, long, long, long, long, char, char, MKTVOL_DATA*, T_CURVE const*);
	    
/* statvar_mc.c */
int     Fix3_DoubleArrayFloorIdx(double const*, int, double);
int     Fix3_DoubleVectSort(double *, int);
int     Fix3_DoubleQuadraticInterp(double const *, double const *, int, double, double *);
int     Fix3_Initialize_IR_SIM (IR_SIM *, long, int);
int     Fix3_IndexCovar (double *, double *, double *, long,long,char,char,long,long, char,char,MKTVOL_DATA *mktvol_data,T_CURVE const*);
int     Fix3_LCoef (double[3][3], long, MKTVOL_DATA *);      
int     Fix3_IndexLimits_Classic (double* ,double*, long, long, long, int, char, char,  MKTVOL_DATA*, T_CURVE const*);
int     Fix3_SwapYield_MC (double[MAXNBSTATES],double*, double*, int, long ,long, char, char, IR_SIM  *, MKTVOL_DATA *,T_CURVE const*);            
/* time.c */
int     Fix3_Time_Line (long,int,CRIT_DATE *,char,FIX3_TREE_DATA *);
int     Fix3_Time_Line_Nmr (long,int*,CRIT_DATE **,char,MKTVOL_DATA *,FIX3_TREE_DATA *);

/* tree.c */
int     Fix3_Forward_Curve (double *,double *,double *,double *,T_CURVE const*,long const*,int);
int     Fix3_Tree_Limits_Classic (int *,int *,int *,int *,int **,int **,int ***,int ***,int, FIX3_TREE_DATA *tree_data, MKTVOL_DATA *mktvol_data);
int     Fix3_Build_Tree (T_CURVE const*,MKTVOL_DATA*,FIX3_TREE_DATA *);

/* tree_smd.c */
int     Fix3_Tree_Limits_Smd (int *,int *,int *,int *,int **,int **,int ***,int ***,int, FIX3_TREE_DATA *tree_data, MKTVOL_DATA *mktvol_data);

/* tree_timedep.c */
int     Fix3_Tree_Limits_TimeDep (int *,int *,int *,int *,int **,int **,int ***,int ***,int, FIX3_TREE_DATA *tree_data, MKTVOL_DATA *mktvol_data);

/* util.c */
double  Fix3_GetIndexStep (double const *,int,int,int,int,int,FIX3_TREE_DATA const *);

/* zerobank.c */
int     Fix3_ZbkUpdate(CLAIM_BANK *,int,long,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
double *Fix3_ZbkReadZero(CLAIM_BANK const*,long,int,long,int,FIX3_TREE_DATA const*);
int     Fix3_ZbkParYield_t(double *,double *,CLAIM_BANK *,long,long,int,char,char,double,int,FIX3_TREE_DATA *);
int     Fix3_ZbkAnnuity_t(double *,CLAIM_BANK *,long,long,int,char,char,int,FIX3_TREE_DATA *);

/* zerobank_plus.cpp */
int     Fix3_ZbkDLFromIdx_Plus(long, int, long *, long *, char, int, char, int *, long **, int *, long **);
	
/* zerobinary.c */
int     Fix3_ZeroBinary_t(double **, long, double *, double, char, int, double **, double, double, int, FIX3_TREE_DATA *);

/* zeros.c */
int     Fix3_Zero_t (double *,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);
int     Fix3_Zero_Bank (double **,long *,int *,int,long,long,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* zip.c */
int  Fix3_Zip_t(double *,long,long,double,double *,double,double,double,double *,double *,int,int,int,int,FIX3_DEV_DATA *,FIX3_TREE_DATA *);

/* dllutil.c */
int Fix3_PackMktVolAndTreeData(  
        MKTVOL_DATA*    VolData,        
        FIX3_TREE_DATA* Tree,          
        long const*     BaseDatesL,   
        long const*     VolDatesL,   
        long const*     VolMatsL,   
        double const*   VolsL,     
        /* proper constness requires too many changes in product */
        char **    VolFreqL, 
        char **    VolDCCL,       
        char **    VolTypeL,     
        char **    VolCalibL,   
        double const*   MrParamsL,  
        double const*   SmileParamsL, 
        long const*     TreeParamsL);


int Fix3_UnPackMktVolAndTreeData(  
        MKTVOL_DATA const*      VolData,
        FIX3_TREE_DATA const*   Tree,  
        long*                   BaseDatesL,
        long*                   VolDatesL,
        long*                   VolMatsL,
        double*                 VolsL,      
        char**                  VolFreqL,  
        char**                  VolDCCL,  
        char**                  VolTypeL,
        char**                  VolCalibL,   
        double*                 MrParamsL,  
        double*                 SmileParamsL, 
        long*                   TreeParamsL);
int     Fix3_PackMktVolAndModelData   (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,long *,long *,double *,char **,long *,double *,long *,double *,long *);
int     Fix3_UnPackMktVolAndModelData (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,long *,long *,double *,char **,long *,double *,long *,double *,long *);
int     Fix3_PackModelData_Classic    (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int     Fix3_PackModelData_TimeDep    (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int     Fix3_PackModelData_Smd        (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int     Fix3_PackModelData_Tmx        (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int     Fix3_PackModelData_E2Q        (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int     Fix3_UnPackModelData_Classic  (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int     Fix3_UnPackModelData_TimeDep  (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int     Fix3_UnPackModelData_Smd      (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int     Fix3_UnPackModelData_Tmx      (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);
int     Fix3_UnPackModelData_E2Q      (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *);

/* concatenate FIX3 style deal, market and model params data into one file */
void Fix3_ConcatenateFix3Input(char const* output, char const* deal);


/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif



#ifdef __cplusplus
// Model parameters input - stream version

#include <iosfwd>

void Fix3_ParamInput(MKTVOL_DATA *,FIX3_TREE_DATA *,int,char[6][MAXBUFF], std::istream& is);
void Fix3_ParamInput(MKTVOL_DATA *,FIX3_TREE_DATA *,int,char[6][MAXBUFF]);

#endif



#endif /* _fix123head_h */
