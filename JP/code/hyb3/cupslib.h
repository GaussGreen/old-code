/****************************************************************************/
/*      Function templates for library.                                	    */
/****************************************************************************/
/*      CUPSLIB.h                                                             */
/****************************************************************************/


/*
$Header$
*/

#include "eslhead.h"

#ifndef _CUPSLIB_H
#define _CUPSLIB_H


#include "cupsmodl.h"

#ifdef __cplusplus
extern "C" {
#endif


/* ALLOC.C */
/* SwapRate data */
void    Hyb3_SwapRate_Init    (SWAPRATE_DATA *);
int     Hyb3_SwapRate_Alloc   (SWAPRATE_DATA *, HYB3_TREE_DATA *);
void    Hyb3_SwapRate_Free    (SWAPRATE_DATA *, HYB3_TREE_DATA *);
/* TreeSim data */
int     Hyb3_TreeSim_Init     (TREESIM_DATA *, int);
int     Hyb3_TreeSim_Alloc    (TREESIM_DATA *, HYB3_TREE_DATA *);
void    Hyb3_TreeSim_FreeSR   (TREESIM_DATA *, HYB3_TREE_DATA *);
void    Hyb3_TreeSim_Free     (TREESIM_DATA *, HYB3_TREE_DATA *);
/* Tree timeline   */
void    Hyb3_Tree_Init (HYB3_TREE_DATA *);            
int     Hyb3_Tree_Alloc (HYB3_TREE_DATA *);
int     Hyb3_Tree_Free (HYB3_TREE_DATA *);
/* Tree variables  */
void    Hyb3_Dev_Init (HYB3_DEV_DATA *);
int     Hyb3_Dev_Alloc (HYB3_DEV_DATA *,HYB3_TREE_DATA *); 
int     Hyb3_Dev_Free (HYB3_DEV_DATA *,HYB3_TREE_DATA *);
/* Slice functions */
TSLICE  Hyb3_Alloc_Slice (HYB3_TREE_DATA *,int);      
void    Hyb3_Free_Slice (TSLICE,HYB3_TREE_DATA *,int);
/* Dimension function */
int     Hyb3_Tree_Dim(HYB3_TREE_DATA *);



/* BOND.C */
int     Hyb3_Bond_Price(TSLICE,double,long,double,long,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_FwdBond_t(TSLICE,TSLICE,TSLICE,double,long,double,long,char,double,long,char,long,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *); 

/* CALLTURBO.C */
int    Hyb3_Callturbo_t(TSLICE,TSLICE,TSLICE,int,int,double,int,TSLICE,char,char,int,double,int,double,TSLICE,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *) ;

/* CALLZERO.C */
int     Hyb3_Callzero_t(TSLICE,TSLICE,TSLICE,TSLICE,TSLICE,int,int,int,double,double,double,char,char,char,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);

/* CET.C */
int     Hyb3_Cet_Main(int,int,T_CURVE [2][3],MKTVOL_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Cet_Manager(int,MKTVOL_DATA *,HYB3_TREE_DATA *,MKTVOL_DATA *,HYB3_TREE_DATA *); 
int     Hyb3_Cet_Schedule(long, CET_OUT_DATA *,MKTVOL_DATA *,HYB3_TREE_DATA *); 
int     Hyb3_Cet_Calc(MKTVOL_DATA *,HYB3_TREE_DATA *,CET_OUT_DATA *);
int     Hyb3_Cet_Print_Vol_File(int,CET_OUT_DATA *,MKTVOL_DATA *,int,int);

/* claimbank.c */
void    Hyb3_CbkInit(CLAIM_BANK *);
int     Hyb3_CbkAlloc(CLAIM_BANK *,int,HYB3_TREE_DATA *, int);
int     Hyb3_CbkFree(CLAIM_BANK *,HYB3_TREE_DATA *, int);
int     Hyb3_CbkSizeFromDL(long *,int,long *,int);
int     Hyb3_CbkCalcSize(int,CRIT_DATE *,int,int);
int     Hyb3_CbkDev(CLAIM_BANK *,long,int,int,int, int, HYB3_DEV_DATA *,HYB3_TREE_DATA *);
TSLICE  Hyb3_CbkPopSlice(CLAIM_BANK *);
int     Hyb3_CbkPushSlice(CLAIM_BANK *, TSLICE,long,long);
TSLICE  Hyb3_CbkReadSlice(CLAIM_BANK *,long);
int     Hyb3_CbkGetOffset(CLAIM_BANK *,long,int);

/* CORRELATIONS.C */
int     Hyb3_AssetIrCorr(double *,double,double *,int,double *,double *,
                         int,long,long,char,char,long,T_CURVE *);
int     Hyb3_IR_IR_Corr( double *,double,double *,double *,int,int,double *,
                         double *,double *,double *,int,long,long,long,long,
                         char,char,char,char,long,long,T_CURVE *,T_CURVE *);

/* DERIVE.C */
double  Hyb3_Delta (double *,double *);
double  Hyb3_Gamma (double *,double *);
double  Hyb3_Delta1 (double *,double *,char);
double  Hyb3_Gamma1 (double *,double *,char);

/* DEV.C */
int     Hyb3_Dev(TSLICE,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Ev (double*, int,int,int,HYB3_DEV_DATA   *,HYB3_TREE_DATA  *);

/* DRIFT.C */
int     Hyb3_SimpleFwdCurve(double *,T_CURVE const*,long *,int );
int     Hyb3_Forward_Curve(double *,double *,double *,T_CURVE const*,long const*,int);
int     Hyb3_CUP_Drift(double *[],double ** ,int ,double *,double *[],double *,int);
int     Hyb3_Find_Drift_1D(double *,double *,double *,double,double,double,double,double *,int *,int *,int *,int *,int,double *,double *,HYB3_TREE_DATA *);
int     Hyb3_Find_Drift_2D(double *,double *,double *,double,double,double,double[3],
                             double **,int *,int *,int *,int *,int **,int **,int **,
                             int **,int,double *,double *,HYB3_TREE_DATA *);
int     Hyb3_CUPSAdjust_FXSmile(HYB3_TREE_DATA *, MKTVOL_DATA *, FX_DATA *);
int     Hyb3_FX_MomentMatching(MKTVOL_DATA *, HYB3_TREE_DATA *, double *);


/* FLOATER.C */
int     Hyb3_Floater_Price (TSLICE,int,TSLICE,TSLICE,long,double,double,char,char,double,double,long,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_FwdFloater_t(TSLICE,TSLICE,int,TSLICE,int,TSLICE,long,double,double,long,char,double,char,long,char,char,long,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Collaret_t(TSLICE,TSLICE,TSLICE,int,double,double,double,char,double,char,double,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Floater_Pmt (TSLICE,int,TSLICE,TSLICE,double,double,double,double,int,int,HYB3_TREE_DATA *);

/* FXTOOLS.C */

int     Hyb3_Gfunc(double *,double,double,double,double);
int     Hyb3_GDashfunc(double*,double,double,double,double);
int     Hyb3_GDriftfunc(double*,double*,double,double,double,double);
int     Hyb3_Kfunc(double *,double,double,double,double);
int     Hyb3_KDashTimesX(double *,double,double,double,double);
int     Hyb3_UpdateSmileCache(HYB3_FXSMILE_CACHE *,int,long *,long *,int,int,double,double *,double *,int,double,double,double);
int     Hyb3_BuildFXSmileCache(HYB3_TREE_DATA *,double *,double *);

int     Hyb3_FillGrid_2d(HYB3_TREE_DATA *,TSLICE,TSLICE,TSLICE,TSLICE,double *,double *,double *,int, int);
int     Hyb3_FillGrid_3d(HYB3_TREE_DATA *,TSLICE,TSLICE,TSLICE,TSLICE,double *,double *,double *,int, int);

int     Hyb3_FillGridEqFwd_2d(HYB3_TREE_DATA *,HYB3_DEV_DATA *,int);
int     Hyb3_FillGridEqFwd_3d(HYB3_TREE_DATA *,HYB3_DEV_DATA *,int);

/* FXVOL.C */
int     Hyb3_Get_TreeFxSpotVol(double *,long *,long,double *,long *,int,double *,long *,int,double *,double *,double *,double *,double *,double *,double *,double*,double *,double *,double*,int,int,double ,double *,double*,int,int,long,int,double*,double *,double *);
int     Hyb3_FX_Vol(double *,double *, int, int, double *,double *, double *, double *, double * , double*[2],double*,double*,double[3],double[3],long *,long*,double*,double*,double*,double[3],double[3],long*,int);               
int     Hyb3_Forward_FX(double*,int,double *[2][3],int*,double);
int     Hyb3_FXMidNode(double  *,int ,double  *,double  *, int   *,double *,double, double  *,double  *,double  *);
int     Hyb3_FX_Fwd (double *,double,T_CURVE *,T_CURVE *,long,long);
int     Hyb3_Get_FxVol(T_CURVE *,T_CURVE *,MKTVOL_DATA *,MKTVOL_DATA *,FX_DATA *,int ,double *,int ,long *,double *);
int     Hyb3_MultiFac_FxVol2(double  *,double  *,double  *,double  *,double  *,double  *,double  *,double  *,double  *,int,int,double  *,double  *,double  *,double  *,double  *,double  *,double  *,double  *,double  *,long    *,long    *,int,long    *,int);
int     Hyb3_MultiFac_Spot_FxVol2(double  *,double  *,double  *,double  *,double  *,double  *,double  *,double  *,double  *,double  *,double  *,double  *,int,long    *,long    *,double  *,int, int,double,double  *,double  *,int,int,double *,double*,double*,long*,int);   
int     Hyb3_MultiFac_Spot_FxVol(double  *, double  *, double  *, double  *, double  *, double  *, double  *,double  *,double  *,double  *, double  *, double  *,int , long    *, long    *, double  *, int, int , double  , double  *, double  *, int, int , long *, int );     
int     Hyb3_MultiFac_FxVol(double  *,double  *,double  *, double  *, double  *, double  *,double  *,double  *,double  *,int, int , double  *,double  *, double  *,double  *,double  *,double  *,long *,int ,long *,int);
int     Hyb3_DiscretiseInput(double *,double *,int ,int ,long *,int,long *);

int     Hyb3_For_MultiFac_FxVol(double  *, double  *, double  *, double  *,double  *, double  *,double  *,double  *,  double  *,   int ,int , double  *,  double  *,  double  *,  double  *,  double  *,double  *, long , long,long *,int   );     
int     Hyb3_For_FX_Vol(double *, double *, double  *[6],double *[2], double *, double *, double [3], double [3], long , long , double *, double *, double *,double [3], double [3], long *,	int );              
int     Hyb3_FX_Fwd_Ratio (double *,T_CURVE *,T_CURVE *,long,long);
int     Hyb3_ForwardFX(double *,long,double,long,T_CURVE const*, T_CURVE const*);
int     Hyb3_Get_FxVol2 (T_CURVE *,T_CURVE *,MKTVOL_DATA *,MKTVOL_DATA *,FX_DATA *,int,long,long *,double *,int,long *,double *);
int     Hyb3_Get_FxSpotVol(T_CURVE *, T_CURVE *,MKTVOL_DATA *,MKTVOL_DATA *,FX_DATA *,int,double *);  

int     Hyb3_CalcPerturb(double *,long,double *,double *,double *,double *,double *);
int     Hyb3_CalcVolA_xx(double *,long,double *,double *,double *,double *,double *);

/* KOOPT.C */
int     Hyb3_KoOption_t (TSLICE,int,TSLICE,long,double,double,int,TSLICE,char,char,int,int,HYB3_TREE_DATA * );
int     Hyb3_KiOption_t (TSLICE,TSLICE,int,TSLICE,long,double,double,char,char,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA  *); 
int     Hyb3_KoOption_Am_t (TSLICE,int,TSLICE,long,double,double,int,TSLICE,char,char,int,int,HYB3_TREE_DATA *, int );
int     Hyb3_KiOption_Am_t (TSLICE,TSLICE,int,TSLICE,long,double,double,char,char,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA  *, int ); 
int     Hyb3_Trigger_t(TSLICE,TSLICE,int,TSLICE,double,double,long,int,double,double,char,char,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);

int     Hyb3_MultiKo_t (TSLICE,TSLICE,int,TSLICE,long,double,double,int,TSLICE,char,char,int,int,HYB3_TREE_DATA *);
int     Hyb3_MultiKi_t (TSLICE,TSLICE,TSLICE,int,TSLICE,long,double,double,char,char,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Multi_Smooth(double *,double,double,double,double,double,double,double);

int     Hyb3_KOIndicator( double*,double const*,double,double,char,int,int,HYB3_TREE_DATA const* );

/* LATTICE.C */
int     Hyb3_Lattice(HYB3_DEV_DATA *, int, int, MKTVOL_DATA  *, HYB3_TREE_DATA *);

/* MC.C */
int     Hyb3_MC_CheckDate(HYB3_TREE_DATA *,int,TREESIM_DATA *,int,int);
int     Hyb3_MC_GetSliceValue(double *,HYB3_TREE_DATA *,TREESIM_DATA *,TSLICE,int,int,int);
void    Hyb3_MC_StartTreeSim(TREESIM_DATA *);
int     Hyb3_MC_UpdateTreeSim(int,int,TREESIM_DATA *,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_MC_ProcessSamples(TREESIM_DATA *);
int     Hyb3_MC_CheckSamples(TREESIM_DATA *);
int     IsGreater(const void *, const void *);

/* OPT.C */
int     Hyb3_Option_t(TSLICE,TSLICE,double,double,long,int,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_OptionPlus_t(TSLICE   Opt,double*,double*,double*,TSLICE,double,double,long,int,int,double*,int,int,int,int,HYB3_DEV_DATA*,HYB3_TREE_DATA*);
int     Hyb3_Exch_Option_t(TSLICE,TSLICE,TSLICE,double,double,double,long,int,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Option_Series_t(TSLICE,TSLICE,double,double,long,int,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Option_Simple_t(TSLICE,double,double,int,int,int,HYB3_TREE_DATA*);

/* PARYIELD.C */
int     Hyb3_Par_YieldPlus_t(TSLICE,TSLICE,int,int,TSLICE *,long const*,long,long,long,int,char,char,double,int,HYB3_TREE_DATA const*);
int     Hyb3_Par_Yield_t(TSLICE,int,int,TSLICE *,long const*,long,long,long,int,char,char,double,int,HYB3_TREE_DATA const*);

/* SLICE.C */
int     Hyb3_Init_Slice (double *,int,double,HYB3_TREE_DATA *);
int     Hyb3_SetSlice(TSLICE,int,double,int,HYB3_TREE_DATA *);
int     Hyb3_CopySlice(TSLICE,TSLICE,int,int,HYB3_TREE_DATA *);
int     Hyb3_ExpandSlice(TSLICE,int,TSLICE,int,int,HYB3_TREE_DATA *);
int     Hyb3_AddTwoSlices(TSLICE,int,TSLICE,TSLICE,int,HYB3_TREE_DATA *);
int     Hyb3_SubtractTwoSlices(TSLICE,int,TSLICE,TSLICE,int,HYB3_TREE_DATA *);
int     Hyb3_MultiplyTwoSlices(TSLICE,int,TSLICE,TSLICE,int,HYB3_TREE_DATA *);
int     Hyb3_DivideTwoSlices(TSLICE,int,TSLICE,TSLICE ,int ,HYB3_TREE_DATA   *); 
int     Hyb3_LCombTwoSlices(TSLICE,int,TSLICE,double,TSLICE,double,int,HYB3_TREE_DATA *);
int     Hyb3_MultTwoSlicesAddAll(double *,int,double *,double *,int,HYB3_TREE_DATA *);
int     Hyb3_MaxTwoSlices(TSLICE,int,TSLICE,TSLICE,int,HYB3_TREE_DATA *);
int     Hyb3_MinTwoSlices(TSLICE,int,TSLICE,TSLICE,int,HYB3_TREE_DATA *);
int     Hyb3_AddScalar(TSLICE,int,double,int,HYB3_TREE_DATA *);
int     Hyb3_MultiplyScalar(TSLICE,int,double,int,HYB3_TREE_DATA *);
int     Hyb3_MaxMinOnSlice(TSLICE,int,double,double,int,HYB3_TREE_DATA *);
int     Hyb3_MaxOnSlice(TSLICE,int,double,int,HYB3_TREE_DATA *);
int     Hyb3_MinOnSlice(TSLICE,int,double,int,HYB3_TREE_DATA *);
int     Hyb3_CheckDimRebate(int);
int     Hyb3_CheckDimIndex (int);
double  Hyb3_GetValueAtNode(int,TSLICE,int,int,int,int,HYB3_TREE_DATA *);
int     Hyb3_SmoothStepUp(double*,double*,double,double,int,int,int,HYB3_TREE_DATA *);
int     Hyb3_Node_Offset (int,int,int,int,HYB3_TREE_DATA const*);
int     Hyb3_Slice_Dim(int);

/* STATEPRICES.C */
int     Hyb3_UpdateStatePrices3D(int ,MKTVOL_DATA  *,HYB3_TREE_DATA *,HYB3_DEV_DATA *,TSLICE ,TSLICE );
int     Hyb3_UpdateStatePrices2D(int ,MKTVOL_DATA  *,HYB3_TREE_DATA *,HYB3_DEV_DATA *,TSLICE ,TSLICE );

/* STDINPUT.C */
/* General */
int     Hyb3_FindAndSkipComLineOptional (FILE *,char const *,char const *, char const *);
int     Hyb3_Param_Input (MKTVOL_DATA *,char [MAXINDEX], HYB3_TREE_DATA *,int,char[4][MAXBUFF],char *);
int     Hyb3_Param_Check (int,MKTVOL_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Eq_Input_W_WithSmile(EQ_DATA *,long,char [MAXBUFF],char *,char *,char [MAXBUFF],char [MAXNBDATE][MAXBUFF]);
int     Hyb3_Eq_Check_W(EQ_DATA *);
int     Hyb3_Eq_Input_Dyn_WType4(EQ_DATA *,long,char *,char *,char *,char [MAXBUFF],char [MAXNBDATE][MAXBUFF]);
int     Hyb3_Eq_Check_Dyn_WType4(EQ_DATA *,long); 
int     Hyb3_Eq_Input_Dyn_WType6(EQ_DATA *,long,char *,char,char *,char [MAXBUFF],char [MAXNBDATE][MAXBUFF]);
int     Hyb3_Eq_Check_Dyn_WType6(EQ_DATA *,long,char);
int     Hyb3_Eq_Input_Sta_W (EQ_DATA *,long,char *);
int     Hyb3_Eq_Check_Sta_W (EQ_DATA *);
int     Hyb3_Fx_Input_W_WithSmile(FX_DATA *,char const[MAXBUFF], char const*,char const*,char const[MAXBUFF],char const[MAXNBDATE][MAXBUFF]);
int     Hyb3_Fx_Input_W(FX_DATA   *,char      [MAXBUFF],char      *);
int     Hyb3_Fx_Check_W(FX_DATA *);
int     Hyb3_Correl_Input_WType3(FX_DATA *,char [3][MAXBUFF],char *);
int     Hyb3_Correl_Check_WType3(FX_DATA *);
int     Hyb3_Correl_Input_WType4(EQ_DATA *,char [MAXBUFF],char *);
int     Hyb3_Correl_Check_WType4(EQ_DATA *);
int     Hyb3_Correl_Input_WType6(FX_DATA *,EQ_DATA *,char,char [6][MAXBUFF],char *);
int     Hyb3_Correl_Check_WType6(FX_DATA *,EQ_DATA *);
int     Hyb3_ModelOverwrites_Input_W(FILE *,char *,char[MAXBUFF],char [MAXBUFF],char [MAXNBDATE][MAXBUFF],char [3][MAXBUFF],char OwriteMParamF[5][MAXBUFF],char[5][MAXBUFF],HYB3_TREE_DATA*);
int     Hyb3_ModelOverwritesMultFactCups_Input_W(FILE *,char *,char [MAXBUFF],char [MAXBUFF],char [MAXNBDATE][MAXBUFF],char [3][MAXBUFF],char [5][MAXBUFF],char [5][MAXBUFF],HYB3_TREE_DATA *);

/* STOCK.C */
int     Hyb3_Forward_Stock(double *,long *,double *,long,double,int,double *,long *,char *,int,long *,long *,char,int,long *,double *,T_CURVE const*,long *,double *,double *,int);
int     Hyb3_Eq_Shift_Centre(double *,double *,double *,int);
int     Hyb3_Eq_Spot_Vol (double *,double *,double *,double *,double *,double *,double *,long *,long,int,long *,double *,double,double,int);
int	    Hyb3_Eq_Vol(double *,double *,double *,double *,double *,double [3],long *,double *,double *,double *,double [3],long *,int);
int    	Hyb3_Composite_Eq(double *,double *,double *,double **,long *,T_CURVE const*, T_CURVE const*,double *,double *,double *,int);
int     Hyb3_CUPSEqMidNode(double *,double *,double *,double *,double *,int);
int     Hyb3_Get_TreeEqSpotVol(double *,double *,long *,long,double *,long *,int,double *,long *,int,double *,double *,double *,double *,int,int,double,double *,int,long,int,double *,double *,double *);

/* SWAP.C */
int     Hyb3_FwdSwap_t(TSLICE,TSLICE,TSLICE,int,TSLICE,TSLICE,TSLICE,long,double,double,double,long,char,long,double,double,double,long,char,char,long,double,char,char,long,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);

/* SWAPVOL.C */
double  Hyb3_ExpDecay (double,double);
int     Hyb3_BFactor(double *,long,long,char,char,int,double *,T_CURVE const* crv);
int     Hyb3_SpotVol(double[6][MAXNBDATE],long,int,long *,double *,int *,char,char,long *,long *,double,double,double,int,double *,double *,double *,int,int,T_CURVE const*);
int     Hyb3_Interp_SpotVol(double **,long,int,long *,double[6][MAXNBDATE],int,double *,int,double *,double *,int);
int     Hyb3_IndexVol (double *,long,long,char,char,int,int,long,long *,double[6][MAXNBDATE],double,double,double,int,double *,T_CURVE const*);
int     Hyb3_IndexSpotVol(double *,char,char,long,long,int,double *,double *,T_CURVE const*);

/* TIME.C */
int     Hyb3_Time_Line (long,int,CRIT_DATE *,char,HYB3_TREE_DATA *);


/* TREE.C */
int     Hyb3_Build_Tree (int,T_CURVE [2][3],MKTVOL_DATA *,FX_DATA *,EQ_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Tree_Limits(int *,int *, int *, int *,int *,int *,int **,int **,int ***,int ***,int ****,int ****,double,double,int,int,int,double *[2][3],int *,double *[2][3],double *,double *,double *[6],double *,double *,double *,double *);
int     Hyb3_Tree_Limits_IR1(int *,int *,int *,int *,double *,double *, double *,int,double,double);
int     Hyb3_Tree_Limits_IR2(int *,int *,int *,int *,int **, int **, double **,int , double *, double *,int,double,double * Beta);
int     Hyb3_Free_TreeLimits_IR2(int *,int *,int **, int **, int , int);
int     Hyb3_Tree_Weights(double **,double *[5],int,int,double *[2][3],double *,double *);
int     Hyb3_Tree_Correlations(double **,int,long *,int,int,long [MAXNBDATE+1],double [3][MAXNBDATE+1],int,long [MAXNBDATE],double [3][MAXNBDATE+1]);
int     Hyb3_Transform_SmileParams(double,double,double,double *,double *,double *,char *);
int     Hyb3_Tree_SmileParams(int *,double*,double*,double*,double*,double*,double*,long*,int,double[MAXNBDATE],double[MAXNBDATE],double[MAXNBDATE],int,long[MAXNBDATE],long[MAXNBDATE],long[MAXNBDATE]);
int     Hyb3_Find_IR_StdDev(double*,double,double,double,double*[2],double,int,double*[2][3],double*,double*,long  *,double,double,double,double,double,double);


/* TURBO.C */
int   Hyb3_Turbo_Pmt(TSLICE,TSLICE,double,double,char,TSLICE,double,double,char,double,double,double,double,TSLICE,int,char,int,int,int,HYB3_DEV_DATA   *,HYB3_TREE_DATA  *);   
int   Hyb3_FiTurbo_Pmt(TSLICE,TSLICE,int,double,double,double,double,TSLICE,int,double,double,double,double,TSLICE,double,double,TSLICE,int,int,HYB3_TREE_DATA  *);
int   Hyb3_Turbo_Price(TSLICE,TSLICE,double,double,char,TSLICE,double,double,char,double,double,double,double,TSLICE,int,char,int,double,int,double,TSLICE,int,int,int,HYB3_DEV_DATA   *,HYB3_TREE_DATA  *);
int   Hyb3_FwdTurbo_t(TSLICE,TSLICE,TSLICE,double,double,char,TSLICE,double,double,char,double,double,double,double,TSLICE,int,long,char,char,char,int,double,char,char,long,int,int,int,HYB3_DEV_DATA  *,HYB3_TREE_DATA *) ;
int   Hyb3_BTurbo_Price(TSLICE,TSLICE,double,double,char,TSLICE,double,double,char,double,double,double,double,TSLICE,int,char,int,double,int,double,TSLICE,int,double *,double *,char,int,int,int,HYB3_DEV_DATA *, HYB3_TREE_DATA *);
int   Hyb3_FwdBTurbo_t(TSLICE,TSLICE,TSLICE,double,double,char,TSLICE,double,double,char,int,double,double,double,double,TSLICE,int,double *,double *,char,int,double,char,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int   Hyb3_Redeemer_Price(TSLICE,TSLICE,TSLICE,double,double,double,int,HYB3_TREE_DATA  *); 
int   Hyb3_BTurbo_Flows(TSLICE,TSLICE,double,double,TSLICE,double,double,double,double,double,double,TSLICE,int,char,int,double, TSLICE,int,double,TSLICE,TSLICE,TSLICE,int,double *,double *,char,int,int,int,HYB3_DEV_DATA *, HYB3_TREE_DATA *);
int   Hyb3_FwdBTurbo_Flows(TSLICE,TSLICE,TSLICE,double,double,TSLICE,double,double,int,double,double,double,double,int, double, int, double, TSLICE,TSLICE,int,double *,double *,char,int,double,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int   Hyb3_BTurbo_Flows_Pmt(TSLICE,TSLICE,double,double,TSLICE,double,double,double,double,double,double,TSLICE,double, TSLICE,double,TSLICE,TSLICE,TSLICE,int,double *,double *,char,int,HYB3_TREE_DATA *);

/* UTIL.C */
double  Hyb3_GetIndexStep (TSLICE,int,int,int,int,int,HYB3_TREE_DATA const*);
int     Hyb3_CashAnnuity_t(TSLICE,TSLICE,int,char,int,int,HYB3_TREE_DATA   *tree_data);
int     PrintTree( FILE *fp, HYB3_TREE_DATA , int PrintTreeSizeInfo );
int     PrintFXData( FILE *fp, FX_DATA );
int     PrintEQData( FILE *fp, EQ_DATA );
int     PrintTCurve( FILE *fp, T_CURVE );
int     PrintMktVolData( FILE *fp, MKTVOL_DATA );
int     PrintSlice_2d(HYB3_TREE_DATA *,TSLICE,int,FILE *);
int     PrintSlice_3d(HYB3_TREE_DATA *,TSLICE,int,FILE *);
int     PrintCriticalDates(FILE* stream, int nbCritDate, CRIT_DATE* critDate);
int     Hyb3_DoubleQuadraticInterp(double const *,double const *,int,double,double *);
int     Hyb3_Triangulation4D (double [4][4],int,double [4][4]);

/* zerobank.c */
int     Hyb3_ZbkUpdate(CLAIM_BANK *,int,long,long,int,int, int, int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_FXZbkUpdate(CLAIM_BANK *,int,long,long,int,int, int, int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
double *Hyb3_ZbkReadZero(CLAIM_BANK *,long,int,long, int, int,HYB3_TREE_DATA *);
int     Hyb3_ZbkParYield_t(double *,double *, int,CLAIM_BANK *,long,long,int,char,char,double,int,HYB3_TREE_DATA *);
int     Hyb3_ZbkAnnuity_t(double *,int,CLAIM_BANK *,long,long,int,char,char,int,HYB3_TREE_DATA *);

/* ZEROS.C */
int     Hyb3_Zero_t (TSLICE,long,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_FXZero_t (TSLICE,long,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);
int     Hyb3_Zero_Bank(TSLICE *,long *,int *,int,long,long,int,int,int,int,HYB3_DEV_DATA *,HYB3_TREE_DATA *);

/* STICKY.C */
int     Hyb3_StickyTurboSwap_t(TSLICE *,int,int,TSLICE,TSLICE,double,double,double,double,int,TSLICE,char,int,int,double,double,double**,int,HYB3_TREE_DATA*);

/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif /* CUPSLIB_H */

