/****************************************************************************/
/*      Function templates for library.                                     */
/****************************************************************************/
/*      FIX123HEAD.h                                                        */
/****************************************************************************/

/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/include/fix123head.h,v 1.75 2005/02/04 20:01:48 skuzniar Exp $
*/


/* Use safe muliple inclusion for Jerry Cohen's stream libraries */

#ifndef _fix123head_h
#define _fix123head_h

#include <stdio.h>
#include "fix123.h"


#ifdef __cplusplus
extern "C" {
#endif


/* alloc.c */
void    Tree_Init (TREE_DATA *);
int     Tree_Alloc(TREE_DATA *);
int     Tree_Free (TREE_DATA *);
void    Dev_Init  (DEV_DATA    *dev_data);
int     Dev_Alloc (DEV_DATA *,TREE_DATA *);
int     Dev_Free (DEV_DATA *,TREE_DATA *);
void    MktVol_Init(MKTVOL_DATA *);
void    Opt_Out_Data_Init(OPT_OUT_DATA  *ood);
double  *Alloc_Slice (TREE_DATA *);
int     Free_Slice (double *,TREE_DATA *);
void    *DR_Array (int,int,int);
void    *DR_Matrix (int,int,int,int,int);
void    *DR_Cube (int,int,int,int,int,int,int);
int     Free_DR_Array (void *,int,int,int);
int     Free_DR_Matrix (void *,int,int,int,int,int);
int     Free_DR_Cube (void *,int,int,int,int,int,int,int);

/* amortloan.c */
int     AnnAmortSwap_t (double **,double *,double *,double *,char,char,double,double,double,long,long,long,int,double,double,double **,int,int,int,DEV_DATA *,TREE_DATA *);

/* bas.c */
double  NDensity (double);
double  NormalH (double);
double  Normal_InvH (double);
double  Call_BSQ (double,double,double,double,double);
double  Put_BSQ (double,double,double,double,double);
double  Vega_BSQ (double,double,double,double,double);
double  Binary_BS (double,double,double,double,double,char);
double  Call_BS (double,double,double,double,double);
double  Put_BS (double,double,double,double,double);
double  D_Call_BS (double,double,double,double,double);
double  D_Put_BS (double,double,double,double,double);
double  Gamma_BS (double,double,double,double,double);
double  Vega_BS (double,double,double,double,double);
double  T_Call_BS (double,double,double,double,double);
double  T_Put_BS (double,double,double,double,double);
double  R_Call_BS (double,double,double,double,double);
double  R_Put_BS (double,double,double,double,double);
double  CCEquation_BS2Q (double,double,double,double,double);
int     ConvexityC_BS2Q (double,double,double,double,double *);
double  Option_BS2Q (double,double,double,double,char,double,double,double);
double  ImpVol_BS2Q (double,double,double,double,char,double,double,double,double);

/* bond.c */
int     Bond_t   (double *,double,long,double,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     FwdBond_t(double *,double *,double *,long,double,double,char,long,double,double,long,double,char,long,int,int,int,DEV_DATA *,TREE_DATA *);

/* callzero.c */
int     Callzero_t (double *,double *,double *,double *,double *,long,long,long,double,double,double,char,char,char,int,int,int,DEV_DATA *,TREE_DATA *);

/* cap.c */
int     Cap_t (double *,double *,double *,int,long,double,double,char,double,char,double,int,int,int,DEV_DATA *,TREE_DATA *);
int     Caplet_t (double *,double *,double *,int,long,double,double,char,double,char,double,int,int,int,DEV_DATA *,TREE_DATA *);
int     Collaret_t (double *,double *,double *,long,double,double,double,char,double,char,double,int,int,int,DEV_DATA *,TREE_DATA *);

/* cet.c */
int     Cet_Main(int,T_CURVE *,MKTVOL_DATA *,TREE_DATA *);
int     Cet_Manager(MKTVOL_DATA *,TREE_DATA *,MKTVOL_DATA *,TREE_DATA *); 
int     Cet_Schedule(long, CET_OUT_DATA *,MKTVOL_DATA *,TREE_DATA *,TREE_DATA *); 
int     Cet_Calc(MKTVOL_DATA *,TREE_DATA *,CET_OUT_DATA *);
void    Cet_Tree_Reset(TREE_DATA *);
int     Cet_Print_Vol_File(CET_OUT_DATA *,MKTVOL_DATA *,int,int);

/* cet_new.c */
int     Cet_Main_New(int,T_CURVE *,MKTVOL_DATA *,TREE_DATA *);

/* cet_mr.c */
int     Cet_MR_Main(int,T_CURVE *,MKTVOL_DATA *,TREE_DATA *);


/* chooser.c */
int     Chooser_t (double **,double *,long,int,int,int,int,DEV_DATA *,TREE_DATA *);
int     Survivor_t (double **,double *,double *,long,int,int,double,char,int,int,int,DEV_DATA *,TREE_DATA *);

/* claimbank.c */
void    CbkInit(CLAIM_BANK *);
int     CbkAlloc(CLAIM_BANK *,int,TREE_DATA *);
int     CbkFree(CLAIM_BANK *,TREE_DATA *);
int     CbkSizeFromDL(long *,int,long *,int);
int     CbkCalcSize(int,CRIT_DATE *,int,int);
int     CbkDev(CLAIM_BANK *,long,int,int,int,DEV_DATA *,TREE_DATA *);
double *CbkPopSlice(CLAIM_BANK *);
int     CbkPushSlice(CLAIM_BANK *,double *,long,long);
int     CbkProcessDL(int *,long **,int *,long **);
double *CbkReadSlice(CLAIM_BANK *,long);
int     CbkGetOffset(CLAIM_BANK *,long,int);

/* date.c */
void    Dsplit (long,long *,long *,long *);
int     Isleap (long);
int     Isimm(long);
long    ThirdWed(long, long);
int     Holiday (long);
int     Dateok (long);
long    Y2toy4 (long);
long    Datepack (long,long,long);
long    Y4toy2 (long);
void    Y2date_str (long,char *);
long    eval_date (char *);
long    eval_date2 (char *);
long    Dayofwk (long);
long    Days360 (long,long);
long    Daysact (long,long);
long    Months360 (long,long);
double 	YearFraction (long,long);
long    Nxtday (long,long);
long    Nxtmth (long,long,long);
long    Nxtwkday (long,long);
long    Nxtbusday (long,long);
long    Nxtimm(long, long);
int     AddDateToList(int *,long **,long);
int     SortDateList(int,long *,long *);
int     GetDLOffset(int,long *,long,int);
int         DrDayCountFraction (long,long,char,double *);
int         DrDateFwdAny (long,int,char,char,long *);
int         DateListFromFreq (long,long,char,char,int *,long **);  
EVENT_LIST  *DrNewEventListFromFreq (int,long *,char,char,char,double *,double *,double *,double *,double *);
int         DrDatesInSchedule (int,long *,long,long,char,char);
int         DrSameDateSchedules (int *,long *,long,long,long,char,char,char);
int         DrSameDateSets (int *,long *,int,long *,long);
int         DrSameDateSetsFlows (int *,long *,int,long *,long);
int         DrDatesIn2FreqSchedule (int,long *,long,long,long,char,char,char,int *,long *);
int         DrDatesInSet(int,long *,int,long *);
void        DrFreeEventList (EVENT_LIST  *);
int         hasStub (long,long,char);
int         Date_CheckAndReport(long);
EVENT_LIST * DrNewEventListFromFreqWithInterpType( 
	ESL_INTERP_TYPE,  /**< (I) Interpolation type to use*/
        int,  /**< (I) Nb of dates input directly   */
        long *,    /**< (I) Dates given directly by user */
        char,        /**< (I) Frequency of event           */
        char,        /**< (I) Stub location Front or Back  */
        char,     /**< (I) Y=input dates must be in list*/
        double *,      /**< (I) Set of values for event      */
        double *,      /**< (I) Set of values for event      */
        double *,      /**< (I) Set of values for event      */
        double *,      /**< (I) Set of values for event      */
        double *);     /**< (I) Set of values for event      */
unsigned 
ExpandDateSchedule(unsigned sSize,          /**< (I) size of the input array */
                   long const* sSched,  /**< (I) input array             */
                   unsigned tSize,          /**< (I) size of the output array*/
                   long* tSched,        /**< (I) output array            */
                   char frequency);         /**< (I) frequency-I,D,M,W,Q,S,A */


/* dev.c */
int     Dev (double *,int,int,int,DEV_DATA *,TREE_DATA *);
int     Ev  (double *,int,int,DEV_DATA *,TREE_DATA *);

/* drift.c */
int     Drift_1D (MKTVOL_DATA *,TREE_DATA *);
int     Drift_2D (MKTVOL_DATA *,TREE_DATA *);
int     Drift_3D (MKTVOL_DATA *,TREE_DATA *);

/* error.c */
void	DR_Error (char *format,  ...);

/* warning.c */
void    DR_Warning (int, char *, ...);

/* flexswap.c */
int     Flexswap_t (double **,double **,long *,char *, long,double,int,double *,double,double,double,double,long,double *,int,int,int,DEV_DATA *,TREE_DATA *);

/* floater.c */
int     Floater_t (double *,double *,double *,long,double,double,char,char,double,double,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     FwdFloater_t (double *,double *,double *,double *,long,double,double,long,char,double,char,long,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     CmpFloater_t (double *,double *,double *,long,long,double,double,char,char,double,long,int,int,int,DEV_DATA *,TREE_DATA *);

/* idxam.c */
int     IdxamSwap_t (double **,long,long,long,double *,double *,double,double,double,double,double,double,double *,double *,int,char,char,char,int,double,double,double **,int,int,int,DEV_DATA *,TREE_DATA *);


/* irrloan.c */
int     IRRLoanSwap_t (double **,long,long,double *,double *,double *,double,double,double,double,double,double,double,char,char,int,int,double,double,double **,int,TREE_DATA *);


/* iou.c */
int     Ioucap_t (double **,double *,double *,long,int,double,double,char,char,int,int,int,DEV_DATA *,TREE_DATA *);
int     Iouswap_t(double **,double *,double *,long,int,double,double,double,double **,char,char,char,int,int,int,DEV_DATA *,TREE_DATA *);


/* koopt.c */
int     KoOption_t (double *,double *,long,double,double,double,char,char,int,TREE_DATA *);
int     KiOption_t (double *,double *,double *,long,double,double,char,char,int,int,int,DEV_DATA *,TREE_DATA *);
int     Trigger_t (double *,double *,double *,double,double,long,int,double,double,char,char,int,int,int,DEV_DATA *, TREE_DATA *);
int     Top_t(double *,double *,double *,double,double,long,int,double,double,char,char,int,int,int,DEV_DATA *,TREE_DATA *);

/* kostream.c */
int     KoCap_t (double *,double *,double *,double *,double *,long,double,double,double,char,char,char,long,long,double,long,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     KoSwap_t (double *,double *,double *,double *,double *,double *,long,double,double,double,char,char,char,long,long,double,long,char,long,double,double,long,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     KoOptionVarRebate_t (double *,double *,double *,long,double,double,char,char,int,int,TREE_DATA *);
int     KoFwdSwap_t(double *,double *,double *,double *,double *,double *,double *,long,double,double,double,char,char,char,char,long,long,double,double,long,char,long,double,double,long,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);

/* ladder.c */
int     LadderSwap_t (double **,long,long,double *, double *, double *,double,double,double,double,char,char,int,double *,double *,int,TREE_DATA *);
int     LadderStep_t (double *,double *,double,double,double,double,double,char,int,TREE_DATA *);

/* lattice.c */
int     Lattice (DEV_DATA *,int,int,MKTVOL_DATA *,TREE_DATA *);

/* multiserial.cpp */
int     Multiserial_t(double **, int, double *, double, double, char, int, double **, double, double, int, TREE_DATA *);

/* opt.c */
int     Option_t (double *,double *,double,double,long,int,int,int,int,DEV_DATA *,TREE_DATA *);
int     OptionPlus_t (double *,double *,double *,double *,double *,double,double,long,int,int,double *,int,int,int,DEV_DATA *,TREE_DATA *);

/* paryield.c */
int     Discount_Bond(int, double *, long *, long, long, double *);
int     Par_Yield_Plus(double *, double *, int, double *, long *, long, long, char, int, char, char);
int     Par_Yield (double *,double *,int,double *,long *,long,long,int,char,char);
int     Par_Yield_From_Dates(double  *,double  *,long,long,char,char,char,int,double *,long *,long);
int     Par_Yield_t (double *,int,double **,long *,long,long,long,int,char,char,double,int,int,int,DEV_DATA *,TREE_DATA *);
int     ParYieldRatio(double *,long,long,double,int,long *,double *,long,int,char,char);

/* slice.c */
int     Set_Slice (double *,double,int,TREE_DATA *);
int     Copy_Slice (double *,double *,int,TREE_DATA *);
int     LCombTwoSlices(double *,double *,double,double *,double,int,TREE_DATA *);
int     AddTwoSlices(double *,double *,double *,int,TREE_DATA *);
int     MultiplyTwoSlices(double *,double *,double *,int,TREE_DATA *);
int     MultTwoSlicesAddAll(double *,double *,double *,int,TREE_DATA *);
int     AddScalar(double *,double,int,TREE_DATA *);
int     MultiplyScalar(double *,double,int,TREE_DATA *);
int     MaxMinOnSlice(double *,double,double,int,TREE_DATA *);
int     MaxOnSlice(double *,double,int,TREE_DATA *);
int     MinOnSlice(double *,double,int,TREE_DATA *);
int     Node_Offset (int,int,int,int, const TREE_DATA *);
int     Init_Slice (double *,double,TREE_DATA *);
int     MaxOfTwoSlices (double *,double *,double *,int,TREE_DATA *);
int     MinOfTwoSlices (double *,double *,double *,int,TREE_DATA *);
int     AccrueSlice(double *, double, int, int, TREE_DATA *);
int     DivideTwoSlices (double *,double *,double *,int,TREE_DATA *);
int     SmoothStepUp(
                  double*        Slice,           /**< (O) Slice to be modified */
                  double const*  Index,           /**< (I) Reference slice      */
                  double         Up,              /**< (I) Maximium level       */
                  double         Barrier,         /**< (I) Barrier level        */
		          int            SmoothingOn,     /**< (I) TRUE = smoothing on  */
                  int            t,               /**< (I) Current time period  */
                  TREE_DATA const*tree_data); /**< (I) Tree data            */
/* slice_oper.cpp */
int    KOIndicator(double *indicator, double const* index, double low, double high, char io, int t, TREE_DATA const* tree);

/* statevar.c */
double  StateVar_Interp(int, int, double *, double **, double[MAXNBSTATES][3], double);
int     StateVar_Generate(double **, int, double, double, char);

/* stdinput.c */
int     FindAndSkipComLine (FILE *,char const*,char const*,char const*);
int     Term_Input_W (T_CURVE *,char const*);
int     Term_Check_W (T_CURVE *);
int     BaseVol_Input_W (MKTVOL_DATA *,char const*);
int     SwapVol_Input_W (int *,int *,long **,long **,double ***,char const*);
int     MktVol_Input_W (MKTVOL_DATA *mktvol_data,char *,char *,T_CURVE *,char const*,char const*);
int     MktVol_Check_W (MKTVOL_DATA *);
int     Param_Input (MKTVOL_DATA *,TREE_DATA *,int,char[6][MAXBUFF],char const*);
int     Param_Check (int,MKTVOL_DATA *,TREE_DATA *);
int     MR_Input_W (MKTVOL_DATA *, char *);
int     MR_Input_W1 (MKTVOL_DATA *, char *);
int     MR_Input_W2 (MKTVOL_DATA *, char *);
int     MR_Input_W3 (MKTVOL_DATA *, FILE *, char *);
int     MrExpToDates(int, MKTVOL_DATA *);
int     MRTS_Input_W (MKTVOL_DATA *, FILE *, char *);

/* sticky.c */
int     Sticky_t(double**,long,double*,double*,double,double,double,double,double,double,double,double,double,double,double,double,int,int,int,int,int,double,double,double**,int,int,int,DEV_DATA*,TREE_DATA*); 
int     StickySwap_t(double **,long,long,double *,double *,double *,double *,double,double,double,double,double,double,double,double,double,double,double,double,double,int,int,int,int,int,double,double,double **,int,TREE_DATA *);
int     VolBond_t(double **,long,double *,double *,double,double,double,double,double,double,int,int,double,double,double **,int,TREE_DATA *);
int     VolSwap_t(double **,long,long,double *,double *,double *,double *,double,double,double,double,double,double,double,double,double,double,double,int,int,double,double,double **,int,TREE_DATA *);


/* swap.c */
int     ParSwap_t (double *,double *,double *,double,long,double,char,long,double,long,double,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     AmortSwap_t (double *,double *,double *,double *,double *,double *,long,double,double,double,long,char,long,double,double,long,char,double,char,long,double,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     Swaplet_t (double *,double *,double *,long,double,double,char,double,char,double,int,int,int,DEV_DATA *,TREE_DATA*);
int     FwdSwap_t (double *,double *,double *,double *,double *,double *,double *,long,double,double,double,long,char,long,double,double,long,char,char,long,double,char,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     FwdSwapFlows_t(double *,CLAIM_BANK *,CLAIM_BANK *,double *,double,double,double,long,long,long,char,char,double *,double,double *,double,double *,double,long,long,long,long,char,char,char,double,long,int,TREE_DATA *);

/* swapvol.c */
double	ExpDecay (double,double);
int     BFactor (double *,long,long,char,char,int,double *,double,double,double,int,double *,long *,long);
int     SpotVol (double[6][MAXNBDATE],double[3][MAXNBDATE],long *,int *,long,int,long *,double *,int *,char,char,long *,long *,double,double,double,double,double,double,int,double *,int,long *,double[3][MAXNBDATE],double *,int,int,int,double *,long *,long);
int     Interp_SpotVol (double **,double **,int *,long *,double **,long,int,long *,double[6][MAXNBDATE],int,int,long *,int,long *,double[3][MAXNBDATE],double,double,double,int,double *,double *,int, long);
int     IndexVol (double *,long,long,char,char,int,int,long,long *,double[6][MAXNBDATE],double,double,double,double,double,double,int,int,long *,double[3][MAXNBDATE],double[3][MAXNBDATE],int,double *,long *,long);
int     Filtered_SpotVol (double[6][MAXNBDATE],double[3][MAXNBDATE],long *,int *,long,int,long *,double *,int *,char,char,long *,long *,double,double,double,double,double,double,int,double *,int,long *,double[3][MAXNBDATE],double *,int,int,double *,long *,long);
int     Interp_SpotVol_New (double **,double **,int *,long *,double **,long,int,long *,double[6][MAXNBDATE],int,int,long *,int,long *,double[3][MAXNBDATE],double,double,double,int,double *,double *,int, long);
int     BFactor_New (double *,long,long,char,char,int,int,long *,double[3][MAXNBDATE],double,double,double,int,double *,long *,long);

/* statvar.c */
int     DoubleArrayFloorIdx(double *, int, double);
int     DoubleVectSort(double *, int);
int     DoubleQuadraticInterp(double *, double *, int, double, double *);
int     Initialize_IR_SIM (IR_SIM *, long, int);
int     SwapYield_MC (double[MAXNBSTATES], double *, double *,int,long,long,char,char,IR_SIM *,int,int,long,long *,double[6][MAXNBDATE],double,double,double,double,double,double,int,double *,int,double *,long *,long);
int     IndexCovar (double *, double *, double *, long,long,char,char,long,long, char,char,int,int,long,long *,double[6][MAXNBDATE],double,double,double,int,double *,int,double *,long *,long);
int     LCoef (double[3][3],long,int,int,long ,long *, double[6][MAXNBDATE], int, double *);

/* time.c */
int     Add_To_DateList (int *,CRIT_DATE **,long,int,double,double,double,double,double,long,long,long);
int     Time_Line (long,int,CRIT_DATE *,char,TREE_DATA *);
int     TimeInterp (long,char *,char,int *,long *,double *,double *,double *,double *,double *);
int     Sort_CritDate (int,CRIT_DATE *);
int     DrExtendAmerSch(long,int,int *,long **,long **,char,double **,double **,double **,double **,double **);

/* tree.c */
int     Forward_Curve (double *,double *,double *,long,int,double *,long *,long *,int);
int     Tree_Limits (int *,int *,int *,int *,int **,int **,int ***,int ***,int,int,int,double **,double **,double *,double *, int, long*, double **);
int     Build_Tree (T_CURVE *,MKTVOL_DATA *,TREE_DATA *);
int     Set_Tree_MR(TREE_DATA *,MKTVOL_DATA *);

/* util.c */
double  GetIndexStep (double const*,int,int,int,int,int,TREE_DATA const*);
int     Smooth_Step (double *,double,double,double,double,double);
double  DrSmoothStep(double,double,double,double,double);
double  DrSmoothMax(double,double,double);
int     Conv_Freq (char);
int     Conv_Index (int *,char *,char *,char *,T_CURVE *);
void    pksort2 (int,long *,long *);
double  NR_Poly (double,double *,int);
int     Quadratic_Solve (double *, double *);
int     gaussj (double  M[3][3],int ,double *);
int     polint (double *,double *,int,double,double *,double *);
int     dlinterp (long,double *,long,long,double,double);
int     linterp(double,double *,double,double,double,double);
void    qinterp(double *,double *,double,double *,double);
void    d4interp(double,double *,double *,double *);
int     sqinterp(double,double,double,double,double,double,double,double,double,double,double*);
void    tableinterp (double,double *,double *,double *,int); 
double  AnnPV (double,double,int,int);
int     AnnIRR (double,double,int,int,double *);
int     Conv_Index_Flex(int *, char *, char *, char *, char *, T_CURVE *);

/* util_nr.c */
int     gaussjORI(double **,int,double **,int);
int     MatrixInverse(long,double **,double **);

/* zerobank.c */
int     ZbkUpdate(CLAIM_BANK *,int,long,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     ZbkDLFromIdx(int,long *,long *,int,char,int *,long **,int *,long **);
double *ZbkReadZero(CLAIM_BANK *,long,int,long,int,TREE_DATA *);
int     ZbkOptDates(int,long *,int,long *,int,CRIT_DATE *,long,int *,long **,int *,long **ZbkErs);
int     ZbkParYield_t(double *,double *,CLAIM_BANK *,long,long,int,char,char,double,int,TREE_DATA *);
int     ZbkAnnuity_t(double *,CLAIM_BANK *,long,long,int,char,char,int,TREE_DATA *);

/* zerobank_plus.cpp */
int     ZbkDLFromIdx_Plus(long, int, long *, long *, char, int, char, int *, long **, int *, long **);
	
/* zerobinary.c */
int     ZeroBinary_t(double **, long, double *, double, char, int, double **, double, double, int, TREE_DATA *);

/* zeros.c */
int     Zero_t (double *,long,int,int,int,DEV_DATA *,TREE_DATA *);
int     Zero_Bank (double **,long *,int *,int,long,long,int,int,int,DEV_DATA *,TREE_DATA *);
double  ZeroPrice(long,long,int,long *,double *);

/* zip.c */
int  Zip_t(double *,long,long,double,double *,double,double,double,double *,double *,int,int,int,int,DEV_DATA *,TREE_DATA *);

/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif /* _fix123head_h */
