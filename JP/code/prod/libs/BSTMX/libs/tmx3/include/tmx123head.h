/****************************************************************************/
/*      Function templates for library.                                     */
/****************************************************************************/
/*      TMX123HEAD.h                                                        */
/****************************************************************************/

/* Use safe muliple inclusion for Jerry Cohen's stream libraries */

#ifndef    _tmx123head_h
#define    _tmx123head_h

#include   "tmx123.h"

/* alloc.c */
void        Tree_Init (TREE_DATA *);
int         Tree_Alloc(TREE_DATA *);
int         Tree_Free (TREE_DATA *);
void        Dev_Init  (DEV_DATA    *dev_data);
int         Dev_Alloc (DEV_DATA *,TREE_DATA *);
int         Dev_Free (DEV_DATA *,TREE_DATA *);
void        MktVol_Init(MKTVOL_DATA *);
void        Opt_Out_Data_Init(OPT_OUT_DATA  *ood);
double     *Alloc_Slice (TREE_DATA *);
int         Free_Slice (double *,TREE_DATA *);
void       *DR_Array (int,int,int);
void       *DR_Matrix (int,int,int,int,int);
void       *DR_Cube (int,int,int,int,int,int,int);
int         Free_DR_Array (void *,int,int,int);
int         Free_DR_Matrix (void *,int,int,int,int,int);
int         Free_DR_Cube (void *,int,int,int,int,int,int,int);
int        *Alloc_Slice_Int (TREE_DATA *);
int         Free_Slice_Int (int *,TREE_DATA *);

/* amortloan.c */
int         AnnAmortSwap_t (double **,double *,double *,double *,char,char,double,double,double,long,long,long,int,double,double,double **,int,int,int,DEV_DATA *,TREE_DATA *);

/* bas.c */
double      NDensity (double);
double      NormalH (double);
double      Normal_InvH (double);
double      Call_BSQ (double,double,double,double,double);
double      Put_BSQ (double,double,double,double,double);
double      Vega_BSQ (double,double,double,double,double);
double      Binary_BS (double,double,double,double,double,char);
double      Call_BS (double,double,double,double,double);
double      Put_BS (double,double,double,double,double);
double      D_Call_BS (double,double,double,double,double);
double      D_Put_BS (double,double,double,double,double);
double      Gamma_BS (double,double,double,double,double);
double      Vega_BS (double,double,double,double,double);
double      T_Call_BS (double,double,double,double,double);
double      T_Put_BS (double,double,double,double,double);
double      R_Call_BS (double,double,double,double,double);
double      R_Put_BS (double,double,double,double,double);
double      CCEquation_BS2Q (double,double,double,double,double);
int         Tmx_ConvexityC_BS2Q (double,double,double,double,double *);
double      Option_BS2Q (double,double,double,double,char,double,double,double);
double      ImpVol_BS2Q (double,double,double,double,char,double,double,double,double);
double      Option_BSQ (double,double,double,double,char,double);
double      ImpVol_BSQ (double,double,double,double,char,double,double);

/* bond.c */
int         Bond_t   (double *,double,long,double,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         FwdBond_t(double *,double *,double *,long,double,double,char,long,double,double,long,double,char,long,int,int,int,DEV_DATA *,TREE_DATA *);

/* callzero.c */
int         Callzero_t (double *,double *,double *,double *,double *,long,long,long,double,double,double,char,char,char,int,int,int,DEV_DATA *,TREE_DATA *);

/* cap.c */
int         Cap_t (double *,double *,double *,int,long,double,double,char,double,char,double,int,int,int,DEV_DATA *,TREE_DATA *);
int         Caplet_t (double *,double *,double *,int,long,double,double,char,double,char,double,int,int,int,DEV_DATA *,TREE_DATA *);
int         Collaret_t (double *,double *,double *,long,double,double,double,char,double,char,double,int,int,int,DEV_DATA *,TREE_DATA *);

/* cet.c */
int         Cet_Main (int,T_CURVE *,MKTVOL_DATA *,TREE_DATA *);
int         Cet_Calc (T_CURVE *,MKTVOL_DATA *,TREE_DATA *);
int         Cet_Manager (TREE_DATA *,TREE_DATA *); 
int         Cet_NmrVanl (long,T_CURVE *,MKTVOL_DATA *,int);
int         Cet_Print (int,double *,MKTVOL_DATA *);
int         Cet_Schedule (long,MKTVOL_DATA *,TREE_DATA *);
int         Cet_SwapVols (int,int , double *, double *, double *, double *, MKTVOL_DATA *, TREE_DATA *, DEV_DATA *, int); 
void        Cet_Tree_Reset(TREE_DATA *);
int         Cet_UpdSwapSml(long,MKTVOL_DATA *);
             

/* chooser.c */
int         Chooser_t (double **,double *,long,int,int,int,int,DEV_DATA *,TREE_DATA *);
int         Survivor_t (double **,double *,double *,long,int,int,double,char,int,int,int,DEV_DATA *,TREE_DATA *);

/* claimbank.c */
void        CbkInit(CLAIM_BANK *);
int         CbkAlloc(CLAIM_BANK *,int,TREE_DATA *);
int         CbkFree(CLAIM_BANK *,TREE_DATA *);
int         CbkSizeFromDL(long *,int,long *,int);
int         CbkCalcSize(int,CRIT_DATE *,int,int);
int         CbkDev(CLAIM_BANK *,long,int,int,int,DEV_DATA *,TREE_DATA *);
double     *CbkPopSlice(CLAIM_BANK *);
int         CbkPushSlice(CLAIM_BANK *,double *,long,long);
int         CbkProcessDL(int *,long **,int *,long **);
double     *CbkReadSlice(CLAIM_BANK *,long);
int         CbkGetOffset(CLAIM_BANK *,long,int);

/* date.c */
int         Isimm(long);
long        ThirdWed(long, long);
int         Holiday (long);
long        Y2toy4 (long);
long        Datepack (long,long,long);
long        Y4toy2 (long);
void        Y2date_str (long,char *);
long        eval_date (char *);
long        eval_date2 (char *);
long        Dayofwk (long);
long        Days360 (long,long);
long        Months360 (long,long);
double      YearFraction (long,long);
long        Nxtday (long,long);
long        Nxtmth (long,long,long);
long        Nxtwkday (long,long);
long        Nxtbusday (long,long);
long        Nxtimm(long, long);
int         AddDateToList(int *,long **,long);
int         SortDateList(int,long *,long *);
int         MergeDateLists(int,long *,int *,long **);
int         DrDayCountFraction (long,long,char,double *);
int         DrDateFwdAny (long,int,char,char,long *);
int         DateListFromFreq (long,long,char,char,int *,long **);  
EVENT_LIST *DrNewEventListFromFreqWithInterpType(INTERP_TYPE,int,long *,char,char,char,double *,
                                                 double *,double *,double *,double *);
EVENT_LIST *DrNewEventListFromFreq (int,long *,char,char,char,double *,double *,double *,double *,double *);
int         DrDatesInSchedule (int,long *,long,long,char,char);
int         DrSameDateSchedules (int *,long *,long,long,long,char,char,char);
int         DrSameDateSets (int *,long *,int,long *,long);
int         DrSameDateSetsFlows (int *,long *,int,long *,long);
int         DrDatesIn2FreqSchedule (int,long *,long,long,long,char,char,char,int *,long *);
int         DrDatesInSet(int,long *,int,long *);
void        DrFreeEventList (EVENT_LIST  *);
int         hasStub (long,long,char);
int         Date_CheckAndReport(long);
int         NearestDate(int,long *,long);

/* dev.c */
int         Dev (double *,int,int,int,DEV_DATA *,TREE_DATA *);
int         Ev  (double *,int,int,    DEV_DATA *,TREE_DATA *);

/* drift.c */
int         Drift_1D (MKTVOL_DATA *, TREE_DATA *);


/* flexswap.c */
int         Flexswap_t (double **,double **,long *,char *, long,double,int,double *,double,double,double,double,long,double *,int,int,int,DEV_DATA *,TREE_DATA *);

/* floater.c */
int         Floater_t (double *,double *,double *,long,double,double,char,char,double,double,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         FwdFloater_t (double *,double *,double *,double *,long,double,double,long,char,double,char,long,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         CmpFloater_t (double *,double *,double *,long,long,double,double,char,char,double,long,int,int,int,DEV_DATA *,TREE_DATA *);

/* idxam.c */
int         IdxamSwap_t (double **,long,long,long,double *,double *,double,double,double,double,double,double,double *,double *,int,char,char,char,int,double,double,double **,int,int,int,DEV_DATA *,TREE_DATA *);

/* irrloan.c */
int         IRRLoanSwap_t (double **,long,long,double *,double *,double *,double,double,double,double,double,double,double,char,char,int,int,double,double,double **,int,TREE_DATA *);

/* iou.c */
int         Ioucap_t (double **,double *,double *,long,int,double,double,char,char,int,int,int,DEV_DATA *,TREE_DATA *);
int         Iouswap_t(double **,double *,double *,long,int,double,double,double,double **,char,char,char,int,int,int,DEV_DATA *,TREE_DATA *);

/* koopt.c */
int         KoOption_t (double *,double *,long,double,double,double,char,char,int,TREE_DATA *);
int         KiOption_t (double *,double *,double *,long,double,double,char,char,int,int,int,DEV_DATA *,TREE_DATA *);
int         Trigger_t (double *,double *,double *,double,double,long,int,double,double,char,char,int,int,int,DEV_DATA *, TREE_DATA *);
int         Top_t(double *,double *,double *,double,double,long,int,double,double,char,char,int,int,int,DEV_DATA *,TREE_DATA *);

/* kostream.c */
int         KoCap_t (double *,double *,double *,double *,double *,long,double,double,double,char,char,char,long,long,double,long,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         KoSwap_t (double *,double *,double *,double *,double *,double *,long,double,double,double,char,char,char,long,long,double,long,char,long,double,double,long,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         KoOptionVarRebate_t (double *,double *,double *,long,double,double,char,char,int,int,TREE_DATA *);
int         KoFwdSwap_t(double *,double *,double *,double *,double *,double *,double *,long,double,double,double,char,char,char,char,long,long,double,double,long,char,long,double,double,long,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);

/* ladder.c */
int         LadderSwap_x (double **,long,long,double *, double *, double *,double,double,double,double,char,char,int,double *,double *,int,TREE_DATA *);
int         LadderStep_x (double *,double *,double,double,double,double,double,char,int,TREE_DATA *);

/* lattice.c */
int         Lattice (DEV_DATA *,int,int,MKTVOL_DATA *,TREE_DATA *, T_CURVE *);

/* numer.c */
int         Nmr_Alloc(MKTVOL_DATA *,TREE_DATA *);
int         Nmr_Schedule(long,long,MKTVOL_DATA *);
int         Nmr_Calc(T_CURVE *,MKTVOL_DATA *,TREE_DATA *);
int         Nmr_Interp(double *,int,int,MKTVOL_DATA *,TREE_DATA *, T_CURVE *);
int         TMX_StateVar(double *, double *, long, double, long, long, char, char, MKTVOL_DATA *, T_CURVE *);

/* opt.c */
int         Option_t (double *,double *,double,double,long,int,int,int,int,DEV_DATA *,TREE_DATA *);
int         OptionPlus_t (double *,double *,double *,double *,double *,double,double,long,int,int,double *,int,int,int,DEV_DATA *,TREE_DATA *);

/* paryield.c */
int         Par_Yield_Plus(double *, double *, int, double *, long *, long, long, char, int, char, char);
int         Par_Yield (double *,double *,int,double *,long *,long,long,int,char,char);
int         Par_Yield_From_Dates(double  *,double  *,long,long,char,char,char,int,double *,long *,long);
int         Par_Yield_t (double *,int,double **,long *,long,long,long,int,char,char,double,int,int,int,DEV_DATA *,TREE_DATA *);
int         ParYieldRatio(double *,long,long,double,int,long *,double *,long,int,char,char);

/* slice.c */
int         Set_Slice (double *,double,int,TREE_DATA *);
int         Copy_Slice (double *,double *,int,TREE_DATA *);
int         LCombTwoSlices(double *,double *,double,double *,double,int,TREE_DATA *);
int         AddTwoSlices(double *,double *,double *,int,TREE_DATA *);
int         MultiplyTwoSlices(double *,double *,double *,int,TREE_DATA *);
int         DivideTwoSlices(double *,double *,double *,int,TREE_DATA *);
int         MultTwoSlicesAddAll(double *,double *,double *,int,TREE_DATA *);
int         AddScalar(double *,double,int,TREE_DATA *);
int         MultiplyScalar(double *,double,int,TREE_DATA *);
int         MaxMinOnSlice(double *,double,double,int,TREE_DATA *);
int         MaxOnSlice(double *,double,int,TREE_DATA *);
int         MinOnSlice(double *,double,int,TREE_DATA *);
int         Node_Offset (int,int,int,int,TREE_DATA *);
int         Init_Slice (double *,double,TREE_DATA *);
int         MaxOfTwoSlices (double *,double *,double *,int,TREE_DATA *);
int         MinOfTwoSlices (double *,double *,double *,int,TREE_DATA *);
int         SlicePlusScalar (double *,double *,double,int,TREE_DATA *);
int         SliceTimesScalar (double *,double *,double,int,TREE_DATA *);
int         AccrueSlice(double *, double, int, int, TREE_DATA *);
int         Print_Slice1D (double *,char *,int, TREE_DATA *);

/* stdinput.c */
int         BaseSmile_Input_W (MKTSMILE_DATA *,int,int,long,char *);
int         SwapSmile_Input_W (MKTSMILE_DATA *,int,int,long,int,char *);
int         MktVol_Input_W (MKTVOL_DATA *,char *,T_CURVE *,char[MAXBUFF],char *,char *,char *,char *,char *);
int         MktVol_Check_W (MKTVOL_DATA *);
int         MktSmile_Input_W (MKTSMILE_DATA *,char[MAXBUFF],long,char,int,char *);
int         smileinterp (MKTVOL_DATA *, MKTSMILE_DATA *); 
int         Convert_VoV (MKTVOL_DATA *, T_CURVE *);

/* sticky.c */
int         Sticky_t(double**,long,double*,double*,double,double,double,double,double,double,double,double,double,double,double,double,int,int,int,int,int,double,double,double**,int,int,int,DEV_DATA*,TREE_DATA*); 
int         StickySwap_t(double **,long,long,double *,double *,double *,double *,double,double,double,double,double,double,double,double,double,double,double,double,double,int,int,int,int,int,double,double,double **,int,TREE_DATA *);
int         VolBond_t(double **,long,double *,double *,double,double,double,double,double,double,int,int,double,double,double **,int,TREE_DATA *);
int         VolSwap_t(double **,long,long,double *,double *,double *,double *,double,double,double,double,double,double,double,double,double,double,double,int,int,double,double,double **,int,TREE_DATA *);

/* swap.c */
int         ParSwap_t (double *,double *,double *,double,long,double,char,long,double,long,double,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         AmortSwap_t (double *,double *,double *,double *,double *,double *,long,double,double,double,long,char,long,double,double,long,char,double,char,long,double,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         Swaplet_t (double *,double *,double *,long,double,double,char,double,char,double,int,int,int,DEV_DATA *,TREE_DATA*);
int         FwdSwap_t (double *,double *,double *,double *,double *,double *,double *,long,double,double,double,long,char,long,double,double,long,char,char,long,double,char,char,char,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         FwdSwapFlows_t(double *,CLAIM_BANK *,CLAIM_BANK *,double *,double,double,double,long,long,long,char,char,double *,double,double *,double,double *,double,long,long,long,long,char,char,char,double,long,int,TREE_DATA *);

/* swapvol.c */
double      TMX_ExpDecay (double,double);
int         Interp_SpotVol(double  **,long , int , long *, double [6][MAXNBDATE], long, double  *, int );
int         BFactor (double *,long,long,char,char,double,double,double*, long, int,double *,long *,long);
int         SpotVol (double [6][MAXNBDATE], double [NBVOLPARS][MAXNBDATE], int, long *, int *, long , int , long *, double  [NBVOLPARS][MAXNBDATE], int *, char , char , long *, long *, int , double *, double  *, double *, int , int , double *, long *, long );          
int         IndexVol (double *,long,long,char,char,int,long,long *,double [6][MAXNBDATE],int,double *,int,double *,
                      long *,long);
int         StateVar(double *,double *,double,long,long,char,char,int,MKTVOL_DATA *,T_CURVE *);
 
/* time.c */
int         Add_To_DateList (int *,CRIT_DATE **,long,int,double,double,double,double,double,long,long,long);
int         Time_Line (long,int,CRIT_DATE *,char,TREE_DATA *);
int         TimeInterp (long,char *,char,int *,long *,double *,double *,double *,double *,double *);
int         Sort_CritDate (int,CRIT_DATE *);
int         DrExtendAmerSch(long,int,int *,long **,long **,char,double **,double **,double **,double **,double **);

/* tree.c */
int         Forward_Curve (double *,double *,double *,double *,long,int,double *,long *,long *,int);
int         Tree_Limits (int *,int *,int *,int *,int **,int **,int ***,int ***,int,int,int,double *,double **,double *,double *);
int         Build_Tree (T_CURVE *,MKTVOL_DATA *,TREE_DATA *);

/* util.c */
double      GetIndexStep (double *,int,int,int,int,int,TREE_DATA *);
int         Smooth_Step (double *,double,double,double,double,double);
double      DrSmoothStep(double,double,double,double,double);
double      DrSmoothMax(double,double,double);
int         Conv_Freq (char);
int         Conv_Index (int *,char *,char *,char *,T_CURVE *);
void        pksort2 (int,long *,long *);
double      NR_Poly (double,double *,int);
int         Quadratic_Solve (double *, double *);
int         gaussj (double  M[3][3],int ,double *);
int         polint (double *,double *,int,double,double *,double *);
int         linterp(double,double *,double,double,double,double);
int         linterp2d (double,double,double *,double,double,double,double,double,double,double,double);
void        qinterp(double *,double *,double,double *,double);
void        d4interp(double,double *,double *,double *);
int         sqinterp(double,double,double,double,double,double,double,double,double,double,double*);
void        tableinterp (double,double *,double *,double *,int); 
int         DoubleArrayFloorIdx(double *, int, double);
int         DoubleVectSort(double *, int);
int         DoubleQuadraticInterp(double *, double *, int, double, double *);
void        matrixinterp (double,double,double *,double *,double *,double **,int,int);
double      AnnPV (double,double,int,int);
int         AnnIRR (double,double,int,int,double *);
int         Conv_Index_Flex(int *, char *, char *, char *, char *, T_CURVE *);
int         DrlSplineInterp1d(double *, double *, int , double , double , double ,  double *);      
int         DrlSplineInterp1dInit(double *, double *, int , double , double , double *);        
int         DrlSplineInterp1dInterp(double *, double *, double *, int , double ,double *);  

/* warning.c */
void        DR_Warning (int, char *, ...);

/* zerobank.c */
int         ZbkUpdate(CLAIM_BANK *,int,long,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         ZbkDLFromIdx(int,long *,long *,int,char,int *,long **,int *,long **);
double     *ZbkReadZero(CLAIM_BANK *,long,int,long,int,TREE_DATA *);
int         ZbkOptDates(int,long *,int,long *,int,CRIT_DATE *,long,int *,long **,int *,long **ZbkErs);
int         ZbkAnnuity_t(double *,CLAIM_BANK *,long,long,int,char,char,int,TREE_DATA *);
int         ZbkParYield_t(double *, double *, CLAIM_BANK *, long, long, int, char, char, double, int, TREE_DATA *);

/* zeros.c */
int         Zero_t (double *,long,int,int,int,DEV_DATA *,TREE_DATA *);
int         Zero_Bank (double **,long *,int *,int,long,long,int,int,int,DEV_DATA *,TREE_DATA *);
double      ZeroPrice(long,long,int,long *,double *);

/* zip.c */
int         Zip_t(double *,long,long,double,double *,double,double,double,double *,double *,int,int,int,int,DEV_DATA *,TREE_DATA *);


#if defined(CRXFLOW)
extern void        Dsplit (long,long *,long *,long *);
extern int         Isleap (long);
extern int         Dateok (long);
extern long        Daysact (long,long);
extern int         GetDLOffset(int,long *,long,int);
extern int         Term_Input_W (T_CURVE *,char *);
extern int         Term_Check_W (T_CURVE *);
extern void        DR_Error (char *format, ...);
double      TMX_ExpDecay (double,double);
int         TMX_dlinterp (long,double *,long,long,double,double);
int         TMX_FindAndSkipComLine (FILE *,char *,char *,char *);
int         TMX_Param_Input (MKTVOL_DATA *,TREE_DATA *,int,char[6][MAXBUFF],char *, char *);
int         TMX_Param_Check (int,MKTVOL_DATA *,TREE_DATA *);
#else
void        Dsplit (long,long *,long *,long *);
int         Isleap (long);
int         Dateok (long);
long        Daysact (long,long);
int         GetDLOffset(int,long *,long,int);
int         Term_Input_W (T_CURVE *,char *);
int         Term_Check_W (T_CURVE *);
void        DR_Error (char *format, ...);
int         dlinterp (long,double *,long,long,double,double);
double      ExpDecay (double,double);
int         FindAndSkipComLine (FILE *,char *,char *,char *);
int         Param_Input (MKTVOL_DATA *,TREE_DATA *,int,char[6][MAXBUFF],char *, char *);
int         Param_Check (int,MKTVOL_DATA *,TREE_DATA *);
#endif


#endif /* _tmx123head_h */
