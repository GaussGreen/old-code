/****************************************************************************/
/*      TMX3: Numeraire calculation and utilities                           */
/****************************************************************************/
/*      numer.c                                                             */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bmx123head.h"
#include "q3.h"


static int Nmr_t(double **,double *,double *,int,int,MQDATA *,MKTVOL_DATA *,TREE_DATA *);                  
static int Nmr_Print0(MKTVOL_DATA *);
static int Nmr_Print1(MKTVOL_DATA *,MQDATA *);
static int Nmr_Print2(MKTVOL_DATA *, TREE_DATA *);
static int Nmr_Print3(double *,double *,int,int,double,int,int,MKTVOL_DATA *,TREE_DATA *);
static int Nmr_Stat(int,int,double *,int *,double *,MKTVOL_DATA *,TREE_DATA *);
static int Nmr_Vanl_Init(long,T_CURVE *,MQDATA *,MKTVOL_DATA *,TREE_DATA *);
static int Nmr_Init(MKTVOL_DATA *);

/* New functions added by BM */
static int SortTree_Slices(double *, double *, double *, double *, double *, long *, int , double *, double *, double *, double *, double [2], long , long , TREE_DATA *);
static int SortList(long, double *, long *, long *, long *);   
static int InvMap(double *, double *, long , double[2], long , long, TREE_DATA *);         
static int Yield_t(double **, double [2], double *, double, double, double, double, int , MQDATA *, MKTVOL_DATA *, TREE_DATA *);
static int Yield2_Vanl(double *,double *, double *, double *, double *, long, long, long, long, char , char , T_CURVE *, MQDATA *, MKTVOL_DATA *, long, TREE_DATA *); 


/*****  Nmr_Schedule *****************************************************/
/*
* TMX3: Generates numeraire dates and inserts them into MKTVOL_DATA
* structure. The numeraire dates are regularly spaced, with frequency 
* determined by the index rate frequency. At nmr date compute swap maturity 
 * from adjacent benchmarks.
*/
int Nmr_Schedule (long        ValueDate,       /* (I) Value date          */
                  long        LastProdDate,    /* (I) Last product date   */
                  MKTVOL_DATA *mvd)            /* (I/O) Mkt vol data      */
{
    int     status = FAILURE; /* Error status = FAILURE initially */  

    double  *SwapStT=NULL,*TenorT=NULL;    
    double  Expiry;
    long    TermDate;
    int     intvl, row, NbExpiry, NbNmr, Mat, MatL, MatR;
    int     i, j, s;

    
    /* Choose  benchmarks: last one must expire AFTER LastProdDate */
    TermDate = LastProdDate;
    for (i = 0; i < mvd->NbVol; i++)
    {
        TermDate = MAX(TermDate, mvd->SwapMat[i]);
        if (mvd->SwapSt[i] >= LastProdDate) break;
    }

    /* Number of benchmarks actually used */
    NbExpiry = MIN(i+1, mvd->NbVol);

    /* Compute number of numeraire dates and shift terminal date to multiple 
     * of swap frequency. There must be at least 2 numeraire dates on or after
     * TermDate */
    intvl = 12 / Conv_Freq(mvd->Freq);
    NbNmr = (int) ceil(Months360(ValueDate, TermDate) / (double)intvl) + 2;

    /* Update terminal date and store in MKTVOL_DATA */
    TermDate = Nxtmth(ValueDate, NbNmr * intvl,1L);
    mvd->TermDate = TermDate;

    /* Count value date as well and copy to MKTVOL_DATA */
    NbNmr++;
    mvd->NbNmr = NbNmr;

    if (NbNmr > MAXNBDATE)
    {
        DR_Error ("Nb of numeraire dates %d exceeds maximum of %d. "
                  "(Nmr_Dates)\n", NbNmr, MAXNBDATE);
        goto RETURN;
    }  

    /* Interpolation variables */
    SwapStT  = (double *) DR_Array(DOUBLE,0,NbExpiry-1);
    TenorT   = (double *) DR_Array(DOUBLE,0,NbExpiry-1);
    for (i=0; i<NbExpiry; i++)
    {
        SwapStT[i]  = Daysact(ValueDate, mvd->SwapSt[i])/365.;
        TenorT[i]   = (double) Months360(mvd->SwapSt[i], mvd->SwapMat[i]);
    }

    /* Create MKTVOL_DATA numeraire dates. Nothing to do on TermDate */
    i = NbNmr-1;
    mvd->NmrDate[i]    = TermDate;
    mvd->NmrSwapSt[i]  = 0L; /* TermDate; */
    mvd->NmrSwapMat[i] = 0L; /* TermDate; */

    for (i = 0; i < NbNmr-1; i++) 
    {
        double  NmrDateT;
        double  NmrTenorT;

        /* Numeraire dates */
        mvd->NmrDate[i]   = Nxtmth(ValueDate, i * intvl, 1L);
        mvd->NmrSwapSt[i] = mvd->NmrDate[i];

        /* Interpolate swap info */
        NmrDateT = Daysact(ValueDate, mvd->NmrDate[i])/365.;
        tableinterp(NmrDateT, 
                    &NmrTenorT, 
                    SwapStT,
                    TenorT,
                    NbExpiry);
        NmrTenorT = intvl*floor(NmrTenorT/intvl+0.5);
        
        /* Calculate interpolated swap maturity */
        mvd->NmrSwapMat[i] = Nxtmth(mvd->NmrSwapSt[i], (int)NmrTenorT, 1);
        mvd->NmrSwapMat[i] = MIN(mvd->NmrSwapMat[i],TermDate);
    }

    /* Initialize MKTVOL_DATA */
    if (Nmr_Init(mvd) == FAILURE) goto RETURN;

    status = SUCCESS;

RETURN:

    Nmr_Print0(mvd);

    Free_DR_Array (SwapStT, DOUBLE, 0, NbExpiry-1);
    Free_DR_Array (TenorT,  DOUBLE, 0, NbExpiry-1);

    if (status == FAILURE)
    {
        DR_Error("Nmr_Schedule: Failed!");
    }

    return (status);

}  /* Nmr_Schedule */


/*****  Nmr_Init *****************************************************/
/*
*/
int Nmr_Init (MKTVOL_DATA *mvd)           /* (I/O) Mkt vol data      */
{
    int     n, s;


    for (n=0; n<mvd->NbNmr; n++)
    {
        for (s=0; s<NBSTRIKE; s++)
        {
            mvd->Strike[n][s] = 0;
            mvd->VanlPr[n][s] = 0;
            mvd->TreePr[n][s] = 0;
            mvd->Strike[n][s] = 0;
        }
        for (s=0; s<5; s++)
        {
            mvd->AASta[n][s]   = 0;
            mvd->OUSta[n][s]   = 0;
            mvd->TreeSta[n][s] = 0;
        }
    }

    return SUCCESS;
}


/*****  Nmr_Alloc ********************************************************/
/*
*/
int Nmr_Alloc (MKTVOL_DATA *mvd,             /* (I)    Mkt vol data      */
               TREE_DATA   *tree_data)       /* (I/O)  Tree data         */
{
    int   status =  FAILURE; /* Error status = FAILURE initially */
    int   NbNmr  = mvd->NbNmr;
    long  Area   = 0;
    int i;

    /* Slices Size */
    if (tree_data->NbFactor == 1)
    {
        Area = tree_data->Width[0];
    }
    if (tree_data->NbFactor == 2)
    {
        Area = tree_data->Width[0] * tree_data->Width[1];
    }
    if (tree_data->NbFactor == 3)
    {
        Area = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];
    }


    if (NbNmr < 0) goto RETURN;

    /* Allocate inverse numeraire and dates */
    tree_data->NmrInv  = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 
    tree_data->NmrDate = (long *)   DR_Array(LONG,      0,NbNmr-1);

    /* Allocate annuity and zero slices, and yield mapping arrays */
    tree_data->Ann1    = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 
    tree_data->Ann0    = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 
    tree_data->Zero1   = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 
    tree_data->Zero0   = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 

    tree_data->YLMin   = (int *)    DR_Array(LONG,      0,NbNmr-1);
    tree_data->YLMax   = (int *)    DR_Array(LONG,      0,NbNmr-1);
    tree_data->Yield   = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 
    tree_data->ZAInv   = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 
    tree_data->YCumP   = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 
    tree_data->MappSize= (long *)   DR_Array(LONG,      0,NbNmr-1);  



    if ((tree_data->NmrInv  == NULL) ||
        (tree_data->NmrDate == NULL) ||
        (tree_data->Ann1    == NULL) ||
        (tree_data->Ann0    == NULL) ||
        (tree_data->Zero1   == NULL) ||
        (tree_data->Zero0   == NULL) ||
        (tree_data->YLMin   == NULL) ||
        (tree_data->YLMax   == NULL) ||
        (tree_data->Yield   == NULL) ||
        (tree_data->ZAInv   == NULL) ||
        (tree_data->YCumP   == NULL) ||
        (tree_data->MappSize== NULL)  ) goto RETURN;

    /* Initialize pointers to NULL for safe freeing */
    for (i=0; i<NbNmr; i++)
    {
        tree_data->NmrInv[i] = NULL;
        tree_data->Ann1[i]   = NULL;
        tree_data->Ann0[i]   = NULL;
        tree_data->Zero1[i]  = NULL;
        tree_data->Zero0[i]  = NULL;
        tree_data->Yield[i]  = NULL;
        tree_data->ZAInv[i]  = NULL;
        tree_data->YCumP[i]  = NULL;
    }

    for (i=0; i<NbNmr; i++)
    {
        tree_data->NmrInv[i] = Alloc_Slice(tree_data);
        tree_data->Ann1[i]   = Alloc_Slice(tree_data);
        tree_data->Ann0[i]   = Alloc_Slice(tree_data);
        tree_data->Zero1[i]  = Alloc_Slice(tree_data);
        tree_data->Zero0[i]  = Alloc_Slice(tree_data);
        if ((tree_data->NmrInv[i] == NULL) ||
            (tree_data->Ann1[i]   == NULL) ||
            (tree_data->Ann0[i]   == NULL) ||
            (tree_data->Zero1[i]  == NULL) ||
            (tree_data->Zero0[i]  == NULL)) goto RETURN;

        tree_data->Yield[i] = (double *)DR_Array (DOUBLE,-1,Area);
        tree_data->ZAInv[i] = (double *)DR_Array (DOUBLE,-1,Area);
        tree_data->YCumP[i] = (double *)DR_Array (DOUBLE,-1,Area);
        if ((tree_data->Yield[i] == NULL) ||
            (tree_data->YCumP[i] == NULL) ||
            (tree_data->ZAInv[i] == NULL)) goto RETURN;
    }

    /* Copy numeraire dates into tree structure */
    tree_data->NbNmr = NbNmr;
    for (i=0; i < NbNmr; i++) tree_data->NmrDate[i] = mvd->NmrDate[i];

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("Nmr_Alloc: Failed!");
    }

    return (status);
}


#define Q3_NCK_5Q 2000002769
/*****  Nmr_Vanl **********************************************************/
/*
*       Calculate vanilla smile info for nmr dates.
*/
int Nmr_Vanl (long        ValueDate,       /* (I) Value date              */
              T_CURVE     *t_curve,        /* (I) Zero curve data         */
              MQDATA      *NmrMQ,          /* (I) MultiQ for vanilla      */
              MKTVOL_DATA *mvd,            /* (I/O) Mkt vol data          */
              TREE_DATA   *tree_data)      /* (I) Tree data               */
{
    int     status = FAILURE; /* Error status = FAILURE initially */  
    char    ErrorMsg[MAXBUFF];

    double  *SwapStT=NULL,*TenorT=NULL;    
    double  Expiry;
    long    TermDate;
    int     intvl, row, NbExpiry, NbNmr, Mat, MatL, MatR;
    int     i, j, s;
    int     CvDiff;    /* Diffuse curve number */
    

    TermDate = mvd->TermDate;
    NbExpiry = mvd->NbVol;
    NbNmr    = mvd->NbNmr;

    CvDiff   = tree_data->CvDiff;

    /* Create VNFM vols for nil calibration */
    if (!mvd->CalibFlag)
    {
        for (i=0; i<mvd->NbVol; i++)
        {
            if (IndexVol (&(mvd->Vol[0][i]),
                          mvd->SwapSt[i],
                          mvd->SwapMat[i],
                          mvd->Freq,
                          mvd->DCC,
                          mvd->CalibFlag,
                          mvd->NbVol,
                          mvd->BaseDate,
                          mvd->VolDate,
                          mvd->Aweight,
                          mvd->Bbq,
                          mvd->VolNorm,
                          mvd->VolLogn,
                          tree_data->NbFactor,
                          mvd->Beta,
                          t_curve[CvDiff].NbZero,
                          t_curve[CvDiff].Zero,
                          t_curve[CvDiff].ZeroDate,
                          t_curve[CvDiff].ValueDate) == FAILURE)
            {
                goto RETURN;
            }   
       }
    }

    /* Compute terminal zero */
    mvd->TermZero0 = ZeroPrice (TermDate,
                                ValueDate,
                                t_curve[CvDiff].NbZero,
                                t_curve[CvDiff].ZeroDate,
                                t_curve[CvDiff].Zero);

    /* Interpolate on time */
    SwapStT  = (double *) DR_Array(DOUBLE,0,NbExpiry-1);
    TenorT   = (double *) DR_Array(DOUBLE,0,NbExpiry-1);
    for (i=0; i<NbExpiry; i++)
    {
        SwapStT[i]  = Daysact (ValueDate, mvd->SwapSt[i])/365.;
        TenorT[i]   = (double) Months360(mvd->SwapSt[i], mvd->SwapMat[i]);
    }


    i = 0;
    /* Compute forward from curve  */
    if (Par_Yield_From_Dates (&(mvd->ParYield0[i]),
                              &(mvd->Annuity0[i]),
                              mvd->NmrSwapSt[i],
                              mvd->NmrSwapMat[i],
                              mvd->DCC,
                              mvd->Freq,
                              'F', 
                              t_curve[CvDiff].NbZero,
                              t_curve[CvDiff].Zero,
                              t_curve[CvDiff].ZeroDate,
                              ValueDate) == FAILURE)
    {
        goto RETURN;
    }
    for (s=0; s<NBVOLPARS; s++) mvd->NmrVol[i][s] = 0;

    i = NbNmr-1;
    mvd->ParYield0[i]  = 0;
    for (s=0; s<NBVOLPARS; s++) mvd->NmrVol[i][s] = 0;

    for (i = 1; i < NbNmr-1; i++) 
    {
        double  NmrDateT;
        double  Expiry;
        double  TotVol;
        double  NmrSmile[NBVOLPARS+4]; /* Extra space for NbSig,Nck,dN,tN */
        double  Delta[NBSTRIKE] = {0.05,0.25,0.5,0.75,0.99};

        /* Interpolate swap info on numeraire dates */
        NmrDateT = Daysact(ValueDate, mvd->NmrDate[i])/365.;
        for (s=0; s<NBVOLPARS; s++)
        {
            tableinterp(NmrDateT, 
                        &(mvd->NmrVol[i][s]), 
                        SwapStT,
                        mvd->Vol[s],
                        NbExpiry);
        }

        /* Compute forward from curve  */
        if (Par_Yield_From_Dates (&(mvd->ParYield0[i]),
                                  &(mvd->Annuity0[i]),
                                  mvd->NmrSwapSt[i], 
                                  mvd->NmrSwapMat[i],       
                                  mvd->DCC,           
                                  mvd->Freq,          
                                  'F',  
                                  t_curve[CvDiff].NbZero,        
                                  t_curve[CvDiff].Zero,         
                                  t_curve[CvDiff].ZeroDate,     
                                  ValueDate) == FAILURE)
        {
             goto RETURN;
        }

        /* Compute discount to numeraire date */
        mvd->StZero0[i] = ZeroPrice (mvd->NmrDate[i],
                                     ValueDate,
                                     t_curve[CvDiff].NbZero,
                                     t_curve[CvDiff].ZeroDate,
                                     t_curve[CvDiff].Zero);

        /* Clean up */
        for (s=0; s<NBSTRIKE; s++) mvd->Strike[i][s]=0;
        for (s=0; s<NBSTRIKE; s++) mvd->VanlPr[i][s]=0;
        for (s=0; s<NBSTRIKE; s++) mvd->TreePr[i][s]=0;

        /* Expiry */
        Expiry = Daysact(ValueDate, mvd->NmrDate[i]) / 365.;

        /* ATM price for nmr calib */
        mvd->AtmPrice0[i] = Call_BSQ (mvd->ParYield0[i],
                                      mvd->ParYield0[i],
                                      Expiry,
                                      mvd->NmrVol[i][0],
                                      1.); 

        /* MultiQ: fill numerical tree parameters and 
         * load all into smile vector */
        NmrSmile[NBVOLPARS]   = mvd->NbSigmaMQ;
        NmrSmile[NBVOLPARS+1] = mvd->NckMQ;
        NmrSmile[NBVOLPARS+2] = mvd->DeltaNMQ;
        NmrSmile[NBVOLPARS+3] = mvd->TauNMQ;
        for (s=0; s<NBVOLPARS; s++) NmrSmile[s] = mvd->NmrVol[i][s];

        /* Calibrate internal MultiQ distribution */
        if (Q3SVToMQ (mvd->ParYield0[i],
                      mvd->NmrVol[i][0],
                      Expiry,
                      NmrSmile+1,
                      &(NmrMQ[i])) == FAILURE)
        {
            DR_Error ("Calibration of MultiQ parameters failed "
                      "at date %d.\n", mvd->NmrDate[i]);
            goto RETURN;
        }

        /* Vanila stats */        
        if (mvd->NmrStatFlag)   
        {
            TotVol = mvd->NmrVol[i][0] * sqrt(Expiry);
            for (s=0; s<NBSTRIKE; s++)
            {
                mvd->Strike[i][s] = exp(Normal_InvH(Delta[s])*TotVol)
                                  * mvd->ParYield0[i];
                if (Q3MQPricer (&(NmrMQ[i]),
                                Q3_CALL,
                                mvd->Strike[i][s],
                                &(mvd->VanlPr[i][s])) == FAILURE)
                {
                    goto RETURN;
                }
            }
        }
    } /* for i < NbNmr-1 */

    status = SUCCESS;

RETURN:

    Nmr_Print1(mvd, NmrMQ);

    Free_DR_Array (SwapStT, DOUBLE, 0, NbExpiry-1);
    Free_DR_Array (TenorT,  DOUBLE, 0, NbExpiry-1);

    if (status == FAILURE)
    {
        DR_Error("Nmr_Schedule: Failed!");
    }

    return (status);

}  /* Nmr_Vanl */


/*****  Nmr_Calc  **********************************************************/
/*
*       TMX3: Numeraire calculation at HK dates 
*/
int Nmr_Calc (T_CURVE              *t_curve,       /* (I) Zero curve      */
              MKTVOL_DATA          *mvd,           /* (I) Volatility data */
              TREE_DATA            *tree_data)     /* (O) Tree data       */
{

    DEV_DATA    dev_data;         /* Dev data structure                    */
    CLAIM_BANK  ZeroBank;         /* Zeroes used in numeraire calculation  */

    MQDATA      *NmrMQ = NULL;    /* Local MultiQ structure for vanilla    */

    /* Slices */
    double      *Annuity  = NULL; /* Annuity slice                         */
    double      *LastZero = NULL; /* Longest maturity zero                 */
    double      *newZero  = NULL; /* New zero added to zerobank            */
    double      *NmrInv;          /* This is only pointer, NO ALLOCATION   */ 


    /* Numeraire date info */
    int         NbNmr;            /* Number of numeraire dates             */
    int         NmrIdx;           /* Numeraire date index                  */
    int         NmrFlag;          /* Is current date a numeraire date      */
    int         EDevIdx;          /* Express DEV date index                */

    /* Numerical variables */
    long        TerminalDate;
    long        CurrentDate;
    long        ErDate;
    long        ValueDate;
    int         IndexMat;         /* Swap maturity in months               */
    int         CvDiff;           /* Diffuse curve number                  */
    int         i, s;             /* Node indices                          */
    int         t;                /* Current time point                    */   
    int         T;                /* Last time point                       */
    int         status = FAILURE; /* Error status = FAILURE initially      */

    
    /* Dates and indices */
    T            = tree_data->NbTP;
    NbNmr        = tree_data->NbNmr;
    ValueDate    = t_curve->ValueDate;
    TerminalDate = tree_data->NmrDate[NbNmr-1];

    CvDiff       = tree_data->CvDiff;

    /* Initialize and allocate DEV_DATA structure */
    Dev_Init (&dev_data);
    if (Dev_Alloc (&dev_data, tree_data) == FAILURE) goto RETURN;

    /* Initialize and allocate zero bank */
    CbkInit(&ZeroBank);
    if (CbkAlloc (&ZeroBank, NbNmr, tree_data) == FAILURE) goto RETURN;

    /* Numeraire cannot be interpolated yet */
    tree_data->NmrInterpOn = FALSE;

    /* Numeraire is not needed in lattice */
    dev_data.NmrToCcy = FALSE;
    dev_data.CcyToNmr = FALSE;

    /* Allocate MQ structures */
    NmrMQ = (MQDATA *) calloc (NbNmr, sizeof(MQDATA));
    if (NmrMQ == NULL) goto RETURN;

    /* Calibrate MQ parameters, price vanilla if NmrStatFlag */
    if (Nmr_Vanl(ValueDate,
                 t_curve,
                 NmrMQ,
                 mvd,
                 tree_data) == FAILURE) goto RETURN;

    /* Allocate variables for numeraire calculation */
    if ((Annuity  = Alloc_Slice (tree_data)) == NULL)
    {
        DR_Error("Nmr_Calc: Could not allocate memory.\n");
        goto RETURN;
    }

    /* Step backward through the tree */
    for (t = T; t >= 0; t--)
    {
        /* Current Date */
        CurrentDate = tree_data->TPDate[t];
        
        /* Set flags */
        NmrFlag = tree_data->TPtype[NMREVENT][t];

        /* 'Update' tree */
        if (Lattice (&dev_data,
                     t,
                     T,
                     mvd,
                     tree_data) == FAILURE) goto RETURN;

        /* Ev Zero Bank; do not add any new zero yet */
        if (CbkDev (&ZeroBank,
                    CurrentDate,
                    t,
                    T,
                    CvDiff,
                    &dev_data,
                    tree_data) == FAILURE) goto RETURN;

        /* On numeraire dates, compute numeraire */
        if (NmrFlag)
        {
            NmrIdx = (int)(tree_data->CritDate[NMREVENT][t]).Value[0];
            if (NmrIdx < 0 || NmrIdx >= NbNmr) goto RETURN;

            NmrInv = NULL;

            /* Annuity and Last Zero slices for Nmr mapping */
            if (0 < NmrIdx && NmrIdx < NbNmr-1) 
            {
                IndexMat = Months360 (mvd->NmrSwapSt[NmrIdx],
                                      mvd->NmrSwapMat[NmrIdx]);

                if (ZbkAnnuity_t (Annuity,
                                  &ZeroBank,
                                  CurrentDate,
                                  mvd->NmrSwapSt[NmrIdx],
                                  IndexMat,
                                  mvd->DCC,
                                  mvd->Freq,
                                  t,
                                  tree_data) == FAILURE) goto RETURN;

                LastZero = ZbkReadZero (&ZeroBank,
                                        mvd->NmrSwapMat[NmrIdx],
                                        FALSE,
                                        CurrentDate,
                                        t,
                                        tree_data);
                if (LastZero == NULL) goto RETURN;
            }

            /* Calculate inverse numeraire */
            if (Nmr_t (&NmrInv,
                       Annuity,
                       LastZero,
                       t,
                       NmrIdx,
                       &(NmrMQ[NmrIdx]),
                       mvd,
                       tree_data) == FAILURE) goto RETURN;
           
            /* Add inverse numeraire as new zero to ZeroBank as in ZbkUpdate */
            newZero = CbkPopSlice(&ZeroBank);
            if (newZero == NULL) goto RETURN;

            if (Copy_Slice (newZero,
                            NmrInv, 
                            t,
                            tree_data) == FAILURE) goto RETURN;

            /* Find out earliest use date */
            ErDate = 99999999;
            for (i = 0; i < mvd->NbNmr; i++)
            {
                if (mvd->NmrSwapMat[i] >= CurrentDate)
                {
                    ErDate = MIN (ErDate, mvd->NmrSwapSt[i]);
                }
            }

            if (CbkPushSlice (&ZeroBank,
                              newZero,
                              CurrentDate,
                              ErDate) == FAILURE) goto RETURN;

            if (NmrIdx > 0 && NmrIdx < NbNmr-1)
            {
                /* Store current annuity and last zero in tree_data */
                if (Copy_Slice (tree_data->Ann1[NmrIdx],
                                Annuity,
                                t,
                                tree_data) == FAILURE) goto RETURN;

                LastZero = ZbkReadZero (&ZeroBank,
                                        mvd->NmrSwapMat[NmrIdx],
                                        FALSE,
                                        CurrentDate,
                                        t,
                                        tree_data);
                if (LastZero == NULL) goto RETURN;

                if (Copy_Slice (tree_data->Zero1[NmrIdx],
                                LastZero,
                                t,
                                tree_data) == FAILURE) goto RETURN;

                IndexMat = Months360 (mvd->NmrSwapSt[NmrIdx-1],
                                      mvd->NmrSwapMat[NmrIdx-1]);

                /* Store next annuity and last zero in tree_data */
                if (ZbkAnnuity_t (tree_data->Ann0[NmrIdx],
                                  &ZeroBank,
                                  CurrentDate,
                                  mvd->NmrSwapSt[NmrIdx-1],
                                  IndexMat,
                                  mvd->DCC,
                                  mvd->Freq,
                                  t,
                                  tree_data) == FAILURE) goto RETURN;

                LastZero = ZbkReadZero (&ZeroBank,
                                        mvd->NmrSwapMat[NmrIdx-1],
                                        FALSE,
                                        CurrentDate,
                                        t,
                                        tree_data);
                if (LastZero == NULL) goto RETURN;

                if (Copy_Slice (tree_data->Zero0[NmrIdx],
                                LastZero,
                                t,
                                tree_data) == FAILURE) goto RETURN;

            } /* if NmrIdx > 0 */
        } /* if (NmrFlag) */
    } /* for t */

    /* From now on it is possible to intepolate numeraire */
    tree_data->NmrInterpOn = TRUE;

    status = SUCCESS;

    RETURN:

    Nmr_Print2 (mvd, tree_data);

    if (NmrMQ != NULL) free(NmrMQ);

    Free_Slice (Annuity,   tree_data);
    Dev_Free   (&dev_data, tree_data);
    CbkFree    (&ZeroBank, tree_data); 

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_Calc: Failed!");
    }

    return (status);

}  /* Nmr_Calc */
       


#define   SIGMA_DELTA    0.0000001
#define   MAX_STEPS      20
#define   MAX_ERR        0.10
#define   CUTOFF_DELTA   0.99
/*****  Nmr_t   ***************************************************************/
/*
*       Generate numeraire slice at time t out of vanilla smile information
*       Nmr put directly into tree structure AND output the
*/
int     Nmr_t (  double      **NmrInv,    /* (O) Nmr pointer                */
                 double      *Annuity,    /* (I) Annuity for calib swap     */
                 double      *LastZero,   /* (I) Longest mat in annuity     */ 
                 int         t,           /* (I) Current time point         */
                 int         NmrIdx,      /* (I) Numeraire index            */
                 MQDATA      *mq,         /* (I) Vanilla data structure     */
                 MKTVOL_DATA *mvd,        /* (I) Volatility data structure  */
                 TREE_DATA   *tree_data)  /* (I/)) Tree data structure      */
{
    /* Local slice pointers */
    int     *LevelL;
    double  *AnnuityL;
    double  *LastZeroL; 
    double  *NmrInvL;      
    double  *StPriceL;

    /* Mapped swap yield information */
    double  *YieldL;              /* calibrating swap yield            */
    double  *ZAInvL;              /* dynamic lower bound for yield     */
    double  *YCumP_OUL;           /* OU cum prob up to given level     */ 
    double  *YProb_OU = NULL;     /* OU probability up to given level  */
    double  *YCumP_AA = NULL;     /* AA cum prob up to given level     */
    double  *YProb_AA = NULL;     /* AA probability of given level     */
    double  *YStrk = NULL;        /* corresponding point in Y space    */
    int     *Level = NULL;        /* node Y-level                      */

    /* 2-factor variables  */
    double *StPriceG  = NULL, 
           *AnnuityG  = NULL,
           *LastZeroG = NULL,
           *NmrInvG   = NULL;
    long   MappSize   = 0;


    /* Variables */
    double  ParYield0;
    double  AAtoOUFact; 
    double  AASlope;
    double  AASta[AA_ALL];
    
    /* NR variables */
    double  Call1, 
            Call2, 
            der, 
            diff,
            a = -0.10,  /* Initialization */ 
            target, 
            fwd, 
            HiYield, 
            Expiry, 
            TotVol,  
            TreeVol,
            VanlVol,
            ErrVol;
    long    CurrNode = 0,
            Iter = 0;

   
    /* Variables for binary swaption algorithm */
    double  binProb;
    double  prevPrice, currPrice;
    double  YMin = -9.999999;

    /* Limits and indices */
    int     Top1, Bottom1;                  /* Tree limits (1rst dim)      */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)       */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)       */ 
    int     LMinInit=0, LMaxInit=0;          /* Min/max yield   levels      */
    int     LMin=0, LMax=0;

    int     EDevIdx;                        /* Express DEV index           */
    int     i, j, k;                        /* Node indices                */
    int     offset;                         /* Node offset                 */
    long    Area = 0;                       /* Slices size                 */
    int     l;                              /* Level index                 */
    int     s;                              /* Stat index                  */
    int     status = FAILURE;               /* Error status                */


    *NmrInv = NULL;
    if (NmrIdx < 0 || NmrIdx >= tree_data->NbNmr) 
    {
        return (status);
    }

    /* Set Nmr pointer to current Nmr slice in tree_data */
    *NmrInv = tree_data->NmrInv[NmrIdx];

    
    /* Get state price slice index */
    EDevIdx = GetDLOffset (tree_data->NbEDevDates,
                           tree_data->EDevDate,
                           tree_data->TPDate[t],
                           CbkEXACT);
    if (EDevIdx < 0)
    {
        DR_Error ("State prices not available for current numeraire date %ld\n",
                  tree_data->TPDate[t]);
        return(status);
    }

    /* Slices size */
    if (tree_data->NbFactor == 1)
    {
        Area = tree_data->Width[0];
    }
    if (tree_data->NbFactor == 2)
    {
        Area = tree_data->Width[0] * tree_data->Width[1];
    }
    if (tree_data->NbFactor == 3)
    {
        Area = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];
    }

        
    /* Set mapping pointers to current Nmr slice in tree_data */
    YieldL     = tree_data->Yield[NmrIdx];
    ZAInvL     = tree_data->ZAInv[NmrIdx];
    YCumP_OUL  = tree_data->YCumP[NmrIdx];

    /* Clean up */
    for (s=0; s<AA_ALL; s++) AASta[s] = 0;
    for (s=-1;s<=Area;s++) ZAInvL[s] = -9.999999;

    /* Special Nmr dates: 0 and last */
    if (NmrIdx == tree_data->NbNmr-1) 
    {
        if (Set_Slice (*NmrInv,
                       1.,
                       t,
                       tree_data) == FAILURE)
        {
            goto RETURN;
        }
        /* Yield is undefined; later, will use the previous one;
         * however, at this point it has not been computed yet */
        return (SUCCESS);
    }
    else if (NmrIdx == 0)
    {
        if (Set_Slice (*NmrInv,
                       1. / mvd->TermZero0,
                       t,
                       tree_data) == FAILURE)
        {
            goto RETURN;
        }

        /* Mapping is constant */
        tree_data->YLMin[NmrIdx] = 0;
        tree_data->YLMax[NmrIdx] = 0;

        YieldL[0]    = mvd->ParYield0[0];
        YCumP_OUL[0] = 1.;

        return (SUCCESS);
    }
    
    /* Binary swaption algorithm */
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    /* Local variables */
    ParYield0  = mvd->ParYield0[NmrIdx];
    AAtoOUFact = mvd->TermZero0 / mvd->Annuity0[NmrIdx];

    /* Allocate level slice */
    Level = Alloc_Slice_Int (tree_data);
    if (Level == NULL)
    {
        DR_Error("Nmr_t: Could not allocate memory.\n");
        goto RETURN;
    }
   
    
    /*  1 FACTOR CASE */
    if (tree_data->NbFactor == 1)
    {
        offset    = Node_Offset(1, 0, 0, t, tree_data);
        LevelL    = Level    + offset;
        AnnuityL  = Annuity  + offset;
        LastZeroL = LastZero + offset;
        NmrInvL   = tree_data->NmrInv[NmrIdx]       + offset;
        StPriceL  = tree_data->EDevStPrice[EDevIdx] + offset;

        /* Assign level */
        for (i = Bottom1; i <= Top1; i++)
        {
            LevelL[i] = i - Bottom1;
            if (LevelL[i] < 0 || LevelL[i] > Area)
            {
                DR_Error ("Exceeded level array bounds\n");
                goto RETURN;
            }

            LMinInit = MIN(0,      LevelL[i]);
            LMaxInit = MAX(Area,   LevelL[i]);    

            /* update lower bound for yield */
            ZAInvL[LevelL[i]] = MAX (ZAInvL[LevelL[i]], 
                                     -LastZeroL[i]/AnnuityL[i]);
        }


        YCumP_AA = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        YProb_AA = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        YStrk    = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        YProb_OU = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);


        if ((YCumP_AA == NULL) ||
            (YProb_AA == NULL) ||
            (YStrk    == NULL) ||
            (YProb_OU == NULL)) goto RETURN;        

        for (l = LMinInit; l <= LMaxInit; l++) YProb_AA[l] = 0;

        /* Binaries */
        for (i = Bottom1; i <= Top1; i++)
        {
            YProb_AA[LevelL[i]] += StPriceL[i] * AAtoOUFact * AnnuityL[i];
            YProb_OU[LevelL[i]] += StPriceL[i];
        }

        /* Cumulative */
        YCumP_AA[LMinInit]  = YProb_AA[LMinInit];
        YCumP_OUL[LMinInit] = YProb_OU[LMinInit];
        for (l = LMinInit+1; l <= LMaxInit; l++) 
        {
            YCumP_AA[l]  = YCumP_AA[l-1]  + YProb_AA[l];
            YCumP_OUL[l] = YCumP_OUL[l-1] + YProb_OU[l];
        } 
        if (fabs(YCumP_AA[LMaxInit]-1.)>1e-6)
        {
            DR_Error ("Nmr_t: state prices don't add up to 1 at date %ld!",
                      tree_data->NmrDate[NmrIdx]);
            goto RETURN;
        }
        mvd->TreeSta[NmrIdx][0] = 1. - YCumP_AA[LMaxInit];
        
        /* search for cutoff levels based on cum probability */
        LMin = LMinInit;
        for (l=LMinInit; l<LMaxInit; l++)
        {
            if (YCumP_OUL[l]  < NormalH (-mvd->NbSigmaBinL)) LMin = l+1;
            else break;
        }

        LMax = LMaxInit;
        for (l=LMaxInit; l>LMinInit; l--)
        {
            if (YCumP_OUL[l] > NormalH (mvd->NbSigmaBinR)) LMax = l-1;
            else break;
        }

        /* Compute binary level */
        for (l = LMin; l <= LMax; l++)
        {
            if (Q3MQMap (mq,
                         mq->muMQ + mq->sigMQ * Normal_InvH (YCumP_AA[l]),
                         &YStrk[l]) == FAILURE)
            { 
                 DR_Error ("Could not compute MqMap\n");
                goto RETURN;
            }
        }

        /* Deal with negative yields: cutoff strike */
        for (i=Bottom1;i<=Top1;i++)
        {
            l = LevelL[i];

            if (YStrk[l]<0) 
            {
                YMin = MAX(YMin,-LastZeroL[i]/AnnuityL[i]);
            }
        }
        mvd->TreeSta[NmrIdx][4] = YMin;


        /* Find cutoff at 99 %-delta */
        Expiry  = Daysact(mvd->NmrDate[0], mvd->NmrDate[NmrIdx]) / 365.0;
        TotVol  = mvd->NmrVol[NmrIdx][0] * sqrt(Expiry);
        HiYield = exp(Normal_InvH(CUTOFF_DELTA)*TotVol) * ParYield0;
            

        /* Yield on the given level */
        prevPrice = mq->fwdRate;
        for (l = LMin; l <= LMax-1; l++)
        {
            /* Call at new YStrk point */
            if (Q3MQPricer (mq,
                            Q3_CALL,
                            YStrk[l],
                            &currPrice) == FAILURE)
            {
                goto RETURN;
            }

            /* Add up binary */
            currPrice += YStrk[l] * (YCumP_AA[LMaxInit] - YCumP_AA[l]);

            /* Finally, assign yield value on this level */
            binProb = (l == LMin ? YCumP_AA[l] : YCumP_AA[l] - YCumP_AA[l-1]);
            if (binProb < NMR_CUTOFF)
            {
                YieldL[l] = YieldL[l-1];
            }
            else
            {
                YieldL[l] = (prevPrice - currPrice) / binProb;
            }

            if (YieldL[l] > HiYield)
            {
                break;
            }

            /* Store stuff for reuse at the next level */
            prevPrice = currPrice;

            CurrNode = l;
        } 

        /* Flat left Yield extrapolation  */
        for (l=LMinInit; l<LMin; l++) YieldL[l] = YieldL[LMin];
        
        /* Semi-parametric TMX : right tail extrapolation */
        target  = mvd->VanlPr[NmrIdx][NBSTRIKE-1];              
        do
        {
            /* Initialization */
            Call1 = 0.0;
            Call2 = 0.0; 

            /* 99-delta Call */ 
            for (l = LMin; l < CurrNode; l++)
            {
                Call1      += YProb_AA[l] * MAX(YieldL[l] - mvd->Strike[NmrIdx][NBSTRIKE-1], 0);
            }
            
            for (l = CurrNode; l <= LMaxInit; l++)
            {
                /* Functional form */
                YieldL[l]   = YieldL[CurrNode] * (1.0 + (  exp(a * (l - CurrNode)) - 1.0 )   /a );
                Call1      += YProb_AA[l] * MAX(YieldL[l] - mvd->Strike[NmrIdx][NBSTRIKE-1], 0);
            }

            /* Bias */
            diff = Call1 - target;


            /* Convergence reached */
            if (fabs(diff) < TINY)
                break;

            a += SIGMA_DELTA;

            /* Find Gradient */
            for (l = LMin; l < CurrNode; l++)
            {
                Call2      += YProb_AA[l] * MAX(YieldL[l] - mvd->Strike[NmrIdx][NBSTRIKE-1], 0);
            }
            
            for (l = CurrNode; l <= LMaxInit; l++)
            {
                YieldL[l]   = YieldL[CurrNode] * (1.0 + (  exp(a * (l - CurrNode)) - 1.0 )   /a );
                Call2      += YProb_AA[l] * MAX(YieldL[l] - mvd->Strike[NmrIdx][NBSTRIKE-1], 0);
            }

            /* Gradient */
            der = (Call2 - Call1) / SIGMA_DELTA;

            
            /* Newton-Raphson */
            if (fabs(der) < TINY)
            {
                a += SIGMA_DELTA;
            }
            else
            { 
                a -= diff / der;
            }

            /*  Incremented Iterations */
            Iter++;

        }
        while (Iter < MAX_STEPS);


        /* Now check convergence */
        VanlVol = ImpVol_BSQ(ParYield0,
                             mvd->Strike[NmrIdx][NBSTRIKE-1],
                             Expiry,
                             target,
                             'C',
                             1, /*  1, 0, */
                             mvd->NmrVol[NmrIdx][0]);
        
        TreeVol = ImpVol_BSQ(ParYield0,
                             mvd->Strike[NmrIdx][NBSTRIKE-1],
                             Expiry,
                             ((Iter < MAX_STEPS) ? Call1 : Call2),
                             'C',
                             1, /*  1, 0, */
                             mvd->NmrVol[NmrIdx][0]);
        
        /* Calibration error */
        ErrVol = fabs(TreeVol-VanlVol);

        /* Check ErrVol is less than MAX_ERR */
        if (fabs(ErrVol) > MAX_ERR)
        {
            DR_Error("Nmr_t : Right tail extrapolation failed!");
            goto RETURN;
        }

        /* Now make sure the forward is exactly repriced */
        fwd  = 0.0; 
        for (l=LMinInit; l <= LMaxInit; l++)
        {
            fwd += YProb_AA[l] * YieldL[l];
        }

        /* Adjustment to reprice the forward */
        for (l=LMinInit; l <= LMaxInit; l++)
        {
            YieldL[l] *= ParYield0 / fwd;
        }

        
        /* Record min and max levels */
        tree_data->YLMin[NmrIdx] = LMin;
        tree_data->YLMax[NmrIdx] = LMax;

        /* Assign 0 cumulative probability below LMin */
        YCumP_OUL[LMin-1] = 0.;

        /* Compute AA tree price stats */
        for (i = Bottom1; i <= Top1; i++)
        {
            double AStP = AnnuityL[i] * StPriceL[i];
            
            AASta[AA_FWD] += AStP * YieldL[LevelL[i]];
            AASta[AA_ATM] += AStP * MAX (YieldL[LevelL[i]] - ParYield0, 0);
        }
        AASta[AA_FWD] *= AAtoOUFact;
        AASta[AA_ATM] *= AAtoOUFact;
        /* Store in MKTVOL_DATA */
        mvd->AASta[NmrIdx][AA_FWD] = AASta[AA_FWD];
        mvd->AASta[NmrIdx][AA_ATM] = AASta[AA_ATM];

        /* Calibrate ATM price and store VolC */
        AASlope = mvd->AtmPrice0[NmrIdx] / AASta[AA_ATM]; 
        mvd->AASta[NmrIdx][AA_VOL] = AASlope;
        for (l = LMinInit; l <= LMaxInit; l++)
        {
            YieldL[l] = ParYield0 + AASlope * (YieldL[l] - ParYield0);
        }
  
        /* Compute inverse numeraire slice */
        mvd->TreeSta[NmrIdx][5] = 0.;
        for (i = Bottom1; i <= Top1; i++)
        {
            NmrInvL[i] = LastZeroL[i] + YieldL[LevelL[i]] * AnnuityL[i];

            /* Store prob of negative numeraire in mvd BEFORE adjusting */
            if (NmrInvL[i] < NMR_CUTOFF)
            {
                mvd->TreeSta[NmrIdx][5] += StPriceL[i]; 
                DR_Error ("Negative numeraire at date %ld",
                          tree_data->NmrDate[NmrIdx]);
                goto RETURN; 
            }
        } 

        /* Print debug info */
        if (Nmr_Stat  (t,
                       NmrIdx,
                       Annuity,
                       Level,
                       YieldL,
                       mvd,
                       tree_data) == FAILURE) goto RETURN;

    } /* 1-factor */


    /* 2 FACTOR CASE */
    else if (tree_data->NbFactor == 2)
    {
        
        /* Allocate memory */
        StPriceG    = (double *) DR_Array (DOUBLE, 0, Area);
        AnnuityG    = (double *) DR_Array (DOUBLE, 0, Area);
        LastZeroG   = (double *) DR_Array (DOUBLE, 0, Area);
        NmrInvG     = (double *) DR_Array (DOUBLE, 0, Area);    
        
        /* Project the 2D slices on a 1D grid */
        if ( SortTree_Slices(AnnuityG,      
                             LastZeroG,  
                             NULL,
                             NULL,
                             StPriceG,
                             &MappSize,
                             2,          /* 2 slices to project */
                             Annuity,
                             LastZero,
                             NULL,
                             NULL,
                             tree_data->NmrMap_dir,   /* Mapping direction   */
                             EDevIdx,
                             t,
                             tree_data) == FAILURE)
        {
            goto RETURN;
        }

        /* Store the MappSize for interpolation */
        tree_data->MappSize[NmrIdx] = MappSize;

        /* Assign level */
        for (i = 0; i < MappSize; i++)
        {
            
            ZAInvL[i] = MAX (ZAInvL[i], 
                            -LastZeroG[i]/AnnuityG[i]);
        }

        
        /* Initialize cuotff parameters */      
        YCumP_AA = (double *) DR_Array (DOUBLE, 0, MappSize);
        YProb_AA = (double *) DR_Array (DOUBLE, 0, MappSize);
        YStrk    = (double *) DR_Array (DOUBLE, 0, MappSize);
        YProb_OU = (double *) DR_Array (DOUBLE, 0, MappSize);
        

        if ((YCumP_AA == NULL) ||
            (YProb_AA == NULL) ||
            (YStrk    == NULL) ||
            (YProb_OU == NULL)) goto RETURN;        

        for (l =0; l < MappSize; l++) YProb_AA[l] = 0;
        
        /* Binaries */
        for (i = 0; i < MappSize; i++)
        {
            YProb_AA[i] += StPriceG[i] * AAtoOUFact * AnnuityG[i];
            YProb_OU[i] += StPriceG[i];
        }
        
        /* Cumulative */
        YCumP_AA[0]  = YProb_AA[0];
        YCumP_OUL[0] = YProb_OU[0];
        for (l = 1; l < MappSize; l++) 
        {
            YCumP_AA[l]  = YCumP_AA[l-1]  + YProb_AA[l];
            YCumP_OUL[l] = YCumP_OUL[l-1] + YProb_OU[l];
        } 
        if (fabs(YCumP_AA[MappSize - 1]-1.)>1e-6)
        {
            DR_Error ("Nmr_t: state prices don't add up to 1 at date %ld!",
                      tree_data->NmrDate[NmrIdx]);
            goto RETURN;
        }
        mvd->TreeSta[NmrIdx][0] = 1. - YCumP_AA[MappSize - 1];
        
        /* search for cutoff levels based on cum probability */
        LMin = 0;
        for (l = 0; l < MappSize - 1; l++)
        {
            if (YCumP_OUL[l]  < NormalH (-mvd->NbSigmaBinL)) LMin = l+1;
            else break;
        }

        LMax = MappSize - 1;
        for (l = MappSize-1; l > 0; l--)
        {
            if (YCumP_OUL[l] > NormalH (mvd->NbSigmaBinR)) LMax = l-1;
            else break;
        }

        /* Compute binary level */
        for (l = LMin; l <= LMax; l++)
        {
            if (Q3MQMap (mq,
                         mq->muMQ + mq->sigMQ * Normal_InvH (YCumP_AA[l]),
                         &YStrk[l]) == FAILURE)
            { 
                 DR_Error ("Could not compute MqMap\n");
                goto RETURN;
            }
        }

        /* Deal with negative yields: cutoff strike */
        for (i = 0; i < MappSize ; i++)
        {
            if (YStrk[i]<0) 
            {
                YMin = MAX(YMin,-LastZeroG[i]/AnnuityG[i]);
            }
        }
        mvd->TreeSta[NmrIdx][4] = YMin;

         /* Find cutoff at 95 %-delta */
        Expiry  = Daysact(mvd->NmrDate[0], mvd->NmrDate[NmrIdx]) / 365.0;
        TotVol  = mvd->NmrVol[NmrIdx][0] * sqrt(Expiry);
        HiYield = exp(Normal_InvH(CUTOFF_DELTA)*TotVol) * ParYield0;


        /* Yield on the given level */
        for (l = LMin; l <= LMax-1; l++)
        {
            YieldL[l] = YStrk[l];
            
            if ((l > LMin) && (YieldL[l] < YieldL[l-1]))
            {
                YieldL[l] = YieldL[l-1];
            }

            
            if (YieldL[l] > HiYield)
            {
                break;
            }

            
            CurrNode = l;
        } 

        /* Flat left Yield extrapolation  */
        for (l = 0; l<LMin; l++) YieldL[l] = YieldL[LMin];
        
        /* Semi-parametric TMX : right tail extrapolation */
        target  = mvd->VanlPr[NmrIdx][NBSTRIKE-1];              
        do
        {
            /* Initialization */
            Call1 = 0.0;
            Call2 = 0.0; 

            /* 99-delta Call */ 
            for (l = LMin; l < CurrNode; l++)
            {
                Call1      += YProb_AA[l] * 
                              MAX(YieldL[l] - mvd->Strike[NmrIdx][NBSTRIKE-1], 0);
            }
            
            for (l = CurrNode; l < MappSize; l++)
            {
                /* Functional form */
                YieldL[l]   = YieldL[CurrNode] * 
                              (1.0 + (  exp(a * (l - CurrNode)) - 1.0 )   /a );
                Call1      += YProb_AA[l] * 
                               MAX(YieldL[l] - mvd->Strike[NmrIdx][NBSTRIKE-1], 0);
            }

            /* Bias */
            diff = Call1 - target;


            /* Convergence reached */
            if (fabs(diff) < TINY)
                break;

            a += SIGMA_DELTA;

            /* Find Gradient */
            for (l = LMin; l < CurrNode; l++)
            {
                Call2      += YProb_AA[l] * MAX(YieldL[l] - mvd->Strike[NmrIdx][NBSTRIKE-1], 0);
            }
            
            for (l = CurrNode; l < MappSize; l++)
            {
                YieldL[l]   = YieldL[CurrNode] * 
                              (1.0 + (  exp(a * (l - CurrNode)) - 1.0 )   /a );
                Call2      += YProb_AA[l] * 
                              MAX(YieldL[l] - mvd->Strike[NmrIdx][NBSTRIKE-1], 0);
            }

            /* Gradient */
            der = (Call2 - Call1) / SIGMA_DELTA;

            
            /* Newton-Raphson */
            if (fabs(der) < TINY)
            {
                a += SIGMA_DELTA;
            }
            else
            { 
                a -= diff / der;
            }

            /*  Incremented Iterations */
            Iter++;

        }
        while (Iter < MAX_STEPS);


        /* Now check convergence */
        VanlVol = ImpVol_BSQ(ParYield0,
                             mvd->Strike[NmrIdx][NBSTRIKE-1],
                             Expiry,
                             target,
                             'C',
                             1, /*  1, 0, */
                             mvd->NmrVol[NmrIdx][0]);
        
        TreeVol = ImpVol_BSQ(ParYield0,
                             mvd->Strike[NmrIdx][NBSTRIKE-1],
                             Expiry,
                             ((Iter < MAX_STEPS) ? Call1 : Call2),
                             'C',
                             1, /*  1, 0, */
                             mvd->NmrVol[NmrIdx][0]);
        
        /* Calibration error */
        ErrVol = fabs(TreeVol-VanlVol);

        /* Check ErrVol is less than MAX_ERR */
        if (fabs(ErrVol) > MAX_ERR)
        {
            DR_Error("Nmr_t : Right tail extrapolation failed!");
            goto RETURN;
        }

        /* Now make sure the forward is exactly repriced */
        fwd  = 0.0; 
        for (l = 0; l < MappSize; l++)
        {
            fwd += YProb_AA[l] * YieldL[l];
        }

        /* Adjustment to reprice the forward */
        for (l = 0; l < MappSize; l++)
        {
            YieldL[l] *= ParYield0 / fwd;
        }

        
        /* Record min and max levels */
        tree_data->YLMin[NmrIdx] = LMin;
        tree_data->YLMax[NmrIdx] = LMax;

        /* Assign 0 cumulative probability below LMin */
        YCumP_OUL[-1] = 0.;

        /* Compute AA tree price stats */
        for (i = 0; i < MappSize; i++)
        {
            double AStP = AnnuityG[i] * StPriceG[i];
            
            AASta[AA_FWD] += AStP * YieldL[i];
            AASta[AA_ATM] += AStP * MAX (YieldL[i] - ParYield0, 0);
        }
        AASta[AA_FWD] *= AAtoOUFact;
        AASta[AA_ATM] *= AAtoOUFact;
        /* Store in MKTVOL_DATA */
        mvd->AASta[NmrIdx][AA_FWD] = AASta[AA_FWD];
        mvd->AASta[NmrIdx][AA_ATM] = AASta[AA_ATM];

        /* Calibrate ATM price and store VolC */
        AASlope = mvd->AtmPrice0[NmrIdx] / AASta[AA_ATM]; 
        mvd->AASta[NmrIdx][AA_VOL] = AASlope;
        for (l = 0; l < MappSize; l++)
        {
            YieldL[l] = ParYield0 + AASlope * (YieldL[l] - ParYield0);
        }
  
        /* Compute inverse numeraire slice */
        mvd->TreeSta[NmrIdx][5] = 0.;
        for (i = 0; i < MappSize; i++)
        {
            NmrInvG[i] = LastZeroG[i] + YieldL[i] * AnnuityG[i];

            /* Store prob of negative numeraire in mvd BEFORE adjusting */
            if (NmrInvG[i] < NMR_CUTOFF)
            {
                mvd->TreeSta[NmrIdx][5] += StPriceG[i]; 
                DR_Error ("Negative numeraire at date %ld",
                          tree_data->NmrDate[NmrIdx]);
                goto RETURN; 
            }
        } 

        /* Now finally populate the NmrInv slice */
        if (InvMap(*NmrInv,
                   NmrInvG,
                   MappSize,
                   tree_data->NmrMap_dir,  /* Mapping direction */
                   EDevIdx,
                   t,
                   tree_data) == FAILURE)
        {
            goto RETURN;
        }


        
        /* Print debug info */
        if (Nmr_Stat  (t,
                       NmrIdx,
                       Annuity,
                       NULL,
                       YieldL,
                       mvd,
                       tree_data) == FAILURE) goto RETURN;


    }   
    else if (tree_data->NbFactor == 3)
    {
        goto RETURN;
    }    

    status = SUCCESS;

    RETURN:

    Free_Slice_Int (Level, tree_data);

    if (tree_data->NbFactor == 1)
    {
        Free_DR_Array (YProb_AA, DOUBLE, LMinInit-1, LMaxInit);
        Free_DR_Array (YCumP_AA, DOUBLE, LMinInit-1, LMaxInit);
        Free_DR_Array (YStrk,    DOUBLE, LMinInit-1, LMaxInit);
        Free_DR_Array (YProb_OU, DOUBLE, LMinInit-1, LMaxInit);
    }    
    if (tree_data->NbFactor == 2)
    {
        
        Free_DR_Array (YProb_AA, DOUBLE, 0, MappSize);
        Free_DR_Array (YCumP_AA, DOUBLE, 0, MappSize);
        Free_DR_Array (YStrk,    DOUBLE, 0, MappSize);
        Free_DR_Array (YProb_OU, DOUBLE, 0, MappSize);

        Free_DR_Array(StPriceG, DOUBLE, 0, Area);
        Free_DR_Array(AnnuityG, DOUBLE, 0, Area);
        Free_DR_Array(LastZeroG,DOUBLE, 0, Area);
        Free_DR_Array(NmrInvG,  DOUBLE, 0, Area);
       
    }

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_t: Failed!");
    }

    return (status);

} /* Nmr_t */


/*****  Nmr_Interp ***********************************************************/
/*
 * Compute interpolated numeraire in-between numeraire dates
 */
int     Nmr_Interp (double      *NmrInv,     /* (O) Nmr pointer             */
                    int         NmrIdx,      /* (I) Numeraire index         */
                    int         t,           /* (I) Current time point      */
                    MKTVOL_DATA *mktvol_data,/* (I) Volatility data struct  */
                    TREE_DATA   *tree_data)  /* (I) Tree data structure     */
{
    /* Local pointers */
    double  *Ann1L;        
    double  *Ann0L;
    double  *Zero1L;
    double  *Zero0L;
    double  *YCumP1L;
    double  *Yield1L;
    double  *YCumP0L;
    double  *Yield0L;
    double  *NmrInvL;
    double  *StPriceL;
    int     *LevelL;

    /* Used slices */
    double  *Yield      = NULL;          /* Yield values, current slice     */
    double  *YCumP      = NULL;          /* Probability on given level      */
    double  *YProb      = NULL;          /* Probability on given level      */
    int     *Level      = NULL;          /* Yield level, current slice      */


    int     EDevIdx;                     /* Express DEV date index          */
    double  tFrac;                       /* Fract within interval btw nmr   */
    double  NmrInv0, NmrInvAdj;          /* Exp & adj of inverse nmr        */
    double  probCutoff;
    double  prevProb, nextProb;
    double  prevCumP, nextCumP;   
    
    long    PrevT, NextT;
    double  PrevYield,NextYield;

    /* 2-factor interpolation */
    long    MappSize  = 0;
    double *StPriceG  = NULL, 
           *Ann1G     = NULL,
           *Zero1G    = NULL,
           *Ann0G     = NULL,
           *Zero0G    = NULL,
           *NmrInvG   = NULL;


    int     Top1, Bottom1;               /* Tree limits (1rst dim)          */
    int     *Top2, *Bottom2;             /* Tree limits (2nd dim)           */
    int     **Top3, **Bottom3;           /* Tree limits (3rd dim)           */
    int     OutTop1, OutBottom1;         /* Outer tree limits               */
    int     *OutTop2, *OutBottom2;
    int     **OutTop3, **OutBottom3;

    int     CvDiff;                      /* Diffuse curve number            */

    int     i, j, k;                     /* Node indices                    */
    int     offset;                      /* Node offset                     */

    int     LMinInit=0, LMaxInit=0;
    int     LMin,LMax;

    long    Area = 0;                    /* Slices size                     */

    int     l;                           /* Yield level index               */
    int     status = FAILURE;            /* Error status                    */

    long    CurrentDate = tree_data->TPDate[t];

    
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    CvDiff = tree_data->CvDiff;
 
    /* Allocate level slice */
    Level = Alloc_Slice_Int (tree_data);
    if (Level == NULL)
    {
        DR_Error("Nmr_t: Could not allocate memory.\n");
        goto RETURN;
    }

    /* Compute fraction within interval between numeraires */
    tFrac  = (double)
        Daysact (tree_data->NmrDate[NmrIdx], CurrentDate) /
        Daysact (tree_data->NmrDate[NmrIdx], 
                 tree_data->NmrDate[NmrIdx-1]);

    /* Identify state price slice corresponding to current TP */
    EDevIdx = GetDLOffset (tree_data->NbEDevDates,
                           tree_data->EDevDate,
                           tree_data->TPDate[t],
                           CbkEXACT);
    if (EDevIdx < 0)
    {
        DR_Error ("State prices not available for current date %ld\n", 
                  CurrentDate);
        goto RETURN;
    }
    
    /* Find previous numeraire date time */
    PrevT = GetDLOffset (tree_data->NbTP, 
                         tree_data->TPDate,
                         tree_data->NmrDate[NmrIdx-1],
                         CbkEXACT);
    
    if (tree_data->TPDate[PrevT] != tree_data->NmrDate[NmrIdx-1])
    {
        DR_Error ("State prices not available for next Nmr date %ld\n", 
                  tree_data->TPDate[t]);
        goto RETURN;
    }
    
    /* Find next numeraire date time */
    NextT = GetDLOffset (tree_data->NbTP,        
                         tree_data->TPDate,
                         tree_data->NmrDate[NmrIdx],
                         CbkEXACT);
    if (tree_data->TPDate[NextT] != tree_data->NmrDate[NmrIdx])
    {
        DR_Error ("State prices not available for next Nmr date %ld\n", 
            tree_data->TPDate[t]);
        goto RETURN;
    }

    /* Slices size */
    if (tree_data->NbFactor == 1)
    {
        Area = tree_data->Width[0];
    }
    if (tree_data->NbFactor == 2)
    {
        Area = tree_data->Width[0] * tree_data->Width[1];
    }
    if (tree_data->NbFactor == 3)
    {
        Area = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];
    }
    
    
    /* 1D arrays */
    YCumP1L = tree_data->YCumP[NmrIdx];
    Yield1L = tree_data->Yield[NmrIdx];
    YCumP0L = tree_data->YCumP[NmrIdx-1];
    Yield0L = tree_data->Yield[NmrIdx-1];
     


    /* 1 FACTOR CASE */
    if (tree_data->NbFactor == 1)
    {
        /* Assign local pointers */
        offset   = Node_Offset(1, 0, 0, t, tree_data);
        StPriceL = tree_data->EDevStPrice[EDevIdx] + offset;
        Ann1L    = tree_data->Ann1[NmrIdx]      + offset;
        Zero1L   = tree_data->Zero1[NmrIdx]     + offset;
        Ann0L    = tree_data->Ann0[NmrIdx]      + offset;
        Zero0L   = tree_data->Zero0[NmrIdx]     + offset;
        NmrInvL  = NmrInv + offset;
        LevelL   = Level  + offset;
        
       
        /* Assign level */
        for (i = Bottom1; i <= Top1; i++)
        {
            LevelL[i] = i - Bottom1;
            if (LevelL[i] < 0 || LevelL[i] > Area)
            {
                DR_Error ("Exceeded level array bounds\n");
                goto RETURN;
            }
            
            LMinInit = MIN(0,        LevelL[i]);
            LMaxInit = MAX(Area    , LevelL[i]);
        }
        
        /* Allocate numeraire and yield slices */
        Yield      = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        YCumP      = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        YProb      = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        if ((Yield  == NULL) ||
            (YCumP  == NULL) ||
            (YProb  == NULL))
        {
            DR_Error("Nmr_Interp: Could not allocate memory.\n");
            goto RETURN;
        }
        
        /* Binaries */
        for (i = Bottom1; i <= Top1; i++)
        {
            YProb[LevelL[i]] += StPriceL[i];
        }
        
        /* Cumulative */
        YCumP[LMinInit]   = YProb[LMinInit];
        for (l = LMinInit+1; l <= LMaxInit; l++) 
        {
            YCumP[l] = YCumP[l-1] + YProb[l];
        }
        
        if (fabs(YCumP[LMaxInit]-1.)>1e-6)
        {
            DR_Error ("Nmr_t: state prices don't add up to 1! ");
            goto RETURN;
        }
        
        /* Seek for minimal and maximal cumprobs */
        probCutoff = NormalH (- tree_data->NbSigmaMax);
        LMin = LMinInit;
        for (l=LMinInit; l<=LMaxInit; l++)
        {
            if (YCumP[l] < probCutoff) LMin = l+1;
            else break;
        }
        
        LMax = LMaxInit;
        for (l=LMaxInit; l>=LMinInit; l--)
        {
            if (YCumP[l] > 1-probCutoff) LMax = l-1;
            else break;
        }
                
        for (l=LMin; l<=LMax; l++)
        {
            /* Interpolate Previous and next numeraire array */
            /* Compute fraction within interval between numeraires */
            if (NmrIdx == 1)
            {
                PrevYield = mktvol_data->ParYield0[0];
            }
            else
            {
                tableinterp (YCumP[l],
                             &PrevYield,
                             YCumP0L,
                             Yield0L,
                             tree_data->Top1[PrevT]-tree_data->Bottom1[PrevT]+1);
            }
            
            
            tableinterp (YCumP[l],
                         &NextYield,
                         YCumP1L,
                         Yield1L,
                         tree_data->Top1[NextT]-tree_data->Bottom1[NextT]+1);            
            
            /* interpolated yield */
            Yield[l] = (1.0-tFrac) * NextYield  + 
                        tFrac      * PrevYield;
        }
        
        /* Extend flatly beyond LMin and LMax */
        for (l=LMinInit; l<LMin; l++) Yield[l] = Yield[LMin];
        for (l=LMaxInit; l>LMax; l--) Yield[l] = Yield[LMax];
        
        /* Create inverse numeraire slices */
        NmrInv0 = 0;
        for (i = Bottom1; i <= Top1; i++)
        {
            NmrInvL[i] = 
                    (1.0-tFrac) * (Yield[LevelL[i]] * Ann1L[i] + Zero1L[i]) +
                     tFrac      * (Yield[LevelL[i]] * Ann0L[i] + Zero0L[i]);
            
            /* ensure positivity */
            NmrInvL[i] = MAX (NmrInvL[i], NMR_CUTOFF);
            
            /* update expectation */
            NmrInv0    += NmrInvL[i] * StPriceL[i];
        }
        
        NmrInvAdj = 1/ (tree_data->TermZero[CvDiff][t] * NmrInv0);
        
        /* Readjust numeraire */
        for (i = Bottom1; i <= Top1; i++)
        {
            NmrInvL[i] *= NmrInvAdj;
        }
        
        if (Nmr_Print3 (Yield,
                        YCumP,
                        LMin,
                        LMax,
                        NmrInvAdj,
                        NmrIdx,
                        t,
                        mktvol_data,
                        tree_data) == FAILURE) goto RETURN;

    }

    /* 2 FACTOR CASE */
    else if (tree_data->NbFactor == 2)
    {
        
         /* Allocate memory */
        StPriceG    = (double *) DR_Array (DOUBLE, 0, Area);
        Ann1G       = (double *) DR_Array (DOUBLE, 0, Area);
        Zero1G      = (double *) DR_Array (DOUBLE, 0, Area);
        Ann0G       = (double *) DR_Array (DOUBLE, 0, Area);
        Zero0G      = (double *) DR_Array (DOUBLE, 0, Area);
        NmrInvG     = (double *) DR_Array (DOUBLE, 0, Area);
        

        /* Project the 2D slices on a 1D grid */
        if ( SortTree_Slices(Ann1G,      
                             Ann0G,  
                             Zero1G,
                             Zero0G,
                             StPriceG,
                             &MappSize,
                             4,                    /* 4 Slices to project */
                             tree_data->Ann1[NmrIdx],
                             tree_data->Ann0[NmrIdx],
                             tree_data->Zero1[NmrIdx],
                             tree_data->Zero0[NmrIdx],
                             tree_data->NmrMap_dir,       /* Mapping direction */
                             EDevIdx,
                             t,
                             tree_data) == FAILURE)
        {
            goto RETURN;
        }   

                
        /* Allocate numeraire and yield slices */
        Yield      = (double *) DR_Array (DOUBLE, 0, MappSize);
        YCumP      = (double *) DR_Array (DOUBLE, 0, MappSize);
        YProb      = (double *) DR_Array (DOUBLE, 0, MappSize);
        
        if ((Yield  == NULL) ||
            (YCumP  == NULL) ||
            (YProb  == NULL))
        {
            DR_Error("Nmr_Interp: Could not allocate memory.\n");
            goto RETURN;
        }
        
        /* Binaries */
        for (i = 0; i < MappSize; i++)
        {
            YProb[i] = StPriceG[i];
        }
        
        /* Cumulative */
        YCumP[0]   = YProb[0];
        for (l = 1; l < MappSize; l++) 
        {
            YCumP[l] = YCumP[l-1] + YProb[l];
        }
        
        if (fabs(YCumP[MappSize -1]-1.)>1e-6)
        {
            DR_Error ("Nmr_Interp: state prices don't add up to 1! ");
            goto RETURN;
        }
        
        /* Seek for minimal and maximal cumprobs */
        probCutoff = NormalH (- tree_data->NbSigmaMax);
        
        LMin = 0;
        for (l = 0; l < MappSize; l++)
        {
            if (YCumP[l] < probCutoff) LMin = l+1;
            else break;
        }
        
        LMax = MappSize -1;
        for (l = MappSize -1  ; l>= 0; l--)
        {
            if (YCumP[l] > 1-probCutoff) LMax = l-1;
            else break;
        }


        /* Quadratic interpolation of the yield */
        for (l=LMin; l<=LMax; l++)
        {
            /* Interpolate Previous and next numeraire array */
            /* Compute fraction within interval between numeraires */
            if (NmrIdx == 1)
            {
                PrevYield = mktvol_data->ParYield0[0];
            }
            else
            {
                tableinterp (YCumP[l],
                             &PrevYield,
                             YCumP0L,
                             Yield0L,
                             tree_data->MappSize[NmrIdx-1]);
                
            }
            
            tableinterp (YCumP[l],
                         &NextYield,
                         YCumP1L,
                         Yield1L,
                         tree_data->MappSize[NmrIdx]);            
            
            
            /* interpolated yield */
            Yield[l] = (1.0-tFrac) * NextYield  + 
                        tFrac      * PrevYield;
        }
        
        /* Extend flatly beyond LMin and LMax */
        for (l = 0         ; l < LMin; l++) Yield[l] = Yield[LMin];
        for (l = MappSize-1; l > LMax; l--) Yield[l] = Yield[LMax];


        /* Create inverse numeraire slices */
        NmrInv0 = 0;
        for (i = 0; i < MappSize; i++)
        {
            NmrInvG [i] = 
                    (1.0-tFrac) * (Yield[i] * Ann1G[i] + Zero1G[i]) +
                     tFrac      * (Yield[i] * Ann0G[i] + Zero0G[i]);
            
            /* ensure positivity */
            NmrInvG[i] = MAX (NmrInvG[i], NMR_CUTOFF);
            
            /* update expectation */
            NmrInv0    += NmrInvG[i] * StPriceG[i];
        }
        
        NmrInvAdj = 1/ (tree_data->TermZero[CvDiff][t] * NmrInv0);
        
        /* Readjust numeraire */
        for (i = 0; i < MappSize; i++)
        {
            NmrInvG[i] *= NmrInvAdj;
        }
        
        
        /* Now populate the numeraire slice with the interpplated values */
        if (InvMap(NmrInv,
                   NmrInvG,
                   MappSize,
                   tree_data->NmrMap_dir,
                   EDevIdx,
                   t,
                   tree_data) == FAILURE)
        {
            goto RETURN;
        }


        
        if (Nmr_Print3 (Yield,
                        YCumP,
                        LMin,
                        LMax,
                        NmrInvAdj,
                        NmrIdx,
                        t,
                        mktvol_data,
                        tree_data) == FAILURE) goto RETURN;
        
    }

    status = SUCCESS;

RETURN:

    if (tree_data->NbFactor == 1)
    {
        Free_DR_Array  (Yield, DOUBLE, LMinInit-1, LMaxInit);
        Free_DR_Array  (YCumP, DOUBLE, LMinInit-1, LMaxInit);
        Free_DR_Array  (YProb, DOUBLE, LMinInit-1, LMaxInit);
        Free_Slice_Int (Level, tree_data);
    }

    else
    {
        
        Free_DR_Array  (Yield, DOUBLE, 0, MappSize);
        Free_DR_Array  (YCumP, DOUBLE, 0, MappSize);
        Free_DR_Array  (YProb, DOUBLE, 0, MappSize);
        Free_Slice_Int (Level, tree_data);

        Free_DR_Array  (StPriceG, DOUBLE , 0, Area);
        Free_DR_Array  (Zero1G,   DOUBLE , 0, Area);
        Free_DR_Array  (Ann0G,    DOUBLE , 0, Area);
        Free_DR_Array  (Zero0G,   DOUBLE , 0, Area);
        Free_DR_Array  (NmrInvG,  DOUBLE , 0, Area);
    }

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_Interp: Failed!");
    }

    return (status);

} /* Nmr_Interp */





/*****  Nmr_Stat  **********************************************************/
/*
*       Calculate Nmr stats 
*/
int     Nmr_Stat   (int         t,               /* (I) Current time point */
                    int         NmrIdx,          /* (I) Numeraire index    */
                    double      *Annuity,        /* (I) Annuity slice      */
                    int         *Level,          /* (I) Node level         */
                    double      *YieldL,         /* (I) Yield slice        */
                    MKTVOL_DATA *mvd,            /* (I) Volatility data    */
                    TREE_DATA   *tree_data)      /* (I) Tree data          */
{
    double  *NmrInvL;
    double  *StPriceL;
    int     *LevelL;
    double  *AnnuityL;

    /* 2-factor variables */
    double *AnnuityG = NULL, 
           *StPriceG = NULL, 
           *NmrInvG  = NULL;
    long   MappSize  = 0;
    long   Area      = 0;
         
    double  OUNmr;             /* E{NmrInv} in OU measure     */
    double  AAStaFwd;          /* Fwd in tree                 */
    double  AAStaVar;          /* Variance in tree            */
    double  AAtoOUFact;        /* Change of measure constant  */
    double  TreePr[NBSTRIKE];  /* Options in tree             */
    double  OUFwd;             /* Fwd in OU measure           */
    double  OUATM;             /* ATM in OU measure           */
    double  OUSta[5];          /* Statistic in OU measure     */
    double  TreeSta[5];        /* Statistic in tree           */
    double  ParYield0 = mvd->ParYield0[NmrIdx];

    int     EDevIdx;           /* Express DEV index           */
    int     Top1, Bottom1;
    int     i, s, level;
    int     offset;
    int     status = FAILURE;               /* Error status */


    if (!mvd->NmrStatFlag)  return SUCCESS;

    AAtoOUFact = mvd->TermZero0 / mvd->Annuity0[NmrIdx];

    /* Get state price slice index */
    EDevIdx = GetDLOffset (tree_data->NbEDevDates,
                           tree_data->EDevDate,
                           tree_data->TPDate[t],
                           CbkEXACT);
    if (EDevIdx < 0)
    {
        DR_Error ("State prices not available for current numeraire date %ld\n",
                  tree_data->TPDate[t]);
        return (FAILURE);
    }

    /* Slices size */
    if (tree_data->NbFactor == 1)
    {
        Area = tree_data->Width[0];
    }
    if (tree_data->NbFactor == 2)
    {
        Area = tree_data->Width[0] * tree_data->Width[1];
    }
    if (tree_data->NbFactor == 3)
    {
        Area = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];
    }
    
    /* Clean up */
    AAStaFwd = AAStaVar = 0.;
    for (s=0; s<NBSTRIKE; s++) TreePr[s] = 0;
    for (s=0; s<OU_ALL; s++) OUSta[s] = 0;
    for (s=0; s<4; s++) TreeSta[s] = 0;

    if (tree_data->NbFactor == 1)
    {
        
        /* Tree architecture */
        Top1     = tree_data->Top1[t];
        Bottom1  = tree_data->Bottom1[t];
        
        offset    = Node_Offset(1, 0, 0, t, tree_data);
        LevelL    = Level    + offset;
        AnnuityL  = Annuity  + offset;
        NmrInvL   = tree_data->NmrInv[NmrIdx]       + offset;
        StPriceL  = tree_data->EDevStPrice[EDevIdx] + offset;
        
        /* Calc stats */
        for (i = Bottom1; i <= Top1; i++) 
        {
            double AStP = AnnuityL[i] * StPriceL[i];
            int    l    = LevelL[i];
            
            /* AA measure */
            AAStaFwd += AStP * YieldL[i];
            
            for (s = 0; s < NBSTRIKE; s++) 
            {
                TreePr[s] += AStP * MAX(YieldL[l]-mvd->Strike[NmrIdx][s],0);
            }
            
            /* OU measure */
            OUSta[OU_FWD] += StPriceL[i] * YieldL[l];
            OUSta[OU_NMR] += StPriceL[i] * NmrInvL[i] ;
            
            TreeSta[1] += AStP;
            TreeSta[2] += AStP * ParYield0;
            TreeSta[3] += AStP * (YieldL[l]-ParYield0);
        }
        
        AAStaFwd *= AAtoOUFact;
        
        for (s=0; s<NBSTRIKE; s++)
        {
            mvd->TreePr[NmrIdx][s] = TreePr[s] * AAtoOUFact;
        }
        
        for (i = Bottom1; i <= Top1; i++) 
        {
            double AStP = AnnuityL[i] * StPriceL[i];
            int    l    = LevelL[i];
            
            AAStaVar += AStP * pow (YieldL[i]-AAStaFwd,2.);
            
            OUSta[OU_ATM] += StPriceL[i] * MAX(YieldL[l]-OUSta[OU_FWD],0);
            OUSta[OU_VAR] += StPriceL[i] * pow(YieldL[l]-OUSta[OU_FWD],2.);
        }
        AAStaVar *= AAtoOUFact;
        mvd->AASta[NmrIdx][AA_VAR] = AAStaVar;
        
        mvd->OUSta[NmrIdx][OU_FWD] = OUSta[OU_FWD];             
        mvd->OUSta[NmrIdx][OU_ATM] = OUSta[OU_ATM];             
        mvd->OUSta[NmrIdx][OU_NMR] = OUSta[OU_NMR];
        mvd->OUSta[NmrIdx][OU_VAR] = OUSta[OU_VAR];
        
        for (s=1; s<4; s++) mvd->TreeSta[NmrIdx][s] = TreeSta[s] * AAtoOUFact;
        mvd->TreeSta[NmrIdx][1] -= 1.;
        mvd->TreeSta[NmrIdx][2] -= ParYield0;

    }

    else if (tree_data->NbFactor == 2)
    {
               
        AnnuityG = (double *) DR_Array(DOUBLE, 0, Area);
        StPriceG = (double *) DR_Array(DOUBLE, 0, Area);
        NmrInvG  = (double *) DR_Array(DOUBLE, 0, Area); 


        /* Project the 2D slices on a 1D grid */
        if ( SortTree_Slices(AnnuityG,                        
                             NmrInvG,
                             NULL,
                             NULL,        
                             StPriceG, 
                             &MappSize,
                             2,                    /* 2 Slices to project */
                             Annuity,
                             tree_data->NmrInv[NmrIdx],
                             NULL,
                             NULL,
                             tree_data->NmrMap_dir,  /* Numeraire mapping direction */
                             EDevIdx,
                             t,
                             tree_data) == FAILURE)
        {
            goto RETURN;
        }

        /* Calc stats */
        for (i = 0; i < MappSize; i++) 
        {
            double AStP = AnnuityG[i] * StPriceG[i];            
            
            /* AA measure */
            AAStaFwd += AStP * YieldL[i];
            
            for (s = 0; s < NBSTRIKE; s++) 
            {
                TreePr[s] += AStP * MAX(YieldL[i]-mvd->Strike[NmrIdx][s],0);
            }
            
            /* OU measure */
            OUSta[OU_FWD] += StPriceG[i] * YieldL[i];
            OUSta[OU_NMR] += StPriceG[i] * NmrInvG[i] ;
            
            TreeSta[1] += AStP;
            TreeSta[2] += AStP * ParYield0;
            TreeSta[3] += AStP * (YieldL[i]-ParYield0);
        }
        
        AAStaFwd *= AAtoOUFact;
        
        for (s=0; s<NBSTRIKE; s++)
        {
            mvd->TreePr[NmrIdx][s] = TreePr[s] * AAtoOUFact;
        }
        
        for (i = 0; i < MappSize; i++) 
        {
            double AStP = AnnuityG[i] * StPriceG[i];
                        
            AAStaVar += AStP * pow (YieldL[i]-AAStaFwd,2.);
            
            OUSta[OU_ATM] += StPriceG[i] * MAX(YieldL[i]-OUSta[OU_FWD],0);
            OUSta[OU_VAR] += StPriceG[i] * pow(YieldL[i]-OUSta[OU_FWD],2.);
        }
        AAStaVar *= AAtoOUFact;
        mvd->AASta[NmrIdx][AA_VAR] = AAStaVar;
        
        mvd->OUSta[NmrIdx][OU_FWD] = OUSta[OU_FWD];             
        mvd->OUSta[NmrIdx][OU_ATM] = OUSta[OU_ATM];             
        mvd->OUSta[NmrIdx][OU_NMR] = OUSta[OU_NMR];
        mvd->OUSta[NmrIdx][OU_VAR] = OUSta[OU_VAR];
        
        for (s=1; s<4; s++) mvd->TreeSta[NmrIdx][s] = TreeSta[s] * AAtoOUFact;
        mvd->TreeSta[NmrIdx][1] -= 1.;
        mvd->TreeSta[NmrIdx][2] -= ParYield0;
    }
    
    
    status = SUCCESS;

RETURN:

    if (tree_data->NbFactor == 2)
    {
        Free_DR_Array(AnnuityG, DOUBLE, 0, Area);
        Free_DR_Array(StPriceG, DOUBLE, 0, Area);
        Free_DR_Array(NmrInvG,  DOUBLE, 0, Area);
    }
    
    if (status != SUCCESS)
    {        
        DR_Error("Nmr_Stat: Failed!");
    }

    return (status);

} /* Nmr_Stat */

    
 
/*****  Nmr_Print0  ******************************************************/
/*
*       Print numeraire information 
*/
int     Nmr_Print0 (MKTVOL_DATA *mvd)  /* (I) Volatility data    */
{
    char    FileName[MAXBUFF];
    FILE    *stream = NULL;
    int     status = FAILURE;               /* Error status */
    int     n;


    if (!mvd->Trace) return SUCCESS;

    strcpy (FileName, "NMR.prn");
    if ((stream = fopen (FileName, "w")) == NULL)
    {
        DR_Error ("Could not open file %s.\n", FileName);
        goto RETURN;
    }

    /* Print vanilla date info */
    fprintf (stream, "###### NMR0 DATE INFO ######\n"); 
    fprintf (stream, "# TermDate \n%8ld\n",mvd->TermDate);
    fprintf (stream, "# [idx] NmrDate    SwapSt    SwapMat\n");

    for (n=0; n<mvd->NbNmr; n++)
    {
        fprintf (stream, "  [%3d] %8ld   %8ld  %8ld   \n",
                 n,
                 mvd->NmrDate[n],
                 mvd->NmrSwapSt[n],
                 mvd->NmrSwapMat[n]
                );
    }

    status = SUCCESS;

    RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_Print0: Failed!");
    }

    return (status);

} /* Nmr_Print0 */


/*****  Nmr_Print1  ******************************************************/
/*
*       Print numeraire information 
*/
int     Nmr_Print1(MKTVOL_DATA *mvd,        /* (I) Volatility data */
                   MQDATA      *mq)         /* (I) Tree data       */
{
    char    FileName[MAXBUFF];
    FILE    *stream = NULL;
    int     status = FAILURE;               /* Error status */
    int     n;


    if (!mvd->Trace) return SUCCESS;

    strcpy (FileName, "NMR.prn");
    if ((stream = fopen (FileName, "a")) == NULL)
    {
        DR_Error ("Could not open file %s.\n", FileName);
        goto RETURN;
    }

    /* Print vanila input info */
    fprintf (stream, "###### NMR0 VANILA MVD INFO ######\n"); 
    fprintf (stream, "# TermDate \n%8ld\n",mvd->TermDate);
    fprintf (stream, "# [idx] NmrDate    SwapSt    SwapMat    FwdRate    "
             "vol_____  q___  vvol  bbV_  bbR_  dL__  tL___  dR__  tR___\n");

    for (n=0; n<mvd->NbNmr; n++)
    {
        fprintf (stream, "  [%3d] %8ld   %8ld  %8ld   %8.4f   %8.4f  "
                 "%4.2f  %4.2f  %4.2f  %4.2f  %4.2f  %5.2f  %4.2f  %5.2f\n",
                 n,
                 mvd->NmrDate[n],
                 mvd->NmrSwapSt[n],
                 mvd->NmrSwapMat[n],
                 mvd->ParYield0[n]*100,
                 mvd->NmrVol[n][0]*100,
                 mvd->NmrVol[n][1],
                 mvd->NmrVol[n][2],
                 mvd->NmrVol[n][3],
                 mvd->NmrVol[n][4],
                 mvd->NmrVol[n][5],
                 mvd->NmrVol[n][6],
                 mvd->NmrVol[n][7],
                 mvd->NmrVol[n][8]
                );
    }
    /* Print vanila output info */
    fprintf (stream, "###### NMR0 VANILLA MQ CALIB INFO ######\n"); 
    fprintf (stream, "# [idx] expiry__  fwd_____  sigATM__  muMQ____  sigMQ___"
             "  qL[1]___  qL[0]___  qR[0]___  qR[1]___  qR[2]___"
             "  fwdC____  volC____  \n");

    for (n=0; n<mvd->NbNmr; n++)
    {
        fprintf (stream, "  [%3d] %8.2f  %8.4f  %8.4f  %8.4f  %8.4f  "
                 "%8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  \n",
                 n,
                 mq[n].expiry,
                 100*mq[n].fwdRate,
                 100*mq[n].sigATM,
                 mq[n].muMQ,
                 mq[n].sigMQ,
                 mq[n].qL[1],
                 mq[n].qL[0],
                 mq[n].qR[0],
                 mq[n].qR[1],
                 mq[n].qR[2],
                 mq[n].fwdC,
                 mq[n].volC                 
                );
    }

    status = SUCCESS;

    RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_Print1: Failed!");
    }

    return (status);

} /* Nmr_Print1 */


/*****  Nmr_Print2  ******************************************************/
/*
*       Print numeraire information 
*/
int     Nmr_Print2 (MKTVOL_DATA *mvd,  /* (I) Volatility data    */
                    TREE_DATA   *trd)  /* (I) Tree data          */
{
    char    FileName[MAXBUFF];
    FILE    *stream = NULL;
    int     status = FAILURE;               /* Error status */
    int     n, s;
    double  ErrFwdC=0, ErrVolC=0, ErrNmrC=0;


    if (!mvd->Trace) return SUCCESS;

    strcpy (FileName, "NMR.prn");
    if ((stream = fopen (FileName, "a")) == NULL)
    {
        DR_Error ("Could not open file %s.\n", FileName);
        goto RETURN;
    }

    /* Print nmr calib info */
    fprintf (stream, "###### NMR0 TREE NMR CALIB INFO ######\n"); 
    fprintf (stream, "# [idx] NmrDate    Fwd_____  FwdC______  "
                     "ATMCall_   VolC_____   AAStd__  "
                     "1_NmrInv____  NmrC______   OUStd__  "
                     "TreeSta0__  TreeSta1__  TreeSta2__  TreeSta3__  "
                     "TreeSta4_    TreeSta5__  \n");

    for (n=1; n<mvd->NbNmr-1; n++)
    {
        double   FwdC = mvd->AASta[n][AA_FWD]/mvd->ParYield0[n];
        double   VolC = mvd->AASta[n][AA_VOL];
        double   NmrC = mvd->TermZero0*mvd->OUSta[n][OU_NMR]/mvd->StZero0[n];
        double   TExp = (double)Daysact(mvd->NmrDate[0],mvd->NmrDate[n])/365.;

        ErrFwdC = MAX(ErrFwdC, fabs(FwdC-1.));
        ErrVolC = MAX(ErrVolC, fabs(VolC-1.));
        ErrNmrC = MAX(ErrNmrC, fabs(NmrC-1.));

        fprintf (stream, "  [%3d] %8ld   %8.4f  %10.2e  %8.4f  %10.2e  "
                 "%8.4f  %12.6f  %10.2e  %8.4f  "
                 "%10.2e  %10.2e  %10.2e  %10.2e  %9.4f  %10.2e\n",
                 n,
                 mvd->NmrDate[n],
                 mvd->AASta[n][AA_FWD]*100,
                 FwdC-1.,
                 mvd->AASta[n][AA_VOL]*100,
                 VolC-1.,
                 sqrt (mvd->AASta[n][AA_VAR]/TExp)*100,
                 mvd->TermZero0*mvd->OUSta[n][OU_NMR]*10000,
                 NmrC-1.,
                 sqrt (mvd->OUSta[n][OU_VAR]/TExp)*100,
                 mvd->TreeSta[n][0],
                 mvd->TreeSta[n][1],
                 mvd->TreeSta[n][2],
                 mvd->TreeSta[n][3],
                 mvd->TreeSta[n][4]*100.,
                 mvd->TreeSta[n][5]
                );
    }
    fprintf (stream, "## ErrFwdC___  ErrVolC___  ErrNmrC___\n");
    fprintf (stream, "  %10.2e  %10.2e  %10.2e\n",
             ErrFwdC, ErrVolC, ErrNmrC);

    /* Print mapped yield info */
    fprintf (stream, "##### MAPPED YIELD INFO IN OU MEASURE#####\n");
    for (n=1; n<mvd->NbNmr-1; n++)
    {
        fprintf (stream,"%8d |  Index  ",mvd->NmrDate[n]);
        for (s=trd->YLMin[n];s<=trd->YLMax[n];s++) fprintf (stream,"%9d  ",s);
        fprintf (stream,"\nYMin     |  YCumP  ");
        for (s=trd->YLMin[n];s<=trd->YLMax[n];s++) 
            fprintf (stream,"% 8.6f  ",trd->YCumP[n][s]);
        fprintf (stream,"\n% 8.6f|  Yield  ", mvd->TreeSta[n][4]);
        for (s=trd->YLMin[n];s<=trd->YLMax[n];s++) 
            fprintf (stream,"% 8.6f  ",trd->Yield[n][s]);
        fprintf (stream,"\n            ZAInv  ");
        for (s=trd->YLMin[n];s<=trd->YLMax[n];s++) 
            fprintf (stream,"% 8.6f  ",trd->ZAInv[n][s]);
        fprintf (stream,"\n");
    }

    /* Print smile stat info */
    fprintf (stream, "###### NMR0 SMILE VOL INFO ######\n"); 
    fprintf (stream, "# [idx] NmrDate   Strike/VanlVol/TreeVol   "
             "                             ErrSmile    ErrAtm\n");

    for (n=1; n<mvd->NbNmr-1; n++)
    {
        double   VanlVol[NBSTRIKE];
        double   TreeVol[NBSTRIKE];
        double   Expiry;
        double   ErrVol=0;
        double   ErrAtm=0;

        Expiry = Daysact(mvd->NmrDate[0], mvd->NmrDate[n]) / 365.;
        
        for (s=0; s<NBSTRIKE; s++)
        {
            VanlVol[s] = ImpVol_BSQ (mvd->ParYield0[n],
                                      mvd->Strike[n][s],
                                      Expiry,
                                      mvd->VanlPr[n][s],
                                      'C',
                                      1, /*  1, 0, */ 
                                      mvd->NmrVol[n][0]);
            TreeVol[s] = ImpVol_BSQ (mvd->ParYield0[n],
                                      mvd->Strike[n][s],
                                      Expiry,
                                      mvd->TreePr[n][s],
                                      'C',
                                      1, /* 1, 0, */ 
                                      mvd->NmrVol[n][0]);

            ErrVol = MAX(ErrVol,fabs(TreeVol[s]-VanlVol[s]));
        }
        /* Use the fact that s=2 corresponds to ATM strike */
        ErrAtm = MAX (ErrAtm, fabs (TreeVol[2]-VanlVol[2]));

        fprintf (stream, "  [%3d] %8ld  ", n, mvd->NmrDate[n]);
        for (s=0;s<NBSTRIKE;s++) 
            fprintf(stream,"%8.4f  ",mvd->Strike[n][s]*100);
        fprintf(stream, "%8.2f    %8.2f",ErrVol*100, ErrAtm*100);
        fprintf(stream, "\n");

        fprintf (stream, "                  ");
        for (s=0;s<NBSTRIKE;s++)
            fprintf(stream,"%8.2f  ",VanlVol[s]*100);
        fprintf(stream, "\n");

        fprintf (stream, "                  ");
        for (s=0;s<NBSTRIKE;s++)
            fprintf(stream,"%8.2f  ",TreeVol[s]*100);
        fprintf(stream, "\n");
    }

    fprintf (stream, "###### NMR0 SMILE PRICE INFO ######\n"); 
    fprintf (stream, "# [idx] NmrDate   Strike/VanlPrice/TreePrice \n");

    for (n=1; n<mvd->NbNmr-1; n++)
    {
        fprintf (stream, "  [%3d] %8ld  ", n, mvd->NmrDate[n]);
        for (s=0;s<NBSTRIKE;s++) 
            fprintf(stream,"%8.4f  ",mvd->Strike[n][s]*100);
        fprintf(stream,"\n");

        fprintf (stream, "                  ");
        for (s=0;s<NBSTRIKE;s++)
            fprintf(stream,"%8.2f  ",mvd->VanlPr[n][s]*mvd->Annuity0[n]*10000);
        fprintf(stream,"\n");

        fprintf (stream, "                  ");
        for (s=0;s<NBSTRIKE;s++)
            fprintf(stream,"%8.2f  ",mvd->TreePr[n][s]*mvd->Annuity0[n]*10000);
        fprintf(stream,"\n");
    }


    /* Print convexity info */
    fprintf (stream, "###### NMR0 CONVEXITY INFO ######\n"); 
    fprintf (stream, "# [idx] NmrDate   Fwd    OUFwd   Vol   OUVol\n");

    for (n=1; n<mvd->NbNmr-1; n++)
    {
        double Expiry, OUVol;

        Expiry = Daysact(mvd->NmrDate[0], mvd->NmrDate[n]) / 365.;
        OUVol = ImpVol_BS2Q (mvd->OUSta[n][OU_FWD],
                             mvd->OUSta[n][OU_FWD],
                             Expiry,
                             mvd->OUSta[n][OU_ATM],
                             'C',
                             1, 1, 0, 
                             mvd->NmrVol[n][0]);

        fprintf (stream, "  [%3d] %8ld  %8.4f  %8.4f  %8.4f  %8.4f\n",
                 n,
                 mvd->NmrDate[n],
                 mvd->ParYield0[n]*100,
                 mvd->OUSta[n][OU_FWD]*100,
                 mvd->NmrVol[n][0]*100,
                 OUVol*100
                 );
    }

    /* Print header for statistics on interpolated yield info;
     * actual printout part of Lattice */
    fprintf (stream, "\n\n####### INTERP NMR INFO #######\n");
    fprintf (stream, "#  Node  NmrDate1  CurrDate  NmrDate2    "
             "YExp1     Y_Exp     Y_Exp2    Y_Std1    Y_Std     "
             "Y_Std2    NmrAdj\n");


    status = SUCCESS;

    RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_Print2: Failed!");
    }

    return (status);

} /* Nmr_Print2 */


/*****  Nmr_Print3  ******************************************************/
/*
*       Print mapped yield statistics 
*/
int     Nmr_Print3 (double          *Yield,    /* (I) Mapped yield       */
                    double          *YCumP,    /* (I) Cumulative prob    */
                    int             LMin,      /* (I) Minimum level      */
                    int             LMax,      /* (I) Maximum level      */
                    double          NmrInvAdj, /* (I) Nmr adjustment     */
                    int             NmrIdx,    /* (I) Numeraire index    */
                    int             t,         /* (I) Current TP         */
                    MKTVOL_DATA     *mvd,      /* (I) Volatility data    */
                    TREE_DATA       *trd)      /* (I) Tree data          */
{
    /* Local pointers */
    double  *YCumP1L;
    double  *Yield1L;
    double  *YCumP0L;
    double  *Yield0L;
    double  *StPriceL;
    int     *LevelL;

    /* Expectations (E) and standard deviations (S) */
    double EY=0.,EY1=0.,EY0=0.,SY=0.,SY1=0.,SY0=0.;
    int    LMin1, LMax1;
    int    LMin0, LMax0;
    int    l;
    long   CurrentDate = trd->TPDate[t];

    FILE    *stream = NULL;
    char    FileName[MAXBUFF];
    int     status = FAILURE;


    if (!mvd->Trace) return SUCCESS;

    strcpy (FileName, "NMR.prn");
    if ((stream = fopen (FileName, "a")) == NULL)
    {
        DR_Error ("Could not open file %s.\n", FileName);
        goto RETURN;
    }

    /* 1D arrays */
    LMin1    = trd->YLMin[NmrIdx];
    LMax1    = trd->YLMax[NmrIdx];
    LMin0    = trd->YLMin[NmrIdx-1];
    LMax0    = trd->YLMax[NmrIdx-1];
    YCumP1L  = trd->YCumP[NmrIdx];
    Yield1L  = trd->Yield[NmrIdx];
    YCumP0L  = trd->YCumP[NmrIdx-1];
    Yield0L  = trd->Yield[NmrIdx-1];

    for (l= LMin1; l<=LMax1; l++)
    {
        EY1 += Yield1L[l] * (YCumP1L[l] - YCumP1L[l-1]);
    }
    for (l= LMin1; l<=LMax1; l++)
    {
        SY1 += pow (Yield1L[l]-EY1,2) * (YCumP1L[l] - YCumP1L[l-1]);
    }
    SY1 = sqrt (SY1 * 365./Daysact(trd->NmrDate[0],trd->NmrDate[NmrIdx]));

    for (l= LMin; l<=LMax; l++)
    {
        EY += Yield[l] * (YCumP[l] - YCumP[l-1]);
    }
    for (l= LMin; l<=LMax; l++)
    {
        SY += pow (Yield[l]-EY,2) * (YCumP[l] - YCumP[l-1]);
    }
    SY = sqrt (SY * 365./Daysact(trd->NmrDate[0], CurrentDate));
                                     
    for (l= LMin0; l<=LMax0; l++)
    {
        EY0 += Yield0L[l] * (YCumP0L[l] - YCumP0L[l-1]);
    }
    for (l= LMin0; l<=LMax0; l++)
    {
        SY0 += pow (Yield0L[l]-EY0,2) * (YCumP0L[l] - YCumP0L[l-1]);
    }
    SY0 = sqrt (SY0 * 365./Daysact(trd->NmrDate[0],
                                     trd->NmrDate[NmrIdx-1]));
    
    fprintf (stream,"  [%4d] %ld  %ld  %ld  %8.4f  %8.4f  %8.4f  "
                    "%8.4f  %8.4f  %8.4f    %10.8f\n",
             t, 
             trd->NmrDate[NmrIdx], 
             CurrentDate,
             trd->NmrDate[NmrIdx-1],
             EY1*100,
             EY*100,
             EY0*100,
             SY1*100,
             SY*100,
             SY0*100,
             NmrInvAdj);
    fclose (stream);

    status = SUCCESS;

RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_Print3: Failed!");
    }

    return (status);

} /* Nmr_Print3 */


 

/************************************************************************************
                          NEW ROUTINES ADDED BY BM FOR THE 2D IMPLEMENTATION
************************************************************************************/


 /*****  SortList  **********************************************/
/*
 *      Given a list and its size
 *      sort its contents in 1st: date ascending order, 
 */

int     SortList(long    ListSize,          /* (I) : list size          */
                 double  *Index,            /* (I) : list               */
                 long    *SuppValue,        /* (I) : list index         */
                 long    *SortVal,          /* (O) : sorted index list  */
                 long    *NbVal)            /* (O) : Grid size          */
{
    double  tmpDate;
    long    tmpValue, counter;  /* temp storage */
    int     i, j;
    long    *Val = NULL;


    if (Index == NULL) return (FAILURE);    
   Val = SuppValue;
   

   /* First sort in increasing order */
   for (j = 1 ; j < ListSize; j++)
   {
       tmpDate = Index[j];
       tmpValue = Val[j];
       
       i = j-1;
       while ((i >= 0) && (Index[i] > tmpDate))
       {
           Index[i+1] = Index[i];
           Val[i+1]      = Val[i];
           i--;
       }
       
       Index[i+1] = tmpDate;
       Val[i+1] = tmpValue;
   
   }  /* for j */
   
   
   /* Now sort in stricly increasing order */   
   counter = 0;
   SortVal[0] = Val[0];
   for (j=1; j < ListSize; j++)
   {
       if ( (Index[j] - Index[j-1]) > TINY)
       {
           counter ++;
           SortVal[counter] = Val[j];
       }
   }   
   *NbVal = counter + 1;
   
   return (SUCCESS);
} /* SortList */


/*****  SortTree_Slices  ********************************************************/
/*
*       Given the weights vector, projects weight1 * x + weight2 * y in a 1D grid
        after sorting their values in an increasing order.
        Gibven that all Zeros are monotonic functions of the driver, all the 
        input slices are also sorted in the same fashion                           
*/
int     SortTree_Slices(double         *SliceGrid1,        /* (O) : 1D grid obtained after projecting the slice 1            */
                        double         *SliceGrid2,        /* (O) : 1D grid obtained after projecting the slice 2            */
                        double         *SliceGrid3,        /* (O) : 1D grid obtained after projecting the slice 3            */
                        double         *SliceGrid4,        /* (O) : 1D grid obtained after projecting the slice 4            */                
                        double         *StatePrGrid,       /* (O) : 1D grid obtained after projecting the state prices grid  */    
                        long           *MappSize,          /* (O) : Number of grid points after sorting                      */
                        int            NbSlices,           /* (I) : Number of slices to project                              */ 
                        double         *Slice1,            /* (I) : 2D slice 1                                               */   
                        double         *Slice2,            /* (I) : 2D slice 2                                               */                 
                        double         *Slice3,            /* (I) : 2D slice 3                                               */                 
                        double         *Slice4,            /* (I) : 2D slice 4                                               */                 
                        double         weights[2],         /* (I) : Gaussian drivers weights for the mapping                 */
                        long           EDevIdx,            /* (I) : Dev Index                                                */
                        long           t,                  /* (I) : current time t                                           */
                        TREE_DATA      *tree_data)         /* (I) : tree data structure                                      */
{
    /* Indices */
    long       i , k ,Node, counter, SliceSize;
    int        offset;
    long       Area;

    /* Pointers to the input 2D slices */
    double     *YValuesL;          /* Local variable for YValues  slice */
    double     *XValuesL;          /* Local variable for YValues  slice */  
    double     *StatePrL;          /* Local variable for StatePrL slice */
    double     *Slice1L;           /* Local variable for Slice 1        */
    double     *Slice2L;           /* Local variable for Slice 2        */
    double     *Slice3L;           /* Local variable for Slice 3        */
    double     *Slice4L;           /* Local variable for Slice 4        */


    /* Undimensional grids after slices projection */  
    double     *GaussianDriverG;          /* Unidimensional Grid for the gaussian driver */
    double     *StatePrG;                 /* Unidimensional Grid for StatePrL slice      */
    double     *Slice1G = NULL;            /* Unidimensional Grid for slice 1             */
    double     *Slice2G = NULL;           /* Unidimensional Grid for other 2             */
    double     *Slice3G = NULL;           /* Unidimensional Grid for other 3             */
    double     *Slice4G = NULL;           /* Unidimensional Grid for other 4             */
    
    /* Mapping correspondance between tree and projected grids */
    long       *Mapping, 
               *SortIndex;


    /* Limits */
    int     Top1, Bottom1;              /* Tree limits (1rst dim)           */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)            */
    int     OutTop1, OutBottom1;        /* Outer tree limits                */
    int     *OutTop2, *OutBottom2;
    int     status = FAILURE;
    
    
    /* Slices size */
    Area = tree_data->Width[0] * tree_data->Width[1];

    /* Make sure the number of factors is two */
    if (tree_data->NbFactor != 2)
    {
        DR_Error("SortTree_Slices : the number of factor must be 2 "
                 "when this function is called!");
        return(status);
    }   
    

    /* Check the number of slices is less than 4 */
    if (NbSlices > 4)
    {
        DR_Error("SortTree_Slices : Only supports 4 slices projection !");
        return(FAILURE);
    }

    /* Allocate memory for the unidimensional grids */  
    GaussianDriverG = DR_Array (DOUBLE, 0, Area);
    StatePrG        = DR_Array (DOUBLE, 0, Area);
    Mapping         = DR_Array (LONG,   0, Area);
    SortIndex       = DR_Array (LONG,   0, Area);

    /* Allocate memory */
    Slice1G              = DR_Array (DOUBLE, 0, Area);   
    if (NbSlices >= 2)
        Slice2G          = DR_Array (DOUBLE, 0, Area);
    if (NbSlices >= 3)
        Slice3G          = DR_Array (DOUBLE, 0, Area);
    if (NbSlices >= 4)
        Slice4G          = DR_Array (DOUBLE, 0, Area);

    
    /* Tree limits for convenience */
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];



    /* Project all the 2D slices from the tree on a unidimensional grid */
    counter = 0;


    /* 1 slice to project */
    if (NbSlices == 1)
    {
        for (i = Bottom1; i <= Top1; i++)
        {
            offset     =  Node_Offset(2, i, 0, t, tree_data);
            StatePrL   =  tree_data->EDevStPrice[EDevIdx] + offset;
            YValuesL   =  tree_data->YValues[EDevIdx] + offset;
            XValuesL   =  tree_data->XValues[EDevIdx] + offset;
            Slice1L    =  Slice1  + offset;
            
            for (k = Bottom2[i]; k <= Top2[i]; k++)
            {
                StatePrG[counter]         = StatePrL[k];
                GaussianDriverG[counter]  = weights[0] * XValuesL[k] + weights[1] * YValuesL[k];
                Mapping[counter]          = counter;        
                Slice1G[counter]          = Slice1L[k];     
                counter ++;
            }
        }
    }


    /* 2 slices to project */
    if (NbSlices == 2)
    {
        for (i = Bottom1; i <= Top1; i++)
        {
            offset     =  Node_Offset(2, i, 0, t, tree_data);
            StatePrL   =  tree_data->EDevStPrice[EDevIdx] + offset;
            YValuesL   =  tree_data->YValues[EDevIdx] + offset;
            XValuesL   =  tree_data->XValues[EDevIdx] + offset;
            Slice1L    =  Slice1  + offset;
            Slice2L    =  Slice2  + offset;
            
            
            for (k = Bottom2[i]; k <= Top2[i]; k++)
            {
                StatePrG[counter]         = StatePrL[k];
                GaussianDriverG[counter]  = weights[0] * XValuesL[k] + weights[1] * YValuesL[k];
                Mapping[counter]          = counter;
                Slice1G[counter]          = Slice1L[k];
                Slice2G[counter]          = Slice2L[k];             
                counter ++;
            }
        }
    }

    /* 3 slices to project */
    if (NbSlices == 3)
    {
        for (i = Bottom1; i <= Top1; i++)
        {
            offset     =  Node_Offset(2, i, 0, t, tree_data);
            StatePrL   =  tree_data->EDevStPrice[EDevIdx] + offset;
            YValuesL   =  tree_data->YValues[EDevIdx] + offset;
            XValuesL   =  tree_data->XValues[EDevIdx] + offset;
            Slice1L    =  Slice1  + offset;
            Slice2L    =  Slice2  + offset;
            Slice3L    =  Slice3  + offset;     
            
            for (k = Bottom2[i]; k <= Top2[i]; k++)
            {
                StatePrG[counter]         = StatePrL[k];
                GaussianDriverG[counter]  = weights[0] * XValuesL[k] + weights[1] * YValuesL[k];
                Mapping[counter]          = counter;
                Slice1G[counter]          = Slice1L[k];
                Slice2G[counter]          = Slice2L[k];     
                Slice3G[counter]          = Slice3L[k];             
                counter ++;
            }
        }
    }


    /* 4 slices to project */
    if (NbSlices == 4)
    {
        for (i = Bottom1; i <= Top1; i++)
        {
            offset     =  Node_Offset(2, i, 0, t, tree_data);
            StatePrL   =  tree_data->EDevStPrice[EDevIdx] + offset;
            YValuesL   =  tree_data->YValues[EDevIdx] + offset;
            XValuesL   =  tree_data->XValues[EDevIdx] + offset;
            Slice1L    =  Slice1  + offset;
            Slice2L    =  Slice2  + offset;
            Slice3L    =  Slice3  + offset;     
            Slice4L    =  Slice4  + offset;     
            
            
            for (k = Bottom2[i]; k <= Top2[i]; k++)
            {
                StatePrG[counter]         = StatePrL[k];
                GaussianDriverG[counter]  = weights[0] * XValuesL[k] + weights[1] * YValuesL[k];
                Mapping[counter]          = counter;                
                Slice1G[counter]          = Slice1L[k];
                Slice2G[counter]          = Slice2L[k];     
                Slice3G[counter]          = Slice3L[k];                             
                Slice4G[counter]          = Slice4L[k];                             
                counter ++;
            }
        }
    }   

    /* Sort now the values of Gaussian driver  in ascending order and 
       get the mapping between the new Gaussian driver array and the tree nodes */
    if (SortList(counter,  
                 GaussianDriverG,
                 Mapping,
                 SortIndex,
                 MappSize) == FAILURE)
    {
        DR_Error ("SortTree_Slices: Could not sort the Gaussian driver array !");
        goto RETURN;
    }

    /* Store slices size for reuse */
    SliceSize   = counter;

    
    /* Now populate the projected grids */   
    
    /* 1 slice to project */
    if (NbSlices == 1)
    {
        for (i=0; i <  *MappSize; i++)
        {
            Node            = SortIndex[i];
            SliceGrid1[i]   = Slice1G[Node];
        }
    }

    /* 2 slices to project */
    if (NbSlices == 2)
    {
        for (i=0; i <  *MappSize; i++)
        {
            Node            = SortIndex[i];
            SliceGrid1[i]   = Slice1G[Node];
            SliceGrid2[i]   = Slice2G[Node];

        }
    }

    /* 3 slices to project */
    if (NbSlices == 3)
    {
        for (i=0; i <  *MappSize; i++)
        {
            Node            = SortIndex[i];
            SliceGrid1[i]   = Slice1G[Node];
            SliceGrid2[i]   = Slice2G[Node];
            SliceGrid3[i]   = Slice3G[Node];
        }
    }

    /* 4 slices to project */
    if (NbSlices == 4)
    {
        for (i=0; i <  *MappSize; i++)
        {
            Node            = SortIndex[i];
            SliceGrid1[i]   = Slice1G[Node];
            SliceGrid2[i]   = Slice2G[Node];
            SliceGrid3[i]   = Slice3G[Node];
            SliceGrid4[i]   = Slice4G[Node];
        }
    }



    /* FINALLY PROJECT THE STATE PRICES */

    /* Start counter */
    counter = 0;

    /* Initialization */
    Node = Mapping[0];  
    StatePrGrid[0] = StatePrG[Node];  
    
    /* The State price grid requires more work : addding the probs when Gaussian driver is flat */
    for (i=1;  i < SliceSize; i++)
    {
        if ((GaussianDriverG[i] -  GaussianDriverG[i-1]) > TINY)
        {
            counter ++;
            Node                 = Mapping[i];
            StatePrGrid[counter] = StatePrG[Node];
        }
        else
        {
            Node                  = Mapping[i];
            StatePrGrid[counter] += StatePrG[Node];
        }
    }


    status = SUCCESS;


RETURN:

    /* Free memory */
    Free_DR_Array (GaussianDriverG, DOUBLE,  0, Area);
    Free_DR_Array (StatePrG,        DOUBLE,  0, Area);
    Free_DR_Array (Mapping,         LONG,    0, Area);
    Free_DR_Array (SortIndex,       LONG,    0, Area);
    Free_DR_Array (Slice1G,         DOUBLE,  0, Area);
    
    if (NbSlices >= 2)
        Free_DR_Array (Slice2G,          DOUBLE,  0, Area);
    
    if (NbSlices >= 3)
        Free_DR_Array (Slice3G,          DOUBLE,  0, Area);
    
    if (NbSlices >= 4)
        Free_DR_Array (Slice4G,          DOUBLE,  0, Area);  
    
    return (status);
}






/*****  InvMap  ********************************************************/
/*
*       Buid a tree slice given a undimensional grid
*/
int     InvMap(double         *Slice,             /* (O) : 2D - Tree Slice                                          */
               double         *SliceGrid,         /* (O) : Unidimensional grid                                      */    
               long           ArraySize1,         /* (I) : size of Slice Grid                                       */
               double         weights[2],         /* (I) : Gaussian drivers weights for the mapping                 */
               long           EDevIdx,            /* (I) : Dev Index                                                */
               long           t,                  /* (I) : current time t                                           */
               TREE_DATA      *tree_data)         /* (I) : tree data structure                                      */
{
    /* Indices */
    long       i , k ,Node, counter;
    int        offset;

    /* Slices size */
    long       Area = 0;
    
    /* Slices */
    double     *YValuesL;          /* Local variable for YValues  slice */
    double     *XValuesL;          /* Local variable for YValues  slice */  
    double     *SliceL;            /* Local variable for Slice          */

    /* 1D Grids */  
    double     *GaussianDriverG;          /* Unidimensional Grid for the gaussian driver */
    double     *SliceG;                   /* Unidimensional Grid for Annuity slice       */
    long       *Mapping, *SortIndex;
    long       ArraySize2;
    long       SlicesSize;
    double     PrevVal;


    /* Limits */
    int     Top1, Bottom1;              /* Tree limits (1rst dim)           */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)            */
    int     OutTop1, OutBottom1;        /* Outer tree limits                */
    int     *OutTop2, *OutBottom2;
    int     status = FAILURE;

    
    /* slices size */
    Area = tree_data->Width[0] * tree_data->Width[1];
    
    /* Allocate memory for the unidimensional grids */  
    GaussianDriverG = DR_Array (DOUBLE,   0, Area);
    Mapping         = DR_Array (LONG,     0, Area);
    SortIndex       = DR_Array (LONG,     0, Area);
    SliceG          = DR_Array (DOUBLE,   0, Area);

    /* Tree limits for convenience */
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    
    
    /* Project all the 2D slices from the tree on a unidimensional grid */
    counter = 0;
    for (i = Bottom1; i <= Top1; i++)
    {
        offset     =  Node_Offset(2, i, 0, t, tree_data);
        YValuesL   =  tree_data->YValues[EDevIdx] + offset;    
        XValuesL   =  tree_data->XValues[EDevIdx] + offset;    


        for (k = Bottom2[i]; k <= Top2[i]; k++)
        {
            GaussianDriverG[counter]  = weights[0] * XValuesL[k] + weights[1] * YValuesL[k];
            Mapping[counter]          = counter;
            counter ++;
        }
    }

    
    /* Sort now the values of Gaussian driver  in acsending order and 
       get the mapping between the new Gaussian driver array and the tree nodes */
    if (SortList(counter,  
                 GaussianDriverG,
                 Mapping,
                 SortIndex,
                 &ArraySize2) == FAILURE)
    {
        DR_Error ("InvMap: Could not sort the Gaussian driver array !");
        goto RETURN;
    }


    /* Check the consistency between the input/output arrays 
       size                                                   */ 
    if (ArraySize1 != ArraySize2)
    {
        DR_Error("InvMap: Size not consistent !");
        goto RETURN;
    }


    /* Store the slices size for reuse */
    SlicesSize = counter;

    /* Now populate the sorted array in a 2D projected array */
    counter      = 0;
    Node         = Mapping[0];
    SliceG[Node] = SliceGrid[0];
    PrevVal      = SliceGrid[0];
    for (i=1; i < SlicesSize; i++)
    {
        if ((GaussianDriverG[i] -  GaussianDriverG[i-1]) <= TINY)
        {
            Node          = Mapping[i];
            SliceG[Node]  = PrevVal;
        }
        else
        {
            counter ++;
            Node          = Mapping[i];
            SliceG[Node]  = SliceGrid[counter];
        }

        PrevVal = SliceG[Node];
    }


    /* Check again consistency with input/output arrays size */
    if (counter + 1 != ArraySize1)
    {
        DR_Error("InvMap : Size not consistent !");
        goto RETURN;
    }


    /* Now populate the slice on the tree*/
    counter = 0;
    for (i = Bottom1; i <= Top1; i++)
    {
        offset     =  Node_Offset(2, i, 0, t, tree_data);
        SliceL     =  Slice  + offset;

        for (k = Bottom2[i]; k <= Top2[i]; k++)
        {
            SliceL[k] = SliceG[counter];  
            counter ++;
        }
    }        


    status = SUCCESS;


RETURN:
    
    /* Free memory */
    Free_DR_Array (GaussianDriverG, DOUBLE,  0, Area);
    Free_DR_Array (SliceG,          DOUBLE,  0, Area);
    Free_DR_Array (Mapping,         LONG,    0, Area);
    Free_DR_Array (SortIndex,       LONG,    0, Area);

    return (status);
}


/*****  Yield_t   ***************************************************************/
/*
*       Generate yield slice at time t out of vanilla smile information
*/
int     Yield_t (double      **Yield,     /* (O) Yield pointer              */
                 double      weights[2],  /* (I) Mapping weights            */
                 double      *Annuity,    /* (I) Annuity for calib swap     */
                 double      ParYield0,   /* (I) Forward yield              */
                 double      Annuity0,    /* (I) Forward annuity            */
                 double      TermZero0,   /* (I) Forward last zero          */
                 double      AtmPrice0,   /* (I) ATM forward price          */
                 int         t,           /* (I) Current time point         */
                 MQDATA      *mq,         /* (I) Vanilla data structure     */
                 MKTVOL_DATA *mvd,        /* (I) market vol data structure  */
                 TREE_DATA   *tree_data)  /* (I/)) Tree data structure      */
{

    /* Mapped swap yield information */
    double  *YProb_OU = NULL;     /* OU probability up to given level  */
    double  *YCumP_OU = NULL;     /* OU probability up to given level  */
    double  *YCumP_AA = NULL;     /* AA cum prob up to given level     */
    double  *YProb_AA = NULL;     /* AA probability of given level     */
    double  *YStrk    = NULL;     /* corresponding point in Y space    */
    
    /* 2-factor variables  */
    double *StPriceG  = NULL, 
           *AnnuityG  = NULL,
           *YieldG    = NULL;
    
    /* Variables */
    double  AAtoOUFact; 
    double  AASlope;
    double  A_ATM;

    /* Variables for binary swaption algorithm */
    double  binProb;
    double  prevPrice, currPrice;

    /* Mapping counters */
    long    Area     = 0;   
    long    MappSize = 0;

    /* NR variables */
    double  AtmVol,
            Call1, 
            Call2, 
            der, 
            diff,
            a = -0.10,  /* Initialization */ 
            target, 
            fwd, 
            HiYield, 
            Expiry, 
            TotVol,  
            TreeVol,
            VanlVol,
            ErrVol;

    long    CurrNode = 0,
            Iter = 0;

    /* Limits and indices */
    int     Top1, Bottom1;                  /* Tree limits (1rst dim)      */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)       */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)       */ 
    int     LMin=0, LMax=0;

    int     EDevIdx;                        /* Express DEV index           */
    int     i, j, k;                        /* Node indices                */
    int     offset;                         /* Node offset                 */
    int     l;                              /* Level index                 */
    int     s;                              /* Stat index                  */
    int     status = FAILURE;               /* Error status                */
   
    

    /* Get state price slice index */
    EDevIdx = GetDLOffset (tree_data->NbEDevDates,
                           tree_data->EDevDate,
                           tree_data->TPDate[t],
                           CbkEXACT);
    if (EDevIdx < 0)
    {
        DR_Error ("State prices not available for current numeraire date %ld\n",
                  tree_data->TPDate[t]);
        return(status);
    }


    /* Deterministic case */
    if (t == 0)
    {
        if (Set_Slice (*Yield,
                       ParYield0,
                       t,
                       tree_data) == FAILURE)
        {
            goto RETURN;
        }

        return (SUCCESS);
    }




    /* slices size */
    Area = tree_data->Width[0] * tree_data->Width[1];
    
    /* Allocate memory */
    StPriceG    = (double *) DR_Array (DOUBLE, 0, Area);
    AnnuityG    = (double *) DR_Array (DOUBLE, 0, Area);
    YieldG      = (double *) DR_Array (DOUBLE, 0, Area);  
    
   
    /* Binary swaption algorithm */
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    /* Local variables */
    AAtoOUFact = TermZero0 / Annuity0;


    /* 1 FACTOR CASE */  
    if (tree_data->NbFactor == 1)
    {
        DR_Error("Only 2-factor supported!");
        goto RETURN;
    }
        
    /* 2 FACTOR CASE */
    else if (tree_data->NbFactor == 2)
    {
        /* Project the 2D slices on a 1D grid 
           given the mapping direction        */
        if ( SortTree_Slices(AnnuityG,      
                             NULL,  
                             NULL,
                             NULL,
                             StPriceG,
                             &MappSize,
                             1,                    /* 1 Slice to project */
                             Annuity,
                             NULL,
                             NULL,
                             NULL,
                             weights,
                             EDevIdx,
                             t,
                             tree_data) == FAILURE)
        {
            goto RETURN;
        }

        
        /* Initialize cuotff parameters */      
        YCumP_AA = (double *) DR_Array (DOUBLE, 0, MappSize);
        YProb_AA = (double *) DR_Array (DOUBLE, 0, MappSize);
        YStrk    = (double *) DR_Array (DOUBLE, 0, MappSize);
        YProb_OU = (double *) DR_Array (DOUBLE, 0, MappSize);
        YCumP_OU = (double *) DR_Array (DOUBLE, 0, MappSize);
        

        if ((YCumP_AA == NULL) ||
            (YProb_AA == NULL) ||
            (YStrk    == NULL) ||
            (YCumP_OU == NULL) ||
            (YProb_OU == NULL)) goto RETURN;        

        for (l = 0; l < MappSize; l++) YProb_AA[l] = 0;
        
        /* Binaries */
        for (i = 0; i < MappSize; i++)
        {
            YProb_AA[i] += StPriceG[i] * AAtoOUFact * AnnuityG[i];
            YProb_OU[i] += StPriceG[i];
        }
        
        /* Cumulative */
        YCumP_AA[0]  = YProb_AA[0];
        YCumP_OU[0] = YProb_OU[0];
        for (l = 1; l < MappSize; l++) 
        {
            YCumP_AA[l]  = YCumP_AA[l-1]  + YProb_AA[l];
            YCumP_OU[l] = YCumP_OU[l-1] + YProb_OU[l];
        } 
        if (fabs(YCumP_AA[MappSize - 1]-1.)>1e-4)
        {
            DR_Error ("Yield_t: state prices don't add up to 1 at date %ld!",
                      tree_data->TPDate[t]);
            goto RETURN;
        }

        /* search for cutoff levels based on cum probability */
        LMin = 0;
        for (l = 0; l < MappSize - 1; l++)
        {
            if (YCumP_OU[l]  < NormalH (-mvd->NbSigmaBinL)) LMin = l+1;
            else break;
        }

        LMax = MappSize - 1;
        for (l = MappSize-1; l > 0; l--)
        {
            if (YCumP_OU[l] > NormalH (mvd->NbSigmaBinR)) LMax = l-1;
            else break;
        }

        /* Compute binary level */
        for (l = LMin; l <= LMax; l++)
        {
            if (Q3MQMap (mq,
                         mq->muMQ + mq->sigMQ * Normal_InvH (YCumP_AA[l]),
                         &YStrk[l]) == FAILURE)
            { 
                 DR_Error ("Yield_t : Could not compute MqMap\n");
                goto RETURN;
            }
        }


        /* Find cutoff at 99 %-delta */
        Expiry  = Daysact(tree_data->TPDate[0], tree_data->TPDate[t]) / 365.0;      
        AtmVol  = ImpVol_BSQ (ParYield0,ParYield0, Expiry, AtmPrice0, 'C', 1.0, mq->atmTol);            
        TotVol  = AtmVol * sqrt(Expiry);
        HiYield = exp(Normal_InvH(CUTOFF_DELTA)*TotVol) * ParYield0;
            

        /* Yield on the given level */
        for (l = LMin; l <= LMax-1; l++)
        {
            YieldG[l] = YStrk[l];

            if ((l > LMin) && (YieldG[l] < YieldG[l-1]))
            {
                YieldG[l] = YieldG[l-1];
            }

            if (YieldG[l] > HiYield)
            {
                break;
            }            
            CurrNode = l;
        } 

        /* Flat left Yield extrapolation  */
        for (l = 0; l<LMin; l++) YieldG[l] = YieldG[LMin];
        
        /* Semi-parametric TMX : right tail extrapolation */
        target  = mvd->VanlPr[0][NBSTRIKE-1];                               
        do
        {
            /* Initialization */
            Call1 = 0.0;
            Call2 = 0.0; 

            /* 99-delta Call */ 
            for (l = LMin; l < CurrNode; l++)
            {
                Call1      += YProb_AA[l] * 
                              MAX(YieldG[l] - mvd->Strike[0][NBSTRIKE-1], 0);
            }
            
            for (l = CurrNode; l < MappSize; l++)
            {
                /* Functional form */
                YieldG[l]   = YieldG[CurrNode] * 
                              (1.0 + (  exp(a * (l - CurrNode)) - 1.0 )   /a );
                Call1      += YProb_AA[l] * 
                               MAX(YieldG[l] - mvd->Strike[0][NBSTRIKE-1], 0);
            }

            /* Bias */
            diff = Call1 - target;


            /* Convergence reached */
            if (fabs(diff) < TINY)
                break;

            a += SIGMA_DELTA;

            /* Find Gradient */
            for (l = LMin; l < CurrNode; l++)
            {
                Call2      += YProb_AA[l] * MAX(YieldG[l] - mvd->Strike[0][NBSTRIKE-1], 0);
            }
            
            for (l = CurrNode; l < MappSize; l++)
            {
                YieldG[l]   = YieldG[CurrNode] * 
                              (1.0 + (  exp(a * (l - CurrNode)) - 1.0 )   /a );
                Call2      += YProb_AA[l] * 
                              MAX(YieldG[l] - mvd->Strike[0][NBSTRIKE-1], 0);
            }

            /* Gradient */
            der = (Call2 - Call1) / SIGMA_DELTA;

            
            /* Newton-Raphson */
            if (fabs(der) < TINY)
            {
                a += SIGMA_DELTA;
            }
            else
            { 
                a -= diff / der;
            }

            /*  Incremented Iterations */
            Iter++;

        }
        while (Iter < MAX_STEPS);


        /* Now check convergence */
        VanlVol = ImpVol_BSQ(ParYield0,
                             mvd->Strike[0][NBSTRIKE-1],
                             Expiry,
                             target,
                             'C',
                             1, /*  1, 0, */
                             AtmVol);
        
        TreeVol = ImpVol_BSQ(ParYield0,
                             mvd->Strike[0][NBSTRIKE-1],
                             Expiry,
                             ((Iter < MAX_STEPS) ? Call1 : Call2),
                             'C',
                             1, /*  1, 0, */
                             AtmVol);
        
        /* Calibration error */
        ErrVol = fabs(TreeVol-VanlVol);

        /* Check ErrVol is less than MAX_ERR */
        if (fabs(ErrVol) > MAX_ERR)
        {
            DR_Error("Yield_t : Right tail extrapolation failed!");
            goto RETURN;
        }

        /* Now make sure the forward is exactly repriced */
        fwd  = 0.0; 
        for (l = 0; l < MappSize; l++)
        {
            fwd += YProb_AA[l] * YieldG[l];
        }

        /* Adjustment to reprice the forward */
        for (l = 0; l < MappSize; l++)
        {
            YieldG[l] *= ParYield0 / fwd;
        }

        

        /* Compute AA tree price stats */        
        A_ATM = 0.0;        
        for (i = 0; i < MappSize; i++)
        {
            double AStP = AnnuityG[i] * StPriceG[i];
            A_ATM      += AStP * MAX (YieldG[i] - ParYield0, 0);
        }
        A_ATM *= AAtoOUFact;
        

        /* Calibrate ATM price */
        AASlope = AtmPrice0 / A_ATM; 
        for (l = 0; l < MappSize; l++)
        {
            YieldG[l] = ParYield0 + AASlope * (YieldG[l] - ParYield0);
        }        

        /* Now invert the mapping  */
        if (InvMap(*Yield,
                   YieldG,
                   MappSize,
                   weights,
                   EDevIdx,
                   t,
                   tree_data) == FAILURE)
        {
            goto RETURN;
        }

        /* Calc stats */
        for (s = 0; s < NBSTRIKE; s++) 
        {
            mvd->TreePr[0][s] = 0.0;
        }   
        
        for (i = 0; i < MappSize; i++) 
        {
            double AStP = AnnuityG[i] * StPriceG[i];
            for (s = 0; s < NBSTRIKE; s++) 
            {
                mvd->TreePr[0][s] += AStP * MAX(YieldG[i]-mvd->Strike[0][s],0.0) * AAtoOUFact;
            }
        }
   }    
    else if (tree_data->NbFactor == 3)
    {
        goto RETURN;
    }
    

    status = SUCCESS;

    RETURN:

    Free_DR_Array (YProb_AA, DOUBLE, 0, MappSize);
    Free_DR_Array (YCumP_AA, DOUBLE, 0, MappSize);
    Free_DR_Array (YStrk,    DOUBLE, 0, MappSize);
    Free_DR_Array (YProb_OU, DOUBLE, 0, MappSize);
    Free_DR_Array (YCumP_OU, DOUBLE, 0, MappSize);    
    
    
    Free_DR_Array(StPriceG, DOUBLE, 0, Area);
    Free_DR_Array(AnnuityG, DOUBLE, 0, Area);
    Free_DR_Array(YieldG,   DOUBLE, 0, Area);

    
    if (status != SUCCESS)
    {        
        DR_Error("Yield_t: Failed!");
    }

    return (status);

} /* Yield_t */




/*****  Yield2_Vanl **********************************************************/
/*
*       Calculate vanilla smile info for Yield2 and calibrate MQ strcuture.
*/
int Yield2_Vanl (double      *ParYield0,      /* (I) Par Yield from curve    */
                 double      *Annuity0,       /* (I) Annuity from curve      */  
                 double      *TermZero0,      /* (I) Terminal Zero           */
                 double      *StZero0,        /* (I) Zero to current date    */              
                 double      *AtmPrice0,      /* (I) Atm Price               */  
                 long        ValueDate,       /* (I) Value date              */
                 long        TermDate,        /* (I) Terminal date           */
                 long        SwapSt,          /* (I) Swap start date         */
                 long        SwapMat,         /* (I) Swap Maturity           */
                 char        Freq,            /* (I) Swap Frequency          */
                 char        DCC,             /* (I) Swap DCC                */
                 T_CURVE     *t_curve,        /* (I) Zero curve data         */
                 MQDATA      *NmrMQ,          /* (I) MultiQ for vanilla      */
                 MKTVOL_DATA *mvd,            /* (I/O) Mkt vol data          */
                 long        CurrentDate,     /* (I) Current time            */
                 TREE_DATA   *tree_data)      /* (I) Tree data               */
{


    int     status = FAILURE; /* Error status = FAILURE initially */  
    double  Vol[NBVOLPARS];     
    double  *SwapStT=NULL,*TenorT=NULL;    
    double  Expiry;
    double  TotVol;
    int     NbExpiry;    
    int     i, j, s;
    double  Delta[NBSTRIKE] = {0.05,0.25,0.5,0.75,0.99};
    int     CvDiff;    /* Diffuse curve number */ 
    double  NmrSmile[NBVOLPARS+4]; /* Extra space for NbSig,Nck,dN,tN */
    
    

    NbExpiry = mvd->NbVol;
    CvDiff   = tree_data->CvDiff;



    /* Compute terminal zero */
    *TermZero0 = ZeroPrice (TermDate,
                           ValueDate,
                           t_curve[CvDiff].NbZero,
                           t_curve[CvDiff].ZeroDate,
                           t_curve[CvDiff].Zero);

    /* Interpolate on time */
    SwapStT  = (double *) DR_Array(DOUBLE,0,NbExpiry-1);
    TenorT   = (double *) DR_Array(DOUBLE,0,NbExpiry-1);
    for (i=0; i<NbExpiry; i++)
    {
        SwapStT[i]  = Daysact (ValueDate, mvd->SwapSt[i])/365.;
        TenorT[i]   = (double) Months360(mvd->SwapSt[i], mvd->SwapMat[i]);
    }
    
    /* Expiry */
    Expiry = Daysact(ValueDate, SwapSt) / 365.;
    
        
    /* Interpolate swap info on current date */ 
    for (s=0; s<NBVOLPARS; s++)
    {
        tableinterp(Expiry,
                    &(Vol[s]),
                    SwapStT,
                    mvd->Vol[s],
                    NbExpiry);  
    }
    
    /* Compute forward from curve  */   
    if (Par_Yield_From_Dates (ParYield0,
                              Annuity0,
                              SwapSt, 
                              SwapMat,       
                              DCC,           
                              Freq,          
                              'F',  
                              t_curve[CvDiff].NbZero,        
                              t_curve[CvDiff].Zero,         
                              t_curve[CvDiff].ZeroDate,     
                              ValueDate) == FAILURE)
    {
        goto RETURN;
    }
    
    /* Compute discount to current date  */
    *StZero0 = ZeroPrice (CurrentDate,
                          ValueDate,
                          t_curve[CvDiff].NbZero,
                          t_curve[CvDiff].ZeroDate,
                          t_curve[CvDiff].Zero);
    
    /* ATM price for Yield2 calib */
    *AtmPrice0 = Call_BSQ (*ParYield0,
                           *ParYield0,
                           Expiry,
                           Vol[0],
                           1.); 
    
        
                          
    /* MultiQ: fill numerical tree parameters and 
       * load all into smile vector */
    NmrSmile[NBVOLPARS]   = mvd->NbSigmaMQ;
    NmrSmile[NBVOLPARS+1] = mvd->NckMQ;
    NmrSmile[NBVOLPARS+2] = mvd->DeltaNMQ;
    NmrSmile[NBVOLPARS+3] = mvd->TauNMQ;
    
    for (s=0; s<NBVOLPARS; s++) NmrSmile[s] = Vol[s];

    if (CurrentDate > tree_data->TPDate[0])
    {
        /* Calibrate internal MultiQ distribution */
        if (Q3SVToMQ (*ParYield0,
                      Vol[0],
                      Expiry,
                      NmrSmile+1,
                      NmrMQ) == FAILURE)
        {
            DR_Error ("Calibration of MultiQ parameters failed "
                      "at date %d.\n", CurrentDate);
            goto RETURN;
        }
        
        /* Vanila stats */           
        TotVol = Vol[0] * sqrt(Expiry);
        for (s=0; s<NBSTRIKE; s++)
        {
            mvd->Strike[0][s] = exp(Normal_InvH(Delta[s])*TotVol)
                                * (*ParYield0);
            
            if (Q3MQPricer (NmrMQ,
                            Q3_CALL,
                            mvd->Strike[0][s],
                            &(mvd->VanlPr[0][s])) == FAILURE)
            {
                goto RETURN;
            }
        }
    }



    status = SUCCESS;

RETURN:


    Free_DR_Array (SwapStT, DOUBLE, 0, NbExpiry-1);
    Free_DR_Array (TenorT,  DOUBLE, 0, NbExpiry-1);

    if (status == FAILURE)
    {
        DR_Error("Yield2_Vanl: Failed!");
    }

    return (status);

}  /* Yield2_Vanl */



/*****  Yield_Calc   ***************************************************************/
/*
*       Calculate the Yield2 Slice : Newton-Raphson algorithm
*/
int     Yield_Calc (double      **Yield,    /* (O) : Calibrated Yield2 slice  */     
                    double      *Correl,    /* (O) : Realized correlation     */
                    double      Map_dir[2], /* (O) : gaussian weights for map */
                    double      *Annuity,   /* (I) : Annuity slice            */
                    double      *NmrInv,    /* (I) : Inv nmr slice            */
                    double      *Index0,    /* (I) : first rate mapped        */
                    double      Target,     /* (I) : Correlation target       */
                    long        SwapSt,     /* (I) : Swap start date          */
                    long        IndexMat,   /* (I) : Index Maturity           */
                    char        Freq,       /* (I) : Swap frequency           */
                    char        DCC,        /* (I) : Swap DCC                 */
                    long        TermDate,   /* (I) : Terminal Date            */
                    long        t,          /* (I) : current time             */
                    T_CURVE     *t_curve,   /* (I) : T curve structure        */
                    MKTVOL_DATA *mvd,       /* (I) : Market Vol structure     */
                    TREE_DATA   *tree_data) /* (I) : Tree data structure      */
{   
    int s;

    /* Swap details */
    double      ParYield0, 
                Annuity0, 
                TermZero0, 
                Zero0, 
                AtmPrice0;
    
    /* Annuity Hat slice */
    double      *AnnuityH;

    /* Local array to store the weights */
    double      Map_dirL[2];
    
    /* Value Date */
    long        ValueDate = t_curve[tree_data->CvDiff].ValueDate, SwapMat;        

    /* Newton-Raphson variables */
    double      diff = 1.0,
                Mindiff = 1.0,
                der,
                corr,
                corr_d,
                MinCorr = 1.0;

    /* Dummy variables */
    double     cov, 
               vol1, 
               vol2, 
               fwd1, 
               fwd2;
    
    /* Number of iterations */
    long        Iter = 0;

    /* Printout variables */
    double VanlVol[NBSTRIKE];
    double TreeVol[NBSTRIKE];
    double Expiry;
    double ErrVol=0;
    double ErrAtm=0;
    char    FileName[MAXBUFF];
    FILE    *stream = NULL;

    /* MQ structure */
    MQDATA      *NmrMQ = NULL;              /* Local MultiQ structure for vanilla    */
    int         status = FAILURE;
    
    /* Allocate MQ structures */
    NmrMQ = (MQDATA *) calloc (1, sizeof(MQDATA));
    SwapMat = Nxtmth(SwapSt,IndexMat,1L);

    /* Initialize the Annuity-Hat structure */  
    AnnuityH = Alloc_Slice(tree_data);


    /* Initialize the mapping */
    Map_dir[0] = 0.0;
    Map_dir[1] = 1.0;


    /* Initialize the local variables (used if the NR does not converge after 20 iterations */
    Map_dirL[0] = Map_dir[0];
    Map_dirL[1] = Map_dir[1];

    /* First Calculate the Vanilla information and calibrate the MultiQ structure */
    if (Yield2_Vanl (&ParYield0,      
                     &Annuity0,         
                     &TermZero0,      
                     &Zero0,                         
                     &AtmPrice0,        
                     ValueDate,       
                     TermDate, 
                     SwapSt,          
                     SwapMat,         
                     Freq,            
                     DCC,             
                     t_curve,        
                     NmrMQ,          
                     mvd,            
                     tree_data->TPDate[t],               
                     tree_data) == FAILURE)
    {
        DR_Error("Can't calibrate the vanilla structure!");
        goto RETURN;
    }


    /* Get Annuity-Hat slice for measure-change */ 
    if (MultiplyTwoSlices (AnnuityH,
                           Annuity,
                           NmrInv,
                           t,
                           tree_data) == FAILURE) goto RETURN;


    /* Start Newton-Raphson */
    do
    {
        /* Calibrate the Yield */
        if (Yield_t(Yield,
                    Map_dir,
                    AnnuityH,     
                    ParYield0,   
                    Annuity0,    
                    TermZero0,   
                    AtmPrice0,  
                    t,           
                    NmrMQ,         
                    mvd,        
                    tree_data) == FAILURE)
        {
            DR_Error("Can't calibrate the vanilla structure!");
            goto RETURN;
        }
        
        /* Now calculate and store instantaneous rate cross-correlations */
        if (Covar_Slices (&corr,
                          &cov,
                          &fwd1,
                          &fwd2, 
                          &vol1,
                          &vol2,
                          *Yield,
                          Index0,
                          NmrInv,
                          TermDate,
                          t,
                          &t_curve[tree_data->CvDiff], 
                          tree_data) == FAILURE)
        {
            goto RETURN;
        }       

        
        /* Correlation error */
        diff = corr - Target;

        /* If the calibration error is reduced , store the value of the mapping weights */
        if (fabs(diff) < fabs(Mindiff))
        {
            Map_dirL[0] = Map_dir[0];
            Map_dirL[1] = Map_dir[1];
            MinCorr     = corr;
        }


        if (fabs(diff) < CORR_ERROR)
            break;

        Map_dir[0] += DELTA_CORR;

        /* Calibrate the Yield */
        if (Yield_t(Yield,
                    Map_dir,
                    AnnuityH,     
                    ParYield0,   
                    Annuity0,    
                    TermZero0,   
                    AtmPrice0,  
                    t,           
                    NmrMQ,         
                    mvd,        
                    tree_data) == FAILURE)
        {
            DR_Error("Can't calibrate the vanilla structure!");
            goto RETURN;
        }
        
        /* Now calculate and store instantaneous rate cross-ccorrelations */
        if (Covar_Slices (&corr_d,
                          &cov,
                          &fwd1,
                          &fwd2, 
                          &vol1,
                          &vol2,
                          *Yield,
                          Index0,
                          NmrInv,
                          TermDate,
                          t,
                          &t_curve[tree_data->CvDiff], 
                          tree_data) == FAILURE)
        {
            goto RETURN;
        }       


        Map_dir[0] -= DELTA_CORR;

        
        
        /* NR derivative */     
        der  = (corr_d - corr) / DELTA_CORR;
        
        if ( fabs(der) < TINY || fabs(diff / der) > 1.0 / DELTA_CORR)
        {
            Map_dir[0] -= diff ;
        }
        else
        {
            Map_dir[0] -= diff/ der;
        }

        /* Minimum correlation error so far */
        Mindiff = MIN(fabs(Mindiff), fabs(diff));

        /* Increment number of iterations */
        Iter ++;
    }
    while ((fabs(diff) > CORR_ERROR) && (Iter < MAX_ITERATIONS));


    /* if the NR has not converged ... */
    if (Iter >= MAX_ITERATIONS)
    {
        Map_dir[0] = Map_dirL[0];
        Map_dir[1] = Map_dirL[1];

        /* Calibrate the Yield */
        if (Yield_t(Yield,
                    Map_dir,
                    AnnuityH,     
                    ParYield0,   
                    Annuity0,    
                    TermZero0,   
                    AtmPrice0,  
                    t,           
                    NmrMQ,         
                    mvd,        
                    tree_data) == FAILURE)
        {
            DR_Error("Can't calibrate the vanilla structure!");
            goto RETURN;
        }
    }
    /* Output */
    *Correl = MinCorr;

    /* if t = 0, nothing to print */
    if (t == 0)
    {
        status = SUCCESS;
        goto RETURN;
    }



    /* Print smile info on the secopnd yield */
    strcpy (FileName, "Yield.prn");
    if ((stream = fopen (FileName, "a")) == NULL)
    {
        DR_Error ("Could not open file %s.\n", FileName);
        goto RETURN;
    }
    
    /* Print smile information now  */
    fprintf (stream, "###### MAPPED YIELD SMILE INFO ######\n"); 
    fprintf (stream, "# Date  Yield tenor      Strike/VanlVol/TreeVol   "
             "                       ErrSmile    ErrAtm\n");
    
    Expiry = Daysact(tree_data->TPDate[0], tree_data->TPDate[t]) / 365.;
    
    for (s=0; s<NBSTRIKE; s++)
    {
        VanlVol[s] = ImpVol_BSQ (ParYield0,
                                 mvd->Strike[0][s],
                                 Expiry,
                                 mvd->VanlPr[0][s],
                                 'C',
                                 1, /*  1, 0, */ 
                                 0.20);  
        
        TreeVol[s] = ImpVol_BSQ (ParYield0,
                                 mvd->Strike[0][s],
                                 Expiry,
                                 mvd->TreePr[0][s],
                                 'C',
                                 1, /* 1, 0, */
                                 0.20);
        
        ErrVol = MAX(ErrVol,fabs(TreeVol[s]-VanlVol[s]));
    }
    
    /* Use the fact that s=2 corresponds to ATM strike */
    ErrAtm = MAX (ErrAtm, fabs (TreeVol[2]-VanlVol[2]));
    
    fprintf (stream, "%8ld  %8ld", tree_data->TPDate[t], (long) (IndexMat / 12.0));
    for (s=0;s<NBSTRIKE;s++) 
        fprintf(stream,"%8.4f  ",mvd->Strike[0][s]*100);
    fprintf(stream, "%8.2f    %8.2f",ErrVol*100, ErrAtm*100);
    fprintf(stream, "\n");
    
    fprintf (stream, "                  ");
    for (s=0;s<NBSTRIKE;s++)
        fprintf(stream,"%8.2f  ",VanlVol[s]*100);
    fprintf(stream, "\n");
    
    fprintf (stream, "                  ");
    for (s=0;s<NBSTRIKE;s++)
        fprintf(stream,"%8.2f  ",TreeVol[s]*100);
    fprintf(stream, "\n");  
    
    status = SUCCESS;

RETURN:

    if (NmrMQ != NULL) free(NmrMQ);
    Free_Slice(AnnuityH, tree_data);

    if (stream != NULL)
    {
        fclose (stream);
    }
    
    if (status != SUCCESS)
    {        
        DR_Error("Yield_Calc: Failed!");
    }

    return (status);
}



/*****  Covar_Slices  ********************************************************/
/*
*       Calculates the cov between 2 slices
*/
int Covar_Slices (double     *Corr,                /* Correlation between 2 slices */
                  double     *Covar,               /* Covar between 2 slices       */
                  double     *Fwd1,                /* expectation of first slice   */
                  double     *Fwd2,                /* ex[pectation of second slice */
                  double     *Vol1,                /* Vol of first rate            */
                  double     *Vol2,                /* Vol of second rate           */
                  double     *Slice1,              /* First slice                  */
                  double     *Slice2,              /* second slice                 */
                  double     *InvNmr,              /* Inverse numeraire            */
                  long       TerminalDate,         /* Tree terminal date           */
                  long       t,                    /* current time in the tree     */
                  T_CURVE    *t_curve,             /* t curve structure data       */
                  TREE_DATA  *tree_data)
{
    long   i, k;
    double correl,
           cov    = 0.0,
           exp1   = 0.0,
           exp2   = 0.0,
           var1   = 0.0,
           var2   = 0.0;
    long   EDevIdx;
    double *StatePr, *Slice1L, *Slice2L, *NmrInvL;  /* Local pointers */
    double ZerotoT, ZerotoTerm;
    
    int    status = FAILURE;


    /* Get state price slice index */
    EDevIdx = GetDLOffset (tree_data->NbEDevDates,
                           tree_data->EDevDate,
                           tree_data->TPDate[t],
                           CbkEXACT);

    if ( EDevIdx < 0)
        goto RETURN;

    /* Zero to current date */
    ZerotoT = ZeroPrice (tree_data->TPDate[t],
                         t_curve->ValueDate,
                         t_curve->NbZero,
                         t_curve->ZeroDate,
                         t_curve->Zero);

    /* Zero to terminal date */
    ZerotoTerm = ZeroPrice (TerminalDate,
                            t_curve->ValueDate,
                            t_curve->NbZero,
                            t_curve->ZeroDate,
                            t_curve->Zero);

    for (i = tree_data->Bottom1[t]; i <= tree_data->Top1[t]; i++)
    {
        Slice1L    =  Slice1                          + Node_Offset(2, i, 0, t, tree_data);;
        StatePr    =  tree_data->EDevStPrice[EDevIdx] + Node_Offset(2, i, 0, t, tree_data);;
        Slice2L    =  Slice2                          + Node_Offset(2, i, 0, t, tree_data);;    
        NmrInvL    =  InvNmr                          + Node_Offset(2, i, 0, t, tree_data);; 

        for (k = tree_data->Bottom2[t][i]; k <= tree_data->Top2[t][i]; k++)
        {
            exp1 += Slice1L[k] * StatePr[k] * NmrInvL[k] * ZerotoTerm / ZerotoT              ;
            exp2 += Slice2L[k] * StatePr[k] * NmrInvL[k] * ZerotoTerm / ZerotoT              ;
            cov  += Slice1L[k] * StatePr[k] * Slice2L[k] * NmrInvL[k] * ZerotoTerm / ZerotoT ;
            var1 += Slice1L[k] * StatePr[k] * Slice1L[k] * NmrInvL[k] * ZerotoTerm / ZerotoT ;
            var2 += Slice2L[k] * StatePr[k] * Slice2L[k] * NmrInvL[k] * ZerotoTerm / ZerotoT ;          
        }
    }

    cov    -= exp1 * exp2             ;
    var1   -= exp1 * exp1             ;
    var2   -= exp2 * exp2             ;
    correl = cov/ (sqrt(var1 * var2)) ;

    *Corr   = correl;
    *Covar  = cov;
    *Fwd1   = exp1;
    *Fwd2   = exp2;
    *Vol1   = sqrt(var1);
    *Vol2   = sqrt(var2);


    status = SUCCESS;


RETURN:
    return(status);

}  /* Covar_Slices */













