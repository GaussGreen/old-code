/****************************************************************************/
/*      TMX3: Numeraire calculation and utilities                           */
/****************************************************************************/
/*      numer.c                                                             */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fix123head.h"
#include "q3_tmx.h"


static int Nmr_t(double **,double *,int,int,MQDATA *,MKTVOL_DATA *,FIX3_TREE_DATA *);                  
static int Nmr_Print0(MKTVOL_DATA *);
static int Nmr_Print1(MKTVOL_DATA *,MQDATA *);
static int Nmr_Vanl(MQDATA *,MKTVOL_DATA *, FIX3_TREE_DATA *, T_CURVE const*);
static int Nmr_Init(MKTVOL_DATA *);


#define  AA_FWD            0
#define  AA_ATM            1
#define  AA_ALL            4


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
    int         i;
    long        TermDate, 
                LibMatDate;  
    int         intvl, NbExpiry, NbNmr; 
    CRIT_DATE   *CritDate = NULL;             /* Critical date list               */
    int         NbCritDate = 0;               /* Number of critical dates         */   
    
    int     status = FAILURE;                 /* Error status = FAILURE initially */  


    /* Allocate memory for Crit Date structure */
    CritDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
    if (CritDate == NULL)
    {
        DR_Error("Nmr_Schedule: unable to allocate memory "
                 "for critical dates array !");
        goto RETURN;
    }


    /* Initialize mktvol_data structure first */
    if (Nmr_Init(mvd) == FAILURE) goto RETURN;


    /* Choose  benchmarks: last one must expire AFTER LastProdDate */
    TermDate = LastProdDate;
    for (i = 0; i < mvd->NbVol; i++)
    {
        TermDate = MAX(TermDate, mvd->SwapMat[i]);
        if (mvd->SwapSt[i] >= LastProdDate) break;
    }

    /* Number of benchmarks actually used */
    NbExpiry = MIN(i+1, mvd->NbVol);


    /*Overwite NbVol with number of effective benchmarks */
    mvd->NbVol = NbExpiry;

    /* Terminal date is obtained by shitfting the lastprod/vol date by two mapping intervals */
    /* 2 periods to cater for 31st Mar maturities where the mat could fall in the last Nmr   */
    /* interval if we only add one - we cannot interpolate in the last interval              */
    if (Conv_Freq(mvd->Freq) == 0)
    {
        DR_Error("Internal error - result of Conv_Freq(char Freq) function (with argument %c) = 0, "
                 "- Invalid/incorrect Freq argument", mvd->Freq);
        goto RETURN;
    }

    intvl         = 12 / Conv_Freq(mvd->Freq);
    TermDate      = Nxtmth(TermDate, 2*intvl,1L);
    mvd->TermDate = TermDate;

    
    /* Always add value date to critical date list */
    if (Add_To_DateList (&NbCritDate,
                         &CritDate,              
                         ValueDate,
                         NBCRITDATE,    /* No specific type */
                         0, 0, 0, 0, 0,              
                         0, 0, 0) == FAILURE)
    {
        goto RETURN;
    }

    
    /* Add Term date to critical date list */
    if (Add_To_DateList (&NbCritDate,
                         &CritDate,              
                         TermDate,
                         NBCRITDATE,    /* No specific type */
                         0, 0, 0, 0, 0,              
                         0, 0, 0) == FAILURE)
    {
        goto RETURN;
    }


    
    /* Add Benchmark dates now */
    for (i=0; i < NbExpiry; i++)
    {
        if (Add_To_DateList (&NbCritDate,
                             &CritDate,
                             mvd->SwapSt[i],
                             NBCRITDATE,
                             0, 0, 0, 0, 0,
                             0, 0, 0) == FAILURE) 
        {
            goto RETURN;
        }
    }

    /* Add offsets of Value Date */
    LibMatDate = Nxtmth(ValueDate, intvl,1L);
    
    while (LibMatDate < mvd->SwapSt[0])     
    {
        /* To handle the case where tenor is 1 day - which    */
        /* could potentially result in dccfrac = 0 if DCC = 3 */
        /* & MatDate = 31st. Assumes intvl > 1 day            */
        if ( Daysact(LibMatDate, mvd->SwapSt[0]) <= 2)
            LibMatDate = Nxtday(LibMatDate, -1L);
        

        if (Add_To_DateList (&NbCritDate,
                             &CritDate,
                             LibMatDate,
                             NBCRITDATE,
                             0, 0, 0, 0, 0,
                             0, 0, 0) == FAILURE)       
        {       
            goto RETURN;            
        }
        
        LibMatDate = Nxtmth(LibMatDate, intvl,1L);
    }

    /* Add Offsets of Benchmark Dates */
    for (i=0; i < NbExpiry - 1; i++)
    {
        LibMatDate = mvd->SwapSt[i];

        while (LibMatDate < mvd->SwapSt[i+1]) 
        {
            /* To handle the case where tenor is 1 day - which    */        
            /* could potentially result in dccfrac = 0 if DCC = 3 */        
            /* & MatDate = 31st. Assumes intvl > 1 day            */
            if ( Daysact(LibMatDate, mvd->SwapSt[i+1]) <= 2)
                LibMatDate = Nxtday(LibMatDate, -1L);

            if (Add_To_DateList (&NbCritDate,
                                 &CritDate,
                                 LibMatDate,
                                 NBCRITDATE,
                                 0, 0, 0, 0, 0,
                                 0, 0, 0) == FAILURE)
            {
                goto RETURN;
            }

            LibMatDate = Nxtmth(LibMatDate, intvl,1L);
        }
    }

    /* Add offset of last benchmark date till the tree last date */
    LibMatDate = mvd->SwapSt[NbExpiry - 1];

    while (LibMatDate < TermDate)
    {
        /* To handle the case where tenor is 1 day - which    */
        /* could potentially result in dccfrac = 0 if DCC = 3 */
        /* & MatDate = 31st. Assumes intvl > 1 day            */
        if ( Daysact(LibMatDate, TermDate) <= 2)
            LibMatDate = Nxtday(LibMatDate, -1L);

        if (Add_To_DateList (&NbCritDate,
                             &CritDate,
                             LibMatDate,
                             NBCRITDATE,
                             0, 0, 0, 0, 0,
                             0, 0, 0) == FAILURE)
        {
            goto RETURN;
        }
        LibMatDate = Nxtmth(LibMatDate, intvl,1L);
    }


        
    /* Sort the CritDate */
    if (Sort_CritDate ( NbCritDate,
                        CritDate) == FAILURE)
    {
        goto RETURN;
    }

    /* Copy numeraire dates into mvd structure */
    mvd->NmrDate[0] = ValueDate;
    NbNmr   = 1; 
    for (i=1; i < NbCritDate; i++)
    {
        if (CritDate[i].CritDate != CritDate[i-1].CritDate)
        {
            mvd->NmrDate[NbNmr] = CritDate[i].CritDate;         
            
            /* Increment Number of numeraire dates */
            NbNmr++;
        }
    }

    /* Store Nbr of Nmr dates */
    mvd->NbNmr = NbNmr;

    /* Check the numeraire schedule is consitent with Terminal date */  
    if (mvd->NbNmr > MAXNBDATE)
    {
        DR_Error ("Nb of numeraire dates %d exceeds maximum of %d. "
                  "(Nmr_Dates)\n", NbNmr, MAXNBDATE);
        goto RETURN;
    }
    if (mvd->NmrDate[mvd->NbNmr-1] != mvd->TermDate)
    {
        DR_Error ("Nmr_Schedule: Last Nmr Date is bot consistent with Term date");
        goto RETURN;
    }

    status = SUCCESS;

RETURN:

    Nmr_Print0(mvd);

    Free_DR_Array (CritDate, CRITDATE, 0, NbCritDate-1);

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
    int     i;

    /* Initialize Flags & dates  */
    for (i=0; i < MAXNBDATE; i++)
    {
        mvd->NmrDate[i]         = 999L;
    }

    return SUCCESS;
}


/*****  Nmr_Alloc ********************************************************/
/*
*/
int Nmr_Alloc (MKTVOL_DATA      *mvd,        /* (I/O)  Mkt vol data     */
               FIX3_TREE_DATA   *tree_data)  /* (I/O)  Tree data        */
{
    int   status =  FAILURE; /* Error status = FAILURE initially */
    int   NbNmr  = mvd->NbNmr;
    long  Area   = 0;
    int i;


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

    if (NbNmr < 0) goto RETURN;

    /* Allocate inverse numeraire and dates */
    tree_data->NmrInv  = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 

    /* Allocate annuity and zero slices, and yield mapping arrays */
    tree_data->Libor    = (double **)DR_Array(DOUBLE_PTR,0,NbNmr-1); 
    tree_data->LastZero = Fix3_Alloc_Slice(tree_data);

    if ((tree_data->NmrInv   == NULL) ||
        (tree_data->LastZero == NULL) ||
        (tree_data->Libor    == NULL) ) goto RETURN;

    /* Initialize pointers to NULL for safe freeing */
    for (i=0; i<NbNmr; i++)
    {
        tree_data->NmrInv[i]  = NULL;
        tree_data->Libor[i]   = NULL;
    }

    for (i=0; i<NbNmr; i++)
    {
        tree_data->NmrInv[i]  = Fix3_Alloc_Slice(tree_data);

        if (tree_data->NmrInv[i]  == NULL) goto RETURN;

        tree_data->Libor[i] = (double *)DR_Array (DOUBLE,-1,Area);
        if (tree_data->Libor[i] == NULL) goto RETURN;
    }

    /* Copy Nb of Nmrs to tree structure */
    tree_data->NbNmr = NbNmr;

    /* Initialise Libor */
    for (i = 0; i < NbNmr; i++)
    {
        mvd->AtmLiborPrice0[i] = 0.0;
        mvd->ParLibor0[i]      = 0.0;
        mvd->AnnuityLibor0[i]  = 0.0;
    }

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("Nmr_Alloc: Failed!");
    }

    return (status);
}

#define TMX_6Q 6
/*****  Nmr_Vanl **********************************************************/
/*
*    Maps an existing Libor SV-MQ distribution to a MQ distribution & 
*    Calcs Libor Vanilla stats based on that mapping
*/
int Nmr_Vanl (MQDATA         *NmrMQ,       /* (O) MultiQ for vanilla      */
              MKTVOL_DATA    *mvd,         /* (I/O) Mkt vol data          */
              FIX3_TREE_DATA *tree_data,   /* (I) Tree data structure     */
              T_CURVE const  *t_curve)     /* (I) T_curve data structure  */
{
    int     NbNmr, i, s;
    int     CvDiff;
    int     status = FAILURE;         /* Error status = FAILURE initially */  

    /* MultiQ variables */
    double  NckMQ;
    long key; 
    long    i3;

    /* Initialisations       */
    NbNmr     = mvd->NbNmr;
    NckMQ     = mvd->NckMQ;
    CvDiff    = tree_data->CvDiff;
    

    /* 
     * TMX always uses MultiQ with 6Q ! 
     * Change Numerical Key configuration in MultiQ
     */
    key   = (long) (NckMQ / 10000.);
    i3    = key % 10; 
    NckMQ = NckMQ + (TMX_6Q - i3) * 10000.;  


   /* Compute terminal zero      */
    mvd->TermZero0 = GetZeroPrice(mvd->TermDate, &t_curve[CvDiff]);

    /* Map each Nmr Date */
    for (i = 1; i < NbNmr - 1; i++) 
    {
        double  Expiry, ParLibor0, NmrLiborAtm;
        double  dccfrac;
        double  NmrSmile[NBVOLPARS+4]; /* Extra space for NbSig,Nck,dN,tN */
        
        /* Get Par Libor */     
        if (ParYieldFromDates (&(mvd->ParLibor0[i]),
                               &(mvd->AnnuityLibor0[i]),
                               mvd->NmrDate[i],
                               mvd->NmrDate[i+1],
                               mvd->DCC,
                               mvd->Freq,
                               'F',
                               &t_curve[CvDiff]) == FAILURE)
        {
            goto RETURN;
        }


        /* Initialisations  */
        Expiry      = Daysact(GetZeroCurveBaseDate(&t_curve[CvDiff]), mvd->NmrDate[i]) / 365.;
        ParLibor0   = mvd->ParLibor0[i];
        NmrLiborAtm = mvd->NmrLibVol[0][i];
        
        if (DrDayCountFraction(mvd->NmrDate[i], 
                               mvd->NmrDate[i+1], 
                               mvd->DCC, 
                               &dccfrac) == FAILURE)
        {
            DR_Error("Nmr_Vanl: Failed to compute cash accrual!");
            goto RETURN;
        }
        
        for (s=0; s<NBVOLPARS; s++)
        {
            NmrSmile[s] = mvd->NmrLibVol[s][i];
        }

        mvd->AtmLiborPrice0[i] = Call_BSQ (ParLibor0,
                                           ParLibor0,
                                           Expiry,
                                           NmrLiborAtm,
                                           1.);

        /* MultiQ: fill numerical tree parameters and 
         * load all into smile vector */
        NmrSmile[NBVOLPARS]   = mvd->NbSigmaMQ;
        NmrSmile[NBVOLPARS+1] = NckMQ;
        NmrSmile[NBVOLPARS+2] = Q3TMX_6Q_DELTAN_DEFAULT;
        NmrSmile[NBVOLPARS+3] = 1.0 / (dccfrac * ParLibor0);
        
        /* Calibrate internal MultiQ distribution */
        if (Q3TMXSVToMQ (ParLibor0,
                      NmrLiborAtm,
                      Expiry,
                      NmrSmile+1,
                     &(NmrMQ[i])) == FAILURE)
        {
            DR_Error ("Calibration of MultiQ parameters failed "
                          "at date %ld.\n", mvd->NmrDate[i]);
            DR_Error ("ATM Libor vol: %.2f\n", 100 * NmrLiborAtm);
            DR_Error ("Libor smile  : %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",
                      NmrSmile[1],NmrSmile[2],NmrSmile[3],NmrSmile[4],
                      NmrSmile[5],NmrSmile[6],NmrSmile[7],NmrSmile[8]); 
            goto RETURN;
        }        
 
    } /* for i < NbNmr-1 */
    
    status = SUCCESS;

RETURN:

    Nmr_Print1(mvd, NmrMQ);

    if (status == FAILURE)
    {
        DR_Error("Nmr_Vanl: Failed!");
    }

    return (status);

}  /* Nmr_Vanl */


/*****  Nmr_Calc  **********************************************************/
/*
*       TMX3: Numeraire calculation at HK dates - done in terms of Libor now
*/
int Nmr_Calc (T_CURVE const        *t_curve,       /* (I) Zero curve       */
              MKTVOL_DATA          *mvd,           /* (I) Volatility data  */
              FIX3_TREE_DATA       *tree_data)     /* (O) Tree data        */
{
    FIX3_DEV_DATA dev_data;       /* Dev data structure                    */
    CLAIM_BANK  ZeroBank;         /* Zeroes used in numeraire calculation  */

    MQDATA      *NmrMQ = NULL;    /* Local MultiQ structure for vanilla    */

    /* Slices */
    double      *FirstZero = NULL; /* Shortest maturity zero               */
    double      *newZero   = NULL; /* New zero added to zerobank           */
    double      *NmrInv;           /* This is only pointer, NO ALLOCATION  */ 

    /* Numeraire date info */
    int         NbNmr;            /* Number of numeraire dates             */
    int         NmrIdx;           /* Numeraire date index                  */
    int         NmrFlag;          /* Is current date a numeraire date      */

    int         intvl; 

    /* Numerical variables */
    long        TerminalDate;
    long        CurrentDate;
    long        ErDate;
    int         CvDiff;           /* Diffuse curve number                  */
    int         t;                /* Current time point                    */   
    int         T;                /* Last time point                       */
    int         status = FAILURE; /* Error status = FAILURE initially      */
    
    /* Dates and indices */
    T            = tree_data->NbTP;
    NbNmr        = mvd->NbNmr;
    TerminalDate = mvd->NmrDate[NbNmr-1];

    CvDiff       = tree_data->CvDiff;

    intvl = 12 / Conv_Freq(mvd->Freq);

    /* Initialize and allocate DEV_DATA structure */
    Fix3_Dev_Init (&dev_data);
    if (Fix3_Dev_Alloc (&dev_data, tree_data) == FAILURE) goto RETURN;

    /* Initialize and allocate zero bank */
    Fix3_CbkInit(&ZeroBank);
    if (Fix3_CbkAlloc (&ZeroBank, NbNmr, tree_data) == FAILURE) goto RETURN;

    /* Numeraire cannot be interpolated yet */
    tree_data->NmrInterpOn = FALSE;

    /* Numeraire is not needed in lattice */
    dev_data.NmrToCcy = FALSE;
    dev_data.CcyToNmr = FALSE;

    /* Allocate MQ structures */
    NmrMQ = (MQDATA *) calloc (NbNmr, sizeof(MQDATA));
    if (NmrMQ == NULL) goto RETURN;

    /* Calibrate MQ parameters, price vanilla if NmrStatFlag */
    if (Nmr_Vanl(NmrMQ,
                 mvd,
                 tree_data,
                 t_curve) == FAILURE) goto RETURN;

    /* Allocate variables for numeraire calculation */
    /* Step backward through the tree */
    for (t = T; t >= 0; t--)
    {
        /* Current Date */
        CurrentDate = tree_data->TPDate[t];
        
        /* Set flags */
        NmrFlag = tree_data->TPtype[NMREVENT][t];

        /* 'Update' tree */
        if (Fix3_Lattice_Tmx (
                     &dev_data,
                     t,
                     T,
                     mvd,
                     tree_data) == FAILURE) goto RETURN;

        /* Ev Zero Bank; do not add any new zero yet */
        if (Fix3_CbkDev (
                    &ZeroBank,
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
                FirstZero = Fix3_ZbkReadZero (
                                         &ZeroBank,
                                         mvd->NmrDate[NmrIdx+1],
                                         FALSE,
                                         CurrentDate,
                                         t,
                                         tree_data);

                if (FirstZero == NULL) goto RETURN;
            }

            /* Calculate inverse numeraire */
            if (Nmr_t (&NmrInv,
                       FirstZero,
                       t,
                       NmrIdx,
                       &(NmrMQ[NmrIdx]),
                       mvd,
                       tree_data) == FAILURE) goto RETURN;

            /* Add inverse numeraire as new zero to ZeroBank as in ZbkUpdate */
            newZero = Fix3_CbkPopSlice(&ZeroBank);

            if (newZero == NULL) goto RETURN;

            if (Fix3_Copy_Slice (newZero,
                                 NmrInv, 
                                 t,
                                 tree_data) == FAILURE) goto RETURN;


            /* Find out earliest use of numeraire */
            ErDate = mvd->NmrDate[NmrIdx];
            if (NmrIdx > 0 )
                ErDate = mvd->NmrDate[NmrIdx-1];



            if (Fix3_CbkPushSlice (&ZeroBank,
                                   newZero,
                                   CurrentDate,
                                   ErDate) == FAILURE) goto RETURN;

        } /* if (NmrFlag) */

    } /* for t */

    /* From now on it is possible to intepolate numeraire */
    tree_data->NmrInterpOn = TRUE;

    status = SUCCESS;

    RETURN:

    if (NmrMQ != NULL) free(NmrMQ);

    Fix3_Dev_Free   (&dev_data, tree_data);
    Fix3_CbkFree    (&ZeroBank, tree_data); 

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_Calc: Failed!");
    }

    return (status);

}  /* Nmr_Calc */


/*****  Nmr_t   ***************************************************************/
/*
*       Generate numeraire slice at time t out of vanilla smile information
*       Nmr put directly into tree structure AND output the
*/
int     Nmr_t (  double      **NmrInv,       /* (O) Nmr pointer               */
                 double      *FirstZero,     /* (I) Shortest mat in annuity   */ 
                 int         t,              /* (I) Current time point        */
                 int         NmrIdx,         /* (I) Numeraire index           */
                 MQDATA      *mq,            /* (I) Vanilla data structure    */
                 MKTVOL_DATA *mvd,           /* (I) Volatility data structure */
                 FIX3_TREE_DATA *tree_data)  /* (I/)) Tree data structure     */
{
    /* Local slice pointers */
    int     *LevelL;
    double  *FirstZeroL; 
    double  *NmrInvL;      
    double  *StPriceL;
    double  *AuxSliceL;

    /* Mapped swap yield information */
    double  *LiborL;              /* calibrating libor                 */
    double  *YCumP_OU = NULL;     /* OU cum prob up to given level     */ 
    double  *YProb_OU = NULL;     /* OU probability up to given level  */
    double  *YCumP_AA = NULL;     /* AA cum prob up to given level     */
    double  *YProb_AA = NULL;     /* AA probability of given level     */
    double  *YStrk    = NULL;     /* corresponding point in Y space    */
    int     *Level    = NULL;     /* node Y-level                      */

    /* Variables */
    double  ParLibor0;
    double  AAtoOUFact;
    double  DCCFrac;
    double  AASlope;
    double  AASta[AA_ALL];
    double  Step = 0.0;
    double  *AuxSlice = NULL;     /* Aux slice for smoothing calcs     */

    /* Variables for binary swaption algorithm */
    double  binProb;
    double  prevPrice, currPrice;

    /* Limits and indices */
    int     Top1, Bottom1;                  /* Tree limits (1rst dim)      */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)       */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)       */ 
    int     LMinInit=0, LMaxInit=0;          /* Min/max yield   levels      */
    int     LMin=0, LMax=0;

    int     EDevIdx;                        /* Express DEV index           */
    int     i;                              /* Node indices                */
    int     offset;                         /* Node offset                 */
    long    Area = 0;                       /* Slices size                 */
    int     l;                              /* Level index                 */
    int     s;                              /* Stat index                  */
    int     SmoothingOn;
    int     status = FAILURE;               /* Error status                */
 
    *NmrInv = NULL;
    if (NmrIdx < 0 || NmrIdx >= tree_data->NbNmr) goto RETURN;

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

    /* Set mapping pointers to current Nmr slice in tree_data */
    LiborL = tree_data->Libor[NmrIdx];

    /* Clean up */
    for (s=0; s<AA_ALL; s++) AASta[s] = 0;

    /* Special Nmr dates: 0 and last */
    if (NmrIdx == tree_data->NbNmr-1) 
    {
        if (Fix3_Set_Slice (*NmrInv,
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
        if (Fix3_Set_Slice (*NmrInv,
                            1. / mvd->TermZero0,
                            t,
                            tree_data) == FAILURE)
        {
            goto RETURN;
        }

        /* Mapping is constant */
        LiborL[0]    = mvd->ParLibor0[0];

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
    ParLibor0  = mvd->ParLibor0[NmrIdx];
    AAtoOUFact = mvd->TermZero0 / mvd->AnnuityLibor0[NmrIdx];

    if (DrDayCountFraction(mvd->NmrDate[NmrIdx], 
                           mvd->NmrDate[NmrIdx+1], 
                           mvd->DCC, 
                           &DCCFrac) == FAILURE)
    {
        DR_Error("Nmr_t: Failed to compute cash accrual!");
        goto RETURN;
    }

    /* Smoothing flag */
    SmoothingOn = (mvd->SmoothingFlag == 'Y');

    AuxSlice = Fix3_Alloc_Slice (tree_data);
    
    if (AuxSlice == NULL)
    {
        DR_Error ("Nmr_t: could not allocate memory!");
        goto RETURN;        
    }

    /* Allocate level slice */
    Level = Fix3_Alloc_Slice_Int (tree_data);

    if (Level == NULL)
    {
        DR_Error("Nmr_t: Could not allocate memory.\n");
        goto RETURN;
    }

    if (tree_data->NbFactor == 1)
    {
        offset     = Fix3_Node_Offset(1, 0, 0, t, tree_data);
        LevelL     = Level    + offset;
        FirstZeroL = FirstZero + offset;
        NmrInvL    = tree_data->NmrInv[NmrIdx]       + offset;
        StPriceL   = tree_data->EDevStPrice[EDevIdx] + offset;
        AuxSliceL  = AuxSlice + offset;

        /* Assign level */
        for (i = Bottom1; i <= Top1; i++)
        {
            LevelL[i] = i - Bottom1;
            if (LevelL[i] < 0 || LevelL[i] > Area)
            {
                DR_Error ("Exceeded level array bounds\n");
                goto RETURN;
            }

            LMinInit = MIN(0       , LevelL[i]);
            LMaxInit = MAX(Area    , LevelL[i]);
        }

        YCumP_AA = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        YCumP_OU = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        YProb_AA = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        YStrk    = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);
        YProb_OU = (double *) DR_Array (DOUBLE, LMinInit-1, LMaxInit);

        if ((YCumP_AA == NULL) ||
            (YCumP_OU == NULL) ||
            (YProb_AA == NULL) ||
            (YStrk    == NULL) ||
            (YProb_OU == NULL)) goto RETURN;

        for (l = LMinInit; l <= LMaxInit; l++) YProb_AA[l] = 0;

        /* Binaries */
        for (i = Bottom1; i <= Top1; i++)
        {
            YProb_AA[LevelL[i]] += StPriceL[i] * AAtoOUFact * DCCFrac * FirstZeroL[i];
            YProb_OU[LevelL[i]] += StPriceL[i];
        }

        /* Cumulative */
        YCumP_AA[LMinInit]  = YProb_AA[LMinInit];
        YCumP_OU[LMinInit]  = YProb_OU[LMinInit];

        for (l = LMinInit+1; l <= LMaxInit; l++) 
        {
            YCumP_AA[l]  = YCumP_AA[l-1] + YProb_AA[l];
            YCumP_OU[l]  = YCumP_OU[l-1] + YProb_OU[l];
        }
        
        if (fabs(YCumP_AA[LMaxInit]-1.)>1e-6)
        {
            DR_Error ("Nmr_t: state prices don't add up to 1 at date %ld!",
                      mvd->NmrDate[NmrIdx]);
            goto RETURN;
        }

        /* search for cutoff levels based on cum probability */
        LMin = LMinInit;
        
        for (l=LMinInit; l<LMaxInit; l++)
        {
            if (YCumP_OU[l]  < NormalH (-tree_data->NbSigmaMax)) LMin = l+1;
            else break;
        }

        LMax = LMaxInit;
        
        for (l=LMaxInit; l>LMinInit; l--)
        {
            if (YCumP_OU[l] > NormalH (tree_data->NbSigmaMax)) LMax = l-1;
            else break;
        }

        /* Compute binary level */
        for (l = LMin; l <= LMax; l++)
        {
            if (Q3TMXMQMap (mq,
                         mq->muMQ + mq->sigMQ * Normal_InvH (YCumP_AA[l]),
                         &YStrk[l]) == FAILURE)
            { 
                 DR_Error ("Could not compute MqMap\n");
                goto RETURN;
            }
        }

        /* Yield on the given level */
        prevPrice = mq->fwdRate;

        for (l = LMin; l <= LMax-1; l++)
        {
            /* Call at new YStrk point */
            if (Q3TMXMQPricer (mq,
                            Q3TMX_CALL,
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
                LiborL[l] = LiborL[l-1]; 
            }
            else
            {
                LiborL[l] = (prevPrice - currPrice) / binProb;
            }

            /* Store stuff for reuse at the next level */
            prevPrice = currPrice;

        } /* for l */

        /* Libor on top level */
        currPrice    = 0.0;
        LiborL[LMax] = (prevPrice - currPrice) / (YCumP_AA[LMaxInit] - YCumP_AA[LMax-1]);


        /* Flat left Yield extrapolation  */
        for (l=LMinInit; l<LMin; l++) LiborL[l] = LiborL[LMin];
        for (l=LMaxInit; l>LMax; l--) LiborL[l] = LiborL[LMax];

        /* Compute AA tree price stats */
        for (i = Bottom1; i <= Top1; i++)
        {
            l            = LevelL[i];
            AuxSliceL[i] = DCCFrac * FirstZeroL[i] * (LiborL[l] - ParLibor0);
        }        

        for (i = Bottom1; i <= Top1; i++)
        {
            l              = LevelL[i];
            AASta[AA_FWD] += YProb_AA[l] * LiborL[l];

            if (SmoothingOn)            
            {            
                Step = Fix3_GetIndexStep(AuxSlice,1,i,0,0,t,tree_data);                
            }

            AASta[AA_ATM] += YProb_OU[l] * DrSmoothMax(AuxSliceL[i], 0.0, Step);
        }        

        AASta[AA_ATM] *= AAtoOUFact;

        /* Calibrate ATM price */
        AASlope = mvd->AtmLiborPrice0[NmrIdx] / AASta[AA_ATM]; 
        
        for (i = Bottom1; i <= Top1; i++)
        {
            l         = LevelL[i];
            LiborL[l] = ParLibor0 + AASlope * (LiborL[l] - ParLibor0);
        }    

  
        /* Compute inverse numeraire slice */
        for (i = Bottom1; i <= Top1; i++)
        {
            l          = LevelL[i];
            NmrInvL[i] = FirstZeroL[i] * (1.0 + DCCFrac * LiborL[LevelL[i]]);

            /* Store prob of negative numeraire in mvd BEFORE adjusting */
            if (NmrInvL[i] < NMR_CUTOFF)
            {
                DR_Error ("Negative numeraire at date %ld",
                          mvd->NmrDate[NmrIdx]);
                goto RETURN; 
            }
        } 

    } /* 1-factor */
    else if (tree_data->NbFactor == 2)
    {
        goto RETURN;
    }
    else if (tree_data->NbFactor == 3)
    {
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    Fix3_Free_Slice (AuxSlice, tree_data);
    Fix3_Free_Slice_Int (Level, tree_data);

    Free_DR_Array (YProb_AA, DOUBLE, LMinInit-1, LMaxInit);
    Free_DR_Array (YCumP_OU, DOUBLE, LMinInit-1, LMaxInit);
    Free_DR_Array (YCumP_AA, DOUBLE, LMinInit-1, LMaxInit);
    Free_DR_Array (YStrk,    DOUBLE, LMinInit-1, LMaxInit);
    Free_DR_Array (YProb_OU, DOUBLE, LMinInit-1, LMaxInit);

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_t: Failed!");
    }

    return (status);

} /* Nmr_t */



/*****  Nmr_Interp ***********************************************************/
/*
 * Compute interpolated numeraire in-between numeraire dates Ti-1 <= t <=Ti,
 * i being the current Nmr Idx:
 */
int     Nmr_Interp (double         *NmrInv,     /* (O) Nmr pointer            */
                    int            NmrIdx,      /* (I) Numeraire index        */
                    int            t,           /* (I) Current time point     */
                    MKTVOL_DATA    *mktvol_data,/* (I) Volatility data struct */
                    FIX3_TREE_DATA *tree_data)  /* (I) Tree data structure    */
{
    /* Local pointers */
    double  *LastZeroL;  
    double  *NmrInvL;
    double  *StPriceL;

    /* Used slices */
    double  *InterpZero = NULL;          /* Interp Zero to Nxt Crit Date    */
    double  *NxtOU      = NULL;
    double  *PrevOU     = NULL;
    double  *NxtZero    = NULL;
    double  *PrevZero   = NULL;    
    
    /* Interpolation variables */
    int     EDevIdx;                     /* Express DEV date index          */
    double  NmrInv0, NmrInvAdj;          /* Exp & adj of inverse nmr        */  
    long    NxtCritDate;
    long    NextT, PrevT;
    double  xt;
    double  Taut;
    double  PrevZeroDcf, NxtZeroDcf, InterpZeroDcf;
    int     CvDiff = tree_data->CvDiff;

    /* Tree limits */
    int     Top1, Bottom1;               /* Tree limits (1rst dim)          */
    int     *Top2, *Bottom2;             /* Tree limits (2nd dim)           */
    int     **Top3, **Bottom3;           /* Tree limits (3rd dim)           */
    int     NxtTop1, NxtBottom1;
    int     PrevTop1, PrevBottom1;
    
    int     i, l;                        /* Node indices                    */
    int     offset;                      /* Node offset                     */
    long    Area = 0;                    /* Slices size                     */
    int     status = FAILURE;            /* Error status                    */
    long    CurrentDate = tree_data->TPDate[t];

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
    
     /* Allocate numeraire and libor slices */
    InterpZero = (double *) DR_Array (DOUBLE, 0, Area);
    NxtOU      = (double *) DR_Array (DOUBLE, 0, Area);
    PrevOU     = (double *) DR_Array (DOUBLE, 0, Area);
    PrevZero   = (double *) DR_Array (DOUBLE, 0, Area);
    NxtZero    = (double *) DR_Array (DOUBLE, 0, Area);

    if ((InterpZero == NULL) ||
        (NxtOU      == NULL) ||
        (PrevOU     == NULL) ||
        (PrevZero   == NULL) ||
        (NxtZero    == NULL) )
    {
        DR_Error("Nmr_Interp: Could not allocate memory.\n");
        goto RETURN;
    }

    /* Initialisations        */    
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    /* Check numeraire index consistency */
    if (NmrIdx <1 || NmrIdx >= mktvol_data->NbNmr - 1)
    {
        DR_Error("Nmr_Interp: Can not perform numeraire interpolation.\n");
        goto RETURN;
    }

    /* Next Critical date in the tree */
    NxtCritDate = tree_data->NxtCritDate;
        
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

    /* Compute fraction within interval between numeraires */
    Taut  = (double) Daysact (CurrentDate, mktvol_data->NmrDate[NmrIdx]) 
          / (double) Daysact (mktvol_data->NmrDate[NmrIdx-1], mktvol_data->NmrDate[NmrIdx]);

    if (tree_data->NbFactor == 1)
    {
        /* Assign local pointers */
        offset    = Fix3_Node_Offset(1, 0, 0, t, tree_data);
        StPriceL  = tree_data->EDevStPrice[EDevIdx] + offset;
        LastZeroL = tree_data->LastZero + offset;
        NmrInvL   = NmrInv + offset;        
        
        /*
         * NEXT NUMERAIRE DATE INFO PROCESSING *
         */

        /* Find next numeraire date time */
        NextT = GetDLOffset (tree_data->NbTP, 
                             tree_data->TPDate,
                             mktvol_data->NmrDate[NmrIdx],
                             CbkEXACT);
        
        if (tree_data->TPDate[NextT] != mktvol_data->NmrDate[NmrIdx])
        {
            DR_Error ("State prices not available for next Nmr date %ld\n", 
                tree_data->TPDate[t]);
            goto RETURN;
        }
        
        /* Generate OU process at Ti (next Nmr Date) for interpolation      */
        NxtBottom1 = tree_data->Bottom1[NextT];
        NxtTop1    = tree_data->Top1[NextT];   
        
        for (i = NxtBottom1; i <= NxtTop1; i++)
        {
            l        = i - NxtBottom1;
            NxtOU[l] = i * tree_data->Aweight[0][NextT-1] * 
                       sqrt(JUMPCOEFF * tree_data->LengthJ[NextT-1]);
        }

       /* 
        * Get next Zero
        */
        if (DrDayCountFraction(CurrentDate, 
                               NxtCritDate, 
                               mktvol_data->DCC, 
                               &InterpZeroDcf) == FAILURE)
        {
            DR_Error("Nmr_Interp: Failed to compute cash accrual!");
            goto RETURN;            
        }

        if (DrDayCountFraction(mktvol_data->NmrDate[NmrIdx], 
                               mktvol_data->NmrDate[NmrIdx+1], 
                               mktvol_data->DCC, 
                               &NxtZeroDcf) == FAILURE)
        {
            DR_Error("Nmr_Interp: Failed to compute cash accrual!");
            goto RETURN;
        }

        for (l=0; l < Area; l++)
        {
            /* Deduce nxt zero fron mapped libor */
            NxtZero[l] = 1.0 / (1.0 + NxtZeroDcf * tree_data->Libor[NmrIdx][l]);

            /* Flat forward assumption */
            NxtZero[l] = pow(NxtZero[l], InterpZeroDcf/NxtZeroDcf);
        }

        /*
         * PREVIOUS NUMERAIRE DATE INFIO PROCESSING *
         */     
        
        /* Find previous numeraire date time */
        PrevT = GetDLOffset (tree_data->NbTP, 
                             tree_data->TPDate,
                             mktvol_data->NmrDate[NmrIdx-1],
                             CbkEXACT);
        
        if (tree_data->TPDate[PrevT] != mktvol_data->NmrDate[NmrIdx-1])
        {
            DR_Error ("State prices not available for next Nmr date %ld\n", 
                tree_data->TPDate[t]);
            goto RETURN;
        }
        
        /* Generate OU process at Ti (previous Nmr Date) for interpolation      */
        PrevBottom1 = tree_data->Bottom1[PrevT];
        PrevTop1    = tree_data->Top1[PrevT];   
        
        for (i = PrevBottom1; i <= PrevTop1; i++)
        {
            l        = i - PrevBottom1;
            PrevOU[l] = i * tree_data->Aweight[0][PrevT-1] * 
                sqrt(JUMPCOEFF * tree_data->LengthJ[PrevT-1]);
        }

         /* 
        * Get previous Zero
        */
        if (DrDayCountFraction(mktvol_data->NmrDate[NmrIdx-1], 
                               mktvol_data->NmrDate[NmrIdx], 
                               mktvol_data->DCC, 
                               &PrevZeroDcf) == FAILURE)
        {
            DR_Error("Nmr_Interp: Failed to compute cash accrual!");
            goto RETURN;
        }

        for (l=0; l < Area; l++)
        {
            /* Deduce previous zero from mapped libor */
            PrevZero[l] = 1.0 / (1.0 + PrevZeroDcf * tree_data->Libor[NmrIdx-1][l]);

            /* Flat forward assumption */
            PrevZero[l] = pow(PrevZero[l], InterpZeroDcf/PrevZeroDcf);
        }

        /* Interpolate Zero in X space */       
        for (i = Bottom1; i <= Top1; i++)
        {
            double PrevInterpZero;
            double NxtInterpZero;
            
            
            l = i - Bottom1;
            
            /* compute xt for this state   */
            xt = i * tree_data->Aweight[0][t-1] * 
                sqrt(JUMPCOEFF * tree_data->LengthJ[t-1]);
            
            /* interpolate Zero            */
            if (NmrIdx == 1)
            {
                tableinterp (xt,
                            &(NxtInterpZero),
                            NxtOU,
                            NxtZero,                   
                            NxtTop1 - NxtBottom1 + 1);
                
                InterpZero[l] = NxtInterpZero;
            }
            else
            {
                /* Interpolation from previous nmr zero */
                tableinterp (xt,
                            &(PrevInterpZero),
                            PrevOU,
                            PrevZero,
                            PrevTop1 - PrevBottom1 + 1);


                /* Interpolation from nxt nmr zero */
                tableinterp (xt,
                             &(NxtInterpZero),
                             NxtOU,
                             NxtZero,
                             NxtTop1 - NxtBottom1 + 1);
                
                InterpZero[l] = Taut * PrevInterpZero + (1.0 - Taut) * NxtInterpZero;               
            }

        } /* for i >= Bottom1 */        
                
        /* Create inverse numeraire slice */
        NmrInv0 = 0;
        
        for (i = Bottom1; i <= Top1; i++)
        {   
            l          = i - Bottom1;
            NmrInvL[i] = LastZeroL[i] / InterpZero[l];

            /* Ensure positivity */
            NmrInvL[i] = MAX(NmrInvL[i], NMR_CUTOFF);
            
            /* update expectation */
            NmrInv0   += NmrInvL[i] * StPriceL[i];
        }
        
        NmrInvAdj = 1/ (tree_data->TermZero[CvDiff][t] * NmrInv0);

        /* Check adjusment is small */
        if (fabs(NmrInvAdj - 1.0) > 0.05)
        {
            DR_Error("Nmr_Interp: adjusment to numeraire is out of range!");
            goto RETURN;
        }
        
        
        /* Adjust numeraire */
        for (i = Bottom1; i <= Top1; i++)
        {
            NmrInvL[i] *= NmrInvAdj;
        }

    }
    else if (tree_data->NbFactor == 2)
    {
        ;
    }
    else if (tree_data->NbFactor == 3)
    {
        ;
    }

    status = SUCCESS;

RETURN:

    Free_DR_Array  (InterpZero, DOUBLE, 0, Area);
    Free_DR_Array  (NxtOU,      DOUBLE, 0, Area);
    Free_DR_Array  (PrevOU,     DOUBLE, 0, Area);
    Free_DR_Array  (PrevZero,   DOUBLE, 0, Area);
    Free_DR_Array  (NxtZero,    DOUBLE, 0, Area);

    if (status != SUCCESS)
    {        
        DR_Error("Nmr_Interp: Failed!");
    }

    return (status);

} /* Nmr_Interp */


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
    int     SwapIdx = 0;

    if (mvd->TraceFlag == 'N') return SUCCESS;

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
        if (mvd->NmrDate[n] == mvd->SwapSt[SwapIdx])
        {
            
            fprintf (stream, "  [%3d] %8ld   %8ld  %8ld   \n",
                                n,
                                mvd->NmrDate[n],
                                mvd->SwapSt[SwapIdx],
                                mvd->SwapMat[SwapIdx]);

            if (SwapIdx < mvd->NbVol - 1)
                SwapIdx ++;
        }
        else
        {
            fprintf (stream, "  [%3d] %8ld     \n",
                                n,
                                mvd->NmrDate[n]);
        }

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

    if (mvd->TraceFlag == 'N') return SUCCESS;

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

    for (n=1; n<mvd->NbNmr-1; n++)
    {
        fprintf (stream, "  [%3d] %8ld   %8ld  %8ld   %8.4f   %8.4f  "
                 "%4.2f  %4.2f  %4.2f  %4.2f  %4.2f  %5.2f  %4.2f  %5.2f\n",
                 n,
                 mvd->NmrDate[n],
                 mvd->NmrDate[n],
                 mvd->NmrDate[n+1],
                 mvd->ParLibor0[n]*100,
                 mvd->NmrLibVol[0][n]*100,
                 mvd->NmrLibVol[1][n],
                 mvd->NmrLibVol[2][n],
                 mvd->NmrLibVol[3][n],
                 mvd->NmrLibVol[4][n],
                 mvd->NmrLibVol[5][n],
                 mvd->NmrLibVol[6][n],
                 mvd->NmrLibVol[7][n],
                 mvd->NmrLibVol[8][n]
                );
    }

    /* Print vanila output info */
    fprintf (stream, "###### NMR0 VANILLA MQ CALIB INFO ######\n"); 
    fprintf (stream, "# [idx] expiry__  fwd_____  sigATM__  muMQ____  sigMQ___"
             "  qL[2]___  qL[1]___  qL[0]___  qR[0]___  qR[1]___  qR[2]___"
             "  fwdC____  volC____  \n");

    for (n=1; n<mvd->NbNmr-1; n++)
    {
        fprintf (stream, "  [%3d] %8.2f  %8.4f  %8.4f  %8.4f  %8.4f  "
                 "%8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  \n",
                 n,
                 mq[n].expiry,
                 100*mq[n].fwdRate,
                 100*mq[n].sigATM,
                 mq[n].muMQ,
                 mq[n].sigMQ,
                 mq[n].qL[2],
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








    
    








    
    








