#include <fix123head.h>    
#include <dll_util.h>    
#include <math.h>    

/* Pack Market Volatility Data */  
/* Task : Packs up MKTVOL_DATA and FIX3_TREE_DATA structures.    
*/    
int Fix3_PackMktVolAndTreeData(  
        MKTVOL_DATA*        VolData,         /* (O) Vol data structure          */
        FIX3_TREE_DATA*     Tree,            /* (O) Tree data structure         */
        long const*         BaseDatesL,      /* (I) Base dates                  */
        long const*         VolDatesL,       /* (I) Vol dates                   */
        long const*         VolMatsL,        /* (I) Vol underlying maturity date*/
        double const*       VolsL,           /* (I) Vols                        */
        char **             VolFreqL,        /* (I) Vol underlying frequency    */
        char **             VolDCCL,         /* (I) Vol underlying DCC          */
        char **             VolTypeL,        /* (I) Vol type ('N'orm, 'L'og)    */
        char **             VolCalibL,       /* (I) Vol calib flags             */
        double const*       MrParamsL,       /* (I) Model MR parameters         */
        double const*       SmileParamsL,    /* (I) Model smile parameters      */
        long const*         TreeParamsL)     /* (I) Tree parameters             */
{    
static  char    routine[] = "PackMktVolAndTreeData";

   int  i;
   int  NumFact;
   int  MrParamsSize;
   double norm;

   ASSERT_OR_FAIL ((int)BaseDatesL[0] == (L_BASEDATE_SIZE-1));

   /* Check inputs */
   ASSERT_OR_FAIL ((int)MrParamsL[0]    >  3);
   ASSERT_OR_FAIL ((int)SmileParamsL[0] == (L_SMILEPARAM_SIZE-1));
   ASSERT_OR_FAIL ((int)TreeParamsL[0]  == (L_TREEPARAM_SIZE-1));
   ASSERT_OR_FAIL ((int)VolDatesL[0]    == (int)VolMatsL[0]);
   ASSERT_OR_FAIL ((int)VolDatesL[0]    == (int)VolsL[0]);
   ASSERT_OR_FAIL ((int)VolCalibL[0]    == (L_CALIBFLAG_SIZE-1));

   NumFact  = (int)MrParamsL[1];

   if (NumFact < 0 ||
       NumFact > 3)
   {
        DR_Error("%s: invalid number of factors (%d).\n",
                 routine, NumFact);
        return (FAILURE);
   }

   VolData->Alpha[1] = -999.;
   VolData->Alpha[2] = -999.;
   VolData->Beta[1]  = -999.;
   VolData->Beta[2]  = -999.;
   VolData->Rho[0]   = -999.;
   VolData->Rho[1]   = -999.;
   VolData->Rho[2]   = -999.;


   MrParamsSize = 2 + 2*NumFact + NumFact*(NumFact-1)/2;
   if ((int)MrParamsL[0] != MrParamsSize) 
   {
        DR_Error("%s: expects mr parameter array length %d+2 with "
                 "%d factors (got total %d).\n", 
                 routine, MrParamsSize-2, NumFact, (int)MrParamsL[0]);
        return (FAILURE);
   }
   
   VolData->BaseDate = DRDate(BaseDatesL[L_VolBaseDate]);   

   VolData->NbVol    = (int)VolDatesL[0] ;    

   COPY_FROM_DATES_L (VolDatesL, MAXNBDATE, VolData->VolDate);
   COPY_FROM_DATES_L (VolDatesL, MAXNBDATE, VolData->SwapSt);
   COPY_FROM_DATES_L (VolMatsL,  MAXNBDATE, VolData->SwapMat);
   COPY_FROM_ARRAY_L (VolsL,     MAXNBDATE, VolData->Vol);

   for ( i = 0; i < (int)VolDatesL[0]; i++ )
   {
       VolData->VolUsed[i] = TRUE;
   }

    /* Initialize model-dependent interface */
    if (Fix3_Model_Interface_Init (VolData->ModelChoice) == FAILURE) return (FAILURE);

   VolData->Freq   = (char)toupper( VolFreqL[1][0] );    
   VolData->DCC    = (char)toupper( VolDCCL[1][0] );  
  
   /* Vol unit */
   VolData->VolUnit =  
      ( toupper( VolTypeL[1][0] ) == 'L' ) ? 1 : 0;

   /* Smile */
   VolData->QLeft     = 1.0 - SmileParamsL[L_QLeft]; 
   VolData->QRight    = 1.0 - SmileParamsL[L_QRight];
   VolData->FwdShift  = SmileParamsL[L_FwdShift];
   VolData->CetNbIter = (int)SmileParamsL[L_CetNbIter];


   for (i=0; i<NumFact; i++) 
   {
        VolData->Beta[i]  = MrParamsL[2+i];
        VolData->Alpha[i] = MrParamsL[2+NumFact+i];
   }
   for (i=0; i<NumFact*(NumFact-1)/2; i++)
   {
        VolData->Rho[i] = MrParamsL[2+2*NumFact+i];
   }

   VolData->Bbq       = 1.0 - MrParamsL[MrParamsSize];

   norm = 0.;
   for (i = 0; i < NumFact; i++) 
        norm += VolData->Alpha[i] * VolData->Alpha[i]; 
   norm = sqrt(norm);
   if (fabs(norm) < ERROR)
   {      
       DR_Error("Total alpha is too small !");
       return (FAILURE);
   }       
   if (IS_EQUAL(VolData->Bbq,1))
   {
       VolData->VolNorm = 0.;
       VolData->VolLogn = norm;
   }
   else
   if (IS_EQUAL(VolData->Bbq,0))
   {
       VolData->VolNorm = norm;
       VolData->VolLogn = 0.;
   }
   else
   {
       DR_Error("Bbq parameter must be either 0 or 1 !");
       return (FAILURE);
   }       

   VolData->SkipFlag =   
      ( toupper( VolCalibL[L_VolSkipFlag][0] ) == 'Y' ) ? TRUE : FALSE;
   VolData->CalibFlag =  
      ( toupper( VolCalibL[L_VolCalibFlag][0] ) == 'Y' ) ? TRUE : FALSE;
   VolData->FilterSpotVolFlag =  
      ( toupper( VolCalibL[L_VolFilterSpotVolFlag][0] ) == 'Y' ) ? TRUE : FALSE;
   VolData->SmoothingFlag =  toupper( VolCalibL[L_VolSmoothingFlag][0] );


   Tree->NbFactor    = NumFact;
   Tree->Ppy         = (int)TreeParamsL[L_Ppy];
   Tree->NbSigmaMax  = (int)TreeParamsL[L_NbSigmaMax];
   Tree->CvDiff      = (int)TreeParamsL[L_CvDiff];
   Tree->CvDisc      = (int)TreeParamsL[L_CvDisc];
  
    if (Tree->CvDiff == 0)
    {
        Tree->CvIdx1 = 1;
        Tree->CvIdx2 = 2;
    }
    else if (Tree->CvDiff == 1)
    {
        Tree->CvIdx1 = 0;
        Tree->CvIdx2 = 2;
    }
    else if (Tree->CvDiff == 2)
    {
        Tree->CvIdx1 = 0;
        Tree->CvIdx2 = 1;
    }
    else
    {
        DR_Error ("Diffuse curve must be 0,1,2");
        return (FAILURE);
    }

    return( SUCCESS );    

}    




/* Unpack Market Volatility Data */  
/* Task : Unpacks MKTVOL_DATA and FIX3_TREE_DATA structures.    
*/    
int Fix3_UnPackMktVolAndTreeData(  
        MKTVOL_DATA const*      VolData,/* (I) Vol data structure          */
        FIX3_TREE_DATA const*   Tree,   /* (I) Tree data structure         */
        long          *BaseDatesL,      /* (O) Base dates                  */
        long          *VolDatesL,       /* (O) Vol dates                   */
        long          *VolMatsL,        /* (O) Vol underlying maturity date*/
        double        *VolsL,           /* (O) Vols                        */
        char          **VolFreqL,       /* (O) Vol underlying frequency    */
        char          **VolDCCL,        /* (O) Vol underlying DCC          */
        char          **VolTypeL,       /* (O) Vol type ('N'orm, 'L'og)    */
        char          **VolCalibL,      /* (O) Vol calib flags             */
        double        *MrParamsL,       /* (O) Model MR parameters         */
        double        *SmileParamsL,    /* (O) Model smile parameters      */
        long          *TreeParamsL)     /* (O) Tree parameters             */
{    
   int  i;
   int  NumFact;
   int  MrParamsSize;


   /* Volatility */
   BaseDatesL[L_VolBaseDate] = YMDDateFromIRDate(VolData->BaseDate);


   COPY_TO_DATES_L (VolDatesL, VolData->NbVol, VolData->VolDate);
   COPY_TO_DATES_L (VolDatesL, VolData->NbVol, VolData->SwapSt);
   COPY_TO_DATES_L (VolMatsL,  VolData->NbVol, VolData->SwapMat);
   COPY_TO_ARRAY_L (VolsL,     VolData->NbVol, VolData->Vol);


   VolFreqL[0] = (char *) 1;
   sprintf(VolFreqL[1], "%c", VolData->Freq);
   VolDCCL[0]  = (char *) 1;
   sprintf(VolDCCL[1], "%c", VolData->DCC);


   /* Smile */
   SmileParamsL[0]           = (double)(L_SMILEPARAM_SIZE - 1);
   SmileParamsL[L_QLeft]     = 1.0 - VolData->QLeft;
   SmileParamsL[L_QRight]    = 1.0 - VolData->QRight;
   SmileParamsL[L_FwdShift]  = VolData->FwdShift;
   SmileParamsL[L_CetNbIter] = (double)VolData->CetNbIter;

   /* Model parameter */
   NumFact = Tree->NbFactor;
   MrParamsSize = 2 + 2*NumFact + NumFact*(NumFact-1)/2;

   MrParamsL[0] = (double)MrParamsSize;
   MrParamsL[1] = (double)NumFact;
   
   for (i=0; i<NumFact; i++) 
   {
        MrParamsL[2+i]         = VolData->Beta[i];
        MrParamsL[2+NumFact+i] = VolData->Alpha[i];
   }
   for (i=0; i<NumFact*(NumFact-1)/2; i++)
   {
        MrParamsL[2+2*NumFact+i] = VolData->Rho[i];
   }

   MrParamsL[MrParamsSize] = 1.0 - VolData->Bbq;

   /* Vol type- normal or lognormal */
   sprintf(VolTypeL[1],  "%c", (VolData->VolUnit) ?'L':'N'); 

   /* Vol calib */
   VolCalibL[0] = (char *) (L_CALIBFLAG_SIZE - 1);
   sprintf(VolCalibL[L_VolSkipFlag],         "%c", (VolData->SkipFlag) ?'Y':'N');   
   sprintf(VolCalibL[L_VolCalibFlag],        "%c", (VolData->CalibFlag) ?'Y':'N'); 
   sprintf(VolCalibL[L_VolFilterSpotVolFlag],"%c", (VolData->FilterSpotVolFlag)?'Y':'N'); 
   sprintf(VolCalibL[L_VolSmoothingFlag],    "%c", VolData->SmoothingFlag);

   TreeParamsL[0]            = (int)(L_TREEPARAM_SIZE -1);
   TreeParamsL[L_Ppy]        = Tree->Ppy;
   TreeParamsL[L_NbSigmaMax] = Tree->NbSigmaMax;
   TreeParamsL[L_CvDiff]     = Tree->CvDiff;
   TreeParamsL[L_CvDisc]     = Tree->CvDisc;
  
   return( SUCCESS );    
}    


/* Pack Market Volatility Data */  
/* Task : Packs up MKTVOL_DATA and FIX3_TREE_DATA structures.    
*/    
int Fix3_PackMktVolAndModelData(  
        MKTVOL_DATA        *VolData,         /* (O) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (O) Tree data structure         */
        long               *BaseDatesL,      /* (I) Base dates                  */
        long               *VolDatesL,       /* (I) Vol dates                   */
        long               *VolMatsL,        /* (I) Vol underlying maturity date*/
        double             *VolsL,           /* (I) Volatilities                */
        char               **VolCalibL,      /* (I) Vol calibration flags       */
        long               *MrDatesL,        /* (I) Mean reversion dates        */
        double             *MrParamsL,       /* (I) Model MR parameters         */
        long               *SmileDatesL,     /* (I) Smile benchmark dates       */
        double             *SmileParamsL,    /* (I) Model smile parameters      */
        long               *NumParamsL)      /* (I) Numerical parameters        */
{    

    static  char    routine[] = "Fix3_PackMktVolAndTreeData";
    char ModelChoiceString[MAXBUFF];
    int  i;

    long BaseDates[L_BASEDATE_SIZE-1];


    ASSERT_OR_FAIL ((int)BaseDatesL[0] == (L_BASEDATE_SIZE-1));
    COPY_FROM_DATES_L (BaseDatesL, L_BASEDATE_SIZE, BaseDates);

    ASSERT_OR_FAIL ((int)VolDatesL[0]  == (int)VolMatsL[0]);
    ASSERT_OR_FAIL ((int)VolDatesL[0]  == (int)VolsL[0]);
    ASSERT_OR_FAIL ((int)VolCalibL[0]  == (VOLCALIB_SIZE-1));

    VolData->NbVol    = (int)VolDatesL[0];
    VolData->BaseDate = BaseDates[L_VolBaseDate-1];   

    COPY_FROM_DATES_L (VolDatesL, MAXNBDATE, VolData->VolDate);
    COPY_FROM_DATES_L (VolDatesL, MAXNBDATE, VolData->SwapSt);
    COPY_FROM_DATES_L (VolMatsL,  MAXNBDATE, VolData->SwapMat);
    COPY_FROM_ARRAY_L (VolsL,     MAXNBDATE, VolData->Vol);

    for ( i = 0; i < (int)VolDatesL[0]; i++ ) VolData->VolUsed[i] = TRUE;

    strcpy (ModelChoiceString, VolCalibL[VOL_ModelChoice]);

    /* Identify model string */
    if      (!strcmp(ModelChoiceString,"ORIGINAL") || !strcmp(ModelChoiceString,"original"))   
            VolData->ModelChoice = FIX3_ORIGINAL;
    else if (!strcmp(ModelChoiceString,"CLASSIC")  || !strcmp(ModelChoiceString,"classic"))   
            VolData->ModelChoice = FIX3_CLASSIC;
    else if (!strcmp(ModelChoiceString,"TIMEDEP")  || !strcmp(ModelChoiceString,"timedep"))    
            VolData->ModelChoice = FIX3_TIMEDEP;
    else if (!strcmp(ModelChoiceString,"2Q")       || !strcmp(ModelChoiceString,"2q")) 
            VolData->ModelChoice = FIX3_TIMEDEP;
    else if (!strcmp(ModelChoiceString,"SMD")      || !strcmp(ModelChoiceString,"smd"))  
            VolData->ModelChoice = FIX3_SMD;
    else if (!strcmp(ModelChoiceString,"TMX")      || !strcmp(ModelChoiceString,"tmx"))    
            VolData->ModelChoice = FIX3_TMX;
    else if (!strcmp(ModelChoiceString,"E2Q")      || !strcmp(ModelChoiceString,"e2q"))    
            VolData->ModelChoice = FIX3_E2Q;
    else 
            VolData->ModelChoice = -999;

    VolData->IsNmrModel = (VolData->ModelChoice == FIX3_TMX);

    /* Initialize model-dependent interface */
    if (Fix3_Model_Interface_Init (VolData->ModelChoice) == FAILURE) return (FAILURE);

    /* Volatility calibration flags */
    VolData->SkipFlag          = (toupper (VolCalibL[VOL_SkipFlag][0])=='Y') ? TRUE : FALSE;
    VolData->CalibFlag         = (toupper (VolCalibL[VOL_CalibFlag][0])=='Y') ? TRUE : FALSE;
    VolData->SmoothingFlag     = toupper (VolCalibL[VOL_SmoothingFlag][0]);
    VolData->Freq              = toupper (VolCalibL[VOL_Freq][0]);    
    VolData->DCC               = toupper (VolCalibL[VOL_DCC][0]);  
    VolData->VolUnit           = (toupper (VolCalibL[VOL_Type][0]) == 'L' ? 1 : 0);
    VolData->Bbq               = (toupper (VolCalibL[VOL_Backbone][0]) == 'L' ? 0 : 1);


    /* Numerical parameters */
    Tree->Ppy          = (int) NumParamsL[NUM_Ppy];
    Tree->NbSigmaMax   = (int) NumParamsL[NUM_NbSigmaMax];
    VolData->CetNbIter = (int) NumParamsL[NUM_CetNbIter];
    Tree->CvDiff       = (int) NumParamsL[NUM_CvDiff];
    Tree->CvDisc       = (int) NumParamsL[NUM_CvDisc];

    /* Additional numerics for TMX */
    if (VolData->ModelChoice == FIX3_TMX)
    {
        VolData->NbSigmaMQ = (int) NumParamsL[NUM_NbSigmaMQ];
        VolData->NckMQ     = (int) NumParamsL[NUM_NckMQ];
    }

    /* FIXME: IRX curves do not support linear interpolation. */
    EslSetZeroInterpolation((int)NumParamsL[NUM_InterpType] ? ESL_INTERP_FLATFWD : ESL_INTERP_LINEAR);

    /* Assign "other" curves */
    /* Set indices for index curves according to diff curve */ 
    if (Tree->CvDiff == 0)
    {
        Tree->CvIdx1 = 1;
        Tree->CvIdx2 = 2;
    }
    else if (Tree->CvDiff == 1)
    {
        Tree->CvIdx1 = 0;
        Tree->CvIdx2 = 2;
    }
    else if (Tree->CvDiff == 2)
    {
        Tree->CvIdx1 = 0;
        Tree->CvIdx2 = 1;
    }
    else
    {
        DR_Error ("Diffuse curve must be 0,1,2");
        return (FAILURE);
    }

  
    /* Model parameters: Vnfm, smile, numerics */
    return (Fix3_PackModelData (VolData,
                                Tree,
                                MrDatesL,
                                MrParamsL,
                                SmileDatesL,
                                SmileParamsL));
}   


int Fix3_PackModelData_Classic (  
        MKTVOL_DATA        *VolData,         /* (O) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (O) Tree data structure         */
        long               *MrDatesL,        /* (I) Mean reversion dates        */
        double             *MrParamsL,       /* (I) Model MR parameters         */
        long               *SmileDatesL,     /* (I) Smile benchmark dates       */
        double             *SmileParamsL)    /* (I) Model smile parameters      */
{    
    int  i;
    int  NbFactor;
    int  MrParamsSize;
    double norm;
    static char routine[] = "Fix3_PackModelData_Classic";

    /* Check inputs */
    ASSERT_OR_FAIL ((int)MrParamsL[0] >= 2);
    ASSERT_OR_FAIL ((int)SmileParamsL[0] == (CLASSIC_SMILEPARAM_SIZE-1));

    NbFactor = (int) MrParamsL[1];
    if (NbFactor < 0 || NbFactor > 3)
    {
        DR_Error("%s: invalid number of factors (%d).\n", routine, NbFactor);
        return (FAILURE);
    }
    /* Copy into tree and mktvol data as well */
    Tree->NbFactor    = NbFactor;
    VolData->NbFactor = NbFactor;

    /* Check size */
    MrParamsSize = 2*NbFactor + NbFactor*(NbFactor-1)/2;
    if ((int)MrParamsL[0] != MrParamsSize+1) 
    {
        DR_Error("%s: expects mr parameter array length %d+1 with "
                 "%d factors (got total %d).\n", 
                 routine, MrParamsSize, NbFactor, (int)MrParamsL[0]);
        return (FAILURE);
    }

    /* Initialization -- needed in 1F and 2F case */
    VolData->Alpha[0] = -999;
    VolData->Alpha[1] = -999.;
    VolData->Alpha[2] = -999.;
    VolData->Beta[0]  = -999.;
    VolData->Beta[1]  = -999.;
    VolData->Beta[2]  = -999.;
    VolData->Rho[0]   = -999.;
    VolData->Rho[1]   = -999.;
    VolData->Rho[2]   = -999.;

    /* Initialize TD data -- not used */
    VolData->NbTDInp = 1;
    VolData->TDInpDate[0] = VolData->BaseDate;

    /* Read Vnfm parameters */
    for (i=0; i<NbFactor; i++) 
    {
         VolData->Beta[i]  = MrParamsL[2+i];
         VolData->Alpha[i] = MrParamsL[2+NbFactor+i];
    }
    for (i=0; i<NbFactor*(NbFactor-1)/2; i++)
    {
         VolData->Rho[i] = MrParamsL[2+2*NbFactor+i];
    }

    /* Assign volatility scale */
    norm = 0.;
    for (i = 0; i < NbFactor; i++) 
        norm += VolData->Alpha[i] * VolData->Alpha[i]; 
    norm = sqrt(norm);
    if (fabs(norm) < ERROR)
    {      
        DR_Error("Total alpha is too small !");
        return (FAILURE);
    }       
    if (IS_EQUAL(VolData->Bbq,1))
    {
        VolData->VolNorm = 0.;
        VolData->VolLogn = norm;
    }
    else if (IS_EQUAL(VolData->Bbq,0))
    {
        VolData->VolNorm = norm;
        VolData->VolLogn = 0.;
    }
    else
    {
        DR_Error("Bbq parameter must be either 0 or 1 !");
        return (FAILURE);
    }       

    /* Smile */
    VolData->QLeft     = 1.0 - SmileParamsL[CLS_QLeft]; 
    VolData->QRight    = 1.0 - SmileParamsL[CLS_QRight];
    VolData->FwdShift  = SmileParamsL[CLS_FwdShift];


    return (SUCCESS);
} /* Fix3_PackModelData_Classic */



int Fix3_PackModelData_TimeDep (  
        MKTVOL_DATA        *VolData,         /* (O) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (O) Tree data structure         */
        long               *MrDatesL,        /* (I) Mean reversion dates        */
        double             *MrParamsL,       /* (I) Model MR parameters         */
        long               *SmileDatesL,     /* (I) Smile benchmark dates       */
        double             *SmileParamsL)    /* (I) Model smile parameters      */
{    
    int  i,j;
    int  idx = 0;  /* skip entry #0 in counted array */
    int  NbFactor;
    int  MrParamsSize;
    double norm;
    static char routine[] = "Fix3_PackModelData_TimeDep";

    /* Check inputs */
    ASSERT_OR_FAIL ((int)MrParamsL[0] >= 2);
 

    NbFactor = (int) MrParamsL[1];
    if (NbFactor < 0 || NbFactor > 3)
    {
        DR_Error("%s: invalid number of factors (%d).\n", routine, NbFactor);
        return (FAILURE);
    }
    /* Copy into tree and mktvol data as well */
    Tree->NbFactor    = NbFactor;
    VolData->NbFactor = NbFactor;

    /* Check size */
    MrParamsSize = 2*NbFactor + NbFactor*(NbFactor-1)/2;
    if ((int)MrParamsL[0] != MrParamsSize * (int)MrDatesL[0] + 1)
    {
        DR_Error("%s: expects mr parameter array length %d for "
                 "%d factors and %d dates (got total %d).\n", 
                 routine, MrParamsSize * (int)MrDatesL[0] + 1,
                 NbFactor, (int)MrDatesL[0],(int)MrParamsL[0]);
        return (FAILURE);
    }

    /* Read smile dates */
    VolData->NbSmileDates = (int) SmileDatesL[0];
    COPY_FROM_DATES_L (SmileDatesL, MAXNBDATE, VolData->SmileDate);

    if ((int)SmileParamsL[0] != 3 * (int)SmileDatesL[0])
    {
        DR_Error("%s: expects mr parameter array length %d "
                 "%d dates (got total %d).\n", 
                 routine, 3 * (int)SmileDatesL[0] ,
                 (int)SmileDatesL[0],(int)SmileParamsL[0]);
        return (FAILURE);
    }
   


    /* Get time-dependent smile parameters*/
    for (i = 0; i < VolData->NbSmileDates; i++)
    {
        VolData->QLeftTD[i] = 1.0 - SmileParamsL[++idx]; 
    }
    for (i = 0; i < VolData->NbSmileDates; i++)
    {
        VolData->QRightTD[i] = 1.0 - SmileParamsL[++idx]; 
    }
    for (i = 0; i < VolData->NbSmileDates; i++)
    {
        VolData->FwdShiftTD[i] = SmileParamsL[++idx]; 
    }


    /* Read Vnfm dates */
    VolData->NbTDInp = (int) MrDatesL[0];
    COPY_FROM_DATES_L (MrDatesL, MAXNBDATE, VolData->TDInpDate);

    idx = 1;
    /* Get time-dependent mean reversions */
    for (i = 0; i < NbFactor; i++)
    {
        for (j = 0; j < VolData->NbTDInp; j++)
        {
            VolData->BetaTD[i][j] = MrParamsL[++idx];
        }
    }
    for (i = 0; i < NbFactor; i++)
    {
        for (j = 0; j < VolData->NbTDInp; j++)
        {
            VolData->AlphaTD[i][j] = MrParamsL[++idx];
        }
    }
    for (i = 0; i < NbFactor*(NbFactor-1)/2 ; i++)
    {
        for (j = 0; j < VolData->NbTDInp; j++)
        {
            VolData->RhoTD[i][j] = MrParamsL[++idx];
        }
    }

    /* Initialization */
    VolData->Alpha[0] = -999;
    VolData->Alpha[1] = -999.;
    VolData->Alpha[2] = -999.;
    VolData->Beta[0]  = -999.;
    VolData->Beta[1]  = -999.;
    VolData->Beta[2]  = -999.;
    VolData->Rho[0]   = -999.;
    VolData->Rho[1]   = -999.;
    VolData->Rho[2]   = -999.;

    /* Assign volatility scale */
    norm = 1.;
    if (IS_EQUAL(VolData->Bbq,1))
    {
        VolData->VolNorm = 0.;
        VolData->VolLogn = norm;
    }
    else if (IS_EQUAL(VolData->Bbq,0))
    {
        VolData->VolNorm = norm;
        VolData->VolLogn = 0.;
    }
    else
    {
        DR_Error("Bbq parameter must be either 0 or 1 !");
        return (FAILURE);
    }      

    


    return (SUCCESS);
} /* Fix3_PackModelData_TimeDep */


int Fix3_PackModelData_Smd (  
        MKTVOL_DATA        *VolData,         /* (O) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (O) Tree data structure         */
        long               *MrDatesL,        /* (I) Mean reversion dates        */
        double             *MrParamsL,       /* (I) Model MR parameters         */
        long               *SmileDatesL,     /* (I) Smile benchmark dates       */
        double             *SmileParamsL)    /* (I) Model smile parameters      */
{    
    int  i;
    int  NbFactor;
    int  MrParamsSize;
    double norm;
    static char routine[] = "Fix3_PackModelData_Smd";

    /* Check inputs */
    ASSERT_OR_FAIL ((int)MrParamsL[0] >= 2);
    ASSERT_OR_FAIL ((int)SmileParamsL[0] == (SMD_SMILEPARAM_SIZE-1));

    NbFactor = (int) MrParamsL[1];
    if (NbFactor < 0 || NbFactor > 3)
    {
        DR_Error("%s: invalid number of factors (%d).\n", routine, NbFactor);
        return (FAILURE);
    }
    /* Copy into tree and mktvol data as well */
    Tree->NbFactor    = NbFactor;
    VolData->NbFactor = NbFactor;

    /* Check size */
    MrParamsSize = 2*NbFactor + NbFactor*(NbFactor-1)/2;
    if ((int)MrParamsL[0] != MrParamsSize+1) 
    {
        DR_Error("%s: expects mr parameter array length %d+1 with "
                 "%d factors (got total %d).\n", 
                 routine, MrParamsSize, NbFactor, (int)MrParamsL[0]);
        return (FAILURE);
    }

    /* Initialization -- needed in 1F and 2F case */
    VolData->Alpha[0] = -999;
    VolData->Alpha[1] = -999.;
    VolData->Alpha[2] = -999.;
    VolData->Beta[0]  = -999.;
    VolData->Beta[1]  = -999.;
    VolData->Beta[2]  = -999.;
    VolData->Rho[0]   = -999.;
    VolData->Rho[1]   = -999.;
    VolData->Rho[2]   = -999.;


    /* Read Vnfm parameters */
    VolData->NbTDInp = 1;
    for (i=0; i<NbFactor; i++) 
    {
         VolData->Beta[i]  = MrParamsL[2+i];
         VolData->Alpha[i] = MrParamsL[2+NbFactor+i];
    }
    for (i=0; i<NbFactor*(NbFactor-1)/2; i++)
    {
         VolData->Rho[i] = MrParamsL[2+2*NbFactor+i];
    }

    /* Assign volatility scale */
    norm = 1.; 
    if (IS_EQUAL(VolData->Bbq,1))
    {
        VolData->VolNorm = 0.;
        VolData->VolLogn = norm;
    }
    else if (IS_EQUAL(VolData->Bbq,0))
    {
        VolData->VolNorm = norm;
        VolData->VolLogn = 0.;
    }
    else
    {
        DR_Error("Bbq parameter must be either 0 or 1 !");
        return (FAILURE);
    }      

    /* Smile */
    VolData->QLeft     = 1.0 - SmileParamsL[SMD_QLeft]; 
    VolData->QRight    = 1.0 - SmileParamsL[SMD_QRight];
    VolData->FwdShift  = SmileParamsL[SMD_FwdShift];
    VolData->Afac      = SmileParamsL[SMD_Afac];
    VolData->Bfac      = SmileParamsL[SMD_Bfac];
    VolData->Cfac      = SmileParamsL[SMD_Cfac];
    VolData->Dfac      = SmileParamsL[SMD_Dfac];


    return (SUCCESS);
} /* Fix3_PackModelData_Smd */



int Fix3_PackModelData_Tmx (  
        MKTVOL_DATA        *VolData,         /* (O) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (O) Tree data structure         */
        long               *MrDatesL,        /* (I) Mean reversion dates        */
        double             *MrParamsL,       /* (I) Model MR parameters         */
        long               *SmileDatesL,     /* (I) Smile benchmark dates       */
        double             *SmileParamsL)    /* (I) Model smile parameters      */
{
    static char routine[] = "Fix3_PackModelData_Tmx";

    int  i, j, l, idx;
    int  NbFactor;
    int  MrParamsSize;
    double norm;

    /* temporary storage for vol / smile */
    double TmpVols[MAXNBDATE * (NBVOLPARS-1)];

    /* Check inputs */
    ASSERT_OR_FAIL ((int)MrParamsL[0] >= 2);
    ASSERT_OR_FAIL ((int)SmileDatesL[0] <= VolData->NbVol);
    ASSERT_OR_FAIL ((int)SmileParamsL[0] == (int)SmileDatesL[0] * (NBVOLPARS-1));

    NbFactor = (int) MrParamsL[1];
    if (NbFactor < 0 || NbFactor > 3)
    {
        DR_Error("%s: invalid number of factors (%d).\n", routine, NbFactor);
        return (FAILURE);
    }
    /* Copy into tree and mktvol data as well */
    Tree->NbFactor    = NbFactor;
    VolData->NbFactor = NbFactor;

    /* Check size */
    MrParamsSize = 2*NbFactor + NbFactor*(NbFactor-1)/2;
    if ((int)MrParamsL[0] != MrParamsSize+1) 
    {
        DR_Error("%s: expects mr parameter array length %d+1 with "
                 "%d factors (got total %d).\n", 
                 routine, MrParamsSize, NbFactor, (int)MrParamsL[0]);
        return (FAILURE);
    }

    /* Initialization -- needed in 1F and 2F case */
    VolData->Alpha[0] = -999;
    VolData->Alpha[1] = -999.;
    VolData->Alpha[2] = -999.;
    VolData->Beta[0]  = -999.;
    VolData->Beta[1]  = -999.;
    VolData->Beta[2]  = -999.;
    VolData->Rho[0]   = -999.;
    VolData->Rho[1]   = -999.;
    VolData->Rho[2]   = -999.;

    /* Initialize TD data -- not used */
    VolData->NbTDInp = 1;
    VolData->TDInpDate[0] = VolData->BaseDate;

    /* Read Vnfm parameters */
    for (i=0; i<NbFactor; i++) 
    {
         VolData->Beta[i]  = MrParamsL[2+i];
         VolData->Alpha[i] = MrParamsL[2+NbFactor+i];
    }
    for (i=0; i<NbFactor*(NbFactor-1)/2; i++)
    {
         VolData->Rho[i] = MrParamsL[2+2*NbFactor+i];
    }

    /* Assign volatility scale */
    norm = 0.;
    for (i = 0; i < NbFactor; i++) 
        norm += VolData->Alpha[i] * VolData->Alpha[i]; 
    norm = sqrt(norm);
    if (fabs(norm) < ERROR)
    {      
        DR_Error("Total alpha is too small !");
        return (FAILURE);
    }       
    if (IS_EQUAL(VolData->Bbq,1))
    {
        VolData->VolNorm = 0.;
        VolData->VolLogn = norm;
    }
    else if (IS_EQUAL(VolData->Bbq,0))
    {
        VolData->VolNorm = norm;
        VolData->VolLogn = 0.;
    }
    else
    {
        DR_Error("Bbq parameter must be either 0 or 1 !");
        return (FAILURE);
    }       

    /* Smile */   
    for (i = 0; i < VolData->NbVol; i++)
    {
        VolData->SmlLiqDate[i] = 0;
    }

    for (i = 0; i < (int) SmileDatesL[0]; i++)
    {
        /* Check that smile dates are a subset of ATM volatility dates */
        idx = GetDLOffset (VolData->NbVol,
                           VolData->VolDate,
                           SmileDatesL[i+1],
                           CbkEXACT);
        if (idx < 0)
        {
            DR_Error ("Smile date %ld is not among ATM volatility dates",
                      SmileDatesL[i+1]);
            return (FAILURE);
        }
        VolData->SmlLiqDate[idx] = 1;
    }

    /* Copy smile parameters */
    COPY_FROM_ARRAY_L (SmileParamsL, MAXNBDATE * (NBVOLPARS-1),
                   TmpVols);
    l = 0;   /* index for TmpVols vector */

    for (j = 1; j < NBVOLPARS; j++)
    {
        for (i = 0; i < VolData->NbVol; i++)
        {
            if (VolData->SmlLiqDate[i] == 1)
            {
                VolData->Smile[j][i] = TmpVols[l++];
            }
        }
    }

    return (SUCCESS);
} /* Fix3_PackModelData_Tmx */



int Fix3_PackModelData_E2Q (  
        MKTVOL_DATA        *VolData,         /* (O) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (O) Tree data structure         */
        long               *MrDatesL,        /* (I) Mean reversion dates        */
        double             *MrParamsL,       /* (I) Model MR parameters         */
        long               *SmileDatesL,     /* (I) Smile benchmark dates       */
        double             *SmileParamsL)    /* (I) Model smile parameters      */
{    
    int  i;
    int  NbFactor;
    int  MrParamsSize;
    double norm;
    static char routine[] = "Fix3_PackModelData_E2Q";

    /* Check inputs */
    ASSERT_OR_FAIL ((int)MrParamsL[0] >= 2);
    ASSERT_OR_FAIL ((int)SmileParamsL[0] == (E2Q_SMILEPARAM_SIZE-1));

    NbFactor = (int) MrParamsL[1];
    if (NbFactor < 0 || NbFactor > 3)
    {
        DR_Error("%s: invalid number of factors (%d).\n", routine, NbFactor);
        return (FAILURE);
    }
    /* Copy into tree and mktvol data as well */
    Tree->NbFactor    = NbFactor;
    VolData->NbFactor = NbFactor;

    /* Check size */
    MrParamsSize = 2*NbFactor + NbFactor*(NbFactor-1)/2;
    if ((int)MrParamsL[0] != MrParamsSize+1) 
    {
        DR_Error("%s: expects mr parameter array length %d+1 with "
                 "%d factors (got total %d).\n", 
                 routine, MrParamsSize, NbFactor, (int)MrParamsL[0]);
        return (FAILURE);
    }

    /* Initialization -- needed in 1F and 2F case */
    VolData->Alpha[0] = -999;
    VolData->Alpha[1] = -999.;
    VolData->Alpha[2] = -999.;
    VolData->Beta[0]  = -999.;
    VolData->Beta[1]  = -999.;
    VolData->Beta[2]  = -999.;
    VolData->Rho[0]   = -999.;
    VolData->Rho[1]   = -999.;
    VolData->Rho[2]   = -999.;

    /* Initialize TD data -- not used */
    VolData->NbTDInp = 1;
    VolData->TDInpDate[0] = VolData->BaseDate;

    /* Read Vnfm parameters */
    for (i=0; i<NbFactor; i++) 
    {
         VolData->Beta[i]  = MrParamsL[2+i];
         VolData->Alpha[i] = MrParamsL[2+NbFactor+i];
    }
    for (i=0; i<NbFactor*(NbFactor-1)/2; i++)
    {
         VolData->Rho[i] = MrParamsL[2+2*NbFactor+i];
    }

    /* Assign volatility scale */
    norm = 0.;
    for (i = 0; i < NbFactor; i++) 
        norm += VolData->Alpha[i] * VolData->Alpha[i]; 
    norm = sqrt(norm);
    if (fabs(norm) < ERROR)
    {      
        DR_Error("Total alpha is too small !");
        return (FAILURE);
    }       
    if (IS_EQUAL(VolData->Bbq,1))
    {
        VolData->VolNorm = 0.;
        VolData->VolLogn = norm;
    }
    else if (IS_EQUAL(VolData->Bbq,0))
    {
        VolData->VolNorm = norm;
        VolData->VolLogn = 0.;
    }
    else
    {
        DR_Error("Bbq parameter must be either 0 or 1 !");
        return (FAILURE);
    }       

    /* Smile */
    VolData->QLeft     = 1.0 - SmileParamsL[E2Q_QLeft]; 
    VolData->QRight    = 1.0 - SmileParamsL[E2Q_QRight];
    VolData->FwdShift  = SmileParamsL[E2Q_FwdShift];
    VolData->Amap      = SmileParamsL[E2Q_Amap];
    VolData->Bmap      = SmileParamsL[E2Q_Bmap];

    return (SUCCESS);
} /* Fix3_PackModelData_E2Q */




/* Unpack Market Volatility Data */  
/* Task : Unpacks MKTVOL_DATA and FIX3_TREE_DATA structures.    
*/    
int Fix3_UnPackMktVolAndModelData(  
        MKTVOL_DATA        *VolData,         /* (I) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (I) Tree data structure         */
        long               *BaseDatesL,      /* (O) Base dates                  */
        long               *VolDatesL,       /* (O) Vol dates                   */
        long               *VolMatsL,        /* (O) Vol underlying maturity date*/
        double             *VolsL,           /* (O) Volatilities                */
        char               **VolCalibL,      /* (O) Vol calibration flags       */
        long               *MrDatesL,        /* (O) Mean reversion dates        */
        double             *MrParamsL,       /* (O) Model MR parameters         */
        long               *SmileDatesL,     /* (O) Smile benchmark dates       */
        double             *SmileParamsL,    /* (O) Model smile parameters      */
        long               *NumParamsL)      /* (O) Numerical parameters        */
{    
    
    static  char    routine[] = "Fix3_UnPackMktVolAndModelData";

    /* Volatility */
    BaseDatesL[L_VolBaseDate] = VolData->BaseDate;

    COPY_TO_DATES_L (VolDatesL, VolData->NbVol, VolData->VolDate);
    COPY_TO_DATES_L (VolDatesL, VolData->NbVol, VolData->SwapSt);
    COPY_TO_DATES_L (VolMatsL,  VolData->NbVol, VolData->SwapMat);
    COPY_TO_ARRAY_L (VolsL,     VolData->NbVol, VolData->Vol);

    /* Model choice */
    switch (VolData->ModelChoice)
    {
    case (FIX3_ORIGINAL): sprintf (VolCalibL[VOL_ModelChoice],"ORIGINAL"); break;
    case (FIX3_CLASSIC):  sprintf (VolCalibL[VOL_ModelChoice],"CLASSIC");  break;
    case (FIX3_TIMEDEP):  sprintf (VolCalibL[VOL_ModelChoice],"TIMEDEP");  break;
    case (FIX3_SMD):      sprintf (VolCalibL[VOL_ModelChoice],"SMD");      break;
    case (FIX3_TMX):      sprintf (VolCalibL[VOL_ModelChoice],"TMX");      break;
    case (FIX3_E2Q):      sprintf (VolCalibL[VOL_ModelChoice],"E2Q");      break;
    default:              sprintf (VolCalibL[VOL_ModelChoice],"XXXXXXX");
    }

    /* Volatility calibration flags */
    VolCalibL[0] = (char *) (VOLCALIB_SIZE - 1);

    sprintf (VolCalibL[VOL_SkipFlag],         "%c", (VolData->SkipFlag) ?'Y':'N');   
    sprintf (VolCalibL[VOL_CalibFlag],        "%c", (VolData->CalibFlag) ?'Y':'N'); 
    sprintf (VolCalibL[VOL_SmoothingFlag],    "%c", VolData->SmoothingFlag);
    sprintf (VolCalibL[VOL_Freq],             "%c", VolData->Freq);   
    sprintf (VolCalibL[VOL_DCC],              "%c", VolData->DCC);  
    sprintf (VolCalibL[VOL_Type],             "%c", (VolData->VolUnit == 1 ? 'L' : 'N'));  
    sprintf (VolCalibL[VOL_Backbone],         "%c", (IS_EQUAL(VolData->Bbq, 0) ? 'L' : 'N'));

    /* Numerical parameters */
    NumParamsL[0]              = NUMERICS_SIZE-1;
    NumParamsL[NUM_Ppy]        = Tree->Ppy;
    NumParamsL[NUM_NbSigmaMax] = Tree->NbSigmaMax;
    NumParamsL[NUM_CetNbIter]  = VolData->CetNbIter;  
    NumParamsL[NUM_CvDiff]     = Tree->CvDiff;
    NumParamsL[NUM_CvDisc]     = Tree->CvDisc;
    NumParamsL[NUM_InterpType] = EslGetZeroInterpolation();
    NumParamsL[NUM_NbSigmaMQ]  = (long) VolData->NbSigmaMQ;
    NumParamsL[NUM_NckMQ]      = (long) VolData->NckMQ;

    /* Model parameters: Vnfm, smile, numerics */
    return (Fix3_UnPackModelData (VolData,
                                  Tree,
                                  MrDatesL,
                                  MrParamsL,
                                  SmileDatesL,
                                  SmileParamsL));
}    


int Fix3_UnPackModelData_Classic (  
        MKTVOL_DATA        *VolData,         /* (I) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (I) Tree data structure         */
        long               *MrDatesL,        /* (O) Mean reversion dates        */
        double             *MrParamsL,       /* (O) Model MR parameters         */
        long               *SmileDatesL,     /* (O) Smile benchmark dates       */
        double             *SmileParamsL)    /* (O) Model smile parameters      */
{    
    int  i;
    int  NbFactor;
    int  MrParamsSize;
    static char routine[] = "Fix3_UnPackModelData_Classic";

   /* Vnfm parameters */
    NbFactor = VolData->NbFactor;
    MrParamsSize = 2*NbFactor + NbFactor*(NbFactor-1)/2;

    MrParamsL[0] = (double)MrParamsSize+1;
    MrParamsL[1] = NbFactor;
   
    for (i=0; i<NbFactor; i++) 
    {
         MrParamsL[2+i]          = VolData->Beta[i];
         MrParamsL[2+NbFactor+i] = VolData->Alpha[i];
    }
    for (i=0; i<NbFactor*(NbFactor-1)/2; i++)
    {
         MrParamsL[2+2*NbFactor+i] = VolData->Rho[i];
    }

    /* Smile */
    SmileParamsL[0]           = (double)(CLASSIC_SMILEPARAM_SIZE - 1);
    SmileParamsL[CLS_QLeft]     = 1.0 - VolData->QLeft;
    SmileParamsL[CLS_QRight]    = 1.0 - VolData->QRight;
    SmileParamsL[CLS_FwdShift]  = VolData->FwdShift;

    /* Unused outputs */
    MrDatesL[0]    = 0;
    SmileDatesL[0] = 0;

    return (SUCCESS);
}   /* Fix3_UnPackModelData_Classic */



int Fix3_UnPackModelData_TimeDep (  
        MKTVOL_DATA        *VolData,         /* (I) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (I) Tree data structure         */
        long               *MrDatesL,        /* (O) Mean reversion dates        */
        double             *MrParamsL,       /* (O) Model MR parameters         */
        long               *SmileDatesL,     /* (O) Smile benchmark dates       */
        double             *SmileParamsL)    /* (O) Model smile parameters      */
{    
    int  i,j;
    int  idx = 1;
    int  NbFactor;
    int  MrParamsSize;
    static char routine[] = "Fix3_UnPackModelData_TimeDep";

    /* Vnfm parameters */
    NbFactor = VolData->NbFactor;
    MrParamsSize = (2*NbFactor + NbFactor*(NbFactor-1)/2) * VolData->NbTDInp;
    MrParamsL[0] = (double)MrParamsSize+1;
    MrParamsL[1] = NbFactor;
   
    /* Time dependent parameters */
    COPY_TO_DATES_L (MrDatesL, VolData->NbTDInp, VolData->TDInpDate);

    for (i = 0; i < NbFactor; i++)
    {
        for (j = 0; j < VolData->NbTDInp; j++)
        {
            MrParamsL[++idx] = VolData->BetaTD[i][j];
        }
    }
    for (i = 0; i < NbFactor; i++)
    {
        for (j = 0; j < VolData->NbTDInp; j++)
        {
            MrParamsL[++idx] = VolData->AlphaTD[i][j];
        }
    }
    for (i = 0; i < NbFactor*(NbFactor-1)/2; i++)
    {
        for (j = 0; j < VolData->NbTDInp; j++)
        {
            MrParamsL[++idx] = VolData->RhoTD[i][j];
        }
    }

    /* Smile dates and parameters*/
    SmileParamsL[0] = 3 * VolData->NbSmileDates ;
    SmileDatesL[0]  = VolData->NbSmileDates;      
    

    idx = 0; 
    COPY_TO_DATES_L (SmileDatesL, VolData->NbSmileDates, VolData->SmileDate);
    for (i = 0; i < VolData->NbSmileDates; i++)
    {
        SmileParamsL[++idx] = 1.0 - VolData->QLeftTD[i];
    }
    for (i = 0; i < VolData->NbSmileDates; i++)
    {
        SmileParamsL[++idx] = 1.0 - VolData->QRightTD[i];
    }
      for (i = 0; i < VolData->NbSmileDates; i++)
    {
        SmileParamsL[++idx] = VolData->FwdShiftTD[i];
    }

   
    return (SUCCESS);
}   /* Fix3_UnPackModelData_TimeDep */



int Fix3_UnPackModelData_Smd (  
        MKTVOL_DATA        *VolData,         /* (I) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (I) Tree data structure         */
        long               *MrDatesL,        /* (O) Mean reversion dates        */
        double             *MrParamsL,       /* (O) Model MR parameters         */
        long               *SmileDatesL,     /* (O) Smile benchmark dates       */
        double             *SmileParamsL)    /* (O) Model smile parameters      */
{    
    int  i;
    int  NbFactor;
    int  MrParamsSize;
    static char routine[] = "Fix3_UnPackModelData_Smd";

   /* Vnfm parameters */
    NbFactor = VolData->NbFactor;
    MrParamsSize = 2*NbFactor + NbFactor*(NbFactor-1)/2;

    MrParamsL[0] = (double)MrParamsSize+1;
    MrParamsL[1] = NbFactor;
   
    for (i=0; i<NbFactor; i++) 
    {
         MrParamsL[2+i]          = VolData->Beta[i];
         MrParamsL[2+NbFactor+i] = VolData->Alpha[i];
    }
    for (i=0; i<NbFactor*(NbFactor-1)/2; i++)
    {
         MrParamsL[2+2*NbFactor+i] = VolData->Rho[i];
    }

    /* Smile */
    SmileParamsL[0]               = (double)(SMD_SMILEPARAM_SIZE - 1);
    SmileParamsL[SMD_QLeft]     = 1.0 - VolData->QLeft;
    SmileParamsL[SMD_QRight]    = 1.0 - VolData->QRight;
    SmileParamsL[SMD_FwdShift]  = VolData->FwdShift;
    SmileParamsL[SMD_Afac]          = VolData->Afac;
    SmileParamsL[SMD_Bfac]          = VolData->Bfac;
    SmileParamsL[SMD_Cfac]          = VolData->Cfac;
    SmileParamsL[SMD_Dfac]          = VolData->Dfac;

    /* Unused outputs */
    MrDatesL[0]    = 0;
    SmileDatesL[0] = 0;

    return (SUCCESS);
}   /* Fix3_UnPackModelData_Smd */


int Fix3_UnPackModelData_Tmx (  
        MKTVOL_DATA        *VolData,         /* (I) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (I) Tree data structure         */
        long               *MrDatesL,        /* (O) Mean reversion dates        */
        double             *MrParamsL,       /* (O) Model MR parameters         */
        long               *SmileDatesL,     /* (O) Smile benchmark dates       */
        double             *SmileParamsL)    /* (O) Model smile parameters      */
{    
    int  i, j, l;
    int  NbFactor;
    int  MrParamsSize;
    int  NbSmileDates;
    static char routine[] = "Fix3_UnPackModelData_Tmx";

   /* Vnfm parameters */
    NbFactor = VolData->NbFactor;
    MrParamsSize = 2*NbFactor + NbFactor*(NbFactor-1)/2;

    MrParamsL[0] = (double)MrParamsSize+1;
    MrParamsL[1] = NbFactor;
   
    for (i=0; i<NbFactor; i++) 
    {
         MrParamsL[2+i]          = VolData->Beta[i];
         MrParamsL[2+NbFactor+i] = VolData->Alpha[i];
    }
    for (i=0; i<NbFactor*(NbFactor-1)/2; i++)
    {
         MrParamsL[2+2*NbFactor+i] = VolData->Rho[i];
    }

    /* Smile */
    NbSmileDates = 0;
    for (i = 0; i < VolData->NbVol; i++)
    {
        if (VolData->SmlLiqDate[i] == 1) NbSmileDates++;
    }

    SmileDatesL[0] = NbSmileDates;
    l = 1;
    for (i = 0; i < VolData->NbVol; i++)
    {
        if (VolData->SmlLiqDate[i] == 1) 
            SmileDatesL[l++] = VolData->VolDate[i];
    }

    SmileParamsL[0] = NbSmileDates * (NBVOLPARS-1);
    l = 1;
    for (j = 1; j < NBVOLPARS; j++)
    {
        for (i = 0; i < VolData->NbVol; i++)
        {
            if (VolData->SmlLiqDate[i] == 1) 
                SmileParamsL[l++] = VolData->Smile[j][i];
        }
    }

    /* Unused outputs */
    MrDatesL[0]    = 0;

    return (SUCCESS);
}   /* Fix3_UnPackModelData_Tmx */


int Fix3_UnPackModelData_E2Q (  
        MKTVOL_DATA        *VolData,         /* (I) Vol data structure          */
        FIX3_TREE_DATA     *Tree,            /* (I) Tree data structure         */
        long               *MrDatesL,        /* (O) Mean reversion dates        */
        double             *MrParamsL,       /* (O) Model MR parameters         */
        long               *SmileDatesL,     /* (O) Smile benchmark dates       */
        double             *SmileParamsL)    /* (O) Model smile parameters      */
{    
    int  i;
    int  NbFactor;
    int  MrParamsSize;
    static char routine[] = "Fix3_UnPackModelData_E2Q";

   /* Vnfm parameters */
    NbFactor = VolData->NbFactor;
    MrParamsSize = 2*NbFactor + NbFactor*(NbFactor-1)/2;

    MrParamsL[0] = (double)MrParamsSize+1;
    MrParamsL[1] = NbFactor;
   
    for (i=0; i<NbFactor; i++) 
    {
         MrParamsL[2+i]          = VolData->Beta[i];
         MrParamsL[2+NbFactor+i] = VolData->Alpha[i];
    }
    for (i=0; i<NbFactor*(NbFactor-1)/2; i++)
    {
         MrParamsL[2+2*NbFactor+i] = VolData->Rho[i];
    }

    /* Smile */
    SmileParamsL[0]           = (double)(E2Q_SMILEPARAM_SIZE - 1);
    SmileParamsL[E2Q_QLeft]     = 1.0 - VolData->QLeft;
    SmileParamsL[E2Q_QRight]    = 1.0 - VolData->QRight;
    SmileParamsL[E2Q_FwdShift]  = VolData->FwdShift;
    SmileParamsL[E2Q_Amap]      = VolData->Amap;
    SmileParamsL[E2Q_Bmap]      = VolData->Bmap;

    /* Unused outputs */
    MrDatesL[0]    = 0;
    SmileDatesL[0] = 0;

    return (SUCCESS);
}   /* Fix3_UnPackModelData_E2Q */
