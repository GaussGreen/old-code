#include <fix123head.h>    
#include <dll_util.h>    
#include <math.h>    

/* Pack Market Volitility Data */  
/* Task : Packs up MKTVOL_DATA and FIX3_TREE_DATA structures.    
*/    
int Fix3_PackMktVolAndTreeData(  
        MKTVOL_DATA*        VolData,         /* (O) Vol data structure          */
        FIX3_TREE_DATA*     Tree,            /* (O) Tree data structure         */
        long const*         BaseDatesL,      /* (I) Base dates                  */
        long const*         VolDatesL,       /* (I) Vol dates                   */
        long const*         VolMatsL,        /* (I) Vol underlying maturity date*/
        double const*       VolsL,           /* (I) Vols                        */
        char **        VolFreqL,       /* (I) Vol underlying frequency    */
        char **        VolDCCL,        /* (I) Vol underlying DCC          */
        char **        VolTypeL,       /* (I) Vol type ('N'orm, 'L'og)    */
        char **        VolCalibL,      /* (I) Vol calib flags             */
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
   
   VolData->BaseDate = LDate2ADate(BaseDatesL[L_VolBaseDate]);   

   VolData->NbVol    = (int)VolDatesL[0] ;    

   COPY_FROM_DATES_L (VolDatesL, MAXNBDATE, VolData->VolDate);
   COPY_FROM_DATES_L (VolDatesL, MAXNBDATE, VolData->SwapSt);
   COPY_FROM_DATES_L (VolMatsL,  MAXNBDATE, VolData->SwapMat);
   COPY_FROM_ARRAY_L (VolsL,     MAXNBDATE, VolData->Vol);

   for ( i = 0; i < (int)VolDatesL[0]; i++ )
   {
       VolData->VolUsed[i] = TRUE;
   }

   VolData->Freq  = (char)toupper( VolFreqL[1][0] );    
   VolData->DCC   = (char)toupper( VolDCCL[1][0] );    

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
  
   return( SUCCESS );    

}    




/* Unpack Market Volitility Data */  
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
   BaseDatesL[L_VolBaseDate] = ADate2LDate(VolData->BaseDate);


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
