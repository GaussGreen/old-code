#include "dll_util.h"

#include "esl_alloc.h"
#include "esl_date.h"

#ifdef ESL_NEW_CURVE
#include "dll_util_new.c"
#else

/* Create counted array of strings for regression test */
char **Alloc_Strng_L (long NbStrings,    /* size of string array      */
                      long MaxStrLen)    /* length of each string     */
{
    int  i;
    char **StrngL;

    StrngL = (char **) DR_Array (CHAR_PTR, 0, NbStrings);
    if (StrngL == NULL) 
    {
        DR_Error ("Failed to allocate CHAR_PTR array");
        return (NULL);
    }

    /* do not allocate 0-th pointer */
    for (i=1;i<=NbStrings;i++)
    {
        StrngL[i] = (char *) DR_Array (CHAR, 0, MaxStrLen);
        if (StrngL[i] == NULL)
        {
            DR_Error ("Failed to allocate CHAR array");
            return (NULL);
        }
    }

    /* instead, use the address to store number of strings */
    StrngL[0] = (char *) NbStrings;

    return (StrngL);
}

/* Free counted array of strings; make sure not to free 0th string! */
int Free_Strng_L (char **StrngL,         /* string array to be freed */        
                  long NbStrings,        /* size of string array     */
                  long MaxStrLen)        /* length of each string    */
{
    int i;

    /* do not free 0th pointer, which is not allocated! */
    for (i=1; i <= NbStrings; i++)
    {
        if (StrngL[i] != NULL)
        {
            Free_DR_Array (StrngL[i], CHAR, 0, MaxStrLen);
        }
    }

    if (StrngL != NULL)
    {
        Free_DR_Array (StrngL, CHAR_PTR, 0, NbStrings);
    }

    return (SUCCESS);
}

  
/* Pack Zero-Curves */  

/* Task  : Packs up a single T_CURVE structure.  
 */  
int PackSingleTCurve( T_CURVE*      crv,  
                      long          Today,
                      long          ValueDate,   
                      char const*   SwapFreq, 
                      char const*   DCC, 
                      char const*   MMB, 
                      long          SpotDays,
                      long const*   ZeroDatesL,
                      double const* ZeroRatesL )
{  
   ASSERT_OR_FAIL ((int)ZeroDatesL[0] == (int)ZeroRatesL[0]);

   crv->Today       = Today ;  
   crv->ValueDate   = ValueDate ;  
   crv->NbZero      = (int)ZeroDatesL[0] ;
   
   crv->SwapFreq   = toupper( SwapFreq[0] ) ; 
   if ( DCC[0] == '0' ) 
   { 
      strcpy( crv->SwapDCC, "360" ) ; 
   } 
   else if ( DCC[0] == '3' ) 
   { 
      strcpy( crv->SwapDCC, "ACT" ) ; 
   } 
   else if ( DCC[0] == '5' ) 
   { 
      strcpy( crv->SwapDCC, "365" ) ; 
   } 
   else 
   { 
      DR_Error( "Value Of Zero DCC not valid" ) ; 
      return( FAILURE ) ; 
   }
   if ( MMB[0] == '0' ) 
   { 
      crv->MMB = 360 ; 
   } 
   else if ( MMB[0] == '5' ) 
   { 
      crv->MMB = 365 ; 
   } 
   else 
   { 
      DR_Error( "Value Of Zero MMB not valid" ) ; 
      return( FAILURE ) ; 
   } 
 
   crv->SpotDays = SpotDays ; 
        
   COPY_FROM_DATES_L (ZeroDatesL, MAXNBDATE, crv->ZeroDate);
   COPY_FROM_ARRAY_L (ZeroRatesL, MAXNBDATE, crv->Zero);

   return( SUCCESS ) ;  
}  
  

/* Task : Packs up all three zero curves.  
*/  
int PackTCurves(  T_CURVE    *Curves,  
                  long       *BaseDatesL,  
                  char       **ZeroFrequenciesL,
                  char       **ZeroDCCsL,
                  char       **ZeroMMBsL,
                  long       *ZeroSpotDaysL,
                  long       *Zero0DatesL,
                  double     *Zero0RatesL,
                  long       *Zero1DatesL,
                  double     *Zero1RatesL,
                  long       *Zero2DatesL,
                  double     *Zero2RatesL )
{ 
   ASSERT_OR_FAIL ((int)BaseDatesL[0] == (L_BASEDATE_SIZE-1));

   if (PackSingleTCurve(&Curves[0],  
                         BaseDatesL[L_Today],  
                         BaseDatesL[L_ValueDate0], 
                         ZeroFrequenciesL[1], 
                         ZeroDCCsL[1],
                         ZeroMMBsL[1],       
                         ZeroSpotDaysL[1],  
                         Zero0DatesL,
                         Zero0RatesL) != SUCCESS ||
       PackSingleTCurve(&Curves[1],    
                         BaseDatesL[L_Today],    
                         BaseDatesL[L_ValueDate1],   
                         ZeroFrequenciesL[2], 
                         ZeroDCCsL[2],
                         ZeroMMBsL[2],         
                         ZeroSpotDaysL[2],    
                         Zero1DatesL,
                         Zero1RatesL) != SUCCESS ||    
       PackSingleTCurve(&Curves[2],    
                         BaseDatesL[L_Today],    
                         BaseDatesL[L_ValueDate2],   
                         ZeroFrequenciesL[3],   
                         ZeroDCCsL[3],
                         ZeroMMBsL[3],
                         ZeroSpotDaysL[3],    
                         Zero2DatesL,
                         Zero2RatesL) != SUCCESS)    
      return FAILURE;    

   return SUCCESS;
}  
  

/* Unpack Zero-Curves */  
/* Task  : Unpacks up a single T_CURVE structure.  
 */  
int UnPackSingleTCurve( 
                      T_CURVE const* Curve,  
                      long       *ValueDate,   
                      char       *SwapFreq, 
                      char       *DCC, 
                      char       *MMB, 
                      long       *SpotDays,
                      long       *ZeroDatesL,
                      double     *ZeroRatesL )
{  
   *ValueDate = Curve->ValueDate;  
   *SwapFreq  = Curve->SwapFreq;

   if (toupper(Curve->SwapDCC[0]) == 'A')
   {
        strcpy(DCC, "3");
   }
   else if (!strcmp (Curve->SwapDCC, "360"))
   {
        strcpy(DCC, "0");
   }
   else if (!strcmp (Curve->SwapDCC, "365"))
   {
        strcpy(DCC, "5");
   }
   else 
   { 
        DR_Error( "Value Of Zero DCC not valid" ) ; 
        return( FAILURE ) ; 
   }

   if ( Curve->MMB == 360 )
   { 
        strcpy(MMB, "0");
   } 
   else if ( Curve->MMB == 365 ) 
   { 
        strcpy(MMB, "5");
   } 
   else 
   { 
      DR_Error( "Value Of Zero MMB not valid" ) ; 
      return( FAILURE ) ; 
   } 
 
   *SpotDays = Curve->SpotDays; 
        
   ZeroDatesL[0] = (long)Curve->NbZero;
   ZeroRatesL[0] = (double)Curve->NbZero;
   
   COPY_TO_DATES_L (ZeroDatesL, Curve->NbZero, Curve->ZeroDate);
   COPY_TO_ARRAY_L (ZeroRatesL, Curve->NbZero, Curve->Zero);

   return( SUCCESS ) ;  
}  
  

/* Task : Unpacks up all three zero curves.  
*/  
int UnPackTCurves(  
                  T_CURVE const* Curves,  
                  long       *BaseDatesL,  
                  char       **ZeroFrequenciesL,
                  char       **ZeroDCCsL,
                  char       **ZeroMMBsL,
                  long       *ZeroSpotDaysL,
                  long       *Zero0DatesL,
                  double     *Zero0RatesL,
                  long       *Zero1DatesL,
                  double     *Zero1RatesL,
                  long       *Zero2DatesL,
                  double     *Zero2RatesL )
{ 
   BaseDatesL[0] = (int)(L_BASEDATE_SIZE - 1);
   BaseDatesL[L_Today] = Curves[0].Today;  

   if ( UnPackSingleTCurve( &Curves[0],  
                            &(BaseDatesL[L_ValueDate0]),
                            ZeroFrequenciesL[1],
                            ZeroDCCsL[1],
                            ZeroMMBsL[1],   
                            &(ZeroSpotDaysL[1]),
                            Zero0DatesL,
                            Zero0RatesL  
                         )  == FAILURE )  
   {  
      return( FAILURE ) ;  
   }  
   if ( UnPackSingleTCurve( &Curves[1],    
                            &(BaseDatesL[L_ValueDate1]),
                            ZeroFrequenciesL[2],
                            ZeroDCCsL[2],
                            ZeroMMBsL[2], 
                            &(ZeroSpotDaysL[2]),
                            Zero1DatesL,
                            Zero1RatesL 
                         )  == FAILURE )    
   {    
      return( FAILURE ) ;    
   }    
   if ( UnPackSingleTCurve( &Curves[2],    
                            &(BaseDatesL[L_ValueDate2]),
                            ZeroFrequenciesL[3],
                            ZeroDCCsL[3],
                            ZeroMMBsL[3],   
                            &(ZeroSpotDaysL[3]),
                            Zero2DatesL,
                            Zero2RatesL    
                         )  == FAILURE )    
   {    
      return( FAILURE ) ;    
   }    

   return( SUCCESS ) ;
}  
  

#endif
