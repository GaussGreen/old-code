/****************************************************************************/
/*      Construct the FX KO option series time line.                        */
/****************************************************************************/
/*      TIME15.C                                                            */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "template15.h"





/*****  Fxkoseries_Schedule **********************************************/
/*
 *       Sets up the time line for the fx option series.
 *
 *       In this  product, the critical  dates are  simply the  exercise
 *       dates.
 *
 *       This routine deals with the detailed date generation from input
 *       parameters and the subsequent addition of dates to the timeline
 *       critical date list.
 */

int   Fxkoseries_Schedule 
          (long               ValueDate,      /* (I) Value date             */
           FXKOSERIES_DATA   *fxkoseries_data,/* (I) FX series deal data    */
           HYB3_TREE_DATA         *tree_data)      /* (O) Tree data              */
{






    CRIT_DATE   *CritDate = NULL;        /* Critical date list             */
    EVENT_LIST  *ExerEventList = NULL;
    EVENT_LIST  *KoEvList = NULL;


    int         NbCritDate = 0;          /* Number of critical dates       */

    int         i;                       /* General index counter          */
    int         status = FAILURE;        /* Status = FAILURE initially     */






    /* Allocate empty critical date list and then add dates successively.  */
    CritDate = (CRIT_DATE *)DR_Array(CRITDATE, 0, 0);
    if (CritDate == NULL)
    {
        DR_Error("Unable to allocate memory "
                "for critical dates array !(Fxkoseries_Schedule)");
        goto RETURN;
        
    }  /* if */

    
    /*  Always add value date to critical date list.  */
    if (Add_To_DateList (&NbCritDate,
                         &CritDate,              
                         ValueDate,
                         NBCRITDATE,             /* No specific type */
                         0, 0, 0, 0, 0,              
                         0, 0, 0) == FAILURE)
    {
        goto RETURN;
        
    }  /* if */





    /*************************************/
    /*  TYPE  0: KNOCK OUT DATES         */ 
    /*************************************/

    if (fxkoseries_data->KoFreq == 'N')
    {               
        tree_data->CritType[0] = 'L';    /* Linear */
    }           
    else    
    {                                                             
        tree_data->CritType[0] = 'D';    /* Discrete */
    }
    tree_data->NbZeros[0] = 0; /* No zeros associated */
   

    /* Generate knock out/in event list  */
    KoEvList = DrNewEventListFromFreq 
                        (fxkoseries_data->NbKoDates,
                         fxkoseries_data->KoDates,
                         fxkoseries_data->KoFreq,
                         'N',         /* Stubs not allowed           */
                         'Y',         /* Input dates must be in list */
                         fxkoseries_data->LoBarrier,
                         fxkoseries_data->HiBarrier,
                         fxkoseries_data->Rebate, 
                         NULL, NULL);
    
    if (KoEvList == NULL)
    {
        goto RETURN;
       
    }
           
        
    /*   Drop dates which fall before value date, redefining the start and  */
    /*   end of the  exer window and interp'ing for the  new start and end. */
    if (TimeInterp (ValueDate,                                      
                    "Knock-out", 
                    tree_data->CritType[0],
                    &(KoEvList->NbEntries),
                    KoEvList->Dates,
                    KoEvList->Curve[0], 
                    KoEvList->Curve[1],
                    KoEvList->Curve[2],
                    NULL, NULL) == FAILURE)
    {
        goto RETURN;
        
    }

   
    /*   Add dates to crit datelist */
    for (i = 0; i < KoEvList->NbEntries; i++)          
    {   

        /* Only keep knock out/in dates falling after value date */
        if (KoEvList->Dates[i] >= ValueDate)   
        {
            if (Add_To_DateList(&NbCritDate,  /* Notification */
                                &CritDate,
                                KoEvList->Dates[i],
                                0,
                                KoEvList->Curve[0][i],
                                KoEvList->Curve[1][i],
                                KoEvList->Curve[2][i], 
                                0, 0, 
                                0, 0, 0) == FAILURE)        
            {
               goto RETURN;
                
            }  /* if */

        } /* If */

    }  /* for i */

    

    /*************************************/
    /*  TYPE  1: EXERCISE DATES          */ 
    /*************************************/
    tree_data->CritType[1] = 'D';
    tree_data->NbZeros[1]  =  0;


    /* Process the option inputs to generate the event list */
    /* which will be a temporary EVENT_LIST structure       */
    ExerEventList = DrNewEventListFromFreq
                       (fxkoseries_data->NbExer,
                        fxkoseries_data->Exer,
                        fxkoseries_data->ExerFreq,
                        'N',  /* No stubs allowed            */
                        'Y',  /* Input dates must be in list */
                        fxkoseries_data->Strike,
                        NULL, NULL, NULL, NULL);

    if (ExerEventList == NULL)
    {
        goto RETURN;
       
    }


    if (TimeInterp(ValueDate,
                   "Exer",
                   tree_data->CritType[1],
                   &(ExerEventList->NbEntries),
                   ExerEventList->Dates,
                   ExerEventList->Curve[0], 
                   NULL, NULL, NULL, NULL) == FAILURE)
    {
        goto RETURN;
    }


    for (i = 0; i < ExerEventList->NbEntries; i++)          
    {
        if (Add_To_DateList(&NbCritDate,
                            &CritDate,
                            ExerEventList->Dates[i],
                            1,
                            ExerEventList->Curve[0][i], 0.0, 0.0, 0.0, 0.0,
                            0L, 0L, 0L) == FAILURE)
        {
            goto RETURN;
        }
    } /* For i */


    


    /********************************/
    /*  TYPES 2 TO NBCRITDATE:  N/A */ 
    /********************************/       
                    
    for (i = 2; i < NBCRITDATE; i++)
    {
        tree_data->CritType[i] = 'D';
        tree_data->NbZeros[i]  =  0;
    }




    /* 
     *   Finally construct the time line using 'I'ncreasing time steps.
     */
    if (Hyb3_Time_Line (ValueDate,                      
                   NbCritDate,
                   CritDate,
                   'I',                            
                   tree_data) == FAILURE)
    {
        goto RETURN;
        
    }  /* if */




    
    status = SUCCESS;


  RETURN:                            

    Free_DR_Array(CritDate, CRITDATE, 0, NbCritDate-1);
    DrFreeEventList(ExerEventList);
    DrFreeEventList(KoEvList);
        
    if (status ==FAILURE)
    {
        DR_Error("Fxkoseries_Schedule: Failed.\n");
    
    }  /* if */
    
    return (status);




}  /* End of Fxkoseries_Schedule */
