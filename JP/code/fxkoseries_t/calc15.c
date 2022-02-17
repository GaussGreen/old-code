/****************************************************************************/
/*      Main routine for foreign exchange option series: calculation going  */
/*      backward in the tree.                                               */
/****************************************************************************/
/*      CALC15.C                                                            */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "template15.h"





/*****  Calc_Fxkoseries  ***************************************************/
/*
*       Main calculation routine: discounted expected value of cash-flows
*       going backward in the tree.
*/
int     Calc_Fxkoseries
    (MKTVOL_DATA       *mktvol_data,           /* (I) Market vol data       */
     FXKOSERIES_DATA   *fxkoseries_data,       /* (I) Deal data             */
     HYB3_TREE_DATA         *tree_data,             /* (I) Tree data             */
     OPT_OUT_DATA      *opt_out_data)          /* (O) Output data           */
{




    HYB3_DEV_DATA     dev_data;

    TSLICE   FxKoSeries = NULL; /* Variable containing fx option series    */ 
    TSLICE   FxSeries   = NULL; 

    int
            CoPInteger,        /* Call or put in integer format            */
            ExerciseFlag,      /* TRUE if current date is an exercise date */
            KoFlag;

    long
            LastExerciseDate,  /* End of the tree timeline                 */
            CurrentDate,       /* Current timeline date                    */
            ValueDate;      

    double
            Strike,            /* Strike of equity option                  */
            LoBarrier,
            HiBarrier,
            Rebate,
            Notional;          /* Notional in payment currency             */

    int
            IdxCurveF,
            DCurveF,
            IdxCurveD,
            DCurveD,

            T,                /* Total number of period in the hybrids tree*/
            t,                /* Current time period                       */
            Offset3D,
            status = FAILURE; /* Error status = FAILURE initially          */

    



    /* initialise dev structure */
    Hyb3_Dev_Init (&dev_data);

    /* Total size of tree timeline */     
    T   = tree_data->NbTP;
    ValueDate = tree_data->TPDate[0];

    /* Assigment of discount and index curves in the engine */
    IdxCurveF = tree_data->CvDiff[0];
    DCurveF   = tree_data->CvDisc[0];
    IdxCurveD = tree_data->CvDiff[1];
    DCurveD   = tree_data->CvDisc[1];

   
        
    /* Allocate tree variables */
    FxKoSeries = Hyb3_Alloc_Slice(tree_data, 3);
    FxSeries = Hyb3_Alloc_Slice(tree_data, 3);
    if (FxKoSeries == NULL ||
        FxSeries   == NULL)
    {
        goto RETURN;
    }
      
    if (Hyb3_Dev_Alloc(&dev_data, tree_data) == FAILURE)
    {
        goto RETURN;
    }


    Notional = fxkoseries_data->Notional;
    CoPInteger = fxkoseries_data->CoP == 'C'?  1 : -1;
    LastExerciseDate = fxkoseries_data->Exer[fxkoseries_data->NbExer - 1];

        
       
    for (t = T; t >= 0; t--)
    {     

        
        CurrentDate  = tree_data->TPDate[t];
        KoFlag     = tree_data->TPType[0][t];
        LoBarrier  = tree_data->CritDate[0][t].Value[0];
        HiBarrier  = tree_data->CritDate[0][t].Value[1];
        Rebate     = tree_data->CritDate[0][t].Value[2];


        ExerciseFlag = tree_data->TPType[1][t];
        Strike       = tree_data->CritDate[1][t].Value[0];



        /*  Update tree */
        if (Hyb3_Lattice(&dev_data,
                    t,
                    T,
                    mktvol_data,
                    tree_data) == FAILURE)
        {
            goto RETURN;
        }  
            
        if (CurrentDate <= LastExerciseDate)
        {
            
            if (Hyb3_Option_Series_t(FxSeries,
                                dev_data.FxSpot,
                                Notional,
                                Strike,
                                ExerciseFlag,
                                CoPInteger,
                                t,
                                T,
                                DCurveD,
                                DISC_3D_CUPS,
                                &dev_data,
                                tree_data) == FAILURE)
            {
                goto RETURN;
            }
        
            
            if (fxkoseries_data->KnockIoO == 'I')
            {

                if (Hyb3_KiOption_t(FxKoSeries,
                               FxSeries,
                               3,
                               dev_data.FxSpot,
                               KoFlag,
                               LoBarrier,
                               HiBarrier,
                               fxkoseries_data->IoO,
                               fxkoseries_data->SmoothFlag,
                               t,
                               T,
                               DCurveD,
                               DISC_3D_CUPS,
                               &dev_data,
                               tree_data) == FAILURE)
                {
                   goto RETURN;  
                }  
            }
            else
            {
                if (Hyb3_KoOption_t(FxSeries,
                               3,
                               dev_data.FxSpot,
                               KoFlag,
                               LoBarrier,
                               HiBarrier,
                               0,
                               (void *)&Rebate,
                               fxkoseries_data->IoO,
                               fxkoseries_data->SmoothFlag,
                               t,
                               DISC_3D_CUPS,
                               tree_data) == FAILURE)
                {
                    goto RETURN;
                } 
            
            } 

        }

    }  /* End of loop backwards in tree */
                     
 
    /* Prepare to address TSLICES */
    Offset3D = Hyb3_Node_Offset(3,0,0,0,tree_data);


    if (fxkoseries_data->KnockIoO == 'I')
    {
        opt_out_data->Option = ((double *)FxKoSeries + Offset3D)[0];
        opt_out_data->Price[0] = ((double *)FxSeries + Offset3D)[0];
    }
    else
    {
        opt_out_data->Option = ((double *)FxSeries + Offset3D)[0];
    }

    
        

    /* Adjust for long or short the option */
    if (fxkoseries_data->LoS == 'S')
    {
        opt_out_data->Option = -opt_out_data->Option;
        opt_out_data->Price[0] = -opt_out_data->Price[0];
    }
        
        

    status = SUCCESS;

    RETURN:


    Hyb3_Free_Slice(FxSeries,tree_data, 3);
    Hyb3_Free_Slice(FxKoSeries, tree_data, 3);
    Hyb3_Dev_Free(&dev_data, tree_data);

                                    
    return (status);


}  /* Calc_Fxkoseries */
