/****************************************************************************/
/*                                                                          */
/*      SubMain Routine (Ladder)                                            */
/*                                                                          */
/****************************************************************************/
/*      MAIN39.c                                                            */
/****************************************************************************/

/*
$Header$
*/
    


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fix123head.h" 
#include "template39.h"

        
/*****  main  ***************************************************************/
/**
*       Main subroutine.
*/
int     main (void)
{

    MKTVOL_DATA         mktvol_data;    /* Structure of vol data            */
    LADDER_DATA         ladder_data;    /* Structure of ladder data         */
    T_CURVE             t_curve[3];     /* Structure of zero curve data     */
    FIX3_TREE_DATA      tree_data;      /* Structure of tree data           */
    OPT_OUT_DATA        opt_out_data;   /* Structure of output data         */
    int                 status = FAILURE;        /* Status */                
    int                 i;
    char                ErrorMsg[MAXBUFF];
    FILE                *fp = NULL;
   

    ladder_data.State = NULL;
	
    printf ("Calibrating ...\n");

    /* 
     *  Initialize tree data structure to NULL.
     */
    Fix3_Tree_Init (&tree_data);
    MktVol_Init(&mktvol_data);
    Opt_Out_Data_Init(&opt_out_data);
    Ladder_Data_Init(&ladder_data);

    /* 
    *   I/O manager 
    */
    if (Ladder_Manager (t_curve,
                        &mktvol_data,
                        &ladder_data,
                        &tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    if (Ladder_PreProcess (t_curve,
                           &mktvol_data,
                           (t_curve[0]).ValueDate, 
                           &ladder_data,
                           &tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* 
    *   Build time line.
    */
    if (Ladder_Schedule ((t_curve[0]).ValueDate,         
                         t_curve,
                         &mktvol_data,
                         &ladder_data,
                         &tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /*
    *   CET 
    */
    if (Fix3_Cet_Main (TRUE,          /* output to CET.prn */
                  t_curve,
                  &mktvol_data,
                  &tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /*
    *   Build tree.
    */
    if (Fix3_Build_Tree (t_curve,
                    &mktvol_data,
                    &tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }
    
    /* 
    *   Print debug information in an ascii file. 
    */
    if (Print_Ladder (t_curve,                                        
                      &mktvol_data,
                      &ladder_data,
                      &tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

#ifdef DEBUG
   Print_Flows_Deal(&mktvol_data,&ladder_data,&tree_data); 
#endif
 

    printf ("Calculating option price ...\n\n");

    /*
    *   Main calculation routine.
    */
    if (Calc_Ladder (&mktvol_data,                                  
                     &ladder_data,
                     &tree_data,
                     &opt_out_data) == FAILURE)
    {
        DR_Error ("Could not calculate price of ladder! (main)");
        goto FREE_MEM_AND_RETURN;
    }
    
    printf ("Option price:      %16.4f \n", opt_out_data.Option);
    printf ("Current notional:  %16.4f \n", ladder_data.InitOuts);

    if (ladder_data.CalcStats == 'Y')
    {
        printf ("Exer Prob:         %12.4f \n", opt_out_data.prob_out_data[EXER_EVENT].TotalEventProb);
        printf ("Exer Time:         %12.4f \n", opt_out_data.prob_out_data[EXER_EVENT].ExpEventTime);
        printf ("Exer Time s.d.:    %12.4f \n", opt_out_data.prob_out_data[EXER_EVENT].StdEventTime);
        printf ("Fugit:             %12.4f \n", opt_out_data.prob_out_data[EXER_EVENT].Fugit);

        printf ("\nExercise probabilities - can be a partial list\n");

        for (i=0; i<opt_out_data.prob_out_data[EXER_EVENT].Count; i++)
            printf ("%ld   %12.6f\n", 
			    YMDDateFromIRDate(opt_out_data.prob_out_data[EXER_EVENT].EventDate[i]), 
			    opt_out_data.prob_out_data[EXER_EVENT].EventProb[i]);

	printProbabilityInFile (&opt_out_data);
    }

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    /* 
    *   Print price in output file for wrapper.
    */
    fp = fopen ("price", "w");                      
    
    if (fp == NULL)
    {
        sprintf (ErrorMsg, "Could not open file price! (main)");                 
        DR_Error (ErrorMsg);
    }
    else if (status == SUCCESS)
    {
        fprintf (fp, "%f", opt_out_data.Option);
        fclose (fp);
    }
    else
    {
        fclose (fp);
    }

    Fix3_Tree_Free (&tree_data);                                             
    Ladder_Data_Free(&ladder_data);
	
	if (ladder_data.State != NULL)
	{
             Free_DR_Matrix(ladder_data.State, DOUBLE,
                          ladder_data.FirstResetI -1 , (ladder_data.EndResetI - 1L),
                          0, (ladder_data.NbStates - 1L));
	}
	
    DestroyZeroCurve(&t_curve[0]);
    DestroyZeroCurve(&t_curve[1]);
    DestroyZeroCurve(&t_curve[2]);

    return (status);

}  /* main */

