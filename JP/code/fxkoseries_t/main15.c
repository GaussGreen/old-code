/****************************************************************************/
/*      SubMain Routine                                                     */
/****************************************************************************/
/*      MAIN15.C                                                            */
/****************************************************************************/
/*                                                                          */
/*    This file contains the main() routine for pricing FX knock-out option */
/*    series.                                                               */
/*                                                                          */
/****************************************************************************/


/*
$Header$
*/

#include         <stdio.h>
#include         <stdlib.h>
#include         "template15.h"



        
          


/*****  main  ***************************************************************/
/*
*       Main subroutine.
*/
int     main (void)
{


    T_CURVE          t_curve[2][3];      /* Zero curve data          */
    MKTVOL_DATA      mktvol_data[2];     /* Market volatility data     */
    FX_DATA          fx_data;            /* Structure of fx data     */
    
    HYB3_TREE_DATA        tree_data;          /* Structure of tree data   */

    FXKOSERIES_DATA  fxkoseries_data;   
    OPT_OUT_DATA     opt_out_data;       /* Structure of output data */

    FILE             *fp = NULL;         /* File pointer for output  */

    int              status = FAILURE;   /* Error status = FAILURE   */
        
    



    ZeroInterpTypeFlag = 0;        /* 0=Linear Zero Cpn; 1=Flat Fwd */
    ZeroInterpTypeFlagStub = 0;    /* 0=Linear Stub; 1=Flat Stub    */
    
    Hyb3_Tree_Init(&tree_data);
    tree_data.TreeType = TTYPE_FX2IR;

    printf ("Calibrating ...\n");


    if (Fxkoseries_Manager(t_curve,
                           mktvol_data,
                           &fx_data,
                           &fxkoseries_data,
                           &tree_data) == FAILURE)
    {                    
        goto RETURN;   
    }                                            

                                                                                   

    if (Fxkoseries_Schedule(fx_data.ValueDate,  /* Set up the product timeline */
                            &fxkoseries_data,
                            &tree_data) == FAILURE)
    {
        goto RETURN;
  
    }
                                            


    if (Hyb3_Build_Tree(TRUE,
                   t_curve,             /* Build the tree (i.e. populate the */
                   mktvol_data,         /* timeline,  calibrate  drifts  and */
                   &fx_data,            /* populate the timeline with deter- */
                   NULL,                /* ministic values)                  */
                   &tree_data) == FAILURE)
    {
        goto RETURN;
        
    }
       
    

    if (Print_Fxkoseries(t_curve,   /* Print term struct and deal info */
                         &tree_data) == FAILURE)
    {
        goto RETURN;
        
    }



    printf ("Calculating FX KO series price ...\n\n");



    if (Calc_Fxkoseries(mktvol_data,  /* Calculate option price */
                        &fxkoseries_data,
                        &tree_data,
                       &opt_out_data) == FAILURE)
    {
        DR_Error ("Could not calculate FX KO series price! (main)");
        goto RETURN;
        
    }  /* if */                                                


    /* Print price to the command line */
    printf("Calculating FX knock-");
    if (fxkoseries_data.KnockIoO == 'I')
        printf("in");
    else
        printf("out");

    printf(" price ...\n\n");


    printf ("Price:     %20.6f \n\n", opt_out_data.Option);
    if (fxkoseries_data.KnockIoO == 'I')
        printf("Fx Series: %20.6f \n\n", opt_out_data.Price[0]);

   
                 

    status = SUCCESS;

    RETURN:



    /* Print price in output file for wrapper */
    fp = fopen ("price", "w");                      
    
    if (fp == NULL)
    {         
        DR_Error("Could not open file price! (main)\n");
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


    Hyb3_Tree_Free (&tree_data);   /* Free memory allocated for the tree */

    return (status);

}  /* main */
