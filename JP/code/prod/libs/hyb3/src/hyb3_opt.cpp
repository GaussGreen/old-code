/****************************************************************************/
/*      Calculation of option price in the lattice.                         */
/****************************************************************************/
/*      OPT.c                                                               */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"


/*****  Hyb3_Option_t  ***********************************************************/
/**
*       Standard option.
*/
int     Hyb3_Option_t(TSLICE      Opt,             /**< (I/O) Option prices         */
                 TSLICE      Under,           /**< (I) Underlying prices       */
                 double      Notional,        /**< (I) Notional                */
                 double      Strike,          /**< (I) Strike                  */
                 long        ExerFlag,        /**< (I) Exercise flag           */
                 int         CoP,             /**< (I) =1 for call, -1 for put */
                 int         t,               /**< (I) Current time point      */
                 int         T,               /**< (I) Total number of point   */
                 int         DCurve,          /**< (I) Discount curve          */
                 int         DMode,           /**< (I) Dim of DEV and slices   */
                 HYB3_DEV_DATA    *dev_data,       /**< (I) Hyb3_Dev data structure      */
                 HYB3_TREE_DATA   *tree_data)      /**< (I) Tree data structure     */
{

    double  *OptL; 
    double  *UnderL;

    int	                                
            Top1,    Bottom1,       /* Tree limits (1rst dim) */
           *Top2,   *Bottom2,       /* Tree limits (2nd dim)  */
          **Top3,  **Bottom3,       /* Tree limits (3rd dim)  */
            i, j, k,                /* Node indices           */
            offset,
            sliceDim,
            status = FAILURE;       /* Error status	          */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];
    sliceDim = Hyb3_Slice_Dim( DMode);



    if (Hyb3_Dev(Opt,
            t,
            T,
            DCurve,
            DMode,
            dev_data,
            tree_data) == FAILURE)
    {
        goto RETURN;
                    
    }  /* if */


    if (ExerFlag)
    {
 
        switch (sliceDim)
        {
        case 1:
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            OptL   = (double *)Opt + offset;
            UnderL = (double *)Under + offset;


            /* 
            *   If this is the last exercise date OptL has been initialized to
            *   zero so the exercise decision is Max (intrinsic value, 0).
            *   If this is not the last exercise date we have an American 
            *   exercise decision Max (intrinsic value, live option).
            */

            for (i = Bottom1; i <= Top1; i ++)                                      
            {
                OptL[i] = MAX (Notional * CoP * (UnderL[i] - Strike), OptL[i]);

            }  /* for i */
            break;

        case 2:
            for (i = Bottom1; i <= Top1; i ++)
            {

                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                OptL   = (double *)Opt + offset;
                UnderL = (double *)Under + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    OptL[j] = MAX (Notional * CoP * 
                         (UnderL[j] - Strike), OptL[j]);

                }  /* for j */	
            }
            break;
        case 3:
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    OptL   = (double *)Opt + offset;
                    UnderL = (double *)Under + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        OptL[k] = MAX (Notional * CoP * 
                             (UnderL[k] - Strike), OptL[k]);

                    }  /* for k */
                }
            }
            break;

        default:
            DR_Error("Unknown Dev Mode (Hyb3_Option_t)! \n");
            break;
        }  /*switch slice dim*/

    }  /* if */    

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Option_t */
/*****************************************************************************/
/**
   Standard option plus the following stats:
      - prob of exercise
      - expected time to exercise
      - expected time to exercise squared
  
   Only perform calculations for stats info if slice is not NULL             */

int     Hyb3_OptionPlus_t(TSLICE   Opt,             /**< (I/O) Option prices        */
                     double   *ExerProb,       /**< (I/O) exer prob            */
                     double   *ExerTime,       /**< (I/O) time to exer         */
                     double   *ExerTimeSqr,    /**< (I/O) time to exer squared */
                     TSLICE   Under,           /**< (I) Underlying prices      */
                     double      Notional,     /**< (I) Notional                */
                     double      Strike,       /**< (I) Strike                  */
                     long        ExerFlag,     /**< (I) Exercise flag           */
                     int         CoP,          /**< (I) =1 for call, -1 for put */
                     int         SmoothingOn,  /**< (I) TRUE = smoothing on     */
                     double      *AuxSlice,    /**< (I) Aux Slice for smoothing */
                     int         t,            /**< (I) Current time point      */
                     int         T,            /**< (I) Total number of point   */
                     int         DCurve,       /**< (I) Discount curve          */
                     int         DMode,        /**< (I) Dim of DEV and slices   */
                     HYB3_DEV_DATA    *dev_data,    /**< (I) Hyb3_Dev data structure      */
                     HYB3_TREE_DATA   *tree_data)   /**< (I) Tree data structure     */
{

    double	*OptL; 
    double	*UnderL;
    double  *AuxSliceL   = NULL;
    double  *ExerProbL   = NULL;
    double  *ExerTimeL   = NULL;
    double  *ExerTimeSqrL = NULL;


    double  Step;                   /* StepSize for smoothing */
    
    double  CurrTimeInYrs;        /* for ExerTime stats     */
    double  CurrTimeInYrsSqr;       /* for ExerTimeSqr stats  */
    int	                                
            Top1,    Bottom1,       /* Tree limits (1rst dim) */
           *Top2,   *Bottom2,       /* Tree limits (2nd dim)  */
          **Top3,  **Bottom3,       /* Tree limits (3rd dim)  */
            i, j, k,                /* Node indices           */
            offset,
            sliceDim,
            status = FAILURE;	    /* Error status	          */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];
    sliceDim = Hyb3_Slice_Dim( DMode);


    if (Hyb3_Dev(Opt,
            t,
            T,
            DCurve,
            DMode,
            dev_data,
            tree_data) == FAILURE)
    {
        goto RETURN;
                    
    }  /* if */

    if (ExerProb != NULL)
    {
       if(Hyb3_Ev(ExerProb,t,T,DMode,dev_data,tree_data) == FAILURE) goto RETURN;
    }

    if (ExerTime != NULL)
    {
       if(Hyb3_Ev(ExerTime,t,T,DMode,dev_data,tree_data) == FAILURE) goto RETURN;
    }

    if (ExerTimeSqr != NULL)
    {
       if(Hyb3_Ev(ExerTimeSqr,t,T,DMode,dev_data,tree_data) == FAILURE) goto RETURN;
    }


    if (ExerFlag)
    {   

        CurrTimeInYrs    = Daysact(tree_data->TPDate[0],
                           tree_data->TPDate[t])/365.0;
        CurrTimeInYrsSqr = CurrTimeInYrs * CurrTimeInYrs;
        Step              = 0.0;            /* initialise to smoothing off*/
        

        /* Store "Notional * CoP *(Underlying - Strike) - Opt" in AuxSlice */
        if (Hyb3_LCombTwoSlices(AuxSlice,
                           DMode,
                           Under,
                           (Notional * CoP),
                           Opt,
                           -1.0,
                           t,
                           tree_data) == FAILURE) goto RETURN;

        if (Hyb3_AddScalar (AuxSlice,
                       DMode,
                      (-1.0 * Notional * CoP * Strike),
                      t,
                      tree_data) == FAILURE) goto RETURN;

        switch (sliceDim)
        {
        case 1:
        
            offset    = Hyb3_Node_Offset(1,0,0,t,tree_data);
            OptL      = (double *)Opt + offset;
            AuxSliceL = AuxSlice      + offset;
            UnderL =    (double *)Under + offset;

            if (ExerProb    != NULL) ExerProbL    = ExerProb     + offset;
            if (ExerTime    != NULL) ExerTimeL    = ExerTime     + offset;
            if (ExerTimeSqr != NULL) ExerTimeSqrL = ExerTimeSqr + offset; 


            /* 
            *   If this is the last exercise date OptL has been initialized to
            *   zero so the exercise decision is Max (intrinsic value, 0).
            *   If this is not the last exercise date we have an American 
            *   exercise decision Max (intrinsic value, live option).
            */

            for (i = Bottom1; i <= Top1; i ++)
            {   

                if (SmoothingOn)
                {
                    Step = Hyb3_GetIndexStep(AuxSlice,1,i,0,0,t,tree_data);
                }
                
                if (ExerProb != NULL)
                    ExerProbL[i] = DrSmoothStep(1.0,
                                                ExerProbL[i],
                                                AuxSliceL[i],
                                                0.0,
                                                Step);

                if (ExerTime != NULL)
                    ExerTimeL[i] = DrSmoothStep(CurrTimeInYrs,
                                                ExerTimeL[i],
                                                AuxSliceL[i],
                                                0.0,
                                                Step);
                if (ExerTimeSqr != NULL)
                    ExerTimeSqrL[i] = DrSmoothStep(CurrTimeInYrsSqr,
                                                   ExerTimeSqrL[i],
                                                   AuxSliceL[i],
                                                   0.0,
                                                   Step);


                OptL[i] = DrSmoothMax(AuxSliceL[i],0.0,Step) + OptL[i];

            }  /* for i */	
            break;

        case 2:

            for (i = Bottom1; i <= Top1; i ++)
            {

                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                OptL   = (double *)Opt + offset;
                AuxSliceL = AuxSlice + offset;                

                if (ExerProb    != NULL) ExerProbL    = ExerProb     + offset;
                if (ExerTime    != NULL) ExerTimeL    = ExerTime     + offset;
                if (ExerTimeSqr != NULL) ExerTimeSqrL = ExerTimeSqr  + offset; 

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    if (SmoothingOn)
                    {
                        Step = Hyb3_GetIndexStep(AuxSlice,2,i,j,0,t,tree_data);
                    }

                    if (ExerProb != NULL)
                        ExerProbL[j] = DrSmoothStep(1.0,
                                                ExerProbL[j],
                                                AuxSliceL[j],
                                                0.0,
                                                Step);

                    if (ExerTime != NULL)
                        ExerTimeL[j] = DrSmoothStep(CurrTimeInYrs,
                                                ExerTimeL[j],
                                                AuxSliceL[j],
                                                0.0,
                                                Step);
                    if (ExerTimeSqr != NULL)
                        ExerTimeSqrL[j] = DrSmoothStep(CurrTimeInYrsSqr,
                                                   ExerTimeSqrL[j],
                                                   AuxSliceL[j],
                                                   0.0,
                                                   Step);


                    OptL[j] = DrSmoothMax(AuxSliceL[j],0.0,Step) + OptL[j];
                    

                }  /* for j */	
            }
            break;

        case 3:
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    OptL   = (double *)Opt + offset;
                    AuxSliceL = AuxSlice + offset;

                    if (ExerProb    != NULL) ExerProbL    = ExerProb   +offset;
                    if (ExerTime    != NULL) ExerTimeL    = ExerTime   +offset;
                    if (ExerTimeSqr != NULL) ExerTimeSqrL = ExerTimeSqr+offset;


                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        if (SmoothingOn)
                        {
                            Step = Hyb3_GetIndexStep(AuxSlice,3,i,j,k,t,tree_data);
                        }

                        if (ExerProb != NULL)
                            ExerProbL[k] = DrSmoothStep(1.0,
                                                        ExerProbL[k],
                                                        AuxSliceL[k],
                                                        0.0,
                                                        Step);

                        if (ExerTime != NULL)
                            ExerTimeL[k] = DrSmoothStep(CurrTimeInYrs,
                                                        ExerTimeL[k],
                                                        AuxSliceL[k],
                                                        0.0,
                                                        Step);
                        if (ExerTimeSqr != NULL)
                            ExerTimeSqrL[k] = DrSmoothStep(CurrTimeInYrsSqr,
                                                           ExerTimeSqrL[k],
                                                           AuxSliceL[k],
                                                           0.0,
                                                           Step);

                        OptL[k] = DrSmoothMax(AuxSliceL[k],0.0,Step) + OptL[k];

                    }  /* for k */
                }
            }
            break;

        default:
            DR_Error("Unknown Dev Mode (Hyb3_Option_t)! \n");
            break;
        }  /* switch slice dim*/
    }  /* if */    

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Option_t */



/*****  Hyb3_Exch_Option_t  *******************************************************/
/**
*       Exchange option.
*/
int     Hyb3_Exch_Option_t(TSLICE      Opt,        /**< (I/O) Option prices         */
                      TSLICE      Under1,     /**< (I) First underlying        */
                      TSLICE      Under2,     /**< (I) Second underlying       */
                      double      Notional1,  /**< (I) First notional          */
                      double      Notional2,  /**< (I) Second notional         */
                      double      Strike,     /**< (I) Strike                  */
                      long        ExerFlag,   /**< (I) Exercise flag           */
                      int         CoP,        /**< (I) =1 for call, -1 for put */
                      int         t,	      /**< (I) Current time period     */
                      int         T,	      /**< (I) Total number of period  */
                      int         DCurve,     /**< (I) Discount curve          */
                      int         DMode,      /**< (I) Dim of DEV and of slices*/
                      HYB3_DEV_DATA    *dev_data,  /**< (I) Hyb3_Dev data structure      */
                      HYB3_TREE_DATA   *tree_data) /**< (I) Tree data structure     */
{


    double    *OptL;
    double    *Under1L;
    double    *Under2L;
    
    int	   
                                 
        Top1,   Bottom1,   /* Tree limits (1rst dim) */
       *Top2,  *Bottom2,   /* Tree limits (2nd dim)  */
      **Top3, **Bottom3,   /* Tree limits (3rd dim)  */
        i, j, k,           /* Node indices           */
        offset,
        sliceDim,

        status = FAILURE;  /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];
    sliceDim = Hyb3_Slice_Dim( DMode);


    if (Hyb3_Dev(Opt,
            t,
            T,
            DCurve,
            DMode,
            dev_data,
            tree_data) == FAILURE)
    {
        goto RETURN;
                    
    }  /* if */


    if (ExerFlag)
    {
        switch (sliceDim)
        {
        case 1:
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            OptL    = (double *)Opt    + offset;
            Under1L = (double *)Under1 + offset;
            Under2L = (double *)Under2 + offset;


            for (i = Bottom1; i <= Top1; i ++)
            {
                OptL[i] = MAX (CoP * (Notional1 * Under1L[i] 
                                 - Notional2 * Under2L[i] - Strike), OptL[i]);

            }  /* for i */
            break;
        
        case 2:
        
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                OptL    = (double *)Opt    + offset;
                Under1L = (double *)Under1 + offset;
                Under2L = (double *)Under2 + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    OptL[j] = MAX (CoP * (Notional1 * Under1L[j] 
                           - Notional2 * Under2L[j] - Strike), OptL[j]);

                }  /* for j */
            }
        break;

        case 3:
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    OptL    = (double *)Opt    + offset;
                    Under1L = (double *)Under1 + offset;
                    Under2L = (double *)Under2 + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        OptL[k] = MAX (CoP*(Notional1*Under1L[k] 
                     - Notional2 * Under2L[k] - Strike), OptL[k]);

                    }  /* for k */
                }
            }
            break;

        default:
            DR_Error("Unknown Dev Mode (Hyb3_Option_t)! \n");
            break;
        } /* Switch slice Dim */
    }  /* if */    

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Exch_Option_t */



/*****  Hyb3_Option_Series_t  ****************************************************/
/**
*       Series of European options.
*/
int     Hyb3_Option_Series_t(TSLICE      Opt,       /**< (I/O) Option prices        */
                        TSLICE      Under,     /**< (I) Underlying prices      */
                        double      Notional,  /**< (I) Notional               */
                        double      Strike,    /**< (I) Strike 	             */
                        long        ExerFlag,  /**< (I) Exercise flag          */
                        int         CoP,       /**< (I) =1 for call, -1  put   */
                        int         t,	       /**< (I) Current time period    */
                        int         T,	       /**< (I) Total number of period */
                        int         DCurve,    /**< (I) Discount curve         */
                        int         DMode,     /**< (I) Dim of DEV and slices  */
                        HYB3_DEV_DATA    *dev_data, /**< (I) Hyb3_Dev data structure     */
                        HYB3_TREE_DATA   *tree_data)/**< (I) Tree data structure    */
{


    double    *OptL;
    double    *UnderL;

    int	 
                                   
         Top1,  Bottom1,     /* Tree limits (1rst dim) */
        *Top2, *Bottom2,     /* Tree limits (2nd dim)  */
       **Top3,**Bottom3,     /* Tree limits (3rd dim)  */
        i, j, k,             /* Node indices           */
        offset,
        sliceDim,

        status = FAILURE;	 /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];
    sliceDim = Hyb3_Slice_Dim( DMode);


    if (Hyb3_Dev(Opt,
            t,
            T,
            DCurve,
            DMode,
            dev_data,
            tree_data) == FAILURE)
    {
        goto RETURN;
                    
    }  /* if */


    if (ExerFlag)
    {
        switch (sliceDim)
        {
        case 1:
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            OptL   = (double *) Opt   + offset;
            UnderL = (double *) Under + offset;


            for (i = Bottom1; i <= Top1; i ++)                                      
            {
                OptL[i] += MAX (Notional * CoP * (UnderL[i] - Strike), 0.);

            }  /* for i */
            break;

        case 2:
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                OptL   = (double *) Opt   + offset;
                UnderL = (double *) Under + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    OptL[j] += MAX (Notional * CoP * (UnderL[j] - Strike), 0.);

                }  /* for j */	
            }
            break;

        case 3:
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    OptL   = (double *) Opt   + offset;
                    UnderL = (double *) Under + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        OptL[k] += MAX (Notional * CoP * (UnderL[k] - Strike), 0.);

                    }  /* for k */
                }
            }
            break;

        default:
            DR_Error("Unknown Dev Mode (Hyb3_Option_t)! \n");
            break;
        }  /* switch sliceDim */
    }  /* if */    

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Option_Series_t */



/*****  Hyb3_Option_Simple_t *****************************/
/**
*       Series of European options.
*/
int     Hyb3_Option_Simple_t(TSLICE      Opt,       /**< (I/O) Option prices        */                        
                        double      Notional,  /**< (I) Notional               */
                        double      Strike,    /**< (I) Strike 	             */                        
                        int         CoP,       /**< (I) =1 for call, -1  put   */
                        int         t,	       /**< (I) Current time period    */
                        int         DMode,     /**< (I) Dim of DEV and slices  */
                        HYB3_TREE_DATA   *tree_data)/**< (I) Tree data structure    */
{


    double    *OptL;
   

    int	 
                                   
         Top1,  Bottom1,     /* Tree limits (1rst dim) */
        *Top2, *Bottom2,     /* Tree limits (2nd dim)  */
       **Top3,**Bottom3,     /* Tree limits (3rd dim)  */
        i, j, k,             /* Node indices           */
        offset,
        sliceDim,

        status = FAILURE;	 /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];
    sliceDim = Hyb3_Slice_Dim( DMode);


    
    
    switch (sliceDim)
    {
    case 1:
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        OptL   = (double *) Opt   + offset;
        


        for (i = Bottom1; i <= Top1; i ++)                                      
        {
            OptL[i] = MAX (Notional * CoP * (OptL[i] - Strike), 0.);

        }  /* for i */
        break;
    
    case 2:
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            OptL   = (double *) Opt   + offset;
        

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                OptL[j] = MAX (Notional * CoP * (OptL[j] - Strike), 0.);

            }  /* for j */	
        }
        break;

    case 3:
        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                OptL   = (double *) Opt   + offset;
                

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                {
                    OptL[k] = MAX (Notional * CoP * (OptL[k] - Strike), 0.);

                }  /* for k */
            }
        }
        break;

    default:
        DR_Error("Unknown Dev Mode (Hyb3_Option_t)! \n");
        break;
    }  /* switch sliceDim */

    

    
    status = SUCCESS;

    return (status);

}  /* Hyb3_Option_Series_t */

