/****************************************************************************/
/*      Time slice manipulation.                                            */
/****************************************************************************/
/*      SLICE.c                                                             */
/****************************************************************************/


/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"


int     Fix3_ShowOffset ( 
                    int         t,              /**< (I) Current time point  */
                    FIX3_TREE_DATA const* tree_data) /**< (I) Tree data structure */
{

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j;                              /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);
            printf("t=%d off=%d\n", t, offset);
            break;
        }
        case 2:
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);
                printf("t=%d i=%d off=%d\n", t, i, offset);
            }  /* for i */

            break;
        }
        case 3:
        {        
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);
                    printf("t=%d i=%d j=%d off=%d\n", t, i, j, offset);
                }  /* for j */

            break;
        }
        default:
        {
            DR_Error("Fix3_ShowOffset: invalid slice dimension (must be 1,2 or 3)!");
            return(FAILURE);        
        }
    } /* End of switch() */

    return (SUCCESS);

}



/*****  Fix3_Set_Slice  *********************************************************/
/**
*       Set a time slice to a given value 
*/
int     Fix3_Set_Slice ( 
                    double      *Slice,         /**< (O) Slice prices        */
                    double      sliceValue,     /**< (I) Value to be set to  */
                    int         t,              /**< (I) Current time point  */
                    FIX3_TREE_DATA const* tree_data) /**< (I) Tree data structure */
{
    double  *SliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;
    
            for (i = Bottom1; i <= Top1; i ++)
            {
                SliceL[i] = sliceValue;
            }

            break;
        }
        case 2:
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    SliceL[j] = sliceValue;
                }
            }  /* for i */

            break;
        }
        case 3:
        {        
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        SliceL[k] = sliceValue;
                    }
                }  /* for j */

            break;
        }
        default:
        {
            DR_Error("Fix3_Set_Slice: invalid slice dimension (must be 1,2 or 3)!");
            return(FAILURE);        
        }
    } /* End of switch() */

    return (SUCCESS);

}  /* Fix3_Set_Slice */



/*****  Fix3_Copy_Slice  *********************************************************/
/**
*       Copy a time slice from a given time slice 
*/
int     Fix3_Copy_Slice(
                   double      *CopySlice,   /**< (O) Slice prices            */
                   double      *OrgSlice,    /**< (I) Original to be repli-ed */
                   int         t,            /**< (I) Current time point      */
                   FIX3_TREE_DATA const* tree_data) /**< (I) Tree data structure     */
{

    double  *CopySliceL;
    double  *OrgSliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            CopySliceL = CopySlice + offset;
            OrgSliceL  = OrgSlice  + offset;
    
            for (i = Bottom1; i <= Top1; i ++)
            {
                CopySliceL[i] = OrgSliceL[i];
            }

            break;
        }
        case 2:
        {        
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                CopySliceL = CopySlice + offset;
                OrgSliceL  = OrgSlice  + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    CopySliceL[j] = OrgSliceL[j];
                }
            }  /* for i */

            break;
        }
        case 3:
        {        
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    CopySliceL = CopySlice + offset;
                    OrgSliceL  = OrgSlice  + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        CopySliceL[k] = OrgSliceL[k];
                    }    
                }  /* for j */

            break;
        }
        default:
        {
            DR_Error("Fix3_Copy_Slice: invalid slice dimension "
                        "(must be 1,2 or 3)!");
            return(FAILURE);        
        }
    } /* End of switch() */
    
    return (SUCCESS);

}  /* Fix3_Copy_Slice */


/*****   Fix3_LCombTwoSlices  ****************************************************/
/**
 *  Linear combination of two slices onto a CombSlice. All five slice must 
 *  be of same dimension and memory for SumSlice must be pre-allocated.
 *
 *    CombSlice = a1 * Slice1 + a2 * Slice2
 *
 */
int Fix3_LCombTwoSlices(
                   double        *CombSlice,  /**< (O) Slice to be modified   */
                   double        *Slice1,     /**< (I) First one to add       */
                   double         a1,         /**< (I) Coeff of first slice   */
                   double        *Slice2,     /**< (I) Second one to add      */
                   double         a2,         /**< (I) Coeff of second slice  */
                   int            t,          /**< (I) Current time point     */
                   FIX3_TREE_DATA const*tree_data) /**< (I) Tree data         */
{

    double  *CombSliceL;
    double  *Slice1L;
    double  *Slice2L;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            CombSliceL = CombSlice + offset;
            Slice1L    = Slice1    + offset;
            Slice2L    = Slice2    + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                CombSliceL[i] = a1*Slice1L[i] + a2*Slice2L[i];
            }

            break;
        }
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                CombSliceL = CombSlice + offset;
                Slice1L    = Slice1    + offset;
                Slice2L    = Slice2    + offset;
                
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    CombSliceL[j] = a1*Slice1L[j] + a2*Slice2L[j];
                }
            }

            break;
        }
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    CombSliceL = CombSlice + offset;
                    Slice1L    = Slice1    + offset;
                    Slice2L    = Slice2    + offset;
                    
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        CombSliceL[k] =a1*Slice1L[k] + a2*Slice2L[k];
                    }
                }
            }

            break;
        }
        default:
        {
            DR_Error("Fix3_LCombTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  Fix3_LCombTwoSlices  */



/*****   Fix3_AddTwoSlices  *******************************************************/
/**
 *  Add two slices into a SumSlice. All slice must be of same
 *  dimension and memory for SumSlice must be pre-allocated.
 *
 */
int Fix3_AddTwoSlices(
                 double        *Slice,        /**< (O) Slice to be modified */
                 double        *Slice1,       /**< (I) First one to add     */
                 double        *Slice2,       /**< (I) Second one to add    */
                 int            t,            /**< (I) Current time point   */
                 FIX3_TREE_DATA const*tree_data) /**< (I) Tree data         */
{

    double  *SliceL;
    double  *Slice1L;
    double  *Slice2L;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL  = Slice  + offset;
            Slice1L = Slice1 + offset;
            Slice2L = Slice2 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = Slice1L[i] + Slice2L[i];
            }

            break;
        }
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL  = Slice  + offset;
                Slice1L = Slice1 + offset;
                Slice2L = Slice2 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = Slice1L[j] + Slice2L[j];
                }
            }

            break;
        }
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL  = Slice  + offset;
                    Slice1L = Slice1 + offset;
                    Slice2L = Slice2 + offset;
        
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = Slice1L[k] + Slice2L[k];
                    }
                }
            }

            break;
        }
        default:
        {
            DR_Error("Fix3_AddTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  Fix3_AddTwoSlices  */



/*****   Fix3_MultiplyTwoSlices  *************************************************/
/**
 *  Multipy two slices. All slices must be of same dimension and memory
 *  must be pre-allocated.
 *
 */
int Fix3_MultiplyTwoSlices(
                      double      *Slice,     /**< (O) Slice to be modified */
                      double      *Slice1,    /*<* (I) First one to add     */
                      double      *Slice2,    /*<* (I) Second one to add    */
                      int         t,          /*<* (I) Current time point   */
                      FIX3_TREE_DATA const*tree_data) /*<* (I) Tree data    */
{

    double  *SliceL;
    double  *Slice1L;
    double  *Slice2L;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL  = Slice  + offset;
            Slice1L = Slice1 + offset;
            Slice2L = Slice2 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = Slice1L[i] * Slice2L[i];
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL  = Slice  + offset;
                Slice1L = Slice1 + offset;
                Slice2L = Slice2 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = Slice1L[j] * Slice2L[j];
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL  = Slice  + offset;
                    Slice1L = Slice1 + offset;
                    Slice2L = Slice2 + offset;
        
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = Slice1L[k] * Slice2L[k];
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_MultiplyTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  Fix3_MultiplyTwoSlices  */


/*****   Fix3_MultTwoSlicesAddAll  *********************************************/
/**
 *  Multipy two slices and perform addition of all nodes.  This operation
 *  is meaningful in the context of the "express DEV" facility.
 *  All slices must be of same dimension and memorym ust be pre-allocated.
 *
 */
int Fix3_MultTwoSlicesAddAll(
                        double      *Value,     /**< (O) Return value         */
                        double      *Slice1,    /**< (I) First one to add     */
                        double      *Slice2,    /**< (I) Second one to add    */
                        int         t,          /**< (I) Current time point   */
                        FIX3_TREE_DATA const*tree_data) /**< (I) Tree data            */
{
    double  *Slice1L;
    double  *Slice2L;

    double   ValueL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  

    ValueL = 0.0;

    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            
            Slice1L = Slice1 + offset;
            Slice2L = Slice2 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                ValueL += Slice1L[i] * Slice2L[i];
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                Slice1L = Slice1 + offset;
                Slice2L = Slice2 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    ValueL += Slice1L[j] * Slice2L[j];
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    Slice1L = Slice1 + offset;
                    Slice2L = Slice2 + offset;
        
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        ValueL += Slice1L[k] * Slice2L[k];
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_MultiplyTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */


    /* Output value */
    *Value = ValueL;

    return (SUCCESS);

}  /*  Fix3_MultTwoSlicesAddAll  */




/*****  Fix3_AddScalar  ***********************************************************/
/**
 *  Add a scalar to a slice of dimension up to three.
 */
int    Fix3_AddScalar(
                 double        *Slice,         /**< (O) Slice to be modified */
                 double         Scalar,        /**< (I) Scalar to be added   */
                 int            t,             /**< (I) Current time point   */
                 FIX3_TREE_DATA const*tree_data) /* (I) Tree data            */
{
    double  *SliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  
    
    
    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] += Scalar;
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] += Scalar;
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] += Scalar;
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_AddScalar: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */
    
    return (SUCCESS);
    
}  /*  Fix3_AddScalar  */



/*****  Fix3_MultiplyScalar  ******************************************************/
/**
 *  Multiply a slice of dimension up to three by a scalar.
 */
int    Fix3_MultiplyScalar(
                      double       *Slice,     /**< (O) Slice to be modified */
                      double        Scalar,    /**< (I) Scalar to be mult    */
                      int           t,         /**< (I) Current time point   */
                      FIX3_TREE_DATA const* tree_data) /**< (I) Tree data    */
{
    double  *SliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  
    

    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] *= Scalar;
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] *= Scalar;
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] *= Scalar;
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_MultiplyScalar: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  Fix3_MultiplyScalar  */



/*****  Fix3_AccrueSlice  ******************************************************/
/**
 *  Multiply a slice of dimension up to three by a scalar.
 */
int    Fix3_AccrueSlice(
                   double       *Slice,     /**< (O) Slice to be modified        */
                   double        DayCount,  /**< (I) Scalar to be mult           */
                   int           isSimple,  /**< (I) If true, simple compounding */
                   int           t,         /**< (I) Current time point          */
                   FIX3_TREE_DATA const*tree_data) /**< (I) Tree data            */
{
    double  *SliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  
    

    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = ACC_FN(SliceL[i], DayCount, isSimple);
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
	                SliceL[j] = ACC_FN(SliceL[j], DayCount, isSimple);
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
		                SliceL[k] = ACC_FN(SliceL[k], DayCount, isSimple);
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_AccrueSlice: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  Fix3_MultiplyScalar  */



/*****  Fix3_MaxMinOnSlice  ******************************************************/
/**
 *  Apply a maximum and a minimum to the values in a slice.
 */
int    Fix3_MaxMinOnSlice(
                     double        *Slice,     /**< (O) Slice to be modified */
                     double         Min,       /**< (I) Minimum level        */
                     double         Max,       /**< (I) Maximum level        */
                     int            t,         /**< (I) Current time point   */
                     FIX3_TREE_DATA const*tree_data) /**< (I) Tree data      */
{
    double  *SliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  
    

    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MAX(MIN(SliceL[i],Max),Min);
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = MAX(MIN(SliceL[j],Max),Min);
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MAX(MIN(SliceL[k],Max),Min);
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_MaxMinOnSlice: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }                    
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  Fix3_MaxMinOnSlice  */



/*****  Fix3_MaxOnSlice  *********************************************************/
/**
 *  Apply a maximum to the values in a slice.
 */
int    Fix3_MaxOnSlice(
                  double        *Slice,        /**< (O) Slice to be modified */
                  double         Max,          /**< (I) Maximum level        */
                  int            t,            /**< (I) Current time period  */
                  FIX3_TREE_DATA const*tree_data) /**< (I) Tree data         */
{
    double  *SliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  
    

    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MIN(SliceL[i],Max);
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = MIN(SliceL[j],Max);
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MIN(SliceL[k],Max);
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_MaxOnSlice: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  Fix3_MaxOnSlice  */


/*****  Fix3_MinOnSlice  *********************************************************/
/**
 *  Apply a minimum to the values in a slice.
 */
int    Fix3_MinOnSlice(
                  double        *Slice,        /**< (O) Slice to be modified */
                  double         Min,          /**< (I) Minimum level        */
                  int            t,            /**< (I) Current time period  */
                  FIX3_TREE_DATA const*tree_data) /**< (I) Tree data         */
{
     return Fix3_SmoothMinOnSlice(
                        Slice,
                        Min,
                        FALSE,
                        t,
                        tree_data);
}

/*****  Fix3_SmoothMinOnSlice  *********************************************************/
/**
 *  Apply a minimum to the values in a slice.
 */
int    Fix3_SmoothMinOnSlice(
                  double        *Slice,        /**< (O) Slice to be modified */
                  double         Min,          /**< (I) Minimum level        */
		  int            SmoothingOn,  /**< (I) TRUE = smoothing on  */
                  int            t,            /**< (I) Current time period  */
                  FIX3_TREE_DATA const*tree_data) /**< (I) Tree data         */
{
    double  *SliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim)  */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)   */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)   */

    int     i, j, k;                           /* Node indices            */
    int     offset;                            /* Node offset             */
    double  Step;                              /* step size for smoothing */

    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  
    Step = 0.0;  /* initialise to smoothing off */
    
    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                if (SmoothingOn)
                {
                    Step = Fix3_GetIndexStep(Slice,1,i,0,0,t,tree_data);
                }

                SliceL[i] = DrSmoothMax(SliceL[i],Min,Step);
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    if (SmoothingOn)
                    {
                        Step = Fix3_GetIndexStep(Slice,2,i,j,0,t,tree_data);
                    }

                    SliceL[j] = DrSmoothMax(SliceL[j],Min,Step);
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        if (SmoothingOn)
                        {
                            Step = Fix3_GetIndexStep(Slice,3,i,j,k,t,tree_data);
                        }

                        SliceL[k] = DrSmoothMax(SliceL[k],Min,Step);
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_SmoothMinOnSlice: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  Fix3_SmoothMinOnSlice  */



/*****  Fix3_StepUp  *********************************************************/
/**
 *  In DR spirit - I just added another useless 100 lines of code 
 */
int    Fix3_StepUp(
                  double*        Slice,           /**< (O) Slice to be modified */
                  double const*  Index,           /**< (I) Reference slice      */
                  double         Up,              /**< (I) Maximium level       */
                  double         Barrier,         /**< (I) Barrier level        */
		          int            SmoothingOn,     /**< (I) TRUE = smoothing on  */
                  int            t,               /**< (I) Current time period  */
                  FIX3_TREE_DATA const*tree_data) /**< (I) Tree data            */
{
    double  *SliceL;
    double const*IndexL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim)  */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)   */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)   */

    int     i, j, k;                           /* Node indices            */
    int     offset;                            /* Node offset             */
    double  Step;                              /* step size for smoothing */

    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  
    Step = 0.0;  /* initialise to smoothing off */
    
    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;
            IndexL = Index + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                if (SmoothingOn)
                    Step = Fix3_GetIndexStep(Index,1,i,0,0,t,tree_data);
                SliceL[i] = DrSmoothStep(Up, SliceL[i], IndexL[i], Barrier, Step);
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;
                IndexL = Index + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    if (SmoothingOn)
                        Step = Fix3_GetIndexStep(Index,2,i,j,0,t,tree_data);
                    SliceL[j] = DrSmoothStep(Up, SliceL[j], IndexL[j], Barrier, Step);
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;
                    IndexL = Index + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        if (SmoothingOn)
                            Step = Fix3_GetIndexStep(Index,3,i,j,k,t,tree_data);
                        SliceL[k] = DrSmoothStep(Up, SliceL[k], IndexL[k], Barrier, Step);
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_StepUp: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }
    } /* End of switch() */
    
    return (SUCCESS);

}



/*****  Fix3_Node_Offset  ********************************************************/
/**
 *      Calculate array offsets. There are three components:
 *
 *      1)  node indices run into [-halfwidth, halfwidth] in each dimension 
 *          whereas memory is allocated as [0, 2 * halfwidth]. This gives an 
 *          offset of halfwidth.
 *      2)  memory is allocated linearly in dimension 2 and 3. This gives an
 *          offset of i * width[1] in dimension 2 and i * width[1] * width[2]
 *          + j * width[1] in dimension 2.
 *      3)  the axis of the ellipse has a slope. The coordinate of the axis
 *          is given by (bottom2[i]+top2[i])/2 at index [i] in dimension 2 and
 *          (bottom3[i][j]+top3[i][j])/2 at index [i][j] in dimension 3.
 */
int     Fix3_Node_Offset (   
                        int         Dim,        /**< Dimension          */
                        int         i,          /**< Node indices       */
                        int         j,          /**<                    */
                        int         t,          /**< Current time point */
                        FIX3_TREE_DATA const*tree_data) /**< Tree data  */
{
    int     Offset;       /* Returned offset */

    if (Dim == 3)
    {
        Offset = (i + tree_data->HalfWidth[0]) 
                    * tree_data->Width[1] * tree_data->Width[2]
            +(j-((tree_data->OutTop2[t][i]+tree_data->OutBottom2[t][i])>>1)
                    + tree_data->HalfWidth[1]) * tree_data->Width[2]
         -((tree_data->OutTop3[t][i][j]+tree_data->OutBottom3[t][i][j])>>1)
              + tree_data->HalfWidth[2];
    }
    else if (Dim == 2)
    {
        Offset = (i + tree_data->HalfWidth[0]) * tree_data->Width[1]
               -((tree_data->OutTop2[t][i]+tree_data->OutBottom2[t][i])>>1)
                + tree_data->HalfWidth[1];
    }
    else
    {
        Offset = tree_data->HalfWidth[0];
    }

    return (Offset);

}  /*  Fix3_Node_Offset  */



/*****  Fix3_Init_Slice  ********************************************************/
/**
*       Initalise the whole slice to a given value
*/
int  Fix3_Init_Slice (
                 double      *Slice,     /**< (I/O) slice to be init'd  */
                 double       InitValue, /**< (I)   init value          */
                 FIX3_TREE_DATA const*tree_data) /* (I)   Tree data structure */
{
    int   status = FAILURE;

    long  ArraySize = 0L;
    long  i;

    if ((Slice == NULL) || (tree_data == NULL)) goto RETURN;

    /* calc the size of the slice array */
    switch (tree_data->NbFactor)
    {
    case 1:
        ArraySize = tree_data->Width[0];
        break;
    case 2:
        ArraySize = tree_data->Width[0] * 
                    tree_data->Width[1];
        break;
    case 3:
        ArraySize = tree_data->Width[0] * 
                    tree_data->Width[1] * 
                    tree_data->Width[2];
        break;
    default:
        goto RETURN;
    }

    /* set the value */
    for (i=0L; i<ArraySize; i++) Slice[i] = InitValue;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("Fix3_Init_Slice: failed.");
    }
    return(status);

}  /* Fix3_Init_Slice */



/*****   Fix3_MaxOfTwoSlices  ******************************************************/
/**
 * Take the max of two Slice
 */
int Fix3_MaxOfTwoSlices(
                   double        *Slice,        /**< (O) Slice to be modified */
                   double        *Slice1,       /**< (I) First one to add     */
                   double        *Slice2,       /**< (I) Second one to add    */
                   int            t,            /**< (I) Current time point   */
                   FIX3_TREE_DATA const*tree_data) /**< (I) Tree data         */
{

    double  *SliceL;
    double  *Slice1L;
    double  *Slice2L;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL  = Slice  + offset;
            Slice1L = Slice1 + offset;
            Slice2L = Slice2 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MAX(Slice1L[i], Slice2L[i]);
            }

            break;
        }
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL  = Slice  + offset;
                Slice1L = Slice1 + offset;
                Slice2L = Slice2 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = MAX(Slice1L[j], Slice2L[j]);
                }
            }

            break;
        }
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL  = Slice  + offset;
                    Slice1L = Slice1 + offset;
                    Slice2L = Slice2 + offset;
        
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MAX(Slice1L[k], Slice2L[k]);
                    }
                }
            }

            break;
        }
        default:
        {
            DR_Error("Fix3_AddTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  Fix3_MaxOfTwoSlices  */


/*****   Fix3_SmoothMaxOfTwoSlices  ******************************************************/
/**
 * Take the smooth max between two slices using DrSmoothMax
 */
int Fix3_SmoothMaxOfTwoSlices(
                   double        *Slice,        /**< (O) Slice to be modified */
                   double        *Slice1,       /**< (I) First one to add     */
                   double        *Slice2,       /**< (I) Second one to add    */
		   int            SmoothingOn,  /**< (I) TRUE = smoothing on  */
                   int            t,            /**< (I) Current time point   */
                   FIX3_TREE_DATA const*tree_data)  /**< (I) Tree data        */
{

    double  *SliceL;
    double  *Slice1L;
    double  *Slice2L;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim)  */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)   */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)   */

    int     i, j, k;                           /* Node indices            */
    int     offset;                            /* Node offset             */
    double  Step1, Step2, Step;                /* step size for smoothing */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  
    Step    = 0.0;                    /* initialise to smoothing off */

    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL  = Slice  + offset;
            Slice1L = Slice1 + offset;
            Slice2L = Slice2 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                if (SmoothingOn)
                {
                    Step1 = Fix3_GetIndexStep(Slice1,1,i,0,0,t,tree_data);
                    Step2 = Fix3_GetIndexStep(Slice2,1,i,0,0,t,tree_data);
                    Step  = MAX(Step1,Step2);
                }

                SliceL[i] = DrSmoothMax(Slice1L[i], Slice2L[i],Step);
            }

            break;
        }
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL  = Slice  + offset;
                Slice1L = Slice1 + offset;
                Slice2L = Slice2 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    if (SmoothingOn)
                    {
                        Step1 = Fix3_GetIndexStep(Slice1,2,i,j,0,t,tree_data);
                        Step2 = Fix3_GetIndexStep(Slice2,2,i,j,0,t,tree_data);
                        Step  = MAX(Step1,Step2);
                    }

                    SliceL[j] = DrSmoothMax(Slice1L[j], Slice2L[j],Step);
                }
            }

            break;
        }
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL  = Slice  + offset;
                    Slice1L = Slice1 + offset;
                    Slice2L = Slice2 + offset;
        
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        if (SmoothingOn)
                        {
                            Step1 = Fix3_GetIndexStep(Slice1,3,i,j,k,t,tree_data);
                            Step2 = Fix3_GetIndexStep(Slice2,3,i,j,k,t,tree_data);
                            Step  = MAX(Step1,Step2);
                        }

                        SliceL[k] = DrSmoothMax(Slice1L[k], Slice2L[k],Step);
                    }
                }
            }

            break;
        }
        default:
        {
            DR_Error("Fix3_AddTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  Fix3_SmoothMaxOfTwoSlices  */


/*****   Fix3_MinOfTwoSlices  ******************************************************/
/*
 */
int Fix3_MinOfTwoSlices(
                   double        *Slice,        /**< (O) Slice to be modified */
                   double        *Slice1,       /**< (I) First one to add     */
                   double        *Slice2,       /**< (I) Second one to add    */
                   int            t,            /**< (I) Current time point   */
                   FIX3_TREE_DATA const*tree_data) /**< (I) Tree data         */
{
    double  *SliceL;
    double  *Slice1L;
    double  *Slice2L;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL  = Slice  + offset;
            Slice1L = Slice1 + offset;
            Slice2L = Slice2 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MIN(Slice1L[i], Slice2L[i]);
            }

            break;
        }
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL  = Slice  + offset;
                Slice1L = Slice1 + offset;
                Slice2L = Slice2 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = MIN(Slice1L[j], Slice2L[j]);
                }
            }

            break;
        }
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL  = Slice  + offset;
                    Slice1L = Slice1 + offset;
                    Slice2L = Slice2 + offset;
        
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MIN(Slice1L[k], Slice2L[k]);
                    }
                }
            }

            break;
        }
        default:
        {
            DR_Error("Fix3_AddTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  Fix3_MinOfTwoSlices  */


/*****   Fix3_SmoothMinOfTwoSlices  ******************************************************/
/**
 * Take the smooth min between two slices using DrSmoothMax
 */
int Fix3_SmoothMinOfTwoSlices(
                   double        *Slice,        /**< (O) Slice to be modified */
                   double        *Slice1,       /**< (I) First one to add     */
                   double        *Slice2,       /**< (I) Second one to add    */
		   int            SmoothingOn,  /**< (I) TRUE = smoothing on  */
                   int            t,            /**< (I) Current time point   */
                   FIX3_TREE_DATA const*tree_data)  /**< (I) Tree data        */
{

    double  *SliceL;
    double  *Slice1L;
    double  *Slice2L;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim)  */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)   */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)   */

    int     i, j, k;                           /* Node indices            */
    int     offset;                            /* Node offset             */
    double  Step1, Step2, Step;                /* step size for smoothing */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  
    Step    = 0.0;                    /* initialise to smoothing off */

    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL  = Slice  + offset;
            Slice1L = Slice1 + offset;
            Slice2L = Slice2 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                if (SmoothingOn)
                {
                    Step1 = Fix3_GetIndexStep(Slice1,1,i,0,0,t,tree_data);
                    Step2 = Fix3_GetIndexStep(Slice2,1,i,0,0,t,tree_data);
                    Step  = MAX(Step1,Step2);
                }

                SliceL[i] = -1*DrSmoothMax(-1*Slice1L[i], -1*Slice2L[i],Step);
            }

            break;
        }
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL  = Slice  + offset;
                Slice1L = Slice1 + offset;
                Slice2L = Slice2 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    if (SmoothingOn)
                    {
                        Step1 = Fix3_GetIndexStep(Slice1,2,i,j,0,t,tree_data);
                        Step2 = Fix3_GetIndexStep(Slice2,2,i,j,0,t,tree_data);
                        Step  = MAX(Step1,Step2);
                    }

                    SliceL[j] = -1*DrSmoothMax(-1*Slice1L[j], -1*Slice2L[j],Step);
                }
            }

            break;
        }
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL  = Slice  + offset;
                    Slice1L = Slice1 + offset;
                    Slice2L = Slice2 + offset;
        
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        if (SmoothingOn)
                        {
                            Step1 = Fix3_GetIndexStep(Slice1,3,i,j,k,t,tree_data);
                            Step2 = Fix3_GetIndexStep(Slice2,3,i,j,k,t,tree_data);
                            Step  = MAX(Step1,Step2);
                        }

                        SliceL[k] = -1*DrSmoothMax(-1*Slice1L[k], -1*Slice2L[k],Step);
                    }
                }
            }

            break;
        }
        default:
        {
            DR_Error("Fix3_AddTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

} /*  Fix3_SmoothMinOfTwoSlices  */


/*****   Fix3_DivideTwoSlices  *************************************************/
/**
 *  Divide two slices. All slices must be of same dimension and memory
 *  must be pre-allocated. The caller must check that the second slice
 *  is non-zero.
 */
int Fix3_DivideTwoSlices(
                    double      *Slice,     /**< (O) Slice to be modified */
                    double      *Slice1,    /**< (I) Numerator            */
                    double      *Slice2,    /**< (I) Denominator          */
                    int         t,          /**< (I) Current time point   */
                    FIX3_TREE_DATA const*tree_data) /**< (I) Tree data    */
{

    double  *SliceL;
    double  *Slice1L;
    double  *Slice2L;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL  = Slice  + offset;
            Slice1L = Slice1 + offset;
            Slice2L = Slice2 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = Slice1L[i] / Slice2L[i];
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL  = Slice  + offset;
                Slice1L = Slice1 + offset;
                Slice2L = Slice2 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = Slice1L[j] / Slice2L[j];
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL  = Slice  + offset;
                    Slice1L = Slice1 + offset;
                    Slice2L = Slice2 + offset;
        
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = Slice1L[k] / Slice2L[k];
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Fix3_MultiplyTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  Fix3_DivideTwoSlices  */



/*****   SliceTimesScalar  *************************************************/
/*
 *  Multiply slice by scalar and put in another slice
 *
 */
int Fix3_SliceTimesScalar (double      *Slice, /* (O) Slice to be modified */
                           double      *Slice1,/* (I) Input slice          */
                           double      scalar, /* (I) Factor               */
                           int         t,      /* (I) Current time point   */
                           FIX3_TREE_DATA const* tree_data) /* Tree data   */
{

    double  *SliceL;
    double  *Slice1L;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];                 
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];                 
    Bottom3 = tree_data->Bottom3[t];                  


    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL  = Slice  + offset;
            Slice1L = Slice1 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = Slice1L[i] * scalar;
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL  = Slice  + offset;
                Slice1L = Slice1 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = Slice1L[j] * scalar;
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL  = Slice  + offset;
                    Slice1L = Slice1 + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = Slice1L[k] * scalar;
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("SliceTimesScalar: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /* Fix3_SliceTimesScalar */



int     Fix3_Dump_Slice (double const* Slice,         /**< (O) Slice prices        */
                         int         t,              /**< (I) Current time point  */
                         FIX3_TREE_DATA const* tree_data) /**< (I) Tree data structure */
{
    double const*SliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */

    int     i, j, k;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


printf("t=%d d=%ld\n", t, tree_data->TPDate[t]);
    switch (tree_data->NbFactor)
    {
        case 1:
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;
    
            for (i = Bottom1; i <= Top1; i ++)
            {
                printf("%.6f ", SliceL[i]);
            }
printf("\n");  

            break;
        }
        case 2:
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    printf("%.6f ", SliceL[j]);
                }
printf("\n");  
            }  /* for i */
printf("\n");  

            break;
        }
        case 3:
        {        
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        printf("%.6f ", SliceL[k]);
                    }
printf("\n");  
                }  /* for j */
printf("\n");  
            }
printf("\n");  

            break;
        }
        default:
        {
            DR_Error("Fix3_Dump_Slice: invalid slice dimension (must be 1,2 or 3)!");
            return(FAILURE);        
        }
    } /* End of switch() */

    return (SUCCESS);
}



