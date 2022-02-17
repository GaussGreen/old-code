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
#include "tmx123head.h"



/*****  Set_Slice  *********************************************************/
/*
*       Set a time slice to a given value 
*/
int     Set_Slice ( double      *Slice,         /* (O) Slice prices        */
                    double      sliceValue,     /* (I) Value to be set to  */
                    int         t,              /* (I) Current time point  */
                    TREE_DATA   *tree_data)     /* (I) Tree data structure */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("Set_Slice: invalid slice dimension (must be 1,2 or 3)!");
            return(FAILURE);        
        }
    } /* End of switch() */

    return (SUCCESS);

}  /* Set_Slice */



/*****  Copy_Slice  *********************************************************/
/*
*       Copy a time slice from a given time slice 
*/
int     Copy_Slice(double      *CopySlice,   /* (O) Slice prices            */
                   double      *OrgSlice,    /* (I) Original to be repli-ed */
                   int         t,            /* (I) Current time point      */
                   TREE_DATA   *tree_data)   /* (I) Tree data structure     */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("Copy_Slice: invalid slice dimension "
                        "(must be 1,2 or 3)!");
            return(FAILURE);        
        }
    } /* End of switch() */
    
    return (SUCCESS);

}  /* Copy_Slice */


/*****   LCombTwoSlices  ****************************************************/
/*
 *  Linear combination of two slices onto a CombSlice. All five slice must 
 *  be of same dimension and memory for SumSlice must be pre-allocated.
 *
 *    CombSlice = a1 * Slice1 + a2 * Slice2
 *
 */
int LCombTwoSlices(double        *CombSlice,  /* (O) Slice to be modified   */
                   double        *Slice1,     /* (I) First one to add       */
                   double         a1,         /* (I) Coeff of first slice   */
                   double        *Slice2,     /* (I) Second one to add      */
                   double         a2,         /* (I) Coeff of second slice  */
                   int            t,          /* (I) Current time point     */
                   TREE_DATA     *tree_data)  /* (I) Tree data              */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("LCombTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  LCombTwoSlices  */



/*****   AddTwoSlices  *******************************************************/
/*
 *  Add two slices into a SumSlice. All slice must be of same
 *  dimension and memory for SumSlice must be pre-allocated.
 *
 */
int AddTwoSlices(double        *Slice,        /* (O) Slice to be modified */
                 double        *Slice1,       /* (I) First one to add     */
                 double        *Slice2,       /* (I) Second one to add    */
                 int            t,            /* (I) Current time point   */
                 TREE_DATA     *tree_data)    /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("AddTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  AddTwoSlices  */



/*****   MultiplyTwoSlices  *************************************************/
/*
 *  Multipy two slices. All slices must be of same dimension and memory
 *  must be pre-allocated.
 *
 */
int MultiplyTwoSlices(double      *Slice,     /* (O) Slice to be modified */
                      double      *Slice1,    /* (I) First one to add     */
                      double      *Slice2,    /* (I) Second one to add    */
                      int         t,          /* (I) Current time point   */
                      TREE_DATA   *tree_data) /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("MultiplyTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  MultiplyTwoSlices  */


/*****   MultTwoSlicesAddAll  *********************************************/
/*
 *  Multipy two slices and perform addition of all nodes.  This operation
 *  is meaningful in the context of the "express DEV" facility.
 *  All slices must be of same dimension and memorym ust be pre-allocated.
 *
 */
int MultTwoSlicesAddAll(double      *Value,     /* (O) Return value         */
                        double      *Slice1,    /* (I) First one to add     */
                        double      *Slice2,    /* (I) Second one to add    */
                        int         t,          /* (I) Current time point   */
                        TREE_DATA   *tree_data) /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

            
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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("MultiplyTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */


    /* Output value */
    *Value = ValueL;

    return (SUCCESS);

}  /*  MultTwoSlicesAddAll  */




/*****  AddScalar  ***********************************************************/
/*
 *  Add a scalar to a slice of dimension up to three.
 */
int    AddScalar(double        *Slice,         /* (O) Slice to be modified */
                 double         Scalar,        /* (I) Scalar to be added   */
                 int            t,             /* (I) Current time point   */
                 TREE_DATA     *tree_data)     /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("AddScalar: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */
    
    return (SUCCESS);
    
}  /*  AddScalar  */



/*****  MultiplyScalar  ******************************************************/
/*
 *  Multiply a slice of dimension up to three by a scalar.
 */
int    MultiplyScalar(double       *Slice,     /* (O) Slice to be modified */
                      double        Scalar,    /* (I) Scalar to be mult    */
                      int           t,         /* (I) Current time point   */
                      TREE_DATA    *tree_data) /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("MultiplyScalar: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  MultiplyScalar  */



/*****  AccrueSlice  ******************************************************/
/*
 *  Multiply a slice of dimension up to three by a scalar.
 */
int    AccrueSlice(double       *Slice,     /* (O) Slice to be modified        */
                   double        DayCount,  /* (I) Scalar to be mult           */
				   int			 isSimple,  /* (I) If true, simple compounding */
                   int           t,         /* (I) Current time point          */
                   TREE_DATA    *tree_data) /* (I) Tree data                   */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("AccrueSlice: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  MultiplyScalar  */



/*****  MaxMinOnSlice  ******************************************************/
/*
 *  Apply a maximum and a minimum to the values in a slice.
 */
int    MaxMinOnSlice(double        *Slice,     /* (O) Slice to be modified */
                     double         Min,       /* (I) Minimum level        */
                     double         Max,       /* (I) Maximum level        */
                     int            t,         /* (I) Current time point   */
                     TREE_DATA     *tree_data) /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("MaxMinOnSlice: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }                    
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  MaxMinOnSlice  */



/*****  MaxOnSlice  *********************************************************/
/*
 *  Apply a maximum to the values in a slice.
 */
int    MaxOnSlice(double        *Slice,        /* (O) Slice to be modified */
                  double         Max,          /* (I) Maximum level        */
                  int            t,            /* (I) Current time period  */
                  TREE_DATA     *tree_data)    /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("MaxOnSlice: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  MaxOnSlice  */



/*****  MinOnSlice  *********************************************************/
/*
 *  Apply a minimum to the values in a slice.
 */
int    MinOnSlice(double        *Slice,        /* (O) Slice to be modified */
                  double         Min,          /* (I) Minimum level        */
                  int            t,            /* (I) Current time period  */
                  TREE_DATA     *tree_data)    /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MAX(SliceL[i],Min);
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = MAX(SliceL[j],Min);
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
                    offset = Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MAX(SliceL[k],Min);
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("MinOnSlice: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */
    
    return (SUCCESS);

}  /*  MinOnSlice  */



/*****  Node_Offset  ********************************************************/
/*
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
int     Node_Offset (   int         Dim,        /* Dimension          */
                        int         i,          /* Node indices       */
                        int         j,          /*                    */
                        int         t,          /* Current time point */
                        TREE_DATA  *tree_data)  /* Tree data          */
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

}  /*  Node_Offset  */



/*****  Init_Slice  ********************************************************/
/*
*       Initalise the whole slice to a given value
*/
int  Init_Slice (double      *Slice,     /* (I/O) slice to be init'd  */
                 double       InitValue, /* (I)   init value          */
                 TREE_DATA   *tree_data) /* (I)   Tree data structure */
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
        DR_Error("Init_Slice: failed.");
    }
    return(status);

}  /* Init_Slice */



/*****   MaxOfTwoSlices  ******************************************************/
/*
 */
int MaxOfTwoSlices(double        *Slice,        /* (O) Slice to be modified */
                   double        *Slice1,       /* (I) First one to add     */
                   double        *Slice2,       /* (I) Second one to add    */
                   int            t,            /* (I) Current time point   */
                   TREE_DATA     *tree_data)    /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("MaxOfTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  MaxOfTwoSlices  */



/*****   MinOfTwoSlices  ******************************************************/
/*
 */
int MinOfTwoSlices(double        *Slice,        /* (O) Slice to be modified */
                   double        *Slice1,       /* (I) First one to add     */
                   double        *Slice2,       /* (I) Second one to add    */
                   int            t,            /* (I) Current time point   */
                   TREE_DATA     *tree_data)    /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("MinOfTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /*  MinOfTwoSlices  */


/*****   SlicePlusScalar  *************************************************/
/*
 *  Add scalar to slice and put in another slice
 *
 */
int SlicePlusScalar (double      *Slice,     /* (O) Slice to be modified */
                     double      *Slice1,    /* (I) Input slice          */
                     double      scalar,     /* (I) Spread               */
                     int         t,          /* (I) Current time point   */
                     TREE_DATA   *tree_data) /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

            SliceL  = Slice  + offset;
            Slice1L = Slice1 + offset;
        
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = Slice1L[i] + scalar;
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                SliceL  = Slice  + offset;
                Slice1L = Slice1 + offset;
        
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = Slice1L[j] + scalar;
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
                    offset = Node_Offset(3, i, j, t, tree_data);

                    SliceL  = Slice  + offset;
                    Slice1L = Slice1 + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = Slice1L[k] + scalar;
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("SlicePlusScalar: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /* SlicePlusScalar */


/*****   SliceTimesScalar  *************************************************/
/*
 *  Multiply slice by scalar and put in another slice
 *
 */
int SliceTimesScalar (double      *Slice,     /* (O) Slice to be modified */
                      double      *Slice1,    /* (I) Input slice          */
                      double      scalar,     /* (I) Factor               */
                      int         t,          /* (I) Current time point   */
                      TREE_DATA   *tree_data) /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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

}  /* SliceTimesScalar */


/*****   DivideTwoSlices  *************************************************/
/*
 *  Divide two slices. All slices must be of same dimension and memory
 *  must be pre-allocated. Caller must check that 2nd slice is non-zero
 *
 */
int DivideTwoSlices(double      *Slice,     /* (O) Slice to be modified */
                    double      *Slice1,    /* (I) Numerator            */
                    double      *Slice2,    /* (I) Denominator          */
                    int         t,          /* (I) Current time point   */
                    TREE_DATA   *tree_data) /* (I) Tree data            */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

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
                offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("DivideTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */

    return (SUCCESS);

}  /* DivideTwoSlices  */


/*****  Print_Slice1D  *******************************************************/
/*
*       Print slice to slice.prn
*/
int     Print_Slice1D (double      *Slice,      /* (O) Slice prices        */
                       char        *name,       /* (I) Slice name          */
                       int         t,           /* (I) Current time point  */
                       TREE_DATA   *tree_data)  /* (I) Tree data structure */
{

    double  *SliceL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim) */

    int     i;                                 /* Node indices           */
    int     offset;                            /* Node offset            */

    FILE    *fp;


    fp = fopen ("slice.prn", "a");
    if (fp == NULL) return (FAILURE);

    Top1    = tree_data->Top1[t];                 
    Bottom1 = tree_data->Bottom1[t];
    offset  = Node_Offset(1, 0, 0, t, tree_data);

    SliceL = Slice + offset;
    
    fprintf (fp, "Current date: %ld\n", tree_data->TPDate[t]);

    fprintf (fp, "Node:        ");
    for (i = Bottom1; i <= Top1; i ++)
    {
        fprintf (fp, "%3d        ", i);
    }

    fprintf (fp, "\n%12s:", name); 
    for (i = Bottom1; i <= Top1; i ++)
    {
        fprintf (fp, "%10.6f ", SliceL[i]);
    }

    fprintf (fp, "\n\n");
    fclose (fp);

    return (SUCCESS);

} /* Print_Slice1D */



/*****  CoVarTwoSlices  ****************************************************/
/*
*       Computes statistics of two slices, given the probabilities 
*       Returns: means, std devs, correlation 
*/
int StatsTwoSlices (double      *Mean1,         /* (O) Mean first slice    */
                    double      *StDev1,        /* (O) StDev first slice   */
                    double      *Mean2,         /* (O) Mean second slice   */
                    double      *StDev2,        /* (O) StDev second slice  */
                    double      *Corr,          /* (O) Correlation         */
                    double      *Slice1,        /* (I) First slice         */
                    double      *Slice2,        /* (I) Second slice        */
                    double      *Prob,          /* (I) Probabilities       */
                    int         t,              /* (I) Current time point  */
                    TREE_DATA   *tree_data)     /* (I) Tree data structure */
{
    double E1, E2, V1, V2, C;
    int    status = FAILURE;

    double *AuxSlice = Alloc_Slice (tree_data);
    if (AuxSlice == NULL) goto RETURN;


    if ((MultTwoSlicesAddAll (&E1,Slice1,Prob,t,tree_data) == FAILURE) ||
        (MultTwoSlicesAddAll (&E2,Slice1,Prob,t,tree_data) == FAILURE))
    {
        goto RETURN;
    }

    if ((MultiplyTwoSlices(AuxSlice,Slice1,Slice1,t,tree_data) == FAILURE) ||
        (MultTwoSlicesAddAll (&V1,AuxSlice,Prob,t,tree_data)   == FAILURE))
    {
        goto RETURN;
    }
    V1 -= (E1 * E1);
    V1 = sqrt (V1);

    if ((MultiplyTwoSlices(AuxSlice,Slice2,Slice2,t,tree_data) == FAILURE) ||
        (MultTwoSlicesAddAll (&V2,AuxSlice,Prob,t,tree_data)   == FAILURE))
    {
        goto RETURN;
    }
    V2 -= (E2 * E2);
    V2 = sqrt (V2);

    if ((MultiplyTwoSlices(AuxSlice,Slice1,Slice2,t,tree_data) == FAILURE) ||
        (MultTwoSlicesAddAll (&C,AuxSlice,Prob,t,tree_data)    == FAILURE))
    {
        goto RETURN;
    }
    C -= (E1 * E2);
    C /= (V1 * V2);    
    
    *Mean1  = E1;
    *Mean2  = E2;
    *StDev1 = V1;
    *StDev2 = V2;
    *Corr   = C;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error ("CoVarTwoSlices failed");
    }

    Free_Slice (AuxSlice, tree_data);

    return (status);

} /* CovarTwoSlices */
