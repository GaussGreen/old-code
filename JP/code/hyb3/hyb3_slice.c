/****************************************************************************/
/*      Time slice utility routines.                                        */
/****************************************************************************/
/*      SLICE.C                                                             */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"

/*****  Hyb3_Init_Slice  ********************************************************/
/*
*       Initalise the whole slice to a given value
*/
int  Hyb3_Init_Slice (double      *Slice,     /* (I/O) Slice to be init'd   */
                 int          dimension, /* (I)   Dim of desired slice */
                 double       InitValue, /* (I)   Init value           */
                 HYB3_TREE_DATA   *tree_data) /* (I)   Tree data structure  */
{
    int   status = FAILURE;

    long  ArraySize = 0L;
    long  i;

    if ((Slice == NULL) || (tree_data == NULL)) goto RETURN;

    /* calc the size of the slice array */
    switch (dimension)
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
        DR_Error("Hyb3_Init_Slice: failed.");
    }
    return(status);

}  /* Hyb3_Init_Slice */



/*****   Hyb3_SetSlice  ***********************************************************/
/*
 *    Set a slice of dimension up to three to a given value.
 */
int    Hyb3_SetSlice(TSLICE         Slice,         /* (O) Slice to be modified   */
                int            SliceDim,      /* (I) Dimension of above     */
                double         Value,         /* (I) Scalar to set slice to */
                int            t,             /* (I) Current time period    */
                HYB3_TREE_DATA     *tree_data)     /* (I) Tree data              */
{


    double
         *SliceL;


    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */


        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:
            SliceL = (double *)Slice + Hyb3_Node_Offset(1,0,0,t,tree_data);
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = Value;
            }
            break;

        case 2:
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL = (double *)Slice + Hyb3_Node_Offset(2,i,0,t,tree_data);
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = Value;
                }
            }
            break;

        case 3:
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL = (double *)Slice + Hyb3_Node_Offset(3,i,j,t,tree_data);
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = Value;
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_SetSlice)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb3_SetSlice  */




/*****  Hyb3_Copy_Slice  *********************************************************/
/*
*       Copy a time slice from a given time slice 
*/
int     Hyb3_CopySlice(TSLICE      Hyb3_CopySlice,    /* (O) Slice prices            */
                  TSLICE      OrgSlice,     /* (I) Original to be copied   */
                  int         SliceDim,     /* (I) Dimension of above      */
                  int         t,            /* (I) Current time point      */
                  HYB3_TREE_DATA   *tree_data)   /* (I) Tree data structure     */
{

    double  *CopySliceL;
    double  *OrgSliceL;

    int     Top1, Bottom1;                     /* Tree limits (1st dim)  */
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


    switch (SliceDim)
    {
        case 1:
        {
            offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

            CopySliceL = (double *)Hyb3_CopySlice + offset;
            OrgSliceL  = (double *)OrgSlice  + offset;
    
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
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                CopySliceL = (double *)Hyb3_CopySlice + offset;
                OrgSliceL  = (double *)OrgSlice  + offset;
    
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
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                    CopySliceL = (double *)Hyb3_CopySlice + offset;
                    OrgSliceL  = (double *)OrgSlice  + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        CopySliceL[k] = OrgSliceL[k];
                    }    
                }  /* for j */

            break;
        }
        default:
        {
            DR_Error("Hyb3_Copy_Slice: invalid slice dimension "
                        "(must be 1,2 or 3)!");
            return(FAILURE);        
        }
    } /* End of switch() */
    
    return (SUCCESS);

}  /* Hyb3_CopySlice */





/*****   Hyb3_ExpandSlice  ********************************************************/
/*
 *    Expand a given slice onto a larger dimension. Expansion can be 1->2,
 *  1->3 or 2->3.
 *
 *  Memory for the expanded slice MUST BE PRE-ALLOCATED!
 *
 */
int    Hyb3_ExpandSlice(TSLICE        NewSlice,    /* (O) Expanded slice         */
                   int           NewSliceDim, /* (I) Dimension of above     */
                   TSLICE        Slice,       /* (I) Input slice            */
                   int           SliceDim,    /* (I) Dimension of above     */
                   int           t,           /* (I) Current time period    */
                   HYB3_TREE_DATA    *tree_data)   /* (I) Tree data              */
{


    double
          *SliceL,
     
          *NewSliceL;

    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */


        i, j, k;                       /* Node indices                 */
       

    if(NewSliceDim == SliceDim)         /* Copy the slice across */
        return Hyb3_CopySlice(NewSlice,        /* (O) Slice prices            */
                         Slice,            /* (I) Original to be copied   */
                         SliceDim,        /* (I) Dimension of above      */
                         t,            /* (I) Current time point      */
                         tree_data);   /* (I) Tree data structure     */
        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    if ((SliceDim == 1) && (NewSliceDim == 2))
    {
        
        SliceL = (double *)Slice + Hyb3_Node_Offset(1,0,0,t,tree_data);

        for (i=Bottom1; i<=Top1; i++)
        {
            NewSliceL = (double *)NewSlice + Hyb3_Node_Offset(2,i,0,t,tree_data);

            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                NewSliceL[j] = SliceL[i];
            }
        }
    }
    else if ((SliceDim == 1) && (NewSliceDim == 3))
    {

        SliceL = (double *)Slice + Hyb3_Node_Offset(1,0,0,t,tree_data);

        for (i=Bottom1; i<=Top1; i++)
        {
            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                NewSliceL = (double *)NewSlice + Hyb3_Node_Offset(3,i,j,t,tree_data);

                for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                {
                    NewSliceL[k] = SliceL[i];

                }
            }
        }
    }
    else if ((SliceDim == 2) && (NewSliceDim == 3))
    {
        
        for (i=Bottom1; i<=Top1; i++)
        {

            SliceL = (double *)Slice + Hyb3_Node_Offset(2,i,0,t,tree_data);
            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                NewSliceL = (double *)NewSlice + Hyb3_Node_Offset(3,i,j,t,tree_data);
                for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                {
                    NewSliceL[k] = SliceL[j];

                }
            }
        }
    }
    else
    {
        DR_Error("Invalid slice expansion (must be 1->2, 2->3 or 1->3)."
                    " (Hyb3_ExpandSlice)\n");
        return(FAILURE);

    } 
     
            
    return (SUCCESS);



}  /*  Hyb3_ExpandSlice  */





/*****   Hyb3_AddTwoSlices  *******************************************************/
/*
 *    Add up to four slices onto a SumSlice. All five slice must be of same
 *  dimension and memory for SumSlice must be pre-allocated.
 *
 */
int Hyb3_AddTwoSlices(TSLICE         Slice,        /* (O) Slice to be modified   */
                 int            SliceDim,     /* (I) Dimension of above     */
                 TSLICE         Slice1,       /* (I) First one to add       */
                 TSLICE         Slice2,       /* (I) Second one to add      */
                 int            t,            /* (I) Current time period    */
                 HYB3_TREE_DATA     *tree_data)    /* (I) Tree data              */
{



    double
          *SumSliceL,
        
          *Slice1L,

          *Slice2L;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            SumSliceL = (double *)Slice  + offset;
            Slice1L   = (double *)Slice1 + offset;
            Slice2L   = (double *)Slice2 + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                SumSliceL[i] = Slice1L[i] + Slice2L[i];
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                SumSliceL = (double *)Slice  + offset;  
                Slice1L   = (double *)Slice1 + offset;  
                Slice2L   = (double *)Slice2 + offset;                                          

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SumSliceL[j] = Slice1L[j] + Slice2L[j];
                }
            }
            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {

                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    SumSliceL = (double *)Slice  + offset;      
                    Slice1L   = (double *)Slice1 + offset;      
                    Slice2L   = (double *)Slice2 + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SumSliceL[k] =  Slice1L[k]+Slice2L[k];
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_AddTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb3_AddTwoSlices  */


/*****   Hyb3_SubtractTwoSlices  **************************************************/
/*
 *    Subtract Slice2 from Slice1 (ie. Slice = Slice1-Slice2)  
 *    All slices must be of same dimension and memory preallocated
 *
 */
int Hyb3_SubtractTwoSlices(TSLICE     Slice,  /* (O) Slice to be modified   */
                 int            SliceDim,     /* (I) Dimension of above     */
                 TSLICE         Slice1,       /* (I) Slice to subtract from */
                 TSLICE         Slice2,       /* (I) Slice to subtract      */
                 int            t,            /* (I) Current time period    */
                 HYB3_TREE_DATA *tree_data)   /* (I) Tree data              */
{



    double
          *SubSliceL,
        
          *Slice1L,

          *Slice2L;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            SubSliceL = (double *)Slice  + offset;
            Slice1L   = (double *)Slice1 + offset;
            Slice2L   = (double *)Slice2 + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                SubSliceL[i] = Slice1L[i] - Slice2L[i];
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                SubSliceL = (double *)Slice  + offset;  
                Slice1L   = (double *)Slice1 + offset;  
                Slice2L   = (double *)Slice2 + offset;                                          

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SubSliceL[j] = Slice1L[j] - Slice2L[j];
                }
            }
            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {

                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    SubSliceL = (double *)Slice  + offset;      
                    Slice1L   = (double *)Slice1 + offset;      
                    Slice2L   = (double *)Slice2 + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SubSliceL[k] =  Slice1L[k] - Slice2L[k];
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_SubtractTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb3_SubtractTwoSlices  */




/*****   Hyb3_MultiplyTwoSlices  *************************************************/
/*
 *    Multipy two slices. All slices must be of same dimension and memory
 *  must be pre-allocated.
 *
 */
int Hyb3_MultiplyTwoSlices(TSLICE       Slice,     /* (O) Slice to be modified   */
                      int          SliceDim,  /* (I) Dimension of above     */
                      TSLICE       Slice1,    /* (I) First one to add       */
                      TSLICE       Slice2,    /* (I) Second one to add      */
                      int          t,         /* (I) Current time period    */
                      HYB3_TREE_DATA   *tree_data) /* (I) Tree data              */
{
    double
          *ProdSliceL,
        
          *Slice1L,

          *Slice2L;


    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */

        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            ProdSliceL = (double *)Slice  + offset;
            Slice1L    = (double *)Slice1 + offset;
            Slice2L    = (double *)Slice2 + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                ProdSliceL[i] = Slice1L[i] * Slice2L[i];
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                ProdSliceL = (double *)Slice  + offset;  
                Slice1L    = (double *)Slice1 + offset;  
                Slice2L    = (double *)Slice2 + offset;                                          

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    ProdSliceL[j] = Slice1L[j] * Slice2L[j];
                }
            }
            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {

                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    ProdSliceL = (double *)Slice  + offset;      
                    Slice1L    = (double *)Slice1 + offset;      
                    Slice2L    = (double *)Slice2 + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        ProdSliceL[k] =  Slice1L[k] * Slice2L[k];
    
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_MultiplyTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb3_MultiplyTwoSlices  */


/*****   Hyb3_DivideTwoSlices  *************************************************/
/*
 *    Divide two slices. All slices must be of same dimension and memory
 *  must be pre-allocated.
 *
 */
int Hyb3_DivideTwoSlices(TSLICE       Slice,     /* (O) Slice to be modified   */
                    int          SliceDim,  /* (I) Dimension of above     */
                    TSLICE       Slice1,    /* (I) First one to add       */
                    TSLICE       Slice2,    /* (I) Second one to add      */
                    int          t,         /* (I) Current time period    */
                    HYB3_TREE_DATA   *tree_data) /* (I) Tree data              */
{


 

    double
          *ProdSliceL,
        
          *Slice1L,

          *Slice2L;



    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */

        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            ProdSliceL = (double *)Slice  + offset;
            Slice1L    = (double *)Slice1 + offset;
            Slice2L    = (double *)Slice2 + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                ProdSliceL[i] = Slice1L[i] / (double) Slice2L[i];
            
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                ProdSliceL = (double *)Slice  + offset;  
                Slice1L    = (double *)Slice1 + offset;  
                Slice2L    = (double *)Slice2 + offset;                                          

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    ProdSliceL[j] = Slice1L[j] / (double) Slice2L[j];
                }
            }
            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {

                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    ProdSliceL = (double *)Slice  + offset;      
                    Slice1L    = (double *)Slice1 + offset;      
                    Slice2L    = (double *)Slice2 + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        ProdSliceL[k] =  Slice1L[k] / (double) Slice2L[k];

                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_MultiplyTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /* Hyb3_DivideTwoSlices  */


/*****   Hyb3_LCombTwoSlices  ****************************************************/
/*
 *    Linear combination of two slices onto a CombSlice. All five slice must 
 *  be of same dimension and memory for SumSlice must be pre-allocated.
 *
 *    CombSlice = a1 * Slice1 + a2 * Slice2
 *
 */
int Hyb3_LCombTwoSlices(TSLICE         CombSlice,  /* (O) Slice to be modified   */
                   int            SliceDim,   /* (I) Dimension of above     */
                   TSLICE         Slice1,     /* (I) First one to add       */
                   double         a1,         /* (I) Coeff of first slice   */
                   TSLICE         Slice2,     /* (I) Second one to add      */
                   double         a2,         /* (I) Coeff of second slice  */
                   int            t,          /* (I) Current time period    */
                   HYB3_TREE_DATA     *tree_data)  /* (I) Tree data              */
{


    double
          *CombSliceL,
          *Slice1L,
          *Slice2L;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */

        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            CombSliceL = (double *)CombSlice  + offset;
            Slice1L    = (double *)Slice1 + offset;
            Slice2L    = (double *)Slice2 + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                CombSliceL[i] = a1 * Slice1L[i] + a2 * Slice2L[i];
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                CombSliceL = (double *)CombSlice  + offset;  
                Slice1L    = (double *)Slice1 + offset;  
                Slice2L    = (double *)Slice2 + offset;                                          

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    CombSliceL[j] = a1 * Slice1L[j] + a2 * Slice2L[j];
                }
            }
            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {

                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    CombSliceL = (double *)CombSlice  + offset;      
                    Slice1L    = (double *)Slice1 + offset;      
                    Slice2L    = (double *)Slice2 + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        CombSliceL[k] =  a1* Slice1L[k] + a2 * Slice2L[k];
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_MultiplyTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
       

    return (SUCCESS);


}  /*  Hyb3_LCombTwoSlices  */


/*****   Hyb3_MultTwoSlicesAddAll  *********************************************/
/*
 *  Multipy two slices and perform addition of all nodes.  This operation
 *  is meaningful in the context of the "express DEV" facility.
 *  All slices must be of same dimension and memorym ust be pre-allocated.
 *
 */
int Hyb3_MultTwoSlicesAddAll(double      *Value,     /* (O) Return value         */
                        int          SliceDim,
                        double      *Slice1,    /* (I) First one to add     */
                        double      *Slice2,    /* (I) Second one to add    */
                        int         t,          /* (I) Current time point   */
                        HYB3_TREE_DATA   *tree_data) /* (I) Tree data            */
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

    switch (SliceDim)
    {
        case 1:
        {
            offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

            
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
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

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
            DR_Error("Hyb3_MultiplyTwoSlices: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }        
    } /* End of switch() */


    /* Output value */
    *Value = ValueL;

    return (SUCCESS);

}  /*  Hyb3_MultTwoSlicesAddAll  */




/*****   Hyb3_MaxTwoSlices  *******************************************************/
/*
 *     Make up an output slice by taking the max of two input slices.  All 
 *   input slices must be of same dimension and memory for MaxSlice must 
 *   be pre-allocated.
 *
 */
int Hyb3_MaxTwoSlices(TSLICE         Slice,        /* (O) Slice to be modified   */
                 int            SliceDim,     /* (I) Dimension of above     */
                 TSLICE         Slice1,       /* (I) First input slice      */
                 TSLICE         Slice2,       /* (I) Second input slice     */
                 int            t,            /* (I) Current time period    */
                 HYB3_TREE_DATA     *tree_data)    /* (I) Tree data              */
{



    double
          *MaxSliceL,
        
          *Slice1L,

          *Slice2L;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            MaxSliceL = (double *)Slice  + offset;
            Slice1L   = (double *)Slice1 + offset;
            Slice2L   = (double *)Slice2 + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                MaxSliceL[i] = MAX(Slice1L[i],Slice2L[i]);
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                MaxSliceL = (double *)Slice  + offset;  
                Slice1L   = (double *)Slice1 + offset;  
                Slice2L   = (double *)Slice2 + offset;                                          

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    MaxSliceL[j] = MAX(Slice1L[j],Slice2L[j]);
                }
            }
            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {

                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    MaxSliceL = (double *)Slice  + offset;      
                    Slice1L   = (double *)Slice1 + offset;      
                    Slice2L   = (double *)Slice2 + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        MaxSliceL[k] =  MAX(Slice1L[k],Slice2L[k]);
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_MaxTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb3_MaxTwoSlices  */





/*****   Hyb3_MinTwoSlices  *******************************************************/
/*
 *     Make up an output slice by taking the min of two input slices.  All 
 *   input slices must be of same dimension and memory for MinSlice must 
 *   be pre-allocated.
 *
 */
int Hyb3_MinTwoSlices(TSLICE         Slice,        /* (O) Slice to be modified   */
                 int            SliceDim,     /* (I) Dimension of above     */
                 TSLICE         Slice1,       /* (I) First input slice      */
                 TSLICE         Slice2,       /* (I) Second input slice     */
                 int            t,            /* (I) Current time period    */
                 HYB3_TREE_DATA     *tree_data)    /* (I) Tree data              */
{



    double
          *MinSliceL,
        
          *Slice1L,

          *Slice2L;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            MinSliceL = (double *)Slice  + offset;
            Slice1L   = (double *)Slice1 + offset;
            Slice2L   = (double *)Slice2 + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                MinSliceL[i] = MIN(Slice1L[i],Slice2L[i]);
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                MinSliceL = (double *)Slice  + offset;  
                Slice1L   = (double *)Slice1 + offset;  
                Slice2L   = (double *)Slice2 + offset;                                          

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    MinSliceL[j] = MIN(Slice1L[j],Slice2L[j]);
                }
            }
            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {

                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    MinSliceL = (double *)Slice  + offset;      
                    Slice1L   = (double *)Slice1 + offset;      
                    Slice2L   = (double *)Slice2 + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        MinSliceL[k] =  MIN(Slice1L[k],Slice2L[k]);
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_MaxTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb3_MinTwoSlices  */







/*****  Hyb3_AddScalar  ***********************************************************/
/*
 *    Add a scalar to a slice of dimension up to three.
 */
int    Hyb3_AddScalar(TSLICE         Slice,         /* (O) Slice to be modified  */
                 int            SliceDim,      /* (I) Dimension of above    */
                 double         Scalar,        /* (I) Scalar to be added    */
                 int            t,             /* (I) Current time period   */
                 HYB3_TREE_DATA     *tree_data)     /* (I) Tree data             */
{




    double
          *SliceL;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            SliceL = (double *)Slice  + offset;
            

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = SliceL[i] + Scalar;
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                SliceL = (double *)Slice  + offset;                                  

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = SliceL[j] + Scalar;
                }
            }
            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {

                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    SliceL = (double *)Slice  + offset;      
       
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] =  SliceL[k] + Scalar;
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_AddScalar)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb3_AddScalar  */







/*****  Hyb3_MultiplyScalar  ******************************************************/
/*
 *    Multiply a slice of dimension up to three by a scalar.
 */
int    Hyb3_MultiplyScalar(TSLICE        Slice,     /* (O) Slice to be modified  */
                      int           SliceDim,  /* (I) Dimension of above    */
                      double        Scalar,    /* (I) Scalar to be mult     */
                      int           t,         /* (I) Current time period   */
                      HYB3_TREE_DATA    *tree_data) /* (I) Tree data             */
{


 

    double
          *SliceL;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            SliceL = (double *)Slice  + offset;
            

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = SliceL[i] * Scalar;
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                SliceL = (double *)Slice  + offset;                                  

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = SliceL[j] * Scalar;
                }
            }
            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {

                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    SliceL = (double *)Slice  + offset;      
       
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] =  SliceL[k] * Scalar;
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_MultiplyScalar)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);



}  /*  Hyb3_MultiplyScalar  */







/*****  Hyb3_MaxMinOnSlice  ******************************************************/
/*
 *    Apply a maximum and a minimum to the values in a slice.
 */
int    Hyb3_MaxMinOnSlice(TSLICE         Slice,     /* (O) Slice to be modified  */
                     int            SliceDim,  /* (I) Dimension of above    */
                     double         Min,       /* (I) Minimum level         */
                     double         Max,       /* (I) Maximum level         */
                     int            t,         /* (I) Current time period   */
                     HYB3_TREE_DATA     *tree_data) /* (I) Tree data             */
{


    double
          *SliceL;


    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */

        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            SliceL = (double *)Slice  + offset;
            
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MAXMIN(SliceL[i],Max,Min);
            }
            break;

        case 2:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                SliceL = (double *)Slice  + offset;         

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = MAXMIN(SliceL[j],Max,Min);
                }
            }
            break;

        case 3:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    SliceL = (double *)Slice  + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MAXMIN(SliceL[k],Max,Min);
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_MaxMinOnSlice)\n");
            return(FAILURE);

    } /* End of switch() */


    return (SUCCESS);


}  /*  Hyb3_MaxMinOnSlice  */







/*****  Hyb3_MaxOnSlice  *********************************************************/
/*
 *    Apply a maximum to the values in a slice.
 */
int    Hyb3_MaxOnSlice(TSLICE         Slice,        /* (O) Slice to be modified  */
                  int            SliceDim,     /* (I) Dimension of above    */
                  double         Max,          /* (I) Maximum level         */
                  int            t,            /* (I) Current time period   */
                  HYB3_TREE_DATA     *tree_data)    /* (I) Tree data             */
{


    double
          *SliceL;


    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */

        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            SliceL = (double *)Slice  + offset;
            
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MIN(SliceL[i],Max);
            }
            break;

        case 2:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                SliceL = (double *)Slice  + offset;     
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = MIN(SliceL[j],Max);
                }
            }
            break;

        case 3:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    SliceL = (double *)Slice  + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MIN(SliceL[k],Max);
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_MaxOnSlice)\n");
            return(FAILURE);

    } /* End of switch() */


    return (SUCCESS);


}  /*  Hyb3_MaxOnSlice  */






/*****  Hyb3_MinOnSlice  *********************************************************/
/*
 *    Apply a minimum to the values in a slice.
 */
int    Hyb3_MinOnSlice(TSLICE         Slice,        /* (O) Slice to be modified  */
                  int            SliceDim,     /* (I) Dimension of above    */
                  double         Min,          /* (I) Minimum level         */
                  int            t,            /* (I) Current time period   */
                  HYB3_TREE_DATA     *tree_data)    /* (I) Tree data             */
{


    double
          *SliceL;


    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */


        i, j, k,                       /* Node indices                 */
        offset;
       

        
        Top1    = tree_data->Top1[t];                    
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];                    
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];                    
        Bottom3 = tree_data->Bottom3[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            SliceL = (double *)Slice  + offset;
            
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MAX(SliceL[i],Min);
            }
            break;

        case 2:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                SliceL = (double *)Slice  + offset; 
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    SliceL[j] = MAX(SliceL[j],Min);
                }
            }
            break;

        case 3:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);    
                    SliceL = (double *)Slice  + offset;    
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MAX(SliceL[k],Min);
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2 or 3)."
                    " (Hyb3_MinOnSlice)\n");
            return(FAILURE);

    } /* End of switch() */


    return (SUCCESS);


}  /*  Hyb3_MinOnSlice  */





/*****  Hyb3_CheckDimRebate  *****************************************************/
/*                                                                           
 *      The rebate can be a constant (dimension 0), a 1-D variable or a 2-D
 *      variable.
 *
 */
int   Hyb3_CheckDimRebate(int     Dim)

{
    if (Dim != 0 && Dim != 1 && Dim != 2)
    {
        DR_Error("Dimension of rebate in KO must be 0, 1 or 2.\n");
        return(FAILURE);
    }
    return(SUCCESS);
   
}



/*****  Hyb3_CheckDimIndex  *****************************************************/
/*                                                                           
 *      The index can be either a 1-D variable or a 2-D variable.
 *
 */
int   Hyb3_CheckDimIndex(int     Dim)

{
    if (Dim != 1 && Dim != 2)
    {
        DR_Error("Dimension of index for KO must be 1 or 2.\n");
        return(FAILURE);
    }
    return(SUCCESS);
   
}





/*****  Hyb3_GetValueAtNode  *****************************************************/
/*                                                                           
 *      Obtains the 'local' value of the rebate, be it a constant rebate, a
 *      1-D variable or a 2-D variable.   This function hides the dimension
 *      generality from  the  caller, but it has no means to check that the 
 *      adequate  amount of space has been allocated under the void * being
 *      passed.
 *
 */
double   Hyb3_GetValueAtNode(int          Dim,
                        TSLICE       ValuePtr,
                        int          i,
                        int          j,
                        int          k,
                        int          t,
                        HYB3_TREE_DATA   *tree_data)

{
    
    
    double    *SliceL;
    double     Value;
   

    /* NB: Dim must be checked beforehand by Hyb3_CheckDimRebate()  */
    
    switch (Dim)
    {
        case 0:
            Value = ((double *)ValuePtr)[0];
            break;
            
        case 1:
            SliceL = (double *)ValuePtr + Hyb3_Node_Offset(1,0,0,t,tree_data);
            Value  = SliceL[i];
            break;
            
        case 2:
            SliceL = (double *)ValuePtr + Hyb3_Node_Offset(2,i,0,t,tree_data);
            Value  = SliceL[j];
            break;

        default:
            SliceL = (double *)ValuePtr + Hyb3_Node_Offset(3,i,j,t,tree_data);
            Value  = SliceL[k];

    }

    return(Value);
}

 
/*****  Hyb3_SmoothStepUp  *********************************************************/
/**
 *  Smoothing a payoff depending on the dimension
 */
int    Hyb3_SmoothStepUp(
                  double*        Slice,           /**< (O) Slice to be modified */
                  double*        Index,           /**< (I) Reference slice      */
                  double         Up,              /**< (I) Maximium level       */
                  double         Barrier,         /**< (I) Barrier level        */
                  int            SmoothingOn,     /**< (I) TRUE = smoothing on  */
                  int            t,               /**< (I) Current time period  */
                  int            tree_dim,        /**< (I) tree dimension       */
                  HYB3_TREE_DATA *tree_data)      /**< (I) Tree data            */
{
    double  *SliceL;
    double  *IndexL;

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
    
    switch (tree_dim)
    {
        case 1:
        {
            offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL = Slice + offset;
            IndexL = Index + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                if (SmoothingOn)
                    Step = Hyb3_GetIndexStep(Index,1,i,0,0,t,tree_data);
                SliceL[i] = DrSmoothStep(Up, SliceL[i], IndexL[i], Barrier, Step);
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                SliceL = Slice + offset;
                IndexL = Index + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    if (SmoothingOn)
                        Step = Hyb3_GetIndexStep(Index,2,i,j,0,t,tree_data);
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
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                    SliceL = Slice + offset;
                    IndexL = Index + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        if (SmoothingOn)
                            Step = Hyb3_GetIndexStep(Index,3,i,j,k,t,tree_data);
                        SliceL[k] = DrSmoothStep(Up, SliceL[k], IndexL[k], Barrier, Step);
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Hyb3_SmoothStepUp: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }
    } /* End of switch() */
    
    return (SUCCESS);

}


/*****  Hyb3_Node_Offset  ********************************************************/
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
int     Hyb3_Node_Offset (   int         Dim,        /* Dimension          */
                        int         i,          /* Node indices       */
                        int         j,          /*                    */
                        int         t,          /* Current time point */
                        HYB3_TREE_DATA const* tree_data)  /* Tree data          */
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

}  /*  Hyb3_Node_Offset  */

/*****  Hyb3_Slice_Dim  *********************************************************/
/**
*        Returns the slice dimension corresponding to a given dev mode.
*/
int Hyb3_Slice_Dim(int DMode)
{


    if (DMode == DISC_1D_NOCUPS)     return 1;
    if (DMode == DISC_2D_CUPS)     return 2;
    if (DMode == DISC_3D_CUPS)   return 3;
    if (DMode == DISC_2D_NOCUPS)   return 2;
    if (DMode == DISC_2D_1IR2F_NOCUPS)   return 2;
    if (DMode == DISC_3D_2IR2F1D_CUPS)  return 3;


    return -1; /* not a valid Dmode type */
}  /* Hyb3_Free_Slice */

