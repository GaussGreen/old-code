/****************************************************************************/
/*      Time slice utility routines.                                        */
/****************************************************************************/
/*      HYB4_SLICE.C                                                        */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hyb4_lib.h"

int Hyb4_Offset_Calc(int             Dim,
                     HYB4_TREE_DATA *tree_data)
{
    int    iMin;
    int    iMax;

    int    jMin;
    int    jMax;
    
    int    kMin;
    int    kMax;

    int    LMin;
    int    LMax;

    int    y0, y1, y2, y3;
    int    m0, m1, m2, m3;

    int    t;
    int    T = tree_data->NbTP;
    int    xT = tree_data->xT;
    int    i, j, k;

    for (t=0; t<=T; ++t)
    {
        y0 = 0;
        y1 = 0;
        y2 = 0;
        y3 = 0;

        iMin = tree_data->iMin[t];
        iMax = tree_data->iMax[t];

        y0 -= iMin;

        tree_data->NodeOffset0[t] = y0;

        y0 += iMax + 1;

        if (Dim > 1)
        {
            for (i=iMin; i<=iMax; ++i)
            {
                jMin = tree_data->jMin[t][i];
                jMax = tree_data->jMax[t][i];

                y1 -= jMin;

                tree_data->NodeOffset1[t][i] = y1;

                y1 += jMax + 1;
            
                if (Dim > 2)
                {
                    for (j=jMin; j<=jMax; ++j)
                    {
                        kMin = tree_data->kMin[t][i][j];
                        kMax = tree_data->kMax[t][i][j];

                        y2 -= kMin;

                        tree_data->NodeOffset2[t][i][j] = y2;

                        y2 += kMax + 1;
            
                        if ((Dim > 3) && (t <= xT))
                        {
                            for (k=kMin; k<=kMax; ++k)
                            {
                                LMin = tree_data->LMin[t][i][j][k];
                                LMax = tree_data->LMax[t][i][j][k];

                                y3 -= LMin;

                                tree_data->NodeOffset3[t][i][j][k] = y3;

                                y3 += LMax + 1;
                            }
                        }
                    }
                }
            }
        }
    }

    m0 = 0;
    m1 = 0;
    m2 = 0;
    m3 = 0;

    for (t=0; t<=T; ++t)
    {
        if (Dim > 0)
        {
            iMax = tree_data->iMax[t];

            m0 = tree_data->NodeOffset0[t] + iMax;
            
            tree_data->MaxIndex[0] = MAX(tree_data->MaxIndex[0], m0);
        }

        if (Dim > 1)
        {
            jMax = tree_data->jMax[t][iMax];
            
            m1 = tree_data->NodeOffset1[t][iMax] + jMax;
            
            tree_data->MaxIndex[1] = MAX(tree_data->MaxIndex[1], m1);
        }

        if (Dim > 2)
        {
            kMax = tree_data->kMax[t][iMax][jMax];
            
            m2 = tree_data->NodeOffset2[t][iMax][jMax] + kMax;
            
            tree_data->MaxIndex[2] = MAX(tree_data->MaxIndex[2], m2);
        }

        if ((Dim > 3) && (t <= xT))
        {
            LMax = tree_data->LMax[t][iMax][jMax][kMax];

            m3 = tree_data->NodeOffset3[t][iMax][jMax][kMax] + LMax;

            tree_data->MaxIndex[3] = MAX(tree_data->MaxIndex[3], m3);
        }
    }

    return SUCCESS;
}

/*****  Hyb4_Init_Slice  ********************************************************/
/*
*       Initalise the whole slice to a given value
*/
int  Hyb4_Init_Slice (double           *Slice,     /* (I/O) Slice to be init'd   */
                      int               dimension, /* (I)   Dim of desired slice */
                      double            InitValue, /* (I)   Init value           */
                      HYB4_TREE_DATA   *tree_data) /* (I)   Tree data structure  */
{
    int   status = FAILURE;

    int   MaxIndex;
    int   i;

    if ((Slice == NULL) || (tree_data == NULL))
    {
        goto RETURN;
    }

    if ((dimension < 1) ||
        (dimension > 4))
    {
        DR_Error("Bad slice dimension. (Hyb4_Init_Slice)\n");
        goto RETURN;
    }

    MaxIndex = tree_data->MaxIndex[dimension-1];

    for (i=0; i<=MaxIndex; ++i)
    {
        Slice[i] = InitValue;
    }

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("Hyb4_Init_Slice: failed.");
    }

    return status;

}  /* Hyb4_Init_Slice */

/*****   Hyb4_SetSlice  ***********************************************************/
/*
 *    Set a slice of dimension up to four to a given value.
 */
int    Hyb4_SetSlice(TSLICE         Slice,         /* (O) Slice to be modified   */
                int            SliceDim,      /* (I) Dimension of above     */
                double         Value,         /* (I) Scalar to set slice to */
                int            t,             /* (I) Current time period    */
                HYB4_TREE_DATA   *tree_data)     /* (I) Tree data              */
{


    double
         *SliceL;


    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (4rd dim) */


        i, j, k, L;                    /* Node indices                 */
       

        
    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];



    switch (SliceDim)
    {
        case 1:
            SliceL = (double *)Slice + tree_data->NodeOffset0[t];
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = Value;
            }
            break;

        case 2:
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL = (double *)Slice + tree_data->NodeOffset1[t][i];
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
                    SliceL = (double *)Slice + tree_data->NodeOffset2[t][i][j];
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = Value;
                    }
                }
            }
            break;

        case 4:
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL = (double *)Slice + tree_data->NodeOffset3[t][i][j][k];

                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            SliceL[L] = Value;
                        }
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2, 3, or 4)."
                    " (Hyb4_SetSlice)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb4_SetSlice  */




/*****  Hyb4_Copy_Slice  *********************************************************/
/*
*       Copy a time slice from a given time slice 
*/
int     Hyb4_CopySlice(TSLICE      CopySlice,    /* (O) Slice prices            */
                  TSLICE      OrgSlice,     /* (I) Original to be copied   */
                  int         SliceDim,     /* (I) Dimension of above      */
                  int         t,            /* (I) Current time point      */
                  HYB4_TREE_DATA *tree_data)   /* (I) Tree data structure     */
{

    double  *CopySliceL;
    double  *OrgSliceL;

    int     Top1, Bottom1;                     /* Tree limits (1st dim)  */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)  */
    int     ***Top4, ***Bottom4;                 /* Tree limits (3rd dim)  */

    int     i, j, k, L;                           /* Node indices           */
    int     offset;                            /* Node offset            */


    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];

    switch (SliceDim)
    {
        case 1:
        {
            offset = tree_data->NodeOffset0[t];

            CopySliceL = (double *)CopySlice + offset;
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
                offset = tree_data->NodeOffset1[t][i];

                CopySliceL = (double *)CopySlice + offset;
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
                    offset = tree_data->NodeOffset2[t][i][j];

                    CopySliceL = (double *)CopySlice + offset;
                    OrgSliceL  = (double *)OrgSlice  + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        CopySliceL[k] = OrgSliceL[k];
                    }    
                }  /* for j */

            break;
        }
        case 4:
        {        
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];

                        CopySliceL = (double *)CopySlice + offset;
                        OrgSliceL  = (double *)OrgSlice  + offset;

                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                             CopySliceL[L] = OrgSliceL[L];
                        }
                    }    
                }
            }

            break;
        }
        default:
        {
            DR_Error("Hyb4_CopySlice: invalid slice dimension "
                        "(must be 1,2, 3, or 4)!");
            return(FAILURE);        
        }
    } /* End of switch() */
    
    return (SUCCESS);

}  /* Hyb4_CopySlice */





/*****   Hyb4_ExpandSlice  ********************************************************/
/*
 *    Expand a given slice onto a larger dimension. Expansion can be 1->2,
 *  1->3 or 2->3, or 1->4, or 2->4, or 3->4
 *
 *  Memory for the expanded slice MUST BE PRE-ALLOCATED!
 *
 */
int    Hyb4_ExpandSlice(TSLICE        NewSlice,    /* (O) Expanded slice         */
                   int           NewSliceDim, /* (I) Dimension of above     */
                   TSLICE        Slice,       /* (I) Input slice            */
                   int           SliceDim,    /* (I) Dimension of above     */
                   int           t,           /* (I) Current time period    */
                   HYB4_TREE_DATA  *tree_data)   /* (I) Tree data              */
{


    double
          *SliceL,
     
          *NewSliceL;

    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (4rd dim) */


        i, j, k, L;                       /* Node indices                 */
       

    if(NewSliceDim == SliceDim)         /* Copy the slice across */

        return Hyb4_CopySlice(NewSlice,        /* (O) Slice prices            */
                         Slice,            /* (I) Original to be copied   */
                         SliceDim,        /* (I) Dimension of above      */
                         t,            /* (I) Current time point      */
                         tree_data);   /* (I) Tree data structure     */
        
    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];



    if ((SliceDim == 1) && (NewSliceDim == 2))
    {
        
        SliceL = (double *)Slice + tree_data->NodeOffset0[t];

        for (i=Bottom1; i<=Top1; i++)
        {
            NewSliceL = (double *)NewSlice + tree_data->NodeOffset1[t][i];

            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                NewSliceL[j] = SliceL[i];
            }
        }
    }
    else if ((SliceDim == 1) && (NewSliceDim == 3))
    {

        SliceL = (double *)Slice + tree_data->NodeOffset0[t];

        for (i=Bottom1; i<=Top1; i++)
        {
            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                NewSliceL = (double *)NewSlice + tree_data->NodeOffset2[t][i][j];

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

            SliceL = (double *)Slice + tree_data->NodeOffset1[t][i];
            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                NewSliceL = (double *)NewSlice + tree_data->NodeOffset2[t][i][j];
                for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                {
                    NewSliceL[k] = SliceL[j];

                }
            }
        }
    }
    else if ((SliceDim == 1) && (NewSliceDim == 4))
    {
        SliceL = (double *)Slice + tree_data->NodeOffset0[t];
            
        for (i=Bottom1; i<=Top1; i++)
        {
            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                {
                    NewSliceL = (double *)NewSlice + tree_data->NodeOffset3[t][i][j][k];
                
                    for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                    {
                        NewSliceL[L] = SliceL[i];
                    }
                }
            }
        }
    }
    else if ((SliceDim == 2) && (NewSliceDim == 4))
    {
        for (i=Bottom1; i<=Top1; i++)
        {
            SliceL = (double *)Slice + tree_data->NodeOffset1[t][i];
            
            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                {
                    NewSliceL = (double *)NewSlice + tree_data->NodeOffset3[t][i][j][k];
                
                    for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                    {
                        NewSliceL[L] = SliceL[j];
                    }
                }
            }
        }
    }
    else if ((SliceDim == 3) && (NewSliceDim == 4))
    {
        for (i=Bottom1; i<=Top1; i++)
        {
            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                SliceL = (double *)Slice + tree_data->NodeOffset2[t][i][j];
            
                for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                {
                    NewSliceL = (double *)NewSlice + tree_data->NodeOffset3[t][i][j][k];
                
                    for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                    {
                        NewSliceL[L] = SliceL[k];
                    }
                }
            }
        }
    }
    else
    {
        DR_Error("Invalid slice expansion (must be 1->2, 2->3, 1->3, 1->4, 2->4, 3->4)."
                    " (Hyb4_ExpandSlice)\n");
        return(FAILURE);

    } 
     
            
    return (SUCCESS);



}  /*  Hyb4_ExpandSlice  */





/*****   Hyb4_AddTwoSlices  *******************************************************/
/*
 */
int Hyb4_AddTwoSlices(TSLICE         Slice,        /* (O) Slice to be modified   */
                 int            SliceDim,     /* (I) Dimension of above     */
                 TSLICE         Slice1,       /* (I) First one to add       */
                 TSLICE         Slice2,       /* (I) Second one to add      */
                 int            t,            /* (I) Current time period    */
                 HYB4_TREE_DATA   *tree_data)    /* (I) Tree data              */
{



    double
          *SumSliceL,
        
          *Slice1L,

          *Slice2L;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k, L;                    /* Node indices                 */
       

        
    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];



    switch (SliceDim)
    {
        case 1:

            offset = tree_data->NodeOffset0[t];    
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
               
                offset = tree_data->NodeOffset1[t][i];    
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

                    offset = tree_data->NodeOffset2[t][i][j];    
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

        case 4:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];    

                        SumSliceL = (double *)Slice  + offset;      
                        Slice1L   = (double *)Slice1 + offset;      
                        Slice2L   = (double *)Slice2 + offset;      
    
                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            SumSliceL[L] =  Slice1L[L]+Slice2L[L];
                        }
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2, 3 or 4)."
                    " (Hyb4_AddTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb4_AddTwoSlices  */

/*****   Hyb4_MultiplyTwoSlices  *******************************************************/

int Hyb4_MultiplyTwoSlices(TSLICE         Slice,        /* (O) Slice to be modified   */
                 int            SliceDim,     /* (I) Dimension of above     */
                 TSLICE         Slice1,       /* (I) First one to add       */
                 TSLICE         Slice2,       /* (I) Second one to add      */
                 int            t,            /* (I) Current time period    */
                 HYB4_TREE_DATA   *tree_data)    /* (I) Tree data              */
{



    double
          *ProdSliceL,
        
          *Slice1L,

          *Slice2L;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k, L;                    /* Node indices                 */
       

        
    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];



    switch (SliceDim)
    {
        case 1:

            offset = tree_data->NodeOffset0[t];    
            ProdSliceL = (double *)Slice  + offset;
            Slice1L   = (double *)Slice1 + offset;
            Slice2L   = (double *)Slice2 + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                ProdSliceL[i] = Slice1L[i] * Slice2L[i];
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = tree_data->NodeOffset1[t][i];    
                ProdSliceL = (double *)Slice  + offset;  
                Slice1L   = (double *)Slice1 + offset;  
                Slice2L   = (double *)Slice2 + offset;                                          

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

                    offset = tree_data->NodeOffset2[t][i][j];    
                    ProdSliceL = (double *)Slice  + offset;      
                    Slice1L   = (double *)Slice1 + offset;      
                    Slice2L   = (double *)Slice2 + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        ProdSliceL[k] =  Slice1L[k] * Slice2L[k];
                    }
                }
            }
            break;

        case 4:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];    

                        ProdSliceL = (double *)Slice  + offset;      
                        Slice1L   = (double *)Slice1 + offset;      
                        Slice2L   = (double *)Slice2 + offset;      
    
                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            ProdSliceL[L] =  Slice1L[L] * Slice2L[L];
                        }
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2, 3 or 4)."
                    " (Hyb4_MultiplyTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb4_MultiplyTwoSlices  */

/*****   Hyb4_InvertSlice  *******************************************************/

int Hyb4_InvertSlice(TSLICE          Slice,        /* (O) Slice to be modified   */
                     int             SliceDim,     /* (I) Dimension of above     */
                     int             t,            /* (I) Current time period    */
                     HYB4_TREE_DATA *tree_data)    /* (I) Tree data              */
{
    double

          *SliceL;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k, L;                    /* Node indices                 */
       
    int status = FAILURE;
        
    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];

    switch (SliceDim)
    {
        case 1:

            offset = tree_data->NodeOffset0[t];    

            SliceL   = (double *)Slice + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                if (fabs(SliceL[i]) < TINY)
                {
                    goto done;
                }

                SliceL[i] = 1 / SliceL[i];
            }

            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = tree_data->NodeOffset1[t][i];    

                SliceL   = (double *)Slice + offset;                                          

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    if (fabs(SliceL[j]) < TINY)
                    {
                        goto done;
                    }

                    SliceL[j] = 1 / SliceL[j];
                }
            }

            break;

        case 3:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = tree_data->NodeOffset2[t][i][j];    

                    SliceL = (double *)Slice + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        if (fabs(SliceL[k]) < TINY)
                        {
                            goto done;
                        }

                        SliceL[k] = 1 / SliceL[k];
                    }
                }
            }

            break;

        case 4:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];    

                        SliceL = (double *)Slice  + offset;      
    
                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            if (fabs(SliceL[L]) < TINY)
                            {
                                goto done;
                            }

                            SliceL[L] = 1 / SliceL[L];
                        }
                    }
                }
            }

            break;

        default:

            DR_Error("Invalid slice dimension (must be 1, 2, 3 or 4)."
                    " (Hyb4_InvertSlice)\n");
            
            goto done;

    } /* End of switch() */

    status = SUCCESS;

done:

    return status;

}  /*  Hyb4_InvertSlice  */

/*****   Hyb4_LCombTwoSlices  ****************************************************/
/*
 *    Linear combination of two slices onto a CombSlice. All five slice must 
 *  be of same dimension and memory for SumSlice must be pre-allocated.
 *
 *    CombSlice = a1 * Slice1 + a2 * Slice2
 *
 */
int Hyb4_LCombTwoSlices(TSLICE         CombSlice,  /* (O) Slice to be modified   */
                   int            SliceDim,   /* (I) Dimension of above     */
                   TSLICE         Slice1,     /* (I) First one to add       */
                   double         a1,         /* (I) Coeff of first slice   */
                   TSLICE         Slice2,     /* (I) Second one to add      */
                   double         a2,         /* (I) Coeff of second slice  */
                   int            t,          /* (I) Current time period    */
                   HYB4_TREE_DATA   *tree_data)  /* (I) Tree data              */
{


    double
          *CombSliceL,
          *Slice1L,
          *Slice2L;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (4rd dim) */

        offset,
        i, j, k, L;                       /* Node indices                 */
       

        
    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];



    switch (SliceDim)
    {
        case 1:

            offset = tree_data->NodeOffset0[t];
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
               
                offset = tree_data->NodeOffset1[t][i];
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

                    offset = tree_data->NodeOffset2[t][i][j];
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

        case 4:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];

                        CombSliceL = (double *)CombSlice  + offset;      
                        Slice1L    = (double *)Slice1 + offset;      
                        Slice2L    = (double *)Slice2 + offset;      

                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            CombSliceL[L] =  a1* Slice1L[L] + a2 * Slice2L[L];
                        }
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2, 3 or 4)."
                    " (Hyb4_MultiplyTwoSlices)\n");
            return(FAILURE);

    } /* End of switch() */
       

    return (SUCCESS);


}  /*  Hyb4_LCombTwoSlices  */

/*****  Hyb4_AddScalar  ***********************************************************/
/*
 *    Add a scalar to a slice of dimension up to three.
 */
int    Hyb4_AddScalar(TSLICE         Slice,         /* (O) Slice to be modified  */
                 int            SliceDim,      /* (I) Dimension of above    */
                 double         Scalar,        /* (I) Scalar to be added    */
                 int            t,             /* (I) Current time period   */
                 HYB4_TREE_DATA   *tree_data)     /* (I) Tree data             */
{




    double
          *SliceL;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k, L;                       /* Node indices                 */
       

        
    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];


    switch (SliceDim)
    {
        case 1:

            offset = tree_data->NodeOffset0[t];
            SliceL = (double *)Slice  + offset;
            

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = SliceL[i] + Scalar;
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = tree_data->NodeOffset1[t][i];
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

                    offset = tree_data->NodeOffset2[t][i][j];
                    SliceL = (double *)Slice  + offset;      
       
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] =  SliceL[k] + Scalar;
                    }
                }
            }
            break;

        case 4:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];
                        SliceL = (double *)Slice  + offset;      
       
                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            SliceL[L] =  SliceL[L] + Scalar;
                        }
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2, 3, or 4)."
                    " (Hyb4_AddScalar)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}  /*  Hyb4_AddScalar  */


int Hyb4_MultiplyScalar(TSLICE         Slice,         /* (O) Slice to be modified  */
                 int            SliceDim,      /* (I) Dimension of above    */
                 double         Scalar,        /* (I) Scalar to be added    */
                 int            t,             /* (I) Current time period   */
                 HYB4_TREE_DATA   *tree_data)     /* (I) Tree data             */
{




    double
          *SliceL;

    int

          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (3rd dim) */


        offset,
        i, j, k, L;                       /* Node indices                 */
       

        
    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];


    switch (SliceDim)
    {
        case 1:

            offset = tree_data->NodeOffset0[t];
            SliceL = (double *)Slice  + offset;
            

            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = SliceL[i] * Scalar;
            }
            break;

        case 2:

            for (i=Bottom1; i<=Top1; i++)
            {
               
                offset = tree_data->NodeOffset1[t][i];
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

                    offset = tree_data->NodeOffset2[t][i][j];
                    SliceL = (double *)Slice  + offset;      
       
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] =  SliceL[k] * Scalar;
                    }
                }
            }
            break;

        case 4:

            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];
                        SliceL = (double *)Slice  + offset;      
       
                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            SliceL[L] =  SliceL[L] * Scalar;
                        }
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2, 3 or 4)."
                    " (Hyb4_AddScalar)\n");
            return(FAILURE);

    } /* End of switch() */
            

    return (SUCCESS);


}

/*****  Hyb4_GetValueAtNode  *****************************************************/
/*                                                                           
 *
 */
double   Hyb4_GetValueAtNode(int          Dim,
                        TSLICE       ValuePtr,
                        int          i,
                        int          j,
                        int          k,
                        int          L,
                        int          t,
                        HYB4_TREE_DATA   *tree_data)

{
    
    
    double    *SliceL;
    double     Value;
   
    switch (Dim)
    {
        case 0:
            Value = ((double *)ValuePtr)[0];
            break;
            
        case 1:
            SliceL = (double *)ValuePtr + tree_data->NodeOffset0[t];
            Value  = SliceL[i];
            break;
            
        case 2:
            SliceL = (double *)ValuePtr + tree_data->NodeOffset1[t][i];
            Value  = SliceL[j];
            break;

        case 3:
            SliceL = (double *)ValuePtr + tree_data->NodeOffset2[t][i][j];
            Value  = SliceL[k];
            break;

        default:
            SliceL = (double *)ValuePtr + tree_data->NodeOffset3[t][i][j][k];
            Value  = SliceL[L];

    }

    return(Value);
}

 
/*****  Hyb4_SmoothStepUp  *********************************************************/
/**
 *  Smoothing a payoff depending on the dimension
 */
int    Hyb4_SmoothStepUp(
                  double*        Slice,           /**< (O) Slice to be modified */
                  double*        Index,           /**< (I) Reference slice      */
                  double         Up,              /**< (I) Maximium level       */
                  double         Barrier,         /**< (I) Barrier level        */
                  int            SmoothingOn,     /**< (I) TRUE = smoothing on  */
                  int            t,               /**< (I) Current time period  */
                  int            tree_dim,        /**< (I) tree dimension       */
                  HYB4_TREE_DATA    *tree_data)      /**< (I) Tree data            */
{
    double  *SliceL;
    double  *IndexL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim)  */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)   */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)   */
    int     ***Top4, ***Bottom4;               /* Tree limits (4th dim)   */

    int     i, j, k, L;                        /* Node indices            */
    int     offset;                            /* Node offset             */
    double  Step;                              /* step size for smoothing */


    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];

    Step = 0.0;  /* initialise to smoothing off */
    
    switch (tree_dim)
    {
        case 1:
        {
            offset = tree_data->NodeOffset0[t];

            SliceL = Slice + offset;
            IndexL = Index + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                if (SmoothingOn)
                    Step = Hyb4_GetIndexStep(Index,1,i,0,0,0,t,tree_data);
                SliceL[i] = DrSmoothStep(Up, SliceL[i], IndexL[i], Barrier, Step);
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = tree_data->NodeOffset1[t][i];

                SliceL = Slice + offset;
                IndexL = Index + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    if (SmoothingOn)
                        Step = Hyb4_GetIndexStep(Index,2,i,j,0,0,t,tree_data);
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
                    offset = tree_data->NodeOffset2[t][i][j];

                    SliceL = Slice + offset;
                    IndexL = Index + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        if (SmoothingOn)
                            Step = Hyb4_GetIndexStep(Index,3,i,j,k,0,t,tree_data);
                        SliceL[k] = DrSmoothStep(Up, SliceL[k], IndexL[k], Barrier, Step);
                    }
                }
            }

            break;
        }        
        case 4:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];

                        SliceL = Slice + offset;
                        IndexL = Index + offset;

                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            if (SmoothingOn)
                                Step = Hyb4_GetIndexStep(Index,3,i,j,k,L,t,tree_data);
                            SliceL[L] = DrSmoothStep(Up, SliceL[L], IndexL[L], Barrier, Step);
                        }
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Hyb4_SmoothStepUp: invalid slice dimension "
                        "(must be 1, 2, 3 or 4)!");
            return(FAILURE);
        }
    } /* End of switch() */
    
    return (SUCCESS);

}

int Hyb4_MaxMinOnSlice(TSLICE         Slice,     /* (O) Slice to be modified  */
                     int            SliceDim,  /* (I) Dimension of above    */
                     double         Min,       /* (I) Minimum level         */
                     double         Max,       /* (I) Maximum level         */
                     int            t,         /* (I) Current time period   */
                     HYB4_TREE_DATA     *tree_data) /* (I) Tree data             */
{


    double
          *SliceL;


    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (3rd dim) */

        offset,
        i, j, k, L;                       /* Node indices                 */
       

        
        Top1    = tree_data->iMax[t];                    
        Bottom1 = tree_data->iMin[t];
        Top2    = tree_data->jMax[t];                    
        Bottom2 = tree_data->jMin[t];
        Top3    = tree_data->kMax[t];                    
        Bottom3 = tree_data->kMin[t];                    
        Top4    = tree_data->LMax[t];                    
        Bottom4 = tree_data->LMin[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = tree_data->NodeOffset0[t];
            SliceL = (double *)Slice  + offset;
            
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MAXMIN(SliceL[i],Max,Min);
            }
            break;

        case 2:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = tree_data->NodeOffset1[t][i];
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
                    offset = tree_data->NodeOffset2[t][i][j];
                    SliceL = (double *)Slice  + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MAXMIN(SliceL[k],Max,Min);
                    }
                }
            }
            break;

        case 4:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];
                        SliceL = (double *)Slice  + offset;      

                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            SliceL[L] = MAXMIN(SliceL[L],Max,Min);
                        }
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2, 3 or 4)."
                    " (Hyb4_MaxMinOnSlice)\n");
            return(FAILURE);

    } /* End of switch() */


    return (SUCCESS);


}  /*  Hyb4_MaxMinOnSlice  */

int Hyb4_MinOnSlice(TSLICE          Slice,     /* (O) Slice to be modified  */
                    int             SliceDim,  /* (I) Dimension of above    */
                    double          Min,       /* (I) Minimum level         */
                    int             t,         /* (I) Current time period   */
                    HYB4_TREE_DATA *tree_data) /* (I) Tree data             */
{


    double
          *SliceL;


    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */
       ***Top4,  ***Bottom4,           /* Limits of the tree (3rd dim) */

        offset,
        i, j, k, L;                       /* Node indices                 */
       

        
        Top1    = tree_data->iMax[t];                    
        Bottom1 = tree_data->iMin[t];
        Top2    = tree_data->jMax[t];                    
        Bottom2 = tree_data->jMin[t];
        Top3    = tree_data->kMax[t];                    
        Bottom3 = tree_data->kMin[t];                    
        Top4    = tree_data->LMax[t];                    
        Bottom4 = tree_data->LMin[t];                    


    switch (SliceDim)
    {
        case 1:

            offset = tree_data->NodeOffset0[t];
            SliceL = (double *)Slice  + offset;
            
            for (i=Bottom1; i<=Top1; i++)
            {
                SliceL[i] = MAX(SliceL[i],Min);
            }
            break;

        case 2:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = tree_data->NodeOffset1[t][i];
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
                    offset = tree_data->NodeOffset2[t][i][j];
                    SliceL = (double *)Slice  + offset;      

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        SliceL[k] = MAX(SliceL[k],Min);
                    }
                }
            }
            break;

        case 4:
            
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];
                        SliceL = (double *)Slice  + offset;      

                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            SliceL[L] = MAX(SliceL[L],Min);
                        }
                    }
                }
            }
            break;

        default:
            DR_Error("Invalid slice dimension (must be 1, 2, 3 or 4)."
                    " (Hyb4_MinOnSlice)\n");
            return(FAILURE);

    }


    return SUCCESS;


}  /*  Hyb4_MinOnSlice  */
