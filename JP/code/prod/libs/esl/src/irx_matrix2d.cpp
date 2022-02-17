/*------------------------------------------------------------------
C SOURCE FILE:  matrix2d.c

CREATED BY:     06/24/93 K. Varikooty

PURPOSE:        Matrix Operations

$Header: /home/alibcvs/alibcvs/.cvsroot-alib/utils/src/matrix2d.c,v 1.39 2004/12/22 14:15:25 cbauer Exp $
---------------------------------------------------------------------- 
Proprietary software, whether developed for Morgan by in-house
staff or external contractors, must not be copied or removed from Morgan
premises without the approval of a Senior Vice President and Audit.

This proprietary software has been developed strictly for Morgan's own
internal use.  Any use or misuse, intentional or otherwise which contradicts
or places this policy in jeopardy is strictly forbidden.  Also, any actions or 
inactions resulting in a breach of this goal will be dealt with harshly.

Do Not Distribute this to anyone outside the Analytics Library
Group without authorization from its management. 

Copyright 1995 J.P. Morgan & Co. Incorporated.   All rights reserved.
-------------------------------------------------------------------------  */

#include "irx/matrix2d.h"

#include <stdio.h>

#include "irx/error.h"
#include "irx/macros.h"             /* NEW */
#include "irx/memutils.h"

/*-----------------------------------------------------------------------
FUNCTION:      irxMatrixGetValue

CREATED BY:   06/24/93 K. Varikooty

DESCRIPTION:  Returns the element of a matrix
---------------------------------------------------------------------- */
int   irxMatrixGetValue
   (const IrxTMatrix2D *theMatrix,               /* (I) */
    int        idx1,                    /* (I) */
    int        idx2,                    /* (I) */
    double     *theValue                /* (O) */
   )
{
    static char routine[] = "irxMatrixGetValue";
    int         status    = FAILURE;

    if ((theMatrix == NULL) ||
        (theValue == NULL))
    {
        irxError("%s: NULL inputs.\n", routine);
        goto RETURN;
    }

    if ((idx1 < 0) || (idx1 >= theMatrix->numDim1) ||
        (idx2 < 0) || (idx2 >= theMatrix->numDim2))
    {
        irxError("%s: Element requested (%d,%d) is outside matrix "
                  "range %dx%d.\n",
                  routine, idx1+1, idx2+1,
                  theMatrix->numDim1, theMatrix->numDim2);
        goto RETURN;
    }

    status = SUCCESS;
    *theValue = theMatrix->data[idx1][idx2];

 RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return (status);
}

/*-----------------------------------------------------------------------
FUNCTION:      irxMatrixCopy

CREATED BY:   11/28/94 D. Gallager

DESCRIPTION:  Allocates a new matrix, and copies a matrix values into it.
---------------------------------------------------------------------- */
IrxTMatrix2D *   irxMatrixCopy(
    const IrxTMatrix2D *matrix)                    /* (I) */
{
    static char routine[] = "irxMatrixCopy";
    int         status    = FAILURE;

    IrxTMatrix2D *newMat = NULL;

    if (matrix == NULL)
    {
        irxError("%s: Cannot copy NULL matrix.\n", routine);
        goto RETURN;
    }

    newMat = irxMatrixNew(matrix->numDim1, matrix->numDim2, 
                          matrix->data);
    if (newMat == NULL)
        goto RETURN;

    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
    {
        irxMatrixFree(newMat);
        newMat = NULL;
        irxErrorFailure(routine);
    }
    return(newMat);
}

/*f
***************************************************************************
** FUNCTION: irxMatrixCopyData
** AUTHOR:   Peter Taylor (7 October 1998)
**
** Sets the values of a matrix from another matrix. Checks that the
** dimensions of the two matrices are the same.
***************************************************************************
*/
int  irxMatrixCopyData
(const IrxTMatrix2D *src,   /* (I) */
 IrxTMatrix2D *dst    /* (O) */
)
{
    static char routine[] = "irxMatrixCopyData";

    int idx1; /* Iterator for the first dimension */
    int idx2; /* Iterator for the second dimension */

    if (src->numDim1 != dst->numDim1 ||
        src->numDim2 != dst->numDim2)
    {
        irxError ("%s: Matrix size mismatch. Source = [%d,%d]. "
                   "Target = [%d,%d].\n",
                   routine,
                   src->numDim1, src->numDim2,
                   dst->numDim1, dst->numDim2);
        return irxErrorFailure (routine);
    }

    for (idx1 = 0; idx1 < src->numDim1; ++idx1)
    {
        for (idx2 = 0; idx2 < src->numDim2; idx2++)
        {
            dst->data[idx1][idx2] = src->data[idx1][idx2];
        }
    }
    return (SUCCESS);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrix2DTo1DNew

CREATED BY:   08/18/99 C. Irwin

DESCRIPTION:  Creates a 1D array from 2D data.
---------------------------------------------------------------------- */
double * irxMatrix2DTo1DNew(
    const IrxTMatrix2D *matrix)   /* (I) Matrix */
{
    static char routine[] = "irxMatrix2DTo1DNew";
    int         status    = FAILURE;

    long    numRows;
    long    numCols;
    long    i;
    long    j;
    double *matrixData = NULL; /* 1D array to return */

    if (matrix == NULL)
    {
        irxError("%s: Cannot process NULL inputs.\n", routine);
        goto RETURN;
    }

    numRows = matrix->numDim1;
    numCols = matrix->numDim2;

    matrixData = NEW_ARRAY(double, numRows * numCols);
    if (matrixData == NULL)
        goto RETURN;

    for (i=0; i < numRows; i++) 
    {
        for (j=0; j < numCols; j++)
        {
            matrixData[numCols*i + j] = matrix->data[i][j];
        }
    }


    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
    {
        FREE(matrixData);
        matrixData = NULL;
        irxErrorFailure(routine);
    }
    return (matrixData);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrix1DTo2DNew

CREATED BY:   08/18/99 C. Irwin

DESCRIPTION:  Creates a 2D matrix from 1D data.
---------------------------------------------------------------------- */
IrxTMatrix2D * irxMatrix1DTo2DNew(
    long    numRows,     /* (I) Number of rows */
    long    numCols,     /* (I) Number of columns */
    double *matrixData)  /* (I) Matrix data */
{
    static char routine[] = "irxMatrix1DTo2DNew";
    int         status    = FAILURE;

    IrxTMatrix2D *retval = NULL;

    long        i;
    long        j; 
    
    if (matrixData == NULL)
    {
        irxError("%s: Cannot process NULL inputs.\n", routine);
        goto RETURN;
    }

    retval = irxMatrixNewEmpty (numRows, numCols); 
    if (retval == NULL) 
        goto RETURN; 

    /* 
    ** Default assumption is that matrix will be "rowized" 
    ** by Excel. 
    */
    
    for (i=0; i < numRows; i++) 
    {
        for (j=0; j < numCols; j++)
        {
            retval->data[i][j]  = matrixData[numCols*i + j]; 
        }
    }

    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
    {
        irxMatrixFree(retval);
        retval = NULL;
        irxErrorFailure(routine);
    }
    return (retval);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixDiagonalNew

CREATED BY:   8/23/99 Charles Irwin

DESCRIPTION:  Creates a diagonal matrix with data as its diagnonal elements
---------------------------------------------------------------------- */
IrxTMatrix2D * irxMatrixDiagonalNew(
    int     n,        /* (I) Size of diagonal */
    double *data)     /* (I) Diagonal elements */
{
    static char routine[] = "irxMatrixDiagonalNew";
    int         status    = FAILURE;

    IrxTMatrix2D *retVal = NULL;
    int i;
    int j;

    retVal = irxMatrixNewEmpty(n, n);
    if (retVal == NULL)
        goto RETURN;

    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            if (i != j)
            {
                retVal->data[i][j] = 0.0;
            }
            else
            {
                retVal->data[i][i] = data[i];
            }
        }
    }

    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
    {
        irxMatrixFree(retVal);
        retVal = NULL;
        irxErrorFailure(routine);
    }
    return (retVal);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixNew

CREATED BY:   06/24/93 K. Varikooty

DESCRIPTION:  Creates a matrix with data as its elements
---------------------------------------------------------------------- */
IrxTMatrix2D *   irxMatrixNew
    (int      numDim1,                   /* (I) */
     int      numDim2,                   /* (I) */
     double **data                       /* (I) */
   )
{
    static char routine[]="irxMatrixNew";
    int    idx1;
    int    idx2;
    
    IrxTMatrix2D *mtx = NULL;
    double **mtxDatap;          /* Pointer to data in matrix */

    /* Allocate a new empty matrix.
     */
    mtx = irxMatrixNewEmpty(numDim1, numDim2);
    if (mtx == NULL)
        goto RETURN;              /* Failed */


    /* Insert data in matrix.
     */
    for (idx1 = 0, mtxDatap = mtx->data;
         idx1 < numDim1; idx1++)
    {
        for (idx2 = 0; idx2 < numDim2; idx2++)
        {
            mtxDatap[idx1][idx2] = data[idx1][idx2];
        }
    }

  RETURN:
    if (mtx == NULL)
        irxError("%s: Failed.\n", routine);

    return(mtx);
}


/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixNewEmpty

CREATED BY:   06/24/93 K. Varikooty

DESCRIPTION:  Creates an empty matrix
---------------------------------------------------------------------- */
IrxTMatrix2D *   irxMatrixNewEmpty
    (int  numDim1,                       /*  (I) */
     int  numDim2                        /*  (I) */
     )
{
    IrxTMatrix2D *mtx =  NULL;
    static char routine[]="irxMatrixNewEmpty";
    int status = FAILURE;       /* Until proven successful */
    
    if (numDim1 < 1 || numDim2 < 1)
    {
        irxError("%s: Number dim1 (%d) or num dim2 (%d) out of range.\n",
                  routine, numDim1,  numDim2);
        goto RETURN;               /* Failed */
    }
    
    mtx = NEW(IrxTMatrix2D);
    if (mtx == NULL)
        goto RETURN;

    mtx->numDim1 = numDim1;
    mtx->numDim2 = numDim2;
    mtx->data = (double **)irxArray2DNew((int)numDim1, (int)numDim2, 
                                         sizeof(double));
    if (mtx->data == NULL)
        goto RETURN;
    status = SUCCESS;

  RETURN:
    if (status == FAILURE)
    {
        irxError("%s: Failed.\n", routine);
        irxMatrixFree(mtx);
        return (NULL);
    }
    
    return(mtx);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixFree

CREATED BY:   06/24/93 K. Varikooty

DESCRIPTION:  Frees the memory associated with a matrix
---------------------------------------------------------------------- */

void   irxMatrixFree
    (IrxTMatrix2D *mtx               /* (I) */
     )
{
    if (mtx != NULL)
        FREE((void **)mtx->data);
    FREE(mtx);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixAdd

CREATED BY:   9/15/99 Charles Irwin

DESCRIPTION:  Adds two matrices together.
---------------------------------------------------------------------- */
IrxTMatrix2D * irxMatrixAdd(
    const IrxTMatrix2D *A,
    const IrxTMatrix2D *B)
{
    static char routine[] = "irxMatrixAdd";
    int         status    = FAILURE;

    IrxTMatrix2D *C = NULL;
    int i,j;

    if ((A == NULL) || (B == NULL))
    {
        irxError("%s: NULL inputs.\n", routine);
        goto RETURN;
    }

    if ((A->numDim1 != B->numDim1) ||
        (A->numDim2 != B->numDim2))
    {
        irxError("%s: Matrix sizes (%dx%d) and (%dx%d) do not match.\n",
                  routine, A->numDim1, A->numDim2,
                  B->numDim1, B->numDim2);
        goto RETURN;
    }

    C = irxMatrixNewEmpty(A->numDim1, A->numDim2);
    if (C == NULL)
        goto RETURN;

    for (i=0; i<A->numDim1; i++)
        for (j=0; j<A->numDim2; j++)
        {
            C->data[i][j] = A->data[i][j] + B->data[i][j];
        }

    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
    {
        irxMatrixFree(C);
        C = NULL;
        irxErrorFailure(routine);
    }
    return (C);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixScalarMultiply

CREATED BY:   8/23/99 Charles Irwin

DESCRIPTION:  Multiplies a scalar times a matrix.
---------------------------------------------------------------------- */
IrxTMatrix2D *   irxMatrixScalarMultiply(
    double     scalar,            /* (I) Scalar */
    const IrxTMatrix2D *matrix)            /* (I) Matrix */
{
    static char routine[] = "irxMatrixScalarMultiply";
    int         status    = FAILURE;

    IrxTMatrix2D *retVal = NULL;
    int i;
    int j;

    if (matrix == NULL)
    {
        irxError("%s: Cannot process NULL inputs.\n", routine);
        goto RETURN;
    }

    retVal = irxMatrixCopy(matrix);
    if (retVal == NULL)
        goto RETURN;

    for (i=0; i<retVal->numDim1; i++)
    {
        for (j=0; j<retVal->numDim2; j++)
        {
            retVal->data[i][j] = scalar * matrix->data[i][j];
        }
    }

    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
    {
        irxMatrixFree(retVal);
        retVal = NULL;
        irxErrorFailure(routine);
    }
    return (retVal);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixKroneckerMultiply

CREATED BY:   1/5/00 Charles Irwin

DESCRIPTION:  Multiplies matrices a and b and returns the TENSOR product
---------------------------------------------------------------------- */
IrxTMatrix2D * irxMatrixKroneckerMultiply(
    const IrxTMatrix2D *a,
    const IrxTMatrix2D *b)
{
    static char routine[] = "irxMatrixKroneckerMultiply";
    int         status    = FAILURE;

    IrxTMatrix2D *c = NULL;

    int m,n,o,p,i,j,k,l;
    
    m = a->numDim1;
    n = a->numDim2;
    o = b->numDim1;
    p = b->numDim2;

    c = irxMatrixNewEmpty(m*o, n*p);
    if (c == NULL)
        goto RETURN;

    for (i=0; i<m; i++)
    {
        for (j=0; j<n; j++)
        {
            for (k=0; k<o; k++)
            {
                for (l=0; l<p; l++)
                {
                    c->data[i*o + k][j*p + l] = a->data[i][j]*b->data[k][l];
                }
            }
        }
    }

    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
    {
        irxMatrixFree(c);
        c = NULL;
        irxErrorFailure(routine);
    }
    return c;

}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixMultiply

CREATED BY:   06/24/93 K. Varikooty/D. Gallager

DESCRIPTION:  Multiplies matrix a and b and returns the product
---------------------------------------------------------------------- */

IrxTMatrix2D *   irxMatrixMultiply
    (const IrxTMatrix2D *a,                /* (I) */
     const IrxTMatrix2D *b                 /* (I) */
   )
{
    static char routine[]="irxMatrixMultiply";
    int status = FAILURE;
    
    IrxTMatrix2D *c = NULL;
    double   sum;                        /* Sum for dot product */
    int      idx1;
    int      idx2;
    int      idxSum;
    int      numDim1;
    int      numDim2;
    int      numCom;
                                         /* = # cols in A = # rows in B */
    double  **prodData;                  /* Product matrix values */
    double  **aData;                     /* A matrix values */
    double  **bData;                     /* B matrix values */
    double  *aDataTmp;                   /* Tmp for A */

#ifdef DEBUG
   irxMatrixPrint(a, " A ");
   irxMatrixPrint(b, " B ");
#endif   
   
   /* Error checking.
    */
   if (a->numDim1 < 1 || a->numDim2 < 1 ||
       b->numDim1 < 1 || b->numDim2 < 1)
   {
       irxError("%s: Matrices must have dim at least 1 by 1.\n", routine);
       goto RETURN;                       /* Failed */
   }
       
   if (a->numDim2 != b->numDim1)
   {
       irxError("%s: # Cols matrix A (%d) != # Rows matrix B (%d)\n",
                 routine, a->numDim2, b->numDim1);
       goto RETURN;                       /* Failed */
   }


   c = irxMatrixNewEmpty(a->numDim1,b->numDim2);
   if (c == (IrxTMatrix2D *)NULL)
       goto RETURN;                       /* Failed */


   /* Preassign for efficiency.
    */
   numCom   = a->numDim2;
   numDim1  = c->numDim1;
   numDim2  = c->numDim2;
   prodData = c->data;               
   aData    = a->data;
   bData    = b->data;

   /* For each row of product matrix...
    */   
   for (idx1 = 0; idx1 < numDim1; idx1++, aData++)
   {
      /* For each column of product matrix...
       */
      for (idx2 = 0; idx2 < numDim2; idx2++)
      {
          sum = 0.;             /* Initialize next element */
          aDataTmp = *aData;

          /* For each element in the dot product...
           */
          for (idxSum = 0; idxSum < numCom; idxSum++)
          {
               sum += *aDataTmp++ * bData[idxSum][idx2];
          }
          prodData[idx1][idx2] = sum;
      }                                 /* for idx2 */
   }                                    /* for idx1 */
   status = SUCCESS;

#ifdef DEBUG
   irxMatrixPrint(c, " C ");
#endif   
   
RETURN:
    if (status == FAILURE)
    {
        irxMatrixFree(c);
        c = NULL;
        irxError("%s: Failed.\n", routine);
    }
    return (c);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixMultiplyVectors

CREATED BY:   11/28/94 D. Gallager

DESCRIPTION:  Multiplies two identically dimensioned
              vectors and places outputs in third vector.
              Note that the vectors have to be the same length,
              but not vary in the same dimension. 
---------------------------------------------------------------------- */

int   irxMatrixMultiplyVectors
  (const IrxTMatrix2D *a,                          /* (I) */
   const IrxTMatrix2D *b,                          /* (I) */
   IrxTMatrix2D *c                           /* (O) */
   )
{
    static char routine[]="irxMatrixMultiplyVectors";
    int status  = FAILURE;              /* Until proven successful */
    IrxTBool aDim1;                     /* If Dim1 varies for matrices a+b */
    IrxTBool cDim1;                     /* If Dim1 varies for matrix c */
    int      dimAB;                     /* Dimension of vector a+b */
    int      dimC;                      /* Dimension of vector c */
    int      idx;

    if (a->numDim1 != b->numDim1 || a->numDim2 != b->numDim2)
    {
        irxError("%s: Matrices a,b must have the same dimensions.\n",
                  routine);
        goto RETURN;                      /* Failed */
    }

    if ((a->numDim1 != 1 && a->numDim2 != 1) ||
        (c->numDim1 != 1 && c->numDim2 != 1))
    {
        irxError("%s: Matrices a,b,c must be vectors.\n",
                  routine);
        goto RETURN;                      /* Failed */
    }

    cDim1 = c->numDim1 == 1 ? FALSE : TRUE;
    aDim1 = a->numDim1 == 1 ? FALSE : TRUE;
    dimAB = aDim1 ? a->numDim1 : a->numDim2;
    dimC = cDim1 ? c->numDim1 : c->numDim2;
    if (dimC != dimAB)
    {
        irxError("%s: Length of all vectors must be the same.\n", routine);
        goto RETURN;
    }

    for (idx = 0; idx < dimAB; idx++)
    {
        c->data [cDim1?idx:0] [cDim1?0:idx] = 
            a->data [aDim1?idx:0] [aDim1?0:idx]  
                * b->data [aDim1?idx:0] [aDim1?0:idx];
    }

    status = SUCCESS;

RETURN:
    if (status == FAILURE)
    {
        irxError("%s: Failed.\n", routine);
    }
    return (status);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixTranspose

CREATED BY:   06/24/93 K. Varikooty

DESCRIPTION:  Returns the transpose of matrix a
---------------------------------------------------------------------- */

IrxTMatrix2D *   irxMatrixTranspose
    (const IrxTMatrix2D *mtx               /* (I) Input matrix*/
     )
{
    static char routine[]="irxMatrixTranspose";
    IrxTMatrix2D *tMtx;                     /* Transposed matrix */
    double **tData;                      /* Transposed data */
    double **oData;                      /* Original data */
    int     idx1;                        /* Transposed row index */
    int     idx2;                        /* Transposed column index */
    int     dim1;                        /* # transposed dim1  = # orig dim2*/
    int     dim2;                        /* # transposed dim2 = # orig rows */
    
    tMtx = irxMatrixNewEmpty(mtx->numDim2, mtx->numDim1);
    if (tMtx == NULL)
        goto RETURN;                       /* Failed */
    
    tData = tMtx->data;                  /* Transposed data */
    oData = mtx->data;                   /* Original data */
    
    dim1  = tMtx->numDim1;               /* Preset for efficiency */
    dim2  = tMtx->numDim2;               /* Preset for efficiency */
    
    for (idx1 = 0; idx1 < dim1; idx1++)
    {
        for (idx2 = 0; idx2 < dim2; idx2++)
        {
            tData[idx1][idx2] = oData[idx2][idx1];
        }
    }

 RETURN:
    if (tMtx == NULL)
        irxError("%s: Failed.\n", routine);

   return(tMtx);
}




/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixPrint

CREATED BY:   06/24/93 K. Varikooty

DESCRIPTION:  prints out a matrix of name  matrixName
---------------------------------------------------------------------- */

void   irxMatrixPrint(
   const IrxTMatrix2D *aMatrix,                  /* (I) */
   char *matrixName                     /* (I) */
   )
{
   int  idx2, idx1;

   if (aMatrix == NULL)
   {
       irxError("Matrix %s\n", matrixName);
       irxError("NULL\n");
   }
   else
   {
       irxError("Matrix %s Size: %dx%d\n", matrixName, aMatrix->numDim1,
                 aMatrix->numDim2);

       for (idx1 = 0; idx1 < aMatrix->numDim1; idx1++)
       {
           for (idx2 = 0; idx2 < aMatrix->numDim2; idx2++)
           {
               irxError("%f\t", aMatrix->data[idx1][idx2]);
           }
           irxError("\n");
       }
   }
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixMultiply3

CREATED BY:   06/24/93 K. Varikooty/D. Gallager

DESCRIPTION:  Computes the product of 3 matrices

---------------------------------------------------------------------- */

IrxTMatrix2D *   irxMatrixMultiply3(
   const IrxTMatrix2D *matrix1,                    /* (I) */
   const IrxTMatrix2D *matrix2,                    /* (I) */
   const IrxTMatrix2D *matrix3                     /* (I) */
   )
{
   IrxTMatrix2D *matrix12 = (IrxTMatrix2D *)NULL; /* Product of 1 & 2 */
   IrxTMatrix2D *matrix123 = (IrxTMatrix2D *)NULL; /* Product of 1,2,3 */

   static char routine[]="irxMatrixMultiply3";
   int status = FAILURE;                    /* Until proven successful */


   matrix12 = irxMatrixMultiply(matrix1, matrix2);
   if (matrix12 == (IrxTMatrix2D *)NULL)
       goto RETURN;                       /* Failed */

   matrix123 = irxMatrixMultiply(matrix12, matrix3);
   if (matrix123 == (IrxTMatrix2D *)NULL)
       goto RETURN;                       /* Failed */
   status = SUCCESS;


 RETURN:
   if (status == FAILURE)
       irxError("%s:  Failed.\n", routine);
   irxMatrixFree(matrix12);
   return(matrix123);
}


/*-----------------------------------------------------------------------
FUNCTION:    irxMatrix2DWriteWrap

CREATED BY:  D. Gallager

DESCRIPTION:     Writes a IrxTMatrix2D to the wrapper interface.
   Unless transpose is specified, dim2 varies fastest in the output
   range. Typically the software which sits between the wrapper and 
   the spreadsheet will read/write data to the spreadsheet row by row. 
   Thus typically, unless transpose is set, the dim1 index specifies
   the row, and the dim2 index specifies the column.
   Returns SUCCESS/FAILURE.
---------------------------------------------------------------------- */
int irxMatrix2DWriteWrap
   (const IrxTMatrix2D  *matrix,         /* (I) Matrix */
    char       *matName,        /* (I) Name of matrix, for err msgs */
    IrxTBool    transpose,      /* (I) T=perform transpose */
    double     *wrapOut)        /* (O) Wrapper output array */
{
    static char routine[] = "irxMatrix2DWriteWrap";
    int         status = FAILURE;
    int         i, j;
    int         dim1, dim2;


    dim1 = matrix->numDim1;
    dim2 = matrix->numDim2;

    /* Make sure we have enough room.
     */
    if (wrapOut[0] < dim1*dim2)
    {
        irxError("%s: Only %d output cells provided for output %s matrix.\n"
                  "\t%d = %dx%d are required.\n",
                  routine, (int)wrapOut[0], matName,
                  dim1*dim2, dim1, dim2);
        goto RETURN;
    }

    /* Indicate how many vols are being returned if more than enough
     * space was provided.
     */
    wrapOut[0] = dim1*dim2;

    /* Write matrix values. This could be made more efficient.
     */
    for (i=0; i<dim1; i++) 
        for (j=0; j<dim2; j++) 
        {
            if (transpose)
                wrapOut[j*dim1 + i + 1] =  matrix->data[i][j];
            else
                wrapOut[i*dim2 + j + 1] =  matrix->data[i][j];
        }

            /* made it through OK */
    status = SUCCESS;

RETURN:
    if (status != SUCCESS) 
    {
        irxError("%s: Failed.\n", routine);
    }
    return(status);
}

/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixInnerProduct

CREATED BY:   Chris Bauer (2nd July 2004)

DESCRIPTION:  Computes the inner product AxMxC where A and C are
              vectors and M is a matrix. 

---------------------------------------------------------------------- */

int  irxMatrixInnerProduct
(double    *vectorA,              /* (I) LHS vector A */
 long       sizeOfA,              /* (I) Size of vector A */ 
 const IrxTMatrix2D *matrixB,              /* (I) Matrix B */
 double    *vectorC,              /* (I) RHS vector C */
 long       sizeOfC,              /* (I) Size of vector A */
 double    *result                /* (O) Inner product AxBxC */
)
{
    static char routine[]="irxMatrixInnerProduct";
    int         status = FAILURE;
    double      sum = 0.0;
    double    **matrix = NULL;
    int         i;
    int         j;

    REQUIRE(vectorA != NULL);
    REQUIRE(matrixB != NULL);
    REQUIRE(vectorC != NULL);
    REQUIRE(result != NULL);

    REQUIRE(sizeOfA == matrixB->numDim1);
    REQUIRE(sizeOfC == matrixB->numDim2);

    matrix = matrixB->data;
    for (i = 0; i < sizeOfA; i++)
    {
        for (j = 0; j < sizeOfC; j++)
        {
            sum += vectorA[i] * matrix[i][j] * vectorC[j];
        }
    }

    status = SUCCESS;
    *result = sum;

 RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return(status);
}



/*-----------------------------------------------------------------------
FUNCTION:     irxMatrixTrace

CREATED BY:   Chris Bauer (22nd December 2004)

DESCRIPTION:  Computes the trace of a matrix

---------------------------------------------------------------------- */

int  irxMatrixTrace
(const IrxTMatrix2D *matrix,               /* (I) Matrix */
 double    *result                /* (O) Trace of matrix */
)
{
    static char routine[]="irxMatrixTrace";
    int         status = FAILURE;
    double      sum = 0.0;
    double    **data = NULL;
    int         size;
    int         i;

    REQUIRE (matrix != NULL);
    REQUIRE (result != NULL);
    REQUIRE (matrix->numDim1 == matrix->numDim2);

    size = matrix->numDim1;
    data = matrix->data;

    for (i = 0; i < size; i++)
    {
        sum += data[i][i];
    }

    status = SUCCESS;
    *result = sum;

 RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return(status);
}



