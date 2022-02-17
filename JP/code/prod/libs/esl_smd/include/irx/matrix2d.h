/*---------------------------------------------------------------------
HEADER FILE:    matrix2d.h

CREATED BY:     K. Varikooty

PURPOSE:        External interface for irxMatrix routines.

$Header: /home/alibcvs/alibcvs/.cvsroot-alib/utils/include/matrix2d.h,v 1.36 2004/12/22 14:14:56 cbauer Exp $
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
#ifndef IRX_MATRIX2D_H
#define IRX_MATRIX2D_H

#include "cgeneral.h"

/*t
 * Defines a two-dimensional matrix with data of type double.
 */
typedef struct _IrxTMatrix2D
{
    int      numDim1;                   /* Num elems for 1st index */
    int      numDim2;                   /* Num elems for 2nd index */
    double **data;                      /* Data in matrix */
} IrxTMatrix2D;



#ifdef __cplusplus
extern "C"
{
#endif


/*f
 * Returns the element of a 1x1 matrix
 */
extern int   irxMatrixGetValue
   (const IrxTMatrix2D *theMatrix,               /* (I) */
    int        idx1,                    /* (I) */
    int        idx2,                    /* (I) */
    double     *theValue                /* (O) */
   );

/*f
 * Creates an empty matrix
 */
extern IrxTMatrix2D *   irxMatrixNewEmpty(
   int numRows,      /*  (I) */
   int numColumns    /*  (I) */
   );

/*f
 * Creates a 1D array from a 2D matrix.
 */
double * irxMatrix2DTo1DNew(
    const IrxTMatrix2D *matrix);  /* (I) Matrix */

/*f
 * Creates a 2D matrix from a 1D array.
 */
IrxTMatrix2D * irxMatrix1DTo2DNew(
    long    numRows,     /* (I) Number of rows */
    long    numCols,     /* (I) Number of columns */
    double *matrixData); /* (I) Matrix data */


/*f
 * Creates a diagonal matrix with data as its diagnonal elements
 */
IrxTMatrix2D * irxMatrixDiagonalNew(
    int     n,        /* (I) Size of diagonal */
    double *data);    /* (I) Diagonal elements */

/*f
 * Creates a matrix with data as its elements
 */
extern IrxTMatrix2D *   irxMatrixNew(
   int numRows,      /* (I) */
   int numColumns,   /* (I) */
   double **data     /* (I) */
   );

#define freeMatrix irxMatrixFree /* Backward compatibility */

/*f
 * Frees the memory associated with a matrix
 */
extern void   irxMatrixFree(
  IrxTMatrix2D *aMatrix   /* (I) */
  );

#define prinIrxTMatrix2D irxMatrixPrint /* Backward compatibility */
/*f
 * prints out a matrix of name  matrixName
 */
extern void   irxMatrixPrint(
   const IrxTMatrix2D *aMatrix, /* (I) */
   char *matrixName  /* (I) */
   );

/*f
 * Adds two matrices together.
 */
IrxTMatrix2D * irxMatrixAdd(
    const IrxTMatrix2D *A,
    const IrxTMatrix2D *B);


/*f
 * Multiplies a scalar times a matrix.
 */
IrxTMatrix2D *   irxMatrixScalarMultiply(
    double     scalar,            /* (I) Scalar */
    const IrxTMatrix2D *matrix);           /* (I) Matrix */

/*f
 * Computes the tensor product of two matrices.
 */
IrxTMatrix2D * irxMatrixKroneckerMultiply(
    const IrxTMatrix2D *a,
    const IrxTMatrix2D *b);

/*f
 * Multiplies matrix a and b and returns the product
 */
extern IrxTMatrix2D *   irxMatrixMultiply(
   const IrxTMatrix2D *a,       /* (I) */
   const IrxTMatrix2D *b        /* (I) */
   );

/*f
 * Computes the product of 3 matrices
 */
extern IrxTMatrix2D *   irxMatrixMultiply3(
   const IrxTMatrix2D *matrix1,                    /* (I) */
   const IrxTMatrix2D *matrix2,                    /* (I) */
   const IrxTMatrix2D *matrix3                     /* (I) */
   );

/*f
 * Multiplies two identically dimensioned
 * vectors and places outputs in third vector.
 * Note that the vectors have to be the same length,
 * but not vary in the same dimension. 
 */
extern int   irxMatrixMultiplyVectors
  (const IrxTMatrix2D *a,                          /* (I) */
   const IrxTMatrix2D *b,                          /* (I) */
   IrxTMatrix2D *c                           /* (O) */
   );

/*f
 * Returns the transpose of matrix a
 */
extern IrxTMatrix2D *   irxMatrixTranspose(
    const IrxTMatrix2D *a        /* (I) */
   );

/*f
 * Allocates a new matrix, and copies a matrix values into it.
 */
extern IrxTMatrix2D *   irxMatrixCopy
  (const IrxTMatrix2D *matrix    /* (I) */
   );

/*f
 * Sets the values of a matrix from another matrix. Checks that the
 * dimensions of the two matrices are the same.
 */
extern int  irxMatrixCopyData
  (const IrxTMatrix2D *src,   /* (I) */
   IrxTMatrix2D *dst    /* (O) */
   );


/*f
 * Writes a IrxTMatrix2D to the wrapper interface.
 * Unless transpose is specified, dim2 varies fastest in the output
 * range. Typically the software which sits between the wrapper and 
 * the spreadsheet will read/write data to the spreadsheet row by row. 
 * Thus typically, unless transpose is set, the dim1 index specifies
 * the row, and the dim2 index specifies the column.
 * Returns SUCCESS/FAILURE.
 */
extern int irxMatrix2DWriteWrap
   (const IrxTMatrix2D  *matrix,         /* (I) Matrix */
    char       *matName,        /* (I) Name of matrix, for err msgs */
    IrxTBool    transpose,      /* (I) T=perform transpose */
    double     *wrapOut);       /* (O) Wrapper output array */

/*f
***************************************************************************
** This function computes the inner product AxMxC where A and C are
** vectors and M is a matrix. 
***************************************************************************
*/
int  irxMatrixInnerProduct
(double    *vectorA,              /* (I) LHS vector A */
 long       sizeOfA,              /* (I) Size of vector A */ 
 const IrxTMatrix2D *matrixB,              /* (I) Matrix B */
 double    *vectorC,              /* (I) RHS vector C */
 long       sizeOfC,              /* (I) Size of vector A */
 double    *result                /* (O) Inner product AxBxC */
);


/*f
***************************************************************************
** This function computes the trace of a matrix
***************************************************************************
*/
int  irxMatrixTrace
(const IrxTMatrix2D *matrix,               /* (I) Matrix */
 double    *result                /* (O) Trace of matrix */
);


#ifdef __cplusplus
}
#endif

#endif    /* IRX_MATRIX2D_H */
