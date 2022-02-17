#ifndef _DRLMATRIXO_H 
#define _DRLMATRIXO_H

/* These are in ~clibdev/clutil. Lib is libclutil.a/clutil.lib */
#include "varapi.h"            /* TVar */ 
#include "objman.h"            /* TObjectType, TObjectManager */
#include "objapi.h"            /* TEncBuffer, TDecBuffer */

#include "bastypes.h"             /* TMatrix2D ... */ 

extern  TObjectType MAT_OBJECT_TYPE; 

/*------------------------------------------------------------------- 
 *                        TYPEDEFS 
 *-------------------------------------------------------------------  
 */


/*------------------------------------------------------------------- 
 *                  PROTOTYPES FOR THE METHODS 
 *------------------------------------------------------------------- 
 */

GTO_EXPORT(TMatrix2D *)  MatrixNew (
    long         numRows,      /* (I)  */ 
    long         numCols,      
    double     * data);        /*      */ 
/*
   GtoMatrixNew 
*/ 


GTO_EXPORT(TMatrix2D *)  MatrixResize (
    long         numRows,      /* (I)  */ 
    long         numCols,      
    long         numNewRows,      
    double     * data);        /*      */ 


GTO_EXPORT(TMatrix2D *) MatrixMakeCopy (
        TMatrix2D  * copyMe);

/*
   GtoMatrixNew (copyMe->numDim1, copyMe->numDim2, copyMe->data); 
*/ 

GTO_EXPORT(void)        MatrixDelete (
        TMatrix2D     * deleteMe); 

GTO_EXPORT(void)         MatrixPrint (
        TMatrix2D     * printMe);

/*
   GtoMatrixFree
*/ 

GTO_EXPORT(TVar *)      MatrixGetField (
        TMatrix2D     * matrix, 
        char          * filedName); 


/* 
    You have to write 
*/ 

GTO_EXPORT(int)         MatrixRegister ( 
        TObjectManager * om); 


GTO_EXPORT(int)         MatrixEncode (
        TEncBuffer *eb,
        char       *objName,
        TMatrix2D  *mat);

GTO_EXPORT(TMatrix2D*) MatrixDecode (
        TDecBuffer *db,
        char *objName);



#endif /* _DRLMATRIXO_H */
