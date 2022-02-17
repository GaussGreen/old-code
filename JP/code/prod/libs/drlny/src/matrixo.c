#include <string.h>
#include "macros.h"    /* FREE, FREE_ARRAY */
#include "cerror.h"    /* GtoErrMsg */ 
#include "gtomat.h"    /* GtoMatrixNew */ 
#include "classreg.h"  /* GtoRegisterClasses */

#include "drlmatrixo.h"


/*-------------------------------------------------------------------
 *                CLASS REGISTRATION RECORD 
 *-------------------------------------------------------------------
 */ 

TObjectType MAT_OBJECT_TYPE =
{ "MAT",                 		/* class short name */ 
  NULL,                  	      	/* no base class */
  (TObjectFreeFunc*)    MatrixDelete, 	/* destructor */
  NULL,                  		/* no base class to get */ 
  (TObjectGetFunc*)     MatrixGetField, /* get method */ 
  NULL,                  		/* no coerceTo function  */ 
  NULL,                  		/* no coerceFrom */  
  NULL,                  		/* no put function */ 
  (TObjectEncodeFunc*)  MatrixEncode,   /* encode function */
  (TObjectDecodeFunc*)  MatrixDecode,   /* decode function */
  (TObjectCopyFunc*)    MatrixMakeCopy 	/* copy function */
};


GTO_EXPORT(int)         MatrixRegister ( 
        TObjectManager *om) 
{ 
    static char         routine[]="MatrixRegister"; 
    int                 status=FAILURE; 
    
    if (om IS NULL) 
    {
        GtoErrMsg ("%s: input object manager is NULL\n",routine); 
        goto     done; 
    } 
    
    if (GtoRegisterClasses(om) ISNT SUCCESS ||
        GtoRegisterObjectType (om, &MAT_OBJECT_TYPE) ISNT SUCCESS)
    {
        goto done; /* failure */
    }


    status=SUCCESS; 
    
done: 
    if (status IS FAILURE) 
        GtoErrMsg ("%s: Failed to register class %s\n", 
                   routine, "TMatrix2D"); 
    return (status) ;
}


GTO_EXPORT(int) MatrixEncode (
     TEncBuffer *eb, char *objName, TMatrix2D *mat)
{
    return SUCCESS;
}



GTO_EXPORT(TMatrix2D*) MatrixDecode (
    TDecBuffer *db, char *objName)
{
    return NULL;
}



GTO_EXPORT(void) MatrixPrint (TMatrix2D* mat)
{
    int i, j;
    static char         routine[]="MatrixPrint"; 
    GtoErrMsg ("%s:Matrix has fields:numDim1:%dnumDim2:%d\n",
		routine, mat->numDim1, mat->numDim2);
    for ( i = 0; i < mat->numDim1; ++i)
         for (j= 0; j < mat->numDim2; ++j)
             GtoErrMsg ("element [%d][%d] :  %f\n", i,j, mat->data[i][j]);
}



/*------------------------------------------------------------------
	Convert 1-D array to 2-D matrix
*/

GTO_EXPORT(TMatrix2D *) MatrixNew (
    long         numRows,     
    long         numCols,      
    double      *data)       
{ 
    static char         routine[]="MatrixNew"; 
    int status = FAILURE;

    TMatrix2D  *retval = NULL; 
    long         i, j; 
    
    retval   = GtoMatrixNewEmpty (numRows, numCols); 
    if (retval IS NULL) 
        goto done; 

    /* 
    ** Default assumption is that matrix will be "rowized" 
    ** by Excel. 
    */
    
    for (i=0; i < numRows; i++) 
        for (j=0; j < numCols; j++)
            retval->data[i][j]  = data[numCols*i + j]; 
    
    status = SUCCESS;

  done: 
    if (status IS FAILURE) 
        GtoErrMsg ("%s: Failed to construct class %s\n", 
                   routine, "TMatrix2D"); 
    return (retval); 
}




/*---------------------------------------------------------------------
	Convert to a matrix with row number = numNewRows, and column number
	to be determined so that it will accommadate the original matrix
	size, numRows x numCols.
*/
GTO_EXPORT(TMatrix2D *) MatrixResize (
    long         numRows,    
    long         numCols,   
    long	 numNewRows,
    double     * data)     
{ 
    static char         routine[]="MatrixResize"; 
    int status = FAILURE;
    TMatrix2D  * retval = NULL; 
    long         i, j, i_old, j_old; 
    long	 numNewCols;

    if (numNewRows >= numRows)
     	numNewCols = numCols;
    else
	numNewCols = numCols*(long)ceil((double)numRows/(double)numNewRows);
    
    retval = GtoMatrixNewEmpty (numNewRows, numNewCols); 
    if (retval IS NULL) 
        goto done; 

    /* 
    ** Default assumption is that matrix will be "rowized" 
    ** by Excel. 
    */
    
    if (numNewRows >= numRows) {
    	for (i=0; i < numRows; i++) 
            for (j=0; j < numCols; j++)
            retval->data[i][j]  = data[numCols*i + j]; 
    }
    else { 
	for (i=0; i < numNewRows; i++) 
            for (j=0; j < numNewCols; j++){
	        i_old = numNewRows*(long)floor((double)j/(double)numCols) + i;
	        j_old = j - numCols*(long)floor((double)j/(double)numCols);
                retval->data[i][j] = (i_old >= numRows)?
				      0:
				      data[numCols*i_old + j_old]; 
		}
    }	
    status = SUCCESS;

  done: 
    if (status IS FAILURE) 
        GtoErrMsg ("%s: Failed to construct class %s\n", 
                   routine, "TMatrix2D"); 
    return (retval); 
}



GTO_EXPORT(TMatrix2D *) MatrixMakeCopy (
        TMatrix2D  * copyMe)
{ 
    static char routine[]="MatrixMakeCopy";

    TMatrix2D  * retval = NULL; 

    retval =    GtoMatrixNew (copyMe->numDim1, copyMe->numDim2, copyMe->data);
    if (retval IS (TMatrix2D*)NULL)
        GtoErrMsg("%s: Failed.\n", routine);

    return (retval); 
} 


GTO_EXPORT(void)        MatrixDelete (
        TMatrix2D     * deleteMe)
{     
    GtoMatrixFree (deleteMe); 
} 



/*-----------------------------------------------------------------
	Get specified column or row of 2-D matrix
*/

GTO_EXPORT(TVar*)  MatrixGetField (
        TMatrix2D     *matrix, 
        char          *fieldName)
{ 
    static char         routine[]="MatrixGetField"; 
    int     status = FAILURE;
    TVar    *retval = NULL; 

    char    *pNum = NULL; 
    double  *data = NULL;
    long     colNum, rowNum; 
    long     k; 
    long     numField;

    char    *fieldStr[] = {
		"numrows",
    		"numcols",
    		"rowi,where i is the row number",
    		"coli,where i is the col number"
	    };

    numField = 4;
    
    /**
     ** If the field name is an empty string, then return all the field-names. 
     **/ 
    if (fieldName[0] IS '\0') {
	retval = GtoNewVarStringArray (numField, fieldStr);
        if (retval IS NULL) 
            goto done; 
    }
    else 
    { 
        if (strcmp (fieldName, fieldStr[0]) == 0)
        { 
            retval = GtoNewVarLong (matrix->numDim1); 
            if (retval IS NULL) 
                goto done; 
        } 
        else if (strcmp (fieldName, fieldStr[1]) == 0)
        { 
            retval = GtoNewVarLong (matrix->numDim2); 
            if (retval IS NULL) 
                goto done; 
        } 
        else if (strncmp(fieldName, "col", 3) == 0) 
        {
	    pNum = fieldName + 3;
            /**
             ** Parse the column number from the 
             ** input field-name
             **/
            if (pNum IS NULL) 
                goto done;	 /* failed */  
            colNum = atoi (pNum) - 1;
	    if (colNum > matrix->numDim2 - 1 ||
		colNum < 0)
	    {
		status = FAILURE;
                GtoErrMsg("The input column number %d is outside the "
			  "maximum range [1, %d].\n", 
			  colNum+1, matrix->numDim2);
		goto done;
	    }
                        
            data = NEW_ARRAY(double, matrix->numDim1); 
            if (data IS NULL) 
                goto done; /* failed */  

            for (k=0; k < matrix->numDim1; k++)
                data[k] = matrix->data[k][colNum]; 
            
            retval = GtoNewVarDoubleArray (matrix->numDim1, data); 

	    FREE(data);

            if (retval IS NULL) 
                goto done; 
	    
        } 
        else if (strncmp(fieldName, "row", 3) == 0) 
        {
            pNum = fieldName + 3; 

            /**
             ** Parse the column number from the 
             ** input field-name
             **/
            if (pNum IS NULL) 
                goto done;	 /* failed */  

            rowNum = atoi (pNum) - 1;

	    if (rowNum > matrix->numDim1 - 1 ||
		rowNum < 0)
	    {
		status = FAILURE;
                GtoErrMsg("The input row number %d is outside the "
			  "maximum range [1, %d].\n", 
			  rowNum+1, matrix->numDim1);
		goto done;
	    }
                        
            data = NEW_ARRAY(double, matrix->numDim2); 
            if (data IS NULL) 
                goto done; /* failed */  

            for (k=0; k < matrix->numDim2; k++)
                data[k] = matrix->data[rowNum][k]; 
            
            retval = GtoNewVarDoubleArray (matrix->numDim2, data); 

	    FREE(data);

            if (retval IS NULL) 
                goto done; 

        } 
    } 
    status = SUCCESS;

    done:
         if (status IS FAILURE)
            GtoErrMsg("%s: Failed.\n", routine);

         return retval;

} 

