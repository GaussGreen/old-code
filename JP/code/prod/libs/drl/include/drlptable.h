/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	PTABLE - Tables of Variable Type
 * File:	drlptable.h
 * Function:	Wrapper utility routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlptable_H
#define	_drlptable_H
#include "drlstd.h"	/* wrapper type def (DVType, FloatL, etc.) */

#include "drlsmat.h"	/* DSwopMat */

/*t-@CDOC(idxn=DVTable,catn=structdef)-----------------------------------
 * Pivot table data structure for elements of arbitrary type.
 */

typedef	struct	{
	int	numFields;		/* number of fields (columns) */
	DVType	*fieldType;		/* field type [numFields] */
	char	**fieldName;		/* field name [numFields] */
	int	*numFieldValues;	/* nb different val [numFields] */
	int	*numFieldValuesMax;	/* max size (realloc) [numFields] */
	void	**fieldValues;		/* field value [numFields][numValues] */

	int	numData;		/* current number of elements (rows) */
	int	numDataMax;		/* max size (for realloc) */
	void	**dataValue;		/* data [numFields][numData] */

	void	**tmpData;		/* working space value [numFields] */

} DVTable;
/*e*/

extern	DVTable*	DrlDVTableNew(int numFields, DVType *fieldType);
extern	DVTable*	DrlDVTableNewFpRead(FILE *fp, int numFields,
				DVType *fieldType);
extern	DVTable*	DrlDVTableNewFileRead(char *fnam, int numFields,
				DVType *fieldTypes);
extern	int		DrlDVTableFree(DVTable* that);
extern	int		DrlDVTableFpWrite(DVTable *that, FILE *fp);

extern	int		DrlDVTableFindFieldName(DVTable *that,
				char *nameField);

extern	int		DrlDVTableExtractMatrix(DVTable *that,
				int indexFieldX, int indexFieldY,
				int indexFieldP, void *valFieldP,
				int indexFieldXY,
				int *numItemsX, int *numItemsY,
				DDate baseDate,
				double ***matrix);
extern	int		DrlDVTableExtractMatrixN(DVTable *that,	
				char *nameFieldX, char *nameFieldY,
				char *nameFieldP, void *valFieldP,
				char *nameFieldXY,
				int *numItemsX, int *numItemsY,
				DDate baseDate,
				double ***matrix);
extern	int		DrlDVTableExtractDSwopMat(DVTable *that,
				char *nameFieldX, char *nameFieldY,
				char *nameFieldP, void *valFieldP,
				char *nameFieldXY,
				DDate baseDate,
        			long freq,
				DSwopMat **mat);

/*t-@CDOC(idxn=TCTable,catn=structdef)-----------------------------------
 * Structure for table of string elements.
 */

typedef	struct	{
	int	numRows;
	int	numCols;
	int	numRowsMax;
	int	numColsMax;
	char	***cells;		/* [row][vol] */
} DCSheet;


#define	DCSheet_END	INT_MAX
#define	DCSheet_START	0



extern	DCSheet* DrlDCSheetNewEmpty(
			int numRowsMax,		/* (I) max number of rows */
			int numColsMax);	/* (I) max number of cols */
extern	int	DrlDCSheetFree(DCSheet* that);

extern	int	DrlDCSheetSet(DCSheet* that,
			int row, int col, const char* value);
extern	char*	DrlDCSheetGet(DCSheet* that, int row, int col);
extern	int	DrlDCSheetGetVType(
			DCSheet* that,		/* (I) table */
			int row,		/* (I) row number */
			int col, 		/* (I) column number */
			DVType fieldType,	/* (I) element type */
			void *ptr);		/* (O) */

extern	int	DrlDCSheetFpRead(
			DCSheet **that,		/* (O) output sheet */
			FILE *fp);		/* (I) file pointer */
extern	int	DrlDCSheetFileRead(
			DCSheet **that,		/* (O) output sheet */
			const char *fnam);	/* (I) file name */

extern	int	DrlDCSheetFpWrite(
			DCSheet *that,		/* (I) table */
			FILE *fp);		/* (I) NULL=errlog */










#endif	/*_drlptable_H*/



