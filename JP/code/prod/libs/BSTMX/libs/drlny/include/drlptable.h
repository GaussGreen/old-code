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
#include "drlstd.h"	/* wrapper type def (TVType, FloatL, etc.) */

#include "drlsmat.h"	/* TSwaptionMatrix2D */

/*t-@CDOC(idxn=TVTable,catn=structdef)-----------------------------------
 * Pivot table data structure for elements of arbitrary type.
 */

typedef	struct	{
	int	numFields;		/* number of fields (columns) */
	TVType	*fieldType;		/* field type [numFields] */
	char	**fieldName;		/* field name [numFields] */
	int	*numFieldValues;	/* nb different val [numFields] */
	int	*numFieldValuesMax;	/* max size (realloc) [numFields] */
	void	**fieldValues;		/* field value [numFields][numValues] */

	int	numData;		/* current number of elements (rows) */
	int	numDataMax;		/* max size (for realloc) */
	void	**dataValue;		/* data [numFields][numData] */

	void	**tmpData;		/* working space value [numFields] */

} TVTable;
/*e*/

extern	TVTable*	DrlTVTableNew(int numFields, TVType *fieldType);
extern	TVTable*	DrlTVTableNewFpRead(FILE *fp, int numFields,
				TVType *fieldType);
extern	TVTable*	DrlTVTableNewFileRead(char *fnam, int numFields,
				TVType *fieldTypes);
extern	int		DrlTVTableFree(TVTable* that);
extern	int		DrlTVTableFpWrite(TVTable *that, FILE *fp);

extern	int		DrlTVTableFindFieldName(TVTable *that,
				char *nameField);

extern	int		DrlTVTableExtractMatrix(TVTable *that,
				int indexFieldX, int indexFieldY,
				int indexFieldP, void *valFieldP,
				int indexFieldXY,
				int *numItemsX, int *numItemsY,
				TDate baseDate,
				double ***matrix);
extern	int		DrlTVTableExtractMatrixN(TVTable *that,	
				char *nameFieldX, char *nameFieldY,
				char *nameFieldP, void *valFieldP,
				char *nameFieldXY,
				int *numItemsX, int *numItemsY,
				TDate baseDate,
				double ***matrix);
extern	int		DrlTVTableExtractTSwaptionMatrix2D(TVTable *that,
				char *nameFieldX, char *nameFieldY,
				char *nameFieldP, void *valFieldP,
				char *nameFieldXY,
				TDate baseDate,
        			long freq,
				TSwaptionMatrix2D **mat);

/*t-@CDOC(idxn=TCTable,catn=structdef)-----------------------------------
 * Structure for table of string elements.
 */

typedef	struct	{
	int	numRows;
	int	numCols;
	int	numRowsMax;
	int	numColsMax;
	char	***cells;		/* [row][vol] */
} TCSheet;


#define	TCSheet_END	INT_MAX
#define	TCSheet_START	0



extern	TCSheet* DrlTCSheetNewEmpty(
			int numRowsMax,		/* (I) max number of rows */
			int numColsMax);	/* (I) max number of cols */
extern	int	DrlTCSheetFree(TCSheet* that);

extern	int	DrlTCSheetSet(TCSheet* that,
			int row, int col, const char* value);
extern	char*	DrlTCSheetGet(TCSheet* that, int row, int col);
extern	int	DrlTCSheetGetVType(
			TCSheet* that,		/* (I) table */
			int row,		/* (I) row number */
			int col, 		/* (I) column number */
			TVType fieldType,	/* (I) element type */
			void *ptr);		/* (O) */

extern	int	DrlTCSheetFpRead(
			TCSheet **that,		/* (O) output sheet */
			FILE *fp);		/* (I) file pointer */
extern	int	DrlTCSheetFileRead(
			TCSheet **that,		/* (O) output sheet */
			const char *fnam);	/* (I) file name */

extern	int	DrlTCSheetFpWrite(
			TCSheet *that,		/* (I) table */
			FILE *fp);		/* (I) NULL=errlog */










#endif	/*_drlptable_H*/



