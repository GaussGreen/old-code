/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	vtsparh.h
 * Function:	
 * Author:	C. Daher
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_vtsparh_H
#define	_vtsparh_H
#include "drlstd.h"	/* Definitions and error handling */


#ifdef __cplusplus
extern "C" {
#endif


/*
 * The following is included in the parser source code
 * only.
 */


#define INSPARMAXS	32	/* maximum number of char in symbol */
#define	MAXSTRLEN	256	/* maximum string length  */

/* init/clean the lexer */
extern  int	vtsparLexInit(const char *s);
extern  int	vtsparLexClean();

/*
 * Definition of a symbol.
 */

typedef	struct {
	char		name[INSPARMAXS];	/* symbol */
	int		setFlag;		/* T=symbol is set */
	long		type;			/* symbol type */

	int		ival;
	double		dval;			/* double value */
	long		lval;			/* long value */
	TDate 		dtval;			/* date value */
	TDateInterval	itval;			/* date interval value */
	TDayCount	dcval;			/* day count value */

	void		*ptr;			/* symbol value */
} VTSSym;


/* Symbol table */

typedef	struct {
	int		numSym;			/* current # of symbols */
	int		numSymMax;		/* maximum # of symbols */
	VTSSym		*table;			/* array of symbols */
} VTSSymTable;


int	vtsParCopyReturn(VTSSym *symb);
int	vtsParCopyReturnConstant(double value);
int	vtsParGetArg(VTSSym *symbO, int idx);
int	vtsParLetDouble(VTSSym *symbO, double value);
int	vtsParUnaryOper(VTSSym *symbO, VTSSym *symb1, char *oper);
int	vtsParScalarOper(VTSSym *symbO, VTSSym *symb1, double value, char *oper);
int	vtsParBinaryOper(VTSSym *symbO, VTSSym *symb1, VTSSym *symb2, char *oper);
int	vtsParTernaryOper(VTSSym *symbO, VTSSym *symb1, VTSSym *symb2,
			VTSSym *symb3, char *oper);






#ifdef __cplusplus
};
#endif


#endif	/*_vtsparh_H*/

