/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	FPARS - Formula Parser
 * File:	drlgpars.h
 * Function:	header
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlgpars_H
#define	_drlgpars_H

#include "drlstd.h"

/*
 *  Formula parser
 */
extern	DLL_EXPORT(int)	DrlGParEval(
	int n,
	double *x,
	char *form,
	double *retVal);


/*
 * The following is included in the parser source code
 * only.
 */
#if defined(_DRLGPAR_SRC)
extern	int		_gParEval(int n, double *x, char *form,
			double *retVal);



#define GPARNSYMS	32	/* maximum number of symbols */
#define GPARMAXS	32	/* maximum number of char in symbol */

typedef	struct {
	char	name[GPARMAXS];
	int	setFlag;
	double	(*funcptr)();
	double	value;
} TgparSymbol;


typedef	struct {
	int		numSym;
	TgparSymbol	table[GPARNSYMS];
} TgparSymbolTable;


extern	int		_gparSymbolTableClear(void);
extern	TgparSymbol*	_gparSymbolTableLookup(char *s);

#endif	/*UNIX*/


#endif	/* _drlgpars_H */
