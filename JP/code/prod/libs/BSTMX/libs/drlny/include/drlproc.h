/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	PROC - Process Control
 * File:	drlproc.h
 * Function:	Processes and system functions.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlproc_H
#define	_drlproc_H

#include "drlstd.h"

/* System Info */
extern	DLL_EXPORT(char*)	DrlSystemInfo(char *s);
extern	DLL_EXPORT(const char*)	DrlStrerror();

/* Eleapsed time */
extern	DLL_EXPORT(int)		DrlChrono(char *what, int chronoNo, ...);


#endif	/* _drl_H */
