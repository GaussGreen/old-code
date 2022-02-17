/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	A platform independent getopt
 * File:	drlgetop.h
 * Function:	header
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlgetop_H
#define	_drlgetop_H

#if !defined(UNIX) || defined(INCLUDE_DRL_GETOPT)
#if defined(UNIX)
#include <unistd.h>
#endif

/*
 * On non-unix platforms, getopt() is not available.
 */
extern	char	*drl_optarg;
extern	int	drl_optind;
extern	int	drl_opterr;
extern	int	drl_getopt(int argc, char *const *argv,
				const char *shortopts);
/* Redefines */
#define	optarg	drl_optarg
#define	optind	drl_optind
#define	opterr	drl_opterr
#define	getopt	drl_getopt

#else
#include <unistd.h>
#endif



#endif	/* _drlgetop_H */
