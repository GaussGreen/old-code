/************************************************************************
 * Module:	DRL - ERR
 * File:	drlerr.h
 * Function:	
 * Revision:	$Header$
 ************************************************************************/
#ifndef _drlerr_H
#define	_drlerr_H

/*
 * Error log: default is ON, prints to stdout.
 */
extern	int	DrlErrOn(int boolean);
extern	void	DrlErrMsg(char *fmt, ...);
extern	int	DrlErrCallbackSet(int (*cbfunc)(const char *string));

/*
 * Trace: default level is 0 (off), prints to stdout.
 */
extern	long	DrlTraceLevelGet(void);
extern	void	DrlTraceLevelSet(long level);
extern	int	DrlTraceLevelCallbackSet(int (*cbset)(long), long (*cbget)());

/* Backward compatibility macro */
#define	DRL_IF_LOGGING(statement)	\
	if (DrlTraceLevelGet() > 0) { statement; }


#endif



