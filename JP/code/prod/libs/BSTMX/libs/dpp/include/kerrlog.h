/***************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	kerrlog.h
 * Function:	Error amd logging messages management.
 * Author:	Christian Daher
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kerrlog_H
#define	_kerrlog_H

#include "kplatdep.h"
#include "kexcept.h"

//==============================================================
//			OVERVIEW 
// 
// Enables and configure error and logging for Dpp modules.
// 
// There are two logging channels for the Dpp library:
// a standard logging to report progress of the program and
// an error logging to report errors. 
// 
// The standard logging is accessed through the routines
// DppLogOs (for an ostream) or DppLogFp (for a FILE pointer).
// There is a global variable that defines the logging level
// whose value is set by "logLevel" and
// can be retrieved by the routine DppLogLevel.
// 
//
//==============================================================


//==============================================================
//
// 	Error Messages Functionalty
//
//==============================================================


/**
 * Returns the current error ostream (default standard output).
 */
ostream&	DppErrOs();


/**
 * Useful macro.
 */

# undef  dppErr
# define dppErr		(DppErrOs())


extern "C" {
/**
 * Prints a message to the error message stream.
 * Defined with C linkage so it can be called from C.
 */
extern	void	DppErrMsg(char *fmt, ...);
};



//==============================================================
//
// 	Logging functionalty
//
//==============================================================


/**
 * Returns the standard logging level.
 */

int	DppLogLevel();

//
// Useful macro for logging
//
#define PG_IF_LOGGING(loggingLevel, printStatement)   \
	{if (GtoLoggingGet() >= loggingLevel) \
	{ printStatement; }}


/**
 * Returns the logging ostream (default standard output).
 */
ostream&	DppLogOs();


/**
 * Returns the logging FILE pointer (default standard output).
 */
FILE*		DppLogFp();


/**
 * Sets the logging to be redirected to a file.
 */
void		DppLogMsgSetFile(const char* fnam);


/**
 * Close the logging file if opened.
 */
void 		DppLogMsgCloseFile();




#define	DPP_ERR_MSG_OFF		0
#define	DPP_ERR_MSG_COUT	1
#define	DPP_ERR_MSG_GTO		2
#define	DPP_ERR_MSG_BUFFER	3

/**
 * Set the error logging method. 
 * The error logging uses the following schemes:
 * 1. stdout or cout as default
 * 2. Use GtoErrMsg() to write to a specified file, which is set by 
 *    DppSetToUseGtoErrMsg.
 * 3. Writing to a string buffer, initiated by calling DppErrMsgBufferInit.
 */
void		DppErrMsgSet(int msgType);


/**
 * Return the error logging status.
 */
int 		DppErrMsgStatus();


/**
 * Sets the error logging to use GtoErrMsg.
 */
void		DppSetToUseGtoErrMsg(const char* fnam, TBoolean append);


/**
 * Initializes the redirection of error messages to an internal buffer
 * and also captures all ALIB error messages.
 * The error messages are retrieved by the functions
 * DppErrMsgBufferSize and DppErrMsgBufferRetrieve
 */

extern	void	DppErrMsgBufferInit(void);


/**
 * Retrieves error messages that are written to the  global buffer.
 * Returns FAILURE if outBufLen is less than the length
 * of the accumulated error message string
 */ 
extern	int	DppErrMsgBufferRetrieve(
	long	outBufLen,	// (I) length of the input buffer.
	char	*outBuffer);	// (I/O) output buffer.

/**
 * Return the size of error buffer.
 */
extern	long	DppErrMsgBufferSize(void);




# undef  dppLog
# define dppLog		(DppLogOs())
# undef  dppFpLog
# define dppFpLog	(DppLogFp())



//==============================================================
//
// Logging level predefined constant
//
//==============================================================


#define	PG_LOG_CALCSENS_D	7
#define	PG_LOG_CALCSENS		7
#define	PG_LOG_PRIPF		10
#define	PG_LOG_PRI		11
#define	PG_LOG_PRI_D		18
#define	PG_LOG_VRP		2
#define	PG_LOG_VRP_D		3

#define	PG_LOG_MKTD		4
#define	PG_LOG_SML		13

#define K_NAME_LEN      	32              // max char in name


#endif


