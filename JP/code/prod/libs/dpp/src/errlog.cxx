/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher, D. Liu
 * Revision:	$Header$
 ************************************************************************/
#if defined(_WIN32)
# include <cstdio>
#endif
#include <errno.h>
#include "kstdinc.h"
#include "kutilios.h"

#if defined(WIN_NT)
# include <io.h>
# include <fcntl.h>
# include <sys\stat.h>
# include <sys\types.h>
#endif


extern	"C" {
#include "cgeneral.h"
#include "cerror.h"
#include "bastypes.h"
#include "macros.h"             /* MAX */
#include "ldate.h"              /* GtoDayCountFraction */
#include "ratelist.h"		/* TRateList */
#include "tcurve.h"
#include "date_sup.h"

#include "drlio.h"		/* FScanStruct */
#include "drlstr.h"		/* StringLineVarScan */
#include "drltime.h"		/* TDatePrint */
#include "drloptio.h"		/* Black */
};


// Callback used for the capture of ALIB error messages.

extern "C" {
	PRIVATE TBoolean errMsgCallBackFunc(
				char *message, void *callBackData);
};



//---------------------------------------------------------------
// Logger class:
// Contains two separate logging and error message facilities.
//
//---------------------------------------------------------------

class KDppLogger {
public:
	KDppLogger();
	~KDppLogger();

//private:
			// log messages 
	ostream	*mLogStream;
	FILE	*mLogFILE;
	bool	mLogSet;
	int	mLogLevelOffset;
			// Error messages
	ostream	*mErrStream;
	FILE	*mErrFILE;
	int	mErrSet;
};

//---------------------------------------------------------------
// Global logger object: automatically intialized and cleaned.

	KDppLogger	dppLogger;


//---------------------------------------------------------------
// By default, redirects log to stdout and err to stderr.
//

KDppLogger::KDppLogger()
{
	// Initialize to stdout by default
	setvbuf(stdout, NULL, _IONBF, 0);

	mLogFILE = stdout;
	mLogStream = &cout;
	mLogSet = false;
	mLogLevelOffset = 0;

	mErrFILE = stdout;
	mErrStream = &cout;
	mErrSet = DPP_ERR_MSG_OFF;

	cout.sync_with_stdio();
	cerr.sync_with_stdio();

	// Register ALIB callback (only once)
	GtoErrMsgAddCallback(
		(TErrCallBackFunc*)errMsgCallBackFunc, 
		FALSE,
		NULL); 

	// Turn ALIB error message on
	GtoErrMsgOn();

}



//---------------------------------------------------------------

KDppLogger::~KDppLogger()
{
	if (mLogSet)  {
		delete mLogStream;
		mLogStream = NULL;
		fclose(mLogFILE);
		mLogFILE = NULL;
	}

	if (mErrSet == DPP_ERR_MSG_BUFFER) {
		delete mErrStream;
		mErrStream = NULL;
		fclose(mErrFILE);
		mErrFILE = NULL;
	}
}


//---------------------------------------------------------------

int	DppLogLevel()
{
	return GtoLoggingGet() + dppLogger.mLogLevelOffset;
}




//---------------------------------------------------------------


void
DppLogMsgSetFile(const char* fnam)
{
static	char	routine[] = "DppLogMsgSetFile";

    try {

	// Backwards compatibility ...
	if ((fnam == NULL) || (!strcmp(fnam, "stdout"))) {
		return;
	}

	// Close if opened
	if (dppLogger.mLogSet)  {
		delete dppLogger.mLogStream;
		dppLogger.mLogStream = NULL;
		fclose(dppLogger.mLogFILE);
		dppLogger.mLogFILE = NULL;
	}


	// Open fd
/*
	fd = _open(fnam, _O_WRONLY | _O_CREAT | _O_TRUNC,
		_S_IREAD | _S_IWRITE);
	if (fd == -1) {
		cerr << strerror(errno) << endl;
		cerr << routine << ": _open failed." << endl;
		throw KFailure();
	}
*/
	//
	dppLogger.mLogSet = true;

	// Open a FILE for logging
	dppLogger.mLogFILE = fopen(fnam, "w");
	if (dppLogger.mLogFILE == NULL) {
		cerr << strerror(errno) << endl;
		cerr << routine << ": new FILE failed." << endl;
		throw KFailure();
	}
	setvbuf(dppLogger.mLogFILE, NULL, _IONBF, 0);


	// Open a stream for logging
/*
#if defined(UNIX) && defined(SOLARIS)
	dppLogger.mLogStream = new ofstream((int) fileno(dppLogger.mLogFILE));
#else
*/
	dppLogger.mLogStream = new ofstream(fnam, ios_base::out | ios_base::app);

	if (dppLogger.mLogStream == NULL) {
		cerr << strerror(errno) << endl;
		cerr << routine << ": new fstream failed." << endl;
		throw KFailure();
	}


    }
    catch (...) {
	cerr << routine << ": failed on `" <<
		(fnam ? fnam : "NULL") << "'." << endl;
	throw KFailure();
    }
}


//---------------------------------------------------------------


void
DppLogMsgCloseFile()
{
	// Close if opened
	if (dppLogger.mLogSet)  {
		delete dppLogger.mLogStream;
		dppLogger.mLogStream = NULL;
		fclose(dppLogger.mLogFILE);
		dppLogger.mLogFILE = NULL;
	}

	dppLogger.mLogSet = false;
}


//---------------------------------------------------------------
//	Error logging status
//

int
DppErrMsgStatus()
{
	return dppLogger.mErrSet;

}


//---------------------------------------------------------------
//
void
DppErrMsgSet(int msgType)
{
	dppLogger.mErrSet = msgType;

}



//---------------------------------------------------------------
//
void
DppSetToUseGtoErrMsg(const char* fnam, TBoolean append)
{
static	char	routine[] = "DppSetToUseGtoErrMsg";

    try {

	char fn[256];	

	strcpy(fn, fnam);

	// Turn on error log
	//
	dppLogger.mErrSet = DPP_ERR_MSG_GTO;

	IF_FAILED_THROW(GtoErrMsgFileName(fn, append));

    }
    catch (...) {
	cerr << routine << ": failed on `" <<
		(fnam ? fnam : "NULL") << "'." << endl;
    }
}






//---------------------------------------------------------------


ostream& DppLogOs()
{
	return *dppLogger.mLogStream;
}


//---------------------------------------------------------------


ostream& DppErrOs()
{
	return *dppLogger.mErrStream;
}


//---------------------------------------------------------------

FILE*	DppLogFp()
{
	return dppLogger.mLogFILE;
}


//---------------------------------------------------------------
// Error logging with C linkage (to be called from C).

extern "C" {
void
DppErrMsg(char *fmt, ...)
{
	va_list	ap;
static	char	tmpBuf[1024];


	va_start(ap, fmt);
	vsprintf(tmpBuf, fmt, ap);
	va_end(ap);

	dppErr << tmpBuf;

	return;
}
};



//===============================================================
//
//	Error Retrieving From Buffer
//
//===============================================================


//---------------------------------------------------------------
// Initializes the redirection of error messages to an internal buffer
// and also captures all ALIB error messages.
//

void
DppErrMsgBufferInit (void)
{
static  char    routine[] = "DppErrMsgBufferInit";


	// Free if allocated
	if (dppLogger.mErrSet == DPP_ERR_MSG_BUFFER)
	{
		delete dppLogger.mErrStream; 
		dppLogger.mErrStream = NULL;

	}

	// Create a string stream
	dppLogger.mErrSet = DPP_ERR_MSG_BUFFER;
	dppLogger.mErrStream = new ostringstream();
	if (dppLogger.mErrStream == NULL) {
		cerr << strerror(errno) << endl;
		cerr << routine << ": new string stream failed." << endl;
		throw KFailure();
	}

	return;
}


//---------------------------------------------------------------

long
DppErrMsgBufferSize(void)
{
	// Add one char for '\0' 
	long strSize  =
	   	(dynamic_cast<ostringstream*>(dppLogger.mErrStream))->str().size() + 1;

	return strSize;
}


//---------------------------------------------------------------


int
DppErrMsgBufferRetrieve(
	long	outBufLen,	// (I) length of the input buffer.
	char	*outBuffer)	// (I/O) output buffer.
{
static  char    routine[] = "DppErrMsgBufferRetrieve";
	static int status = FAILURE;
 
	long   		globStrlen;

	ostringstream	*currErrStream = NULL;
 
  try{

	if (dppLogger.mErrStream == NULL)
	{
		cerr << "Error buffer is null." << endl;
		return status;

	}

	//
	// Cast to ostringstream, so can work with string class.
	//
	currErrStream = dynamic_cast<ostringstream*>(dppLogger.mErrStream);

	// String length
	globStrlen = currErrStream->str().size();

	if (outBuffer IS NULL ||
	    outBufLen < globStrlen)
	{
		cerr << routine << ": Either output buffer size ("
		     << outBufLen << ") is smaller than required ("
		     << globStrlen << "), or output buffer is NULL." 
		     << endl;
		return (status);
	}
 
	// Copy the result, including '\0'.
	strncpy (outBuffer, currErrStream->str().c_str(), globStrlen+1);

	status = SUCCESS;

	return (status);
 
    }
    catch (...) {
        cerr << format("%s: failed (uncaught).\n", routine);
        return (status);
   }
}




//---------------------------------------------------------------
// Callback used for the capture of ALIB error messages.

extern "C" {
PRIVATE TBoolean
errMsgCallBackFunc(
	char *message, 		// (I) input error message 
	void *callBackData) 	// (I/O) call back data 
{ 
	// Print to error stream
	(*dppLogger.mErrStream) << message;

	// Exit with normal status
	return (TRUE); 
}

}; // extern "C"



