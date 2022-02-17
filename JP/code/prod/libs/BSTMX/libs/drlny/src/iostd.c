/************************************************************************
 * Module:	DRL - IO
 * Function:	I/O utilities
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlsys.h"		/* platform dependent system calls */
#include "drlstd.h"		/* platform compatibility */

#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "drlio.h"		/* prototype consistency */


#if defined(_WINDLL) && !(defined(WIN32) || defined(WIN32))
# define stderr	((FILE*)NULL)
#endif

/*f-------------------------------------------------------------
 * I/O: fprintf (admits NULL file pointer).
 *
 * <br><br>
 * Works just as a <i> fprintf</i>, but if the first argument <i> fp</i>
 * is NULL, the string is written to the error log.
 */

DLL_EXPORT(int)
DrlFPrintf(FILE *fp, char *fmt,...)
{
	va_list	ap;
	char	buf[2048], buf2[256];

	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);

	if (fp != (FILE*)NULL) {
		fputs(buf, fp);
	} else {
		char	*p = buf, *q;
		int	i, imax=255;
		while (*p != '\0') {
			q = buf2;
			i = 0;
			while ((i++ <= imax) && (*q++ = *p++));
			p--; q--;
			*q = '\0';
			GtoErrMsg("%s", buf2);
		}
	}
	return(0);
}



/*f-------------------------------------------------------------
 * I/O: fprintf directly to file.
 *
 * <br><br>
 * Works just as a <i> fprintf</i>, but the first argument <i> fnam</i>
 * is a filename where the formatted string is appened.
 * If <i> fnam</i> is a NULL or empty string, does nothing.
 * If <i> fnam</i> is the string <i> "stdout"</i>, writes to
 * the standard output.
 */

DLL_EXPORT(int)
DrlFilePrintf(char *fnam, char *fmt,...)
{
	FILE	*fp;
	va_list	ap;
	char	buf[256];

	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);

	if ((fnam == NULL) || (fnam[0] == '\0')) return(0);

	if (!strcmp(fnam, "stdout")) {
		fputs(buf, stdout);
	} else {
	    if ((fp = fopen(fnam, "a")) == NULL) {
		GtoErrMsg("FilePrintf: can't open `%s'\n", fnam);
		return(-1);
	    } else {
		fputs(buf, fp);
		fclose(fp);
	    }
	}
	return(0);
}



/*--------------------------------------------------------------
 * A dummy fscanf
 */

#if defined(CLIB) && (defined(_WINDLL) && !(defined(_WIN32) && defined(WIN32)))
DLL_EXPORT(int)
_dummyFscanf(FILE *fp, char *fmt,...)
{
	GtoErrMsg("fscanf not implemented\n");
	return(0);
}
#endif





#if defined(_WINDLL)
static	char	_stdoutFnam[] = "c:\\stdout.log";
static	char	_stderrFnam[] = "c:\\stderr.log";
static	FILE	*_stdoutFp = (FILE*) NULL,
		*_stderrFp = (FILE*) NULL;
#endif
#if defined(UNIX)
static	int	_stdoutSaveFd = -1,
		_stderrSaveFd = -1;
#endif



#if defined(_WINDLL)
/*f---------------------------------------------------------------------
 * Returns the file pointer to the standard output.
 * This function is available only in WINDLL where it replaces the
 * file pointer <i> stdout</i> (not available in windows) through
 * the macro
 * \begin{verbatim}
 * #if defined(_WINDLL)
 * # define stdout (DrlStdoutFp())
 * #endif
 * \end{verbatim}
 * The file opened is c:\stdout.log, unless it has been
 * changed by the function <i> StdoutFileSet</i> (see below).
 */

DLL_EXPORT(FILE*)
DrlStdoutFp(void)
{
	if (_stdoutFp == NULL) {
		DrlStdoutFileSet(_stdoutFnam, "w");
	}
	return(_stdoutFp);
}
#endif

/*f---------------------------------------------------------------------
 * Redirects the standard output to a file "pathname" to be opened
 * in "mode" (same as in <i> fopen</i>).
 * <br>
 * <br> In UNIX, if <i> pathname</i> is NULL, and 
 * the environment variable STDOUT\_LOG points to a valid pathname,
 * stdout is redirected to that file.
 * If it is not, it is redirected to <i> /dev/null</i>.\\
 * It makes a copy of the <i> stdout</i> file descriptor and returns it 
 * (the original stdout is retrieved by <i> DrlStdoutClose</i>.
 * <br> In Windows, opens a file <i> pathname</i>.
 * <br>
 * Returns -1 iff failed, a nonnegative integer otherwise.
 */

DLL_EXPORT(int)
DrlStdoutFileSet(char *pathname, char *mode)
{
static	char	routine[] ="DrlStdoutFileSet";
#if defined(_WINDLL)
	time_t	caltp;

	/* close if already opened */
	if (_stdoutFp != (FILE*)NULL) {
		fclose(_stdoutFp);
		_stdoutFp = (FILE*)NULL;
	}

	/* open new file */
	if ((_stdoutFp = fopen(pathname, "w")) == (FILE*)NULL) {
		GtoErrMsg("%s: can't open `%s' (%s)\n", routine,
			pathname, DrlStrerror());
		abort();
	}

	/* print current time */
	if (time(&caltp) != -1) {
        	fprintf(stdout, "\t%s", ctime(&caltp));
	}


	return(SUCCESS);

#elif defined(UNIX)
	int	newfd, fd, oflag;
	mode_t	tmode;
	time_t	caltp;
static	char	devNullFnam[] = "/dev/null";

	/* make a copy of stdout */
	if ((newfd = dup(STDOUT_FILENO)) < 0) {
		GtoErrMsg("%s: duplicate failed (%s)\n",
			routine, DrlStrerror());
		return(-1);
	}
	_stdoutSaveFd = newfd;


	/* open a file */
	tmode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
	if (mode[0] == 'w') {
		oflag = O_WRONLY | O_CREAT | O_TRUNC;
	} else if (mode[0] == 'a') {
		oflag = O_WRONLY | O_CREAT | O_APPEND;
	} else {
		GtoErrMsg("%s: bad mode `%s'\n", routine, mode);
		return(-2);
	}

	/* if env variable STDOUT_LOG is set */
	if (pathname == NULL) {
	    if ((pathname = getenv("STDOUT_LOG")) == NULL) {
		pathname = devNullFnam;
	    }
	}


	if ((fd = open(pathname, oflag, tmode)) <= 0) {
		GtoErrMsg("%s: open `%s' failed (%s)\n",
			routine, pathname, DrlStrerror());
		return(fd);
	}

	/* unbuffered I/O */
	setvbuf(stdout, NULL, 0, _IONBF);

	/* copy */
	if (dup2(fd, STDOUT_FILENO) < 0) {
		GtoErrMsg("%s: duplicate failed (%s)\n",
			routine, DrlStrerror());
		return(fd);
	}
	if (close(fd) != 0) {
		GtoErrMsg("%s: close failed (%s)\n",
			routine, DrlStrerror());
		return(-1);
	}


	/*if (fcntl(newfd, F_SETFL, O_SYNC) < 0) {
	    GtoErrMsg("%s: fcntl failed (%s)\n", routine,
				DrlStrerror());
	    return(-1);
        }*/


	/* print current time */
	if (time(&caltp) != -1) {
        	fprintf(stdout, "\t%s", ctime(&caltp));
	}


	return(newfd);
#else
	GtoErrMsg("%s: not implemented.\n", routine);
	return(-1);
#endif
}



/*f---------------------------------------------------------------------
 * Closes the standard output under WINDLL.
 * Under UNIX, closes the duplicate file descriptor of stdout (if any)
 * and restores the original stdout.
 * {\bf This function  must be called before returning from an
 * addin call.}
 */

DLL_EXPORT(int)
DrlStdoutClose(void)
{
static	char	routine[] = "DrlStdoutClose";
#if defined(_WINDLL)
	/* close if already opened */
	if (_stdoutFp != (FILE*)NULL) {
		fclose(_stdoutFp);
		_stdoutFp = (FILE*)NULL;
	}
#elif defined(UNIX)
	fflush(stdout);
	if (_stdoutSaveFd >= 0) {
	    if (dup2(_stdoutSaveFd, STDOUT_FILENO) < 0) {
		GtoErrMsg("%s: duplicate failed (%s)\n",
			routine, DrlStrerror());
		return(-1);
	    }
	    if (close(_stdoutSaveFd) != 0) {
		GtoErrMsg("%s: close failed (%s)\n",
			routine, DrlStrerror());
		return(-1);
	    }
	    _stdoutSaveFd = -1;
	}
#endif
	return(0);
}


/*f---------------------------------------------------------------------
 * Similar to <i> DrlStdoutFileSet</i>, but for the standard error
 * (in UNIX, the environment variable STDERR\_LOG must be set to redirect).
 */

DLL_EXPORT(int)
DrlStderrFileSet(char *pathname, char *mode)
{
static	char	routine[] ="DrlStderrFileSet";
#if defined(_WINDLL)
	time_t	caltp;

	/* close if already opened */
	if (_stderrFp != (FILE*)NULL) {
		fclose(_stderrFp);
		_stderrFp = (FILE*)NULL;
	}

	/* open new file */
	if ((_stderrFp = fopen(pathname, "w")) == (FILE*)NULL) {
		GtoErrMsg("%s: can't open `%s' (%s)\n", routine,
			pathname, strerror(errno));
		abort();
	}

	/* print current time */
	if (time(&caltp) != -1) {
        	fprintf(stderr, "\t%s", ctime(&caltp));
	}


	return(SUCCESS);

#elif defined(UNIX)
	int	newfd, fd, oflag;
	mode_t	tmode;
	time_t	caltp;
static	char	devNullFnam[] = "/dev/null";

	/* make a copy of stderr */
	if ((newfd = dup(STDERR_FILENO)) < 0) {
		GtoErrMsg("%s: duplicate failed (%s)\n",
			routine, strerror(errno));
		return(-1);
	}
	_stderrSaveFd = newfd;


	/* open a file */
	tmode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
	if (mode[0] == 'w') {
		oflag = O_WRONLY | O_CREAT | O_TRUNC;
	} else if (mode[0] == 'a') {
		oflag = O_WRONLY | O_CREAT | O_APPEND;
	} else {
		GtoErrMsg("%s: bad mode `%s'\n", routine, mode);
		return(-2);
	}

	/* if env variable STDERR_LOG is set */
	if (pathname == NULL) {
	    if ((pathname = getenv("STDERR_LOG")) == NULL) {
		pathname = devNullFnam;
	    }
	}


	if ((fd = open(pathname, oflag, tmode)) <= 0) {
		GtoErrMsg("%s: open `%s' failed (%s)\n",
			routine, pathname, strerror(errno));
		return(fd);
	}

	/* unbuffered I/O */
	setvbuf(stderr, NULL, 0, _IONBF);

	/* copy */
	if (dup2(fd, STDERR_FILENO) < 0) {
		GtoErrMsg("%s: duplicate failed (%s)\n",
			routine, strerror(errno));
		return(fd);
	}
	if (close(fd) != 0) {
		GtoErrMsg("%s: close failed (%s)\n",
			routine, strerror(errno));
		return(-1);
	}


	/*if (fcntl(newfd, F_SETFL, O_SYNC) < 0) {
	    GtoErrMsg("%s: fcntl failed (%s)\n", routine, strerror(errno));
	    return(-1);
        }*/


	/* print current time */
	if (time(&caltp) != -1) {
        	fprintf(stderr, "\t%s", ctime(&caltp));
	}


	return(newfd);
#else
	GtoErrMsg("%s: not implemented.\n", routine);
	return(-1);
#endif
}



/*f---------------------------------------------------------------------
 * Similar to <i> DrlStdoutClose</i> but for stderr.
 */

DLL_EXPORT(int)
DrlStderrClose(void)
{
static	char	routine[] = "DrlStderrClose";
#if defined(_WINDLL)
	/* close if already opened */
	if (_stderrFp != (FILE*)NULL) {
		fclose(_stderrFp);
		_stderrFp = (FILE*)NULL;
	}
#elif defined(UNIX)
	if (_stderrSaveFd >= 0) {
	    if (dup2(_stderrSaveFd, STDERR_FILENO) < 0) {
		GtoErrMsg("%s: duplicate failed (%s)\n",
			routine, strerror(errno));
		return(-1);
	    }
	    if (close(_stderrSaveFd) != 0) {
		GtoErrMsg("%s: close failed (%s)\n",
			routine, strerror(errno));
		return(-1);
	    }
	    _stderrSaveFd = -1;
	}
#endif
	return(0);
}


#if !defined(CLIB)
static	char	*errLogFnam[] = "error.log";
static	int	errLogSet = 1;

/*--------------------------------------------------------------
 * Only if CLIb not available
 */

DLL_EXPORT(int)
DrlPrintErrMsg(char *fmt,...)
{
	FILE	*fp;
	va_list	ap;
	char	buf[256];
static	char	fnam[] = "error.log";

	if (errLogSet != FALSE) 
		return(SUCCESS);

	if ((fp = fopen(fnam, "a")) != NULL) {
		va_start(ap, fmt);
		vsprintf(buf, fmt, ap);
		va_end(ap);
		fputs(buf, fp);
		fclose(fp);
		return(SUCCESS);
	} else {
		return(FAILURE);
	}
}

/*--------------------------------------------------------------
 * Only if CLIb not available
 */

DLL_EXPORT(int)
PrintErrMsgSet(int errMsgStatus)
{
	if (errLogSet != 0) {
		errLogSet = 0;
		return(1);
	} else {
		errLogSet = 1;
		return(0);
	}
}

#endif


