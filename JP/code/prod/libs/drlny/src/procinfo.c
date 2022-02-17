/************************************************************************
 * Module:	DRL
 * Submodule:	PROC
 * File:	
 * Function:	Processes Management
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlsys.h"		/* system calls */
#include "drlstd.h"		/* platform compatibility */

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#ifndef	_WIN32
# include <errno.h>
#endif

#include "drlmem.h"
/*#include "drlio.h"*/

#include "drlproc.h"		/* Prototype consistency */


#if defined(UNIX)
static	char*		ElpsTimePrint(char *s, clock_t start,
				struct tms *tmsstart,
				clock_t end,
				struct tms *tmsend);
# define NMAX_CHRONO	12
# if defined(_USETMS)
static	struct tms	chronoTmsstart[NMAX_CHRONO],
			chronoTmsend[NMAX_CHRONO];
# endif /*_USETMS*/
static	double		chronoElps[NMAX_CHRONO];
static	clock_t		chronoStart[NMAX_CHRONO],
			chronoEnd[NMAX_CHRONO];
#elif !defined(_WINDLL)
# define NMAX_CHRONO	2
static	double		chronoElps[NMAX_CHRONO];
static	clock_t		chronoStart[NMAX_CHRONO],
			chronoEnd[NMAX_CHRONO];
#endif



/*f-----------------------------------------------------
 * System : compute CPU time.
 * 
 * <br><br>
 * Routine to compute CPU time used by a process.\\
 * <b> Example:</b>
 * <pre>
 *        int chronoNo = 0;
 *        DrlChrono("START", chronoNo);
 *        ...
 *        DrlChrono("STOP", chronoNo);
 *        ...
 *        DrlChrono("RESUME", chronoNo);
 *        ...
 *        DrlChrono("STOP", chronoNo);
 *        ...
 *        printf("Real time: %d\n", DrlChrono("REAL", chronoNo));
 *        ...
 *        DrlChrono("SPRINT", chronoNo, (char*) buf);
 *        printf("%s\n", buf);
 *        ...
 *        DrlChrono("FPRINT", chronoNo, (FILE*) fp);
 *        ...
 * </pre>
 */

DLL_EXPORT(int)
DrlChrono(char *what, int chronoNo, ...)
{
static	char	routine[] = "DrlChrono";
#if !defined(_WINDLL)
	int	i;
	char	uWhat[32],
		*s;
	FILE	*fp;
	va_list	ap;

#ifdef	_SKIP
static	long	clktck  = 0 ;

	/* get the clock tick */
	if (clktck == 0) {
		if ((clktck = sysconf(_SC_CLK_TCK)) < 0) {
			GtoErrMsg("%s: clock tick N/A (sysconf error)\n",
				routine);
			return(-1);
		}
		/*printf("clktck = %ld\n", clktck);*/
	}
#endif


	/* convert argument to uppercase */
	for (i = 0; i<= strlen(what); i++)
		uWhat[i] = (char) toupper((int) what[i]);

	/* check chrono # is OK */
	if ((chronoNo < 0) || (chronoNo >= NMAX_CHRONO)) return(2);

	/* check action to perfom */
	if (strcmp(uWhat, "START") == 0) {
	    /*
	     * START
	     */
#if defined(_USETMS)
	    if ((chronoStart[chronoNo] =
			times(&chronoTmsstart[chronoNo])) == -1)
			return(-1);
#endif	/*_USETMS*/
	    if ((chronoStart[chronoNo] = clock()) == -1)
			return(-2);
	    chronoElps[chronoNo] = 0.;

	} else if (strcmp(uWhat, "STOP") == 0) {
	    /*
	     * STOP
	     */
#if defined(_USETMS)
	    if ((chronoEnd[chronoNo] =
			times(&chronoTmsend[chronoNo])) == -1)
			return(-1);
#endif /*_USETMS*/
	    if ((chronoEnd[chronoNo] = clock()) == -1)
			return(-2);

	    chronoElps[chronoNo] += ((double) chronoEnd[chronoNo]
				- (double) chronoStart[chronoNo])
				/ (double) CLOCKS_PER_SEC;

	} else if (strcmp(uWhat, "RESUME") == 0) {
	    /*
	     * RESUME
	     */
#if defined(_USETMS)
	    if ((chronoStart[chronoNo] =
			times(&chronoTmsstart[chronoNo])) == -1)
			return(-1);
#endif /*_USETMS*/
	    if ((chronoStart[chronoNo] = clock()) == -1)
			return(-2);

	} else if (strcmp(uWhat, "REAL") == 0) {
	    /*
	     * REAL
	     */
	    return((int) chronoElps[chronoNo]);
	} else if (strcmp(uWhat, "SPRINT") == 0) {
	    /*
	     * SPRINT
	     */
	    /* Print real time in a string */
	    va_start(ap, chronoNo);
	    s = (char*) va_arg(ap, char*);
	    va_end(ap);
	    sprintf(s, "Real: %7.2f", chronoElps[chronoNo]);
	} else if (strcmp(uWhat, "FPRINT") == 0) {
	    /*
	     * FPRINT
	     */
	    /* Print real time in a file pointer */
	    va_start(ap, chronoNo);
	    fp = (FILE*) va_arg(ap, FILE*);
	    va_end(ap);
	    if (fp != NULL) {
		fprintf(fp, "Real: %7.2f\n", chronoElps[chronoNo]);
	    } else {
		GtoErrMsg("Real: %7.2f\n", chronoElps[chronoNo]);
	    }
	} else {
	    GtoErrMsg("%s: bad flag `%s'\n", routine, uWhat);
	    return(-1);
	}
	return(0);

#else /*defined(_WINDLL)*/
	GtoErrMsg("%s: routine not available\n", routine);
	return(-1);
#endif
}


/*------------------------------------------------------
 * clock_t	start, end;
 * struct tms	tmsstart, tmsend;
 * ....  assert((start = times(&tmsstart)) != -1);
 * ....  assert((end = times(&tmsend)) != -1);
 * ....  ElpsTime(...)
 */


#if defined(_USETMS)
static	char*
ElpsTimePrint(char *s,
	clock_t start, struct tms *tmsstart,
	clock_t end, struct tms *tmsend)
{
static	char	buf[255], *tmp ;
static	long	clktck  = 0 ;

	tmp = (s != NULL ? s : buf);

	if (clktck == 0)
		if ((clktck = sysconf(_SC_CLK_TCK)) < 0) {
			strcat(tmp, "N/A (sysconf error)");
			return(tmp);
		}
	/*printf("clktck = %ld\n", clktck);*/

	sprintf(tmp,
	    "Real: %7.2f User: %7.2f  Sys: %7.2f  "
	    "CUser: %7.2f CSys: %7.2f",
	    ((double) end - (double) start) / ((double) CLOCKS_PER_SEC),
	    (tmsend->tms_utime - tmsstart->tms_utime) / ((double) clktck),
	    (tmsend->tms_stime - tmsstart->tms_stime) / ((double) clktck),
	    (tmsend->tms_cutime - tmsstart->tms_cutime) / ((double) clktck),
	    (tmsend->tms_cstime - tmsstart->tms_cstime) / ((double) clktck));

	return(tmp);

}
#endif	/*_USETMS*/





/*f-------------------------------------------------------------
 * System : information printout.
 *
 * <br><br>
 * Prints in the string <i> s</i> a one line summary of the system
 * (hostname, operating system, etc).
 */

DLL_EXPORT(char*)
DrlSystemInfo(char *s)
{
static	char	buf[255];
#if defined(UNIX)
	struct utsname	theName;

	s = (s != NULL ? s : buf);

	if (uname(&theName) == -1) {
		strcat(s, "N/A");
	} else {
		sprintf(s, "%s %s R%s V%s",
			theName.sysname,
			theName.nodename,
			theName.release,
			theName.version);
	}
	return(s);
#elif defined(_WINDLL) && !(defined(_WIN32) || defined(WIN32))
	s = (s != NULL ? s : buf);
	return("DOS");
#elif defined(_WIN32) || defined(WIN32)
	s = (s != NULL ? s : buf);
	return("WIN32");
#else 
	s = (s != NULL ? s : buf);
	return("N/A");
#endif
}


/*f-------------------------------------------------------------
 * System : error char string.
 *
 * <br><br>
 * Returns a string corresponding to the last system error.
 */

DLL_EXPORT(const char*)
DrlStrerror()
{
#ifdef	_WIN32
	return _strerror(NULL);
#else
	return strerror(errno);
#endif
}


