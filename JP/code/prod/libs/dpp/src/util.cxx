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



/*f---------------------------------------------------------------------
 * Executes a command.
 */

void
DppExecCmd(const char *fmt,...)
{
static	char	routine[]= "DppExecCmd";
	va_list	ap;
static	char	tmpBuf[1024];

	va_start(ap, fmt);
	vsprintf(tmpBuf, fmt, ap);
	va_end(ap);

	//dppLog << routine << ": exec `" << tmpBuf << '.' << endl);

	if (system(tmpBuf) == -1) {
		throw KFailure("%s: exec `%s' failed (%s).\n",
			routine, tmpBuf, strerror(errno));
	}
}

