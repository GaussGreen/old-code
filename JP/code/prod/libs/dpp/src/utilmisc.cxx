/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:	$Header$
 ************************************************************************/
#define	_dpp_SRC
#include <errno.h>
#include "kstdinc.h"


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


#undef	__DEBUG__
#define	__DEBUG__








