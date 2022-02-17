/*
***************************************************************************
** SOURCE FILE: DateListMakeRegular.c
**
** Tests the binary search routines.
**
** $Header: /nasdev/export2/home/drdev/cvsadmin/cvs/credit-hybrids/test/drivers/src/DateListMakeRegular.c,v 1.1 2004/12/06 20:10:43 ptaylor Exp $
***************************************************************************
*/

#include <assert.h>
#include <string.h>

#include "cxmacros.h"
#include "datelist.h"
#include <alib/cerror.h>
#include <alib/convert.h>
#include <alib/dtlist.h>

/* we use main.inc and hence we get inputs as regular C-arrays instead
   of as counted arrays */
#include <alib/main.inc>

/* define the input file */

ARG_DEF g_ArgInfoArray[]=
{
    {"startDate:",       AT_DATE},
    {"endDate:",         AT_DATE},
    {"interval:",        AT_STRING},
    {"stubType:",        AT_STRING}
};

long g_NumArgs = sizeof(g_ArgInfoArray)/sizeof(ARG_DEF);

void callRoutine(FILE *fp)
{
    static char routine[] = "DateListMakeRegular";
    int         status    = FAILURE;

    TDate  startDate   = *(TDate*)TestArgCacheGet(0);
    TDate  endDate     = *(TDate*)TestArgCacheGet(1);
    char*  strInterval = *(char**)TestArgCacheGet(2);
    char*  strStubType = *(char**)TestArgCacheGet(3);

    TDateInterval interval;
    CxTStubType   stubType;

    TDateList *dl = NULL;

    if (GtoStringToDateInterval (strInterval, "interval",
                                 &interval) != SUCCESS)
        goto done; /* failure */

    if (CxStubTypeFromString (strStubType, &stubType) != SUCCESS)
        goto done; /* failure */

    dl = CxDateListMakeRegular (startDate, endDate, &interval, stubType);
    if (dl == NULL)
        goto done; /* failure */

    GtoPrintDateList (dl, routine);
    
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    GtoFreeDateList (dl);
}



