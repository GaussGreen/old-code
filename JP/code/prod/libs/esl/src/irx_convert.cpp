/*
***************************************************************************
** convert.c
**
** Various useful conversions routines which don't really have anywhere
** else to go.
***************************************************************************
*/
#include <ctype.h>                      /* toupper */
#include <string.h>                     /* For strcmp, strtok */

#include "irx/convert.h"                    /* Prototype consistency */
#include "irx/cgeneral.h"                   /* TBoolean */
#include "irx/error.h"                      /* irxError */
#include "irx/macros.h"                     /* GTO_WRAP_CHECK_SCALAR */

/*
 * Formats a TDate. Can be called eight times from the same print
 * statement, but not more. Format is DD-MMM-YYYY.
 *
 * This function is not thread safe!
 */
char* irxFormatDate(IrxTDate date) /* (I) */
{
    static int ibuf;
#define MAX_STR_LEN 16
#define MAX_AT_ONCE 8                  /* Must be a power of 2 */
    static char format[MAX_AT_ONCE][MAX_STR_LEN];
    ibuf = (ibuf+1)&(MAX_AT_ONCE-1); /* Toggle buffers */

    return irxDateFormat (date, "DD-MMM-YYYY", &(format[ibuf][0]));
}


/*
***************************************************************************
** Extracts an array from a structure array.
**
** For example, suppose you have a structure X which contains a date and
** a double, e.g. TDate aDate and double aDouble.
**
** Then given an array of structure X (not of pointers to structure X),
** then to extract the array of dates you would do something like this:
**
** size_t arraySize;
** X structArray[arraySize];
** TDate *dateArray;
** double *doubleArray;
**
** dateArray = (TDate*) GtoArrayFromStructArray (sizeof(TDate),
**                                               sizeof(X),
**                                               &(structArray[0].aDate),
**                                               arraySize);
**
** doubleArray = (double*) GtoArrayFromStructArray (sizeof(double),
**                                                  sizeof(X),
**                                                  &(structArray[0].aDouble),
**                                                  arraySize);
***************************************************************************
*/
void* irxArrayFromStructArray
(size_t elemSize,    /* (I) Size of elements of returned array: use sizeof */
 size_t structSize,  /* (I) Size of structure: use sizeof                  */
 void  *base,        /* (I) Address of field in first item in struct array */
 size_t numElems     /* (I) Size of the array                              */
)
{
    size_t idx;
    void  *data = NULL; /* to be returned */

    if (numElems > 0)
    {
        data = irxMemAlloc (elemSize * numElems);
        if (data == NULL)
            return NULL;
    }

    for (idx = 0; idx < numElems; ++idx)
    {
        memcpy ((char*)data + elemSize*idx, 
                (char*)base + structSize*idx,
                elemSize);
    }

    return data;
}


