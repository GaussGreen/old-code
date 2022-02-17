/*
***************************************************************************
** SOURCE FILE: strutils.c
**
** Original authors include Doug Gallager, Bruce Broder, Peter Taylor and
** Simon Meldrum, although a lot has been removed!
**
** Various string utility routines.
***************************************************************************
*/

#include "irx/strutils.h"

#include <assert.h>
#include <ctype.h>                      /* toupper */
#include <string.h>
#include <stdarg.h>

#include "irx/error.h"
#include "irx/memutils.h"
#include "irx/macros.h"

/*
 * case indepedent string comparison
 *
 * in the system library this seems to have a different name on different
 * platforms
 *
 * using the system library version would almost certainly be quicker
 */
int irxStrncmpi(const char *str1, const char *str2, int len)
{
    while(len--)
    {
        int c1 = toupper(*str1);
        int c2 = toupper(*str2);
        if (c1 > c2) return 1;
        if (c1 < c2) return -1;
        if (c1 == 0 && c2 == 0) return 0;
        ++str1;
        ++str2;
    }
    return 0;
}

/*
 * case independent string comparison
 *
 * in the system library this seems to have a different name on different
 * platforms
 *
 * using the system library version would almost certainly be quicker
 */
int irxStrcmpi (const char *p1, const char *p2)
{
    int c1;
    int c2;
    while (toupper(*p1) == toupper(*p2))
    {
        if (*p1 == '\0') return 0;
        ++p1;
        ++p2;
    }
    /* on first different character */
    c1 = toupper(*p1);
    c2 = toupper(*p2);
    if (c1 > c2) return 1;
    if (c1 < c2) return -1;
    assert (c1 != c2); /*NOTREACHED*/
    return 0;
}

/*
***************************************************************************
** FUNCTION: irxStringConcatenate
** AUTHOR:   Peter Taylor (August 2001)
**
** Concatenates multiple strings into a single string, allocating the
** space for the output string. NULL input strings are ignored.
***************************************************************************
*/
char* irxStringConcatenate(
    size_t numStrings,   /* (I) Number of strings to concatenate */
    ...)                 /* (I) Extra arguments are all strings */
{
    static char routine[] = "irxStringConcatenate";

    va_list     ap;

    size_t      idx;
    size_t      len;
    char       *out = NULL;
    char      **strings = NULL;

    va_start (ap, numStrings);

    strings = NEW_ARRAY (char*, numStrings);
    if (strings == NULL)
        goto RETURN; /* failure */

    len = 0;
    for (idx = 0; idx < numStrings; ++idx)
    {
        strings[idx] = va_arg (ap, char*);
        if (strings[idx] != NULL)
            len += strlen(strings[idx]);
    }

    out = NEW_ARRAY (char, len+1);

    for (idx = 0; idx < numStrings; ++idx)
    {
        if (strings[idx] != NULL)
            strcat (out, strings[idx]);
    }
    
 RETURN:

    FREE (strings);
    va_end (ap);

    if (out == NULL)
        irxErrorFailure (routine);

    return out;
}


/*------------------------------------------------------------------------
FUNCTION:       irxStringToUpper

By:             Bruce Broder, 2/97

Description:    Convert a string to upper case, allocating memory for
                the output.
------------------------------------------------------------------------*/
char* irxStringToUpper( 
   const char *in)                /* (I) input string to make upper case */
{
   static char routine[] = "irxStringToUpper";
   char *out=NULL;
   char *outPtr;

   out = STRDUP( in );

   if (out==NULL)
   {
       irxError("%s Failed.\n", routine);
   }
   else
   {
       outPtr = out;
       while( (*outPtr = toupper(*outPtr)) )
           outPtr++;
   }

   return(out);
}

/*------------------------------------------------------------------------
FUNCTION:       irxStringToLower

By:             Bruce Broder, 2/97

Description:    Convert a string to lower case, allocating memory for
                the output.
------------------------------------------------------------------------*/
char* irxStringToLower( 
   const char *in)                /* (I) input string to make lower case */
{
    static char routine[] = "irxStringToLower";
    char *out;
    char *outPtr;
    
    out = STRDUP( in );
    
    if (out==NULL)
    {
        irxError("%s Failed.\n", routine);
    }
    else
    {
        outPtr = out;
        while( (*outPtr = tolower(*outPtr)) )
            outPtr++;
    }
    
    return(out);
}



/*
***************************************************************************
** FUNCTION: irxStringParser
** AUTHOR:   Peter Taylor (30 December 1998)
**
** Simple string parser.
**
** The idea is that the input string is parsed into its constituent
** components using a number of delimiters.
**
** If a delimiter is missing in the input string, then this is not an
** error. Instead all remaining constituent strings are returned as NULL.
**
** You need to provide the same number of pointers to char* as the number
** delimiters provided.
**
** The input string is modified - all the delimiters found are replaced
** with the '\0' to terminate the C-string, and the pointers returned are
** within the input string. No new memory is assigned.
**
** Therefore if you wish to preserve the input string, then you should
** take a copy of it before calling this function.
**
** Example:
**   Suppose you have a comma delimited string with three constituents.
**   Your code would be as follows (ignoring error checking):
**
**   char *str1;
**   char *str2;
**   irxStringParser (str, ",,", &str1, &str2);
**
**
** Suppose str = "abc,def,ghi"
**
** Then after this call we would have:
** 
**   str = "abc"
**   str1 = "def"
**   str2 = "ghi"
**
***************************************************************************
*/
int irxStringParser
(char *str,        /* (I/O) Input string - modified by this function */
 char *delimiters, /* (I) Delimiters */
 ...               /* (O) Need to provide strlen(delimiters) char** */
)
{
    static char routine[] = "irxStringParser";
    int         status    = FAILURE;

    va_list     ap;
    int         len;
    int         idx;
    char      **subString;
    char       *buf;        /* Current pointer */
    char       *delimiter;

    va_start (ap, delimiters);

    if (str == NULL || delimiters == NULL)
    {
        irxError ("%s: NULL inputs.\n", routine);
        goto RETURN; /* failure */
    }

    len = strlen(delimiters);

    buf = str;
    for (idx = 0; idx < len; ++idx)
    {
        subString = va_arg (ap, char**);
        if (buf != NULL)
        {
            delimiter = strchr (buf, delimiters[idx]);
            if (delimiter == NULL)
            {
                buf = NULL;
            }
            else
            {
                *delimiter = '\0'; /* Terminate previous component */
                buf = delimiter + 1;
            }
        }
        *subString = buf;
    }

    status = SUCCESS;

RETURN:

    va_end (ap);

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}

/*f
***************************************************************************
** Splits a string into components given a separator. This function
** destroys the input string by replacing the separators with the zero
** character.
**
** For example, bond.type with a separator of '.' would be split into
** "bond" and "type".
**
** Returns the number of components, and an array of strings. You must free
** the array of strings using FREE() or irxFreeSafe(). Do not use
** FREE_PTR_ARRAY to free.
**
** Note that the split array returned is an array of (numItems+1) with the
** final item being NULL. This means that you can iterate this array either
** using numItems or by testing the string pointer.
***************************************************************************
*/
int irxStringSplit
(char     *str,       /* (I) Input string */
 char      separator, /* (I) Separator - a character not a string */
 size_t   *numItems,  /* (O) Number of items found */
 char   ***split      /* (O) Array of strings found within input. User must
                           FREE() this once they have finished with it. */
)
{
    static char routine[] = "irxStringSplit";
    int         status    = FAILURE;

    size_t      myNumItems = 0;
    size_t      len;
    char        c;
    char      **mySplit = NULL;
    char       *myString;
    char       *ptr;
    size_t      pos;

    if (split != NULL) *split = NULL;
    if (numItems != NULL) *numItems = 0;
    if (str == NULL || numItems == NULL || split == NULL)
    {
        irxError ("%s: NULL inputs\n", routine);
        goto RETURN; /* failure */
    }

    len = 0;
    myNumItems = 1;
    ptr = str;
    while ((c = *ptr) != '\0')
    {
        ++len;
        if (c == separator) 
            ++myNumItems;
        ++ptr;
    }

    /*
    ** mySplit allocates all the memory that we need.
    ** This consists of the array of char*, plus a copy of the input string.
    ** Within the input string, we will replace instances of separator with
    ** the end of the string.
    */
    mySplit  = (char**)(irxMemAlloc(sizeof(char*)*(myNumItems+1)+len+1));
    if (mySplit == NULL)
        goto RETURN; /* failure */
    myString = ((char*)mySplit) + sizeof(char*)*(myNumItems+1);
    strcpy (myString, str);

    ptr = myString;
    pos = 0;

    mySplit[pos] = ptr;
    while (*ptr != '\0')
    {
        if (*ptr == separator)
        {
            *ptr = '\0';
            ++pos;
            mySplit[pos] = ptr+1;
        }
        ++ptr;
    }
    
    assert (pos+1 == myNumItems);
    
    *split    = mySplit;
    *numItems = myNumItems;
    mySplit   = NULL;
    status    = SUCCESS;

 RETURN:

    FREE(mySplit);
    if (status != SUCCESS)
        irxErrorFailure (routine);
    return status;
}


