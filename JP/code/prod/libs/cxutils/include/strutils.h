/*
***************************************************************************
** HEADER FILE: strutils.h
**
** Some string utility routines.
**
** $Header$
***************************************************************************
*/

#ifndef CX_STRUTILS_H
#define CX_STRUTILS_H

#include "defines.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*f
***************************************************************************
** Case independent string comparison (fixed length)
**
** In the system library this seems to have a different name on different
** platforms
**
** Using the system library version would almost certainly be quicker.
***************************************************************************
*/
int CxStrncmpi(const char *str1, const char *str2, int len);

/*f
***************************************************************************
** Case independent string comparison
**
** In the system library this seems to have a different name on different
** platforms
**
** Using the system library version would almost certainly be quicker.
***************************************************************************
*/
int CxStrcmpi (const char *str1, const char *str2);


/*f
***************************************************************************
** Concatenates multiple strings into a single string, allocating the
** space for the output string. NULL input strings are ignored.
***************************************************************************
*/
char* CxStringConcatenate(
    size_t numStrings,   /* (I) Number of strings to concatenate */
    ...);                /* (I) Extra arguments are all strings */

/*f
***************************************************************************
** Convert a string to upper case, allocating memory for the output.
***************************************************************************
*/
char* CxStringToUpper( 
    const char *in);               /* (I) input string to make upper case */

/*f
***************************************************************************
** Convert a string to lower case, allocating memory for the output.
***************************************************************************
*/
char* CxStringToLower( 
    const char *in);               /* (I) input string to make lower case */

/*f
***************************************************************************
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
**   CxStringParser (str, ",,", &str1, &str2);
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
int CxStringParser
(char *str,        /* (I/O) Input string - modified by this function */
 char *delimiters, /* (I) Delimiters */
 ...               /* (O) Need to provide strlen(delimiters) char** */
);

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
** the array of strings using FREE() or CxFreeSafe(). Do not use
** FREE_PTR_ARRAY to free.
**
** Note that the split array returned is an array of (numItems+1) with the
** final item being NULL. This means that you can iterate this array either
** using numItems or by testing the string pointer.
***************************************************************************
*/
int CxStringSplit
(char     *str,       /* (I) Input string */
 char      separator, /* (I) Separator - a character not a string */
 size_t   *numItems,  /* (O) Number of items found */
 char   ***split      /* (O) Array of strings found within input. User must
                           FREE() this once they have finished with it. */
);

#ifdef __cplusplus
}
#endif

#endif
