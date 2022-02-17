/*
***************************************************************************
** SOURCE FILE: fileutil.c
**
** CREATED BY:  Peter Taylor (10 October 1997)
**
** Various file manipulation utilities.
**
** $Header: /home/alibcvs/alibcvs/.cvsroot-alib/utils/src/fileutil.c,v 1.5 2004/09/07 12:28:26 ptaylor Exp $
***************************************************************************
*/

#include "irx/fileutil.h"

#include <string.h>     /* strcpy */

#include "irx/error.h"     /* irxErrMsg */
#include "irx/macros.h"     /* NEW_ARRAY etc */
#include "irx/strutils.h"   /* irxStringDuplicate */

/*
***************************************************************************
** FUNCTION: irxFileReadAll
**
** AUTHOR:   Peter Taylor (10 October 1997)
**
** Reads a whole text file and returns a character string. This string must
** be free'd by the user.
***************************************************************************
*/
char * irxFileReadAll
(char *filename /* (I) File name to read */
)
{
    static char routine[] = "irxFileReadAll";
    int         status    = FAILURE;

    char        *str=NULL;       /* To be returned */

    FILE        *fp = NULL;

/* Open the input file */
    fp = fopen (filename, "r");
    if (fp == NULL)
        goto RETURN; /* failure */

/*
** Use irxOpenFileReadAll to read the file
*/
    str = irxOpenFileReadAll (fp, filename, NULL);
    if (str == NULL)
        goto RETURN; /* failure */

    status = SUCCESS;

RETURN:

    FCLOSE(fp);

    if (status != SUCCESS)
    {
        FREE_ARRAY(str);
        str = NULL;
        irxError("%s: Failed.\n", routine);
    }

    return str;
}








/*f
***************************************************************************
** FUNCTION: irxOpenFileReadAll
** AUTHOR:   Peter Taylor (August 2001)
**
** Reads a whole text file which has been previously opened and returns
** a character string. This string must be free'd by the user, and the
** file closed by the caller.
**
** Optionally you can provide the text that has been read in so far and
** this will be attached to the beginning of the return string. The idea
** behind this is to allow you to detect the file type by just reading
** the beginning of the file, and then follow-up by reading the entire
** file.
***************************************************************************
*/
char* irxOpenFileReadAll
(FILE *fp,         /* (I) File to read */
 char *filename,   /* (I) For error messages */
 char *initialText /* (I) Text that has been read so far - can be NULL */
)
{
    static char routine[] = "irxOpenFileReadAll";
    int         status    = FAILURE;

    char        *str=NULL;       /* To be returned */

    char         buf[512];
    int          bufsize = sizeof(buf);
    long         allocLen;
    long         len;
    int          buflen; /* Actual number of bytes read */
    int          initlen; /* Size of initial text */

    if (initialText != NULL)
        initlen = strlen(initialText);
    else
        initlen = 0;

/*
** Read the file in chunks into the fixed size buffer.
*/
    buflen = fread (buf, 1, bufsize-1, fp);
    if (buflen < 1)
        goto RETURN; /* failure */

    buf[buflen] = '\0';
    allocLen = bufsize+initlen;
    str = NEW_ARRAY(char, allocLen+1);
    if (str == NULL)
        goto RETURN; /* failure */

    /*
    ** Join together the initial text, the last character from file and
    ** the buffer just read from file.
    */
    if (initlen > 0)
    {
        strcpy (str, initialText);
    }
    strcpy (str+initlen, buf);

    len = strlen(str);

    while (buflen = fread (buf, 1, bufsize-1, fp), 
           buflen > 0)
    {
        if ((len+buflen) > allocLen)
        {
            char *tmp = NULL;

            allocLen *= 2;

            /*
            ** Since allocLen will always be increasing by at least bufsize,
            ** it should never be the case that (len+buflen) > allocLen.
            ** It will only be the case when allocLen has overflowed.
            */
            if ((len+buflen) > allocLen)
            {
                irxError ("%s: File %s is too large to read as one string.\n",
                           routine, filename);
                goto RETURN; /* failure */
            }
            
            tmp = NEW_ARRAY (char, allocLen+1);
            if (tmp == NULL)
                goto RETURN; /* failure */

            strcpy (tmp, str);
            FREE_ARRAY (str);
            str = tmp;
        }

        buf[buflen] = '\0';
        strcpy (str+len, buf);
        len += strlen(buf);
    }

    status = SUCCESS;

RETURN:

    if (status != SUCCESS)
    {
        FREE_ARRAY(str);
        str = NULL;
        irxError("%s: Failed.\n", routine);
    }

    return str;
}

/*
***************************************************************************
** Returns the point at which a file name splits into dirname and basename.
** If NULL is returned, it means that there is no split and the file name
** is already a base name.
***************************************************************************
*/
static char* splitFileNamePosition(char* fullFileName)
{
    char       *c1;
    char       *c2;

    c1 = strrchr (fullFileName, '\\');
    c2 = strrchr (fullFileName, '/');

    if (c1 == NULL && c2 == NULL)
        return NULL;
    if (c1 == NULL)
        return c2;
    if (c2 == NULL)
        return c1;
    if (c1 < c2)
        return c2;
    return c1;
}


/*f
***************************************************************************
** Returns the basename of a file from the full file path.
**
** Will accept either Unix separators ("/") or Windows separators ("\") in
** the file name.
**
** You should free the resulting string using FREE.
***************************************************************************
*/
char* irxBasename (char* fullFileName /* (I) */)
{
    static char routine[] = "irxBasename";
    
    char       *split;
    char       *basename;

    if (fullFileName == NULL)
    {
        irxError ("%s: NULL inputs\n", routine);
        return NULL;
    }

    split = splitFileNamePosition (fullFileName);
    if (split == NULL)
        basename = fullFileName;
    else
        basename = split+1;

    return STRDUP (basename);
}
            
/*f
***************************************************************************
** Returns the directory of a file from the full file path.
**
** Will accept either Unix separators ("/") or Windows separators ("\") in
** the file name. Does not change the separators within the file name.
**
** If there is no directory in the name, then returns a string of length
** zero rather than NULL - a NULL return is always an error case.
**
** You should free the resulting string using FREE.
***************************************************************************
*/
char* irxDirname (char* fullFileName /* (I) */)
{
    static char routine[] = "irxDirname";
    
    char       *dirname;
    char       *split;

    if (fullFileName == NULL)
    {
        irxError ("%s: NULL inputs\n", routine);
        return NULL;
    }

    dirname = STRDUP (fullFileName);
    if (dirname == NULL)
        return NULL; /* an error case */

    split = splitFileNamePosition (dirname);
    if (split == NULL)
        dirname[0] = '\0';
    else
        *split = '\0';

    return dirname;
}






