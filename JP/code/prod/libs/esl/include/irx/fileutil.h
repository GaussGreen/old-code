/*
***************************************************************************
** HEADER FILE: fileutil.h
**
** CREATED BY:  Peter Taylor (10 October 1997)
**
** Various file manipulation utilities.
**
** $Header: /home/alibcvs/alibcvs/.cvsroot-alib/utils/include/fileutil.h,v 1.5 2004/09/07 12:28:26 ptaylor Exp $
***************************************************************************
*/

#ifndef IRX_FILEUTIL_H
#define IRX_FILEUTIL_H

#include "cgeneral.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*f
***************************************************************************
** Reads a whole text file and returns a character string. This string must
** be free'd by the user.
***************************************************************************
*/
extern char * irxFileReadAll
(char *fileName /* (I) File name to read */
);



/*f
***************************************************************************
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
);


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
char* irxBasename (char* fullFileName /* (I) */);

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
char* irxDirname (char* fullFileName /* (I) */);

#ifdef __cplusplus
}
#endif

#endif    /* IRX_FILEUTIL_H */










