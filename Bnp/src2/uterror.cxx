/* ==========================================================================
   FILE_NAME:	uterror.cxx

   PURPOSE:     return errors as strings...
   ========================================================================== */

#include "utallhdr.h"

/* ------------------------------------------------------------------------
        Works like sprintf        , except the string is defined internally.
        Can be used to return an error message with extra arguments...
        See Kerninghah & Ritchie page 174 for the code inspiration.
   ------------------------------------------------------------------------ */
static char _the_serror_string[256];

Err serror(char *format, ...)

{
  va_list args;
  static char temp[256];

  va_start(args, format);
  vsprintf(temp, format, args);
  va_end(args);

  strcpy(_the_serror_string, temp);

  return _the_serror_string;
}
