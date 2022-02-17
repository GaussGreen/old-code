/* -------------------------------------------------------------------------
        FILE_NAME:	utmessage.cxx
   ------------------------------------------------------------------------ */

#include "utallhdr.h"

/* ------------------------------------------------------------------------
        This function allows the user to return a message either to a file or
        to the interface (linking against another function) without
        interupting the program.

        On ALPHA:
        It works like vfprintf:
                prints to the default file if it has been defined

        On PC:
    It works like vsprintf        , and
                sends the string to the Trace Window
    It works like vfprintf        , and
                prints to the default file if it has been defined

        To set a message FILE        , do:
                smessage("SETFILE"        , FILE *file_path_and_name);
   ------------------------------------------------------------------------ */

#ifdef _WINDOWS
#include "traceapi.h"
#include "windows.h"
#endif

void smessage(String message, ...) {
  static char SRT_MESSAGE_OUTPUT_STRING[SORT_MESS_BUF_SIZE];
  static FILE *fout = NULL;
  static String sout = NULL;
  va_list args;

/* For a PC environment        , sets the string */
#ifdef _WINDOWS
  sout = (String)SRT_MESSAGE_OUTPUT_STRING;
#endif

  va_start(args, message);

  if (!strcmp(message, "SETFILE"))
  /* Set default error file if required */
  {
    fout = (FILE *)(va_arg(args, FILE *));
    va_end(args);
  }

  if (fout)
  /* Prints the message to the message file if defined */
  {
    fprintf(fout, "\n");
    vfprintf(fout, message, args);
    fprintf(fout, "\n");
    va_end(args);
  }

  if (sout)
  /* Prints the message to the message string if defined */
  {
    vsprintf(sout, message, args);
#ifdef _WINDOWS
    EMTraceMsg(sout);
#endif
  }
}