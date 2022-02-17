/* -------------------------------------------------------------------------
        FILE_NAME:	utmessage.h
   ------------------------------------------------------------------------ */

#ifndef UTMESSAGE_H
#define UTMESSAGE_H

#define SORT_MESS_BUF_SIZE 300

/* ------------------------------------------------------------------------
        This function allows the user to return a message either to a file or
        to the interface (linking against another function) without
        interupting the program.

        On ALPHA:
        It works like vfprintf:
                prints to the default file if it has been defined

        On PC:
    It works like vsprintf  , and
                sends the string to the Trace Window
    It works like vfprintf  , and
                prints to the default file if it has been defined

        To set a message FILE  , do:
                smessage("SETFILE"  , FILE *file_path_and_name);
   ------------------------------------------------------------------------ */

void smessage(char *message, ...);

#endif
