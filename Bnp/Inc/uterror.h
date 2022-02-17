/* ======================================================
   FILENAME :  uterror.h

   PURPOSE  :  Functions to handle and return errors
   ====================================================== */
/*
        serror(fmt  ,...) works like printf; it prints to a file  , (Default
   NULL) and also to a string  , which it returns.

        To change default string  , error("SETSTRING"  ,String
   new_default_string); To change default FILE  , error("SETFILE"  ,FILE
   *new_default_file);

        To return an error:
        return serror("problem: rate %d cannot be negative\n"  ,i);
*/

#ifndef SRT_H_UTLERR_H
#define SRT_H_UTLERR_H

#define SORT_ERR_BUF_SIZE 300

typedef char *Err;
typedef char *SrtErr;
// dnewo'mwedenwoewnewnoipewpoewp[wew[pwemjepewop

Err serror(char *format, ...);

#endif
