/* -------------------------------------------------------------------------

   File: CCdebug.c
   Path: /home/julia/projets/dev/Common/libdebug/SCCS/s.CCdebug.c
   Description: debug library source file
   Created: 96/09/17
   Author: Jacques WERNERT
   Modified: 98/11/05 13:59:45
   Last maintained by: Jacques WERNERT
   Revision: 1.15

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#include <CCcommon.h>
SCCS_ID (debug_c_SccsId, "@(#)CCdebug.c	1.15, modified 98/11/05");

#define _debug_c
#include "CCdebug.h"

/*
 **      Function name : debug_mode
 **
 **      Description : set/unset/get debug mode
 **      Input : FDEBUG_CMD_TYPE type
 **      Output : always return current mode
*/
static int debug_mode (cmd)
{
#define FDEBUG_ON  1
#define FDEBUG_OFF 0

#ifdef DEBUG
  static int mode = FDEBUG_ON;
#else
  static int mode = FDEBUG_OFF;
#endif
  switch (cmd) {
  case FDEBUG_SET:
    mode = FDEBUG_ON;
    break;
  case FDEBUG_UNSET:
    mode = FDEBUG_OFF; 
    break;
  case FDEBUG_GET:
    break;	
  }
  return (mode==FDEBUG_ON);
#undef FDEBUG_ON
#undef FDEBUG_OFF
}

/*
 **      Function name : DEB_fdebug_parameters
 **
 **      Description : client hidden function. Sets line, file and debug mode
 **      Input : FDEBUG_CMD_TYPE type 
 **      Output : always return current mode
*/
#ifdef __STDC__
int DEB_fdebug_parameters (FDEBUG_CMD_TYPE cmd, ...)
#else
int DEB_fdebug_parameters (cmd, va_alist)
     FDEBUG_CMD_TYPE cmd;
     va_dcl
#endif
{
  va_list args;
  int status_code = 1;
  static char *file;
  static int line;
 
  CC_VA_START (args, cmd); 

  switch (cmd) {
  case FDEBUG_SET:
    file = va_arg(args, char *);
    line = va_arg(args, int);
    if (va_arg(args, int) == FDEBUG_SET)
      debug_mode (FDEBUG_SET);
    break;
    
  case FDEBUG_GET:
    *(va_arg(args, char **)) = file;
    *(va_arg(args, int *)) = line;
    break;
    
  default:
    status_code = 0;
  }
  va_end(args);
  return(status_code);
}

#ifdef __STDC__
int DEB_the_fdebug(FILE *stream, char *format, ...)
#else
int DEB_the_fdebug(stream, format, va_alist)
     FILE *stream;
     char *format;
     va_dcl
#endif
{
  va_list args;
  int ret_code;

  CC_VA_START (args, format);

  ret_code = DEB_the_vfdebug (stream, format, args);

  va_end(args);
	
  return (ret_code);
}

/*
 **      Function name : DEB_fprint_file_line
 **
 **      Description : 
 **      Input : flux
 **      Output : 'fprintf' return code
*/
int DEB_fprint_file_line (stream)
     FILE *stream;
{
  char *file;
  int line;

  /*
   * Impression du fichier et de la ligne concernee
   */    
  DEB_fdebug_parameters (FDEBUG_GET, &file, &line);
  return (fprintf (stream, "[%s@%i] ", file, line));
}

/*
 **      Function name : DEB_the_vfdebug
 **
 **      Description : client function. The 'debug' function
 **      Input : similar to 'fprintf' but prefix each line by the date, file
 **			and line in which this function was called
 **      Output : last internal fprintf call / similar to 'fprintf' return code
*/
int DEB_the_vfdebug (stream, format, args)
     FILE *stream;
     char *format;
     va_list args;
{
  time_t timestamp;
  char *heurelocale;
  int slen, error;
  
  error = 0;

  /*
   * Si on est en mode debug ...
   */
  if (debug_mode (FDEBUG_GET)) {

    /*
     * Ajout d'un entete
     */
    DEB_ASSERT (time(&timestamp));
    
    heurelocale = ctime(&timestamp);
    DEB_ASSERT(heurelocale);
    ((slen=strlen(heurelocale)) > 0)?(heurelocale[slen-1]='\0'):(heurelocale[0]='\0');
    fprintf (stream, "[%s] ", heurelocale);	

    /*
     * Impression du fichier et de la ligne concernee
     */    
    DEB_fprint_file_line (stream);

    vfprintf (stream, format, args);
    error = fprintf (stream, "\n");
    fflush(stream);

  }
 
  
  return(error);
}

/*
 **      Function name : DEB_flog
 **
 **      Description : client function. The 'log' function
 **      Input : similar to 'fprintf' buf prefix each line by the current date
 **      Output : last internal fprintf call / similar to 'fprintf' return code
 */
#ifdef __STDC__
int DEB_flog(FILE *stream, char *format, ...)
#else
int DEB_flog(stream, format, va_alist)
     FILE *stream;
     char *format;
     va_dcl
#endif
{
  va_list args;
  time_t timestamp;
  char *heurelocale;
  int slen, error;
  
  CC_VA_START (args, format);

  DEB_ASSERT(time(&timestamp));
  
  heurelocale = ctime(&timestamp);
  DEB_ASSERT(heurelocale);
  ((slen=strlen(heurelocale)) > 0)?(heurelocale[slen-1]='\0'):(heurelocale[0]='\0');
  fprintf (stream, "[%s] ", heurelocale);	
  
  vfprintf (stream, format, args);
  error = fprintf (stream, "\n");
  fflush (stream);
  va_end(args);
  
  return(error);
}

/*
 **      Function name : DEB_fdebug_init
 **
 **      Description : client function. The 'init' function for 'fdebug'
 **      Input : 'argc' and 'argv' which are the 'main' arguments plus
 **              a null terminated list of string in which each element
 **              is a command line argument which tells the library to
 **              turn on the debug mode.
 **      Output : the debug mode
*/
#ifdef __STDC__
int DEB_fdebug_init (int argc, char *argv[], ...)
#else
int DEB_fdebug_init (argc, argv, va_alist)
     int argc;
     char *argv[];
     va_dcl
#endif
{
  va_list args;
  char *arguments;
  int i;
  
  CC_VA_START (args, argv);

  while ((arguments=va_arg(args, char *)))
    for (i=1; i < argc; i++)
      if (!strcmp(arguments, argv[i]))
	debug_mode (FDEBUG_SET);
  
  va_end(args);
  return(debug_mode(FDEBUG_GET));
}

/* EOF CCdebug.c */
