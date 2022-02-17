/* -------------------------------------------------------------------------

   File: CCmessage.c
   Path: /home/julia/projets/dev/Common/libmessage/SCCS/s.CCmessage.c
   Description: Librairie de gestion des messages d'erreurs
   Created: 96/10/10
   Author: Jacques WERNERT
   Modified: 99/04/28 12:12:51
   Last maintained by: Jacques WERNERT
   Revision: 1.24

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#include <CCcommon.h>
SCCS_ID (CCmessage_c_SccsId, "@(#)CCmessage.c	1.24, modified 99/04/28");

#define _CCmessage_c
#include "CCmessage.h"

int MSG_verbose = MSG_DISABLED;
int MSG_stderr_quiet = MSG_DISABLED;    /* quiet sur stderr */

#ifdef DEBUG
  int MSG_debug = MSG_ENABLED;
#else
  int MSG_debug = MSG_DISABLED;
#endif DEBUG

// FIXMEFRED: mig.vc8 (25/05/2007 10:44:57): stderr is not a constant
static FILE *stdlog /* = stderr */;
static char *message_entete = "";
static int allocated_message_entete = 0;

/*
 * Indique si la fin de ligne doit etre prefixee
 * d'un \r
 */
int MSG_carriage_return = MSG_DISABLED;

/*
 **      Function name: MSG_print_eol
 **
 **      Description: termine une impression / permet de gerer le \r et flush
 **      Input: flux ouvert en ecriture
 **      Output: id a fprintf
 */
int MSG_print_eol (flux)
     FILE *flux;
{

  int ret;

  if (MSG_carriage_return == MSG_ENABLED) {
    ret = fprintf (flux, "\r");
  }
  ret = fprintf (flux, "\n");

  return (ret);

} /* MSG_print_eol */

/*
 **      Function name: MSG_set_header
 **
 **      Description: positionne/recupere la chaine de description
 **      Input: chaine qui sera dupliquee ou NULL
 **      Output: la chaine utilisee
 */
char *MSG_set_header (entete)
     char *entete;
{
  if (entete) {
    if (allocated_message_entete) free (message_entete);
    message_entete = (char *)strdup (entete);
    allocated_message_entete = 1;
  }

  return (message_entete);

} /* MSG_set_header */

/*
 **      Function name: MSG_set_stdlog
 **
 **      Description: declare le flux de log
 **      Input: flux ouvert en ecriture
 **      Output: none
 */
void MSG_set_stdlog (new_stdlog)
    FILE *new_stdlog;
{
  if (new_stdlog)
     stdlog = new_stdlog;
}

/*
 **      Function name: MSG_the_printf_message
 **
 **      Description: la fonction d'impression, interface
 **      Input: type, format, args
 **      Output: id que fprintf
 */
#ifdef __STDC__
int MSG_the_printf_message (int message_type, const char *format, ...)
#else
int MSG_the_printf_message (message_type, format, va_alist)
     char *format;
     va_dcl
#endif __STDC__
{
  va_list args;
  int ret_code;

  CC_VA_START (args, format);
  ret_code = MSG_the_vprintf_message (message_type, format, args);
  va_end (args);

  return (ret_code);

} /* MSG_the_printf_message */

/*
 **      Function name: MSG_the_vprintf_message
 **
 **      Description: implementation de la fonction d'impression
 **      Input: type, format, va_list
 **      Output: id a fprintf
 */
#ifdef __STDC__
int MSG_the_vprintf_message (int message_type, const char *format, va_list args)
#else
int MSG_the_vprintf_message (message_type, format, args)
     char *format;
     va_list args;
#endif __STDC__
{

  char *heurelocale = NULL;
  static char buf[DAT_STRING_SIZE];

#define GET_HEURE_LOCALE heurelocale ? heurelocale : (heurelocale = DAT_gmt_to_mystring (DAT_now, buf, DAT_STRING_SIZE))

  switch (message_type) {

  case MSG_TRACE:

    /* 
     * mode debug: ecriture sur stdlog
     * mode verbose: ecriture sur stdlog si non quiet sur stderr
     */
    if ((MSG_debug) || 
        ((MSG_verbose) && (! MSG_must_be_quiet_on_stderr (stdlog)))
	) {
      fprintf (stdlog, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_TRACE_STR, 
	       message_entete);
      if (MSG_debug)
        DEB_fprint_file_line (stdlog);
      
      vfprintf (stdlog, format, args);
      MSG_print_eol (stdlog);
    }
    
    /*
     * ecriture du stderr en mode debug
     *    si les flux de log et d'erreur sont differents
     */
    if ((MSG_debug) && (stdlog != stderr)) {
      fprintf (stderr, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_TRACE_STR, 
	       message_entete);
      DEB_fprint_file_line (stderr);
      vfprintf (stderr, format, args);
      MSG_print_eol (stderr);
    }
    break;
      
  case MSG_INFO:

    /* 
     * messages d'informations appli
     * ecriture sur stdlog
     *   si il ne faut pas causer sur stderr (debug non compris)
     */
    if ((MSG_debug) || (! MSG_must_be_quiet_on_stderr (stdlog))) {

      fprintf (stdlog, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_INFO_STR,
	       message_entete);
      
      if (MSG_debug)
	DEB_fprint_file_line (stdlog);
      
      vfprintf (stdlog, format, args);
      MSG_print_eol (stdlog);
    }
    
    /*
     * ecriture du stderr
     *   si pas deja fait sur stdlog
     *   et si il ne faut pas causer sur stderr (debug non compris)
     */
    if ( (stdlog != stderr) && ((MSG_debug) || (! MSG_stderr_quiet)) ) {
      fprintf (stderr, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_INFO_STR,
	       message_entete);
      
      if (MSG_debug)
        DEB_fprint_file_line (stderr);
      
      vfprintf (stderr, format, args);
      MSG_print_eol (stderr);
    }
    
    break;
    
  case MSG_ERROR:
  case MSG_ERROR_DETAIL:
    
    /*
     * messages d'erreurs non fatales
     * ecriture sur stdlog
     *   si il ne faut pas causer sur stderr (debug non compris)
     */
    if ((MSG_debug) || (! MSG_must_be_quiet_on_stderr (stdlog))) {

      fprintf (stdlog, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_ERROR_STR,
	       message_entete);
    
      if (MSG_debug)
	DEB_fprint_file_line (stdlog);
    
      vfprintf (stdlog, format, args);
      MSG_print_eol (stdlog);
    
    }

    /*
     * ecriture du stderr
     *   si pas deja fait sur stdlog
     *   et si il ne faut pas causer sur stderr (debug non compris)
     */
    if ( (stdlog != stderr) && ((MSG_debug) || (! MSG_stderr_quiet)) ) {

      fprintf (stderr, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_ERROR_STR,
	       message_entete);	
	
      if (MSG_debug)
	DEB_fprint_file_line (stderr);
      
      vfprintf (stderr, format, args);
      MSG_print_eol (stderr);
    }
    
    break;
    
  case MSG_FATAL:

    /*
     * erreurs fatales
     * ecriture sur stdlog
     *   si il ne faut pas causer sur stderr (debug non compris)
     */
    if ((MSG_debug) || (! MSG_must_be_quiet_on_stderr (stdlog))) {
      
      fprintf (stdlog, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_FATAL_STR,
	       message_entete);	
      
      if (MSG_debug)
	DEB_fprint_file_line (stdlog);
      
      vfprintf (stdlog, format, args);
      MSG_print_eol (stdlog);
    
    }

    /*
     * ecriture du stderr
     *   si pas deja fait sur stdlog
     *   et si il ne faut pas causer sur stderr (debug non compris)
     */
    if ( (stdlog != stderr) && ((MSG_debug) || (! MSG_stderr_quiet)) ) {

      fprintf (stderr, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_FATAL_STR,
	       message_entete);	
      
      if (MSG_debug)
        DEB_fprint_file_line (stderr);
      
      vfprintf (stderr, format, args);
      MSG_print_eol (stderr);
    }

    break;
    
  default:

    /*
     * erreur fatale, code inconnu
     * ecriture sur stdlog
     *   si il ne faut pas causer sur stderr (debug non compris)
     */
    if ((MSG_debug) || (! MSG_must_be_quiet_on_stderr (stdlog))) {

      fprintf (stdlog, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_FATAL_STR,
	       message_entete);	
      
      if (MSG_debug)
	DEB_fprint_file_line (stdlog);
      
      vfprintf (stdlog, format, args);
      MSG_print_eol (stdlog);
    }

    /*
     * ecriture du stderr
     *   si pas deja fait sur stdlog
     *   et si il ne faut pas causer sur stderr (debug non compris)
     */
    if ( (stdlog != stderr) && ((MSG_debug) || (! MSG_stderr_quiet)) ) {

      fprintf (stderr, "[%s] [%s] [%s] ", GET_HEURE_LOCALE, MSG_FATAL_STR,
               message_entete);	
      
      if (MSG_debug)
        DEB_fprint_file_line (stderr);
      
      vfprintf (stderr, format, args);
      MSG_print_eol (stderr);

    }

    CC_ASSERT (0, MSG_F_MIS_INCONSISTENT_VALUE);

  } /* switch */

  return(1);

} /* MSG_the_vprintf_message */

/* EOF CCmessage.c */
