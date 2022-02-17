/* -------------------------------------------------------------------------

   File: %M%
   Path: %P%
   Description: utilitaires pour les signaux
   Created: 97/03/19
   Author: Jacques WERNERT
   Modified: %E% %U%
   Last maintained by: Jacques WERNERT
   Revision: %I%

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#include <CCcommon.h>
SCCS_ID (sigutil_c_SccsId, "%W%, modified %E%");

#include "CCsigutil.h"

/*
 **      Function name: COM_get_signal_str
 **
 **      Description: affectete une chaine imprimable a un numero de signal
 **      Input: signal
 **      Output: chaine STATIQUE
 */
char *COM_get_signal_str (signal)
     int signal;
{
  char *ret = DEFAULT_NULL_STRING;
 
  switch (signal) {
    
  case SIGINT:
    ret = SIGINT_STR;
    break;

#ifndef WIN32
  case SIGQUIT:
    ret = SIGQUIT_STR;
    break;
#endif /* WIN32 */
    
  case SIGFPE:
    ret = SIGFPE_STR;
    break;

#ifndef WIN32
  case SIGKILL:
    ret = SIGKILL_STR;
    break;
#endif /* WIN32 */

#ifndef WIN32
  case SIGBUS:
    ret = SIGBUS_STR;
    break;
#endif /* WIN32 */

  case SIGSEGV:
    ret = SIGSEGV_STR;
    break;
    
#ifndef WIN32
  case SIGSYS:
    ret = SIGSYS_STR;
    break;
#endif /* WIN32 */
    
#ifndef WIN32
  case SIGPIPE:
    ret = SIGPIPE_STR;
    break;
#endif /* WIN32 */
    
#ifndef WIN32
  case SIGALRM:
    ret = SIGALRM_STR;
    break;
#endif /* WIN32 */
    
  case SIGTERM:
    ret = SIGTERM_STR;
    break;

#ifndef WIN32
  case SIGCHLD:
    ret = SIGCHLD_STR;
    break;
#endif /* WIN32 */
    
#ifndef WIN32
  case SIGSTOP:
    ret = SIGSTOP_STR;
    break;
#endif /* WIN32 */
    
#ifndef WIN32
  case SIGCONT:
    ret = SIGCONT_STR;
    break;
#endif /* WIN32 */
    
#ifndef WIN32
  case SIGTTIN:
    ret = SIGTTIN_STR;
    break;
#endif /* WIN32 */
    
#ifndef WIN32
  case SIGTTOU:
    ret = SIGTTOU_STR;
    break;
#endif /* WIN32 */
    
  default:
    ret = SIGUNKNOWN_STR;
    
  } /* switch */
  
  return (ret);

} /* COM_get_signal_str */

/* EOF %M% */

