/* -------------------------------------------------------------------------

   File: CCsigutil.h
   Path: /home/julia/projets/dev/Common/libcommon/SCCS/s.CCsigutil.h
   Description: utilitaires de manipulation des signaux
   Created: 97/03/19
   Author: Jacques WERNERT
   Modified: 00/06/19 17:32:43
   Last maintained by: Jacques WERNERT
   Revision: 1.6

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#ifndef _sigutil_h
#define _sigutil_h

#include <CCcommon.h>
SCCS_ID (sigutil_h_SccsId, "@(#)CCsigutil.h	1.6, modified 00/06/19");

#include <CCmessage.h>
#include <sys/types.h>
#include <signal.h>

/*
 *******************
 * Partie publique *
 *******************
 */
CCEXTERN_FUNCTION (char *COM_get_signal_str, (int signal));

/*
 *****************
 * Partie privee *
 *****************
 */
#ifdef _sigutil_c

#endif _sigutil_c

#endif _sigutil_h

/* EOF CCsigutil.h */

