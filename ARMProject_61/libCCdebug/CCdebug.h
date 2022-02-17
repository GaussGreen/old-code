/* -------------------------------------------------------------------------

   File: CCdebug.h
   Path: /home/julia/projets/dev/Common/libdebug/SCCS/s.CCdebug.h
   Description: debug library header file
   Created: 96/09/18
   Author: Jacques WERNERT
   Modified: 99/07/15 17:11:55
   Last maintained by: Jacques WERNERT
   Revision: 1.10

   -------------------------------------------------------------------------

   Note:
	  3 fonctions sont definies dans cette librairie:

      o flog. Cette fonction dont les parametres sont identiques a 'fprintf'
        permet d'ecrire un message prefixe par l'heure courante.
        Le retour chariot est implicite.

      o fdebug. Cette fonction permet, en mode 'debug' d'imprimer un message
        a la 'flog' avec, en plus, le nom de fichier et la ligne en
        cours. Elle prend les memes types d'arguments que 'fprintf'
        Le retour chariot est implicite.

      o fdebug_init. Cette fonction permet de basculer en mode debug par
        parametre de ligne de commande. Cette fonction recoit en
        'argc' et 'argv' du 'main' ainsi que toutes les valeurs
        possibles du parametre de ligne de commande pour le mode
        'debug'. La liste doit etre 'NULL terminated'.

   Exemple:
        #include <stdio.h>
        #include <CCdebug.h>

        #ifdef __STDC__
          int main(int argc, char *argv[])
        #else
          int main(argc, argv)
                int argc;
                char *argv[];
        #endif
        {
            fdebug_init(argc, argv, "-debug", "-d", NULL);
            flog(stderr, "Un essai");
            flog(stderr, "Essai no %i", 2);
            fdebug(stderr, "Un essai");
            fdebug(stderr, "Essai no %i", 2);
            return(0);
        }

        % flog_test
        [Wed Sep 18 08:46:47 1996] Un essai
        [Wed Sep 18 08:46:47 1996] Essai no 2

        % flog_test -d
        [Wed Sep 18 08:47:08 1996] Un essai
        [Wed Sep 18 08:47:08 1996] Essai no 2
        [Wed Sep 18 08:47:08 1996] [flog_test.c@15] Un essai
        [Wed Sep 18 08:47:08 1996] [flog_test.c@16] Essai no 2

   Compilation:
        Si le flag DEBUG est positionne, le mode 'debug' pour la fonction
        'fdebug' est force.

   ------------------------------------------------------------------------- */

#ifndef _debug_h
# define _debug_h

#include <CCcommon.h>
SCCS_ID (debug_h_SccsId, "@(#)CCdebug.h	1.10, modified 99/07/15");

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>

#define DEB_ASSERT(x) if (!(x)) { fprintf(stderr, "CCdebug.h library. Fatal Error\n"); assert((x)); }

#define DEB_fdebug DEB_the_fdebug
#define DEB_vfdebug DEB_the_vfdebug

CCEXTERN_FUNCTION (int DEB_flog,
		 (FILE *stream, char *format, ...));
CCEXTERN_FUNCTION (int DEB_fdebug_init,
		 (int argc, char *argv[], ...));
CCEXTERN_FUNCTION (int DEB_fdebug,
		 (FILE *stream, char *format, ...));
CCEXTERN_FUNCTION (int DEB_vfdebug,
		 (FILE *stream, char *format, va_list args));
CCEXTERN_FUNCTION (int DEB_fprint_file_line,
                 (FILE *stream));

#undef DEB_fdebug
#undef DEB_vfdebug

/* -------------------------------------------------------------------------
   The following is private, my eyes only !
   ------------------------------------------------------------------------- */
  typedef enum { FDEBUG_SET, FDEBUG_GET, FDEBUG_UNSET } FDEBUG_CMD_TYPE;

CCEXTERN_FUNCTION (int DEB_fdebug_parameters,
		   (FDEBUG_CMD_TYPE cmd, ...));

# ifdef DEBUG
#   define DEB_fdebug DEB_fdebug_parameters(FDEBUG_SET, __FILE__, __LINE__, FDEBUG_SET); DEB_the_fdebug
#   define DEB_vfdebug DEB_fdebug_parameters(FDEBUG_SET, __FILE__, __LINE__, FDEBUG_SET); DEB_the_vfdebug
# else
#   define DEB_fdebug DEB_fdebug_parameters(FDEBUG_SET, __FILE__, __LINE__, FDEBUG_GET); DEB_the_fdebug
#   define DEB_vfdebug DEB_fdebug_parameters(FDEBUG_SET, __FILE__, __LINE__, FDEBUG_GET); DEB_the_vfdebug
# endif /* DEBUG */

#ifdef _debug_c

CCSTATIC_FUNCTION (int debug_mode,
                   (int cmd));

#endif _debug_c

#endif _debug_h

/* EOF CCdebug.h */

