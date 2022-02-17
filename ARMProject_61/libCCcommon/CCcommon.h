/* -------------------------------------------------------------------------

   File: CCcommon.h
   Path: /home/julia/projets/dev/Common/libcommon/SCCS/s.CCcommon.h
   Description: Header de fichier commun
   Created: 96/12/10
   Author: Jacques WERNERT
   Modified: 99/08/24 12:03:46
   Last maintained by: Jacques WERNERT
   Revision: 1.40

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#ifndef _common_h
#define _common_h

#pragma warning(disable : 4273)

#include <sys/types.h>
#include <assert.h>

/*
 * macro speciale anti "-Wall"
 *    elle permet, lorsque le compilateur GNU CC est utilise avec l'option
 *    -Wall (Warning all) d'eviter d'avoir un message du type:
 *    "warning: unused variable `pouetpouet'". Il suffit de placer un
 *    USE(pouetpouet) pour eviter ce message apres la declaration.
 *    Elle peut etre utilisee en particulier pour les SccsId ou RcsId
 */
#if defined(__GNUC__) && defined (__STDC__)
#  define CC_USE(var) static void *use_##var = (&use_##var, (void *) &var)
#endif

/*
 * macro de definition des SccsID
 */
#if ! (defined(NO_WHATLIST) || defined(lint))
#  if defined(__GNUC__) && defined(__STDC__)
#    define SCCS_ID(file, value) static char file[] = value; CC_USE(file)
#  else
#    define SCCS_ID(file, value) static char file[] = value 
#  endif
#else
#  define SCCS_ID(file, value)  
#endif

SCCS_ID (comm_h_SccsID, "@(#)CCcommon.h	1.40, modified 99/08/24");

#ifdef __STDC__
#  include <stdarg.h>
#else
#  include <varargs.h>
#endif

#include <CCfdecl.h>
#ifndef WIN32
#include <CCmotifcst.h>
#endif

/*
 * macro d'init des fonctions a argument variable
 */
#ifdef __STDC__
#  define CC_VA_START(args, last_parameter) va_start(args, last_parameter)
#else
#  define CC_VA_START(args, last_var) va_start(args)
#endif __STDC__

#ifdef CC_COMPAT
#  define VA_START(args, last_parameter)  CC_VA_START ((args), (last_parameter))
#endif /* CC_COMPAT */

/*
 * messages d'erreurs extraits de MSG_message.h
 */
#define COM_F_MI_CONTACT_SUPPORT "Erreur Fatale: %s\n\
Contactez votre support technique avec un descriptif des actions et \
les lignes affichees\n"
#define COM_F_MI_NO_DESCRIPTION  "(pas de description)"
#define COM_F_MI_ABW             "array bounds write"
/*
 * fonction de test bas-niveau
 * si la condition est fausse, stoppe le process en affichant "msg", le nom de 
 * fichier et la ligne ou l'evenement s'est produit
 */
#define CC_ASSERT(cond, msg) if (!(cond)) { \
                             fprintf (stderr, COM_F_MI_CONTACT_SUPPORT, \
				      msg?msg:COM_F_MI_NO_DESCRIPTION);\
                             assert ((cond)); \
                          }

#ifdef CC_COMPAT
#define ASSERT(cond, msg) CC_ASSERT ((cond), (msg))
#endif /* CC_COMPAT */

/*
 * fonction de copie de chaine 'securisee' a base de strncpy
 */
#ifdef DEBUG
#  define CC_STRCPY(dest, src, len) strncpy(dest, src, len); \
                                CC_ASSERT(dest[len-1] == '\0', COM_F_MI_ABW)
#  define CC_STRCAT(dest, src, len) strncat(dest, src, len); \
                                CC_ASSERT(dest[len-1] == '\0', COM_F_MI_ABW)
#else
#  define CC_STRCPY(dest, src, len) strncpy(dest, src, len); \
                                 dest[len-1] = '\0'
#  define CC_STRCAT(dest, src, len) strncat(dest, src, len - strlen (dest) - 1); \
                                 dest[len-1] = '\0'
#endif

#ifdef CC_COMPAT
#  define STRCPY(dest, src, len) CC_STRCPY ((dest), (src), (len))
#  define STRCAT(dest, src, len) CC_STRCAT ((dest), (src), (len))
#endif

#define CC_SAFE_STR(str) ((str) ? (str) : "(null)")

/*
 * acces au messages d'erreurs
 */

#include <errno.h>

#if !(defined (WIN32) && defined(_MT))
extern int errno;             /* pas tjrs exporte dans errno.h */
#endif

// FIXMEFRED: mig.vc8 (25/05/2007 10:42:59): sys_errlist and sys_nerr are defined in <stdlib.h>, and not in <errno.h>
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	#include <stdlib.h>
#else
	extern char *sys_errlist[];   /* pas tjrs exporte dans errno.h */
#endif

#define ERRNO_GET_MESSAGE     (sys_errlist[errno])

/*
 * Macro pour dimensionner une chaine destinee a recevoir un numerique
 */
#define CHARS_TO_PRINT_BYTES(nb_of_bytes) (3 * nb_of_bytes)
#define CHARS_TO_PRINT_SIGNED_BYTES(nb_of_bytes) (3 * nb_of_bytes + 1)

/*
 * gestion des codes fournis par l'appel systeme wait
 *                  poids Fort,   Faible
 * bien termine: 	[ exit ] [ 0 ]
 * terminaison syst: 	[ 0 ] 	 [ (core ?) signal ]
 */ 
#define NO_PROCESS_NUMBER               (-1)
#define WAIT_NORMAL_RET(status)		(WAIT_ABNORMAL_RET (status) == 0)
#define WAIT_ABNORMAL_RET(status)	((status) & 0377)
#define WAIT_STATUS_CODE(status)	((status) >> 8)
#define WAIT_SIGNAL_CODE(status)	((status) & 0177)
#define WAIT_CORE_DUMP(status)		(((status) & 0200) ? 1 : 0)

#define DEFAULT_NULL_STRING             ""

/*
 * Certains systemes ne definissent meme pas les "protos" des fonctions de la libc et 
 *  des appels systemes.
 */
#ifdef WIN32
# ifdef __STDC__
# define strdup _strdup
# else
  char *strdup ();
# endif
#endif

/*
 * Definitions de constantes non definies sous Win32
 */

#ifdef WIN32
#define MAXNAMELEN	256

#define rint(x)	(((double)(x) - (int)(x) <= 0.5) ? floor ((x)) : ceil ((x)))

#endif	/* WIN32 */

#define CC_MIN(a, b)  ((a) < (b) ? (a) : (b))
#define CC_MAX(a, b)  ((a) > (b) ? (a) : (b))

#ifdef NEED_SYS_PROTO

  int wait ();
  int gethostname ();
  int fprintf ();
  int vfprintf ();
  int fflush ();
  int sscanf ();
  int fclose();
  char *getwd ();
  int lockf ();
  void setbuf ();
  int setlinebuf ();
  time_t time ();
  int puts ();
  int printf ();
  int semop ();
  int semctl ();
  int strncmp ();
  int fprintf ();
  key_t ftok ();
  int shmctl ();
  int shmget ();
  int semget ();
  char *shmat ();
  int shmdt ();
  void perror ();
  int strftime ();
  char *re_comp ();
  int re_exec ();
  int select ();
  void bzero ();
  int tolower();
#endif

#endif _common_h

/* EOF CCcommon.h */
