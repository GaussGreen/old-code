/* -------------------------------------------------------------------------

   File: CCmessage.h
   Path: /home/julia/projets/dev/Common/libmessage/SCCS/s.CCmessage.h
   Description: base de messages
   Created: 96/10/13
   Author: Jacques WERNERT
   Modified: 99/07/15 17:13:34
   Last maintained by: Jacques WERNERT
   Revision: 1.66

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#ifndef _CCmessage_h
#define _CCmessage_h

#include <CCcommon.h>
SCCS_ID (CCmessage_h_SccsId, "@(#)CCmessage.h	1.66, modified 99/07/15");

#include <stdio.h>
#include <sys/types.h>
#include <time.h>

#include <CCdebug.h>
#include <CCdate.h>

extern int MSG_verbose;
extern int MSG_debug;
extern int MSG_stderr_quiet;
extern int MSG_carriage_return;

#ifdef DEBUG
# define MSG_printf_message DEB_fdebug_parameters(FDEBUG_SET, __FILE__, __LINE__, FDEBUG_SET); MSG_the_printf_message
# define MSG_vprintf_message DEB_fdebug_parameters(FDEBUG_SET, __FILE__, __LINE__, FDEBUG_SET); MSG_the_vprintf_message
#elif defined MSG_FORMAT_TEST
#  define MSG_printf_message  fprintf
#  define MSG_vprintf_message vfprintf
#else
#  define MSG_printf_message DEB_fdebug_parameters(FDEBUG_SET, __FILE__, __LINE__, FDEBUG_GET); MSG_the_printf_message
#  define MSG_vprintf_message DEB_fdebug_parameters(FDEBUG_SET, __FILE__, __LINE__, FDEBUG_GET); MSG_the_vprintf_message
# endif

typedef enum { MSG_DISABLED = 0, MSG_ENABLED, MSG_MODE_MAX_TYPE } MSG_MODES;

#define MSG_enable_verbose()	MSG_verbose = MSG_ENABLED;
#define MSG_disable_verbose()	MSG_verbose = MSG_DISABLED;
#define MSG_enable_debug()	MSG_debug = MSG_ENABLED;
#define MSG_disable_debug()	MSG_debug = MSG_DISABLED;
#define MSG_set_stderr_quiet()	MSG_stderr_quiet = MSG_ENABLED;
#define MSG_unset_stderr_quiet() MSG_stderr_quiet = MSG_DISABLED;
#define MSG_enable_carriage_return() MSG_carriage_return = MSG_ENABLED;
#define MSG_disable_carriage_return() MSG_carriage_return = MSG_DISABLED;

/*
 * Base de messages. Format:
 *
 * format: MSG_<type>_<module>_<message>
 * 
 * ou type est {F|I|E|T}: Fatal, Informatif, Erreur, Trace
 * et module est {UNX|SYB|DB|TIB|MIS} : UNiX, SYBase, DataBase, TIB, 
 */

/*
 * codes des niveaux d'erreurs
 */
#define MSG_TRACE           40      /* debug/verbose */
#define MSG_INFO            30      /* infos console */
#define MSG_ERROR           20      /* erreurs non fatales */
#define MSG_ERROR_DETAIL    25      /* erreurs non fatales */
#define MSG_FATAL           10      /* erreurs fatales */

/* 
 * libelle des niveaux d'erreurs 
 */
#define MSG_TRACE_STR         "Trace"
#define MSG_INFO_STR          "Info"
#define MSG_ERROR_STR         "Erreur"
#define MSG_FATAL_STR         "Fatal"


/*
 * messages d'erreurs systemes (UNX)
 */
#define MSG_F_UNX_MEMORY                    "memoire insuffisante: %s"
#define MSG_F_UNX_PARAMETER                 "retour invalide"
#define MSG_F_UNX_ABORT                     "Program stoppe"
#define MSG_F_UNX_FORK_FAILED               "appel systeme fork impossible: %s"
#define MSG_F_UNX_DUP2_FAILED               "appel systeme dup2 impossible: %s"
#define MSG_F_UNX_EXEC_FAILED               "appel systeme exec impossible: %s"
#define MSG_F_UNX_EXEC_NAMED_FAILED         "appel systeme exec sur le fichier %s impossible: %s"
#define MSG_F_UNX_PIPE_FAILED               "appel systeme pipe impossible: %s"
#define MSG_F_UNX_WAIT_FAILED               "appel systeme wait impossible: %s"
#define MSG_F_UNX_STAT_FAILED               "appel systeme stat sur l'entree %s impossible: %s"
#define MSG_F_UNX_LOCKF_FAILED              "appel systeme lockf impossible sur le fichier '%s': %s"
#define MSG_F_UNX_LOCKF_ULOCK_FAILED        "appel systeme lockf, impossible de deposer un verrou sur '%s': %s"
#define MSG_F_UNX_LOCKF_TEST_FAILED         "appel systeme lockf, impossible de recuperer les informations de verrou de '%s': %s"
#define MSG_F_UNX_FCNTL_FAILED              "appel systeme fcntl impossible: %s"
#define MSG_F_UNX_FCNTL_GETLK_FAILED        "appel systeme fcntl, impossible de recuperer les informations de verrou de '%s': %s"
#define MSG_F_UNX_FCNTL_SETLK_FAILED        "appel systeme fcntl, impossible de (de)poser un verrou sur '%s': %s"
#define MSG_E_UNX_ALREADY_LOCKED            "il existe deja un verrou sur le fichier '%s'"
#define MSG_E_UNX_GETPWUID_FAILED           "appel systeme getpwuid, impossible de recuperer le repertoire 'home' de l'utilisateur"
#define MSG_E_UNX_GETWD_FAILED              "appel systeme getwd, impossible de recuperer le repertoire courant" 
#define MSG_E_UNX_SELECT_ERROR              "erreur sur appel systeme select: %s"

#define MSG_E_UNX_DATA_ON_UNMANAGED_FD      "select indique des donnees sur un descripteur non gere"


/*
 *
 */
#define MSG_I_APP_INIT_START               "initialisation en cours ..."
#define MSG_I_APP_INIT_OK                  "initialisation terminee"
#define MSG_F_APP_INIT_FAILED              "initialisation avortee"
#define MSG_I_APP_END                      "fin de l'application"
#define MSG_I_APP_END_WITH_ERROR           "fin avec erreur de l'application"
#define MSG_I_APP_CHILD_NORMAL_END         "fin volontaire du fils %d avec code de retour %d"
#define MSG_I_APP_CHILD_SIGNAL_END         "fin par signal %d \"%s\" du fils %d%s"
#define MSG_I_APP_END_WTH_CORE(status)     (WAIT_CORE_DUMP (status) ? " (core genere)" : "")
#define MSG_I_APP_SIGNAL_TRAP              "reception du signal %d \"%s\""

#define SIGINT_STR          "interruption (ctrl-C)"
#define SIGQUIT_STR         "quitter (ctrl-\\)"
#define SIGFPE_STR          "erreur mathematique (floating point exception)"
#define SIGKILL_STR         "fin irremediable"
#define SIGBUS_STR          "adresse incorrecte (bus error)" 
#define SIGSEGV_STR         "adresse (segmentation violation)"
#define SIGSYS_STR          "(bad argument to system call)"
#define SIGPIPE_STR         "ecriture pipe sans lecteur"
#define SIGALRM_STR         "delai a expire (alarm)"
#define SIGTERM_STR         "arret standard (software termination)"
#define SIGCHLD_STR        "mort d'un fils"
#define SIGSTOP_STR         "stopper"
#define SIGCONT_STR         "continuer"
#define SIGTTIN_STR         "tentative de lecture sur terminal en tache de fond"
#define SIGTTOU_STR         "tentative d'ecriture sur terminal en table de fond"
#define SIGUNKNOWN_STR      "pas de description"

/*
 * messages d'erreurs sybase (SYB) et base de donnees (DB)
 */
#define MSG_I_SYB_RECONNECTION              "re-connexion au serveur Sybase en cours"
#define MSG_E_SYB_ERROR                     "erreur Serveur-SQL: %d, etat %d, gravite %d: %s"
#define MSG_E_SYB_CONNECTION_LOST           "connexion au SQL-serveur perdue"
#define MSG_E_SYB_LEAVING_CONNECTION        "fermeture de la connexion au SQL-serveur"
#define MSG_F_SYB_SQL_FAILED                "l'execution de la requete SQL a echoue"
#define MSG_F_SYB_DBNEXTROW_FAILED          "erreur a la recuperation de l'enregistrement suivant"
#define MSG_F_SYB_LIBRARY_ERROR             "erreur DB-lib: %d, gravite %d: %s"
#define MSG_F_SYB_INIT_FAILED               "erreur a l'initialisation Sybase"
#define MSG_F_DB_CNX_PB                     "probleme de connexion au serveur %s/%s/%s"
#define MSG_I_DB_CNX_OK                     "connexion au serveur %s ok"
#define MSG_F_DB_USE_PB                     "probleme de connexion a la base %s"
#define MSG_F_DB_FATAL_ERROR                "erreur fatale SQL-serveur"
#define MSG_F_DB_SYS_ERROR                  "erreur systeme SQL-serveur: %s"
#define MSG_F_SYB_DBRESULTS_FAILED          "Dbresults failed"
#define MSG_E_SYB_DBDBG                     "dbdbg trace: %s"
#define MSG_E_SYB_NOT_ENOUGH_RESOURCES      "probleme de ressources systeme"
#define MSG_F_SYB_BAD_USER_PROGRAM          "probleme applicatif"
#define MSG_E_SYB_SQL_SERVER_ERROR_SRVNAME  "Serveur %s,"
#define MSG_E_SYB_SQL_SERVER_ERROR_PROCNAME "Procedure %s,"
#define MSG_E_SYB_SQL_SERVER_ERROR_LINE     "Ligne %d,"
#define MSG_E_SYB_SQL_SERVER_USER_ERROR     "Erreur utilisateur en provenance du SQL-Server"
#define MSG_E_SYB_SQL_SERVER_SYSTEM_ERROR   "Erreur systeme en provenance du SQL-Server. Contacter l'Administrateur Systeme"
#define MSG_E_SYB_SQL_PROC_BAD_RETURN	    "Valeur de retour de la procedure incorrecte (%d)"

/*
 * messages d'erreurs des fonctions Tib (TIB) (y compris l'API signal, 
 * io et timer)
 */
#define MSG_I_TIB_UNSUBSCRIBE_OK            "arret de l'abonnement Tib ok pour le marche %s"
#define MSG_I_TIB_SUBSCRIBE_OK              "abonnement Tib reussi pour le sujet %s"
#define MSG_E_TIB_UNKNOWN_CLASS             "classe inconnue: %s"
#define MSG_E_TIB_SUBSCRIBE_FAILED          "erreur d'abonnement Tib pour le sujet %s: %s"
#define MSG_E_TIB_GET_FIELD_ID_FAILED       "erreur de recuperation du field id du champ %s de la classe %s: %s"
#define MSG_E_TIB_INCONSISTENT_DISCIPLINE   "type d'abonnement inconnu: %s"
#define MSG_I_TIB_NO_CI_SERVER              "ne peut se connecter au ciServer"

#define MSG_F_TIB_FORM                      "erreur au chargement des formes"
#define MSG_F_TIB_FORMCLASS_INIT_FAILED     "erreur a l'initialisation des \"FormClass\" TIB"
#define MSG_F_TIB_SUBJECT_INIT_FAILED       "erreur a l'initialisation de la liste des abonnements TIB"
#define MSG_F_TIB_PUBLISH_INIT_FAILED       "erreur a l'initialisation de la publication TIB : %s"
#define MSG_F_TIB_CREATE_CLASS_FROM_FILE_FAILED  "erreur a l'initialisation fichier de la classe %s"
#define MSG_E_TIB_FIELD_SET_DATA_FAILED     "erreur a l'initialisation de la valeur: %s"
#define MSG_E_TIB_FORM_CREATE_FAILED        "erreur a la creation de la form: %s"
#define MSG_E_TIB_PUBLISH_FAILED            "erreur a la publication : %s"
#define MSG_F_TIB_SUBSCRIBE_FAILED          "erreur a l'initialisation des abonnements TIB"
#define MSG_I_TIB_RE_SUBSCRIBE              "re-initialisation des abonnements TIB"
#define MSG_I_TIB_RE_SUBSCRIBE_END          "fin de la re-initialisation des abonnements TIB"
#define MSG_F_TIB_INIT_ADDIO                "erreur d'initialisation des callback sur flux"
#define MSG_F_TIB_INIT_ADDSIG               "erreur a l'initialisation des \"handler\" de signaux"
#define MSG_F_TIB_ADDIO_FAILED              "erreur a l'ajout d'un lecteur sur flux (AddIo)"
#define MSG_F_TIB_REMOVEIO_FAILED           "erreur a la suppression d'un lecteur sur flux (RemoveIo)"
#define MSG_F_TIB_ADDTIMER_FAILED           "erreur a l'ajout d'un 'timer' (AddTimer)"
#define MSG_F_TIB_REMOVETIMER_FAILED        "erreur a la suppression d'un 'timer' (RemoveTimer)"
#define MSG_F_TIB_MAIN_LOOP_ERROR           "erreur de la \"main loop\""
#define MSG_E_TIB_FORMAT_NULL_MSG           "reception d'un message Tib NULL sous le sujet %s"
#define MSG_E_TIB_UNKNOWN_FORMAT_MSG        "reception d'un message Tib de format inconnu sous le sujet %s (%d)"

#define MSG_E_RPC_UDP_SERVICE_CREATION_FAILED   "impossible de creer un service UDP pour le serveur RPC 0x%x version %d"
#define MSG_E_RPC_TCP_SERVICE_CREATION_FAILED   "impossible de creer un service TCP pour le serveur RPC 0x%x version %d"
#define MSG_E_RPC_UDP_CANT_REGISTER             "impossible d'enregistrer le serveur 0x%x version %d aupres du portmapper pour le service UDP"
#define MSG_E_RPC_TCP_CANT_REGISTER             "impossible d'enregistrer le serveur 0x%x version %d aupres du portmapper pour le service TCP"
#define MSG_E_RPC_NETPATH_CREATION_FAILED       "impossible de creer le serveur 0x%x version %d aupres du 'netpath'"

/*
 * misc: les X-files des messages d'erreurs, les messages generaux
 *       et les messages d'assert
 */
#define MSG_I_MIS_INIT_OK                    "initialisation terminee"
#define MSG_F_MIS_UNKNOWN_ERROR              "erreur inconnue"
#define MSG_F_MIS_INCONSISTENT_VALUE         "valeur inconsistante"
#define MSG_F_MIS_INCONSISTENT_INT_VALUE     "valeur %d inconsistante"
#define MSG_F_MIS_ABW                        "ecriture en dehors des bornes du tableau (array bounds write)"
#define MSG_F_MIS_INIT_DATA_FILE_FAILED      "erreur a l'initialisation des fichiers de donnees"
#define MSG_F_MIS_INCONSISTENT_CONSTANT      "probleme de definition de constante"
#define MSG_F_MIS_INCONSISTENT_RETURN        "valeur retournee invalide"
#define MSG_I_MIS_WELCOME_TO_WITH_REV	     "** %s version %s **"
#define MSG_I_MIS_COMMAND_LINE_HELP          "Options de la ligne de commande:"
/*
 * messages exportes dans MSG_comm.h
 *
 * #define MSG_F_MIS_CONTACT_SUPPORT            
 * #define MSG_F_MIS_NO_DESCRIPTION
 */

#define MSG_get_entete() MSG_set_header (NULL)

// FIXMEFRED: mig.vc8 (25/05/2007 10:38:53):missing return type
CCEXTERN_FUNCTION (int MSG_the_printf_message,
		 (int message_type, const char *format, ...));
CCEXTERN_FUNCTION (int MSG_the_vprintf_message,
		 (int message_type, const char *format, va_list args));
CCEXTERN_FUNCTION (char *MSG_set_header,
		 (char *entete));
#define MSG_message_set_stdlog MSG_set_stdlog  /* backward compatibility */
CCEXTERN_FUNCTION (void MSG_set_stdlog, 
                 (FILE *new_stdlog));
CCEXTERN_FUNCTION (char *MSG_ctime, (time_t timestamp));
CCEXTERN_FUNCTION (int MSG_print_eol, (FILE *flux));

#ifdef _CCmessage_c

#define MSG_must_be_quiet_on_stderr(file) ((MSG_stderr_quiet) && (file == stderr))

#endif _CCmessage_c

#endif /* _CCmessage_h */

/* EOF CCmessage.h */
