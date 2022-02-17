/* lsconfig.h */

/*===================================================================*/
/*  lsconfig.h contains macro definitions which define the code      */
/*  configuration for the lsgrg system build and which define        */
/*  a number of hardcoded limits and lengths                         */
/*                                                                   */
/*  ------------------------------ Change Log  --------------------  */
/*  ** 2/01 jcp ** add definition for LEN_RETRY_HISTORY              */
/*  ** 4/02 jcp ** add definition for LSGRG_TIME_LIMIT               */
/*===================================================================*/
/*---------------------------------------------------------------------*/
/*   OQMS_GAMS specifies that oqms is running under gams               */
/*---------------------------------------------------------------------*/
/*  #define OQMS_GAMS  */
/*---------------------------------------------------------------------*/
/*  defining USER_TERMINATION_ENABLED causes the typedef's for         */
/*  P_GCOMP and P_PARSH to declare functions returning long ints       */
/*  ( the defaults are 'void') so that if either function returns      */
/*  a zero value, the problem run is terminated with a return status   */
/*  of _LSGRG_USER_TERMINATION                                         */
/*---------------------------------------------------------------------*/
/* #define USER_TERMINATION_ENABLED */
/*---------------------------------------------------------------------*/
/*  include definitions, interface code, and library for */
/*  lsgrgc/optquest global optimization system           */

/* #define  LSGRG_OPTQUEST */
/*---------------------------------------------------------------------*/
/*  if we are using the gams interface, define LSGRG_GAMS to activate  */
/*  relevant code                                                      */
/*---------------------------------------------------------------------*/
/*  #define  LSGRG_GAMS */
/*---------------------------------------------------------------------*/
/* include calls to  lsgrg_report() */

/*  #define USER_REPORT  */
/*---------------------------------------------------------------------*/
/* enable calls to clock and timestats */

#define TIMING_ENABLED
/*--------------------------------------------------------------------*/
/* enable sprintf writes to lsgrg_msgbuffer   */

/* #define IO_ENABLED */
//#define IO_ENABLED

/*----------------------------------------------------------------------*/
/* enable fprintf writes from lsgrg_msgbuffer */
/*  and directly from echo routines           */
/*  FILE_IO_ENABLED implies IO_ENABLED        */

/* #define FILE_IO_ENABLED */
//#define FILE_IO_ENABLED

/*----------------------------------------------------------------------*/
/*  DEBUG_GLOBALS controls whether echo_globals will dump all  */
/*  global variables when called                               */

/*  #define DEBUG_GLOBALS */

/*----------------------------------------------------------------------*/
/* flush file buffers on each write to help debugging */

/* #define FLUSH_ALWAYS */

/*---------------------------------------------------------------------*/
/*  defining LSGRG_TESTHARNESS includes code for interface to          */
/*  development testharness                                            */
/*---------------------------------------------------------------------*/

/* #define LSGRG_TESTHARNESS */

/*=========================================================================*/
/*    ------------  size and length parameters ---------                   */
/*                                                                         */
/*  MAXOPTIONS         =  size of table of user settable options           */
/*  MAXALLOC           =  size of allocation table                         */
/*  TITLELENGTH        =  length of problem title string                   */
/*  MSGBUFFER_LENGTH   =  size of memory buffer for formatted writes       */
/*  MAXVARS            =  limit of nbr of vars                             */
/*  MAXROWS            =  limit on nbr of rows                             */
/*  NAMELENGTH         =  maximum allowable length of var/row names        */
/*  MIN_LGRAD          =  minimum allowable length of sparse jacobian      */
/*  J_GROWTH_ALLOWANCE =  allowance for extra jacobian nonzeros            */
/*  MINLBINV           =  minimum allocation for basis inverse             */
/*  BINVFACTOR         =  rough estimate of basis inverse fillin density   */
/*  PLINFY             =  default value for +infinity                      */
/*  PLZERO             =  smallest nonzero value                           */
/*  LEN_RETRY_HISTORY  =  length of retry_history vector.
/*-------------------------------------------------------------------------*/
#define MAXOPTIONS                 50
#define MAXALLOC                  100
#define TITLELENGTH                75
#define MSGBUFFER_LENGTH        1024l
#define MAXVARS                10000l
#define MAXROWS                10000l
#define NAMELENGTH                 8l
#define MIN_LGRAD               3000l
#define J_GROWTH_ALLOWANCE       500l
#define MINLBINV                2000l
#define BINVFACTOR               0.05
#define PLINFY                  1.e30
#define PLZERO                 1.e-30
#define LEN_RETRY_HISTORY          15
#define LSGRG_TIME_LIMIT         1000
#define LSGRG_ITERATION_LIMIT   10000

