#include <StdAfx.h>
#include "lsinfo.h"
#define IOINFO  &(_info->io_info)
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)

/*  ** 2/01 jcp **                                  */
/*  enumeration constant for set of possible fixups */
/*  this goes in grgsub.c or lsinfo.h               */
enum fixup_codes { RESTART = 1, SCALE = 2, PERTURB_POINT = 3,
                   RESIZE_BINV=4, RESIZE_JAC=5,EXIT = 6};

/***********************************************************************
 **       FILE:              GRGSUB.C                                  *
 **       AUTHOR:            John Plummer   March 1998                 *
 **       LAST UPDATE:                                                 *
 **       LAST CHANGE:       Feb 2001                                  *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **                                                                    *
 **  This is a complete rewrite of the lsgrgc callable interface       *
 **                                                                    *
 **  ROUTINE             DESCRIPTION                      LAST UPDATE  *
 ** ---------            -------------                   ------------- *
 **  grgsub()            lsgrg callable interface            jul 1998  *
 **                                                                    *
 **  lsgrg_setvarname()  sets var names                      oct 1998  *
 **  lsgrg_setrowname()  sets row names                      oct 1998  *
 **                                                                    *
 **                                                                    *
 **  lsgrg_errorexit() called to abort lsgrg                           *
 **                    execution on fatal errors                       *
 **                    during solution process                         *
 **                    calls longjmp() to                              *
 **                    grgsub to return to caller                      *
 **                                                                    *
 ** lsgrg_initialize_grgtest()                                         *
 **                    initializes test harness fields                 *
 **                                                                    *
 ** lsgrg_get_terminationmsg()                                         *
 **                    returns pointer to test harness termination     *
 **                    string in case user wants it                    *
 **                                                                    *
 **   lsgrg_xfer_alg_term_code()                                       *
 **                    maps algorithm termination code into            *
 **                    more detailed set defined for shell             *
 **                                                                    *
 **  void lsgrg_set_termin_msg()                                       *
 **                    copies alg termination msg into                 *
 **                    test harness termination string for             *
 **                    query by user                                   *
 **                                                                    *
 **  void lsgrg_set_termin_code()                                      *
 **                    copies alg termination string into test         *
 **                    harness s_termin for use in test harness        *
 **                    statistics.                                     *
 ***********************************************************************
 * */
/***********************************************************************
 ***********************************************************************
 ***                         LSGRG2                                  ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
*/
int LSGRGCLASS grgsub (LsgrgInfo *user_info,
            long nvars_in, long nfun_in, long nobj_in,long maximize,
            long lvars[],
            double blvar[],double buvar[],double blcon[],double bucon[],
            double xx[],double fcns[],double rmults[],long nonbas[],
            double redgr[],long inbind[],long *nbind,long *nnonb,
            P_GCOMP p_user_gcomp,P_PARSH p_user_parsh, long nnz_user)
{
/*         ***************************************************
 *         ***************************************************
 *         ***                 LSGRG2                      ***
 *         ***     COPYRIGHT  1991, 1992, 1993, 1994, 1998 ***
 *         ***              LEON S. LASDON                 ***
 *         ***                    &                        ***
 *         ***              ALLAN D. WAREN                 ***
 *         ***                    &                        ***
 *         ***               STUART SMITH                  ***
 *         ***                    &                        ***
 *         ***              John C. Plummer                ***
 *         ***************************************************
 *         ***************************************************
 *
 *     THE ARGUMENTS ARE DIVIDED INTO 5 SECTIONS:
 *
 *     1. Required Scalar Inputs:
 *         . info                  pointer to LsgrgInfo structure
 *                                   returned by call to LsgrgInitialize
 *         . nvars_in              nbr of problem variables
 *         . nfun_in               nbr of problem functions(including obj)
 *         . nobj_in               index of objective function
 *         . maximize              logical for sense of optimization
 *                        0 = minimize
 *                        1 = maximize
 *                           transfered to global maxim
 *         . nnz_user              number of nonzero jacobian elements
 *                                 if user is supplying derivatives
 *                                 (set to 0 otherwise)
 *     2. Required Array Inputs:
 *        NOTE:  all array inputs and outputs assume that values start
 *               at array element 1 (except for lvars as below)
 *
 *         . lvars                 set of purely linear variables
 *                                 lvars[0] = nbr of linear variables(required)
 *                                 if(lvars[0]>0) then
 *                                  lvars[k] = index of kth linear var
 *
 *         . blvar                 variable lower bounds
 *         . buvar                 variable upper bounds
 *         . blcon                 function lower bounds
 *         . bucon                 function upper bounds
 *         . xx                    variable initial values
 *                                 (final values on output)
 *     3. lsgrg outputs:
 *         . fcns[]                final function values
 *         . rmults[]              final multipliers
 *         . redgr[]               final reduced gradients
 *         . nonbas[]              indices of nonbasic vars
 *         . inbind[]              indices of binding constraints
 *         . xx[]                  final variable values
 *         . nbind                 nbr of binding constraints
 *         . nnonb                 nbr of nonbasic variables
 *
 *     TERMINATION CODES (returned as grgsub result)
 *     -----------------------------------------------
 *
 *      The lsgrgc 3.0 termination codes are parametrized with an
 *        enum constant defined in lscodes.h
 *
 *  Change Log
 *  ----------
 *
 *
 *  11/00 jcp fixed bug in loop which loads up the reduced gradients
 *   for return to the user. loop was 1=0;i<=_info.n, should have been
 *   i=1;i<=_info.n
 *
 *   2/01 jcp  add retry loop to grgitn call with call to
 *     lsgrg_fixup_strategy to set fixup strategy
 *     add retry_history vector
 *
 *
 *   4/01    jcp  complete setup/run/shutdown logic and add
 *                code to copy point and bounds on run only
 *                . add user functions to specify setup/run/shutdown
 *
 *                .  grgsub needs to copy its local LsgrgInfo structure
 *                   back to the users structure at exit of setup phase if
 *                   doing setup/run/shutdown logic to save allocation
 *                   pointers and values set during setup
 *                .  add additional pointers for row/col to LsgrgInfo
 *                   so that we can store locally allocated addresses
 *                   between calls
 *
 *
 *  4/27/02 jcp  remove local handling of row and var names
 *   to fix memory leak problems.  also ensure that _LSGRG_ALLOCATION_FAILURE
 *   is returned for any memory allocation failure
 *
 *  NOTE:  the copy of the user's structure to the local structure
 *  will copy the pointers allocated by set_rownames and set_varnames
 *
 *  ** 7/26/03 jcp ** if this is a fullrun, reset maxr to -1
 *  so that if the user uses the same LsgrgInfo structure for
 *  multiple runs, maxr will be reset to the actual nbr of vars
 *  #### set the value in the user structure AND the local
 *   just in case we are doing a copyback
 *  **** doing this is probably dangerous (for a user) ****

  **********************************************************************/
/*---------------------------------------------------------------------*/
/* NOTE: lsgrg_set_defaults is called by LsgrgInitialize, so debug     */
/*  flags (and other values set there) are available here at entry     */
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
/*         local declarations                                          */
/*  create local LsgrgInfo structure                                   */
/*---------------------------------------------------------------------*/
   static nruns=0; /* **fixme** take out after testing */
    LsgrgInfo _info;
    int grgitn_return;
    long i,j,nindex,ninf;
    double tstart,temp;
    long ibmap[10];  /* dummy arrays for v1 grgitn call */
    double zdummy[10];
    int retry_count,max_retry_count=3,fixup_code;
    char error_classification[100];
/*    char ** user_varnames, **user_rownames; */
/*    char *p1,*p2; */
    int debug = user_info->dbg.grgsub;
/*---------------------------------------------------------------------*/
/*   each entry to grgsub gets a new jmp_buf for error terminations    */
/*---------------------------------------------------------------------*/
    jmp_buf lsexit_jmpbuf;
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*  if all three run mode flags are 0, do a simple setup/run/shutdown  */
/*  setup_done will track whether a setup has been done before a       */
/*  run or shutdown is called                                          */
/*---------------------------------------------------------------------*/
    int lsgrg_fullrun;
/*---------------------------------------------------------------------*/
/*  .copy the user's LsgrgInfo structure to our local copy so that     */
/*   all access past this point is to the local copy                   */
/*  . set msgbuffer pointer in io_info substructure to point to        */
/*    allocated msgbuffer                                              */
/*  . set pointer in LsgrgInfo struct to point to local jmp_buf        */
/*---------------------------------------------------------------------*/
     _info = *user_info;
     _info.io_info.lsgrg_msgbuffer      = _info.lsgrg_msgbuffer;
     _info.lsexit_jmpbuf = &lsexit_jmpbuf;
        lsgrg_initialize_grgtest(&_info);  /* init test harness fields */

/*---------------------------------------------------------------------*/
/*  if all three run mode flags are 0, do a simple setup/run/shutdown  */
/*  setup_done will track whether a setup has been done before a       */
/*  run or shutdown is called                                          */
/*  . check for invalid combination of modes                           */
/*---------------------------------------------------------------------*/
#ifdef IO_ENABLED
    if(debug) {
       sprintf(_info.lsgrg_msgbuffer,
             "\n .... entry grgsub ...."
              "\n .. _info.lsgrg_setup = %d _info.lsgrg_run = %d"
              " _info.lsgrg_shutdown = %d",_info.lsgrg_setup,
              _info.lsgrg_run,_info.lsgrg_shutdown);
        lsgrg_screen_msg(&_info.io_info,_info.lsgrg_msgbuffer);
    }
#endif
        i = _info.lsgrg_setup+_info.lsgrg_run+_info.lsgrg_shutdown;
        if(i !=0 && i !=1) {
#ifdef IO_ENABLED
         sprintf(_info.lsgrg_msgbuffer,
           "\n Error: Inconsistent Setup/Run/Shutdown Flags "
           " lsgrg_setup = %d lsgrg_run = %d lsgrg_shutdown = %d",
             _info.lsgrg_setup,_info.lsgrg_run,_info.lsgrg_shutdown);
          lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
            return _LSGRG_BAD_COMMAND;
        }
/*  check for disabled copyback of info structure  */

        if( _info.lsgrg_setup|| _info.lsgrg_run || _info.lsgrg_shutdown) {
            if(_info.disable_copyback) {
#ifdef IO_ENABLED
         sprintf(_info.lsgrg_msgbuffer,
           "\n Error: grgsub called with Setup/Run/Shutdown requested,"
           " but \n  "
           " LsgrgInfo structure copyback disabled. Unallowable "
           " Configuration");
          lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
            return _LSGRG_BAD_COMMAND;
            }
        }
        lsgrg_fullrun =
        _info.lsgrg_setup==0 && _info.lsgrg_run==0 &&
        _info.lsgrg_shutdown==0;
/*---------------------------------------------------------------------*/
/*  6/01 jcp modify io assignments                                     */
/*  .initialize default output units if file-io is enabled             */
/*                                                                     */
/*     if user has set ioerr, send error msgs to that destination      */
/*     only.  set ioerr_to_ioterm and ioerr_to_ioout to 0              */
/*                                                                     */
/*     if ioout==NUll, default all output to screen including error    */
/*       messages.  if user has not set ioterm, disable screen output  */
/*       to avoid duplication                                          */
/*     if user set ioout, set ioerr_to_ioout to send error msgs to     */
/*        ioout.  if user has set ioterm set ioerr_to_ioterm to        */
/*        send error messages there as well.                           */
/*     if the user has not set ioterm, but has set ioout, set ioterm   */
/*        to stdout and set ioerr_to_ioterm to 1 to send error messages*/
/*        to screen                                                    */
/*     NOTE:                                                           */
/*       the logicals error_output and screen_output override          */
/*       destination assignments and are handled by the io functions   */
/*       in grgio.c                                                    */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*---------------------------------------------------------------------*/
#ifdef FILE_IO_ENABLED
     if(_info.io_info.lsgrg_ioerr != NULL) {
          _info.io_info.ioerr_to_ioterm = 0;
          _info.io_info.ioerr_to_ioout  = 0;
     }
     if(_info.io_info.lsgrg_ioout==NULL)  {
         _info.io_info.lsgrg_ioout = stdout; /* all output to screen */
         _info.io_info.ioerr_to_ioterm = 0;  /* no duplicate msgs */
         _info.io_info.ioerr_to_ioout  = 1;  /* send errmsgs to screen */
     }
     else {
          _info.io_info.ioerr_to_ioout = 1;
     }
     if(_info.io_info.lsgrg_ioterm == NULL)  {
        if(_info.io_info.lsgrg_ioout != stdout)  /* dont send ioout and ioterm */
           _info.io_info.lsgrg_ioterm = stdout; /* to stdout */
        if(_info.io_info.ioout_set_by_user== 1)
           _info.io_info.ioerr_to_ioterm = 1;
     }
     else
         _info.io_info.ioerr_to_ioterm = 1;

#endif
/*---------------------------------------------------------------------*/
          _info.lsgrg_return_status = _LSGRG_STATUS_NOT_SET;

    if(debug) lsgrg_echo_globals(&_info,"Entry grgsub");
/*---------------------------------------------------------------------*/
/*   initialize members which replace statics in alg routines          */
/*   do this on each call so that algorithm is restarted correctly     */
/*---------------------------------------------------------------------*/
    _info.firstcall.newton =  0;
    _info.firstcall.redobj =  0;
    _info.firstcall.ph0log = TRUE;

    _info.cgStatic.initcg = 1;
    _info.cgStatic.itncg  = 0;

    _info.direcStatic.nsupvm = -1;
    _info.direcStatic.msgvm  =  0;
    _info.direcStatic.rtinsb = 4.0e0;

    _info.ph1objStatic.sinf0 = 0.0e0;
    _info.ph1objStatic.true0 = 0.0e0;

    _info.ph0logStatic.termheader = 0;
    _info.ph0logStatic.iprhd3     = 0;
    _info.ph0logStatic.iprhld     = 0;
/*=====================================================================*/
/*      end code block executed on each entry to grgsub                */
/*=====================================================================*/
/*---------------------------------------------------------------------*/
/*   begin lsgrg2c setup code                                          */
/*---------------------------------------------------------------------*/
    if(debug) lsgrg_screen_msg(&_info.io_info,"at if(info.lsgrg_setup)");
    if(_info.lsgrg_setup || lsgrg_fullrun) {
       if(debug) lsgrg_screen_msg(&_info.io_info,"initiating problem setup");

      _info.maxim = maximize; /* transfer max/min flag from user arg */
/*---------------------------------------------------------------------*/
/*   get initial time and initialize system parms                      */
/* 1.  if no user calls to lsgrg_setparameter(), which will cause      */
/*      parm defaults to be set, call lsgrg_set_defaults() here        */
/* **fixme** revise timing results                                     */
/*---------------------------------------------------------------------*/
        lsgrg_timestats( &_info,T_INIT, (double) 0.0);

        _info.getsiz = FALSE;
/*---------------------------------------------------------------------*/
/*   rough check on input problem dimensions                           */
/*---------------------------------------------------------------------*/
        if( nvars_in <  1 || nvars_in > MAXVARS ) {
#ifdef IO_ENABLED
           sprintf( _info.lsgrg_msgbuffer,
            "\n ERROR - nvars_in < 1 or > %d  (%5ld)",
              MAXVARS,nvars_in );
           lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
            _info.lsgrg_return_status = _LSGRG_DIMENSION_ERROR;
        }
        if( nfun_in < 1 || nfun_in > MAXROWS) {
#ifdef IO_ENABLED
           sprintf( _info.lsgrg_msgbuffer,
            "\n ERROR - nfun_in < 1 or > %d  (%5ld)",
              MAXROWS,nfun_in );
           lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
            _info.lsgrg_return_status = _LSGRG_DIMENSION_ERROR;
        }
        if(_info.lsgrg_return_status) return _info.lsgrg_return_status;

        if( nobj_in < 1 || nobj_in > nfun_in ) {
#ifdef IO_ENABLED
           sprintf( _info.lsgrg_msgbuffer,
            "\n ERROR - NNOBJ=(%5ld)  < 1 OR > NFUN (%5ld)",
              nobj_in,nfun_in );
           lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
           _info.lsgrg_return_status =_LSGRG_BAD_NOBJ;
           return _info.lsgrg_return_status;
        }

        if( lvars[0] < 0 || lvars[0] > nvars_in )  {
#ifdef IO_ENABLED
           sprintf( _info.lsgrg_msgbuffer,
            "\n ERROR - NUMBER OF LINEAR VARIABLES (%10ld) \n"
            "           IS < 0 OR > NVARS(%d) ",
            lvars[0], nvars_in);
            lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
            return _LSGRG_LINEAR_VARS_ERROR;
        }

        _info.nvars  = nvars_in;   /* transfer problem dimensions */
        _info.nrows  = nfun_in;
        _info.nlin   = lvars[0];
        _info.nnlin  = _info.nvars - _info.nlin;
        _info.nobj   = nobj_in;
        _info.nnlequ = _info.nrows;

        _info.ngcomp   = 0;   /* initialize gcomp call counter */
        _info.nparsh   = 0;   /* initialize parsh call counter */
/*---------------------------------------------------------------------*/
/*   set up problem: allocate memory, check user parms, etc.           */
/*   . assign user-supplied gcomp pointer to lsgrg_user_gcomp          */
/*   . lsgrg_memsetup will allocate lsgrg arrays and will return 0     */
/*     if any allocation fails. also calls lsgrg_setupj to initialize  */
/*     jacobian data sructures                                         */
/* NOTE:  all model array initialization (status,bounds,etc.) MUST     */
/*        be called  from memsetup after the memory has been allocated.*/
/*        Also must be called before setupj is called in memsetup to   */
/*        get nonzeros at initial point                                */
/*---------------------------------------------------------------------*/
    if(p_user_gcomp == NULL) {
#ifdef IO_ENABLED
           sprintf( _info.lsgrg_msgbuffer,
            "\n ERROR - User Supplied Pointer to Function Evaluation"
            "\n         Routine is NULL");
            lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
            _info.lsgrg_return_status = _LSGRG_MISSING_GCOMP;
            return _info.lsgrg_return_status;
     }
     _info.lsgrg_user_gcomp = p_user_gcomp;

     if(_info.kderiv == 2) {
        if(p_user_parsh == NULL) {
#ifdef IO_ENABLED
           sprintf( _info.lsgrg_msgbuffer,
            "\n ERROR - User-Supplied Derivatives Specified (kderiv=2)"
            "\n         But Pointer to Derivative Function is NULL");
            lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
            _info.lsgrg_return_status = _LSGRG_MISSING_PARSH;
            return _info.lsgrg_return_status;
        }
        _info.lsgrg_user_parsh = p_user_parsh;
     }
     if(_info.dbg.grgsub)
       lsgrg_echobnd(&_info,xx, blvar, buvar, blcon, bucon,  lvars,
                "in grgsub before call to memsetup");
/* maxr is set to -1 in set_defaults(), if user has not speciried */
/* a value, set to default(n) here */

        if(_info.maxr == -1)
           _info.maxr = _info.nvars;
/*--------------------------------------------------------------------*/
/*  now set all global parameters which are computed from problem     */
/*    dependent values                                                */
/*                                                                    */
/*   memsetup calls lsgrg_setbounds which may set lsgrg_return_status */
/*   if memsetup returns 0, an allocation failed so return            */
/*   _LSGRG_INSFMEMORY        , otherwise, if return_status is < 0    */
/*   a setup failure occured, so return that code                     */
/*   . the allocation routines free all allocated memory (except     */
/*       char arrays) on an allocation failure                        */
/*--------------------------------------------------------------------*/
        lsgrg_setcomputedparms(&_info);

        i = lsgrg_memsetup( &_info,
              xx,blvar, buvar, blcon, bucon, lvars,nnz_user);
        if(i == 0) return _LSGRG_INSFMEMORY;
        if(i < 0 ) return _info.lsgrg_return_status;

        lsgrg_set_invert_tolerances(&_info);
        if(_info.inprnt) {
           lsgrg_echoparms(&_info);
           lsgrg_tablin(&_info);
        }
        if(_info.dbg.grgsub) lsgrg_printgrad(&_info,
                                  1,"Jacobian at Initial Point");
        _info.lsgrg_setup_done = 1;

#ifdef IO_ENABLED
       if(_info.dbg.grgsub) {
          sprintf(_info.lsgrg_msgbuffer,"\n ..... grgsub setup done ....");
          lsgrg_msg(&_info.io_info,_info.lsgrg_msgbuffer);
       }
#endif
        _info.lsgrg_return_status = _LSGRG_SETUP_SUCCESS;

/* if this is a setup call only, return and copyback userinfo struct */

       if(_info.lsgrg_setup && !lsgrg_fullrun) {
           *user_info = _info; /* store back info values if this is */
                               /* a setup-only call                 */
       if(_info.ipr > 0)
          lsgrg_msg(&_info.io_info,
          "\n ... exiting grgsub problem setup ...");
           if(debug)
              lsgrg_screen_msg(&_info.io_info,
                 "setup info stored to user_info -- exiting setup code");
       return _info.lsgrg_return_status;
       }

     }    /*  end of lsgrg_setup code */
/*---------------------------------------------------------------------*/
/*     code for lsgrg_run  get initial time                            */
/*---------------------------------------------------------------------*/
/* zzzzz test basis inv reallocation */
/* _info.lbinv = 5;  */

        tstart = lsgrg_timer();

#ifdef IO_ENABLED
    if(debug) {
      sprintf(_info.lsgrg_msgbuffer,
         "\n ...at run code: if(_info.lsgrg_run)"
         "\n  _info.lsgrg_run = %d lsgrg_fullrun = %d",
        _info.lsgrg_run,lsgrg_fullrun);
      lsgrg_screen_msg(&_info.io_info,_info.lsgrg_msgbuffer);
    }

#endif
    if(_info.lsgrg_run ||  lsgrg_fullrun) {
        if(_info.ipr > 1)
          lsgrg_msg(&_info.io_info,
          "\n  Beginning Solution run \n");
          if(debug) lsgrg_screen_msg(&_info.io_info,
             "\n .....entry run code .... ");

       if(_info.lsgrg_setup_done==0) { /* no prior setup call */

#ifdef IO_ENABLED
         sprintf(_info.lsgrg_msgbuffer,
           "\n Error: LSGRGR Run Requested Before Setup Call");
          lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
            return _LSGRG_BAD_COMMAND;
        }
/*---------------------------------------------------------------------*/
/*  10/00 jcp if this is a run call following a setup call, copy       */
/*    saved locally allocated var/row name pointers back onto the      */
/*    ones used by tablin/outres   **fixme** check on this             */
/*                                                                     */
/*   1.  check that nvars_in and nfun_in have not changed              */
/*   2.  check that lvars and user_nnz have not changed                */
/*   3.  call lsgrg_setbounds to copy initial x values for this run    */
/*       from xx to x, set bounds values for this run, and project     */
/*       x values onto bounds (and perform some consistency checks)    */
/*                                                                     */
/*  ** 1/3/02 jcp **  reset all fixed variable status codes to         */
/*    nonfixed so that consistency check in setbnd will not            */
/*    detect a problem with variable count                             */
/*                                                                     */
/*                                                                     */
/*---------------------------------------------------------------------*/
        if(_info.lsgrg_run ==1 && (lsgrg_fullrun==0) ) { /* test for */
                                                         /* run only */
         if(debug) lsgrg_screen_msg(&_info.io_info,
         "... checking structure for run-only call");
/*  ** 1/3/02 jcp ** */
          for(i=1;i<=_info.nvars;++i) {
              if(_info.ivstat[i]==V_LINEAR_FIXED) _info.ivstat[i] = V_LINEAR;
              if(_info.ivstat[i]==V_NONLINEAR_FIXED) _info.ivstat[i] = V_NONLINEAR;
          }
#ifdef IO_ENABLED
         if(debug) {
           sprintf(_info.lsgrg_msgbuffer,
           "\n ---- checking problem structure on run call -----\n"
            " _info.nvars = %5d nvars_in = %5d\n"
            " _info.nrows = %5d nfun_in  = %5d\n lvars[0] = %d",
             _info.nvars,nvars_in,_info.nrows,nfun_in,lvars[0]);
          lsgrg_msg(&_info.io_info,_info.lsgrg_msgbuffer);
         }
#endif
/*  check lvars correspondence */
            j = 0;   /* use j as error counter */
            if(_info.nvars != nvars_in || _info.nrows != nfun_in) ++j;
            for(i=1;i<=_info.nvars;++i)   /* use x0 as temp array for */
               _info.x0[i] = V_NONLINEAR; /* lvars check */
            for(i=1;i<=lvars[0];++i) _info.x0[lvars[i]] = V_LINEAR;

#ifdef IO_ENABLED
         if(debug) {
             sprintf(_info.lsgrg_msgbuffer,
                "\n ..... checking variable status .....");
          lsgrg_msg(&_info.io_info,_info.lsgrg_msgbuffer);
        }
#endif
            for(i=1;i<=_info.nvars;++i)  {
#ifdef IO_ENABLED
          if(debug) {
             sprintf(_info.lsgrg_msgbuffer,
                "\n var %3d _info.x0 = %2d _info.ivstat = %2d",i,
                   (int) _info.x0[i],_info.ivstat[i]);
            lsgrg_msg(&_info.io_info,_info.lsgrg_msgbuffer);
          }
#endif
               if( (int) _info.x0[i] !=  _info.ivstat[i]) ++j;
           }
            if(j > 0) {
#ifdef IO_ENABLED
         sprintf(_info.lsgrg_msgbuffer,
           "\n Error: LSGR2C Run Requested with problem structure"
           "\n        different from setup. Check nvars_in, nfun_in"
           "\n        and lvars");
          lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
               return  _LSGRG_PROBLEM_STRUCTURE;
            }
            lsgrg_setbounds(&_info,blvar,buvar,blcon,bucon,xx,
                            lvars,_info.ivstat);

        } /* end setup for run-only call */
/*====================================================================*/
/*  top of algorithm retry loop -- do until fixup_code = EXIT         */
/*====================================================================*/
        fixup_code = EXIT;
        retry_count = 0;
        for(i=0; i < LEN_RETRY_HISTORY; ++i) /* zero out retry_history */
         _info.retry_history[0] = 0;
do {

/*  make initial gcomp call to set initial function values */

        lsgrg_gcomp1(&_info, _info.g, _info.x );
        if(_info.dbg.grgsub)
           gcomp_check(&_info,"at initial call in grgsub",_info.n, _info.x, _info.mp1,_info.g, 1);

        for( i = 1; i <= _info.nrows; i++ )  /* store initial function values */
                _info.g0[i] = _info.g[i];
        for( i = 1; i <= _info.nvars; i++ )  /* store initial variable values */
                _info.x0[i] = _info.x[i];

/*  echo system parameters and initial var/con status if requested */
/*  ** 8/01 jcp do not repeat tablin and report if this is a retry */
/*    fixup_code is EXIT on 1st trip through loop , another value  */
/*    otherwise since test is at bottom                            */
/*    also dont call echoparms since problem structure hasn't      */
/*    changed since setup                                          */

        if( _info.inprnt && fixup_code == EXIT) {
           if(!lsgrg_fullrun) {
   /*        lsgrg_echoparms(&_info);   */
             lsgrg_tablin(&_info);
           }
#ifdef USER_REPORT
          if(fixup_code == EXIT )
           lsgrg_report( _info.g, _info.x, _info.mp1, _info.n,
                   _info.usernames,_info.lsgrg_varnames,
                   _info.lsgrg_rownames,_info.x0);
#endif
        }
/*------------------------------------------------------------------*/
/*  implant test harness values                                     */
/*------------------------------------------------------------------*/
   _info.grgtest.ncols    = nvars_in;
   _info.grgtest.nrows    = nfun_in;
   strcpy(_info.grgtest.s_alg,"lsgrgc3");
/*------------------------------------------------------------------*/
/* function header for grgitn call                                  */
/*------------------------------------------------------------------*/
/*

int   grgitn(long int ibmap[], double zmem[], double grad[],
         long int ihag[], long int iheg[], long int ihegl[], double r[],
         double alb[], double ub[], double x[], double gradf[], double g[],
         double v[], double d[], double u[], double gbest[], double xbest[],
         double ascale[], long int inbv[], long int iub[], double rowb[],
         double colb[], long int ibc[], long int ibv[], long int istat[],
         double xstat[], long int ivstat[], double xb1[], double xb2[],
         double xb3[], double gradfp[], double dbnd[], double gg[], double rr[],
         double y[], long int icols[], long int iobjpt[], long int icand[],
         long int inlin[], long int irank[], double cdnum[], long int *tablst,
         double albo[], double ubo[], long int lksame[], double *ycg,
         double *scg, double *cgscr, long int ibvbct[], double paij[],
         long int iprow[], long int ipcol[], long int ipmap[])

*/
/*-------------------------------------------------------------------*/
/*   The translated algorithm routines use 0-based arrays whereas    */
/*   the shell assumes 1-based arrays.  Increment all array pointers */
/*   by 1.                                                           */
/*   NOTE:  the pointers are global, ( so DO NOT use increment       */
/*   expressions, ( pass (p+1) to grgitn                             */
/*                                                                   */
/*  NOTE: use iheg0 and ihegl0 which have their pointer values       */
/*   decremented by 1 for v1 routines                                */
/*  NOTE: if a fatal error is encountered during the algorithm       */
/*        a longjump will be executed to the point below where       */
/*        setjmp is set. setjmp returns 0 on set call, longjump      */
/*        causes a nonzero value                                     */
/*  ** 2/01 jcp ** wrap a dowhile loop around the entire             */
/*   grgitn call logic to facilitate retrys on certain error         */
/*   conditions                                                      */
/*                                                                   */
/* ** 5/01 jcp ** there is a call to lsgrg_put_globals executed      */
/*  prior to each entry to grgitn and a call to lsgrg_get_globals    */
/*  following the return to move values between the shell values     */
/*  and the equivalent values in the structures which derive from    */
/*  the old fortran common blocks.                                   */
/*                                                                   */
/*  colmax was added to the values xfered back in get_globals since  */
/*  it is updated as the jacobian grows and, if the jacobian is      */
/*  reallocated, the current value of colmax needs to be retained    */
/*  so that the next call to put_globals transfers the latest        */
/*  colmax value, rather than the one computed in the setup phase.   */
/*                                                                   */
/*  NOTE: a call to get_globals was added at the top of the if       */
/*   branch for error terminations to cause the transfer             */
/*                                                                   */
/*  NOTE:  there may be other such values as we shake down the       */
/*    retry mechanism                                                */
/*                                                                   */
/*  NOTE2:  we may be able to get away with simply conditioning      */
/*    the call to put_globals on retry_count so that the transfer    */
/*    is not done on algorithm reentries????                         */
/*                                                                   */
/*                                                                   */
/*-------------------------------------------------------------------*/
   fixup_code = EXIT;  /* exit at bottom if algorithm terminates */
                       /* normally */
      _info.abort = setjmp(lsexit_jmpbuf);
      if(!_info.abort) {

              /* copy globals from shell to alg variables*/
      lsgrg_put_globals(&_info);

      setobj( &_info,((_info.ihag)+1),(_info.iheg0), (_info.iobjpt)+1);

      grgitn_return =  grgitn(&_info,
                           (ibmap), (zdummy), ((_info.grad)+1),
                ((_info.ihag)+1),  (_info.iheg0),  (_info.ihegl0),
                ((_info.r)+1),
                ((_info.alb)+1),   ((_info.ub)+1),     ((_info.x)+1),
                ((_info.gradf)+1), ((_info.g)+1),
                ((_info.v)+1),     ((_info.d)+1),      ((_info.u)+1),
                ((_info.gbest)+1), ((_info.xbest)+1),
                ((_info.ascale)+1),((_info.inbv)+1),   ((_info.iub)+1),
                ((_info.rowb)+1),
                ((_info.colb)+1),  ((_info.ibc)+1),    ((_info.ibv)+1),
                ((_info.istat)+1),
                ((_info.xstat)+1), ((_info.ivstat)+1), ((_info.xb1)+1),
                ((_info.xb2)+1),
                ((_info.xb3)+1),   ((_info.gradfp)+1), ((_info.dbnd)+1),
                ((_info.gg)+1),    ((_info.rr)+1),
                ((_info.y)+1),     ((_info.icols)+1),  ((_info.iobjpt)+1),
                ((_info.icand)+1),
                ((_info.inlin)+1), ((_info.irank)+1),  ((_info.cdnum)+1),
                ((_info.tablst)+1),
                ((_info.albo)+1),  ((_info.ubo)+1),    ((_info.lksame)+1),
                ((_info.ycg)+1),
                ((_info.scg)+1),   ((_info.cgscr)+1),  ((_info.ibvbct)+1),
                ((_info.paij)+1),
                ((_info.iprow)+1), ((_info.ipcol)+1),  ((_info.ipmap)+1) );

       lsgrg_get_globals(&_info);  /* xfer info back from alg globals */

/* map algorithm termination code into shell code and set status */
       _info.lsgrg_return_status =
            lsgrg_xfer_alg_term_code(&_info,grgitn_return);

/* if maximization and an lp, the obj fcn value is not unnegated */
/* in the alg code, so do it here                                */
    if(maximize && lvars[0]==nvars_in)
        _info.g[nobj_in] = -_info.g[nobj_in];
/*-------------------------------------------------------------------*/
/*    record run statistics for test harness                         */
/*-------------------------------------------------------------------*/
     _info.grgtest.itns    = _info.nsear;
     _info.grgtest.obj     = _info.g[_info.nobj];
     _info.grgtest.ninf    = _info.ninf;
     _info.grgtest.ngcomps = _info.ngcomp;
     _info.grgtest.sinf    = lsgrg_sinf(nfun_in, nobj_in,_info.g,
                           blcon,bucon);
/*-------------------------------------------------------------------*/
/*  if timing enabled, get finish time and report if otprnt          */
/*-------------------------------------------------------------------*/
        if( _info.otprnt) {
           lsgrg_outres(&_info,_info.lsgrg_return_status);

#ifdef TIMING_ENABLED
        lsgrg_timestats(&_info,T_TOTAL,(lsgrg_timer() - tstart) );
        if(_info.otprnt) lsgrg_timestats(&_info,T_PRINT,(double)0.0);
#endif

#ifdef USER_REPORT
           lsgrg_report( _info.g, _info.x, _info.mp1, _info.n, _info.usernames,
                        _info.lsgrg_varnames, _info.lsgrg_rownames, _info.x0);
#endif
        }

/*    return final values  */

        for( i = 1; i <= _info.n;   i++ )     xx[i]   = _info.x[i];
        for( i = 1; i <= _info.mp1; i++ )     fcns[i] = _info.g[i];

/*  7/01 jcp add capability to return mults unpacked */
/* RETURN MULTIPLIERS OF BINDING CONSTRAINTS */
        if(_info.return_mults_unpacked) {
           for( i=1; i <= _info.mp1; ++i) rmults[i] = 0.0;
           for( i=1; i <= _info.nbc; ++i) {
               j = _info.ibc[i];
               inbind[i] = j ;
               rmults[j] = _info.u[j];
           }
        } else {
             for( i = 1; i <= _info.nbc; i++ ) {
                 j = _info.ibc[i];
                 inbind[i] = j ;
                 rmults[i] = _info.u[j];
             }
          }
/*   RETURN REDUCED GRADIENT OF NONBASIC VARIABLES */

        for( i = 1; i <= _info.n; i++ ) {
                nonbas[i] = 0;
                redgr[i] = 0.0;
        }
        if(_info.return_mults_unpacked) {
           for( i = 1; i <= _info.n; i++ ) {
               nindex = _info.inbv[i];
              if( nindex < _info.nvars ) {
                nonbas[nindex] = nindex;
                redgr[nindex] = _info.gradf[i];
              }
            }
         } else {
                   j = 1;
                   for( i = 1; i <= _info.n; i++ ) {
                      nindex = _info.inbv[i];
                      if( nindex < _info.nvars ) {
                          nonbas[j] = nindex;
                          redgr[j] = _info.gradf[i];
                         ++j;
                      }
                   }
           }
/*  set scalar return arguments */
        *nnonb = j;
        *nbind = _info.nbc;
        fixup_code = EXIT;  /* set value to exit enveloping dowhile */

    } /* end of if(!_info.abort) 'normal' algorithm termination  */

    else {
/*------------------------------------------------------------------*/
/*  algorithm terminated with call to lsgrg_errorexit               */
/*  call fixup_strategy to set fixup_code, then take specified      */
/*  action                                                          */
/* ** 2/01 jcp ** eventually parametrize the perturbation values    */
/*                                                                  */
/* ** 5/01 jcp ** call get_globals to make sure that appropriate    */
/*  values are moved back to the shell so that the call to          */
/*  put_globals preceeding the grgitn call puts them back           */
/*                                                                  */
/* ** 7/01 jcp ** if fixups are disabled by user, just exit         */
/*------------------------------------------------------------------*/
       lsgrg_get_globals(&_info);  /* xfer info back from alg globals */

              fixup_code = lsgrg_fixup_strategy(&_info,_info.abort,
                                     &retry_count,max_retry_count,
                                     error_classification,
                                     _info.retry_history);
          if(_info.disable_retries) fixup_code = EXIT;
#ifdef IO_ENABLED
         if(_info.ipr > 2 ) {
           sprintf(_info.lsgrg_msgbuffer,
             "\n Algorithm Returned Error Status Code = %d"
             "\n     [%s]",_info.abort,error_classification);
            lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
          }
         else if(_info.ipr > 1) {
                  sprintf(_info.lsgrg_msgbuffer,
                     "\n Algorithm Returned Fatal Error "
                    "\n     [%s]\n",error_classification);
                  lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
          }
#endif
              switch(fixup_code) {
             case EXIT:      /* fall out of bottom of loop */
                    break;
             case RESTART:  /*  restart at current point */
                    break;
             case  SCALE:  /*   switch to scaling if not already on */
                   if(_info.iscale != 1) _info.iscale = 1;
                   break;
             case PERTURB_POINT:  /* perturb point slightly */
                for(i=1;i<=nvars_in;++i) {
                     if(i%2==0)
                        temp = _info.x[i] * 1.05;
                     else
                        temp = _info.x[i] * -1.05;
                     if(temp < blvar[i]) temp = _info.x[i] * 1.05;
                     if(temp > buvar[i]) temp = _info.x[i] * -1.05;
                     if(temp < blvar[i] || temp > buvar[i])
                          temp = blvar[i] + (buvar[i] - blvar[i]) / 2.0;
                     _info.x[i] = temp;
                 }
                 break;
          case RESIZE_BINV:
              if(lsgrg_resize_binv(&_info)==0)
                   fixup_code = EXIT;
                   _info.abort = _LSGRG_INSFMEMORY;
                break;
          case RESIZE_JAC:
                if(lsgrg_resize_jacobian(&_info,_info.iheg0)==0)
                   fixup_code = EXIT;
                   _info.abort = _LSGRG_INSFMEMORY;
                break;
          default:
#ifdef IO_ENABLED
         sprintf(_info.lsgrg_msgbuffer,
           "\n Internal Error: fixup_strategy returned unknown "
            "fixup code = %d",fixup_code);
          lsgrg_error_msg(&_info.io_info,_info.lsgrg_msgbuffer);
#endif
          fixup_code = EXIT;
          break;
      }  /* end of switch on fixup_code */

           _info.lsgrg_return_status = _info.abort; /* last error code */

/*-------------------------------------------------------------------*/
/* return last point visited                                         */
/* if fixup_code is EXIT and termination was not user initiated      */
/* call gcomp to set the functions to correspond to this point       */
/* and set whether this point is feasible or not                     */
/*-------------------------------------------------------------------*/
        for( i = 1; i <= nvars_in; i++ )     xx[i]   = _info.x[i];
        if(fixup_code==EXIT && _info.abort != _LSGRG_USER_TERMINATION)
              lsgrg_gcomp1(&_info, fcns, xx );
        else
            for(i=1; i <= nfun_in; ++i) fcns[i] = _info.g[i];

              ninf = 0;
              for(i=1;i<=nfun_in;++i) {
                  if(fcns[i] < blcon[i]) temp = blcon[i] - fcns[i];
                  else
                     if(fcns[i] > bucon[i]) temp = fcns[i] - bucon[i];
                         else continue;
                  if( temp/(1.0e0 + fabs(fcns[i])) > _info.limits.epnewt)
                       ++ninf;
               }
              _info.last_point_feasible = ninf > 0;
/*
  if we are exiting with an error condition, copy last_point_feasible
  and retry_history to the users info_structure
*/
              if(fixup_code == EXIT) {
                 user_info->last_point_feasible = _info.last_point_feasible;
                 for(i=0; i< LEN_RETRY_HISTORY;++i)
                   user_info->retry_history[i] = _info.retry_history[i];
              }

         } /* end if !info.abort  else block */

       } while(fixup_code != EXIT); /* loop until 'normal' termination */
                                    /* or fixup options exhausted      */

/*  if this is a run-only call, return with last status code set */
/*  comment this out and let logic below do copyback and return */
      /* xxxxx  if(!lsgrg_fullrun) return _info.lsgrg_return_status; */

     } /* end of lsgrg_run code */
/*-------------------------------------------------------------------*/
/*  ** 7/26/03 jcp ** if this is a fullrun, reset maxr to -1         */
/*  so that if the user uses the same LsgrgInfo structure for        */
/*  multiple runs, maxr will be reset to the actual nbr of vars      */
/*  #### set the value in the user structure AND the local           */
/*   just in case we are doing a copyback                            */
/*  **** doing this is probably dangerous (for a user) ****          */
/*                                                                   */
/*-------------------------------------------------------------------*/
    if(lsgrg_fullrun) {
       _info.maxr = -1;
       user_info->maxr = -1;
   }
/*-------------------------------------------------------------------*/
/*       lsgrg_shutdown code                                         */
/* 11/00 jcp free the varname/rowname storage allocated from         */
/*    grgsub. then copy the user's pointers back on top of the       */
/*    lsgrg_varnames/lsgrg_rownames pointers to allow LsgrgBye       */
/*    to free that storage.                                          */
/* 4/27/02 jcp ** delete additional var/row name storage ptrs        */
/*-------------------------------------------------------------------*/
     if(_info.lsgrg_shutdown || lsgrg_fullrun ) {
        lsgrg_free_resources(&_info);
        _info.lsgrg_setup_done = 0;  /* unset setup status */


        if(_info.lsgrg_return_status == _LSGRG_STATUS_NOT_SET)
           _info.lsgrg_return_status = _LSGRG_SHUTDOWN_SUCCESS;
     }  /*  end of shutdown code */
/*-------------------------------------------------------------------*/
/*      END OF GRGSUB                                                */
/*  set return code into  grgtest results structure                  */
/*  copy local grgtest structure back to user structure              */
/*  /* 7/01 jcp ** if copyback not disabled, copy local info struct  */
/*    back to users structure                                        */
/*  NOTE:  copyback can only be disabled for a fullrun.  code above  */
/*  at entry checks for disabled copyback in a multirun setup and    */
/*  exits with an error message                                      */
/*-------------------------------------------------------------------*/
        _info.grgtest.info = _info.lsgrg_return_status;

  if(!_info.disable_copyback) {
         user_info->grgtest = _info.grgtest;
         *user_info  = _info;
   }
        return _info.lsgrg_return_status;
} /* end of function grgsub */

void LSGRGCLASS lsgrg_echobnd(LsgrgInfo *info,
                    double xx[],double blvar[], double buvar[],
                    double blcon[],double bucon[], long lvars[],
                    char *msg)
{
#ifdef IO_ENABLED
   int i;
   sprintf(info->lsgrg_msgbuffer,
     "\n............\n bounds and x echo: msg = %s",msg);
   lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
   sprintf(info->lsgrg_msgbuffer,
    "\n lvars[0] = %d",lvars[0]);
   lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
   for(i=1;i<=info->nvars; ++i) {
      sprintf(info->lsgrg_msgbuffer,
       "\n x[%5d] = %15.8e blvar = %15.8e buvar = %15.8e",
         i,xx[i],blvar[i],buvar[i]);
       lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
      if(lvars[0]) {
        sprintf(info->lsgrg_msgbuffer,
        " \n      lvars = %d", lvars[i]);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
      }
   }
   for(i=1;i<=info->nrows;++i) {
      sprintf(info->lsgrg_msgbuffer,
      "\n   row = %5d blcon = %15.8e bucon = %15.8e",
      i,blcon[i],bucon[i]);
     lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
   }
   sprintf(info->lsgrg_msgbuffer,"\n..........\n");
#endif
    return;
}
void LSGRGCLASS lsgrg_errorexit(jmp_buf lsexit_jmpbuf,int errorcode)
{
     longjmp(lsexit_jmpbuf,errorcode);
}
/* **fixme** break here to file grgmisc.c */
void LSGRGCLASS lsgrg_initialize_grgtest(LsgrgInfo *info)
{
/* .... initialize grgtest global info fields at start of run ...*/
   info->grgtest.obj     = 0.0;
   info->grgtest.sinf    = 0.0;
   info->grgtest.time    = 0.0;
   info->grgtest.ncols   = 0  ;
   info->grgtest.nrows   = 0  ;
   info->grgtest.itns    = 0  ;
   info->grgtest.ninf    = 0  ;
   info->grgtest.ngcomps = 0  ;
   info->grgtest.info    = -99;
   strcpy( info->grgtest.s_alg,        "<xxx>");
   strcpy( info->grgtest.s_termin,     "<xxx>");
   strcpy( info->grgtest.s_termination,"<xxx>");
   strcpy( info->grgtest.s_datetime,   "<xxx>");
   return;
}
char* LSGRGCLASS lsgrg_get_terminationmsg(LsgrgInfo *info)
{   return info->grgtest.s_termination; }

int LSGRGCLASS lsgrg_xfer_alg_term_code(LsgrgInfo *shell,int info)
{
/*--------------------------------------------------------------------*/
/*  lsgrg_xfer_alg_term_code() matches the algorithm termination code */
/*  to the corresponding shell termination code.  The algorithm codes */
/*  do not exactly correspond to the more detailed codes in the       */
/*  lsgrg_termination_codes enum, so this routine avoids the need to  */
/*  repeatedly re-align the codes, since the algorithm codes are      */
/*  parameterized as well.  If we put more or more detailed codes     */
/*  into the algorithm, simply extend the switch stmt below           */
/*--------------------------------------------------------------------*/
     switch(info) {
       case _LSGRGALG_KTC:         return _LSGRG_KTC;
       case _LSGRGALG_FC:
           if(shell->ninf == 0) return _LSGRG_FRACTCHG;
           return _LSGRG_INFEASIBLE_FRACTCHG;

       case _LSGRGALG_ALLREMEDIES:
           if(shell->ninf == 0) return _LSGRG_ALLREMEDIES;
           return _LSGRG_INFEASIBLE_ALLREMEDIES;

       case _LSGRGALG_ITN_LIMIT:
           if(shell->ninf == 0) return _LSGRG_ITERATIONS;
           return _LSGRG_INFEASIBLE_ITERATIONS;

       case _LSGRGALG_UNBOUNDED:   return _LSGRG_UNBOUNDED ;
       case _LSGRGALG_INF:         return _LSGRG_INFEASIBLE;
       case _LSGRGALG_RUNTIME_ERR: return _LSGRG_OTHER_RUNTIME;
       case _LSGRGALG_USER_TERMIN: return _LSGRG_USER_TERMINATION;
       case _LSGRGALG_TIME_LIMIT:  return _LSGRG_TIME_LIMIT_EXCEEDED;
    };
    return shell->lsgrg_return_status;
}
void LSGRGCLASS lsgrg_set_termin_msg(LsgrgInfo *info,char *buffer)
{
      strcpy(info->grgtest.s_termination,buffer);
}
void LSGRGCLASS lsgrg_set_termin_code(LsgrgInfo *info,char *string)
{
     strcpy(info->grgtest.s_termin,string);
}
double LSGRGCLASS lsgrg_sinf(long m, long nobj,double g[], double blcon[],
                  double bucon[])
{
    long i; double sinf=0.0;
    for(i=1;i<=m;++i) {
        if(i==nobj) continue;
        if( (g[i] - bucon[i]) > 0.0) sinf += (g[i] - bucon[i]);
        else
           if((blcon[i] - g[i]) > 0.0 ) sinf += (blcon[i] - g[i]);
    }
    return sinf;
}
int LSGRGCLASS lsgrg_fixup_strategy(LsgrgInfo *_info,int return_code,
                        int *retry_count, int max_retry_count,
                        char *msg, int  retry_history[])
{
      char nn[10];
      int i;
/*---------------------------------------------------------------*/
/*  set msg to the text of the error code and add this error code*/
/*  to the list in retry_history  so that the history is saved   */
/*  **fixme** for now, use a fixed length array of length        */
/*  LEN_RETRY_HISTORY.  if this length is exceeded, simply shift */
/*  codes up and add this code to the end.  we shouldn't need    */
/*  a lot of codes here (default=14) so a dynamic table is not   */
/*  needed                                                       */
/*                                                               */
/*  return error code name in msg for message issuance in caller */
/*                                                               */
/*                                                               */
/*---------------------------------------------------------------*/
       msg[0] = '\0';
       if(++retry_history[0] >= LEN_RETRY_HISTORY - 1) {
          for(i=1;i< LEN_RETRY_HISTORY-1;++i)
              retry_history[i-1] = retry_history[i];
          retry_history[LEN_RETRY_HISTORY-1] = return_code;
       }
       else
          retry_history[retry_history[0]] = return_code;

       lsgrg_int_to_string((*retry_count+1),nn);
       switch(return_code) {
/*---------------------------------------------------------------*/
/*  The codes below involve singularities, bad conditioning,     */
/*  and other algorithm terminations which might be the result   */
/*  of numerical difficulties and/or particularly unsiutable     */
/*  points.  we initiate a sequence of retries which include     */
/*      cold restart from current point, switch to scaling (if   */
/*  not already on), restart from a perturbed point,             */
/*  or anything else that we can think of                        */
/*                                                               */
/*  To add fixups, add fixup name to enum fixup_codes, then add  */
/*  switch statement below
/*---------------------------------------------------------------*/
/* --- set string with name of error condition for return to */
/*     grgsub and inclusion in msg (if io is enabled)        */

  case _LSGRG_REDOBJ_CONSTVIOL:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_REDOBJ_CONSTVIOL");
  case _LSGRG_REDGRA_NB_LE_0:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_REDGRA_NB_LE_0 ");
  case _LSGRG_CHUZQ_BADPIVOT:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_CHUZQ_BADPIVOT ");
  case _LSGRG_XPIVOT_BASIS_ILLCOND:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_XPIVOT_BASIS_ILLCOND");
  case _LSGRG_XPIVOT_BASIS_SING:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_XPIVOT_BASIS_SING");
  case _LSGRG_PH0FAC_INVERT_FAILURE:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_PH0FAC_INVERT_FAILURE");
  case _LSGRG_PH0PIV_XPIVOT_FAILURE:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_PH0PIV_XPIVOT_FAILURE");
  case _LSGRG_CONSBS_BASIC_SLACK:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_CONSBS_BASIC_SLACK");
  case _LSGRG_CONSBS_NOINVERT_SEARCH:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_CONSBS_NOINVERT_SEARCH");
  case _LSGRG_CONSBS_NB_TOOBIG:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_CONSBS_NB_TOOBIG");
/*
  retry strategy for stuck problems or problems with possible
  numeric difficulties
*/
        switch (++(*retry_count)) {
        case 1:  /*  restart from current point */
           strcat(msg,"\n restarting algorithm from current point "
                         "-- retry count= ");
           strcat(msg,nn);
           return RESTART;
        case 2: /* switch to scaling if not already on */
                /*  if already on, fall through to next fixup */
           if(!_info->iscale) {
              strcat(msg,"\n restarting algorithm from current point "
                            "with scaling  -- retry count = ");
              strcat(msg,nn);
              return SCALE;
           }
         case 3: /* restart from perturbed point */
             strcat(msg,"\n restarting algorithm from preturbed point "
                            "-- retry count = ");
              strcat(msg,nn);
              return PERTURB_POINT;
          default: /* exit after all remedies exhausted */
              strcat(msg,"\n all remedies exhausted -- exiting --"
                       " retries = ");
              lsgrg_int_to_string((*retry_count-1),nn);
              strcat(msg,nn);
                          strcat(msg,"\n");
              return EXIT;
        }
        return EXIT;
/*---------------------------------------------------------------*/
/*  column length errors, most of these should never occur since */
/*  the jacobian growth code      should cause a reset of colmax */
/*  when a column grows in length                                */
/*  if they are encountered, they should be treated as internal  */
/*  shouldNeverOccur errors                                      */
/*  for now, just exit  >>>> am I right about this<<<< ????      */
/*---------------------------------------------------------------*/
  case _LSGRG_XDOT_COLLEN:
      if(msg[0]=='\0')
      strcpy(msg,"_LSGRG_XDOT_COLLEN");
  case _LSGRG_XSAXPY_COLLEN:
      if(msg[0]=='\0')
      strcpy(msg,"_LSGRG_XSAXPY_COLLEN");
  case _LSGRG_XPIVOT_COLLEN:
      if(msg[0]=='\0')
      strcpy(msg,"_LSGRG_XPIVOT_COLLEN");
  case _LSGRG_CONDNM_BAD_COLLEN:
      if(msg[0]=='\0')
      strcpy(msg,"_LSGRG_CONDNM_BAD_COLLEN");

      strcat(msg,"\n  Internal Error -- Terminating Execution");
        return EXIT;
/*---------------------------------------------------------------*/
/*  insufficient memory errors.  Return RESIZE_BINV for          */
/*  basis inverse storage exceeded and RESIZE_JAC for jacobian   */
/*  growth allowance exceeded.                                   */
/*---------------------------------------------------------------*/
  case _LSGRG_GETBAS_INSFMEM: /* basis inverse out of space */
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_GETBAS_INSFMEM");
  case _LSGRG_XPIVOT_INSFMEM: /* basis inverse out of space */
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_XPIVOT_INSFMEM");
      strcat(msg,": Basis Inverse Space Allocation Exceeded");
      return RESIZE_BINV;

  case _LSGRG_JAC_OVERFLOW:   /* jacobian out of space      */
       if(msg[0]=='\0') {
          strcpy(msg,"_LSGRG_JAC_OVERFLOW");
          strcat(msg,": Jacobian Growth Allowance Exceeded");
        }
        return RESIZE_JAC;
/*---------------------------------------------------------------*/
/*  I think that these errors result from internal consistency   */
/*  checks and (1) should never occur, (2) are indicative of a   */
/*  bug if they do.  all we can do is exit and examine           */
/*---------------------------------------------------------------*/
  case _LSGRG_CONSBS_BASIS_STRUCTURE:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_CONSBS_BASIS_STRUCTURE");
  case _LSGRG_CONSBS_REINVERT_BC:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_CONSBS_REINVERT_BC");
  case _LSGRG_PH0PIV_BAD_INDEX:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_PH0PIV_BAD_INDEX");
  case _LSGRG_PH0PIV_BAD_ICOLS1:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_PH0PIV_BAD_ICOLS1");
  case _LSGRG_PH0PIV_BAD_ICOLS2:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_PH0PIV_BAD_ICOLS2");
  case _LSGRG_CONDNM_BAD_NBLOCKS:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_CONDNM_BAD_NBLOCKS");
  case _LSGRG_DIREC_UPDATE_ERR:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_DIREC_UPDATE_ERR");
  case _LSGRG_CONSBS_JPIV_0:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_CONSBS_JPIV_0");
  case _LSGRG_XPIVOT_OTHER_ERR:
      if(msg[0]=='\0')
         strcpy(msg,"_LSGRG_XPIVOT_OTHER_ERR");

         strcat(msg,"\n   Internal Error -- Terminating Execution");
         return EXIT;
/*---------------------------------------------------------------*/
/*  user termination should cause an exit.  for now,             */
/*  there is no code which returns the _LSGRG_OTHER_RUNTIME      */
/*  code. I put it in as a catchall if needed                    */
/*---------------------------------------------------------------*/
  case _LSGRG_OTHER_RUNTIME:
      if(msg[0]=='\0')
      strcpy(msg,"_LSGRG_OTHER_RUNTIME");
  case _LSGRG_USER_TERMINATION:
      if(msg[0]=='\0')
      strcpy(msg,"_LSGRG_USER_TERMINATION"
         ": Execution Terminated By User");

         return EXIT;

  default:  /* unknown or return code added to codes.h but */
            /* not present in this list.  return exit      */
            /* and add to this list                        */
#ifdef IO_ENABLED
         sprintf(LSGRG_MSGBUFFER,
          "\n Error: lsgrg_fixup_strategy detected an unknown"
          " termination code = %d",return_code);
          lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
         strcpy(msg," Error in lsgrg_fixup_strategy"
            ": Unknown Error Code = ");
         lsgrg_int_to_string(return_code,nn);
         strcat(msg,nn);

         return EXIT;
   } /* end of switch */
}
int LSGRGCLASS lsgrg_get_retry_history(LsgrgInfo *_info,
                         int retry_history[], int *len_retry_history)
{
    int i,nbr_retrys;
/*------------------------------------------------------------------*/
/*     return retry history in users array                          */
/*  nbr_retrys is real retry count, len_retry_history is nbr of     */
/*  stored retry codes                                              */
/*------------------------------------------------------------------*/
    nbr_retrys = _info->retry_history[0];
    if(nbr_retrys < LEN_RETRY_HISTORY-1)
          *len_retry_history = nbr_retrys;
    else
      *len_retry_history = LEN_RETRY_HISTORY-1;
    for(i=1;i<*len_retry_history;++i)
         retry_history[i] = _info->retry_history[i];
    return nbr_retrys;
}
void LSGRGCLASS LsgrgProblemSetup(LsgrgInfo *info)
{
    info->lsgrg_setup = 1; info->lsgrg_run = 0; info->lsgrg_shutdown = 0;
}
void LSGRGCLASS LsgrgProblemRun(LsgrgInfo *info)
{
    info->lsgrg_setup = 0; info->lsgrg_run = 1; info->lsgrg_shutdown = 0;
}
void LSGRGCLASS LsgrgProblemShutdown(LsgrgInfo *info)
{
    info->lsgrg_setup = 0; info->lsgrg_run = 0; info->lsgrg_shutdown = 1;
}
int LSGRGCLASS lsgrg_resize_jacobian(LsgrgInfo *_info,long *iheg) {
/*-------------------------------------------------------------------*/
/*  This routine doubles the jacobian growth allowance               */
/*  then the size of the jacobian (grad/ihag) is increased to        */
/*         (original allocation + original growth allowance)  +      */
/*  (new jacobian growth allowance.  the old values are copied to    */
/*   the new arrays                                                  */
/*  use the nostats version of allocate_darray so that an additional */
/*  pointer is not added to the pointer table                        */
/*-------------------------------------------------------------------*/
    double *pold,*pnew;
    long *pold1,*pnew1;
    long ncols = _info->n + _info->mp1,i;
    long nnz = iheg[ncols]-1,new_length;



    long old_growth_allowance = _info->jacobian_growth_allowance;
/*  increase growth allowance */
    _info->jacobian_growth_allowance = (long) (_info->jacobian_growth_allowance *
          _info->jacobian_growth_realloc_factor) ;
    new_length = nnz + _info->jacobian_growth_allowance;

#ifdef IO_ENABLED
    if(_info->ipr > 0) {
       sprintf(LSGRG_MSGBUFFER,
       "\n Reallocating Jacobian Arrays (grad/ihag)"
       "\n   Old length(including new nonzeros) = %10d, new length = %10d"
       "\n   Jacobian Growth Allowance increased from %10d to %10d\n",
        nnz,new_length,old_growth_allowance,_info->jacobian_growth_allowance);
        lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
     }
#endif
/*---------------------------------------------------------------------*/
/*  reallocate grad                                                    */
/*  copy from 0 to <= nnz since grad+1 is passed to grgitn and         */
/*    the [1] element is the [0] element inside of grgitn              */
/* NOTE: replace the entry in the ptable pointer table with the new    */
/*    pointer so that the free_resources call will free this allocation*/
/*  if the new allocation fails, the old pointer is not touched        */
/*  so that the free_resources call will free that memory              */
/*---------------------------------------------------------------------*/
/*
         printf(" pointer tbl loc for grad = %d"
                 "\n  pointer tbl loc for ihag = %d",_info->tbl_loc_grad,
                 _info->tbl_loc_ihag);
*/
      pold = _info->grad;
      if( (pnew = allocate_darray_nostats(_info,"grad",
                   new_length))==NULL) {
                   lsgrg_free_resources(_info);
                   return 0;
      }

    (_info->ptable)[_info->tbl_loc_grad] = NULL;
    (_info->ptable)[_info->tbl_loc_grad] = pnew;
     _info->grad = pnew;
     for(i=0;i<=nnz;++i) pnew[i] = pold[i];
     free(pold);

/*---------------------------------------------------------------------*/
/*  reallocate ihag                                                    */
/*  copy from 0 to <= nnz since iheg+1 is passed to grgitn and         */
/*    the [1] element is the [0] element inside of grgitn              */
/*  if the new allocation fails, the old pointer is not touched        */
/*  so that the free_resources call will free that memory              */
/*---------------------------------------------------------------------*/
      pold1 = _info->ihag;

      if( (pnew1 = allocate_larray_nostats(_info,"ihag",
                   new_length))==NULL)  {
                   lsgrg_free_resources(_info);
                   return 0;
      }
    (_info->ptable)[_info->tbl_loc_ihag] = NULL;
    (_info->ptable)[_info->tbl_loc_ihag] = pnew1;
     _info->ihag = pnew1;
     for(i=0;i<=nnz;++i) pnew1[i] = pold1[i];
     free(pold1);

     return 1;
}
 int LSGRGCLASS lsgrg_resize_binv(LsgrgInfo *_info) {
/*-------------------------------------------------------------------*/
/*  this routine increases the size of the basis inverse by a factor */
/*  binv_realloc_factor set in lsgrg_set_defaults                    */
/*  use the nostats version of allocate_darray so that an additional */
/*  pointer is not added to the pointer table                        */
/* NOTE: replace the entry in the ptable pointer table with the new    */
/*    pointer so that the free_resources call will free this allocation*/
/*  if the new allocation fails, the old pointer is not touched        */
/*  so that the free_resources call will free that memory              */
/*-------------------------------------------------------------------*/
    int lbinv_old = _info->lbinv;
    double *pnew,*pold = _info->inv_binv;

    _info->lbinv = (long) (_info->lbinv * _info->binv_realloc_factor);
      if( (pnew = allocate_darray_nostats(_info,"inv_binv",
                   _info->lbinv))==NULL) {
                   lsgrg_free_resources(_info);
                   return 0;
      }

    (_info->ptable)[_info->tbl_loc_inv_binv] = NULL;
    (_info->ptable)[_info->tbl_loc_inv_binv] = pnew;
     _info->inv_binv = pnew;
     free(pold);
#ifdef IO_ENABLED
     if(_info->ipr > 0) {
        sprintf(LSGRG_MSGBUFFER,
          "\n Reallocating Storage for Basis Inverse from %8d dwords"
          " to %8d dwords",lbinv_old,_info->lbinv);
        lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
     }
#endif
     return 1;
}
void LSGRGCLASS lsgrg_int_to_string(int n, char s[]){
/* converts positive integer to string */
         int i=0,j,t;
         do {
               s[i++] = n % 10 + '0';
             } while( (n /= 10) > 0);
         s[i] = '\0';
         j = i-1;
         n = i;
         for(i = 0; i < n/2; ++i) {
              t = s[i];
              s[i] = s[j];
              s[j--] = t;
         }
}
void LSGRGCLASS lsgrg_disable_retries(LsgrgInfo *info) {
      info->disable_retries = 1;
}
void LSGRGCLASS lsgrg_enable_retries(LsgrgInfo *info) {
      info->disable_retries = 0;
}
void LSGRGCLASS lsgrg_return_mults_unpacked(LsgrgInfo *info) {
      info->return_mults_unpacked = 1;
}
void LSGRGCLASS lsgrg_return_mults_packed(LsgrgInfo *info) {
      info->return_mults_unpacked = 0;
}
void LSGRGCLASS lsgrg_disable_copyback(LsgrgInfo *info) {
      info->disable_copyback = 1;
}
void LSGRGCLASS lsgrg_enable_copyback(LsgrgInfo *info) {
      info->disable_copyback = 0;
}
#include "lsinfo.h"
#define IOINFO  &(_info->io_info)
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)
/***********************************************************************
 **                 FILE:        GRGSEAR FORTRAN                       *
 **                 AUTHOR:      STUART SMITH                          *
 **                 LAST UPDATE: 27 FEB  1992                          *
 **        **************************************************          *
 **        **************************************************          *
 **        ***               LSGRG2C                      ***          *
 **        ***          COPYRIGHT   1991,1998             ***          *
 **        ***            LEON S. LASDON                  ***          *
 **        ***                 &                          ***          *
 **        ***            ALLAN D. WAREN                  ***          *
 **        ***                 &                          ***          *
 **        ***            STUART SMITH                    ***          *
 **        ***                 &                          ***          *
 **        ***         John C. Plummer                    ***          *
 **        **************************************************          *
 **        **************************************************          *
 **                                                                    *
 **                                                                    *
 **      This file contains the following routines for LSGRG2, a       *
 **  large scale implementation of GRG2:                               *
 **                                                                    *
 **  ROUTINE           DESCRIPTION                        LAST UPDATE  *
 ** ---------         -------------                      ------------- *
 **  PH1OBJ     Computes the PHASE I objective            14 AUG 1989  *
 **  QUAD       Compute quadratic initial estimates       05 MAY 1988  *
 **  TANG       Computes the tangent vector               17 AUG 1989  *
 **  REDGRA     Computes the Reduced Gradient             27 FEB 1992  *
 **  GETCB      Allows direct acces to the individual     04 MAR 1994  *
 **             Obj Function Gradient elements                         *
 **  SEARCH     Performs the one-dimensional line search  27 FEB 1992  *
 ***********************************************************************
 *
 * */
void /*FUNCTION*/ LSGRGCLASS ph1obj(LsgrgInfo *_info,long int ibc[], double g[], double ub[],
         double alb[])
{
long int i, i_, j;
double bnd, factor, gj, sinf, t;

/* replace statics with alg members */

    #define sinf0  _info->ph1objStatic.sinf0
    #define true0  _info->ph1objStatic.true0


/*
   static double sinf0 = 0.0e0;
   static double true0 = 0.0e0;
*/


        /* .....................................................................
         *     THIS SUBROUTINE COMPUTES PHASE 1 OBJECTIVE = SUM OF ABSOLUTE
         *     VALUES OF CONSTRAINT VIOLATIONS AND STORES AS G(NOBJ). TRUE
         *     OBJECTIVE SAVED AS TRUOBJ
         * ..................................................................... */

        if( _info->nintbk.ipr >= 5 )
                lsgrg_msg(IOINFO," ****ENTERING SUBROUTINE: PH1OBJ\n" );


        factor = 100.0e0;
        _info->nintbk.ninf = 0;
        sinf = 0.0e0;
        if( _info->bind.nnbc == 0 )
                goto L_300;


        for( i = 1; i <= _info->bind.nnbc; i++ ){
                i_ = i - 1;
                j = ibc[_info->bind.nbc + i_];
                gj = g[j - 1];
                bnd = alb[_info->dimen.n + j - 1];
                t = bnd - gj;
                if( t <= _info->limits.epnewt*(1.0e0 + fabs( bnd )) )
                        goto L_100;
                _info->nintbk.ninf = _info->nintbk.ninf + 1;
                sinf = sinf + t;
                goto L_200;
L_100:
                bnd = ub[_info->dimen.n + j - 1];
                t = gj - bnd;
                if( t <= _info->limits.epnewt*(1.0e0 + fabs( bnd )) )
                        goto L_200;
                _info->nintbk.ninf = _info->nintbk.ninf + 1;
                sinf = sinf + t;
L_200:
                ;
        }

L_300:
        _info->bestbk.truobj = g[_info->nintbk.nobj - 1];
        if( _info->nintbk.ninf == 0 )
                return;
        t = _info->bestbk.truobj;
        if( _info->optblk.maxim )
                t = -t;

        /*     If in middle of search, do not recompute PHMULT
         * */
        if( _info->ph1bk.initph == 0 )
                goto L_400;

        /*     If starting PHASE 1, compute PHMULT from scratch
         * */
        if( _info->ph1bk.initph == 2 )
                goto L_350;
        _info->bestbk.truobj = _info->redph.trubst;
        t = _info->bestbk.truobj;
        if( _info->optblk.maxim )
                t = -t;

        /*     Right after CONSBS call, see if SINF has changed enough to restart
         * */
        if( _info->ph1bk.phmult == 0.0e0 )
                goto L_350;
        if( sinf0 > factor*sinf )
                goto L_350;
        if( true0*_info->bestbk.truobj < 0.0e0 )
                goto L_350;
        if( fabs( true0 ) < factor*fabs( _info->bestbk.truobj ) && fabs( _info->bestbk.truobj ) <
         factor*fabs( true0 ) )
                goto L_400;
L_350:
        ;
        _info->ph1bk.phmult = 0.0e0;
        sinf0 = sinf;
        true0 = _info->bestbk.truobj;
        _info->logblk.restrt = TRUE;
        if( fabs( _info->bestbk.truobj ) < 0.01e0 )
                goto L_400;
        _info->ph1bk.phmult = fabs( _info->ph1bk.ph1eps*sinf/_info->bestbk.truobj );
        if( _info->optblk.maxim )
                _info->ph1bk.phmult = -_info->ph1bk.phmult;
L_400:
#ifdef FILE_IO_ENABLED
        if( _info->nintbk.ipr3 >= 2 )
                {
                sprintf(LSGRG_MSGBUFFER, " NEW PHMULT =%12.5e  SUM OF INFEASIBILITIES = %15.8e  TRUOBJ =%15.8e\n",
                 _info->ph1bk.phmult, sinf, t );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

        if( _info->nintbk.ninf != 0 )
                g[_info->nintbk.nobj - 1] = sinf + _info->ph1bk.phmult*_info->bestbk.truobj;
        /*@@@@ IF (PH1EPS.NE.0.0D0.AND.IPR3.GE.2) WRITE (IOOUT,510) SINF,T,NINF */

#ifdef IO_ENABLED
        if( (_info->nintbk.ninf != 0) && (_info->nintbk.ipr >= 5) )
                {
                sprintf(LSGRG_MSGBUFFER, " SINF =%15.8e TRUOBJ =%15.8e NINF = %6ld\n",
                 sinf, t, _info->nintbk.ninf );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
        if( _info->nintbk.ipr >= 5 )
                lsgrg_msg(IOINFO, " ****EXITING SUBROUTINE: PH1OBJ\n" );


#endif

#ifdef sinf0
   #undef sinf0
#endif
#ifdef true0
   #undef true0
#endif
        return;

        /*     End of PH1OBJ
         * */
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS lsquad(LsgrgInfo *_info,long int ibv[], double x[], double v[], double xb1[],
         double xb2[], double xb3[])
{
long int i, i_, ii;
double aa, t2, t3, t32, ta, w1, w2, w3;

        if( _info->quadbk.icon == 2 )
                goto L_20;
        aa = powi(_info->bestbk.step/_info->quadbk.a2, 2);
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                ii = ibv[i_];
                x[ii - 1] = xb1[i_] + _info->bestbk.step*v[i_] + aa*(xb2[i_] - xb1[i_] -
                 _info->quadbk.a2*v[i_]);
        }
        goto L_40;
L_20:
        ;
        t2 = _info->quadbk.a2 - _info->quadbk.a1;
        t3 = _info->quadbk.a3 - _info->quadbk.a1;
        t32 = _info->quadbk.a3 - _info->quadbk.a2;
        ta = _info->bestbk.step - _info->quadbk.a1;
        aa = ta*ta;
        w1 = 1.0e0 - (ta*(t2 + t3) - aa)/(t2*t3);
        w2 = (ta*t3 - aa)/(t2*t32);
        w3 = (aa - ta*t2)/(t3*t32);
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                ii = ibv[i_];
                x[ii - 1] = w1*xb1[i_] + w2*xb2[i_] + w3*xb3[i_];
        }
L_40:
        return;

        /*     End of QUAD
         * */
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS tang(LsgrgInfo *_info,long int ibmap[], double zmem[], double grad[],
         long int ihag[], long int iheg[], long int inbv[], double v[],
         double d[], double rhs[])
{
long int i, i_, j, j_, jj;

        /************************************************************************
         *     Computes Tangent Vector V = - BINV * JACOBIAN * DIRECTION
         *         BINV is the basis inverse
         *         Jacobian has I,J element = Partial I'th binding constraint
         *         with respect to J'th nonbasic variable
         ************************************************************************
         * */

        /*      --------------------------
         *     | Compute RHS = GRAD * D
         *      --------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                rhs[i_] = 0.0e0;
        }

        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;
                jj = inbv[j_];
                xsaxpy(_info, grad, ihag, iheg, rhs, d[j_], _info->nintbk.nb, jj );
        }

        /*      --------------------------
         *     | Compute V = -BINV * RHS
         *      --------------------------
         * */
        xftran(_info, ibmap, zmem, rhs, _info->nintbk.nb );

        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                v[i_] = -rhs[i_];
        }
        v[_info->nintbk.nobj - 1] = _info->slpobj.slope;
        /*xxx      IF (MAXIM) V(NOBJ)=-V(NOBJ)
         * */
        return;

        /*     *** End of TANG ***
         * */
} /* end of function */



void /*FUNCTION*/ LSGRGCLASS redgra(LsgrgInfo *_info,double redgr[], long int ibmap[], double zmem[],
         double grad[], long int ihag[], long int iheg[], double u[],
         long int ibv[], long int inbv[], double g[], double alb[], double ub[],
         long int iobjpt[])
{
long int i, i_, ii, k, k_;
double bnd, f, gk, rdot;

        if( _info->nintbk.ipr >= 6 )
                lsgrg_msg(IOINFO, "\n ****ENTERING SUBROUTINE: REDGRA \n" );


        if( _info->nintbk.nb <= 0 ){
                lsgrg_error_msg(IOINFO, "ERROR...REDGRA: NB LE ZERO\n" );
         /*     xerror( 0, iounit.ioerr ); */
                lsgrg_errorexit(JUMPBUF,_LSGRG_REDGRA_NB_LE_0);
                return;
        }

        /*C     ----------------------------------------------------------
         *C    | Store the appropriate CB in U
         *     | If PHASE-1, compute Gradient of sum of infeasibilities
         *      ----------------------------------------------------------
         *
         *
         * */
        if( _info->nintbk.ninf == 0 ){
                for( k = 1; k <= _info->nintbk.nb; k++ ){
                        k_ = k - 1;
                        u[k_] = getcb(_info, ibv[k_], iobjpt, _info->dimen.n, _info->memory.lgrad,
                         grad );
                }
                u[_info->nintbk.nobj - 1] = 1.0e0;
                if( _info->optblk.maxim ){
                        for( k = 1; k <= _info->nintbk.nb; k++ ){
                                k_ = k - 1;
                                u[k_] = -u[k_];
                        }
                }
        } else{
                for( k = 1; k <= _info->dimen.mp1; k++ ){
                        k_ = k - 1;
                        gk = g[k_];
                        bnd = alb[_info->dimen.n + k_];
                        if( gk < (bnd - _info->limits.epnewt*(1.0e0 + fabs( bnd ))) )
                                goto L_20;
                        bnd = ub[_info->dimen.n + k_];
                        if( gk > (bnd + _info->limits.epnewt*(1.0e0 + fabs( bnd ))) )
                                goto L_30;
                        u[k_] = 0.0e0;
                        goto L_40;
L_20:
                        u[k_] = -1.0e0;
                        goto L_40;
L_30:
                        u[k_] = 1.0e0;
L_40:
                        ;
                }
                if( _info->ph1bk.phmult != 0.0e0 ){
                        if( _info->optblk.maxim ){
                                f = -_info->ph1bk.phmult;
                        } else{
                                f = _info->ph1bk.phmult;
                        }
                        for( i = 1; i <= _info->nintbk.nb; i++ ){
                                i_ = i - 1;
                                u[i_] = u[i_] + f*getcb(_info, ibv[i_], iobjpt, _info->dimen.n,
                                 _info->memory.lgrad, grad );
                        }
                }
        }

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 6 )
                {
                lsgrg_msg(IOINFO,  "   THE VECTOR CB IS: \n " );
                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.3e", u[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,  "\n" );
                }
#endif

        /*      ------------------------------------------------
         *     | Computes Lagrange Multipliers U = CB*BINV
         *      ------------------------------------------------
         * */
        xbtran(_info, ibmap, zmem, u, _info->nintbk.nb );

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 6 )
                {
                lsgrg_msg(IOINFO,  "   THE LAGRANGE MULTIPLIERS:\n " );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", u[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,  "\n" );
                }
#endif

        /*      ---------------------------------------------
         *     | Finally, compute REDGRA using multipliers.
         *     | Note:  Cnb is still stored in GRAD, so no need
         *     |        to access it with IOBJPT
         *      ---------------------------------------------
         * */
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                ii = inbv[i_];
                xdot(_info, grad, ihag, iheg, u, _info->nintbk.nb, ii, &rdot );
                redgr[i_] = -rdot;
        }

        if( _info->nintbk.ipr >= 6 )
            lsgrg_msg(IOINFO, " ****EXITING SUBROUTINE: REDGRA\n" );


        return;


        /*     *** End of REDGRA *** */
} /* end of function */


double /*FUNCTION*/ LSGRGCLASS getcb(LsgrgInfo *_info,long int ivar, long int iobjpt[], long int n,
         long int lgrad, double grad[])
{
double getcb_v;

        /********************************************************************
         *     This function returns the OBJ GRADIENT element for "IVAR"
         ********************************************************************
         * */



        /*     *  Start of GETCB * */
        if( ivar > n ){
                getcb_v = 0.0e0;
        } else if( iobjpt[ivar - 1] == 0 ){
                getcb_v = 0.0e0;
        } else{
                getcb_v = grad[iobjpt[ivar - 1] - 1];
        }
        return( getcb_v );

        /*     ** END OF GETCB */
} /* end of function */



void /*FUNCTION*/ LSGRGCLASS search(LsgrgInfo *_info,long int ibmap[], double zmem[], double grad[],
         long int ihag[], long int iheg[], double d[], double x[], double g[],
         long int ibv[], double v[], double xb1[], double xb2[], double xb3[],
         long int inbv[], double alb[], double ub[], double xstat[], double gbest[],
         double xbest[], long int ibc[], double rowb[], double colb[],
         double row[], double ascale[])
{
LOGICAL32 failp;
long int i, i_, ii, ipt, j, k, lmquit, maxcut, maxdub, ncut,/*  ncut_, */
     /*     ndub_, */ next, nsbbnd;
double bt, df, f2, f3, fa, fb, fc, fd, ftemp, objint, pctchg, sa,
         sb, sc, sd, si, smax_, steplg, stpmxt, t2, t3, tmp, /*tmp2, tmp3,*/
         ts, ts2, xi;

        /********************************************************************
         *     This subroutine performs the one-dimensional minimization.
         ********************************************************************
         * */

        ipt = _info->nintbk.ipr;
        if( _info->nintbk.ipr >= 5 )
            lsgrg_msg(IOINFO, "\n SEARCH ENTERED \n" );

        _info->misc.nsear = _info->misc.nsear + 1;
        /*      --------------------------------------------------
         *     | Initialization for quadratic extrapolation scheme
         *      -------------------------------------------------- */
        _info->quadbk.icon = 0;
        _info->quadbk.a1 = 0.0e0;
        if( _info->nintbk.nb == 0 )
                goto L_20;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                ii = ibv[i_];
                xb1[i_] = x[ii - 1];
        }
L_20:
        ;

        /*      ----------------
         *     | Initialization
         *      ---------------- */
        _info->newsrc.iter = 0;
        _info->initbk.lastcl = 1;
        _info->srchlg.fail = FALSE;
        _info->srchlg.uncon = TRUE;
        _info->srchlg.succes = TRUE;
        _info->srchlg.mxstep = FALSE;
        _info->srchlg.jstfes = FALSE;
        _info->supblk.basbnd = FALSE;
        _info->bestbk.step = 0.0e0;
        lmquit = _info->limits.itlim/2;
        if( _info->quadbk.iquad == 1 )
                lmquit = 4;
        maxcut = 10;
        maxdub = 30;
        sa = 0.0e0;
        sb = _info->bestbk.stpbst;
        fa = g[_info->nintbk.nobj - 1];
        objint = fa;
        ftemp = -_info->tols.plinfy;
        _info->bestbk.stpbst = 0.0e0;
        _info->redph.trubst = _info->bestbk.truobj;
        _info->bestbk.objbst = g[_info->nintbk.nobj - 1];
        _info->redser.ninfb = _info->nintbk.ninf;

        if( _info->nintbk.nb == 0 )
                goto L_40;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                xbest[i_] = x[j - 1];
        }
L_40:
        ;
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                gbest[i_] = g[i_];
        }

        smax_ = _info->tols.eps;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                ts = fabs( d[i_] );
                if( smax_ < ts )
                        smax_ = ts;
        }

        /*      --------------------------------------------------------------
         *     | Compute STEPMX = Max step until first super basic hits bound
         *     |         STEPLG = Max step until last super basic hits bound
         *      -------------------------------------------------------------- */
        _info->bestbk.stepmx = _info->tols.plinfy;
        steplg = -_info->tols.plinfy;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                si = d[i_];
                if( fabs( si ) < _info->tols.tolz )
                        goto L_90;
                j = inbv[i_];
                bt = ub[j - 1];
                if( si < 0.0e0 )
                        bt = alb[j - 1];
                if( fabs( bt ) > 1.0e20 )
                        bt = sign( 1.0e20, bt );
                ts = (bt - x[j - 1])/si;
                if( ts > _info->bestbk.stepmx )
                        goto L_85;
                _info->misc.jqq = i;
                _info->bestbk.stepmx = ts;
L_85:
                steplg = fmax( steplg, ts );
L_90:
                ;
        }
        stpmxt = _info->bestbk.stepmx;

        /*      ------------------------------------------
         *     | If projecting superbasics, reset STEPMX
         *      ------------------------------------------ */
        if( _info->optblk.multsb )
                _info->bestbk.stepmx = steplg;

        /*      -----------------------------
         *     | Determine initial step size
         *      ----------------------------- */
 /* L_95: unreferenced??? */

        ts = _info->epscom.eps3/smax_;
        if( sb < ts )
                sb = ts;
        pctchg = 0.05e0;
        ts = _info->bestbk.stepmx;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                j = inbv[i_];
                xi = x[j - 1];
                xstat[i_] = xi;
                xi = fabs( xi );
                si = fabs( d[i_] );
                if( si < _info->tols.tolz )
                        goto L_120;
                if( xi < 1.0e0 )
                        goto L_100;
                ts2 = pctchg*xi/si;
                goto L_110;
L_100:
                ts2 = pctchg/si;
L_110:
                if( ts > ts2 )
                        ts = ts2;
L_120:
                ;
        }

        /*      ------------------------------------------------------------
         *     | Set SB to TS unless previous search was unconstrained and
         *     | SB is smaller than TS.
         *      ------------------------------------------------------------ */
        if( !_info->srchlg.unconp || sb > ts )
                sb = ts;
        if( sb > _info->bestbk.stepmx )
                sb = _info->bestbk.stepmx;
        tmp = g[_info->nintbk.nobj - 1];
        if( _info->optblk.maxim && _info->nintbk.ninf == 0 )
                tmp = -tmp;
        /*     IF (IPR.GT.3) WRITE (IOOUT,310) TMP,SB,STEPMX */

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 4 )
                {
                sprintf(LSGRG_MSGBUFFER, " OBJECTIVE = %13.6e INITIAL STEP = %13.6e STEP TO BND = %13.6e\n",
                 tmp, sb, stpmxt );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

        /*      ----------------------------------------------------------
         *     | This loop computes FB and cuts back stepsize if FB > FA.
         *      ---------------------------------------------------------- */
        if( _info->nintbk.ninf == 0 && _info->sernew.edf <= 0.0e0 )
                _info->sernew.edf = fabs( _info->slpobj.slope )/2.0;
        for( ncut = 1; ncut <= maxcut; ncut++ ){
         /*        ncut_ = ncut - 1; unreferenced??? */
                failp = _info->srchlg.fail;
                _info->bestbk.step = sb;


                redobj(_info, ibmap, zmem, grad, ihag, iheg, x, g, ibv, v, xb1,
                 xb2, xb3, inbv, alb, ub, xstat, d, gbest, xbest, ibc, rowb,
                 colb, row, ascale );
                if( _info->glberr.abort )
                        return;
                if( _info->srchlg.fail )
                        goto L_125;


                fb = g[_info->nintbk.nobj - 1];
                tmp = fb;
                if( _info->optblk.maxim && _info->nintbk.ninf == 0 )
                        tmp = -tmp;
                if( _info->srchlg.jstfes )
                        goto L_200;
                if( fb <= fa + _info->tols.eps )
                        goto L_150;
                ftemp = fb;
                sc = _info->bestbk.step;
                goto L_130;
L_125:
                _info->counts.nstepc = _info->counts.nstepc + 1;

                /*          -------------------
                 *         | Reduce stepsize
                 *          ------------------- */
L_130:
                sb = _info->bestbk.step/(ipow(2, ncut));
        }
        goto L_240;

        /*      -------------------------------------------------------------
         *     | Step reduction phase completed--Have found a better point.
         *     | Consider all possible cases.
         *      -------------------------------------------------------------- */
L_150:
        ;
        if( _info->misc.lv != 0 )
                goto L_210;
        if( _info->srchlg.mxstep )
                goto L_220;
        if( failp )
                goto L_190;
        fc = ftemp;

        /*      ----------------------------------------
         *     | Begin quadratic interpolation block
         *      ---------------------------------------- */
        if( _info->quadbk.iquad == 0 || _info->nintbk.nb == 0 )
                goto L_2160;
        _info->quadbk.icon = 1;
        _info->quadbk.a2 = sb;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                ii = ibv[i_];
                xb2[i_] = x[ii - 1];
        }
        next = 3;
L_2160:
        ;

        /*      --------------------------------------------
         *     | Interpolate if a bracket has been found.
         *      -------------------------------------------- */
        if( fc > fb + _info->tols.eps )
                goto L_170;

        /*      -----------------------
         *     | Step doubling phase
         *      ----------------------- */
        for( _info->counts.ndub = 1; _info->counts.ndub <= maxdub; _info->counts.ndub++ ){
          /*       ndub_ = _info->counts.ndub - 1; unreferenced??? */

                /*          ---------------------------------------------
                 *         | Quit search if Newton failure anticipated.
                 *          --------------------------------------------- */
                if( _info->newsrc.iter >= lmquit )
                        goto L_230;
                sc = sb + sb;
                _info->bestbk.step = sc;


                redobj(_info, ibmap, zmem, grad, ihag, iheg, x, g, ibv, v, xb1,
                 xb2, xb3, inbv, alb, ub, xstat, d, gbest, xbest, ibc, rowb,
                 colb, row, ascale );
                if( _info->glberr.abort )
                        return;


                if( _info->srchlg.fail )
                        goto L_180;
                sc = _info->bestbk.step;
                fc = g[_info->nintbk.nobj - 1];
                tmp = fc;
                if( _info->optblk.maxim && _info->nintbk.ninf == 0 )
                        tmp = -tmp;
                if( _info->srchlg.jstfes )
                        goto L_200;

                /*      --------------------------------------
                 *     | Begin quadratic interpolation block
                 *      -------------------------------------- */
                if( _info->quadbk.iquad == 0 || _info->nintbk.nb == 0 )
                        goto L_2260;
                _info->quadbk.icon = 2;
                switch( next ){
                        case 1: goto L_2200;
                        case 2: goto L_2220;
                        case 3: goto L_2240;
                        }
L_2200:
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        ii = ibv[i_];
                        xb1[i_] = x[ii - 1];
                }
                _info->quadbk.a1 = sc;
                next = 2;
                goto L_2260;
L_2220:
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        ii = ibv[i_];
                        xb2[i_] = x[ii - 1];
                }
                _info->quadbk.a2 = sc;
                next = 3;
                goto L_2260;
L_2240:
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        ii = ibv[i_];
                        xb3[i_] = x[ii - 1];
                }
                _info->quadbk.a3 = sc;
                next = 1;
L_2260:
                ;

                /*      -----------------------------------------
                 *     | Interpolate if a bracket has been found.
                 *      ----------------------------------------- */
                if( fc > fb + _info->tols.eps )
                        goto L_170;
                if( _info->nintbk.ninf > 0 && fc >= (fb - _info->epscom.eps0) )
                        goto L_170;
                if( _info->misc.lv != 0 )
                        goto L_210;
                if( _info->srchlg.mxstep )
                        goto L_220;

                /*          ---------------------------------------
                 *         | Move 3 point pattern one step ahead.
                 *          --------------------------------------- */
                fa = fb;
                sa = sb;
                fb = fc;
                sb = sc;

        }
        _info->counts.ndub = maxdub;
        goto L_250;

        /*      -----------------------
         *     | Interpolation phase
         *      ----------------------- */
L_170:
        ;
        t2 = sb - sa;
        t3 = sc - sa;
        f2 = (fb - fa)*t3;
        f3 = (fc - fa)*t2;
        if( fabs( f2 - f3 ) < _info->tols.plzero )
                goto L_260;

        /*      -----------------------------------------
         *     | SD is minimum point for quadratic fit.
         *      ----------------------------------------- */
        sd = sa + 0.5e0*(t3*f2 - t2*f3)/(f2 - f3);
        if( sd <= sa || sd >= sc )
                goto L_260;

#ifdef IO_ENABLED
        if( _info->nintbk.ipr < 3 )
                goto L_178;
        lsgrg_msg(IOINFO,  "     QUADRATIC INTERPOLATION\n" );
        lsgrg_msg(IOINFO,  "     A          B          C            FA         FB         FC         D\n" );
        tmp = fa;
        tmp2 = fb;
        tmp3 = fc;
        if( !_info->optblk.maxim || _info->nintbk.ninf != 0 )
                goto L_175;
        tmp = -tmp;
        tmp2 = -tmp2;
        tmp3 = -tmp3;
L_175:
        sprintf(LSGRG_MSGBUFFER, "%11.4e%11.4e%11.4e %11.4e%11.4e%11.4e %11.4e\n",
         sa, sb, sc, tmp, tmp2, tmp3, sd );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
L_178:
        ;
#endif
        if( _info->nintbk.ipr > 3 )
            lsgrg_msg(IOINFO,"\n         QUADRATIC INTERPOLATION\n" );


        /*      ---------------------------------
         *     | Compute objective at SD point.
         *      --------------------------------- */
        _info->bestbk.step = sd;


        redobj(_info, ibmap, zmem, grad, ihag, iheg, x, g, ibv, v, xb1, xb2,
         xb3, inbv, alb, ub, xstat, d, gbest, xbest, ibc, rowb, colb,
         row, ascale );
        if( _info->glberr.abort )
                return;



        if( _info->srchlg.fail )
                goto L_180;
        fd = g[_info->nintbk.nobj - 1];
        if( _info->srchlg.jstfes )
                goto L_200;
        if( _info->misc.lv != 0 && fd < fb )
                goto L_210;
        if( fd > fb )
                _info->initbk.lastcl = 0;
        goto L_270;
L_180:
        ;

        /*      ----------------------------------------------------
         *     | Quit because Newton failed and better point found.
         *      ---------------------------------------------------- */
        if( _info->nintbk.ipr > 2 )
                lsgrg_msg(IOINFO, "\n         NEWTON FAILURE\n" );

        _info->initbk.lastcl = 0;
        goto L_270;
L_190:
        ;
        if( _info->nintbk.ipr > 2 )
            lsgrg_msg(IOINFO, "\n         EARLIER NEWTON FAILURE\n" );

        goto L_270;
L_200:
        ;

        /*      ----------------------------
         *     | Have just become FEASIBLE
         *      ---------------------------- */
        if( _info->nintbk.ipr > 1 )
           lsgrg_msg(IOINFO,
              "         ALL VIOLATED CONSTRAINTS SATISFIED.  "
              "NOW BEGIN TO OPTIMIZE TRUE OBJECTIVE \n" );

        _info->nintbk.ipr = ipt;
        return;
L_210:
        ;

        /*      ---------------------------
         *     | Basic variable hit bound
         *      --------------------------- */
        _info->srchlg.uncon = FALSE;
        _info->counts.nbs = _info->counts.nbs + 1;
        i = ibv[_info->misc.lv - 1];
        if( i <= _info->dimen.n )
                _info->supblk.basbnd = TRUE;
        if( _info->nintbk.ipr < 2 )
                goto L_270;
        k = i - _info->dimen.n;
#ifdef IO_ENABLED
        if( i <= _info->dimen.n )
                {
                sprintf(LSGRG_MSGBUFFER, "         BASIC VARIABLE #%4ld HIT BOUND\n",
                 i );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
        if( i > _info->dimen.n )
                {
                sprintf(LSGRG_MSGBUFFER, "         CONSTRAINT #%4ld NOW AT BOUND \n",
                 k );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        goto L_270;
L_220:
        ;

        /*      -------------------------------
         *     | Superbasic variable hit bound
         *      ------------------------------- */
        _info->srchlg.uncon = FALSE;
#ifdef IO_ENABLED
        if( (!_info->optblk.multsb) && (_info->nintbk.ipr > 1) )
                {
                sprintf(LSGRG_MSGBUFFER, "         SUPERBASIC VARIABLE #%4ld HIT BOUND\n",
                 inbv[_info->misc.jqq - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        goto L_270;
L_230:
        ;

        /*      ----------------------
         *     | Newton took too long
         *      ---------------------- */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 2 )
                {
                lsgrg_msg(IOINFO,  "         ANTICIPATED NEWTON FAILURE\n" );
                }
#endif
        goto L_270;
L_240:
        ;


        _info->srchlg.succes = FALSE;
        _info->initbk.lastcl = 0;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 1 )
                {
                sprintf(LSGRG_MSGBUFFER, "         NO OBJECTIVE IMPROVEMENT AFTER%4ld STEPSIZE REDUCTIONS\n",
                 maxcut );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        goto L_270;
L_250:
        ;

        /*      ------------------------------
         *     | Step size doubled NDUB times
         *      ------------------------------ */
        _info->srchlg.unbd = TRUE;
        goto L_270;
L_260:
        ;

        /*      --------------------------------------
         *     | Quadratic interpolation out of range
         *      -------------------------------------- */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 2 )
                {
                lsgrg_msg(IOINFO,  " QUADRATIC INTERPOLATION ABANDONED \n" );
                }
#endif
        goto L_270;
L_270:
        ;
        _info->srchlg.unconp = _info->srchlg.uncon;

        /*      -------------------------------------------
         *     | Pick up best point encountered and return.
         *      ------------------------------------------- */
        _info->bestbk.step = _info->bestbk.stpbst;
        _info->bestbk.stepmx = stpmxt;
        _info->bestbk.truobj = _info->redph.trubst;
        _info->nintbk.ninf = _info->redser.ninfb;

        /*      ---------------------------------------------------------
         *     | Compute the new superbasic values
         *     | If projecting superbasics, count the number of
         *     | superbasics which reached a bound during the search
         *      --------------------------------------------------------- */
        nsbbnd = 0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                j = inbv[i_];
                x[j - 1] = xstat[i_] + _info->bestbk.step*d[i_];
                if( !_info->optblk.multsb )
                        goto L_280;
                if( x[j - 1] >= ub[j - 1] - _info->limits.epboun*(1.0e0 + fabs( ub[j - 1] )) ){
                        if( x[j - 1] >= ub[j - 1] )
                                x[j - 1] = ub[j - 1];
                        nsbbnd = nsbbnd + 1;
                } else if( x[j - 1] <= alb[j - 1] + _info->limits.epboun*(1.0e0 +
                 fabs( alb[j - 1] )) ){
                        if( x[j - 1] <= alb[j - 1] )
                                x[j - 1] = alb[j - 1];
                        nsbbnd = nsbbnd + 1;
                }
L_280:
                ;
        }
#ifdef IO_ENABLED
        if( (_info->nintbk.ipr > 1 && nsbbnd != 0) && _info->optblk.multsb )
                {
                sprintf(LSGRG_MSGBUFFER, "         %4ld SUPERBASICS HIT BOUNDS\n",
                 nsbbnd );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

        if( _info->nintbk.nb == 0 )
                goto L_302;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                x[j - 1] = xbest[i_];
        }
L_302:
        ;
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                ts = gbest[i_];
                g[i_] = ts;
                x[_info->dimen.n + i_] = ts;
        }
        if( _info->nintbk.ninf == 0 ){
                df = fabs( g[_info->nintbk.nobj - 1] - objint );
#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 4 )
                        {
                        sprintf(LSGRG_MSGBUFFER, "   ESTIMATED DF  = %16.8e  ACTUAL DF = %16.8e\n",
                         _info->sernew.edf, df );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                _info->sernew.edf = 0.5000*df + 0.5000*_info->sernew.edf;
        }
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                lsgrg_msg(IOINFO,  " SEARCH COMPLETED \n" );
                }
#endif
        _info->nintbk.ipr = ipt;
        return;


        /*     End of SEARCH
         * */
} /* end of function */

#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"

#define IOINFO  &_info->io_info
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)
/***********************************************************************
 **                 FILE:        GRGNEWT  FORTRAN                      *
 **                 AUTHOR:      STUART SMITH                          *
 **                 LAST UPDATE:  Nov 1998                             *
 **                                                                    *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **                                                                    *
 **                                                                    *
 **      This file contains the following routines for LSGRG2, a       *
 **  large scale implementation of GRG2:                               *
 **                                                                    *
 **  ROUTINE           DESCRIPTION                        LAST UPDATE  *
 ** ---------         -------------                      ------------- *
 **  NEWTON     Performs Newton's method                  22 MAY 1996  *
 **  REDOBJ     Computes the Reduced Objective Function   23 OCT 1992  *
 **  CALFUN     Computes the problem function values      12 JUL 1993  *
 **  **change log**                                                    *
 **   10/98 jcp  fresh conversion from fortran.  cputime calls are     *
 **      replaced by calls to lsgrg_timer() (new shell routine)        *
 **      all timing calls wrapped around gcomp calls are commented out *
 **      since gcomp shell has timing calls there. timing calls for    *
 **      newton gcomp calls are left in since newton time exclusive    *
 **      of gcomp time is kept.                                        *
 **                                                                    *
 **                                                                    *
 ***********************************************************************
 *
 * */
void /*FUNCTION*/ LSGRGCLASS newton(LsgrgInfo *_info,long int ibmap[], double zmem[], double x[],
         long int inbv[], double g[], double xstat[], double d[], long int ibc[],
         long int ibv[], double rowb[], double colb[], double error[],
         LOGICAL32 *bviol, double alb[], double ub[], double ascale[])
{
long int i, i_, iii,  ip, isave, isv, j, j_,
         jr, k, k_, kcnt, lll, lvv, nbctmp;
double aaaa, delmax, delta, delta1, delta2, denom, ernorm, err, oldnrm,
         relerr, rr, tfinal, tgc, tin, tstart, tt, xlb, xnum, xub, xx;

/* replace statics with _info members */

   #define icheck _info->firstcall.newton

/*    static long icheck = 0; */

        /*.......................................................................
         *
         *     *** START OF ROUTINE: NEWTON ***
         *C */
        /* 4312    FORMAT( ' CORNER CHECK -- XNUM = ',D14.7,' DENOM = ',D14.7) */

        if( _info->nintbk.ipr >= 5 )
            lsgrg_msg(IOINFO," NEWTON ENTERED\n" );

        tgc = 0.0e0;
        tin = lsgrg_timer();
/*      cputime( &tin, &irc ); */

        /*      --------------------------------------------------------------
         *C    | FAIL is nonconvergence flag = FALSE CONVERGED =TRUE DIDNT
         *     | BVIOL is basic variable bound violation flag
         *      -------------------------------------------------------------- */
        _info->srchlg.fail = FALSE;
        *bviol = FALSE;
        oldnrm = _info->tols.plinfy;
        _info->counts.ncalls = _info->counts.ncalls + 1;


        /*C     --------------------------------------------------------
         *C    | If backing up to a constraint, we must make that
         *C    | constraint temporarily binding so that call to GCOMPX
         *C    | will evaluate it.  It's value is needed to set the slack
         *C    | during thel Newton iteration.
         *C     --------------------------------------------------------
         * */

        nbctmp = _info->bind.nbc;
        isave = 0;
        jr = 0;
        if( _info->misc.lv != 0 ){
                jr = ibv[_info->misc.lv - 1];
                if( jr <= _info->dimen.n ){
                        x[jr - 1] = _info->rednew.xb;
                } else if( _info->bind.nbc < _info->dimen.m ){
                        nbctmp = _info->bind.nbc + 1;
                        isave = ibc[nbctmp - 1];
                        /*@@@@@      IBC(NBCTMP) = JR - N */
                }
        }

        /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
        icheck = icheck + 1;
#ifdef IO_ENABLED
        if( icheck == 1 && _info->gmode.jdbg < 0 )
                {
                sprintf(LSGRG_MSGBUFFER,
                 "\n\n --- FIRST ENTRY NEWTON, HRDBND = %1c JDBG = %6ld\n",
                 TorF(_info->optblk.hrdbnd), _info->gmode.jdbg );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 6 )
                {
                lsgrg_msg(IOINFO,
                 "     BASIC VARIABLES IN BINDING CONSTRAINT ORDER:\n " );
                for( i = 1; i <= _info->bind.nbc; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", x[ibv[ibc[i - 1] -
                         1] - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
         }
#endif
        for( iii = 1; iii <= _info->limits.itlim; iii++ ){
                /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                 *XX   IF(.NOT.HRDBND) GO TO 7
                 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                 *   WILL DELETE THE NEXT FEW LINES AFTER DEBUG */

#ifdef IO_ENABLED
                if( _info->gmode.jdbg < 0 ){
                        lsgrg_msg(IOINFO,
                        " \n\n ------ ENTRY NEWTON BASIC VAR CHECK -------\n" );
                        for( i = 1; i <= _info->nintbk.nb; i++ ){
                                i_ = i - 1;
                                j = ibv[i_];
                                sprintf(LSGRG_MSGBUFFER, " BASIC VAR %4ld X(%3ld) = %14.7g\n",
                                 i, j, x[j - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                }
#endif
                /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
                lvv = 0;
                /*xxx  DELMAX = EPNEWT */
                delmax = _info->limits.epboun;
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        j = ibv[i_];
                        if( j > _info->dimen.n )
                                goto L_3;
                        delta1 = (alb[j - 1] - x[j - 1])/(1.0e0 + fabs( alb[j - 1] ));
                        delta2 = (x[j - 1] - ub[j - 1])/(1.0e0 + fabs( ub[j - 1] ));
                        delta1 = fmax( delta1, delta2 );
                        /*ZZZZZZZZZ--  DEBUG WRITE FOR DEVELOPMENT --- ZZZZZZZZZZZZZZZZ
                         *&&&&& IF(DELTA1.GT.EPBOUN.AND.JDBG.GT.0) */
#ifdef IO_ENABLED
                        if( delta1 > _info->limits.epboun && _info->gmode.jdbg < 0 )
                                {
                                sprintf(LSGRG_MSGBUFFER, " -- BASIC VAR %4ld X(%4ld) VIOLATES BOUND: \n X(J)  ALB(J), UB(J), VIOLATION = %14.7g %14.7g %14.7g %14.7g \n",
                                 i, j, x[j - 1], alb[j - 1], ub[j - 1], delta1 );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                        /*ZZZZZZZZZZZZ */
                        if( delta1 <= delmax )
                                goto L_3;
                        delmax = delta1;
                        lvv = i;
L_3:
                        ;
                }
                if( lvv == 0 )
                        goto L_7;
#ifdef IO_ENABLED
                if( _info->gmode.jdbg < 0 )
                        {
                        sprintf(LSGRG_MSGBUFFER,
                         "\nMAX BOUND VIOLATION: BASIC VAR %4ld BY %14.7g\n VAR = %5ld NEWTON ITERATION # = %5ld\n",
                         lvv, delmax, ibv[lvv - 1], iii - 1 );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                /*!!!! IF(.NOT.HRDBND) GO TO 7 */
                if( _info->optblk.hrdbnd ){
                        _info->misc.lv = lvv;
                        *bviol = TRUE;
                        goto L_140;
                        /*!!!!    RETURN */
                }
                /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
L_7:
                _info->newsrc.iter = iii - 1;
                if( _info->misc.lv == 0 )
                        goto L_20;

                /*      ------------------------------------------------
                 *     | Modify Superbasics if in BACKUP PHASE
                 *      ------------------------------------------------ */
                for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                        i_ = i - 1;
                        j = inbv[i_];
                        x[j - 1] = xstat[i_] + _info->rednew.xstep*d[i_];
                        if( !_info->optblk.multsb )
                                goto L_10;
                        if( x[j - 1] > ub[j - 1] ){
                                x[j - 1] = ub[j - 1];
                        } else if( x[j - 1] < alb[j - 1] ){
                                x[j - 1] = alb[j - 1];
                        }
L_10:
                        ;
                }
                /*CXXXXXXXXXXXXXXXXXXXXXXX */
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 6 )
                        {
                        lsgrg_msg(IOINFO,
                         "\n FOLLOWING MODIFICATION OF SUPERBASICS"
                         " (BACKUP PHASE) X = \n " );
                        for( j = 1; j <= _info->dimen.npmp1; j++ ){
                                sprintf(LSGRG_MSGBUFFER, "%13.6e", x[j - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                        lsgrg_msg(IOINFO, "\n" );
                        }
#endif
                /*CXXXXXXXXXXXXXXXXXXXXXXXX */

                /*          ---------------------------------------------------
                 *         | Evaluate (binding) constraints
                 *         | If no binding constraints, skip GCOMP call
                 *          --------------------------------------------------- */
L_20:
                ;
                if( _info->optblk.subset && nbctmp == 0 )
                        goto L_36;

                /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                 *     ERROR CHECK -- CHECK IF ANY VARIABLE ARE OUTSIDE THEIR BOUNDS
                 *                    WHEN HRDBND IS TRUE
                 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                 * */
                if( _info->optblk.hrdbnd ){
                        /*        DELMAX = EPBOUN */
                        delmax = _info->limits.epboun;
                        lll = 0;
                        aaaa = 0.1*_info->tols.plinfy;
                        /*!!!!    DO 333 I=1,NB */
                        for( i = 1; i <= _info->dimen.n; i++ ){
                      /*          i_ = i - 1;*/
                                /*!!!!       J = IBV(I) */
                                j = i;
                                if( j > _info->dimen.n )
                                        goto L_333;
                                xub = ub[j - 1];
                                xlb = alb[j - 1];
                                delta1 = (xlb - x[j - 1])/(1.0e0 + fabs( xlb ));
                                delta2 = (x[j - 1] - xub)/(1.0e0 + fabs( xub ));
                                /*!!!!       DELTA1 = (XLB - X(J))
                                 *!!!!       DELTA2 = (X(J)   - XUB) */
                                delta1 = fmax( delta1, delta2 );
                                if( delta1 <= delmax )
                                        goto L_333;
                                delmax = delta1;
                                lll = i;
                                k = j;
L_333:
                                ;
                        }
                        if( lll != 0 ){
                                xub = ub[k - 1];
                                xlb = alb[k - 1];
                                if( _info->scal.scaled ){
                                        xx = x[k - 1]*ascale[k - 1];
                                        if( xub < 0.1*_info->tols.plinfy )
                                                xub = xub*ascale[k - 1];
                                        if( xlb > -0.1*_info->tols.plinfy )
                                                xlb = xlb*ascale[k - 1];
                                }
#ifdef IO_ENABLED
                                sprintf( LSGRG_MSGBUFFER,
                                 "\n***NEWTON ERROR AT GCOMP CALL: HRDBND IS"
                                 " TRUE\nMAX BOUND VIOLATION:  BASIC VAR"
                                 " %4ld BY %14.7g\nVAR: %4ld X: %14.7e UB: %14.7e LB: %14.7e\n",
                                 lll, delmax, k, xx, xub, xlb );
                                lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                        }
                }
                /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                 *
                 *          --------------------------
                 *         | Unscale the variables
                 *          --------------------------
                 * */
                if( _info->scal.scaled ){
                        for( i = 1; i <= _info->dimen.n; i++ ){
                                i_ = i - 1;
                                x[i_] = x[i_]*ascale[i_];
                        }
                }

                /*          -------------------------
                 *         |  Evaluate the functions
                 *          -------------------------
                 * */
       /*       cputime( &tstart, &irc ); */
                tstart = lsgrg_timer();

                if( _info->optblk.subset ){
                        if( isave != 0 )
                                ibc[nbctmp - 1] = jr - _info->dimen.n;
                        lsgrg_gcompx(_info, g, x, &nbctmp, ibc );
                        if( isave != 0 )
                                ibc[nbctmp - 1] = isave;
                } else{
                        lsgrg_gcomp0( _info,g, x );
                        if( _info->optblk.maxim )
                                g[_info->nintbk.nobj - 1] = -g[_info->nintbk.nobj - 1];
                }
                if( _info->pardat.kderiv == 3 && _info->gmsblk.gmserr ){
                    lsgrg_error_msg(IOINFO,
                      " *** GAMS INTERP ERROR -  FORCING A NEWTON FAILURE\n" );
                        goto L_130;
                }
      /*        cputime( &tfinal, &irc ); */
                tfinal = lsgrg_timer();

                tt = tfinal - tstart;
                tgc = tgc + tt;
                _info->gctime.tgcomp = _info->gctime.tgcomp + tt;
                _info->counts.nftn = _info->counts.nftn + 1;


                /*          -------------------------------------------------------
                 *         | Rescale the variables and the true binding constraints
                 *          -------------------------------------------------------
                 * */
                if( _info->scal.scaled ){
                        for( i = 1; i <= _info->dimen.n; i++ ){
                                i_ = i - 1;
                                x[i_] = x[i_]/ascale[i_];
                        }
                        for( j = 1; j <= _info->bind.nbc; j++ ){
                                j_ = j - 1;
                                i = ibc[j_];
                                g[i - 1] = g[i - 1]*ascale[_info->dimen.n + i - 1];
                        }
                }

                /*          -------------------------------------------
                 *         | Set non-binding constraints ERROR to 0
                 *          -------------------------------------------
                 * */
L_36:
                if( _info->bind.nnbc == 0 )
                        goto L_47;
                for( i = 1; i <= _info->bind.nnbc; i++ ){
                        i_ = i - 1;
                        error[ibc[_info->bind.nbc + i_] - 1] = 0.0e0;
                }
                /*@@@@@    IF (ISAVE .NE. 0) ERROR(ISAVE) = 0.0D0
                 *
                 *      ---------------------------------------
                 *     | Compute NORM of equation ERROR
                 *      ---------------------------------------
                 * */
L_47:
                if( _info->misc.lv == 0 )
                        goto L_30;
                if( jr > _info->dimen.n ){
                        x[jr - 1] = g[jr - _info->dimen.n - 1];
                        if( _info->scal.scaled )
                                x[jr - 1] = x[jr - 1]*ascale[jr - 1];
                }
                xnum = _info->rednew.xb - x[jr - 1];
L_30:
                ;
                ernorm = 0.0e0;
                if( _info->misc.lv != 0 )
                        ernorm = fabs( xnum )/(1.0e0 + fabs( _info->rednew.xb ));
                if( _info->bind.nbc == 0 )
                        goto L_50;

                /*          -----------------------------------------
                 *         | Compute ERROR for binding constraints
                 *          -----------------------------------------
                 * */
                for( k = 1; k <= _info->bind.nbc; k++ ){
                        k_ = k - 1;
                        i = ibc[k_];
                        err = g[i - 1] - x[_info->dimen.n + i - 1];
                        relerr = fabs( err )/(1.0e0 + fabs( x[_info->dimen.n + i - 1] ));
                        if( relerr > ernorm )
                                ernorm = relerr;
                        error[i - 1] = err;
                }
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 6 )
                    {
                        lsgrg_msg(IOINFO,
                         " ERROR ARRAY FOR BINDING CONSTRAINTS:\n " );
                        for( i = 1; i <= _info->bind.nbc; i++ ){
                                sprintf(LSGRG_MSGBUFFER, "%13.6e", error[ibc[i - 1] -
                                 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                        lsgrg_msg(IOINFO, "\n" );
                    }
#endif

                /*      -----------------------------
                 *     | Test for convergence
                 *      -----------------------------
                 * */
L_50:
                ;
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " ERNORM = %13.6e\n", ernorm );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( ernorm < _info->limits.epnewt )
                        goto L_140;
                if( ernorm > oldnrm || _info->newsrc.iter == _info->limits.itlim )
                        goto L_130;
                ip = (_info->limits.itlim - _info->newsrc.iter)/2;


                /*          ------------------------------------------
                 *         | Test for predicted failure based on ERNORM
                 *          ------------------------------------------
                 * */
                if( _info->newsrc.iter > 0 ){
                        aaaa = ernorm*powi(ernorm/oldnrm, ip);
                        if( aaaa > _info->limits.epnewt ){
#ifdef IO_ENABLED
                                if( _info->nintbk.ipr >= 5 )
                                        {
                                        sprintf(LSGRG_MSGBUFFER, " PREDICTED NEWTON FAILURE -- ERNORM = %14.7e OLDNRM = %14.7e\n 0.5*ITNS TO GO = %4ld PRED ERRNORM = %14.7e\n",
                                         ernorm, oldnrm, ip, aaaa );
                                        lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                        }
#endif
                                goto L_130;
                        }
                }

                /*          --------------------------
                 *         | Proceed with Newton
                 *          -------------------------- */
                if( _info->misc.lv == 0 )
                        goto L_90;
                /*          ---------------------------------
                 *         | Computations for BACKUP phase
                 *          --------------------------------- */
                denom = -_info->rednew.corner;
                if( _info->nintbk.nb == 0 )
                        goto L_70;
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        rr = rowb[i_];
                        xnum = xnum + rr*error[i_];
                        denom = denom + rr*colb[i_];
                }
L_70:
                ;
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " XNUM = %15.8e DENOM = %15.8e\n",
                         xnum, denom );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( fabs( denom ) < _info->tols.tolz )
                        goto L_130;
                /*XXXX     IF (DABS(DENOM).LT.TOLZ) THEN
                 *XXXX         WRITE(IOOUT,4312) XNUM,DENOM
                 *XCC          GO TO 130
                 *XXX      ENDIF */
                delta = xnum/denom;
                _info->rednew.xstep = _info->rednew.xstep - delta;
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " DELTA = %15.8e XSTEP = %15.8e\n",
                         delta, _info->rednew.xstep );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( _info->nintbk.nb == 0 )
                        goto L_120;
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        error[i_] = error[i_] - delta*colb[i_];
                }
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 6 )
                   {
                        lsgrg_msg(IOINFO,
                         " ERROR ARRAY FOR BINDING CONSTRAINTS:\n " );
                        for( i = 1; i <= _info->nintbk.nb; i++ ){
                                sprintf(LSGRG_MSGBUFFER, "%13.6e", error[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                        lsgrg_msg(IOINFO, "\n" );
                    }
#endif

                /*C         ------------------------------
                 *C        | Compute Newton correction
                 *C         ------------------------------
                 * */
L_90:
                ;
                xftran(_info, ibmap, zmem, error, _info->nintbk.nb );

#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 6 )
                    {
                        lsgrg_msg(IOINFO,
                         " NEWTON CORRECTION FOR BINDING CONSTRAINTS:\n " );
                        for( i = 1; i <= _info->bind.nbc; i++ ){
                                sprintf(LSGRG_MSGBUFFER, "%13.6e", error[ibc[i - 1] -
                                 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                        lsgrg_msg(IOINFO, "\n" );
                    }
#endif

                /* chg 05/22/96 SHS - Force basic slacks to stay within their bounds.  This
                 *         asssures that when NEWTON converges, the binding constraints
                 *         will be no more than EPNEWT away from their bounds.
                 * */
                for( j = 1; j <= _info->bind.nbc; j++ ){
                        j_ = j - 1;
                        i = ibc[j_];
                        k = ibv[i - 1];
                        x[k - 1] = x[k - 1] - error[i - 1];
                        if( k > _info->dimen.n ){
                                if( x[k - 1] > ub[k - 1] ){
                                        x[k - 1] = ub[k - 1];
                                } else if( x[k - 1] < alb[k - 1] ){
                                        x[k - 1] = alb[k - 1];
                                }
                        }
                }
                if( _info->misc.lv != 0 && jr <= _info->dimen.n )
                        x[jr - 1] = _info->rednew.xb;
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 6 )
                    {
                     lsgrg_msg(IOINFO,
                      "     BASIC VARIABLES IN BINDING CONSTRAINT ORDER:\n " );
                     for( i = 1; i <= _info->bind.nbc; i++ ){
                             sprintf(LSGRG_MSGBUFFER, "%13.6e", x[ibv[ibc[i - 1] -
                              1] - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                     lsgrg_msg(IOINFO, "\n" );
                    }
#endif
                oldnrm = ernorm;
                _info->counts.nit = _info->counts.nit + 1;
L_120:
                ;
        }

        /*      ------------------------------
         *     | Failure --  No Convergence
         *      ------------------------------ */
L_130:
        _info->srchlg.fail = TRUE;
        _info->counts.nnfail = _info->counts.nnfail + 1;
        kcnt = 1;
        goto L_695;
L_140:
        ;
        kcnt = 0;

        /*      ----------------------------------------------------------
         *     | Evaluate and scale Non-Binding constraints at final point
         *      ----------------------------------------------------------
         * */
        if( _info->bind.nnbc == 0 )
                goto L_690;
        if( _info->optblk.subset ){

                /*          --------------------------
                 *         | Unscale the variables
                 *          --------------------------
                 * */
                if( _info->scal.scaled ){
                        for( i = 1; i <= _info->dimen.n; i++ ){
                                i_ = i - 1;
                                x[i_] = x[i_]*ascale[i_];
                        }
                }

                /*@@@@@    IF (ISAVE .NE. 0) IBC(NBC+1) = ISAVE */
        /*      cputime( &tstart, &irc ); */
                tstart = lsgrg_timer();

                lsgrg_gcompx( _info,g, x, &_info->bind.nnbc, &ibc[_info->bind.nbc] );
          /*    cputime( &tfinal, &irc ); */
                tfinal = lsgrg_timer();

                tt = tfinal - tstart;
                tgc = tgc + tt;
                _info->gctime.tgcomp = _info->gctime.tgcomp + tt;
                if( _info->optblk.maxim )
                        g[_info->nintbk.nobj - 1] = -g[_info->nintbk.nobj - 1];
                /*???????  NFTN=NFTN+1
                 *
                 *          -----------------------
                 *         | Rescale the variables
                 *          -----------------------
                 * */
                if( _info->scal.scaled ){
                        for( i = 1; i <= _info->dimen.n; i++ ){
                                i_ = i - 1;
                                x[i_] = x[i_]/ascale[i_];
                        }
                }
        }


        /*      -----------------------------------------------------
         *     | Rescale the  functions  and set the slack variables
         *     | of non-binding constraints.
         *      -----------------------------------------------------
         * */
        if( _info->scal.scaled ){
                for( j = 1; j <= _info->bind.nnbc; j++ ){
                        j_ = j - 1;
                        i = ibc[_info->bind.nbc + j_];
                        g[i - 1] = g[i - 1]*ascale[_info->dimen.n + i - 1];
                }
        }

        for( j = 1; j <= _info->bind.nnbc; j++ ){
                j_ = j - 1;
                i = ibc[_info->bind.nbc + j_];
                x[_info->dimen.n + i - 1] = g[i - 1];
        }
L_690:
        ;

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 6 )
            {
                lsgrg_msg(IOINFO, "\n     FINAL X = \n " );
                for( i = 1; i <= _info->dimen.npmp1; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", x[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
            }
        if( _info->nintbk.ipr >= 5 )
                {
                lsgrg_msg(IOINFO, " NEWTON COMPLETED\n" );
                }
#endif
L_695:
        ;
        /*@@@@ IF (ISAVE .NE. 0) IBC(NBC+1) = ISAVE
         *
         *      -------------------------
         *     | Update failure counts
         *      -------------------------
         * */
        _info->nwtcnt.inwtpt = _info->nwtcnt.inwtpt + 1;
        if( _info->nwtcnt.inwtpt > 20 )
                _info->nwtcnt.inwtpt = 1;
        if( _info->nwtcnt.inwtk < 20 ){
                _info->nwtcnt.inwtk = _info->nwtcnt.inwtk + 1;
                isv = 0;
        } else{
                isv = _info->nwtcnt.inwtlk[_info->nwtcnt.inwtpt - 1];
        }
        _info->nwtcnt.inwtfk = _info->nwtcnt.inwtfk + kcnt - isv;
        _info->nwtcnt.inwtlk[_info->nwtcnt.inwtpt - 1] = kcnt;

   /*   cputime( &tfinal, &irc ); */
        tfinal = lsgrg_timer();

        tt = tfinal - tin;
        _info->nwtim.tnewt = _info->nwtim.tnewt + tt;
        _info->nwtim.tnewtx = _info->nwtim.tnewtx + tt - tgc;
        return;

        /*C
         *C     END OF NEWTON
         *C */

  #undef icheck

} /* end of function */




void /*FUNCTION*/ LSGRGCLASS redobj(LsgrgInfo *_info,long int ibmap[], double zmem[], double grad[],
         long int ihag[], long int iheg[], double x[], double g[], long int ibv[],
         double v[], double xb1[], double xb2[], double xb3[], long int inbv[],
         double alb[], double ub[], double xstat[], double d[], double gbest[],
         double xbest[], long int ibc[], double rowb[], double colb[],
         double row[], double ascale[])
{
LOGICAL32 bviol;
long int i, i_, j, j_, jbnd, jr, k, kk, knb, l, lv1, lvv, nj,
         nviol;
double alpha, b1, b2, d1, d2, denom, est, gj, ratio, reltol, t, tb,
         thet, thmin, tmp, tol, tst, tstnum, xbestj, xbi,
         xdiff, xj;

/* replace statics with _info members */

      #define icallr _info->firstcall.redobj

/*    static long icallr = 0; */

        /*..................................................................
         *
         *     **** START OF ROUTINE: REDOBJ ****
         *
         *
         *C-----------------------------------------------------------------
         *C  NEW CODE FOR HARD BOUNDS ON BASIC VARAIBLES 12/86 JC PLUMMER
         *C
         *C  HRDBND = .TRUE. ==>  NEWTON WILL NEVER CALL GCOMP WITH BASIC
         *C  VAR WHICH VIOLATES ITS BOUND.  NEWTON WILL TERMINATE AND
         *C  SIGNAL REDOBJ TO INITIATE BACKUP PHASE ON THAT VAR
         *C  NEWTON WILL RETURN BVIOL = .TRUE., LV = INDEX OF BASIC
         *C  VAR WITH LARGEST BOUND VIOLATION (LV COMES THROUGH COMMON)
         *C------------------------------------------------------------------ */
        _info->gmode.jdbg = 1;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                lsgrg_msg(IOINFO, " REDOBJ ENTERED\n" );
                }
        if( icallr == 0 && _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " ... FIRST ENTRY REDOBJ, EPS = %10.3e PLINFY = %10.3e\n",
                 _info->tols.eps, _info->tols.plinfy );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( icallr == 0 )
                icallr = icallr + 1;

        _info->srchlg.fail = FALSE;
        _info->rednew.xb = 0.0e0;
        _info->misc.lv = 0;
        lv1 = 0;
        nviol = 0;
        _info->srchlg.mxstep = _info->bestbk.step >= _info->bestbk.stepmx;
        if( _info->srchlg.mxstep )
                _info->bestbk.step = _info->bestbk.stepmx;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                j = inbv[i_];
                x[j - 1] = xstat[i_] + _info->bestbk.step*d[i_];
                if( !_info->optblk.multsb )
                        goto L_10;
                if( x[j - 1] > ub[j - 1] ){
                        x[j - 1] = ub[j - 1];
                } else if( x[j - 1] < alb[j - 1] ){
                        x[j - 1] = alb[j - 1];
                }
L_10:
                ;
        }

        if( _info->nintbk.nb != 0 )
                goto L_20;
        /*     10     NO BINDING CONSTRAINTS
         *     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I */
        if( _info->scal.scaled ){
                for( i = 1; i <= _info->dimen.n; i++ ){
                        i_ = i - 1;
                        x[i_] = x[i_]*ascale[i_];
                }
        }

 /*     cputime( &tstart, &irc ); */
        lsgrg_gcomp0(_info, g, x );

/*      cputime( &tfinal, &irc );
        _info->gctime.tgcomp = _info->gctime.tgcomp + tfinal - tstart;
*/
        if( _info->optblk.maxim )
                g[_info->nintbk.nobj - 1] = -g[_info->nintbk.nobj - 1];
        _info->counts.nftn = _info->counts.nftn + 1;
        if( _info->scal.scaled ){
                for( i = 1; i <= _info->dimen.n; i++ ){
                        i_ = i - 1;
                        x[i_] = x[i_]/ascale[i_];
                }
                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        i_ = i - 1;
                        g[i_] = g[i_]*ascale[_info->dimen.n + i_];
                }
        }
        goto L_60;
        /*     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I */
L_20:
        ;
        if( _info->quadbk.icon == 0 )
                goto L_30;
        /*     OBTAIN INITIAL ESTIMATE OF BASIC VARIABLES BY QUADRATIC EXTRAPOLAT
         *     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I */
        lsquad(_info, ibv, x, v, xb1, xb2, xb3 );
        /*     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I */
        goto L_50;
L_30:
        ;
        /*----------------------------------------------------------------------
         *     LINEAR EXTRAPOLATION
         *     RETURN HERE IF BOUND VIOLATION DETECTED IN NEWTON TO REESTIMATE
         *     BASIC VARIABLES WITH REDUCED STEP
         *---------------------------------------------------------------------- */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                k = ibv[i_];
                x[k - 1] = xbest[i_] + (_info->bestbk.step - _info->bestbk.stpbst)*v[i_];
        }
L_50:
        ;
        /*     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I
         *
         *     Make sure basic variables are between bounds
         * CHG 5/24/96 SHS - Make sure basic slacks are always between bounds,
         *                   regardless of HRDBND setting.
         *
         * ------ New basic bound logic is below ------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                if( (j > _info->dimen.n) || _info->optblk.hrdbnd ){
                        if( x[j - 1] < alb[j - 1] ){
                                x[j - 1] = alb[j - 1];
                        } else if( x[j - 1] > ub[j - 1] ){
                                x[j - 1] = ub[j - 1];
                        }
                }
        }
        /*  --- End of New logic - Old logic is commented out below -----
         *XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         *xxx  IF (HRDBND) THEN
         *XXX     DO 44 I=1,NB
         *XXX        J = IBV(I)
         *XXX        IF ( J .GT. N) GO TO 44
         *XXX        IF (X(J) .LT. ALB(J)) THEN
         *XXX            X(J) = ALB(J)
         *XXX        ELSEIF (X(J) .GT. UB(J)) THEN
         *XXX            X(J) = UB(J)
         *XXX        ENDIF
         *X44    CONTINUE
         *xxx  ENDIF
         *XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         * -----  End of old logic -----------------------------------
         *
         *     BVIOL = .FALSE.
         * */
        newton(_info, ibmap, zmem, x, inbv, g, xstat, d, ibc, ibv, rowb, colb,
         row, &bviol, alb, ub, ascale );


        if( !_info->optblk.hrdbnd )
                goto L_58;
#ifdef IO_ENABLED
        if( _info->gmode.jdbg < 0 ){
                lsgrg_msg(IOINFO,
                 "\n ------ RETURN FROM NEWTON -- BASIC VARS -----\n" );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        j = ibv[i_];
                        sprintf(LSGRG_MSGBUFFER, " BASIC VAR %4ld X(%3ld) = %14.7g\n",
                         i, j, x[j - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
        }
#endif
        /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
        if( !bviol )
                goto L_58;
        nviol = nviol + 1;
        if( nviol > 1 )
                goto L_57;
        /*C
         *C  BOUND VIOLATION IN NEWTON -- REESTIMATE BASIC VARS FOR
         *C  BACKUP PHASE REENTRY.  RESTRICT NEW ESTIMATES TO
         *C  BE WITHIN VAR BOUNDS
         *C */
        alpha = _info->tols.plinfy/2.0;
        lvv = 0;
        tstnum = _info->tols.eps*_info->tols.plinfy;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                if( j > _info->dimen.n )
                        goto L_55;
                /*!!!!        TST = UB(J) + EPBOUN*(1.0D0+DABS(UB(J))) */
                tst = ub[j - 1];
                /*!!!!        IF(V(I).LT.0.0) TST = ALB(J)-EPBOUN*(1.0D0+DABS(ALB(J))) */
                if( v[i_] < 0.0 )
                        tst = alb[j - 1];
                ratio = _info->tols.plinfy;
                xdiff = fabs( tst - xbest[i_] );
                if( fabs( v[i_] ) > _info->tols.eps && xdiff <= tstnum )
                        ratio = (tst - xbest[i_])/v[i_];
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 6 )
                        {
                        sprintf(LSGRG_MSGBUFFER, ".. I = %4ld J = %4ld V(I) = %14.7g TST = %14.7g\n    XBEST(I) = %14.7g RATIO = %14.7g ALPHA = %14.7g\n",
                         i, j, v[i_], tst, xbest[i_], ratio, alpha );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( ratio >= alpha )
                        goto L_55;
                alpha = ratio;
                lvv = i;
L_55:
                ;
        }
        if( lvv != 0 )
                _info->misc.lv = lvv;
        xdiff = 0.01*(_info->bestbk.step - _info->bestbk.stpbst);
        if( alpha <= xdiff )
                alpha = xdiff;
        _info->rednew.xstep = alpha + _info->bestbk.stpbst;
        if( _info->rednew.xstep > (_info->bestbk.step - xdiff) )
                _info->rednew.xstep = _info->bestbk.step - xdiff;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 1 && lvv == 0 )
                {
                sprintf(LSGRG_MSGBUFFER,
                "\nRATIO TEST FOR BACKUP STEP FAILED, LV SET TO  BASIC VAR# %4ld X(%4ld)\n",
                 _info->misc.lv, ibv[_info->misc.lv - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
/* L_555: unreferenced ??? */
        ;
#ifdef IO_ENABLED
        if( _info->gmode.jdbg < 0 )
                {
                sprintf(LSGRG_MSGBUFFER, "\n-- BASIC VAR BOUND VIOLATION IN NEWTON, XSTEP  SET TO %14.7g LV = %4ld\n",
                 _info->rednew.xstep, _info->misc.lv );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                if( j > _info->dimen.n )
                        goto L_556;
                if( i == _info->misc.lv )
                        goto L_556;
                est = xbest[i_] + (_info->rednew.xstep - _info->bestbk.stpbst)*v[i_];
                if( est < alb[j - 1] )
                        est = alb[j - 1];
                if( est > ub[j - 1] )
                        est = ub[j - 1];
                x[j - 1] = est;
L_556:
                ;
        }
        j = ibv[_info->misc.lv - 1];
        if( x[j - 1] < alb[j - 1] )
                x[j - 1] = alb[j - 1];
        if( x[j - 1] > ub[j - 1] )
                x[j - 1] = ub[j - 1];
        /*CZZZZZZZZZ */
#ifdef IO_ENABLED
        if( _info->gmode.jdbg < 0 ){
                lsgrg_msg(IOINFO,
                 "\n ------ ECHO REDOBJ NEW BASI VAR ESTIMATES ----\n" );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        j = ibv[i_];
                        sprintf(LSGRG_MSGBUFFER, " BASIC VAR %4ld X(%3ld) = %14.7g\n",
                         i, j, x[j - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
        }
#endif
        /*CZZZZZZZZZ */
        goto L_167;
        /*C */
L_57:
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 1 )
                {
                lsgrg_msg(IOINFO,
                 "\nBASIC VARIABLE VIOLATED BOUND DURING REDOBJ\n"
                 "    BACKUP PHASE -- NEWTON FAILURE ASSUMED\n" );
                }
#endif

        _info->srchlg.fail = TRUE;
        _info->counts.nnfail = _info->counts.nnfail + 1;
        /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         *     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I */
L_58:
        if( _info->srchlg.fail )
                goto L_330;

L_60:
        knb = _info->nintbk.nb;
/* L_80:  unreferenced??? */
        if( knb == 0 ){
#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 2 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " STEP = %15.6e  OBJ = %15.6e  NEWTON ITERS = %4ld\n",
                         _info->bestbk.step, g[_info->nintbk.nobj - 1], _info->newsrc.iter );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                goto L_370;
        }
        if( _info->misc.lv == 0 )
                goto L_90;
        if( _info->rednew.xstep >= _info->bestbk.step || _info->rednew.xstep <= _info->bestbk.stpbst )
                goto L_320;
        _info->bestbk.step = _info->rednew.xstep;
L_90:
        if( _info->nintbk.ninf == 0 )
                goto L_120;

        /*     COME HERE IF INFEASIBLE
         *
         *
         *                      -BEGIN CHECKING FOR CHANGES IN STATUS OF BASIC
         *     VARIABLES
         * */
        _info->misc.lv = 0;
        thmin = _info->tols.plinfy;
        for( i = 1; i <= knb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                xj = x[j - 1];
                xbestj = xbest[i_];
                denom = xj - xbestj;
                if( fabs( denom ) <= _info->tols.eps )
                        goto L_110;
                tol = _info->limits.epboun;
                if( j > _info->dimen.n )
                        tol = _info->limits.epnewt;
                /*CXXX     IF ( J .GT. N) GO TO 110 */
                tb = alb[j - 1];
                reltol = tol*(1.0e0 + fabs( tb ));
                t = tb - reltol;
                jbnd = 0;
                if( xbestj > tb + reltol )
                        jbnd = 1;
                if( (xbestj >= t) && (xj < t) )
                        goto L_100;
                tb = ub[j - 1];
                reltol = tol*(1.0e0 + fabs( tb ));
                t = tb + reltol;
                jbnd = 0;
                if( xbestj < tb - reltol )
                        jbnd = 2;
                if( xbestj <= t && xj > t )
                        goto L_100;
                goto L_110;
L_100:
                if( jbnd != 0 )
                        t = tb;
                thet = (t - xbestj)/denom;
                if( thet >= thmin )
                        goto L_110;
                thmin = thet;
                _info->misc.lv = i;
L_110:
                ;
        }

        /*      TEST IF BACKUP PHASE NEEDED */
        if( _info->misc.lv != 0 )
                goto L_160;
        /*      NO BACKUP PHASE NEEDED  CHECK IF FEASIBLE
         *
         * */
        ph1obj(_info, ibc, g, ub, alb );


        if( _info->nintbk.ninf == 0 )
                goto L_340;
        /*    STILL INFEASIBLE */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 2 )
                {
                sprintf(LSGRG_MSGBUFFER, " STEP = %15.6e  OBJ = %15.6e  NEWTON ITERS = %4ld\n",
                 _info->bestbk.step, g[_info->nintbk.nobj - 1], _info->newsrc.iter );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        goto L_350;
        /*----------------------------------------------------------------------
         *     WE WERE FEASIBLE BEFORE NEWTON.  CHECK BOUNDS ON BASICS TO SEE IF
         *     WE ARE STILL FEASIBLE.
         *---------------------------------------------------------------------- */
L_120:
        _info->misc.lv = 0;
        thmin = _info->tols.plinfy;
        for( i = 1; i <= knb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                xj = x[j - 1];
                tol = _info->limits.epboun;
                if( j > _info->dimen.n )
                        tol = _info->limits.epnewt;
                xbestj = xbest[i_];
                denom = xj - xbestj;
                if( fabs( denom ) <= _info->tols.eps )
                        goto L_150;
                if( (alb[j - 1] - xj) <= (tol*(1.0e0 + fabs( alb[j - 1] ))) )
                        goto L_130;
                t = alb[j - 1] - xbestj;
                goto L_140;
L_130:
                if( (xj - ub[j - 1]) <= (tol*(1.0e0 + fabs( ub[j - 1] ))) )
                        goto L_150;
                t = ub[j - 1] - xbestj;
L_140:
                thet = t/denom;
                if( thet >= thmin )
                        goto L_150;
                thmin = thet;
                _info->misc.lv = i;
L_150:
                ;
        }
L_160:
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " LV = %5ld THMIN = %13.6e\n",
                 _info->misc.lv, thmin );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( thmin < 0.0e0 )
                goto L_320;
        if( _info->nintbk.ninf == 0 && _info->optblk.maxim )
                g[_info->nintbk.nobj - 1] = -g[_info->nintbk.nobj - 1];
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 2 )
                {
                sprintf(LSGRG_MSGBUFFER, " STEP = %15.6e  OBJ = %15.6e  NEWTON ITERS = %4ld\n",
                 _info->bestbk.step, g[_info->nintbk.nobj - 1], _info->newsrc.iter );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( _info->nintbk.ninf == 0 && _info->optblk.maxim )
                g[_info->nintbk.nobj - 1] = -g[_info->nintbk.nobj - 1];
        if( _info->misc.lv == 0 )
                goto L_350;
        if( _info->nintbk.ipr < 3 )
                goto L_167;
        i = ibv[_info->misc.lv - 1];
#ifdef IO_ENABLED
        if( i > _info->dimen.n ){
                i = i - _info->dimen.n;
                sprintf(LSGRG_MSGBUFFER,
                 "                              CONSTRAINT # %4ld VIOLATED BOUND\n",
                 i );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        } else{
                sprintf(LSGRG_MSGBUFFER,
                 "                              BASIC VARIABLE  X(%4ld)  VIOLATED BOUND\n",
                 i );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif
L_167:
        ;
        /*         **************************************************************
         *     BOUND VIOLATED--START BACK UP PHASE
         *         **************************************************************
         *
         *     SET XB = BOUND NEAREST X(LV)
         * */
        jr = ibv[_info->misc.lv - 1];
        b1 = alb[jr - 1];
        b2 = ub[jr - 1];
        d1 = x[jr - 1] - b1;
        d2 = x[jr - 1] - b2;
        _info->rednew.xb = b1;
        if( fabs( d1 ) > fabs( d2 ) )
                _info->rednew.xb = b2;
        if( !bviol )
                _info->rednew.xstep = _info->bestbk.stpbst + thmin*(_info->bestbk.step - _info->bestbk.stpbst);
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 2 )
                {
                sprintf(LSGRG_MSGBUFFER, " REDOBJ BACKING UP.  OBJ = %15.8e LV = %5ld XSTEP = %15.8e\n",
                 g[_info->nintbk.nobj - 1], _info->misc.lv, _info->rednew.xstep );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( _info->nintbk.nb == 0 )
                goto L_210;

        /*      --------------------------------------------------------
         *     | Compute COLB=B2*D, column part of Jacobian border
         *      --------------------------------------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                colb[i_] = 0.0e0;
        }

        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;
                kk = inbv[j_];
                xsaxpy(_info, grad, ihag, iheg, colb, d[j_], _info->nintbk.nb, kk );
        }

        /*      ------------------------------------------------------
         *     | Do row border and corner element calculations
         *      ------------------------------------------------------
         * */
        if( _info->misc.lv > _info->nintbk.nb )
                goto L_210;

        /*      -----------------------------------------------
         *     | Case of basic variable violating bound
         *      -----------------------------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                row[i_] = 0.0e0;
        }

        row[_info->misc.lv - 1] = 1.0e0;
        _info->rednew.corner = 0.0e0;
        goto L_270;

        /*      ---------------------------------------
         *     | Case of constraint violating bound
         *     | SHOULD NOT HAPPEN !!!!
         *      ---------------------------------------
         * */
L_210:
        ;
        j = ibv[_info->misc.lv - 1];
        l = j - _info->dimen.n;
#ifdef IO_ENABLED
        sprintf( LSGRG_MSGBUFFER,
        "\nREDOBJ...ERROR: CONSTRAINT # %6ld VIOLATED A BOUND BASIC VAR"
         " = %7ld\n",
         l, _info->misc.lv );
        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
        lsgrg_errorexit(JUMPBUF,_LSGRG_REDOBJ_CONSTVIOL);
        return;

        if( _info->nintbk.nb == 0 )
                goto L_240;

        /*       --------------------------------
         *      | Get basic row L of Jacobian
         *       --------------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                rowb[i_] = 0.0e0;
        }
        rowb[l - 1] = 1.0e0;

        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                k = ibv[i_];
                xdot(_info, grad, ihag, iheg, rowb, _info->nintbk.nb, k, &tmp );
                row[i_] = tmp;
        }
L_240:
        _info->rednew.corner = 0.0e0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                j = inbv[i_];
                xdot(_info, grad, ihag, iheg, rowb, _info->nintbk.nb, j, &tmp );
                _info->rednew.corner = _info->rednew.corner + tmp*d[i_];
        }

        /*      --------------------------------
         *     | Compute ROWB = ROW * BINV
         *      --------------------------------
         * */
        if( _info->nintbk.nb == 0 )
                goto L_310;

L_270:
        ;


        xbtran(_info, ibmap, zmem, row, _info->nintbk.nb );

        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                rowb[i_] = row[i_];
        }

        /*C????????????????????????????????????????????????????????????????? */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                lsgrg_msg(IOINFO, " ROW BORDER FOLLOWED BY CORNER = \n " );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%15.8e", rowb[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                sprintf(LSGRG_MSGBUFFER, "%15.8e\n", _info->rednew.corner );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
        if( bviol )
                goto L_310;
        /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

        /*      ----------------------------------------
         *     | Compute estimates of basic variables
         *      ---------------------------------------- */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                xbi = xbest[i_];
                x[j - 1] = xbi + (x[j - 1] - xbi)*thmin;
        }
L_310:
        lv1 = _info->misc.lv;
        goto L_50;

        /*     BACKUP FAILED
         * */
L_320:
        ;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 1 )
                {
                lsgrg_msg(IOINFO,
                 " REDOBJ CANNOT BACKUP TO FIRST VIOLATED CONSTRAINT\n" );
                }
#endif

        /*      ----------------------
         *     | Newton call FAILED
         *      ---------------------- */
L_330:
        ;
        _info->srchlg.fail = TRUE;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 1 )
                {
                sprintf(LSGRG_MSGBUFFER, " NEWTON FAILED TO CONVERGE.  NO.ITER. = %4ld STEP = %12.5e\n",
                 _info->newsrc.iter, _info->bestbk.step );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        goto L_440;

        /*      ----------------------
         *     | Just became feasible
         *      ---------------------- */
L_340:
        ;
        _info->srchlg.jstfes = TRUE;
        g[_info->nintbk.nobj - 1] = _info->bestbk.truobj;
        goto L_440;
        /*      ----------------------
         *     | Normal termination
         *      ---------------------- */
L_350:
        ;
        if( _info->misc.lv > 0 || _info->bind.nnbc == 0 )
                goto L_370;

        /*      ---------------------------------------
         *     | Check for new binding constraints
         *      --------------------------------------- */
        for( i = 1; i <= _info->bind.nnbc; i++ ){
                i_ = i - 1;
                j = ibc[_info->bind.nbc + i_];
                gj = g[j - 1];
                nj = _info->dimen.n + j;
                if( (fabs( gj - alb[nj - 1] ) < (_info->limits.epnewt*(1.0e0 + fabs( alb[nj - 1] )))) ||
                 (fabs( gj - ub[nj - 1] ) < (_info->limits.epnewt*(1.0e0 + fabs( ub[nj - 1] )))) )
                        lv1 = j;
        }

        /*      ---------------------
         *     | Update best values
         *      --------------------- */
L_370:
        ;
        if( g[_info->nintbk.nobj - 1] > _info->bestbk.objbst )
                goto L_430;
        if( _info->nintbk.nb == 0 )
                goto L_390;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                xbest[i_] = x[j - 1];
        }

L_390:
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                gbest[i_] = g[i_];
        }
        _info->bestbk.stpbst = _info->bestbk.step;
        _info->bestbk.objbst = g[_info->nintbk.nobj - 1];
        _info->redph.trubst = _info->bestbk.truobj;
        _info->redser.ninfb = _info->nintbk.ninf;
L_430:
        _info->misc.lv = lv1;
L_440:
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                lsgrg_msg(IOINFO, " REDOBJ COMPLETED\n" );
                }
#endif
        return;


        /*     End of REDOBJ
         * */

  #undef icallr

} /* end of function */


void /*FUNCTION*/ LSGRGCLASS calfun(LsgrgInfo *_info,double g[], double x[], double ascale[])
{
long int i, i_;

        /*.......................................................................
         *
         *     *** START OF ROUTINE: CALFUN ***
         *
         *
         *
         *      --------------------------
         *     | Unscale the variables
         *      --------------------------
         * */
        if( _info->scal.scaled ){
                for( i = 1; i <= _info->dimen.n; i++ ){
                        i_ = i - 1;
                        x[i_] = x[i_]*ascale[i_];
                }
        }

        /*      -------------------------
         *     |  Evaluate the functions
         *      ------------------------- */


  /*    cputime( &tstart, &irc ); */

        lsgrg_gcomp0(_info, g, x );
        if( _info->optblk.maxim )
                g[_info->nintbk.nobj - 1] = -g[_info->nintbk.nobj - 1];

/*      cputime( &tfinal, &irc );
        _info->gctime.tgcomp = _info->gctime.tgcomp + tfinal - tstart;
*/
        _info->counts.nftn = _info->counts.nftn + 1;


        /*      -------------------------------------------------
         *     | Rescale the variables and the constraints
         *      -------------------------------------------------
         * */
        if( _info->scal.scaled ){
                for( i = 1; i <= _info->dimen.n; i++ ){
                        i_ = i - 1;
                        x[i_] = x[i_]/ascale[i_];
                }
                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        i_ = i - 1;
                        g[i_] = g[i_]*ascale[_info->dimen.n + i_];
                }
        }
        return;
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS chkfun(LsgrgInfo *_info,double g[], double gg[], double x[], double ascale[],
         char *msg)
{
//long int i, i_;

        /*.......................................................................
         *
         *     *** START OF ROUTINE: CHKFUN ***
         *
         *      -------------------------
         *     |  EVALUATE THE FUNCTIONS
         *      -------------------------
         * */
        calfun(_info, gg, x, ascale );

#ifdef IO_ENABLED

        sprintf(LSGRG_MSGBUFFER, "\n  CHECKING G IN ....%s\n", msg );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                sprintf(LSGRG_MSGBUFFER, " %3ld:  G: %15.7e  GG: %15.7e  ERR: %15.7e\n",
                 i, g[i_], gg[i_], fabs( gg[i_] - g[i_] ) );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif
        return;
} /* end of function */

#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#define LSGRGCLASS

/*
#include "lsgrgver.h"
#ifdef LSGRGCPP
    #include "lsgrgc3.h"
#endif
*/

/*
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
*/
#ifndef LSGRGCPP
   #include <stddef.h>
   #include <stdlib.h>
   #include <math.h>


#define         fmax(a,b)       (double)( (a) > (b) ? (a) : (b) )

#define          max(a,b)               ( (a) > (b) ? (a) : (b) ) 
/*-------------------------------------------------------------------*/
/*  NOTE:  the invert subsystem is treated as a separate module      */
/*     from the rest of lsgrgc and does not include the system .h    */
/*     file.  the prototypes are included here for internal calls    */
/*  ** 6/99 jcp ** condition prototypes and includes on whether      */
/*     this is the c++ version                                       */
/*-------------------------------------------------------------------*/
void lsdetm(float*,long*,long,long*,double[],long);
void lslubt(long,long,double[],short*,long*,short*,double[],
        double[],long[]);
void lslucc(long,long,short[],long[],short[],long*,long*);
void lslucr(long,long,short[],double[],long[],short[],long*,
        long*);
void lslufc(double[],short*,long,long,long,long*,short*,
        double[],long*,long,long[],double[],long*);
void lsluft(long,long,double[],short*,long*,short*,double[],
        double[],long*,long*,long[]);
void lslupv(long,long,short*,double[],long*,short*,double[],
        long,long[],double[],long*);
void lsortr(long,long,double[],short[],long[],short[],short[]);

#endif
/***********************************************************************
 **                 FILE:        GRGINVRT  FORTRAN                     *
 **                 AUTHOR:      LINUS SCHRAGE                         *
 **                 LAST UPDATE: 05 NOV   1990                         *
 **                                                                    *
 **      This file contains the following routines for LSGRG2, a       *
 **  large scale implementation of GRG2. These routines create a       *
 **  sparse LU factorization of the basis,and either update or use     *
 **  it in some way.                                                    *
 **                                                                    *
 **        **************************************************          *
 **        **************************************************          *
 **        ***               LSGRG2                       ***          *
 **        ***          COPYRIGHT   1991                  ***          *
 **        ***            LEON S. LASDON                  ***          *
 **        ***                 &                          ***          *
 **        ***            ALLAN D. WAREN                  ***          *
 **        ***                 &                          ***          *
 **        ***            STUART SMITH                    ***          *
 **        **************************************************          *
 **        **************************************************          *
 **                                                                    *
 **                                                                    *
 **  ROUTINE           DESCRIPTION                        LAST UPDATE  *
 ** ---------         -------------                      ------------- *
 **  LSLUFC       Creates L and U factors                 05 NOV 1990  *
 **  LSLUFT       Does an FTRAN using these factors       05 NOV 1990  *
 **  LSLUBT       Does a BTRAN                            05 NOV 1990  *
 **  LSLUCC       Compresses file of positive integers    05 NOV 1990  *
 **  LSLUCR       Compresses file of pos.int.and reals    05 NOV 1990  *
 **  LSORTR       Construct and sort a row list starting  05 NOV 1990  *
 **               with (I,J,A) triplets                                *
 **  LSLUPV       Does a pivot                            05 NOV 1990  *
 **  LSDETM       Computes determinant of inverse         05 NOV 1990  *
 **                                                                    *
 ***********************************************************************
 * */
void /*FUNCTION*/ LSGRGCLASS lslufc(double a[], short int *ind, long int nzero,
         long int nzdim, long int n, long int *ip, short int *lnrc, double w[],
         long int *kdropd, long int keepsq, long int iluprm[], double dluprm[],
         long int *kalinv)
{
#define IND(I_,J_)      (*(ind+(I_)*(nzdim)+(J_)))
#define IP(I_,J_)       (*(ip+(I_)*(n)+(J_)))
#define LNRC(I_,J_)     (*(lnrc+(I_)*(n)+(J_)))
long int i, i_, idum, idum_, ii, ii_, il, in, ip0, ipu, ipu_, ipv,
         ipv1, ipv_, ipvt, ir, ir_, j, j_, jpvt, jr, k, k1, k2, k_, kc,
         kc_, kci, kcl, kcost, kj, kk, kk_, kl, klc, kn, knp, knp_, kp,
         kp1, kp2, kpc, kpr2, kpri, kprl, kr, kri, krl, ks, ks_, l, l_,
         lcol, lnrc4, lntmp, lnu, lowstl, lrow, lstinl, mkost, mt[2],
         mxcmp, nc, nc_, ndrpd[2], nfill, nfillr, nonstk, nps, nspike,
         ntcmp, nvinl, nzc, nzcr, nzcr_, nzinl, nzinu, nzr, nzs;
float dak;
double am, amax, au, auab, dpiv, g, gromx, u, zerotl;
static long maxps = 4;


        /*   INVERT A SPARSE MATRIX INTO FACTORED LU FORM
         *
         *  INPUTS:
         *    N = NO OF ROWS AND COLS IN MATRIX TO BE INVERTED
         *    NZERO = NO. OF NONZEROES IN THE MATRIX
         *    NZDIM = LENGTH OF THE A() AND IND(,) VECTORS
         *    A() = VECTOR STORING THE NONZEROES
         *    IND( K, 1) = ROW IN WHICH KTH NONZERO APPEARS
         *    IND( K, 2) = COL IN WHICH KTH NONZERO APPEARS
         *      THIS SUBROUTINE WANTS A(), IND(,) TO BE IN COLUMN ORDER,
         *      IF THIS IS INCONVENIENT, INSERT AN EXTRA 0.0 ANYWHERE IN
         *        THE MATRIX AND IT WILL BE PROCESSED CORRECTLY.
         *    W() = A WORK VECTOR
         *    KEEPSQ = 1 IMPLIES TRY TO USE EXISTING PIVOT SEQUENCE
         *               SUPPLIED IN LNRC( , 3/4), E.G., FROM PREVIOUS INVERT
         *           = 0 IF FIND NEW GOOD SEQUENCE;
         *           KEEPSQ = 1 IS CURRENTLY DISREGARDED, I.E., NOT IMPLEMENTED;
         *    DLUPRM( 1) = RELATIVE ROW SIZE THRESHOLD FOR BEING A PIVOT
         *                 E.G., 0.1
         *    DLUPRM( 2) = TOLERANCE FOR DISREGARDING A NUMBER CLOSE TO ZERO
         *                 E.G., .1E-15
         *
         *  OUTPUTS:
         *    KDROPD = NO. ROWS/COLS DROPPED, THESE APPEAR AT BOTTOM OF U
         *    DLUPRM( 3) = LARGEST COEF IN ABS VALUE FOUND IN A() ON INPUT
         *    DLUPRM( 4) = LARGEST COEF IN ABS VALUE IN INVERSE; IF THIS
         *                 IS MUCH LARGER THAN DLUPRM( 3), E.G., BY A FACTOR
         *                 OF 10.E15, YOU MAY WISH TO REINVERT WITH A
         *                 LARGER DLUPRM( 1).
         *    A() = ELEMENTS OF U AND L IN FACTORED FORM(L/U)
         *    IP( I, 1/2) = START OF ROW/COL I IN U.
         *    LNRC( I, 1/2) = NUMBER OF NON-ZEROS IN ROW/COL I OF U.
         *    LNRC( I, 3/4) = ROW/COL PIVOTED IN POSITION I, FROM TOP OF U.
         *    LNRC( I, 5) = POSITION IN U OF ROW I
         *    IND( K, 1) = ROW IN WHICH KTH NONZERO APPEARS IN COL LISTS
         *    IND( K, 2) = COL IN WHICH KTH NONZERO APPEARS IN ROW LISTS
         *    ILUPRM( 7) = NO. ELEMENTS IN U.
         *    ILUPRM( 6) = NO. ELEMENTS IN L, NOTE STORED DOWNWARD FROM NZDIM
         *    ILUPRM( 9) = NO. OF ELEMENTS IN ROW FILE OF U, INCLUDING EMPTIES
         *    ILUPRM( 8) = NO. OF ELEMENTS IN COL FILE OF U, INCLUDING EMPTIES
         *    ILUPRM( 14) = NO. SINGLETON COLS IN BASIS
         *    ILUPRM( 15) = NO. TRIANGULAR COLS IN BASIS
         *    ILUPRM( 16) = NO. SPIKES IN INVERSE
         *    ILUPRM( 17) = NO. COLS IN L
         *    ILUPRM( 18) = NO. TRIANGULAR ROWS IN INVERSE
         *    KALINV = 0 IF NO ERROR OCCURRED ,
         *             9 IF RAN OUT OF SPACE
         *             10 IF A COL HAD > 1 ENTRY IN A ROW,
         *             11 IF MATRIX HAD < 1 ROWS.
         *
         *
         *  REFERENCES:
         *   ALGORITHMS:
         *    SUHL, U. H. AND L. AITTONIEMI(1987), "COMPUTING SPARSE
         *    LU-FACTORIZATIONS FOR LARGE-SCALE LINEAR PROGRAMMING BASES",
         *    FREIE UNIVERSITAET BERLIN.
         *
         *    REID, J. K.(1982), "A SPARSITY-EXPLOITING VARIANT OF THE BARTELS-
         *    GOLUB DECOMPOSITION FOR LINEAR PROGRAMMING BASES", MATH. PROG.,
         *    VOL. 24, PP. 55-69.
         *
         *    REID, J. K,(1976), "FORTRAN SUBROUTINES FOR HANDLING SPARSE LINEAR
         *    PROGRAMMING BASES", REPORT AERE-R.8269, HARWELL.
         *
         *    MARKOWITZ, H. M.(1957), "THE ELIMINATION FORM OF THE INVERSE AND
         *    ITS APPLICATION TO LINEAR PROGRAMMING", MAN. SCI., VOL. 3,
         *    PP. 255-269
         *
         *   DATA STRUCTURES:
         *    GUSTAVSON, F. G.(1972), "SOME BASIC TECHNIQUES FOR SOLVING SPARSE
         *    SYSTEMS OF LINEAR EQUATIONS" IN SPARSE MATRICES AND THEIR
         *    APPLICATIONS ED. BY D. J. ROSE AND R. A. WILLOUGHBY, PLENUM PRESS,
         *    NY, PP. 41-52.
         *
         *    DUFF, I. S., AND J. K. REID(1979), "SOME DESIGN FEATURES OF A
         *    SPARSE MATRIX CODE", TRANS. MATH. SOFT., VOL. 5, NO. 1, PP. 18-35.
         * */
        /*  LIST POINTERS FOR DROPPED ROWS/COLS, AND NO. DROPPPED */
        /*  MAXPS = SUHL'S P PARAMETER FOR MAX MARKOWITZ PASSES ALLOWED */

        /*  DATE 25 MAY 1990
         *
         *  GET OUT IF INPUT BLATANTLY BAD */
        if( n < 1 )
                goto L_8800;
        *kalinv = 0;

        /*  MXCMP = MAXIMUM NUMBER  COMPRESSES BEFORE GIVING OUT OF SPACE
         *   ERROR RETURN */
        mxcmp = max( 25, n/8 );
        ntcmp = 0;

        /*   IF USER FORGOT TO SET DLUPRM( 1), SET TO REASONABLE VALUE */
        if( dluprm[0] > 1. )
                dluprm[0] = .0625;
        if( dluprm[0] < 0.0 )
                dluprm[0] = .0625;
        u = dluprm[0];

        /*  G = LARGEST VALUE ENCOUNTERED */
        g = 0.;

        /*  INITIALIZE LISTS OF DROPPED ROWS */
        mt[0] = 0;
        mt[1] = 0;
        /*  INVERSE STATS */
        iluprm[15] = 0;
        iluprm[16] = 0;
        /*  NO ROWS/COLS HAVE BEEN DROPPED AS YET */
        ndrpd[0] = 0;
        ndrpd[1] = 0;
        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                w[i_] = 0.0;
        }
        for( j = 1; j <= 6; j++ ){
                j_ = j - 1;
                for( i = 1; i <= n; i++ ){
                        i_ = i - 1;
                        LNRC(j_,i_) = 0;
                }
        }

        /*  COUNT NONZEROES IN EACH ROW AND COL, AFTER DELETING SMALL ONES */
        lnu = nzero;
        zerotl = dluprm[1];
        for( k = lnu; k >= 1; k-- ){
                k_ = k - 1;
                dak = (float)  fabs( a[k_] ); /* warning */
                if( dak <= zerotl )
                        goto L_1040;
                if( dak > g )
                        g = dak;
                i = IND(0,k_);
                LNRC(0,i - 1) = (short) (LNRC(0,i - 1) + 1);
                i = IND(1,k_);
                LNRC(1,i - 1) = (short) (LNRC(1,i - 1) + 1);
                goto L_1050;
                /*  OVERWRITE THIS SMALL ONE */
L_1040:
                a[k_] = a[lnu - 1];
                IND(0,k_) = IND(0,lnu - 1);
                IND(1,k_) = IND(1,lnu - 1);
                lnu = lnu - 1;
L_1050:
                ;
        }

        dluprm[2] = g;
        nzinl = 0;
        nzinu = lnu;
        lrow = lnu;
        lcol = lnu;

        /*  CONSTRUCT THE ROW LISTS OF COLS IN EACH ROW */
        lsortr( n, nzinu, a, &IND(1,0), ip, &IND(0,0), &LNRC(0,0) );

        /*  CONSTRUCT THE COL LISTS OF ROW INDICES IN EACH COL
         *   IF NZERO = LNU THEN COL LISTS OF ROW INDICES ARE STILL INTACT */
        if( lnu == nzero )
                goto L_1170;
        /*  WE HAVE TO DO IT THE HARD WAY
         *   INITIALIZE IP( I, 2) TO POINT 1 BEYOND THE END OF COL I */
        k = 1;
        for( ir = 1; ir <= n; ir++ ){
                ir_ = ir - 1;
                k = k + LNRC(1,ir_);
                IP(1,ir_) = k;
        }
        /*  FILL EACH COL BACKWARDS FROM TOP SO IP( J, 2) IS LEFT POINTING
         *   TO START OF COL WHEN DONE */
        krl = nzinu;
        /*  LOOP OVER EACH ROW */
        for( ir = n; ir >= 1; ir-- ){
                ir_ = ir - 1;
                kri = IP(0,ir_);
                if( kri > krl )
                        goto L_1155;
                for( k = kri; k <= krl; k++ ){
                        k_ = k - 1;
                        j = IND(1,k_);
                        kr = IP(1,j - 1) - 1;
                        IP(1,j - 1) = kr;
                        IND(0,kr - 1) = (short) ir;
                }
L_1155:
                krl = kri - 1;
        }
        /*  DONE CONSTRUCTING COL LISTS THE HARD WAY */
        goto L_1184;

        /*  CONSTRUCT COL LISTS THE EASY WAY, JUST HAVE TO SET UP IP( J, 2) */
L_1170:
        k = 1;
        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                IP(1,i_) = k;
                k = k + LNRC(1,i_);
        }

        /*  CODE TO IDENTIFY THE TRIANGULAR COLS OF U ABOVE THE BUMP;
         *    IN THIS SECTION, USAGE OF VARIABLES IS:
         *   LNRC( I, 1) = NO. NONZEROES IN ROW I
         *   LNRC( J, 2) = CURRENT NO. NONZEROES IN COL J
         *   LNRC( J, 3) = THE STACK OF ID'S OF SINGLETON COLS
         *   LNRC( J, 4) = ORIGINAL NO. OF NONZEROES IN COL J
         *   LNRC( I, 5) = NEGATIVE OF PIVOT POSITION OF ROW I
         *   LNRC( J, 6) = NEGATIVE OF PIVOT POSITION OF COL J
         *   LNRC( J, 8) = SUCCESSOR OF J IF COL J HAD TO BE DROPPED
         *
         *  PUT ALL SINGLETON COLS ON THE STACK */
L_1184:
        nonstk = 0;
        /*  DO IT BACKWARDS SO COL 1 DEFINITELY GETS SELECTED */
        for( j = n; j >= 1; j-- ){
                j_ = j - 1;
                LNRC(3,j_) = LNRC(1,j_);
                if( LNRC(1,j_) != 1 )
                        goto L_3730;
                nonstk = nonstk + 1;
                LNRC(2,nonstk - 1) = (short) j;
L_3730:
                ;
        }
        /*  RECORD SINGLETON COLS */
        iluprm[13] = nonstk;

        ipv = 0;
        /*  START LOOP: IDENTIFY TOP PART OF U BEFORE BUMP
         *   IS THERE SOMETHING ON THE STACK? */
L_3745:
        if( nonstk <= 0 )
                goto L_3800;
        /*  YES, TAKE TOP ELEMENT OFF THE STACK; JPVT = PIVOT COL */
        jpvt = LNRC(2,nonstk - 1);
        LNRC(2,nonstk - 1) = 0;
        nonstk = nonstk - 1;

        /*  IS THIS COL EMPTY? */
        if( LNRC(1,jpvt - 1) <= 0 )
                goto L_3745;
        /*   NO, SEARCH COL LIST JPVT TO FIND THE PIVOT ROW... */
        kci = IP(1,jpvt - 1);
        kcl = kci + LNRC(3,jpvt - 1) - 1;
        for( k = kci; k <= kcl; k++ ){
                k_ = k - 1;
                /*  IS THIS THE PIVOT ROW */
                ip0 = IND(0,k_);
                /*  WAS THIS ROW ALREADY PIVOTED? */
                if( LNRC(4,ip0 - 1) >= 0 )
                        ipvt = ip0;
                /*   .. AND EMPTY THE LIST */
                IND(0,k_) = 0;
        }
        LNRC(1,jpvt - 1) = 0;

        /*  IPVT IS THE PIVOT ROW, JPVT = PIVOT COL; RECORD SEQUENCE NO. OF PIVOT */
        ipv = ipv + 1;
        LNRC(4,ipvt - 1) = (short) -ipv;
        LNRC(5,jpvt - 1) = (short) -ipv;
        /*  GET READY TO */
        kp = IP(0,ipvt - 1);
        kl = kp + LNRC(0,ipvt - 1) - 1;
        /*   LOOP OVER THE PIVOT ROW */
        for( k = kp; k <= kl; k++ ){
                k_ = k - 1;
                j = IND(1,k_);
                if( j == jpvt )
                        kr = k;
                /*  DECREASE THE COUNT FOR ALL COLS IN THIS ROW */
                LNRC(1,j - 1) = (short) (LNRC(1,j - 1) - 1);
                /*  PUT THIS COL ON THE STACK IF IT HIT ONE */
                if( LNRC(1,j - 1) != 1 )
                        goto L_3790;
                nonstk = nonstk + 1;
                LNRC(2,nonstk - 1) = (short) j;
L_3790:
                ;
        }
        /*  BRING PIVOT TO FRONT OF ITS ROW */
        au = a[kr - 1];
        a[kr - 1] = a[kp - 1];
        a[kp - 1] = au;
        IND(1,kr - 1) = IND(1,kp - 1);
        IND(1,kp - 1) = (short) jpvt;
        goto L_3745;
        /*  END OF LOOP TO IDENTIFY TOP PART OF U
         *
         *  NO. SPIKES THIS FAR */
L_3800:
        nspike = 0;
        /*  NO. TRIANGULAR COLS */
        iluprm[14] = ipv;
        /*  NO. VECTORS IN L */
        nvinl = 0;
        /*  POSN OF LAST COL IN L */
        lstinl = n;

        /*  DELETE PIVOTED ROWS FROM COL LISTS OF UNPIVOTED COLS
         *    WHERE ANY COLS UNPIVOTED? */
        if( ipv == n )
                goto L_5000;
        /*   YES, TRACK THEM DOWN */
        for( j = 1; j <= n; j++ ){
                j_ = j - 1;
                lnrc4 = LNRC(3,j_);
                LNRC(3,j_) = 0;
                if( LNRC(5,j_) < 0 )
                        goto L_3890;
                kci = IP(1,j_);
                kcl = kci + lnrc4 - 1;
                if( kci > kcl )
                        goto L_3890;
                k1 = kci;
                for( k = kci; k <= kcl; k++ ){
                        k_ = k - 1;
                        i = IND(0,k_);
                        IND(0,k_) = 0;
                        /*  HAS ROW I ALREADY BEEN PIVOTED? */
                        if( LNRC(4,i - 1) < 0 )
                                goto L_3880;
                        /*   NO, COPY THIS ELEMENT DOWN TO THE GOOD PART */
                        IND(0,k1 - 1) = (short) i;
                        k1 = k1 + 1;
L_3880:
                        ;
                }
L_3890:
                ;
        }

        /*  PUT UNPIVOTED ROWS AND COLS IN LINKED LISTS..
         *   ACCORDING TO NO. NON-ZEROS */
        for( l = 1; l <= 2; l++ ){
                l_ = l - 1;
                /*  DO IN REVERSE ORDER SO COL 1 GETS SELECTED FIRST */
                for( i = n; i >= 1; i-- ){
                        i_ = i - 1;
                        /*  WAS THIS VECTOR ALREADY PIVOTED? */
                        if( LNRC(l_ + 4,i_) < 0 )
                                goto L_3970;
                        lntmp = LNRC(l_,i_);
                        /*  IS THIS VECTOR EMPTY? */
                        if( lntmp > 0 )
                                goto L_3960;
                        /*  YES, PUT IR ON THE LIST OF DROPPED ROWS/COLS */
                        LNRC(l_ + 6,i_) = (short) mt[l_];
                        mt[l_] = i;
                        ndrpd[l_] = ndrpd[l_] + 1;
                        goto L_3970;
                        /*   NO, GET CURRENT 1ST ON LIST */
L_3960:
                        in = LNRC(l_ + 2,lntmp - 1);
                        LNRC(l_ + 2,lntmp - 1) = (short) i;
                        LNRC(l_ + 6,i_) = (short) in;
                        LNRC(l_ + 4,i_) = 0;
                        if( in != 0 )
                                LNRC(l_ + 4,in - 1) = (short) i;
L_3970:
                        ;
                }
        }


        /*  DURING PROCESSING OF BUMP..
         *    LNRC( I, 3/4) = FIRST ROW/COL WITH I NONZEROES, 0 IF NONE
         *    LNRC( I, 5/6) = ROW/COL PRIOR TO I ON ITS LIST, 0 IF NONE,
         *                   IF I HAS NOT YET BEEN PIVOTED.
         *                - NEGATIVE OF PIVOT POSITION IF I HAS BEEN PIVOTED
         *    LNRC( I, 7/8) = ROW/COL AFTER I ON ITS LIST, 0 IF NONE
         *
         *  NOW GET READY TO PIVOT NON-TRIANGULAR COLS
         *   RESUME WHERE WE LEFT OFF */
        ipv1 = ipv + 1;
        /*  LOWSTL = ( LOWEST POSITION OCCUPIED BY L ) - 1 */
        lowstl = nzdim - nzinl;

        for( ipv = ipv1; ipv <= n; ipv++ ){
                ipv_ = ipv - 1;
                /*  GET OUT IF ALL COLS HAVE BEEN EITHER PIVOTED OR DROPPED. */
                if( ipv + (ndrpd[0] + ndrpd[1])/2 > n )
                        goto L_5000;

                /*  FIND PIVOT IN ROW IPVT AND COL JPVT
                 *   MKOST IS MARKOWITZ COST OF CHEAPEST PIVOT FOUND SO FAR,
                 *     WHICH IS IN ROW IPVT AND COLUMN JPVT. */
                mkost = n*n;
                kcost = 0;
                /*  NPS = NO. ROW/COL LENGTHS CHECKED THUS FAR */
                nps = 0;
                /*  LOOP OVER LENGTH OF COL/ROW BEING SEARCHED */
                for( nzcr = 1; nzcr <= n; nzcr++ ){
                        nzcr_ = nzcr - 1;
                        kp2 = nzcr - 1;
                        kp2 = kp2*kp2;
                        if( mkost <= kp2 )
                                goto L_4400;
                        j = LNRC(3,nzcr_);
                        /*  SEARCH COLUMNS WITH NZCR NON-ZEROS. */
                        for( idum = 1; idum <= n; idum++ ){
                                idum_ = idum - 1;
                                if( j <= 0 )
                                        goto L_4300;
                                kci = IP(1,j - 1);
                                kcl = kci + LNRC(1,j - 1) - 1;

                                /*   PREPROCESS THE LIST OF ROW IDS TO PUT SHORT ROWS FIRST
                                 *    DO NOT BOTHER IF ONLY 1 ELEMENT */
                                if( nzcr <= 1 )
                                        goto L_4245;
                                /*  NZS = SHORTEST ROW FOUND THUS FAR IN THIS COL LIST */
                                nzs = n;
                                for( k = kci; k <= kcl; k++ ){
                                        k_ = k - 1;
                                        i = IND(0,k_);
                                        nzr = LNRC(0,i - 1);
                                        if( nzr > nzs )
                                                goto L_4240;
                                        if( nzr == nzs )
                                                goto L_4230;
                                        /*  HAVE WE ALREADY CHECKED THIS ROW? */
                                        if( nzr < nzcr )
                                                goto L_4240;
                                        /*  NO */
                                        kp1 = kci;
                                        nzs = nzr;
L_4230:
                                        IND(0,k_) = IND(0,kp1 - 1);
                                        IND(0,kp1 - 1) = (short) i;
                                        /*  GET OUT IF THIS LENGTH IS GOOD ENOUGH */
                                        if( nzs == nzcr )
                                                goto L_4245;
                                        kp1 = kp1 + 1;
L_4240:
                                        ;
                                }

                                /*  LOOP OVER ROWS IN THIS COL TO FIND BEST */
L_4245:
                                for( k = kci; k <= kcl; k++ ){
                                        k_ = k - 1;
                                        i = IND(0,k_);
                                        /*  IF IT HAS ONLY ONE ROW... */
                                        if( nzcr > 1 )
                                                goto L_4247;
                                        /*   WE TAKE THIS COL JPVT AND ROW I */
                                        jpvt = j;
                                        goto L_4270;
L_4247:
                                        nzr = LNRC(0,i - 1);
                                        kcost = (nzcr - 1)*(nzr - 1);
                                        if( kcost >= mkost )
                                                goto L_4280;
                                        /*  HAVE WE ALREADY CHECKED THIS ROW? */
                                        if( nzr < nzcr )
                                                goto L_4280;
                                        /*  NO, FIND COL IN ROW I WITH NZCR NONZEROES WITH LARGEST COEF DPIV */
                                        dpiv = 0.0;
                                        /*  FIND LARGEST ELEMENT IN ROW, AMAX */
                                        amax = dpiv;
                                        k1 = IP(0,i - 1);
                                        k2 = LNRC(0,i - 1) + k1 - 1;
                                        for( kk = k1; kk <= k2; kk++ ){
                                                kk_ = kk - 1;
                                                auab = fabs( a[kk_] );
                                                /*  COULD THIS ONE BE BETTER THAN BIRD IN HAND? */
                                                if( auab <= dpiv )
                                                        goto L_4260;
                                                /*   YES, IS IT THE OVERALL MAX IN ROW? */
                                                if( auab > amax )
                                                        amax = auab;
                                                /*  IS THIS COL SHORT ENOUGH? */
                                                if( LNRC(1,IND(1,kk_) - 1) != nzcr )
                                                        goto L_4260;
                                                /*   YES, MAKE IT THE NEW CANDIDATE */
                                                dpiv = auab;
                                                kj = kk;
                                                /*  SHOULD WE CHECK SIZE EARLY HERE? */
L_4260:
                                                ;
                                        }
                                        /*  IS CANDIDATE PIVOT BIG ENOUGH IN ITS ROW? */
                                        if( dpiv < amax*u )
                                                goto L_4280;
                                        /*   YES */
                                        jpvt = IND(1,kj - 1);
L_4270:
                                        mkost = kcost;
                                        ipvt = i;
                                        if( mkost <= kp2 )
                                                goto L_4400;
L_4280:
                                        ;
                                }
                                /*  GET NEXT COL WITH NZCR NONZEROES */
                                j = LNRC(7,j - 1);
                        }
                        nps = nps + 1;
                        /*  SHOULD WE BE CONTENT WITH WHAT WE HAVE? */
                        if( nps < maxps )
                                goto L_4300;
                        if( mkost < n*n )
                                goto L_4400;

                        /*  SEARCH ROWS WITH NZCR NON-ZEROS. */
L_4300:
                        i = LNRC(2,nzcr_);
                        for( idum = 1; idum <= n; idum++ ){
                                idum_ = idum - 1;
                                if( i <= 0 )
                                        goto L_4380;
                                amax = 0.;
                                kp = IP(0,i - 1);
                                kl = kp + LNRC(0,i - 1) - 1;
                                /*  FIND LARGEST COEF IN THE ROW */
                                for( k = kp; k <= kl; k++ ){
                                        k_ = k - 1;
                                        amax = fmax( fabs( a[k_] ), amax );
                                }
                                au = amax*u;
                                for( k = kp; k <= kl; k++ ){
                                        k_ = k - 1;
                                        /*  IS PROSPECTIVE PIVOT BIG ENOUGH IN ITS ROW? */
                                        if( fabs( a[k_] ) < au )
                                                goto L_4360;
                                        j = IND(1,k_);
                                        kcost = (nzcr - 1)*(LNRC(1,j - 1) - 1);
                                        if( kcost >= mkost )
                                                goto L_4360;
                                        ipvt = i;
                                        jpvt = j;
                                        mkost = kcost;
                                        /*  IS MKOST .LE. ( NZCR - 1) * ( NZCR - 1)? */
                                        if( mkost <= kp2 )
                                                goto L_4400;
L_4360:
                                        ;
                                }
                                /*  GET NEXT ROW WITH NZCR NONZEROES */
                                i = LNRC(6,i - 1);
                        }
                        nps = nps + 1;
L_4380:
                        ;
                }

                /*  NOW DO A PIVOT IN ROW IPVT, COL JPVT
                 *  REMOVE ROWS/COLS WITH COEFS IN PIVOT ROW/COL FROM LINKED LISTS */
L_4400:
                kpri = IP(1,jpvt - 1);
                kprl = LNRC(1,jpvt - 1) + kpri - 1;
                for( l = 1; l <= 2; l++ ){
                        l_ = l - 1;
                        for( k = kpri; k <= kprl; k++ ){
                                k_ = k - 1;
                                i = IND(l_,k_);
                                il = LNRC(l_ + 4,i - 1);
                                in = LNRC(l_ + 6,i - 1);
                                if( il == 0 )
                                        goto L_4420;
                                LNRC(l_ + 6,il - 1) = (short) in;
                                goto L_4430;
L_4420:
                                lntmp = LNRC(l_,i - 1);
                                LNRC(l_ + 2,lntmp - 1) = (short) in;
L_4430:
                                if( in > 0 )
                                        LNRC(l_ + 4,in - 1) = (short) il;
                        }
                        kpri = IP(0,ipvt - 1);
                        kprl = kpri + LNRC(0,ipvt - 1) - 1;
                }

                /*  RECORD SEQUENCE NO. OF PIVOT OF ROW IPVT AND COL JPVT */
                LNRC(4,ipvt - 1) = (short) -ipv;
                LNRC(5,jpvt - 1) = (short) -ipv;
                /*  REMOVE PIVOT ROW IPVT FROM COLUMN LISTS AND FIND JPVT IN ROW IPVT. */
                for( k = kpri; k <= kprl; k++ ){
                        k_ = k - 1;
                        j = IND(1,k_);
                        if( j == jpvt )
                                kr = k;
                        /*  FIND PIVOT ROW IN COL J AND DELETE IT(DO WE NEED TO DO THIS NOW?) */
                        kpc = IP(1,j - 1);
                        LNRC(1,j - 1) = (short) (LNRC(1,j - 1) - 1);
                        klc = kpc + LNRC(1,j - 1);
                        for( kc = kpc; kc <= klc; kc++ ){
                                kc_ = kc - 1;
                                if( ipvt == IND(0,kc_) )
                                        goto L_4470;
                        }
L_4470:
                        IND(0,kc - 1) = IND(0,klc - 1);
                        IND(0,klc - 1) = 0;
                }
                /*  MOVE PIVOT TO FRONT OF ROW IPVT */
                dpiv = a[kr - 1];
                a[kr - 1] = a[kpri - 1];
                a[kpri - 1] = dpiv;
                IND(1,kr - 1) = IND(1,kpri - 1);
                IND(1,kpri - 1) = (short) jpvt;
                /*  KPR2 POINTS TO START OF NONPIVOT COEFS IN PIVOT ROW IPVT */
                kpr2 = kpri + 1;

                /*  PREPARE TO SUBTRACT PIVOT ROW FROM OTHER ROWS IN PIVOT COL */
                nzc = LNRC(1,jpvt - 1);
                if( nzc == 0 )
                        goto L_4850;
                /*  BUMP NO. COLS IN L */
                nvinl = nvinl + 1;
                nfill = kprl - kpri;
                /*  RECORD POSN OF LAST COL IN L */
                lstinl = ipv;
                if( nfill <= 0 )
                        goto L_4500;
                /*  INCREMENT NO. SPIKES */
                nspike = nspike + 1;
                /*  LOAD PIVOT ROW, SANS PIVOT ELEMENT, INTO W() */
                for( k = kpr2; k <= kprl; k++ ){
                        k_ = k - 1;
                        w[IND(1,k_) - 1] = a[k_];
                }

                /*  LOOP ON NON-ZEROS IN PIVOT COLUMN TO ADD TO L AND ELIMINATE FROM U */
L_4500:
                for( nc = 1; nc <= nzc; nc++ ){
                        nc_ = nc - 1;
                        kc = IP(1,jpvt - 1) + nc - 1;
                        ir = IND(0,kc - 1);
                        /*  IN NONPIVOT ROW IR, FIND ELEMENT IN PIVOT COL JPVT */
                        kr = IP(0,ir - 1);
                        krl = kr + LNRC(0,ir - 1) - 1;
                        for( knp = kr; knp <= krl; knp++ ){
                                knp_ = knp - 1;
                                if( jpvt == IND(1,knp_) )
                                        goto L_4530;
                        }
                        /*  COMPUTE THE NEXT TERM IN L */
L_4530:
                        am = -a[knp - 1]/dpiv;
                        /*  OVERWRITE PIV COL COEF WITH LAST COEF IN THE ROW IR */
                        a[knp - 1] = a[krl - 1];
                        IND(1,knp - 1) = IND(1,krl - 1);
                        IND(1,krl - 1) = 0;
                        krl = krl - 1;
                        LNRC(0,ir - 1) = (short) (LNRC(0,ir - 1) - 1);

                        /*  PREPARE TO UPDATE EXISTING ELEMENTS IN ROW IR
                         *  GROMX = LARGEST ELEMENT FOUND IN ROW */
                        gromx = 0.0;
                        /*  K1 = LOCATION OF BIGGEST ELEMENT IN NON-PIVOT ROW IR */
                        k1 = 0;
                        /*  NFILLR = NO. FILL ELEMENTS REMAINING */
                        nfillr = nfill;
                        if( kr > krl )
                                goto L_4650;
                        /*  LOOP OVER EXISTING ELEMENTS IN ROW IR TO UPDATE THEM IN PLACE */
                        kj = krl;
                        for( ks = kj; ks >= kr; ks-- ){
                                ks_ = ks - 1;
                                j = IND(1,ks_);
                                au = w[j - 1];
                                if( au != 0.0 )
                                        goto L_4540;
                                auab = fabs( a[ks_] );
                                goto L_4550;
L_4540:
                                nfillr = nfillr - 1;
                                /*  THIS COL CAN NO LONGER PRODUCE A FILL */
                                w[j - 1] = 0.0;
                                au = a[ks_] + am*au;
                                auab = fabs( au );
                                if( auab <= zerotl )
                                        goto L_4560;
                                a[ks_] = au;
                                /*  SHOULD THIS ELEMENT GO FIRST? */
L_4550:
                                if( auab <= gromx )
                                        goto L_4590;
                                /*  YES */
                                gromx = auab;
                                k1 = ks;
                                goto L_4590;

                                /*  REMOVE AN ELEMENT FROM U */
L_4560:
                                nzinu = nzinu - 1;
                                /*  OVERWRITE IT IN ROW FILE */
                                a[ks_] = a[krl - 1];
                                IND(1,ks_) = IND(1,krl - 1);
                                IND(1,krl - 1) = 0;
                                if( k1 == krl )
                                        k1 = ks;
                                krl = krl - 1;
                                /*  GET READY TO DELETE FROM COL FILE */
                                k = IP(1,j - 1);
                                kl = k + LNRC(1,j - 1) - 1;
                                LNRC(1,j - 1) = (short) (kl - k);
                                /*  FIND ELEMENT IN COL FILE */
                                for( kk = k; kk <= kl; kk++ ){
                                        kk_ = kk - 1;
                                        if( IND(0,kk_) == ir )
                                                goto L_4580;
                                }
                                /*  AND REMOVE IT BY OVERWRITING IT */
L_4580:
                                IND(0,kk - 1) = IND(0,kl - 1);
                                IND(0,kl - 1) = 0;
L_4590:
                                ;
                        }

                        /*  PREPARE TO DO ANY FILL ELEMENTS, ARE THERE ANY? */
L_4650:
                        if( kpr2 > kprl )
                                goto L_4800;
                        /*  CAN WE LEAVE ROW IR IN PLACE BECAUSE FILL <= 1? */
                        if( nfillr <= 1 )
                                goto L_4700;
                        /*  NO, MAKE K1 RELATIVE TO ROW START KR, TEMPORARILY */
                        k1 = k1 - kr;
                        /*   DO WE NEED TO COMPRESS TO MAKE ROOM FOR NEW ROW? */
                        if( lrow + LNRC(0,ir - 1) + nfillr <= lowstl )
                                goto L_4670;
                        /*  YES, COMPRESS IF IT IS WORTHWHILE */
                        if( ntcmp >= mxcmp || nzinu + LNRC(0,ir - 1) + nfillr >
                         lowstl )
                                goto L_8655;
                        /*  MUST UPDATE LNRC( IR) FOR LSLUCR */
                        LNRC(0,ir - 1) = (short) (krl - kr + 1);
                        lslucr( n, nzdim, &IND(1,0), a, ip, lnrc, &lrow, &ntcmp );
                        kpri = IP(0,ipvt - 1);
                        kpr2 = kpri + 1;
                        kprl = kpri + LNRC(0,ipvt - 1) - 1;
                        kr = IP(0,ir - 1);
                        krl = kr + LNRC(0,ir - 1) - 1;

                        /*  TRANSFER ROW IR TO HIGH MEMORY FOR ROOM TO GROW */
L_4670:
                        IP(0,ir - 1) = lrow + 1;
                        if( kr > krl )
                                goto L_4690;
                        for( k = kr; k <= krl; k++ ){
                                k_ = k - 1;
                                lrow = lrow + 1;
                                a[lrow - 1] = a[k_];
                                IND(1,lrow - 1) = IND(1,k_);
                                IND(1,k_) = 0;
                        }
L_4690:
                        krl = lrow;
                        kr = IP(0,ir - 1);
                        /*  RESTORE K1 TO ABSOLUTE FORM */
                        if( k1 > 0 )
                                k1 = kr + k1;

                        /*  LOOP OVER PIVOT ROW IPVT TO CREATE THE FILLS */
L_4700:
                        for( ks = kpr2; ks <= kprl; ks++ ){
                                ks_ = ks - 1;
                                j = IND(1,ks_);
                                au = w[j - 1];
                                /*  THIS COL ALREADY PROCESSED? */
                                if( au != 0.0 )
                                        goto L_4710;
                                /*   YES, PREPARE IT FOR THE NEXT ROW */
                                w[j - 1] = a[ks_];
                                goto L_4790;
L_4710:
                                au = am*au;
                                auab = fabs( au );
                                if( auab <= zerotl )
                                        goto L_4790;
                                krl = krl + 1;
                                /*  SHOULD THIS ELEMENT GO FIRST? */
                                if( auab <= gromx )
                                        goto L_4720;
                                /*  YES */
                                gromx = auab;
                                k1 = krl;
L_4720:
                                a[krl - 1] = au;
                                IND(1,krl - 1) = (short) j;
                                nzinu = nzinu + 1;

                                /*  PUT NEW COEF IN COL LIST */
                                lntmp = LNRC(1,j - 1);
                                kci = IP(1,j - 1);
                                kcl = kci + lntmp - 1;
                                /*  TRY TO PLACE NEW COEF AT END OF EXISTING COL LIST */
                                if( kcl != lcol )
                                        goto L_4730;
                                if( lcol >= lowstl )
                                        goto L_4740;
                                goto L_4770;
L_4730:
                                if( IND(0,kcl) != 0 )
                                        goto L_4740;
                                IND(0,kcl) = (short) ir;
                                goto L_4780;

                                /*  COL LIST MUST BE MOVED TO END; IS THERE ROOM? */
L_4740:
                                if( lcol + lntmp + 1 < lowstl )
                                        goto L_4750;
                                /*  NO, WILL COMPRESSION HELP? */
                                if( ntcmp >= mxcmp || nzinu + lntmp + 1 >= lowstl )
                                        goto L_8655;
                                /*  YES, DO COMPRESSION OF COL LISTS */
                                lslucc( n, nzdim, ind, &IP(1,0), &LNRC(1,0), &lcol,
                                 &ntcmp );
                                kci = IP(1,j - 1);
                                kcl = kci + lntmp - 1;
                                /*  NOW MOVE THIS COL LIST TO END */
L_4750:
                                IP(1,j - 1) = lcol + 1;
                                for( kk = kci; kk <= kcl; kk++ ){
                                        kk_ = kk - 1;
                                        lcol = lcol + 1;
                                        IND(0,lcol - 1) = IND(0,kk_);
                                        IND(0,kk_) = 0;
                                }

                                /*  ADD NEW ENTRY AT END OF COL LIST */
L_4770:
                                lcol = lcol + 1;
                                IND(0,lcol - 1) = (short) ir;
L_4780:
                                LNRC(1,j - 1) = (short) (lntmp + 1);
L_4790:
                                ;
                        }
                        /*  IF UPDATED ROW WAS PUT AT END, UPDATE LROW */
                        if( krl > lrow )
                                lrow = krl;

                        /*  MOVE BIGGEST COEF TO FRONT OF ROW IR, IF THERE WAS ONE */
L_4800:
                        g = fmax( g, gromx );
                        if( k1 == 0 )
                                goto L_4815;
                        auab = a[kr - 1];
                        a[kr - 1] = a[k1 - 1];
                        a[k1 - 1] = auab;
                        kk = IND(1,kr - 1);
                        IND(1,kr - 1) = IND(1,k1 - 1);
                        IND(1,k1 - 1) = (short) kk;

L_4815:
                        LNRC(0,ir - 1) = (short) (krl - kr + 1);

                        /*  STORE MULTIPLIER DOWNWARD IN L; IS THERE ROOM?
                         *   WILL IT RUN INTO COL STORAGE OF U */
                        if( lcol + 1 <= lowstl )
                                goto L_4820;
                        /*  COMPRESS COL FILE IF IT WILL HELP */
                        if( ntcmp >= mxcmp || nzinu + 1 > lowstl )
                                goto L_8655;
                        lslucc( n, nzdim, ind, &IP(1,0), &LNRC(1,0), &lcol, &ntcmp );
                        /*   WILL IT RUN INTO ROW STORAGE OF U */
L_4820:
                        if( lrow + 1 < lowstl )
                                goto L_4825;
                        /*  COMPRESS ROW FILE IF IT WILL HELP */
                        if( ntcmp >= mxcmp || nzinu + 1 > lowstl )
                                goto L_8655;
                        lslucr( n, nzdim, &IND(1,0), a, ip, lnrc, &lrow, &ntcmp );
                        kpri = IP(0,ipvt - 1);
                        kpr2 = kpri + 1;
                        kprl = kpri + LNRC(0,ipvt - 1) - 1;
L_4825:
                        k = lowstl;
                        nzinl = nzinl + 1;
                        lowstl = lowstl - 1;
                        a[k - 1] = am;
                        IND(0,k - 1) = (short) ipvt;
                        IND(1,k - 1) = (short) ir;
                        /*  COEF AT ( IR, JPVT) IS NO LONGER IN U */
                        nzinu = nzinu - 1;
                }

                /*  ZERO OUT THE W() VECTOR */
                if( nfill <= 0 )
                        goto L_4850;
                for( k = kpr2; k <= kprl; k++ ){
                        k_ = k - 1;
                        w[IND(1,k_) - 1] = 0.0;
                }

                /*  PUT ROWS/COLS WITH COEFS IN PIVOT COL/ROW BACK ON LINKED LISTS
                 *    WITH SAME NO. OF NON-ZEROS
                 *   THEY WERE PREVIOUSLY REMOVED FROM THEIR OLD HOMES
                 *  GET THE START OF THE PIVOT COLUMN */
L_4850:
                k1 = IP(1,jpvt - 1);
                /*   AND ITS END */
                k2 = k1 + LNRC(1,jpvt - 1) - 1;
                /*  PIVOT COL IS NO LONGER IN THE RUNNING */
                LNRC(1,jpvt - 1) = 0;
                for( l = 1; l <= 2; l++ ){
                        l_ = l - 1;
                        if( k2 < k1 )
                                goto L_4880;
                        for( k = k1; k <= k2; k++ ){
                                k_ = k - 1;
                                ir = IND(l_,k_);
                                if( l == 1 )
                                        IND(l_,k_) = 0;
                                lntmp = LNRC(l_,ir - 1);
                                /*  DID THIS VECTOR BECOME EMPTY? */
                                if( lntmp > 0 )
                                        goto L_4860;
                                /*  YES, PUT IR ON THE LIST OF DROPPED ROWS/COLS */
                                LNRC(l_ + 6,ir - 1) = (short) mt[l_];
                                mt[l_] = ir;
                                ndrpd[l_] = ndrpd[l_] + 1;
                                goto L_4870;

                                /*  NO, GET ID OF FIRST ENTRY IN ITS NEW HOME */
L_4860:
                                in = LNRC(l_ + 2,lntmp - 1);
                                /*  PUT IR BEFORE OLD FIRST ENTRY */
                                LNRC(l_ + 6,ir - 1) = (short) in;
                                /*  IR HAS NO PREDECESSOR BECAUSE IT IS FIRST */
                                LNRC(l_ + 4,ir - 1) = 0;
                                /*   POINT TO THE FIRST ENTRY ON THIS LIST, I.E., IR */
                                LNRC(l_ + 2,lntmp - 1) = (short) ir;
                                /*  IF THE LIST HAS A SECOND ELEMENT, POINT IT BACK TO IR */
                                if( in != 0 )
                                        LNRC(l_ + 4,in - 1) = (short) ir;
L_4870:
                                ;
                        }

                        /*  NOW GET ROW START */
L_4880:
                        k1 = IP(0,ipvt - 1) + 1;
                        /*  AND ITS END */
                        k2 = k1 + LNRC(0,ipvt - 1) - 2;
                }
        }
        /*  END OF MAIN LOOP
         *
         *  WERE ANY ROWS DROPPED? */
L_5000:
        if( ndrpd[0] <= 0 )
                goto L_5500;
        /*   YES, PUT THEM AT END OF PIVOT ORDER, IPV = NEXT PIVOT POSN
         *  DO WE HAVE ENOUGH ROOM IN THE ROW FILE? */
        if( lrow + ndrpd[0] < lowstl )
                goto L_5050;
        /*   NOT YET, TRY TO RECLAIM SOME */
        lslucr( n, nzdim, &IND(1,0), a, ip, lnrc, &lrow, &ntcmp );
        /*  NOW IS THERE ENOUGH ROOM? */
        if( lrow + ndrpd[0] >= lowstl )
                goto L_8655;
        /*  YES
         *  IS THERE ENOUGH ROOM IN THE COL FILE? */
L_5050:
        if( lcol + ndrpd[0] < lowstl )
                goto L_5060;
        /*   NOT YET, TRY TO RECLAIM SOME */
        lslucc( n, nzdim, ind, &IP(1,0), &LNRC(1,0), &lcol, &ntcmp );
        /*  NOW IS THERE ENOUGH? */
        if( lcol + ndrpd[0] >= lowstl )
                goto L_8655;
L_5060:
        ir = mt[0];
        jr = mt[1];
        /*  LOOP OVER DROPPED ROWS AND COLS */
        for( ipu = ipv; ipu <= n; ipu++ ){
                ipu_ = ipu - 1;
                /*  ADD A 1.0 TO U
                 *   IN THE ROW FILE... */
                lrow = lrow + 1;
                IP(0,ir - 1) = lrow;
                LNRC(0,ir - 1) = 1;
                IND(1,lrow - 1) = (short) jr;
                a[lrow - 1] = 1.0;
                /*   AND IN THE COL FILE */
                lcol = lcol + 1;
                IP(1,jr - 1) = lcol;
                LNRC(1,jr - 1) = 1;
                IND(0,lcol - 1) = (short) ir;
                nzinu = nzinu + 1;
                /*  GIVE PIVOT POSITIONS */
                LNRC(4,ir - 1) = (short) -ipu;
                LNRC(5,jr - 1) = (short) -ipu;
                /*  GET NEXT ROW AND COL */
                ir = LNRC(6,ir - 1);
                jr = LNRC(7,jr - 1);
        }

        /*  RETURN NO. OF DROPPED COLUMNS */
L_5500:
        *kdropd = ndrpd[0];

        /*  SET:
         *    LNRC( I, 3) = ROW IN ITH POSITION FROM TOP OF U
         *    LNRC( I, 4) = COL IN ITH POSITION FROM TOP OF U
         *    LNRC( I, 5) = POSITION FROM TOP OF U OF ITH ROW */
        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                /*  J = POSITION OF ROW I IN SEQUENCE */
                j = -LNRC(4,i_);
                LNRC(4,i_) = (short) j;
                /*  SO JTH PIVOT IS ROW I */
                LNRC(2,j - 1) = (short) i;
                /*  SAME THING FOR COLUMNS */
                j = -LNRC(5,i_);
                LNRC(3,j - 1) = (short) i;
                LNRC(1,i_) = 0;
        }

        /*  CONSTRUCT COLUMN FILE FOR U
         *   PUT LENGTH OF EACH COL IN LNRC(.,2) */
        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                kp = IP(0,i_);
                kl = LNRC(0,i_) + kp - 1;
                for( k = kp; k <= kl; k++ ){
                        k_ = k - 1;
                        j = IND(1,k_);
                        LNRC(1,j - 1) = (short) (LNRC(1,j - 1) + 1);
                }
        }

        /*  POINT IP( I, 2) TO THE END OF COL I */
        k = 1;
        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                k = k + LNRC(1,i_);
                IP(1,i_) = k;
        }

        /*  NOW PUT THE NONZEROES INTO PLACE AND
         *   AND DECREMENT IP( J, 2) TO THE START OF COL J */
        lcol = k - 1;
        for( ii = 1; ii <= n; ii++ ){
                ii_ = ii - 1;
                i = LNRC(2,ii_);
                kp = IP(0,i - 1);
                kl = LNRC(0,i - 1) + kp - 1;
                for( k = kp; k <= kl; k++ ){
                        k_ = k - 1;
                        j = IND(1,k_);
                        kn = IP(1,j - 1) - 1;
                        IP(1,j - 1) = kn;
                        IND(0,kn - 1) = (short) i;
                }
        }

        /*  RETURN INVERSE SHAPE */
        iluprm[5] = nzinl;
        iluprm[11] = nzinl;
        iluprm[6] = nzinu;
        iluprm[12] = nzinu;
        iluprm[7] = lcol;
        iluprm[8] = lrow;
        iluprm[9] = ntcmp;
        iluprm[10] = mxcmp;
        /*  RETURN STATS */
        iluprm[15] = nspike;
        iluprm[16] = nvinl;
        iluprm[17] = n - lstinl;
        dluprm[0] = u;
        dluprm[3] = g;
        goto L_9100;

        /*  VARIOUS ERROR RETURNS
         *  RAN OUT OF SPACE */
L_8655:
        *kalinv = 9;
        goto L_9100;

        /*  COL HAD MORE THAN ONE ENTRY IN A ROW */
/* L_8700: unreferenced??? */
        *kalinv = 10;
        goto L_9100;

        /*  ASKED TO INVERT A MATRIX WITH <= 0 ROWS */
L_8800:
        *kalinv = 11;

L_9100:
        return;
#undef  LNRC
#undef  IND
#undef  IP
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS lsluft(long int m, long int nzdim, double a[], short int *ind,
         long int *ip, short int *lnrc, double w[], double b[], long int *nzinw,
         long int *nzinb, long int iluprm[])
{
#define IND(I_,J_)      (*(ind+(I_)*(nzdim)+(J_)))
#define IP(I_,J_)       (*(ip+(I_)*(m)+(J_)))
#define LNRC(I_,J_)     (*(lnrc+(I_)*(m)+(J_)))
long int i, i_, ii, ii_, j, k, k_, kri, krl, l2, lntmp, lowstl, nzint;
double dtmp;


        /*  DO AN FTRAN
         *
         *  INPUTS:
         *    W() = VECTOR TO BE PRE-MULTIPLIED BY BASIS INVERSE
         *    NZINW = NO. OF NONZEROES IN W(),
         *      IF INCONVENIENT TO PROVIDE THIS, SET NZINW = 0 AND
         *      IT WILL BE OK, JUST SLOWER
         *    LNRC( ,7) = LIST OF ROWS OF W() CONTAINING A NONZERO.
         *      IF INCONVENIENT TO PROVIDE THIS, SET NZINW = 0
         *    M = NO. OF ROWS/COLS
         *    A() = NONZEROES OF LU INVERSE
         *    NZDIM = LENGTH OF THE A() VECTOR
         *    ILUPRM( 6) = NO. OF ELEMENTS IN L
         *    IP( I, 1), IP( I, 2) POINT TO START OF ROW/COLUMN I OF U.
         *    LNRC( I, 1), LNRC( I, 2) ARE LENGTHS OF ROW/COL I OF U.
         *    LNRC(.,3), LNRC(.,4) HOLD ROW/COL NUMBERS IN PIVOTAL ORDER.
         *  OUTPUTS:
         *    W() = ( L INVERSE) * W()
         *    B() = ( U INVERSE) * ( L INVERSE) * W()
         *    LNRC( , 7) = LIST OF INDICES IN W() WITH NONZEROES
         *    LNRC( , 8) = LIST OF INDICES IN B() WITH NONZEROES
         *    ILUPRM( 19) = NO. NONZEROES IN W()
         *    NZINW = NO. OF NONZEROES IN W()
         *    ILUPRM( 20) = NO. NONZEROES IN B()
         *    NZINB = NO. OF NONZEROES IN B()
         * */

        /*  DATE 31 DEC 1989
         *
         *  FIRST SET W() = L INVERSE * W()
         *   ARE THERE ANY ELEMENTS IN L? */
        if( iluprm[5] <= 0 )
                goto L_1800;
        /*   YES, GET THEIR LOW POINT */
        lowstl = nzdim - iluprm[5] + 1;
        for( k = nzdim; k >= lowstl; k-- ){
                k_ = k - 1;
                /*  GET THE PIVOT ROW */
                i = IND(0,k_);
                if( w[i - 1] == 0. )
                        goto L_1700;
                j = IND(1,k_);
                w[j - 1] = w[j - 1] + w[i - 1]*a[k_];
L_1700:
                ;
        }
        goto L_1900;

L_1800:
        if( *nzinw > 0 )
                goto L_2000;

        /*  WE CONSTRUCT THE LIST OF NONZEROES IN W() */
L_1900:
        nzint = 0;
        for( i = 1; i <= m; i++ ){
                i_ = i - 1;
                if( w[i_] == 0.0 )
                        goto L_1990;
                nzint = nzint + 1;
                LNRC(6,nzint - 1) = (short) i;
L_1990:
                ;
        }
        *nzinw = nzint;
        /*#############################
         *   CHANGED 7/14/90 -- SHS
         *      MUST SET ILUPRM(19)
         *#############################
         * */
L_2000:
        iluprm[18] = *nzinw;
        /*  NOW SET B() =  U INVERSE * W()
         *2000 NZINT = 0 */
        nzint = 0;
        for( ii = m; ii >= 1; ii-- ){
                ii_ = ii - 1;
                /*  GET THE ROW */
                i = LNRC(2,ii_);
                dtmp = w[i - 1];
                kri = IP(0,i - 1);
                /*  DID LOWER PIVOTS IN U AFFECT THIS ROW */
                if( kri > 0 )
                        goto L_2700;
                /*  YES.. */
                kri = -kri;
                IP(0,i - 1) = kri;
                lntmp = LNRC(0,i - 1);
                krl = kri - 1 + lntmp;
                /*  LOOP OVER LOWER PIVOTS IN U */
                for( k = kri + 1; k <= krl; k++ ){
                        k_ = k - 1;
                        j = IND(1,k_);
                        dtmp = dtmp - b[j - 1]*a[k_];
                }
                /*  DONE ACCOUNTING FOR LOWER PIVOTS IN U
                 *  GET THE COLUMN THAT IS BASIC IN THIS ROW */
L_2700:
                j = IND(1,kri - 1);
                if( dtmp != 0.0 )
                        goto L_2750;
                b[j - 1] = 0.0;
                goto L_2900;
L_2750:
                b[j - 1] = dtmp/a[kri - 1];
                /*  ADD THIS ONE TO THE LIST OF NONZEROES IN B */
                nzint = nzint + 1;
                LNRC(7,nzint - 1) = (short) j;

                l2 = LNRC(1,j - 1);
                if( l2 == 1 )
                        goto L_2900;
                kri = IP(1,j - 1);
                krl = l2 + kri - 1;
                /*  MARK ALL HIGHER ROWS THAT THIS VAR AFFECTS */
                for( k = kri + 1; k <= krl; k++ ){
                        k_ = k - 1;
                        i = IND(0,k_);
                        IP(0,i - 1) = -labs( IP(0,i - 1) );
                }
L_2900:
                ;
        }

        /*  STORE NO. NONZEROES IN THE RESULT */
        iluprm[19] = nzint;
        *nzinb = nzint;
        return;
#undef  LNRC
#undef  IP
#undef  IND
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS lslubt(long int n, long int nzdim, double anv[],
         short int *ind, long int *ip, short int *lnrc, double w[], double b[],
         long int iluprm[])
{
#define IND(I_,J_)      (*(ind+(I_)*(nzdim)+(J_)))
#define IP(I_,J_)       (*(ip+(I_)*(n)+(J_)))
#define LNRC(I_,J_)     (*(lnrc+(I_)*(n)+(J_)))
long int i, ii, ii_, j, k, k_, kri, krl;
double dtmp;

        /*  DO A BTRAN, I.E., MULTIPLY B() BY INVERSE OF TRANSPOSE OF U */

        /*  INPUTS:
         *    N = NO. OF ROWS/COLS
         *    W() = THE VECTOR TO BE TRAN-ED
         *    ANV() = NONZEROES OF LU INVERSE
         *    NZDIM = LENGTH OF THE ANV() VECTOR
         *    IP( I,1), IP( I,2) POINT TO START OF ROW/COLUMN I OF U.
         *    LNRC( I,1), LNRC( I,2) ARE LENGTHS OF ROW/COL I OF U.
         *    LNRC(.,3), LNRC(.,4) HOLD ROW/COL NUMBERS IN PIVOTAL ORDER.
         *    ILUPRM( 6) = NO. OF ELEMENTS IN L
         *  OUTPUTS:
         *    B() = W() * ( U INVERSE) * ( L INVERSE)
         *    W() = IS MODIFIED
         * */

        /*  DATE 31 DEC 1989
         *
         *  SET B() = W() * U INVERSE */
        for( ii = 1; ii <= n; ii++ ){
                ii_ = ii - 1;
                j = LNRC(2,ii_);
                i = LNRC(3,ii_);
                dtmp = w[i - 1];
                if( dtmp != 0.0 )
                        goto L_1300;
                b[j - 1] = 0.0;
                goto L_1500;
L_1300:
                kri = IP(0,j - 1);
                dtmp = dtmp/anv[kri - 1];
                b[j - 1] = dtmp;
                if( LNRC(0,j - 1) == 1 )
                        goto L_1500;
                krl = LNRC(0,j - 1) + kri - 1;
                for( k = kri + 1; k <= krl; k++ ){
                        k_ = k - 1;
                        i = IND(1,k_);
                        w[i - 1] = w[i - 1] - dtmp*anv[k_];
                }
L_1500:
                ;
        }

        /*  SET B() = B() * L INVERSE */
        kri = nzdim - iluprm[5] + 1;
        if( kri > nzdim )
                goto L_9000;
        for( k = kri; k <= nzdim; k++ ){
                k_ = k - 1;
                dtmp = b[IND(1,k_) - 1];
                if( dtmp == 0. )
                        goto L_1950;
                i = IND(0,k_);
                b[i - 1] = b[i - 1] + dtmp*anv[k_];
L_1950:
                ;
        }

L_9000:
        return;
#undef  IP
#undef  LNRC
#undef  IND
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS lslucc(long int n, long int nzdim, short int irc[],
         long int ip[], short int lnrc[], long int *last, long int *ntcmp)
{
long int iprve, j, j_, k, k_, kn, lntmp;


        /*  COMPRESS FILE OF POSITIVE INTEGERS.
         *
         *  INPUTS:
         *    N = NO. OF VECTORS, SOME OF 0 LENGTH
         *    IP( J) = START IN IRC() OF JTH VECTOR
         *    LNRC( J) = LENGTH OF JTH VECTOR, POSSIBLY WITH SIGN REVERSED
         *    IRC() = VECTOR TO BE COMPRESSED,
         *            NOTE! EMPTY ELS MUST = 0
         *    NZDIM = LENGTH OF IRC()
         *    NTCMP = NO. OF COMPRESSES DONE THIS FAR
         *    LAST = LAST OR HIGHEST POSITIVE ELEMENT IN IRC()
         *  OUTPUTS:
         *    LAST = LAST NONZERO ELEMENT IN COMPRESSED VECTOR
         *    IRC() = COMPRESSED SO NO INTERVENING ZEROES,
         *               ORDER OF NONZEROES IS UNCHANGED
         *    IP( J) = START IN IRC() OF JTH VECTOR
         *
         * */

        /*  DATE 20 MAY 1990
         *
         * */
        *ntcmp = *ntcmp + 1;
        /*  FLAG THE END OF EACH VECTOR WITH NEGATIVE OF ITS ID.... */
        for( j = 1; j <= n; j++ ){
                j_ = j - 1;
                lntmp = lnrc[j_];
                if( lntmp <= 0 )
                        goto L_1500;
                k = ip[j_] + lntmp - 1;
                /*  BUT 1ST STORE CONTENTS IN IP( J) */
                ip[j_] = irc[k - 1];
                irc[k - 1] = (short) -j;
L_1500:
                ;
        }

        /*  KN = HIGHEST FILLED POSITION THUS FAR IN COMPRESSED FILE */
        kn = 0;
        /*  IPRVE = WHERE PREVIOUS VECTOR ENDED IN COMPRESSED FILE */
        iprve = 0;
        /*  LOOP TO BRING NONZEROES DOWN... */
        for( k = 1; k <= *last; k++ ){
                k_ = k - 1;
                if( irc[k_] == 0 )
                        goto L_1900;
                kn = kn + 1;
                /*  LAST NONZERO IN CURRENT VECTOR? */
                if( irc[k_] >= 0 )
                        goto L_1800;
                /*  YES,  RESET IRC( K) AND REPOINT IP( J) */
                j = -irc[k_];
                irc[k_] = (short) ip[j - 1];
                ip[j - 1] = iprve + 1;
                iprve = kn;
L_1800:
                irc[kn - 1] = irc[k_];
L_1900:
                ;
        }
        *last = kn;
        return;
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS lslucr(long int n, long int nzdim, short int irc[],
         double anv[], long int ip[], short int lnrc[], long int *last,
         long int *ntcmp)
{
long int iprve, j, j_, k, k_, kn, lntmp;


        /*  COMPRESS FILE OF POSITIVE INTEGERS AND REALS.
         *
         *  INPUTS:
         *    N = NO. OF VECTORS, SOME PERHAPS OF 0 LENGTH
         *    IP( J) = START IN IRC() OF JTH VECTOR
         *    LNRC( J) = LENGTH OF JTH VECTOR
         *    IRC() = VECTOR TO BE COMPRESSED,
         *            NOTE!  EMPTY ELS MUST = 0
         *    ANV() = VECTOR OF REALS PARALLEL TO IRC()
         *    NZDIM = LENGTH OF IRC()
         *    NTCMP = NO. OF COMPRESSES DONE THUS FAR
         *    LAST = LAST OR HIGHEST NONZERO ELEMENT IN ANV()/IRC()
         *  OUTPUTS:
         *    LAST = LAST NONZERO ELEMENT IN COMPRESSED VECTOR
         *    IRC() = COMPRESSED SO NO INTERVENING ZEROES
         *    ANV() = COMPRESSED IN PARALLEL WITH IRC,
         *              ORDER OF NONZEROES IS UNCHANGED
         *    IP( J) = START IN RN() OF JTH VECTOR
         *
         * */

        /*  DATE 21 MAY 1990
         * */
        *ntcmp = *ntcmp + 1;
        /*  FLAG THE END OF EACH VECTOR WITH NEGATIVE OF ITS ID.... */
        for( j = 1; j <= n; j++ ){
                j_ = j - 1;
                lntmp = lnrc[j_];
                if( lntmp <= 0 )
                        goto L_1000;
                k = ip[j_] + lntmp - 1;
                /*  BUT 1ST STORE CONTENTS IN IP( J) */
                ip[j_] = irc[k - 1];
                irc[k - 1] = (short) -j;
L_1000:
                ;
        }

        /*  KN = HIGHEST FILLED POSITION THUS FAR IN COMPRESSED FILE */
        kn = 0;
        /*  IPRVE = WHERE PREVIOUS VECTOR ENDED IN COMPRESSED FILE */
        iprve = 0;
        /*  LOOP TO BRING NONZEROES DOWN... */
        for( k = 1; k <= *last; k++ ){
                k_ = k - 1;
                if( irc[k_] == 0 )
                        goto L_2100;
                kn = kn + 1;
                anv[kn - 1] = anv[k_];
                /*  LAST NONZERO IN THIS ROW? */
                if( irc[k_] >= 0 )
                        goto L_2000;
                /*  YES,  RESET IRC( K) AND REPOINT IP( J) */
                j = -irc[k_];
                irc[k_] = (short) ip[j - 1];
                ip[j - 1] = iprve + 1;
                iprve = kn;
L_2000:
                irc[kn - 1] = irc[k_];
L_2100:
                ;
        }
        *last = iprve;
        return;
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS lsortr(long int nc, long int nzero, double anv[],
         short int jcn[], long int iptr[], short int irn[], short int lnrc[])
{
long int i, i_, it1, it2, jt1, jt2, k, k1, k1_, k_, loc;
double at1, at2;

        /*  CONSTRUCT/SORT A ROW LIST OF NONZEROES FROM (I, J, A) TRIPLETS
         *   NOTE, IPTR CAN BE INTEGER*2 IF NO. NONZEROES < 32K */
        /*   THE FOLLOWING INTEGER*2 DECLARATION ASSUMES NC < 32K */
        /*c
         *jcp  DIMENSION JCN( NZERO), IRN( NZERO),  ANV( NZERO)
         *jcp  DIMENSION LNRC( NC), IPTR( NC) */

        /*  INPUTS:
         *    NC = NO. OF ROWS( AND COLS)
         *    NZERO = NO. OF NONZEROES IN ANV(), JCN(), AND IRN()
         *    ANV( K) = THE KTH NONZERO
         *    IRN( K) = ROW OF THE KTH NONZERO
         *    JCN( K) = COL OF THE KTH NONZERO
         *    LNRC( J) = NO. NONZEROES IN THE JTH ROW
         *  OUTPUTS:
         *    IPTR( I) = LOCN IN ANV() AND JCN() OF FIRST NONZERO IN ROW I
         *    ANV() = REORDERED BY ROW
         *    JCN() = REORDERED BY ROW,
         *         NOTE IRN IS UNALTERED
         * */

        /*  DATE 23 NOV 1989
         *
         *  POINT IPTR( I) TO 1 PAST WHERE LAST ENTRY IN ROW I WILL BE AT END */
        k = 1;
        for( i = 1; i <= nc; i++ ){
                i_ = i - 1;
                k = k + lnrc[i_];
                iptr[i_] = k;
        }

        /*  DURING 3900 LOOP INTERPRETATION OF JCN( K)  IS:
         *          > 0 MEANS K NEEDS TO BE PUT IN PLACE
         *          = 0 MEANS LOCN K IS EMPTY
         *          < 0 MEANS K IS ALREADY IN PLACE */
        for( k = 1; k <= nzero; k++ ){
                k_ = k - 1;
                /*  GET COL NO. OF THIS ELEMENT */
                jt1 = jcn[k_];
                /*  IS THIS GUY ALREADY IN PLACE? */
                if( jt1 < 0 )
                        goto L_3800;
                /*  NO,  FIND OUT HIS ROW */
                it1 = irn[k_];
                /*  REMOVE FROM CURRENT PLACE */
                at1 = anv[k_];
                /*  MARK THIS PLACE EMPTY */
                jcn[k_] = 0;

                /*  DO MUSICAL CHAIRS UNTIL WE FIND A VACANT CELL */
                for( k1 = 1; k1 <= nzero; k1++ ){
                        k1_ = k1 - 1;
                        /*  LOC IS WHERE CURRENT HOMELESS ELEMENT GOES */
                        loc = iptr[it1 - 1] - 1;
                        iptr[it1 - 1] = loc;
                        /*  EVICT CURRENT INHABITANTS */
                        jt2 = jcn[loc - 1];
                        at2 = anv[loc - 1];
                        it2 = irn[loc - 1];
                        /*  OLD HOMELESS GETS STORED */
                        anv[loc - 1] = at1;
                        jcn[loc - 1] = (short) -jt1;
                        /*  IF LOC WAS VACANT WE ARE DONE WITH THIS ROUND OF MUSICAL CHAIRS */
                        if( jt2 == 0 )
                                goto L_3800;
                        /*  NOT VACANT, NOW WE HAVE A NEW HOMELESS
                         *   THIS LOOP IS UNROLLED ONCE
                         *  LOC IS WHERE CURRENT HOMELESS ELEMENT GOES */
                        loc = iptr[it2 - 1] - 1;
                        iptr[it2 - 1] = loc;
                        /*  EVICT CURRENT INHABITANTS */
                        jt1 = jcn[loc - 1];
                        at1 = anv[loc - 1];
                        it1 = irn[loc - 1];
                        /*  OLD HOMELESS GETS STORED */
                        anv[loc - 1] = at2;
                        jcn[loc - 1] = (short) -jt2;
                        /*  IF LOC WAS VACANT WE ARE DONE WITH THIS ROUND OF MUSICAL CHAIRS */
                        if( jt1 == 0 )
                                goto L_3800;
                }

                /*  WE WILL NOT LOOK THIS LOW AGAIN, SO PUT IN CORRECT SIGN */
L_3800:
                jcn[k_] = (short) -jcn[k_];
        }
        return;
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS lslupv(long int n, long int nzdim, short int *ind,
         double anv[], long int *ip, short int *lnrc, double w[], long int ipvt,
         long int iluprm[], double dluprm[], long int *kalinv)
{
#define IND(I_,J_)      (*(ind+(I_)*(nzdim)+(J_)))
#define IP(I_,J_)       (*(ip+(I_)*(n)+(J_)))
#define LNRC(I_,J_)     (*(lnrc+(I_)*(n)+(J_)))
long int i, i_, ii, ii_, ij, ij_, im, in, ins, ipvro, ipvtl, ir, is,
         j, jpvcl, k, k1, k1_, k_, kci, kcl, kj, kk, kk_, kl, km, km_,
         knp, knp_, kp, kpl, kq, kr, kri, krl, ks, ks_, l, lbmp, lbmp1,
         lbmp2, lbmpt, lcol, lntmp, lowstl, lrow, m, m1, mxcmp, nfill,
         ntcmp, nzinl, nzinu, nzinv;
double au, dtmp, dtmp1, g, u, zerotl;



        /*  DO A PIVOT ...
         *
         *  INPUTS:
         *    W() = ENTERING COLUMN MULTIPLIED BY L INVERSE,
         *        IT IS DESTROYED.
         *    IPVT = COLUMN POSITION BEING REPLACED
         *          (THE PIVOT ROW TO THE OUTSIDE WORLD)
         *    ILUPRM( 19) = NO. NONZEROES IN W()
         *    LNRC( I, 7) = LIST OF ROWS WITH NONZERO IN W(),
         *       BOTH ILUPRM( 19) AND LNRC(,7) WERE PREPARED BY LSLUFT
         *    N = NO. OF ROWS/COLS IN U
         *    NZDIM = TOTAL SPACE AVAILABLE
         *    ILUPRM( 6) = NO. NONZEROES IN L INVERSE
         *    ILUPRM( 7) = NO. NONZEROES IN U
         *    ILUPRM( 8) = HIGHEST CELL USED FOR COL STORAGE IN IND(,2)
         *    ILUPRM( 9) = HIGHEST CELL USED FOR ROW STORAGE IN IND(,1),ANV()
         *    ILUPRM( 10) = NO. COMPRESSES DONE SINCE START OF CURRENT INVERSE
         *    ILUPRM( 11) = COMPRESSES ALLOWED BEFORE GIVING UP
         *    LNRC( I, 3/4) = ROW/COL IN POSITION I OF U INVERSE, TOP TO BOTTOM
         *    LNRC( I, 5) = POSITION IN U OF ROW I
         *    IND( K, 1) = ROW IN WHICH NONZERO K APPEARS
         *    IND( K, 2) = COL IN WHICH NONZERO K APPEARS
         *    IP( I, 1/2) = LOCATION OF FIRST NONZERO IN ROW/COL I OF U INVERSE
         *    DLUPRM( 1) = RELATIVE PIVOT SIZE TOLERANCE, E.G., .1
         *    DLUPRM( 2) = SET TO ZERO TOLERANCE, E.G., .1E-16
         *    DLUPRM( 4) = LARGEST ELEMENT IN ABSOLUTE VALUE IN INVERSE
         *
         *  OUTPUTS:
         *    DLUPRM( 4) = LARGEST ELEMENT IN ABSOLUTE VALUE IN INVERSE
         *    LNRC( I, 1/2) = NO. NONZEROES IN ROW/COL I
         *    LNRC( I, 3/4) = ROW/COL IN POSITION I OF U INVERSE, TOP TO BOTTOM
         *    LNRC( I, 5) = POSITION IN U OF ROW I
         *    IND( K, 1) = ROW IN WHICH NONZERO K APPEARS
         *    IND( K, 2) = COL IN WHICH NONZERO K APPEARS
         *    IP( I, 1/2) = LOCATION OF FIRST NONZERO IN ROW/COL I OF U INVERSE
         *    KALINV = 0  IF NO ERRORS
         *           = 7  IF MATRIX BECAME SINGULAR
         *           = 9  IF RAN OUT OF SPACE
         * */

        /*   REFERENCES:
         *    REID, J. K.(1982), "A SPARSITY-EXPLOITING VARIANT OF THE BARTELS-
         *    GOLUB DECOMPOSITION FOR LINEAR PROGRAMMING BASES", MATH. PROG.,
         *    VOL. 24, PP. 55-69.
         *
         *    REID, J. K,(1976), "FORTRAN SUBROUTINES FOR HANDLING SPARSE LINEAR
         *    PROGRAMMING BASES", REPORT AERE-R.8269, HARWELL.
         *
         *  DATE 4 JAN 1990
         *
         *  USE LOCAL VARIABLES TO STORE VARIOUS INVERSE PARAMETERS */
        lcol = iluprm[7];
        lrow = iluprm[8];
        nzinv = iluprm[18];
        u = dluprm[0];
        zerotl = dluprm[1];
        g = dluprm[3];
        *kalinv = 0;

        nzinl = iluprm[5];
        /*  HIGHEST CELL NOT OCCUPIED BY L */
        lowstl = nzdim - nzinl;
        ipvtl = ipvt;
        /*  ALLOW 20 COMPRESSES BEFORE ERROR EXIT BECAUSE OF SPACE */
        ntcmp = iluprm[9];
        mxcmp = ntcmp + 20;
        /*  DELETE OLD COLUMN FROM U */
        nzinu = iluprm[6] - LNRC(1,ipvtl - 1);
        /*  GET START OF OLD COLUMN */
        kci = IP(1,ipvtl - 1);
        /*   AND THE ROW IN WHICH IT IS BASIC */
        im = IND(0,kci - 1);
        /*   AND THE LENGTH OF OLD COLUMN - 1 */
        kcl = kci + LNRC(1,ipvtl - 1) - 1;
        /*   STRIKE IT OUT */
        LNRC(1,ipvtl - 1) = 0;
        /*  LOOP OVER NONZEROES OF OLD COLUMN */
        for( k = kci; k <= kcl; k++ ){
                k_ = k - 1;
                /*   GET NEXT ROW NO. IN OLD COLUMN */
                i = IND(0,k_);
                /*   STRIKE IT OUT */
                IND(0,k_) = 0;
                /*  GET START OF THIS ROW */
                kri = IP(0,i - 1);
                /*  NO. OF NONZEROES IN THIS ROW DECREASES BY 1 */
                lntmp = LNRC(0,i - 1) - 1;
                LNRC(0,i - 1) = (short) lntmp;
                /*  GET END OF THIS ROW */
                krl = kri + lntmp;
                /*   LOOP OVER THIS ROW, LOOKING FOR COL IPVTL */
                for( km = kri; km <= krl; km++ ){
                        km_ = km - 1;
                        /*    WE FOUND IT */
                        if( IND(1,km_) == ipvtl )
                                goto L_1580;
                }

                /*  OVERWRITE WITH LAST COL IN ROW */
L_1580:
                anv[km - 1] = anv[krl - 1];
                IND(1,km - 1) = IND(1,krl - 1);
                IND(1,krl - 1) = 0;
        }
        /*  END OF: DELETE OLD COLUMN
         *
         *  NOW INSERT NEW COLUMN INTO U
         *   LBMP WILL = LOWEST ROW IN BUMP WITH NONZERO IN W() */
        lbmp = 0;
        /*  SET M = PIVOT NO. OF PIVOT ROW */
        m = LNRC(4,im - 1);
        /*  NOW LOOP OVER NONZEROES OF W() */
        for( ii = 1; ii <= iluprm[18]; ii++ ){
                ii_ = ii - 1;
                /*  GET THE ROW */
                i = LNRC(6,ii_);
                if( fabs( w[i - 1] ) <= zerotl )
                        goto L_1690;
                if( lbmp < LNRC(4,i - 1) )
                        lbmp = LNRC(4,i - 1);
                /*  ANOTHER NONZERO IN U */
                nzinu = nzinu + 1;

                /*  APPEND A NONZERO TO COL IPVTL
                 *   GOT ENOUGH SPACE TO STORE NEW COL? */
                if( lcol < lowstl )
                        goto L_1630;
                /*  WILL COMPRESSING COLUMN FILE HELP? */
                if( ntcmp >= mxcmp || nzinu >= lowstl )
                        goto L_9200;
                lslucc( n, nzdim, ind, &IP(1,0), &LNRC(1,0), &lcol, &ntcmp );

                /*  APPEND THIS NONZERO */
L_1630:
                lcol = lcol + 1;
                lntmp = LNRC(1,ipvtl - 1);
                /*  POINT TO START IF THIS IS THE FIRST */
                if( lntmp == 0 )
                        IP(1,ipvtl - 1) = lcol;
                /*  ADJUST LENGTH OF THIS COL IN U */
                LNRC(1,ipvtl - 1) = (short) (lntmp + 1);
                /*   RECORD THE ROW OF THIS NONZERO */
                IND(0,lcol - 1) = (short) i;

                /*  APPEND A NONZERO TO ROW I
                 *  GET THE LENGTH OF THIS ROW */
                lntmp = LNRC(0,i - 1);
                /*  GET THE END OF THIS ROW + 1 */
                krl = IP(0,i - 1) + lntmp;
                if( krl > lrow )
                        goto L_1640;
                /*   IS THERE AN EMPTY SLOT AT END OF THIS ROW? */
                if( IND(1,krl - 1) == 0 )
                        goto L_1680;

                /*  NO, NEW ROW LIST WILL BE WRITTEN FOR THIS ROW.
                 *   DO WE HAVE ENOUGH ROOM IN THE ROW FILE? */
L_1640:
                if( lrow + lntmp < lowstl )
                        goto L_1650;
                /*   NO, COMPRESS ROW FILE IF IT WILL DO ANY GOOD. */
                if( ntcmp >= mxcmp || nzinu + lntmp >= lowstl )
                        goto L_9200;
                lslucr( n, nzdim, &IND(1,0), anv, ip, lnrc, &lrow, &ntcmp );

                /*  GET THE START OF THIS ROW */
L_1650:
                kp = IP(0,i - 1);
                /*  REWRITE IT AT THE END OF THE ROW FILE */
                IP(0,i - 1) = lrow + 1;
                if( lntmp == 0 )
                        goto L_1670;
                krl = kp + lntmp - 1;
                for( k = kp; k <= krl; k++ ){
                        k_ = k - 1;
                        lrow = lrow + 1;
                        anv[lrow - 1] = anv[k_];
                        IND(1,lrow - 1) = IND(1,k_);
                        IND(1,k_) = 0;
                }

L_1670:
                lrow = lrow + 1;
                krl = lrow;
                /* PLACE NEW ELEMENT AT END OF ROW.
                 *  NO. OF NONZEROES IN THIS ROW */
L_1680:
                LNRC(0,i - 1) = (short) (lntmp + 1);
                /*  STORE THE COEFFICIENT */
                anv[krl - 1] = w[i - 1];
                /*   AS WELL AS ITS COLUMN NO. */
                IND(1,krl - 1) = (short) ipvtl;
L_1690:
                w[i - 1] = 0.;
        }

        /*  NEW COLUMN HAS BEEN INSERTED
         *
         *  DID WE GET A SINGULAR MATRIX? */
        if( (LNRC(0,im - 1) == 0 || LNRC(1,ipvtl - 1) == 0) || m > lbmp )
                goto L_9100;

        /*  FIND SINGLETON COLS THAT CAN BE PUT ABOVE BUMP IN U
         *   ZERO OUT LNRC( , 7) TO BE USED AS A MARKER OF COLS IN BUMP */
        for( ii = m; ii <= lbmp; ii++ ){
                ii_ = ii - 1;
                j = LNRC(3,ii_);
                LNRC(6,j - 1) = 0;
        }

        /*  NON-SINGLETONS ARE MARKED WITH LNRC( J, 7) = 1
         *   ONLY LNRC( .,3) IS UPDATED
         *        LNRC(.,4)  IS USED TO HOLD ROWS STILL IN BUMP */
        ins = m;
        m1 = m;
        LNRC(6,ipvtl - 1) = 1;

        /*  LOOP OVER THE COLS IN THE BUMP FROM TOP TO BOTTOM */
        for( ii = m; ii <= lbmp; ii++ ){
                ii_ = ii - 1;
                i = LNRC(2,ii_);
                j = LNRC(3,ii_);
                /*  IF THIS COL HAD NO NONZERO IN EARLIER ROW OF BUMP */
                if( LNRC(6,j - 1) != 0 )
                        goto L_1750;
                /*  THEN IT IS A SINGLETON.  PUT ITS ROW ABOVE BUMP */
                LNRC(2,m1 - 1) = (short) i;
                m1 = m1 + 1;
                goto L_1770;
                /*   COL IS NOT A BUMP SINGLETON SO IT STAYS IN BUMP */
L_1750:
                kri = IP(0,i - 1);
                krl = kri + LNRC(0,i - 1) - 1;

                /*  MARK THE COLUMNS NOT SINGLETONS */
                for( k = kri; k <= krl; k++ ){
                        k_ = k - 1;
                        j = IND(1,k_);
                        LNRC(6,j - 1) = 1;
                }

                LNRC(3,ins - 1) = (short) i;
                ins = ins + 1;
L_1770:
                ;
        }

        /* PLACE NON-SINGLETONS IN NEW POSITION. */
        ij = m + 1;
        if( m1 >= lbmp )
                goto L_1790;
        lbmp1 = lbmp - 1;
        for( ii = m1; ii <= lbmp1; ii++ ){
                ii_ = ii - 1;
                LNRC(2,ii_) = LNRC(3,ij - 1);
                ij = ij + 1;
        }

        /*  PLACE SPIKE AT END. */
L_1790:
        LNRC(2,lbmp - 1) = (short) im;

        /*  FIND SINGLETON ROWS, EXCLUDING SPIKE, THAT CAN BE PUT BELOW BUMP
         *   NON-SINGLETONS ARE MARKED WITH LNRC( I, 7) = 2
         *   ONLY LNRC( , 3) IS UPDATED
         *        LNRC( , 4) STORES ROWS STILL IN BUMP
         *   ZERO OUT LNRC( , 7) TO MARK ROWS */
        for( ii = m1; ii <= lbmp; ii++ ){
                ii_ = ii - 1;
                i = LNRC(2,ii_);
                LNRC(6,i - 1) = 0;
        }

        lbmp1 = lbmp;
        lbmpt = lbmp;
        LNRC(6,im - 1) = 2;
        j = ipvtl;

        /*  PROCESS BUMP FROM BOTTOM TO TOP */
        for( ii = lbmp; ii >= m1; ii-- ){
                ii_ = ii - 1;
                i = LNRC(2,ii_);
                if( LNRC(6,i - 1) != 2 )
                        goto L_1870;
                /*  GOT A NONSINGLETON ROW */
                k = IP(0,i - 1);
                if( ii != lbmp )
                        j = IND(1,k - 1);
                kci = IP(1,j - 1);
                kcl = kci + LNRC(1,j - 1) - 1;
                LNRC(3,lbmpt - 1) = (short) i;
                lbmpt = lbmpt - 1;
                /*  LOOP OVER THIS COL, ANY ROW IN IT IS NOT A SINGLETON */
                for( k = kci; k <= kcl; k++ ){
                        k_ = k - 1;
                        i = IND(0,k_);
                        LNRC(6,i - 1) = 2;
                }
                goto L_1880;
                /*  GOT A SINGLETON ROW, PUT IT LOW IN BUMP */
L_1870:
                LNRC(2,lbmp1 - 1) = (short) i;
                lbmp1 = lbmp1 - 1;
L_1880:
                ;
        }

        /*  PUT BUMP ROWS BELOW SINGLETON COLS, ABOVE SINGLETON ROWS AND
         *   MARK WITH LNRC( I, 7) = 3 */
        for( ii = m1; ii <= lbmp1; ii++ ){
                ii_ = ii - 1;
                lbmpt = lbmpt + 1;
                i = LNRC(3,lbmpt - 1);
                LNRC(6,i - 1) = 3;
                LNRC(2,ii_) = (short) i;
        }

        /*  LOOK FOR SINGLETON ENTERING COLUMN
         *   LOOP OVER BUMP FROM TOP TO BOTTOM */
        for( ii = m1; ii <= lbmp1; ii++ ){
                ii_ = ii - 1;
                kci = IP(1,ipvtl - 1);
                kcl = kci + LNRC(1,ipvtl - 1) - 1;
                is = 0;
                for( k = kci; k <= kcl; k++ ){
                        k_ = k - 1;
                        l = IND(0,k_);
                        if( LNRC(6,l - 1) != 3 )
                                goto L_1910;
                        if( is != 0 )
                                goto L_1950;
                        is = 1;
                        i = l;
                        knp = k;
L_1910:
                        ;
                }

                if( is == 0 )
                        goto L_9100;
                /*  COL IPVTL IS A SINGLETON IN THE BUMP IN ROW I */
                IND(0,knp - 1) = IND(0,kci - 1);
                IND(0,kci - 1) = (short) i;
                kci = IP(0,i - 1);
                for( k = kci; k <= nzdim; k++ ){
                        k_ = k - 1;
                        if( IND(1,k_) == ipvtl )
                                goto L_1930;
                }
                /*  MOVE IPVTL TO THE FRONT OF ITS ROW */
L_1930:
                dtmp = anv[kci - 1];
                anv[kci - 1] = anv[k - 1];
                anv[k - 1] = dtmp;
                IND(1,k - 1) = IND(1,kci - 1);
                IND(1,kci - 1) = (short) ipvtl;
                ipvtl = IND(1,k - 1);
                LNRC(3,ii_) = (short) i;
                LNRC(6,i - 1) = 2;
        }

        ii = lbmp1;
        goto L_1970;
L_1950:
        in = m1;
        for( ij = ii; ij <= lbmp1; ij++ ){
                ij_ = ij - 1;
                LNRC(3,ij_) = LNRC(2,in - 1);
                in = in + 1;
        }
L_1970:
        lbmp2 = lbmp1 - 1;
        if( m1 == lbmp1 )
                goto L_5880;
        for( i = m1; i <= lbmp2; i++ ){
                i_ = i - 1;
                LNRC(2,i_) = LNRC(3,i_);
        }
        m1 = ii;
        if( m1 == lbmp1 )
                goto L_5880;

        /*  NOW PROCESS NON-REDUCIABLE BUMP, FROM M1 DOWN TO LBMP1 */
        ir = LNRC(2,lbmp1 - 1);
        for( ii = m1; ii <= lbmp1; ii++ ){
                ii_ = ii - 1;
                ipvro = LNRC(2,ii_);
                kp = IP(0,ipvro - 1);
                kr = IP(0,ir - 1);
                jpvcl = IND(1,kp - 1);
                if( ii == lbmp1 )
                        jpvcl = ipvtl;
                /*  IN NON-PIVOT ROW IR,
                 *   FIND PIVOT COL COEF... */
                krl = kr + LNRC(0,ir - 1) - 1;
                for( knp = kr; knp <= krl; knp++ ){
                        knp_ = knp - 1;
                        if( jpvcl == IND(1,knp_) )
                                goto L_2940;
                }
                if( ii < lbmp1 )
                        goto L_4900;
                goto L_9100;

                /*  BRING PIV COL COEF TO FRONT OF ROW */
L_2940:
                dtmp = anv[knp - 1];
                dtmp1 = anv[kr - 1];
                anv[knp - 1] = dtmp1;
                anv[kr - 1] = dtmp;
                IND(1,knp - 1) = IND(1,kr - 1);
                IND(1,kr - 1) = (short) jpvcl;
                /*  DECIDE WHICH ELEMENT IS THE PIVOT */
                if( ii == lbmp1 )
                        goto L_2950;
                /*  IS ANV( KP) DEFINITELY TOO SMALL? */
                if( fabs( anv[kp - 1] ) < u*fabs( dtmp ) )
                        goto L_2950;
                /*   NO, IS DTMP DEFINITELY TOO SMALL? */
                if( fabs( dtmp ) < u*fabs( anv[kp - 1] ) )
                        goto L_2970;
                /*   IS DTMP TOO SMALL RELATIVE TO CURRENT PIVOT IN ROW IR? */
                if( fabs( dtmp ) < u*dtmp1 )
                        goto L_2970;
                /*   NO, GO WITH THE SHORTER */
                if( LNRC(0,ipvro - 1) <= LNRC(0,ir - 1) )
                        goto L_2970;
                /*  SWITCH PIVOTS */
L_2950:
                LNRC(2,lbmp1 - 1) = (short) ipvro;
                LNRC(2,ii_) = (short) ir;
                ir = ipvro;
                ipvro = LNRC(2,ii_);
                k = kr;
                kr = kp;
                kp = k;
                kj = IP(1,jpvcl - 1);
                /*  BRING PIVOT TO FRONT OF COLUMN */
                for( k = kj; k <= nzdim; k++ ){
                        k_ = k - 1;
                        if( IND(0,k_) == ipvro )
                                goto L_2960;
                }
L_2960:
                IND(0,k - 1) = IND(0,kj - 1);
                IND(0,kj - 1) = (short) ipvro;

                /*  STAY WITH CURRENT PIVOT IN THE COLUMN */
L_2970:
                if( anv[kp - 1] == 0. )
                        goto L_9100;
                if( ii == lbmp1 )
                        goto L_4900;
                /*  KP = 1ST ELEMENT IN PIVOT ROW;  KR = ELEMENT TO BE ELIMINATED */
                dtmp = -anv[kr - 1]/anv[kp - 1];

                /*  IS THERE CLEARLY ENOUGH ROOM AT THE TOP FOR NEW ROW. */
                if( lrow + LNRC(0,ir - 1) + LNRC(0,ipvro - 1) <= lowstl )
                        goto L_2980;
                /*   NO, WILL COMPRESSION HELP? */
                if( ntcmp >= mxcmp || nzinu + LNRC(0,ir - 1) + LNRC(0,ipvro - 1) >
                 lowstl )
                        goto L_9200;
                /*   YES, DO A ROW COMPRESSION */
                lslucr( n, nzdim, &IND(1,0), anv, ip, lnrc, &lrow, &ntcmp );
                kp = IP(0,ipvro - 1);
                kr = IP(0,ir - 1);

L_2980:
                krl = kr + LNRC(0,ir - 1) - 1;
                kq = kp + 1;
                /*  COMPUTE MAX POSSIBLE NO. OF FILL ELEMENTS */
                nfill = LNRC(0,ipvro - 1) - 1;
                kpl = kp + nfill;
                if( nfill <= 0 )
                        goto L_3132;
                /*  LOAD PIVOT ROW, EXCLUDING EL IN PIV COL, INTO W */
                for( k = kq; k <= kpl; k++ ){
                        k_ = k - 1;
                        j = IND(1,k_);
                        w[j - 1] = anv[k_];
                }

                /*  ANY ELEMENTS LEFT IN ROW IR? */
                if( kr >= krl )
                        goto L_3125;
                /*  YES,
                 *  LOOP OVER NONPIVOT ROW IR TO COMPUTE NONFILL ELEMENTS */
                kj = krl;
                for( ks = kj; ks >= (kr + 1); ks-- ){
                        ks_ = ks - 1;
                        j = IND(1,ks_);
                        au = w[j - 1];
                        if( au == 0.0 )
                                goto L_3110;
                        nfill = nfill - 1;
                        au = anv[ks_] + dtmp*au;
                        /* IF ELEMENT IS VERY SMALL REMOVE IT FROM U. */
                        if( fabs( au ) <= zerotl )
                                goto L_3070;
                        g = fmax( g, fabs( au ) );
                        anv[ks_] = au;
                        goto L_3100;

                        /*  DELETE ELEMENT KS FROM ROW IR */
L_3070:
                        nzinu = nzinu - 1;
                        /* REMOVE ELEMENT FROM COL FILE. */
                        k = IP(1,j - 1);
                        kl = k + LNRC(1,j - 1) - 1;
                        LNRC(1,j - 1) = (short) (kl - k);
                        for( kk = k; kk <= kl; kk++ ){
                                kk_ = kk - 1;
                                if( IND(0,kk_) == ir )
                                        goto L_3090;
                        }
L_3090:
                        IND(0,kk - 1) = IND(0,kl - 1);
                        IND(0,kl - 1) = 0;
                        /*  DELETE ELEMENT FROM ROW FILE BY OVERWRITING IT */
                        anv[ks_] = anv[krl - 1];
                        IND(1,ks_) = IND(1,krl - 1);
                        IND(1,krl - 1) = 0;
                        krl = krl - 1;

L_3100:
                        w[j - 1] = 0.;
L_3110:
                        ;
                }

                /*  ANY FILLS LEFT TO DO? */
L_3125:
                if( nfill != 1 )
                        goto L_3130;
                /*  MAKE ROOM AT THE TOP FOR THE ONE NONZERO */
                anv[kr - 1] = anv[krl - 1];
                IND(1,kr - 1) = IND(1,krl - 1);
                IND(1,krl - 1) = 0;
                krl = krl - 1;
                goto L_3180;
L_3130:
                if( nfill > 0 )
                        goto L_3140;
                /*  NO FILLS TO BE GENERATED, ZERO OUT OLD FIRST ELEMENT */
L_3132:
                IND(1,kr - 1) = 0;
                /*   AND REPOINT START */
                IP(0,ir - 1) = kr + 1;
                goto L_3300;

                /*   MUST TRANSFER ROW TO MORE SPACE */
L_3140:
                IP(0,ir - 1) = lrow + 1;
                IND(1,kr - 1) = 0;
                kr = kr + 1;
                if( kr > krl )
                        goto L_3170;
                for( k = kr; k <= krl; k++ ){
                        k_ = k - 1;
                        lrow = lrow + 1;
                        anv[lrow - 1] = anv[k_];
                        IND(1,lrow - 1) = IND(1,k_);
                        IND(1,k_) = 0;
                }
L_3170:
                krl = lrow;

                /*  LOOP OVER PIVOT ROW TO GENERATE FILLS */
L_3180:
                if( kq > kpl )
                        goto L_3300;
                for( k1 = kq; k1 <= kpl; k1++ ){
                        k1_ = k1 - 1;
                        j = IND(1,k1_);
                        au = dtmp*w[j - 1];
                        if( fabs( au ) <= zerotl )
                                goto L_3289;
                        krl = krl + 1;
                        anv[krl - 1] = au;
                        IND(1,krl - 1) = (short) j;
                        nzinu = nzinu + 1;

                        /*  INSERT NEW NONZERO IN COLUMN FILE */
                        kci = IP(1,j - 1);
                        lntmp = LNRC(1,j - 1);
                        kcl = kci + lntmp - 1;
                        /*  TRY TO PUT IN EMPTY SLOT AT END */
                        if( kcl != lcol )
                                goto L_3220;
                        if( lcol >= lowstl )
                                goto L_3240;
                        lcol = lcol + 1;
                        goto L_3230;
L_3220:
                        if( IND(0,kcl) != 0 )
                                goto L_3240;
L_3230:
                        IND(0,kcl) = (short) ir;
                        goto L_3280;
                        /*  IS THERE ROOM AT TOP FOR NEW ENTRY? */
L_3240:
                        if( lcol + lntmp + 1 < lowstl )
                                goto L_3260;
                        /*   NO, WILL COMPRESSING COL FILE HELP? */
                        if( ntcmp >= mxcmp || nzinu + lntmp + 1 >= lowstl )
                                goto L_9200;
                        lslucc( n, nzdim, ind, &IP(1,0), &LNRC(1,0), &lcol, &ntcmp );
                        kci = IP(1,j - 1);
                        kcl = kci + lntmp - 1;
                        /*  COPY COL TO END OF COL SPACE */
L_3260:
                        IP(1,j - 1) = lcol + 1;
                        for( kk = kci; kk <= kcl; kk++ ){
                                kk_ = kk - 1;
                                lcol = lcol + 1;
                                IND(0,lcol - 1) = IND(0,kk_);
                                IND(0,kk_) = 0;
                        }
                        /*  AND FINALLY, ADD THE NEW NONZERO */
                        lcol = lcol + 1;
                        IND(0,lcol - 1) = (short) ir;
L_3280:
                        g = fmax( g, fabs( au ) );
                        LNRC(1,j - 1) = (short) (lntmp + 1);
L_3289:
                        w[j - 1] = 0.;
                }

L_3300:
                if( lrow < krl )
                        lrow = krl;
                /*  SET LENGTH OF ROW IR */
/* L_4810: unreferenced??? */
                LNRC(0,ir - 1) = (short) (krl - IP(0,ir - 1) + 1);

                /*  IS THERE ENOUGH ROOM TO STORE ELEMENT IN L? */
                if( lcol + 1 <= lowstl )
                        goto L_4830;
                /*   NO, WILL A COL COMPRESSION HELP? */
                if( ntcmp >= mxcmp )
                        goto L_9200;
                lslucc( n, nzdim, ind, &IP(1,0), &LNRC(1,0), &lcol, &ntcmp );

                /*  ADD ELEMENT TO L */
L_4830:
                anv[lowstl - 1] = dtmp;
                IND(0,lowstl - 1) = (short) ipvro;
                IND(1,lowstl - 1) = (short) ir;
                nzinl = nzinl + 1;
                lowstl = lowstl - 1;
                /*  REMOVE FROM PIVOT COL */
                lntmp = LNRC(1,jpvcl - 1) - 1;
                LNRC(1,jpvcl - 1) = (short) lntmp;
                kci = IP(1,jpvcl - 1);
                kcl = kci + lntmp;
                for( k = kci; k <= kcl; k++ ){
                        k_ = k - 1;
                        if( IND(0,k_) == ir )
                                goto L_4880;
                }

L_4880:
                IND(0,k - 1) = IND(0,kcl - 1);
                IND(0,kcl - 1) = 0;
                nzinu = nzinu - 1;
L_4900:
                ;
        }

        /*  SET:
         *    LNRC( II, 4) = COL IN PIVOT POSITION II
         *    LNRC( I, 5) = PIVOT POSITION OF ROW I */
L_5880:
        for( ii = m; ii <= lbmp; ii++ ){
                ii_ = ii - 1;
                i = LNRC(2,ii_);
                LNRC(4,i - 1) = (short) ii;
                k = IP(0,i - 1);
                j = IND(1,k - 1);
                LNRC(3,ii_) = (short) j;
        }

        /*  PUT LOCAL VARS BACK INTO INVERSE DATA OBJECT */
        iluprm[5] = nzinl;
        iluprm[6] = nzinu;
        iluprm[7] = lcol;
        iluprm[8] = lrow;
        iluprm[9] = ntcmp;
        iluprm[10] = mxcmp;
        dluprm[3] = g;
        goto L_9999;

        /*  ERROR RETURN CODES
         *   MATRIX BECAME SINGULAR */
L_9100:
        *kalinv = 7;
        goto L_9999;
        /*   RAN OUT OF SPACE */
L_9200:
        *kalinv = 9;
L_9999:
        return;
#undef  LNRC
#undef  IND
#undef  IP
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS lsdetm(float *determ, long int *idtrm, long int m,
         long int *ip, double anv[], long int ndim)
{
#define IP(I_,J_)       (*(ip+(I_)*(m)+(J_)))
long int i, i_;


        /*  COMPUTE DETERMINANT OF CURRENT INVERSE
         *
         *  INPUTS:
         *    M = NO. OF ROWS
         *    IP( I, 1) = ADDRESS OF FIRST(AND PIVOT) ELEMENT IN ROW I
         *    ANV() = VECTOR OF NONZEROES IN U
         *  OUTPUTS:
         *    DETERM = MANTISSA OF DETERMINANT
         *    IDTRM = EXPONENT OF 10 OF DETRMINANT
         *
         *  DATE 26 MAY 1990
         * */
        *determ = 1.0;
        for( i = 1; i <= m; i++ ){
                i_ = i - 1;
                *determ = (float)( *determ*anv[IP(0,i_) - 1]); /* warning */
L_3150:
                if( fabs( *determ ) < 100. )
                        goto L_3200;
                *idtrm = *idtrm + 2;
                *determ = *determ/ (float) 100.; /* warning */
                goto L_3150;
L_3175:
                *determ = *determ*(float)10.; /* warning */
                *idtrm = *idtrm - 1;
L_3200:
                if( fabs( *determ ) < .1 )
                        goto L_3175;
        }
        return;
#undef  IP
} /* end of function */

#include "lsinfo.h"

#define IOINFO  &(_info->io_info)
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)
/***********************************************************************
 **                 FILE:        QPINTR FORTRAN                        *
 **                 AUTHOR:      STUART SMITH                          *
 **                 LAST UPDATE: 04 MAR   1994                         *
 **                                                                    *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **      This file contains the following routines for LSGRG2, a       *
 **  large scale implementation of GRG2:                               *
 **                                                                    *
 **                                                                    *
 **  ROUTINE             DESCRIPTION                      LAST UPDATE  *
 ** ---------           -------------                    ------------- *
 **  XBTRAN     Performs a BTRAN of a vector and BINV     29 JAN 1991  *
 **  XDOT       Takes a dot product of a vector and a     07 JUN 1993  *
 **             column of the packed GRAD array                        *
 **  XSAXPY     Performs a SAXPY operation where X is     07 JUN 1993  *
 **             a column of the packed GRAD array                      *
 **  XFACT      Calls XFACTX to factor the basis matrix   04 MAR 1994  *
 **  XFTRAN     Performs a FTRAN of a vector and BINV     29 JAN 1991  *
 **  GETBAS     Retrieves the basis matrix from GRAD      04 MAR 1994  *
 **  XPIVOT     Pivots a new column into the basis        04 MAR 1994  *
 **  XERROR     A general error stop routine              29 JAN 1991  *
 **                                                                    *
 **  NOTES:    AS OF 6/7/93  Slacks are handled seperately and not     *
 **            stored in the GRAD array.                               *
 **                                                                    *
 **            As of 3/4/94  Obj Row Gradient Coefficients are not     *
 **            stored in the basis matrix to facilitate problems with  *
 **            all linear constraints.                                 *
 ***********************************************************************
 *
 * */
void /*FUNCTION*/ LSGRGCLASS xbtran(LsgrgInfo *_info,long int ibmap[], double zmem[], double row[],
         long int nb)
{
long int i, i_ /* , ka, kind, kip, klnrc, kw */;

/*----------------------------------------------------------------*/
/*  working storage for invert subsystem                          */
/*  define local pointers to alias global ones                    */
/*  lnrc and ind are defined as short in invert subsystem         */
/*  so cast global long pointers to short here                    */
/*----------------------------------------------------------------*/
  /*  DoubleArray   w1   = _info->inv_w1; */
  /*  DoubleArray   w2   = _info->inv_w2; */
      DoubleArray   w    = _info->inv_w;
      DoubleArray   a    = _info->inv_binv;
      IntegerArray  ip   = _info->inv_ip;
      Integer2Array lnrc = (short*) _info->inv_lnrc;
      Integer2Array ind  = (short*) _info->inv_ind;
/*----------------------------------------------------------------*/

        /*...................................................................
         *.................. COMMON DELCARATIONS ............................
         *...................................................................
         * */




        /*.....................................................................
         *................... ARGUMENT DECLARATION ............................
         *.....................................................................
         * */

        /*....................................................................
         *.....................LOCAL DECLARATIONS.............................
         *....................................................................
         * */

        /*......................................................................
         *     *****FUNCTIONS
         *     NONE
         *     *****SUBROUTINES CALLED
         *     LSLUBT,XERROR
         *     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****PURPOSE:
         *     THIS SUBROUTINE IMPLEMENTS THE BTRAN OR BACKWARD TRANSFORMATION
         *     ROW*(BASIS INVERSE)
         *
         *     *****ARGUMENT DESCRIPTION:
         *     SEE XMP DICTIONARY FOR ARGUMENTS NOT DESCRIBED HERE.
         *     ON INPUT:
         *     IBMAP    - POINTS TO BASIS INVERSE DATA STUCTURE MAP
         *     ZMEM     - ALLOCATED STORAGE
         *     NB       - ORDER OF BASIS MATRIX
         *     ROW      -  IS THE ROW TO BE BTRAN'ED.
         *
         *     ON OUTPUT:
         *     ROW      - CONTAINS (INPUT ROW)*(BASIS INVERSE)
         *
         *     *****APPLICATION AND USAGE RESTRICTIONS:
         *     THIS ROUTINE ASSUMES DATA STRUCTURES USED BY LSINVRT ROUTINES
         *     FOR LU DECOMPOSITION.
         *
         *     *****ALGORITHM NOTES:
         *
         *
         *     DATA STRUCTURE FOR THE PROBLEM DATA:
         *
         *     IBMAP(3)      POINTS TO THE LNRC INTEGER*2 ARRAY.
         *     IBMAP(4)      POINTS TO THE IP  INTEGER ARRAY.
         *     IBMAP(5)      POINTS TO THE W REAL*8 ARRAY.
         *     IBMAP(6)      POINTS TO THE IND INTEGER*2 ARRAY.
         *     IBMAP(7)      POINTS TO  A ARRAY.
         *
         *     *****REFERENCES:
         *     'FORTRAN SUBROUTINES FOR HANDLING SPARSE LINEAR PROGRAMMING BASES'
         *     J. K. REID, REPORT AERE-R8269, COMPUTER SCIENCE AND SYSTEMS
         *     DIVISION, AERE HARWELL, OXFORDSHIRE,ENGLAND.  JANUARY, 1976.
         *     'A SPARSITY-EXPLOITING VARIANT OF THE BARTELS-GOLUB
         *     DECOMPOSITION FOR LINEAR PROGRAMMING BASES', J. K. REID,
         *     MATHEMATICAL PROGRAMMING, VOL. 24, NO. 1, SEPTEMBER 1982.
         *     *****HISTORY:
         *     WRITTEN BY STUART H. SMITH, DEPARTMENT OF MECHANICAL ENGINEERING
         *                                 UNIVERSITY OF TEXAS, AUSTIN, TXI 78712
         *     DATE LAST MODIFIED:  JAN     31, 1988
         *     *****GENERAL
         *     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *     *****BODY OF PROGRAM  (XBTRAN)
         *
         *C----------------------------------------------------------------
         *C    GET WORK ARRAYS
         *C----------------------------------------------------------------
         *
         * */
/*
        klnrc = ibmap[2];
        kip = ibmap[3];
        kw = ibmap[4];
        kind = ibmap[5];
        ka = ibmap[6];
*/

        /*      -----------------------
         *     : PUT THE ROW INTO W
         *      -----------------------
         * */
        for( i = 1; i <= nb; i++ ){
                i_ = i - 1;
      /*          zmem[kw + i_ - 1] = row[i_]; */
                   w[i_] = row[i_];
        }


        /*      ---------------------------
         *     : CALL THE LSINVRT ROUTINE
         *      ---------------------------
         * */
/*
        lslubt( nb, _info->memory.lbinv, &zmem[ka - 1], (short*)&zmem[kind - 1],
         (long*)&zmem[kip - 1], (short*)&zmem[klnrc - 1], &zmem[kw - 1],
         row, _info->lsinvt.iluprm );
*/

        lslubt( nb, _info->memory.lbinv, a            , ind                    ,
         ip                   ,lnrc                     , w            ,
         row, _info->lsinvt.iluprm );

        /********
         *     NO ERROR CONDITIONS ARE SET BY LSLUBT
         *********
         * */
        return;
        /*     :::::::::: LAST CARD OF XBTRAN  :::::::::: */
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS xdot(LsgrgInfo *_info,double grad[], long int ihag[], long int iheg[],
         double row[], long int nb, long int icol, double *zdot)
{
long int collen, /* error , */  ifirst, ik, ik_, ikx, ilast;
double temp, zabs;



        /*....................................................................
         *........................COMMON DECLARATIONS ........................
         *....................................................................
         *
         * */







        /*......................................................................
         *.................... ARGUMENT DECLARATIONS ...........................
         *......................................................................
         * */

        /*......................................................................
         *.................... LOCAL DECLARATIONS ..............................
         *......................................................................
         * */

        /*......................................................................
         *     *****FUNCTIONS
         *     DABS
         *     *****SUBROUTINES CALLED
         *     XERROR
         *     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****PURPOSE:
         *     THE PURPOSE OF THIS ROUTINE IS TO PERFORM
         *     AN INNER PRODUCT BETWEEN A ROW VECTOR AND A
         *     PACKED MATRIX COLUMN.
         *     *****ARGUMENT DESCRIPTION:
         *     ON INPUT:
         *     GRAD        PACKED JACOBIAN MATRIX
         *     ICOL        COLUMN NUMBER IN PACKED MATRIX
         *     ROW         IS THE ROW VECTOR.
         *
         *     ON OUTPUT:
         *     ZDOT        CONTAINS THE INNER PRODUCT BETWEEN
         *                 ROW AND THE "ICOL" COLUMN IN THE PACKED MATRIX
         *
         *     *****APPLICATION AND USAGE RESTRICTIONS:
         *     *****ALGORITHM NOTES:
         *     *****REFERENCES:
         *     *****HISTORY:
         *     WRITTEN BY STUART H. SMITH,  DEPARTMENT OF MECHANICAL ENGINEERING
         *                UNIVERSITY OF TEXAS, AUSTIN, TX 85721.
         *     DATE LAST MODIFIED: APRIL 06, 1988
         *     *****GENERAL
         *     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****BODY OF PROGRAM  (XDOT)
         *
         * */
        *zdot = 0.0;
        zabs = 0.0;

        if( icol > _info->dimen.n ){
                *zdot = -row[icol - _info->dimen.n - 1];
                zabs = fabs( *zdot );
        } else{
                ifirst = iheg[icol - 1];
                ilast = iheg[icol] - 1;
                collen = ilast - ifirst + 1;
                if( collen >= 1 && collen <= _info->cmax.colmax ){
                        for( ik = ifirst; ik <= ilast; ik++ ){
                                ik_ = ik - 1;
                                ikx = ihag[ik_];
                                temp = row[ikx - 1]*grad[ik_];
                                *zdot = *zdot + temp;
                                zabs = zabs + fabs( temp );
                        }
                } else if( collen != 0 ){
#ifdef IO_ENABLED
                        sprintf( LSGRG_MSGBUFFER,
                         " XDOT...ERROR: COLUMN%6ld HAS LENGTH%6ld BUT"
                         " COLMAX = %6ld\n",
                         icol, collen, _info->cmax.colmax );
                         lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
             /*         error = 0; */
            /*          xerror( error, iounit.ioerr ); */
                        lsgrg_errorexit(JUMPBUF,_LSGRG_XDOT_COLLEN);
                }
        }

        if( fabs( *zdot ) < _info->tols.eps || fabs( *zdot ) < _info->tols.eps*zabs )
                *zdot = 0.0;


        return;
        /* :::::::::::::::::::::::::::::::::::  LAST CARD OF XDOT  :::::::::: */
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS xsaxpy(LsgrgInfo *_info,double grad[], long int ihag[], long int iheg[],
         double col[], double a, long int nb, long int icol)
{
long int collen, /* error,*/  ifirst, ik, ik_, ikx, ilast;



        /*......................................................................
         *..................... COMMON DELCARATIONS ............................
         *......................................................................
         * */





        /*......................................................................
         *...................... ARGUMENT DECLARATIONS .........................
         *......................................................................
         * */

        /*......................................................................
         *...................... LOCAL DECLARATIONS ............................
         *......................................................................
         * */

        /*......................................................................
         *     *****FUNCTIONS
         *     NONE
         *     *****SUBROUTINES CALLED
         *     XERROR
         *     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****PURPOSE:
         *     THE PURPOSE OF THIS ROUTINE IS TO PERFORM
         *     THE VECTOR EQUATION Y = A * X + Y WHERE A IS A CONSTANT AND X IS A
         *     PACKED MATRIX COLUMN.
         *     *****ARGUMENT DESCRIPTION:
         *     SEE XMP DICTIONARY FOR ARGUMENTS NOT DESCRIBED HERE.
         *     ON INPUT:
         *     GRAD        IS THE PACKED JACOBIAN MATRIX
         *     ICOL        IS THE NUMBER OF THE PACKED COLUMN IN GRAD
         *     COL         IS THE COL VECTOR.
         *
         *     ON OUTPUT:
         *     COL         CONTAINS THE UPDATED VECTOR A * X + Y
         *
         *     *****APPLICATION AND USAGE RESTRICTIONS:
         *     *****ALGORITHM NOTES:
         *     *****REFERENCES:
         *     *****HISTORY:
         *     WRITTEN BY STUART H. SMITH,  DEPARTMENT OF MECHANICAL ENGINEERING
         *                UNIVERSITY OF TEXAS, AUSTIN, TX 85721.
         *     DATE LAST MODIFIED: APRIL 06, 1988
         *     *****GENERAL
         *     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****BODY OF PROGRAM  (XSAXPY)
         * */
        if( icol > _info->dimen.n ){
                col[icol - _info->dimen.n - 1] = col[icol - _info->dimen.n - 1] - a;
        } else{
                ifirst = iheg[icol - 1];
                ilast = iheg[icol] - 1;
                collen = ilast - ifirst + 1;
                if( collen >= 1 && collen <= _info->cmax.colmax ){
                        for( ik = ifirst; ik <= ilast; ik++ ){
                                ik_ = ik - 1;
                                ikx = ihag[ik_];
                                col[ikx - 1] = col[ikx - 1] + a*grad[ik_];
                        }
                } else if( collen != 0 ){
#ifdef IO_ENABLED
                        sprintf( LSGRG_MSGBUFFER,
                         " XSAXPY...ERROR: COLUMN%6ld HAS LENGTH%6ld"
                         " BUT COLMAX = \n",
                         icol, collen );
                        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
            /*          error = 0; */
            /*          xerror( error, iounit.ioerr ); */
                        lsgrg_errorexit(JUMPBUF,_LSGRG_XSAXPY_COLLEN);
                }
        }
/* L_999: unreferenced??? */
        return;
        /* :::::::::::::::::::::::::::::::::::  LAST CARD OF XSAXPY  :::::::::: */
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS xfact(LsgrgInfo *_info,double grad[], long int ihag[], long int iheg[],
         long int ibmap[], double zmem[], long int ibv[], LOGICAL32 bschng,
         long int *rcode, double cdnum[], long int irank[], LOGICAL32 getcd)
{
long int i, iflag, /* ka,  */ kdropd, keepsq, /* kind, kip, klnrc, kw, kw1,
         kw2 , */ nz;


/*----------------------------------------------------------------*/
/*  working storage for invert subsystem                          */
/*  define local pointers to alias global ones                    */
/*  lnrc and ind are defined as short in invert subsystem         */
/*  so cast global long pointers to short here                    */
/*----------------------------------------------------------------*/
      DoubleArray   w1   = _info->inv_w1;
      DoubleArray   w2   = _info->inv_w2;
      DoubleArray   w    = _info->inv_w;
      DoubleArray   a    = _info->inv_binv;
      IntegerArray  ip   = _info->inv_ip;
      Integer2Array lnrc = (short*) _info->inv_lnrc;
      Integer2Array ind  = (short*) _info->inv_ind;
/*----------------------------------------------------------------*/
        /*......................................................................
         *.................... ARGUMENT DECLARATIONS ...........................
         *......................................................................
         *
         *     DOUBLE PRECISION ZMEM(LMEM),CDNUM(NB)
         *     INTEGER  IBMAP(LBMAP),IHAG(LGRAD),IHEG(NP1),IRANK(NB)
         *     INTEGER  IBV(NB),RCODE
         *     LOGICAL BSCHNG,GETCD
         *     DOUBLE PRECISION GRAD(LGRAD)
         *
         * */
        /*......................................................................
         *.................... LOCAL DECLARATIONS ..............................
         *......................................................................
         * */

        /*......................................................................
         *     *****SUBROUTINES CALLED
         *     GETBAS,CONDNM,XERROR
         *     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****PURPOSE:
         *     THIS SUBROUTINE IS AN INTERFACE ROUTINE WHICH ALLOWS ACCESS TO
         *     THE LSINVERT DATA STRUCTURES.  FIRST, GETBAS IS CALLED TO PLACE
         *     THE CURRENT BASIS INTO THE PROPER DATA STRUCURES.  THEN LSLUFC IS
         *     CALLED TO FACTOR THE BASIS.  TH OPTION TO SELECT NEW PIVOTS OR TO
         *     USE PREVIOUS PIVOTS IS DECIDED BY THE LOGICAL VARIABLE "BSCHNG"
         *
         *     *****COMMON DESCRIPTION:
         *
         *     *****ARGUMENT DESCRIPTION:
         *     ON INPUT:
         *
         *     IBMAP    - BASIS INVERSE WORK ARRAY MAP
         *     NB       - NUMBER OF BASIC VARIABLES. (in COMMON)
         *     BSCHNG   - LOGICAL INDICATING WHETHER THE BASIS HAS CHANGED
         *     GRAD     - PACKED JACOBIAN
         *     IHAG     - ROW NUMBERS FOR GRAD
         *     IHEG     - COLUMN POINTER ARRAY
         *     IBV      - LIST OF BASIC VARIABLE INDECIES
         *     GETCD    - LOGICAL INDICATING WHETHER TO CALCULATE BASIS CONDITION
         *
         *     ON OUTPUT:
         *     CDNUM    - CONDITION NUMBER ESTIMATE
         *     IRANK    - RANK OF EACH DIAGONAL BLOCK
         *     RCODE    - IS A RETURN CODE FOR XFACT:
         *                 1 MEANS EVERYTHING OK;
         *                 2 MEANS THAT THE BASIS IS SINGULAR;
         *                 3 MEANS STORAGE OVERFLOW.
         *                 4 MEANS OTHER ERROR CONDITIONS
         *                 5 MEANS BASIS IS ILL-CONDITIONED
         *
         *     THE CURRENT BASIS HAS BEEN FACTORED BY LSINVERT
         *
         *     *****APPLICATION AND USAGE RESTRICTIONS:
         *     THIS SUBROUTINE IS COMPATIBLE WITH THE LSINVERT
         *     ROUTINES FOR HANDLING AN LU FACTORIZATION OF THE BASIS MATRIX.
         *     THE 'BODY OF PROGRAM' MUST BE RE-WRITTEN IF SOME OTHER INVERSE
         *     REPRESENTATION IS USED.
         *
         *     *****ALGORITHM NOTES:
         *     DATA STRUCTURE FOR THE PROBLEM DATA:
         *
         *     IBMAP( 1)  -   POINTS TO THE W1 WORK ARRAY (COND)
         *     IBMAP( 2)  -   POINTS TO THE W2 WORK ARRAY (COND)
         *     IBMAP( 3)  -   POINTS TO THE LRNC INTEGER*2 ARRAY.
         *     IBMAP( 4)  -   POINTS TO THE IP INTEGER ARRAY.
         *     IBMAP( 5)  -   POINTS TO THE W REAL ARRAY.
         *     IBMAP( 6)  -   POINTS TO THE IND ARRAY.
         *     IBMAP( 7)  -   POINTS TO THE A (BINV) ARRAY.
         *
         *
         *     *****REFERENCES:
         *     'FORTRAN SUBROUTINES FOR HANDLING SPARSE LINEAR PROGRAMMING BASES'
         *     BY J. K. REID, REPORT AERE-R8269, COMPUTER SCIENCE AND SYSTEMS
         *     DIVISION, AERE HARWELL, OXFORDSHIRE, ENGLAND. JANUARY, 1976.
         *     'A SPARSITY-EXPLOITING VARIANT OF THE BARTELS-GOLUB
         *     DECOMPOSITION FOR LINEAR PROGRAMMING BASES', J. K. REID,
         *     MATHEMATICAL PROGRAMMING, VOL. 24, NO. 1, SEPTEMBER 1982.
         *     *****HISTORY:
         *     WRITTEN BY:  STUART H. SMITH  DEPARTMENT OF MECHANICAL ENGINEERING
         *                  UNIVERSITY OF TEXAS, AUSTIN, TEXAS 78712.
         *     DATE LAST MODIFIED:  SEPT  12 1988
         *     *****GENERAL
         *     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *     *****BODY OF PROGRAM  (XFACT)
         * */

        /*C-----------------------------------------------------
         *C    GET POINTERS TO WORK ARRAYS
         *C----------------------------------------------------- */
/*
        kw1 = ibmap[0];
        kw2 = ibmap[1];
        klnrc = ibmap[2];
        kip = ibmap[3];
        kw = ibmap[4];
        kind = ibmap[5];
        ka = ibmap[6];
*/
        nz = _info->nzerob.nzbas;

        /*C     ----------------------------------------------
         *C     | CALL GETBAS TO FILL IN DATA STRUCURES
         *C     ----------------------------------------------
         * */
#ifdef IO_ENABLED
        if( _info->dbug.debug || _info->nintbk.ipr > 5 )
                {
                sprintf(LSGRG_MSGBUFFER,
                 "\n\n***ENTERING SUBROUTINE: XFACT  -  NB = %5ld  BSCHNG = %3c\n",
                 _info->nintbk.nb, TorF(bschng) );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
/*
        getbas(_info, grad, ihag, iheg, &zmem[ka - 1], (short*)&zmem[kind - 1],
         ibv, _info->dimen.n, &_info->nzerob.nzbas );
*/
        getbas(_info, grad, ihag, iheg, a            , ind                    ,
         ibv, _info->dimen.n, &_info->nzerob.nzbas );


        /*      ------------------------------------------
         *     | CALL LSLUFC TO FACTOR BASIS
         *      ------------------------------------------
         * */
        keepsq = 0;
        if( !bschng ){
                keepsq = 1;
        }

/*
        lslufc( &zmem[ka - 1], (short*)&zmem[kind - 1], _info->nzerob.nzbas,
         _info->memory.lbinv, _info->nintbk.nb, (long*)&zmem[kip - 1], (short*)&zmem[klnrc - 1],
         &zmem[kw - 1], &kdropd, keepsq, _info->lsinvt.iluprm, _info->lsinvt.dluprm,
         &iflag );
*/

        lslufc( a            , ind                    , _info->nzerob.nzbas,
         _info->memory.lbinv, _info->nintbk.nb, ip                   , lnrc               ,
         w            , &kdropd, keepsq, _info->lsinvt.iluprm, _info->lsinvt.dluprm,
         &iflag );

        if( kdropd != 0 )
                iflag = 7;

#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 5 ){
                sprintf(LSGRG_MSGBUFFER, " XFACT... RETURN FROM LSLUFC...IFLAG = %4ld\n",
                 iflag );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                sprintf(LSGRG_MSGBUFFER, " FACTORIZATION COMPLETE -- STABILITY COMMON BLOCK PARAMETERS ARE:\n W(1) = %15.8e\tW(2) = %15.8e\n",
                 _info->lsinvt.dluprm[2], _info->lsinvt.dluprm[3] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        if( getcd && iflag == 0 ){
/*
                condnm(_info, _info->nintbk.nb, (short*)&zmem[kind - 1], &zmem[ka - 1],
                 _info->memory.lbinv, (short*)&zmem[klnrc - 1], (long*)&zmem[kip - 1],
                 _info->lsinvt.iluprm, &zmem[kw - 1], irank, cdnum, grad, ihag, iheg,
                 ibv, &zmem[kw1 - 1], &zmem[kw2 - 1] );
*/
                condnm(_info, _info->nintbk.nb, ind                    , a            ,
                 _info->memory.lbinv, lnrc                    , ip                   ,
                 _info->lsinvt.iluprm, w            , irank, cdnum, grad, ihag, iheg,
                 ibv, w1            , w2             );

                if( _info->zblck.condmx > _info->zcond.cndtol )
                        iflag = -20;
        }


        /*C---------------------------------------------------
         *C    ERROR RETURNS
         *C---------------------------------------------------
         * */
        i = labs( iflag );

        /*#######################
         *    CHANGED 7/12/90  --- INVERT ROUTINES RETURN 0 IF OK  (SHS)
         *#######################
         * */
        if( iflag == 0 ){
                *rcode = 1;
        } else if( i == 7 ){
           if( _info->nintbk.ipr > 2 )
               lsgrg_error_msg(IOINFO,
                 "\n XFACT... CURRENT BASIS MATRIX IS SINGULAR \n" );

                *rcode = 2;
        } else if( iflag == -20 ){
                *rcode = 5;
        } else if( iflag == 9 ){
                lsgrg_msg(IOINFO,
                 "\n XFACT... INSUFFICIENT SPACE FOR BASIS INVERSE \n" );
                *rcode = 3;
        } else{
                *rcode = 4;
        }

        if( _info->nintbk.ipr > 4 )
            lsgrg_msg(IOINFO," ***EXITING SUBROUTINE: XFACT\n" );

        return;
        /* :::::::::::::::::::::::::::::::::::  LAST CARD OF XFACT  :::::::::: */
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS xftran(LsgrgInfo *_info,long int ibmap[], double zmem[], double column[],
         long int nb)
{
long int i, i_, /* ka, kind, kip, klnrc, kw, */ nzinb, nzinw;

/*----------------------------------------------------------------*/
/*  working storage for invert subsystem                          */
/*  define local pointers to alias global ones                    */
/*  lnrc and ind are defined as short in invert subsystem         */
/*  so cast global long pointers to short here                    */
/*----------------------------------------------------------------*/
  /*  DoubleArray   w1   = _info->inv_w1; */
  /*  DoubleArray   w2   = _info->inv_w2; */
      DoubleArray   w    = _info->inv_w;
      DoubleArray   a    = _info->inv_binv;
      IntegerArray  ip   = _info->inv_ip;
      Integer2Array lnrc = (short*) _info->inv_lnrc;
      Integer2Array ind  = (short*) _info->inv_ind;
/*----------------------------------------------------------------*/
        /*......................................................................
         *     *****FUNCTIONS
         *     NONE
         *     *****SUBROUTINES CALLED
         *     LSLUFT,XERROR
         *     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****PURPOSE:
         *     THIS SUBROUTINE IMPLEMENTS THE FTRAN OR FORWARD TRANSFORMATION
         *     (BASIS INVERSE)*COLUMN
         *
         *     *****ARGUMENT DESCRIPTION:
         *     SEE XMP DICTIONARY FOR ARGUMENTS NOT DESCRIBED HERE.
         *     ON INPUT:
         *     IBMAP    - POINTS TO BASIS INVERSE DATA STUCTURE MAP
         *     ZMEM     - ALLOCATED STORAGE
         *     NB       - ORDER OF BASIS MATRIX
         *     COLUMN   - IS THE COLUMN TO BE FTRAN'ED.
         *
         *     ON OUTPUT:
         *     COLUMN   - CONTAINS (BASIS INVERSE)*(INPUT COLUMN)
         *
         *     *****APPLICATION AND USAGE RESTRICTIONS:
         *     THIS ROUTINE ASSUMES DATA STRUCTURES USED BY THE LSINVERT
         *     ROUTINES FOR LU DECOMPOSITION.
         *
         *     *****ALGORITHM NOTES:
         *
         *     DATA STRUCTURE FOR THE PROBLEM DATA:
         *
         *     IBMAP(3)      POINTS TO THE LNRC INTEGER*2 ARRAY.
         *     IBMAP(4)      POINTS TO THE IP INTEGER ARRAY.
         *     IBMAP(5)      POINTS TO THE W REAL*8ARRAY.
         *     IBMAP(6)      POINTS TO THE IND INTEGER*2 ARRAY.
         *     IBMAP(7)      POINTS TO  A ARRAY.
         *
         *     *****REFERENCES:
         *     'FORTRAN SUBROUTINES FOR HANDLING SPARSE LINEAR PROGRAMMING BASES'
         *     J. K. REID, REPORT AERE-R8269, COMPUTER SCIENCE AND SYSTEMS
         *     DIVISION, AERE HARWELL, OXFORDSHIRE,ENGLAND.  JANUARY, 1976.
         *     'A SPARSITY-EXPLOITING VARIANT OF THE BARTELS-GOLUB
         *     DECOMPOSITION FOR LINEAR PROGRAMMING BASES', J. K. REID,
         *     MATHEMATICAL PROGRAMMING, VOL. 24, NO. 1, SEPTEMBER 1982.
         *     *****HISTORY:
         *     WRITTEN BY STUART H. SMITH, DEPARTMENT OF MECHANICAL ENGINEERING
         *                                 UNIVERSITY OF TEXAS, AUSTIN, TXI 78712
         *     DATE LAST MODIFIED:  JAN  31, 1990
         *     *****GENERAL
         *     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *     *****BODY OF PROGRAM  (XFTRAN)
         *
         *C----------------------
         *C    GET WORK ARRAYS
         *C----------------------
         *
         *
         * */
/*
        klnrc = ibmap[2];
        kip = ibmap[3];
        kw = ibmap[4];
        kind = ibmap[5];
        ka = ibmap[6];
*/

        /*      ------------------------------------------
         *     : COPY COLUMN INTO THE W ARRAY
         *     | LET LSLUFT COUNT THE NUMBER OF NONZEROS
         *      ------------------------------------------
         * */
        for( i = 1; i <= nb; i++ ){
                i_ = i - 1;
          /*    zmem[kw + i_ - 1] = column[i_]; */
                   w[i_]           = column[i_];

        }
        nzinw = 0;

        /*      ---------------------------------
         *     : CALL THE LSINVERT FTRAN ROUTINE
         *      ---------------------------------
         * */
/*
        lsluft( nb, _info->memory.lbinv, &zmem[ka - 1], (short*)&zmem[kind - 1],
         (long*)&zmem[kip - 1], (short*)&zmem[klnrc - 1], &zmem[kw - 1],
         column, &nzinw, &nzinb, _info->lsinvt.iluprm );
*/
        lsluft( nb, _info->memory.lbinv, a            , ind                    ,
         ip                   , lnrc                    , w            ,
         column, &nzinw, &nzinb, _info->lsinvt.iluprm );

        /*******
         *     NO ERROR CONDITIONS ARE SET BY LSLUFT
         ********
         * */
        return;
        /*     :::::::::: LAST CARD OF XFTRAN  :::::::::: */
} /* end of function */


/*C
 * */
void /*FUNCTION*/ LSGRGCLASS getbas(LsgrgInfo *_info,double grad[], long int ihag[], long int iheg[],
         double binv[], short int *ind, long int ibv[], long int n, long int *nzbas)
{
#define IND(I_,J_)      (*(ind+(I_)*(_info->memory.lbinv)+(J_)))
long int collen, /* error, */  ictr, ictr_, iend, istart, j, j_, jcol, nz;

        /*......................................................................
         *     *****FUNCTIONS
         *     MIN0
         *     *****SUBROUTINES CALLED
         *     XERROR
         *
         *     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****PURPOSE:
         *     THIS ROUTINE RETRIEVES THE BASIS FROM THE PACKED GRAD ARRAY AND
         *     PLACES IT INTO THE DATA STRUCURES REQUIRED BY THE LSINVERT
         *     ROUTINES IN PREPARATION FOR A CALL TO LSLUFC .
         *     NOTE:  THE OBJECTIVE ROW COEFFICIENTS ARE NOT STORED FOR THE
         *            STRUCTURAL VARIABLES.  ONLY THE SLACK ELEMENT IS STORED.
         *
         *     *****ARGUMENT DESCRIPTION:
         *     SEE XMP DICTIONARY FOR ARGUMENTS NOT DESCRIBED HERE.
         *     ON INPUT:
         *     GRAD        IS PACKED JACONIAN ARRAY
         *     IHAG,IHEG   POINTER ARRAYS FOR PACKED STRUCTURE
         *     IBV         BASIC VARIABLE INDICIES
         *
         *     ON OUTPUT:
         *     BINV        BASIS IN PACKED FORM
         *     IND(*,1)    ROW INDICES FOR PACKED BASIS
         *     IND(*,2)    COLUMN INDICES FOR PACKED BASIS
         *     NZBAS       NUMBER OF NONZEROS IN CURRENT BASIS
         *
         *     *****APPLICATION AND USAGE RESTRICTIONS:
         *     THIS SUBROUTINE IS COMPATIBLE WITH THE LSINVERT ROUTINES
         *     FOR HANDLING AN LU FACTORIZATION OF THE BASIS MATRIX.
         *
         *     *****ALGORITHM NOTES:
         *     *****REFERENCES:
         *     'FORTRAN SUBROUTINES FOR HANDLING SPARSE LINEAR PROGRAMMING BASES'
         *     J. K. REID, REPORT AERE-R8269, COMPUTER SCIENCE AND SYSTEMS
         *     DIVISION, AERE HARWELL, OXFORDSHIRE,ENGLAND.  JANUARY, 1976.
         *     'A SPARSITY-EXPLOITING VARIANT OF THE BARTELS-GOLUB
         *     DECOMPOSITION FOR LINEAR PROGRAMMING BASES', J. K. REID,
         *     MATHEMATICAL PROGRAMMING, VOL. 24, NO. 1, SEPTEMBER 1982.
         *     *****HISTORY:
         *     WRITTEN BY: STUART H. SMITH, DEPARTMENT OF MECHANICAL ENGINEERING
         *                 UNIVERSITY OF TEXAS, AUSTIN, TEXAS 78712.
         *     DATE LAST MODIFIED: MARCH  04, 1994
         *
         *         3/4/94  -- REMOVED OBJECTIVE ROW COEFFICIENTS FROM BASIS
         *                    THE OBJECTIVE ROW IS SIMPLY A UNIT ROW.
         *     7/4/96 -- jcp modified to work with expanded jacobian
         *        (merge with changes to process new nonzero entries
         *         in getgr)
         *     *****GENERAL
         *     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****BODY OF PROGRAM  (GETBAS)
         * */
        /*--------------------------------------------------------
         *     GET START TIME
         *--------------------------------------------------------
         *
         *
         *
         *------------------------------------------------------------------
         *     FIND MINIMUM ARRAY LENGTH TO CHECK FOR OVERFLOW
         *------------------------------------------------------------------
         * */
#ifdef IO_ENABLED
       if( _info->dbug.debug ) {
            sprintf(LSGRG_MSGBUFFER,
             "\n*** ENTERING SUBROUTINE: GETBAS  -  NB = %6ld\n",
                 _info->nintbk.nb );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
           }
#endif

        nz = 0;

        if( _info->dbug.debug )
            lsgrg_msg(IOINFO,"\n     THE BASIS MATRIX IS:\n\n" );

        for( j = 1; j <= _info->nintbk.nb; j++ ){
                j_ = j - 1;
                jcol = ibv[j_];
                if( jcol > n ){
                        nz = nz + 1;
                        if( nz > _info->memory.lbinv )
                                goto L_140;
                        binv[nz - 1] = -1.0e0;
                        IND(0,nz - 1) = (short) (jcol - n);
                        IND(1,nz - 1) = (short)  j;
#ifdef IO_ENABLED
                        if( _info->dbug.debug )
                                {
                                sprintf(LSGRG_MSGBUFFER,
                                 "\n ...COLUMN = %6ld  VARIABLE = %6ld  SLACK #  %6ld\n",
                                 j, jcol, jcol - n );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                } else{
                        istart = iheg[jcol - 1];
                        iend = iheg[jcol] - 1;
                        collen = iend - istart + 1;
#ifdef IO_ENABLED
                        if( _info->dbug.debug )
                                {
                                sprintf(LSGRG_MSGBUFFER, "\n ...COLUMN = %6ld  VARIABLE = %6ld  LENGTH = %6ld\n",
                                 j, jcol, collen );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                        if( collen <= 0 || collen > _info->cmax.colmax )
                                goto L_130;
                        for( ictr = istart; ictr <= iend; ictr++ ){
                                ictr_ = ictr - 1;
                                if( ihag[ictr_] == _info->nintbk.nobj )
                                        goto L_110;     /* Skip Obj Row */
                                nz = nz + 1;
                                if( nz > _info->memory.lbinv )
                                        goto L_140;
                                binv[nz - 1] = grad[ictr_];
                                IND(0,nz - 1) = (short) ihag[ictr_];
                                IND(1,nz - 1) = (short) j;
#ifdef IO_ENABLED
                                if( _info->dbug.debug )
                                        {
                                        sprintf(LSGRG_MSGBUFFER, " NZ = %5ld  ROW = %5d     COEFF = %13.6g\n",
                                         nz, IND(0,nz - 1), binv[nz - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                        }
#endif
L_110:
                                ;
                        }
                }

        }

#ifdef IO_ENABLED
        if( _info->dbug.debug )
                {
                sprintf(LSGRG_MSGBUFFER,
                 "\n***EXITING SUBROUTINE: GETBAS  -  NZBAS = %6ld\n",
                 nz );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        *nzbas = nz;


        return;

L_130:
#ifdef IO_ENABLED
        sprintf( LSGRG_MSGBUFFER,
           "\n GETBAS...ERROR: THE BASIC COLUMN IN POSITION%6ld HAS LENGTH%6ld"
           "      (ACTUAL VARIABLE NUMBER = %6ld)\n",
           j, collen, ibv[j - 1] );
         lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
        goto L_150;

L_140:
#ifdef IO_ENABLED
        lsgrg_error_msg(IOINFO,
          "\n GETBAS...INSUFFICIENT MEMORY FOR BASIS MATRIX"
          ", INITIATING REALLOCATION  \n" );
#endif

L_150:
        ;
/*      error = 0; */
/*      xerror( error, iounit.ioerr ); */
                        lsgrg_errorexit(JUMPBUF,_LSGRG_GETBAS_INSFMEM);
        return;


        /* :::::::::::::::::::::::::::::::::::  LAST CARD OF GETBAS  :::::::::: */
#undef  IND
} /* end of function */



void /*FUNCTION*/ LSGRGCLASS xpivot(LsgrgInfo *_info,double grad[], long int ihag[], long int iheg[],
         long int ibmap[], double zmem[], long int ibv[], long int ient,
         long int lv, double gcol[], long int *rcode, double cdnum[],
         long int irank[], LOGICAL32 getcd)
{
long int collen, /* error,  */ i, i_, iend, iflag, ilv, ilvlen, irow, istart,
         /* ka, kind, kip, klnrc, kw, kw1, kw2, */  nzinb, nzinw;

/*----------------------------------------------------------------*/
/*  working storage for invert subsystem                          */
/*  define local pointers to alias global ones                    */
/*  lnrc and ind are defined as short in invert subsystem         */
/*  so cast global long pointers to short here                    */
/*----------------------------------------------------------------*/
      DoubleArray   w1   = _info->inv_w1;
      DoubleArray   w2   = _info->inv_w2;
      DoubleArray   w    = _info->inv_w;
      DoubleArray   a    = _info->inv_binv;
      IntegerArray  ip   = _info->inv_ip;
      Integer2Array lnrc = (short*) _info->inv_lnrc;
      Integer2Array ind  = (short*) _info->inv_ind;
/*----------------------------------------------------------------*/
        /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *......................................................................
         *     *****FUNCTIONS
         *     NONE
         *     *****SUBROUTINES CALLED
         *     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****PURPOSE:
         *     THE PURPOSE OF THIS ROUTINE IS TO REPLACE BASIC VARIABLE LV
         *     BY PIVOTING VARIABLE IENT INTO THE BASIS. FIRST, THE ENTERING
         *     COLUMN IS FTRAN'ED, THEN LSLUPV IS CALLED TO UPDATE THE BASIS
         *     INVERSE REPRESENTATION.
         *
         *
         *     *****ARGUMENT DESCRIPTION:
         *     ON INPUT:
         *
         *     GRAD     - PACKED JACOBIAN
         *     IHAG     - ROW NUMBERS FOR GRADF
         *     IHEG     - COLUMN POINTER ARRAY
         *     BINV     - PACKED BASIS INVERSE REPRESENTATION
         *     IBMAP    - BASIS INVERSE WORK ARRAY MAP
         *     ZMEM     - PARTIONED WORKING STORAGE
         *     IBV      - LIST OF BASIC VARIABLE INDECIES
         *     IENT     - ENTERING VARIABLE NUMBER
         *     LV       - LEAVING BASIS VARIABLE
         *     NB       - NUMBER OF BASIC VARIABLES.  (in COMMON)
         *     ON OUTPUT:
         *     NZBAS    - NUMBER OF NON-ZEROES IN BASIS MATRIX (in COMMON)
         *     RCODE    - IS A RETURN CODE FOR XPIVOT
         *                 1 MEANS EVERYTHING OK;
         *                 2 MEANS THAT THE BASIS IS SINGULAR;
         *                 3 MEANS STORAGE OVERFLOW.
         *                 4 MEANS OTHER ERROR CONDITIONS
         *                 5 MEANS BASIS IS ILL-CONDITIONED
         *
         *
         *     THE CURRENT BASIS HAS BEEN RE-FACTORED
         *
         *     *****APPLICATION AND USAGE RESTRICTIONS:
         *     THIS SUBROUTINE IS COMPATIBLE WITH THE LSINVERT ROUTINES
         *     FOR HANDLING AN LU FACTORIZATION OF THE BASIS MATRIX.  THE
         *     'BODY OF PROGRAM' MUST BE RE-WRITTEN IF SOME OTHER INVERSE
         *     REPRESENTATION IS USED.
         *
         *     *****ALGORITHM NOTES:
         *     DATA STRUCTURE FOR THE PROBLEM DATA:
         *
         *     IBMAP(3)  -   POINTS TO THE LNRC INTEGER*2 ARRAY.
         *     IBMAP(4)  -   POINTS TO THE IP INTEGER ARRAY.
         *     IBMAP(5)  -   POINTS TO THE W REAL*8 ARRAY.
         *     IBMAP(6)  -   POINTS TO THE IND INTEGER*2 ARRAY.
         *     IBMAP(7)  -   POINTS TO THE A REAL*8 ARRAY.
         *
         *     THE WORK ARRAY GCOL IS USED TO UNPACK THE ENTERING COLUMN
         *
         *
         *     *****REFERENCES:
         *     'FORTRAN SUBROUTINES FOR HANDLING SPARSE LINEAR PROGRAMMING BASES'
         *     BY J. K. REID, REPORT AERE-R8269, COMPUTER SCIENCE AND SYSTEMS
         *     DIVISION, AERE HARWELL, OXFORDSHIRE, ENGLAND. JANUARY, 1976.
         *     'A SPARSITY-EXPLOITING VARIANT OF THE BARTELS-GOLUB
         *     DECOMPOSITION FOR LINEAR PROGRAMMING BASES', J. K. REID,
         *     MATHEMATICAL PROGRAMMING, VOL. 24, NO. 1, SEPTEMBER 1982.
         *     *****HISTORY:
         *     WRITTEN BY:  STUART H. SMITH  DEPARTMENT OF MECHANICAL ENGINEERING
         *                  UNIVERSITY OF TEXAS, AUSTIN, TEXAS 78712.
         *     DATE LAST MODIFIED:  JAN 31 1990
         *     *****GENERAL
         *::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         *
         *     *****BODY OF PROGRAM  (XPIVOT)
         *
         *C     --------------------------
         *C    | CHECK ENTERING COLUMN
         *C     --------------------------
         * */
#ifdef IO_ENABLED
        if( _info->dbug.debug )
                {
                sprintf(LSGRG_MSGBUFFER,
                 "\n ****ENTERING SUBROUTINE: XPIVOT   -- GETCD = %2c  ENTERING = %6ld LEAVING = %6ld\n",
                 TorF(getcd), ient, lv );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

        /*C     -------------------
         *C    | GET WORK ARRAYS
         *C     -------------------
         *
         * */
/*
        kw1 = ibmap[0];
        kw2 = ibmap[1];
        klnrc = ibmap[2];
        kip = ibmap[3];
        kw = ibmap[4];
        kind = ibmap[5];
        ka = ibmap[6];
*/

        /*      -----------------------------------------------------------
         *     : COPY COLUMN INTO THE W ARRAY
         *     | LET LSLUFT COUNT THE NUMBER OF NONZEROS
         *     | SKip the objective row, since it isn't stored in basis
         *      -----------------------------------------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
         /*     zmem[kw + i_ - 1] = 0.0e0;  */
                     w[i_]        = 0.0e0;
        }
        nzinw = 0;
        if( ient <= _info->dimen.n ){
                istart = iheg[ient - 1];
                iend = iheg[ient] - 1;
                collen = iend - istart + 1;
                if( collen <= 0 || collen > _info->cmax.colmax )
                        goto L_300;
                for( i = istart; i <= iend; i++ ){
                        i_ = i - 1;
                        irow = ihag[i_];
                        if( irow != _info->nintbk.nobj )
                          /*    zmem[kw + irow - 2] = grad[i_]; */
                                   w[irow - 1]      = grad[i_];
                }
        } else{
         /*     zmem[kw + ient - _info->dimen.n - 2] = -1.0e0; */
                   w[ient - _info->dimen.n - 1]      = -1.0e0;
                collen = 1;
        }
 /*     prtdarray("w","w before call to lsluft",w,_info->nintbk.nb); */

        /*C     ---------------------------------------------------------
         *C    |    REPLACE LEAVING COLUMN WITH ENTERING COLUMN IN IBV
         *C     ---------------------------------------------------------
         * */

        ilv = ibv[lv - 1];
        if( ilv > _info->dimen.n ){
                ilvlen = 1;
        } else{
                ilvlen = iheg[ilv] - iheg[ilv - 1];
        }
        ibv[lv - 1] = ient;
        _info->nzerob.nzbas = _info->nzerob.nzbas + collen - ilvlen;

        /*      ---------------------------------
         *     : CALL THE LSINVERT FTRAN ROUTINE
         *      ---------------------------------
         * */
/*
        lsluft( _info->nintbk.nb, _info->memory.lbinv, &zmem[ka - 1], (short*)&zmem[kind - 1],
         (long*)&zmem[kip - 1], (short*)&zmem[klnrc - 1], &zmem[kw - 1],
         gcol, &nzinw, &nzinb, _info->lsinvt.iluprm );
*/

        lsluft( _info->nintbk.nb, _info->memory.lbinv, a            , ind                    ,
         ip                   , lnrc                    , w            ,
         gcol, &nzinw, &nzinb, _info->lsinvt.iluprm );

 /*     prtdarray("w","w after call to lsluft",w,_info->nintbk.nb); */
        /*C------------------------------------------------------------
         *C    CALL LSLUPV TO UPDATE THE INVERSE
         *C------------------------------------------------------------
         * */

        *rcode = 1;

/*
        lslupv( _info->nintbk.nb, _info->memory.lbinv, (short*)&zmem[kind - 1], &zmem[ka - 1],
         (long*)&zmem[kip - 1], (short*)&zmem[klnrc - 1], &zmem[kw - 1],
         lv, _info->lsinvt.iluprm, _info->lsinvt.dluprm, &iflag );
*/

        lslupv( _info->nintbk.nb, _info->memory.lbinv, ind                    , a            ,
         ip                   , lnrc                    , w            ,
         lv, _info->lsinvt.iluprm, _info->lsinvt.dluprm, &iflag );


        /*     IF(DEBUG) WRITE(IOOUT,1220) IFLAG
         * */
        if( getcd && iflag == 0 ){
/*
                condnm(_info, _info->nintbk.nb, (short*)&zmem[kind - 1], &zmem[ka - 1],
                 _info->memory.lbinv, (short*)&zmem[klnrc - 1], (long*)&zmem[kip - 1],
                 _info->lsinvt.iluprm, &zmem[kw - 1], irank, cdnum, grad, ihag, iheg,
                 ibv, &zmem[kw1 - 1], &zmem[kw2 - 1] );
*/

                condnm(_info, _info->nintbk.nb, ind                    , a            ,
                 _info->memory.lbinv, lnrc                    , ip                   ,
                 _info->lsinvt.iluprm, w            , irank, cdnum, grad, ihag, iheg,
                 ibv, w1            , w2             );

                if( _info->zblck.condmx > _info->zcond.cndtol )
                        iflag = -20;
        }


        /*C---------------------------------------------------
         *C    ERROR RETURNS
         *C---------------------------------------------------
         * */
        i = labs( iflag );

        /*#############################
         *   CHANGED 7/12/90 --- INVERT ROUTINES NOW RETURN 0 IF OK (SHS)
         *#############################
         *
         *     IF (IFLAG .EQ. 1) THEN */
        if( iflag == 0 ){
                *rcode = 1;
        } else if( i == 7 ){
                lsgrg_error_msg(IOINFO,
                  "\n XPIVOT... CURRENT BASIS MATRIX IS SINGULAR \n" );
                *rcode = 2;
        } else if( iflag == -20 ){
                *rcode = 5;
        } else if( iflag == 9 ){
                lsgrg_error_msg(IOINFO,
                  "\n XPIVOT... INSUFFICIENT SPACE FOR BASIS INVERSE \n" );
                *rcode = 3;
        } else{
                *rcode = 4;
        }

        /*     IF(RCODE .EQ. 1)RETURN
         *
         *     ** PUT RETURN CODES HERE ***
         *     ** OVERFLOW AND SINGULAR *** */
        return;

        /*---------------------------------------------
         *     COLUMN LENGTH ERROR
         *----------------------------------------------
         * */
L_300:
        ;
#ifdef IO_ENABLED
        sprintf( LSGRG_MSGBUFFER,
          "\n XPIVOT...ERROR: THE ENTERING COLUMN IN POSITION%6ld"
          " HAS LENGTH%6ld\n",
         ient, collen );
        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
/*      error = 0; */
/*      xerror( error, iounit.ioerr ); */
                        lsgrg_errorexit(JUMPBUF,_LSGRG_XPIVOT_COLLEN);

        return;


        /* :::::::::::::::::::::::::::::::::::  LAST CARD OF XPIVOT  :::::::::: */
} /* end of function */
void LSGRGCLASS prtdarray(LsgrgInfo *_info,char *name,char *msg, double x[], long n)
{
#ifdef IO_ENABLED
    long i;
    sprintf(LSGRG_MSGBUFFER,
     "\n ... print of array : %s",msg);
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
    for(i=0;i< n;++i)
        sprintf(LSGRG_MSGBUFFER,"\n %s[%5d] = %14.7e",name,i,x[i]);
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
    lsgrg_msg(IOINFO,"\n.....................");
#endif
    return;
}
#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"

#define IOINFO  &(_info->io_info)
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)
/***********************************************************************
 **                 FILE:        GRGCONSB  FORTRAN                     *
 **                 AUTHOR:      STUART SMITH                          *
 **                 CREATED:      5 MARCH 1988                         *
 **                 LAST UPDATE: 06 MAY   1998                         *
 **                                                                    *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **      This file contains the following routines for lsgrg2, a       *
 **  large scale implementation of GRG2:                               *
 **                                                                    *
 **  ROUTINE           DESCRIPTION                        LAST UPDATE  *
 ** ---------         -------------                      ------------- *
 **  CHUZQ      Performs a degenerate pivot               04 JAN 1990  *
 **  CHUZR      Checks for degenercy                      22 MAY 1996  *
 **  CONSBS     Constructs the basis matrix               06 MAY 1998  *
 **  GETGR      Computes the Jacobian in sparse form      10 JUN 1993  *
 **  GRDUMP     Prints the sparse Jacobian structure      29 JAN 1991  *
 **  PARSHC     Central difference partial derivatives    16 MAR 1988  *
 **  PARSHF     Forward difference partial derivatives    16 MAR 1988  *
 **  PRSCAL     Scales the problem data                   29 JUN 1993  *
 **  SCLJAC     Scales the Jacobian matrix                29 JUN 1993  *
 **  M2SCAL     Calculates scale factors                  10 JUN 1993  *
 **  CONDNM     Estimates cond. number of basis           29 JAN 1991  *
 **  CHKELT     Computes the largest Jacobian element     29 JUB 1993  *
 **             and checks for reasonable size                         *
 **                                                                    *
 ** REVISIONS:                                                         *
 ** 05/22/96 SHS - Add argument ICOLS to CHUZR.  CHUZR will set ICOLS*
 **    to determine the bound status of basic variables for use by     *
 **    perturbation logic.                                             *
 **    See notes in each subroutine
 **
 ** 06/28/96 JCP -- Add code to GETGR to handle new nonzeros
 **    by storing revised column and shifting rest of jacobian to right
 **    (modifications to VARMEM to preallocate expansion space)
 **    added new routine MODJAC
 **
 ** 05/06/98 SHS -- Add code to handle new Analytical PARSH.  The new
 **    PARSH returns 3 arrays:
 **        PAIJ  - Array of nonzero elements
 **        IPROW - Row index of each element (< 0 implies a linear element)
 **        IPCOL - Col index of each element
 **
 **    The array IPMAP is used to map each element in PAIJ into the
 **    correct position into the GRAD array.
 **
 **
 ***********************************************************************
 *
 * */
void /*FUNCTION*/ LSGRGCLASS chuzq(LsgrgInfo *_info,double grad[], long int ihag[], long int iheg[],
         long int ibmap[], double zmem[], double v[], double d[], long int inbv[],
         double alb[], double ub[], double x[], long int ibv[], long int iub[],
         long int icand[], long int ivstat[], long int *tablst, double cdnum[],
         long int irank[])
{
#define TABLST(I_,J_)   (*(tablst+(I_)*(_info->chq.maxtab)+(J_)))
long int i, i_, icol, ii, ipchq, j, j_, jq, jq1, jq2, jq3, k, k_,
         rcode;
double d1, d2, dmax, pivmax, pivot, pvt, sum, tol, xj;

        /* ....................................................................
         *
         *     Choose maximum pivot from columnl of BINV
         *     D(1),...,D(NSUPER) will hold the pivot choices
         * */
        ipchq = _info->nintbk.ipr;

#ifdef IO_ENABLED
        if( ipchq >= 5 ){
               lsgrg_msg(IOINFO, "\n CHUZQ ENTERED \n" );
                sprintf(LSGRG_MSGBUFFER, " LEAVING VARIABLE IS BASIC VARIABLE NO.  %5ld VAR# %5ld\n",
                 _info->misc.lv, ibv[_info->misc.lv - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                sprintf(LSGRG_MSGBUFFER, " NSUPER = %6ld\n", _info->nintbk.nsuper );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
               lsgrg_msg(IOINFO, " INBV IS" );
                for( i = 1; i <= _info->dimen.n; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", inbv[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n IUB IS " );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                for( i = 1; i <= _info->dimen.n; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", iub[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
               lsgrg_msg(IOINFO, "\n IBV IS " );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibv[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
               lsgrg_msg(IOINFO, "\n ICAND IS" );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                for( i = 1; i <= (_info->nintbk.nb + _info->dimen.n); i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", icand[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                sprintf(LSGRG_MSGBUFFER, "\n CHUZQ... LTAB = %5ld  LPTR = %5ld",
                 _info->chq.ltab, _info->chq.lptr );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
               lsgrg_msg(IOINFO, "\n TABLST... \n VARS:\t" );
                for( i = 1; i <= _info->chq.ltab; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", TABLST(0,i - 1) );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                if( _info->chq.ltab > 0 )
                        {
                       lsgrg_msg(IOINFO, " ROWS:    " );
                        for( i = 1; i <= _info->chq.ltab; i++ ){
                                sprintf(LSGRG_MSGBUFFER, "%4ld", TABLST(1,i - 1) );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                        lsgrg_msg(IOINFO,"\n" );
                        }
        }
#endif

        /*C     --------------------------------
         *C    | Compute Pivot Row in BINV
         *C     --------------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                v[i_] = 0.0e0;
        }
        v[_info->misc.lv - 1] = 1.0e0;
        xbtran(_info, ibmap, zmem, v, _info->nintbk.nb );

#ifdef IO_ENABLED
        if( ipchq >= 6 )
                {
               lsgrg_msg(IOINFO, " ROW OF BINV....\n " );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", v[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                }

#endif
        /*C     ---------------------------------------------
         *C    | Compute updated Pivot Row for SuperBasics
         *C    | JQ = Index of Superbasic with largest pivot
         *C    |      which is NOT on the Tabu List.
         *C     ---------------------------------------------
         * */
        jq = 0;
        pivmax = 0.0e0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                ii = inbv[i_];
                d[i_] = 0.0e0;
                if( ii > _info->dimen.n )
                        goto L_43;
                if( ivstat[ii - 1] < 0 )
                        goto L_50;
L_43:
                if( _info->chq.ltab == 0 )
                        goto L_47;
                for( j = 1; j <= _info->chq.ltab; j++ ){
                        j_ = j - 1;
                        if( (TABLST(0,j_) == ii) && (TABLST(1,j_) == _info->misc.lv) )
                                goto L_50;
                }
L_47:
                xdot(_info, grad, ihag, iheg, v, _info->nintbk.nb, ii, &sum );
                d[i_] = sum;
                if( pivmax >= fabs( sum ) )
                        goto L_50;
                pivmax = fabs( sum );
                jq = i;
L_50:
                ;
        }
#ifdef IO_ENABLED
        if( ipchq >= 5 ){
               lsgrg_msg(IOINFO,
                " UPDATED SUPERBASIC PIVOT ROW IS .... \n " );
                for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", d[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                sprintf(LSGRG_MSGBUFFER, " CHUZQ JQ = %5ld LARGEST PIVOT = %13.6e\n",
                 jq, pivmax );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        /*      -------------------------------------------------------
         *     | If all pivot elements too small,  go look at
         *     | Nonbasic Variables.
         *      -------------------------------------------------------
         * */
        if( jq == 0 || pivmax < _info->limits.epspiv )
                goto L_150;

        /*C     ------------------------------------------------
         *C    | Choose one away from its bounds if possible.
         *C    | JQ2 = Index of Superbasic farthest from bound and
         *C    |       with "large enough" pivot
         *C     ------------------------------------------------
         * */
        tol = fmax( _info->pvtblk.xpvpct*pivmax, _info->limits.epspiv );

L_152:
        dmax = 0.0e0;
        pvt = 0.0e0;
#ifdef IO_ENABLED
        if( ipchq > 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " CHECKING SUPERBASIC PIVOTS -- TOL = %14.7e\n",
                 tol );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

        jq2 = 0;
        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;

#ifdef IO_ENABLED
                if( ipchq > 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " SUPERBASIC NO.  %4ld    ROW ELEM = %15.8e\n",
                         j, d[j_] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( fabs( d[j_] ) < tol )
                        goto L_60;
                k = inbv[j_];
                xj = x[k - 1];
                d1 = fabs( xj - alb[k - 1] )/(1.0e0 + fabs( alb[k - 1] ));
                d2 = fabs( ub[k - 1] - xj )/(1.0e0 + fabs( ub[k - 1] ));
                if( d1 > d2 )
                        d1 = d2;

#ifdef IO_ENABLED
                if( ipchq > 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " VARIABLE NO.    %4ld    DBOUND =   %15.8e\n",
                         k, d1 );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( dmax < d1 ){
                        pvt = fabs( d[j_] );
                        dmax = d1;
                        jq2 = j;
                }
L_60:
                ;
        }

        /*     -----------------------------------
         *    | Found a pivot away its bound
         *     -----------------------------------
         * */
        if( jq2 > 0 && pvt >= _info->limits.epspiv ){
                jq = jq2;
                goto L_70;
        }

        /*     -----------------------------------------------------------
         *    | Look for one away from bound and at least EPSPIV
         *     -----------------------------------------------------------
         * */
        if( tol > _info->limits.epspiv ){
                tol = _info->limits.epspiv;
                goto L_152;
        }

        /*     ------------------------------------------------------------
         *    | If pivot is large enough, just use the Super Basic with
         *    | the largest pivot (Held in JQ).
         *     ------------------------------------------------------------
         * */
        if( jq > 0 )
                goto L_70;


        /*    -----------------------------------------------------------------
         *   | Choose a Nonbasic Variable  - No Super Basic could be chosen.
         *   | Find Nonbasic Variable with largest pivot element.
         *   | JQ1 = Index of Nonbasic with largest pivot element
         *    -----------------------------------------------------------------
         * */
L_150:
        ;
        jq1 = 0;
        jq = 0;
        dmax = 0.0e0;

#ifdef IO_ENABLED
        if( ipchq >= 5 )
                {
               lsgrg_msg(IOINFO,
                " CHECKING NON-BASIC PIVOTS.....     \n" );
                }
#endif

        for( i = _info->nintbk.nsuper + 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                ii = inbv[i_];
                if( _info->chq.ltab == 0 )
                        goto L_147;
                for( j = 1; j <= _info->chq.ltab; j++ ){
                        j_ = j - 1;
                        if( (TABLST(0,j_) == ii) && (TABLST(1,j_) == _info->misc.lv) )
                                goto L_151;
                }
L_147:
                xdot(_info,grad, ihag, iheg, v, _info->nintbk.nb, ii, &sum );

#ifdef IO_ENABLED
                if( ipchq >= 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " NONBASIC NO.  %4ld  VARIABLE NO.  %4ld    ROW ELEM = %15.8e\n",
                         i, ii, sum );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                pivot = fabs( sum );
                if( dmax >= pivot )
                        goto L_151;
                dmax = pivot;
                jq1 = i;
L_151:
                ;
        }
        if( dmax > _info->limits.epspiv )
                jq = jq1;

        /*      -----------------------------------------------------------
         *     | Pivot has been chosen, JQ should be > 0
         *     | If JQ = 0, use a variable in the TABU LIST if possible
         *      -----------------------------------------------------------
         * */
        if( jq > 0 )
                goto L_70;

#ifdef IO_ENABLED
        if( ipchq >= 5 )
                {
               lsgrg_msg(IOINFO,
                 " CHECKING TABU LIST PIVOTS.....     \n" );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( _info->chq.ltab == 0 )
                goto L_67;
        for( i = 1; i <= _info->chq.ltab; i++ ){
                i_ = i - 1;
                if( TABLST(1,i_) != _info->misc.lv )
                        goto L_66;
                ii = TABLST(0,i_);
                xdot(_info, grad, ihag, iheg, v, _info->nintbk.nb, ii, &sum );
#ifdef IO_ENABLED
                if( ipchq >= 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER, "  VARIABLE NO.  %4ld    ROW ELEM = %15.8e\n",
                         ii, sum );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( fabs( sum ) >= _info->limits.epspiv ){
                        for( k = 1; k <= _info->dimen.n; k++ ){
                                k_ = k - 1;
                                if( inbv[k_] != ii )
                                        goto L_68;
                                jq = k;
                                goto L_70;
L_68:
                                ;
                        }
                }
L_66:
                ;
        }

        /*     -----------------------------------
         *    | All pivots < EPSPIV -- So ABORT
         *     -----------------------------------
         * */
L_67:
#ifdef IO_ENABLED
        sprintf( LSGRG_MSGBUFFER, " ***ERROR...CHUZQ -- COULD NOT FIND A PIVOT > EPSPIV IN ROW %4ld  WHEN ATTEMPTING TO PIVOT OUT VARIABLE %5ld\n",
         _info->misc.lv, ibv[_info->misc.lv - 1] );
        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
    /*  xerror( 0, iounit.ioerr );                */
        lsgrg_errorexit(JUMPBUF,_LSGRG_CHUZQ_BADPIVOT);
        return;


        /*C     ----------------------
         *C    | Perform the pivot
         *C     ----------------------
         * */
L_70:
        icol = inbv[jq - 1];

        jq2 = ibv[_info->misc.lv - 1];

#ifdef IO_ENABLED
        if( ipchq >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " VARIABLE%4ld LEAVING BASIS  - BASIC VAR NO.  %4ld\n VARIABLE%4ld ENTERING BASIS - NON-BASIC NO. %4ld\n",
                 jq2, _info->misc.lv, icol, jq );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        jq3 = icand[_info->misc.lv - 1];

        xpivot(_info, grad, ihag, iheg, ibmap, zmem, ibv, icol, _info->misc.lv, v,
         &rcode, cdnum, irank, TRUE );

        if( rcode == 1 )
                goto L_163;

#ifdef IO_ENABLED
        if(ipchq > 1) {
           sprintf( LSGRG_MSGBUFFER,
            " CHUZQ....WARNING: RETURN FROM XPIVOT - RCODE = %5ld\n",
            rcode );
           lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
           lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);

           if( rcode == 5 ){
               sprintf( LSGRG_MSGBUFFER,
                " CHUZQ....WARNING: UPDATED BASIS MATRIX IS ILL-CONDITIONED\n" );
               lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
               lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
           }
        }
#endif
        if( rcode == 2 ) {
#ifdef IO_ENABLED
            sprintf( LSGRG_MSGBUFFER,
             " CHUZQ....ERROR: UPDATED BASIS MATRIX IS SINGULAR \n" );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
            lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
    /*          xerror( 0, iounit.ioerr );         */
            lsgrg_errorexit(JUMPBUF,_LSGRG_XPIVOT_BASIS_SING);
        }
                return;

        /*      ----------------------------------------------------------
         *     | Update index sets of Basic and Nonbasic variables and IUB
         *      ----------------------------------------------------------
         *C */
L_163:
        ;
        icand[_info->misc.lv - 1] = icol;
        if( jq < _info->nintbk.nsuper )
                goto L_167;
        if( jq == _info->nintbk.nsuper )
                goto L_168;

        /*      -------------------------------
         *     | Entering column is Nonbasic
         *      ------------------------------- */
        _info->nintbk.nsuper = _info->nintbk.nsuper + 1;
        if( jq == _info->nintbk.nsuper )
                goto L_168;
        k = jq - 1;
        for( j = _info->nintbk.nsuper; j <= k; j++ ){
      /*        j_ = j - 1; */
                i = _info->nintbk.nsuper + k - j;
                ii = i + _info->nintbk.nb;
                icand[ii] = icand[ii - 1];
                iub[i] = iub[i - 1];
                inbv[i] = inbv[i - 1];
        }
        goto L_168;

        /*      --------------------------------
         *     | Entering column is Superbasic
         *      -------------------------------- */
L_167:
        k = _info->nintbk.nsuper - 1;
        for( i = jq; i <= k; i++ ){
                i_ = i - 1;
                ii = i + _info->nintbk.nb;
                icand[ii - 1] = icand[ii];
                iub[i_] = iub[i_ + 1];
                inbv[i_] = inbv[i_ + 1];
        }

        /*      ---------------------------------------------
         *     | Now store leaving variable in proper place
         *      --------------------------------------------- */
L_168:
        ;
        iub[_info->nintbk.nsuper - 1] = 0;
        icand[_info->nintbk.nb + _info->nintbk.nsuper - 1] = jq3;
        tol = _info->limits.epboun;
        if( fabs( x[jq2 - 1] - ub[jq2 - 1] ) <= (tol*(1.0e0 + fabs( ub[jq2 - 1] ))) ){
                /*XXXXX      X(JQ2)=UB(JQ2) */
                iub[_info->nintbk.nsuper - 1] = 1;
        } else if( fabs( x[jq2 - 1] - alb[jq2 - 1] ) <= (tol*(1.0e0 +
         fabs( alb[jq2 - 1] ))) ){
                /*XXXX       X(JQ2)=ALB(JQ2) */
                iub[_info->nintbk.nsuper - 1] = -1;
        }
        inbv[_info->nintbk.nsuper - 1] = jq2;

        /*      --------------------------
         *     |  Update the Tabu list
         *      --------------------------
         * */
        if( _info->chq.maxtab <= 0 )
                goto L_4323;
        if( _info->chq.ltab < _info->chq.maxtab )
                _info->chq.ltab = _info->chq.ltab + 1;
        _info->chq.lptr = _info->chq.lptr + 1;
        if( _info->chq.lptr > _info->chq.maxtab )
                _info->chq.lptr = 1;
        TABLST(0,_info->chq.lptr - 1) = jq2;
        TABLST(1,_info->chq.lptr - 1) = _info->misc.lv;
L_4323:
        if( iub[_info->nintbk.nsuper - 1] != 0 )
                _info->nintbk.nsuper = _info->nintbk.nsuper - 1;

#ifdef IO_ENABLED
        if( ipchq >= 4 ){
                sprintf(LSGRG_MSGBUFFER, " CHUZQ... LTAB = %5ld  LPTR = %5ld",
                 _info->chq.ltab, _info->chq.lptr );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
               lsgrg_msg(IOINFO, "\n TABLST... \n VARS:\t" );
                for( i = 1; i <= _info->chq.ltab; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", TABLST(0,i - 1) );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
               lsgrg_msg(IOINFO, "\n ROWS:    " );
                for( i = 1; i <= _info->chq.ltab; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", TABLST(1,i - 1) );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                sprintf(LSGRG_MSGBUFFER, "\n NSUPER = %6ld\n", _info->nintbk.nsuper );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
               lsgrg_msg(IOINFO, " INBV IS" );
                for( i = 1; i <= _info->dimen.n; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", inbv[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
               lsgrg_msg(IOINFO, "\n IUB IS " );
                for( i = 1; i <= _info->dimen.n; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", iub[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
               lsgrg_msg(IOINFO, "\n IBV IS " );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibv[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
               lsgrg_msg(IOINFO, "\n ICAND IS" );
                for( i = 1; i <= (_info->nintbk.nb + _info->dimen.n); i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", icand[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
        }
        if( ipchq >= 5 )
                {
               lsgrg_msg(IOINFO, " CHUZQ COMPLETED \n" );
                }
#endif
        return;



        /*     End of CHUZQ
         * */
#undef  TABLST
} /* end of function */



void /*FUNCTION*/ LSGRGCLASS chuzr(LsgrgInfo *_info,double alb[], double ub[], double x[], double dx[],
         double g[], long int ibv[], long int icols[], long int ibc[],
         long int *jp, double *theta, LOGICAL32 *move)
{
        /* 05/01/02 jcp ** comment out unrefernced vars */
long int /*i,*/ ibnd, j, j_, jr, k;
double bndjr, d, d1, /* d2,*/  pertbn, pivot, psi, t, tmax, tr, xjr, ztol;

        /*SHS  SUBROUTINE CHUZR (ALB,UB,X,DX,G,IBV,IBC,JP,THETA,MOVE)
         *
         *
         * SHS 05/01/96 - Changed CHUZR to not consider basic variables away from
         *                bound as degenerate.  For a variable to be selected as
         *                degenerate, the  variable must be at a bound and the step
         *                along the tangent vector must be small (i.e.  < TOLX).
         * LSL 06/14/96 - Changed tolerance for skipping variables with tiny
         *                tangent vectors to sqr root of eps. Without this,h20 and h16
         *                do not solve due to small pivot and degen probs. Prev tol
         *                ,eps, was too small.
         * SHS 05/22/96 - Add ICOLS to argument list.  ICOLS(J) will be set in CHUZR
         *                to determine the bound status of basic variable "J".
         *                Possible values are:
         *
         *            0  -> Variable can be moved in tanget direction
         *           -1  -> Variable at lower bound and tanget points into bound
         *            1  -> Variable at upper bound and tangent points into bound.
         *
         *-------------------------------------------------------------------------------
         *
         * */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
               lsgrg_msg(IOINFO, " CHUZR ENTERED \n" );
                }
#endif
        *move = TRUE;
        *theta = _info->tols.plinfy;
        psi = _info->tols.plinfy;
        /*CXX  PERTBN = EPNEWT */
        pertbn = _info->limits.epboun;
        ztol = pow(_info->tols.eps, 0.5e0);
        *jp = 0;
        /* SHS 05/01/96 - Compute D1 = Distance from bound for basic variables
         * SHS 05/22/96 - Set value for ICOLS to flag degenerate basic variables.
         * */
        for( j = 1; j <= _info->nintbk.nb; j++ ){
                j_ = j - 1;
                icols[j_] = 0;
                if( j == _info->nintbk.nobj )
                        goto L_100;
                t = -dx[j_];
                if( fabs( t ) <= ztol )
                        goto L_100;
                k = ibv[j_];
                if( _info->nintbk.ninf > 0 && k > _info->dimen.n )
                        goto L_100;
                if( t >= 0.0e0 ){
                        d = x[k - 1] - alb[k - 1] + pertbn*(1.0e0 + fabs( alb[k - 1] ));
                        d1 = (x[k - 1] - alb[k - 1])/(1.0e0 + fabs( alb[k - 1] ));
                        ibnd = -1;
                } else{
                        d = x[k - 1] - ub[k - 1] - pertbn*(1.0e0 + fabs( ub[k - 1] ));
                        d1 = (ub[k - 1] - x[k - 1])/(1.0e0 + fabs( ub[k - 1] ));
                        ibnd = 1;
                }
                if( fabs( d ) > 1.0e20 )
                        d = sign( 1.0e20, d );

                /*  SHS 5/01/96 - Only look at variables at bounds
                 *  SHS 5/22/96 - Set ICOLS - Should we consider size of tangent when
                 *                setting ICOLS?
                 * */
                if( d1 <= _info->limits.epboun ){
                        icols[j_] = ibnd;
                        t = d/t;
                        if( t < psi ){
                                psi = t;
                                *jp = j;
                        }
                }
L_100:
                ;
        }
        if( *jp == 0 )
                goto L_160;

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " PSI =%13.6e", psi );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
               lsgrg_msg(IOINFO, " DX IS\n " );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", dx[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                }
#endif
        jr = ibv[*jp - 1];

        /*      -----------------------------
         *     | SECOND PASS OF HARRIS
         *      -----------------------------
         *
         *!!!  IF (PSI .LE. EPS) THEN
         *!!!     JR = IBV(JP)
         *!!!     THETA = 0.0
         *!!!     GO TO 2000
         *!!!  ENDIF */
        tmax = 0.0e0;
        /* SHS 05/01/96 - Compute D1 = Distance from bound for basic variables */
        for( j = 1; j <= _info->nintbk.nb; j++ ){
                j_ = j - 1;
                if( j == _info->nintbk.nobj )
                        goto L_150;
                t = -dx[j_];
                if( fabs( t ) < ztol )
                        goto L_150;
                k = ibv[j_];
                if( _info->nintbk.ninf > 0 && k > _info->dimen.n )
                        goto L_150;
                if( t >= 0.0e0 ){
                        d = x[k - 1] - alb[k - 1];
                        d1 = (x[k - 1] - alb[k - 1])/(1.0e0 + fabs( alb[k - 1] ));
                } else{
                        d = x[k - 1] - ub[k - 1];
                        d1 = (ub[k - 1] - x[k - 1])/(1.0e0 + fabs( ub[k - 1] ));
                }
                if( fabs( d ) > 1.0e20 )
                        d = sign( 1.0e20, d );
                tr = d/t;
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " K =%5ld TR =%13.6e\n", k,
                         tr );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif

                /* SHS 5/01/96  Skip variables away from bounds
                 * */
                if( d1 <= _info->limits.epboun ){
                        if( (tr <= psi) && (fabs( t ) > tmax) ){
                                tmax = fabs( t );
                                *jp = j;
                                *theta = tr;
#ifdef IO_ENABLED
                                if( _info->nintbk.ipr >= 5 )
                                        {
                                        sprintf(LSGRG_MSGBUFFER, " JP =%5ld TMAX = %13.6e\n",
                                         *jp, tmax );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                        }
#endif
                        }
                }
L_150:
                ;
        }
        jr = ibv[*jp - 1];
        if( jr <= _info->dimen.n )
                goto L_2000;
        bndjr = alb[jr - 1];
        pivot = dx[*jp - 1];
        if( pivot > 0.0e0 )
                bndjr = ub[jr - 1];
        xjr = x[jr - 1];
        if( jr > _info->dimen.n )
                xjr = g[jr - _info->dimen.n - 1];
        *theta = (bndjr - xjr)/pivot;
L_2000:
        ;
        *move = *theta > _info->tols.tolx;
L_160:
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
               lsgrg_msg(IOINFO, " CHUZR COMPLETED \n" );
                }
#endif
        return;

        /*     End of CHUZR
         * */
} /* end of function */



/*    chg 5/6/98 SHS - Add parsh output arrays to argument list */

void /*FUNCTION*/ LSGRGCLASS consbs(LsgrgInfo *_info,long int ibmap[], double zmem[], double grad[],
         long int ihag[], long int iheg[], long int ihegl[], long int ibv[],
         double x[], double alb[], double ub[], long int inbv[], long int iub[],
         long int ibc[], double g[], long int inlin[], long int ivstat[],
         long int istat[], double dbnd[], double gg[], double pvrow[],
         long int icols[], long int iwork[], double gcol[], long int icand[],
         double cdnum[], long int irank[], double ascale[], long int ibvbct[],
         long int *nbfree, long int iobjpt[], double paij[], long int iprow[],
         long int ipcol[], long int ipmap[])
{
LOGICAL32 all, bschng, degprb, dropbv, first, incpvt, jreset, lastpv,
         newbas, reinvt;
long int i, i_, ibctr, icsupe, ictrb, ii, ipcb3, ipconb, /* irf, */ irow,
        /*  irs, */ islack, j, j_, jcol, jj, jpiv, k, nbcp, nbp, nchg, npnb,
         rcode,lsgrg_return_status;
double bnd, d1, d2, dd, dmax, eltmax, epnorm, sum, tol, umold, zlast,
         zstart, ztemp, ztime;

        /*.....................................................................
         * ----------------------------------------------------------------
         *      Subroutine CONSBS constructs a basis from the Jacobian of
         *      binding constraints using a modified complete pivoting procedure.
         *
         *
         *     INPUT VARIABLES
         *      ISTAT  ... Status array for constraints
         *
         *     OUTPUT VARIABLES
         *      NB     ... No. of basic variables
         *      IBV    ... Index set of Basic Variables
         *      INBV   ... Index set of Non-basic Variables
         *      IBC    ... Index set of binding constraints
         *      INBC   ... Index set of Nonbinding constraints
         *      BINV   ... Inverse of basis
         *      IUB    ... Bound indicator array for Nonbasics
         *
         *
         *
         *C    INITIALIZATIONS
         *C     NB   -  Number of Basic Variables (Always MP1 after first call)
         *C     NBP  -  Number of Basic Varibales in previous basis (MP1)
         *C     NBCP -  Number of binding constraints in previous call
         *C
         *C     NOTE --- NBP may be zero only on first call to CONSBS since
         *C              NB = MP1 after the first call.
         *C----------------------------------------------------------------- */
        /*XXXX2 46H CONDITION NUMBERS OF DIAGONAL BLOCKS ARE ... ,/(1X,10D13.6)) */

        /*      -------------------------
         *     | Save start time
         *      ------------------------- */
  /*    cputime( &zstart, &irs ); */
        zstart = lsgrg_timer();
        zlast = zstart;

        ipconb = _info->nintbk.ipr;
        ipcb3 = _info->nintbk.ipr3;

#ifdef IO_ENABLED
        if( ipconb >= 5 || _info->dbug.debug )
                {
               lsgrg_msg(IOINFO,
                " ***** ENTERING SUBROUTINE: CONSBS\n" );
                }
        if( ipconb >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " SMODE = %2c BASBND = %2c NEWPT = %2c JSTFES = %2c\n",
                 TorF(_info->cbmode.smode), TorF(_info->supblk.basbnd), TorF(_info->cbmode.newpt),
                 TorF(_info->srchlg.jstfes) );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        _info->cbscnt.ncbs = _info->cbscnt.ncbs + 1;
        first = _info->cbscnt.ncbs == 1;
        nbp = _info->nintbk.nb;
        nbcp = _info->bind.nbc;
        epnorm = _info->limits.epspiv;

        /*      --------------------------------------------------------------
         *     |               Set the logicals
         *     |    BSCHNG - Indicates whether basis has changed or not
         *      --------------------------------------------------------------
         * */
        bschng = ((_info->cbmode.smode || _info->supblk.basbnd) || first) || _info->cbmode.phase0;
        _info->supblk.sbchng = !bschng;
        reinvt = !((_info->cbmode.newpt || _info->cbmode.smode) || _info->cbmode.phase0);

        /*      -----------------------------------------------------
         *     | Set the logicals to see if XPVPCT needs to be changed
         *      -----------------------------------------------------
         * */
        if( _info->degcnt.idegk >= 20 )
                _info->degcnt.redpvt = _info->degcnt.redpvt || (_info->degcnt.idegfk >= _info->degcnt.idegk/
                 2);
        incpvt = _info->nwtcnt.inwtk >= 20 && (_info->nwtcnt.inwtfk >= _info->nwtcnt.inwtk/
         3);
        degprb = _info->degcnt.redpvt;
        /*C
         *      ----------------------------------------------------
         *     |  Compute Jacobian of the functions, scale the
         *     |  nonlinear elements if needed, and find the
         *     |  largest Jacobian element.
         *     |  NOTE:  GETGR just returns if all linear
         *      ----------------------------------------------------
         *
         *     chg 5/6/98 SHS -- Add PARSH output arrays to argument list of GETGR */
        if( _info->cbmode.compgr ){
                if( _info->lincnt.nnlin > 0 ){
                        _info->bestbk.objbst = g[_info->nintbk.nobj - 1];
                        jreset = FALSE;
                        getgr(_info, x, grad, ihag, iheg, ihegl, inlin, g, gg, (double*)iwork,
                         gcol, ub, alb, ascale, &jreset, iobjpt, paij, iprow,
                         ipcol, ipmap );
                        if( _info->glberr.abort )
                                return;
                        /*jcp        if(jreset) smode = .true. */
                        if( _info->nintbk.ninf != 0 )
                                g[_info->nintbk.nobj - 1] = _info->bestbk.objbst;
                        if( _info->scal.iscale > 0 )
                                scljac(_info, grad, ihag, iheg, ihegl, inlin, ascale, FALSE );
                }
                all = first && _info->scal.iscale == 0;
                if( _info->lincnt.nnlin > 0 || all )
                        chkelt(_info, grad, ihag, iheg, ihegl, inlin, all );
                if( _info->glberr.abort )
                        return;

                if( _info->ztol.bigelt < _info->limits.epspiv ){
                        bschng = TRUE;
#ifdef IO_ENABLED
                        if( ipconb > 2 )
                                {
                                sprintf(LSGRG_MSGBUFFER, " CONSBS...ERROR:  LARGEST ELEMENT IN GRAD ARRAY IS %16.8e\n  A NONSINGULAR BASIS MATRIX CAN NOT BE CHOSEN--USING AN ALL SLACK BASIS\n",
                                 _info->ztol.bigelt );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
                        sprintf(LSGRG_MSGBUFFER,
                         " CONSBS...ERROR:  LARGEST ELEMENT IN GRAD"
                         " ARRAY IS %16.8e\n  A NONSINGULAR BASIS MATRIX"
                         " CAN NOT BE CHOSEN--USING AN ALL SLACK BASIS\n",
                         _info->ztol.bigelt );
                   lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                }
        }

        /*      --------------------------------------
         *     |  Jacobian is now computed
         *      --------------------------------------
         *
         * */
        if( nbcp > 0 ){
                for( i = 1; i <= nbcp; i++ ){
                        i_ = i - 1;
                        iwork[i_] = ibc[i_];
                }
        }

        /*      -----------------------------------------------------------------
         *     | This loop determines index sets of binding and nonbinding
         *     | constraints IBC. Sets slacks of binding constraints.
         *      -----------------------------------------------------------------
         * */
        _info->bind.nbc = 0;
        _info->bind.nnbc = 0;
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                ii = istat[i_];
                islack = _info->dimen.n + i;

                /*     If ignored constraint or objective, skip
                 * */
                if( ii <= 0 )
                        goto L_40;

                if( ii == 2 )
                        goto L_20;

                /*     Equality Constraint
                 * */
                bnd = alb[islack - 1];
                if( fabs( g[i_] - bnd ) <= (_info->limits.epnewt*(1.0e0 + fabs( bnd ))) )
                        goto L_30;
                if( _info->cbmode.phase0 )
                        goto L_30;
                goto L_40;

                /*     Inequality Constraint - Check both bounds
                 * */
L_20:
                bnd = alb[islack - 1];
                if( fabs( g[i_] - bnd ) <= (_info->limits.epnewt*(1.0e0 + fabs( bnd ))) )
                        goto L_30;
                if( !_info->cbmode.phase0 )
                        goto L_25;
                if( g[i_] <= (bnd - _info->limits.epnewt*(1.0e0 + fabs( bnd ))) )
                        goto L_30;
                /*???         IF (G(I) .LE. (BND-EPNEWT*(1.0D0+DABS(BND)))) GO TO 40 */
L_25:
                bnd = ub[islack - 1];
                if( fabs( g[i_] - bnd ) <= (_info->limits.epnewt*(1.0e0 + fabs( bnd ))) )
                        goto L_30;
                if( !_info->cbmode.phase0 )
                        goto L_40;
                if( g[i_] >= (bnd + _info->limits.epnewt*(1.0e0 + fabs( bnd ))) )
                        goto L_30;
                goto L_40;

                /*         Infeasible or binding const, Set slack to closest bound
                 * */
L_30:
                _info->bind.nbc = _info->bind.nbc + 1;
                ibc[_info->bind.nbc - 1] = i;
                x[islack - 1] = bnd;
                goto L_50;

                /*         A feasible, non-binding constraint
                 * */
L_40:
                ibc[_info->dimen.mp1 - _info->bind.nnbc - 1] = i;
                _info->bind.nnbc = _info->bind.nnbc + 1;
                x[islack - 1] = g[i_];
L_50:
                ;
        }

        if( _info->bind.nbc > _info->dimen.nbmax )
                goto L_500;

#ifdef IO_ENABLED
        if( ipconb <= 4 )
                goto L_70;
        if( _info->bind.nbc > 0 )
                {
               lsgrg_msg(IOINFO, " BINDING CONSTRAINTS ARE " );
                for( i = 1; i <= _info->bind.nbc; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibc[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                }
        if( _info->bind.nnbc > 0 )
                {
               lsgrg_msg(IOINFO, " NONBINDING CONSTRAINTS ARE " );
                for( i = 1; i <= _info->bind.nnbc; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibc[_info->dimen.mp1 - i] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                }
L_70:
        ;
#endif

        /*C     -------------------------------------------
         *C    | Check if the basis must be changed
         *C    | RUPDT is logical to inhibit reseting of
         *C    | HESSIAN (R) when true
         *C     -------------------------------------------
         * */
        if( _info->bind.nbc != nbcp ){
                bschng = TRUE;
                _info->supblk.sbchng = FALSE;
        }
        /*C   -----------------------------------------------------
         *C  | DO 74 Loop checks for change in set of binding
         *C  | constraints, sets RUPDT to .FALSE. and, also forces
         *C  | BSCHNG to .TRUE.  RUPDT also depends upon a check for
         *C  | changes in the  set of basic variables, DO 75 stores
         *C  | current set in IWORK for later comparison.
         *C   ----------------------------------------------------- */
        if( (!_info->supblk.sbchng) || _info->bind.nbc == 0 )
                goto L_77;
        for( i = 1; i <= _info->bind.nbc; i++ ){
                i_ = i - 1;
                if( ibc[i_] != iwork[i_] )
                        _info->supblk.sbchng = FALSE;
        }
        bschng = bschng || (!_info->supblk.sbchng);
L_77:
        ;
#ifdef IO_ENABLED
        if( bschng && ipconb > 5 )
                {
                lsgrg_msg(IOINFO, " BASIS CHANGE\n" );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }

        if( ipconb > 6 ){
                sprintf(LSGRG_MSGBUFFER, " .... CONSBS DO 74  LOOP .....\n NBC = %5ld NBCP = %5ld RUPDT = %1c\n",
                 _info->bind.nbc, nbcp, TorF(_info->supblk.sbchng) );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
               lsgrg_msg(IOINFO,
                   " IBC   = " );
                for( i = 1; i <= _info->bind.nbc; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibc[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
               lsgrg_msg(IOINFO,
                   " IWORK = " );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                for( i = 1; i <= nbcp; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", iwork[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
        }
#endif

        /*C     ---------------------------------------------------------
         *     | Reset number of Basic Variables and other constants
         *     | so that basis size is MP1.
         *C     ---------------------------------------------------------
         *
         *      ----------------------------------------
         *C    |If same basis, go and Re-invert
         *      ----------------------------------------
         * */
        if( (bschng && reinvt) && (!_info->cbmode.pertub) ){
#ifdef IO_ENABLED
                sprintf( LSGRG_MSGBUFFER,
                 " ***ERROR...CONSBS --- RE-INVERT WAS TRUE BUT"
                 " BINDING CONSTRAINTS HAVE CHANGED\n" );
                lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
    /*           xerror( 0, iounit.ioerr );         */
        lsgrg_errorexit(JUMPBUF,_LSGRG_CONSBS_REINVERT_BC);
                return;
        }

        npnb = _info->dimen.npmp1;
        _info->nintbk.nb = _info->dimen.mp1;

        /*      --------------------------------------------------
         *C    | Save current set of Basic Variables
         *C     -------------------------------------------------
         * */
        if( nbp > 0 ){
                for( i = 1; i <= nbp; i++ ){
                        i_ = i - 1;
                        iwork[i_] = ibv[i_];
                }
        }

        /*C     ---------------------------------------------------------------------
         *     |                More Initializations:
         *     | ICOLS is used to indicate those variables that are at bounds
         *     | and those which have already been made Basic:
         *     |    ICOLS(J) =  0  --> Variable J at bound
         *     |    ICOLS(J) =  1  --> Variable J *not* at bound and *not* chosen
         *     |    ICOLS(J) = -1  --> Variable J fixed or already basic
         *     |    ICOLS(J) = -2  --> Variable J must *not* be made basic
         *C     ---------------------------------------------------------------------
         *
         *      -----------------------------------------------------
         *     | LOOP 78 - Come here if a basis was Ill-Conditioned
         *      -----------------------------------------------------
         * */
L_78:
        ;

        /*     --------------------------------------------
         *    | Reduce or increase XPVTPCT if necessary
         *     --------------------------------------------
         * */
        if( _info->pivots.fixpiv )
                goto L_1200;
        if( _info->degcnt.redpvt ){
                if( _info->pivots.ipvtpt > 2 ){
                        _info->pivots.ipvtpt = _info->pivots.ipvtpt - 1;
                        _info->pvtblk.xpvpct = _info->pivots.pivtol[_info->pivots.ipvtpt - 1];
#ifdef IO_ENABLED
                        if( ipconb > 0 )
                                {
                                sprintf( LSGRG_MSGBUFFER,
                                 "\n\n  DEGENERACY PROBLEMS -- PIVOT PERCENTAGE"
                                 " HAS BEEN REDUCED TO %12.4e\n",
                                 _info->pvtblk.xpvpct );
                                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                } else{
#ifdef IO_ENABLED
                        if( ipconb > 0 )
                                {
                                sprintf( LSGRG_MSGBUFFER,
                                 "\n\n  DEGENERACY PROBLEMS -- PIVOT PERCENTAGE"
                                 " IS AT ITS LOWER BOUND OF %12.4e\n",
                                 _info->pvtblk.xpvpct );
                                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                }
                _info->degcnt.redpvt = FALSE;
                _info->degcnt.idegk = 0;
                _info->degcnt.idegfk = 0;
                _info->degcnt.idegpt = 0;
        } else if( incpvt ){
                if( _info->pivots.ipvtpt < 7 ){
                        _info->pivots.ipvtpt = _info->pivots.ipvtpt + 1;
                        _info->pvtblk.xpvpct = _info->pivots.pivtol[_info->pivots.ipvtpt - 1];
#ifdef IO_ENABLED
                        if( ipconb > 0 )
                                {
                                sprintf( LSGRG_MSGBUFFER,
                                 "\n\n  NEWTON FAILURES -- PIVOT PERCENTAGE"
                                 " HAS BEEN INCREASED TO %12.4e\n",
                                 _info->pvtblk.xpvpct );
                                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                } else{
#ifdef IO_ENABLED
                        if( ipconb > 0 )
                                {
                                sprintf( LSGRG_MSGBUFFER,
                                "\n\n  NEWTON FAILURES -- PIVOT PERCENTAGE"
                                " IS AT ITS UPPER BOUND OF %12.4e\n",
                                 _info->pvtblk.xpvpct );
                                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                }
                _info->nwtcnt.inwtk = 0;
                _info->nwtcnt.inwtfk = 0;
                _info->nwtcnt.inwtpt = 0;
        }
        incpvt = FALSE;
        _info->degcnt.redpvt = FALSE;

L_1200:
        ;
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                icols[i_] = 0;
                if( ivstat[i_] < 0 )
                        icols[i_] = -1;

                /*      ----------------------------------------------
                 *     | Set distances of variables from nearest bound
                 *      ----------------------------------------------
                 * */
                if( alb[i_] <= -_info->tols.plinfy ){
                        d1 = _info->tols.plinfy;
                } else{
                        d1 = (x[i_] - alb[i_])/(1.0e0 + fabs( alb[i_] ));
                }
                if( ub[i_] >= _info->tols.plinfy ){
                        d2 = _info->tols.plinfy;
                } else{
                        d2 = (ub[i_] - x[i_])/(1.0e0 + fabs( ub[i_] ));
                }
                dd = fmin( d1, d2 );
                if( dd <= _info->limits.epboun )
                        dd = 0.0e0;
                if( dd > 1.0e-6*_info->tols.plinfy )
                        dd = _info->tols.plinfy*1.0e-6;
                dbnd[i_] = dd;
                if( dd > 0.0e0 )
                        icols[i_] = 1;
        }
#ifdef IO_ENABLED
        if( ipconb >= 6 )
                {
                lsgrg_msg(IOINFO, " DBND =" );
                for( i = 1; i <= _info->dimen.n; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%15.7e", dbnd[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                }
#endif


/* L_79: unreferenced??? */
        ;
        nchg = 0;
        if( !bschng ){
                _info->cbscnt.nreinv = _info->cbscnt.nreinv + 1;
#ifdef IO_ENABLED
                if( ipconb > 3 )
                        {
                        lsgrg_msg(IOINFO,
                         " ...CONSBS - ENTERING RE-INVERSION MODE ....\n" );
                        }
#endif
                lastpv = TRUE;
                newbas = FALSE;
                goto L_205;
        } else if( _info->cbmode.smode ){
                _info->cbscnt.nser = _info->cbscnt.nser + 1;
#ifdef IO_ENABLED
                if( ipconb > 4 )
                        {
                        lsgrg_msg(IOINFO,
                         " ...CONSBS - ENTERING SEARCH MODE ....\n" );
                        }
#endif
        } else{
                _info->cbscnt.nrser = _info->cbscnt.nrser + 1;
#ifdef IO_ENABLED
                if( ipconb > 4 )
                        {
                        lsgrg_msg(IOINFO,
                         " ...CONSBS - ENTERING RESTRICTED SEARCH MODE ....\n" );
                        }
#endif
        }

        /*      -----------------------------------------------------------
         *     |       **** SEARCH MODE STARTS HERE ****
         *     |
         *C    | Set up an initial basis:
         *     |   If in Search Mode use an all slack basis
         *     |   If in Restricted Search Mode, use old basis and
         *     |     mark all the active set changes
         *     |
         *     | ICOLS(N+1) .. ICOLS(N+NCHG) hold the rows which must
         *     |     be pivoted in
         *     |
         *     | NEWBAS - Used to indicate whether the starting Search
         *     |     Mode basis is the same as previous basis
         *      -----------------------------------------------------------
         *
         *
         *XXXX CALL CPUTIME(ZLAST,IRF) */
#ifdef IO_ENABLED
        if( ipconb >= 6 )
                {
                lsgrg_msg(IOINFO,
                 "  IBV FROM PREVIOUS ITERATION....\n " );
                for( i = 1; i <= nbp; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibv[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                }
#endif

        if( (_info->lincnt.nfix == _info->dimen.n || _info->bind.nbc == 0) || _info->ztol.bigelt <
         _info->limits.epspiv ){
                newbas = TRUE;
                nchg = 0;
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        ibv[i_] = _info->dimen.n + i;
                }
        } else if( _info->cbmode.smode ){
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        j = _info->dimen.n + i;
                        ibv[i_] = j;
                        dbnd[j - 1] = 0.0e0;
                }
                for( i = 1; i <= _info->bind.nbc; i++ ){
                        i_ = i - 1;
                        icols[_info->dimen.n + i_] = ibc[i_];
                }
                nchg = _info->bind.nbc;
                newbas = TRUE;
        } else{
                newbas = first;
                if( _info->bind.nnbc > 0 ){
                        for( i = 0; i <= (_info->bind.nnbc - 1); i++ ){
                                i_ = i - 1;
                                irow = ibc[_info->dimen.mp1 - i - 1];
                                islack = irow + _info->dimen.n;
                                if( ibv[irow - 1] != islack ){
                                        ibv[irow - 1] = islack;
                                        newbas = TRUE;
                                }
                        }
                }
                nchg = 0;
                for( i = 1; i <= _info->bind.nbc; i++ ){
                        i_ = i - 1;
                        irow = ibc[i_];
                        jcol = ibv[irow - 1];
                        if( jcol > _info->dimen.n ){
                                islack = _info->dimen.n + irow;
                                if( jcol != islack ){
                                        ibv[irow - 1] = islack;
                                        newbas = TRUE;
                                        dropbv = FALSE;
                                } else{
                                        dropbv = TRUE;
                                }
                        } else{
                                dropbv = (icols[jcol - 1] == 0) && (degprb || ibvbct[jcol - 1] <
                                 _info->pvtblk.ibvblm);
                        }

                        if( dropbv ){
                                nchg = nchg + 1;
                                icols[_info->dimen.n + nchg - 1] = irow;
                        } else if( jcol <= _info->dimen.n ){
                                icols[jcol - 1] = -1;
                        }
                }
        }


        /*      ---------------------------------------------------------------
         *     | If no binding constraints or no active changes, go factor basis
         *     |
         *     | LASTPV - Indicates that no more pivots have to be performed
         *      ---------------------------------------------------------------
         * */
        lastpv = _info->bind.nbc == 0 || nchg == 0;
        if( lastpv )
                goto L_205;
        /*C
         * */
#ifdef IO_ENABLED
        if( ipconb >= 6 ){
                sprintf(LSGRG_MSGBUFFER,
                 "  ENTERING SEARCH MODE  -- NUM OF CHANGES = %6ld\n",
                 nchg );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                lsgrg_msg(IOINFO,
                 " ORIGINAL IBV AT START OF SEARCH MODE....\n " );
                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibv[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                lsgrg_msg(IOINFO,
                 "  ICOLS.... (0 INDICATES AT BOUND)\n" );
                for( j = 1; j <= _info->dimen.n; j++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", icols[j - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                lsgrg_msg(IOINFO, "  PIVOT ROWS = \t" );
                for( j = _info->dimen.n + 1; j <= (_info->dimen.n + nchg); j++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", icols[j - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
        }
#endif
        epnorm = _info->limits.epspiv*_info->ztol.bigelt;
        if( epnorm > _info->limits.epspiv )
                epnorm = _info->limits.epspiv;

L_205:
#ifdef IO_ENABLED
        if( ipconb > 4 )
                {
                sprintf(LSGRG_MSGBUFFER, " NUMBER OF PIVOTS TO CHOOSE = %4ld\n",
                 nchg );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

        if( reinvt ){
                newbas = TRUE;
                lastpv = TRUE;
                umold = _info->lsinvt.dluprm[0];
                _info->lsinvt.dluprm[0] = 0.50e0;
        }
        if( !newbas ){
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        i_ = i - 1;
                        if( ibv[i_] == iwork[i_] )
                                goto L_4562;
#ifdef IO_ENABLED
                        sprintf( LSGRG_MSGBUFFER,
                        " CONSBS...ERROR: BASIS HAS CHANGED BEFORE INITIAL"
                        " XFACT CALL BUT BSCHNG WAS FALSE.\nBASIC VAR #%4ld"
                        " IS VAR #%4ld BUT WAS #%4ld\n",
                         i, ibv[i_], iwork[i_] );
                        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
    /*                   xerror( 0, iounit.ioerr );  */
        lsgrg_errorexit(JUMPBUF,_LSGRG_CONSBS_BASIS_STRUCTURE);
                        return;
L_4562:
                        ;
                }
        }

        /*      -----------------------------
         *C    | Factor the original basis
         *c  ** 8/96 jcp **  force reinvert if jacobian expanded
         *c    jreset is set in getgr/modjac
         *      -----------------------------
         * */
        if( (newbas || jreset) || !_info->equblk.lincon ){
                if( jreset )
                        newbas = TRUE;
    /*          cputime( &ztime, &irf ); */
                ztime = lsgrg_timer();

                xfact(_info, grad, ihag, iheg, ibmap, zmem, ibv, newbas, &rcode,
                 cdnum, irank, lastpv );
     /*         cputime( &ztemp, &irf ); */
                ztemp = lsgrg_timer();
                _info->temp.tfact = _info->temp.tfact + ztemp - ztime;
                if( reinvt )
                        _info->lsinvt.dluprm[0] = umold;


                /*         ------------------------------------------
                 *        : Check the return code from factorization
                 *         ------------------------------------------
                 * */
                if( rcode != 1 ){
                        /**********************
                         *   CHANGE THIS LATER
                         *        MAY WANT TO ALLOW IN ILL-CONDITIONED ORIGINAL FACTOR
                         ********************** */
                        if( !_info->cbmode.smode && rcode != 3 ){
#ifdef IO_ENABLED
                          if( ipconb >= 2 )
                                  {
                                  sprintf(LSGRG_MSGBUFFER,
                                   " CONSBS...WARNING: COULD NOT INVERT"
                                   " ORIGINAL BASIS MATRIX...FORCING CONSBS"
                                   " INTO SEARCH MODE \n RCODE FROM"
                                   " XFACT = %3ld\n",
                                   rcode );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                   }
#endif
                                _info->cbmode.smode = TRUE;
                                newbas = TRUE;
                                bschng = TRUE;
                                _info->cbscnt.nsing = _info->cbscnt.nsing + 1;
                                goto L_78;
                        } else{
#ifdef IO_ENABLED
                            sprintf( LSGRG_MSGBUFFER,
                             " CONSBS...ERROR: COULD NOT INVERT ORIGINAL"
                             " BASIS MATRIX WHILE IN SEARCH MODE \n RCODE"
                             " FROM XFACT = %3ld\n",
                                 rcode );
                            lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
    /*                           xerror( 0, iounit.ioerr ); */
        lsgrg_errorexit(JUMPBUF,_LSGRG_CONSBS_NOINVERT_SEARCH);
                                return;
                        }
                }
        }

        /*      ----------------------------------------------------------------
         *     | Start main pivot loop. -- Choose a basic variable for each
         *     | constraint marked above. "JPIV" will be basic variable for
         *     | row "IROW"
         *      ----------------------------------------------------------------
         * */
        if( lastpv )
                goto L_430;
  /*    cputime( &zlast, &irf ); */
        zlast = lsgrg_timer();
        for( i = 1; i <= nchg; i++ ){
                i_ = i - 1;
                jpiv = 0;
                irow = icols[_info->dimen.n + i_];

                /*        -------------------------------------------------------------
                 *       | Row "IROW" of BINV = (UNIT VECTOR) * BINV    (Held in GG)
                 *        -------------------------------------------------------------
                 * */
                for( j = 1; j <= _info->nintbk.nb; j++ ){
                        j_ = j - 1;
                        gg[j_] = 0.0e0;
                }
                gg[irow - 1] = 1.0e0;

                xbtran(_info, ibmap, zmem, gg, _info->nintbk.nb );

#ifdef IO_ENABLED
                if( ipconb >= 6 )
                        {
                        sprintf(LSGRG_MSGBUFFER, "  ROW %5ld", irow );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        lsgrg_msg(IOINFO, " OF BINV IS:\n " );
                        for( j = 1; j <= _info->nintbk.nb; j++ ){
                                sprintf(LSGRG_MSGBUFFER, "%12.6e", gg[j - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                        sprintf(LSGRG_MSGBUFFER, " \n" );
                        }
#endif


                /*        ----------------------------------------------------
                 *       | Find largest pivot element from variables
                 *       | away from bounds
                 *        ----------------------------------------------------
                 * */
                eltmax = 0.0e0;
                for( j = 1; j <= _info->dimen.n; j++ ){
                        j_ = j - 1;
                        pvrow[j_] = 0.0e0;
                        if( icols[j_] != 1 )
                                goto L_230;
                        xdot(_info, grad, ihag, iheg, gg, _info->nintbk.nb, j, &sum );
                        pvrow[j_] = sum;
                        sum = fabs( sum );
                        if( eltmax < sum )
                                eltmax = sum;
L_230:
                        ;
                }

#ifdef IO_ENABLED
                if( ipconb >= 6 )
                        {
                        sprintf(LSGRG_MSGBUFFER, "  UPDATED PIVOT ROW %5ld",
                         irow );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        lsgrg_msg(IOINFO, " IS:\n " );
                        for( j = 1; j <= _info->dimen.n; j++ ){
                                sprintf(LSGRG_MSGBUFFER, "%12.6e", pvrow[j - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                        sprintf(LSGRG_MSGBUFFER, " \n" );
                        }
#endif

                /*         -------------------------------------------------------
                 *        | If all pivot elements too small, go look at variables
                 *        | at bounds.
                 *         -------------------------------------------------------
                 * */
                if( eltmax < epnorm )
                        goto L_250;

                /*         -----------------------------------------------------------
                 *        | Find variable farthest from bound and with pivot element
                 *        | at least XPVPCT*ELTMAX.  Note, there will always be one with
                 *        | ELTMAX.
                 *         -----------------------------------------------------------
                 * */
                tol = fmax( _info->pvtblk.xpvpct*eltmax, epnorm );
                dmax = 0.0e0;
                for( j = 1; j <= _info->dimen.n; j++ ){
                        j_ = j - 1;
                        if( fabs( pvrow[j_] ) < tol )
                                goto L_240;
                        if( dmax >= dbnd[j_] )
                                goto L_240;
                        dmax = dbnd[j_];
                        jpiv = j;
L_240:
                        ;
                }

                /*         ----------------------------------------------------
                 *        | Found a basic coulmn, so go pivot it into basis
                 *         ----------------------------------------------------
                 * */
                goto L_270;

                /*         -----------------------------------------------------
                 *        | Come here to check variables at bounds if needed.
                 *        | Look for variable at bound with largest pivot element.
                 *         -----------------------------------------------------
                 * */
L_250:
                ;
#ifdef IO_ENABLED
                if( ipconb >= 6 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " LOOKING AT VARS AT BOUNDS IN ROW %5ld ELTMAX OF %15.8e  IS LESS THAN EPNORM = %14.7e\n",
                         irow, eltmax, epnorm );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                eltmax = 0.0e0;
                for( j = 1; j <= _info->dimen.n; j++ ){
                        j_ = j - 1;
                        if( icols[j_] != 0 )
                                goto L_260;
                        xdot(_info, grad, ihag, iheg, gg, _info->nintbk.nb, j, &sum );
                        pvrow[j_] = sum;
                        sum = fabs( sum );
                        if( eltmax > sum )
                                goto L_260;
                        eltmax = sum;
                        jpiv = j;
L_260:
                        ;
                }

                /*        --------------------------------------
                 *       | If pivot is large enough, use it.
                 *        --------------------------------------
                 * */
                if( eltmax > epnorm )
                        goto L_270;

                /*        --------------------------------------------------
                 *       | NO pivot could be chosen,  use slack at bound.
                 *        --------------------------------------------------
                 * */
                islack = irow + _info->dimen.n;
#ifdef IO_ENABLED
                if( ipconb >= 6 )
                        {
                        sprintf(LSGRG_MSGBUFFER,
                         " NO PIVOTS > ESPPIV  ---  PIVOT ON SLACK"
                         "   %5ld\n VARIABLE NO. %5ld\n",
                         irow, islack );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( ibv[irow - 1] == islack )
                        goto L_210;
                for( j = 1; j <= _info->dimen.mp1; j++ ){
                        j_ = j - 1;
                        if( ibv[j_] == islack )
                                goto L_313;
                }
                goto L_314;
L_313:
                ;
                if( _info->cbmode.smode ){
#ifdef IO_ENABLED
                        sprintf( LSGRG_MSGBUFFER,
                         " CONSBS...ERROR: THE SLACK FOR ROW %4ld IS"
                         " ALREADY BASIC IN COL %4ld\t COULD NOT CHOOSE"
                         " A PIVOT FOR ROW %4ld IN SEARCH MODE\n",
                         irow, j, irow );
                        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
    /*                   xerror( 0, iounit.ioerr );  */
        lsgrg_errorexit(JUMPBUF,_LSGRG_CONSBS_BASIC_SLACK);
                        return;
                } else{
#ifdef IO_ENABLED
                        if( ipconb >= 2 )
                                {
                                sprintf(LSGRG_MSGBUFFER,
                                 " CONSBS...COULD NOT CHOOSE A PIVOT FOR"
                                 " ROW %4ld\n          AND SLACK IS ALREADY"
                                 " BASIC IN COL %4ld  --- FORCING CONSBS"
                                 " INTO SEARCH MODE\n",
                                 irow, j );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                        _info->cbmode.smode = TRUE;
                        newbas = TRUE;
                        bschng = TRUE;
                        _info->cbscnt.nsing = _info->cbscnt.nsing + 1;
                        incpvt = TRUE;
                        _info->degcnt.redpvt = FALSE;
                        goto L_78;
                }

L_314:
                jpiv = islack;


                /*        ------------------------------------------------------
                 *       |  Pivot has been chosen --- JPIV should be > 0
                 *        ------------------------------------------------------
                 * */
L_270:
                ;

                if( jpiv <= 0 ){
#ifdef IO_ENABLED
                    sprintf( LSGRG_MSGBUFFER,
                     " CONSBS...ERROR: JPIV CAN NOT BE ZERO, BUT IS \n" );
                    lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
    /*                   xerror( 0, iounit.ioerr );  */
        lsgrg_errorexit(JUMPBUF,_LSGRG_CONSBS_JPIV_0);
                        return;
                }
#ifdef IO_ENABLED
                if( ipconb >= 6 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " PIVOT ROW IS %6ld  MAXIMUM ELEMENT = %14.7e\n",
                         irow, eltmax );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                if( ipconb >= 6 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " PIVOT COLUMN IS %6ld  PIVOT ELEMENT = %14.7e  DISTANCE FROM BOUND = %13.7e\n",
                         jpiv, pvrow[jpiv - 1], dbnd[jpiv - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif


                /*          -------------------------
                 *         : Update basis factors
                 *          -------------------------
                 * */
                lastpv = i == nchg;
                xpivot(_info, grad, ihag, iheg, ibmap, zmem, ibv, jpiv, irow, gg,
                 &rcode, cdnum, irank, lastpv );

                /*         --------------------------
                 *        |  Check return code
                 *  ** 8/01 jcp use return_status as flag to avoid
                 *    redundant error msgs
                 *         --------------------------
                 * */
                lsgrg_return_status = _LSGRG_STATUS_NOT_SET;
                switch( rcode ){
                        case 1: goto L_211;
                        case 2:
                           lsgrg_return_status = _LSGRG_XPIVOT_BASIS_SING;
                           goto L_212;
                        case 3:
                           lsgrg_return_status = _LSGRG_XPIVOT_INSFMEM;
                           goto L_214;
                        case 4:
                           lsgrg_return_status = _LSGRG_XPIVOT_OTHER_ERR;
                           goto L_215;
                        case 5:
                           lsgrg_return_status = _LSGRG_XPIVOT_BASIS_ILLCOND;
                           goto L_213;
                        }

                /*         ----------------------------------------------------------
                 *        | Singular Matrix -- Throw out chosen variable and try again
                 *        |                    Force a Search Mode.
                 *         ---------------------------------------------------------- */
L_212:
                ;
                if( _info->cbmode.smode ){
#ifdef IO_ENABLED
                  sprintf( LSGRG_MSGBUFFER,
                   " CONSBS...ERROR: CHOOSING COLUMN %6ld CREATES A"
                   " SINGULAR BASIS WHILE IN SEARCH MODE \n",
                         jpiv );
                  lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                        goto L_215;
                } else{
#ifdef IO_ENABLED
                   sprintf( LSGRG_MSGBUFFER,
                    " CONSBS...WARNING: CHOOSING COLUMN %6ld CREATES A"
                    " SINGULAR BASIS---FORCING CONSBS INTO SEARCH MODE \n",
                         jpiv );
                   lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                        _info->cbmode.smode = TRUE;
                        newbas = TRUE;
                        bschng = TRUE;
                        _info->cbscnt.nsing = _info->cbscnt.nsing + 1;
                        incpvt = TRUE;
                        _info->degcnt.redpvt = FALSE;
                        goto L_78;
                }

                /*          -----------------------------------------------
                 *         : Ill-Conditioned Matrix -- Force a Search Mode
                 *          -----------------------------------------------
                 * */
L_213:
                ;

                if( lastpv ){
                        if( _info->cbmode.smode ){

                                /*                 -----------------------------------------------------
                                 *                : Just print ERROR and continue if inverted matrix
                                 *                : was Ill-Conditioned in Search Mode
                                 *                 -----------------------------------------------------
                                 * */
#ifdef IO_ENABLED
                     if(ipconb > 1) {
                           sprintf( LSGRG_MSGBUFFER,
                            " CONSBS...WARNING: SEARCH MODE BASIS MATRIX IS"
                            " ILL-CONDITIONED \n\n" );
                           lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
                      }
#endif
                                goto L_211;
                        } else{
#ifdef IO_ENABLED
                                if( ipconb > 1 )
                                        {
                                        sprintf(LSGRG_MSGBUFFER,
                                         " CONSBS...CURRENT BASIS IS ILL"
                                         " CONDITIONED -- FORCING BASIS"
                                         " CHANGE \n" );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                        }
#endif
                                _info->cbscnt.nicond = _info->cbscnt.nicond + 1;
                                _info->cbmode.smode = TRUE;
                                newbas = TRUE;
                                bschng = TRUE;
                                goto L_78;
                        }
                } else{
#ifdef IO_ENABLED
               if(ipconb > 1) {
                  sprintf( LSGRG_MSGBUFFER,
                   " CONSBS...ERROR: XPIVOT CALCULATED A CONDITION"
                   " NUMBER WHEN LASTPV WAS FALSE\n" );
                  lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
                        goto L_215;
                }

                /*         ----------------------------------------------------------
                 *        | Storage Overflow -- Print an ERROR message and STOP
                 *  ** jcp 8/01  this msg is redundant, dont print on ipr=1
                 *         ---------------------------------------------------------- */
L_214:
                ;
#ifdef IO_ENABLED
            if(ipconb > 1) {
                sprintf( LSGRG_MSGBUFFER,
                 " CONSBS...ERROR: BASIS INVERSE STORAGE OVERFLOW  \n" );
                lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
            }
#endif

                /*         ------------------------
                 *        | Other ERROR returns
                 *  ** 8/01 jcp  dont print generic  return code msg
                 *  if basis inverse overflow or low ipr
                 *         ------------------------
                 * */
L_215:
                ;
#ifdef IO_ENABLED
              if(lsgrg_return_status != _LSGRG_XPIVOT_INSFMEM
                   || ipconb > 2 ) {
                  sprintf( LSGRG_MSGBUFFER,
                    " CONSBS...ERROR: RETURN FROM XPIVOT - RCODE = %5ld\n",
                    rcode );
                   lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
              }
#endif
    /*           xerror( 0, iounit.ioerr );        */
                 lsgrg_errorexit(JUMPBUF,lsgrg_return_status);
                return;




                /*         ------------------------------
                 *        | Inversion was successfull
                 *         ------------------------------
                 * */
L_211:
                ;
                if( jpiv <= _info->dimen.n )
                        icols[jpiv - 1] = -1;

L_210:
                ;
        }



#ifdef IO_ENABLED
        if( ipconb > 3 )
                {
               lsgrg_msg(IOINFO,
                " ...CONSBS - EXITING  SEARCH MODE ....\n" );
                }
#endif
  /*    cputime( &ztime, &irf ); */
        ztime = lsgrg_timer();
        _info->temp.tsear = _info->temp.tsear + ztime - zlast;


        /*      End of basis selection for binding constraints
         *
         *      **** END OF SEARCH MODE ***
         *
         *###################################################################### */

L_430:
        ;

        /*C   -------------------------------------------------------------
         *C  | Check for changes in set of basic variables.
         *C  | Set RUPDT to .FALSE. if changes have taken place
         *C   -------------------------------------------------------------
         *C */
        if( !bschng )
                goto L_498;
        if( !_info->supblk.sbchng )
                goto L_498;
        for( j = 1; j <= npnb; j++ ){
                j_ = j - 1;
                icols[j_] = 0;
        }
        for( j = 1; j <= _info->nintbk.nb; j++ ){
                j_ = j - 1;
                icols[ibv[j_] - 1] = 1;
        }
        for( j = 1; j <= nbp; j++ ){
                j_ = j - 1;
                jj = iwork[j_];
                if( icols[jj - 1] == 1 )
                        goto L_497;
                _info->supblk.sbchng = FALSE;
                goto L_498;
L_497:
                ;
        }


L_498:
        ;

        /*     -----------------------------------------------------------
         *    | Note -- Only change INBV if RUPDT is FALSE; All other
         *    |         changes in Non-Basics is done in DIREC
         *     -----------------------------------------------------------
         * */
        if( _info->supblk.sbchng )
                goto L_492;

        /*C     -------------------------------------------------------------
         *     | Set up index set of nonbasic variables, Superbasics first and
         *     | Nonbasic slacks last
         *C     -------------------------------------------------------------
         * */
        for( i = 1; i <= npnb; i++ ){
                i_ = i - 1;
                icols[i_] = 0;
        }
        if( _info->nintbk.nb == 0 )
                goto L_448;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                j = ibv[i_];
                icand[i_] = j;
                icols[j - 1] = 1;
        }

L_448:
        ;

        /*      ---------------------
         *     | Do Nonbasic slacks
         *      --------------------- */
        ibctr = _info->dimen.n;
        ictrb = _info->dimen.n + _info->nintbk.nb;
        if( _info->nintbk.nb == 0 )
                goto L_470;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                k = _info->dimen.n + i;
                if( icols[k - 1] == 1 )
                        goto L_460;
                icand[ictrb - 1] = k;
                ictrb = ictrb - 1;
                inbv[ibctr - 1] = k;
                ibctr = ibctr - 1;
L_460:
                ;
        }
        /*      --------------------------
         *     | Handle rest of Nonbasics
         *      -------------------------- */
L_470:
        ;
        icsupe = _info->nintbk.nb;
        _info->nintbk.nsuper = 0;
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                if( icols[i_] == 1 )
                        goto L_490;
                if( dbnd[i_] > 0.0e0 )
                        goto L_480;
                inbv[ibctr - 1] = i;
                ibctr = ibctr - 1;
                icand[ictrb - 1] = i;
                ictrb = ictrb - 1;
                goto L_490;
L_480:
                ;
                /*           ---------------------
                 *          | Superbasic Variables
                 *           --------------------- */
                icsupe = icsupe + 1;
                icand[icsupe - 1] = i;
                _info->nintbk.nsuper = _info->nintbk.nsuper + 1;
                inbv[_info->nintbk.nsuper - 1] = i;
L_490:
                ;
        }

L_492:
        ;

        /*      ----------------------------------
         *C    | Check for variables at bounds
         *      ----------------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                jcol = ibv[i_];
                if( jcol > _info->dimen.n )
                        goto L_505;
                if( dbnd[jcol - 1] != 0.0e0 || degprb ){
                        ibvbct[jcol - 1] = 0;
                } else{
                        ibvbct[jcol - 1] = ibvbct[jcol - 1] + 1;
                }
L_505:
                ;
        }


        /*      ----------------------------------------------------------
         *     |  Now set the bound status of Nonbasic Variables
         *     |  and count the number of Nonbasics not at bounds
         *      ----------------------------------------------------------
         * */
        *nbfree = 0;
        tol = _info->limits.epboun;
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                iub[i_] = 0;
                k = inbv[i_];
                if( x[k - 1] >= (ub[k - 1] - (tol*(1.0e0 + fabs( ub[k - 1] )))) ){
                        iub[i_] = 1;
                } else if( x[k - 1] <= (alb[k - 1] + (tol*(1.0e0 + fabs( alb[k - 1] )))) ){
                        iub[i_] = -1;
                } else{
                        *nbfree = *nbfree + 1;
                }
                if( k <= _info->dimen.n )
                        ibvbct[k - 1] = 0;
        }

#ifdef IO_ENABLED
        if( ipconb > 6 ){
                sprintf(LSGRG_MSGBUFFER, " .... CONSBS DO 497 LOOP .....\n NB = %5ld NBP = %5ld RUPDT = %1c\n",
                 _info->nintbk.nb, nbp, TorF(_info->supblk.sbchng) );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                sprintf(LSGRG_MSGBUFFER, " %s", "IBV   = " );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibv[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
               lsgrg_msg(IOINFO, " IWORK = " );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                for( i = 1; i <= nbp; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", iwork[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
        }
        if( ipcb3 <= 3 && ipconb <= 4 )
                goto L_9987;
        sprintf(LSGRG_MSGBUFFER, " NUMBER OF SUPERBASICS = %4ld\n",
         _info->nintbk.nsuper );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
       lsgrg_msg(IOINFO, " INBV IS" );
        for( i = 1; i <= _info->dimen.n; i++ ){
                sprintf(LSGRG_MSGBUFFER, "%4ld", inbv[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
        lsgrg_msg(IOINFO,"\n" );
       lsgrg_msg(IOINFO, " IUB IS " );
        for( i = 1; i <= _info->dimen.n; i++ ){
                sprintf(LSGRG_MSGBUFFER, "%4ld", iub[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
        lsgrg_msg(IOINFO,"\n" );
        if( _info->nintbk.nb > 0 )
                {
               lsgrg_msg(IOINFO, " IBV IS" );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibv[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
                }

L_9987:
        ;

        if( ipconb >= 5 || _info->dbug.debug )
                {
               sprintf(LSGRG_MSGBUFFER,
                " ***** EXITING  SUBROUITNE: CONSBS\n" );
                }
#endif

        /*      -----------------------------------
         *     | Reset SMODE and save PRSMOD
         *      -----------------------------------
         * */
        _info->cbmode.prsmod = _info->cbmode.smode;
        _info->cbmode.smode = FALSE;
        _info->degcnt.redpvt = FALSE;
        _info->supblk.basbnd = FALSE;

        /*      ------------------
         *     | Get final time
         *      ------------------
         * */
  /*    cputime( &ztime, &irf ); */
        ztime = lsgrg_timer();
        _info->contim.cbtime = _info->contim.cbtime + ztime - zstart;

        return;
        /*      -----------------------------------------
         *     | Number of binding constraints too large
         *      ----------------------------------------- */
L_500:
#ifdef IO_ENABLED
        sprintf( LSGRG_MSGBUFFER,
         "\n\n CONSBS -- ERROR: LIMIT ON NUMBER OF BINDING"
         " CONSTRAINTS (%5ld) EXCEEDED  nb = %5ld\n",
         _info->dimen.nbmax,_info->nintbk.nb );
        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
    /*   xerror( 0, iounit.ioerr );                */
        lsgrg_errorexit(JUMPBUF,_LSGRG_CONSBS_NB_TOOBIG);
        return;


        /*     End of CONSBS
         * */
} /* end of function */


/**DECK GETGR */

/*     chg 5/6/98 SHS - Add PARSH output arrays to argument list */

void /*FUNCTION*/ LSGRGCLASS getgr(LsgrgInfo *_info,double x[], double grad[], long int ihag[],
         long int iheg[], long int ihegl[], long int inlin[], double g[],
         double gplus[], double gminus[], double gcol[], double ub[],
         double alb[], double ascale[], LOGICAL32 *jreset, long int iobjpt[],
         double paij[], long int iprow[], long int ipcol[], long int ipmap[])
{
long int i, i_, iparsh, /* irf,*/  irow, /* irs,*/ j, j_, jnl, k1, k2, /* maxnnz, */
         nelts, nnewnz, nz, nztot, ldummy;
double cmax0, gvabs, gval, upbnd, zertol, zfinal, zstart;




        /* ....................................................................
         *
         *     *** START OF GETGR ***
         *
         *      -----------------------------------------
         *     |  If no nonlinear elements, just return
         *      ----------------------------------------- */
        if( _info->lincnt.nnlin == 0 )
                return;


        /*      ------------------
         *     | Get start time
         *      ------------------
         *
         *     change 4/9/98 shs
         *     Use maxgro as parameter to modjac so that limit is decremented
         *     maxgro == space remaining for expansion in grad
         *xxx  maxnnz = maxgro */
  /*    cputime( &zstart, &irs ); */
        zstart = lsgrg_timer();
        _info->itnggr.ngetgr = _info->itnggr.ngetgr + 1;

#ifdef IO_ENABLED
        if( _info->dbug.debug )
                {
                sprintf(LSGRG_MSGBUFFER,
                 "\n ***** ENTERING SUBROUTINE: GETGR  - NUMBER OF NONLINEAR VARIABLES = %6ld\n\n\n\t *** DUMP OF NONZERO PATTERN ***\n\n",
                 _info->lincnt.nnlin );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        _info->counts.ngrad = _info->counts.ngrad + 1;


        /*      --------------------------------------
         *     | If data is scaled, unscale X and G
         *      --------------------------------------
         * */

        if( _info->scal.scaled ){
                for( i = 1; i <= _info->dimen.n; i++ ){
                        i_ = i - 1;
                        x[i_] = x[i_]*ascale[i_];
                }
                if( _info->pardat.kderiv <= 1 ){
                        for( i = 1; i <= _info->dimen.mp1; i++ ){
                                i_ = i - 1;
                                g[i_] = g[i_]/ascale[_info->dimen.n + i_];
                        }
                        _info->bestbk.truobj = _info->bestbk.truobj/ascale[_info->dimen.n + _info->nintbk.nobj - 1];
                }
        }


        /*      ------------------------------------------------------
         *     | KDERIV = 2 --> Use the Analytical Derivatives
         *     |                Store the derivatives directly in GRAD
         *     |
         *     | KDERIV = 3 --> GAMS interpreted derivatives
         *     |  10/98 jcp ifdef out parsh_gams call until/if interface
         *     |         reimplemented for c version
         *      ------------------------------------------------------
         *
         *
         *     chg 5/6/98 SHS - Add  Analatical PARSH call for kderiv=2 */
        if( _info->pardat.kderiv == 3 ){
#ifdef GAMS_INTERFACE
                parsh_gams( g, x, &_info->dimen.mp1, &_info->dimen.n, grad );
#endif
        } else if( _info->pardat.kderiv == 2 ){
/* call analytic parsh wrapper for 0-based arrays */

               lsgrg_parsh0( _info, x, _info->dimen.n, grad,
                 paij, iprow, ipcol,
                 ipmap, _info->nzerog.nzgrad );

                if( _info->glberr.abort )
                        return;
        } else{
                /*         -------------------------------------------
                 *        | Must use finite differencing
                 *        | Set branch for derivative option
                 *         ------------------------------------------
                 * */
                if( _info->pardat.kderiv == 0 )
                        iparsh = 21;
                if( _info->pardat.kderiv == 1 )
                        iparsh = 22;

                /*         ------------------------------------------------------
                 *        | If scaling then unscaling causes small deriv. errors
                 *        | so adjust ZERTOL
                 *         ------------------------------------------------------
                 *
                 *   chg 3/31/98 adjust zertol if recomputing scale factors to account
                 *     for finite differencing errors */

                zertol = _info->ztol.xajtol;
                if( _info->scal.scaled || _info->scal.newscl ){
                        zertol = zertol/_info->stepbk.pstep;
                        cmax0 = 0.0e0;
                        for( i = 1; i <= _info->dimen.mp1; i++ ){
                                i_ = i - 1;
                                cmax0 = fmax( cmax0, fabs( g[i_] ) );
                        }
                        zertol = zertol*cmax0;
                }

                /*         --------------------------------------------------
                 *        | Call the appropriate PARSH to get elelments
                 *        | Compute the Jacobian for each nonlinear column
                 *         --------------------------------------------------
                 * */
                nztot = 0;
                for( j = 1; j <= _info->lincnt.nnlin; j++ ){
                        j_ = j - 1;
                        jnl = inlin[j_];
                        switch( iparsh ){
                                case 21: goto L_21;
                                case 22: goto L_22;
                                }
L_21:
                        upbnd = ub[jnl - 1];
                        if( _info->scal.scaled )
                                upbnd = upbnd*ascale[jnl - 1];
                        lsgrg_parshf0( _info,x, g, gcol, jnl, gplus, upbnd );
                        goto L_25;
L_22:
           /*  lsgrgc3 parshc does not have the 'g' argument */
           /*           parshc( x, g, gcol, jnl, gplus, gminus ); */
                        lsgrg_parshc0( _info, x,  gcol, jnl, gplus, gminus );
L_25:
                        ;

                        /*            -----------------------------
                         *C          | Store predicted nonzeroes
                         *C           -----------------------------
                         *c   ** 6/96 jcp  .. provide check for new nonzeros by storing copy
                         *c                   of gcol in gminus
                         *C */
                        for( i = 1; i <= _info->dimen.mp1; i++ ){
                                i_ = i - 1;
                                gminus[i_] = gcol[i_];
                        }
                        /*c */
                        k2 = iheg[jnl] - 1;
                        k1 = iheg[jnl - 1];
                        nelts = k2 - k1 + 1;
                        nz = 0;
                        if( nelts <= 0 )
                                goto L_45;
                        for( i = k1; i <= k2; i++ ){
                                i_ = i - 1;
                                irow = ihag[i_];
                                gval = gcol[irow - 1];
                                gvabs = fabs( gval );
                                grad[i_] = gval;
                                if( gvabs >= zertol )
                                        nz = nz + 1;
                                /*XXXXX         IF (IROW .EQ. NOBJ .AND. MAXIM) GRAD(I) = -GVAL */
                                gcol[irow - 1] = 0.0e0;
                        }
L_45:
                        ;
                        nztot = nztot + nz;

                        /*C           -------------------------------------------
                         *C          | Check if any new nonzeroes were created
                         *C           -------------------------------------------
                         * */
                        nnewnz = 0;
                        for( i = 1; i <= _info->dimen.mp1; i++ ){
                                i_ = i - 1;
                                if( fabs( gcol[i_] ) < zertol )
                                        goto L_5035;
                                nnewnz = nnewnz + 1;
#ifdef IO_ENABLED
                                if( _info->nintbk.ipr >  1 )
                                    {
                                    sprintf(LSGRG_MSGBUFFER,
                                     "\n In GETGR,New Nonzero Detected In"
                                     " Col = %6ld   Row = %6ld\n"
                                     "   Value = %14.6e ZERTOL = %14.6e\n",
                                    jnl, i, gcol[i_], zertol );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                    }
#endif
L_5035:
                                ;
                        }
                        /*4059       FORMAT(/,' IN GETGR AFTER DO 5035 LOOP NNEWNZ = ',I3,/)
                         *           WRITE(IOOUT,4059)NNEWNZ
                         *c
                         *c   ** 6/96 jcp *** call modjac to modify jacobian
                         *c
                         *c     change 4/9/98 shs
                         *c     use maxgro as parameter into modjac, replacing maxnnz */
  /* jcp 10/98 replace ioout arg with long dummy to remove */
  /* FILE * arg */
/* 05/01/02 jcp **  set ldummy to 0 to quiet uninitalized var warning */
                                                ldummy = 0;
                        if( nnewnz > 0 ){
                                modjac(_info, ldummy      , nnewnz, jnl, gminus, gcol, iheg,
                                 ihag, grad, &_info->jgrow.maxgro, _info->dimen.n, zertol, _info->dimen.mp1,
                                 _info->nintbk.ipr, _info->dbug.debug, &_info->cmax.colmax, ihegl, iobjpt );
                                /*xxx *      maxnnz,n,zertol,mp1,ipr,debug,colmax,ihegl,iobjpt) */
                                *jreset = TRUE;
                                _info->nzerog.nzgrad = _info->nzerog.nzgrad + nnewnz;
                                _info->memory.lgrad = _info->memory.lgrad + nnewnz;
                                _info->nzerog.nznlin = _info->nzerog.nznlin + nnewnz;
                        }

                }

                /*C        --------------------------------------
                 *C       | If new nonzeros, just STOP for now
                 *C        --------------------------------------
                 *c
                 * */
#ifdef IO_ENABLED
                if( _info->dbug.debug )
                        {
                        sprintf(LSGRG_MSGBUFFER, "\n ***** NONZERO NONLINEAR ELEMENTS:   PREDICTED = %6ld     ACTUAL = %6ld\n\n ***** EXITING  SUBROUTINE: GETGR \n",
                         _info->nzerog.nznlin, nztot );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
        }

        /*      ------------------------------------
         *     |  Jacobian has now been computed.
         *      ------------------------------------ */


        /*      ------------------------------------------
         *     | If data was scaled,  Re-Scale  X and G
         *      ------------------------------------------ */

        /*     CALL FDCHECK(X,GRAD,IHAG,IHEG,G,GPLUS,GMINUS,GCOL,
         *    1                   NZGRAD,NP1) */
        if( _info->scal.scaled ){

                for( i = 1; i <= _info->dimen.n; i++ ){
                        i_ = i - 1;
                        x[i_] = x[i_]/ascale[i_];
                }

                if( _info->pardat.kderiv <= 1 ){
                        for( i = 1; i <= _info->dimen.mp1; i++ ){
                                i_ = i - 1;
                                g[i_] = g[i_]*ascale[_info->dimen.n + i_];
                        }
                        _info->bestbk.truobj = _info->bestbk.truobj*ascale[_info->dimen.n + _info->nintbk.nobj - 1];
                }

        }

        /*C     --------------------------------------
         *C    | Print GRAD array if IPN4 < 0
         *C     --------------------------------------
         * */
#ifdef IO_ENABLED
        if( _info->mngrg.ipn4 < 0 || _info->dbug.debug ){
                sprintf(LSGRG_MSGBUFFER, "\n\nGRAD ARRAY IS \n" );
                grdump(_info, grad, ihag, iheg, ihegl, ihag, _info->dimen.mp1, _info->dimen.n,
                 _info->nzerog.nzgrad, _info->memory.np1, ldummy, "GRAD ARRAY IN GETGR"
                  );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        /*      -------------------------------
         *     | Compute Total PARSH time
         *      -------------------------------
         * */
        _info->cbmode.compgr = FALSE;
        /*XXX  NEWPT=.FALSE. */
  /*    cputime( &zfinal, &irf ); */
        zfinal = lsgrg_timer();
        _info->gradtm.tgrad = _info->gradtm.tgrad + zfinal - zstart;

        return;




        /*     *** End of GETGR ***
         * */
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS grdump(LsgrgInfo *_info,double a[], long int ha[], long int he[],
         long int hel[], long int jtype[], long int nb, long int n, long int nz,
         long int np1, long int ioout, char *msg)
{
//long int i, i_, iend, irow, istart, j, j_;
//double coeff;

        /*----------------------------------------------------------------------
         *             *** START OF GRDUMP ***
         * */
#ifdef IO_ENABLED
       sprintf(LSGRG_MSGBUFFER,
        "\n\n==================================================================================================================================\n\n"
"                                              SPARSE JACOBIAN DUMP\n\n                     CAL"
"==> %s\n\n          NUMBER OF ROWS = \t%5ld\n          "
"NUMBER OF VARIABLES = \t%5ld\nNUMBER OF NON-ZEROES = \t%5ld\n",
         msg, nb, n, nz );
       lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);


       lsgrg_msg(IOINFO, "\n\n THE COLUMN POINTERS HE AND HEL ARE: \n\n\n" );
       lsgrg_msg(IOINFO, " J = " );
        for( j = 1; j <= n; j++ ){
                sprintf(LSGRG_MSGBUFFER, "%5ld HE(J) = %5ld HEL(J) = %5ld", j,
                 he[j - 1], hel[j - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
        lsgrg_msg(IOINFO,"\n" );
        sprintf(LSGRG_MSGBUFFER, " J = %5ld HE(J) = %5ld\n", np1, he[np1 - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        for( j = 1; j <= n; j++ ){
                j_ = j - 1;
                istart = he[j_];
                iend = he[j_ + 1] - 1;
                if( iend < istart ){
                        sprintf(LSGRG_MSGBUFFER, "\n   COL = \t%5ld**** HAS ZERO LENGTH ***\n\n",
                         j );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        goto L_21;
                }
                for( i = istart; i <= iend; i++ ){
                        i_ = i - 1;
                        irow = ha[i_];
                        coeff = a[i_];
                        sprintf(LSGRG_MSGBUFFER, "   ROW = \t%5ld COL = \t%5ld COEFF = \t%14.7g\tJTYP = %5ld\n",
                         irow, j, coeff, jtype[i_] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
L_21:
                ;
        }
        /*----------------------------------------------------- */
       lsgrg_msg(IOINFO,
        "\n\n==============================================================================\n" );
#endif
        return;
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS prscal(LsgrgInfo *_info,long int istat[], long int iheg[], long int ihegl[],
         long int ihag[], double grad[], double ascale[], long int ivstat[],
         long int inlin[], double g[], double gplus[], double gminus[],
         double gcol[], double x[], double alb[], double ub[], double albo[],
         double ubo[], long int iobjpt[], double paij[], long int iprow[],
         long int ipcol[], long int ipmap[])
{
LOGICAL32 jreset;
long int i, i1, i2, i_, ir, j, j_;
double bplus;



        /***********************************************************************
         **                       SUBROUTINE PRSCAL                            *
         **     Scales the problem data, including variables, funtions, bounds,*
         **  slacks, and the Jacobian.  M2SCAL is called to calculate the scale*
         **  factors.                                                          *
         ***********************************************************************
         *
         * */


        /*.......................................................................
         *
         *     **** START OF PRSCAL ****
         * */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 3 )
                {
                sprintf(LSGRG_MSGBUFFER, "\n\n*** PRSCAL ENTERED --- SCALED = %3c\n",
                 TorF(_info->scal.scaled) );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( _info->scal.scaled && _info->lincnt.nnlin == 0 )
                return;
        bplus = 0.10e0*_info->tols.plinfy;

        /*      ----------------------------------------------------------
         *     | If problem was previously scaled, unscale the original
         *     | bounds, the variables, the slacks, and the Jacobian.
         *     | Note that only the linear Jacobian elements are unscaled if
         *     | the Jacobian is UPDATED (COMPGR=TRUE) since the nonlinear
         *     | elements will be changed by GETGR.
         *      ----------------------------------------------------------
         * */

        if( _info->scal.scaled ){
#ifdef IO_ENABLED
               lsgrg_msg(IOINFO,
                 " *** UNSCALING PREVIOUSLY SCALED DATA....\n" );
#endif
                for( j = 1; j <= _info->dimen.n; j++ ){
                        j_ = j - 1;
                        if( alb[j_] > -bplus )
                                alb[j_] = alb[j_]*ascale[j_];
                        if( ub[j_] < bplus )
                                ub[j_] = ub[j_]*ascale[j_];
                        if( albo[j_] > -bplus )
                                albo[j_] = albo[j_]*ascale[j_];
                        if( ubo[j_] < bplus )
                                ubo[j_] = ubo[j_]*ascale[j_];
                        x[j_] = x[j_]*ascale[j_];
                        if( _info->cbmode.compgr ){
                                i1 = ihegl[j_];
                        } else{
                                i1 = iheg[j_];
                        }
                        i2 = iheg[j_ + 1] - 1;
                        if( i2 < i1 )
                                goto L_50;
                        for( i = i1; i <= i2; i++ ){
                                i_ = i - 1;
                                ir = ihag[i_];
                                grad[i_] = grad[i_]/(ascale[_info->dimen.n + ir - 1]*ascale[j_]);
                        }
L_50:
                        ;
                }
                for( j = _info->dimen.n + 1; j <= _info->dimen.npmp1; j++ ){
                        j_ = j - 1;
                        g[j_ - _info->dimen.n] = g[j_ - _info->dimen.n]/ascale[j_];
                        x[j_] = x[j_]/ascale[j_];
                        if( alb[j_] > -bplus )
                                alb[j_] = alb[j_]/ascale[j_];
                        if( ub[j_] < bplus )
                                ub[j_] = ub[j_]/ascale[j_];
                        if( albo[j_] > -bplus )
                                albo[j_] = albo[j_]/ascale[j_];
                        if( ubo[j_] < bplus )
                                ubo[j_] = ubo[j_]/ascale[j_];
                }
                _info->bestbk.truobj = _info->bestbk.truobj/ascale[_info->dimen.n + _info->nintbk.nobj - 1];
                _info->scal.scaled = FALSE;
        }

        /*      ---------------------------------------
         *     |  Update the Jacobian if necessary.
         *      ---------------------------------------
         * */
        if( _info->cbmode.compgr ){
                jreset = FALSE;
                getgr(_info, x, grad, ihag, iheg, ihegl, inlin, g, gplus, gminus,
                 gcol, ub, alb, ascale, &jreset, iobjpt, paij, iprow, ipcol,
                 ipmap );
        }

        /*      ----------------------------------------
         *     |  Compute the new scale factors
         *      ----------------------------------------
         * */
        m2scal(_info, _info->dimen.mp1, _info->dimen.n, istat, iheg, ihegl, ihag, grad, ascale,
         ivstat, gcol, gplus );

        /*      -------------------------------------------------------
         *     | Scale the Jacobian elements and find largest element
         *      -------------------------------------------------------
         * */
        scljac(_info, grad, ihag, iheg, ihegl, inlin, ascale, TRUE );

        chkelt(_info, grad, ihag, iheg, ihegl, inlin, TRUE );
        if( _info->glberr.abort )
                return;

        /*      -----------------------------------------------------------
         *     | Scale bounds, slacks, and the variables.  If bounds have
         *     | been perturbed, then update as needed.
         *      ----------------------------------------------------------- */
        for( j = 1; j <= _info->dimen.n; j++ ){
                j_ = j - 1;
                if( alb[j_] > -bplus )
                        alb[j_] = alb[j_]/ascale[j_];
                if( ub[j_] < bplus )
                        ub[j_] = ub[j_]/ascale[j_];
                if( albo[j_] > -bplus )
                        albo[j_] = albo[j_]/ascale[j_];
                if( ubo[j_] < bplus )
                        ubo[j_] = ubo[j_]/ascale[j_];
                x[j_] = x[j_]/ascale[j_];
        }
        for( j = _info->dimen.n + 1; j <= _info->dimen.npmp1; j++ ){
                j_ = j - 1;
                if( alb[j_] > -bplus )
                        alb[j_] = alb[j_]*ascale[j_];
                if( ub[j_] < bplus )
                        ub[j_] = ub[j_]*ascale[j_];
                if( albo[j_] > -bplus )
                        albo[j_] = albo[j_]*ascale[j_];
                if( ubo[j_] < bplus )
                        ubo[j_] = ubo[j_]*ascale[j_];
                x[j_] = x[j_]*ascale[j_];
                g[j_ - _info->dimen.n] = g[j_ - _info->dimen.n]*ascale[j_];
        }
        /*c    set TRUOBJ to newly scaled objective */
        _info->bestbk.truobj = _info->bestbk.truobj*ascale[_info->dimen.n + _info->nintbk.nobj - 1];
        /*XX   TRUOBJ = G(NOBJ) */
        _info->scal.scaled = TRUE;
        _info->scal.newscl = FALSE;
        return;



        /*     END OF PRSCAL */
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS scljac(LsgrgInfo *_info,double grad[], long int ihag[], long int iheg[],
         long int ihegl[], long int inlin[], double ascale[], LOGICAL32 sclall)
{
long int i, i1, i2, i_, ir, j, j_, jj, jj_;
/* double bplus; */



        /***********************************************************************
         **                       SUBROUTINE SCLJAC                            *
         **     Scales the elements in the Jacobian matrix.  All elements are  *
         **  scaled if SCLALL is TRUE, and only nonlinear elements if FALSE.   *
         ***********************************************************************
         *
         * */


        /*.......................................................................
         *
         *     **** START OF SCLJAC ****
         *
         *
         *      -----------------------------
         *     |  Scale the linear elements
         *      ----------------------------- */
        if( sclall ){
                for( j = 1; j <= _info->dimen.n; j++ ){
                        j_ = j - 1;
                        i1 = ihegl[j_];
                        i2 = iheg[j_ + 1] - 1;
                        for( i = i1; i <= i2; i++ ){
                                i_ = i - 1;
                                ir = ihag[i_];
                                grad[i_] = grad[i_]*(ascale[_info->dimen.n + ir - 1]*ascale[j_]);
                        }
                }
        }

        /*      ---------------------------------------
         *     |  Now, Scale the nonlinear elements
         *      --------------------------------------- */
        for( jj = 1; jj <= _info->lincnt.nnlin; jj++ ){
                jj_ = jj - 1;
                j = inlin[jj_];
                i1 = iheg[j - 1];
                i2 = ihegl[j - 1] - 1;
                for( i = i1; i <= i2; i++ ){
                        i_ = i - 1;
                        ir = ihag[i_];
                        grad[i_] = grad[i_]*(ascale[_info->dimen.n + ir - 1]*ascale[j - 1]);
                }
        }
        return;
        /*     *** End of SCLJAC *** */
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS m2scal(LsgrgInfo *_info,long int m, long int n, long int istat[],
         long int he[], long int hel[], long int ha[], double a[], double ascale[],
         long int ivstat[], double rmin[], double rmax[])
{
long int i, i1, i2, i_, ilin, ir, j, j_, jlin, k, k_, maxk, mn, npass;
double ac, amax, amin, ar, aratio, bplus, cmax, cmin, cratio, damp,
         one, scltol, sratio, zero;



        /***********************************************************************
         **  M2SCAL  -  Calculates row and column scale factors                *
         ***********************************************************************
         * */

        /*.......................................................................
         *
         *
         *
         *     M2SCAL computes scale factors ASCALE from A.
         *
         *     In Phase 1, an iterative method based on geometric means is
         *     used to compute scales from A alone. This procedure is derived
         *     from a routine written by Robert Fourer, 1979. The main steps
         *     are:
         *
         *        (1) COMPUTE ARATIO = MAX(I1,I2,J)  A(I1,J) / A(I2,J)
         *        (2) DIVIDE EACH ROW I BY
         *               ( MIN(J) A(I,J) * MAX(J) A(I,J) ) ** 1/2
         *        (3) DIVIDE EACH COLUMN J BY
         *               ( MIN(I) A(I,J) * MAX(I) A(I,J) ) ** 1/2
         *        (4) COMPUTE SRATIO AS IN (1)
         *        (5) IF SRATIO .LT. SCLTOL * ARATIO, REPEAT FROM STEP (1)
         *
         *     Free rows (ISTAT=0) and fixed columns (BL=BU) are not used
         *     in defining scales for other rows and columns.  In Phase 2,
         *     free rows are scaled by their largest element (so their
         *     PI values will not affect PINORM).  The OBJ ROW is then
         *     reset to have SCALE = 1.0.  During the subsequent optimization
         *     its size is accounted for by PINORM.
         *
         *     INITIAL VERSION  MARCH 1981.
         *     FEB 1983.  For convenience, everything is now DOUBLE-PRECISION.
         *     AUG 1983.  Revised to damp the effect of very small elements.
         *                On each pass, a new row or column scale will not be
         *                smaller than  DSQRT(DAMP)  Times the largest (SCALED)
         *                element in that row or column.
         *     ------------------------------------------------------------------
         *
         *
         *     **** START OF M2SCAL ****
         * */
        mn = m + n;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 3 )
                {
               lsgrg_msg(IOINFO,
                "\n\n SCALING\n -------\n             MIN ELEM    MAX ELEM       MAX COL RATIO\n\n" );
                }
#endif
        zero = 0.0e0;
        bplus = 0.1*_info->tols.plinfy;
        one = 1.0e0;
        damp = 1.0e-4;
        scltol = 0.99e0;
        /*      ---------------------------------
         *     | Initialize scale factors to 1
         *      ---------------------------------
         * */
        for( j = 1; j <= mn; j++ ){
                j_ = j - 1;
                ascale[j_] = one;
        }

        aratio = bplus;
        ilin = 1;
        jlin = 1;


        /*      ****************
         *          MAIN LOOP
         *      **************** */
        maxk = 11;

        for( k = 1; k <= maxk; k++ ){
                k_ = k - 1;

                /*        Find the largest column ratio.
                 *        Also set new column scales (Except on PASS 0).
                 * */
                npass = k - 1;
                amin = bplus;
                amax = zero;
                sratio = one;

                for( j = jlin; j <= n; j++ ){
                        j_ = j - 1;
                        if( ivstat[j_] < 0 )
                                goto L_250;
                        i1 = he[j_];
                        i2 = he[j_ + 1] - 1;
                        cmin = bplus;
                        cmax = zero;

                        /*XXXX       IF (I2 .LT. I1) GO TO 232 */
                        for( i = i1; i <= i2; i++ ){
                                i_ = i - 1;
                                ir = ha[i_];
                                /*-I            IR     = IABS(IR) */
                                if( istat[ir - 1] == 0 )
                                        goto L_230;
                                ar = fabs( a[i_] );
                                if( ar <= _info->ztol.xajtol )
                                        goto L_230;
                                ar = ar/ascale[n + ir - 1];
                                cmin = fmin( cmin, ar );
                                cmax = fmax( cmax, ar );
L_230:
                                ;
                        }

/* L_232: unreferenced???*/
                        ac = fmax( cmin, damp*cmax );
                        ac = sqrt( ac )*sqrt( cmax );
                        if( cmax <= _info->ztol.xajtol )
                                ac = one;
                        if( npass > 0 )
                                ascale[j_] = ac;
                        amin = fmin( amin, cmin/ascale[j_] );
                        amax = fmax( amax, cmax/ascale[j_] );
                        cratio = cmax/cmin;
                        sratio = fmax( sratio, cratio );
L_250:
                        ;
                }

#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 3 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " AFTER%3ld%12.2e%12.2e%20.2f\n",
                         npass, amin, amax, sratio );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( sratio >= aratio*scltol )
                        goto L_600;
                if( k == maxk )
                        goto L_600;
                aratio = sratio;


                /*        Set new row scales for next pass.
                 * */
                if( ilin > m )
                        goto L_500;

                for( i = ilin; i <= m; i++ ){
                        i_ = i - 1;
                        rmin[i_] = bplus;
                        rmax[i_] = zero;
                }

                for( j = 1; j <= n; j++ ){
                        j_ = j - 1;
                        if( ivstat[j_] < 0 )
                                goto L_350;
                        i1 = he[j_];
                        i2 = he[j_ + 1] - 1;
                        /*XXXX       IF (I2 .LT. I1) GO TO 350 */
                        for( i = i1; i <= i2; i++ ){
                                i_ = i - 1;
                                ir = ha[i_];
                                /*-I            IR     = IABS(IR) */
                                if( ir < ilin )
                                        goto L_330;
                                ar = fabs( a[i_] );
                                if( ar <= _info->ztol.xajtol )
                                        goto L_330;
                                ar = ar/ascale[j_];
                                rmin[ir - 1] = fmin( rmin[ir - 1], ar );
                                rmax[ir - 1] = fmax( rmax[ir - 1], ar );
L_330:
                                ;
                        }
L_350:
                        ;
                }

                for( i = ilin; i <= m; i++ ){
                        i_ = i - 1;
                        ar = fmax( rmin[i_], damp*rmax[i_] );
                        ascale[n + i_] = sqrt( ar )*sqrt( rmax[i_] );
                }
L_500:
                ;
        }

        /*     ********************
         *       END OF MAIN LOOP
         *     ********************
         *
         *     Deal with empty rows and free rows.
         *     Invert the row and column scales so that only multiplications
         *     need be performed when scaling the Jacobian Matrix.
         * */
L_600:
        if( ilin > m )
                goto L_650;

        for( i = ilin; i <= m; i++ ){
                i_ = i - 1;
                ar = ascale[n + i_];
                if( istat[i_] == 0 )
                        ar = rmax[i_];
                if( ar <= _info->ztol.xajtol )
                        ar = one;
                ascale[n + i_] = one/ar;
        }

        /*X650 IF (NOBJ .NE. 0) ASCALE(N+NOBJ) = ONE */
L_650:
        ;
        _info->sclobj.objscl = ascale[n + _info->nintbk.nobj - 1];

        /*     Invert the column scales
         * */
        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                ascale[i_] = one/ascale[i_];
        }
#ifdef IO_ENABLED
        if( ilin <= m && _info->nintbk.ipr >= 3 ){
               lsgrg_msg(IOINFO,
                 "\n\n\n ROW SCALES  R(I)         (SCALED AIJ)  =  (ORIGINAL AIJ)  *  R(I)  *  C(J)\n ----------------\n\n" );
                for( i = ilin; i <= m; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%5ld) %16.8e", i, ascale[n + i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
               lsgrg_msg(IOINFO,
                 "\n\n\n COLUMN SCALES  C(J)\n -------------------\n\n" );
                for( j = jlin; j <= n; j++ ){
                        sprintf(LSGRG_MSGBUFFER, "%5ld) %16.8e", j, ascale[j - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
        }
#endif
        return;


        /*     End of M2SCAL */
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS condnm(LsgrgInfo *_info,long int n, short int *ind, double a[], long int nzdim,
         short int *lnrc, long int *ip, long int iluprm[], double w[],
         long int irank[], double cnum[], double grad[], long int ihag[],
         long int iheg[], long int ibv[], double w1[], double w2[])
{
#define IND(I_,J_)      (*(ind+(I_)*(nzdim)+(J_)))
#define LNRC(I_,J_)     (*(lnrc+(I_)*(n)+(J_)))
#define IP(I_,J_)       (*(ip+(I_)*(n)+(J_)))
LOGICAL32 len1;
long int collen, error, i, i_, ictr, ictr_, iend, ifirst, ii, ii_,
         ilast, /* ipcd,*/  irow, is, is_, istart, j, j2, j3, jcol, jj, jj_,
         jpiv, num, nz;
double abspiv, ainvn, anorm, colnrm, sm, sp, wii, wsm1, xnorm, ynorm,
         zm, zp;


        /*......................................................................
         * CONDNM uses the factors produced by the invert routines  to
         *     estimate the condition number of each diagonal block.
         *
         * WE NOW DESCRIBE THE ARGUMENT LIST FOR CONDNM
         * N  IS AN INTEGER VARIABLE SET TO THE ORDER OF THE MATRIX. IT IS NOT
         *     ALTERED BY THE SUBROUTINE.
         * ICN IS AN INTEGER*2  ARRAY OF LENGTH LICN. ENTRIES IDISP(1) TO
         *     IDISP(2) SHOULD BE UNCHANGED SINCE THE LAST CALL TO MA30A/AD. IF
         *     THE MATRIX HAS MORE THAN ONE DIAGONAL BLOCK, THEN COLUMN INDICES
         *     CORRESPONDING TO NON-ZEROS IN SUB-DIAGONAL BLOCKS OF PAQ MUST
         *     APPEAR IN POSITIONS 1 TO IDISP(1)-1. FOR THE SAME ROW THOSE
         *     ENTRIES MUST BE CONTIGUOUS, WITH THOSE IN ROW I PRECEDING THOSE
         *     IN ROW I+1 (I=1,...,N-1) AND NO WASTED SPACE BETWEEN ROWS.
         *     ENTRIES MAY BE IN ANY ORDER WITHIN EACH ROW. IT IS NOT ALTERED
         *     BY CONDNM.
         * A  IS A REAL/DOUBLE PRECISION ARRAY OF LENGTH LICN.  ENTRIES
         *     IDISP(1) TO IDISP(2) SHOULD BE UNCHANGED SINCE THE LAST CALL TO
         *     MA30A/AD OR MA30B/BD.  IF THE MATRIX HAS MORE THAN ONE DIAGONAL
         *     BLOCK, THEN THE VALUES OF THE NON-ZEROS IN SUB-DIAGONAL BLOCKS
         *     MUST BE IN POSITIONS 1 TO IDISP(1)-1 IN THE ORDER GIVEN BY ICN.
         *     IT IS NOT ALTERED BY MA30C/CD.
         * LICN  IS AN INTEGER VARIABLE SET TO THE SIZE OF ARRAYS ICN AND A.
         *     IT IS NOT ALTERED BY CONDNM.
         * LENR,LENRL ARE INTEGER*2  ARRAYS OF LENGTH N WHICH SHOULD BE
         *     UNCHANGED SINCE THE LAST CALL TO MA30A/AD. THEY ARE NOT ALTERED
         *     BY MA30C/CD.
         * IP,IQ ARE INTEGER*2 ARRAYS OF LENGTH N WHICH SHOULD BE UNCHANGED
         *     SINCE THE LAST CALL TO MA30A/AD. THEY ARE NOT ALTERED BY
         *     MA30C/CD.
         * W  IS A REAL/DOUBLE PRECISION ARRAY OF LENGTH N WHICH IS USED AS
         *     WORKSPACE BY MA30C/CD.
         *
         *......................................................................
         *
         * */

        /*...................................................................
         *
         * */
        /*     *** START OF CONDNM *** */




#ifdef IO_ENABLED
        if( _info->dbug.debug )
                {
               lsgrg_msg(IOINFO,
                " ENTERING SUBROUTINE COND....\n" );
                }
#endif

        /*     IFIRST holds the index of the first row in the current block.
         * */
        ifirst = 1;
        ilast = n;

        /*     NBLCK counts the number of diagonal blocks
         * */
        _info->zblck.nblck = 0;

        /*     NOW SOLVE UT Z = B
         * */
        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                w[i_] = 0.0e0;
        }

        /*     Start of WHILE  - Loop for each block
         * */
L_710:
        ;
        if( ifirst > ilast )
                goto L_500;

        _info->zblck.nblck = _info->zblck.nblck + 1;
/* L_730: unreferenced??? */
        irank[_info->zblck.nblck - 1] = ilast - ifirst + 1;

        /*        J1 points to the position of the beginning of row I (LT part)
         *        JPIV points to the pivot in current row.
         *
         *C******************************************************
         *      WRITE(IOOUT,6767) NBLCK, IRANK(NBLCK)
         *6767  FORMAT( // '  THE LU FACTORS FOR BLOCK # ',I3,'  RANK = ',I3,/)
         *      WRITE(IOOUT,6768)(IP(I),I=IFIRST,ILAST)
         *6768  FORMAT('  IP: ',(T8,15I4))
         *      WRITE(IOOUT,6769)(IQ(I),I=IFIRST,ILAST)
         *6769  FORMAT('  IQ: ',(T8,15I4))
         *      J1 = IBLOCK
         *      DO 7878 II=IFIRST,ILAST
         *      JPIV = J1 + LENRL(II)
         *      J2 = JPIV - 1
         *      J3 = J1 + LENR(II) - 1
         *      WRITE(IOOUT,6780) II,IP(II),LENR(II),LENRL(II)
         *6780  FORMAT(' *** ROW: ',I3, ' ACTUAL ROW:',I3,' LENGTH:',I3,
         *     1       ' LENGTH OF L:', I3)
         *      WRITE(IOOUT,6781) (A(I),I=J1,J2)
         *6781  FORMAT(' ROW OF L: ',(T10,5E15.8))
         *      WRITE(IOOUT,6782) (ICN(I),I=J1,J2)
         *6782  FORMAT(' COL IND. :',(T10,5I15))
         *      WRITE(IOOUT,6783) (A(I),I=JPIV,J3)
         *6783  FORMAT(' ROW OF U: ',(T10,5E15.8))
         *      WRITE(IOOUT,6784) (ICN(I),I=JPIV,J3)
         *6784  FORMAT(' COL IND. :',(T10,5I15))
         *      J1 = J3+1
         *7878  CONTINUE
         *        Do First row separately
         * */
        irow = LNRC(3,ifirst - 1);
        jcol = LNRC(2,ifirst - 1);
        jpiv = IP(0,jcol - 1);
        wsm1 = 1.0e0/a[jpiv - 1];
        w[irow - 1] = wsm1;
        if( ifirst == ilast )
                goto L_610;
        len1 = LNRC(0,jcol - 1) == 1;
        if( !len1 ){
                j2 = jpiv + 1;
                j3 = LNRC(0,jcol - 1) + jpiv - 1;
        }

        /*        For each row, Compute W(I) and eliminate
         * */
        for( is = ifirst + 1; is <= ilast; is++ ){
                is_ = is - 1;

                /*        Update the PJ (S-1)
                 * */
                if( len1 )
                        goto L_750;
                for( jj = j2; jj <= j3; jj++ ){
                        jj_ = jj - 1;
                        j = IND(1,jj_);
                        w[j - 1] = w[j - 1] + a[jj_]*wsm1;
                }

                /*           Do comparison to find W(IROW)
                 * */
L_750:
                irow = LNRC(3,is_);
                jcol = LNRC(2,is_);
                jpiv = IP(0,jcol - 1);

                zp = 1.0e0 - w[irow - 1];
                zm = -1.0e0 - w[irow - 1];
                sp = fabs( zp );
                sm = fabs( zm );
                zp = zp/a[jpiv - 1];
                zm = zm/a[jpiv - 1];

                for( i = is + 1; i <= ilast; i++ ){
                        i_ = i - 1;
                        ii = LNRC(3,i_);
                        w1[ii - 1] = w[ii - 1];
                        w2[ii - 1] = w[ii - 1];
                }

                len1 = LNRC(0,jcol - 1) == 1;
                if( len1 )
                        goto L_770;
                j2 = jpiv + 1;
                j3 = LNRC(0,jcol - 1) + jpiv - 1;
                for( jj = j2; jj <= j3; jj++ ){
                        jj_ = jj - 1;
                        j = IND(1,jj_);
                        w1[j - 1] = w1[j - 1] + a[jj_]*zp;
                        w2[j - 1] = w2[j - 1] + a[jj_]*zm;
                }

L_770:
                ;
                for( i = is + 1; i <= ilast; i++ ){
                        i_ = i - 1;
                        ii = LNRC(3,i_);
                        sp = sp + fabs( w1[ii - 1] );
                        sm = sm + fabs( w2[ii - 1] );
                }

                if( sp >= sm ){
                        wsm1 = zp;
                } else{
                        wsm1 = zm;
                }

                w[irow - 1] = wsm1;
        }


        /*        Back Substitution.
         *        This loop does the back substitution on the rows of the block I
         *        The reverse order doing it simultaneously on the L TRANSPOSE
         *        part of the Diagonal Blocks and the Off-Diagonal Blocks.
         *
         *        Set W() = W() * L INVERSE
         * */
        istart = nzdim - iluprm[5] + 1;
        if( istart > nzdim )
                goto L_300;
        for( ii = istart; ii <= nzdim; ii++ ){
                ii_ = ii - 1;
                wii = w[IND(1,ii_) - 1];
                if( wii == 0.0e0 )
                        goto L_270;
                i = IND(0,ii_);
                w[i - 1] = w[i - 1] + wii*a[ii_];
L_270:
                ;
        }


        /*        Now compute the 1-NORM of the vector W (X)
         * */
L_300:
        ;
        xnorm = 0.0e0;
        for( i = ifirst; i <= ilast; i++ ){
                i_ = i - 1;
                xnorm = xnorm + fabs( w[i_] );
        }


        /*        We now solve   A * W1 = W.
         * */
        nz = 0;
        lsluft( n, nzdim, a, ind, ip, lnrc, w, w1, &nz, &ii, iluprm );

        /*        Find the 1-NORM of W1 (Y) and the estimate of NORM OF AINV
         * */
        ynorm = 0.0e0;
        for( i = ifirst; i <= ilast; i++ ){
                i_ = i - 1;
                ynorm = ynorm + fabs( w1[i_] );
        }

        ainvn = ynorm/xnorm;

        /*        Now find the 1-NORM of this block
         * */
        anorm = 0.0e0;
        for( i = ifirst; i <= ilast; i++ ){
                i_ = i - 1;
                jcol = ibv[i_];
                if( jcol >= _info->memory.np1 ){
                        colnrm = 1.0e0;
                } else{
                        istart = iheg[jcol - 1];
                        iend = iheg[jcol] - 1;
                        collen = iend - istart + 1;
                        if( collen <= 0 || collen > _info->cmax.colmax )
                                goto L_430;
                        colnrm = 0.0e0;
                        for( ictr = istart; ictr <= iend; ictr++ ){
                                ictr_ = ictr - 1;
                                colnrm = colnrm + fabs( grad[ictr_] );
                        }
                }
                if( colnrm > anorm )
                        anorm = colnrm;
        }

        /*        Now compute the Condition Number estimate
         * */
        cnum[_info->zblck.nblck - 1] = anorm*ainvn;
#ifdef IO_ENABLED
        if( _info->dbug.debug ){
                sprintf(LSGRG_MSGBUFFER, " COND NUMBER RESULTS FOR BLOCK # %5ld  RANK = %5ld\n",
                 _info->zblck.nblck, irank[_info->zblck.nblck - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                sprintf(LSGRG_MSGBUFFER, " XNORM = %15.8e  YNORM = %15.8e\n AINVN = %15.8e  ANORM = %15.8e  COND NUM = %15.8e\n",
                 xnorm, ynorm, ainvn, anorm, cnum[_info->zblck.nblck - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        goto L_620;

        /*        Come here if block is Singular
         * */
/* L_600: */
        cnum[_info->zblck.nblck - 1] = _info->tols.plinfy;
#ifdef IO_ENABLED
        if( _info->dbug.debug )
                {
                sprintf(LSGRG_MSGBUFFER, " BLOCK NUMBER %5ld IS SINGULAR\n",
                 _info->zblck.nblck );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        goto L_620;

        /*        Come here for a 1 X 1 block
         * */
L_610:
        abspiv = fabs( a[jpiv - 1] );
        if( abspiv >= _info->limits.epspiv ){
                cnum[_info->zblck.nblck - 1] = 1.0e0;
#ifdef IO_ENABLED
                if( _info->dbug.debug )
                        {
                        sprintf(LSGRG_MSGBUFFER, " BLOCK NUMBER %5ld IS A 1 X 1 BLOCK PIVOT = %15.8e\n",
                         _info->zblck.nblck, a[jpiv - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
        } else{
#ifdef IO_ENABLED
                if( _info->dbug.debug )
                        {
                        sprintf(LSGRG_MSGBUFFER,
                         " WARNING... THE 1 X 1 BLOCK # %5ld HAD A PIVOT"
                         " LESS THAN EPSPIV  -  PIVOT = %14.8e\n",
                         _info->zblck.nblck, a[jpiv - 1] );
                   lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                /*CXX        CNUM(NBLCK) = PLINFY */
                cnum[_info->zblck.nblck - 1] = 1.0e0;
        }
        anorm = abspiv;
        ainvn = 1.0e0/abspiv;



        /*        Now reset IBLOCK and IFIRST for the next block
         * */
L_620:
        ifirst = ilast + 1;
        goto L_710;
        /*     END  OF WHILE LOOP AT 710
         *
         *     Check to make sure the number of blocks is correct
         * */
L_500:
        ;
        num = 1;
        if( _info->zblck.nblck != num ){
#ifdef IO_ENABLED
                sprintf( LSGRG_MSGBUFFER,
                 " COND....ERROR: THE NUMBER OF BLOCKS FOUND IS DIFFERENT"
                 " THAN IN MA28\n BLOCKS FOUND = %5ld EXPECTED = %5ld\n",
                 _info->zblck.nblck, num );
                 lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                error = 0;
    /*           xerror( error, iounit.ioerr );    */
        lsgrg_errorexit(JUMPBUF,_LSGRG_CONDNM_BAD_NBLOCKS);
                return;
        }

#ifdef IO_ENABLED
        if( _info->dbug.debug ){
                sprintf(LSGRG_MSGBUFFER, " NUMBER OF DIAGONAL BLOCKS IS %6ld\n",
                 _info->zblck.nblck );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
               lsgrg_msg(IOINFO,
                " IRANK IS ....\n" );
                for( i = 1; i <= _info->zblck.nblck; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%5ld", irank[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
               lsgrg_msg(IOINFO,
                 " COND NUMBERS ARE ....\n" );
                for( i = 1; i <= _info->zblck.nblck; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", cnum[i - 1] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,"\n" );
               lsgrg_msg(IOINFO,
                " EXITING SUBROUTINE COND....\n" );
        }
#endif

        /*        Find the largest Condition Number */

        _info->zblck.condmx = 0.0e0;
        for( i = 1; i <= _info->zblck.nblck; i++ ){
                i_ = i - 1;
                zp = cnum[i_];
                if( zp > _info->zblck.condmx )
                        _info->zblck.condmx = zp;
#ifdef IO_ENABLED
                if( zp > _info->zcond.cndtol && _info->nintbk.ipr >= 4 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " DIAGONAL BLOCK # %5ld IS ILL-COND -- COND NUM = %15.8e\n",
                         i, zp );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
        }





        return;

L_430:
#ifdef IO_ENABLED
        sprintf( LSGRG_MSGBUFFER,
         " COND.....ERROR: THE BASIC COLUMN IN POSITION%6ld HAS"
         " LENGTH%6ld      (ACTUAL VARIABLE NUMBER = %6ld)\n",
         i, collen, ibv[i - 1] );
        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
        error = 0;
    /*   xerror( error, iounit.ioerr );            */
        lsgrg_errorexit(JUMPBUF,_LSGRG_CONDNM_BAD_COLLEN);
        return;


        /* ::::::::::::::::::::::::::::::::::: END  OF CONDNM    ::::::::::
         *
         * */
#undef  LNRC
#undef  IND
#undef  IP
} /* end of function */



void /*FUNCTION*/ LSGRGCLASS chkelt(LsgrgInfo *_info,double grad[], long int ihag[], long int iheg[],
         long int ihegl[], long int inlin[], LOGICAL32 chkall)
{
long int i, i1, i2, i_, icolbig, irowbig, isvl, isvnl, j, j_, jj,
         jj_, jsvl, jsvnl;



        /***********************************************************************
         **                       SUBROUTINE CHKELT                            *
         **     Computes the largest Jacobian element and then checks it to    *
         **  insure that it is within a reasonable range.  If extremely large  *
         **  elements are found, LSGRG2 uses this as an indication that errors *
         **  are present in the Jacobian computation and ABORTS
         ***********************************************************************
         *
         * */

        /*....................................................................... */

        /*     **** START OF CHKELT ****
         *
         *
         *      ------------------------------------
         *     |  Find largest  Nonlinear Element
         *      ------------------------------------ */
        _info->ztol.bigelt = 0.0e0;
        for( jj = 1; jj <= _info->lincnt.nnlin; jj++ ){
                jj_ = jj - 1;
                j = inlin[jj_];
                i1 = iheg[j - 1];
                i2 = ihegl[j - 1] - 1;
                for( i = i1; i <= i2; i++ ){
                        i_ = i - 1;
                        if( fabs( grad[i_] ) > _info->ztol.bigelt ){
                                _info->ztol.bigelt = fabs( grad[i_] );
                                isvnl = ihag[i_];
                                jsvnl = j;
                        }
                }
        }

        /*      -------------------------------------------
         *     |  Find largest Linear element if desired
         *      ------------------------------------------- */
        if( chkall ){
                _info->ztol.biglin = 0.0e0;
                for( j = 1; j <= _info->dimen.n; j++ ){
                        j_ = j - 1;
                        i1 = ihegl[j_];
                        i2 = iheg[j_ + 1] - 1;
                        for( i = i1; i <= i2; i++ ){
                                i_ = i - 1;
                                if( fabs( grad[i_] ) > _info->ztol.biglin ){
                                        _info->ztol.biglin = fabs( grad[i_] );
                                        isvl = ihag[i_];
                                        jsvl = j;
                                }
                        }
                }
        }
        irowbig = isvnl;
        icolbig = jsvnl;
        if( _info->ztol.biglin > _info->ztol.bigelt ){
                irowbig = isvl;
                icolbig = jsvl;
                _info->ztol.bigelt = _info->ztol.biglin;
        }
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 4 )
                {
                sprintf(LSGRG_MSGBUFFER, " THE LARGEST JACOBIAN ELEMENT IS %15.8e\n",
                 _info->ztol.bigelt );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }

        if( fabs( _info->ztol.bigelt ) < _info->limits.epspiv ){
                if( _info->nintbk.ipr >= 1 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " CHKELT... ERROR: MAX JACOBIAN ELEMENT = %16.8e IS LESS THAN %12.6e CHECK DERIVATIVE COMPUTATION\n",
                         _info->ztol.bigelt, _info->limits.epspiv );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                /*XXX     CALL XERROR(0,IOERR) */
        } else if( fabs( _info->ztol.bigelt ) > 1.0e12 ){
                if( _info->nintbk.ipr >= 1 ){
                        sprintf(LSGRG_MSGBUFFER, " CHKELT... ERROR: MAX JACOBIAN ELEMENT = %16.8e IS LARGER THAN %12.6e CHECK DERIVATIVE COMPUTATION\n",
                         _info->ztol.bigelt, 1.0e12 );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        sprintf(LSGRG_MSGBUFFER, " CHEKELT...MAX ELT IS IN ROW %6ld AND IN COLUMN %6ld\n",
                         irowbig, icolbig );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                /*XXX     CALL XERROR(0,IOERR) */
        }
#endif
        return;


        /*     *** END OF CHKELT *** */
} /* end of function */
#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"

#define IOINFO  &(_info->io_info)
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)
/***********************************************************************
 **                 FILE:        GRGDIREC FORTRAN                      *
 **                 AUTHOR:      STUART SMITH                          *
 **                 LAST UPDATE: 27 FEB 1992                           *
 **                                                                    *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **      This file contains the following routines for LSGRG2, A       *
 **  large scale implementation of GRG2:                               *
 **                                                                    *
 **                                                                    *
 **
 **  ROUTINE             DESCRIPTION                      LAST UPDATE  *
 ** ---------           -------------                    ------------- *
 **  ADDCOL     Adds a unit column to approx. Hessian     04 DEC 1989  *
 **  ADDCG      Adds a Superbasic for CG method           12 MAR 1990  *
 **  CG         Conjugate Gradient search direction       27 FEB 1992  *
 **  COMDFP     Quasi-Newton search direction             25 APR 1989  *
 **  DELCOL     Deletes a column from approx. Hessian     27 FEB 1992  *
 **  DIREC      Search direction driving routine          27 FEB 1992  *
 ***********************************************************************
 * */
void /*FUNCTION*/ LSGRGCLASS addcol(LsgrgInfo *_info,double r[], double y[], long int ny, double diag)
{
long int i, i_, k, nsprev;


        nsprev = _info->nintbk.nsuper - 1;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                y[i_] = 0.0e0;
        }

        /*      -------------------------
         *     | Insert new column of R
         *      ------------------------ */
        y[_info->nintbk.nsuper - 1] = diag;
        k = nsprev*_info->nintbk.nsuper/2;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                r[k + i_] = y[i_];
        }
        return;

        /*     End of ADDCOL
         * */
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS cg(LsgrgInfo *_info,double gradfp[], double gradf[], double d[],
         long int *msgcg, long int *update, double *y, double *s, double *work)
{
#define Y(I_,J_)        (*(y+(I_)*(_info->dimen.n)+(J_)))
#define S(I_,J_)        (*(s+(I_)*(_info->dimen.n)+(J_)))
#define WORK(I_,J_)     (*(work+(I_)*(_info->cgbk.mcgm1-(0)+1)+(J_)))
long int i, i_, ii, ii_, j, j_, k, k_, newest;
double alphai, betaj, cgbeta, g1, g1n, g2, g2n, gamma, gamtmp, gtd,
         gty, gyd, rho, sj, std, yj, ytd, yts, yty;

/* replace statics by alg */

   #define initcg _info->cgStatic.initcg
   #define itncg  _info->cgStatic.itncg

/*
   static long initcg = 1;
   static long itncg = 0;
*/

        /* ................................................
         *      -------------------------------------------
         *     | Conjugate Gradient method on Superbasics
         *      ------------------------------------------- */
        if( _info->nintbk.ipr >= 5 )
                {
                lsgrg_msg(IOINFO, " CG ENTERED \n" );
                }
        if( initcg != 0 )
                *msgcg = 1;
        initcg = 0;
        _info->hesblk.gamma0 = 1.0e0;
        _info->hesblk.dirsc = 1.0e0;
        if( *msgcg != 0 )
                {
                lsgrg_msg(IOINFO,
                 " HESSIAN IS TOO LARGE FOR VARIABLE METRIC--SWITCH TO"
                 " CONJUGATE GRADIENTS \n" );
                }
        if( (!_info->logblk.restrt) && (itncg <= _info->nintbk.nsuper) )
                goto L_20;

        /*      ----------
         *     | Restart
         *      ---------- */
L_9:
        ;
        if( _info->nintbk.ipr > 3 )
                {
                lsgrg_msg(IOINFO,
                   "  CONJUGATE GRADIENT DIRECTION RESET TO -GRADF \n" );
                }
        *update = 1;
        itncg = 0;
  /*    cgbeta = 0.0e0;  unused assignment */

        /*       --------------------------------------------
         *      | These are pointers for circular memory list
         *       --------------------------------------------- */
        _info->cgbk.icgcnt = 0;
        _info->cgbk.icgptr = 0;

        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                d[i_] = -gradf[i_];
        }
        if( *msgcg <= 0 )
                goto L_210;
        switch( _info->cgbk.modcg ){
                case 1: goto L_11;
                case 2: goto L_12;
                case 3: goto L_13;
                case 4: goto L_14;
                case 5: goto L_15;
                case 6: goto L_16;
                }
L_11:
        lsgrg_msg(IOINFO, "\n FLETCHER-REEVES DIRECTION WILL BE USED \n" );
        goto L_210;
L_12:
        lsgrg_msg(IOINFO, "\n POLAK-RIBIERE DIRECTION WILL BE USED \n" );
        goto L_210;
L_13:
        lsgrg_msg(IOINFO, "\n PERRY (MARCH 76) DIRECTION WILL BE USED \n" );
        goto L_210;
L_14:
        lsgrg_msg(IOINFO, "\n DFP DIRECTION WILL BE USED \n" );
        goto L_210;
L_15:
        lsgrg_msg(IOINFO, "\n COMPLEMENTARY DFP DIRECTION WILL BE USED \n" );
        goto L_210;
L_16:
        lsgrg_msg(IOINFO, "\n LIMITED MEMORY BFGS DIRECTION WILL BE USED  \n" );
        goto L_210;

L_20:
        *update = 3;
        switch( _info->cgbk.modcg ){
                case 1: goto L_30;
                case 2: goto L_60;
                case 3: goto L_90;
                case 4: goto L_110;
                case 5: goto L_120;
                case 6: goto L_400;
                }

        /*      -----------------------
         *     | Limited memory BFGS
         *      ----------------------- */
L_400:
        if( *msgcg > 0 )
             {
                lsgrg_msg(IOINFO, " LIMITED MEMORY BFGS DIRECTION WILL BE USED  \n" );
             }

        /*      ---------------------------------------------------------
         *     | Set correct index into Circular List:
         *     |    ICGPTR -- Position of first (oldest) data
         *     |    NEWEST -- Position where current data is to be stored
         *      -----------------------------------------------------------
         * */
        if( _info->cgbk.icgcnt < _info->cgbk.memcg ){
                newest = _info->cgbk.icgcnt;
                _info->cgbk.icgcnt = _info->cgbk.icgcnt + 1;
        } else{
                newest = _info->cgbk.icgptr;
                _info->cgbk.icgptr = (_info->cgbk.icgptr + 1)%_info->cgbk.memcg;
        }

        /*      ----------------------------------------
         *     | Calculate Y, S, YTY, and YTS
         *     | WORK(*,1) holds the coefficients RHO
         *      ---------------------------------------- */
        yts = 0.0e0;
        yty = 0.0e0;
        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;
                sj = _info->bestbk.step*d[j_];
                yj = gradf[j_] - gradfp[j_];
                Y(newest,j_) = yj;
                S(newest,j_) = sj;
                yts = yts + yj*sj;
                yty = yty + yj*yj;
        }
        if( fabs( yts ) <= _info->tols.eps ){
#ifdef IO_ENABLED
                sprintf(LSGRG_MSGBUFFER, " DIREC ERROR .... CG -divide by 0 -- yts = %13.6e\n",
                 yts );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                goto L_9;
        }

        WORK(0,newest) = 1.0e0/yts;

        /*      --------------------------------------------
         *     | Calculate H*G using the recursion formula
         *     | WORK(I,2) will hold ALPHA(I)
         *      -------------------------------------------- */
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                d[i_] = -gradf[i_];
        }

        for( i = _info->cgbk.icgcnt - 1; i >= 0; i-- ){
       /*       i_ = i - 1;  unreferenced??? */
                j = (_info->cgbk.icgptr + i)%_info->cgbk.memcg;
                std = 0.0e0;
                for( k = 1; k <= _info->nintbk.nsuper; k++ ){
                        k_ = k - 1;
                        std = std + S(j,k_)*d[k_];
                }
                alphai = WORK(0,j)*std;
                WORK(1,i) = alphai;
                for( ii = 1; ii <= _info->nintbk.nsuper; ii++ ){
                        ii_ = ii - 1;
                        d[ii_] = d[ii_] - alphai*Y(j,ii_);
                }
        }

        /*      -----------------------------
         *     | Automatic Hessian scaling
         *      ----------------------------- */
        if( _info->cgbk.hscale ){
                gamtmp = yts/yty;
                if( gamtmp >= 1.0e-7 && gamtmp <= 1.0e7 ){
                        _info->hesblk.gamma0 = gamtmp;
                        _info->hesblk.dirsc = _info->hesblk.gamma0;
                        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                                i_ = i - 1;
                                d[i_] = _info->hesblk.gamma0*d[i_];
                        }
                }
        }

        for( i = 0; i <= (_info->cgbk.icgcnt - 1); i++ ){
      /*        i_ = i - 1;  unreferenced??? */
                j = (_info->cgbk.icgptr + i)%_info->cgbk.memcg;
                ytd = 0.0e0;
                for( k = 1; k <= _info->nintbk.nsuper; k++ ){
                        k_ = k - 1;
                        ytd = ytd + Y(j,k_)*d[k_];
                }
                betaj = WORK(0,j)*ytd;
                for( ii = 1; ii <= _info->nintbk.nsuper; ii++ ){
                        ii_ = ii - 1;
                        d[ii_] = d[ii_] + S(j,ii_)*(WORK(1,i) - betaj);
                }
        }

        goto L_210;


        /*      -------------------
         *     | FLETCHER-REEVES
         *      ------------------- */
L_30:
        if( *msgcg > 0 )
             {
                lsgrg_msg(IOINFO, " FLETCHER-REEVES DIRECTION WILL BE USED \n" );
             }
        g1n = 0.0e0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                g1n = g1n + powi(gradfp[i_], 2);
        }
        if( g1n <= _info->tols.eps )
                goto L_9;
        g2n = 0.0e0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                g2n = g2n + powi(gradf[i_], 2);
        }
        cgbeta = g2n/g1n;
        if( _info->nintbk.ipr < 5 )
                goto L_190;

#ifdef IO_ENABLED
        sprintf(LSGRG_MSGBUFFER, " CGBETA =%15.8e\n", cgbeta );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        lsgrg_msg(IOINFO, " GRADF =" );
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                sprintf(LSGRG_MSGBUFFER, "%15.8e", gradf[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
        lsgrg_msg(IOINFO, "\n" );
        lsgrg_msg(IOINFO, " GRADFP=" );
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                sprintf(LSGRG_MSGBUFFER, "%15.8e", gradfp[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
        lsgrg_msg(IOINFO, "\n" );
#endif
        goto L_190;

        /*      ----------------
         *     | POLAK-RIBIERE
         *      ---------------- */
L_60:
        if( *msgcg > 0 )
            {
               lsgrg_msg(IOINFO, "\n POLAK-RIBIERE DIRECTION WILL BE USED \n" );
            }
        g1n = 0.0e0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                g1n = g1n + powi(gradfp[i_], 2);
        }
        if( g1n <= _info->tols.eps )
                goto L_9;
        gty = 0.0e0;
        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;
                g1 = gradfp[j_];
                g2 = gradf[j_];
                gty = gty + g2*(g2 - g1);
        }
        cgbeta = gty/g1n;
        goto L_190;

        /*      ---------
         *     | PERRY
         *      --------- */
L_90:
        if( *msgcg > 0 )
            {
              lsgrg_msg(IOINFO, "\n PERRY (MARCH 76) DIRECTION WILL BE USED \n" );
            }
        gyd = 0.0e0;
        ytd = 0.0e0;
        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;
                yj = gradf[j_] - gradfp[j_];
                gyd = gyd + gradf[j_]*(yj - _info->bestbk.step*d[j_]);
                ytd = ytd + yj*d[j_];
        }
        if( fabs( ytd ) <= _info->tols.eps )
                goto L_9;
        cgbeta = gyd/ytd;
        goto L_190;

        /*      -----
         *     | DFP
         *      ----- */
L_110:
        if( *msgcg > 0 )
            {
               lsgrg_msg(IOINFO, "\n DFP DIRECTION WILL BE USED \n" );
            }
        goto L_130;

        /*      ----------------------
         *     | COMPLEMENTARY DFP
         *      ---------------------- */
L_120:
        if( *msgcg > 0 )
            {
               lsgrg_msg(IOINFO, " COMPLEMENTARY DFP DIRECTION WILL BE USED \n" );
            }
L_130:
        gty = 0.0e0;
        ytd = 0.0e0;
        yty = 0.0e0;
        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;
                g1 = gradfp[j_];
                g2 = gradf[j_];
                yj = g2 - g1;
                gty = gty + g2*yj;
                ytd = ytd + yj*d[j_];
                yty = yty + powi(yj, 2);
        }
        if( fabs( ytd ) <= _info->tols.eps )
                goto L_9;
        if( fabs( yty ) <= _info->tols.eps )
                goto L_9;
        gtd = 0.0e0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                gtd = gtd + gradf[i_]*d[i_];
        }
        if( _info->cgbk.modcg == 5 )
                goto L_160;

        cgbeta = -_info->bestbk.step*gtd/ytd;
        gamma = gty/yty;
        goto L_170;

L_160:
        rho = _info->bestbk.step + yty/ytd;
        cgbeta = (gty - rho*gtd)/ytd;
        gamma = gtd/ytd;

L_170:
        itncg = itncg + 1;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " CGBETA =%14.7e GAMMA =%14.7e YTD =%14.7e YTY =%14.7e GTD =%14.7e GTY =%14.7e\n",
                 cgbeta, gamma, ytd, yty, gtd, gty );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;
                d[j_] = -gradf[j_] + cgbeta*d[j_] + gamma*(gradf[j_] - gradfp[j_]);
        }
        goto L_210;

        /*      ---------------------------------------
         *     | Set up next CG-Type search direction
         *      --------------------------------------- */
L_190:
        itncg = itncg + 1;
        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;
                d[j_] = -gradf[j_] + cgbeta*d[j_];
        }
L_210:
        *msgcg = 0;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                lsgrg_msg(IOINFO, " D =" );
                for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%15.8e", d[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
                }
        if( _info->nintbk.ipr >= 5 )
                {
                lsgrg_msg(IOINFO, " CG COMPLETED \n" );
                }
#endif
        return;




        /*     End of CG
         * */
#undef  WORK
#undef  S
#undef  Y

#ifdef initcg
    #undef initcg
#endif
#ifdef itncg
    #undef itncg
#endif

} /* end of function */

void /*FUNCTION*/ LSGRGCLASS comdfp(LsgrgInfo *_info,double r[], double y[], double gg[], long int ny,
         double theta, double *yty, double ytp, double *gtp, long int *kadd,
         long int *ksub)
{
long int i, k;
double c1, c2, gamin;

        /*      ----------------------------------------------
         *     | Modify R using Complementary DFP formula
         *      ---------------------------------------------- */
        if( fabs( theta ) <= _info->tols.eps )
                return;
        if( fabs( ytp ) < _info->tols.eps )
                return;
        if( fabs( *gtp ) < _info->tols.eps )
                *gtp = -_info->tols.eps;
        if( fabs( *yty ) < _info->tols.eps )
                *yty = _info->tols.eps;

        /*      -----------------------------------------------
         *     | Scale the Identity if last update was a RESET
         *      ----------------------------------------------- */
        _info->hesblk.gamma0 = 1.0e0;
        gamin = *yty/(theta*ytp);
        if( ((gamin >= 1.0e-7 && gamin <= 1.0e7) && _info->logblk.resetp) &&
         _info->cgbk.hscale ){
                _info->hesblk.gamma0 = 1.0e0/gamin;
                gamin = sqrt( gamin );
                k = 0;
                for( i = 1; i <= ny; i++ ){
              /*        i_ = i - 1;  unreferenced??? */
                        k = k + i;
                        r[k - 1] = gamin;
                }
        }

        c1 = 1.0e0/(theta*ytp);
        c2 = 1.0e0/ *gtp;
        r1mod(_info, r, y, ny, c1, kadd );
        r1mod(_info, r, gg, ny, c2, ksub );

        return;

        /*     End of COMDFP
         * */
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS delcol(LsgrgInfo *_info,double r[], long int jq)
{
long int i, i1,  j, j1, j2, j_, k, k1, ns;
double cs, d, sn, t1, t2;

        /*      ------------------------------------------------
         *     | Delete the JQ'th column of upper triangular R
         *      ------------------------------------------------ */
        ns = _info->nintbk.nsuper;
        if( jq >= ns )
                goto L_60;
        k = jq*(jq + 1)/2;
        i1 = jq + 1;
        for( i = i1; i <= ns; i++ ){
        /*      i_ = i - 1; unreferenced??? */
                k = k + i;
                t1 = r[k - 2];
                t2 = r[k - 1];
                d = sqrt( t1*t1 + t2*t2 );
                r[k - 2] = d;
                if( i == ns )
                        goto L_20;
                cs = t1/d;
                sn = t2/d;
                j1 = i + 1;
                k1 = k + i;
                for( j = j1; j <= ns; j++ ){
          /*            j_ = j - 1; unreferenced??? */
                        t1 = r[k1 - 2];
                        t2 = r[k1 - 1];
                        r[k1 - 2] = cs*t1 + sn*t2;
                        r[k1 - 1] = sn*t1 - cs*t2;
                        k1 = k1 + j;
                }
L_20:
                k1 = i - 1;
                j2 = k - i;
                j1 = j2 - i + 2;
                for( j = j1; j <= j2; j++ ){
                        j_ = j - 1;
                        r[j_] = r[j_ + k1];
                }
        }
L_60:
        _info->nintbk.nsuper = _info->nintbk.nsuper - 1;
        return;


        /*     End of DELCOL
         * */
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS direc(LsgrgInfo *_info,double gradf[], double r[], double gradfp[],
         long int iub[], long int inbv[], double d[], long int istat[],
         double y[], long int ivstat[], long int ibv[], long int lksame[],
         double *ycg, double *scg, double *work)
{
#define YCG(I_,J_)      (*(ycg+(I_)*(_info->dimen.n)+(J_)))
#define SCG(I_,J_)      (*(scg+(I_)*(_info->dimen.n)+(J_)))
#define WORK(I_,J_)     (*(work+(I_)*(_info->cgbk.mcgm1-(0)+1)+(J_)))
LOGICAL32 first, modr, newcg, savres, sbdel;
/* note: compiler flagged iiir as unneeded, assignments,but no use */
/* 1/3/01 jcp take msgcg out of local defs -- added to statics */
long int i, i_, ii, ii_, iii, /* iiir, */ imax, iubi, j, j_, k, k_, kadd,
         kk,  ksub, /* msgcg,*/ nsb, nsbo;
double atst, di, dmax, dmin, gfi, gfpi, gfract, gmax, gtp, gtp1, p,
         rmean, rmnsq, t, tmax, told, tst, yi, ytp, yty;

/* replace statics by alg members */

    #define nsupvm _info->direcStatic.nsupvm
    #define msgvm  _info->direcStatic.msgvm
    #define rtinsb _info->direcStatic.rtinsb
    #define msgcg  _info->direcStatic.msgcg

/*
   static long nsupvm = -1;
   static long msgvm = 0;
   static double rtinsb = 4.0e0;
*/
        modr = FALSE; /* 1/3/01 jcp uninitialized and messed up test */
                      /* at l_990 when cg used */
        if( _info->nintbk.ipr >= 5 )
                lsgrg_msg(IOINFO,"\n DIREC ENTERED \n" );

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " NSUPER = %5ld MAXH = %5ld NSEAR = %5ld STEP = %13.7e\n",
                 _info->nintbk.nsuper, _info->misc.maxh, _info->misc.nsear, _info->bestbk.step );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        sbdel = FALSE;
        _info->dirgrg.cond = 0.0e0;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                lksame[ibv[i_] - 1] = 0;
        }

#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " DROP = %1c VARMET = %1c CONJGR = %1c JSTFES = %1c\n",
                 TorF(_info->logblk.drop), TorF(_info->logblk.varmet), TorF(_info->logblk.conjgr),
                 TorF(_info->srchlg.jstfes) );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( _info->srchlg.jstfes )
                _info->logblk.restrt = TRUE;
        if( _info->logblk.drop )
                goto L_80;
        /*C */
        _info->logblk.move = _info->bestbk.step > _info->tols.eps;

        /*      ----------------------------------------------------------------
         *     | UPDATE is flag indicating type of direction formula used:
         *     |     UPDATE  = 1   -->  No update of Hessian or Neg. Gradient
         *     |     UPDATE  = 2   -->  QN Hessian update
         *     |     UPDATE  = 3   -->  CG  update
         *      ----------------------------------------------------------------
         * */
L_5:
        _info->dirgrg.update = 1;
        savres = FALSE;
        _info->dfblk.dfail = TRUE;
        first = TRUE;
        rmean = 1.0e0;
        _info->hesblk.dirsc = 1.0e0;
        _info->slpobj.slope = 0.0e0;
        _info->logblk.varmet = _info->nintbk.nsuper <= _info->misc.maxh;
        _info->logblk.conjgr = !_info->logblk.varmet;
        if( _info->nintbk.nsuper == 0 ){
                _info->logblk.restrt = TRUE;
                goto L_146;
                /*XXXX    GO TO 500 */
        }
        if( _info->logblk.conjgr )
                goto L_80;

        /*      -----------------------------------------------
         *     | Complementary DFP Variable Metric Algorithm
         *      ----------------------------------------------- */
        if( _info->nintbk.ipr > 1 )
                msgcg = -1;
#ifdef IO_ENABLED
        if( msgvm > 0 )
                {
                sprintf(LSGRG_MSGBUFFER, " HESSIAN SIZE IS LESS THAN %5ld -- SWITCH TO VARIABLE METRIC ROUTINE \n",
                 _info->misc.maxh );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        msgvm = 0;
        if( _info->misc.nsear == 0 ){
#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 3 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " NSEAR = %4ld\n", _info->misc.nsear );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                goto L_60;
        }
        if( !_info->supblk.sbchng ){
#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 3 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " RUPDT = %1c\n", TorF(_info->supblk.sbchng) );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                goto L_60;
        }
        if( _info->logblk.restrt ){
#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 3 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " RESTRT = %1c\n", TorF(_info->logblk.restrt) );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                goto L_60;
        }
        if( !_info->logblk.move ){
#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 3 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " MOVE = %1c\n", TorF(_info->logblk.move) );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                goto L_60;
        }

        /*      ---------------------------------------------------
         *     | MSGCG causes CG subroutine to print on next call
         *      --------------------------------------------------- */
        if( _info->nintbk.ipr3 >= 0 )
                msgcg = 1;
        modr = FALSE;

        /*      ------------------------------------------------------------------
         *     | Use COMDFP update if no more than one superbasic has hit a bound
         *     |            ( STEP  .LE . STEPMX)
         *      ------------------------------------------------------------------ */
        if( _info->bestbk.step > _info->bestbk.stepmx )
                goto L_70;
        ytp = 0.0e0;
        gtp = 0.0e0;
        gtp1 = 0.0e0;
        yty = 0.0e0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                p = d[i_];
                yi = gradf[i_] - gradfp[i_];
                y[i_] = yi;
                ytp = yi*p + ytp;
                yty = yi*yi + yty;
                gtp = gradfp[i_]*p + gtp;
                gtp1 = gradf[i_]*p + gtp1;
        }

        /*      --------------------------------------------------------------------
         *     | Use COMDFP update only if GTP1/GTP < 0.9(say).  (Note that GTP<0)
         *      -------------------------------------------------------------------- */
        gfract = 0.99e0;
        modr = gtp1 > gfract*gtp;
        modr = modr && _info->bestbk.step <= _info->bestbk.stepmx;
        savres = _info->logblk.resetp && (!modr);
#ifdef IO_ENABLED
        if( !modr && _info->nintbk.ipr >= 3 )
                {
                sprintf(LSGRG_MSGBUFFER, " MODR FALSE, SKIP UPDATE OF HESSIAN \n GTP1 =%13.6e GTP =%13.6e\n",
                 gtp1, gtp );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
        if( _info->nintbk.ipr > 6 )
                {
                lsgrg_msg(IOINFO, " Y = " );
                for( i = 1; i <= _info->dimen.n; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%14.7e", y[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, " \n" );
                }
#endif
        goto L_70;

        /*      --------------------------------
         *     | Reset Hessian to the Identity
         *      -------------------------------- */
L_60:
        ;
        resetr(_info, r, _info->dimen.nrtot, &_info->dirgrg.cond );
        _info->dirgrg.update = 1;
    /*    iiir = 1; */
        goto L_80;

        /*      ----------------------------------------------
         *     | Hessian update when no variable hits bound
         *      ---------------------------------------------- */
L_70:
        ;
        if( !modr )
                goto L_80;
        kadd = 4;
        ksub = 4;
        comdfp(_info, r, y, gradfp, _info->nintbk.nsuper, _info->bestbk.step, &yty, ytp, &gtp,
         &kadd, &ksub );
    /*     iiir = 0; */
        _info->dirgrg.update = 2;

        /*      --------------------------------------------
         *     | Now, delete any Superbasics at bounds
         *      -------------------------------------------- */
L_80:
        ;
        if( sbdel )
                goto L_85;
        sbdel = TRUE;
        nsbo = _info->nintbk.nsuper;
        nsb = _info->nintbk.nsuper;
        for( kk = 1; kk <= nsbo; kk++ ){
       /*       kk_ = kk - 1; unreferenced??? */
                i = nsbo + 1 - kk;
                if( iub[i - 1] == 0 )
                        goto L_390;

                /*          -----------------------------------
                 *         | Make this Superbasic a Nonbasic
                 *          ----------------------------------- */
                iii = nsb;
                nsb = nsb - 1;
                if( _info->logblk.varmet )
                        delcol(_info, r, i );
                if( i > nsb )
                        goto L_390;
                j = inbv[i - 1];
                di = d[i - 1];
                gfi = gradf[i - 1];
                gfpi = gradfp[i - 1];
                iubi = iub[i - 1];
                for( k = i; k <= nsb; k++ ){
                        k_ = k - 1;
                        inbv[k_] = inbv[k_ + 1];
                        d[k_] = d[k_ + 1];
                        gradf[k_] = gradf[k_ + 1];
                        gradfp[k_] = gradfp[k_ + 1];
                        iub[k_] = iub[k_ + 1];
                }
                inbv[iii - 1] = j;
                d[iii - 1] = di;
                gradf[iii - 1] = gfi;
                gradfp[iii - 1] = gfpi;
                iub[iii - 1] = iubi;
L_390:
                ;
        }

        _info->nintbk.nsuper = nsb;
        if( _info->logblk.drop )
                goto L_500;
        if( _info->nintbk.nsuper == 0 )
                goto L_146;
        newcg = _info->logblk.conjgr && _info->nintbk.nsuper < nsbo;
        _info->logblk.restrt = _info->logblk.restrt || newcg;
        if( newcg && _info->nintbk.nsuper <= _info->misc.maxh )
                goto L_5;
L_85:
        ;
        if( _info->logblk.conjgr )
                goto L_110;

        /*      -----------------------------
         *     | Compute search direction,D
         *      ----------------------------- */
        for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                j_ = j - 1;
                d[j_] = -gradf[j_];
        }
        rtrsol(_info, r, d, _info->nintbk.nsuper );

        /*      -------------------------------------------------
         *     | Compute lower bound on condition number of R
         *      ------------------------------------------------- */
        dmin = _info->tols.plinfy;
        dmax = 0.0e0;
        k = 0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
        /*      i_ = i - 1;  unreferenced */
                k = k + i;
                t = fabs( r[k - 1] );
                if( dmin > t )
                        dmin = t;
                if( dmax < t )
                        dmax = t;
        }
        _info->dirgrg.cond = _info->tols.plinfy;
        if( dmin >= _info->tols.eps ){
                _info->dirgrg.cond = powi(dmax/dmin, 2);
                rmnsq = dmax*dmin;
                rmean = sqrt( rmnsq );
                _info->hesblk.dirsc = 1.0e0/rmnsq;
        }
#ifdef IO_ENABLED
        if( _info->nintbk.ipr3 >= 5 ){
                k = _info->nintbk.nsuper*(_info->nintbk.nsuper + 1)/2;
                lsgrg_msg(IOINFO, " R IS" );
                for( i = 1; i <= k; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%15.7e", r[i - 1] );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
        }
#endif
        goto L_120;

        /*      -----------------------------
         *     | Conjugate Gradient method
         *      ----------------------------- */
L_110:
        ;
        if( _info->nintbk.ipr > 1 )
                msgvm = 1;
        if( ((!_info->srchlg.uncon || (!_info->supblk.sbchng)) || (!_info->logblk.move)) ||
         _info->bestbk.step > _info->bestbk.stepmx )
                _info->logblk.restrt = TRUE;
        cg(_info, gradfp, gradf, d, &msgcg, &_info->dirgrg.update, ycg, scg, work );

        /*      -----------------------------------
         *     | Check if direction is downhill
         *      ----------------------------------- */
L_120:
        ;
        _info->slpobj.slope = 0.0e0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                _info->slpobj.slope = _info->slpobj.slope + d[i_]*gradf[i_];
        }
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, "  ORIGINAL SLOPE = %16.8e\n",
                 _info->slpobj.slope );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( _info->slpobj.slope < -_info->tols.eps )
                goto L_145;

        /*      ---------------------------
         *     | Bad Direction.  Reset
         *      --------------------------- */
        if( _info->logblk.restrt ){
                if( _info->nintbk.ipr > 0 )
                    {
                        lsgrg_msg(IOINFO,
                           "\nDIRECTION VECTOR < EPS .. TRY DROPPING A CONSTRAINT\n" );
                    }
                goto L_146;
        }
        if( _info->nintbk.ipr > 1 )
                {
                lsgrg_msg(IOINFO, "\n DIRECTION NOT DOWNHILL, RESET \n" );
                }
        _info->logblk.restrt = TRUE;
        if( _info->logblk.varmet )
                goto L_60;
        goto L_110;

        /*      ------------------------------------------------------------
         *     | Decide if any variables should be released from bounds
         *      ------------------------------------------------------------ */
L_145:
        ;
        if( !_info->logblk.move )
                goto L_220;
L_146:
        ;
L_149:
        told = 0.0e0;
        j = _info->nintbk.nsuper + 1;
        gmax = 0.0e0;
        imax = 0;
        if( j > _info->dimen.n )
                goto L_220;
        for( ii = j; ii <= _info->dimen.n; ii++ ){
                ii_ = ii - 1;
                i = inbv[ii_];
                if( i > _info->dimen.n )
                        goto L_403;

                /*         -----------------------
                 *        | Skip fixed variables
                 *         ----------------------- */
                if( ivstat[i - 1] < 0 )
                        goto L_430;
                goto L_405;

                /*         ----------------------------------
                 *        | Slack variable processed below
                 *         ---------------------------------- */
L_403:
                k = i - _info->dimen.n;
                if( istat[k - 1] == 1 )
                        goto L_430;

                /*         -------------------------------------------
                 *        | Regular variables and inequality slacks
                 *         ------------------------------------------- */
L_405:
                ;
                tst = gradf[ii_];
                atst = fabs( tst );
                if( iub[ii_] == 1 )
                        goto L_410;

                /*         ------------------------------
                 *        | Variable at lower bound
                 *         ------------------------------ */
                if( tst >= -told )
                        goto L_427;
                if( atst > gmax )
                        goto L_420;
                goto L_425;

                /*         -------------------------
                 *        | Variable at upper bound
                 *         ------------------------- */
L_410:
                ;
                if( tst <= told )
                        goto L_427;
                if( atst <= gmax )
                        goto L_425;

                /*         -------------------
                 *        | Save largest
                 *         ------------------- */
L_420:
                ;
                gmax = atst;
                imax = ii;
L_425:
                if( first )
                        lksame[i - 1] = lksame[i - 1] + 1;
                goto L_430;
L_427:
                lksame[i - 1] = 0;

L_430:
                ;
        }

        /*      --------------------------------------
         *     | Test for insertion into Super Basics
         *      -------------------------------------- */
        if( imax == 0 )
                goto L_220;
        first = FALSE;
        i = inbv[imax - 1];
        tst = (powi(gmax*lksame[i - 1], 2))*rtinsb*_info->hesblk.dirsc;
        if( tst < fabs( _info->slpobj.slope ) )
                goto L_220;

        /*      --------------------------------
         *     | Make this variable Superbasic
         *      -------------------------------- */
        tmax = gradf[imax - 1];
        iubi = iub[imax - 1];
        _info->nintbk.nsuper = _info->nintbk.nsuper + 1;
        _info->logblk.varmet = _info->nintbk.nsuper <= _info->misc.maxh;
        _info->logblk.conjgr = !_info->logblk.varmet;
        inbv[imax - 1] = inbv[_info->nintbk.nsuper - 1];
        inbv[_info->nintbk.nsuper - 1] = i;
        gradf[imax - 1] = gradf[_info->nintbk.nsuper - 1];
        gradf[_info->nintbk.nsuper - 1] = tmax;
        iub[imax - 1] = iub[_info->nintbk.nsuper - 1];
        iub[_info->nintbk.nsuper - 1] = iubi;
        lksame[i - 1] = 0;
        if( _info->logblk.varmet ){
                addcol(_info, r, y, _info->dimen.n, rmean );
        } else if( _info->nintbk.nsuper - 1 == _info->misc.maxh ){
                _info->logblk.restrt = TRUE;
                goto L_110;
        } else if( _info->cgbk.modcg == 6 ){
                for( ii = 0; ii <= _info->cgbk.mcgm1; ii++ ){
             /*          ii_ = ii - 1; unreferenced */
                        YCG(ii,_info->nintbk.nsuper - 1) = 0.0e0;
                        SCG(ii,_info->nintbk.nsuper - 1) = 0.0e0;
                }
        }
        d[_info->nintbk.nsuper - 1] = -tmax*_info->hesblk.dirsc;
        _info->slpobj.slope = _info->slpobj.slope + tmax*d[_info->nintbk.nsuper - 1];
        _info->dfblk.dfail = FALSE;
        goto L_149;

        /*      ----------------------------------
         *     |  Release all possible Nonbasics
         *      ---------------------------------- */
L_500:
        ;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 1 )
                {
                sprintf(LSGRG_MSGBUFFER, " DIREC...DROPPING ALL NONBASICS - DROP = %2c  NSUPER = %5ld\n",
                 TorF(_info->logblk.drop), _info->nintbk.nsuper );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

        told = 0.0e0;

        j = _info->nintbk.nsuper + 1;
        if( j > _info->dimen.n )
                goto L_540;
        for( ii = j; ii <= _info->dimen.n; ii++ ){
                ii_ = ii - 1;
                i = inbv[ii_];
                if( i > _info->dimen.n )
                        goto L_503;

                /*         ------------------------
                 *        | Skip fixed variables
                 *         ------------------------ */
                if( ivstat[i - 1] < 0 )
                        goto L_530;
                goto L_505;

                /*         --------------------------------
                 *        | Slack variable processed below
                 *         -------------------------------- */
L_503:
                k = i - _info->dimen.n;
                if( istat[k - 1] == 1 )
                        goto L_530;

                /*         -------------------------------------------
                 *        | Regular variables and inequality slacks
                 *         ------------------------------------------- */
L_505:
                ;
                tst = gradf[ii_];
                if( iub[ii_] == 1 )
                        goto L_510;

                /*         ---------------------------
                 *        | Variable at lower bound
                 *         --------------------------- */
                if( tst >= -told )
                        goto L_530;
                if( i > _info->dimen.n && tst >= -10.0e1*told )
                        goto L_530;
                goto L_520;

                /*         ---------------------------
                 *        | Variable at upper bound
                 *         --------------------------- */
L_510:
                ;
                if( tst <= told )
                        goto L_530;
                if( i > _info->dimen.n && tst <= 10.0e1*told )
                        goto L_530;

                /*         -----------------------------------
                 *        | Make this variable Superbasic
                 *         ----------------------------------- */
L_520:
                ;
                _info->nintbk.nsuper = _info->nintbk.nsuper + 1;
                _info->logblk.varmet = _info->nintbk.nsuper <= _info->misc.maxh;
                inbv[ii_] = inbv[_info->nintbk.nsuper - 1];
                inbv[_info->nintbk.nsuper - 1] = i;
                gradf[ii_] = gradf[_info->nintbk.nsuper - 1];
                gradf[_info->nintbk.nsuper - 1] = tst;
                iubi = iub[ii_];
                iub[ii_] = iub[_info->nintbk.nsuper - 1];
                iub[_info->nintbk.nsuper - 1] = iubi;
                if( _info->logblk.varmet )
                        addcol(_info, r, y, _info->dimen.n, rmean );
                d[_info->nintbk.nsuper - 1] = -gradf[_info->nintbk.nsuper - 1];
                _info->dfblk.dfail = FALSE;

L_530:
                ;
        }
        _info->logblk.conjgr = !_info->logblk.varmet;
L_540:
        ;
        if( _info->dfblk.dfail && _info->logblk.restrt )
                goto L_557;

        if( _info->dfblk.dfail && _info->nintbk.ipr > 0 )
                {
                lsgrg_msg(IOINFO,
                 " COULD NOT DROP ANY CONSTRAINT. TRY -VE GRADIENT DIRECTION.\n" );
                }

        /*      --------------------------------------------------
         *     | Update direction vector and calculate new SLOPE
         *      -------------------------------------------------- */
        _info->logblk.restrt = TRUE;
        _info->dirgrg.update = 1;
        if( _info->logblk.varmet ){
                resetr(_info, r, _info->dimen.nrtot, &_info->dirgrg.cond );
                for( j = 1; j <= _info->nintbk.nsuper; j++ ){
                        j_ = j - 1;
                        d[j_] = -gradf[j_];
                }
        } else{
                cg(_info, gradfp, gradf, d, &msgcg, &_info->dirgrg.update, ycg, scg, work );
        }

        _info->slpobj.slope = 0.0e0;
        for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                i_ = i - 1;
                _info->slpobj.slope = _info->slpobj.slope + d[i_]*gradf[i_];
        }

        if( _info->slpobj.slope < -_info->tols.eps )
                goto L_220;

        /*      -----------------------------------
         *     | Negative Gradient is not downhill
         *      ----------------------------------- */
L_557:
        if( _info->nintbk.ipr > 0 )
            {
               lsgrg_msg(IOINFO,
                  "\nNEGATIVE GRADIENT DIRECTION NOT DOWNHILL.  CHECK DERIVATIVES AND/OR TOLERANCES. \n" );
            }

        _info->dfblk.dfail = TRUE;
        goto L_999;

L_220:
        ;
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                gradfp[i_] = gradf[i_];
        }
        if( _info->logblk.varmet )
                nsupvm = _info->nintbk.nsuper;
        _info->dirgrg.nsupp = _info->nintbk.nsuper;
        _info->dfblk.dfail = FALSE;
        /*     WRITE(IOOUT,8890) NSUPER */
        goto L_999;

        /*      ----------------------------------------------------------
         *     | Direction Vector < EPS -- Try dropping a constraint.
         *      ---------------------------------------------------------- */
/* L_235: unreferenced??? */
        _info->logblk.drop = TRUE;
        if( _info->nintbk.ipr > 0 )
                {
                lsgrg_msg(IOINFO,
                   "\nDIRECTION VECTOR < EPS .. TRY DROPPING A CONSTRAINT\n" );
                }
        goto L_500;
L_999:
        ;

        _info->logblk.resetp = _info->logblk.restrt || savres;
        if( (_info->logblk.restrt && _info->dirgrg.update != 1) || ((!_info->logblk.restrt &&
         _info->dirgrg.update == 1) && modr) ){
#ifdef IO_ENABLED
                sprintf(LSGRG_MSGBUFFER, " ***DIREC...ERROR: RESTRT = %2c BUT UPDATE = %3ld\n",
                 TorF(_info->logblk.restrt), _info->dirgrg.update );
                 lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
        /*      xerror( 0, iounit.ioout ); */
                lsgrg_errorexit(JUMPBUF,_LSGRG_DIREC_UPDATE_ERR);
                return;
        }
        _info->logblk.drop = FALSE;
        _info->logblk.restrt = FALSE;

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, "  FINAL SLOPE = %16.8e\n", _info->slpobj.slope );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
        if( _info->nintbk.ipr > 4 )
                {
                lsgrg_msg(IOINFO, " DIREC COMPLETED \n" );
                }
#endif
        return;
        /*----------------------------------------------------------------------
         *
         *     End of DIREC
         * */
#undef  WORK
#undef  SCG
#undef  YCG
    #undef nsupvm
    #undef msgvm
    #undef rtinsb
    #undef msgcg

} /* end of function */


void /*FUNCTION*/ LSGRGCLASS resetr(LsgrgInfo *_info,double r[], long int nrtot, double *cond)
{
long int i, i_, k, ncol;

        /* .................................................
         *
         *   -----------------------------------------------
         *  | Reset the Cholesky factor of the Hessian
         *   ----------------------------------------------- */
        if( _info->nintbk.nsuper == 0 )
                return;
        if( _info->nintbk.ipr > 3 )
                {
                lsgrg_msg(IOINFO,
                   " CHOLESKY FACTOR OF HESSIAN RESET TO I. \n" );
                }
        *cond = 1.0e0;
        ncol = _info->nintbk.nsuper;
        if( _info->misc.maxh < ncol )
                ncol = _info->misc.maxh;
        k = ncol*(ncol + 1)/2;
        for( i = 1; i <= k; i++ ){
                i_ = i - 1;
                r[i_] = 0.0e0;
        }
        k = 0;
        for( i = 1; i <= ncol; i++ ){
         /*     i_ = i - 1; unreferenced??? */
                k = k + i;
                r[k - 1] = 1.0e0;
        }

        _info->logblk.restrt = TRUE;
        return;

        /*     End of RESETR
         * */
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS rtrsol(LsgrgInfo *_info,double r[], double y[], long int ny)
{
long int i, i1, i_, ii,  j, j_, k;
double s, t;

        /* ..............................................................
         *      -------------------------------------------------------
         *     | Solve R'R*Y = Y (R triangular, stored by columns)
         *      ------------------------------------------------------- */
        y[0] = y[0]/r[0];
        k = 1;
        if( ny <= 1 )
                goto L_50;
        for( i = 2; i <= ny; i++ ){
                i_ = i - 1;
                s = y[i_];
                i1 = i - 1;
                for( j = 1; j <= i1; j++ ){
                        j_ = j - 1;
                        k = k + 1;
                        s = s - r[k - 1]*y[j_];
                }
                k = k + 1;
                y[i_] = s/r[k - 1];
        }

L_50:
        for( ii = 1; ii <= ny; ii++ ){
       /*       ii_ = ii - 1; unreferenced???*/
                i = ny + 1 - ii;
                t = y[i - 1]/r[k - 1];
                y[i - 1] = t;
                if( i <= 1 )
                        goto L_80;
                k = k - i;
                i1 = i - 1;
                for( j = 1; j <= i1; j++ ){
                        j_ = j - 1;
                        y[j_] = y[j_] - r[k + j_]*t;
                }
L_80:
                ;
        }
        return;

        /*     End of RTRSOL
         * */
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS r1mod(LsgrgInfo *_info,double r[], double y[], long int ny, double sigma,
         long int *k)
{
LOGICAL32 posdef;
long int i, i_;
double s, t, tol;

        /* ...................................................
         *      --------------------------------------------
         *     | Modify R such that R'R := R'R + SIGMA*YY'
         *      -------------------------------------------- */
        tol = _info->epscom.eps0;
        s = 0.0e0;
        t = sqrt( fabs( sigma ) );
        for( i = 1; i <= ny; i++ ){
                i_ = i - 1;
                s = powi(y[i_], 2) + s;
                y[i_] = y[i_]*t;
        }
        s = sigma*s;
        if( fabs( s ) <= tol )
                return;
        if( s <= tol )
                goto L_200;
        r1add(_info, r, y, ny );
        return;

L_200:
        r1sub(_info, r, y, ny, tol, &posdef );
        if( !posdef )
                *k = -*k;
        return;

        /*     End of R1MOD
         * */
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS r1add(LsgrgInfo *_info,double r[], double y[], long int ny)
{
long int j, j1, j_, k, k0, k1, k_, l;
double cs, d, sn, t1, t2;


        /*      -----------------------------------------------------------------
         *     | Modify R such that R'R := R_R + YY' where R is upper-triangular,
         *     | stored by columns.
         *      ----------------------------------------------------------------- */
        k1 = 1;
        for( k = 1; k <= ny; k++ ){
                k_ = k - 1;
                k0 = k1;
                t1 = r[k1 - 1];
                t2 = y[k_];
                d = sqrt( t1*t1 + t2*t2 );
                cs = t1/d;
                sn = t2/d;
                j1 = k1 + k;
                k1 = j1 + 1;
                if( fabs( sn ) <= _info->epscom.eps0 )
                        goto L_100;
                r[k0 - 1] = d;
                l = k + 1;
                if( l > ny )
                        goto L_100;
                for( j = l; j <= ny; j++ ){
                        j_ = j - 1;
                        t1 = r[j1 - 1];
                        t2 = y[j_];
                        r[j1 - 1] = cs*t1 + sn*t2;
                        y[j_] = sn*t1 - cs*t2;
                        j1 = j1 + j;
                }
L_100:
                ;
        }

/* L_900: unreferenced??? */
        return;

        /*     End of R1ADD
         * */
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS r1sub(LsgrgInfo *_info,double r[], double y[], long int ny, double tol,
         LOGICAL32 *posdef)
{
long int i, i1, i_, ii,  j, j_, k, l;
double cs, d, s, sn, t, t1, t2, u;

        /*     ----------------------------------------------------------------
         *    | Modify R such that R'R := R'R - YY' where R is upper-triangular,
         *    | stored by columns.
         *    | See Saunders, Stanford Technical Report STAN-CS-72-252, Chapter 7.
         *     ----------------------------------------------------------------- */
        /*      -----------------------
         *     | First solve R'P = Y
         *      ----------------------- */
        t = y[0]/r[0];
        d = t*t;
        y[0] = t;
        k = 1;
        if( ny <= 1 )
                goto L_50;
        for( i = 2; i <= ny; i++ ){
                i_ = i - 1;
                s = y[i_];
                i1 = i - 1;
                for( j = 1; j <= i1; j++ ){
                        j_ = j - 1;
                        k = k + 1;
                        s = s - r[k - 1]*y[j_];
                }
                k = k + 1;
                t = s/r[k - 1];
                d = t*t + d;
                y[i_] = t;
        }

        /*      -------------------------------------------
         *     | See if new R will be positive definite.
         *      ------------------------------------------- */
L_50:
        d = 1.0e0 - d;
        *posdef = d > _info->epscom.eps0;
        if( !*posdef )
                return;
        s = sqrt( d );

        /*      -----------------------------------------------
         *     | Perform backward sweep of plane rotations.
         *      ----------------------------------------------- */
        for( ii = 1; ii <= ny; ii++ ){
       /*       ii_ = ii - 1; unreferenced??? */
                i = ny + 1 - ii;
                u = s;
                t = y[i - 1];
                d = t*t + d;
                s = sqrt( d );
                cs = u/s;
                sn = t/s;
                y[i - 1] = 0.0e0;
                l = k;
                k = k - i;
                if( fabs( sn ) <= _info->epscom.eps0 )
                        goto L_80;
                for( j = i; j <= ny; j++ ){
                        j_ = j - 1;
                        t1 = y[j_];
                        t2 = r[l - 1];
                        y[j_] = cs*t1 + sn*t2;
                        r[l - 1] = sn*t1 - cs*t2;
                        l = l + j;
                }
L_80:
                ;
        }
        return;

        /*     End of R1SUB
         * */
} /* end of function */

#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"

#define IOINFO  &_info->io_info
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)
/***********************************************************************
 **                 FILE:        GRGPHAS0 FORTRAN                      *
 **                 AUTHOR:      STUART SMITH                          *
 **                 LAST UPDATE: 24 MAR  1994                          *
 **                                                                    *
 **                                                                    *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **      This file contains the following routines for LSGRG2, a       *
 **  large scale implementation of GRG2:                               *
 **                                                                    *
 **  ROUTINE           DESCRIPTION                        LAST UPDATE  *
 ** ---------         -------------                      ------------- *
 **  PHAS0     Attempts to find a feasible point         24 MAR 1994   *
 **            using Newton's method to solve the                      *
 **            nonlinear equations                                     *
 **  PH0FAC    Performs an LU factorization of the       08 OCT 1993   *
 **            basis matrix at a new point                             *
 **                                                                    *
 **  PH0PIV    Performs a pivot operation to update the  24 MAR 1994   *
 **            basis when a basic variable reaches a                   *
 **            bound or after slow Newton's convergence.               *
 **                                                                    *
 **  SUMRES    Computes the constraint residuals and     24 MAR 1994   *
 **            L-0, L-1, and L-2 norms.                                *
 **                                                                    *
 **  DROPRW    Updates appropriate data when a           24 MAR 1994   *
 **            constraint is dropped from consideration.               *
 **                                                                    *
 **  PH0LOG    Prints a PHASE-0 iteration log            08 OCT 1993   *
 **                                                                    *
 **  **change log**                                                    *
 **  10/98 jcp fresh translation from fortran, cputime calls replaced  *
 **       with calls to new timer lsgrg_timer()                        *
 **                                                                    *
 **  *fixme** move these parm translations to a higher level           *
 ***********************************************************************
 * */
        /* PARAMETER translations */
#define ASMALL  1.0e-3
#define DELTAF  0.25e0
#define DELTAM  0.01e0
#define MXDEG   25
#define MXSMAL  3
#define REDSTP  0.80e0
#define STPMIN  1.0e-5
        /* end of PARAMETER translations */

void /*FUNCTION*/ LSGRGCLASS phas0(LsgrgInfo *_info,long int ibmap[], double zmem[], double grad[],
         long int ihag[], long int iheg[], long int ihegl[], long int ibv[],
         double x[], double alb[], double ub[], long int inbv[], long int iub[],
         long int ibc[], double g[], long int inlin[], long int ivstat[],
         long int istat[], double x0[], double gg[], double pvrow[], long int icols[],
         long int iwork[], double err0[], long int icand[], double cdnum[],
         long int irank[], double ascale[], long int ibvbct[], double delx[],
         double error[], LOGICAL32 *infes, long int icstat[], long int *linect,
         long int iobjpt[], double paij[], long int iprow[], long int ipcol[],
         long int ipmap[])
{
LOGICAL32 bschng, degen, newbas, singlr, slowc;
long int i, i_,   j, j_, jbnd, jmax, jmaxa, k, lv, lvbnd,
         lvv, nbfree, ndegen, ninfa, ninfao, ninfd, ninfdo, nloops, nslow,
         nsmall;

/*  7/26/03 jcp ** add iteration limit */
      long maxIterations;

double aerr, alpha, derr, dist, dmaxx, dxx, err, errmx, errmxo, sinf,
         sinfd, sinfdo, sinfo, t, tfinal, tstart;



        /*C.....................................................................
         *C...................... LOCAL PARAMETERS . ...........................
         *C.....................................................................
         *     STPMIN --  Smallest allowable step size. Steps less that STEPMN
         *                are considered degenerate.
         *     ASMALL --  Used to monitor consecutive small steps.  If the step
         *                is less than ASMALL for MXSMAL iterations, the basis is
         *                changed to attempt to allow larger steps.
         *     DELTAF --  Used to monitor slow Newton's convergence. Convergence
         *                is considered fast if the norm of the residuals satisfies
         *                SINFNEW  <=  (1 - DELTAF*ALPHA)*SINFOLD
         *                ALPHA = current step size.
         *    DELTAM  --  Used to force a minimum relative decrease in the residuals.
         *                If the for the current step, ALPHA,
         *                     SINFNEW  >=  (1 - DELTAM*ALPHA)*SINFOLD
         *                ALPHA will be reduced .
         *    REDSTP  --  The step reduction factor.  When necessary, the step is
         *                reduced by  ALPHA = ALPHA*REDSTP.
         *    MXDEG   --  The maximum number of consecutive degenerate steps.  After
         *                MXDEG consecutive degenerate steps, a constraint will be
         *                dropped.
         *   MXSMAL   -- Maximum number of consecutive steps less than ASMAL.  After
         *               MXSMAL such steps, a basis change is forced.
         *
         *   NOTE:       DELTAM and DELTAF must satisfy:
         *                    0 <= DELTAM <= DELTAF <= 1
         *.............................................................................
         * */


        /*.......................................................................
         *
         *     *** START OF ROUTINE: PHASE0 ***
         *C
         * */
/*      cputime( &tstart, &irc ); */
        tstart = lsgrg_timer();

        if( _info->nintbk.ipr >= 3 )
                {
                lsgrg_msg(IOINFO,"\nPHASE0 ENTERED\n" );
                }
        nloops = 0;
/*  ** 7/26/03 jcp ** move init of alpha to here.  the call   */
/*  to ph0log at 999 flags it as uninitialized (apparently    */
/*  when the transfer to 999 occurs before computations begin */
/*  I think that newbas and degen also have this problem      */

        alpha = _info->tols.plinfy;
        newbas = 0;
        degen = 0;

        _info->nph0.ninf0 = 0;
        _info->iters.nph0it = 0;
        _info->cbmode.phase0 = TRUE;
        *infes = FALSE;
        _info->nph0.ndrop = 0;

/*--------------------------------------------------------------------*/
/*  ** 7/26/03 jcp ** add an interation limit to prevent runaway      */
/*   execution.  **fixme** add time limit when I figure out what a    */
/*   reasonable approach would be and add user options for both       */
/*--------------------------------------------------------------------*/
         maxIterations = _info->limser / 10;
         if(maxIterations < 100) maxIterations = 100;
/*--------------------------------------------------------------------*/


        /*      --------------------------------------------------------------
         *     |  ICSTAT is used to parition the contraints into 3 classes:
         *     |    ICSTAT(I)= 1 --> Constraint "I" is infeasible or binding
         *     |    ICSTAT(I)= 0 --> Constraint "I" is feasible and *not* binding
         *     |    ICSTAT(I)=-1 --> Constraint "I" is dropped.
         *      --------------------------------------------------------------
         * */
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                icstat[i_] = 0;
        }

        /*      ------------------------------------------------------------
         *     |  Outer most loop -- Come here to restart entire process
         *      ------------------------------------------------------------
         * */
L_1:
        nloops = nloops + 1;
        nsmall = 0;
        ndegen = 0;
        sinfo = 0.0e0;
        sinfdo = 0.0e0;
        errmxo = 0.0e0;
        lv = 0;
        lvv = 0;
        lvbnd = 0;

        /*      --------------------------------
         *     |  Set the slack variables
         *      --------------------------------
         * */
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                if( g[i_] < alb[_info->dimen.n + i_] ){
                        x[_info->dimen.n + i_] = alb[_info->dimen.n + i_];
                } else if( g[i_] > ub[_info->dimen.n + i_] ){
                        x[_info->dimen.n + i_] = ub[_info->dimen.n + i_];
                } else{
                        x[_info->dimen.n + i_] = g[i_];
                }
        }

        /*      ------------------------------------------
         *     |  Compute the initial residuals and norms
         *      ------------------------------------------
         * */
        sumres(_info, g, x, alb, ub, error, icstat, &errmx, &sinf, &sinfd, &ninfa,
         &ninfd );
        if( nloops == 1 )
                _info->nph0.ninf0 = ninfa + ninfd;

        /*      ----------------------------
         *     | IF FEASIBLE, THEN RETURN
         *      ----------------------------
         * */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 2 )
                {
                sprintf(LSGRG_MSGBUFFER, " SELECTING ORIGINAL BASIS --- NUM IGNORED CONSTS = %5ld\n NINF = %5ld SINF =  %14.7e MAXERR= %14.7e\n",
                 _info->nph0.ndrop, ninfa, sinf, errmx );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( ninfa == 0 )
                goto L_999;

        /*      ------------------------------------------------------
         *     | Call CONSBS to select and invert the original basis
         *     | and to determine the feasible/infeasible constraints
         *      ------------------------------------------------------
         *
         *     write(ioout,*) ' PHASE0 - nzgrad = ',nzgrad */
        consbs(_info, ibmap, zmem, grad, ihag, iheg, ihegl, ibv, x, alb, ub,
         inbv, iub, ibc, g, inlin, ivstat, istat, x0, gg, pvrow, icols,
         iwork, err0, icand, cdnum, irank, ascale, ibvbct, &nbfree, iobjpt,
         paij, iprow, ipcol, ipmap );
        if( _info->glberr.abort )
                return;
        newbas = TRUE;


        /*----------------------------------------------------------------------
         *     INBV holds the current set of Non-Basic variables
         *     ICOLS holds index of variables in either INBV or IBV
         *         ICOLS(I) = -K  --->  IBV(K) = I
         *         ICOLS(I) =  K   ---> INBV(K) =I
         *     Set status for each variable in INBV: (Set in CONSBS)
         *         IUB(J) = -1    --->  VAR J AT LOWER BOUND
         *         IUB(J) =  0    --->  VAR J FREE
         *         IUB(J) =  1    --->  VAR J AT UPPER BOUND
         *----------------------------------------------------------------------
         * */
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                icols[ibv[i_] - 1] = -i;
        }
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                icols[inbv[i_] - 1] = i;
        }
        for( i = 1; i <= _info->bind.nbc; i++ ){
                i_ = i - 1;
                if( icstat[ibc[i_] - 1] >= 0 )
                        icstat[ibc[i_] - 1] = 1;
        }
        for( i = _info->bind.nbc + 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                if( icstat[ibc[i_] - 1] >= 0 )
                        icstat[ibc[i_] - 1] = 0;
        }

        /*      -----------------------------------------------------
         *     | LOOP -- Come here to compute a new Newton correction
         *      ----------------------------------------------------- */
L_2:
        ;

        /*C     --------------------------------------------------
         *C    | Compute Newton correction.  Needed to check for
         *c    | for critical basic variables.
         *C     --------------------------------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                delx[i_] = error[i_];
        }
        xftran(_info, ibmap, zmem, delx, _info->nintbk.nb );

        /*      ----------------------------------------
         *     | Check for a blocking basic variable
         *      ----------------------------------------
         * */
L_4:
        ;
        if( _info->nph0.ndrop == _info->dimen.m )
                goto L_999;
        if( lv == 0 )
                goto L_3;
        if( delx[lv - 1]*lvbnd >= 0 )
                goto L_3;
        if( icstat[lv - 1] == -1 && ibv[lv - 1] == _info->dimen.n + lv )
                goto L_3;
        if( icstat[lv - 1] == 0 && ibv[lv - 1] == _info->dimen.n + lv ){
                icstat[lv - 1] = 1;
#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER,
                         " REPLACING SLACK FOR NEW ACTIVE CONSTRAINT %5ld\n",
                         lv );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
        }

        /*      ----------------------------------
         *     |  Perform a basis change pivot
         *      ----------------------------------
         * */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, "  REPLACING BLOCKING VARIABLE # %5ld  IN ROW # %5ld\n",
                 ibv[lv - 1], lv );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        ph0piv(_info, grad, ihag, iheg, ibmap, zmem, err0, pvrow, inbv, alb,
         ub, x, ibv, icols, iub, ivstat, cdnum, irank, lv, delx, infes,
         icstat, &singlr );
        if( _info->glberr.abort )
                return;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, "  PIVOTING COMPLETED IN ROW %5ld NEW BASIC VARIABLE = %4ld\n",
                 lv, ibv[lv - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( icstat[lv - 1] == 1 && ibv[lv - 1] == _info->dimen.n + lv ){
                icstat[lv - 1] = 0;

#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 3 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " CONSTRAINT # %5ld HAS BECOME INACTIVE\n",
                         lv );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
        }
        if( singlr )
                goto L_1;
        if( *infes ){
                if( icstat[lv - 1] >= 0 ){

#ifdef IO_ENABLED
                        if( _info->nintbk.ipr >= 2 )
                                {
                                sprintf(LSGRG_MSGBUFFER, "  NO BASIS CHANGE POSSIBLE FOR ROW -- %4ldBASIC VARIABLE = %4ld DROPPING THIS ROW\n",
                                 lv, ibv[lv - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                        droprw(_info, lv, icstat, istat, &_info->nph0.ndrop, &ninfa, &ninfd,
                         &sinf, &sinfd, error, g, x, alb, ub, &errmx );
                        if( _info->nph0.ndrop == _info->dimen.m )
                                goto L_999;
                        goto L_4;       /*  SHOULD WE GO TO 1 INSTEAD??????? */
                } else{
#ifdef IO_ENABLED
                        sprintf(LSGRG_MSGBUFFER, "  NO BASIS CHANGE POSSIBLE FOR DROPPED ROW -- %4ldBASIC VARIABLE = %4ld PHASE-0 TERMINATED \n",
                         lv, ibv[lv - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                        goto L_999;
                }
        } else{
                newbas = TRUE;
        }

        /*C     --------------------------------------------------
         *C    | Compute Newton correction for new basis
         *C     --------------------------------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                delx[i_] = error[i_];
        }
        xftran(_info, ibmap, zmem, delx, _info->nintbk.nb );

        /*      --------------------------------------------
         *     | Compute step length to first basic bound
         *      --------------------------------------------
         * */
L_3:
        alpha = _info->tols.plinfy;
        lv = 0;
        lvv = 0;
        lvbnd = 0;
        for( j = 1; j <= _info->nintbk.nb; j++ ){
                j_ = j - 1;
                t = delx[j_];
                if( fabs( t ) <= _info->tols.eps )
                        goto L_110;
                k = ibv[j_];
                if( (k == _info->dimen.n + j) && (icstat[j_] == -1) )
                        goto L_110;
                if( t < 0.0e0 )
                        goto L_115;
                /*XXXX     DIST=X(K)-ALB(K) + EPBOUN*(1.0D0+DABS(ALB(K))) */
                dist = x[k - 1] - alb[k - 1];
                jbnd = -1;
                goto L_116;
                /*115      DIST=X(K)-UB(K) - EPBOUN*(1.0D0+DABS(UB(K))) */
L_115:
                dist = x[k - 1] - ub[k - 1];
                jbnd = 1;
L_116:
                ;
                if( fabs( dist ) > 1.0e20 )
                        dist = sign( 1.0e20, dist );
                t = dist/t;
                if( alpha <= t )
                        goto L_110;
                alpha = t;
                lv = j;
                lvbnd = jbnd;
L_110:
                ;
        }
        _info->iters.nph0it = _info->iters.nph0it + 1;

/*-----------------------------------------------------------------*/
/*  ** 7/26/03 jcp ** add test for iteration limit exceeded and msg*/
/*-----------------------------------------------------------------*/
        if(_info->iters.nph0it > maxIterations) {
#ifdef IO_ENABLED
           if( _info->nintbk.ipr > 0 )
                {
                sprintf(LSGRG_MSGBUFFER,
                 "\n Phase0 Iteration Limit of %5d Exceeded"
                 " -- Phase0 Terminated",maxIterations);
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
            goto L_999;
         }
/*-----------------------------------------------------------------*/

        /*      -----------------------------------
         *     | Check here for degenerate steps
         *      -----------------------------------
         * */
        degen = alpha <= STPMIN;
        if( degen ){
                if( alpha <= _info->tols.eps )
                        alpha = 0.0e0;

#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 3 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " DEGENERATE STEP -- LV = %4ld BASIC VAR = %4ld STEP = %13.6e\n",
                         lv, ibv[lv - 1], alpha );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                ndegen = ndegen + 1;
                if( ndegen >= MXDEG ){

#ifdef IO_ENABLED
                        if( _info->nintbk.ipr >= 2 )
                                {
                                sprintf(LSGRG_MSGBUFFER, " DEGENERATE FOR %5ld CONSECUTIVE STEPS --  DROPPING CONSTRAINT # %5ld\n",
                                 ndegen, lv );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                        ndegen = 0;
                        droprw(_info, lv, icstat, istat, &_info->nph0.ndrop, &ninfa, &ninfd,
                         &sinf, &sinfd, error, g, x, alb, ub, &errmx );
                        if( _info->nph0.ndrop == _info->dimen.m )
                                goto L_999;
                        if( ibv[lv - 1] == _info->dimen.n + lv )
                                goto L_3;
                }
                ph0log(_info, _info->dimen.n, _info->dimen.mp1, _info->iters.nph0it, _info->nph0.ndrop, _info->nintbk.nb,
                 ninfa, &_info->nintbk.ipr, &_info->nintbk.ipr3, linect, newbas, degen,
                 sinf, sinfd, errmx, alpha, lv );
                goto L_4;
        } else{
                ndegen = 0;
                if( alpha >= 1.0e0 ){
                        alpha = 1.0e0;
                        lv = 0;
                }
        }
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 6 )
                {
                sprintf(LSGRG_MSGBUFFER, " PHASE0 ITER# %5ld STEP= %14.7e LV=%5ld\n",
                 _info->iters.nph0it, alpha, lv );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
        if( _info->nintbk.ipr >= 6 )
                {
                sprintf(LSGRG_MSGBUFFER, " PHASE0 ITER # %5ld", _info->iters.nph0it );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                lsgrg_msg(IOINFO, " -- NEWTON CORRECTION DELX:\n " );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", delx[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
                }
#endif

        /*      ----------------------------------------------------
         *     |  Save the error and the base point in case we
         *     |  have to reduce the step size.
         *      ----------------------------------------------------
         * */
        nslow = 0;
        for( i = 1; i <= _info->dimen.npmp1; i++ ){
                i_ = i - 1;
                x0[i_] = x[i_];
        }

        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                gg[i_] = g[i_];
                err0[i_] = error[i_];
        }

        sinfo = sinf;
        sinfdo = sinfd;
        errmxo = errmx;
        ninfdo = ninfd;
        ninfao = ninfa;

        /*      -------------------------------------------------------------
         *     | Take a step of length ALPHA.  Come back here to reduce step
         *      -------------------------------------------------------------
         * */
L_333:
        ;
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                k = ibv[i_];
                x[k - 1] = x0[k - 1] - alpha*delx[i_];
        }

        /*      ---------------------------------------------
         *     |  Compute the new residuals
         *      ---------------------------------------------
         * */
        calfun(_info, g, x, ascale );

        /*      --------------------------------------------------------------
         *     | Compute norm of constraint error and new error for constraints
         *      ---------------------------------------------------------------
         * */
        sumres(_info, g, x, alb, ub, error, icstat, &errmx, &sinf, &sinfd, &ninfa,
         &ninfd );

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 6 )
                {
                sprintf(LSGRG_MSGBUFFER, "  AFTER NEW RESIDUALS--- SINF= %15.7e ERRMX = %15.7e\n  SINFD= %15.7e TOTAL = %14.6e\n",
                 sinf, errmx, sinfd, sinf + sinfd );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( ninfa == 0 )
                goto L_999;


        /*      ------------------------------------
         *     | Check for convergence problems
         *      ------------------------------------
         *
         *      ------------------------------------------------------------
         *     |  Drop slow convergent constraints only if basis is current
         *      ------------------------------------------------------------
         * */
        bschng = FALSE;
        slowc = sinf >= (1.0e0 - DELTAF*alpha)*sinfo || (nsmall >= MXSMAL);

        if( (sinf >= (1.0e0 - DELTAM*alpha)*sinfo) && newbas ){
                /*      -----------------------------------------------------------
                 *     | Problem -- The residuals increased and Jacobian is current.
                 *     |   Reduce step size if possible.  If step is already too
                 *     |   small then drop a constraint from consideration.
                 *      ----------------------------------------------------------- */

#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 2 )
                        {
                        sprintf(LSGRG_MSGBUFFER, "  RESIDUALS INCREASED FOR CURRENT STEP --   SINF= %14.7e SINFO= %14.7e\n",
                         sinf, sinfo );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                nslow = nslow + 1;

                if( alpha > STPMIN ){
                        alpha = alpha*REDSTP;

#ifdef IO_ENABLED
                        if( _info->nintbk.ipr >= 2 )
                                {
                                sprintf(LSGRG_MSGBUFFER, "  REDUCING STEP SIZE ALPHA TO %14.6e\n",
                                 alpha );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                        lv = 0;
                        goto L_333;
                } else{
                        /*         ----------------------------------------------------------
                         *        | Drop the infeasible constraint with the largest
                         *        | absolute increase in its residual.  Make the appropriate
                         *        | slack the new basic variable for that row.
                         *         ---------------------------------------------------------- */
                        derr = 0.0e0;
                        aerr = 0.0e0;
                        jmax = 0;
                        jmaxa = 0;
                        for( i = 1; i <= _info->dimen.mp1; i++ ){
                                i_ = i - 1;
                                if( icstat[i_] < 0 )
                                        goto L_404;
                                err = fabs( error[i_] );
                                if( err <= _info->limits.epnewt*(1.0e0 + fabs( x[_info->dimen.n + i_] )) )
                                        goto L_404;
                                if( err > aerr ){
                                        aerr = err;
                                        jmaxa = i;
                                }
                                dist = err - fabs( err0[i_] );
                                if( dist <= 0 )
                                        goto L_404;
                                if( dist > derr ){
                                        derr = dist;
                                        jmax = i;
                                }
L_404:
                                ;
                        }
                        if( derr <= _info->limits.epnewt )
                                jmax = jmaxa;
                        if( jmax > 0 ){
                                droprw(_info, jmax, icstat, istat, &_info->nph0.ndrop, &ninfa,
                                 &ninfd, &sinf, &sinfd, error, g, x, alb, ub, &errmx );
                                if( _info->nph0.ndrop == _info->dimen.m )
                                        goto L_999;
                                sumres(_info, gg, x0, alb, ub, err0, icstat, &errmxo, &sinfo,
                                 &sinfdo, &ninfao, &ninfdo );

#ifdef IO_ENABLED
                                if( _info->nintbk.ipr >= 1 )
                                        {
                                        sprintf(LSGRG_MSGBUFFER, " PHASE0 ITER #: %5ld  CONSTRAINT %5ld DROPPED -- BASIC VAR IS%5ld\n",
                                         _info->iters.nph0it, jmax, ibv[jmax - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                        }
#endif
                        }
                }
        } else if( slowc && newbas ){
                /*      -----------------------------------------------------------
                 *     |  Residuals decreased, but slowly.  Make a basis change to
                 *     |  try to increase the rate of convergence by removing the
                 *     |  non-slack basic variable that changes the most.
                 *      ----------------------------------------------------------- */
                bschng = TRUE;

                dmaxx = 0.e0;
                lvv = 0;
                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        i_ = i - 1;
                        k = ibv[i_];
                        if( k <= _info->dimen.n && icstat[i_] >= 0 ){
                                dxx = fabs( delx[i_] );
                                if( dxx > dmaxx ){
                                        dmaxx = dxx;
                                        lvv = i;
                                }
                        }
                }

        } else if( (_info->nph0.ndrop > 0 && sinfd > sinfdo) && newbas ){
                /*      -----------------------------------------------------------
                 *     |  The residuals of the dropped constraints increased.  Make
                 *     |  a basis change to try to avoid this.  Remove a basic
                 *     |  variable by looking at the change in dropped slacks
                 *     |  relative to each basic variable.
                 *      -----------------------------------------------------------
                 * */
                dmaxx = 0.e0;
                lvv = 0;
                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        i_ = i - 1;
                        k = ibv[i_];
                        if( k <= _info->dimen.n && icstat[i_] >= 0 ){
                                dxx = 0.e0;
                                for( j = iheg[k - 1]; j <= (iheg[k] - 1); j++ ){
                                        j_ = j - 1;
                                        if( icstat[ihag[j_] - 1] < 0 )
                                                dxx = dxx + fabs( grad[j_] );
                                }
                                dxx = dxx*alpha*fabs( delx[i_] );
                                if( dxx > dmaxx ){
                                        dmaxx = dxx;
                                        lvv = i;
                                }
                        }
                }
                bschng = TRUE;
        }

        /*      --------------------------------------------------
         *     | If the new point is worse, restore the old point
         *      --------------------------------------------------
         * */
        if( sinf > sinfo ){

                if( _info->nintbk.ipr >= 2 )
                        {
                        lsgrg_msg(IOINFO,
                         " SINF INCREASED --- RESTORING ORIGINAL POINT \n" );
                        }
                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        i_ = i - 1;
                        g[i_] = gg[i_];
                        error[i_] = err0[i_];
                }
                for( i = 1; i <= _info->dimen.npmp1; i++ ){
                        i_ = i - 1;
                        x[i_] = x0[i_];
                }
                sinf = sinfo;
                sinfd = sinfdo;
                errmx = errmxo;
                ninfa = ninfao;
                ninfd = ninfdo;
        }
        ph0log(_info, _info->dimen.n, _info->dimen.mp1, _info->iters.nph0it, _info->nph0.ndrop, _info->nintbk.nb,
         ninfa, &_info->nintbk.ipr, &_info->nintbk.ipr3, linect, newbas, degen, sinf,
         sinfd, errmx, alpha, lv );

        /*      ----------------------------------------------------
         *     | New basis inversion if convergence was slow
         *     |     or if a variable hit a bound.
         *      ----------------------------------------------------
         * */
        if( alpha <= ASMALL ){
                nsmall = nsmall + 1;
        } else{
                nsmall = 0;
        }
        newbas = ((bschng || slowc) || lv != 0) || nslow > 0;
        if( newbas && !_info->equblk.lincon ){
                ph0fac(_info, x, grad, ihag, iheg, ihegl, inlin, ibmap, zmem, ibv,
                 g, gg, (double*)iwork, err0, ub, alb, ascale, _info->ztol.bigelt,
                 cdnum, irank, FALSE, &singlr, iobjpt, paij, iprow, ipcol,
                 ipmap );
                if( _info->glberr.abort )
                        return;
                if( singlr )
                        goto L_1;
        }

        /*      ---------------------------------
         *     | Check for forced basis change
         *      ---------------------------------
         * */
        if( lvv > 0 ){

#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 2 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " CONVERGENCE PROBLEMS -- FORCING A BASIS CHANGE IN ROW %4ldOLD BASIC VARIABLE = %4ld\n",
                         lvv, ibv[lvv - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                ph0piv(_info, grad, ihag, iheg, ibmap, zmem, err0, pvrow, inbv,
                 alb, ub, x, ibv, icols, iub, ivstat, cdnum, irank, lvv, delx,
                 infes, icstat, &singlr );
                if( _info->glberr.abort )
                        return;
                if( singlr )
                        goto L_1;

#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 3 )
                        {
                        sprintf(LSGRG_MSGBUFFER, "  PIVOTING COMPLETED IN ROW %5ld NEW BASIC VARIABLE = %4ld\n",
                         lvv, ibv[lvv - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( icstat[lvv - 1] == 1 && ibv[lvv - 1] == _info->dimen.n + lvv ){
                        icstat[lvv - 1] = 0;

#ifdef IO_ENABLED
                        if( _info->nintbk.ipr >= 5 )
                                {
                                sprintf(LSGRG_MSGBUFFER, " CONSTRAINT # %5ld HAS BECOME INACTIVE\n",
                                 lvv );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                }
#endif
                }

#ifdef IO_ENABLED
                if( *infes && _info->nintbk.ipr >= 5 )
                        {
                        sprintf(LSGRG_MSGBUFFER, "  NO BASIS CHANGE POSSIBLE FOR ROW -- %4ldKEEPING OLD BASIC VARIABLE = %4ld\n",
                         lvv, ibv[lvv - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
        }

        /*      -----------------------------------------
         *     |  Go back and perform a new iteration
         *      -----------------------------------------
         * */
        goto L_2;

        /*      -------------------------------
         *     | Come here when converged
         *      ------------------------------- */
L_999:
        ;
        ph0log(_info, _info->dimen.n, _info->dimen.mp1, _info->iters.nph0it, _info->nph0.ndrop, _info->nintbk.nb,
         ninfa, &_info->nintbk.ipr, &_info->nintbk.ipr3, linect, newbas, degen, sinf,
         sinfd, errmx, alpha, lv );


#ifdef IO_ENABLED
        if( _info->nintbk.ipr == 1 ) {
            lsgrg_msg(IOINFO,"\n Phase0 Iterations Completed");
             if( _info->io_info.screen_output_enabled)
                lsgrg_screen_msg(IOINFO,"\n Phase0 Iterations Completed");
        }
        if( _info->nintbk.ipr >= 2 )
                {
                sprintf(LSGRG_MSGBUFFER, " PHASE0 COMPLETED .... NEWTON ITERS: %5ld\n NUM CONSTRAINTS DROPPED: %5ld NUM CONSTRAINTS CONVERGED: %5ld\n",
                 _info->iters.nph0it, _info->nph0.ndrop, _info->dimen.m - _info->nph0.ndrop );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
        if( _info->nintbk.ipr3 >= 4 )
                {
                lsgrg_msg(IOINFO, "\n     FINAL X = \n " );
                for( i = 1; i <= _info->dimen.npmp1; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", x[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
                }
#endif
        *infes = FALSE;
        _info->cbmode.newpt = TRUE;
        _info->cbmode.compgr = TRUE;
        _info->cbmode.phase0 = FALSE;
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                istat[i_] = labs( istat[i_] );
        }
/*      cputime( &tfinal, &irc ); */
        tfinal = lsgrg_timer();

        _info->ph0tim.tph0 = tfinal - tstart;
        return;


        /*C     End of PHASE0
         *C */
} /* end of function */



/*    *************************************************************
 *    *             Produces an LU factorization of the current
 *    *             basis given in IBV.  "GETGR" is called to
 *    *             compute the Jacobian and then the factorization
 *    *             is performed.
 *    ************************************************************* */

void /*FUNCTION*/ LSGRGCLASS ph0fac(LsgrgInfo *_info,double x[], double grad[], long int ihag[],
         long int iheg[], long int ihegl[], long int inlin[], long int ibmap[],
         double zmem[], long int ibv[], double g[], double gplus[], double gminus[],
         double gcol[], double ub[], double alb[], double ascale[], double bigelt,
         double cdnum[], long int irank[], LOGICAL32 bschng, LOGICAL32 *singlr,
         long int iobjpt[], double paij[], long int iprow[], long int ipcol[],
         long int ipmap[])
{
LOGICAL32 jreset;
long int  rcode;
double ztemp, ztime;



        /*......................................................................
         *
         *     *** START OF PH0FAC ***
         *
         *      ------------------------------------------
         *     |  Compute Jacobian at the current point
         *      ------------------------------------------
         * */
        jreset = FALSE;
        getgr(_info, x, grad, ihag, iheg, ihegl, inlin, g, gplus, gminus, gcol,
         ub, alb, ascale, &jreset, iobjpt, paij, iprow, ipcol, ipmap );

        if( bigelt < _info->limits.epspiv ){

#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 2 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " PH0FAC...ERROR:  LARGEST ELEMENT IN GRAD ARRAY IS %16.8e\n  A NONSINGULAR BASIS MATRIX CAN NOT BE CHOSEN--USING AN ALL SLACK BASIS\n",
                         bigelt );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
#ifdef IO_ENABLED
                sprintf( LSGRG_MSGBUFFER,
                 " PH0FAC...ERROR:  LARGEST ELEMENT IN GRAD ARRAY IS"
                 " %16.8e\n  A NONSINGULAR BASIS MATRIX CAN NOT BE"
                 " CHOSEN--USING AN ALL SLACK BASIS\n",
                 bigelt );
                 lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
        }

        /*      --------------------------------------------
         *     | Perform L/U factorization of basis matrix
         *      --------------------------------------------
         * */
/*      cputime( &ztime, &irf ); */
        ztime = lsgrg_timer();

        xfact(_info, grad, ihag, iheg, ibmap, zmem, ibv, bschng, &rcode, cdnum,
         irank, TRUE );
/*      cputime( &ztemp, &irf ); */
        ztemp = lsgrg_timer();

        _info->temp.tfact = _info->temp.tfact + ztemp - ztime;


        /*      ------------------------------------------
         *     : Check the return code from factorization
         *      ------------------------------------------
         * */
        *singlr = FALSE;
        if( rcode == 2 ){

#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 2 )
                        {
                        sprintf(LSGRG_MSGBUFFER, " PH0FAC...ERROR: CURRENT BASIS MATRIX IS SINGULAR  RCODE FROM XFACT = %3ld\n",
                         rcode );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                /*XXXX     NSING = NSING + 1 */
                *singlr = TRUE;
                return;
        } else if( rcode == 5 ){
                lsgrg_error_msg(IOINFO,
                  " PH0FAC...CURRENT BASIS IS ILL CONDITIONED --"
                  " PROCEEDING WITH NEWTON METHOD \n" );
        } else if( rcode != 1 ){
#ifdef IO_ENABLED
                sprintf( LSGRG_MSGBUFFER,
                  " PH0FAC...ERROR: COULD NOT INVERT CURRENT BASIS RCODE"
                  " FROM XFACT = %3ld\n",
                  rcode );
                lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
      /*        xerror( 0, iounit.ioerr ); */
                lsgrg_errorexit(JUMPBUF,_LSGRG_PH0FAC_INVERT_FAILURE);
                return;
        }
        return;
} /* end of function */



void /*FUNCTION*/ LSGRGCLASS ph0piv(LsgrgInfo *_info,double grad[], long int ihag[], long int iheg[],
         long int ibmap[], double zmem[], double v[], double d[], long int inbv[],
         double alb[], double ub[], double x[], long int ibv[], long int icols[],
         long int iub[], long int ivstat[], double cdnum[], long int irank[],
         long int lv, double delx[], LOGICAL32 *infes, long int icstat[],
         LOGICAL32 *singlr)
{
long int i, i_, icol, ii, in, iprp, j, j_, jmax, jq, jq2, k, rcode;
double d1, delxlv, dxj, pivmax, pvt, smax_, step, sum, tol, xj;

        iprp = _info->nintbk.ipr;
        if( iprp >= 5 )
                lsgrg_msg(IOINFO," PH0PIV ENTERED \n" );

#ifdef IO_ENABLED
        if( iprp >= 5 ){
                lsgrg_msg(IOINFO, " PH0PIV ENTERED \n" );
                sprintf(LSGRG_MSGBUFFER, " LEAVING VARIABLE IS BASIC VARIABLE NO.  %5ld VAR# %5ld\n",
                 lv, ibv[lv - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                lsgrg_msg(IOINFO, " INBV IS" );
                for( i = 1; i <= _info->dimen.n; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", inbv[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
                lsgrg_msg(IOINFO, " IUB IS " );
                for( i = 1; i <= _info->dimen.n; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", iub[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
                lsgrg_msg(IOINFO, " IBV IS " );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%4ld", ibv[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
        }
#endif
        *singlr = FALSE;

        /*C     -----------------------------------------------------
         *C    | Compute pivot row of leaving variable in BINV
         *     | V(1),...,V(NB) will hold the row
         *C     -----------------------------------------------------
         * */
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                v[i_] = 0.0e0;
        }
        v[lv - 1] = 1.0e0;
        xbtran(_info, ibmap, zmem, v, _info->nintbk.nb );

#ifdef IO_ENABLED
        if( iprp >= 6 )
                {
                lsgrg_msg(IOINFO, " PIVOT ROW OF BINV....\n " );
                for( i = 1; i <= _info->nintbk.nb; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", v[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
                }
#endif

        /*C     -----------------------------------------------------
         *C    | Compute updated pivot row
         *     | D(1),...,D(N) will hold the pivot choices
         *C     -----------------------------------------------------
         * */
        jmax = 0;
        pivmax = 0.0e0;
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                ii = inbv[i_];
                d[i_] = 0.0e0;
                if( ii > _info->dimen.n )
                        goto L_47;
                if( ivstat[ii - 1] < 0 )
                        goto L_50;
L_47:
                xdot(_info, grad, ihag, iheg, v, _info->nintbk.nb, ii, &sum );
                d[i_] = sum;
                if( pivmax >= fabs( sum ) )
                        goto L_50;
                pivmax = fabs( sum );
                jmax = i;
L_50:
                ;
        }

#ifdef IO_ENABLED
        if( iprp >= 5 ){
                lsgrg_msg(IOINFO, " UPDATED PIVOT ROW OF NONBASICS .... \n " );
                for( i = 1; i <= _info->dimen.n; i++ ){
                        sprintf(LSGRG_MSGBUFFER, "%13.6e", d[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, "\n" );
                sprintf(LSGRG_MSGBUFFER, " PH0PIV: JMAX = %5ld LARGEST PIVOT = %13.6e\n",
                 jmax, pivmax );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        /*      -------------------------------------------------------
         *     | If all pivot elements too small,  return INFEASIBLE
         *      -------------------------------------------------------
         * */
        if( pivmax < _info->limits.epspiv )
                goto L_67;

        /*C     --------------------------------------------------
         *C    | Choose the entering variable with largest step
         *C     --------------------------------------------------
         * */
        if( icstat[lv - 1] == -1 ){
                in = icols[lv + _info->dimen.n - 1];
                if( in <= 0 ){
#ifdef IO_ENABLED
                        sprintf( LSGRG_MSGBUFFER,
                          "  PH0PIV...ERROR -- INDEX ERROR IN ICOLS FOR"
                          " SLACK #%5ld VAR # %5ld IS BASIC IN ROW %5ld\n",
                          lv, lv + _info->dimen.n, -in );
                        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
            /*          xerror( 0, iounit.ioout ); */
                        lsgrg_errorexit(JUMPBUF,_LSGRG_PH0PIV_BAD_INDEX);
                        return;
                }
#ifdef IO_ENABLED
                if(iprp > 0)
                   sprintf(LSGRG_MSGBUFFER,
                   " PH0PIV.. CONS# %5ld DROPPED -- MAKING SLACK %5ldBASIC\n",
                    lv, lv + _info->dimen.n );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                jq = in;
                goto L_65;
        }

        /*     TOL=DMAX1(XPVPCT*PIVMAX,EPSPIV) */
        tol = _info->limits.epspiv;
        delxlv = delx[lv - 1];

/* L_152: unreferenced??? */
        smax_ = 0.0e0;
        pivmax = 0.0e0;
        pvt = 0.0e0;

#ifdef IO_ENABLED
        if( iprp > 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " CHECKING SUPERBASIC PIVOTS -- TOL = %14.7e\n",
                 tol );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        jq = 0;
        for( j = 1; j <= _info->dimen.n; j++ ){
                j_ = j - 1;

#ifdef IO_ENABLED
                if( _info->dbug.debug )
                        {
                        sprintf(LSGRG_MSGBUFFER, " SUPERBASIC NO.  %4ld    ROW ELEM = %15.8e\n",
                         j, d[j_] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( fabs( d[j_] ) < tol )
                        goto L_60;
                k = inbv[j_];
                if( k > _info->dimen.n && k != lv + _info->dimen.n )
                        goto L_60;
                xj = x[k - 1];
                dxj = delxlv/d[j_];
                d1 = 0.0e0;
                if( dxj < 0 ){
                        if( iub[j_] != 1 )
                                d1 = x[k - 1] - ub[k - 1] - _info->limits.epboun*(1.0e0 +
                                 fabs( ub[k - 1] ));
                } else{
                        if( iub[j_] != -1 )
                                d1 = x[k - 1] - alb[k - 1] + _info->limits.epboun*(1.0e0 +
                                 fabs( alb[k - 1] ));
                }
                step = d1/dxj;

#ifdef IO_ENABLED
                if( _info->dbug.debug )
                        {
                        sprintf(LSGRG_MSGBUFFER, " VARIABLE NO.    %4ld    DBOUND =   %15.8e\n",
                         k, d1 );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                if( step >= 1.0e0 ){
                        smax_ = 1.0e0;
                        pvt = fabs( d[j_] );
                        if( pvt > pivmax ){
                                pivmax = pvt;
                                jq = j;
                        }
                } else if( step > smax_ ){
                        smax_ = step;
                        jq = j;
                }
L_60:
                ;
        }

        /*     -------------------------------------------------------
         *    | Could not find a suitable pivot, so return INFEASIBLE
         *     -------------------------------------------------------
         * */
        if( jq == 0 )
                goto L_67;


        /*C     --------------------------------------------------------
         *C    |         Pivot new variable into basis
         *C    | "ICOL" enters basis and "JQ2" leaves in row "LV"
         *C     --------------------------------------------------------
         * */
L_65:
        icol = inbv[jq - 1];
        jq2 = ibv[lv - 1];

#ifdef IO_ENABLED
        if( iprp >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " VARIABLE%4ld LEAVING BASIS  - BASIC VAR NO.  %4ld\n VARIABLE%4ld ENTERING BASIS - NON-BASIC NO. %4ld\n",
                 jq2, lv, icol, jq );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

        xpivot(_info, grad, ihag, iheg, ibmap, zmem, ibv, icol, lv, v, &rcode,
         cdnum, irank, TRUE );

        if( rcode == 1 )
                goto L_163;

#ifdef IO_ENABLED
        sprintf( LSGRG_MSGBUFFER,
           " PH0PIV...ERROR: RETURN FROM XPIVOT - RCODE = %5ld\n",
             rcode );
        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
        if( rcode == 5 ){
             lsgrg_error_msg(IOINFO,
               " PH0PIV...ERROR: UPDATED BASIS MATRIX IS  ILL-CONDITIONED\n" );
        } else{
                if( rcode == 2 ){
                    lsgrg_error_msg(IOINFO,
                       " PH0PIV...ERROR: UPDATED BASIS MATRIX IS SINGULAR \n" );
                    *singlr = TRUE;
                } else{
#ifdef IO_ENABLED
                     sprintf( LSGRG_MSGBUFFER,
                       " PH0PIV...ERROR: RETURN FROM XPIVOT - RCODE = %5ld\n",
                          rcode );
                     lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
              /*        xerror( 0, iounit.ioerr ); */
                lsgrg_errorexit(JUMPBUF,_LSGRG_PH0PIV_XPIVOT_FAILURE);
                }
                return;
        }

        /*      ------------------------------------------------------------
         *     | Update index sets of basic and nonbasic variables and IUB
         *      ------------------------------------------------------------
         *
         *C */
L_163:
        ;
        if( icols[jq2 - 1] != -lv ){
#ifdef IO_ENABLED
             sprintf( LSGRG_MSGBUFFER,
               " PH0PIV...ERROR: ERROR IN ICOLS INDEX WHEN PIVOTING IN"
               " ROW %5ld\n                 LEAVING VARIABLE IS %5ld AND"
               " ICOLS(%5ld) IS %5ld\n",
                   lv, jq2, jq2, icols[jq2 - 1] );
             lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
          /*    xerror( 0, iounit.ioerr ); */
                lsgrg_errorexit(JUMPBUF,_LSGRG_PH0PIV_BAD_ICOLS1);
                return;
        }
        if( icols[icol - 1] != jq ){
#ifdef IO_ENABLED
            sprintf( LSGRG_MSGBUFFER,
              " PH0PIV...ERROR: ERROR IN ICOLS INDEX FOR NONBASIC # %5ld\n                 NONBASIC VAR IS %5ld BUT ICOLS(%5ld) IS%5ld\n",
                 jq, icol, icol, icols[icol - 1] );
            lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
        /*      xerror( 0, iounit.ioerr ); */
                lsgrg_errorexit(JUMPBUF,_LSGRG_PH0PIV_BAD_ICOLS2);
                return;
        }
        icols[icol - 1] = -lv;
        icols[jq2 - 1] = jq;
        inbv[jq - 1] = jq2;
        goto L_150;


        /*----------------------------------------------------------------
         *    No pivot > EPSPIV could be found -- so STOP
         *----------------------------------------------------------------
         * */
L_67:
#ifdef IO_ENABLED
        sprintf( LSGRG_MSGBUFFER,
          " PH0PIV...ERROR -- COULD NOT FIND A PIVOT > EPSPIV IN ROW"
          " %4ld  WHEN ATTEMPTING TO PIVOT OUT VARIABLE %5ld\n",
             lv, ibv[lv - 1] );
        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
        sprintf( LSGRG_MSGBUFFER,
         " PH0PIV ... MAXIMUM PIVOT ELEMENT = %14.7e IS LESS THAN EPSPIV = \n",
            pivmax );
        lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif

        *infes = TRUE;

L_150:
        ;

        if( iprp >= 5 )
                {
                lsgrg_msg(IOINFO, " PH0PIV COMPLETED \n" );
                }
        return;



        /*     End of PH0PIV
         * */
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS sumres(LsgrgInfo *_info,double g[], double x[], double alb[], double ub[],
         double error[], long int icstat[], double *errmx, double *sinf,
         double *sinfd, long int *ninfa, long int *ninfd)
{
long int i, i_, isl;
double abserr, err, relerr, xval;


        /* ....................................................................
         *
         *      -----------------------------------------------------------
         *     |  Compute the residuals, and various norms
         *      -----------------------------------------------------------
         * */
        if( _info->dbug.debug )
                lsgrg_msg(IOINFO," ** ENTERING SUMRES ....\n" );

        *errmx = 0.0e0;
        *sinf = 0.0e0;
        *sinfd = 0.0e0;
        *ninfa = 0;
        *ninfd = 0;
        x[_info->dimen.n + _info->nintbk.nobj - 1] = g[_info->nintbk.nobj - 1];
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                isl = _info->dimen.n + i;
                if( g[i_] < alb[isl - 1] ){
                        /*???               X(ISL)=ALB(ISL) */
                        xval = alb[isl - 1];
                } else if( g[i_] > ub[isl - 1] ){
                        /*???               X(ISL)=UB(ISL) */
                        xval = ub[isl - 1];
                } else{
                        xval = g[i_];
                }
                /*???         ERR=G(I)-X(ISL) */
                err = g[i_] - xval;
                abserr = fabs( err );
                relerr = abserr/(1.0e0 + fabs( xval ));
                if( icstat[i_] < 0 ){
                        /*XXX        SINFD = SINFD + ABSERR */
                        if( relerr > _info->limits.epnewt ){
                                *ninfd = *ninfd + 1;
                                *sinfd = *sinfd + abserr;
                        }
                        error[i_] = 0.0e0;
                        x[isl - 1] = g[i_];
                } else{
                        /*XXXX       SINF = SINF + ABSERR */
                        *errmx = fmax( *errmx, relerr );
                        if( relerr > _info->limits.epnewt ){
                                *ninfa = *ninfa + 1;
                                *sinf = *sinf + abserr;
                        }
                        error[i_] = err;
                }

#ifdef IO_ENABLED
                if( _info->dbug.debug )
                        {
                        sprintf(LSGRG_MSGBUFFER,
                         " CHECKING CONS # %5ld ICSTAT= %5ld\n G= %14.6e"
                         " LB=%14.6e UB=%14.6e\nABSERR= %14.6e SINF ="
                         " %14.6e\n",
                          i, icstat[i_], g[i_], alb[isl - 1], ub[isl - 1],
                           abserr,
                         *sinf );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
        }
        return;
} /* end of function */


void /*FUNCTION*/ LSGRGCLASS droprw(LsgrgInfo *_info,long int irow, long int icstat[], long int istat[],
         long int *ndrop, long int *ninfa, long int *ninfd, double *sinf,
         double *sinfd, double error[], double g[], double x[], double alb[],
         double ub[], double *errmx)
{

         /*
         *      ---------------------------------------------------
         *     |  Set Dropped Flags, and Re-Compute the residuals
         *      ---------------------------------------------------
         * */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, " ** ENTERING DROPRW .... DROPPING ROW: %5ld\n",
                 irow );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( icstat[irow - 1] < 0 ){
#ifdef IO_ENABLED
                sprintf( LSGRG_MSGBUFFER,
                  " ***DROPRW ERROR .... CONSTRAINT %5ld HAS ALREADY"
                  " BEEN DROPPED. \n",
                    irow );
                lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                return;
        }
        *ndrop = *ndrop + 1;
        if( *ndrop == _info->dimen.m )
                return;
        icstat[irow - 1] = -1;
        istat[irow - 1] = -istat[irow - 1];

        /*      ----------------------------------------------------------
         *     |  For now, let SUMRES do the work.  Consider recalculating
         *     |  these variables directly.
         *      ---------------------------------------------------------- */
        sumres(_info, g, x, alb, ub, error, icstat, errmx, sinf, sinfd, ninfa,
         ninfd );

        return;
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS ph0log(LsgrgInfo *_info,long int n, long int mp1, long int nph0it,
         long int ndrop, long int nb, long int ninf, long int *ipr, long int *ipr3,
         long int *linect, LOGICAL32 newbas, LOGICAL32 degen, double sinf,
         double sinfd, double errmx, double step, long int ilv)
{
/*   ** 3/99 jcp ** ifdef out entire routine  IO_ENABLED not set   */

#ifdef IO_ENABLED

   char t = 'T';
   char f = ' ';
   char dgflag;
   long int k, nitr;
   char *columnheader=
        "\n\n"
        "  Itn   Ndrop  Ninf   Sinf       Sinfd       Max        Step"
        "      Lv   New   Deg\n"
        "  No.   Cons   Cons                         Resid       Size"
        "      Row  Basis Itr\n" ;

/* replace statics with alg members */

   #define termheader _info->ph0logStatic.termheader
   #define iprhd3     _info->ph0logStatic.iprhd3
   #define iprhld     _info->ph0logStatic.iprhld
   #define first      _info->firstcall.ph0log
/*
   static int termheader  = 0;
   static long int iprhd3, iprhld;
   static LOGICAL32 first = TRUE;
*/


        if( first ){
                iprhld = *ipr;
                iprhd3 = *ipr3;
                first = FALSE;
                if(_info->nintbk.ipr > 0 ) {
                   lsgrg_msg(IOINFO,"\n Beginning Phase0 Iterations");
                   if( _info->io_info.screen_output_enabled)
                     lsgrg_screen_msg(IOINFO,"\n Beginning Phase0 Iterations");
                }
        }
        nitr = nph0it - 1;
        if( _info->mngrg.iper != 0 )
                goto L_330;
/* L_3600: unreferenced??? */
        if( *ipr < 1 )
                goto L_380;
        *linect = *linect + 1;
        /*XXX  IF (LINECT.LT.48.AND.IPR3.EQ.0) GO TO 340 */
        if( *linect < 48 )
                goto L_340;
        if( nitr == 0 )
                {
                lsgrg_msg(IOINFO, "\n" );
                }
        if( *ipr3 == 0 && nitr > 1 )
                {
                lsgrg_msg(IOINFO, "\n" );
                }
        lsgrg_msg(IOINFO,columnheader);
        *linect = 0;
        goto L_340;
L_330:
        ;
        k = nitr/_info->mngrg.iper*_info->mngrg.iper;
        if( k != nitr && k != nitr - 1 )
                goto L_340;
        if( nitr == 0 ){
                lsgrg_msg(IOINFO,columnheader);
                /*XXXX    IF (IOTERM .GT. 0) WRITE (IOTERM,770) */
        }
        if( nitr < 2 )
                goto L_340;
        if( k == nitr - 1 ){
                lsgrg_msg(IOINFO,columnheader);
                /*XXXXX   IF (IOTERM .GT. 0) WRITE (IOTERM,770) */
        }
        *ipr = iprhld;
        *ipr3 = iprhd3;
L_340:
        ;
        dgflag = f;
        if( degen )
                dgflag = t;
        sprintf(LSGRG_MSGBUFFER,
         " %4ld %5ld %5ld %11.2e %11.2e %10.2e %10.2e %5ld %3c    %c\n",
         nph0it, ndrop, ninf, sinf, sinfd, errmx, step, ilv, TorF(newbas),
         dgflag );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);


        if( _info->mngrg.iper == 0 )
                goto L_380;
        if( k != nitr - 1 )
                goto L_380;
        if( nitr != 1 ){
                lsgrg_msg(IOINFO,columnheader);
                /*XX      IF (IOTERM .GT. 0) WRITE(IOTERM,770) */
        }
        *ipr = 1;
        *ipr3 = 0;
L_380:
        ;


        if( _info->io_info.screen_output_enabled  && _info->nintbk.ipr > 0){
                if( (nitr%24) == 0 || termheader ==0 ) {
                      sprintf( LSGRG_MSGBUFFER,"%s",columnheader);
                      termheader = 1;
                      lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
                 }

                sprintf( LSGRG_MSGBUFFER,
         " %4ld %5ld %5ld %11.2e %11.2e %10.2e %10.2e %5ld %3c    %c\n",
 /*                " %3ld %6ld  %5ld%12.4e%12.4e %9.2e%10.2e%5ld %4c%5c\n", */
                 nph0it, ndrop, ninf, sinf, sinfd, errmx, step, ilv,
                 TorF(newbas), dgflag );
                lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
        }

   #undef termheader
   #undef iprhd3
   #undef iprhld
   #undef first


#endif
        return;
        /*----------------------------------------------------------------------
         *
         *     End of PH0LOG
         * */
} /* end of function */

#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"

#define IOINFO  &(_info->io_info)
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)
/***********************************************************************
 **                 FILE:              GRGSETUP.C                      *
 **                 AUTHOR:            John C. Plummer                 *
 **        Rewrite of Original Translation to C (ALLAN D Waren)        *
 **        ca. 1992 from fortran source by Stuart Smith                *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **                                                                    *
 **   ... reorganized and rewritten for lsgrgc 2.0 j.c. plummer  3/98  *
 **   .   full dynamic allocation (no partitioning)                    *
 **   .   initial call to lsgrg_setupj() to get actual nonzero         *
 **        count                                                       *
 **  ** 11/98 jcp ** add translations of anajac and fdcheck from       *
 **    fortran version for lsgrgc 3.0 algorithm update. integrate into *
 **    setupj                                                          *
 **                                                                    *
 **                                                                    *
 **                                                                    *
 **  ROUTINE         DESCRIPTION                        LAST UPDATE    *
 ** ---------       -------------                      -------------   *
 **  lsgrg_setupj()                                        mar 1998    *
 **                computes jacobian pattern, sets up                  *
 **                data structures, implants linear                    *
 **                coefficients. counts coefficients                   *
 **                                                                    *
 **  lsgrg_setbounds()                                     mar 1998    *
 **                copies user specified row/var bounds                *
 **                into internal lsgrg arrays and checks               *
 **                consistency.  also assigns status for               *
 **                user specified linear variables                     *
 **                                                                    *
 **  lsgrg_memsetup()                                      apr 1998    *
 **                                                                    *
 **                allocates memory for lsgrg arrays                   *
 **                makes an initial call to lsgrg_setupj to obtain     *
 **                count of nonzeros at initial point for initial      *
 **                allocation of grad/ha                               *
 **                                                                    *
 ** lsgrg_initialize_model()                                           *
 **                assigns default bounds and status       mar 1998    *
 **                values for variables and constraints                *
 **                                                                    *
 **  lsgrg_set_invert_tolerances()                        mar 1998     *
 **                                                                    *
 **                sets some invert tolerances which                   *
 **                depend on user options (e.g. xajtol)                *
 **                                                                    *
 **  lsgrg_machine_precision()                           may 1998      *
 **                computes machine precision                          *
 **                                                                    *
 **  anajac()      builds jacobian from user-supplied                  *
 **                pattern and derivatives               nov 1998      *
 **                                                                    *
 **  fdcheck()     checks forward differences against                  *
 **                centrals                              nov 1998      *
 ***********************************************************************/

int LSGRGCLASS lsgrg_setupj(LsgrgInfo *_info,
           int build,double xtemp[],
           double g[], double gplus[],double gminus[],double gcol[],
           long length_grad, long *nnz_out,long *nnz_total_out,
           long *nnznl_out, long *colmax_out,
           long ipcol[], long iprow[], long ipmap[], double  paij[] )
{
   long int  i, icol,    irow, j, ne, colmax, nnz,nnznl;
   double   gtabs, gtemp;
/*-------------------------------------------------------------------*/
/*  lsgrg_setupj computes jacobian sparsity pattern from user's      */
/*  initial point, builds sparse matrix structures and implants      */
/*  constant coefficients                                            */
/*                                                                   */
/*  if build == 0, do not store structure.  just count nonzeros and  */
/*     return that count in nnz_out for memory allocation purposes   */
/*     2nd call with build==1 will actually build structure          */
/*                                                                   */
/*  return 0 of we run out of memory.  should never occur            */
/*-------------------------------------------------------------------*/
#ifdef IO_ENABLED
     if( _info->dbg.setupj ) {
       sprintf( LSGRG_MSGBUFFER,
        "\n... Entry lsgrg_setupj: build = %d "
        " KDERIV = %3ld   AIJTOL = %13.6e\n"
        " n = %d  nobj = %d mp1 = %d pstep = %13.6e\n",
             build,_info->kderiv, _info->xajtol,_info->n,_info->nobj,_info->mp1,_info->pstep );
        lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif
/*---------------------------------------------------------------------*/
/*        SET INITIAL X VECTOR FOR COEFFICIENT                         */
/*        EVALUATION to user initial point                             */
/* (assume that initial values have already been projected onto bounds)*/
/*---------------------------------------------------------------------*/
        for( i = 1; i <= _info->n; i++ )  xtemp[i] = _info->x[i];

        lsgrg_gcomp1(_info, g, xtemp );     /* get functions as initial pt */

        if( _info->maxim ) g[_info->nobj] = -g[_info->nobj];

#ifdef IO_ENABLED
        if( _info->dbg.setupj ) {
                sprintf( LSGRG_MSGBUFFER,
                 "\n... Initial Point for Jacobian Struture " );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                for( i = 1; i <= _info->n; i++ ) {
                   LSGRG_MSGBUFFER[0] = ' ';
                   if(i%2==1)  LSGRG_MSGBUFFER[0] = '\n';
                   sprintf( LSGRG_MSGBUFFER+1,
                         " x[%6d] = %20.13e", i,xtemp[i] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                sprintf( LSGRG_MSGBUFFER,
                 "\n... Initial g for Jacobian Struture " );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);

                for( i = 1; i <= _info->mp1; i++ ) {
                   sprintf( LSGRG_MSGBUFFER,
                      "\n g[%6d] = %13.6g",i,  g[i] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
        }
#endif
/*----------------------------------------------------------------------*/
/*     CALL APPROPRIATE PARSH TO GET ELELMENTS                          */
/*     COLMAX -- USED TO HOLD LENGTH OF LONGEST COLUMN                  */
/*               IT IS AT LEAST 1, BECAUSE OF SLACKS                    */
/* **fixme**   iheg0,ihegl0 are for 0-based alg routines, expunge later */
/*             or expunge iheg/ihegl
/*----------------------------------------------------------------------*/
        nnz        = 0; /* structural nonzeros */
        nnznl      = 0; /* nonlinear nonzeros ?? */
        colmax     = 1;
/*--------------------------------------------------------------------*/
/*  if user analytic derivatives are specified, initial call to parsh */
/*  was executed in grgmem to get nnz confirmation and derivatives at */
/*  initial point.  Here, call anajac to generate jacobian structures */
/*  NOTE: anajac was translated from fortran and uses 0-based arrays. */
/*  shift pointers for all array arguments                            */
/*  NOTE2: anajac needs an integer work array named icount of length n*/
/*      pass in the basis inverse array and cast it to (long *)       */
/*------------------------------------------------------------------- */
        if(_info->kderiv == USER_ANALYTIC) {
            anajac(_info,
                    (xtemp+1), ((_info->grad)+1), ((_info->ihag)+1),
                    ((_info->iheg)+1),((_info->ihegl)+1), (paij+1), (iprow+1),
                    ((ipcol+1)), (ipmap+1), (long *)_info->inv_binv,
                    &colmax, &nnznl, _info->nzgrad, _info->n);
            nnz = _info->nzgrad;
            for(i=1; i<=_info->n; ++i) {
               _info->iheg0[i-1] = _info->iheg[i];
               _info->ihegl0[i-1] = _info->ihegl[i];
            }
            if(_info->dbg.setupj)
               lsgrg_check_anajac(_info,
                                  _info->n, _info->nzgrad,_info->grad,
                                  _info->ihag, _info->iheg,_info->iheg0,
                                  _info->ihegl,_info->ihegl0,iprow, ipcol,
                                  ipmap, paij);
        }
        else {
/*------------------------------------------------------------------- */
/*  compute jacobian pattern at initial point via finite differences  */
/*------------------------------------------------------------------- */
          ++_info->nparsh;
          for( j = 1; j <= _info->n; j++ )  { /* start loop over all columns */

             if(build)    {
                _info->iheg[j] = nnz + 1;
                _info->iheg0[j-1] = _info->iheg[j];

#ifdef IO_ENABLED
                if( _info->dbg.setupj ) {
                   sprintf( LSGRG_MSGBUFFER,
                    "\n  ... Column %d Derivatives: IHEG(%6ld) = %7ld\n",
                              j, j, _info->iheg[j] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
             }
             if(_info->kderiv == FORWARD_DIFFERENCES) {
                  lsgrg_parshf1(_info, xtemp, g, gcol, j, gplus, _info->ub[j] );
             }

             if(_info->kderiv == CENTRAL_DIFFERENCES) {
                  lsgrg_parshc1(_info, xtemp, gcol, j, gplus, gminus );
             }
/*------------------------------------------------------------------- */
/*    build this column in grad, return 0 if we run out of space      */
/*    which should never occur since we have made an initial call     */
/*    to count the nonzeros                                           */
/*------------------------------------------------------------------- */
                icol = 0;
                for( irow = 1; irow <= _info->mp1; irow++ ) {
                   gtemp = gcol[irow];
                   gtabs = fabs( gtemp );
                   if( gtabs >=  _info->xajtol ) {
                       ++icol;
                       ++nnz;
                       if(build) {
                           if( nnz > length_grad) {
#ifdef IO_ENABLED
                          sprintf( LSGRG_MSGBUFFER,
                          "\n Error: number of nonzeros > lgrad [%d]"
                          " Should Never Occur ",length_grad);
                          lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                             _info->lsgrg_return_status = _LSGRG_INSFMEMORY;
                             return 0;
                          }

                           _info->grad[nnz] = gtemp;
/* take these stmts out
                           if( (irow == _info->nobj) && (_info->maxim) )
                                   _info->grad[nnz] = -gtemp;
*/
                           _info->ihag[nnz] = irow;
                       } /* end if build */
#ifdef IO_ENABLED
                       if( _info->dbg.setupj ) {
                          sprintf( LSGRG_MSGBUFFER,
                          "\n NONZERO # %6ld ROW = %6ld  ELEMENT = %13.6e",
                              nnz, irow, gtemp );
                          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                  } /* end if(gtabs */
                } /* end of row loop this column */

                if( icol > colmax ) colmax = icol;
                if( labs( _info->ivstat[j] ) == 1 ) nnznl += icol;
        } /* end column loop  */
        if(_info->dbg.setupj) lsgrg_msg(IOINFO,"\n .. end of column loop ");

    }  /* end of finite difference pattern detection */
/*------------------------------------------------------------------*/
/*    SET UP SLACK COLUMNS AT END OF GRAD ARRAY                     */
/*------------------------------------------------------------------*/
    ne = nnz;
    for( i = 1; i <= _info->mp1; i++ ) {
        ++ne;
        if(build) {
           _info->iheg[_info->n + i]      = ne;
           _info->iheg0[_info->n + i - 1] = ne;
           _info->grad[ne]                = -1.00e0;
           _info->ihag[ne]                = i;
        }
    }
    if(!build) {   /* exit if we are only counting nonzeros */
       *nnz_out = nnz;
       return 1;
    }
    if(_info->dbg.setupj) lsgrg_msg(IOINFO,"\n... after setup of slack columns..");
    _info->iheg[_info->n+i] = ne + 1;
    _info->iheg0[_info->n+i-1] = ne + 1;
/*------------------------------------------------------------------*/
/*  If we are doing finite differences
/*    STORE THE LINEAR ELEMENT POINTERS IN IHEGL.                   */
/*    THIS IS NOT REALLY NEEDED FOR FINITE DIFFERENCING, BUT IS     */
/*    USED SO THAT GAMS INTERFACE LOOKS THE SAME                    */
/*------------------------------------------------------------------*/
   if(_info->kderiv == FORWARD_DIFFERENCES ||
     _info->kderiv == CENTRAL_DIFFERENCES )  {

        if(_info->dbg.setupj) lsgrg_msg(IOINFO,"\n... starting setup of ihegl....");
        for( i = 1; i <= _info->n; i++ ) {
           if( labs( _info->ivstat[i] ) == 2 ) {
             _info->ihegl[i]    = _info->iheg[i];
             _info->ihegl0[i-1] = _info->iheg[i];
           }
           else {
             _info->ihegl[i]    = _info->iheg[i + 1];
             _info->ihegl0[i-1] = _info->iheg[i + 1];
           }
        }
        if(_info->dbg.setupj) lsgrg_msg(IOINFO,"\n... after  setup of ihegl....");
   } /* end finite difference ihegl setup */

        *nnz_out = nnz;
        *nnznl_out = nnznl;
        *nnz_total_out = ne;
        *colmax_out = colmax;
        return 1;

} /* end of lsgrg_setupj */

/*....................................................................*/
int  LSGRGCLASS lsgrg_setbounds(LsgrgInfo *_info,
                     double blvar[],double buvar[],
                     double blcon[],double bucon[],double xx[],
                     long lvars [], long ivstat[])
{
long int i,  indx, error_count = 0,j;
double bl, bu, plinfy = _info->plinfy;
/*---------------------------------------------------------------------*/
/*   check consistency of user-specified bounds                        */
/*---------------------------------------------------------------------*/
  /* lsgrg_echobnd(xx, blvar, buvar, blcon, bucon,  lvars,
               "at entry setgbounds"); */
        for( i = 1; i <= _info->nvars; i++ ) {

            if( buvar[i] < blvar[i] )  {
                ++error_count;

#ifdef IO_ENABLED
              sprintf( LSGRG_MSGBUFFER,
               "\n Error: Var [%d] \n "
               "Upper Bound = %15.8e  <  Lower Bound = %15.8e",
                 i, buvar[i], blvar[i] );
              lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
            }
        }  /* end var check */

        for( i = 1; i <= _info->nrows; i++ )  {
           if( bucon[i] < blcon[i] && i != _info->nobj ) {
               ++error_count;

#ifdef IO_ENABLED
              sprintf( LSGRG_MSGBUFFER,
               "\n Error: Row [%d] \n "
               "Upper Bound = %15.8e  <  Lower Bound = %15.8e",
                 i, bucon[i], blcon[i] );
              lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
           }
        } /* end row check */
        if(error_count) {
           _info->lsgrg_return_status = _LSGRG_BOUNDS_ERROR;
           return 0;
         }
/*---------------------------------------------------------------------*/
/*   set status code for linear variables  (ivstat = 2)                */
/*---------------------------------------------------------------------*/
        if( lvars[0] ) {
            for( i = 1; i <= lvars[0]  ; i++ )  {
                if( lvars[i] < 1 || lvars[i] > _info->nvars ) {

#ifdef IO_ENABLED
                   sprintf(LSGRG_MSGBUFFER,
                     "\n Error: lvars[%d] = %d < 1 or > nvars [%d]",
                     i,lvars[i],_info->nvars);
                   lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                   _info->lsgrg_return_status = _LSGRG_LINEAR_VARS_ERROR;
                   return 0;
                }
                ivstat[lvars[i]] = V_LINEAR;
             }
        }
/*---------------------------------------------------------------------*/
/*    copy bounds from user arrays to lsgrg arrays                     */
/*     SET STATUS OF FIXED VARIABLES (IVSTAT(I) = -1 ,-2   VAR I FIXED)*/
/*     ALSO SET UP LIST OF NONLINEAR INDICES IN INLIN                  */
/*---------------------------------------------------------------------*/
        indx = 0;
        _info->nfix = 0;
        for( i = 1; i <= _info->n; i++ )  {   /* variable bounds */
                bu = buvar[i];
                bl = blvar[i];
                if( bl < -1.0e-2*plinfy ) bl = -plinfy;
                if( bu >  1.0e-2*plinfy ) bu = plinfy;
                _info->alb[i] = bl;
                 _info->ub[i] = bu;
                if( ivstat[i] == V_NONLINEAR )  {
                        ++indx ;
                        _info->inlin[indx] = i;
                }
                if( bl == bu ) {
                        ivstat[i] = -ivstat[i];
                        ++_info->nfix;
                }
        }
/*---------------------------------------------------------------------*/
/*   check for gross error in lvars indx should be count of nl vars    */
/*---------------------------------------------------------------------*/
        if( (indx + lvars[0] != _info->n) || (indx != _info->nnlin) )  {

#ifdef IO_ENABLED
           sprintf( LSGRG_MSGBUFFER, "\nERROR: ERROR IN TOTAL NUMBER OF "
              "VARIABLES\n N = %6ld NLIN = %6ld NNLIN = %6ld"
              "\n Count of nonlinear variables = %6ld",
                  _info->n, _info->nlin, _info->nnlin, indx );
            lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
            _info->lsgrg_return_status = _LSGRG_LINEAR_VARS_ERROR;
            return 0;
        }
/*---------------------------------------------------------------------*/
        for( i = 1; i <= _info->mp1; i++ )  {      /* constraint bounds */
           if( i == _info->nobj ) continue;
              _info->alb[_info->n + i]  = blcon[i];
               _info->ub[_info->n + i]  = bucon[i];

              if( _info->alb[_info->n + i] < -1.0e-2*plinfy ) _info->alb[_info->n + i] = -plinfy;
              if(  _info->ub[_info->n + i] >  1.0e-2*plinfy ) _info->ub[_info->n + i] = plinfy;
              if( _info->alb[_info->n + i] <= -plinfy && _info->ub[_info->n + i] >= plinfy )
                      _info->istat[i] = 0;
              if( _info->alb[_info->n + i] == _info->ub[_info->n + i] )
                      _info->istat[i] = 1;
        }

        for( i = 1; i <= _info->n; i++ )   /*  assign initial var values */
                _info->x[i] = xx[i];
/*---------------------------------------------------------------------*/
/*  check for vars out of bounds and project onto closest bound        */
/*  if that is the case                                                */
/*---------------------------------------------------------------------*/
        for( i = 1; i <= _info->n; i++ ) {

           if( _info->x[i] > _info->ub[i] )  {


#ifdef IO_ENABLED
             if( _info->ipr >= 1 )  {
                sprintf( LSGRG_MSGBUFFER,
                 "\n x[%d] Initial Value of %15.8e Changed to "
                 " Upper Bound = %15.8e",i,_info->x[i],_info->ub[i]);
                lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
             }
#endif
             _info->x[i] = _info->ub[i];
           }  /* close check on upper bound */

           if(_info->x[i] < _info->alb[i] ) {

#ifdef IO_ENABLED
             if( _info->ipr >= 1 )  {
                sprintf( LSGRG_MSGBUFFER,
                 "\n x[%d] Initial Value of %15.8e Changed to "
                 " Lower Bound = %15.8e",i,_info->x[i],_info->alb[i]);
                lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
             }
#endif
             _info->x[i] = _info->alb[i];
           } /* close check on lower bound */

        _info->x0[i] = _info->x[i];     /* save initial value after any adjustment */
        }

          _info->istat[_info->nobj] = 0;       /* status and bounds for obj row */
        _info->alb[_info->n + _info->nobj] = -plinfy;
         _info->ub[_info->n + _info->nobj] = plinfy;
/*-------------------------------------------------------------------*/
/*  now do a brute force consistency check on the linear var set     */
/*-------------------------------------------------------------------*/
     for(i=1;i<=_info->n;++i) {
         if(abs(ivstat[i]) == V_LINEAR) {
            for(j=1;j<=_info->nlin && (lvars[j] !=i);++j);
            if(j > _info->nlin) {
#ifdef IO_ENABLED
               sprintf(LSGRG_MSGBUFFER,
                 "\n Internal Error: Variable %d with linear ivstat"
                 " not found in lvars",i);
               lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
               _info->lsgrg_return_status = _LSGRG_LINEAR_VARS_ERROR;
               return 0;
            }
         }
     }
     for(i=1;i<=_info->nlin;++i) {
         if(abs(ivstat[lvars[i]]) != V_LINEAR) {
#ifdef IO_ENABLED
               sprintf(LSGRG_MSGBUFFER,
                 "\n Internal Error: Lvars[%d] = %d but ivstat[%d] not"
                 " = linear(%d)",i,lvars[i],ivstat[lvars[i]],V_LINEAR);
               lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
               _info->lsgrg_return_status = _LSGRG_LINEAR_VARS_ERROR;
               return 0;
         }
     }

/*  end of checks */
        return 1;

} /* end of function lsgrg_setbounds */

void LSGRGCLASS lsgrg_initialize_model(LsgrgInfo *_info,
                                       long nvars, long nrows)
{
/*-------------------------------------------------------------------*/
/*  lsgrg_initialize model sets variable and constraint bounds and   */
/*  status vectors to default values                                 */
/*                                                                     */
/*     INITIALIZE BOUNDS ET AL. TO DEFAULT VALUES                      */
/*                                                                     */
/*     IVSTAT IS USED TO INDICATE VARIBLE TYPE:                        */
/*          ABS(IVSTAT(I)) = 1  ==>  VARIABLE "I" IS NONLINEAR         */
/*          ABS(IVSTAT(I)) = 2  ==>  VARIABLE "I" IS LINEAR            */
/*                                                                     */
/*     NOTE: IF IVSTAT(I)  < 0  ==>  VARIABLE "I" IS FIXED             */
/*                                                                     */
/*     DEFAULT VARIABLE STATUS IS NONLINEAR, NOT FIXED (IVSTAT(I) = 1) */
/*---------------------------------------------------------------------*/
long i;
double plinfy = _info->plinfy;
/*-----------------------------------------------------------------*/
/*  initialize lsgrg bounds and status arrays to default values    */
/*-----------------------------------------------------------------*/
        for( i = 1; i <= nvars; i++ ) {
           _info->ivstat[i]  = V_NONLINEAR;  /* default is nonlinear var */
             _info->x[i]     = 0.0;
           _info->alb[i]     = -plinfy;    /*  default is unbounded */
            _info->ub[i]     =  plinfy;
        }


/*   DEFAULT CONSTRAINT STAUS IS INEQUALITY GREATER THAN 0 */

        for( i = 1; i <= nrows; i++ )  {
                  _info->istat[i] = 2;
                _info->alb[_info->n + i] = 0.0;
                 _info->ub[_info->n + i] = plinfy;
        }

     return ;
} /* end function lsgrg_initialize_model  */

/*....................................................................*/
void LSGRGCLASS lsgrg_set_invert_tolerances(LsgrgInfo *_info)
{
/*     SET GRGINVRT TOLERANCES IN DLUPRM   */

        _info->dluprm[0L] = 1.0e-1;
        _info->dluprm[1L] = _info->xajtol;
        return;
}    /* end function lsgrg_set_invert_tolerances */

double LSGRGCLASS lsgrg_machine_precision()
{
   double epsx, eps1;

        epsx = 1.0e0;
        eps1 = 1.1;
        while(eps1 > 1.0e0 ) {
           epsx *= 0.5e0;
           eps1 = epsx + 1.0e0;
        }
        return epsx;
}  /* end function lsgrg_machine_precision */


int  /*FUNCTION*/ LSGRGCLASS anajac(LsgrgInfo *_info,
         double x[], double grad[], long int ihag[],
         long int iheg[], long int ihegl[], double paij[], long int iprow[],
         long int ipcol[], long int ipmap[], long int icount[], long int *colmax,
         long int *nznl, long int nnz_user, long n)
{
long int i, i_, indx, nnz;


        /*.......................................................................
         *
         *                      *** ANAJAC ***
         *
         *     *****PURPOSE:
         *     THIS SUBROUTINE SETUPS THE DATA STRUCUTURES FOR THE JACOBIAN WHEN USING
         *     ANALYTICAL PARSH (KDERIV=2).
         *
         *
         *     *****ARGUMENT DESCRIPTION:
         *     ON INPUT:
         *
         *     X        - Array of variable values
         *     NZ       - THE NUMBER OF NON-ZERO JACOBIAN ELEMENTS.
         *     NP1      - THE NUMBER OF COLUMNS+1.
         *
         *     ON OUTPUT:
         *     GRAD     - CONTAINS VALUES FOR NON-ZERO ELEMENTS PACKED BY COLUMNS
         *     IHAG     - ROW INDECES OF EACH ELEMENT
         *     IHEG     - START OF COLUMN POINTERS
         *     IHEGL    - FIRST LINEAR ELEMENT IN EACH COLUMN
         *     IAMAP    - MAPS ELEMENTS IN INPUT ARRAY AIJ TO CORRESPONDING ELEMENT IN GRAD
         *     NZNL     - NUMBER OF NONLINEAR NONZERO JACOBIAN ELEMENTS
         *
         *     WORK ARRAYS PASSED TO PARSH:
         *     AIJ      - ARRAY PASSED TO PARSH TO STORE NON-ZERO ELEMENTS
         *     IROW     - ROW INDEX FOR EACH ELEMENT OF AIJ
         *     ICOL     - COLUMN INDEX FOR EACH ELEMENT OF AIJ
         *
         *     LOCAL WORK ARRAYS:
         *     ICOUNT   - USED TO COUNT NUMBER OF NONZEROS IN EACH COLUMN
         *
         *
         *     *****REFERENCES:
         *
         *     *****HISTORY:
         *     WRITTEN BY:  STUART H. SMITH
         *                  AUSTIN, TX 78728  (512-990-5296)
         *
         *     DATE LAST MODIFIED:  28 APR 1998
         *
         *  UPDATE LOG:
         * --------------
         * ** 11/98 jcp **
         *      translation to C and mating to lsgrgc3 shell protocols
         *      deleted arg np1, never used
         *
         *.......................................................................
         * */

        /*................................................................
         *................ COMMON DECLARATIONS ...........................
         *................................................................
         *
         * */




        /*.................................................................
         *................ ARGUMENT DECLARATIONS ..........................
         *.................................................................
         * */

        /*.................................................................
         *................. LOCAL DECLARATIONS ............................
         *.................................................................
         * */


        /*.................................................................
         *
         *     *** START OF ANAJAC ****
         *
         *     First, call user supplied PARSH to compute derivitives at
         *      the initial point.
         *     NOTE: We are wasting an evaulation here, should set COMPGR later
         *   NOTE2: parsh assumes 1-based arrays, so adjust pointers
         *          before call
         */
        ++_info->nparsh;
        //_info->lsgrg_user_parsh( (x-1), n, (paij-1),
        //                        (iprow-1), (ipcol-1), &nnz );
		//
		// A.G.
		//
        LSGRGSolver* pObject = (LSGRGSolver*)_info->m_pOwner;
        P_PARSH pfnparsh = _info->lsgrg_user_parsh;
        (pObject->*pfnparsh)( (x-1), n, (paij-1), (iprow-1), (ipcol-1), &nnz );


/* check that parsh returns same nbr of nonzeros as specified in grgsub */

      if(nnz != nnz_user) {
#ifdef IO_ENABLED
      sprintf(LSGRG_MSGBUFFER,
        "\n Error: Nbr of Jacobian Nonzeros Returned From User PARSH"
        "[%d] \n Does Not Match Value in Call to GRGSUB [%d]",
        nnz,nnz_user);
      lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
      _info->lsgrg_return_status = _LSGRG_BAD_USER_NNZ;
      return 0;
    }

        /*     Count the number of nonzeroes in each column
         *     Track the number of nonlinear elements */

        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                icount[i_] = 0;
        }

        *nznl = 0;
        for( i = 1; i <= nnz_user; i++ ){
                i_ = i - 1;
                icount[ipcol[i_] - 1] = icount[ipcol[i_] - 1] + 1;
                if( iprow[i_] > 0 )
                        *nznl = *nznl + 1;
        }

        /*     Set up START of Column pointers
         *     Initialize linear pointers to start of next column */

        iheg[0] = 1;
        *colmax = 0;
        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                iheg[i_ + 1] = iheg[i_] + icount[i_];
                ihegl[i_] = iheg[i_ + 1];
                if( icount[i_] > *colmax )
                        *colmax = icount[i_];
        }

        /*     Now, Store Nonlinear and Linear elements and set up mapping array
         *     NOTE: We even map linear elements, although they will not be changed later
         *     Use ICOUNT to hold place of last nonlinear element stored
         *     IHEGL holds place of last linear element stored */

        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                icount[i_] = iheg[i_] - 1;
        }

        for( i = 1; i <= nnz_user; i++ ){
                i_ = i - 1;

                /*        Store linear elements from bottom up ; nonlinears from top down */
                if( iprow[i_] < 0 ){
                        ihegl[ipcol[i_] - 1] = ihegl[ipcol[i_] - 1] - 1;
                        indx = ihegl[ipcol[i_] - 1];
                } else{
                        icount[ipcol[i_] - 1] = icount[ipcol[i_] - 1] + 1;
                        indx = icount[ipcol[i_] - 1];
                }

                /*        Store the element and row and set pointer */
                grad[indx - 1] = paij[i_];
                ihag[indx - 1] = labs( iprow[i_] );
                ipmap[i_] = indx;
        }

        /*     DEBUG Sanity check, just to make sure we didn't miss anything */
        for( i = 1; i <= n; i++ ){
                i_ = i - 1;
                if( icount[i_] + 1 != ihegl[i_] ){
#ifdef IO_ENABLED
                        sprintf( LSGRG_MSGBUFFER,
                          "\n ANAJAC ERROR: COLUMN COUNTS ARE INCORRECT"
                          " FOR COLUMN %8ld\n LAST NONLINEAR INDEX = %7ld"
                          " FIRST LINEAR INDEX = %8ld\n",
                           i, icount[i_], ihegl[i_] );
                         lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
                         _info->lsgrg_return_status = _LSGRG_ANAJAC_BAD_COL;
                         return 0;
                }
        }

        return 1;

} /* end of function */

void /*FUNCTION*/ LSGRGCLASS fdcheck(LsgrgInfo *_info,
         double x[], double grad[], long int ihag[],
         long int iheg[], double gplus[], double gminus[],
         double gcol[] )
{
//long int i, i_, j, j_;
//double diff;

        /*.......................................................................
         *
         *                      *** FDCHECK ***
         *
         *     *****PURPOSE:
         *     THIS SUBROUTINE CHECKS THE USER SUPPLIED ANALYTICAL DERIVATIVES AGAINST
         *     FINITE DIFFERENCING WHEN KDERIV=2.
         *
         *     ****NOTES:
         *     Finite differencing is performed using central differencing.  An error
         *     is flagged if the relative difference is greater than 100*PSTEP, where
         *     PSTEP is the finite difference step size.
         *
         *
         *     *****ARGUMENT DESCRIPTION:
         *     ON INPUT:
         *
         *     X        - Array of variable values
         *     UB       - Array of upper bounds for variables
         *     GRAD     - CONTAINS VALUES FOR NON-ZERO ELEMENTS PACKED BY COLUMNS
         *     IHAG     - ROW INDECES OF EACH ELEMENT
         *     IHEG     - START OF COLUMN POINTERS
         *     NZ       - Number of nonzeroes in grad array
         *     NP1      - THE NUMBER OF COLUMNS+1.
         *
         *     WORK ARRAYS:
         *     G       -  Work array for function calls in forward/central differencing
         *     GPLUS   -  Work array for function calls in forward/central differencing
         *     GMINUS  -  Work array for function calls in forward/central differencing
         *     GCOL    -  Used to return derivatives for a specific variable
         *
         *
         *     *****REFERENCES:
         *
         *     *****HISTORY:
         *     WRITTEN BY:  STUART H. SMITH
         *                  AUSTIN, TX 78728  (512-990-5296)
         *
         *     DATE LAST MODIFIED:  20 May 1998
         *
         *  UPDATE LOG:
         * --------------
         * ** 11/98 jcp **
         *  translation to c and mating to lsgrgc3 shell protocols
         *  lsgrgc3 parshc0 does not need 'g' as argument, deleted
         *  deleted args 'nz' and 'np1' never used
         *.......................................................................
         * */

        /*.................................................................
         *
         *     *** START OF FDCHECK **** */

#ifdef IO_ENABLED

        /*     IF (IPR .GT. 0 ) THEN */
        lsgrg_msg(IOINFO,
         "\n FDCHECK:  Checking Analytical Derivatives using Central Differencing \n\n" );

        /*     ENDIF
         *
         *      ----------------------------------------------
         *     | Calculate constraints at initial value
         *      ----------------------------------------------
         *
         *     CALL GCOMP_LSGRG(G,X)
         *
         *      --------------------------------------------------
         *     | Use Central Differencing to
         *     | Compute the Jacobian for each nonlinear column
         *     | Check each element to a relative tolerance of PSTEP
         *      --------------------------------------------------
         * */
        for( j = 1; j <= _info->dimen.n; j++ ){
                j_ = j - 1;

                lsgrg_parshc0(_info, x, gcol, j, gplus, gminus );
                ++_info->nparsh;

                /*         -----------------------------
                 *C       | Check user predicted nonzeroes
                 *C        ----------------------------- */

                if( iheg[j_ + 1] <= iheg[j_] )
                        goto L_45;
                for( i = iheg[j_]; i <= (iheg[j_ + 1] - 1); i++ ){
                     i_ = i - 1;
                     diff = fabs( gcol[ihag[i_] - 1] - grad[i_] );
                     if( diff > (100.e0*_info->stepbk.pstep*(1.0e0 + fabs( grad[i_] ))) ){
                         _info->drvchk.nfddif = (_info->drvchk.nfddif) + 1;
                          sprintf(LSGRG_MSGBUFFER, " ROW: %5ld COL: %5ld ANA: %14.7e FD: %14.7e ERROR: %14.7e\n",
                           ihag[i_], j, grad[i_], gcol[ihag[i_] - 1], diff/(1.0e0 +
                           fabs( grad[i_] )) );
                          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                        gcol[ihag[i_] - 1] = 0.0e0;
                }
L_45:
                ;


                /*C        -------------------------------------------
                 *C       | Check if any unexpected nonzeroes were created
                 *C        -------------------------------------------
                 * */
                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        i_ = i - 1;
                        if( fabs( gcol[i_] ) > _info->stepbk.pstep ){
                                sprintf(LSGRG_MSGBUFFER, " ROW: %5ld COL: %5ld ANA: %14.7e FD: %14.7e ERROR: %14.7e\n",
                                 i, j, 0.0e0, gcol[i_], fabs( gcol[i_] ) );
                                _info->drvchk.nfddif = _info->drvchk.nfddif + 1;
                                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
                }

        }

        sprintf(LSGRG_MSGBUFFER,
         "\n FDCHECK: Number of Derivative errors found = %6ld\n",
         _info->drvchk.nfddif );
        lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
#endif

        return;
} /* end of function */
void LSGRGCLASS lsgrg_echo_lvars(LsgrgInfo *_info,
              long nvars, long ivstat[], long lvars[],
                      long nnlin, long nlin)
{
/*  consistency echo for linear variable specs */
#ifdef IO_ENABLED
   int i;
   sprintf(LSGRG_MSGBUFFER,
       "\n\n.... Linear Variables Specified: lvars[0] = %d ...\n"
       " nvars = %d nnlin = %d nlin = %d\n",
        lvars[0], nvars, nnlin, nlin);
   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
   for(i=1;i<=lvars[0];++i){
      sprintf(LSGRG_MSGBUFFER,
         "\n  linear var %5d  is var %5d, ivstat = %d",
           i,lvars[i],ivstat[lvars[i]]);
      lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
   }
#endif
   return;
}
void LSGRGCLASS lsgrg_check_anajac(LsgrgInfo *_info,
                        long n, long nzgrad,double *grad,
                        long *ihag, long *iheg, long *iheg0,
                        long *ihegl,long *ihegl0,long *iprow,
                        long *ipcol,long *ipmap, double *paij)
{
#ifdef IO_ENABLED
    int i,k,k1,k2,kk;
    sprintf(LSGRG_MSGBUFFER,
       "\n ... entry lsgrg_check_anajac: n = %d, nzgrad = %d\n",
         n,nzgrad);
    lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
    for(i=1; i<=nzgrad; ++i) {
        sprintf(LSGRG_MSGBUFFER,
          "\n [%5d] iprow = %5d, ipcol = %5d, paij = %12.6e, ipmap = %5d"
          "\n        ihag[ipmap[i]] = %5d, grad[ipmap[i]] = %12.6e",
          i,iprow[i], ipcol[i], paij[i], ipmap[i],ihag[ipmap[i]],
          grad[ipmap[i]] );
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
    }
    lsgrg_msg(IOINFO,"\n... check of iheg,ihegl.....");
    for(i=1; i<=n;++i) {
        k1 = iheg[i];
        k2 = iheg[i+1] - 1;
        sprintf(LSGRG_MSGBUFFER,
          "\n i = %d"
          "\n  iheg[i]    = %5d, iheg[i+1]-1 = %5d"
          "\n  iheg0[i-1] = %5d  ihegl[i]    = %5d"
          "\n                    ihegl0[i-1] = %5d\n",
          i,k1,k2,ihegl[i],ihegl0[i], iheg0[i-1], ihegl0[i-1]);
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        for(k=k1;k<=k2;++k) {
           for(kk=1;kk<=nzgrad;++kk)
               if(ipmap[kk] == k) break;

           if(kk > nzgrad) kk = -9999;
           sprintf(LSGRG_MSGBUFFER,
            "\n k = %5d ihag[k] = %5d kk = %d",
               k,ihag[k],kk);
           lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
           if(kk != -9999)
              sprintf(LSGRG_MSGBUFFER,
                " iprow[kk] = %5d",iprow[kk]);
       }
    }
#endif
   return;
}
#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"
/*
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
*/
/*------------------------------------------------------------------*/
/*   grgiface.c  Author: John C. Plummer    10/98                   */
/*                                                                  */
/*   grgiface contains the routines which interface the lsgrgc      */
/*   shell with the algorithm routines.  Mainly, these routines     */
/*   copy values from the 'global' structure defined for the shell  */
/*   to and from the 'global' structures (derived from the          */
/*   original fortran common blocks)                                */
/*   in the algorithm routines                                      */
/*                                                                  */
/*  . lsgrg_get_globals() will also transfer the algorithm timing   */
/*     stats via calls to shell routine lsgrg_timestats()           */
/*                                                                  */
/*  . the global gcomp and parsh counts kept by the algorithm are   */
/*    not copied since the shell gcomp and parsh functions keep     */
/*    their own counts and the old counts are redundant             */
/*------------------------------------------------------------------*/
void LSGRGCLASS lsgrg_put_globals(LsgrgInfo *shell)
{
  int i;
/*------------------------------------------------------------------*/
/* lsgrg_put_globals() copies global values from shell to alg       */
/* and initializes some counters and timing stats in the algorithm  */
/* routines                                                         */
/*------------------------------------------------------------------*/
/* NOTE:  optblk.subset controls gcompx calls for individual        */
/*        function access. not implemented this version             */
   shell->optblk.subset = 0;

   shell->dbug.debug = FALSE;                /* set algorithm debug flag */
   if(shell->dbg.algorithm) shell->dbug.debug = TRUE;

   shell->epscom.eps0   =    shell->eps0;  shell->epscom.eps3   =    shell->eps3;
   shell->epscom.eps1   =    shell->eps1;  shell->epscom.eps4   =    shell->eps4;
   shell->epscom.eps2   =    shell->eps2;  shell->epscom.eps5   =    shell->eps5;

   shell->limits.epboun  =   shell->epboun;  shell->limits.epnewt  =   shell->epnewt;
   shell->limits.epspiv  =   shell->epspiv;

   shell->ingrg.epinit   =   shell->epinit;  shell->ingrg.eplast   =   shell->eplast;
   shell->ingrg.epdeg    =   shell->epdeg;

   shell->tols.eps       =   shell->eps;     shell->tols.plinfy    = shell->plinfy;
   shell->tols.tolx      =   shell->tolx;    shell->tols.plzero    = shell->plzero;
   shell->tols.tolz      =   shell->tolz;

   shell->mngrg.epstop   =   shell->epstop;

   shell->zcond.cndtol   =   shell->cndtol;
   shell->sernew.edf     =   shell->edf;

   shell->stepbk.funpr   =   shell->funpr;   shell->stepbk.pstep   =   shell->pstep;

   shell->ph1bk.phmult   =   shell->phmult;  shell->ph1bk.ph1eps   =   shell->ph1eps;

   shell->ztol.xajtol    =   shell->xajtol;
   shell->pvtblk.xpvpct  =   shell->xpvpct;

    for(i=1; i <= 4; ++i) shell->lsinvt.dluprm[i-1] =  shell->dluprm[i];
    for(i=1; i <= 7; ++i) shell->pivots.pivtol[i-1] =  shell->pivtol[i];

    shell->bestbk.truobj = shell->truobj;

  /* longs */

      shell->mngrg.limser  =  shell->limser;  shell->mngrg.ipn4    =  shell->ipn4;
      shell->mngrg.nstop   =  shell->nstop;   shell->mngrg.ipn5    =  shell->ipn5;
      shell->mngrg.iper    =  shell->iper;    shell->mngrg.ipn6    =  shell->ipn6;

      shell->quadbk.iquad  =  shell->iquad;
      shell->cmax.colmax   =  shell->colmax;
      shell->pvtblk.ibvblm =  shell->ibvblm;
      shell->degn.idglim   =  shell->idglim;

      shell->nintbk.ipr    =  shell->ipr;     shell->nintbk.nsuper =  shell->nsuper;
      shell->nintbk.ipr3   =  shell->ipr3;    shell->nintbk.nb     =  shell->nb;
      shell->nintbk.nobj   =  shell->nobj;

      shell->pivots.ipvtpt =  shell->ipvtpt;

      shell->scal.iscale   =  shell->iscale;
      shell->scal.isclag   =  shell->isclag;

      shell->newsrc.iter   =  shell->iter;

      shell->limits.itlim  =  shell->itlim;
      shell->pardat.kderiv =  shell->kderiv;

      shell->memory.lbinv  =  shell->lbinv;    shell->memory.lirn   =  shell->lirn;
      shell->memory.lgrad  =  shell->lgrad;    shell->memory.lmem   =  shell->lmem;
     /*  np1 needed in grdump set to n+1??? */
      shell->memory.np1    =  (shell->n)+1;

      shell->dimen.m       =  shell->m;        shell->dimen.nnbmax  =  shell->nnbmax;
      shell->dimen.n       =  shell->n;        shell->dimen.mpnbmx  =  shell->mpnbmx;
      shell->dimen.mp1     =  shell->mp1;      shell->dimen.npnbmx  =  shell->npnbmx;
      shell->dimen.npmp1   =  shell->npmp1;    shell->dimen.nrtot   =  shell->nrtot;
      shell->dimen.nbmax   =  shell->nbmax;


      shell->cgbk.mcgm1    =  shell->mcgm1;    shell->cgbk.memcg    =  shell->memcg;
      shell->cgbk.modcg    =  shell->modcg;

      shell->lincnt.nfix   =  shell->nfix;     shell->lincnt.nlin   =  shell->nlin;
      shell->lincnt.nnlin  =  shell->nnlin;

      shell->equblk.nnlequ =  shell->nnlequ;

      shell->nzerog.nzgrad =  shell->nzgrad;   shell->nzerog.nzlin  =  shell->nzlin;
      shell->nzerog.nznlin =  shell->nznlin;

      shell->jgrow.maxgro  =  shell->jacobian_growth_allowance;

      shell->bind.nbc      =  shell->nbc;
      shell->nintbk.ninf   =  shell->ninf;
      shell->chq.maxtab    =  shell->maxtab;

   shell->counts.ncalls =  shell->nnewton;       shell->counts.nstepc =  shell->nnewton_stepc;
   shell->counts.nit    =  shell->nnewton_itns;  shell->counts.nftn   =  shell->ngcomp;
   shell->counts.nnfail =  shell->nnewton_fail;
   shell->counts.nbs    =  shell->nbs;

   shell->cbscnt.ncbs   =  shell->ncbs;     shell->cbscnt.nicond =  shell->nicond;
   shell->cbscnt.nser   =  shell->nser;     shell->cbscnt.nsing  =  shell->nsing;
   shell->cbscnt.nrser  =  shell->nrser;    shell->cbscnt.ifrcnt =  shell->ifrcnt;
   shell->cbscnt.nreinv =  shell->nreinv;   shell->cbscnt.nsame  =  shell->nsame;

   shell->misc.nsear  =  shell->nsear;    shell->misc.maxh     =  shell->maxh;

   shell->iters.nph1ls  =  shell->nph1ls;  shell->iters.ndeg    =  shell->ndeg;

  /* ints */

   shell->pivots.fixpiv   =    shell->fixpiv;
   shell->optblk.hrdbnd   =    shell->hrdbnd;
   shell->cgbk.hscale     =    shell->hscale;

   shell->optblk.maxim    =    shell->maxim; shell->optblk.multsb   =    shell->multsb;
   shell->optblk.useph0   =    shell->useph0;
   shell->optblk.gfeas    =    shell->gfeas;

#ifdef FILE_IO_ENABLED
     shell->io_info.lsgrg_ioout  = shell->io_info.lsgrg_ioout;
     shell->io_info.lsgrg_ioerr  = shell->io_info.lsgrg_ioerr;
     shell->io_info.lsgrg_ioterm = shell->io_info.lsgrg_ioterm;
#endif

} /* end of lsgrg_put_globals */

void LSGRGCLASS lsgrg_get_globals(LsgrgInfo *shell)
{
/*------------------------------------------------------------------*/
/* lsgrg_get_globals() copies global values from alg to shell       */
/*------------------------------------------------------------------*/
    shell->nparsh        +=         shell->counts.ngrad;
    shell->nbc            =         shell->bind.nbc;
    shell->ninf           =         shell->nintbk.ninf;
    shell->nnewton        =         shell->counts.ncalls;
    shell->nnewton_stepc  =         shell->counts.nstepc;
    shell->nnewton_itns   =         shell->counts.nit;
    shell->nnewton_fail   =         shell->counts.nnfail;
    shell->nbs            =         shell->counts.nbs;

    shell->ncbs           =         shell->cbscnt.ncbs;
    shell->nicond         =         shell->cbscnt.nicond;
    shell->nser           =         shell->cbscnt.nser;
    shell->nsing          =         shell->cbscnt.nsing;
    shell->nrser          =         shell->cbscnt.nrser;
    shell->ifrcnt         =         shell->cbscnt.ifrcnt;
    shell->nreinv         =         shell->cbscnt.nreinv;
    shell->nsame          =         shell->cbscnt.nsame;

    shell->nsear          =         shell->misc.nsear;
    shell->nph1ls         =         shell->iters.nph1ls;
    shell->ndeg           =         shell->iters.ndeg;
    shell->nb             =         shell->nintbk.nb;
    shell->nsuper         =         shell->nintbk.nsuper;
/* ** 5/01 jcp ** transfer colmax back to shell so that  */
/*  correct value is maintained between grgitn entries   */
/*  when jacobian growth allowance is exceeded and grad  */
/*  is reallocated                                       */

       shell->colmax      =         shell->cmax.colmax;


/* transfer timing stats kept by alg routines */
  lsgrg_timestats(shell,T_CONSBS   ,shell->contim.cbtime);
  lsgrg_timestats(shell,T_PHASE0   ,shell->ph0tim.tph0);
  lsgrg_timestats(shell,T_PHASE1   ,shell->ph0tim.tph1);
  lsgrg_timestats(shell,T_NEWT     ,shell->nwtim.tnewt);
  lsgrg_timestats(shell,T_NEWTONLY ,shell->nwtim.tnewtx);
  lsgrg_timestats(shell,T_SEARCH   ,shell->temp.tsear);
  lsgrg_timestats(shell,T_INVERT   ,shell->temp.tfact);

   return;
}
#include "lsinfo.h"
#ifdef IO_ENABLED
   #define IOINFO &(info->io_info)
#endif
/***********************************************************************
 **                 FILE:              GRGOUTR.C                       *
 **                 AUTHOR:            STUART SMITH                    *
 **    Original     TRANSLATION TO C:  ALLAN D. WAREN                  *
 **                 LAST CHANGE:       Oct 1998                        *
 **                 Rewritten for lsgrgc 2.0 John Plummer 4/98         *
 **                                                                    *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **  ** 3/98 John Plummer ** moved echog3 to this file to put all      *
 **                          output stuff here and renamed to          *
 **                          grgechoparms. also updated tablin and     *
 **                          outres to use row/col name conventions    *
 **                          of revised user interface                 *
 **                                                                    *
 **  ** 4/98 jcp  reworked io to conform to revised io formats         *
 **        .  fprinf's ifdefed on FILE_IO_ENABLED                      *
 **        .  revised for 1-based arrays                               *
 **                                                                    *
 **      THIS FILE CONTAINS THE FOLLOWING ROUTINES FOR LSGRG2, A       *
 **  LARGE SCALE IMPLEMENTATION OF GRG2:                               *
 **                                                                    *
 **  ROUTINE             DESCRIPTION                      LAST UPDATE  *
 ** ---------           -------------                    ------------- *
 **  lsgrg_outres   Prints out final results                 may 1998  *
 **  lsgrg_tablin   Prints fcns/var into at initial point    may 1998  *
 **  lsgrg_echoparms    echos problem parameters             mar 1998  *
 **                  echos global parameters                 mar 1998  *
 **  lsgrg_printgrad  prints jacobian structure at           may 1998  *
 **                   current point                                    *
 **  lsgrg_printgrad0 prints jacobian structure assuming               *
 **                   0-based arrays                                   *
 **                                                                    *
 **  lsgrg_prt_iarray dumps an integer array                           *
 **                                                                    *
 **                                                                    *
 ***********************************************************************
 * */

void LSGRGCLASS lsgrg_outres(LsgrgInfo *info, int return_code)
{
/*-----------------------------------------------------------------*/
/*    all problem data needed here is global                       */
/*  lsgrg_varnames and lsgrg_rownames start at index 1             */
/*  no code here if io and fileio not enabled                      */
/*-----------------------------------------------------------------*/
#ifdef IO_ENABLED
   double plinfy = info->plinfy;
   char bnd[9], nearb[3], status[9], tname[9];
   long int i, k, ni, nparsh_setup;
   double bound, slack, slackl, slacku, ts, xncall, xnit;
   char aobj[9]     = "  OBJ   ";
   char ignore[9]   = "IGNORED ";
   char lowbnd[9]   = "LOWERBND";
   char uppbnd[9]   = "UPPERBND";
   char equal[9]    = "EQUALITY";
   char free[9]     = "  FREE  ";
   char nobnd[9]    = "NO BOUND";
   char supbas[9]   = "SUPBASIC";
   char nonbas[9]   = "NONBASIC";
   char basic[9]    = " BASIC  ";
   char viol[9]     = "VIOLATED";
   char blank[9]    = "        ";
   char fixed[9]    = " fixed  ";
   char tcode[8][33]      ={  "Kuhn-Tucker Conditions Satisfied",
                           "Fractional Change in Objective  ",
                           "All Remedies Failed             ",
                           "Line Search Limit Exceeded      ",
                           "Solution Unbounded              ",
                           "Feasible Point Not Found        ",
                           "GRG Runtime Error - No Solution ",
                           "User Termination at Feasible pt ",
                            };
    char blank2[3]   = "  ";
    char lower[3]    = ":l";
    char upper[3]    = ":u";
    double big       = 1.0e20;
    double zero      = 0.0e0;
     int newpage;
     int print_row_col_info /* ,print_run_stats,print_summary */;
/*-------------------------------------------------------------------*/
/*     have already changed sign of g(nobj) in grgitn                */
/*     don't negate gradf and u -- already done in grgitn            */
/*-------------------------------------------------------------------*/
    print_row_col_info = info->ipr3 > -1 ||
            info->outres_print_level == OUTRES_PRINT_ALL ||
            info->outres_print_level == OUTRES_PRINT_ROWCOL;


   if(print_row_col_info) {  /*  detailed output */

        for( i = 1; i <= info->mp1; i++ )
                info->gg[i] = plinfy;
/*--------------------------------------------------------------------*/
/*       mark the lang. mults of the non-binding constraints          */
/*       they will not be printed                                     */
/*--------------------------------------------------------------------*/
        if( info->nbc > 0 )  {
           for( i = 1; i <= info->nbc; i++ ) {
                k = info->ibc[i];
                info->gg[k] = info->u[k];
           }
        }
/*--------------------------------------------------------------------*/
/*                            output rows                             */
/*--------------------------------------------------------------------*/
        newpage = 1;
        lsgrg_output_headers(info,FINAL_RESULTS,FUNCTIONS,newpage);
        for( i = 1; i <= info->mp1; i++ ) {
             if( i%40==0 )
                 lsgrg_output_headers(info,FINAL_RESULTS,FUNCTIONS,newpage);

              ni = info->n + i;
              strcpy( status, free );
              slacku = info->ub[ni] - info->g[i];
              slackl = info->g[i] - info->alb[ni];
              strcpy( nearb, upper );
              slack = slacku;
              bound = fabs( info->ub[ni] );
              if( info->lsgrg_rownames != NULL )
                  strncpy( tname, info->lsgrg_rownames[i], 8 );
              else
                  strcpy( tname, blank );

                if( slackl < slacku )  {
                    slack = slackl;
                    bound = fabs( info->alb[ni] );
                    strcpy( nearb, lower );
                }

                if( info->alb[ni] == info->ub[ni] ) {     /* equality row */
                    strcpy( status, equal );
                    strcpy( nearb, blank2 );
                    if( slack < (-info->epnewt*(info->rtf + bound)) )
                            strcpy( status, viol );
                }
                else {     /* non equality row */
                       if( slack < (-info->epnewt*(info->rtf + bound)) )
                           strcpy( status, viol );
                       else {
                             if(slack <= (info->epnewt*(info->rtf+bound)) ) {
                                if( slack == slackl )
                                    strcpy( status, lowbnd );
                                if( slack == slacku )
                                    strcpy( status, uppbnd );
                             }
                       }
                       if( info->istat[i] == 0 ) {
                          strcpy( status, ignore );
                          strcpy( nearb, blank2 );
                       }
                       if( i == info->nobj ) {
                           strcpy( status, aobj );
                           sprintf(info->lsgrg_msgbuffer, "\n"
/*                           " %5ld %8.8s %12.5e %12.5e %8.8s\n", */
                             " %5ld %8.8s %12.4e %12.4e %8.8s\n",
                              i, tname, info->g0[i], info->g[i], status );
                           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
                           continue;
                       }
                } /* end non-equality row else */

                sprintf(info->lsgrg_msgbuffer, "\n"
/*                  " %5ld %8.8s %12.5e %12.5e %8.8s %10.3e%2.2s", */
                    " %5ld %8.8s %12.4e %12.4e %8.8s %10.2e%2.2s",
                        i, tname, info->g0[i], info->g[i], status, slack,
                        nearb);
                lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

                if( info->gg[i] != plinfy ) {
                   sprintf(info->lsgrg_msgbuffer," %12.4e",info->gg[i] );
                   lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
                }
        }  /*  end of row loop */
/*--------------------------------------------------------------------*/
/*                            output columns                          */
/*                                                                    */
/*  icols is used to temporarily store the basic/superbasic           */
/*  status of the variables                                           */
/*  0 = nonbasic, 1 = superbasic, 2 = basic  3 = fixed                */
/*  gradfp is used to temporarily store reduced cost                  */
/*--------------------------------------------------------------------*/
        newpage = 0;
        if( info->mp1 + info->n > 34 ) newpage = 1;
        lsgrg_output_headers(info,FINAL_RESULTS,VARIABLES,newpage);

        for( i = 1; i <= info->n; i++ ) {
             info->icols[i] = 0;
             if( info->ivstat[i] < 0 )
                 info->icols[i] = 3;
             info->gradfp[i] = plinfy;
        }

/*  flag superbasics   */

        if( info->nsuper >  0 ) {
           for( i = 1; i <= info->nsuper; i++ ) {
                k = info->inbv[i];
                if( k <= info->n && info->iub[i] == 0 )
                        info->icols[k] = 1;
                }
        }
/*  flag basics */

        if( info->nb > 0 ) {
            for( i = 1l; i <= info->nb; i++ ) {
                 k = info->ibv[i];
                if( k <= info->n )
                    info->icols[k] = 2;
            }
        }
/*  extract reduced gradient components */

        for( i = 1; i <= info->n; i++ ) {
             k = info->inbv[i];
             if( k <= info->n )
                 info->gradfp[k] = info->gradf[i];
        }

/*  determine status and output line for each variable */


        for( i = 1; i <= info->n; i++ ) {
             if( i%40 ==0 )  {
                 newpage = 1;
                 lsgrg_output_headers(info,FINAL_RESULTS,VARIABLES,newpage);
             }

             strcpy( status, nonbas );
             if( info->icols[i] == 1 )
                  strcpy( status, supbas );
             if( info->icols[i] == 2 )
                  strcpy( status, basic );
             if( info->icols[i] == 3 )
                  strcpy( status, fixed );

             slackl = info->x[i] - info->alb[i];
             slacku = info->ub[i] - info->x[i];
             slack = slackl;
             strcpy( nearb, lower );
             bound = fabs( info->alb[i] );
             if( slackl >= slacku )  {
                 slack = slacku;
                 strcpy( nearb, upper );
                 bound = fabs( info->ub[i] );
             }
             strcpy( bnd, blank );
             if( info->alb[i] == info->ub[i] )
                  strcpy( bnd, fixed );

             if( slack > big )
                 strcpy( bnd, nobnd );

             if( fabs( slack ) <= (info->epboun*(info->rtf + bound)) ) {
                 if( slack == slackl )
                     strcpy( bnd, lowbnd );
                 if( slack == slacku )
                     strcpy( bnd, uppbnd );
             }
             if( info->lsgrg_varnames != NULL )
                 strncpy( tname, info->lsgrg_varnames[i], 8 );
             else
                 strcpy( tname, blank );

/*  print output line */

             sprintf(info->lsgrg_msgbuffer,
              "%4ld %8.8s%12.4e%12.4e %8.8s" ,
                 i, tname, info->x0[i], info->x[i], status );
             lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

                if( info->gradfp[i] != plinfy ) {
                   sprintf(info->lsgrg_msgbuffer,
                     "   %8.8s %12.4e \n", bnd ,info->gradfp[i]);
                   lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
                }
                else {
                   sprintf(info->lsgrg_msgbuffer,
                          "%10.2e%2.2s\n",
                           slack, nearb );
                   lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
                }

        } /* end of column print loop */

     if(info->outres_print_level == OUTRES_PRINT_ROWCOL) return;
/*---------------------------------------------------------- */
        sprintf(info->lsgrg_msgbuffer,
            "\n\n\n\nRun Statistics \n %s\n", info->title);
        lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

    }  /*  close detailed print (if  print_row_col_info) */
/*---------------------------------------------------------- */
/*  end of detailed print                                    */
/* ** 3/20/02 jcp ** if we are only printing row/col info    */
/*   return
/*---------------------------------------------------------- */

     if( info->ipr3 <= -1 ) {
         sprintf(info->lsgrg_msgbuffer, "\nTermination Code:   %32.32s\n",
           tcode[return_code] );
         lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
         sprintf(info->lsgrg_msgbuffer, "\nFinal Objective Value =  %15.8e\n",
             info->g[info->nobj] );
         lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
     }
     ts = zero;
     if( info->nnewton > 0 ) {
         xnit = info->nnewton_itns;
         xncall = info->nnewton;
         ts = xnit/xncall;
     }
/*  info->ngcomp is total gcomp calls. back out gradient calls */
/*   there are 2 parsh calls for full variable set in memsetup     */

     nparsh_setup = 2;
     i = info->ngcomp;
     if(info->kderiv < 2)
        i -= (info->nnlin*(info->nparsh-nparsh_setup))
                           + (info->n*nparsh_setup) ;


    if( info->optblk.useph0 )  {
       sprintf(info->lsgrg_msgbuffer,"\n"
    ". Number of Phase-0 Newton Iterations                             = %6ld\n"
    "   Initial Number of Phase-0 Infeasibilities                      = %6ld\n"
    "   Number of Phase-0 Constraints Which Failed to Converge         = %6ld\n",
                 info->iters.nph0it, info->nph0.ninf0, info->nph0.ndrop);
       lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
    }
  sprintf(info->lsgrg_msgbuffer,
    ". Number of One-Dimensional Searches                              = %6ld\n"
    "   Number of Phase-1 One-Dimensional Searches                     = %6ld\n"
    "   Number of Degenerate Iterations                                = %6ld\n",
                 info->nsear, info->nph1ls, info->ndeg );
  lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
  sprintf(info->lsgrg_msgbuffer,"\n"
    " Function Calls =%5ld  Gradient Calls = %5ld\n"
    "    Total Function Calls (inc. For Gradient)                      = %6ld\n",
          i, info->nparsh, info->ngcomp );
  lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
  sprintf(info->lsgrg_msgbuffer,
    "    Nbr Gradient Calls During Setup Phase                         = %6ld\n",
         nparsh_setup);
  lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

  sprintf(info->lsgrg_msgbuffer,"\n"
    ". Newton Calls = %5ld  Newton Iterations = %5ld  Average        = %8.2f\n",
              info->nnewton, info->nnewton_itns, ts );
  lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
  sprintf(info->lsgrg_msgbuffer,
    "   Number of Times Newton Failed to Converge                      = %6ld\n "
    "   Times Stepsize Cut Back Due to Newton Failure                 = %6ld\n",
              info->nnewton_fail, info->nnewton_stepc );
  lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
  sprintf(info->lsgrg_msgbuffer,
    ". Number of Times Basic Variable Violated a Bound                 = %6ld\n",
      info->nbs );
  lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
  sprintf(info->lsgrg_msgbuffer,
    ". Number of Consbs calls                                          = %6ld\n"
    "    Number of Times Consbs Used Search Mode                       = %6ld\n"
    "    Number of Times Consbs Used Restricted Search Mode            = %6ld\n"
    "    Number of Times Consbs Used Re-inversion Mode                 = %6ld\n",
          info->ncbs, info->nser, info->nrser, info->nreinv );
  lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
  sprintf(info->lsgrg_msgbuffer,
    "    Number of Times an Ill-Conditioned Basis Forced a Search Mode = %6ld\n"
    "    Number of Times an Inversion Failure Forced a Search Mode     = %6ld\n"
    "    Number of Times Fractional Change Forced Search Mode          = %6ld\n"
    "    Number of Times CONSBS Chose Previous Basis in Search Mode    = %6ld\n",
              info->nicond, info->nsing, info->ifrcnt, info->nsame );
  lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
#endif
/* ------------------------------------------------------------ */
        return;

} /* end of function lsgrg_outres */


void LSGRGCLASS lsgrg_tablin(LsgrgInfo *info)
{

#ifdef IO_ENABLED

   char type[5], status[5],  name[9];
   long int i, ni;
   double gi, plinfy = info->plinfy;
/*--------------------------------------------------------------------*/
/*     print initial row values                                       */
/*  lsgrg_varnames and lsgrg_rownames start at index 1                */
/*--------------------------------------------------------------------*/
   lsgrg_output_headers(info,INITIAL_VALUES,FUNCTIONS,0);

    for( i = 1; i <= info->mp1; i++ ) {
        ni = info->n + i;
        gi = info->g[i];
        strcpy( name, "   " );
        if(info->lsgrg_rownames != NULL) {
           strncpy( name,info->lsgrg_rownames[i], 8 );
           lsgrg_formatname(name,8);
        }

/*  assign row type string */
         strcpy(type," ");
        if( info->alb[ni] == -plinfy && info->ub[ni] == plinfy )
            strcpy( type, "free" );
        else {
           if( info->alb[ni] == info->ub[ni] )
               strcpy( type, " eq " );
        }
        if(strcmp(type," ")==0) {
           if( info->alb[ni] != -plinfy && info->ub[ni] != plinfy )
               strcpy( type, "rnge" );
           else if( info->alb[ni] != -plinfy )
                    strcpy( type, " ge " );
                else if( info->ub[ni] != plinfy )
                            strcpy( type, " le " );
        }
        if( i == info->nobj )
            strcpy( type,"obj " );

/*  assign row status strings */

       strcpy( status, "    " );
       if( fabs( gi - info->ub[ni] ) < info->epnewt )
           strcpy( status, " ul " );
       if( fabs( info->alb[ni] - gi ) < info->epnewt )
           strcpy( status, " ll " );
       if( fabs( gi - info->ub[ni] ) < info->epnewt &&
           fabs( info->alb[ni] - gi ) < info->epnewt   )
              strcpy( status, " eq " );
       if( gi > info->ub[ni] + info->epnewt || gi < info->alb[ni] - info->epnewt )
           strcpy( status, "viol" );

/*  print row detail */

       sprintf(info->lsgrg_msgbuffer, "%4ld  %8s   %4s    %4s %14.6e",
          i, name, status, type, info->g[i] );
       lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

/*  print lower bound */
       if( info->alb[ni] <= -plinfy)
          lsgrg_msg(IOINFO,"     <none>   ");
       else {
          sprintf(info->lsgrg_msgbuffer, " %14.6e",
            info->alb[ni]);
          lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
       }

/*   print upper bound */
       if( info->ub[ni] >= plinfy )
          lsgrg_msg(IOINFO,"     <none>   \n");
       else {
          sprintf(info->lsgrg_msgbuffer, " %14.6e\n",
            info->ub[ni]);
          lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
       }

    } /* end of row print loop */
/*--------------------------------------------------------------------*/
/*     print initial var values                                       */
/*  lsgrg_varnames and lsgrg_rownames start at index 1                */
/*--------------------------------------------------------------------*/
   lsgrg_output_headers(info,INITIAL_VALUES,VARIABLES,0);

   strcpy(type, "    ");
   for( i = 1; i <= info->n; i++ ) {
       if( labs( info->ivstat[i] ) == 1 )
            strcpy( type, "nlin" );
       else
            strcpy( type, "lin " );
       strcpy( name, "    " );
       if( info->lsgrg_varnames != NULL ) {
           strncpy( name, info->lsgrg_varnames[i], 8 );
           lsgrg_formatname(name,8);
       }

       strcpy( status, "    " );
       if( info->ivstat[i] < 0 || info->alb[i] == info->ub[i] )
           strcpy( status, " fx " );
       else
          if( info->alb[i] == info->x[i] )
              strcpy( status, " ll " );
          else
              if( info->ub[i] == info->x[i] )
                  strcpy( status, " ul " );
          else
              if( info->alb[i] == -plinfy && info->ub[i] == plinfy )
                strcpy( status, "free" );

       sprintf(info->lsgrg_msgbuffer,
        "%4ld  %8s   %4s   %4s %14.6e",
                 i, name, status, type, info->x[i] );
       lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

/* print lower bound */
       if( info->alb[i] <= -plinfy)
          lsgrg_msg(IOINFO,"     <none>   ");
       else {
          sprintf(info->lsgrg_msgbuffer, " %14.6e",
            info->alb[i]);
          lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
       }
/* print upper bound */
       if( info->ub[i] >= plinfy )
          lsgrg_msg(IOINFO,"     <none>   \n");
       else {
          sprintf(info->lsgrg_msgbuffer, " %14.6e\n",
            info->ub[i]);
          lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
       }

   } /* end of variable print loop */

   lsgrg_msg(IOINFO,"\n" );
#endif
   return;
} /* end of function lsgrg_tablin  */

void  LSGRGCLASS lsgrg_echoparms(LsgrgInfo *info)
{
/*--------------------------------------------------------------------*/
/*   prints problem parameter values                                  */
/*--------------------------------------------------------------------*/
#ifdef IO_ENABLED

   char hscind[18], oname[11];
   char sense[10];
   char *cgtype[7]={" ","Fletcher-Reeves",
                        "Polak-Ribiere",
                        "Perry",
                        "1-Step DFP",
                        "1-Step BFGS",
                        "Limited Memory BFGS"};

   char *derivs[4]={"Forward Differences",
                        "Central Differences",
                        "Analytic (user parsh)",
                        "GAMS Interpreter" };

   char *iests[3]={"Tangent Vectors",
                      "Quadratic Interpolation"};
   char blank8[9] = "        ";
   char releasedate[9] = {RELEASEDATE};
   char *months[] = { " ","January","February","March","April","May",
                      "June","July","August","September","October",
                      "November","December"};
   char s_release[20], s_mm[3]="xx",*p,*p2;
   int i,mm;
    char isDefault[]="(default)";
    char isNotDefault[] = " ";

/*--------------------------------------------------------------------*/
/*         build release date string                                  */
/*--------------------------------------------------------------------*/
        p = s_mm;
        if(releasedate[0] != '0') *p++ = releasedate[0];
        *p++ = releasedate[1];
        mm = atoi(s_mm);
        if(mm > 0 && mm < 13) {
           strcpy(s_release,months[mm]);
           p = s_release;
           while (*p) p++;
           *p++ = ' ';
           if(releasedate[2] != '0')
               *p++ = releasedate[2];
           *p++ = releasedate[3];
           *p++ = ',';
           *p++ = ' ';
           for(i=4;i<=7;++i)
               *p++ = releasedate[i];
           *p = '\0';
        }
        else
           strcpy(s_release,"XXXXXXXXX");
/*--------------------------------------------------------------------*/
/*         echo problem specifications                                */
/*--------------------------------------------------------------------*/
        strcpy( oname, "-unnamed-" );
        strcpy( hscind, "Will Not be Used " );
        if( info->hscale )
                strcpy( hscind, "Will Be Used     " );
        if( info->usernames ) {
           if( strcmp(info->lsgrg_rownames[info->nobj],blank8) != 0 )
               strcpy( oname, info->lsgrg_rownames[info->nobj] );
        }

     lsgrg_msg(IOINFO,"\n\n");
     for(i=1;i<=78;++i) lsgrg_msg(IOINFO,"-");
     sprintf(info->lsgrg_msgbuffer,"\n"
       "                         Lsgrgc 3.0  Release of %s\n"
       "                         Problem Parameters and Options\n"
       "                         ------------------------------\n"
       "Problem Name: %s\n", s_release,info->title );
     lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

     sprintf(info->lsgrg_msgbuffer,
       ". Number of Variables                             = %7ld\n"
       "    Number of Fixed Variables                     = %7ld\n"
       "    Number of Linear Variables                    = %7ld\n"
       "    Number of Nonlinear Variables                 = %7ld\n"
       ". Number of Functions                             = %7ld\n",
           info->nvars, info->nfix, info->nlin, info->nnlin, info->nrows );
     lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

      sprintf(info->lsgrg_msgbuffer,
       "    Number of Nonlinear Functions                 = %7ld\n"
       "    Number of Linear Functions                    = %7ld\n\n",
                 info->nnlequ, (info->nrows - info->nnlequ) );
      lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      sprintf(info->lsgrg_msgbuffer,
       ". Total Number of Nonzero Jacobian Elements       = %7ld\n"
       "    Number of Nonzero Nonlinear Jacobian Elements = %7ld\n"
       "    Number of Nonzero Linear Jacobian Elements    = %7ld\n\n",
                 info->nzgrad, info->nznlin, info->nzlin );
      lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      sprintf(info->lsgrg_msgbuffer,
       ". Space Reserved for Hessian Has Dimension        = %7ld\n"
       ". Limit on Binding Constraints                    = %7ld\n"
       ". Length of Chuzq Tabu List                       = %7ld\n\n",
                 info->maxr, info->maxb, info->maxtbu );
      lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

      strcpy(sense,"Minimized");
      if( info->maxim ) strcpy(sense,"Maximized");

      sprintf(info->lsgrg_msgbuffer,
        ". Objective Row is row %5ld [%s]  (%s)\n"
        ". Intitial Estimates Will be Obtained Via %s\n\n",
                 info->nobj, oname, sense, iests[info->iquad] );
      lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      sprintf(info->lsgrg_msgbuffer,
        ". Derivative Computation Method: %s\n",
                 derivs[info->kderiv] );
      lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      if( info->kderiv < 2 ) {
          sprintf(info->lsgrg_msgbuffer,
            "    Finite Difference Perturbation: %10.3e\n",info->pstep );
      lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      }
      if( info->iscale == 1 ) {
          sprintf(info->lsgrg_msgbuffer,
           ". The Jacobian Matrix Will Be Scaled at the Initial Point\n" );
       lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      }
      if( info->isclag != 0l ) {
          sprintf(info->lsgrg_msgbuffer,
           ". The Jacobian Matrix Will be Rescaled After Every"
           " %5ld Line Searches\n", info->isclag );
          lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      }
      if( info->multsb )
            lsgrg_msg(IOINFO,
             ". Superbasics Will be Projected Onto Their Bounds"
             " During the Line Search \n" );
      else
            lsgrg_msg(IOINFO,
             ". The Line Search Will be Terminated If a Superbasic"
             " Variable Reaches a Bound \n" );

      lsgrg_msg(IOINFO,". Basic Variables ");
      if( info->hrdbnd )
          lsgrg_msg(IOINFO,"Will *NOT* ");
      else
          lsgrg_msg(IOINFO,"*Will* ");
      lsgrg_msg(IOINFO,
         "Be Allowed to Violate Their Bounds During\n"
          "     the Line Search \n" );

      if( info->fixpiv ) {
          sprintf(info->lsgrg_msgbuffer,
              ". Threshold Pivot Tolerance Will be Fixed at the Value %12.6e\n",
                 info->xpvpct );
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      }
      if( info->maxh == 0 )  {
          sprintf(info->lsgrg_msgbuffer,
            ". Conjugate Gradient Method #%1ld: %s Will be Used\n",
                  info->modcg, cgtype[info->modcg] );
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
          if( info->modcg == 6 ) {
              sprintf(info->lsgrg_msgbuffer,
               "   Limited Memory Length = %5ld\n   Automatic Hessian"
               " Scaling %17.17s\n", info->memcg, hscind );
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
          }
        }
        else {
               sprintf(info->lsgrg_msgbuffer,
                 ". BFGS Will be Used if Nbr Superbasics =< %4ld\n",info->maxh );
               lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
               if( info->maxh < info->nvars ) {
                  sprintf(info->lsgrg_msgbuffer,
                   "  Conjugate Gradient Method #%1ld: %s"
                   " Will be Used Otherwise\n",
                       info->modcg, cgtype[info->modcg] );
                  lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
                  if( info->modcg == 6 ) {
                      sprintf(info->lsgrg_msgbuffer,
                        "   Limited Memory Length = %5ld\n   Automatic"
                        " Hessian Scaling %s\n",
                          info->memcg, hscind );
                      lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
                  }
               }
        }
        lsgrg_msg(IOINFO,"\n. Phase 0 ");
        if(info->useph0)
            lsgrg_msg(IOINFO," Will Be Used");
        else
            lsgrg_msg(IOINFO," Will not be Used\n");

        lsgrg_msg(IOINFO,
         "\nTolerances, Limits And Print Options"
         "\n ------------------------------------\n" );
        if(info->maxtime == LSGRG_TIME_LIMIT) p = isDefault;
        else p = isNotDefault;
        if(info->limser == LSGRG_ITERATION_LIMIT) p2 = isDefault;
        else p = isNotDefault;
        sprintf(info->lsgrg_msgbuffer,
          "\n Iteration Limit         = %20d %s"
          "\n Time Limit (sec)        = %20.2f %s\n",
        info->limser,p2,info->maxtime,p);
        lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

        sprintf(info->lsgrg_msgbuffer,
          "\n Machine Precision (eps) = %20.14e  (computed)\n",
          info->eps);
        sprintf(info->lsgrg_msgbuffer,
          " EPNEWT =%12.4e EPINIT =%12.4e EPSTOP =%12.4e\n"
          " EPBOUN =%12.4e EPDEG  =%12.4e EPPIV  =%12.4e\n"
          " PH1EPS =%12.4e CNDTOL =%12.4e PSTEP  =%12.4e\n"
          " FUNPR  =%12.4e PIVPCT =%12.4e AIJTOL =%12.4e\n"
          " PLINFY =%12.4e\n\n"
          " NSTOP  =%5ld  ITLIM =%5ld  IDGLIM =%5ld  SEARCH =%10ld\n"
          " BVBND= %5ld\n"
          " IPR =%4ld  PN4 =%4ld  PN5 =%4ld\n PN6 =%4ld  PER =%4ld\n",
              info->epnewt, info->epinit, info->epstop, info->epboun,
              info->epdeg,  info->epspiv, info->ph1eps,
              info->cndtol, info->pstep,  info->funpr, info->xpvpct,
              info->xajtol, info->plinfy, info->nstop,
              info->itlim,  info->idglim, info->limser, info->ibvblm,
               info->ipr, info->ipn4, info->ipn5, info->ipn6, info->iper );
         lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

        lsgrg_memstats(info);
        sprintf(info->lsgrg_msgbuffer, "\n",
          "     Basis Factors         %10d               %10d\n",
                  info->lbinv,(info->lbinv*sizeof(double))  );
        lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

        if( info->dbg.globals || info->ipr >= 5 )
            lsgrg_echo_globals(info,"After Problem Setup");
#endif
        return;

} /* end of function lsgrg_echoparms */

void LSGRGCLASS lsgrg_echo_globals(LsgrgInfo *info,char *msg)
{
/*-------------------------------------------------------------------*/
/*   for development/debug purposes, echo global parameters          */
/*-------------------------------------------------------------------*/
#ifdef IO_ENABLED


int i;

   if(! info->dbg.globals) return;
   sprintf(info->lsgrg_msgbuffer,
    "\n\n                    LSGRGC Shell Global Variables\n "
    "\n  Called From: %s\n",msg     );
   lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

   for(i=1;i<=78;++i) lsgrg_msg(IOINFO,"=");
   lsgrg_msg(IOINFO," ..... #define constants .....");
   lsgrg_lprint(info,"MAXOPTIONS"           ,MAXOPTIONS        );
   lsgrg_lprint(info,"TITLELENGTH"          ,TITLELENGTH       );
   lsgrg_lprint(info,"MAXVARS"              ,MAXVARS);
   lsgrg_lprint(info,"MAXROWS"              ,MAXROWS);
   lsgrg_lprint(info,"NAMELENGTH"           ,NAMELENGTH);
   lsgrg_lprint(info,"MIN_LGRAD"            ,MIN_LGRAD);
   lsgrg_lprint(info,"J_GROWTH_ALLOWANCE"   ,J_GROWTH_ALLOWANCE);
   lsgrg_lprint(info,"MINLBINV"             ,MINLBINV);
   lsgrg_lprint(info,"MSGBUFFER_LENGTH"     ,MSGBUFFER_LENGTH  );
   lsgrg_dprint(info,"BINVFACTOR"           ,BINVFACTOR);
   lsgrg_dprint(info,"PLINFY"               ,PLINFY    );
   lsgrg_dprint(info,"PLZERO"               ,PLZERO    );

   lsgrg_msg(IOINFO,"\n..... global variables .....");
   lsgrg_lprint(info,"usernames"            ,info->usernames );
    lsgrg_lprint(info,"inprnt"         ,              info->inprnt);
    lsgrg_lprint(info,"otprnt"         ,              info->otprnt);
    lsgrg_lprint(info,"nbc"            ,                 info->nbc);
    lsgrg_lprint(info,"ifrcnt"         ,              info->ifrcnt);
    lsgrg_lprint(info,"ncbs"           ,                info->ncbs);
    lsgrg_lprint(info,"nicond"         ,              info->nicond);
    lsgrg_lprint(info,"nreinv"         ,              info->nreinv);
    lsgrg_lprint(info,"nrser"          ,               info->nrser);
    lsgrg_lprint(info,"nsame"          ,               info->nsame);
    lsgrg_lprint(info,"nser"           ,               info->nser);
    lsgrg_lprint(info,"nsing"          ,               info->nsing);
    lsgrg_lprint(info,"mcgm1"          ,               info->mcgm1);
    lsgrg_lprint(info,"memcg"          ,               info->memcg);
    lsgrg_lprint(info,"modcg"          ,               info->modcg);
    lsgrg_iprint(info,"hscale"         ,              info->hscale);
    lsgrg_lprint(info,"maxtab"         ,              info->maxtab);
    lsgrg_lprint(info,"colmax"         ,              info->colmax);
    lsgrg_lprint(info,"nbs"            ,                 info->nbs);
    lsgrg_lprint(info,"nnewton"        ,              info->nnewton);
    lsgrg_lprint(info,"ngcomp"         ,                info->ngcomp);
    lsgrg_lprint(info,"nparsh"         ,               info->nparsh);
    lsgrg_lprint(info,"nnewton_itns"   ,               info->nnewton_itns);
    lsgrg_lprint(info,"idglim"         ,              info->idglim);
    lsgrg_lprint(info,"m"              ,                   info->m);
    lsgrg_lprint(info,"mp1"            ,                 info->mp1);
    lsgrg_lprint(info,"mpnbmx"         ,              info->mpnbmx);
    lsgrg_lprint(info,"n"              ,                   info->n);
    lsgrg_lprint(info,"nbmax"          ,               info->nbmax);
    lsgrg_lprint(info,"nnbmax"         ,              info->nnbmax);
    lsgrg_lprint(info,"npmp1"          ,               info->npmp1);
    lsgrg_lprint(info,"npnbmx"         ,              info->npnbmx);
    lsgrg_lprint(info,"nrtot"          ,               info->nrtot);
    lsgrg_dprint(info,"eps0"           ,                info->eps0);
    lsgrg_dprint(info,"eps1"           ,                info->eps1);
    lsgrg_dprint(info,"eps2"           ,                info->eps2);
    lsgrg_dprint(info,"eps3"           ,                info->eps3);
    lsgrg_dprint(info,"eps4"           ,                info->eps4);
    lsgrg_dprint(info,"eps5"           ,                info->eps5);
    lsgrg_lprint(info,"nnlequ"         ,              info->nnlequ);
    lsgrg_iprint(info,"hrdbnd"         ,              info->hrdbnd);
    lsgrg_dprint(info,"epinit"         ,              info->epinit);
    lsgrg_dprint(info,"eplast"         ,              info->eplast);
    lsgrg_dprint(info,"epdeg"          ,               info->epdeg);
    lsgrg_lprint(info,"ndeg"           ,                info->ndeg);
    lsgrg_dprint(info,"epboun"         ,              info->epboun);
    lsgrg_dprint(info,"epnewt"         ,              info->epnewt);
    lsgrg_dprint(info,"epspiv"         ,              info->epspiv);
    lsgrg_dprint(info,"rtf"            ,                 info->rtf);
    lsgrg_lprint(info,"itlim"          ,               info->itlim);
    lsgrg_lprint(info,"nfix"           ,                info->nfix);
    lsgrg_lprint(info,"nlin"           ,                info->nlin);
    lsgrg_lprint(info,"nnlin"          ,               info->nnlin);
    lsgrg_iprint(info,"maxim"          ,               info->maxim);
    lsgrg_lprint(info,"lgrad"          ,               info->lgrad);
    lsgrg_lprint(info,"lirn"           ,                info->lirn);
    lsgrg_lprint(info,"lmem"           ,                info->lmem);
    lsgrg_lprint(info,"npmp2"          ,               info->npmp2);
    lsgrg_lprint(info,"maxh"           ,                info->maxh);
    lsgrg_lprint(info,"nsear"          ,               info->nsear);
    lsgrg_dprint(info,"epstop"         ,              info->epstop);
    lsgrg_lprint(info,"iper"           ,                info->iper);
    lsgrg_lprint(info,"ipn4"           ,                info->ipn4);
    lsgrg_lprint(info,"ipn5"           ,                info->ipn5);
    lsgrg_lprint(info,"ipn6"           ,                info->ipn6);
    lsgrg_lprint(info,"limser"         ,              info->limser);
    lsgrg_lprint(info,"nstop"          ,               info->nstop);
    lsgrg_lprint(info,"iter"           ,                info->iter);
    lsgrg_lprint(info,"ipr"            ,                 info->ipr);
    lsgrg_lprint(info,"ipr3"           ,                info->ipr3);
    lsgrg_lprint(info,"nb"             ,                  info->nb);
    lsgrg_lprint(info,"ninf"           ,                info->ninf);
    lsgrg_lprint(info,"nobj"           ,                info->nobj);
    lsgrg_lprint(info,"nsuper"         ,              info->nsuper);
    lsgrg_lprint(info,"nzbinv"         ,              info->nzbinv);
    lsgrg_lprint(info,"nzgrad"         ,              info->nzgrad);
    lsgrg_lprint(info,"nzlin"          ,               info->nzlin);
    lsgrg_lprint(info,"nznlin"         ,              info->nznlin);
    lsgrg_lprint(info,"kderiv"         ,              info->kderiv);
    lsgrg_dprint(info,"phmult"         ,              info->phmult);
    lsgrg_dprint(info,"ph1eps"         ,              info->ph1eps);
    lsgrg_lprint(info,"nph1ls"         ,              info->nph1ls);
    lsgrg_lprint(info,"ipvtpt"         ,              info->ipvtpt);
    lsgrg_iprint(info,"fixpiv"         ,              info->fixpiv);
    lsgrg_dprint(info,"xpvpct"         ,              info->xpvpct);
    lsgrg_lprint(info,"ibvblm"         ,              info->ibvblm);
    lsgrg_lprint(info,"iquad"          ,               info->iquad);
    lsgrg_lprint(info,"iscale"         ,              info->iscale);
    lsgrg_lprint(info,"isclag"         ,              info->isclag);
    lsgrg_dprint(info,"edf"            ,                 info->edf);
    lsgrg_iprint(info,"galloc"         ,                info->galloc);
    lsgrg_iprint(info,"getsiz"         ,                info->getsiz);
    lsgrg_lprint(info,"maxb"           ,                info->maxb);
    lsgrg_lprint(info,"maxcg"          ,               info->maxcg);
    lsgrg_lprint(info,"maxr"           ,                info->maxr);
    lsgrg_lprint(info,"maxtbu"         ,              info->maxtbu);
    lsgrg_lprint(info,"nrows"          ,               info->nrows);
    lsgrg_lprint(info,"nvars"          ,               info->nvars);
    lsgrg_iprint(info,"multsb"         ,              info->multsb);
    lsgrg_dprint(info,"funpr"          ,               info->funpr);
    lsgrg_dprint(info,"pstep"          ,               info->pstep);
    lsgrg_dprint(info,"eps"            ,                 info->eps);
    lsgrg_dprint(info,"plinfy"         ,                info->plinfy);
    lsgrg_dprint(info,"plzero"         ,               info->plzero);
    lsgrg_dprint(info,"tolx"           ,                info->tolx);
    lsgrg_dprint(info,"tolz"           ,                info->tolz);
    lsgrg_dprint(info,"cndtol"         ,              info->cndtol);
    lsgrg_dprint(info,"xajtol"         ,              info->xajtol);
    lsgrg_dprint(info,"jacobian_growth_allowance",
                              info->jacobian_growth_allowance);

    lsgrg_msg(IOINFO,"\n"
      " -------------------- End of Global Parameter Dump -----\n");

#endif
        return;


} /* end of function lsgrg_echo_globals*/

void LSGRGCLASS lsgrg_printgrad(LsgrgInfo *info,int print_slacks, char *msg)
{

#ifdef IO_ENABLED

long  i, k, istart, iend, irow, nnz ;
char name[20];
/*---------------------------------------------------------------------*/
/*   lsgrg_printgrad prints the jacobian at the current point          */
/*   code here assumes 1-based arrays                                  */
/*   lsgrg uses column start pointers                                  */
/*---------------------------------------------------------------------*/
     lsgrg_msg(IOINFO,"\n ");
     for(i=1;i<=78;++i) lsgrg_msg(IOINFO,"=");
     sprintf(info->lsgrg_msgbuffer,
     "\n  ........ Jacobian Structure at Current Point ......\n"
     " %s",msg);
     lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

      lsgrg_lprint(info,"nvars",info->nvars); /* print relevant global scalars */
      lsgrg_lprint(info,"nrows",info->nrows);
      lsgrg_lprint(info,"n",info->n);
      lsgrg_lprint(info,"m",info->m);
      lsgrg_lprint(info,"mp1",info->mp1);
      lsgrg_lprint(info,"nzgrad",info->nzgrad);
      lsgrg_lprint(info,"nzlin", info->nzlin );
      lsgrg_lprint(info,"nznlin",info->nznlin);
      lsgrg_msg(IOINFO,"\n");

      lsgrg_msg(IOINFO,"\n   ...... Strucural Columns ......\n");
      for(i=1;i<=info->n+info->mp1;++i) {
          if(i > info->nvars && !print_slacks) break;
          if(i == (info->nvars+1)) lsgrg_msg(IOINFO,
                     "\n   ...... Slack Columns ......\n");

        strcpy(name," ");
        if(info->usernames && i <= info->nvars )
            strncpy(name,info->lsgrg_varnames[i],NAMELENGTH);
        sprintf(info->lsgrg_msgbuffer,"\n"
        " Column %d (%20.14e) %10.10s iheg[i] = %d iheg[i+1]-1 = %d",
        i,(info->x)[i],name,(info->iheg)[i],((info->iheg)[i+1]-1));
        lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

        istart = info->iheg[i];
        iend   = info->iheg[i+1] - 1;
        nnz = iend - istart + 1;
        sprintf(info->lsgrg_msgbuffer,"\n  Column has %d nonzeros ",nnz);
        lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
        if(nnz) {
           for(k=istart;k<=iend;++k) {
              irow = info->ihag[k];
              strcpy(name," ");
              if(info->usernames) strncpy(name,info->lsgrg_rownames[irow],NAMELENGTH);
              sprintf(info->lsgrg_msgbuffer,"\n"
                "      k = %8d row[%7d] %10.10s  coeff = %20.14e",
                 k,irow,name,info->grad[k]);
              lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
           } /* end coefficients loop */
         }
      } /* end column loop */

     lsgrg_msg(IOINFO,"\n  ........ End Jacobian Structure Output ......\n");
     for(i=1;i<=78;++i) lsgrg_msg(IOINFO,"=");

#endif
   return;
} /* end of function lsgrg_printgrad */
void LSGRGCLASS lsgrg_printgrad0(LsgrgInfo *info,
                      int print_slacks, char *msg,long n,
                      long nvars,long m,long nrows,long mp1,
                      long nzgrad, long nzlin, long nznlin,
                      long iheg[], long ihag[] , double grad[] ,
                      double x[])

{

#ifdef IO_ENABLED

long  i, k, istart, iend, irow, nnz ;
char name[20];
/*---------------------------------------------------------------------*/
/*   lsgrg_printgrad0 prints the jacobian at the current point         */
/*   code here assumes 0-based arrays                                  */
/*   lsgrg uses column start pointers                                  */
/*  rownames and varnames are 1-based                                  */
/*---------------------------------------------------------------------*/
     lsgrg_msg(IOINFO,"\n ");
     for(i=1;i<=78;++i) lsgrg_msg(IOINFO,"=");
     sprintf(info->lsgrg_msgbuffer,
     "\n  ........ Jacobian Structure at Current Point ......\n"
     "   ....... using 0-based arrays and iheg/ihegl .....\n"
     " %s",msg);
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

      lsgrg_lprint(info,"nvars",nvars); /* print relevant global scalars */
      lsgrg_lprint(info,"nrows",nrows);
      lsgrg_lprint(info,"n",n);
      lsgrg_lprint(info,"m",m);
      lsgrg_lprint(info,"mp1",mp1);
      lsgrg_lprint(info,"nzgrad",nzgrad);
      lsgrg_lprint(info,"nzlin", nzlin );
      lsgrg_lprint(info,"nznlin",nznlin);
      lsgrg_msg(IOINFO,"\n");

      lsgrg_msg(IOINFO,"\n   ...... Strucural Columns ......\n");
      for(i=1;i<=nvars+mp1;++i) {
          if(i > nvars && !print_slacks) break;
          if(i == (nvars+1)) lsgrg_msg(IOINFO,
                     "\n   ...... Slack Columns ......\n");

        strcpy(name," ");
        if(info->usernames && i <= nvars )
            strncpy(name,info->lsgrg_varnames[i],NAMELENGTH);
        sprintf(info->lsgrg_msgbuffer,"\n"
        " Column %d (%20.14e) %10.10s iheg[i-1] = %d iheg[i]-1 = %d",
        i-1,x[i-1],name,iheg[i-1],iheg[i]-1);
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);

        istart = iheg[i-1];
        iend   = iheg[i]-1;
        nnz = iend - istart + 1;
        sprintf(info->lsgrg_msgbuffer,"\n  Column has %d nonzeros ",nnz);
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
        if(nnz) {
           for(k=istart;k<=iend;++k) {
              irow = ihag[k-1];
              strcpy(name," ");
              if(info->usernames) strncpy(name,info->lsgrg_rownames[irow],NAMELENGTH);
              sprintf(info->lsgrg_msgbuffer,"\n"
                "      k = %8d(%8d) row[%7d] %10.10s  coeff = %20.14e",
                 k,(k-1),irow,name,grad[k-1]);
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
           } /* end coefficients loop */
         }
      } /* end column loop */

     lsgrg_msg(IOINFO,"\n\n  ........ End Jacobian Structure Output ......\n");
     for(i=1;i<=78;++i) lsgrg_msg(IOINFO,"=");

#endif
   return;
} /* end of function lsgrg_printgrad0*/
void LSGRGCLASS lsgrg_output_headers(LsgrgInfo *info,
                      int output,int section, int newpage)
{
#ifdef IO_ENABLED

   if(output == INITIAL_VALUES) {
      if(section == FUNCTIONS ) {
       sprintf(info->lsgrg_msgbuffer,
       "\n\n\n\n                              "
       "OUTPUT OF INITIAL VALUES \n\n "
       "Problem: %s\n\n SECTION 1 -- FUNCTIONS \n",
           info->title );
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
       lsgrg_msg(IOINFO,
            "      Function                    "
           "Initial        Lower         Upper\n "
           "No.   Name     Status  Type       Value         Limit"
           "         Limit\n" );
      }

      if(section == VARIABLES ) {
         sprintf(info->lsgrg_msgbuffer, "\n\n\n\n"
         "Problem: %s"
         "\n\nSECTION 2 -- VARIABLES \n"
         "      Variable                    Initial        Lower"
         "         Upper\n"
         " No.   Name     Status  Type       Value         Limit"
         "         Limit\n",info->title );
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
     }
   } /* end initial values headers */

   if(output == FINAL_RESULTS) {
      if(newpage)  {
         sprintf(info->lsgrg_msgbuffer,
            "\n\f               FINAL RESULTS\n %s\n",info->title );
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      }

      if(section == FUNCTIONS)  {
         sprintf(info->lsgrg_msgbuffer,"\n\nSECTION 1 -- FUNCTIONS   \n\n"
"                                                    Distance\n"
"                  Initial       Final                  From        Lagrange\n"
"   No.   Name      Value        Value       Status    Nearest     Multiplier\n"
"                                                       Bound\n" );
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      }
/*
 xxxxx xxxxxxxx xxxxxxxxxxxx xxxxxxxxxxxx  xxxxxxxx xxxxxxxxxxxx xxxxxxxxxxxx
*/
      if(section == VARIABLES ) {
         sprintf(info->lsgrg_msgbuffer, "\n\nSECTION 2 -- VARIABLES\n\n"
"                                                   Distance\n"
"                Initial      Final                   From      Reduced\n"
" No.   Name      Value       Value     Status       Nearest    Gradient\n"
"                                                    Bound\n" );
           lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      }
/*
 xxxxx xxxxxx xxxxxxxxxxxx xxxxxxxxxxxx xxxxxxxx xxxxxxxxxxxx xxxxxxxxxxxx
*/
    } /* end final results headers */

#endif
  return;
}
void LSGRGCLASS lsgrg_formatname(char *name,int length)
{
/* pad names out to full length to allow %s fields in printf formats */
  int i;
  for(i=::strlen(name)+1;i< length-1;++i) name[i] = ' ';
  name[length-1] = '\0';
  return;
}
void LSGRGCLASS lsgrg_prt_larray(LsgrgInfo *info,
                             char *msg, long x[], long n)
{
#ifdef IO_ENABLED
       long i;
       sprintf(info->lsgrg_msgbuffer,
         "\n\n ... Dump of long array : %s",msg);
       lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
       for(i=0; i<=n; ++i) {
         sprintf(info->lsgrg_msgbuffer,"\n [%5d] = %ld",i,x[i]);
         lsgrg_msg(IOINFO,info->lsgrg_msgbuffer);
      }
       lsgrg_msg(IOINFO,"\n-----------------------------------");
#endif
       return;
}
#ifdef IO_ENABLED
  #undef IOINFO
#endif
#include "lsinfo.h"

#define IOINFO  &_info->io_info
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)

/***********************************************************************/
/*                  File:   grgparsh.c                                 */
/*                  Author: John C. Plummer                            */
/*                                                                     */
/***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
*/
/*  This file contains the finite difference derivative routines       */
/*  for lsgrgc together with the gcomp shells                          */
/*                                                                     */
/*   lsgrg_gcomp1()         calling shell for function evaluations     */
/*                          for lsgrgc      Assumes 1-based g,x        */
/*                          (shell does timing calls for gcomp)        */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*   lsgrg_gcomp0()         calling shell for function evaluations     */
/*                          called from algorithm routines which       */
/*                          assume  0-based g,x.     Calling shell     */
/*                          (grgsub) adjusted pointers passed into     */
/*                           grgitn.  lsgrg_gcomp0()     adjusts       */
/*                           g,x pointers accordingly to work with     */
/*                           1-based gcomp.                            */
/*                          (shell does timing calls for gcomp)        */
/*                                                                     */
/*   lsgrg_parshc1()        computes central difference derivatives    */
/*                            assuming 1-based arrays                  */
/*   lsgrg_parshc0()        computes central difference derivatives    */
/*                          for v1 algorithm caller                    */
/*                          adjusts pointers back to 1-based arrays    */
/*                          0-based arrays                             */
/*                          since v1 assumes 0-based arrays            */
/*                          calls 1-based lsgrg_gcomp()                */
/*                                                                     */
/*   lsgrg_parshf1()        computes forward difference derivatives    */
/*                            assuming 1-based arrays                  */
/*   lsgrg_parshf0()        computes forward difference derivatives    */
/*                          for   caller using 0-based arrays          */
/*                          adjusts pointers back to 1-based arrays    */
/*                          0-based arrays                             */
/*                          calls 1-based lsgrg_gcomp()                */
/*                                                                     */
/*   lsgrg_parsh0()         shell for call to user analytic parsh      */
/*                          for 0-based algorithm.  calls user parsh   */
/*                          and copies elements to grad                */
/*                                                                     */
/*  ** 9/99 jcp ** if  USER_TERMINATION_ENABLED, check return codes    */
/*   from gcomp/parsh and call lsgrg_errorexit if 0                    */
/*                                                                     */
/*---------------------------------------------------------------------*/
void LSGRGCLASS lsgrg_gcomp1(LsgrgInfo *_info,double g[], double x[])
{
/*  gcomp shell for function calls                   */
/*  assumes that callee uses 1-based arrays, so pass */
/*  g,x pointers without adjustment                  */

/*  do timing statistics if timing is enabled        */

#ifdef TIMING_ENABLED
     double t;
     t = lsgrg_timer();
#endif

#ifdef USER_TERMINATION_ENABLED
      if(!_info->lsgrg_user_gcomp(g,x) )   /* call user gcomp function */
          lsgrg_errorexit(JUMPBUF,_LSGRG_USER_TERMINATION);
#else
      //_info->lsgrg_user_gcomp(g,x);  /* call user gcomp function */
	  //
	  // A.G. Dereference the member function
	  //
	  LSGRGSolver* pObject = (LSGRGSolver*)_info->m_pOwner;
	  P_GCOMP pfncomp = _info->lsgrg_user_gcomp;
      (pObject->*pfncomp)(g,x);  /* call user gcomp function */

#endif
      ++(_info->ngcomp);             /* update global call counter */

#ifdef TIMING_ENABLED
        lsgrg_timestats(_info,T_GCOMP,(lsgrg_timer() - t)  );
#endif
}
void LSGRGCLASS lsgrg_gcomp0(LsgrgInfo *_info,double g[], double x[])
{
/*  gcomp shell for calls from lsgrgc algorithm routines which       */
/*  assume  that caller uses 0-based arrays, and that                */
/*  (grgsub) adjusted pointers at call to grgitn, so caller to       */
/*  lsgrg_gcomp_alg1() has g[0] where grgsub has g[1].  decrement    */
/*  g,x pointers since user gcomp assumes 1-based arrays             */

/*  do timing statistics if timing is enabled        */

#ifdef TIMING_ENABLED
     double t;
     t = lsgrg_timer();
#endif

#ifdef USER_TERMINATION_ENABLED
      if(!_info->lsgrg_user_gcomp(--g,--x) )   /* call user gcomp function */
          lsgrg_errorexit(JUMPBUF,_LSGRG_USER_TERMINATION);
#else
      //_info->lsgrg_user_gcomp(--g,--x);  /* call user gcomp function */
	  //
	  // A.G. Dereference the user function
	  //
      LSGRGSolver* pObject = (LSGRGSolver*)_info->m_pOwner;
      P_GCOMP pfncomp = _info->lsgrg_user_gcomp;
      (pObject->*pfncomp)(--g,--x);  /* call user gcomp function */
#endif
      ++(_info->ngcomp);         /* update global call counter */

#ifdef TIMING_ENABLED
        lsgrg_timestats(_info,T_GCOMP,(lsgrg_timer() - t)  );
#endif
}  /* end of lsgrg_gcomp0 */
/*========================================================================*/
void LSGRGCLASS lsgrg_parshf1(LsgrgInfo *_info,double x[], double g[], double gcol[],
                 long col_index,double  gg[],double xub)
{
long  i ;
double dx, tmp, ts;
/*---------------------------------------------------------------------*/
/*  computes finite difference derivatives by forward differences      */
/*  gcomp calls assume that input pointers do 1-based arrays           */
/*---------------------------------------------------------------------*/
#ifdef TIMING_ENABLED
       double t;
       t = lsgrg_timer();
#endif
#ifdef IO_ENABLED
        if(_info->dbg.parshf || _info->ipr > 5) {
           sprintf(LSGRG_MSGBUFFER,
             "\n... entry parshf index = %d nobj = %d pstep = %15.8e",
            col_index,_info->nobj,_info->pstep);
           lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
           for(i=1;i<=_info->n;++i)  {
             sprintf(LSGRG_MSGBUFFER,"\n x[%5d] = %15.8e",i,x[i]);
             lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
           }
         }
#endif

        tmp = g[_info->nobj];
        if( _info->nintbk.ninf != 0 && (!(_info->cbmode.phase0)) )
                       tmp = _info->bestbk.truobj;
        if( _info->maxim )       tmp = -tmp;

        ts = x[col_index];
        dx = (1.e0 + fabs( ts )) *_info->pstep;
        if( xub - ts < dx )    dx = -dx;
        x[col_index] = ts + dx;

#ifdef IO_ENABLED
        if(_info->dbg.parshf) {
             sprintf(LSGRG_MSGBUFFER,"\n dx = %15.8e",
             dx);
             lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        lsgrg_gcomp1(_info, gg, x );

        for( i = 1; i <= _info->mp1; i++ ) {
                gcol[i] = (gg[i] - g[i])/dx;

#ifdef IO_ENABLED
            if(_info->dbg.parshf) {
                sprintf(LSGRG_MSGBUFFER,"\n gg[%5d] = %15.8e"
                " g[%5d] = %15.8e gcol = %15.8e",
                 i,gg[i],i,g[i],gcol[i]);
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
           }
#endif

        }
        gcol[_info->nobj] = (gg[_info->nobj] - tmp)/dx;

        x[col_index] = ts;

#ifdef TIMING_ENABLED
       lsgrg_timestats(_info,T_PARSH,(lsgrg_timer()-t) );
#endif

#ifdef IO_ENABLED
        if(_info->dbg.parshf || _info->ipr > 5) {
           for(i=1;i<=_info->mp1;++i) {
             sprintf(LSGRG_MSGBUFFER,"\n gcol[%5d] = %15.8e",i,gcol[i]);
             lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
           }
           sprintf(LSGRG_MSGBUFFER,"\n... exit parshf...");
           lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
         }
#endif
        return;

} /* end of function lsgrg_parshf1 */

void LSGRGCLASS lsgrg_parshf0(LsgrgInfo *_info,double x[], double g[], double gcol[],
                 long col_index,double  gg[],double xub)
{
/*---------------------------------------------------------------------*/
/*  computes finite difference derivatives by forward differences      */
/*  assumes that caller is using 0-based array references              */
/*  Decrement passed in pointers     so that 1-based references        */
/*  are correct.  then call the 1-based parshf which calls the         */
/*  1-based gcomp                                                      */
/*---------------------------------------------------------------------*/
     lsgrg_parshf1(_info,--x, --g , --gcol, col_index, --gg, xub);

        return;

} /* end of function lsgrg_parshf0 */
/*========================================================================*/

void LSGRGCLASS lsgrg_parshc1(LsgrgInfo *_info,double x[], double gradcol[],
                 long col_index, double gplus[], double gminus[])
{
long i ;
double dx, ts;
/*---------------------------------------------------------------------*/
/*  computes finite difference derivatives by central differences      */
/*   **fixme**  add code to reflect off bounds                         */
/*---------------------------------------------------------------------*/
#ifdef TIMING_ENABLED
       double t;
       t = lsgrg_timer();
#endif

        ts = x[col_index];
        dx = (1.e0 + fabs( ts ))*_info->pstep;
        x[col_index] = ts + dx;

        lsgrg_gcomp1(_info, gplus, x );

        x[col_index] = ts - dx;
        lsgrg_gcomp1(_info, gminus, x );

        for( i = 1; i <= _info->mp1; i++ )
            gradcol[i] = (gplus[i] - gminus[i])/(2.0e0*dx);

        x[col_index] = ts;

#ifdef TIMING_ENABLED
       lsgrg_timestats(_info,T_PARSH,(lsgrg_timer()-t) );
#endif

        return;

} /* end of function lsgrg_parshc1*/

void LSGRGCLASS lsgrg_parshc0(LsgrgInfo *_info,double x[], double gradcol[],
                 long col_index, double gplus[], double gminus[])
{
/*---------------------------------------------------------------------*/
/*  computes finite difference derivatives by central differences      */
/*  assumes that caller uses 0-based arrays                            */
/*     decrement those pointers here so that 1-based references        */
/*  are correct.  then call 1-based parshc which calls 1-based         */
/*  gcomp                                                              */
/*---------------------------------------------------------------------*/

    lsgrg_parshc1(_info,--x,--gradcol,col_index,--gplus,--gminus);

} /* end of function lsgrg_parshc0*/


/*========================================================================*/
void LSGRGCLASS lsgrg_parsh0(LsgrgInfo *_info,double x[],long n,double grad[],
                double paij[],long iprow[], long ipcol[],
                long ipmap[],long nzgrad)
{
long nnz, i;
/*---------------------------------------------------------------------*/
/*  lsgrg_parsh0 is a shell which calls the user analytic parsh        */
/*  function via the supplied pointer.  user routine is assumed to     */
/*  user 1-based arrays, so adjust pointers here                       */
/*---------------------------------------------------------------------*/
#ifdef TIMING_ENABLED
       double t;
       t = lsgrg_timer();
#endif
    --x; --grad; --paij; --iprow; --ipcol; --ipmap;

#ifdef USER_TERMINATION_ENABLED
    if(!_info->lsgrg_user_parsh( x, n, paij, iprow, ipcol, &nnz    ) )
       lsgrg_errorexit(JUMPBUF,_LSGRG_USER_TERMINATION);
#else

    //_info->lsgrg_user_parsh( x, n, paij, iprow, ipcol, &nnz    );
	//
	// A.G.
	//
	LSGRGSolver* pObject = (LSGRGSolver*)_info->m_pOwner;
	P_PARSH pfnparsh = _info->lsgrg_user_parsh;
	(pObject->*pfnparsh)( x, n, paij, iprow, ipcol, &nnz);
#endif

    if(nnz != nzgrad  ) {
#ifdef IO_ENABLED
      sprintf(LSGRG_MSGBUFFER,
        "\n Error: Nbr of Jacobian Nonzeros Returned From User PARSH"
        "[%d] \n Does Not Match Value in Call to GRGSUB [%d]",
        nnz,nzgrad  );
      lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
      lsgrg_errorexit(JUMPBUF, _LSGRG_BAD_USER_NNZ);
    }

      for(i=1;i<=nnz;++i) grad[ipmap[i]] = paij[i];
      ++_info->nparsh;


#ifdef TIMING_ENABLED
       lsgrg_timestats(_info,T_PARSH,(lsgrg_timer()-t) );
#endif

} /* end of function lsgrg_parsh0 */

/*========================================================================*/
void LSGRGCLASS gcomp_check(LsgrgInfo *_info,char *msg,long n, double x[], long mp1,
                 double g[], int onebased)
{
//    long i;
    return;
#ifdef IO_ENABLED
    sprintf(LSGRG_MSGBUFFER,
     "\n\n ... entry gcomp_check: msg = %s",msg);
    lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
   lsgrg_msg(IOINFO,"\n ....... dump of x values .... ");
   if(!onebased) {
      for(i=0; i<n ;++i) {
          sprintf(LSGRG_MSGBUFFER,
           "\n %5ld  %8.8s %12.5e",i,_info->lsgrg_varnames[i+1],x[i]);
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
      }
   }
   else {
      for(i=1; i<=n ;++i) {
          sprintf(LSGRG_MSGBUFFER,
           "\n %5ld  %8.8s %12.5e",i,_info->lsgrg_varnames[i],x[i]);
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
      }
   }
   lsgrg_msg(IOINFO, "\n ....... dump of g values .... ");
   if(!onebased)
      for(i=0; i<mp1 ;++i) {
          sprintf(LSGRG_MSGBUFFER,
           "\n %5ld  %8.8s %12.5e",i,_info->lsgrg_rownames[i+1],g[i]);
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
      }
   else
      for(i=1; i<=mp1 ;++i) {
          sprintf(LSGRG_MSGBUFFER,
           "\n %5ld  %8.8s %12.5e",i,_info->lsgrg_rownames[i],g[i]);
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
      }
   lsgrg_msg(IOINFO,"\n");
   for(i=1;i<=78;++i) lsgrg_msg(IOINFO,"-");
   return;
#endif
}
void LSGRGCLASS lsgrg_gcompx(LsgrgInfo *_info,double g[], double x[],
                             long *nbc, long ibc[])
{
/*  shell for gcomp with individual function access */

   if(_info->lsgrg_user_gcomp_fcns != NULL )
   {
      // _info->lsgrg_user_gcomp_fcns(g,x,nbc,ibc);
	  //
	  // A.G. Dereference the user functions
	  //
      LSGRGSolver* pObject = (LSGRGSolver*)_info->m_pOwner;
      P_GCOMPX pfncompx = _info->lsgrg_user_gcomp_fcns;
      (pObject->*pfncompx)(g,x,nbc,ibc);
   }
   else  {
#ifdef IO_ENABLED
       sprintf(LSGRG_MSGBUFFER,
          "\n Error: Attempt to Call GCOMPX (Individual Function Access)"
          "\n        But No Function Supplied (Pointer is NULL) " );
       lsgrg_error_msg(IOINFO,LSGRG_MSGBUFFER);
#endif
       lsgrg_errorexit(JUMPBUF,_LSGRG_MISSING_GCOMPX);
     }
}
#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"
/*
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
*/
/*--------------------------------------------------------------------*/
/*   FILE:  grgopt.c   Author John C. Plummer                         */
/*                                                                    */
/*  This file contains the functions which obtain option settings     */
/*   from the user, verify, and post those values to the global       */
/*   parameters in the LsgrgInfo* sructure passed as an arg           */
/*                                                                    */
/*                                                                    */
/*  lsgrg_setparameter()  called by user prior to grgsub to set       */
/*                        optional parameters.  1st call sets up      */
/*                        user option name. checks validity and       */
/*                        posts user value to global parameter        */
/*                        overriding default or previous value        */
/*                        returns 0 on invalid value or unknown       */
/*                        name. 1 otherwise                           */
/*                                                                    */
/*  lsgrg_post_option() posts user specified parameter values to      */
/*                      global parameter variables                    */
/*                                                                    */
/*  lsgrg_check_options() checks global parameters for                */
/*                       invalid values to catch bad default values   */
/*                                                                    */
/*  lsgrg_setcomputedparms() sets values for parameters computed from */
/*                        other parameters after user values          */
/*                        have been posted                            */
/*   lsgrg_upcase()       converts a string to upper case             */
/*                                                                    */
/*   lsgrg_lcase()        converts a string to lower case             */
/*                                                                    */
/*   functions to set/recall frequently used values and options       */
/*   ----------------------------------------------------------       */
/*                                                                    */
/*   lsgrg_get_plinfy()   returns default value for +infinity         */
/*                        to user to aid in bounds setting            */
/*                                                                    */
/*   lsgrg_set_title()    copys user supplied string to title         */
/*                        global                                      */
/*   lsgrg_set_fderivatives()    set forward difference derivatives   */
/*   lsgrg_set_cderivatives()    set central difference derivatives   */
/*   lsgrg_set_userderivatives() set user supplied derivatives        */
/*   lsgrg_set_printlevel()      set lsgrg print level                */
/*   lsgrg_set_inprnt()          toggle pre-solution output           */
/*   lsgrg_set_otprnt()          toggle post-solution output          */
/*   lsgrg_set_scalingOn()       turns on problem scaling             */
/*   lsgrg_set_scalingOff()      turns off problem scaling            */
/*   lsgrg_set_error_output()    toggle error output                  */
/*   lsgrg_set_screen_output()   toggle screen output                 */
/*                                                                    */
/* ** 4/27/02 jcp ** add new parameter for maxtime and function       */
/*   lsgrg_set_maxtime(double t) to set it                            */
/*                                                                    */
/*--------------------------------------------------------------------*/
int LSGRGCLASS lsgrg_setparameter(LsgrgInfo *info,char *name,
               long lvalue, double dvalue)
{
     char name_in[20];
     int i;
/*-----------------------------------------------------------------*/
/*  lsgrg_setparameter stores a user specified value for an lsgrg  */
/*  parameter in the  lsgrg_option_name   structure . on the       */
/*  first call, the option names are filled in according to the    */
/*  coding defined in the enum lsgrg_user_options                  */
/*                                                                 */
/*-----------------------------------------------------------------*/
   if(info->nbr_user_options_set == 0) {
      strcpy(info->lsgrg_option_name[KDERIV].name,"kderiv");
      strcpy(info->lsgrg_option_name[MAXBAS].name,"maxbas");
      strcpy(info->lsgrg_option_name[MAXHES].name,"maxhes");
      strcpy(info->lsgrg_option_name[EPNEWT].name,"epnewt");
      strcpy(info->lsgrg_option_name[EPINIT].name,"epinit");
      strcpy(info->lsgrg_option_name[EPSTOP].name,"epstop");
      strcpy(info->lsgrg_option_name[EPSPIV].name,"epspiv");
      strcpy(info->lsgrg_option_name[PH1EPS].name,"ph1eps");
      strcpy(info->lsgrg_option_name[NSTOP ].name,"nstop");
      strcpy(info->lsgrg_option_name[ITLIM ].name,"itlim" );
/* ** 4/27/02 jcp ** new option 'maxtime'   */
      strcpy(info->lsgrg_option_name[MAXTIME].name,"maxtime" );

      strcpy(info->lsgrg_option_name[LIMSER].name,"limser");
      strcpy(info->lsgrg_option_name[IPR   ].name,"ipr"   );
      strcpy(info->lsgrg_option_name[IPN4  ].name,"ipn4"  );
      strcpy(info->lsgrg_option_name[IPN5  ].name,"ipn5"  );
      strcpy(info->lsgrg_option_name[IPN6  ].name,"ipn6"  );
      strcpy(info->lsgrg_option_name[IPER  ].name,"iper"  );
      strcpy(info->lsgrg_option_name[PSTEP ].name,"pstep" );
      strcpy(info->lsgrg_option_name[IQUAD ].name,"iquad" );
      strcpy(info->lsgrg_option_name[MODGG ].name,"modcg" );
      strcpy(info->lsgrg_option_name[AIJTOL].name,"aijtol");
      strcpy(info->lsgrg_option_name[PIVPCT].name,"pivpct");
      strcpy(info->lsgrg_option_name[MXTABU].name,"mxtabu");
      strcpy(info->lsgrg_option_name[FUNPR ].name,"funpr" );
      strcpy(info->lsgrg_option_name[CONDTL].name,"condtl");
      strcpy(info->lsgrg_option_name[IDEGLM].name,"ideglm");
      strcpy(info->lsgrg_option_name[EPBOUN].name,"epboun");
      strcpy(info->lsgrg_option_name[EPDEG ].name,"epdeg" );
      strcpy(info->lsgrg_option_name[ISCAL ].name,"iscale" );
      strcpy(info->lsgrg_option_name[ISCLG ].name,"isclag" );
      strcpy(info->lsgrg_option_name[MEMCG ].name,"memcg" );
      strcpy(info->lsgrg_option_name[IBVBLM].name,"ibvblm");
      strcpy(info->lsgrg_option_name[HRDBND].name,"hrdbnd");
      strcpy(info->lsgrg_option_name[FIXPIV].name,"fixpiv");
      strcpy(info->lsgrg_option_name[INPRNT].name,"inprnt");
      strcpy(info->lsgrg_option_name[OTPRNT].name,"otprnt");
      strcpy(info->lsgrg_option_name[GFEAS ].name,"gfeas" );
      strcpy(info->lsgrg_option_name[USEPH0].name,"useph0");
      strcpy(info->lsgrg_option_name[STEEPD].name,"steepest");
      strcpy(info->lsgrg_option_name[STEPLIM].name,"steplimit");
   }
     strncpy(name_in,name,20); /* copy name string to local var */
     lsgrg_lcase(name_in);     /* and set case to lower         */
/*------------------------------------------------------------------*/
/*  search table for option and store value when found. check for   */
/*  option validity and return failure if invalid                   */
/*  3rd arg to post_option is flag to post to global parm           */
/*------------------------------------------------------------------*/
     for(i=LSGRG_OPTION_START+1;i < LSGRG_OPTION_END  ; ++i)
        if(strcmp(info->lsgrg_option_name[i].name,name_in)==0) {
           if( !lsgrg_post_option(info,i,lvalue,dvalue,TRUE) ) return 0;
           else {
                 ++(info->nbr_user_options_set);
                 return 1;
           }
         }
/*------------------------------------------------------------------*/
/*  unknown option name -- return failure and msg if io enabled     */
/*------------------------------------------------------------------*/
#ifdef IO_ENABLED
       sprintf(info->lsgrg_msgbuffer,
        "\n Error(lsgrg_setparameter): Unknown Parameter Name [%s]",
        name_in);
        lsgrg_error_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
     return 0;
}

int LSGRGCLASS lsgrg_post_option(LsgrgInfo *info,
               int option_index,long lvalue, double dvalue,
                      int post)
{
int error;
/*------------------------------------------------------------------*/
/*  check an individual option value for validity                   */
/*   if 'post' is true, post value to corresponding global value    */
/*   this allows use for checking user settings prior to grgsub call*/
/*   and after the fact checking of all parameter values            */
/* **fixme** double check identities of global parms                */
/*------------------------------------------------------------------*/
    error = 0;
    switch(option_index) {

/* ** 4/27/02 jcp ** add maxtime */
       case MAXTIME:
              if( dvalue  < 0.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "A Negative Value has Been Specified for "
                " maxtime = %10.2f",
                dvalue );
#endif
              }
              else if(post) info->maxtime = dvalue;
              break;

       case KDERIV:
              if( lvalue  > 3 || lvalue < 0  ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "kderiv = %d Valid values are 0,1,2,3",
                lvalue );
#endif
              }
              else if(post) info->kderiv = lvalue;
              break;

       case MAXBAS:
              if( lvalue  < 1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "maxbas = %d Valid values are >= 1",
                lvalue );
#endif
              }
              else if(post) info->maxb   = lvalue;
              break;
       case MAXHES:
              if( lvalue  < 0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "maxhes = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->maxr   = lvalue;
              break;
       case EPNEWT:
              if( dvalue  <= 0.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "epnewt = %15.8e  Valid values are > 0.0",
                dvalue );
#endif
              }
              else if(post) info->epnewt = dvalue;
              break;

       case EPINIT:
              if( dvalue  <= 0.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "epinit = %15.8e  Valid values are > 0.0",
                dvalue );
#endif
              }
              else if(post) info->epinit = dvalue;
              break;

       case EPSTOP:
              if( dvalue  <= 0.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "epstop = %15.8e  Valid values are > 0.0",
                dvalue );
#endif
              }
              else if(post) info->epstop = dvalue;
              break;

       case EPSPIV:
              if( dvalue  <= 0.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "epspiv = %15.8e  Valid values are > 0.0",
                dvalue );
#endif
              }
              else if(post) info->epspiv = dvalue;
              break;

       case PH1EPS:
              if( dvalue  <= 0.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "ph1eps = %15.8e  Valid values are > 0.0",
                dvalue );
#endif
              }
              else if(post) info->ph1eps = dvalue;
              break;

       case NSTOP:
              if( lvalue  <  1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "nstop = %d Valid values are >= 1",
                lvalue );
#endif
              }
              else if(post) info->nstop  = lvalue;
              break;

       case ITLIM :
              if( lvalue  <  1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "itlim  = %d Valid values are >= 1",
                lvalue );
#endif
              }
              else if(post) info->itlim  = lvalue;
              break;

       case LIMSER:
              if( lvalue  <  0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "limser = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->limser = lvalue;
              break;

       case IPR   :
              if( lvalue  <  -1 || lvalue > 6 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "ipr    = %d Valid values are -1 through 6",
                lvalue );
#endif
              }
              else if(post) info->ipr    = lvalue;
              break;

       case IPN4  :
              if( lvalue  <  0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "ipn4   = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->ipn4   = lvalue;
              break;

       case IPN5  :
              if( lvalue  <  0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "ipn5   = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->ipn5   = lvalue;
              break;

       case IPN6  :
              if( lvalue  <  0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "ipn6   = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->ipn6   = lvalue;
              break;

       case IPER  :
              if( lvalue  <  0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "iper   = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->iper   = lvalue;
              break;

       case PSTEP :
              if( dvalue  <= 0.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "pstep  = %15.8e  Valid values are > 0.0",
                dvalue );
#endif
              }
              else if(post) info->pstep  = dvalue;
              break;

       case IQUAD :
              if( lvalue  != 0 && lvalue != 1) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "iquad  = %d Valid values are 0 and 1",
                lvalue );
#endif
              }
              else if(post) info->iquad  = lvalue;
              break;

       case MODGG :
              if( lvalue  < 1 || lvalue > 6  ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "modcg  = %d Valid values are 1 through 6",
                lvalue );
#endif
              }
              else if(post) info->modcg  = lvalue;
              break;

       case AIJTOL:
              if( dvalue  <= info->eps ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "aijtol = %15.8e  Valid values >= eps(%15.8e)",
                dvalue );
#endif
              }
              else if(post) info->xajtol = dvalue;
              break;

       case PIVPCT:
              if( dvalue  < 0.0 || dvalue > 1.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "pivpct = %15.8e  Valid values 0.0 through 1.0",
                dvalue );
#endif
              }
              else if(post) info->xpvpct = dvalue;
              break;

       case MXTABU:
              if( lvalue  < 0  ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "mxtabu = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->maxtbu = lvalue;
              break;

       case FUNPR :
              if( dvalue  < 0.0 || dvalue > 1.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "funpr  = %15.8e  Valid values 0.0 through 1.0",
                dvalue );
#endif
              }
              else if(post) info->funpr  = dvalue;
              break;

       case CONDTL:
              if( lvalue  < 0  ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "condtl = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->cndtol = lvalue;
              break;

       case IDEGLM:
              if( lvalue  <= 0  ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "ideglm = %d Valid values are > 0",
                lvalue );
#endif
              }
              else if(post) info->idglim = lvalue;
              break;

       case EPBOUN:
              if( dvalue  <= 0.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "epboun = %15.8e  Valid values are > 0.0",
                dvalue );
#endif
              }
              else if(post) info->epboun = dvalue;
              break;

       case EPDEG :
              if( dvalue  <= 0.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "epdeg  = %15.8e  Valid values are > 0.0",
                dvalue );
#endif
              }
              else if(post) info->epdeg  = dvalue;
              break;

       case ISCAL :
              if( lvalue  != 0 && lvalue != 1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "iscal  = %d Valid values are 0 and 1",
                lvalue );
#endif
              }
              else if(post) info->iscale = lvalue;
              break;

       case ISCLG :
              if( lvalue  < 0  ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "isclag  = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->isclag = lvalue;
              break;

       case MEMCG :
              if( lvalue  < 0  ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "memcg  = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->memcg  = lvalue;
              break;

       case IBVBLM:
              if( lvalue  < 0  ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "ibvblm = %d Valid values are >= 0",
                lvalue );
#endif
              }
              else if(post) info->ibvblm = lvalue;
              break;

       case HRDBND:
              if( lvalue  !=0 && lvalue != 1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "hrdbnd = %d Valid values are 0 and 1",
                lvalue );
#endif
              }
              else if(post) info->hrdbnd = lvalue;
              break;

       case FIXPIV:
              if( lvalue  !=0 && lvalue != 1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "fixpiv = %d Valid values are 0 and 1",
                lvalue );
#endif
              }
              else if(post) info->fixpiv = lvalue;
              break;

       case INPRNT:
              if( lvalue  !=0 && lvalue != 1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "inprnt = %d Valid values are 0 and 1",
                lvalue );
#endif
              }
              else if(post) info->inprnt = lvalue;
              break;

       case OTPRNT:
              if( lvalue  !=0 && lvalue != 1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "otprnt = %d Valid values are 0 and 1",
                lvalue );
#endif
              }
              else if(post) info->otprnt = lvalue;
              break;

       case GFEAS:
              if( lvalue  !=0 && lvalue != 1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "gfeas = %d Valid values are 0 and 1",
                lvalue );
#endif
              }
              else if(post) info->gfeas  = lvalue;
              break;

       case USEPH0:
              if( lvalue  !=0 && lvalue != 1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "useph0 = %d Valid values are 0 and 1",
                lvalue );
#endif
              }
              else if(post) info->useph0 = lvalue;
              break;

       case STEEPD:
              if( lvalue  !=0 && lvalue != 1 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "steepest = %d Valid values are 0 and 1",
                lvalue );
#endif
              }
              else if(post) info->steepest_descent = lvalue;
              break;

       case STEPLIM:
              if( dvalue  <=0.0 || dvalue > 1.0 ) {
                  ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "steplim = %d Valid values > 0 and <= 1",
                dvalue );
#endif
              }
              else if(post) info->steplimit = lvalue;
              break;
/*------------------------------------------------------------------*/
/*  default case  -- unknown option should never occur              */
/*------------------------------------------------------------------*/
       default:
                 ++error;
#ifdef IO_ENABLED
                sprintf( info->lsgrg_msgbuffer,"\n Error: "
                "Internal Error: Unknown lsgrg option index = %d",
                lvalue );
#endif

        };  /* end option case statement */

            if(error) {
               info->lsgrg_return_status = _LSGRG_BAD_OPTION_VALUE;
#ifdef IO_ENABLED
               lsgrg_error_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
                return 0;
            }
            else return 1;
}

int LSGRGCLASS lsgrg_check_options(LsgrgInfo *info)
{
/*---------------------------------------------------------------------*/
/*  lsgrg_check_option() calls lsgrg_post_option() once for each       */
/*  parameter to be checked for invalid values.  this allows checking  */
/*  for invalid default values. called with 'post' = false             */
/*---------------------------------------------------------------------*/
    int error=0;
/* ** 4/27/02 jcp ** add entry for maxtime */
   error += lsgrg_post_option(info,MAXTIME ,(long) 0,info->maxtime,FALSE);

   error += lsgrg_post_option(info,EPNEWT ,(long) 0,info->epnewt,FALSE);
   error += lsgrg_post_option(info,EPINIT ,(long) 0,info->epinit,FALSE);
   error += lsgrg_post_option(info,EPSTOP ,(long) 0,info->epstop,FALSE);
   error += lsgrg_post_option(info,EPSPIV ,(long) 0,info->epspiv,FALSE);
   error += lsgrg_post_option(info,PH1EPS ,(long) 0,info->ph1eps,FALSE);
   error += lsgrg_post_option(info,PSTEP  ,(long) 0,info->pstep ,FALSE);
   error += lsgrg_post_option(info,EPBOUN ,(long) 0,info->epboun,FALSE);
   error += lsgrg_post_option(info,EPDEG  ,(long) 0,info->epdeg ,FALSE);
   error += lsgrg_post_option(info,AIJTOL ,(long) 0,info->xajtol,FALSE);
   error += lsgrg_post_option(info,PIVPCT ,(long) 0,info->xpvpct,FALSE);
   error += lsgrg_post_option(info,FUNPR  ,(long) 0,info->funpr ,FALSE);
   error += lsgrg_post_option(info,CONDTL ,(long) 0,info->cndtol,FALSE);

   error += lsgrg_post_option(info,KDERIV ,info->kderiv, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,MAXBAS ,info->maxb  , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,MAXHES ,info->maxr  , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,NSTOP  ,info->nstop , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,ITLIM  ,info->itlim , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,LIMSER ,info->limser, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,IPR    ,info->ipr   , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,IPN4   ,info->ipn4  , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,IPN5   ,info->ipn5  , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,IPN6   ,info->ipn6  , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,IPER   ,info->iper  , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,IQUAD  ,info->iquad , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,MODGG  ,info->modcg , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,MXTABU ,info->maxtbu, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,IDEGLM ,info->idglim, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,ISCAL  ,info->iscale, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,ISCLG  ,info->isclag, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,MEMCG  ,info->memcg , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,IBVBLM ,info->ibvblm, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,HRDBND ,info->hrdbnd, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,FIXPIV ,info->fixpiv, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,INPRNT ,info->inprnt, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,OTPRNT ,info->otprnt, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,GFEAS  ,info->gfeas , (double) 0.0,FALSE);
   error += lsgrg_post_option(info,USEPH0 ,info->useph0, (double) 0.0,FALSE);
   error += lsgrg_post_option(info,STEEPD ,info->steepest_descent,
                                              (double) 0.0,FALSE);
   error += lsgrg_post_option(info,STEPLIM ,(long) 0,
                                   info->steplimit,FALSE);
   if(error) return 0;
   else return 1;
}
/*----------------------------------------------------------------------*/
/*  routines to allocate and store user-specified var/row names         */
/*----------------------------------------------------------------------*/
int LSGRGCLASS lsgrg_setvarname(LsgrgInfo *info,long nnvars,
                                long index,char *name)
{
/*   long i; */
   long nvars;
  /*  on first call, attempt to allocate memory for var names */
/* ** 5/6/02 jcp ** check index against nnvars value used for */
/*  allocation */

   if(info->dbg.setvarname)
      lsgrg_msg(&(info->io_info),"\n.... Entry lsgrg_setvarname");

   if(info->lsgrg_varnames==NULL)
      if( !lsgrg_alloc_varnames(info,nnvars)) return 0;

  /* place this var name -- check for bad index */

      nvars = info->nvars_varnames;
      if(nvars != nnvars) {
#ifdef IO_ENABLED
        sprintf(info->lsgrg_msgbuffer,
          "\n Error(lsgrg_setvarname): \n"
          "Value of nnvars (%d) has changed since first"
          " call (%d).  Name Rejected",nnvars,nvars);
        lsgrg_error_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
         return 0;
     }

     if(index < 1 || index > nvars) {
#ifdef IO_ENABLED
        sprintf(info->lsgrg_msgbuffer,
          "\n Error(lsgrg_setvarname): \n"
          "Index (%d) of Variable Name [%s] is < 1"
          " Or > nnvars (%d)",index,name,nnvars);
        lsgrg_error_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
        return 0;
     }
     strncpy(info->lsgrg_varnames[index],name,NAMELENGTH);
     if(info->dbg.setvarname)
        lsgrg_msg(&(info->io_info),"\n.... Exit lsgrg_setvarname");
     return 1;
}
int LSGRGCLASS lsgrg_setrowname(LsgrgInfo *info,
                                long nnrows,long index,char *name)
{
/*   long i; */
     long nrows;
  /*  on first call, attempt to allocate memory for row names */
/* ** 5/6/02 jcp ** check index against nnrows value used for */
/*  allocation */

   if(info->dbg.setrowname)
       lsgrg_msg(&(info->io_info),"\n.... Entry lsgrg_setrowname");

     if(info->lsgrg_rownames==NULL)
        if( !lsgrg_alloc_rownames(info, nnrows)) return 0;


      nrows = info->nrows_rownames;
      if(nrows != nnrows) {
#ifdef IO_ENABLED
        sprintf(info->lsgrg_msgbuffer,
          "\n Error(lsgrg_setrowname): \n"
          "Value of nnrows (%d) has changed since first"
          " call (%d).  Name Rejected",nnrows,nrows);
        lsgrg_error_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
         return 0;
     }

  /* place this row name -- check for bad index */
     if(index < 1 || index > nnrows) {
#ifdef IO_ENABLED
        sprintf(info->lsgrg_msgbuffer,
           "\n Error(lsgrg_setrowname):\n"
           " Index (%d) of Variable Name [%s] is < %1"
           " Or > nnrows (%d)",index,name,nnrows);
        lsgrg_error_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
        return 0;
     }
     strncpy(info->lsgrg_rownames[index],name,NAMELENGTH);
    if(info->dbg.setrowname)
        lsgrg_msg(&(info->io_info),"\n.... Exit lsgrg_setrowname");
     return 1;
}

void LSGRGCLASS lsgrg_upcase(char *string)
{
/*---------------------------------------------------------------------*/
/*   this routine converts 'string' to upper case.  commented out code */
/*   is a use a brute force approach provided for the tiny possiblity  */
/*   that this code might run be run on an EBCDIC machine              */
/*---------------------------------------------------------------------*/
    int i;
     for(i=0;string[i];++i)
         if(string[i] >= 'a' && string[i] <= 'z') string[i] -= (char)32;
/*
     int i,j;
     char letters[53]=
          {"aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ"};

     for(i=0;string[i];++i) {
         for(j=0; j<=52; j+=2) {
            if( string[i] == letters[j] ) {
                string[i] = letters[j+1];
                break;
             }
         }
     }
*/
     return;
}
void LSGRGCLASS lsgrg_lcase(char *string)
{
/*   this routine converts 'string' to lower case.                      */

    int i;
     for(i=0;string[i];++i)
         if(string[i] >= 'A' && string[i] <= 'Z') string[i] += (char)32;
/*
     int i,j;
     char letters[53]=
          {"aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ"};

     for(i=0;string[i];++i) {
         for(j=1; j<=51; j+=2) {
            if( string[i] == letters[j] ) {
                string[i] = letters[j-1];
                break;
             }
         }
     }
*/
     return;
}
double LSGRGCLASS lsgrg_get_plinfy(LsgrgInfo *info)
 { return info->plinfy;}

void LSGRGCLASS lsgrg_set_title(LsgrgInfo*info,char *usertitle)
{
    if(::strlen(usertitle) <= TITLELENGTH)
       strcpy(info->title,usertitle);
    else
		strncpy(info->title,usertitle,(::strlen(info->title)) );
    return;
}
void LSGRGCLASS lsgrg_set_fderivatives(LsgrgInfo *opt)
{    lsgrg_setparameter(opt, "kderiv",(long)0, (double)0.0); }

void LSGRGCLASS lsgrg_set_cderivatives(LsgrgInfo *opt)
{    lsgrg_setparameter(opt,"kderiv",(long)1, (double)0.0); }

void LSGRGCLASS lsgrg_set_userderivatives(LsgrgInfo *opt)
{    lsgrg_setparameter(opt,"kderiv",(long)2, (double)0.0); }

int LSGRGCLASS lsgrg_set_printlevel(LsgrgInfo *opt,long ipr)
{   return (lsgrg_setparameter(opt,"ipr",(long)ipr,(double)0.0)); }

int LSGRGCLASS lsgrg_set_inprnt(LsgrgInfo *opt,long OnOff)
{   return (lsgrg_setparameter(opt,"inprnt",(long)OnOff,(double)0.0)); }

int LSGRGCLASS lsgrg_set_otprnt(LsgrgInfo *opt,long OnOff)
{   return (lsgrg_setparameter(opt,"otprnt",(long)OnOff,(double)0.0)); }

void LSGRGCLASS lsgrg_set_scalingOn(LsgrgInfo *opt)
{   lsgrg_setparameter(opt,"iscale",(long)1,(double)0.0); }

void LSGRGCLASS lsgrg_set_scalingOff(LsgrgInfo *opt)
{   lsgrg_setparameter(opt,"iscale",(long)0,(double)0.0); }

void LSGRGCLASS lsgrg_set_error_output(LsgrgInfo *opt,int OnOff)
{    opt->io_info.error_output_enabled = OnOff; }

void LSGRGCLASS lsgrg_set_screen_output(LsgrgInfo *opt,int OnOff)
{    opt->io_info.screen_output_enabled = OnOff; }

void LSGRGCLASS lsgrg_steepest_descent(LsgrgInfo *opt,int OnOff)
{    lsgrg_setparameter(opt,"steepest", (long) OnOff, (double) 0.0); }
void LSGRGCLASS lsgrg_limit_stepsize(LsgrgInfo *opt,double steplim)
{    lsgrg_setparameter(opt,"steplimit",(long) 0, steplim); }
/* ** 4/27/02 jcp ** add function to set maxtime */
void LSGRGCLASS lsgrg_set_maxtime(LsgrgInfo *opt, double t) {
     lsgrg_setparameter(opt,"maxtime",(long) 0, t);
}
#include "lsinfo.h"
/*
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
*/
/*---------------------------------------------------------------------*/
/*  FILE: grginit.c  Author: John C. Plummer                           */
/*                                                                     */
/*   This file contains functions for setting default values for       */
/*   lsgrg global parameters                                           */
/*                                                                     */
/*   LsgrgInitialize()                                                 */
/*     .  Called by user to initialize lsgrgc system.  allocates       */
/*        memory for an LsgrgInfo structure to hold all shell          */
/*        global values, sets default initial values, and returns      */
/*        pointer to the structure. Returns NULL if allocation fails   */
/*                                                                     */
/*   LsgrgExit()                                                       */
/*     .  Called by user after run completion to deallocate resources  */
/*        in LsgrgInfo structure                                       */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*                                                      Last Update    */
/*                                                      -----------    */
/*   lsgrg_set_defaults()                                may 1998      */
/*          sets default values for lsgrg user-modifiable parameters   */
/*                                                                     */
/*                                                                     */
/*   lsgrg_initialize_1time()                              mar 1998    */
/*             .   performs one-time only initializations              */
/*                 computes machine dependent parameters               */
/*                 assigns default values for some tolerances          */
/*                                                                     */
/*             .   called from lsgrg_setdefaults to insure             */
/*                 correct calling sequence                            */
/*  ** 4/27/02 jcp ** set default for run time limit maxtime           */
/*---------------------------------------------------------------------*/
LsgrgInfo*  LSGRGCLASS LsgrgInitialize(void)
{
/*  LsgrgInitialize creates a structure of type LsgrgInfo, initializes */
/*  all values and returns pointer to caller                           */
/*  zero out all previously global values before implanting defaulte   */

    LsgrgInfo * lsgrginfo;
    lsgrginfo = (LsgrgInfo*) malloc( (size_t) (sizeof(LsgrgInfo) ) );
    if(lsgrginfo != NULL){
        lsgrg_zero_shell_globals(lsgrginfo);
        lsgrg_zero_alg_globals(lsgrginfo);
        lsgrg_set_defaults(lsgrginfo);
        lsgrginfo->nbr_user_options_set = 0;
    }

    return lsgrginfo;
}
void LSGRGCLASS LsgrgBye(LsgrgInfo *_info)
{

/*----------------------------------------------------------------*/
/*   free row and column name memory if allocated                 */
/*----------------------------------------------------------------*/
     if(_info->lsgrg_varnames != NULL) lsgrg_free_varnames(_info);
     if(_info->lsgrg_rownames != NULL) lsgrg_free_rownames(_info);
/*----------------------------------------------------------------*/
/*     free memory allocated for structure                        */
/*----------------------------------------------------------------*/
        free(_info);

        return;
}

void LSGRGCLASS lsgrg_set_defaults(LsgrgInfo *_info)
{
  int i;
/*------------------------------------------------------------------*/
/*  set default values for all scalar parameters                    */
/*  call lsgrg_initialize_1time to compute eps2, tolz, tolx         */
/*  which are used below                                            */
/*------------------------------------------------------------------*/
   lsgrg_initialize_1time(_info);
   _info->plinfy = PLINFY;
   _info->plzero = PLZERO;
/* ** 4/27/02 jcp ** NULL pointer table */
    for(i=0; i < MAXALLOC; ++i) _info->ptable[i] = NULL;

/* ** 4/27/02 jcp ** add maxtime */
    _info->maxtime = LSGRG_TIME_LIMIT;

/* ** 4/01 jcp ** add inits for reallocation factors */
    _info->binv_realloc_factor = 1.5;
    _info->jacobian_growth_realloc_factor = 2.0;

/*  4/01 jcp ** zero retry count
     _info->retry_history[0] = 0;
/*  7/01 jcp ** default retry strategy to enabled */
     _info->disable_retries = 0;
/*  7/01 jcp ** default mult return to packed */
     _info->return_mults_unpacked = 0;
/*  7/01 jcp ** default to copyback user structure*/
     _info->disable_copyback = 0;

/*  set io_info substructure values */
/* 6/01 jcp ** set logicals to indicate user set io units */

   _info->io_info.error_output_enabled  = 1;  /* default is on */
   _info->io_info.screen_output_enabled = 0;  /* default is off */
/* 3/02  add logical for ioout */
   _info->io_info.ioout_enabled = 1;  /* default is on */
   _info->io_info.ioout_set_by_user = 0; /* initialize to no user settings*/
   _info->io_info.ioerr_set_by_user = 0;
   _info->io_info.ioterm_set_by_user = 0;
   _info->io_info.ioerr_to_ioterm  = 0;
   _info->io_info.ioerr_to_ioout   = 0;
#ifdef FILE_IO_ENABLED
   _info->io_info.lsgrg_ioout  = NULL;
   _info->io_info.lsgrg_ioerr  = NULL;
   _info->io_info.lsgrg_ioterm = NULL;
#endif

   _info->inprnt = 0;                    _info->otprnt = 0;
/* 3/20/02 jcp ** add initialization of outres print level */
   _info->outres_print_level = OUTRES_PRINT_ALL;

   _info->maxr = -1;  /* place -1 value so that we can put in */
                     /*  true default in grgsub if user does  */
                     /*  specify a value via setparameter()   */

   _info->epnewt = 1.0e-4;               _info->epstop = 1.0e-4;
   _info->epspiv = 1.0e-6;               _info->epinit = 0.0;
   _info->epdeg  = 1.0e-4;               _info->ph1eps = 0.;
   _info->epboun = _info->tolx;           _info->pstep = _info->eps2;
   _info->xajtol = _info->tolz;           _info->xpvpct = 1.0e-01;
   _info->funpr  = 1.0e-08;              _info->cndtol = 1.0e08;

   _info->multsb = TRUE;                 _info->hscale = FALSE;
   _info->fixpiv = FALSE;                _info->hrdbnd = FALSE;
   _info->modcg  = 6;                    _info->nstop  = 3;
   _info->itlim  = 10;                   _info->idglim = 25;
   _info->iscale = 0;                    _info->isclag = 0;
   _info->limser = LSGRG_ITERATION_LIMIT;  _info->ipr    = 1;
   _info->ipn5   = 0;                    _info->ipn4   = 0;
   _info->ipn6   = 0;                    _info->iper   = 0;
   _info->iquad  = 0;                    _info->kderiv = 0;
   _info->ipvtpt = 3;                    _info->ibvblm = 2;

   _info->maxcg  = 3;                    _info->maxtbu = 25;
   _info->useph0 = TRUE;
   _info->gfeas  = FALSE;

   _info->alloctable_entries = -1; /* index for pointer allocation */
                                  /* table                        */
/*--------------------------------------------------------------------*/
/*  for user analytic parsh, initialize structure pointers to null    */
/*  since they are not allocated unless user specified analytic       */
/*--------------------------------------------------------------------*/
    _info->iprow = NULL;     _info->ipcol = NULL;
    _info->paij  = NULL;     _info->ipmap = NULL;
/*--------------------------------------------------------------------*/
/*  Initialize pointers for row/var names                             */
/*--------------------------------------------------------------------*/
    _info->lsgrg_varnames = NULL; _info->lsgrg_rownames = NULL;
        _info->usernames = 0;
/*--------------------------------------------------------------------*/
/*  zero char mem stats for var/row names                             */
/*--------------------------------------------------------------------*/
    lsgrg_zero_memstats(_info);
/*--------------------------------------------------------------------*/
/*   set internal debug flag values                                   */
/*--------------------------------------------------------------------*/
    _info->dbg.grgsub      = 0;  _info->dbg.setvarname  = 0;
    _info->dbg.setrowname  = 0;  _info->dbg.memsetup    = 0;
    _info->dbg.setupj      = 0;  _info->dbg.memalloc    = 0;
    _info->dbg.parshf      = 0;  _info->dbg.parshf0     = 0;
    _info->dbg.parshc      = 0;  _info->dbg.parshc0     = 0;
    _info->dbg.parsh0      = 0;  _info->dbg.globals     = 0;
    _info->dbg.algorithm   = 0;
/*--------------------------------------------------------------------*/
/*   initialize startup/run/shutdown flags                            */
/*--------------------------------------------------------------------*/
    _info->lsgrg_setup      = 0;      _info->lsgrg_run = 0;
    _info->lsgrg_shutdown   = 0;
    _info->lsgrg_setup_done = 0;
/*--------------------------------------------------------------------*/
/*   initialize grgtest testharness fields which need it              */
/*--------------------------------------------------------------------*/
    _info->grgtest.s_termination[0] = '\0';
    _info->grgtest.s_problem[0]     = '\0';
    _info->grgtest.s_datetime[0]    = '\0';
    _info->grgtest.s_termin[0]      = '\0';
    _info->grgtest.s_alg[0]         = '\0';
    _info->grgtest.s_termination[0] = '\0';

    _info->grgtest.obj = 0.0;
    _info->grgtest.sinf = 0.0;
    _info->grgtest.time = 0.0;

    _info->grgtest.ncols   = 0;
    _info->grgtest.nrows   = 0;
    _info->grgtest.itns    = 0;
    _info->grgtest.ninf    = 0;
    _info->grgtest.ngcomps = 0;
    _info->grgtest.info    = 0;


    _info->nvars_varnames = 0;   _info->nrows_rownames = 0;

} /* end of lsgrg_setdefaults */
void LSGRGCLASS lsgrg_initialize_1time(LsgrgInfo *_info)
{
double eps;
/*-------------------------------------------------------------------*/
/*     lsgrg_initialize_1time sets values which need be set only     */
/*     once per execution                                            */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*  SET MACHINE DEPENDENT PARAMETERS                                 */
/*                                                                   */
/*  EPS IS MACHINE ACCURACY FOR DOUBLE PRECISION VARIABLES.          */
/*  PLINFY IS 'PLUS INFINITY' FOR GIVEN MACHINE.                     */
/*  PLZERO IS POSITIVE VALUE NEAREST ZERO WITHOUT                    */
/*   BEING ZERO ON MACHINE.                                          */
/*  ** 3/98 jcp ** move initialization of plinfy to lsgrg0.c to make */
/*    value available to user prior to grgsub call                   */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*        CALCULATE EPS - MACHINE PRECISION                          */
/*-------------------------------------------------------------------*/
        eps = lsgrg_machine_precision();

        _info->eps = eps;
        _info->eps0 = pow(eps, 4.0/5.0);  /* constants derived from machine */
        _info->eps1 = pow(eps, 2.0/3.0);  /* precision                      */
        _info->eps2 = pow(eps, 1.0/2.0);
        _info->eps3 = pow(eps, 1.0/3.0);
        _info->eps4 = pow(eps, 1.0/4.0);
        _info->eps5 = pow(eps, 1.0/5.0);

        _info->tolz = _info->eps2;                /* lsgrg tolerances */
        _info->tolx =  _info->eps2 > 1.0e-6 ?  _info->eps2: 1.0e-6 ;
/*-------------------------------------------------------------------*/
/*   other one-time only initializations
/*-------------------------------------------------------------------*/
        _info->dluprm[1] = 1.0e-1;
        _info->dluprm[2] = 0.0e0;
        _info->dluprm[3] = 0.0e0;
        _info->dluprm[4] = 0.0e0;

        _info->nnlequ = 0;

        _info->pivtol[1] = 1.0e-2;
        _info->pivtol[2] = 5.0e-2;
        _info->pivtol[3] = 1.0e-1;
        _info->pivtol[4] = 5.0e-1;
        _info->pivtol[5] = 8.0e-1;
        _info->pivtol[6] = 9.0e-1;
        _info->pivtol[7] = 9.5e-1;

        _info->ipvtpt = 3;
        _info->rtf = 1.0e0;

        _info->jacobian_growth_allowance = J_GROWTH_ALLOWANCE;

        return;

} /* end of lsgrg_initialize_1time  */

void LSGRGCLASS lsgrg_setcomputedparms(LsgrgInfo *_info)
{
   _info->maxtab = _info->maxtbu;            _info->memcg = _info->maxcg;
   _info->maxh   = _info->maxr;              _info->mcgm1 = _info->memcg - 1;
   _info->eplast = _info->epnewt;
   if( _info->epinit == 0.0 ) _info->epinit = _info->epnewt;
   _info->ipr3 = _info->ipr - 1;

}
void LSGRGCLASS lsgrg_zero_alg_globals(LsgrgInfo *_info)
{
/*-----------------------------------------------------------------*/
/*  this function explicitly zeros all 'global' vars for the       */
/*  algorithm routines.  since these were previously external,     */
/*  they were initialized by default to 0                          */
/*                                                                 */
/*                                                                 */
/*-----------------------------------------------------------------*/
     int i;

     _info->epscom.eps0 = 0;
     _info->epscom.eps1 = 0;
     _info->epscom.eps2 = 0;
     _info->epscom.eps3 = 0;
     _info->epscom.eps4 = 0;
     _info->epscom.eps5 = 0;
     _info->cgbk.modcg = 0;
     _info->cgbk.memcg = 0;
     _info->cgbk.mcgm1 = 0;
     _info->cgbk.icgcnt = 0;
     _info->cgbk.icgptr = 0;
     _info->cgbk.hscale = 0;
     _info->bestbk.stpbst = 0;
     _info->bestbk.objbst = 0;
     _info->bestbk.step = 0;
     _info->bestbk.stepmx = 0;
     _info->bestbk.truobj = 0;
     _info->redph.trubst = 0;
     _info->congrg.nsear0 = 0;
     _info->congrg.lvlast = 0;
     _info->counts.nftn = 0;
     _info->counts.ngrad = 0;
     _info->counts.nminv = 0;
     _info->counts.nnfail = 0;
     _info->counts.ncalls = 0;
     _info->counts.nit = 0;
     _info->counts.nbs = 0;
     _info->counts.nstepc = 0;
     _info->counts.ndub = 0;
     _info->optblk.maxim = 0;
     _info->optblk.hrdbnd = 0;
     _info->optblk.subset = 0;
     _info->optblk.multsb = 0;
     _info->optblk.gfeas = 0;
     _info->optblk.useph0 = 0;
     _info->optblk.warmst = 0;
     _info->optblk.penlty = 0;
     _info->dfblk.dfail = 0;
     _info->dimen.m = 0;
     _info->dimen.n = 0;
     _info->dimen.mp1 = 0;
     _info->dimen.npmp1 = 0;
     _info->dimen.nbmax = 0;
     _info->dimen.nnbmax = 0;
     _info->dimen.npnbmx = 0;
     _info->dimen.mpnbmx = 0;
     _info->dimen.nrtot = 0;
     _info->dirgrg.cond = 0;
     _info->dirgrg.update = 0;
     _info->dirgrg.nsupp = 0;
     _info->zblck.condmx = 0;
     _info->zblck.nblck = 0;
     _info->infbk.info = 0;
     _info->pardat.kderiv = 0;
     _info->ingrg.epinit = 0;
     _info->ingrg.eplast = 0;
     _info->ingrg.epdeg = 0;
     (_info->inout.title)[0] = '\0';

#ifdef FILE_IO_ENABLED
     _info->iounit.ioin = 0;
     _info->iounit.ioout = 0;
     _info->iounit.iodump = 0;
     _info->iounit.ioerr = 0;
     _info->iounit.ioterm = 0;
#endif

     _info->limits.epboun = 0;
     _info->limits.epnewt = 0;
     _info->limits.epspiv = 0;
     _info->limits.itlim = 0;
     _info->logblk.move = 0;
     _info->logblk.restrt = 0;
     _info->logblk.drop = 0;
     _info->logblk.varmet = 0;
     _info->logblk.conjgr = 0;
     _info->logblk.resetp = 0;
     _info->slpobj.slope = 0;
     _info->sernew.edf = 0;
     _info->sernew.trunew = 0;
     _info->mfact.edfper = 0;
     _info->mfact.stpper = 0;
     _info->mfact.rtnmul = 0;
     _info->misc.maxh = 0;
     _info->misc.nsear = 0;
     _info->misc.jp = 0;
     _info->misc.lv = 0;
     _info->misc.jqq = 0;
     _info->mngrg.epstop = 0;
     _info->mngrg.limser = 0;
     _info->mngrg.nstop = 0;
     _info->mngrg.ierr = 0;
     _info->mngrg.ipn4 = 0;
     _info->mngrg.ipn5 = 0;
     _info->mngrg.ipn6 = 0;
     _info->mngrg.iper = 0;
     _info->iters.nitr = 0;
     _info->iters.ndeg = 0;
     _info->iters.nph0it = 0;
     _info->iters.nph1ls = 0;
     _info->degn.idglim = 0;
     _info->nintbk.nb = 0;
     _info->nintbk.nobj = 0;
     _info->nintbk.ninf = 0;
     _info->nintbk.nsuper = 0;
     _info->nintbk.ipr3 = 0;
     _info->nintbk.ncand = 0;
     _info->nintbk.ipr = 0;
     _info->bind.nbc = 0;
     _info->bind.nnbc = 0;
     _info->ph1bk.phmult = 0;
     _info->ph1bk.ph1eps = 0;
     _info->ph1bk.initph = 0;
     _info->srchlg.uncon = 0;
     _info->srchlg.fail = 0;
     _info->srchlg.jstfes = 0;
     _info->srchlg.mxstep = 0;
     _info->srchlg.unbd = 0;
     _info->srchlg.succes = 0;
     _info->srchlg.unconp = 0;
     _info->supblk.sbchng = 0;
     _info->supblk.basbnd = 0;
     _info->cbmode.smode = 0;
     _info->cbmode.prsmod = 0;
     _info->cbmode.newpt = 0;
     _info->cbmode.phase0 = 0;
     _info->cbmode.compgr = 0;
     _info->cbmode.pertub = 0;
     _info->tols.eps = 0;
     _info->tols.plinfy = 0;
     _info->tols.plzero = 0;
     _info->tols.tolx = 0;
     _info->tols.tolz = 0;
     _info->cbscnt.ncbs = 0;
     _info->cbscnt.nser = 0;
     _info->cbscnt.nreinv = 0;
     _info->cbscnt.nsame = 0;
     _info->cbscnt.nrser = 0;
     _info->cbscnt.nicond = 0;
     _info->cbscnt.nsing = 0;
     _info->cbscnt.ifrcnt = 0;
     _info->intggr.ngetgr = 0;
     _info->chq.maxtab = 0;
     _info->chq.ltab = 0;
     _info->chq.lptr = 0;
     _info->memory.lgrad = 0;
     _info->memory.lbinv = 0;
     _info->memory.lirn = 0;
     _info->memory.lbmap = 0;
     _info->memory.lmem = 0;
     _info->memory.np1 = 0;
     _info->lincnt.nlin = 0;
     _info->lincnt.nnlin = 0;
     _info->lincnt.nfix = 0;
     _info->nwtcnt.inwtk = 0;
     _info->nwtcnt.inwtfk = 0;
     _info->nwtcnt.inwtpt = 0;
      for(i=0;i<20;++i)
        (_info->nwtcnt.inwtlk)[i] = 0;
     _info->degcnt.idegk = 0;
     _info->degcnt.idegfk = 0;
     _info->degcnt.idegpt = 0;
      for(i=0;i<20;++i)
        (_info->degcnt.ideglk)[i] = 0;
     _info->degcnt.redpvt = 0;
     _info->sclobj.objscl = 0;
     _info->scal.iscale = 0;
     _info->scal.isclag = 0;
     _info->scal.scaled = 0;
     _info->scal.newscl = 0;
     _info->setin.lastz = 0;
     _info->setin.galloc = 0;
     _info->setin.getsiz = 0;
     _info->glberr.abort = 0;
     _info->dbug.debug = 0;
     _info->gctime.tgcomp = 0;
     _info->ph0tim.tph0 = 0;
     _info->ph0tim.tph1 = 0;
     _info->pvtblk.xpvpct = 0;
     _info->pvtblk.ibvblm = 0;
      for(i=0;i<4;++i)
        (_info->lsinvt.dluprm)[i] = 0;
      for(i=0;i<24;++i)
        (_info->lsinvt.iluprm)[i] = 0;
     _info->rednew.corner = 0;
     _info->rednew.xb = 0;
     _info->rednew.xstep = 0;
     _info->newsrc.iter = 0;
     _info->gmode.hbnd2 = 0;
     _info->gmode.jdbg = 0;
     _info->nwtim.tnewt = 0;
     _info->nwtim.tnewtx = 0;
     _info->gmsblk.gmserr = 0;
     _info->quadbk.a1 = 0;
     _info->quadbk.a2 = 0;
     _info->quadbk.a3 = 0;
     _info->quadbk.icon = 0;
     _info->quadbk.iquad = 0;
     _info->redser.ninfb = 0;
     _info->hesblk.gamma0 = 0;
     _info->hesblk.dirsc = 0;
     _info->initbk.init = 0;
     _info->initbk.lastcl = 0;
     _info->nph0.ninf0 = 0;
     _info->nph0.ndrop = 0;
     _info->equblk.nnlequ = 0;
     _info->equblk.lincon = 0;
     _info->ztol.xajtol = 0;
     _info->ztol.bigelt = 0;
     _info->ztol.biglin = 0;
     _info->temp.tsear = 0;
     _info->temp.tfact = 0;
     _info->contim.cbtime = 0;
     _info->itnggr.ngetgr = 0;
     _info->nzerog.nzgrad = 0;
     _info->nzerog.nznlin = 0;
     _info->nzerog.nzlin = 0;
     _info->jgrow.maxgro = 0;
     _info->cmax.colmax = 0;
     _info->gradtm.tgrad = 0;
     _info->stepbk.pstep = 0;
     _info->stepbk.funpr = 0;
     _info->zcond.cndtol = 0;
      for(i=0;i<7;++i)
        (_info->pivots.pivtol)[i] = 0;
     _info->pivots.ipvtpt = 0;
     _info->pivots.fixpiv = 0;
     _info->nzerob.nzbas = 0;
     _info->drvchk.nfddif = 0;
     _info->drvchk.chkdrv = 0;

     return;
}
void LSGRGCLASS lsgrg_zero_shell_globals(LsgrgInfo *_info)
{
/*   zeros all previously global vars in used in shell */
  int i;
/*  doubles */
  _info->cndtol = 0;
  _info->edf = 0;
  _info->epboun = 0;
  _info->epnewt = 0;
  _info->epspiv = 0;
  _info->rtf = 0;
  _info->epinit = 0;
  _info->eplast = 0;
  _info->epdeg = 0;
  _info->eps = 0;
  _info->tolx = 0;
  _info->tolz = 0;
  _info->eps0 = 0;
  _info->eps1 = 0;
  _info->eps2 = 0;
  _info->eps3 = 0;
  _info->eps4 = 0;
  _info->eps5 = 0;
  _info->epstop = 0;
  _info->funpr = 0;
  _info->pstep = 0;
  _info->phmult = 0;
  _info->ph1eps = 0;
  _info->xajtol = 0;
  _info->xpvpct = 0;
   for(i=0;i<8;++i)
     _info->pivtol[i] = 0;
   for(i=0;i<5;++i)
  _info->dluprm[i] = 0;
  _info->truobj = 0;
  _info->plinfy = 0;
  _info->plzero = 0;

/*      long int  */
  _info->inprnt = 0;
  _info->otprnt = 0;
/* 3/20/02 jcp ** add initialization of outres print level */
   _info->outres_print_level = OUTRES_PRINT_ALL;
  _info->iper = 0;
  _info->ipn4 = 0;
  _info->ipn5 = 0;
  _info->ipn6 = 0;
  _info->limser = 0;
  _info->nstop = 0;
  _info->iquad = 0;
  _info->colmax = 0;
  _info->ibvblm = 0;
  _info->idglim = 0;
  _info->ipr = 0;
  _info->ipr3  = 0;
  _info->ipvtpt = 0;
  _info->iscale = 0;
  _info->isclag = 0;
  _info->iter = 0;
  _info->itlim = 0;
  _info->kderiv = 0;
  _info->lbinv = 0;
  _info->lgrad = 0;
  _info->lirn = 0;
  _info->lmem = 0;
  _info->npmp2 = 0;
  _info->m = 0;
  _info->mp1 = 0;
  _info->mpnbmx = 0;
  _info->n = 0;
  _info->nbmax = 0;
  _info->nnbmax = 0;
  _info->npmp1 = 0;
  _info->npnbmx = 0;
  _info->nrtot = 0;
  _info->maxb = 0;
  _info->maxcg = 0;
  _info->maxr = 0;
  _info->maxtbu = 0;
  _info->nrows = 0;
  _info->nvars = 0;
  _info->maxh = 0;
  _info->mcgm1 = 0;
  _info->memcg = 0;
  _info->modcg = 0;
  _info->nfix = 0;
  _info->nlin = 0;
  _info->nnlin = 0;
  _info->nnlequ = 0;
  _info->nzbinv = 0;
  _info->nzgrad = 0;
  _info->nzlin = 0;
  _info->nznlin = 0;
  _info->nobj = 0;
  _info->ngcomp = 0;
  _info->nparsh = 0;
  _info->nbc = 0;
  _info->ninf = 0;
  _info->maxtab = 0;
  _info->nnewton = 0;
  _info->nnewton_itns = 0;
  _info->nnewton_fail = 0;
  _info->nnewton_stepc = 0;
  _info->ncbs = 0;
  _info->nser = 0;
  _info->nrser = 0;
  _info->nreinv = 0;
  _info->nicond = 0;
  _info->nsing = 0;
  _info->ifrcnt = 0;
  _info->nsame = 0;
  _info->nsuper = 0;
  _info->nb = 0;
  _info->nsear = 0;
  _info->nph1ls = 0;
  _info->ndeg = 0;
  _info->nbs = 0;
  _info->jacobian_growth_allowance = 0;

/*  ints  */
  _info->fixpiv = 0;
  _info->hrdbnd = 0;
  _info->hscale = 0;
  _info->maxim = 0;
  _info->multsb = 0;
  _info->abort = 0;

/*   LOGICAL32 */
  _info->useph0 = 0;
  _info->gfeas = 0;

   return;
}
#include "lsinfo.h"
/*
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
*/
/*-------------------------------------------------------------------*/
/*  File:  grgio.c  contains the io subsystem for lsgrgc             */
/*   author:  John C. Plummer  4/98                                  */
/*                                                                   */
/*  the purpose of FILE_IO_ENABLED is to inhibit inclusion of        */
/*  any FILE* definitions and fprintf calls when they are not        */
/*  desired                                                          */
/*                                                                   */
/*  the purpose of IO_ENABLED is to exclude any i/o declarations     */
/*  including stdio.h with all physical output is to be excluded     */
/*                                                                   */
/*  if IO_ENABLED is defined, calls to lsgrg_output are active       */
/*  if FILE_IO_ENABLED is defined, output is directed to one of      */
/*     LSGRG_LOG_UNIT, LSGRG_ERROR_UNIT, LSGRG_SCREEN                */
/*       (lsgrg_ioout)  (lsgrg_ioerr)     (lsgrg_ioterm)             */
/*                                                                   */
/*  **fixme** for now, no output if FILE_IO_ENABLED is not active    */
/*    windows interface will probably process buffer and pop up      */
/*    a msg box                                                      */
/*                                                                   */
/*  lsgrg_msg() -- writes contents of 'buffer' to LSGRG_LOG_UNIT     */
/*                                                                   */
/*  lsgrg_screen_msg() -- assigns contents of 'buffer' to            */
/*                         LSGRG_SCREEN                              */
/*                                                                   */
/*  lsgrg_error_msg() -- assigns contents of 'buffer' to             */
/*                         LSGRG_ERROR_UNIT                          */
/*                       and screen and log unit if active           */
/*  lsgrg_output() -- writes contents of 'buffer' to assigned        */
/*                         destination                               */
/*                                                                   */
/*   if FILE_IO_ENABLED, the following functions will be included    */
/*                                                                   */
/*   lsgrg_set_ioterm(FILE*)   set_lsgrg_ioout(FILE*)                */
/*   set_lsgrg_ioerr(FILE*)                                          */
/*                                                                   */
/*  if FLUSH_ALWAYS is defined                                       */
/*     file output will be flushed on each write to assist           */
/*     debugging                                                     */
/*-------------------------------------------------------------------*/
void LSGRGCLASS lsgrg_msg(IoInfo *ioinfo,char *buffer)
{
#ifdef IO_ENABLED
     if(ioinfo->ioout_enabled)
        lsgrg_output(ioinfo,buffer,LSGRG_LOG_UNIT);
#endif
   return;
}
void LSGRGCLASS lsgrg_error_msg(IoInfo *ioinfo,char *buffer)
{
/* send error msg to error_unit.  if physical io enabled */
/* if error output has been disabled, return            */
/* 4/26/01 jcp (1) fix bug with error output test       */
/*   (2) if screen output is enabled, error msgs go to */
/*       screen. write them to regular output also     */
/*                                                     */
/* 6/01 jcp if ioerr is assigned to ioout, send        */
/*   errmsg there. if assigned to ioterm, send errmsg  */
/*   there.  otherwise, send to ioerr                  */

#ifdef IO_ENABLED
   int done=0;
     if(ioinfo->error_output_enabled) {
        if(ioinfo->ioerr_to_ioout) {
           lsgrg_msg(ioinfo,buffer);
           done = 1;
        }
        if(ioinfo->ioerr_to_ioterm) {
           lsgrg_screen_msg(ioinfo,buffer);
           done = 1;
        }
        if(!done)
           lsgrg_output(ioinfo,buffer,LSGRG_ERROR_UNIT);
     }
#endif
   return;
}
void LSGRGCLASS lsgrg_screen_msg(IoInfo *ioinfo,char *buffer)
{
#ifdef IO_ENABLED
     if(ioinfo->screen_output_enabled)
        lsgrg_output(ioinfo,buffer,LSGRG_SCREEN);
#endif
   return;
}
void LSGRGCLASS lsgrg_output(IoInfo *io,char *buffer,int destination)
{
#ifdef FILE_IO_ENABLED
      FILE * out; char name[15];
      switch (destination) {
        case  LSGRG_LOG_UNIT:   out = io->lsgrg_ioout;  break;
        case  LSGRG_ERROR_UNIT: out = io->lsgrg_ioerr;  break;
        case  LSGRG_SCREEN:     out = io->lsgrg_ioterm; break;
        default:
            printf("\n Lsgrg internal error: unknown output destination"
              " code %d \n  buffer contents are: \n",destination);
            printf("%s",buffer);
            return;
     }
      if(out !=NULL) {
         fprintf(out,buffer);
#ifdef FLUSH_ALWAYS
         if (out != stdout && out != stderr) fflush(out);
#endif
      }
      else {
            switch(destination) {
             case  LSGRG_LOG_UNIT:   strcpy(name,"lsgrg_ioout");  break;
             case  LSGRG_ERROR_UNIT: strcpy(name,"lsgrg_ioerr");  break;
             case  LSGRG_SCREEN:     strcpy(name,"lsgrg_ioterm"); break;
            }
            printf("\n ERROR(lsgrg_output): Output Requested to [%s] but"
               " FILE pointer is NULL",name);
            printf("\n buffer is\n%s",buffer);
      }
#endif
      return;
}
void LSGRGCLASS lsgrg_lprint(LsgrgInfo *info,char *name,long lval)
{
/*  print line of form name ..... ival  with ival at col 50*/
#ifdef IO_ENABLED
       int ipos,i;
       ipos = sprintf(info->lsgrg_msgbuffer,"\n %s",name);
       info->lsgrg_msgbuffer[ipos] = ' ' ;
       for(i=ipos+1;i<=50;++i)
          sprintf(info->lsgrg_msgbuffer+i,".");
       sprintf((info->lsgrg_msgbuffer)+i-1," %12d",lval);
       lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
     return;
}
void LSGRGCLASS lsgrg_iprint(LsgrgInfo *info,char *name,int  ival)
{
/*  print line of form name ..... ival  with ival at col 50*/
#ifdef IO_ENABLED
       int ipos,i;
       ipos = sprintf(info->lsgrg_msgbuffer,"\n %s",name);
       info->lsgrg_msgbuffer[ipos] = ' ' ;
       for(i=ipos+1;i<=50;++i)
          sprintf(info->lsgrg_msgbuffer+i,".");
       sprintf((info->lsgrg_msgbuffer)+i-1," %12d",ival);
       lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
     return;
}
void LSGRGCLASS lsgrg_dprint(LsgrgInfo *info,char *name,double dval)
{
/*  print line of form name ..... dval  with dval at col 50*/
#ifdef IO_ENABLED
       int ipos,i;
       ipos = sprintf(info->lsgrg_msgbuffer,"\n %s",name);
       info->lsgrg_msgbuffer[ipos] = ' ' ;
       for(i=ipos+1;i<=50;++i)
          sprintf(info->lsgrg_msgbuffer+i,".");
       sprintf((info->lsgrg_msgbuffer)+i-1," %12.6g",dval);
       lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
     return;
}
#ifdef FILE_IO_ENABLED
void LSGRGCLASS lsgrg_set_ioout(LsgrgInfo *info,FILE *file)
{
     info->io_info.lsgrg_ioout = file;
     if(file== NULL) {
         info->io_info.ioout_enabled = 0;
     }
     else
         info->io_info.ioout_enabled = 1;

     info->io_info.ioout_set_by_user = 1;
     return;
}
void LSGRGCLASS lsgrg_set_ioerr(LsgrgInfo *info,FILE *file)
{
     info->io_info.lsgrg_ioerr = file;
     if(file == NULL)
        info->io_info.error_output_enabled = 0;
     else
        info->io_info.error_output_enabled = 1;

     info->io_info.ioerr_set_by_user = 1;
     return;
}
void LSGRGCLASS lsgrg_set_ioterm(LsgrgInfo *info,FILE *file)
{
/*  setting ioterm to a non-NULL value will set screen_output flag to 1 */
     info->io_info.lsgrg_ioterm = file;
     if(file == NULL)
        info->io_info.screen_output_enabled = 0;
     else
        info->io_info.screen_output_enabled = 1;
     info->io_info.ioterm_set_by_user = 1;
     return;
}
#endif
void LSGRGCLASS lsgrg_disable_ioout(LsgrgInfo *info)
 { info->io_info.ioout_enabled = 0; }
void LSGRGCLASS lsgrg_enable_ioout(LsgrgInfo *info)
 { info->io_info.ioout_enabled = 1; }

void LSGRGCLASS lsgrg_disable_ioerr(LsgrgInfo *info)
 { info->io_info.error_output_enabled = 0; }
void LSGRGCLASS lsgrg_enable_ioerr(LsgrgInfo *info)
 {  info->io_info.error_output_enabled = 1; }

void LSGRGCLASS lsgrg_disable_ioterm(LsgrgInfo *info)
 { info->io_info.screen_output_enabled = 0; }
void LSGRGCLASS lsgrg_enable_ioterm(LsgrgInfo *info)
 {  info->io_info.screen_output_enabled = 1; }
void LSGRGCLASS lsgrg_quiet_all_output(LsgrgInfo *info)
 {
    lsgrg_set_otprnt(info,(long) 0);
    lsgrg_set_inprnt(info,(long) 0);
    lsgrg_set_printlevel(info,(long) 0);
 }
#include "lsinfo.h"

#define IOINFO  &(info->io_info)
#define LSGRG_MSGBUFFER info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(info->lsexit_jmpbuf)
/***********************************************************************
 **                 FILE:              grgtime.c                       *
 **                 AUTHOR:            John C. Plummer                 *
 **                 LAST UPDATE:       20 Mar   2002                   *
 **                 LAST CHANGE:                                       *
 **                                                                    *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **   Timing functions for lsgrgc                                      *
 **   macro in grgglobl.h TIMING_ENABLED will cause activation of      *
 **   clock calls                                                      *
 **                                                                    *
 **  lsgrg_timer()  -- new timer call, returns elapsed time as double  *
 **  lsgrg_timestats(int, double)   -- handles timing statistics       *
 **   1st arg T_INIT     ==> initialize timing subsystem               *
 **           T_PRINT    ==> causes output of timing statistics        *
 **           T_GCOMP    ==> adds value of 2nd arg to gcomp time       *
 **           T_CONSBS   ==> adds value of 2nd arg to consbs time      *
 **           T_PARSH    ==> adds to gradient evaluation time          *
 **           T_NEWT     ==> adds to total newton time                 *
 **           T_NEWTONLY ==> adds to newton time without gcomp         *
 **           T_SEARCH   ==> adds to search mode time                  *
 **           T_INVERT   ==> adds to invert time                       *
 **           T_TOTAL    ==> adds to total time                        *
 **                                                                    *
 **                                                                    *
 **                                                                    *
 **                                                                    *
 ***********************************************************************
 * */
double LSGRGCLASS lsgrg_timer()
{
/* if timing is enabled, lsgrg_timer returns elapsed time in seconds */
/* return value is double                                            */

#ifdef TIMING_ENABLED
        return (double) ((double)clock()   / (double) CLOCKS_PER_SEC);
#else
        return 0.0;
#endif

}
void LSGRGCLASS lsgrg_timestats(LsgrgInfo *info,
                                int command, double tvalue)
{
/*  ** 3/20/02 jcp ** if command is to print, check outres_print_level */
/*   to disable timestats if row/col only print is requested.  this is */
/*   for final output from oqgrg                                       */
//double t_total_chk, t_consbs_chk;
      switch (command) {
         case T_INIT:
                info->tim.gcomp  = 0.0; info->tim.consbs  = 0.0; info->tim.parsh  = 0.0;
                info->tim.newton = 0.0; info->tim.newton_nogcomp = 0.0;
                info->tim.search = 0.0;  info->tim.invert = 0.0; info->tim.total = 0.0;
                info->tim.phase0 = 0.0; info->tim.phase1 = 0.0;
                info->grgtest.time = 0.0;
                break;
         case T_GCOMP:     info->tim.gcomp          += tvalue; break;
         case T_CONSBS:    info->tim.consbs         += tvalue; break;
         case T_PARSH:     info->tim.parsh          += tvalue; break;
         case T_NEWT:      info->tim.newton         += tvalue; break;
         case T_NEWTONLY:  info->tim.newton_nogcomp += tvalue; break;
         case T_SEARCH:    info->tim.search         += tvalue; break;
         case T_INVERT:    info->tim.invert         += tvalue; break;
         case T_PHASE0:    info->tim.phase0         += tvalue; break;
         case T_PHASE1:    info->tim.phase1         += tvalue; break;
         case T_TOTAL:     info->tim.total          += tvalue;
                           info->grgtest.time       += tvalue;
                           break;
      }
/*    if(command == T_GCOMP)
         printf("\n gcomp elapsed time = %12.6e",tvalue); */

      if(command != T_PRINT) return;
/*  ** 3/20/02 jcp ** dont print time stats if outres is printing */
/*   only row/col info                                            */
      if(info->outres_print_level == OUTRES_PRINT_ROWCOL) return;
/*---------------------------------------------------------------------*/
/*   print statistics if io is enabled                                 */
/*   fix denominators for percentages if times are too small           */
/*---------------------------------------------------------------------*/
#ifdef IO_ENABLED
        t_total_chk  = info->tim.total;
        t_consbs_chk = info->tim.consbs;
        if(info->tim.total  <= 0.0 ) t_total_chk = 1.0;
        if(info->tim.consbs <= 0.0 ) t_consbs_chk = 1.0;

        sprintf( LSGRG_MSGBUFFER,
         "\n\n Timing Information (All Times In Seconds):\n" );
        lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);

         sprintf( LSGRG_MSGBUFFER,
        "                                                     Time      %%Total"
        "\n Total CPU Time:                             %12.4f    %8.2f%% "
        "\n Total Time Evaluating Functions (GCOMP      %12.4f    %8.2f%% "
        "\n Total Time Evaluating Gradients             %12.4f    %8.2f%% "
        "\n"
        "\n Total CONSBS Time:                          %12.4f    %8.2f%% "
        "\n Total NEWTON Time:                          %12.4f    %8.2f%% "
        "\n Total NEWTON Time (Excluding GCOMP):        %12.4f    %8.2f%% ",
            info->tim.total, 100.0,info->tim.gcomp, (100.0*info->tim.gcomp/t_total_chk),
            info->tim.parsh,  (100.0*info->tim.parsh/t_total_chk),
            info->tim.consbs, (100.0*info->tim.consbs/t_total_chk),
            info->tim.newton, (100.0*info->tim.newton/t_total_chk),
            info->tim.newton_nogcomp , (100.0*info->tim.newton_nogcomp/t_total_chk)  );
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);

         sprintf( LSGRG_MSGBUFFER,"\n"
        "\n Time Spent in Phase 0                       %12.4f    %8.2f%%"
        "\n Time Spent in Phase 1                       %12.4f    %8.2f%%",
            info->tim.phase0,  (100.0*info->tim.phase0/t_total_chk),
            info->tim.phase1, (100.0*info->tim.phase1/t_total_chk) );
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);

        sprintf( LSGRG_MSGBUFFER, "\n"
        "\n Statistics for CONSBS\n" );
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        sprintf( LSGRG_MSGBUFFER,
        "\n Total CONSBS Time:                          %12.4f"
        "\n      Time Spent in Search Mode:             %12.4f    %8.2f%%"
        "\n      Time Spent Inverting Basis:            %12.4f    %8.2f%%",
          info->tim.consbs, info->tim.search, (100.0*info->tim.search/t_consbs_chk),
          info->tim.invert,  (100.0*info->tim.invert/t_consbs_chk) );
          lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
#endif

        return;

} /* end of function lsgrg_timestats*/
#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"

#define IOINFO  &(_info->io_info)
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)

/*
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 */

void /*FUNCTION*/ LSGRGCLASS modjac(LsgrgInfo *_info,long int ldummy, long int nnewnz, long int current_col,
         double current_jac_col[], double *check_col, long int he[],
         long int ha[], double a[], long int *max_new_nnz, long int nbr_cols,
         double zero_tol, long int nbr_rows, long int ipr, LOGICAL32 debug,
         long int *colmax, long int hel[], long int iobjpt[])
{
long int i, i_, irow,  k, k1, k2, k_, last_col_to_print, nnz;

/*--------------------------------------------------------------------------*/
/*    subroutine modjac makes room for additional                           */
/*    nonzeros in the current column and shifts the rest                    */
/*    of the jacobian to the right                                          */
/*  *** change log ***                                                      */
/*   ** 2/99 jcp ** change check_col from store of coefficients to          */
/*    store of indices to allow for expansion only and prevent              */
/*    skipping of previous nonzeros which become zero at the same           */
/*    time new nonzeros are added.                                          */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*                check available space for expansion                       */
/*--------------------------------------------------------------------------*/
        if(ipr > 2 || debug)
           lsgrg_msg(IOINFO,"\n .... entry modjac .....");

        if( nnewnz > *max_new_nnz ){
           if(ipr > 2)
              lsgrg_msg(IOINFO,
            "\n ...Note: Nbr of New Jacobian nonzeros has exceeded the"
            "\n          the growth allowance for nonzeros not present"
            "\n          at the initial point.  Initiating reallocation");
           lsgrg_errorexit(JUMPBUF,_LSGRG_JAC_OVERFLOW);
        }

        last_col_to_print = min( nbr_cols, 100 );  /* limit debug output */

#ifdef IO_ENABLED
        if( debug || ipr > 3 ){
            sprintf(LSGRG_MSGBUFFER,
              "\n ncols = %6ld nrows = %6ld"
              "\n  last_col_to_print = %6ld\n\n",
              nbr_cols, nbr_rows, last_col_to_print );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
        if( debug || ipr > 4 ){
            sprintf(LSGRG_MSGBUFFER,
              "\n .... current column: %6ld\n\n", current_col );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
            for( i = 1; i <= nbr_rows; i++ ){
                i_ = i - 1;
                sprintf(LSGRG_MSGBUFFER, "  row = %6ld gcol = %20.14g\n",
                  i, current_jac_col[i_] );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
            }
            for( i = current_col; i <= last_col_to_print; i++ ){
                i_ = i - 1;
                k1 = he[i_];
                k2 = he[i_ + 1] - 1;
                sprintf(LSGRG_MSGBUFFER,
                  " ..... current configuration of column %6ld"
                  "\n  col start = %6ld col end = %6ld\n",i, k1, k2 );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                for( k = k1; k <= k2; k++ ){
                    k_ = k - 1;
                    irow = ha[k_];
                    sprintf(LSGRG_MSGBUFFER,
                      " %6ld row = %6ld a(k) = %20.14g\n",
                      k, irow, a[k_] );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
            }
        }
#endif
#ifdef IO_ENABLED
        if( ipr > 1 ){
            sprintf(LSGRG_MSGBUFFER,
              "\n %3ld New Nonzeros Detected in Jacobian Column"
              " %6ld\n     Space Remaining for %6ld New Nonzeros",
              nnewnz, current_col, *max_new_nnz );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        *max_new_nnz = *max_new_nnz - nnewnz; /* update growth allowance */
/*------------------------------------------------------------------------*/
/*       move all subsequent columns right by nnewnz elements             */
/*------------------------------------------------------------------------*/
        nnz = he[nbr_cols] - 1;

#ifdef IO_ENABLED
        if( ipr > 1 || debug)  {
            sprintf( LSGRG_MSGBUFFER,
               "\n     Current Number of Nonzeros = %8ld",nnz );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        k1 = he[current_col];
#ifdef IO_ENABLED
        if(debug ) {
           sprintf(LSGRG_MSGBUFFER,"\n.... nnz = %10d, k1 = %10d",nnz,k1);
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif
        for( i = nnz; i >= k1; i-- ) {
                i_ = i - 1;
                a[i_ + nnewnz] = a[i_];
                ha[i_ + nnewnz] = ha[i_];

#ifdef IO_ENABLED
         if(debug )  {
           sprintf(LSGRG_MSGBUFFER,"\n updated ha[%d] = %10d from ha[%d] = %10d",
              (i_+nnewnz),ha[i_+nnewnz], i_,ha[i_]);
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
         }
#endif
        }

        /*  save start/end of current column before updating he */
        k1 = he[current_col - 1];
        k2 = he[current_col] - 1;

        nnz = nnz + nnewnz;
        for( i = current_col + 1; i <= nbr_cols; i++ ){
                i_ = i - 1;

#ifdef IO_ENABLED
                if(debug ) {
                  sprintf(LSGRG_MSGBUFFER,"\n.. updating he[%d] = %d",i_,he[i_]);
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
                he[i_] = he[i_] + nnewnz;
                hel[i_] = hel[i_] + nnewnz;
#ifdef IO_ENABLED
                if(debug )  {
                    sprintf(LSGRG_MSGBUFFER,"\n..    new he[%d] = %d",i_,he[i_]);
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        }
        he[nbr_cols] = nnz + 1;

#ifdef IO_ENABLED
        if( ipr > 1 || debug )  {
           sprintf( LSGRG_MSGBUFFER,
            "\n     New Number of Nonzeros = %8ld", nnz );
           lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

/*------------------------------------------------------------------------*/
/*     place new current column                                           */
/*     first copy all previously allocated nonzeros  into check_col       */
/*------------------------------------------------------------------------*/
        /*        k1 = he[current_col - 1];
                  k2 = he[current_col] - 1; */

#ifdef IO_ENABLED
        if(debug > 1)  {
          sprintf(LSGRG_MSGBUFFER,"\n k1 = %10d ha[k1-1] = %10d k2 = %10d ha[k2-1] = "
            "%10d",k1,ha[k1-1],k2,ha[k2-1]);
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif
/*------------------------------------------------------------------------*/
/*  ** jcp 2/99 ** new column expansion code to ensure that columns only  */
/*     grow, use indices rather than values                               */
/*------------------------------------------------------------------------*/
        for(i=0;i< nbr_rows;++i) check_col[i] = 0.0;
        for(k=k1; k<=k2; ++k) {   /* implant indices of existing nonzeros */
            k_ = k - 1;
            irow = ha[k_];
            check_col[irow-1] = irow;
        }
#ifdef IO_ENABLED
        if(debug) {
          lsgrg_msg(IOINFO,"\n... i_check_col -- existing indices");
          for(i=1;i<=nbr_rows;++i)  {
             irow = (long) check_col[i-1]; /* warning */
            sprintf(LSGRG_MSGBUFFER,"\n      row = %5d check_col = %5d",
                    i,irow);
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
          }
        }
#endif
        for(i=1; i <= nbr_rows; ++i) {  /* implant indices of all nonzeros */
                                        /* in current_jac_col              */
            if( fabs(current_jac_col[i-1]) > zero_tol)
                 check_col[i-1] = i;
        }
#ifdef IO_ENABLED
        if(debug) {
          lsgrg_msg(IOINFO,"\n... i_check_col -- after new indices");
          for(i=1;i<=nbr_rows;++i) {
             irow = (long) check_col[i-1]; /* warning */
            sprintf(LSGRG_MSGBUFFER,"\n     row = %5d check_col = %5d",
                    i,irow);
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
          }
        }
#endif
        k = k1;
        for(i=1;i<=nbr_rows;++i)  {  /* implant indices and coefficients into */
                                     /* a and ha                              */
            if(check_col[i-1] > 0.0)  {
               irow    = (long) check_col[i-1]; /* warning */
               ha[k-1] = irow;
               a[k-1]  = current_jac_col[i-1];
               ++k;
            }
        }
/*------------------------------------------------------------------------*/
/*   ** 2/99 jcp ** comment out old column expansion code                 */
/*  **fixme** expunge after shakedown                                     */
/*------------------------------------------------------------------------*/
/*
        for( k = k1; k <= k2; k++ ){
                k_ = k - 1;
                irow = ha[k_];
                check_col[irow - 1] = current_jac_col[irow - 1];
        }
/*
        /*c
         *c   now use check_col to fill up expanded column
         *c */
/*
        k = k1 - 1;
        for( irow = 1; irow <= nbr_rows; irow++ ){
                irow_ = irow - 1;
                if( fabs( check_col[irow_] ) > zero_tol ){
                        k = k + 1;
                        ha[k - 1] = irow;
                        a[k - 1] = current_jac_col[irow_];
                }
        }
*/
/*------------------------------------------------------------------------*/
/*  ------------------ end old column expansion code --------------       */
/*------------------------------------------------------------------------*/
        /*  get new length of updated column */
        k1 = he[current_col - 1];
        k2 = he[current_col] - 1;
        k = k2 - k1 + 1;

#ifdef IO_ENABLED
        if( ipr > 1 || debug ) {
           sprintf( LSGRG_MSGBUFFER,
            "\n     Updated Length of Column %d = "
            "%8ld\n",current_col,k );
           lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
             sprintf( LSGRG_MSGBUFFER,
              "\n     Updating Colmax: current value = %6ld", *colmax);
             lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif
        *colmax = max( *colmax, k );  /* update colmax */

#ifdef IO_ENABLED
           if(ipr > 1 || debug) {
             sprintf( LSGRG_MSGBUFFER,
              "\n                          New value = %6ld", *colmax);
             lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
           }
        if( ipr >= 1 || debug ) {
           sprintf( LSGRG_MSGBUFFER,
             "\n %3ld New Nonzeros Have Been Detected in Column %5ld",
             nnewnz, current_col );
           lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif
#ifdef IO_ENABLED
        if( debug || ipr > 4 ){
            for( i = current_col; i <= last_col_to_print; i++ ){
                 i_ = i - 1;
                 k1 = he[i_];
                 k2 = he[i_ + 1] - 1;
                 sprintf(LSGRG_MSGBUFFER,
                   " ..... new configuration of column %6ld  col start"
                   " = %6ld col end = %6ld\n",
                   i, k1, k2 );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                 for( k = k1; k <= k2; k++ ){
                     k_ = k - 1;
                     irow = ha[k_];
                     sprintf(LSGRG_MSGBUFFER, " %6ld row = %6ld a(k) = %20.14g\n",
                              k, irow, a[k_] );
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                 }
            }
        }
#endif
/*-----------------------------------------------------------------------*/
/*    call setobj to rebuild objective pointer set                       */
/*-----------------------------------------------------------------------*/
        setobj(_info, ha, he, iobjpt );

/*  write newline since some tranlsated routines use newline at end */
/*  of line rather than beginning                                   */
        if(ipr > 0) lsgrg_msg(IOINFO,"\n");
        return;
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS setobj(LsgrgInfo *_info,long int ihag[], long int iheg[], long int iobjpt[])
{
long int ind, ind_, jcol, jcol_;

        /*.......................................................................
         *
         *                      *** SETOBJ ***
         *
         *     *****Purpose:
         *     This subroutine calculates index pointers for the
         *     objective function gradient elements.  This allows direct
         *     access to the elements in the Jacobian array "GRAD".
         *
         *     *****Argument Description:
         *     On Input:
         *
         *     IHAG     - The array of row indecies for the sparse Jacobian.
         *     IHEG     - The array of start of column pointers.
         *
         *     On Output:
         *
         *     IOBJPT(I) = Location of the nonzero objective function gradient
         *                 element corresponding to variable "XI"
         *                 in the array "GRAD".
         *                 (0 implies that the corresponding element is 0.0)
         *
         *
         *     *****History:
         *     Written by:  Stuart H. Smith
         *                  Krannert School of Management
         *                  Purdue University
         *                  West Lafayette, IN  47901   (317-494-4531)
         *
         *     Date Last Modified:  JAN   24 1994
         *
         *----------------------------------------------------------------------- */
        if( _info->nintbk.ipr > 4 )
            lsgrg_msg(IOINFO," *** ENTERING SUBROUTINE: SETOBJ ....\n\n" );

        /*      --------------------------------------------------------
         *     | For each column, find Objective Element id it exists
         *      --------------------------------------------------------
         * */
        for( jcol = 1; jcol <= _info->dimen.n; jcol++ ){
                jcol_ = jcol - 1;
                iobjpt[jcol_] = 0;
                for( ind = iheg[jcol_]; ind <= (iheg[jcol_ + 1] - 1); ind++ ){
                        ind_ = ind - 1;
                        if( ihag[ind_] == _info->nintbk.nobj ){
                                iobjpt[jcol_] = ind;
                                goto L_100;
                        }
                }
L_100:
                ;
        }


#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 4 ){
               lsgrg_msg(IOINFO,
                " THE OBJECTIVE FUNCTION GRADIENT POINTERS - IOBJPT:\n" );
                for( ind = 1; ind <= _info->dimen.n; ind++ ){
                     sprintf(LSGRG_MSGBUFFER, "%5ld", iobjpt[ind - 1] );
               lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO,
                   "\n *** EXITING SUBROUTINE: SETOBJ...\n\n" );
        }
#endif

        return;

        /*  ** End of SETOBJ *** */
} /* end of function */
void /*FUNCTION*/ LSGRGCLASS prtha(LsgrgInfo *_info,long *ha, char *msg)
{
#ifdef IO_ENABLED
  long i;
  sprintf(LSGRG_MSGBUFFER,"\n.... print of ha: %s\n",msg);
  lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
  for(i=0;i<=100;++i) {
     if(i%4 ==0) sprintf(LSGRG_MSGBUFFER,"\n");
     sprintf(LSGRG_MSGBUFFER," ha[%3d]= %10d",i,ha[i]);
     lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
  }
#endif
  return;
  }
void /*FUNCTION*/ LSGRGCLASS dumpcol(LsgrgInfo *_info,char *msg, long *ha, long icol, long *he)
{
#ifdef IO_ENABLED
   int i,j=0;
   sprintf(LSGRG_MSGBUFFER,
     "\n.. dumpcol: col = %d msg = %s",icol,msg);
   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
   sprintf(LSGRG_MSGBUFFER,
     "\n he[icol-1] = %10d, he[icol]-1 = %10d\n",(he[icol-1]),
       (he[icol]-1));
   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
   for(i=he[icol-1];i<=he[icol]-1;++i) {
     if(++j%4==0) printf("\n");
     sprintf(LSGRG_MSGBUFFER," ha[%3d] = %10d ",i-1,ha[i-1]);
     lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
     if(ha[i-1]<0) exit(0);
  }
   lsgrg_msg(IOINFO,"\n");
#endif
   return;
}
#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"

#define IOINFO  &(_info->io_info)
#define LSGRG_MSGBUFFER _info->io_info.lsgrg_msgbuffer
#define JUMPBUF *(_info->lsexit_jmpbuf)
/***********************************************************************
 **                 FILE:        GRGLSITN FORTRAN                      *
 **                 AUTHOR:      STUART SMITH                          *
 **                 LAST UPDATE: 06 MAY   1998                         *
 **                                                                    *
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
 **                                                                    *
 **      This file contains the following routines for LSGRG2, a       *
 **  large scale implementation of GRG2:                               *
 **                                                                    *
 **  ROUTINE           DESCRIPTION                        LAST UPDATE  *
 ** ---------         -------------                      ------------- *
 **  GRGITN     GRG Driver routine                        22 MAy 1996  *
 **  ITRLOG     Prints an iteration log                   23 FEB 1991  *
 **                                                                    *
 **  REVISIONS:                                                        *
 **   05/22/96 SHS - Changed bound perturbation logic in GRGITN.       *
 **      Now all basic variables (including slacks) at bounds that     *
 **      are forced into the bound along the tangent vector are        *
 **      perturbed outward.                                            *
 **                                                                    *
 *    06/13/96 LSL - Added calls to stats routine for test harness
 *    and computed true KT error for output purposes. Includes change
 *    to itrlog
 *
 *    05/06/98 SHS - Added PAIJ,IPROW,IPCOL, and IPMAP arrays for use
 *    with new Analytical Parsh subroutine
 *  ** 10/98 jcp ** add ipow, powi functions to C version
 *  ** 1/98 jcp** adapted io to allow selective and total stripping
 *   of io code. NOTE: direct calls to lsgrg_msg do not need to be
 *   ifdef'ed since the acutal output in lsgrg_msg is already ifdefed
 *   on IO_ENABLED
 ***********************************************************************
 * */
int  /*FUNCTION*/ LSGRGCLASS grgitn(LsgrgInfo *_info,long int ibmap[], double zmem[], double grad[],
         long int ihag[], long int iheg[], long int ihegl[], double r[],
         double alb[], double ub[], double x[], double gradf[], double g[],
         double v[], double d[], double u[], double gbest[], double xbest[],
         double ascale[], long int inbv[], long int iub[], double rowb[],
         double colb[], long int ibc[], long int ibv[], long int istat[],
         double xstat[], long int ivstat[], double xb1[], double xb2[],
         double xb3[], double gradfp[], double dbnd[], double gg[], double rr[],
         double y[], long int icols[], long int iobjpt[], long int icand[],
         long int inlin[], long int irank[], double cdnum[], long int *tablst,
         double albo[], double ubo[], long int lksame[], double *ycg,
         double *scg, double *cgscr, long int ibvbct[], double paij[],
         long int iprow[], long int ipcol[], long int ipmap[])
{
#define TABLST(I_,J_)   (*(tablst+(I_)*(_info->chq.maxtab)+(J_)))
#define YCG(I_,J_)      (*(ycg+(I_)*(_info->dimen.n)+(J_)))
#define SCG(I_,J_)      (*(scg+(I_)*(_info->dimen.n)+(J_)))
#define CGSCR(I_,J_)    (*(cgscr+(I_)*(_info->cgbk.mcgm1-(0)+1)+(J_)))
LOGICAL32 degen, droped, infes, kt, newepn, reinvt;
long int i, i_, idegct, ii, ii_, inc, iprhd3, iprhld, /* irc , */ isave,
         istop, j, j_, jr, k, linect, msgcg, nbnd, nbndp, nbviol, nfail,
         ninfsv, nj, nloop, nsbsav, nstopm;
double ainf, bnd, /* epdegl, */ epdegx, epnorg, epstpo, err, kterr, objtst,
         phhold, tfinal, theta, tst, tstart;

        /*         **************************************************
         *         **************************************************
         *         ***                LSGRG2                      ***
         *         ***           COPYRIGHT   1991                 ***
         *         ***            LEON S. LASDON                  ***
         *         ***                  &                         ***
         *         ***            ALLAN D. WAREN                  ***
         *         ***                  &                         ***
         *         ***            STUART SMITH                    ***
         *         **************************************************
         *         **************************************************
         * */


        /* .....................................................................
         *
         *     INITIALIZE PERFORMANCE COUNTERS: */

/* ** 11/98 jcp ** do not set nftn=2 or ngrad=1 for setup calls */
/*  new shell already counts these                              */

        _info->counts.ncalls = 0;
        /*             NCALLS =  Number of Newton calls */
        _info->counts.nit = 0;
        /*             NIT    =  Cumulative No. of Newton iterations */
        _info->counts.nftn = 1;
/*      if( !_info->setin.galloc )
                _info->counts.nftn = 2;  */

        /*             NFTN   =  Total no. of GCOMP calls
         *                       (1 IN SETGR, 1 IN TABLIN) */
        _info->counts.ngrad = 0;

/*      if( !_info->setin.galloc )
                _info->counts.ngrad = 1; */

        /*             NGRAD  =  No of PARSH (gradient) calls
         *                       (1 in SETGR) */
        _info->misc.nsear = 0;
        /*             NSEAR  =  No. of one dimensional searches
         *
         *******************************************************************
         *     These next 2 variables are not used.  They are intended for
         *     Drud's ideas of tightening convergence of Newton as you
         *     approach the optimum.
         *
         *XXXX NFESLS=0
         *             NFESLS = Number of times SEARCH could not tighten final
         *                      point to feasibility
         *XXXX NQFAIL=0
         *             NQFAIL = Number of times QFIT failed to fit a convex
         *                      quadratic to data
         ******************************************************************** */
        _info->iters.nitr = -1;
        /*             NITR   =  No. of iterations (including degenerate ones) */
        _info->iters.ndeg = 0;
        /*             NDEG   =  No. of Degenerate iterations */
        _info->iters.nph1ls = 0;
        /*             NPH1LS =  No. of PHASE-1 line searches performed */
        _info->iters.nph0it = 0;
        /*             NPH0IT =  No. of PHASE-0 Newton iterations */
        istop = 0;
        /*             ISTOP  =  No. of consecutive times relative change in
         *                       objective function is less than EPSTOP */
        _info->counts.nbs = 0;
        /*             NBS    =  No. of times basic variable violates bound */
        _info->counts.nnfail = 0;
        /*             NNFAIL =  No. times Newton failed to converge */
        _info->counts.nstepc = 0;
        /*             NSTEPC =  No. times step size cut back when Newton fails
         * */
        _info->cbscnt.ncbs = 0;
        /*             NCBS   =  No. times CONSBS is called */
        _info->intggr.ngetgr = 0;
        /*             NGETGR =  No. time GETGR is called */
        _info->cbscnt.nser = 0;
        /*             NSER   =  No. times CONSBS is called in SEARCH MODE
         *                       (i.e. find a whole new basis) */
        _info->cbscnt.nrser = 0;
        /*             NRSER  =  No. times CONSBS is called in RESTRICTED SEARCH
         *                       MODE(i.e. just alter previous basis) */
        _info->cbscnt.nreinv = 0;
        /*             NREINV =  No. times CONSBS is called in RE-INVERSION MODE
         *                       (i.e. with the same basis to invert) */
        _info->cbscnt.nsing = 0;
        /*             NSING  =  No. times the invert routines found a singular
         *                       basis and forced CONSBS into SEARCH MODE */
        _info->cbscnt.nicond = 0;
        /*             NICOND =  No. times CONSBS found an Ill-Conditioned basis
         * */
        _info->cbscnt.nsame = 0;
        /*             NSAME  =  No. times CONSBS chooses same basis as previous
         *                       iteration chose.
         * */
        _info->cbscnt.ifrcnt = 0;
        /*             IFRCNT =  No. times fractional change forces SEARCH MODE
         *
         *
         *
         * */
#ifdef IO_ENABLED
        if( _info->dbug.debug || _info->nintbk.ipr >= 5 ){
                lsgrg_msg(IOINFO,
                 "\n ***** ENTERING SUBROUTINE: GRGITN\n" );
                lsgrg_msg(IOINFO, "  ISTAT(" );
                for( j = 1; j <= _info->dimen.mp1; j++ ){
                        sprintf(LSGRG_MSGBUFFER, "%2ld) = %3ld", j, istat[j - 1] );
                        lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
                lsgrg_msg(IOINFO, " ISTAT(\n" );
                lsgrg_msg(IOINFO, "  IVSTAT(" );
                for( j = 1; j <= _info->dimen.n; j++ ){
                        sprintf(LSGRG_MSGBUFFER, "%2ld) = %3ld", j, ivstat[j - 1] );
                }
                sprintf(LSGRG_MSGBUFFER, " IVSTAT(\n" );
        }
/* *** comment this out, printgrad0 requires an LsgrgInfo struture ptr

     if(_info->dbug.debug) {
        lsgrg_printgrad0( 1, "Jacobian at entry grgitn ",
                          _info->dimen.n    , _info->dimen.n,_info->dimen.m,_info->dimen.m, _info->dimen.mp1,
                          _info->nzerog.nzgrad, _info->nzerog.nzlin, _info->nzerog.nznlin,
                          iheg, ihag,grad ,x );
         lsgrg_echo_globals("Entry grgitn");
    }
*/
#endif

        _info->cbmode.smode = !_info->optblk.warmst;
        if( _info->optblk.warmst ){
                _info->nintbk.nb = _info->dimen.mp1;
        } else{
                _info->nintbk.nb = 0;
        }
        _info->bind.nbc = 0;
        _info->nintbk.nsuper = 0;
        _info->limits.epnewt = _info->ingrg.epinit;
        epstpo = _info->mngrg.epstop;
        if( _info->ingrg.eplast < _info->ingrg.epinit )
                _info->mngrg.epstop = 10.0e0*_info->mngrg.epstop;
        phhold = _info->ph1bk.ph1eps;
        _info->bestbk.step = 0.0e0;
        nloop = 0;
        _info->nintbk.ipr = _info->nintbk.ipr3 + 1;
        linect = 60;
        iprhld = _info->nintbk.ipr;
        iprhd3 = _info->nintbk.ipr3;
        _info->glberr.abort = FALSE;
        _info->scal.scaled = FALSE;
        epnorg = _info->limits.epnewt;
        _info->sclobj.objscl = 1.0e0;
        for( i = 1; i <= _info->dimen.npmp1; i++ ){
                i_ = i - 1;
                ascale[i_] = 1.0e0;
        }

        /************************************************************
         *    THESE PERCENTAGES ARE TEMPORARILY HARD-CODED
         ************************************************************
         * */
        _info->mfact.edfper = 1.0e-02;
        _info->mfact.stpper = 1.0e-03;
        _info->mfact.rtnmul = 1.0e03;

        /*      -----------------------------------
         *     : Save the original variable bounds
         *      -----------------------------------
         * CHG 05/22/96 SHS - Save original constraint bounds too. */

        /*SHS      DO 7 I=1,N */
        for( i = 1; i <= _info->dimen.npmp1; i++ ){
                i_ = i - 1;
                albo[i_] = alb[i_];
                ubo[i_] = ub[i_];
        }

        /*       ---------------------------------------
         *      | This is return point for EPNEWT loop
         *       ---------------------------------------
         * */
L_10:
        ;
        nloop = nloop + 1;
        _info->degcnt.idegk = 0;
        _info->degcnt.idegfk = 0;
        _info->degcnt.idegpt = 0;
        _info->nwtcnt.inwtk = 0;
        _info->nwtcnt.inwtfk = 0;
        _info->nwtcnt.inwtpt = 0;
        _info->degcnt.redpvt = FALSE;
        newepn = FALSE;
        nbndp = 0;
        epdegx = _info->ingrg.epdeg;
        _info->congrg.nsear0 = _info->misc.nsear;
        if( _info->limits.epnewt <= _info->ingrg.eplast )
                _info->mngrg.epstop = epstpo;
        _info->scal.newscl = _info->scal.iscale != 0 && nloop == 1;
        _info->cbmode.newpt = TRUE;
        _info->srchlg.succes = TRUE;
        _info->cbmode.compgr = TRUE;
        _info->logblk.drop = FALSE;
        droped = FALSE;
        reinvt = FALSE;
        _info->cbmode.pertub = FALSE;
        _info->logblk.move = TRUE;
        _info->logblk.restrt = TRUE;
        _info->logblk.resetp = FALSE;
        _info->logblk.varmet = TRUE;
        _info->logblk.conjgr = FALSE;
        _info->srchlg.unbd = FALSE;
        _info->srchlg.jstfes = FALSE;
        degen = FALSE;
        _info->srchlg.uncon = FALSE;
        _info->supblk.basbnd = FALSE;
        _info->supblk.sbchng = FALSE;
        idegct = 0;
        /*XXX  NSUPER=0
         *XXX  NSUPP=0 */
        _info->dirgrg.nsupp = _info->nintbk.nsuper;
        _info->redph.trubst = g[_info->nintbk.nobj - 1];
        _info->chq.ltab = 0;
        _info->chq.lptr = 0;
        _info->mngrg.ierr = 0;

        nfail = 0;
        msgcg = 1;
        _info->bestbk.stpbst = 1.0e0;
        _info->misc.lv = 0;
        istop = 0;
        objtst = _info->tols.plinfy;
        _info->sernew.edf = 0.0e0;
        _info->bestbk.step = 0.0e0;
        _info->dirgrg.cond = 0.0e0;
        _info->dirgrg.update = 1;
        _info->nintbk.ninf = 0;
        ninfsv = 0;
        /*PEN  NB=0
         *PEN  NBC=0 */
        nstopm = _info->mngrg.nstop - 1;
        if( _info->mngrg.iper != 0 ){
                _info->nintbk.ipr = 1;
                _info->nintbk.ipr3 = 0;
        }
        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                x[_info->dimen.n + i_] = g[i_];
        }
        for( i = 1; i <= _info->dimen.npmp1; i++ ){
                i_ = i - 1;
                lksame[i_] = (long) 0.0e0; /* warning */
        }

        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                ibvbct[i_] = 0;
        }

        if( _info->scal.newscl )
                prscal(_info,istat, iheg, ihegl, ihag, grad, ascale, ivstat, inlin,
                 g, gg, gbest, colb, x, alb, ub, albo, ubo, iobjpt, paij,
                 iprow, ipcol, ipmap );

        /*      ------------------------------------------------------
         *     | Perform a PHASE-0 to find an intial feasible point
         *      ------------------------------------------------------
         * */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 )
                {
                sprintf(LSGRG_MSGBUFFER, "  CALLING PHASE-0 -- USEPH0 = %3c\n",
                 TorF(_info->optblk.useph0) );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        if( _info->optblk.useph0 ){
                phas0(_info,ibmap, zmem, grad, ihag, iheg, ihegl, ibv, x, alb,
                 ub, inbv, iub, ibc, g, inlin, ivstat, istat, dbnd, gg, d,
                 icols, (long*)gbest, colb, icand, cdnum, irank, ascale, ibvbct,
                 rowb, rr, &infes, (long*)xb1, &linect, iobjpt, paij, iprow,
                 ipcol, ipmap );
                if( _info->glberr.abort )
                        goto L_480;
                linect = 60;
        }
  /*    cputime( &tstart, &irc ); */
        tstart = lsgrg_timer();

        /*      ---------------------------------------
         *     |  This is return point for MAIN loop.
         *      ---------------------------------------
         *
         *     ----------------------------------------------------------------- */
L_40:
        ;
        /*     write (ioout,*)'at 40 newscl = ',newscl,' newpt = ',
         *    * newpt,' isclag = ', isclag, ' nitr = ',nitr */
        if( _info->scal.newscl ){
                prscal(_info,istat, iheg, ihegl, ihag, grad, ascale, ivstat, inlin,
                 g, gg, gbest, colb, x, alb, ub, albo, ubo, iobjpt, paij,
                 iprow, ipcol, ipmap );
                _info->cbmode.smode = TRUE;
        }
        /*lsl 3/30/98  newpt is set to true above to force CONSBS to set
         *    REINVT to FALSE, and avoid an error return
         *
         *
         *      --------------------------------------------------
         *     | Compute basis inverse, except when degenerate
         *      --------------------------------------------------
         *
         *     -----------------------------------------------------------------
         *     write(ioout,*)' before consbs newpt = ',newpt
         *
         *     chg 5/6/98 SHS - Add Analytical parsh arrays to CONSBS argument list
         * */
        consbs(_info,ibmap, zmem, grad, ihag, iheg, ihegl, ibv, x, alb, ub,
         inbv, iub, ibc, g, inlin, ivstat, istat, dbnd, gg, d, icols,
         (long*)gbest, colb, icand, cdnum, irank, ascale, ibvbct, &nsbsav,
         iobjpt, paij, iprow, ipcol, ipmap );
/*       printbasisInv("after consbs call 1"); */
        /*     ----------------------------------------------------------------- */
        if( _info->glberr.abort )
                goto L_480;
        if( _info->nintbk.ninf == 0 || _info->ph1bk.ph1eps == 0.0e0 )
                goto L_50;
        _info->ph1bk.initph = 1;
        ph1obj(_info,ibc, g, ub, alb );
        _info->ph1bk.initph = 0;
L_50:
        ;
        if( _info->misc.nsear != _info->congrg.nsear0 )
                goto L_100;

        /*      -----------------------------------------------------------
         *     | Initializations that must be done after first CONSBS call
         *     |   for each value of EPNEWT
         *      ----------------------------------------------------------- */
        _info->ph1bk.initph = 2;
        ph1obj(_info,ibc, g, ub, alb );
        ninfsv = _info->nintbk.ninf;
        _info->srchlg.jstfes = _info->nintbk.ninf == 0 || _info->srchlg.jstfes;
        _info->ph1bk.initph = 0;
        if( _info->nintbk.nb == 0 )
                goto L_70;
        for( i = 1; i <= _info->nintbk.nb; i++ ){
                i_ = i - 1;
                k = ibv[i_];
                xbest[i_] = x[k - 1];
        }

L_70:
        ;

        for( i = 1; i <= _info->dimen.mp1; i++ ){
                i_ = i - 1;
                gbest[i_] = g[i_];
        }
L_100:
        ;
        _info->iters.nitr = _info->iters.nitr + 1;
        if( _info->iters.nitr > 0 ){
                if( _info->iters.nitr == _info->mngrg.ipn4 )
                        _info->nintbk.ipr = 4;
                if( _info->iters.nitr == _info->mngrg.ipn5 )
                        _info->nintbk.ipr = 5;
                if( _info->iters.nitr == _info->mngrg.ipn6 )
                        _info->nintbk.ipr = 6;
                _info->nintbk.ipr3 = _info->nintbk.ipr - 1;
        }

        /*     Compute Reduced Gradient
         *
         *     ----------------------------------------------------------------- */
        redgra(_info,gradf, ibmap, zmem, grad, ihag, iheg, u, ibv, inbv, g,
         alb, ub, iobjpt );
        /*     ----------------------------------------------------------------- */
        if( _info->glberr.abort )
                goto L_480;
        if( _info->nintbk.ipr < 4 )
                goto L_140;
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                xstat[i_] = gradf[i_];
        }
        if( !_info->optblk.maxim || _info->nintbk.ninf != 0 )
                goto L_130;
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                xstat[i_] = -xstat[i_];
        }
L_130:
        if( _info->nintbk.ipr < 0 )
                goto L_140;
#ifdef IO_ENABLED
        lsgrg_msg(IOINFO, " REDUCED GRADIENT IS \n " );
        for( i = 1; i <= _info->dimen.n; i++ ){
                sprintf(LSGRG_MSGBUFFER, "%16.8e", xstat[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
        lsgrg_msg(IOINFO, "\n" );
#endif

L_140:
        ;
        /*     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         *
         *      Check if any of the STOP criteria are satisfied
         *
         *     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         *     Test if KUHN-TUCKER conditions satisfied
         * */
/* L_155: unreferenced??? */
        kt = FALSE;
        kterr = 0.e0;
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                ii = inbv[i_];
                tst = gradf[i_];
                if( ii <= _info->dimen.n )
                        goto L_160;

                /*         Slack variables
                 * */
                if( istat[ii - _info->dimen.n - 1] == 1 )
                        goto L_190;
                goto L_165;

L_160:
                ;
                /*C
                 *C        Skip fixed variables
                 *C */
                if( ivstat[ii - 1] < 0 )
                        goto L_190;

L_165:
                if( iub[i_] == 0 )
                        goto L_180;
                if( iub[i_] == 1 )
                        goto L_170;

                /*  NONBASIC VAR AT LOWER BOUND
                 * */
                err = 0.e0;
                if( tst < 0.e0 )
                        err = -tst;
                if( err > kterr )
                        kterr = err;
                goto L_190;
L_170:
                ;


                /*  NONBASIC VAR AT UPPER BOUND
                 * */
                err = 0.e0;
                if( tst > 0.e0 )
                        err = tst;
                if( err > kterr )
                        kterr = err;
                goto L_190;

                /*  NONBASIC VARIABLE FREE
                 * */
L_180:
                err = fabs( tst );
                if( err > kterr )
                        kterr = err;

L_190:
                ;
        }

        /*     TEST IF KUHN-TUCKER ERROR LESS THAN EPSTOP
         * */
        if( kterr < _info->mngrg.epstop )
                kt = TRUE;
        /*C */
/* L_200: unreferenced??? */
        itrlog(_info,_info->dimen.n, _info->dimen.mp1, _info->nintbk.nobj, _info->iters.nitr, nsbsav, _info->nintbk.nb,
         _info->nintbk.ninf, &_info->nintbk.ipr, iprhld, &_info->nintbk.ipr3, iprhd3, &linect,
         degen, _info->optblk.maxim, _info->dirgrg.update, _info->dirgrg.cond, _info->zblck.condmx,
         g, x, _info->bestbk.step, kterr );
        /*C
         *C    If K-T conditions are satisfied. go check if EPNEWT at final value
         *C */
        if( kt ){
                _info->infbk.info = 0;
                lsgrg_set_termin_code(_info,"KTC");
#ifdef IO_ENABLED
                  sprintf(LSGRG_MSGBUFFER,
                  "\nTERMINATION CRITERION MET.\n KUHN-TUCKER CONDITIONS"
                  " SATISFIED TO WITHIN %12.5e AT CURRENT POINT \n",
                  _info->mngrg.epstop );
                  lsgrg_set_termin_msg(_info,LSGRG_MSGBUFFER);
                  if( _info->nintbk.ipr > 0 ){
                      lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                      if(_info->io_info.screen_output_enabled)
                         lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
                         linect = linect + 1;
                  }
#endif
                goto L_480;
        } else if( _info->srchlg.jstfes && _info->optblk.gfeas ){
                   _info->infbk.info = 7;
                   lsgrg_set_termin_code(_info,"Feas-UsrTerm");
#ifdef IO_ENABLED
                    sprintf(LSGRG_MSGBUFFER,
                      "\nINITIAL FEASIBLE POINT FOUND -- ALGORITHM"
                      " TERMINATED BY USER OPTION\n NUMBER OF PHASE-1"
                      "  LINE SEARCHES = %8ld\n",
                         _info->iters.nph1ls );
                    lsgrg_set_termin_msg(_info,LSGRG_MSGBUFFER);
                    if( _info->nintbk.ipr > 0 ){
                       lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                      if(_info->io_info.screen_output_enabled)
                       lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
                           linect = linect + 1;
                   }
#endif
                goto L_480;
        }
        /*     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         *      Unconditional STOP if number of linear searches > LIMSER
         *     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         * */
        if( _info->misc.nsear >= _info->mngrg.limser ){
                    lsgrg_set_termin_code(_info,"Iterations");
#ifdef IO_ENABLED
                    sprintf(LSGRG_MSGBUFFER,
                      "\nNUMBER OF COMPLETED ONE-DIMENSIONAL SEARCHES"
                      " = LIMSER = %5ld\n ...  OPTIMIZATION TERMINATED.\n",
                     _info->iters.nitr );
                    lsgrg_set_termin_msg(_info,LSGRG_MSGBUFFER);
                    if( _info->nintbk.ipr > 0 ){
                        lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                      if(_info->io_info.screen_output_enabled)
                        lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
                            linect = linect + 1;
                   }
#endif
                _info->infbk.info = 3;
                _info->mngrg.ierr = 11;
                _info->glberr.abort = TRUE;
                goto L_480;
        }
/*=======================================================================*/
/*  ** 4/27/02 jcp ** add termination on time limit exceeded             */
/*=======================================================================*/
        if( (lsgrg_timer() - tstart) >= _info->maxtime ){
                    lsgrg_set_termin_code(_info,"TimeLimit");
#ifdef IO_ENABLED
                    sprintf(LSGRG_MSGBUFFER,
                      "\nTIME LIMIT OF %10.2f (sec) EXCEEDED "
                      " ...  OPTIMIZATION TERMINATED.\n",
                     _info->maxtime );
                    lsgrg_set_termin_msg(_info,LSGRG_MSGBUFFER);
                    if( _info->nintbk.ipr > 0 ){
                        lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                      if(_info->io_info.screen_output_enabled)
                        lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
                            linect = linect + 1;
                   }
#endif
                _info->infbk.info = _LSGRGALG_TIME_LIMIT;
                _info->mngrg.ierr = 11; /* leave same as itn limit?? */
                _info->glberr.abort = TRUE;
                goto L_480;
        }



/*=======================================================================*/

        /*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         *     Checks if relative change in Objective is less than EPSTOP for
         *     NSTOP consecutive iterations.
         *     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
        _info->cbmode.pertub = FALSE;
        if( !_info->srchlg.succes || degen )
                goto L_240;
        /*XXXX IF (.NOT. UNCON) GO TO 230 */
        if( fabs( g[_info->nintbk.nobj - 1] - objtst ) > (1.0e0 + fabs( objtst ))*
         _info->mngrg.epstop )
                goto L_230;

        /*     Fractional change too small. Count how often consecutively.
         *
         *     If infeasible,  and have reduced NINF, then continue
         * */
        if( _info->nintbk.ninf > 0 && _info->nintbk.ninf < ninfsv ){
                ninfsv = _info->nintbk.ninf;
                goto L_230;
        }

        istop = istop + 1;
        if( istop < nstopm )
                goto L_240;

        if( istop == nstopm )
                goto L_220;

        /*     ISTOP > NSTOP so go check if EPNEWT at final value */

         lsgrg_set_termin_code(_info,"FracChange");

#ifdef IO_ENABLED
         sprintf(LSGRG_MSGBUFFER,
            "\nTOTAL FRACTIONAL CHANGE IN OBJECTIVE LESS"
            " THAN %12.5e\n           FOR%4ld CONSECUTIVE ITERATIONS \n",
            _info->mngrg.epstop, istop );
         lsgrg_set_termin_msg(_info,LSGRG_MSGBUFFER);
         if( _info->nintbk.ipr > 0 ) {
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                      if(_info->io_info.screen_output_enabled)
            lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
           linect = linect + 2;
         }
#endif

        _info->mngrg.ierr = 1;
        _info->infbk.info = 1;
        goto L_480;

        /*     Fractional change too small NSTOP times. If infeasible,
         *     but have reduced NINF, then continue
         * */
L_220:
        if( _info->nintbk.ninf > 0 && _info->nintbk.ninf < ninfsv ){
                ninfsv = _info->nintbk.ninf;
                goto L_230;
        }

        /*     Fractional change too small NSTOP times.
         *     Try selecting a new basis and dropping any possible constraints
         * */
        _info->cbmode.smode = TRUE;
        _info->cbmode.compgr = TRUE;
        _info->cbscnt.ifrcnt = _info->cbscnt.ifrcnt + 1;

#ifdef IO_ENABLED
         sprintf(LSGRG_MSGBUFFER,
             "\n FRACTIONAL CHANGE IN OBJECTIVE LESS THAN %12.5e"
             "\n  FOR%4ld CONSECUTIVE ITERATIONS -- FORCING CONSBS"
             " INTO SEARCH MODE \n",
             _info->mngrg.epstop, istop );
         lsgrg_set_termin_msg(_info,LSGRG_MSGBUFFER);
         if( _info->nintbk.ipr > 0 ) {
             lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
             if(_info->io_info.screen_output_enabled)
                 lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
         }
#endif
        /*@@   IF (ISCALE .NE. 0) THEN
         *@@      NEWSCL=.TRUE.
         *@@      WRITE(IOOUT,636)
         *@@      CALL PRSCAL(ISTAT,IHEG,IHEGL,IHAG,GRAD,ASCALE,IVSTAT,
         *@@  *               INLIN,G,GG,GBEST,COLB,X,ALB,UB,ALBO,UBO)
         *@@   ENDIF
         *@@   DROP=.TRUE.
         *@@   GO TO 40 */



        /*     CHG 5/6/98 SHS - Add analytical PARSH arrays to CONSBS argument list
         *     ----------------------------------------------------------------- */
        consbs(_info,ibmap, zmem, grad, ihag, iheg, ihegl, ibv, x, alb, ub,
         inbv, iub, ibc, g, inlin, ivstat, istat, dbnd, gg, d, icols,
         (long*)gbest, colb, icand, cdnum, irank, ascale, ibvbct, &nsbsav,
         iobjpt, paij, iprow, ipcol, ipmap );
/*       printbasisInv("after consbs call 2"); */
        /*     ----------------------------------------------------------------- */
        if( _info->glberr.abort )
                goto L_480;

        /*     ------------------------------------------------------------------ */
        redgra(_info,gradf, ibmap, zmem, grad, ihag, iheg, u, ibv, inbv, g,
         alb, ub, iobjpt );
        /*     ------------------------------------------------------------------
         * */
        if( _info->glberr.abort )
                goto L_480;
        _info->cbmode.smode = FALSE;
        _info->logblk.drop = TRUE;
        goto L_240;
L_230:
        istop = 0;
        _info->cbmode.smode = FALSE;
        objtst = g[_info->nintbk.nobj - 1];
L_240:
        ;

        /*     ******************************************************************
         *
         *     Compute search direction for Superbasics
         * */
L_250:
        ;

        /*     ------------------------------------------------------------------ */
        direc(_info,gradf, r, gradfp, iub, inbv, d, istat, y, ivstat, ibv,
         lksame, ycg, scg, cgscr );
        /*     ------------------------------------------------------------------
         *XXXX NSBSAV = NSUPER
         * */
        if( _info->glberr.abort )
                goto L_480;
        if( _info->dfblk.dfail && droped )
                goto L_420;
        if( _info->dfblk.dfail )
                goto L_4322;

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 4 ) {
            lsgrg_msg(IOINFO, " DIRECTION VECTOR IS \n " );
            for( i = 1; i <= _info->nintbk.nsuper; i++ ){
                    sprintf(LSGRG_MSGBUFFER, "%13.5e", d[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
            }
            lsgrg_msg(IOINFO, "\n" );
        }
#endif

        if( _info->nintbk.nb == 0 )
                goto L_300;

        /*     Compute Tangent Vector V
         *
         *     ------------------------------------------------------------------ */
        tang(_info,ibmap, zmem, grad, ihag, iheg, inbv, v, d, rr );
        /*     ------------------------------------------------------------------
         * */
#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 5 ) {
            lsgrg_msg(IOINFO, " TANGENT VECTOR IS \n " );
            for( i = 1; i <= _info->nintbk.nb; i++ ){
                    sprintf(LSGRG_MSGBUFFER, "%13.5e", v[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
            }
            lsgrg_msg(IOINFO, "\n" );
        }
#endif

        /*     Find JP, index of first basic variable to hit a bound
         *     ICOLS(J) holds bound status of basic variable "J":
         *            0  -> Variable can be moved in tanget direction
         *           -1  -> Variable at lower bound and tanget points into bound
         *            1  -> Variable at upper bound and tangent points into bound.
         *
         *     ------------------------------------------------------------------ */
        chuzr(_info,alb, ub, x, v, g, ibv, icols, ibc, &_info->misc.jp, &theta, &_info->logblk.move );
        /*     ------------------------------------------------------------------
         *
         *      --------------------------------
         *     | Update degenerate counts
         *      --------------------------------
         * */
        inc = 0;
        if( !_info->logblk.move )
                inc = 1;
        _info->degcnt.idegpt = _info->degcnt.idegpt + 1;
        if( _info->degcnt.idegpt > 20 )
                _info->degcnt.idegpt = 1;
        if( _info->degcnt.idegk < 20 ){
                _info->degcnt.idegk = _info->degcnt.idegk + 1;
                isave = 0;
        } else{
                isave = _info->degcnt.ideglk[_info->degcnt.idegpt - 1];
        }
        _info->degcnt.idegfk = _info->degcnt.idegfk + inc - isave;
        _info->degcnt.ideglk[_info->degcnt.idegpt - 1] = inc;

        if( _info->logblk.move )
                goto L_300;

        /*     Degenerate and no MOVE in basics is possible
         * */
        jr = ibv[_info->misc.jp - 1];

#ifdef IO_ENABLED
        if( _info->nintbk.ipr >= 3 ) {
             sprintf(LSGRG_MSGBUFFER,
              "\nBASIS DEGENERATE--VARIABLE  %5ld LEAVING BASIS \n",
             jr );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        _info->misc.lv = _info->misc.jp;
        degen = TRUE;
        _info->iters.ndeg = _info->iters.ndeg + 1;
        idegct = idegct + 1;

        /* -------------------------------------------------------------------------
         * CHG 05/22/96 SHS -Perturb all basic variables (including slacks) at bounds
         *      that CHUZR determines will violate the bound with current tangent
         *      vector. (Determined by ICOLS in CHUZR)
         *
         *    If degenerate for 'IDGLIM" iters then perturb current bounds for
         *      basic variable at bounds. Note that:
         *     1) Fixed structural variables are ignored, but equality constraints
         *        (fixed slacks) are perturbed.
         *     2) Only the bound identified by CHUZR as being violated (if it exists)
         *        is perturbed.
         *     3) For variables satisfying 2), the current bound (which may already
         *        be perturbed) is perturbed outward by EPDEGX (in a relative sense).
         *
         *    If NBNDP exceeds limit (currently hard-coded at 5), then try
         *      a search mode before giving up.
         *---------------------------------------------------------------------
         *
         * CHG 05/22/96 - New perturbation logic follows:
         * */
        if( idegct < _info->degn.idglim )
                goto L_290;

#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 0 ) {
             sprintf( LSGRG_MSGBUFFER,
              "\nDEGENERATE FOR%5ld STEPS.  PROBABLY CYCLING. \n",
             idegct );
             lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif

        idegct = 0;

        /*  NOTE: Limit is hard-coded for now.  Should be changed to a user parameter.
         * */
        if( nbndp < 5 ){
#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 0 ) {
                    sprintf( LSGRG_MSGBUFFER,
                    "\n PERTURBING VARIABLE BOUNDS -- EPDEG = %14.7e\n",
                    epdegx );
                    lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

                nbnd = 0;
                for( j = 1; j <= _info->nintbk.nb; j++ ){
                        j_ = j - 1;
                        i = ibv[j_];

                        /*        Skip fixed strucutral variables
                         * */
                        if( i <= _info->dimen.n ){
                                if( ivstat[i - 1] < 0 )
                                        goto L_298;
                        }

                        if( icols[j_] == -1 ){
                                nbnd = nbnd + 1;
                                bnd = alb[i - 1];
                                alb[i - 1] = alb[i - 1] - epdegx*(1.0 + fabs( alb[i - 1] ));
                                if( alb[i - 1] < -_info->tols.plinfy )
                                        alb[i - 1] = -_info->tols.plinfy;
#ifdef IO_ENABLED
                                if( _info->nintbk.ipr >= 4 )
                                        {
                                   sprintf(LSGRG_MSGBUFFER,
     "\n  PERTURBING LOWER BOUND FOR VARIABLE %5ld OLD = %15.7e NEW = %15.7e\n",
                                         i, bnd, alb[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                        }
#endif
                        } else if( icols[j_] == 1 ){
                                nbnd = nbnd + 1;
                                bnd = ub[i - 1];
                                ub[i - 1] = ub[i - 1] + epdegx*(1.0 + fabs( ub[i - 1] ));
                                if( ub[i - 1] > _info->tols.plinfy )
                                        ub[i - 1] = _info->tols.plinfy;


#ifdef IO_ENABLED
                                if( _info->nintbk.ipr >= 4 )
                                        {
                                        sprintf(LSGRG_MSGBUFFER,
   "\n  PERTURBING UPPER BOUND FOR VARIABLE %5ld OLD = %15.7e NEW = %15.7e\n",
                                         i, bnd, ub[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                                        }

#endif
                        }
L_298:
                        ;
                }

#ifdef IO_ENABLED
                if( _info->nintbk.ipr >= 2 )
                        {
                        sprintf(LSGRG_MSGBUFFER,
         " FINISHED PERTURBING BASIC BOUNDS: %5ld BOUNDS  CHANGED.\n",
                         nbnd );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                nbndp = nbndp + 1;

                /*        Clear Tabu list by resetting pointers.
                 * */
                _info->chq.ltab = 0;
                _info->chq.lptr = 0;
                _info->cbmode.pertub = TRUE;
                goto L_40;
        }

        /* END 05/22/96 - END of new logic
         *
         * CHG 05/22/96 SHS - Comment out old perturbation logic below
         *XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         *      IF (IDEGCT.LT.IDGLIM) GO TO 290
         *      IF (IPR .GT. 0) WRITE (IOOUT,800) IDEGCT
         *      IDEGCT = 0
         *XXXX IF(NBNDP .GE. 2 .OR. EPDEGX .LE. EPS) GO TO 299
         *XXX  IF(EPDEGX .GT. 1.01D-1 .OR. EPDEGX .LE. EPS) GO TO 299
         *     IF (NBNDP .GE. 10) GO TO 299
         *     IF (IPR .GT. 0) WRITE(IOOUT,802) EPDEGX
         *     DO 298 J=1,NB
         *        I=IBV(J)
         *        IF (I .GT. N) GO TO 298
         *        IF (IVSTAT(I) .LT. 0) GO TO 298
         *        IF (UBO(I) .GE. PLINFY) GO TO 297
         *         UB(I)=UBO(I)+EPDEGX*(1+DABS(UBO(I)))
         * 297     IF (ALBO(I) .LE. -PLINFY) GO TO 298
         *         ALB(I)=ALBO(I)-EPDEGX*(1+DABS(ALBO(I)))
         * 298     CONTINUE
         *
         *     NBNDP = NBNDP+1
         *     EPDEGL=EPDEGX
         *     EPDEGX=EPDEGX*10.0D0
         *     LTAB = 0
         *     LPTR = 0
         *????    SMODE=.TRUE.
         *     GO TO 40
         *XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         * END 05/22/96 SHS
         *
         *     Come to 299 if more than 10 perturbations performed.
         *     Try one final search mode before tightening EPNEWT at 475.
         * */
/* L_299: unreferenced??? */
        if( _info->cbmode.prsmod )
                goto L_475;

        if( _info->nintbk.ipr > 0 ) {
            lsgrg_msg(IOINFO,"\n FORCING CONSBS INTO SEARCH MODE \n" );
        }

        /*XX   IDEGCT = 0 */
        _info->degcnt.redpvt = TRUE;
        _info->cbmode.smode = TRUE;
        goto L_40;

        /*  PRINT ITERATION LOG
         * */
L_290:
        ;
        /*----------------------------------------------------------------------
         *     Exchange basic with some Superbasic and update BINV
         *
         *     ------------------------------------------------------------------ */
        chuzq(_info,grad, ihag, iheg, ibmap, zmem, v, d, inbv, alb, ub, x,
         ibv, iub, icand, ivstat, tablst, cdnum, irank );
        /*     ------------------------------------------------------------------
         *
         *     Set logicals for use by DIREC
         * */
        _info->logblk.restrt = TRUE;
        _info->logblk.drop = FALSE;
        _info->srchlg.uncon = FALSE;
        _info->supblk.sbchng = TRUE;
        _info->srchlg.mxstep = TRUE;
        _info->bestbk.step = 0.0;

        /*     Now go to begin new iteration for degenerate case
         * */
        goto L_100;
L_300:
        ;
        degen = FALSE;
        idegct = 0;
        _info->chq.ltab = 0;
        _info->chq.lptr = 0;

        /*     ------------------------------------------------------------------ */
        search(_info,ibmap, zmem, grad, ihag, iheg, d, x, g, ibv, v, xb1, xb2,
         xb3, inbv, alb, ub, xstat, gbest, xbest, ibc, rowb, colb, rr,
         ascale );
        /*     ------------------------------------------------------------------
         *
         *     If absolute value of X's is very small, change to 0 to avoid
         *     UNDERFLOW. */
        for( i = 1; i <= _info->dimen.n; i++ ){
                i_ = i - 1;
                if( fabs( x[i_] ) < _info->tols.eps )
                        x[i_] = 0.0e0;
        }
        if( _info->glberr.abort )
                goto L_480;

        /*xxx  NEWPT = NEWPT .OR. SUCCES */
        _info->cbmode.newpt = _info->srchlg.succes;
        _info->cbmode.compgr = _info->cbmode.newpt || _info->cbmode.compgr;

        if( _info->srchlg.succes )
                goto L_440;

        /*     TROUBLE -- No function decrese
         *
         *
         *      --------------------------------------------------------------
         *     | First,
         *     | try dropping constraints(and Gradient step ) if not done already
         *      --------------------------------------------------------------
         * */
        if( droped )
                goto L_420;
        if( _info->nintbk.ipr > 0 ) {
            lsgrg_msg(IOINFO,
              "\n LINE SEARCH FAILED TO FIND IMPROVED POINT"
               "-- DROPPING CONSTRAINTS AND USING REDUCED GRADIENT\n" );
        }

        _info->logblk.drop = TRUE;
        droped = TRUE;
        /*CXXX SMODE = .TRUE. */
        goto L_250;

        /*      --------------------------------------------------------------
         *     | Already Dropped, try reinverting current basis matrix with
         *     | larger pivot toerances for MA28 routines
         *      --------------------------------------------------------------
         * */
L_420:
        if( reinvt )
                goto L_425;
        if( _info->nintbk.ipr > 0 ) {
            lsgrg_msg(IOINFO,
              "\n LINE SEARCH FAILED TO FIND IMPROVED POINT"
              " -- REINVERTING BASIS MATRIX\n" );
        }

        reinvt = TRUE;
        _info->cbmode.smode = FALSE;
        goto L_40;

        /*      --------------------------------------------------------------
         *     | Already re-inverted, force CONSBS into SEARCH MODE
         *      --------------------------------------------------------------
         * */
L_425:
        ;
        if( _info->cbmode.prsmod )
                goto L_435;

        if( _info->nintbk.ipr > 0 ) {
            lsgrg_msg(IOINFO,
              "\n LINE SEARCH FAILED TO FIND IMPROVED POINT"
              " -- FORCING CONSBS INTO SEARCH MODE\n" );
        }

        _info->cbmode.smode = TRUE;
        goto L_40;

        /*      ------------------------------------------------------------
         *     | Finally, try tightening Newton tolerances
         *     : Go back to solve new problem with new EPNEWT
         *      ------------------------------------------------------------
         * */
L_435:
        ;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 0 )
                {
                sprintf( LSGRG_MSGBUFFER,
                 "\n LINE SEARCH FAILED TO FIND IMPROVED POINT"
                 " -- EPNEWT REDUCED TO %15.7e\n",
                 _info->limits.epnewt*1.0e-01 );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif

L_475:
        _info->limits.epnewt = _info->limits.epnewt*1.0e-01;
        /*@@@@ GO TO 4322 */
        if( _info->limits.epnewt < _info->epscom.eps2 )
                goto L_4322;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 0 )
                {
                sprintf( LSGRG_MSGBUFFER,
                 "\nRESTARTING WITH CONSTRAINT TOLERANCE"
                 " (EPNEWT) SET TO %12.5e\n",
                 _info->limits.epnewt );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        _info->cbmode.smode = FALSE;
        /*XXXX EDF = EDF*1.0D-01
         *CXXX GETFES = .TRUE. */
        if( _info->nintbk.ninf != 0 )
                g[_info->nintbk.nobj - 1] = _info->bestbk.truobj;
        newepn = TRUE;
        goto L_480;

        /*     NO Improvement in linesearch.
         *     and Newton tolerances are at the minimum
         *
         *     All Remedies Failed.
         * */
L_4322:
        _info->mngrg.ierr = 2;
        lsgrg_set_termin_code(_info,"AllRemedies");

#ifdef IO_ENABLED
         sprintf(LSGRG_MSGBUFFER,
         "\nALL REMEDIES HAVE FAILED TO FIND A BETTER POINT."
         "  PROGRAM TERMINATED. \n" );
         lsgrg_set_termin_msg(_info,LSGRG_MSGBUFFER);
        if( _info->nintbk.ipr > 0 ) {
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                      if(_info->io_info.screen_output_enabled)
            lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
           linect = linect + 2;
       }
#endif

        _info->infbk.info = 2;
        _info->glberr.abort = TRUE;
        goto L_480;

        /*     Line search found a better point
         * */
L_440:
        ;
        if( !_info->srchlg.unbd )
                goto L_450;

        /*     UNBOUNDED solution
         * */
        _info->mngrg.ierr = 20;
        lsgrg_set_termin_code(_info,"Unbounded");

#ifdef IO_ENABLED
         sprintf(LSGRG_MSGBUFFER,
         "\nSOLUTION UNBOUNDED--FUNCTION IMPROVING AFTER"
         " DOUBLING STEP %4ld TIMES\n",
         _info->counts.ndub );
         lsgrg_set_termin_msg(_info,LSGRG_MSGBUFFER);
        if( _info->nintbk.ipr > 0 ) {
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                      if(_info->io_info.screen_output_enabled)
            lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
           linect = linect + 1;
        }
#endif
        _info->infbk.info = 4;
        _info->glberr.abort = TRUE;
        goto L_480;
L_450:
        ;
        nfail = 0;
        _info->logblk.restrt = FALSE;
        _info->logblk.drop = FALSE;
        droped = FALSE;
        /*XXXX GETFES = .FALSE. */
        reinvt = FALSE;
        if( _info->srchlg.jstfes ){

    /*          cputime( &tfinal, &irc ); */
                tfinal = lsgrg_timer();

                _info->ph0tim.tph1 = _info->ph0tim.tph1 + tfinal - tstart;
                _info->iters.nph1ls = _info->misc.nsear;

#ifdef IO_ENABLED
                if( (_info->nintbk.ipr > 1) && _info->limits.epnewt != epnorg )
                        {
                        sprintf(LSGRG_MSGBUFFER,
                         "\nCONSTRAINT TOLERANCE HAS BEEN RESET TO"
                         " ITS ORIGINAL VALUE OF %12.5e\n",
                         epnorg );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                        }
#endif
                _info->limits.epnewt = epnorg;
                if( _info->limits.epnewt > _info->ingrg.eplast )
                        _info->mngrg.epstop = 10.e0*epstpo;
                _info->sernew.edf = 0.0e0;
        }

        /*     Check to see if we need to compute new scale factors
         *
         *     Only rescale if problem hasn't been perturbed */

        if( (_info->scal.isclag != 0) && (nbndp == 0) ){
                _info->scal.newscl = (_info->misc.nsear%_info->scal.isclag) == 0;
        } else{
                _info->scal.newscl = FALSE;
        }
        goto L_40;


        /*      ---------------------------------------------
         *     | Get to 480 when ready to terminate.
         *     | Start Checking what happened.
         *      ---------------------------------------------
         * */
L_480:
        ;

        /*      ----------------------------------------------------------
         *     |  Reset to original bounds and check for infeasibilities
         *      ----------------------------------------------------------
         * */
        nbviol = 0;
        if( nbndp > 0 ){

                /*        First, check the variable bounds.
                 * */
                for( i = 1; i <= _info->dimen.n; i++ ){
                        i_ = i - 1;
                        alb[i_] = albo[i_];
                        ub[i_] = ubo[i_];
                        bnd = alb[i_];
                        if( (alb[i_] - x[i_]) > _info->limits.epboun*(1.0e0 + fabs( alb[i_] )) ){
                                x[i_] = alb[i_];
                                nbviol = nbviol + 1;
                        } else if( (x[i_] - ub[i_]) > _info->limits.epboun*(1.0e0 + fabs( ub[i_] )) ){
                                x[i_] = ub[i_];
                                nbviol = nbviol + 1;
                        }
                }

                /*        If point has changed, re-evaluate the functions at new point
                 * */
                if( nbviol > 0 )
                        calfun(_info,g, x, ascale );

                /*        Now, do the constraint bounds
                 * */
                for( ii = 1; ii <= _info->dimen.mp1; ii++ ){
                        ii_ = ii - 1;
                        i = _info->dimen.n + ii;
                        alb[i - 1] = albo[i - 1];
                        ub[i - 1] = ubo[i - 1];
                        bnd = alb[i - 1];
                        if( (alb[i - 1] - g[ii_]) > _info->limits.epnewt*(1.0e0 + fabs( alb[i - 1] )) ){
                                x[i - 1] = alb[i - 1];
                                nbviol = nbviol + 1;
                        } else if( (g[ii_] - ub[i - 1]) > _info->limits.epboun*(1.0e0 +
                         fabs( ub[i - 1] )) ){
                                x[i - 1] = ub[i - 1];
                                nbviol = nbviol + 1;
                        }
                }
        }

        /*      -------------------------------------
         *     | Consider all the possible cases:
         *      -------------------------------------
         *
         *     If aborting, go to normal termination
         * */
        if( _info->glberr.abort )
                goto L_520;

        /*     If tightening constraint tolerance, go back and re-solve */

        if( newepn )
                goto L_10;

        /*     If perturbed solution is not feasible go back and solve the
         *     original problem.
         * */
        if( nbviol > 0 ){
#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 0 ) {
                    sprintf( LSGRG_MSGBUFFER,
                    " PERTURBED PROBLEM HAS %5ld BOUND VIOLATIONS"
                    "  -- RESETTING BOUNDS AND ATTEMPTING TO SOLVE\n",
                    nbviol );
                    lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }

#endif
                goto L_10;
        }



        /*      -----------------------------------------------------------------
         *     | If we terminated feasible, check if EPNEWT is at final value.
         *     | IF not, adjust EPNEWT and go back to re-solve problem.
         *      ------------------------------------------------------------------
         * */
        if( (((_info->limits.epnewt <= _info->ingrg.eplast) || (_info->nintbk.ninf != 0)) ||
         _info->optblk.gfeas) || (_info->lincnt.nfix == _info->dimen.n) )
                goto L_520;

        /*X    IF(INFO.NE.2)
         *X   *  CALL ITRLOG(N,MP1,NOBJ,NSEAR,NSUPER,NB,NINF,IPR,DEGEN,MAXIM,
         *X   *        COND,CONDMX,G,X,GRADF,STEP) */
        _info->limits.epnewt = _info->ingrg.eplast;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 0 ){
            sprintf( LSGRG_MSGBUFFER,
              "\nCONSTRAINT TOLERANCE HAS BEEN TIGHTENED TO ITS"
              " FINAL VALUE OF %12.5e\n",
              _info->limits.epnewt );
            linect = linect + 2;
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif
        /*CXXX EPSTOP = 0.1D0*EPSTOP */
        epnorg = _info->limits.epnewt;
        _info->ph1bk.ph1eps = 0.2e0;
        _info->cbmode.smode = TRUE;
        goto L_10;

        /*      ------------------------------------------------------------
         *     | Almost done - Must check to see if we terminated in PHASE-I
         *      ------------------------------------------------------------
         * */
L_520:
        ;


        /*      PRINT ITERATION LOG
         *
         *X    CALL ITRLOG(N,MP1,NOBJ,NSEAR,NSUPER,NB,NINF,IPR,DEGEN,MAXIM,
         *X   *      COND,CONDMX,G,X,GRADF,STEP) */

        /*      ------------------------------------
         *     | If feasible, we can finally quit!
         *      ------------------------------------ */

        if( _info->nintbk.ninf == 0 )
                goto L_540;

        /*      ------------------------
         *     | Solution Infeasible
         *      -----------------------
         * */
         lsgrg_set_termin_code(_info,"Infeasible");
#ifdef IO_ENABLED
         sprintf(LSGRG_MSGBUFFER,
           "\nFEASIBLE POINT NOT FOUND.  FINAL VALUE OF TRUE OBJECTIVE"
           " = %13.6e\n THE FOLLOWING%4ld CONSTRAINTS WERE IN VIOLATION"
           "  (EPNEWT = %12.5e):\n",
         _info->bestbk.truobj, _info->nintbk.ninf, _info->limits.epnewt );
         lsgrg_set_termin_msg(_info,LSGRG_MSGBUFFER);
         if( _info->nintbk.ipr > 0 ) {
            lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                      if(_info->io_info.screen_output_enabled)
            lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
            linect = linect + 2;
            /*X521 INFO=INFO+10 */
        }
#endif
        _info->infbk.info = 5;
        _info->mngrg.ierr = 9;

        for( j = 1; j <= _info->dimen.mp1; j++ ){
                j_ = j - 1;
                nj = _info->dimen.n + j;
                if( (g[j_] > alb[nj - 1] - _info->limits.epnewt*(1.0e0 + fabs( alb[nj - 1] ))) &&
                 (g[j_] < ub[nj - 1] + _info->limits.epnewt*(1.0e0 + fabs( ub[nj - 1] ))) )
                        goto L_530;
#ifdef IO_ENABLED
                if( _info->nintbk.ipr > 0 ) {
                   sprintf( LSGRG_MSGBUFFER, "%5ld  %13.6e  \n", j, g[j_] );
                   lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
                ;
L_530:
                ;
        }

        g[_info->nintbk.nobj - 1] = _info->bestbk.truobj;

        if( _info->glberr.abort || _info->lincnt.nfix == _info->dimen.n )
                goto L_540;

        /*      ----------------------------------------------------------------
         *     : Terminated in PHASE-1, reduce EPNEWT and solve the new problem
         *      ----------------------------------------------------------------
         *
         *@@@@ GO TO 540 */
        _info->limits.epnewt = _info->limits.epnewt*1.0e-01;
        if( _info->limits.epnewt < _info->epscom.eps2 )
                goto L_540;
#ifdef IO_ENABLED
        if( _info->nintbk.ipr > 0 )
                {
                sprintf( LSGRG_MSGBUFFER,
                 "\nRESTARTING WITH CONSTRAINT TOLERANCE (EPNEWT) SET"
                 " TO %12.5e\n",
                 _info->limits.epnewt );
                 lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
                }
#endif
        _info->cbmode.smode = TRUE;
        goto L_10;

        /*      ---------------------------------------------
         *     | We have tried everything, so terminate.
         *      --------------------------------------------
         * */
L_540:
        ;
        _info->limits.epnewt = epnorg;
        /*XX   IF (EPNEWT .NE. EPLAST) EPSTOP = (  0.01D0    )*EPSTOP */
        if( _info->limits.epnewt != _info->ingrg.eplast )
                _info->mngrg.epstop = epstpo;
        _info->limits.epnewt = _info->ingrg.eplast;
        _info->ph1bk.ph1eps = phhold;

        /*      ---------------------------------------------------------
         *     :  Unscale X and G, and bounds
         *     :  and adjust Reduced Gradient and Multipliers
         *      ---------------------------------------------------------
         * */
        if( _info->scal.iscale != 0 && _info->scal.scaled ){
                ainf = 0.1*_info->tols.plinfy;
                for( i = 1; i <= _info->dimen.n; i++ ){
                        i_ = i - 1;
                        x[i_] = x[i_]*ascale[i_];
                        if( alb[i_] > -ainf )
                                alb[i_] = alb[i_]*ascale[i_];
                        if( ub[i_] < ainf )
                                ub[i_] = ub[i_]*ascale[i_];
                        j = inbv[i_];
                        if( j <= _info->dimen.n )
                                gradf[i_] = gradf[i_]/(_info->sclobj.objscl*ascale[j - 1]);
                }
                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        i_ = i - 1;
                        j = _info->dimen.n + i;
                        g[i_] = g[i_]/ascale[j - 1];
                        if( alb[j - 1] > -ainf )
                                alb[j - 1] = alb[j - 1]/ascale[j - 1];
                        if( ub[j - 1] < ainf )
                                ub[j - 1] = ub[j - 1]/ascale[j - 1];
                        u[i_] = u[i_]*ascale[_info->dimen.n + i_]/_info->sclobj.objscl;
                }

                _info->scal.scaled = FALSE;
        }

        /*      ---------------------------------------------------------------
         *     :  If maximizing, flip sign of Reduced Gradient and Multipliers
         *      ---------------------------------------------------------------
         * */
        if( _info->optblk.maxim && _info->nintbk.ninf == 0 ){
                g[_info->nintbk.nobj - 1] = -g[_info->nintbk.nobj - 1];
                for( i = 1; i <= _info->dimen.n; i++ ){
                        i_ = i - 1;
                        gradf[i_] = -gradf[i_];
                }

                for( i = 1; i <= _info->dimen.mp1; i++ ){
                        i_ = i - 1;
                        u[i_] = -u[i_];
                }
        }

#ifdef IO_ENABLED
        if( _info->dbug.debug || _info->nintbk.ipr >= 5 )
                {
                lsgrg_msg(IOINFO,
                 "\n ***** EXITING  SUBROUTINE: GRGITN\n" );
                }
#endif

        /*     store final kt error to write to stats file
         *
         *     call stats_nlp('converge',i,KTERR,'xxx') */
        return _info->infbk.info;
        /*--------------------------------------------------------------------- */



        /*     END OF GRGITN
         * */
#undef  TABLST
#undef  CGSCR
#undef  SCG
#undef  YCG
} /* end of function */

void /*FUNCTION*/ LSGRGCLASS itrlog(LsgrgInfo *_info,long int n, long int mp1, long int nobj,
         long int nitr, long int nsuper, long int nb, long int ninf, long int *ipr,
         long int iprhld, long int *ipr3, long int iprhd3, long int *linect,
         LOGICAL32 degen, LOGICAL32 maxim, long int update, double cond,
         double condmx, double g[], double x[], double step, double kterr)
{
/*------------------------------------------------------------------------*/
/* ** 1/98 jcp ** notes on io protocol                                    */
/* . all io to main output (iounit.ioout) is ifdefed on  IO_ENABLED  */
/* . screen iteration log is ifdefed on IO_ENABLED so that a screen       */
/*   (or redirected ) iteration log will appear even if file-io is        */
/*   not enabled                                                          */
/*                                                                        */
/*                                                                        */
/*------------------------------------------------------------------------*/

//char dgflag;
//long int i, k;
//double  grnorm,  gtemp;
char t = 'T';
char blank = ' ';
char updtyp[3][4]={"  F"," QN"," CG"};
//char iterlog_header[200];


#ifdef IO_ENABLED
     sprintf(iterlog_header,
         "\nItn   Objective     Bind  Super Inf     Norm   Hessian"
         "  Basis    Up   Step  Deg"
         "\n No.   Function     Cons  Basic Cons   KT Err  Cond No"
         " Cond No  Date  Size  Itr\n" );

        if( _info->mngrg.iper != 0 )
                goto L_330;
/*  L_3600: unreferenced??? */
        if( *ipr < 1 )
                goto L_380;
        *linect = *linect + 1;
        /*XXX  IF (LINECT.LT.48.AND.IPR3.EQ.0) GO TO 340 */
        if( *linect < 48 )
                goto L_340;
        if( nitr == 0 )
                {
                lsgrg_msg(IOINFO, "\n\n\n\n\n" );
                }
        if( *ipr3 == 0 && nitr > 1 )
                {
                lsgrg_msg(IOINFO, "\n\n\n\n\n" );
                }
        lsgrg_msg(IOINFO,iterlog_header);
        *linect = 0;
        goto L_340;
L_330:
        ;
        k = nitr/_info->mngrg.iper*_info->mngrg.iper;
        if( k != nitr && k != nitr - 1 )
                goto L_340;
        if( nitr == 0 ){
                lsgrg_msg(IOINFO,iterlog_header);
                /*XXXX    IF (IOTERM .GT. 0) WRITE (IOTERM,770) */
        }
        if( nitr < 2 )
                goto L_340;
        if( k == nitr - 1 ){
                lsgrg_msg(IOINFO,iterlog_header);
                /*XXXX    IF (IOTERM .GT. 0) WRITE (IOTERM,770) */
        }
        *ipr = iprhld;
        *ipr3 = iprhd3;
L_340:
        ;
        /*340  GRNORM=0.0D0
         *XX   IF (NSUPER.EQ.0) GO TO 360
         *XX   DO 350 I=1,NSUPER
         *XX       IF (DABS(GRADF(I)).GT.GRNORM) GRNORM=DABS(GRADF(I)) */
/* L_350: unreferenced??? */
        ;
/*  L_360: unreferenced??? */
        ;
        /*?????      IF (NSUPER.GT.MAXH) COND=0.0D0 */
        gtemp = g[nobj - 1];
        if( maxim && ninf == 0 )
                gtemp = -gtemp;
        if( ninf == 0 )
                gtemp = gtemp/_info->sclobj.objscl;
        dgflag = blank;
        if( degen )
                dgflag = t;
        sprintf(LSGRG_MSGBUFFER,
  /*     " %3ld%14.6e%5ld%6ld%7ld%10.2e%9.2e%9.2e%3.3s%10.2e %c\n", */
         "%4ld %14.6e%5ld%5ld%5ld  %9.1e %8.1e %8.1e% 3.3s  %8.1e %c\n",
         nitr, gtemp, _info->bind.nbc, nsuper, ninf, kterr, cond, condmx,
          updtyp[update - 1]
         , step, dgflag );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);

        if( *ipr <= 3 )
                goto L_370;
        lsgrg_msg(IOINFO, " G IS\n" );
        for( i = 1; i <= mp1; i++ ){
                sprintf(LSGRG_MSGBUFFER, "%4ld%14.7e", i, g[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
        lsgrg_msg(IOINFO, "  \n" );
        lsgrg_msg(IOINFO, " X IS\n" );
        for( i = 1; i <= n; i++ ){
                sprintf(LSGRG_MSGBUFFER, "%4ld%14.7e", i, x[i - 1] );
                lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
        }
        lsgrg_msg(IOINFO, "  \n" );
L_370:
        ;
        /*xxx  IF (MAXIM.AND.NINF.EQ.0) G(NOBJ)=-G(NOBJ) */
        if( _info->mngrg.iper == 0 )
                goto L_380;
        if( k != nitr - 1 )
                goto L_380;
        if( nitr != 1 ){
                lsgrg_msg(IOINFO,iterlog_header);
                /*XX      IF (IOTERM .GT. 0) WRITE(IOTERM,770) */
        }
        *ipr = 1;
        *ipr3 = 0;
L_380:
        ;

#endif

#ifdef IO_ENABLED
        if( _info->io_info.screen_output_enabled  && *ipr > 0){
                /*xx      IF (MAXIM.AND.NINF.EQ.0) G(NOBJ)=-G(NOBJ) */
           if( (nitr%24) == 0 || nitr == 0) {
               sprintf( LSGRG_MSGBUFFER,iterlog_header);
               lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
           }
           sprintf( LSGRG_MSGBUFFER,
  /*         " %3ld%14.6e%5ld%6ld%7ld%10.2e%9.2e%9.2e%3.3s%10.2e %c\n", */
         "%4ld %14.6e%5ld%5ld%5ld  %9.1e %8.1e %8.1e% 3.3s  %8.1e %c\n",
             nitr, gtemp, _info->bind.nbc, nsuper, ninf, kterr, cond, condmx,
             updtyp[update - 1], step, dgflag );
                /*xx      IF (MAXIM.AND.NINF.EQ.0) G(NOBJ)=-G(NOBJ) */
          lsgrg_screen_msg(IOINFO,LSGRG_MSGBUFFER);
        }
#endif
        return;
        /*----------------------------------------------------------------------
         *
         *     END OF ITRLOG
         * */
} /* end of function */
double LSGRGCLASS powi(double x, long  n)
{
    double result=x;
    while(--n) result *= x;
    return result;
}
long  LSGRGCLASS ipow(long x, long  n)
{
    int result=x;
    while(--n) result *= x;
    return result;
}
void LSGRGCLASS printbasisInv(LsgrgInfo *_info,char *msg)
{
#ifdef IO_ENABLED
    int i;
    Integer2Array ind = (short*) _info->inv_ind;
    Integer2Array lnrc = (short*) _info->inv_lnrc;

    sprintf(LSGRG_MSGBUFFER,
      "\n ... Print of Basis Inverse (1st 50 elems): %s",msg);
    lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
    for(i=0;i<=49;++i)
     {
       sprintf(LSGRG_MSGBUFFER,
          "\n inv_binv[%3d] = %20.14e inv_ind[%3d] = %d",
            i,_info->inv_binv[i],i, ind[i]);
       lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
     }
    for(i=0;i< _info->mp1;++i)
    {
      sprintf(LSGRG_MSGBUFFER,
        "\n [%3d]  w = %14.7e w1 = %14.7e w2 = %14.7e",
         i,_info->inv_w[i],_info->inv_w1[i],_info->inv_w2[i]);
      lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
    }
    for(i=0;i< 9*_info->mp1;++i)
    {
       if(i < 2*_info->mp1)
          sprintf(LSGRG_MSGBUFFER,
             "\n [%3d]  inv_lnrc = %5d  inv_ip = %5d",
                i, lnrc[i], _info->inv_ip[i]);
       else
          sprintf(LSGRG_MSGBUFFER,
             "\n [%3d]  inv_lnrc = %5d",
              i, lnrc[i]);
       lsgrg_msg(IOINFO,LSGRG_MSGBUFFER);
    }

    lsgrg_msg(IOINFO,"\n..... end of basis inverse print ....");
#endif
    return;
}
#undef IOINFO
#undef LSGRG_MSGBUFFER
#undef JUMPBUF
#include "lsinfo.h"
/*
 ***********************************************************************
 ***             File:  grgmem.c                                     ***
 ***             Author: John C. Plummer                             ***
 ***********************************************************************
 ***********************************************************************
 ***                         LSGRG2C                                 ***
 ***           COPYRIGHT  1991, 1992, 1993, 1994, 1998               ***
 ***                      LEON S. LASDON                             ***
 ***                           &                                     ***
 ***                      ALLAN D. WAREN                             ***
 ***                           &                                     ***
 ***                      STUART SMITH                               ***
 ***                           &                                     ***
 ***                      John C. Plummer                            ***
 ***********************************************************************
 ***********************************************************************
*/
int  LSGRGCLASS lsgrg_memsetup(LsgrgInfo *info,
                    double xx[],double blvar[], double buvar[],
                    double blcon[],double bucon[], long lvars[],
                    long nnz_user)
{
long nnz, nnznl,nnz_total,nz,nznl,i;
int build_jacobian;
/*.......................................................................
 *
 *     *****PURPOSE (**3/98 jcp retain original fortran notes ):
 *     THIS SUBROUTINE CALCULATES MEMORY POINTERS FOR VARIABLE LENGTH
 *     (DATA DEPENDENT) ARRAYS. THESE INCLUDE THE ARRAYS USED TO STORE
 *     BOTH THE JACOBIAN AND THE BASIS INVERSE REPRESENTATION.
 *
 *     *****ARGUMENT DESCRIPTION:
 *     ON INPUT:
 *
 *     Z        - THE ARRAY OF AVAILABLE MEMORY TO BE PARTIONED.
 *     IBMAP    - BASIS INVERSE WORK ARRAY MAP
 *     NCORE    - THE TRUE DIMENSION OF HE Z ARRAY.
 *
 *     ON OUTPUT:
 *     POINTERS IN THE COMMON BLOCKS 'DYNAM5', 'DYNAM7', AND 'MEMORY'
 *     ARE CALCULATED IN ORDER TO PARTION MEMORY.
 *
 *     *****APPLICATION AND USAGE RESTRICTIONS:
 *     THIS SUBROUTINE IS CALLED ONCE TO PARTION MEMORY THAT DEPENDS ON
 *     THE NUMBER OF NONZEROES IN THE JACOBIAN MATRIX.  FIRST 'SETGR' IS
 *     CALLED TO CALCULATE THE NONZERO STRUCTURE AND POINTERS TO THE
 *     APPROPRIATE ARRAYS ARE SET.  FINALLY, THE REST OF AVAILABLE MEMORY
 *     IS ALLOCATED TO THE VARIABLE LENGTH ARRAYS BINV AND ICN.
 *
 *     *****ALGORITHM NOTES:
 *     POINTERS TO THE FOLLOWING ARRAYS ARE CALCULATED (THE ARRAYS IND,
 *     ICN AND BINV OVERLAP AND SHARE MEMORY.  ALSO, THE ARRAYS IKEEP AND
 *     IP SHARE THE SAME LOCATION.  THIS IS ALLOWABLE SINCE
 *     THE LA05 DATA STRUCTURES ARE DESTROYED BY THE MA28 ROUTINES) :
 *  ** 4/98 jcp ** arrays ind, icn, and ikeep are not referenced in the
 *                 current C code ??????
 *
 *                  LSINVERT DATA STRUCTURES
 *                  --------------------
 *
 *     LSINVERT       POINTER       TYPE               LENGTH
 *     --------       -------       ----               -------
 *     HEG             KHEG        INTEGER            N + MP1 + 1
 *     HAG             KHAG        INTEGER            LGRAD (=NZGRAD
 *     GRAD            KGRAD       DOUBLE PRECISION   LGRAD (=NZGRAD
 *     W1 (COND)       IBMAP(1)    DOUBLE PRECISION   MP1
 *                   inv_w1
 *     W2 (COND)       IBMAP(2)    DOUBLE PRECISION   MP1
 *                   inv_w2
 *     LNRC            IBMAP(3)    INTEGER*2          9*MP1
 *                   inv_lnrc
 *     IP              IBMAP(4)    INTEGER            2*MP1
 *                   inv_ip
 *     W               IBMAP(5)    DOUBLE PRECISION   MP1
 *                   inv_w
 *     IND             IBMAP(6)    INTEGER*2          2*LBINV
 *                   inv_ind
 *     A (BINV)        IBMAP(7)    DOUBLE PRECISION   LBINV
 *                   inv_binv
 *     *****REFERENCES:
 *
 *     *****HISTORY:
 *     WRITTEN BY:  STUART H. SMITH  DEPARTMENT OF MECHANICAL ENGINEERING
 *                  UNIVERSITY OF TEXAS, AUSTIN, TEXAS 78712.
 *     DATE LAST MODIFIED:  NOV   05 1990
 *       TRANSLATION TO C:  ALLAN D. WAREN
 *            LAST CHANGE:  MAY   11 1994
 *                          mar      1998 jcp
 *
 *    ** 3/12/98 jcp **
 *        removed partitioning of z and ibmap table for explicit
 *        allocation of all arrays. pointers are global. merge
 *        fixmem and varmem into grgsetupmem
 *
 *        Make initial call to setupj to determine nbr jacobian nonzeros
 *        then allocate larger of actual nnz and MIN_LGRAD
 *
 *       initial allocation for binv is:
 *          lbinv = n*mp1*BINVFACTOR;
 *          lbinv = lbinv < MINLBINV ? MINLBINV : lbinv;
 *       where:
 *         BINVFACTOR is rough estimate of inverse fillin density
 *         MINLBINV is minimum lbinv length allowable
 *
 *  ** 5/99 jcp **
 *     . after allocation, zero g array and arrays used for differencing
 *       some users with ignored rows might not place values into g
 *       thereby messing up the differencing procedure
 *
 *  ** 4/27/02 jcp **
 *     . if any allocation fails, call free_resources to free
 *       all allocated memory before exiting
 *
 *.......................................................................
 * */
    if( info->dbg.memsetup || info->ipr > 4 ) {
        lsgrg_msg(&(info->io_info),"\n .... Entry grgmemsetup .....");
    /*    lsgrg_echobnd(xx, blvar, buvar, blcon, bucon,  lvars,
                 "at entry memsetup"); */
     }

/*  compute values for global variables used in dimensions */

        info->n     = info->nvars;
        info->mp1   = info->nrows;
        info->m     = info->mp1 - 1;
        if( info->m == 0 )  info->m = 1;
        info->npmp1 = info->n + info->mp1;

/*     ALWAYS WORK WITH A BASIS OF SIZE MP1 */

        info->nbmax = info->mp1;
        info->maxb  = info->nbmax;
        if( info->maxr > info->n )
                info->maxr = info->n;
        info->nnbmax = info->n;
        if( info->nbmax > info->n )
                info->nnbmax = info->nbmax;
        info->mpnbmx = info->mp1;
        info->npnbmx = info->n + info->nbmax;
        info->nrtot  = info->maxr*(info->maxr + 1L)/2L;

/*--------------------------------------------------------------------*/
/*  for debug output, write out all dimension variables               */
/*--------------------------------------------------------------------*/
#ifdef IO_ENABLED
   if(info->dbg.memsetup || info->ipr > 2) {
    lsgrg_msg(&(info->io_info),"\n ..... Dimensions for Memory Allocation .....");
    sprintf(info->lsgrg_msgbuffer,
    "\n... kderiv = %d",info->kderiv);
    lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    lsgrg_lprint(info,"nvars",info->nvars);
    lsgrg_lprint(info,"n",    info->n    );
    lsgrg_lprint(info,"nrows",info->nrows);
    lsgrg_lprint(info,"mp1",  info->mp1  );
    lsgrg_lprint(info,"m",    info->m    );
    lsgrg_lprint(info,"nbmax",info->nbmax);
    lsgrg_lprint(info,"maxb", info->maxb );
    lsgrg_lprint(info,"maxr", info->maxr );
    lsgrg_lprint(info,"nnbmax",info->nnbmax);
    lsgrg_lprint(info,"mpnbmx",info->mpnbmx);
    lsgrg_lprint(info,"npnbmx",info->npnbmx);
    lsgrg_lprint(info,"nrtot",info->nrtot);
    lsgrg_msg(&(info->io_info),"\n ..... Beginning Memory Allocations ......");
  }
#endif
/*--------------------------------------------------------------------*/
/*  ALLOCATE SPACE FOR ARRAYS AND INITIALIZE ALL MEMORY BYTES TO ZERO */
/*--------------------------------------------------------------------*/
     if( (info->x       = allocate_darray(info,"x",info->npmp1))==NULL)
        return 0;
     if( (info->g       = allocate_darray(info,"g",info->mp1))==NULL)
        return 0;
     if( (info->alb     = allocate_darray(info,"alb",info->npmp1))==NULL)
        return 0;
     if( (info->ub      = allocate_darray(info,"ub",info->npmp1))==NULL)
        return 0;
     if( (info->ascale  = allocate_darray(info,"ascale",info->npmp1))==NULL)
        return 0;
     if( (info->r       = allocate_darray(info,"r",info->nrtot)) ==NULL)
        return 0;
     if( (info->gradf   = allocate_darray(info,"gradf",info->n))==NULL)
        return 0;
     if( (info->v       = allocate_darray(info,"v",info->nbmax))==NULL)
        return 0;
     if( (info->d       = allocate_darray(info,"d",info->n))==NULL)
        return 0;
     if( (info->u       = allocate_darray(info,"u",info->nbmax))==NULL)
        return 0;
     if( (info->gbest   = allocate_darray(info,"gbest",info->mp1))==NULL)
        return 0;
     if( (info->xbest   = allocate_darray(info,"xbest",info->mp1))==NULL)
        return 0;
     if( (info->xb1     = allocate_darray(info,"xb1",info->nbmax))==NULL)
        return 0;
     if( (info->xb2     = allocate_darray(info,"xb2",info->nbmax))==NULL)
        return 0;
      if( (info->xb3     = allocate_darray(info,"xb3",info->nbmax))==NULL)
        return 0;
     if( (info->dbnd    = allocate_darray(info,"dbnd",info->npnbmx))==NULL)
        return 0;
     if( (info->xstat   = allocate_darray(info,"xstat",info->n))==NULL)
        return 0;
     if( (info->gg      = allocate_darray(info,"gg",info->mp1))==NULL)
        return 0;
     if( (info->rr      = allocate_darray(info,"rr",info->nbmax))==NULL)
        return 0;
     if( (info->y       = allocate_darray(info,"y",info->n))==NULL)
        return 0;
     if( (info->gradfp  = allocate_darray(info,"gradfp",info->n))==NULL)
        return 0;
     if( (info->rowb    = allocate_darray(info,"rowb",info->nbmax))==NULL)
        return 0;
     if( (info->colb    = allocate_darray(info,"colb",info->nbmax))==NULL)
        return 0;
     if( (info->x0      = allocate_darray(info,"x0",info->n))==NULL)
        return 0;
     if( (info->g0      = allocate_darray(info,"g0",info->mp1))==NULL)
        return 0;
     if( (info->ycg     = allocate_darray(info,"ycg",info->maxcg*info->n))==NULL)
        return 0;
     if( (info->scg     = allocate_darray(info,"scg",info->maxcg*info->n))==NULL)
        return 0;
     if( (info->cgscr   = allocate_darray(info,"cgscr",info->maxcg*info->n))==NULL)
        return 0;
      if( (info->albo    = allocate_darray(info,"albo",info->npmp1))==NULL)
        return 0;
     if( (info->ubo     = allocate_darray(info,"ubo",info->npmp1))==NULL)
        return 0;
     if( (info->cdnum   = allocate_darray(info,"cdnum",info->nbmax))==NULL)
        return 0;

     if( (info->istat   = allocate_larray(info,"istat",info->mp1))==NULL)
        return 0;
     if( (info->inbv    = allocate_larray(info,"inbv",info->n))==NULL)
        return 0;
     if( (info->iub     = allocate_larray(info,"iub",info->n))==NULL)
        return 0;
     if( (info->ibc     = allocate_larray(info,"ibc",info->mp1))==NULL)
        return 0;
     if( (info->ibv     = allocate_larray(info,"ibv",info->mp1))==NULL)
        return 0;
     if( (info->ivstat  = allocate_larray(info,"ivstat",info->n))==NULL)
        return 0;
     if( (info->icols   = allocate_larray(info,"icols",info->npnbmx))==NULL)
        return 0;
     if( (info->icand   = allocate_larray(info,"icand",info->npnbmx))==NULL)
        return 0;
     if( (info->inlin   = allocate_larray(info,"inlin",info->nnlin))==NULL)
        return 0;
     if( (info->irank   = allocate_larray(info,"irank",info->nbmax))==NULL)
        return 0;
     if( (info->lksame  = allocate_larray(info,"lksame",info->npmp1))==NULL)
        return 0;
     if( (info->ibvbct  = allocate_larray(info,"ibvbct",info->n))==NULL)
        return 0;
     if( (info->tablst  = allocate_larray(info,"tablst",2L*info->maxtbu))==NULL)
        return 0;
     if( (info->iobjpt  = allocate_larray(info,"iobjpt",info->npmp1))==NULL)
        return 0;
/*------------------------------------------------------------------*/
/*    ---- Allocate space for Basis Inverse Arrays -----            */
/*           LNRC, IP, W, W1, AND W2 FIXED ARRAYS                   */
/*------------------------------------------------------------------*/
      if( (info->inv_w1 = allocate_darray(info,"inv_w1",info->mp1))==NULL)
        return 0;
      if( (info->inv_w2 = allocate_darray(info,"inv_w2",info->mp1))==NULL)
        return 0;
      if( (info->inv_lnrc = allocate_larray(info,"inv_lnrc",9*info->mp1))==NULL)
        return 0;
      if( (info->inv_ip = allocate_larray(info,"inv_ip",2*info->mp1))==NULL)
        return 0;
      if( (info->inv_w  = allocate_darray(info,"inv_w",info->mp1))==NULL)
        return 0;

      info->lbinv = (long) (info->n*info->mp1*BINVFACTOR);
      info->lbinv = (info->lbinv < MINLBINV) ? MINLBINV : info->lbinv;
/*  binv is used in anajac as a work array of length n */
/*  insure that it is at least that long               */
      info->lbinv = (info->lbinv < info->n) ? (info->n + 1) : info->lbinv;

      if( (info->inv_ind = allocate_larray(info,"inv_ind",2*info->lbinv))==NULL)
        return 0;
      if( (info->inv_binv = allocate_darray(info,"inv_binv",info->lbinv))==NULL)
        return 0;
      for(i=1;i<=info->lbinv;++i) info->inv_binv[i] = 0.0;
/*  zero out g and associated arrays */

      for(i=1; i<=info->mp1; ++i) {
          info->g[i]     = 0.0;
          info->gbest[i] = 0.0;
          info->gg[i]    = 0.0;
      }
/*-------------------------------------------------------------------*/
/*   . call lsgrg_initialize_model() to set bounds and status vectors  */
/*     to default values                                               */
/*                                                                     */
/*   .call lsgrg_setbounds to store user initial point into x and to   */
/*    and set bounds and status vectors                                */
/*    project point onto bounds prior to initial jacobian evaluation   */
/*    (if necessary)                                                   */
/*                                                                     */
/*-------------------------------------------------------------------*/
      lsgrg_initialize_model(info,info->n,info->mp1); /* set var/fcn defaults */

      if(!(lsgrg_setbounds(info,blvar,buvar,blcon,bucon,xx,lvars,info->ivstat)) )
             return info->lsgrg_return_status;
      if(info->dbg.memsetup || info->ipr > 1 )
          lsgrg_echo_lvars(info,info->n, info->ivstat, lvars,
                      info->nnlin, info->nlin);
/*-------------------------------------------------------------------*/
/*     if user is specifying jacobian and sparsity pattern,          */
/*     allocate memory for info arrays here based on nnz_user        */
/*     value in grgsub arg list                                      */
/*-------------------------------------------------------------------*/
    if( info->kderiv == USER_ANALYTIC ) {
       info->nzgrad = nnz_user;
      if( (info->ipcol = allocate_larray(info,"ipcol",nnz_user))==NULL)
        return 0;
      if( (info->iprow = allocate_larray(info,"iprow",nnz_user))==NULL)
        return 0;
      if( (info->ipmap = allocate_larray(info,"ipmap",nnz_user))==NULL)
        return 0;
      if( (info->paij   = allocate_darray(info,"paij",nnz_user))==NULL)
        return 0;
      nnz = nnz_user;
    }
/*---------------------------------------------------------------------*/
/*  if using finite differences,                                       */
/*  call lsgrg_setupj to count number of nonzeros at initial point     */
/*---------------------------------------------------------------------*/
    else {
       build_jacobian = 0;
       lsgrg_setupj(info,
               build_jacobian, info->xstat,info->g,
               info->gg, info->gbest, info->colb,
               info->lgrad, &nnz, &nnz_total,&nnznl, &info->colmax,
               info->ipcol, info->iprow, info->ipmap, info->paij );
    }
/*---------------------------------------------------------------------*/
/*  nonzero count established                                          */
/*  allocate to that length and call again to build jacobian data      */
/*  structures                                                         */
/*   We need at least nnz+mp1 locations in grad and hag to hold        */
/*   structural nonzeros and slacks.  need n+mp1+1 locations in        */
/*   heg and hegl to hold structurals and slacks                       */
/*                                                                     */
/*                                                                     */
/*---------------------------------------------------------------------*/
    info->lgrad  = (nnz+info->mp1) > MIN_LGRAD ? (nnz+info->mp1): MIN_LGRAD;
    info->lgrad += info->jacobian_growth_allowance;

#ifdef IO_ENABLED
    if(info->dbg.memsetup || info->ipr > 1 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n ... lsgrg_setupj has found %d nonzeros at initial point"
        "\n     nnz+mp1 = %d",nnz,(nnz+info->mp1) );
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
        sprintf(info->lsgrg_msgbuffer,
         "\n  MIN_LGRAD = %d info->jacobian_growth_allowance = %d "
         "\n  Final value of lgrad = %d",MIN_LGRAD,info->jacobian_growth_allowance,
         info->lgrad);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif
/*-------------------------------------------------------------------*/
/*  allocate grad and associated index arrays                        */
/*-------------------------------------------------------------------*/
   if((info->grad = allocate_darray(info,"grad",info->lgrad))==NULL)
        return 0;
   if((info->ihag  = allocate_larray(info,"ihag",info->lgrad))==NULL)
        return 0;
   if((info->iheg  = allocate_larray(info,"iheg",(long)(info->n+info->mp1+1)))==NULL)
        return 0;
   if((info->ihegl = allocate_larray(info,"ihegl",(long)(info->n+info->mp1+1)))==NULL)
        return 0;

/* **fixme** allocate iheg0 and ihegl0 to store 0-based */
/*           subscripts for v1 algorithm routines       */
/*  either ifdef these or just expunge for v2 alg routines */

    if((info->iheg0  = allocate_larray(info,"iheg0",(long)(info->n+info->mp1+1)))==NULL)
        return 0;
   if((info->ihegl0 = allocate_larray(info,"ihegl0",(long)(info->n+info->mp1+1)))==NULL)
        return 0;
/*---------------------------------------------------------------------*/
/*  Now call setupj to build jacobian structure                        */
/*  . nnz is structural nonzeros                                       */
/*  . nnznl is nonlinear nonzeros  ???                                 */
/*  . nnz_total is total nonzeros                                      */
/*  . lgrad is allocated length of grad                                */
/*                                                                     */
/*---------------------------------------------------------------------*/
   if(info->dbg.memsetup) lsgrg_msg(&(info->io_info),"\n ..... calling setupj to with build=1");

   build_jacobian = 1;
   if( !lsgrg_setupj(info,
           build_jacobian, info->xstat,info->g,
           info->gg, info->gbest, info->colb,
           info->lgrad, &nz, &nnz_total, &nznl, &info->colmax,
           info->ipcol, info->iprow, info->ipmap, info->paij) ) return 0;

    info->nzgrad = nz;
    info->nznlin = nznl;
    info->nzlin = info->nzgrad - info->nznlin;

#ifdef IO_ENABLED
    if(info->dbg.memsetup || info->ipr > 1 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n ---------- Final Jacobian Allocation ----------"
        "\n  Space reserved for               %7d elements"
        "\n  Structural Nonzeros    =         %7d"
        "\n  Slacks                 =         %7d"
        "\n  Total Nonzeros         =         %7d"
        "\n\n  Nonlinear Nonzeros     =         %7d"
        "\n  Space for new Nonzeros =         %7d  elements",
        info->lgrad,info->nzgrad,info->mp1,nnz_total,info->nznlin,
         info->jacobian_growth_allowance);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif
/*------------------------------------------------------------------*/
/*     grad and associated arrays now alocated and built            */
/*------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*   end of storage allocation code                                  */
/*-------------------------------------------------------------------*/
#ifdef IO_ENABLED
    if(info->dbg.memsetup || info->ipr > 4 ) {
        sprintf(info->lsgrg_msgbuffer,
             "\n grgmemsetup: lbinv = %d lgrad = %d"
             "\n ... Exiting grgmemsetup ....\n",
                 info->lbinv, info->lgrad );
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif

        return 1;

/*     end of grgmemsetup  */

} /* end of function */
/*================================================================*/
/*     memory allocation functions                                */
/*   dbl_alloc(), long_alloc(), char_alloc(), and charp_alloc()   */
/*   do actual allocation calls.  allocate_darray(),              */
/*   allocate_larray(), allocate_carray(), allocate_pcarray()     */
/*   are shells which trap failures, accumulate stats, and        */
/*   (for long and double) add pointers to table for table driven */
/*   free() calls                                                 */
/*                                                                */
/*   allocate_darray_nostore() and allocate_larray_nostore()      */
/*   allocate doubles and longs but do not accumulate stats or    */
/*   store pointers in table.  these are provided for allocation  */
/*   calls from test harness which will be freed separately       */
/*   and should not count towards algorithm memory stats          */
/*================================================================*/
double* LSGRGCLASS dbl_alloc(long length)
{
    /*  allocate 2 more than specified length to allow 1-n references */
    /*  plus one buffer word                                          */
    length +=2;
    return ( (double*) (malloc((size_t) length*sizeof(double))) );
}
long* LSGRGCLASS long_alloc(long length)
{
    /*  allocate 2 more than specified length to allow 1-n references */
    /*  plus one buffer word                                          */
    length +=2;
    return ( (long*) (malloc((size_t) length*sizeof(long))) );
}
char* LSGRGCLASS char_alloc(long length)
{
    /*  allocate 1 more char for trailing null */
    return ( (char*) (malloc((size_t) ++length*sizeof(char))) );
}
char** LSGRGCLASS charp_alloc(long length)
{
    /*  allocate 2 more than specified length to allow 1-n references */
    /*  plus one buffer word                                          */
    length +=2;
    return ( (char**) (malloc((size_t) 20+length*sizeof(char*))) );

}
/*--------------------------------------------------------------------*/
double* LSGRGCLASS allocate_darray(LsgrgInfo *info,char *name,long length)
{
     double *p;
/*  ** 4/01 jcp set table positions for grad,inv_binv for dynamic */
/*     resizing                                                   */

#ifdef IO_ENABLED
    if(info->dbg.memalloc || info->ipr > 2 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n Allocating Memory for Double array   [%s] length = %d",
        name,length);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif

    p = dbl_alloc(length+20);

    if(p==NULL)  {
       lsgrg_free_resources(info);
       lsgrg_insfmem(info,name, length);
    }
    else {     /* add pointer to table, unwind and exit on overflow */
       if(!lsgrg_alloctable(info,name, LSGRG_MEMALLOC, (void *) p) ) {
            lsgrg_free_resources(info);
            return NULL;
       }
       info->memstats.doublewords += length+1;
       info->memstats.doublebytes += (length+1)*sizeof(double);
       info->memstats.totalbytes  += info->memstats.doublebytes;
/* save table locations of binv and grad for reallocation */
       if(strcmp(name,"inv_binv")==0) info->tbl_loc_inv_binv =
               (info->alloctable_entries-1);
       if(strcmp(name,"grad")==0) info->tbl_loc_grad =
               (info->alloctable_entries-1);

    }
    return  p;
}

long* LSGRGCLASS allocate_larray(LsgrgInfo *info,char *name,long length)
{
    long *p;

/*  ** 4/01 jcp set table positions for ihag  for dynamic         */
/*     resizing                                                   */
#ifdef IO_ENABLED
    if(info->dbg.memalloc || info->ipr > 2 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n Allocating Memory for Long   array   [%s] length = %d",
        name,length);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif
    p = long_alloc(length+20);

    if(p==NULL) {
       lsgrg_free_resources(info);
       lsgrg_insfmem(info,name, length);
    }
    else {     /* add pointer to table, unwind and exit on overflow */
       if(!lsgrg_alloctable(info,name, LSGRG_MEMALLOC, (void *) p) ) {
           lsgrg_free_resources(info);
            return NULL;
       }
       info->memstats.longwords   += length+1;
       info->memstats.longbytes   += (length+1)*sizeof(long);
       info->memstats.totalbytes  += info->memstats.longbytes;
/* save table location  of ihag for reallocation */
       if(strcmp(name,"ihag")==0) info->tbl_loc_ihag =
            (info->alloctable_entries-1);

    }
    return p;
}

char* LSGRGCLASS allocate_carray(LsgrgInfo *info,char *name,long length)
{
/* allocate char array to store character string */
/*  dont add to pointer table because we will use this for */
/*  allocating row/var names which could be a large number */
/*  the free routines for names will handle                */

     char *p;

#ifdef IO_ENABLED
    if(info->dbg.memalloc || info->ipr > 2 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n Allocating Memory for Char   array   [%s] length = %d",
        name,length);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif
    /*  char_alloc() allocates 1 extra for trailing null    */
    p = char_alloc( length);

    if(p==NULL) {
       lsgrg_insfmem(info,name, length);
       lsgrg_free_resources(info);
    }
    else {
       info->memstats.charwords   += length+1;
       info->memstats.charbytes   += (length+1)*sizeof(char);
       info->memstats.totalbytes  += info->memstats.charbytes;
    }

    return p;
}
char** LSGRGCLASS allocate_pcarray(LsgrgInfo *info,char *name,long length)
{
/*  allocates pointer-to-char arrays to anchor string storage */
    char **p;

#ifdef IO_ENABLED
    if(info->dbg.memalloc || info->ipr > 2 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n Allocating Memory for Char*  array   [%s] length = %d",
        name,length);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif

    p = charp_alloc( length);

    if(p==NULL) {
       lsgrg_insfmem(info,name, length);
       lsgrg_free_resources(info);
    }
    else {
       info->memstats.cpointerwords   += length+1;
       info->memstats.cpointerbytes   += (length+1)*sizeof(char*);
       info->memstats.totalbytes      += info->memstats.cpointerbytes;
    }

    return p;
}
/*-------------------------------------------------------------------*/
/* nostats versions of allocate shells, no storage of pointers int   */
/* table and no statistic accumulation                               */
/* ... responsibility of caller to free this memory                  */
/*-------------------------------------------------------------------*/
double* LSGRGCLASS allocate_darray_nostats(LsgrgInfo *info,char *name,long length)
{
/* allocate double array but do not add to pointer table or */
/* accumulate stats                                         */

     double *p;

#ifdef IO_ENABLED
    if(info->dbg.memalloc || info->ipr > 2 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n Allocating Memory for Double array   [%s] length = %d nostats",
        name,length);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif


    if( (p = dbl_alloc( length) ) ==NULL)
         lsgrg_insfmem(info,name, length);

    return  p;
}  /* end of allocate_darray_nostats() */

long* LSGRGCLASS allocate_larray_nostats(LsgrgInfo *info,char *name,long length)
{
/* allocate long array but do not add to pointer table or */
/* accumulate stats                                         */

    long *p;

#ifdef IO_ENABLED
    if(info->dbg.memalloc || info->ipr > 2 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n Allocating Memory for Long   array   [%s] length = %d nostats",
        name,length);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif

    if( (p = long_alloc( length+20) ) ==NULL)
          lsgrg_insfmem(info,name, length);

    return p;
}

char* LSGRGCLASS allocate_carray_nostats(LsgrgInfo *info,char *name,long length)
{
/* allocate char array to store character string */

     char *p;

#ifdef IO_ENABLED
    if(info->dbg.memalloc || info->ipr > 2 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n Allocating Memory for Char   array   [%s] length = %d nostats",
        name,length);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif
    /*  char_alloc() allocates 1 extra for trailing null    */

    if( ( p = char_alloc( length) ) ==NULL)
          lsgrg_insfmem(info,name, length);
    return p;
}

char** LSGRGCLASS allocate_pcarray_nostats(LsgrgInfo *info,char *name,long length)
{
/*  allocates pointer-to-char arrays to anchor string storage */
    char **p;

#ifdef IO_ENABLED
    if(info->dbg.memalloc || info->ipr > 2 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n Allocating Memory for Char*  array   [%s] length = %d nostats",
        name,length);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif

    if( (p = charp_alloc( length) ) ==NULL)
       lsgrg_insfmem(info,name, length);
    return p;
}

/*-------------------------------------------------------------------*/
void LSGRGCLASS lsgrg_insfmem(LsgrgInfo *info,char *arrayname,long length)
{
/*  issue error msg (if enabled), set status */

#ifdef IO_ENABLED
  sprintf(info->lsgrg_msgbuffer,
   "\n Error: Unable to Allocate Memory for array [%s] length = %d",
     arrayname,length);
  lsgrg_error_msg(&(info->io_info),info->lsgrg_msgbuffer);
  if( strcmp(arrayname,"r")==0) {
     sprintf(info->lsgrg_msgbuffer,"\n"
     "\n r is the approximate Hessian for Quasi-Newton Methods"
     "\n and is allocated to (maxhes*(maxhes+1))/2 which can "
     "\n be very large for large problems.  Specify a small "
     "\n number for maxhes to switch to a conjuate gradient "
     "\n method.");
     lsgrg_error_msg(&(info->io_info),info->lsgrg_msgbuffer);
   }
#endif

   info->lsgrg_return_status = _LSGRG_INSFMEMORY;
   return;
}
void LSGRGCLASS lsgrg_memstats(LsgrgInfo *info)
{
#ifdef IO_ENABLED
   sprintf(info->lsgrg_msgbuffer,
     "\n LSGRGC Memory Allocation Statistics:\n"
     "    Data Type             Nbr Words           Nbr Bytes\n"
     "     Double                %10d               %10d\n"
     "     Long                  %10d               %10d\n"
     "     Char                  %10d               %10d\n"
     "     Pointer to Char       %10d               %10d\n"
     "     TOTAL                                          %10d\n",
     info->memstats.doublewords,  info->memstats.doublebytes,
     info->memstats.longwords,    info->memstats.longbytes,
     info->memstats.charwords,    info->memstats.charbytes,
     info->memstats.cpointerwords,info->memstats.cpointerbytes,
     info->memstats.totalbytes);
     lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
   return;
}
void LSGRGCLASS lsgrg_free_resources(LsgrgInfo *info)
{
/*  free memory acquired in lsgrg_memsetup                    */
/*  doubles and longs are stored in table in lsgrg_alloctbl() */

     void *pdummy=NULL;
     lsgrg_alloctable(info," ",LSGRG_MEMFREE, pdummy);

/*----------------------------------------------------------------*/
/*   free row and column name memory if allocated                 */
/*----------------------------------------------------------------*/
     if(info->lsgrg_rownames != NULL) lsgrg_free_rownames(info);
     if(info->lsgrg_varnames != NULL) lsgrg_free_varnames(info);

     lsgrg_zero_memstats(info);
}
void LSGRGCLASS lsgrg_zero_memstats(LsgrgInfo *info)
{
  info->memstats.doublewords   = 0; /* zero out memory allocation stats */
  info->memstats.doublebytes   = 0;
  info->memstats.longwords     = 0;
  info->memstats.longbytes     = 0;

/*  if no previous calls to setvarnames/setrownames or just freed    */
/*  zero char stats                                                  */

  if(info->lsgrg_rownames==NULL && info->lsgrg_varnames==NULL) {
       info->memstats.charbytes     = 0;
       info->memstats.charwords     = 0;
       info->memstats.cpointerwords = 0;
       info->memstats.cpointerbytes = 0;
    }
    info->memstats.totalbytes    = 0;
}

int LSGRGCLASS lsgrg_alloctable(LsgrgInfo *info,char *name, int cmd, void *p)
{
/*------------------------------------------------------------------*/
/*  alloctable stores pointers for allocated arrays in a table      */
/*  when cmd = LSGRG_MEMALLOC, then issues free() calls for all     */
/*  table entries when cmd = LSGRG_MEMFREE                          */
/*  add ptablechk to check for pointer corruption                   */
/*  remove ptablechk                                                */
/*------------------------------------------------------------------*/
int  i;

        if(cmd == LSGRG_MEMALLOC) {
           if(++(info->alloctable_entries) >= MAXALLOC) {  /* check for table overflow */
#ifdef IO_ENABLED
             sprintf(info->lsgrg_msgbuffer,
               "\n Internal Error: Size of array allocation table"
               " exceeds limit of %d.\n   Increase MAXALLOC in lsconfig.h\n",
                 MAXALLOC);
              lsgrg_error_msg(&(info->io_info),info->lsgrg_msgbuffer);
#endif
              info->lsgrg_return_status = _LSGRG_ALLOCTBL_OVERFLOW;
              return 0;
           }

#ifdef IO_ENABLED
    if(info->dbg.memalloc || info->ipr > 2 ) {
        sprintf(info->lsgrg_msgbuffer,
        "\n Pointer for Array [%s] is table entry %d",
        name,info->alloctable_entries);
        lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
    }
#endif
           info->ptable[info->alloctable_entries] = p;
        }
        if(cmd == LSGRG_MEMFREE) {
            for(i=0;i<=info->alloctable_entries;++i){
#ifdef IO_ENABLED
               if(info->dbg.memalloc ) {
                   sprintf(info->lsgrg_msgbuffer,
                   "\n Freeing Pointer number  %d",i);
                   lsgrg_msg(&(info->io_info),info->lsgrg_msgbuffer);
               }
#endif
             if(info->ptable[i] != NULL) free(info->ptable[i]);
            }
            info->alloctable_entries = -1;
        }  /* end code to free pointers in table */
         return 1;
}

int LSGRGCLASS lsgrg_alloc_varnames(LsgrgInfo *info,long nnvars)
{
/*   allocates memory to store variable names */
/* ** 4/27/02 jcp ** if any failure, free any memory allocated here */
/*   alloc_carray will free anything else   */
     long i;
     if( (info->lsgrg_varnames =
          allocate_pcarray(info,"lsgrg_varnames",nnvars)) == NULL)
            return 0;

     for(i=1;i<= nnvars;++i) info->lsgrg_varnames[i] = NULL;
     for(i=1;i<= nnvars;++i) {
         if( (info->lsgrg_varnames[i] =
            allocate_carray(info,"lsgrg_varnames[i]",NAMELENGTH)) == NULL){
              lsgrg_free_varnames(info);
              return 0;
         }
         strcpy(info->lsgrg_varnames[i]," ");
     }
     info->nvars_varnames = nnvars;
     return 1;
}
int LSGRGCLASS lsgrg_alloc_rownames(LsgrgInfo *info,long nnrows)
{
/*   allocates memory to store row names */
/* ** 4/27/02 jcp ** if any failure, free any memory allocated here */
/*   alloc_carray will free anything else   */
     long i;
     if( (info->lsgrg_rownames =
          allocate_pcarray(info,"lsgrg_rownames",nnrows)) == NULL)
            return 0;

     for(i=1;i<= nnrows;++i) info->lsgrg_rownames[i] = NULL;
     for(i=1;i<= nnrows;++i) {
         if( (info->lsgrg_rownames[i] =
            allocate_carray(info,"lsgrg_rownames[i]",NAMELENGTH)) == NULL){
              lsgrg_free_rownames(info);
              return 0;
         }
         strcpy(info->lsgrg_rownames[i]," ");
     }
     info->nrows_rownames = nnrows;
     return 1;
}
int lsgrg_free_rownames(LsgrgInfo *info)
{
      long i;
      for(i=1; i<= info->nrows_rownames; ++i)  /* free name string storage */
           if(info->lsgrg_rownames[i] !=NULL)
                free(info->lsgrg_rownames[i]);
      if(info->lsgrg_rownames != NULL)
          free(info->lsgrg_rownames);          /* free pointer array storage */
      info->lsgrg_rownames = NULL;
      info->usernames = 0;
      info->nrows_rownames = 0;
      return 1;
}
int lsgrg_free_varnames(LsgrgInfo *info)
{
     long i;
     for(i=1; i<= info->nvars_varnames; ++i)  /* free name string storage */
          if(info->lsgrg_varnames != NULL) free(info->lsgrg_varnames[i]);
     if(info->lsgrg_varnames != NULL)
        free(info->lsgrg_varnames);           /* free pointer array storage */
     info->lsgrg_varnames = NULL;
     info->usernames = 0;
     info->nvars_varnames = 0;
     return 1;
}

