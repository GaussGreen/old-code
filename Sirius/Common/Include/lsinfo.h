/***********************************************************************/
/*   File:         lsinfo.h                                            */
/*   Author:       John C. Plummer                                     */
/*   Last Update:  20 Mar 2002                                         */
/*                                                                     */
/*   This header file contains 'global' structure definitions          */
/*   prototype declarations for lsgrgc v3.0 ANSI C version             */
/***********************************************************************/

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
/* guard against multiple inclusion of this file */

#ifndef LSINFO_INCLUDED
  #define LSINFO_INCLUDED

#include <stdlib.h>

#include <string.h>
#include <math.h>
#include <setjmp.h>
#include <ctype.h>

#ifdef max
  #undef max
#endif

#ifdef min
  #undef min
#endif


#define RELEASEDATE "11012001"
/*--------------------------------------------------------------------*/
/* -----------------------  Change Log  ------------------------------*/
/*  ** 2/01 jcp ** add 'retry_history' array to LsgrgInfo             */
/*  ** 5/01 jcp ** add prototype for lsgrg_msg2 -- screen and file    */
/*  ** 3/02 jcp ** remove optquest/oqgrg stuff except for file ptrs   */
/*    in ioinfo.  everything will be defined in oqgrg.h/oqms.h        */
/*  . leave typedefs for oqgrg info structures and gams info ptr      */
/*  ** 4/02 jcp ** add 'maxtime' to lsgrginfo                         */
/*    add prototype for lsgrg_set_maxtime(LsgrgInfo *lsinfo,double t) */
/*                                                                    */
/*  ** 5/01/02 jcp ** to clear warnings and strange external refs     */
/*  if max and min are defined by stdlib, undef them and define       */
/*  copies here.                                                      */
/*                                                                    */
/*                                                                    */
/*--------------------------------------------------------------------*/


/* include all other lsgrg2c header files here */

#include "lsconfig.h"   /*  code configurations and macro definitions */

//
// A.G. Because we redefine the typedef's, do not include this until after
//		we define LSGRGSolver class
//
//#include "lstypes.h"    /*  typedefs for user gcomp/parsh functions   */
#include "solve.h"
#include "lscodes.h"    /*  enum for lsgrg2c termination codes        */

/*--------------------------------------------------------------------*/
/* inclusions conditional on macros in lsconfig.h                     */
/*--------------------------------------------------------------------*/
/*  LSGRGCLASS is macro class name in  C++ version of code; set to    */
/*     null string  for ANSI C version.  if LSGRGCPP is defined,      */
/*     define                                                         */
/*--------------------------------------------------------------------*/

#ifdef LSGRGCPP
       extern Lsgrg LSGRG;
       #define LSGRGOBJECT LSGRG.
       #define LSGRGCLASS Lsgrg::
#else
       #define LSGRGOBJECT
       #define LSGRGCLASS
#endif

#ifdef IO_ENABLED   /* do not include stdio.h unless io is enabled */
  #include <stdio.h>
#endif

#ifdef TIMING_ENABLED  /* do not include time.h unless timing is enabled */
   #include <time.h>
#endif
/*-------------------------------------------------------------------------*/
/*  define the typedef for the pointer to Optquest_grg_info_t              */
/*    structure which is included in the LsgrgInfo structure def           */
/*  define the typedef for the pointer to linear_row_info_t               */
/*    which is used for oqgrg gams interface                               */
/* NOTE:  these typedefs are hardcoded to refer to structure tags          */
/*        in oqgrg.h                                                       */
/*-------------------------------------------------------------------------*/
#ifdef LSGRG_OPTQUEST
        typedef struct Optquest_grg_info_t *OqGrgInfoP;
        typedef struct linear_row_info_t *linrowinfoP;
        typedef struct OqMsInfo_t        *OqMsInfoP;
#endif
/*-------------------------------------------------------------------------*/
/* io_info_t contains io stuff which needs to be threaded down to          */
/*  io routines char*lsgrg_msgbuffer is a duplicate pointer to the         */
/*  lsgrg_msgbuffer in the shell                                           */
/*  ** 12/99 jcp ** updated to contain io units for optquest inteface      */
/*  ** 5/01 jcp ** add logicals for output destinations to io_info_t  */
/*-------------------------------------------------------------------------*/
typedef   struct io_info_t {   /* structure for io unit stuff */
          int error_output_enabled;  /* logicals for enabling error msgs */
          int screen_output_enabled;
          int ioout_enabled;
          int ioout_set_by_user;
          int ioerr_set_by_user;
          int ioterm_set_by_user;
          int ioerr_to_ioout;
          int ioerr_to_ioterm;
          char * lsgrg_msgbuffer;
#ifdef FILE_IO_ENABLED
          FILE * lsgrg_iodump; /* file pointers if file io enabled */
          FILE * lsgrg_ioerr;
          FILE * lsgrg_ioin;
          FILE * lsgrg_ioout;
          FILE * lsgrg_ioterm;
#ifdef LSGRG_OPTQUEST
          FILE * oqgrg_testlog;  /* file pointers for oqgrg */
          FILE * oqgrg_localsfile;
          FILE * oqgrg_probout;
          FILE * oqgrg_statistics;
          FILE * oqgrg_trialpts;
          FILE * oqgrg_poplist;
#endif

#endif
} IoInfo;
/*---------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*   misc definitions and typedef's                                        */
/*-------------------------------------------------------------------------*/
typedef double*      DoubleArray;
typedef long*        IntegerArray;
typedef short int*   Integer2Array;

typedef long LOGICAL32;
typedef char* STRING;

#ifdef FILE_IO_ENABLED
   typedef FILE * IOUNIT;
#endif
/*-------------------------------------------------------------------------*/
/*    ---- misc macro definitions ----                                     */
/*-------------------------------------------------------------------------*/
#ifndef TRUE
        #define TRUE 1
        #define FALSE 0
#endif

/* ---- needed by algorithm routines -----*/

#define TorF(l) ( l ? 'T' : 'F' )
#define         fmax(a,b)       (double)( (a) > (b) ? (a) : (b) )

/* 05/01/02 jcp comment out these definitions to quiet compiler */
/*  warnings about redefined macros. these duplicate stdlib */

#define          max(a,b)               ( (a) > (b) ? (a) : (b) )
#define          min(a,b)               ( (a) < (b) ? (a) : (b) )

#define         fmin(a,b)       (double)( (a) < (b) ? (a) : (b) )
#define sign(a,b)       (double)( (b) < 0 ? -fabs(a) : fabs(a) )

/*-------------------------------------------------------------------------*/
/*  struct LsgrgInfo contains entries for all user-specifiable             */
/*   options as well as everything needed by the shell                     */
/* . the user defines a pointer to an LsgrgInfo struct and passes it       */
/*    to lsgrg_initialize() which allocates the memory and fills in the    */
/*    initial values                                                       */
/*-------------------------------------------------------------------------*/
struct LsgrgGlobals_t {

/* 4/27/02 jcp ** add maxtime, time limit for lsgrg run */
    double maxtime;

/* 3/20/02 jcp ** add print level for outres to allow printing */
/*  of selected sections from oqgrg */
    long outres_print_level;

/* 7/01 jcp ** logical for pack/unpack return values of multipliers */
    int return_mults_unpacked;

    jmp_buf *lsexit_jmpbuf;  /* pointer to setjmp/longjmp buffer */

	//
	// A.G. 
	//
	LSGRGSolver* m_pOwner;

    P_GCOMP  lsgrg_user_gcomp;  /* pointers to evaluation functions */
    P_PARSH  lsgrg_user_parsh;
    P_GCOMPX lsgrg_user_gcomp_fcns;

    IoInfo io_info;               /* io flags and units */

/*  global status code, msg buffer */

       long lsgrg_return_status;
       char lsgrg_msgbuffer[MSGBUFFER_LENGTH];

/* ** new 2/01 jcp ** retry history vector */
/*  and flag for whether last point is feasible */
       int retry_history[LEN_RETRY_HISTORY];
       int last_point_feasible;
/* ** 7/01 jcp ** flag for disabling retries */
       int disable_retries;
/* ** 7/01 jcp ** flag for disabling copyback of info structures */
       int disable_copyback;
/* ** end new */

/*-------------------------------------------------------------------------*/
/*  for lsgrgc/optquest interface, pointer to OqGrgInfo structure */
/*  ifdef'ed here because typedef is not seen if lsgrgoq.h is not */
/*  included                                                      */
/*-------------------------------------------------------------------------*/

#ifdef LSGRG_OPTQUEST
    OqGrgInfoP OqGrg;
    OqMsInfoP  OqMsInfo;  /* 2nd oqms info structure */
/* for fixing extra objvar in gams interface;  */
/* moved here from oqgrginfo on 4/7/01, ZU     */
       int gams_iobvar;

#ifdef LSGRG_GAMS
/*------------------------------------------------------------------*/
/* pointer to linear row information used by oqgrg when run with    */
/*   GAMS interface                                                 */
/*------------------------------------------------------------------*/
       linrowinfoP linrowinfo;

/*  end ifdef LSGRG_GAMS */
#endif
/*  end ifdef LSGRG_OPTQUEST */
#endif

/*  ** 4/01 jcp ** add variables to hold reallocation factors */
/*   for binv and jacobian growth allowance                   */

    double binv_realloc_factor,jacobian_growth_realloc_factor;

    int lsgrg_setup, lsgrg_run, lsgrg_shutdown, lsgrg_setup_done;
    int usernames, nbr_user_options_set;
    long galloc, getsiz;
    char ** lsgrg_varnames, ** lsgrg_rownames;
    char ** lsgrg_varnames_user, ** lsgrg_rownames_user;

    long nvars_varnames,nrows_rownames;  /* dimensions from setrowname */

      double    cndtol,  edf, epboun, epnewt, epspiv, rtf,
                epinit, eplast, epdeg, eps, tolx, tolz,
                eps0, eps1, eps2, eps3, eps4, eps5, epstop,
                funpr, pstep, phmult, ph1eps, xajtol, xpvpct,
                pivtol[8], dluprm[5],truobj, plinfy, plzero;

      long int  inprnt, otprnt, iper, ipn4, ipn5, ipn6, limser, nstop,
                iquad, colmax, ibvblm, idglim, ipr, ipr3 , ipvtpt,
                iscale, isclag, iter, itlim, kderiv, lbinv, lgrad,
                lirn, lmem, npmp2, m, mp1, mpnbmx, n, nbmax, nnbmax,
                npmp1, npnbmx, nrtot, maxb, maxcg, maxr, maxtbu,
                nrows, nvars, maxh,mcgm1, memcg, modcg, nfix, nlin,
                nnlin, nnlequ, nzbinv, nzgrad, nzlin, nznlin,
                nobj, ngcomp, nparsh, nbc, ninf, maxtab, nnewton,
                nnewton_itns, nnewton_fail, nnewton_stepc,
                ncbs, nser, nrser, nreinv, nicond, nsing, ifrcnt,
                nsame, nsuper, nb, nsear, nph1ls, ndeg, nbs,
                jacobian_growth_allowance;

      int       fixpiv, hrdbnd, hscale, maxim, multsb, abort;

      char title[TITLELENGTH];
      char problemName[TITLELENGTH]; /* problem name for gams/oqgrg */
      LOGICAL32 useph0, gfeas;
/* 9/99 jcp new options for reduced step size */
      LOGICAL32 steepest_descent;  /* force steepest descent direction */
      double  steplimit;        /* restrict step to step_limit fraction */
                                /* of step determined by linesearch     */
/*----------------------------------------------------------------------*/
/* Pointers for Array Allocations                                       */
/*  arrays prefixed by inv_ are for invert subsystem                    */
/*----------------------------------------------------------------------*/
      DoubleArray   x,g,alb,ub,ascale,r,gradf,v,d,u;
      DoubleArray   gbest,xbest,xb1,xb2,xb3,dbnd;
      DoubleArray   xstat,gg,rr,y,gradfp,rowb,colb,x0,g0;
      DoubleArray   ycg,scg,cgscr,albo,ubo,cdnum, grad;
      DoubleArray   inv_w1, inv_w2, inv_w, inv_binv;
      DoubleArray   paij;    /* user parsh coefficients */

      IntegerArray  istat,inbv,iub,ibc,ibv,ivstat;
      IntegerArray  icols,icand,inlin,irank,lksame;
      IntegerArray  ibvbct,tablst ;
      IntegerArray  ihag, iheg, ihegl;
      IntegerArray  inv_lnrc, inv_ip, inv_ind;
      IntegerArray  iprow, ipcol, ipmap;  /* user parsh jacobian map */
      IntegerArray  iobjpt;

  /* ** 7/98 jcp ** pointers for 0-based algorithm routines */
             IntegerArray  iheg0, ihegl0;

     struct lsgrgDebug {   /* structure for routine-level debug flags */
       int grgsub, setvarname,setrowname,memsetup, setupj,
           memalloc, parshf, parshf0, parshc, parshc0, parsh0,
           globals,algorithm;
      } dbg;

/* ** 4/98 jcp ** structure for memory allocation stats */
     struct _lsrgmemorystats_t {
        long doublewords,doublebytes,longwords,longbytes,
             charbytes,
             charwords,cpointerwords,cpointerbytes, totalbytes;

        } memstats ;

/*  table for storage of memory allocation pointers */
    void *ptable[MAXALLOC],*ptablechk[MAXALLOC];
    int  alloctable_entries,tbl_loc_inv_binv,tbl_loc_grad,
           tbl_loc_ihag;
/* **3/98 jcp ** structures for user-settable parameters */

       struct lsgrg_option_name_tbl {
           char name[20];
        } lsgrg_option_name[MAXOPTIONS];
/*-------------------------------------------------------------------------*/
/*   accumulators for timing stats                                         */
/*-------------------------------------------------------------------------*/
    struct timestats_t {
           double gcomp, consbs, parsh, newton,newton_nogcomp,
              search, invert, phase0, phase1, total;
    } tim;
/*-------------------------------------------------------------------------*/
/*   statistics accumulators for test harness                              */
/*-------------------------------------------------------------------------*/
     struct grgtest_global {
           double obj,sinf,time;
           long ncols, nrows,
                itns, ninf, ngcomps;
           int info;
           char s_problem[20];
           char s_termination[150];
           char s_datetime[50];
           char s_termin[20],s_alg[20];
    } grgtest;
                             /* pointer to grgsub arg structure  */
                             /* so that test harness info can be */
                             /* passed back                      */
    struct grgtest_global* p_grgtest;

/* flags to replace static vars in alg routines */

   struct firstcall_flags {
          long newton,redobj;
          LOGICAL32 ph0log;
   }  firstcall;
   struct CgStatic {
      long initcg,itncg;
   }  cgStatic;
   struct DirecStatic { /* 1/3/00 jcp add msgcg */
          long nsupvm,msgvm,msgcg;
          double rtinsb;
   }  direcStatic;
   struct Ph1objStatic {
         double sinf0, true0;
   }  ph1objStatic;
   struct Ph0logStatic {
      int termheader;
      long iprhd3, iprhld;
   } ph0logStatic;
/*---------------------------------------------------------------------*/
/*   structures for previously global vars in algorithm                */
/*---------------------------------------------------------------------*/
       struct t_epscom {
               double eps0, eps1, eps2, eps3, eps4, eps5;
               }       epscom;
       struct t_cgbk {
               long int modcg, memcg, mcgm1, icgcnt, icgptr;
               LOGICAL32 hscale;
               }       cgbk;
       struct t_bestbk {
               double stpbst, objbst, step, stepmx, truobj;
               }       bestbk;
       struct t_redph {
               double trubst;
               }       redph;
       struct t_congrg {
               long int nsear0, lvlast;
               }       congrg;
       struct t_counts {
               long int nftn, ngrad, nminv, nnfail, ncalls, nit, nbs, nstepc,
                ndub;
               }       counts;
       struct t_optblk {
               LOGICAL32 maxim, hrdbnd, subset, multsb, gfeas, useph0, warmst,
                penlty;
               }       optblk;
       struct t_dfblk {
               LOGICAL32 dfail;
               }       dfblk;
       struct t_dimen {
               long int m, n, mp1, npmp1, nbmax, nnbmax, npnbmx, mpnbmx, nrtot;
               }       dimen;
       struct t_dirgrg {
               double cond;
               long int update, nsupp;
               }       dirgrg;
       struct t_zblck {
               double condmx;
               long int nblck;
               }       zblck;
       struct t_infbk {
               long int info;
               }       infbk;
       struct t_pardat {
               long int kderiv;
               }       pardat;
       struct t_ingrg {
               double epinit, eplast, epdeg;
               }       ingrg;
       struct t_inout {
               char title[77];
               }       inout;
#ifdef FILE_IO_ENABLED
       struct t_iounit {
               IOUNIT   ioin, ioout, iodump, ioerr, ioterm;
               }       iounit;
#endif
       struct t_limits {
               double epboun, epnewt, epspiv;
               long int itlim;
               }       limits;
       struct t_logblk {
               LOGICAL32 move, restrt, drop, varmet, conjgr, resetp;
               }       logblk;
       struct t_slpobj {
               double slope;
               }       slpobj;
       struct t_sernew {
               double edf;
               LOGICAL32 trunew;
               }       sernew;
       struct t_mfact {
               double edfper, stpper, rtnmul;
               }       mfact;
       struct t_misc {
               long int maxh, nsear, jp, lv, jqq;
               }       misc;
       struct t_mngrg {
               double epstop;
               long int limser, nstop, ierr, ipn4, ipn5, ipn6, iper;
               }       mngrg;
       struct t_iters {
               long int nitr, ndeg, nph0it, nph1ls;
               }       iters;
       struct t_degn {
               long int idglim;
               }       degn;
       struct t_nintbk {
               long int nb, nobj, ninf, nsuper, ipr3, ncand, ipr;
               }       nintbk;
       struct t_bind {
               long int nbc, nnbc;
               }       bind;
       struct t_ph1bk {
               double phmult, ph1eps;
               long int initph;
               }       ph1bk;
       struct t_srchlg {
               LOGICAL32 uncon, fail, jstfes, mxstep, unbd, succes, unconp;
               }       srchlg;
       struct t_supblk {
               LOGICAL32 sbchng, basbnd;
               }       supblk;

       struct t_cbmode {
               LOGICAL32 smode, prsmod, newpt, phase0, compgr, pertub;
               }       cbmode;
       struct t_tols {
               double eps, plinfy, plzero, tolx, tolz;
               }       tols;
       struct t_cbscnt {
               long int ncbs, nser, nreinv, nsame, nrser, nicond, nsing, ifrcnt;
               }       cbscnt;
       struct t_intggr {
               long int ngetgr;
               }       intggr;
       struct t_chq {
               long int maxtab, ltab, lptr;
               }       chq;
       struct t_memory {
               long int lgrad, lbinv, lirn, lbmap, lmem, np1;
               }       memory;
       struct t_lincnt {
               long int nlin, nnlin, nfix;
               }       lincnt;
       struct t_nwtcnt {
               long int inwtk, inwtfk, inwtpt, inwtlk[20];
               }       nwtcnt;
       struct t_degcnt {
               long int idegk, idegfk, idegpt, ideglk[20];
               LOGICAL32 redpvt;
               }       degcnt;
       struct t_sclobj {
               double objscl;
               }       sclobj;
       struct t_scal {
               long int iscale, isclag;
               LOGICAL32 scaled, newscl;
               }       scal;
       struct t_setin {
               long int lastz;
               LOGICAL32 galloc, getsiz;
               }       setin;
       struct t_glberr {
               LOGICAL32 abort;
               }       glberr;
       struct t_dbug {
               LOGICAL32 debug;
               }       dbug;
       struct t_gctime {
               double tgcomp;
               }       gctime;
       struct t_ph0tim {
               double tph0, tph1;
               }       ph0tim;
       struct t_pvtblk {
               double xpvpct;
               long int ibvblm;
               }       pvtblk;
       struct t_lsinvt {
               double dluprm[4];
               long int iluprm[24];
               }       lsinvt;
       struct t_rednew {
               double corner, xb, xstep;
               }       rednew;
       struct t_newsrc {
               long int iter;
               }       newsrc;
       struct t_gmode {
               LOGICAL32 hbnd2;
               long int jdbg;
               }       gmode;
       struct t_nwtim {
               double tnewt, tnewtx;
               }       nwtim;
       struct t_gmsblk {
               LOGICAL32 gmserr;
               }       gmsblk;
       struct t_quadbk {
              double a1, a2, a3;
              long int icon, iquad;
              }       quadbk;
       struct t_redser {
              long int ninfb;
              }       redser;

       struct t_hesblk {
              double gamma0, dirsc;
              }       hesblk;
       struct t_initbk {
              long int init, lastcl;
              }       initbk;
       struct t_nph0 {
              long int ninf0, ndrop;
              }       nph0;
       struct t_equblk {
              long int nnlequ;
              LOGICAL32 lincon;
              }       equblk;
       struct t_ztol {
              double xajtol, bigelt, biglin;
              }       ztol;
       struct t_temp {
              double tsear, tfact;
              }       temp;
       struct t_contim {
              double cbtime;
              }       contim;
       struct t_itnggr {
              long int ngetgr;
              }       itnggr;
       struct t_nzerog {
              long int nzgrad, nznlin, nzlin;
              }       nzerog;
       struct t_jgrow {
              long int maxgro;
              }       jgrow;
       struct t_cmax {
              long int colmax;
              }       cmax;
       struct t_gradtm {
              double tgrad;
              }       gradtm;
       struct t_stepbk {
              double pstep, funpr;
              }       stepbk;
       struct t_zcond {
              double cndtol;
              }       zcond;
       struct t_pivots {
              double pivtol[7];
              long int ipvtpt;
              LOGICAL32 fixpiv;
              }       pivots;
       struct t_nzerob {
              long int nzbas;
              }       nzerob;
       struct t_drvchk {
              long int nfddif;
              LOGICAL32 chkdrv;
              }       drvchk;
};
typedef struct LsgrgGlobals_t LsgrgInfo;
/*--------------------------------------------------------------------*/
/*   if LSGRG_TESTHARNESS is defined, include that header file        */
/*   this is placed below typedef for LsgrgInfo because a prototype   */
/*   inside references that type                                      */
/*                                                                    */
/*  if we are running the test harness with the optquest interface    */
/*  include the optquest test stuff                                   */
/*                                                                    */
/*--------------------------------------------------------------------*/
#ifdef LSGRG_TESTHARNESS
     #include "grgtest.h"

#ifdef LSGRG_OPTQUEST
       #include "lsgrgoqt.h"
#endif

#endif
/*-------------------------------------------------------------------*/
/*  prototypes for functions in the gams interface which need to     */
/*  be called                                                        */
/*-------------------------------------------------------------------*/
#ifdef LSGRG_GAMS
   int lsgamsGetUserAbortStatus();
#endif
/*====================================================================*/
/*  prototypes for shell functions                                    */
/*====================================================================*/
/*-------------------------------------------------------------------*/
/*  lsgrgc shell  function prototypes                                */
/*-------------------------------------------------------------------*/
/*  per function gcomp not connected this version */
void lsgrg_gcompx(LsgrgInfo *,double[],double[],long*,long[]);

void parsh_gams(double[],double[],long*,long*,double[]);

/*   ----- function and derivitive functions parsh.c  ------- */

void lsgrg_gcomp1(LsgrgInfo *,double g[], double x[]);
void lsgrg_gcomp0(LsgrgInfo *,double g[], double x[]);

/* --- user parsh --- */
void parsh_lsgrg(LsgrgInfo *,double* , long, double* , long* , long* , long*);

void lsgrg_parsh0(LsgrgInfo *,double x[],long n,double grad[],
                double paij[],long iprow[], long ipcol[],
                long ipmap[],long nzgrad);

/* --- finite difference parsh routines */

void lsgrg_parshf1(LsgrgInfo *, double x[], double g[], double gcol[],
                  long j, double gplus[], double xub );
void lsgrg_parshf0(LsgrgInfo *, double x[], double g[], double gcol[],
                        long j, double gplus[], double xub );

void lsgrg_parshc1(LsgrgInfo *, double x[], double gcol[],long  j,
                   double gplus[], double gminus[] );
void lsgrg_parshc0(LsgrgInfo *, double x[], double gcol[],long  j,
                        double gplus[], double gminus[] );

void gcomp_check(LsgrgInfo *,char *msg,long n, double x[], long mp1,
                 double g[], int onebased);
/*
void parshc(double[],double[],double[],long,double[],double[]);
void parshf(double[],double[],double[],long,double[],double);
*/

  /*  prototype for optional user report function */

void lsgrg_report( double g[],double x[], long mp1, long n,
           int usernames,
           char **lsgrg_rownames,char **lsgrg_varnames, double x0[]);

                  /* ---- io functions grgio.c  -----*/

void lsgrg_msg(IoInfo *,char *buffer);
/* new function for concurrent screen and file output */
void lsgrg_msg2(IoInfo *, char *buffer);
void lsgrg_error_msg(IoInfo *,char *buffer);
void lsgrg_screen_msg(IoInfo *,char *buffer);
void lsgrg_output(IoInfo *,char *buffer,int destination);
void lsgrg_lprint(LsgrgInfo *,char *,long);
void lsgrg_dprint(LsgrgInfo *,char *,double);
void lsgrg_iprint(LsgrgInfo *,char *,int);
void lsgrg_disable_ioerr(LsgrgInfo *);
void lsgrg_enable_ioerr(LsgrgInfo *);
void lsgrg_disable_ioterm(LsgrgInfo *info);
void lsgrg_enable_ioterm(LsgrgInfo *info);
void lsgrg_quiet_all_output(LsgrgInfo *info);

#ifdef FILE_IO_ENABLED
   void lsgrg_set_ioout(LsgrgInfo *,FILE *file); /* only active if file io active */
   void lsgrg_set_ioerr(LsgrgInfo *,FILE *file);
   void lsgrg_set_ioterm(LsgrgInfo *,FILE *file);
#endif


                  /* ---- memory functions grgmem.c  -----*/

int     lsgrg_memsetup(LsgrgInfo *,
                       double xx[], double blvar[], double buvar[],
                       double blcon[],double bucon[], long lvars[],
                       long);
void    lsgrg_memstats(LsgrgInfo *);
void    lsgrg_insfmem(LsgrgInfo *,char *, long);
double *dbl_alloc( long);
long   *long_alloc( long);
char   *char_alloc( long);
char   **charp_alloc( long);

double *allocate_darray(LsgrgInfo *,char *,long);
long   *allocate_larray(LsgrgInfo *,char *,long);
char   *allocate_carray(LsgrgInfo *,char *,long);
char **allocate_pcarray(LsgrgInfo *,char *,long);

double *allocate_darray_nostats(LsgrgInfo *,char *,long);
long   *allocate_larray_nostats(LsgrgInfo *,char *,long);
char   *allocate_carray_nostats(LsgrgInfo *,char *,long);
char **allocate_pcarray_nostats(LsgrgInfo *,char *,long);

void    lsgrg_free_resources(LsgrgInfo *);
void lsgrg_zero_memstats(LsgrgInfo *);
int lsgrg_alloctable(LsgrgInfo *,char *, int, void*);
int lsgrg_alloc_rownames(LsgrgInfo *info,long nnrows);
int lsgrg_alloc_varnames(LsgrgInfo *info,long nnvars);
int lsgrg_free_rownames(LsgrgInfo *info);
int lsgrg_free_varnames(LsgrgInfo *info);

         /*   setup functions grgsetup.c  */

int    lsgrg_setupj(LsgrgInfo *,int build_jacobian, double xstat[],
                    double g[], double gg[], double gbest[],
                    double colb[],
                    long lgrad,long  *nnz, long *nnz_total,
                    long *nnznl, long *colmax,
                    long ipcol[], long iprow[], long ipmap[],
                    double paij[] );
void lsgrg_initialize_model(LsgrgInfo *,long nvars, long nrows);
void lsgrg_set_invert_tolerances(LsgrgInfo *);
int  lsgrg_setbounds(LsgrgInfo *,double blvar[],double buvar[],
                     double blcon[],double bucon[],double xx[],
                     long lvars [], long ivstat[]);
int  anajac(LsgrgInfo *,
         double x[], double grad[], long int ihag[],
         long int iheg[], long int ihegl[], double paij[], long int iprow[],
         long int ipcol[], long int ipmap[], long int icount[], long int *colmax,
         long int *nznl, long int nz, long int n);
void lsgrg_check_anajac(LsgrgInfo *,
                        long n, long nzgrad,double *grad,
                        long *ihag, long *iheg, long *iheg0,
                        long *ihegl,long *ihegl0,long *iprow,
                        long *ipcol,long *ipmap, double *paij);
void fdcheck(LsgrgInfo *,
         double x[], double grad[], long int ihag[],
         long int iheg[], double gplus[], double gminus[],
         double gcol[] );
void lsgrg_echo_lvars(LsgrgInfo *,
                      long nvars, long ivstat[], long lvars[],
                      long nnlin, long nlin);

double lsgrg_machine_precision(void);

   /*   timing functions grgtime.c */
double lsgrg_timer(void);
void   lsgrg_timestats(LsgrgInfo *,int, double);

   /*  file grgsub.c */

int grgsub(LsgrgInfo *,long nvars_in, long nfun_in, long nobj_in,long maximize,
            long lvars[],
            double blvar[],double buvar[],double blcon[],double bucon[],
            double xx[],double fcns[],double rmults[],long nonbas[],
            double redgr[],long inbind[],long *nbind,long *nnonb,
            P_GCOMP p_user_gcomp,P_PARSH p_user_parsh, long);
int lsgrg_setrowname(LsgrgInfo *,long nnrows,long index,char *name);
int lsgrg_setvarname(LsgrgInfo *,long nnvars,long index,char *name);
void lsgrg_echobnd(LsgrgInfo *,double xx[],double blvar[], double buvar[],
                    double blcon[],double bucon[], long lvars[],
                    char *msg);
void lsgrg_initialize_grgtest(LsgrgInfo *);
char *lsgrg_get_terminationmsg(LsgrgInfo *);
int lsgrg_xfer_alg_term_code(LsgrgInfo *,int info);
void lsgrg_set_termin_msg(LsgrgInfo *,char *buffer);
void lsgrg_set_termin_code(LsgrgInfo *,char *string);
double lsgrg_sinf(long m, long nobj,double g[], double blcon[],
                  double bucon[]);
int lsgrg_resize_jacobian(LsgrgInfo *_info,long *iheg) ;
int lsgrg_resize_binv(LsgrgInfo *_info) ;

/* **new 2/01 jcp **  */
int lsgrg_get_retry_history(LsgrgInfo *_info,
                         int  retry_history[], int *len_retry_history);
int lsgrg_fixup_strategy(LsgrgInfo *_info,int return_code,
                        int *retry_count, int max_retry_count,
                        char *msg, int retry_history[]);
void lsgrg_int_to_string(int n, char s[]);
void LSGRGCLASS LsgrgProblemSetup(LsgrgInfo *info);
void LSGRGCLASS LsgrgProblemRun(LsgrgInfo *info);
void LSGRGCLASS LsgrgProblemShutdown(LsgrgInfo *info);
/* 7/01 jcp */
void LSGRGCLASS lsgrg_disable_retries(LsgrgInfo *info);
void LSGRGCLASS lsgrg_enable_retries(LsgrgInfo *info) ;
void LSGRGCLASS lsgrg_return_mults_unpacked(LsgrgInfo *info) ;
void LSGRGCLASS lsgrg_return_mults_packed(LsgrgInfo *info) ;
void LSGRGCLASS lsgrg_disable_copyback(LsgrgInfo *info) ;
void LSGRGCLASS lsgrg_enable_copyback(LsgrgInfo *info) ;

/* ** end new */
  /*  echo routines file grgoutr.c */

void lsgrg_outres(LsgrgInfo *,int info);
void lsgrg_tablin(LsgrgInfo *);
void report(double g[], double x[], long mp1, long n, int usernames,
            char **con, char **var, double x0[]);
void lsgrg_echoparms(LsgrgInfo *);
void lsgrg_output_headers(LsgrgInfo *,int,int,int);
void lsgrg_echo_globals(LsgrgInfo *,char *);
void lsgrg_formatname(char *name,int length);
void lsgrg_printgrad(LsgrgInfo *,int print_slacks, char *msg);
void lsgrg_printgrad0(LsgrgInfo *,int print_slacks, char *msg,long n,
                      long nvars,long m,long nrows,long mp1,
                      long nzgrad, long nzlin, long nznlin,
                      long iheg[], long ihag[] , double grad[] ,
                      double x[]);
void lsgrg_prt_larray(LsgrgInfo *,char *msg,long x[], long n);

      /* parameter initializations file grginit.c */

LsgrgInfo *LsgrgInitialize(void);
void LsgrgBye(LsgrgInfo *);
void lsgrg_initialize_1time(LsgrgInfo *);
void lsgrg_set_defaults(LsgrgInfo *);
void lsgrg_put_globals(LsgrgInfo *);
void lsgrg_get_globals(LsgrgInfo *);
void lsgrg_zero_alg_globals(LsgrgInfo *);
void lsgrg_zero_shell_globals(LsgrgInfo *);

     /*  option processing functions grgopt.c */

int lsgrg_setparameter(LsgrgInfo *,char *name,long lvalue,double dvalue);
int lsgrg_post_option(LsgrgInfo *,int index, long lvalue, double dvalue,
                      int post);
int lsgrg_check_options(LsgrgInfo *);
void lsgrg_upcase(char *);
void lsgrg_lcase(char *);
double lsgrg_get_plinfy(LsgrgInfo *);
void lsgrg_set_title(LsgrgInfo *,char *);
void lsgrg_setcomputedparms(LsgrgInfo *);
/*  ** 4/27/02 jcp ** add prototype for function to set maxtime */
void lsgrg_set_maxtime(LsgrgInfo *lsinfo,double t);

void lsgrg_set_fderivatives(LsgrgInfo *);
void lsgrg_set_cderivatives(LsgrgInfo *);
void lsgrg_set_userderivatives(LsgrgInfo *);
int  lsgrg_set_printlevel(LsgrgInfo *,long ipr);
int  lsgrg_set_inprnt(LsgrgInfo *,long OnOff);
int  lsgrg_set_otprnt(LsgrgInfo *,long OnOff);
void lsgrg_set_scalingOn(LsgrgInfo *);
void lsgrg_set_scalingOff(LsgrgInfo *);

void lsgrg_steepest_descent(LsgrgInfo *opt,int OnOff);
void lsgrg_limit_stepsize(LsgrgInfo *opt,double steplim);

/*  lsgrg_errorexit -- contains setjmp/longjmp code           */
/*  for returns from algorithm routines                       */
void lsgrg_errorexit(jmp_buf buf,int exit_status_code);
/*-------------------------------------------------------------------*/
/*  end lsgrgc shell  function prototypes                            */
/*--------------------------------------------------------------------*/
/*====================================================================*/
/*   enums for parametrized values                                    */
/*====================================================================*/
/* lsgrgc termination codes were included from lscodes.h at top       */
/*   Negative values are non-algorithmic terminations                 */
/*   positive are algorithmic terminations                            */
/*--------------------------------------------------------------------*/
/*  maximize/minimize codes */
      enum lsgrg_optimization_sense {
             LSGRG_MINIMIZE=0, LSGRG_MAXIMIZE=1 };

/*  lsgrg algorithm module termination codes */
       enum lsgrg_alg_info_codes {
            _LSGRGALG_KTC          = 0,
            _LSGRGALG_FC           = 1,
            _LSGRGALG_ALLREMEDIES  = 2,
            _LSGRGALG_ITN_LIMIT    = 3,
            _LSGRGALG_UNBOUNDED    = 4,
            _LSGRGALG_INF          = 5,
            _LSGRGALG_RUNTIME_ERR  = 6,
            _LSGRGALG_USER_TERMIN  = 7,
            _LSGRGALG_TIME_LIMIT   = 8
   };
       enum command_mask {
         LSGRG_SETUP = 1, LSGRG_RUN = 2, LSGRG_SHUTDOWN = 4,
         LSGRG_SETUP_RUN  = 3,
         LSGRG_RUN_SHUTDOWN = 6,
         LSGRG_SETUP_RUN_SHUTDOWN = 7 };

         /* ----- io related stuff ----- */

       enum io_units {
         LSGRG_LOG_UNIT = 1, LSGRG_ERROR_UNIT = 2,
         LSGRG_SCREEN = 3 };

        enum lsgrg_derivative_options {
        FORWARD_DIFFERENCES = 0, CENTRAL_DIFFERENCES = 1,
        USER_ANALYTIC = 2 };

       enum lsgrg_timstats_options {
        T_INIT   = 1, T_PRINT  = 2, T_GCOMP = 3,
        T_CONSBS = 4, T_PARSH   = 5, T_NEWT  = 6, T_NEWTONLY = 7,
        T_SEARCH = 8, T_INVERT = 9, T_PHASE0 = 10, T_PHASE1 = 11,
        T_TOTAL = 12 };

       enum lsgrg_output_parms {
       FUNCTIONS = 1, VARIABLES = 2, INITIAL_VALUES = 3,
       FINAL_RESULTS = 4 };


     /*  declarations for parameter setting subsystem */

       enum lsgrg_user_options { LSGRG_OPTION_START = 0,
             KDERIV = 1,
             MAXBAS =  2, MAXHES =  3, EPNEWT =  4, EPINIT =  5,
             EPSTOP =  6, EPSPIV =  7, PH1EPS =  8, NSTOP  =  9,
             ITLIM  = 10, LIMSER = 11, IPR    = 12, IPN4   = 13,
             IPN5   = 14, IPN6   = 15, IPER   = 16, PSTEP  = 17,
             IQUAD  = 18, MODGG  = 19, AIJTOL = 20, PIVPCT = 21,
             MXTABU = 22, FUNPR  = 23, CONDTL = 24, IDEGLM = 25,
             EPBOUN = 26, EPDEG  = 27, ISCAL  = 28, ISCLG  = 29,
             MEMCG  = 30, IBVBLM = 31, HRDBND = 32, FIXPIV = 33,
             INPRNT = 34, OTPRNT = 35, GFEAS  = 36, USEPH0 = 37,
             STEEPD = 38, STEPLIM = 39,
/* ** 4/27/02 jcp ** add MAXTIME */
             MAXTIME = 40,

                LSGRG_OPTION_END = 41   };

       enum mem_tbl_cmds {
              LSGRG_MEMALLOC,LSGRG_MEMFREE };
       enum var_status_codes {
              V_LINEAR_FIXED     = -2,
              V_NONLINEAR_FIXED  = -1,
              V_NONLINEAR        =  1,
              V_LINEAR           =  2   };

       enum outres_print_levels {
              OUTRES_PRINT_ALL = 0,
              OUTRES_PRINT_ROWCOL = 1 };
/*====================================================================*/
/*     prototypes for algorithm routines                              */
/*====================================================================*/
/* --- prototypes for helper functions     -----*/

double   powi(double, long int);  /* raise double to int power */
long int ipow(long int, long int);  /* raise int  to int power */
void printbasisInv(LsgrgInfo *,char *);

/*  ---- prototypes for algorithm functions ----- */


void addcol(LsgrgInfo *_info,double[],double[],long,double);
void calfun(LsgrgInfo *_info,double[],double[],double[]);
void cg(LsgrgInfo *_info,double[],double[],double[],long*,long*,double*,double*,
        double*);
void chkelt(LsgrgInfo *_info,double[],long[],long[],long[],long[],LOGICAL32);
void chkfun(LsgrgInfo *_info,double[],double[],double[],double[],STRING);
void chuzq(LsgrgInfo *_info,double[],long[],long[],long[],double[],double[],double[],
        long[],double[],double[],double[],long[],long[],long[],
        long[],long*,double[],long[]);
void chuzr(LsgrgInfo *_info,double[],double[],double[],double[],double[],long[],long[],
        long[],long*,double*,LOGICAL32*);
void comdfp(LsgrgInfo *_info,double[],double[],double[],long,double,double*,double,
        double*,long*,long*);
void condnm(LsgrgInfo *_info,long,short*,double[],long,short*,long*,long[],
        double[],long[],double[],double[],long[],long[],long[],
        double[],double[]);
void consbs(LsgrgInfo *_info,long[],double[],double[],long[],long[],long[],long[],
        double[],double[],double[],long[],long[],long[],double[],
        long[],long[],long[],double[],double[],double[],long[],
        long[],double[],long[],double[],long[],double[],long[],
        long*,long[],double[],long[],long[],long[]);
void delcol(LsgrgInfo *_info,double[],long);
void direc(LsgrgInfo *_info,double[],double[],double[],long[],long[],double[],long[],
        double[],long[],long[],long[],double*,double*,double*);
void droprw(LsgrgInfo *_info,long,long[],long[],long*,long*,long*,double*,
        double*,double[],double[],double[],double[],double[],double*);
void getbas(LsgrgInfo *_info,double[],long[],long[],double[],short*,long[],long,
        long*);
double getcb(LsgrgInfo *_info,long,long[],long,long,double[]);
void getgr(LsgrgInfo *_info,double[],double[],long[],long[],long[],long[],double[],
        double[],double[],double[],double[],double[],double[],LOGICAL32*,
        long[],double[],long[],long[],long[]);
void grdump(LsgrgInfo *_info,double[],long[],long[],long[],long[],long,long,
        long,long,long,STRING);
int  grgitn(LsgrgInfo *_info,long[],double[],double[],long[],long[],long[],double[],
        double[],double[],double[],double[],double[],double[],double[],
        double[],double[],double[],double[],long[],long[],double[],
        double[],long[],long[],long[],double[],long[],double[],
        double[],double[],double[],double[],double[],double[],double[],
        long[],long[],long[],long[],long[],double[],long*,
        double[],double[],long[],double*,double*,double*,long[],
        double[],long[],long[],long[]);
void itrlog(LsgrgInfo *_info,long,long,long,long,long,long,long,
        long*,long,long*,long,long*,LOGICAL32,LOGICAL32,
        long,double,double,double[],double[],double,double);
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
void m2scal(LsgrgInfo *_info,long,long,long[],long[],long[],long[],double[],
        double[],long[],double[],double[]);
void modjac(LsgrgInfo *_info,long,long,long,double[],double[],long[],long[],
        double[],long*,long,double,long,long,LOGICAL32,
        long*,long[],long[]);
void newton(LsgrgInfo *_info,long[],double[],double[],long[],double[],double[],double[],
        long[],long[],double[],double[],double[],LOGICAL32*,double[],
        double[],double[]);
/*   **fixme** take these out of here
void parsh_gams(double[],double[],long*,long*,double[]);
void parsh_lsgrg(double[],long*,double[],long[],long[],double[],long[],
        long[],long[],long*);
void parshc(double[],double[],double[],long,double[],double[]);
void parshf(double[],double[],double[],long,double[],double);
*/

void ph0fac(LsgrgInfo *_info,double[],double[],long[],long[],long[],long[],long[],
        double[],long[],double[],double[],double[],double[],double[],
        double[],double[],double,double[],long[],LOGICAL32,LOGICAL32*,
        long[],double[],long[],long[],long[]);
void ph0log(LsgrgInfo *_info,long,long,long,long,long,long,long*,
        long*,long*,LOGICAL32,LOGICAL32,double,double,double,
        double,long);
void ph0piv(LsgrgInfo *_info,double[],long[],long[],long[],double[],double[],double[],
        long[],double[],double[],double[],long[],long[],long[],
        long[],double[],long[],long,double[],LOGICAL32*,long[],
        LOGICAL32*);
void ph1obj(LsgrgInfo *_info,long[],double[],double[],double[]);
void phas0(LsgrgInfo *_info,long[],double[],double[],long[],long[],long[],long[],
        double[],double[],double[],long[],long[],long[],double[],
        long[],long[],long[],double[],double[],double[],long[],
        long[],double[],long[],double[],long[],double[],long[],
        double[],double[],LOGICAL32*,long[],long*,long[],double[],
        long[],long[],long[]);
void prscal(LsgrgInfo *_info,long[],long[],long[],long[],double[],double[],long[],
        long[],double[],double[],double[],double[],double[],double[],
        double[],double[],double[],long[],double[],long[],long[],
        long[]);
void lsquad(LsgrgInfo *_info,long[],double[],double[],double[],double[],double[]);
void r1add(LsgrgInfo *_info,double[],double[],long);
void r1mod(LsgrgInfo *_info,double[],double[],long,double,long*);
void r1sub(LsgrgInfo *_info,double[],double[],long,double,LOGICAL32*);
void redgra(LsgrgInfo *_info,double[],long[],double[],double[],long[],long[],double[],
        long[],long[],double[],double[],double[],long[]);
void redobj(LsgrgInfo *_info,long[],double[],double[],long[],long[],double[],double[],
        long[],double[],double[],double[],double[],long[],double[],
        double[],double[],double[],double[],double[],long[],double[],
        double[],double[],double[]);
void resetr(LsgrgInfo *_info,double[],long,double*);
void rtrsol(LsgrgInfo *_info,double[],double[],long);
void scljac(LsgrgInfo *_info,double[],long[],long[],long[],long[],double[],LOGICAL32);
void search(LsgrgInfo *_info,long[],double[],double[],long[],long[],double[],double[],
        double[],long[],double[],double[],double[],double[],long[],
        double[],double[],double[],double[],double[],long[],double[],
        double[],double[],double[]);
void setobj(LsgrgInfo *_info,long[],long[],long[]);
void sumres(LsgrgInfo *_info,double[],double[],double[],double[],double[],long[],double*,
        double*,double*,long*,long*);
void tang(LsgrgInfo *_info,long[],double[],double[],long[],long[],long[],double[],
        double[],double[]);
void xbtran(LsgrgInfo *_info,long[],double[],double[],long);
void xdot(LsgrgInfo *_info,double[],long[],long[],double[],long,long,double*);
void xerror(long,long);
void xfact(LsgrgInfo *_info,double[],long[],long[],long[],double[],long[],LOGICAL32,
        long*,double[],long[],LOGICAL32);
void xftran(LsgrgInfo *_info,long[],double[],double[],long);
void xpivot(LsgrgInfo *_info,double[],long[],long[],long[],double[],long[],long,
        long,double[],long*,double[],long[],LOGICAL32);
void xsaxpy(LsgrgInfo *_info,double[],long[],long[],double[],double,long,long);

void prtha(LsgrgInfo *_info,long *ha, char *msg);
void dumpcol(LsgrgInfo *_info,char *msg, long *ha, long icol, long *he);

void prtdarray(LsgrgInfo *_info,char *name,char *msg, double x[], long n);

#endif
