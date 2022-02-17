#ifndef GRFN_LANG_H
#define GRFN_LANG_H

#define COMLLBUFSZ 32

/******************************************************************************\

    TYPE            :   COMLLType

    DESCRIPTION     :   GRFN Language internal enumerated types

                        COMLLType determines how a particular node in a
                        Command Linked List is treated.

                        CLLNameType is a more general grouping of different
                        sorts of grfn language functions      , i.e. arithmetic
                        functions      , booleans      , etc.

\*******************************************************************************/

typedef enum {
  COMLL_NULL, /* syntax basics */
  COMLL_REAL,
  COMLL_STRING,
  COMLL_VAR,
  COMLL_USERVAR,
  COMLL_CONST,
  COMLL_JUMP,
  COMLL_BRANCH,
  COMLL_ASSIGN,
  COMLL_POP,
  /* 10 */
  COMLL_RESET_STATE,
  COMLL_F_LEVEL, /* financial functions */
  COMLL_F_SWAP,
  COMLL_F_CMS,
  COMLL_F_DF,
  COMLL_F_PAY,
  COMLL_F_PAYFROMPAST,
  COMLL_F_PAYINPAST,
  COMLL_F_FRA,
  COMLL_F_SWAPTION,
  /* 20 */
  COMLL_F_PVRANGE,
  COMLL_F_CVG,
  COMLL_F_ACCINT,
  COMLL_F_YIELD,
  COMLL_F_CLEAN,
  COMLL_A_PLUS, /* Arithmetic functions */
  COMLL_A_TIMES,
  COMLL_A_MINUS,
  COMLL_A_DIVIDE,
  COMLL_A_LOG,
  /* 30 */
  COMLL_A_FLOOR,
  COMLL_A_CEILING,
  COMLL_A_EXP,
  COMLL_A_POW,
  COMLL_A_MAX,
  COMLL_A_MIN,
  COMLL_A_ABS,
  COMLL_A_NORM,
  COMLL_A_IF,
  COMLL_A_BLKSCHLS,
  /* 40 */
  COMLL_A_BLKSCHLSNRM,
  COMLL_B_LT, /* SRT_Boolean (0 or 1) functions */
  COMLL_B_GT,
  COMLL_B_LE,
  COMLL_B_GE,
  COMLL_B_EQ,
  COMLL_B_OR,
  COMLL_B_AND,
  COMLL_C_COLMAX, /* Column functions */
  COMLL_C_COLMIN,
  /* 50 */
  COMLL_C_COLSUM,
  COMLL_C_COLAVG,
  COMLL_C_COLSORT,
  COMLL_C_ROWMAX, /* Row functions */
  COMLL_C_ROWMIN,
  COMLL_C_ROWSUM,
  COMLL_C_ROWAVG,
  COMLL_C_ROWSORT,
  COMLL_C_ROWINTERP,
  COMLL_X_INTERP, /* Functions on auxiliary ranges */
  /* 60 */
  COMLL_T_REF,    /* References to PV vector */
  COMLL_T_CUMREF, /* References to current and cumulative values FWD */
  COMLL_T_CURREF,
  COMLL_T_SET,
  COMLL_T_INCREMENT,
  COMLL_S_UNDERLYING, /* References to underlying */
  COMLL_S_SHORT_RATE,
  COMLL_S_PHI,
  COMLL_S_SIGMA,
  COMLL_S_STATEVAR,
  /* 70 */
  COMLL_S_NUMERAIRE,
  COMLL_S_MMACCOUNT,
  COMLL_S_STOCK,
  COMLL_PRINT, /* Miscellaneous */
  COMLL_ERROR,
  COMLL_HIST,
  COMLL_C_ROWMAXIDX, /* Added by Cyril Godart for Regis Benichou */
  COMLL_C_ROWMINIDX,
  COMLL_CONSTF,
  COMLL_A_BSSABR,          /*for Jerome and end of added by Cyril Godart */
                           /* 80 */
  COMLL_A_NUMDAYS,         /* number of days in a corridor */
  COMLL_F_PAYOFF,          /* payoff of a product */
  COMLL_T_RISKY_REF,       /* Added by Matthieu Autret for GRFNCR*/
  COMLL_T_PVSOURCE_ASSIGN, /* Added by Matthieu Autret for GRFNCR*/
  COMLL_X_PVINTERP,        /* Interp columns */
  COMLL_A_SPREADNUM,
  COMLL_F_PAYN

} COMLLType;

/******************************************************************************\

    TYPE            :   CLLNameType

    DESCRIPTION     :   Command Linked List Name Type subdivided Arithmetic
                        Functions      , SRT_Boolean Functions      , etc ...
                        CLLNameType is ONLY used during compilation of
                        the GRFN language (during parsing)

\******************************************************************************/

typedef enum {
  CLLNFinFunc,
  CLLNArithFunc,
  CLLNCellFunc,
  CLLNVar,
  CLLNUserVar,
  CLLNReal,
  CLLNConst,
  CLLNDateRef,
  CLLNCellRef,
  CLLNPvRef,
  // IMPLEMENTATION OF PVR
  CLLNPvRRef,       // Added by Matthieu Autret
  CLLNSourceAssign, // Added by Matthieu Autret
  CLLNCumRef,
  CLLNCurRef,
  CLLNCvg,
  CLLNI,
  CLLNJ,
  CLLNNow,
  CLLNNxt,
  CLLNPrv,
  CLLNAuxFunc,
  CLLNAuxRef,
  CLLNSample,
  CLLNUtil
} CLLNameType;

/******************************************************************************\

    TYPE            :   GrfnSymbol

    DESCRIPTION     :   Information about GRFN functions and names.
                        Some of this information is also within grf_f_ytab.c
                        and grf_f_lex_yy.c (YACC and LEX files)      ,  which
take care of functions with one letter such as '+' directly.

                        See grf_f_lngsymtab.c for the table of all GrfnSymbols
                        currently understood by GRFN.

\******************************************************************************/

typedef struct {
  String name;          /* Name of symbol                   */
  CLLNameType nametype; /* Type of symbol                   */
  COMLLType type;       /* ENUM of this particular symbol   */
  int ndargs;           /* Number of deterministic args     */
  int nrargs;           /* Number of non-deterministic args */
  int no_opt_args;      /* Number of optional args          */
  int ival;             /* Useful field                     */
  String helpstr;       /* Help string about this argument  */
} GrfnSymbol;

/******************************************************************************\

    STRUCTURE       :   COMLL_PTR      , COMLL_STR

    DESCRIPTION     :   COMmand Linked List.

                        A command like f(g(x)      ,y      ,4+z) * 2 is
represented as a COMLL_PTR as follows

                        4 -> z -> + -> y -> x -> g -> f -> 2 -> * -> NULL

                        Order of evaluation is as in "Reverse Polish".


    NOTE            :   [1] Expressions are simplified as mush as possible.
                            e.g.  2+5 is reduced to 7      , as is d[i]      ,
i+2 (since the symbol i is known)      , max      , pow      , cvg      , etc...

                        [2] Financial functions have their arguments evaluated
                            before hand. Hence swap(s      , e      , c      ,
b) + 2 is represented as:

                                2 -> swap(...) -> +

                                (NOT "2"->"b"->"c"->"e"->"s"->"swap"->"+")

                            i.e. each node in the linked-list contains a place
                            to store the dates of a swap rate.

                        [3] There are several utility functions that operate on
                            a COMLL_PTR. See grf_f_comlst.c for details.

\******************************************************************************/

typedef char comllchar[COMLLBUFSZ];

typedef enum { FWDCELL, BWDCELL } COMLLEvalstatus;

struct COMLL_STRUCT {
  struct COMLL_STRUCT *prev;
  struct COMLL_STRUCT *next;
  struct COMLL_STRUCT *jmp;     /* Used in "if" statements        */
  struct COMLL_STRUCT *nextcol; /* Used in ammc */

  COMLLEvalstatus evalstatus; /* Used in ammc */

  long index; /* Index OF COMLL in the list     */

  COMLLType type; /* Type of this comll             */

  int nargs; /* No of args (if comll=function) */

  void *gdcells; /* Comll points to gd->cells  */

  void *ptr; /* If comll points to something   */

  double dval; /* If comll is of type COMLL_REAL */
               /* Also used for underlying index */

  int ivec[4]; /* Used for index of rows or      */
               /* index of underlying in list    */

  comllchar sval; /* If comll is of type COMLL_STR  */

  long *dfind; /* Indices of df needed by COMLL  */
  long *
      dfindfloat; /* Indices of df needed by COMMLL
                                                                     for
                     floating leg of swap (spreads)
                                                                                                        */

  int dfindlen;      /* Length of dfind                */
  int dfindfloatlen; /* Length of dfindfloat           */

  double *cvg; /* Coverages of times indexed     */
               /* by dfind. From 0 -> dfindlen-1 */
               /* NB. cvg[0] is NOT used         */

  double *cvgfloat; /* Coverages of times indexed     */
                    /* by dfindfloat.
                        /* From 0 -> dfindfloatlen-1      */
                    /* NB. cvg[0] is NOT used         */

  double *spread; /* Spread for the rate/cash rates
                  /* From 0->dfindfloatlen-1        */
                  /* Spreads are stored on end date */
                  /* NB: spread[0] is then not used */

  double long_spread; /* Spread for long ref rates */

  double *dstore; /* Array to store double which */
                  /* are going to be used when */
                  /* evaluating */

  long dstorelen; /* Number of doubles stored in dstore */

  comllchar(*sstore); /* Array to store strings of */
                      /* COMLLBUFSZ characters which */
                      /* are going to be used when */
                      /* evaluating */

  long sstorelen; /* Number of string stores in sstore */
};

typedef struct COMLL_STRUCT COMLL_STR, *COMLL_PTR;

/******************************************************************************\
*                        Prototype Functions for COMLL                         *
\******************************************************************************/

COMLL_PTR comll_node_alloc();
COMLL_PTR comll_free_node(COMLL_PTR a);
COMLL_PTR comll_gototop(COMLL_PTR a);
COMLL_PTR comll_gotobot(COMLL_PTR a);
COMLL_PTR comll_free_list(COMLL_PTR a);
COMLL_PTR comll_insert_after(COMLL_PTR a);
COMLL_PTR comll_insert_before(COMLL_PTR a);
COMLL_PTR comll_join(COMLL_PTR a, COMLL_PTR b);
COMLL_PTR comll_copy(COMLL_PTR a);
void comll_cut(COMLL_PTR a);
int comll_atom(COMLL_PTR a, COMLLType t);
long comll_create_index(COMLL_PTR a);

/******************************************************************************\
*                       Macro functions used in language building              *
\******************************************************************************/

/******************************************************************************\

    NAME            :   test_real (index      , "name"      , comllptr)
    DESCRIPTION     :   Is comll of type COMLL_REAL.
                        If not return an error message.
                        Index is the index of comll as an argument to some
                        function.

\******************************************************************************/

#define test_real(I, N, C)                                                     \
  {                                                                            \
    if (C->type != COMLL_REAL) {                                               \
      return serror("%s #%d in %s", GRERR_NOT_NUMBER, I, N);                   \
    }                                                                          \
  }

/******************************************************************************\

    NAME            :   test_deterministic (index      , "name"      , comllptr)
    DESCRIPTION     :   Is comll of type COMLL_REAL or COMLL_STRING.

                        If not return an error message.
                        Index is the index of comll as an argument to some
                        function.
\******************************************************************************/

#define test_deterministic(I, N, C)                                            \
  {                                                                            \
    if (C->type != COMLL_REAL && C->type != COMLL_STRING) {                    \
      return serror("%s #%d in %s", GRERR_NOT_DETERMINISTIC, I, N);            \
    }                                                                          \
  }

/******************************************************************************\

    NAME            :   CheckAuxRange1(i      , grfndealptr)
    DESCRIPTION     :   Is i a valid index of some auxiliary range in
                        grfndeal.  If not      , return an error message.

\******************************************************************************/

#define CheckAuxRange1(I, GD)                                                  \
  {                                                                            \
    if ((I) < 0 || (I) > GD->auxwidth - 1)                                     \
      return serror("%s:%d", GRERR_BAD_AUX_RNG, (I));                          \
  }

/******************************************************************************\

    NAME            :   CheckAuxRange2 (i      ,j      ,grfndealptr)
    DESCRIPTION     :   Are i and j  both valid indices of auxiliary
                        ranges in grfndeal with the same length.
                        If not      , return an error message.

\******************************************************************************/

#define CheckAuxRange2(I, J, GD)                                               \
  {                                                                            \
    CheckAuxRange1(I, GD);                                                     \
    CheckAuxRange1(J, GD);                                                     \
    if (GD->auxlen[(I)] != GD->auxlen[(J)]) {                                  \
      return serror("%s: lengths %d      ,%d", GRERR_INTERP_SAME_SIZE,         \
                    gd->auxlen[(I)], gd->auxlen[(J)]);                         \
    }                                                                          \
  }

#endif
