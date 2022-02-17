/******************************************************************************\
*
*   MODULE              :   grf_h_types
*
*   DESCRIPTION         :   Header file containing the declaration of GRFN
*                           internal data types.
*
*   NOTE                :   GRFN Language data types are in grf_h_lang.h
*                           Public data types are in grf_h_pubtypes.h
*
\******************************************************************************/

#ifndef GRF_H_TYPES_H
#define GRF_H_TYPES_H

/******************************************************************************\
*                           Import Include Files                               *
\******************************************************************************/

#include "grf_h_lang.h"
#include "grf_h_pubtypes.h"
#include "srt_h_sample.h"

/* ------------------------------------------------------------------------- */

typedef enum {
  GRFN_SS_CPN,
  GRFN_SS_BAL,
  GRFN_SS_VA,
  GRFN_SS_VB,
  GRFN_SS_VC,
  GRFN_SS_VD,
  GRFN_SS_VE,
  GRFN_SS_VF,
  GRFN_SS_VG,
  GRFN_SS_VH,
  GRFN_SS_V0,
  GRFN_SS_V1,
  GRFN_SS_V2,
  GRFN_SS_V3,
  GRFN_SS_V4,
  GRFN_SS_V5,
  GRFN_SS_V6,
  GRFN_SS_V7,
  GRFN_SS_V8,
  GRFN_SS_V9,
  GRFN_SS_V10,
  GRFN_SS_V11,
  GRFN_SS_V12,
  GRFN_SS_V13,
  GRFN_SS_V14,
  GRFN_SS_V15,
  GRFN_SS_V16,
  GRFN_SS_V17,
  GRFN_SS_V18,
  GRFN_SS_V19,
  GRFN_SS_V20,
  GRFN_SS_V21,
  GRFN_SS_V22,
  GRFN_SS_V23,
  GRFN_SS_V24,
  GRFN_SS_V25,
  GRFN_SS_V26,
  GRFN_SS_V27,
  GRFN_SS_V28,
  GRFN_SS_V29,
  GRFN_SS_NUM_PARAM
} GrfnSSParamName;

/******************************************************************************\

    STRUCTURE           :   GrfnSSParam

    DESCRIPTION         :   Contains values of all pre-declared variables in
                            the GRFN language viz. VA  , VB  , etc.

                            During compilation  , COMLL nodes that contain
                            references to these variables are populated with
                            pointers to the locations of VA  , VB  , ... in this
                            structure  , which itself is contained within a
                            GrfnDeal.

    NOTE                :   In order to declare more variableso modify
                            GRFN_SS_NUM_PARAM  , declare the variable in
                            GRFN symbol table grf_f_lngsymtab.c

\******************************************************************************/

typedef double GrfnSSParam[GRFN_SS_NUM_PARAM];
typedef double *GrfnLatticeVec; /* Representation of data being passed
                                   through a lattice */

/******************************************************************************\

    STRUCTURE           :   GrfnEvent

    DESCRIPTION         :   Contains information on a financial event occuring
                            on a particular date described in the GRFN tableau.

                            The action we need to take on an "event date" is
                            contained in a COMLL. Therefore a GRFN event will
                            need to contain information required to evaluate
                            the COMLL. This is performed by grfn_eval_event.

    ELEMENTS            :   t       Summary of grammatical information
                            d       Date in long integer format
                            dflen   Length of dft  , dfd  , df and yp
                            dft     Times in years at which we need to eval df
                            dfd     Times as dates as above
                            df      Discount factors from time t to t+dt
                                    df[i] = df from t to dft[i]
                            yp      Information required for computation of
                                    zero coupon yields once r(t) and other
                                    state variables are known
                                                                        (Y_T_at_t)
                            icom    Command Linked List

                                                        sprdlen Length of spread
                                                                        ( = 0 if
no spread  , = dflen otherwise) spread  Spread of the reference rate to cash
rate  , stored on payment date (end of period)  , for all the df dates
\******************************************************************************/

typedef struct {
  GrfnCellStatus status;
  double t;
  Date d;
  COMLL_PTR icom;
  long dflen[MAXUNDERLYING];
  double *dft[MAXUNDERLYING];
  Date *dfd[MAXUNDERLYING];
  double *df[MAXUNDERLYING];
  void *yp[MAXUNDERLYING];

} GrfnEvent;

/* The following will need to go in a SRT header */

typedef enum {
  GRFNNOMETHOD,
  GRFNMONTECARLO,
  GRFNBACKWARDLATTICE,
  LASTGRFNMETHOD
} GrfnMethod;

/* --------------------------------------------------------------
   TYPE         :  SrtFncWrapper
   DESCRIPTION  :  Struct to wrap a pointer to a function.
                   This is useful for passing pointers to
                                   functions to other machines ---
                                   A function will have a different address on
                                   other machine  , so pointer will be invalid.
                                   So we use the name field to find the new
  address.
  --------------------------------------------------------------*/

typedef struct {
  char name[32];

  Err (*fnc)();

} SrtFncWrapper;

/******************************************************************************\

    STRUCTURE       :   GrfnDeal

    DESCRIPTION     :   Contains the following information:

                        (1) Parsing environment for the creation of GrfnEvents
                            i.e. all data EXCEPT MARKET DATA that the GRFN
                            compilation routines (e.g. grf_f_lang  , etc) will
                            need in order to create GrfnEvents.  This includes
                            the locations of the areas that will be used to
                            keep track of path dependent data  , viz. the place
                            where values of the cells in the tableau are stored.

                        (2) Co-ordinates of the GrfnCell currently being
                            compiled  , and the time of the GrfnEvent currently
                            being created.

                        (3) Values of each cell after evaluation together with
                            values of variables  , etc.

                        (4) Diagnostic information that results from parsing a
                            GRFN tableau  , i.e. whether or not a  particular
                            numerical method has been forced.  (??? Not needed)

\******************************************************************************/

typedef struct {
  Date today;             /* Today                                  */
  Date *event_dates;      /* Days of input financial events         */
  int num_event_dates;    /* Length of event_dates                  */
  int sswidth;            /* Cols in gcells                         */
  int sslength;           /* Rows in gcells                         */
  GrfnCell **gcells;      /* Input spreadsheet                      */
  double **cells;         /* Work area returned by @grfn_cell       */
                          /* Used to keep track of path dependence  */
  int first_unkn_index;   /* Index of first unknown financial event */
  int end_of_day_index;   /* If end-of-day flag is true=1  , else=0   */
  int end_of_day_fixing;  /* If end-of-day fixing flag is true=1  , else=0   */
  int end_of_day_payment; /* Due to misuse of the two previous ones.
                                                                              This
                             one is only use by the GRFN language PAY function
                           */

  GrfnSSParam ssparam;       /* Predeclared vars (va  , vb  ,...)          */
  GrfnSSParam hist_ssparam;  /* Value of pre-declared vars usable after
                                historical events have been evaluated  */
  long sptr_last_index;      /* Index of last element of sptr          */
  GrfnMethod method;         /* Numerical method                       */
  GrfnCellStatus sumstatus;  /* Sum of grammatical information over
                                entire grfn deal                       */
  GrfnCellStatus *rowstatus; /* Sum of grammatical info for strings
                                in different rows of gcells            */
  char domestic_und[32];     /* Name of domestic (P&L) underlying  */
  char disc_curve[32];       /* Name of the associated discount curve */

  Err(*fixing_fct) /* The function to get a fixing from the outside */
      (long, char *, double *); /* (should be of FixingFuncType)  */

  double pv_of_past; /* Cashflows fixed but not yet received   */

  /* The function to return payments and dates to the outside */

  long *pay_dates; /* Dates at which known payments will occur (PAY fct) */
  long *event_pay_dates; /* Event dates on which the pay function is called */
  double *pay_amounts;   /* Amounts to be paid on these dates (PAY fct) */
  long numpay;

  int pass;               /* If 1 parse tableau  , else create even   */
  SRT_Boolean is_history; /* If current date is in the past         */
  SRT_Boolean american;   /* If input string now being evaluated is
                         an auto generated string corresponding
                         to an american cashflow                */
  long amdays;            /* If (american)  , no. of days between
                             generated event time and input event   */
  double amdt;            /* Fraction of a day to add to amdays     */
  int I, J;               /* Row and col of current cell in gcells. */
  Ddate nowdt;            /* If grfn deal is being initialized  , now */
  Ddate prvdt;            /* If grfn deal is being initialized  ,
                             the previous event date                */
  Ddate nxtstpdt;         /* If grfn deal is being initialized  ,
                             the date of the next step              */
  Ddate nxtevdt;          /* If grfn deal is being initialized  ,
                             the date of the next prespecified event*/
  double **aux;           /* Auxiliary information                  */
  long *auxlen;           /* Lengths of arrays on aux i.e.
                              aux[i] is an array running from
                              aux[i][0] to aux[i][auxlen[i]-1]      */
  long auxwidth;          /* Length of auxlen and aux  , i.e. aux and
                              auxlen go from aux[0] and auxlen[0] to
                              aux[auxwidth-1] and auxlen[auxwidth-1]*/
  GrfnRng *grng;          /* Extra named ranges (unused)            */
  int num_grng;           /* Length of grng                         */
  FILE *outfileptr;       /* Debugging information                  */

} GrfnDeal;

/******************************************************************************\
*                       Macro Functions for GRFN Cells                         *
\******************************************************************************/

#define GrfnCGCell(gd) gd->gcells[gd->I][gd->J]
#define GrfnCCell(gd) gd->cells[gd->I][gd->J]

#define GrfnIsStatus(status, mask) (status[mask])
#define GrfnAddStatus(status, mask)                                            \
  { status[mask] = SRT_YES; }

void GrfnResetStatus(GrfnCellStatus status);
void GrfnAggregateStatus(GrfnCellStatus ful_status, GrfnCellStatus cur_status);
void GrfnCopyStatus(GrfnCellStatus new_status, GrfnCellStatus old_status);

/******************************************************************************\
*                       Macro Constants for GRFN Error Messages                *
\******************************************************************************/

#define GRERR_WRONG_NUM_OPT_ARG "Too many optional arguments"
#define GRERR_WRONG_TYPE_OPT_ARG "Optional argument has incorrect type"
#define GRERR_UNKNOWN_SYMBOL "Unknown symbol"
#define GRERR_WRONG_NUM_ARGS "Invalid number of arguments"
#define GRERR_REF_OUT_OF_RANGE "Reference out of tableau range"
#define GRERR_INDEX_OUT_OF_RANGE "Index wanted out of range"
#define GRERR_COL_OUT_OF_RANGE "Bad column reference"
#define GRERR_ROW_OUT_OF_RANGE "Bad row reference"
#define GRERR_ILLEGAL_REF "Illegal reference to future"
#define GRERR_BAD_TYPE "Argument has incorrect data type"
#define GRERR_BAD_START "Start must be >= event date"
#define GRERR_BAD_COMPD "Invalid compounding"
#define GRERR_BAD_CALL_PUT "Invalid call/put type"
#define GRERR_BAD_DIFFUSION_TYPE "Invalid lognormal/normal type"
#define GRERR_BAD_BASIS "Invalid basis"
#define GRERR_BAD_END "Invalid end date"
#define GRERR_BAD_END_OR_NFP "Invalid end or nfp"
#define GRERR_BAD_REF_RATE "Invalid Reference Rate Name"
#define GRERR_BAD_DATE_REF "Reference to date doesn't exist"
#define GRERR_BAD_ASSIGNMENT "Invalid assignment :="
#define GRERR_INTERP_SAME_SIZE "xlength != ylength"
#define GRERR_BAD_AUX_RNG "Auxiliary range doesn't exist"
#define GRERR_NOT_DETERMINISTIC "Invalid or undetermined argument"
#define GRERR_NOT_NUMBER "Argument is not a number"
#define GRERR_METHOD_CONFLICT "Tree and Monte-Carlo don't mix"
#define GRERR_BAD_EVDIM "Invalid event or evdate dimension"
#define GRERR_NO_REFRATE "Need reference rate name (or value)"
#define GRERR_NO_HIST_FNC "Historical data not available"
#define GRERR_MISSING_FIXING "Could not get fixing"
#define GRERR_PAY_LAST_COL "PAY must be in cash-flow column"
#define GRERR_PAY_ALONE "PAY(...) must be used alone in a cell"
#define GRERR_MEMORY_ALLOC "Memory allocation error in "

#endif
