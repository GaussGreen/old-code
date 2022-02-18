#ifndef GRF_H_PUBTYPES_H
#define GRF_H_PUBTYPES_H

#include "uttypes.h"
#include "utstring.h"
#include "uterror.h"
/* ------------------------------------------------------------------ */
/* default size of string in GrfnCell  (equal to the XL one )         */
/* ------------------------------------------------------------------ */

#define GRFN_DEF_ARGBUFSZ 255

/* ------------------------------------------------------------------ */

/* -----------------------------------------------------------------
  TYPE            :GrfnCellType
  DESCRIPTION     :Enumerated type of a Grfn cell
  DEFINITION      :
  ---------------------------------------------------------------------*/

typedef enum GrfnCellType
{
    GRFNDCELL, /* cell contains a number */
    GRFNSCELL, /* cell contains a string */
    GRFNBCELL  /* cell is blank */
} GrfnCellType;

/* ---------------------------------------------------------------------
   TYPE:          : GrfnCellStatusType;
   DESCRIPTION    : An enumerated type to describe all the potential
                    status in a GrfnCellStatus
  ---------------------------------------------------------------------- */
typedef enum
{
    GRFNCSAMERICAN = 0, /* Cell contains AM qualifier                 */
    GRFNCSFUTUREREF,    /* Cell contains reference to following row   */
    GRFNCSPASTREF,      /* Cell contains reference to previous row    */
    GRFNCSVARREF,       /* Cell contains pre-declared variable        */
    GRFNCSPVREF,        /* Cell contains PV[...]                      */
    GRFNCSPAYREF,       /* Cell contains Pay(...)                     */
    GRFNCSLASTSTATUS
} GrfnCellStatusType;

/* -----------------------------------------------------------------
  TYPE            :GrfnCellStatus
  DESCRIPTION     :structure used to represent various bits of information
   about the financial event described by a particular GrfnCell.
   Right now, this is a vector of SRT_Boolean, each element corresponding to
   the StatusType.
  ---------------------------------------------------------------------*/

typedef SRT_Boolean GrfnCellStatus[GRFNCSLASTSTATUS];

/* -----------------------------------------------------------------
  TYPE            :GrfnCell
  DESCRIPTION     :structure to represent what is contained in one cell of a
                   GrfnTableau.
                                   Note there are public Grfn functions to allocate and free
matrices of GrfnCells.  To use GRFN, the user is obliged to populate matrix of
GrfnCells.
  DEFINITION      :
  ---------------------------------------------------------------------*/

typedef struct
{
    GrfnCellType   type;        /* type of cell */
    GrfnCellStatus status;      /* summary of grammatical info (internal)*/
    double         dval;        /* input value of cell of type GRFNDCELL */
    String         sval;        /* input value of cell of type GRFNSCELL */
    SRT_Boolean    str_alloced; /* was sval allocated? */
} GrfnCell;

/* -----------------------------------------------------------------
  TYPE            : GrfnRngType, GrfnRng
  DESCRIPTION     : these do not yet exist; place is provided for them in some
GRFN functions as a sort of expansion slot for when extra inputs more
sophistocated than auxiliary ranges are finally supported.  In GRFN functions
expecting *GrfnRng, a NULL input should for now be used.
 ---------------------------------------------------------------------*/

/* one possible future defintion: */
typedef enum GrfnRngType
{
    GRFNRNGVAR,
    GRFNRNGCONST,
    GRFNRNGCONSTARRAY,
    LASTGRFNRNGTYPE
} GrfnRngType;

typedef struct
{
    GrfnRngType type;
    long        len;
    double*     val;
    double*     init;
    String*     name;
} GrfnRng;

/* ----------------------------------------------------------------------- */

#endif
