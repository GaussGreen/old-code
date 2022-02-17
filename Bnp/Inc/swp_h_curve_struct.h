/* =============================================================================

   FILENAME      : swp_h_curve_struct.h

   PURPOSE       : the Object used to store a Curve in a list

   =============================================================================
 */

#ifndef SWP_H_CURVE_STRUCT_H
#define SWP_H_CURVE_STRUCT_H

#include "srt_h_fwdobj.h"
#include "srt_h_repo_obj.h"
#include "srt_h_types.h"
#include "swp_h_ccy_param.h"
#include "swp_h_cmt_types.h"

/* -----------------------------------------------------------------------------
           THE GENERIC CURVE STORAGE STRUCTURE AND THE CURVE OBJECT ITSELF:
               a SrtCurveDesc (stored in the Double linked list)
   -----------------------------------------------------------------------------
 */

typedef struct {
  char curve_name[SRTBUFSZ];
  char curve_label[SRTBUFSZ]; /* String corresponding to the type */
  SrtCurveType curve_type;
  long curve_ticker;
  char curve_ccy[SRTBUFSZ];

  void *curve_desc;
} SrtCurveDesc, SrtCurveObj, *SrtCurvePtr, *SrtCrvPtr, *SrtCrvObjPtr;

/* ---------------------------------------------------------------------------------
 */

/* ------------------------------- YIELD CURVE
 * ------------------------------------- */

typedef struct {
  SrtYieldCurveType yc_type;
  SrtCcyParam *ccy_param;
  Date clcn_date;
  Date spot_date;

} SrtYCDesc;

/* ------------------------------- FORWARD CURVE
 * ----------------------------------- */

typedef struct {
  SrtDvdObj *dvd_obj;
  SrtRepoObj *repo_obj;
  Date clcn_date;
} SrtFwdDesc;

/* ---------------------------------- CMT CURVE
 * ------------------------------------ */

typedef struct {
  char fwdT_name[SRTBUFSZ];
  CMT_Obj *cmt;
  FwdT_Obj *fwdT;
  CMT_Param *cmtprm;
  Date clcn_date;
  char yc_name[SRTBUFSZ];
} SrtCMTDesc;

/* ====================================================================== */

#endif
