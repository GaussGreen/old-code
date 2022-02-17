/* =============================================================================

   FILENAME      : swp_h_curve_fct.h

   PURPOSE       : functions to access curve information from the object stored
                   in the double linked list :
                   SrtCurve = SrtCurveDesc

   =============================================================================
 */

#ifndef SWP_H_CURVE_FCT_H
#define SWP_H_CURVE_FCT_H

/* -------------------------------------------------------------------------------
   Some Macros to access immediately curve information from the SrtCurve
   -------------------------------------------------------------------------------
 */

/* -------------------------------- Info from any curve
 * --------------------------- */

#define get_curve_name(crvptr) (crvptr->curve_name)

#define get_curve_type(crvptr) (crvptr->curve_type)

#define get_curve_ticker(crvptr) (crvptr->curve_ticker)

#define get_curve_ccy(crvptr) (crvptr->curve_ccy)

#define ISCURVETYPE(crvptr, type) (crvptr->curve_type == (type))

/* ------------------------------- Yield Curves
 * ------------------------------------- */

#define get_clcndate_from_yldcrv(crvptr)                                       \
  (((SrtYCDesc *)(crvptr->curve_desc))->clcn_date)

#define get_ccyparam_from_yldcrv(crvptr)                                       \
  (((SrtYCDesc *)(crvptr->curve_desc))->ccy_param)

#define get_yctype_from_yldcrv(crvptr)                                         \
  (((SrtYCDesc *)(crvptr->curve_desc))->yc_type)

#define get_spotdate_from_yldcrv(crvptr)                                       \
  (((SrtYCDesc *)(crvptr->curve_desc))->spot_date)

#define get_discfunc_from_yldcrv(crvptr)                                       \
  (((SrtYCDesc *)(crvptr->curve_desc))->DiscFunc)

#define get_crvobj_from_yldcrv(crvptr)                                         \
  (((SrtYCDesc *)(crvptr->curve_desc))->crv_object)

#define set_spotdate_from_yldcrv(crvptr, spotdate)                             \
  (((SrtYCDesc *)(crvptr->curve_desc))->spot_date) = spotdate

/* ------------------------------ Forward Curves
 * ------------------------------------ */

#define get_clcndate_from_dvdcrv(crvptr)                                       \
  (((SrtFwdDesc *)(crvptr->curve_desc))->clcn_date)

#define get_dvdobj_from_dvdcrv(crvptr)                                         \
  (((SrtFwdDesc *)(crvptr->curve_desc))->dvd_obj)

#define get_clcndate_from_repocrv(crvptr)                                      \
  (((SrtFwdDesc *)(crvptr->curve_desc))->clcn_date)

#define get_repoobj_from_repocrv(crvptr)                                       \
  (((SrtFwdDesc *)(crvptr->curve_desc))->repo_obj)

/* -------------------------------- Cmt Curves
 * -------------------------------------- */

#define get_clcndate_from_cmtcrv(crvptr)                                       \
  (((SrtCMTDesc *)(crvptr->curve_desc))->clcn_date)

#define get_ycname_from_cmtcrv(crvptr)                                         \
  (((SrtCMTDesc *)(crvptr->curve_desc))->yc_name)

#define get_fwdTname_from_cmtcrv(crvptr)                                       \
  (((SrtCMTDesc *)(crvptr->curve_desc))->fwdT_name)

#define get_cmtobj_from_cmtcrv(crvptr)                                         \
  (((SrtCMTDesc *)(crvptr->curve_desc))->cmt)

#define get_fwdTobj_from_cmtcrv(crvptr)                                        \
  (((SrtCMTDesc *)(crvptr->curve_desc))->fwdT)

#define get_cmtparam_from_cmtcrv(crvptr)                                       \
  (((SrtCMTDesc *)(crvptr->curve_desc))->cmtprm)

/* -------------------------------------------------------------------------------
   Some Functions to access immediately curve information from any SrtCurve
   -------------------------------------------------------------------------------
 */

/* Functions to get today == clcn_date */

Err get_curve_clcndate(SrtCurvePtr crv, long *clcn_date);

long get_clcndate_from_curve(SrtCurvePtr crv);

#define get_today_from_curve get_clcndate_from_curve

/* Functions to get the spot lag */

Err get_curve_spotlag(SrtCurvePtr crv, int *spotlag);

int get_spotlag_from_curve(SrtCurvePtr crv);

/* Functions to get the name of the yield curve asociated to a curve */

Err get_curve_ycname(SrtCurvePtr crv, String *yc_name);

char *get_ycname_from_curve(SrtCurvePtr crv);

#endif
