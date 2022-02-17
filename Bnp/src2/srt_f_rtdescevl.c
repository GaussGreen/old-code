#include "srt_h_all.h"
#include "srt_h_rtdescevl.h"

/*<%%STA>-----------------------------------------------------------------
  MODULE          : srt_f_rtdescevl.c
  AUTHOR          : E.Auld
  DESCRIPTION     : functions to populate SrtRtFnc objects with market info
  DATE            : April 1995

<%%END>---------------------------------------------------------------------*/
/*---------------------------------------------------------------------------
  AMENDMENTS      :
  Reference       :
  Author          :
  Date            :
  Description     :
-----------------------------------------------------------------------------*/

/*
        evaluate a SrtRtFnc from todays yield curve
*/
SrtErr srt_f_rtevlfwd(SrtRtFnc rt, String ycname, double *answer) {
  int i, len;
  double *df;
  SrtCurvePtr yldcrv = lookup_curve(ycname);
  Ddate *dt, today = get_clcndate_from_yldcrv(yldcrv);
  yldcrv = lookup_curve(ycname);
  df = srt_f_rtgetdf(rt);
  len = srt_f_rtgetlen(rt);
  dt = srt_f_rtgetdt(rt);

  /*
   * if this SrtRtFnc has never been evaluated  , it may
   * not have had its df field allocated yet
   */
  if (!df) {
    df = srt_calloc(len, sizeof(double));
  }
  for (i = 0; i < len; i++) {
    df[i] = swp_f_df(today, dt[i], ycname);
  }
  srt_f_rtsetdf(rt, df);
  return srt_f_rtevl(rt, 0.05, answer);
}
