/* ==========================================================================
   FILE_NAME:	srt_h_rtdesc.h

   PURPOSE:     function to deal with rates descriptions...
   ========================================================================== */
#ifndef SRT_H_RTDESC_H
#define SRT_H_RTDESC_H

/*
 * author E.Auld
 * a SrtRtFnc is a structure meant to represent a rate we might want to know
 */

/*-------------------------private structures------------------------*/
typedef enum {
  SRTRTFNCSWP,
  SRTRTFNCLVL,
  SRTRTFNCIMMSWP,
  SRTRTFNCIMMLVL,
  SRTRTFNCDF,
  SRTRTFNCFRA,
  LASTSRTRTFNCTYPE
} SrtRtFncType;

typedef struct _srt_rtfnc_desc {
  SrtRtFncType type;
  Ddate *date;
  double *cvg;
  void *yp; /* hook for hanging interest rate mdl info */
  double *t;
  double *df;
  Ddate evl_date;
  int len;
} SrtRtFncElement, *SrtRtFnc;

/* top level useful functions */
void srt_f_rtprint(FILE *out, SrtRtFnc rt);
Err srt_f_rtevl(SrtRtFnc rt, double short_rate, double *answer);
void srt_f_rtfre(SrtRtFnc rt);
Err srt_f_rtmk(char *rtname, Date start, char *stastr, Date end_or_nfp,
               char *endstr, char *freq, char *basis, SrtRtFnc *ret);

/* lower level for getting at what's inside */
void srt_f_rtsetyp(SrtRtFnc rt, void *ptr);
void *srt_f_rtgetyp(SrtRtFnc rt);
void srt_f_rtsettm(SrtRtFnc rt, double *tm);
double *srt_f_rtgettm(SrtRtFnc rt);
void srt_f_rtsetdf(SrtRtFnc rt, double *df);
double *srt_f_rtgetdf(SrtRtFnc rt);
void srt_f_rtsetdt(SrtRtFnc rt, Ddate *dt);
Ddate *srt_f_rtgetdt(SrtRtFnc rt);
void srt_f_rtsetcvg(SrtRtFnc rt, Ddate *cvg);
Ddate *srt_f_rtgetcvg(SrtRtFnc cvg);
int srt_f_rtgetlen(SrtRtFnc rt);

#endif
