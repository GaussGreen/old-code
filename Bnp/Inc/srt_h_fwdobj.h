/* =======================================================================================

   FILENAME :       srt_h_fwdobj.h

   PURPOSE:         structures and functions to store a forward curve for use in
   Grfn A ForwardCurve is stored as a SrtDvdObj which is a SrtListHdr

   =======================================================================================
 */

#ifndef SRT_H_FWDOBJ_H
#define SRT_H_FWDOBJ_H

/* --------------------------------------------------------------------------------------
   THE FORWARD CURVE TYPE: A Forward Curve is the Header of a SrtList...
   --------------------------------------------------------------------------------------
 */
typedef SrtListHdr SrtDvdObj;

/* ------------------------------------------------------------------------------------
   FUNCTIONS TO CREATE      , INSERT AN ELEMENT IN OR DELETE A SrtDvdObj ==
   SrtList
   ------------------------------------------------------------------------------------
 */

/* Creates a SrtdvdObject with only one  element containing today labelled
 * "_TODAY" */
Err srt_f_dvdobj_create(Date today, char *dvd_name, SrtDvdObj **dvdobj);

/* Delete a SrtdvdObject and everything it contains */
Err srt_f_dvdobj_delete(SrtDvdObj **dvdobj);

/* Insert a point (date      , value) in the forward object */
Err srt_f_dvdobj_insertpoint(SrtDvdObj *dvdobj, Ddate date, double dvd);

/* ----------------------------------------------------------------------------------
         INITIALISATION OF A SrtDvdObj FROM RAW DATA (DATES + VALUES )
   ----------------------------------------------------------------------------------
 */

/* Allocate memory and build a SrtDvdObj from a vector of forwards: date / value
 */
Err srt_f_dvdobj_init(Date today,
                      double **forwards, /* Forward curve [date][value] */
                      int num_cols, int num_rows, char *dvd_name,
                      SrtDvdObj **dvdobj);

/* ----------------------------------------------------------------------------------
   FUINCTIONS TO INTERPOLATE OR EXTRAPOLATE INFORMATION FROM A SrtDvdObj
   ----------------------------------------------------------------------------------
 */

/* Extracts today from the element labelled "_TODAY" in the SrtDvdObj */
Err srt_f_dvdobj_extracttoday(SrtDvdObj *dvdobj, Ddate *today);

Err srt_f_dvdobj_dvd(SrtDvdObj *dvdobj, double start_date, double *dvd);

/* Computes the drift (ratio of log of forwards) from start_date to end_date */
Err srt_f_dvdobj_drift(SrtDvdObj *dvdobj, double start_date, double end_date,
                       double *drift);

/* Computes the implied df (inverse of forward ratio) from start_date to
 * end_date */
Err srt_f_dvdobj_df(SrtDvdObj *dvdobj, double start_end, double end_date,
                    double *df);

#endif
