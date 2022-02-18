/* =================================================================================

   FILENAME: srt_h_lgmdynamics.h

   PURPOSE:  Give all the Equations required to descretise the LGM model:
                                - drift of r1 and r2
                                - expectations of r1 and r2 at t+1 knowing t
                                - local volatilities of r1 and r2
                                - varaince of r1 and r2 at t+1 knowing t
                                (and same for phi's and cross-phi)
   ================================================================================= */

#ifndef SRT_H_LGMDYNAMICS_H
#define SRT_H_LGMYDYNAMICS_H

void srt_f_lgmdrfatsam(SrtStpPtr stp, SrtSample* cur_sam, SrtSample* drift_sam, int index);

#endif