/* =================================================================================

   FILENAME: srt_h_cheydynamics.h

   PURPOSE:  Give all the Equations required to descretise the model:
                                - drift of r1 and r2
                                - expectations of r1 and r2 at t+1 knowing t
                                - local volatilities of r1 and r2
                                - varaince of r1 and r2 at t+1 knowing t
                                (and same for phi's and cross-phi)
   ================================================================================= */

#ifndef SRT_H_CHEYDYNAMICS_H
#define SRT_H_CHEYDYNAMICS_H

/*  =====================================================================
    These two function are defined in srt_f_chefct.c
    ====================================================================== */
void srt_f_chedrfatsam(SrtStpPtr stp, SrtSample* cur_sam, SrtSample* drift_sam, int index);

void srt_f_chevaratsam(SrtStpPtr stp, SrtSample* cur_sam, double* var_at_sam, int index);

#endif