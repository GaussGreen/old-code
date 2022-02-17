/* =================================================================================

   FILENAME: srt_h_cheybetadynamics.h

   PURPOSE:  Give all the functions required to descretise the Cheyette Beta
             model (based on the diffusion equation):
                                - drift of r
                                - expectation of r at t+1 knowing t
                                - local volatility of r
                                - varaince of r at t+1 knowing t (and same for
   phi)
   =================================================================================
 */

#ifndef SRT_H_CHEYDYNAMICS_H
#define SRT_H_CHEYDYNAMICS_H

/* ----------------------------------------------------------------------
   From the short rate  and the power  , returns the local (normal)
   volatility
   ---------------------------------------------------------------------- */

Err srt_f_CheyBeta_local_vol(SrtMdlType mdl_type, SrtIRMTmInf *tminf,
                             SrtSample *cur_sample, int und_index,
                             double *local_vol);

Err srt_f_CheyBeta_drift_at_sam(SrtMdlType mdl_type, SrtStpPtr stp,
                                SrtSample *cur_sam, SrtSample *drift_sam,
                                int und_index);
Err srt_f_CheyBeta_var_at_sam(SrtMdlType mdl_type, SrtStpPtr stp,
                              SrtSample *cur_sam, double *var_at_sam,
                              int und_index);
#endif