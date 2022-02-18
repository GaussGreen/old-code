/* =================================================================================

   FILENAME: srt_h_cheytreefct.h

   PURPOSE:  Useful functions to build a One Factor Cheyette Tree

   ================================================================================= */

#ifndef SRT_H_CHEYTREEFCT_H
#define SRT_H_CHEYTREEFCT_H

void srt_f_chetrerindex(
    SrtStpPtr stp, double* prev_r, SrtTrinTreNdInf* node, SrtGrfnParam* grfnparam, int centering);

void srt_f_chetreprob(
    SrtStpPtr stp, double* prev_r, SrtTrinTreNdInf* node, SrtGrfnParam* grfnparam);

void srt_f_chetredcntvec(
    SrtStpPtr        stp,
    double*          cur,
    SrtTrinTreNdInf* node,
    double***        assets,
    int              node_dim,
    double**         prev_phi);

Err srt_f_cheytreelim(SrtStpPtr stp, SrtGrfnParam* grfnparam, SrtCheTreInf* maxtrinf);

#endif