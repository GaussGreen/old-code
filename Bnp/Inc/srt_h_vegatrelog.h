#ifndef SRT_H_VEGATRELOG_H
#define SRT_H_VEHATRELOG_H

Err srt_f_logtree(SrtUndPtr und, SrtGrfnParam *grfnparam, SrtStpPtr stp,
                  GrfnDeal *gd, EvalEventFct evalcf, SrtIOStruct *iolist,
                  SrtUndInfo *und_info);

void srt_f_trintrexindex(SrtStpPtr stp,  /* current step */
                         double *prev_x, /* state vars at previous step */
                         SrtTrinTreNdInf *node,   /* current node in tree */
                         SrtGrfnParam *grfnparam, /* model parameter */
                         int closestmeth /* CLOSESTINX or CLOSESTINLNX */
);

#endif
