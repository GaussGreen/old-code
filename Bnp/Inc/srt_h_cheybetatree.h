#ifndef SRT_H_CHEYBETATREE_H
#define SRT_H_CHEYBETATREE_H

Err srt_f_CheyBetaTree(SrtUndPtr und,           /* Underlying pointer */
                       SrtGrfnParam *grfnparam, /* Model Parameters */
                       SrtStpPtr stp,           /* Step pointer  */
                       GrfnDeal *gd,            /* GRFN Deal description */
                       EvalEventFct evalcf,     /* Cashflow evalution function*/
                       SrtIOStruct *iolist,     /* List of requests */
                       SrtUndInfo *und_info);   /* Underlying info */

#endif