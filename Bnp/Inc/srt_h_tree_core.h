/* ------------------------------------------------------------------------------
   FILENAME:      srt_h_tree_core.h

   PURPOSE:       THE Core function for a tree discretisation
   ------------------------------------------------------------------------------ */
#ifndef SRT_H_TREE_CORE_H
#define SRT_H_TREE_CORE_H

Err TreeCore(
    SrtUndPtr     und,
    SrtGrfnParam* grfnparam,
    SrtStpPtr     step,
    GrfnDeal*     gd,
    EvalEventFct  evalcf,
    void*         iolist, /* price and shifted prices (greeks)*/
    SrtUndInfo*   und_info);

#endif