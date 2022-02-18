/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins                */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SRT_H_VEGATRELGM.H                                    */
/*                                                                            */
/*      PURPOSE:        Functions to compute lgm tree                         */
/*                                                                            */
/*      AUTHORS:        G.Amblard, E.Auld, K.Chau, J.Malhi, A. Sahuguet       */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     Jasbir S Malhi                                        */
/*                                                                            */
/*      DATE:           March 1995                                            */
/*                                                                            */
/*      REASON:         add greek functions                                   */
/*                                                                            */
/******************************************************************************/

#ifndef VEGATRELGM_H
#define VEGATRELGM_H

/* -------------------------------------------------------------------------- */

/* -> customized version of srt_f_cheinitstp().
This function initializes each step with the values for sigma, tau...	    */

Err srt_f_vegashiftirministp(
    SrtStpPtr   stp,
    SrtUndPtr   und,
    SrtUndInfo* und_info, /* und info */
    double      sigma_bucket_start,
    double      sigma_bucket_end,
    double      sigma_shift,
    int         sigma_shift_type,
    double      tau_bucket_start,
    double      tau_bucket_end,
    double      tau_shift,
    int         tau_shift_type);

/* -------------------------------------------------------------------------- */

/* -> Customized version of function srt_f_srt_f_trelgm1d().
This function fills answers the request through 1 sweep of an Hull & White
tree.									*/

Err srt_f_vegashifttrelgm1d(
    SrtUndPtr     und,
    SrtGrfnParam* grfparam,
    SrtStpPtr     stp,
    GrfnDeal*     gd,
    EvalEventFct  evalcf,
    SrtIOStruct*  iolist,
    SrtUndInfo*   und_info);

#endif
