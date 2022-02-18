/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins                */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SRT_H_GREEK_TOOLS.H                                   */
/*                                                                            */
/*      PURPOSE:                                                              */
/*                                                                            */
/*      AUTHORS:        Arnaud Sahuguet                                       */
/*                                                                            */
/*      DATE:                                                                 */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     Jasbir S Malhi                                        */
/*                                                                            */
/*      DATE:           February 1995                                         */
/*                                                                            */
/*      REASON:         add greek functions                                   */
/*                                                                            */
/******************************************************************************/

#ifndef SRT_H_GREEK_TOOLS_H
#define SRT_H_GREEK_TOOLS_H

/* -------------------------------------------------------------------------- */

Err srt_f_trinfdup(SrtCheTreInf* trinf, SrtCheTreInf** trinf_dup);

/* -------------------------------------------------------------------------- */

SrtStpPtr srt_f_copystp(SrtStpPtr input_dest, SrtStpPtr input_source);

/* -------------------------------------------------------------------------- */

Err srt_f_stpdup(SrtStpPtr stp, SrtStpPtr* stp_dup);

/* -------------------------------------------------------------------------- */

void my_print_Step(SrtStpPtr s, char file_name[], int index);

#endif
