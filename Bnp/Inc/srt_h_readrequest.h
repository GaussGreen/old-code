/*********************************************************************
 *      Name: srt_h_readrequest.h                                    *
 *  Function: Prototypes for request list functions                  *
 * Copyright: (C) Paribas Capital Markets Ltd.                       *
 *-------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                     *
 *      Date: 03/11/95                                               *
 *-------------------------------------------------------------------*
 *    Inputs:                                                        *
 *   Returns:                                                        *
 *   Globals:                                                        *
 *-------------------------------------------------------------------*
 * Modification Record                                               *
 * Date     Inits   Comments                                         *
 * dd/mm/yy                                                          *
 * 03/11/95 FOS     Changed create_request_list to take no arguments *
 *-------------------------------------------------------------------*
 * Modification Record                                               *
 * Date     Inits   Comments                                         *
 * dd/mm/yy                                                          *
 * 15/03/96 DD     Changed create_request_list to take arguments     *
 *********************************************************************/

#ifndef SRT_H_READREQUEST_H
#define SRT_H_READREQUEST_H

/* Build a request list */

Err srt_f_readrequest(String *requests, int *buckets, int nr, TermStruct ts,
                      SrtIOStruct *iolist);

/* read the result of a request */

Err srt_f_readresult(String greek, int bucket, double *ret_val);

Err srt_f_sig_ts_dates(int, Date *, Date *, TermStruct);
Err srt_f_tau_ts_dates(int, Date *, Date *, TermStruct);

Err create_request_list(String request_list_name);
Err destroy_request_list();
SrtIOStruct *get_request_list(void);

#endif
