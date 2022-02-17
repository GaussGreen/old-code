/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT        , Fixed Income 2020 Addins */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SRT_F_REQUEST_LIST                                    */
/*                                                                            */
/*      PURPOSE:        To provide a home for greek requests from 2020        */
/*                                                                            */
/*      AUTHORS:        Jasbir MALHI                  			      */
/*                                                                            */
/*      DATE:           22nd February 1995                                    */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:                                                           */
/*                                                                            */
/*      DATE:                                                                 */
/*                                                                            */
/*                                                                            */
/*      REASON:                                                               */
/*                                                                            */
/*                                                                            */
/******************************************************************************/

/**  Include Statements  ******************************************************/

#include "srt_h_all.h"

/**  Static Declarations  *****************************************************/

static SrtIOStruct *__request_list;

/* -------------------------------------------------------------------------- */

Err create_request_list(String request_list_name) {
  Err err;

  if (err = srt_f_IOstructcreate(&__request_list, request_list_name))
    return err;

  return NULL;
}

Err destroy_request_list() {
  Err err;

  if (err = srt_f_IOstructfree(&__request_list)) {
    return err;
  }

  return NULL;
}

SrtIOStruct *get_request_list(void) { return __request_list; }

/* ========================================================================== */
