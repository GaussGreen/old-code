/* =========================================================================

   FILENAME:    srt_h_und_list.h

   PURPOSE:     Work with Underlyings in lists

   ========================================================================= */

#ifndef SRT_H_UND_LIST_H
#define SRT_H_UND_LIST_H

/* -------------------------------------------------------------------------
   An underlying list is a usual double linked list
   ------------------------------------------------------------------------- */

typedef SrtList SrtUndList, *SrtUndListPtr;

/* -------------------------------------------------------------------------
   Functions to operate on the underlying list
   ------------------------------------------------------------------------- */

/* Allocate memory for the list where all the underlyings will be stored */
Err create_underlying_list(String und_list_name);

/* Destroy the list where all the underlyings are stored, as well as the underlyings */
Err destroy_all_underlyings();

/* An easy way to access the static pointer to the underlying list */
SrtUndListPtr get_underlying_list(void);

/* --------------------------------------------------------------------------
              Functions to operate with an underlying in a list
   -------------------------------------------------------------------------- */

/* Check if a list has go an underlying */
SRT_Boolean srt_f_isundinlist(SrtUndListPtr und_list, String und_name);

/* Get the Underlying (structure) corresponding to the name in a list */
Err srt_f_getundinlist(SrtUndListPtr und_list, String und_name, SrtUndPtr* und);

/* Add an underlying to a list (overwrite the previous one if exist ) */
Err srt_f_addundtolist(
    SrtUndListPtr und_list,
    String        und_name, /* == und_name in SrtUndDesc */
    String        und_lbl,  /* FX_UND, IR_UND, EQ_UND */
    String        und_ccy,
    String        mdl_lbl,
    String        crv_name1, /* == discount curve */
    String        crv_name2, /* == dividend curve */
    String        crv_name3, /* == repo curve */
    TermStruct*   ts,
    double        spot);

/* ---------------------------------------------------------------------------------
                      FUNCTION TO LOOKUP AN UNDERLYING OR DESTROY IT
   --------------------------------------------------------------------------------- */

/* Access an underlying through its name */
SrtUndPtr srt_f_lookup_und(String und_name);

#define lookup_und srt_f_lookup_und

/* Get the FX Underlying  corresponding to the domname and forname in the underlying list used */
SrtUndPtr srt_f_lookup_fxund(String dom_name, String for_name);

/* Destroys an underlying associated to a name */
Err srt_f_destroy_und(String undName);

/* --------------------------------------------------------------------------------- */
/* destructor for undelrying descriptor  (Stanley Mrose: 17.10.2002) */
Err srt_f_unddescfree(SrtUndDesc* spec_desc);

#endif
