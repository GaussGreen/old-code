/* ================================================================================
  
   FILENAME:      swp_f_string_interp.c

   PURPOSE:       string interpretation functions for swaps parameters

   ================================================================================ */

#ifndef SWP_H_STRING_INTERP_H
#define SWP_H_STRING_INTERP_H


Err interp_no_overwrite_flg(String cs, SRT_Boolean *b);
Err translate_no_overwrite_flg(String *cs, SRT_Boolean b);
Err interp_insert_flg(String cs, SRT_Boolean *b);
Err translate_insert_flg(String *cs, SRT_Boolean b);

Err interp_swap_message(String str, Message *mess) ;
Err interp_struct(String str, Message *val) ;

#endif
