#ifndef SRT_H_BNDTYPES_H
#define SRT_H_BNDTYPES_H
/*---------------------------------------------------------------------------
  AMENDMENTS      :
  Reference       :
  Author          :O. Van Eyseren from previous code E. Auld
  Date            :24 Apr 1995
  Description     :for Bond People (do not try to understand  , I do not even
                        known what this is all about...)
                   this file has to be read with srt_h_mdl_types.h  ,
                   from which it is inspired
                   moreover  , tihs is a complement to srt_h_und.h where the
                   other SrtUndDesc are described
-----------------------------------------------------------------------------*/

/**** Bond Auxiliary Parameters ***/

typedef struct {
  double first_coupon;
  Date first_coupon_date;
  double gross_coupon;
  double redemption;
  double face_value;
  int ex_dividend_lag;
  SRT_Boolean dbg;
  double **smile; /* Just for Francois PICARD */
} SrtBndAuxParam;

/*
Err srt_f_bndauxinit(SrtBndPtr bnd);
Err srt_f_bndauxtype(String aux_name  , SrtStrEleType *type);
Err srt_f_bndauxset(SrtBndPtr bnd  , SrtStrEle *ele);
*/
#endif
