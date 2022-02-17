
/***********

        swp_h_cmt_types.h

************/

#ifndef SWP_H_CMT_TYPES_H
#define SWP_H_CMT_TYPES_H

#include "swp_h_cmt_param.h"
/*

typedef struct
{
        int 		  Dwidth   ;
        int 		 *Dlengths ;
        double 		**Dfields  ;

        int 		  Iwidth   ;
        int 		 *Ilengths ;
        long 		**Ifields  ;
}
        Field_List ;
*/

/* field list based types */

/******** Fwd_T_obj --- hold  forward treasury curve ****/

enum FwdT_DFields {
  FwdT_DATE,
  FwdT_TIME,
  FwdT_RATE,
  FwdT_SPREAD,
  FwdT_SWAP_RATE,
  FwdT_DWIDTH
};
enum FwdT_IFields { FwdT_IWIDTH };

typedef struct {
  double FwdT_TODAY;
  double FwdT_SPOT;
  CMTCode FwdT_CMTCODE;
  int FwdT_INTERP_METHOD;
  int FwdT_NUMSPREAD;
  Field_List *data;
} FwdT_Obj;

#define FwdT_field_Dget(fwdTptr, FIELDNAME, index)                             \
  field_Dget(fwdTptr->data, FIELDNAME, index)

#define FwdT_field_Iget(fwdTptr, FIELDNAME, index)                             \
  field_Iget(fwdTptr->data, FIELDNAME, index)

#define FwdT_field_Dateget(fwdTptr, FIELDNAME, index)                          \
  DTOL(field_Dget(zcptr->data, FIELDNAME, index))

#define FwdT_field_Dateset(fwdTptr, FIELDNAME, index, val)                     \
  { field_Dset(fwdTptr->data, FIELDNAME, index, (double)val); }

#define FwdT_Dateget(fwdTptr, FIELDNAME) DTOL((fwdTptr->FIELDNAME))

#define FwdT_Dateset(zcptr, FIELDNAME, val)                                    \
  { (fwdTptr->FIELDNAME) = (double)val; }

#define FwdT_field_Dset(fwdTptr, FIELDNAME, index, val)                        \
  { field_Dset(fwdTptr->data, FIELDNAME, index, val); }

#define FwdT_field_Iset(fwdTptr, FIELDNAME, index, val)                        \
  { field_Dset(fwdTptr->data, FIELDNAME, index, val); }

#define FwdT_field_Dlength(fwdTptr, FIELDNAME)                                 \
  { field_Dlength(fwdTptr->data, FIELDNAME); }

#define FwdT_field_Ilength(fwdTptr, FIELDNAME)                                 \
  { field_Ilength(fwdTptr->data, FIELDNAME); }

#define FwdT_Dget(fwdTptr, FIELDNAME) (fwdTptr->FIELDNAME)

#define FwdT_Iget(fwdTptr, FIELDNAME) (fwdTptr->FIELDNAME)

#define FwdT_Dset(fwdTptr, FIELDNAME, val)                                     \
  { (fwdTptr->FIELDNAME) = (val); }

#define FwdT_Iset(fwdTptr, FIELDNAME, val)                                     \
  { (fwdTptr->FIELDNAME) = (val); }

#define FwdT_field_set_Dlength(fwdTptr, len)                                   \
  { Field_List_set_Dlength(fwdTptr->data, len); }

#define FwdT_field_set_Ilength(fwdTptr, len)                                   \
  { Field_List_set_Ilength(fwdTptr->data, len); }

/******** CMT_obj --- hold CMT-CMS spread curve and dates ****/

enum CMT_DFields {
  CMT_SWAPEND_DATE,
  CMT_SPREAD_LIBOR,
  CMT_SWAP_FIRST_FULL_FIXING,
  CMS_SWAPEND_DATE,
  CMS_SPREAD_LIBOR,
  CMS_CMT_SPREAD_RATE,
  CMS_SWAP_FIRST_FULL_FIXING,
  CMT_DWIDTH
};

enum CMT_IFields {
  CMT_SWAP_NUM_FULL_PERIOD,
  CMT_SWAP_FREQ,
  CMT_SWAP_BASIS_CODE,
  CMS_SWAP_NUM_FULL_PERIOD,
  CMS_SWAP_FREQ,
  CMS_SWAP_BASIS_CODE,

  CMT_IWIDTH
};

typedef struct {
  CMTCode CMT_CODE;
  Ddate CMT_TODAY;
  Ddate CMT_SPOT;
  int CMT_SPOT_EXISTS;
  double CMT_SPOT_SPREAD;
  Ddate CMT_SPOT_DATE;

  int CMT_NUMSWAP;
  int CMT_NUMCASH;

  Field_List *data;
} CMT_Obj;

#define CMT_field_Dget(cmtptr, FIELDNAME, index)                               \
  field_Dget(cmtptr->data, FIELDNAME, index)
#define CMT_field_Iget(cmtptr, FIELDNAME, index)                               \
  field_Iget(cmtptr->data, FIELDNAME, index)

#define CMT_field_Dset(cmtptr, FIELDNAME, index, val)                          \
  { field_Dset(cmtptr->data, FIELDNAME, index, val); }
#define CMT_field_Iset(cmtptr, FIELDNAME, index, val)                          \
  { field_Iset(cmtptr->data, FIELDNAME, index, val); }

#define CMT_field_Dateget(cmtptr, FIELDNAME, index)                            \
  DTOL(field_Dget(cmtptr->data, FIELDNAME, index))
#define CMT_field_Dateset(cmtptr, FIELDNAME, index, val)                       \
  { field_Dset(cmtptr->data, FIELDNAME, index, (double)val); }

#define CMT_Dateget(cmtptr, FIELDNAME) DTOL((cmtptr->FIELDNAME))
#define CMT_Dateset(cmtptr, FIELDNAME, val)                                    \
  { (cmtptr->FIELDNAME) = (double)val; }

#define CMT_field_Dlength(cmtptr, FIELDNAME)                                   \
  field_Dlength(cmtptr->data, FIELDNAME)
#define CMT_field_Ilength(cmtptr, FIELDNAME)                                   \
  field_Ilength(cmtptr->data, FIELDNAME)

#define CMT_Dget(cmtptr, FIELDNAME) (cmtptr->FIELDNAME)
#define CMT_Iget(cmtptr, FIELDNAME) (cmtptr->FIELDNAME)

#define CMT_Dset(cmtptr, FIELDNAME, val)                                       \
  { (cmtptr->FIELDNAME) = (val); }
#define CMT_Iset(cmtptr, FIELDNAME, val)                                       \
  { (cmtptr->FIELDNAME) = (val); }

#define CMT_field_set_Dlength(cmtptr, len)                                     \
  Field_List_set_Dlength(cmtptr->data, len)
#define CMT_field_set_Ilength(cmtptr, len)                                     \
  Field_List_set_Ilength(cmtptr->data, len)

/*********** 	allocation	*******************/

FwdT_Obj *new_FwdT_Obj();
CMT_Obj *new_CMT_Obj();

void free_CMT_Obj(CMT_Obj *cmto);
void free_FwdT_Obj(FwdT_Obj *fwdTo);

#endif
