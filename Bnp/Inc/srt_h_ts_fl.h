/* SRT_H_TS_FL.H */

#ifndef TS_FL_H
#define TS_FL_H

#include "srt_h_all.h"

enum TS_Fields { TIME, TAU, SIG, F, G, H, PSI };

/* -------------------------------------------------------------------------- */

enum Che_TS_DFields {
  CTS_TAU_TIME,
  CTS_TAU,
  CTS_SIG_TIME,
  CTS_SIG,
  CTS_F,
  CTS_G,
  CTS_H,
  CTS_PSI,
  CTS_J,
  CTS_DWIDTH
};

enum Che_TS_IFields { CTS_SIG_DATES, CTS_TAU_DATES, CTS_IWIDTH };

/* -------------------------------------------------------------------------- */

TermStruct *fl2ts(Field_List *fl);

/* -> return a NEW TermStruct (2linked-list) from a Field_List */

/* -------------------------------------------------------------------------- */

#endif
