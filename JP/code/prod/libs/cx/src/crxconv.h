/*
***************************************************************************
** FILENAME: crxconv.h
**
** Converts between CX and CRX data structures.
**
** $Header$
***************************************************************************
*/

#ifndef CX_CRXCONV_H
#define CX_CRXCONV_H

#include "cx.h"
#include <crxflow/include/crcrv.h>

KProtLeg_D* CxProtectionLegConvert
(CxTContingentLeg *cl,  /* (I) Contingent leg       */
 CxTRecoveryCurve *rc,  /* (I) Recovery curve           */
 TDateInterval    *ivl  /* (I) Integration interval */
);

KFeeLeg_D* CxFeeLegConvert
(CxTFeeLeg        *fl,  /* (I) Fee leg              */
 TDateInterval   *ivl  /* (I) Integration interval */
);

#endif

