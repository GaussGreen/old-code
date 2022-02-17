#ifndef _VNFMDATAO_H 
#define _VNFMDATAO_H

/* These are in ~clibdev/clutil. Lib is libclutil.a/clutil.lib */
#include "varapi.h"            /* TVar */ 
#include "objman.h"            /* TObjectType, TObjectManager */
#include "objapi.h"            /* TEncBuffer, TDecBuffer */

#include "bastypes.h"             /* TMatrix2D ... */ 

#include "drlmem.h"
#include "drlvtype.h"                   /* DrlLilStructGet() */
#include "drlsmat.h"
#include "drlts.h"                      /* DrlTCurveWrap() */
#include "drlio.h"                      /* DrlFPrintf() */

#include "vnfmanly.h"
#include "vnfmcali.h"
 
extern  TObjectType VNFMDATA_OBJECT_TYPE; 

/*------------------------------------------------------------------- 
 *                  PROTOTYPES FOR THE METHODS 
 *------------------------------------------------------------------- 
 */

GTO_EXPORT(int)         VnfmObjectRegister ( 
        TObjectManager * om); 

GTO_EXPORT(VnfmData*)  VnfmDataCreate (
	TDateL *refDateL,	/* 'D' (I) reference date */
	TDateL *zcDatesL,       /* 'D' (I) array of zero coupon dates */
        FloatL *zcRatesL,       /* 'F' (I) array of zero coupon rates */	

	FloatL *backboneqL,     /* 'F' (I) back bone Q */
        FloatL *betaL,          /* 'F' (I) array of mr coeff */
        FloatL *alphaL,         /* 'F' (I) array of weight coeff */
        TDateL *dateL,          /* 'D' (I) array of dates */
        FloatL *sigmaL,         /* 'F' (I) volatility arrays */
        FloatL *rhoL);          /* 'F' (I) correlation arrays */



GTO_EXPORT(VnfmData *) VnfmDataMakeCopy (
        VnfmData  *that);


GTO_EXPORT(void)        VnfmDataDelete (
        VnfmData     *that); 

GTO_EXPORT(void)         VnfmDataPrint (
        VnfmData     *that);


GTO_EXPORT(TVar *)      VnfmDataGetField (
        VnfmData     *that, 
        char         *filedName); 


GTO_EXPORT(int)         VnfmDataEncode (
        TEncBuffer *eb,
        char       *objName,
        VnfmData   *that);

GTO_EXPORT(VnfmData*) VnfmDataDecode (
        TDecBuffer *db,
        char *objName);



#endif /* _VNFMDATAO_H */
