/*********************************************************************************
 * IRXUTILIO.H 
 * crx io utils
 *
 ********************************************************************************/
#ifndef _IRX_IRXUTILIO_H_
#define _IRX_IRXUTILIO_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "date.h"
#include "error.h"
#include "irxutils.h"

#define        READ_DATA(type,ptr,str)        \
    { if (irxFScanVType(fp, type, (void*) ptr) != SUCCESS) \
          { irxError("%s: can't read %s.\n", routine, str); \
                           goto RETURN;}}

/* String size */
#define  WRAP_STR_BYTES  128

/* variable type specification */
typedef        long        IrxTVType;
                                       /* Data types */
#define        IRX_NULL_T              ((IrxTVType) 0)
                                                /* C types */
#define        IRX_POINTER_T           ((IrxTVType) 02)

#define        IRX_DOUBLE_T            ((IrxTVType) 10)
#define        IRX_FLOAT_T             ((IrxTVType) 11)
#define        IRX_INT_T               ((IrxTVType) 12)
#define        IRX_LONG_T              ((IrxTVType) 13)
#define        IRX_CHAR_T              ((IrxTVType) 14)        /* single char */
#define        IRX_STRING_T            ((IrxTVType) 20)        /* char pointer (ie char *p) */
#define        IRX_CHAR_ARRAY_T        ((IrxTVType) 21)        /* char array (ie char p[..] */
#define        IRX_TDATE_T             ((IrxTVType) 30)
#define        IRX_TDATEINTERVAL_T     ((IrxTVType) 31)
#define        IRX_TDAYCOUNT_T         ((IrxTVType) 32)
                                        /* Types build form basic types: */
#define        IRX_PERCENT_T           ((IrxTVType) 33)        /* double multilpied by 100 in I/Os */
#define        IRX_CUR_T               ((IrxTVType) 34)        /* currency format I/Os */
#define        IRX_CURK_T              ((IrxTVType) 35)        /* currency format in K I/Os */
#define        IRX_CURM_T              ((IrxTVType) 36)        /* currency format in M I/Os */
#define        IRX_BOOLEAN_T           ((IrxTVType) 37)        /* boolean (TRUE,FALSE) as int */
#define        IRX_CHAR_BLOCK_T        ((IrxTVType) 84)        
#define        IRX_BPOINT_T            ((IrxTVType) 92)        /* IRX_FLOAT_L, but in bp  */

int            irxFScanVType(FILE *fp, IrxTVType type, void *valPtr);
int            irxFPrintVType(FILE *fp, IrxTVType type, void *valPtr);

int            irxVTypeScan(char *str, IrxTVType type, void *p);
char*          irxVTypePrint(char *str, IrxTVType type, void *p);

int            irxVectArrayFpReadV(
                FILE *fp,        /* (I) file pointer */
                int numItems,    /* (I) num elements to be read */
                /* IrxTVType varType, void *lptr,
                 * ...
                 * IrxTVType varType, void *lptr,
                 * IRX_NULL_T (last argument MUST be IRX_NULL_T)
                 */
                ...);


int            irxFScanString(FILE *fp, char *s);
/*f----------------------------------------------------------------------
 * Returns the name of a type "type".
 */
char*          irxVTypeName(IrxTVType type);

/*f----------------------------------------------------------------------
 * Returns SUCCESS if type "type" is a valid C type.
 */
int            irxVTypeCheckValidCType(IrxTVType type);

    
/*f----------------------------------------------------------------------
 * Returns the size of <i> type</i>.
 */
size_t         irxVTypeSizeof(IrxTVType type);


/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif
