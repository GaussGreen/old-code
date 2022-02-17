/* ===========================================================
   FILENAME:      utliststruct.h

   PURPOSE:       THE double linked list used for the storage
                  of any information
   =========================================================== */

/* If you want to add a new class of object in the list  ,
watch lines where '### How to add a new class ###' is present.
An easy example is given to help as comments.
*/

#ifndef UTLISTSTRUCT_H
#define UTLISTSTRUCT_H

#include "utDates.h"
#include "utListStruct.h"

/* --------------------------Object Type Codes ------------------------------ */

enum OBJECT_TYPE_CODE {
  /* Simple Object Types */
  OBJ_INT = 1,
  OBJ_DOUBLE = 2,
  OBJ_STRING = 3,
  OBJ_PTR = 10,

  /* Specific Object Types : TermStructures */
  OBJ_PTR_IRM_TermStruct = 20,
  OBJ_PTR_IRM2f_TermStruct = 21,
  OBJ_PTR_EQ_TermStruct = 22,
  OBJ_PTR_FX_TermStruct = 23,

  /* Specific Object Types : I/O list */
  OBJ_PTR_IO = 30,

  /* Specific Object Types : Curve (Yield  , Spread  , ...) */
  OBJ_PTR_CURVE = 40,

  /* Specific Object Types : Underlying (IRM  , EQ  , FX) */
  OBJ_PTR_UND = 41,

  /* Specific Object Types : Correlation list */
  OBJ_PTR_CorrelStruct = 50,

  /* Specific Object Types : Histograms list */
  OBJ_PTR_HIST = 60

};

/* ----------------------------- SrtObjVal ---------------------------------- */

/* value of an object */
typedef union SrtObjVal {
  int ival;
  double dval;
  char sval[32];
  void *pval;

} SrtObjVal;

/* ------------------------------ OBJECT ------------------------------------ */

/* object itself */
typedef struct SrtObject {
  /* The Value of the object */
  SrtObjVal val;

  /* The Name of the Object for lookups */
  char name[32];

  /* A Key attached to the object to spot its place in the list (time  ,...) */
  double key;

  /* The type of the pval object attached: i  , d  , s  , p  , which p... */
  int type;

  /* The VersionNumber of the object ("ticker") if several time the same one */
  int ticker;

  /* A flag (YES/NO) wether memory has been allocated for *val */
  SRT_Boolean ObjValPtr_allocated;

  /* The function to free the val.pval (if needed) */
  Err (*ObjValPvalFreeFct)(void *);

} SrtObject, SrtListObject;

/* -------------------------------- SrtLst ---------------------------------- */

/* list atom */
typedef struct SrtListAtom {
  struct SrtListAtom *next;
  struct SrtListAtom *previous;
  SrtListObject *element;

} SrtLst, SrtListAtom;

/* -------------------------------- SrtList --------------------------------- */

/* The Hook of the List */
typedef struct SrtList {
  SrtLst *tail;
  SrtLst *head;
  SrtLst *current;
  char name[32];

} SrtList, SrtListHdr;

#endif
