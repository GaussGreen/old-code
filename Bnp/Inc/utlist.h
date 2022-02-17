/* =============================================================================
   
   FILENAME:   	utlist.h

   PURPOSE:     The famous SrtList, that is used everywhere
   
   ============================================================================= */


#ifndef UTLIST_H
#define UTLIST_H

#include "utliststruct.h"
#include "utallhdr.h"


/* -------------------------- Function Prototypes --------------------------- */
/* -------------------------------------------------------------------------- */

/* ------------------------------------------------------------
            TOP LEVEL FUNCTIONS TO WORK ON THE LIST
   ------------------------------------------------------------ */

/* -> return as *PTr a SrtList* structure properly allocated in memory 
			 or NULL*/
Err srt_f_lstcreate(SrtList** Ptr, char* name);


/* -------------------------------------------------------------------------- */

/*-> put object 'obj' into list 'l'*/
Err srt_f_lstset( SrtList *l, SrtObject* obj);


/* ------------------------------------------------------------------------- */

/* -> insert an element in the list with name, key
Its value type is defined by type and Ptr points to its value */

Err srt_f_lstins(
		SrtList*  l,
		char*     name,
		double    key,
		int       type,
		void*     Ptr,
		Err        (*ObjValPvalFreeFct)(void *),
		long*     ticker);


/* -------------------------------------------------------------------------- */

/* -> (*ret) is ptr to SrtObject with name 'name' inside SrtList 'l'
	return NULL if object is found and *ret points to the object
	return sthg else otherwise and *ret = NULL
*/
Err srt_f_lstgetobj(SrtList l, char* name, double key, SrtListObject** ret);

/* -------------------------------------------------------------------------- */
/* Remove an element in the list: free it and re-establish connectivity  */
Err srt_f_lstremobj( SrtList *l, char* objname, double objkey);

/* ------------------------------------------------------------------------- */

/* -> return 1 if object with name 'name' or key number 'key' is found */
SRT_Boolean srt_f_lsthas(SrtList l, char name[], double key);


/* ------------------------------------------------------------------------- */

/* -> srt_free memory used by the list 'l' structure
	the objects of the list are not destroyed if flag is off (0)
	 only the ptr that refers to them
	the objects are destroyed if flag is on (1)
	*/
Err srt_f_lstfree(SrtList *l, int srt_free_obj_flag);


/* -------------------------------------------------------------------------- */

/* ------------------------------------------------------------
          FUNCTIONS TO WORK WITH AN OBJECT OF THE LIST
   ------------------------------------------------------------ */

/* -> return as *Ptr an SrtObject* properly allocated in memory or NULL */
Err srt_f_objset(
			char*        name,
			double       key,
			int          type,
			SrtObjVal*   val,
			SrtObject**  Ptr,
			SRT_Boolean      Ptr_allocated,
			Err         (*ObjPvalFreeFct)(void *));


/* -------------------------------------------------------------------------- */

/* -> return 0 if both objects are equal : same name, same key number */
SRT_Boolean srt_f_objequal(SrtObject* obj1, SrtObject* obj2);


/* -------------------------------------------------------------------------- */

/* -> return the order with which both objects will be put in the list */
int srt_f_objcmp(SrtObject* obj1, SrtObject* obj2);


/* -------------------------------------------------------------------------- */


/* -> srt_free memory used for the object obj */
Err srt_f_objfree(SrtObject* obj);


/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------
          FUNCTIONS TO CREATE AN OBJECTVAL FOR AN OBJECT
   ------------------------------------------------------------ */

/* return a pointer to an object with field pval = ptr */
SrtObjVal* val2obj(void* ptr);

/* ------------------------------------------------------------------------- */

/* ->  return a pointer to an object with field ival = i */
SrtObjVal* i2obj(int i);


/* ------------------------------------------------------------------------- */

/* ->  return a pointer to an object with field dva = d */
SrtObjVal* d2obj(double d);


/*--------------------------------------------------------------------- */

 /* ->  return a pointer to an object with field sval = s
The string s is physically copied into sval */
SrtObjVal* s2obj(char s[]);


/* -------------------------------------------------------------------
SrtObjVal* histdata2obj(SrtHistData* shd);

 ->  return a pointer to an object with field pval = shd 
---------------------------------------------------------------------- */


#endif
