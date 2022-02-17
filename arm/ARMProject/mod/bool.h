/*
 * $Log: bool.h,v $
 * Revision 1.4  2000/10/25 14:20:12  smysona
 * Ajout d'un test STL
 *
 * Revision 1.3  1999/02/15 11:05:33  ypilchen
 * Chgt de UNIX en unix
 *
 * Revision 1.2  1999/02/15 10:28:18  ypilchen
 * Rajout de : ifdef UNIX
 *
 * Revision 1.1  1998/11/17 14:53:48  ypilchen
 * Initial revision
 *
 */

#ifndef _bool_h
#define _bool_h 1

#ifndef _ARM_STL_
#ifdef unix
 
#undef FALSE
#undef TRUE
#undef true
#undef false
 
enum bool { FALSE = 0, false = 0, TRUE = 1, true = 1 };

#endif

 
#endif // unix
#endif // _ARM_STL_

