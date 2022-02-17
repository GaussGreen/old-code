/*
***************************************************************************
** platdep.h
**
** Compiler specific platform dependencies. This file only depends on
** on macros defined by the compiler vendors.
***************************************************************************
*/

#ifndef IRX_PLATDEP_H
#define IRX_PLATDEP_H

#ifdef __GNUC__

/* gcc */

#endif

#ifdef _MSC_VER

/* microsoft */
#define vsnprintf _vsnprintf
#define snprintf  _snprintf

#endif


#endif
