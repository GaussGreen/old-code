/* -------------------------------------------------------------------------

   File: CCfdecl.h
   Path: /home/julia/projets/dev/Common/libcommon/SCCS/s.CCfdecl.h
   Description: macro de definition de fonctions
   Created: 
   Author: Teknekron Software Systems, Inc
   Modified: 99/07/19 18:31:09
   Last maintained by: Jacques WERNERT
   Revision: 1.8

   -------------------------------------------------------------------------

   Note: Copyright 1991, 1992 by Teknekron Software Systems, Inc.
         ALL RIGHTS RESERVED

   ------------------------------------------------------------------------- */

#ifndef _fdecl_h
#define _fdecl_h

#include <CCcommon.h>
SCCS_ID (fdecl_h_SccsId, "@(#)CCfdecl.h	1.8, modified 99/07/19");

/*
 *
 * MACRO		USAGE
 * ---------------	------------------------------------------------
 *
 * CCEXTERN_FUNCTION	Used to declare external functions.
 *
 * CCFUNCTION_TYPE	Used in specifying pointers to functions in
 *			typedefs or casts.
 *
 * CCUNKNOWN_ARGS	Used to specify an unknown argument list.
 *
 *
 * Examples:
 *
 *	CCEXTERN_FUNCTION(int fread, (char *, int, int, FILE *));
 *	CCEXTERN_FUNCTION(void printf, (char *, ...));
 *	CCEXTERN_FUNCTION(int getpid, (void));
 *	CCEXTERN_FUNCTION(void fun_with_unknown_args, (CCUNKNOWN_ARGS));
 *
 *	typedef CCFUNCTION_TYPE(int (*reader_t), (char *, int, int, FILE *));
 *	reader = (CCFUNCTION_TYPE(int (*), (char *, int, int, FILE *))) fread;
 *	reader = (CCFUNCTION_TYPE(int (*), (CCUNKNOWN_ARGS))) fread;
 *
 */


#if defined(__cplusplus)
/*
 * C++ (2.0+)
 */
# define CCEXTERN_FUNCTION(rtn, args) extern "C" { rtn args; }
# define CCSTATIC_FUNCTION(rtn, args) static rtn args
# define CCFUNCTION_TYPE(rtn, args) rtn args
# define CCUNKNOWN_ARGS ...

#elif defined(__STDC__) || defined(_STDC_PROTO_)
/*
 * ANSI C
 */
# define CCEXTERN_FUNCTION(rtn, args) extern rtn args
# define CCSTATIC_FUNCTION(rtn, args) static rtn args
# define CCFUNCTION_TYPE(rtn, args) rtn args
# define CCUNKNOWN_ARGS

#else
/*
 * K&R C
 */
# define CCEXTERN_FUNCTION(rtn, args) extern rtn()
# define CCSTATIC_FUNCTION(rtn, args) static rtn()
# define CCFUNCTION_TYPE(rtn, args) rtn()
# define CCUNKNOWN_ARGS

#endif

#ifdef CC_COMPAT

#define EXTERN_FUNCTION CCEXTERN_FUNCTION
#define STATIC_FUNCTION CCSTATIC_FUNCTION
#define FUNCTION_TYPE 	CCFUNCTION_TYPE
#define UNKNOWN_ARGS 	CCUNKNOWN_ARGS

#endif

#endif 

/* EOF CCfdecl.h */
