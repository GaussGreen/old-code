/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: undef_va_vars.h,v $
 * Revision 1.1  2004/03/23 07:51:19  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file undef_va_vars.h
 *
 *  \brief file to remove the warnings on studio 6.0 basically remove the macro va_start and va_end!
 * C:\Program Files\Microsoft Visual Studio\VC98\INCLUDE\varargs.h(70) : warning C4005: 'va_start' : macro redefinition
 *       C:\Program Files\Microsoft Visual Studio\VC98\INCLUDE\stdarg.h(58) : see previous definition of 'va_start'
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _REMOVE_VA_WARNINGS_H
#define _REMOVE_VA_WARNINGS_H


#ifdef WIN32
	#ifndef _REMOVE_VA_WARNINGS

		#if defined( va_start)
			#undef va_start
		#endif

		#if defined( va_end)
			#undef va_end
		#endif
		
	#endif
#endif

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
