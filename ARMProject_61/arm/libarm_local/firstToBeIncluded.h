/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: firsttobeincluded.h,v $
 * Revision 1.1  2004/03/23 07:51:19  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file firsttobeincluded.h
 *
 *  \brief file to remove certain warnings!
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _FIRSTOBEINCLUDED_H
#define _FIRSTOBEINCLUDED_H

/// File to be included EVERYWHERE
/// disable STL warnings in Windows NT
/// to get back STL warning use the flag _SHOW_STL_WARNING
#ifdef WIN32
	#ifndef _SHOW_STL_WARNING
		#pragma warning(disable : 4786)
		#pragma warning(disable : 4273)
		#pragma warning(disable : 4275)
		#pragma warning(disable : 4503)
		#pragma warning(disable : 4005)
//		#pragma warning(disable : 4530)
	#endif
#endif

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
