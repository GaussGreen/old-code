/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file removeidentifiedwarning.h
 *
 *  \brief files to control warnings
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2004
 */


//////////////////////////////////////////
/// File to be included EVERYWHERE
/// disable STL warnings in Windows NT
/// to get back STL warning use the flag _SHOW_STL_WARNING
///////////////////////////////////////////

#ifdef WIN32
	#ifndef _SHOW_STL_WARNING

/// Compiler Warning (level 1) C4786
/// 'identifier' : identifier was truncated to 'number' characters in the debug information
		#pragma warning(disable : 4786) 

/// Compiler Warning (level 1) C4503
///	'identifier' : decorated name length exceeded, name was truncated
		#pragma warning(disable : 4503) 

/// Compiler Warning (level 2) C4275
/// non DLL-interface classkey 'identifier' used as base for DLL-interface classkey 'identifier' 
		#pragma warning(disable : 4275) 

	#endif
#endif


//-----------------------------------------------------------------------------
/*---- End of file ----*/
