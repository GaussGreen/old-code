/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 * $Log: ostringstream.h,v $
 * Revision 1.1  2003/10/18 16:45:55  ebenhamou
 * Initial revision
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file ostringstream.h
 *
 *  \brief to handle ostringstream uncompatibility on unix
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

 
/*----------------------------------------------------------------------------*/

#ifndef _INGPBASE_OSTRINGSTREAM_H
#define _INGPBASE_OSTRINGSTREAM_H

#include "port.h"


#ifdef unix
	#include "myostringstream.h"
	typedef oStringStream CC_Ostringstream;

#else
	#include <sstream>
	typedef std::ostringstream CC_Ostringstream;
#endif


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
