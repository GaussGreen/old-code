/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 *  \file badalloc.h
 *  \brief Macro Switch on/off the bad allocation handler
 * 
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date December 2005
 */

#ifndef _INGPBASE_BADALLOC_H
#define _INGPBASE_BADALLOC_H


class ARM_BadAlloc
{
public:
// FIXMEFRED: mig.vc8 (22/05/2007 15:58:43):missing return type
	static void BadAllocHandlerWithExp();

	// Switch on the bad allocation handler with exception
	static void SwitchOn();
	// Switch off the bad allocation handler with exception
	static void SwitchOff();
};
	
#endif ///_INGPBASE_BADALLOC_H