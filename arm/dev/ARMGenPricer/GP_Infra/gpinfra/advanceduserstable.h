/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: advanceduserstable.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file advanceduserstable.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPINFRA_ADVANCEDUSERSTABLE_H
#define _INGPINFRA_ADVANCEDUSERSTABLE_H

#include "gpbase/port.h"
#include <string>
CC_USING_NS(std,string)


CC_BEGIN_NAMESPACE( ARM )

extern const string GP_User_Level2_Table[];
extern const size_t GP_User_Level2_Table_Size;

bool UserControl_IsInTheList( const string list[], size_t size, const string& elem );


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

