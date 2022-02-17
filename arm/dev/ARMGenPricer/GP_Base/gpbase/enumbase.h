/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: enumbase.h,v $
 * Revision 1.1  2005/11/21 14:51:19  emezzine
 * Initial revision
 *
 *
 */

/*! \file enumbase.h
 *
 *  \brief 
 *
 *	\author  E.M Ezzine 
 *	\version 1.0
 *	\date November 2005
 */


#ifndef _INGPBASE_ENUMBASE_H
#define _INGPBASE_ENUMBASE_H

#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_InterpolationType
{
	enum InterpolationType
		{
            linear_column_extrapoleCst = 0,
		    linear_row_extrapoleCst,
		    linear_column_row_extrapoleCst_column_row,
		    linear_row_column_extrapoleCst_row_column,
            stepup_right_column,
		    stepup_right_row,
			stepup_left_column,
			stepup_left_row,
            unknown,
		};
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/