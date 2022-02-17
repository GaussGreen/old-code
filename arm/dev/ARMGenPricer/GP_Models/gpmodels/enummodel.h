/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: calibdirection.h,v $
 * Revision 1.1  2005/06/16 14:51:19  emezzine
 * Initial revision
 *
 *
 */

/*! \file enummodel.h
 *
 *  \brief 
 *
 *	\author  E.M Ezzine 
 *	\version 1.0
 *	\date June 2005
 */


#ifndef _INGPMODEL_ENUMMODEL_H
#define _INGPMODEL_ENUMMODEL_H

#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_PricingModelType
{
	enum ModelType
		{
            SFRM1F = 0,
		    SFRM2F,
		    QGM1F,
		    HWM1F,
            HWM2F,
		    QM,
			SBGM,
			HK,
            Unknown,
		};
};

struct ARM_MultiAssetsModelType
{
	enum MultiAssetsModelType
		{
            twoirfx = 0,
		    oneirfx,
		    irfx,
		    unknown,
		};
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/