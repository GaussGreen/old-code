/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: calibdirection.h,v $
 * Revision 1.1  2003/12/02 07:51:19  emezzine, ebenhamou
 * Initial revision
 *
 *
 */

/*! \file calibdirection.h
 *
 *  \brief 
 *
 *	\author  E.M Ezzine E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#ifndef _INGPINFRA_CALIBDIRECTION_H
#define _INGPINFRA_CALIBDIRECTION_H

#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )

    enum ARM_CalibDirection
    {
        CalibDirection_Forward =0,
        CalibDirection_Backward,
        CalibDirection_None
    };

struct ARM_CalibrationTarget
{
    enum TargetType
    {
        PriceTarget=0,
        ImpliedVolatilityTarget,
        UnknownTarget
    };
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/