/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: stringconvert.h,v $
 * Revision 1.1  2003/10/08 16:45:06  ebenhamou
 * Initial revision
 *
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file curveconvert.h
 *
 *  \brief files to convert a curve in a reference value
 *
 *	\author  Richard Guillemot
 *	\version 1.0
 *	\date August 2005
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPBASE_CURVECONVERT_H
#define _INGPBASE_CURVECONVERT_H

/// use our macro for namespace
#include "port.h"
#include "curve.h"
#include "curvetypedef.h"
#include "refvalue.h"

CC_BEGIN_NAMESPACE( ARM )

// Function to convert a curve into a reference value
ARM_ReferenceValue* CurveToRefValue( const ARM_Curve& curve, double asOfDate, double base = 1.0);
ARM_ReferenceValue GPCurveToRefValue( const ARM_Curve& curve, double asOfDate, double base = 1.0 );

// Function to convert a reference value into a curve
ARM_Curve* RefValueToCurve( const ARM_ReferenceValue& refVal, double asOfDate, double base = 1.0);
ARM_Curve RefValueToCurve( ARM_ReferenceValue* refVal, double asOfDate, double base = 1.0, int calcMethod = -1111);

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
