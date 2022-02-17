/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file curvefactory.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */

#ifndef _INGPBASE_CURVEFACTORY_H
#define _INGPBASE_CURVEFACTORY_H

/// this headers has to come first
/// as firsttoinc defines pragma for warning to avoid
#include "port.h"
#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

class ARM_RootObject;

struct ARM_CurveFactory
{
	static ARM_RootObject* CreateGenericCurve(
		const vector<double>& abscisses,
		const vector<double>& ordinates,
		long rowsNb,
		long colsNb,
		bool sortAbscisses,
		const string& interpolatorName,
		bool alwaysMulti);
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
