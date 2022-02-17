/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file curve.h
 *  \brief various utility for the ARM_T_Curve
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date July 2004
 */

#ifndef _INGPBASE_CURVEUTILS_H
#define _INGPBASE_CURVEUTILS_H

#include "port.h"
#include "env.h"			/// to have strict validation in debug mode
#include "curve.h"

CC_BEGIN_NAMESPACE( ARM )

template <typename T, typename U>	T IntegrateStepWise( const ARM_T_Curve<T,U>& curve, T a, T b )
{
	/// eliminate trivial one point case
	if( 1 == curve.size() )
	{
		return  curve.GetOrdinate(0)*(b-a);
	}
	else
	{
		int absB = curve.lower_bound(b)-1;
		double sum = curve.GetOrdinate(absB) * (b-curve.GetAbscisse(absB));
		int absA = curve.lower_bound(a);
	
		/// loop
		for(int i=absB-1; i>=absA; --i )
			sum += curve.GetOrdinate(i) * ( curve.GetAbscisse(i+1)-curve.GetAbscisse(i));

		/// last part
		/// is it the first point?
		if( absA )
		{
			if ( absA < curve.size() )
				sum += curve.GetOrdinate(absA-1) * (curve.GetAbscisse(absA)-a);
			else
				sum += curve.GetOrdinate(absB) * (curve.GetAbscisse(absB)-a);
		}
		else
			sum += curve.GetOrdinate(absA) * (curve.GetAbscisse(absA)-a);
		return sum;
	}

}


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
