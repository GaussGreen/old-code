/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions to convert
 *
 *	\file convertor.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include "gpclosedforms/convertor.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/vanille_bs_interface.h"

CC_BEGIN_NAMESPACE(ARM)

double ConvertFXOptionPriceToStrike( double targetDelta, 
		double fwd,
		double totalVol,
		int callPut, 
		double initValue )
{
	struct FXOptionToInverse : public DoubleToDoubleFunc 
	{
		double itsFwd;
		double itsTotalVol;
		int itsCallPut;
		long itsDeltaType;

		FXOptionToInverse( double fwd,
			double totalVol,
			int callPut,
			long deltaType = 0)
		: 
		    itsFwd(fwd),
			itsTotalVol(totalVol),
			itsCallPut(callPut)
			{}

		virtual double operator() (double strike ) const
		{
			/*switch(itsCallPut)
			case:*/
			double deltaFwd = Export_BlackSholes(0,itsFwd,itsTotalVol,1.0,strike,itsCallPut);
			double premium = Export_BlackSholes(itsFwd,itsTotalVol,1.0,strike,itsCallPut);
			double deltaFwdWpremium = deltaFwd-premium/itsFwd;

			return deltaFwdWpremium;
		}
	};
	
	FXOptionToInverse x(fwd, totalVol, callPut);
	return Inverse(x,Inverse::ALWAYSPOSITIVE)(targetDelta,initValue,initValue/5.0,1e-10);
}

CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
