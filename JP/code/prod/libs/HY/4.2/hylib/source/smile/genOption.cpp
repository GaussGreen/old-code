// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/19/99 Neil Yang
// 
//


#include "cerror.h"
#include "genOption.h"
#include "RootFindBrent.h"




double GenOptionPricer(double  forward,        /* Forward price */
					   BaseOption*  baseOption, /* option info */
					   double  volatility,     /* */
					   double  discountRate,   /* */
					   BaseFunction*  inputProcess)
{

	GtoErrMsgOn();
	double stddev = volatility*sqrt( baseOption->get_yearsToExpiry());
	double boundLo = -INTEGRAL_LIMIT;
	double boundHi =  INTEGRAL_LIMIT;
	double tempResult;

	

	inputProcess->set_stddev(stddev);
	inputProcess->set_time(baseOption->get_yearsToExpiry());

	//inputProcess->set_amplitude(1);
	// integrate to get the match for fwd
	
	CallOption fwd(0,0,0);
	Payoff fwdPayoff(&fwd, inputProcess);
	ForwardSolverFunc  solverFunc(forward, &fwdPayoff);



	/*tempResult= Integral(*inputProcess,         
						 boundLo,          
						 boundHi,           
						 INTEGRAL_MAX_ERROR);  */

	tempResult = RootFindBrent(solverFunc,
							   forward,
							   ROOTLO,
							   ROOTHI);


	// set process to match fwd
	inputProcess->set_amplitude(tempResult);

	//construct payoff

	Payoff  payoffFunc(baseOption,inputProcess);

	return Integral(payoffFunc,         
					boundLo,          
					boundHi);

}







