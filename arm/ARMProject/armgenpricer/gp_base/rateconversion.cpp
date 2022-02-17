/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: rateconversion.cpp,v $
 * Revision 1.1  2003/10/08 16:39:43  ebenhamou
 * Initial revision
 *
 */


/*! \file rateconversion.cpp
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpbase/rateconversion.h"

// Kernel
#include "armdef.h"
#include <cmath>

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////////////////////////
/// function to convert from a decompounding method to another decompounding method
///////////////////////////////////////////////////////////////////////

double ConvertRateToRate(double rate, double term, int fromCompMeth, int toCompMeth)
{
    double newRate;

    switch(fromCompMeth) 
    {
	case K_COMP_ANNUAL :
	case K_COMP_SEMIANNUAL :
	case K_COMP_QUARTERLY :
	case K_COMP_BIMONTHLY :
	case K_COMP_MONTHLY :
        {
            switch (toCompMeth) 
            {
			case K_COMP_ANNUAL :
			case K_COMP_SEMIANNUAL :
			case K_COMP_QUARTERLY :
			case K_COMP_BIMONTHLY :
			case K_COMP_MONTHLY :
                {
                    newRate = (pow(1.0+(rate/double(fromCompMeth)),double(fromCompMeth)/double(toCompMeth))-1.0)* double(toCompMeth);
                    break;
                };
				
			case K_COMP_CONT :
                {
                    newRate = (log(1.0+rate/double(fromCompMeth)))*double(toCompMeth);
                    break;
                };
				
			case K_COMP_PROP :
                {
                    newRate = (pow(1.0+rate/double(fromCompMeth),double(fromCompMeth)*term)-1.0)/term;
                    break;
                };
				
			default :
				newRate = rate;
            }
            break;
        }
		
	case K_COMP_CONT :
        {
            switch (toCompMeth) 
            {
			case K_COMP_ANNUAL :
			case K_COMP_SEMIANNUAL :
			case K_COMP_QUARTERLY :
			case K_COMP_BIMONTHLY :
			case K_COMP_MONTHLY :
				{
                    newRate = (exp(rate/double(toCompMeth)) - 1.0)*double(toCompMeth);
                    break;
				};
			case K_COMP_PROP :
                {
                    newRate = (exp(rate*term) - 1.0)/ term;
                    break;
                };
			default :
				newRate = rate;
            };
            break;
        };
		
	case K_COMP_PROP :
        {
            switch (toCompMeth) 
            {
			case K_COMP_ANNUAL :
			case K_COMP_SEMIANNUAL :
			case K_COMP_QUARTERLY :
			case K_COMP_BIMONTHLY :
			case K_COMP_MONTHLY :
                {
                    newRate = (pow(1.0+rate*term, 1.0/double(term*toCompMeth))- 1.0) *double(toCompMeth);
                    break;
                };
				
			case K_COMP_CONT :
                {
                    newRate = log(1.0+rate*term)/term;
                    break;
                };
			default :
				newRate = rate;
            };
            break;
        };
	default :
		newRate = rate;
    }
    return newRate;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

