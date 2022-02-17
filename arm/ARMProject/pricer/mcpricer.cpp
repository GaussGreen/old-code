/*
 * $Log: mcpricer.cpp,v $
 * Revision 1.1  1998/11/20 11:07:55  nicolasm
 * Initial revision
 *
 */


#include "mcpricer.h"


double ARM_MonteCarloPricer::GetEstimatedSdev(void)
{
	if (!itsIsValidSdev)
				Price();
			
	return itsSdev;
}


void ARM_MonteCarloPricer::SetEstimatedSdev(double new_sdev) 
{
	itsIsValidSdev = 1;

	itsSdev = new_sdev;
}


