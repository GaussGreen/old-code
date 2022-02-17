/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : pricer.cpp                                                   */
/*                                                                            */
/* DESCRIPTION : classe mere de tout les pricers                              */
/*                                                                            */
/* DATE        : Fri Sep 11 1998                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#include "pricer.h"

ARM_Pricer::ARM_Pricer(ARM_Security *sec, ARM_Model *mod)
{
	itsSecurity = sec;

	itsModel    = mod;
}

void ARM_Pricer::StartPricing(void)
{
    SetTime();
}

