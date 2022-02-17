/*
 * $Log: othpricer.h,v $
 * Revision 1.4  2003/06/18 15:14:09  ebenhamou
 * remove include
 *
 * Revision 1.3  2003/06/18 11:32:16  ebenhamou
 * remove STL warning
 *
 */



/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : oldpricer.h                                                  */
/*                                                                            */
/* DESCRIPTION : classe pricers pour compatibilite avec existant              */
/*                                                                            */
/* DATE        : Fri Sep 11 1998                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _OTHPRICER_H
#define _OTHPRICER_H

#include "pricer.h"

class ARM_OtherPricer : public ARM_Pricer
{
    public :
		ARM_OtherPricer(void)
		{
		}

		ARM_OtherPricer(ARM_Security *sec, ARM_Model *mod) : 
                                                  ARM_Pricer(sec, mod)
		{
		};

		double Price(void);

};

#endif
