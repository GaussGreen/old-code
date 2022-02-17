/*
 * $Log: ipricer.h,v $
 * Revision 1.4  2003/06/18 15:08:59  ebenhamou
 * remove include
 *
 * Revision 1.3  2003/06/18 11:27:39  ebenhamou
 * remove STL warning
 *
 */


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : ipricer.h                                                    */
/*                                                                            */
/* DESCRIPTION : Interface pricer classe                                      */             
/*                                                                            */
/* DATE        : Fri Sep 11 1998                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _IFPRICER_H
#define _IFPRICER_H

#include "pricer.h"

#include "security.h"
#include "model.h"

class ARM_IFPricer : public ARM_Object
{
public :
	ARM_IFPricer(void)
	{
		itsPricer = NULL;
	}


	ARM_IFPricer(ARM_Security *sec, ARM_Model *mod);

    ~ARM_IFPricer(void)
    {
        if (itsPricer)
            delete itsPricer;
    }

	double Price(void);
    ARM_Vector* ComputePrices(void);

private:
	ARM_Pricer *itsPricer;

};

#endif
