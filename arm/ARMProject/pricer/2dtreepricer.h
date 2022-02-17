/*
 * $Log: 2dtreepricer.h,v $
 * Revision 1.3  2003/06/18 14:59:41  ebenhamou
 * remove link since included in armglob.h
 *
 * Revision 1.2  2003/06/18 11:13:28  ebenhamou
 * remove STL warning .
 *
 *
 */


 
#ifndef _2DTREEPRICER_H
#define _2DTREEPRICER_H

#include "treepricer.h"


class ARM_2DTreePricer : public ARM_TreePricer
{



public :
	ARM_2DTreePricer(void)
	{
	;
	};

	ARM_2DTreePricer(ARM_Security *sec, ARM_Model *mod);

	double Price(void) ;

	// Services

	ARM_2DTreePricer(const ARM_2DTreePricer &pricer)
    {
        this->BitwiseCopy(&pricer);
    }

	ARM_2DTreePricer &operator = (const ARM_2DTreePricer &pricer)
	{
		 (*this).ARM_TreePricer::operator = (pricer);

        BitwiseCopy(&pricer);

        return (*this);
	}

    
    void BitwiseCopy(const ARM_Object* opricer)
    {
         ;
    }

    void Copy(const ARM_Object* src)
    {    
		
		ARM_TreePricer::Copy(src);

        BitwiseCopy(src);
		
    }

};

#endif

