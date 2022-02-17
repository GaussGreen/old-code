/*
 * $Log
 */ 

#ifndef _TREEPRICER_H
#define _TREEPRICER_H

#include "pricer.h"


class ARM_TreePricer : public ARM_Pricer
{



public :
	ARM_TreePricer(void)
	{
	;
	};

	ARM_TreePricer(ARM_Security *sec, ARM_Model *mod);


	// Services

	ARM_TreePricer(const ARM_TreePricer &pricer)
    {
        this->BitwiseCopy(&pricer);
    }

	ARM_TreePricer &operator = (const ARM_TreePricer &pricer)
	{
		 (*this).ARM_Pricer::operator = (pricer);

        BitwiseCopy(&pricer);

        return (*this);
	}

    
    void BitwiseCopy(const ARM_Object* opricer)
    {
         ;
    }

    void Copy(const ARM_Object* src)
    {    
		
		ARM_Pricer::Copy(src);

        BitwiseCopy(src);
		
    }

};

#endif
