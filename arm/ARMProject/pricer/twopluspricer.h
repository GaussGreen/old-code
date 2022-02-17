/*
 * $Log: twopluspricer.h,v $
 * Revision 1.3  2003/06/18 15:14:45  ebenhamou
 * remove include
 *
 * Revision 1.2  2003/06/18 11:49:56  ebenhamou
 * remove STL warning
 *
 * Revision 1.1  2000/09/28 08:51:18  mab
 * Initial revision
 *
 */ 


/*----------------------------------------------------------------------------*
 
  twopluspricer.h

  Header of the 2+Model Pricer 

*----------------------------------------------------------------------------*/



#ifndef _TWOPLUSPRICER_H
#define _TWOPLUSPRICER_H

#include "treepricer.h"




class ARM_TwoPlusPricer : public ARM_TreePricer
{
    public :


    ARM_TwoPlusPricer(void)
    {
    }

    ARM_TwoPlusPricer(ARM_Security* sec, ARM_Model* mod);



    double Price(void); 



    // Services

    ARM_TwoPlusPricer(const ARM_TwoPlusPricer& pricer)
    {
        this->BitwiseCopy(&pricer);
    }

    ARM_TwoPlusPricer &operator = (const ARM_TwoPlusPricer& pricer)
    {
        (*this).ARM_TreePricer::operator = (pricer);

        BitwiseCopy(&pricer);

        return (*this);
    }

    
    void BitwiseCopy(const ARM_Object* opricer)
    {
    }

    void Copy(const ARM_Object* src)
    {    
        ARM_TreePricer::Copy(src);

        BitwiseCopy(src);
    }

};

#endif
/*------------------------------------------------------------------------------*/
/*---- End Of File ----*/
