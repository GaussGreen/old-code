/*
 * $Log: twopluspricer.cpp,v $
 * Revision 1.1  2000/09/28 08:52:12  mab
 * Initial revision
 *
 */ 




/*----------------------------------------------------------------------------*
 
  twopluspricer.cpp 
 
  2+Model Pricer Methods
 
*----------------------------------------------------------------------------*/



#include "twopluspricer.h"
#include "security.h"




ARM_TwoPlusPricer::ARM_TwoPlusPricer(ARM_Security* sec, ARM_Model* mod)
                                    :ARM_TreePricer( sec, mod)

{
}



double ARM_TwoPlusPricer::Price(void)
{
	try
	{
		itsModel->BeFittedTo(itsSecurity);
	}


	catch(Exception& m)
    {
        m.DebugPrint();
 
        throw Exception(__LINE__, __FILE__, ERR_SWAP_CALC_PB,
                      "Problem in  function ARM_TwoPlusPricer::Price");
    }
	
	try
	{
		itsSecurity->PrepareToPrice(itsModel->GetStartDate());
	}


	catch(Exception& m)
    {
        m.DebugPrint();
 
        throw Exception(__LINE__, __FILE__, ERR_SWAP_CALC_PB,
                      "Problem in  function ARM_TwoPlusPricer::Price");
    }	
				
	return(0.0);
}



/*---------------------------------------------------------------------------*/
/*---- End Of File ----*/
