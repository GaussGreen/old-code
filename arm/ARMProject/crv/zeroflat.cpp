/*----------------------------------------------------------------------------*

 	zeroflat.cpp

*----------------------------------------------------------------------------*/


#include "zeroflat.h"

double ARM_ZeroFlat::DiscountFunction(double yearTerm)
{
        return(exp(-0.01*yearTerm*itsYield));
}


ARM_ZeroFlat::ARM_ZeroFlat(ARM_Date& asOf, double flatRate,
                           ARM_Currency* ccy):ARM_ZeroCurve(asOf)
{
	SetName(ARM_ZERO_FLAT);

	itsYield = flatRate;

    SetCurrencyUnit(ccy);

    GenerateFields();
}



ARM_ZeroFlat::ARM_ZeroFlat(ARM_ZeroFlat& flatCrv):ARM_ZeroCurve(flatCrv)
{
	SetName(ARM_ZERO_FLAT);

    BitwiseCopy(&flatCrv);
}


	
ARM_ZeroFlat& ARM_ZeroFlat::operator = (ARM_ZeroFlat& flatCrv)
{
    (*this).ARM_ZeroCurve::operator = (flatCrv);
		
    BitwiseCopy(&flatCrv);

    return(*this);
}



void ARM_ZeroFlat::Set(ARM_Date& asOf, double flatRate)
{
    SetAsOfDate(asOf);
 
    itsYield = flatRate;

    GenerateFields();
}
 


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
