/*
 * $Log: zerointlight.cpp,v $
 * Revision 1.2  2003/09/25 09:17:32  mab
 * "constification"
 *
 * Revision 1.1  2001/09/14 07:42:29  abizid
 * Initial revision
 *
 */


/*----------------------------------------------------------------------------*
 
    zerointlight.cpp
 
    This file implements the ARM_ZeroLInterpolLight class, a class for 
         computing a ARM_ZeroLInterpol class more fastly.

*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "zerointlight.h"
#include "currency.h"
#include "expt.h"
#include "newton.h"




/*---------------------------------------------------------------------------*/
/*
ARM_ZeroLInterpolLight::ARM_ZeroLInterpolLight(ARM_Date& asOf, 
                                               ARM_Vector* yearTerms,
                                               ARM_Vector* zeroYields, 
                                               int compoundMeth, int lastBucketInt, 
                                               int interpolMeth):
                                         ARM_ZeroLInterpol(asOf, yearTerms, zeroYields, 
                                        compoundMeth, lastBucketInt, interpolMeth)
{
    Init();
}
*/


ARM_ZeroLInterpolLight::ARM_ZeroLInterpolLight(double* yearTerms,
                                               double* zeroYield,
                                               int size):ARM_ZeroLInterpol()
{
    Init();

    SizeVect = size;
    SetTerms(yearTerms);
    SetRates(zeroYield);
}


ARM_ZeroLInterpolLight::ARM_ZeroLInterpolLight(
                                    const ARM_ZeroLInterpolLight& zeroLInterpol)
                   :ARM_ZeroLInterpol(zeroLInterpol)
{
    Init();

    BitwiseCopy(&zeroLInterpol);
}



ARM_ZeroLInterpolLight& ARM_ZeroLInterpolLight::operator = (const ARM_ZeroLInterpolLight&
                                                   zeroLInterpol)
{
    (*this).ARM_ZeroCurve::operator = (zeroLInterpol);

    BitwiseCopy(&zeroLInterpol);

    return(*this);
}



/*----------------------------------------------------------------------------*
   Returns the discount price with maturity yearTerm years from settlement, 
    computed from splines.
   We assume that the interpolation is constant, out of the range of YearTerms
*----------------------------------------------------------------------------*/

double ARM_ZeroLInterpolLight::DiscountFunction(double yearTerm)
{
    double z, intYield;

    if ( yearTerm < 0.0 )
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");
    }

    //  interpolate 
    int lastIndx = SizeVect-1;//SizeVect-1;

    if ( yearTerm < itsTerms[0] - K_DOUBLE_TOL ) 
    {
       intYield = itsZCRates[0]/100.0;
    }
    else if ( yearTerm >= itsTerms[lastIndx] - K_DOUBLE_TOL) 
    {
       intYield = itsZCRates[lastIndx]/100.0;
    }
    else
    {
       intYield = linInterpol(itsTerms, lastIndx, yearTerm, itsZCRates)/100.0;
    }

    z = exp(-yearTerm*intYield);

    return(z);
}

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
