/*----------------------------------------------------------------------------*

    zerovsk.cpp
 
    This file implements the ARM_ZeroVasicek class, a class for computing 
    a ARM_ZeroCurve

*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "zerovsk.h"
#include "expt.h"





ARM_ZeroVasicek::ARM_ZeroVasicek(void)
{
    SetName(ARM_ZERO_VASICEK);
}

 

ARM_ZeroVasicek::ARM_ZeroVasicek(ARM_Date& asOf, 
                                 double* parameters):ARM_ZeroCurve(asOf)
{
    SetName(ARM_ZERO_VASICEK);

    ARM_Vector* param = new ARM_Vector(4, parameters);
 
    SetParameters(param);

    GenerateFields();
}



ARM_ZeroVasicek::ARM_ZeroVasicek(ARM_Date& asOf) : ARM_ZeroCurve(asOf) 
{
   	int i;
 
	double vpar[4];
 
    SetName(ARM_ZERO_VASICEK);
 
    for (i = 0; i < 4; i++)
    {
        vpar[i] = 0.0;
    }
 

    ARM_Vector* param = new ARM_Vector(4, vpar);
 
    SetParameters(param);
 
    GenerateFields();
}



ARM_ZeroVasicek::ARM_ZeroVasicek(ARM_ZeroVasicek& zeroVasi) 
                : ARM_ZeroCurve(zeroVasi)
{
    SetName(ARM_ZERO_VASICEK);

    BitwiseCopy(&zeroVasi);
}



ARM_ZeroVasicek& ARM_ZeroVasicek::operator = (ARM_ZeroVasicek& zeroVasi)
{
    (*this).ARM_ZeroCurve::operator = (zeroVasi);

    BitwiseCopy(&zeroVasi);

    return(*this);
}
    


/*

int ARM_ZeroVasicek::GetParameters(double *vasiCoeffs)
{
    if (vasiCoeffs)
       MEMCPY(vasiCoeffs, GetParameters(), 4 * sizeof(double));

    return(1);
}

*/



/*----------------------------------------------------------------------------*
    Returns the discount price with maturity yearTerm years from settlement, 
    computed from CDC's implementation of vasicek method.
*----------------------------------------------------------------------------*/

double ARM_ZeroVasicek::DiscountFunction(double yearTerm)
{
    double    y, z;


    if ( yearTerm < 0.0 ) 
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");
    }

    y = GetParameters()->Elt(1) - 
        GetParameters()->Elt(2) * g1(GetParameters()->Elt(0), yearTerm) +
        GetParameters()->Elt(3) * g2(GetParameters()->Elt(0), yearTerm) ;

    z = exp(-y*yearTerm);

    return(z);
}



double ARM_ZeroVasicek::D1DiscountFunction(double yearTerm)
{
    double    zp, pas=1e-6;

 
    if ( yearTerm < 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");
    }
 
    zp = (DiscountFunction(yearTerm+yearTerm*pas)
         -DiscountFunction(yearTerm-yearTerm*pas))/(2.*pas*yearTerm);
 
    return(zp);
}



double g1(double a, double theta)
{
    double g1;

    
    if ( theta < K_DOUBLE_TOL )
    {
       g1 = 1.0;
    }
    else
    {
       g1 = (1.0-exp(-a*theta))/(a*theta);
    }
    
    return(g1);
}



double g2(double a, double theta)
{
    double g2;


    if ( theta < K_DOUBLE_TOL )
    {
       g2 = 0.;
    }
    else
    {
       g2 = pow((1.0-exp(-a*theta)),2.0)/(4.0*a*theta);
    }

    return(g2);
}


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
