/*
 * $Log: zerospli.cpp,v $
 * Revision 1.2  2002/05/30 13:35:19  mab
 * Introducing RCS Mark
 *
 */


/*----------------------------------------------------------------------------*

    zerospli.cpp
 
    This file implements the ARM_ZeroSplines class, a class for 
    computing a ARM_ZeroCurve using splines.

*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "zerospli.h"
#include "expt.h"




ARM_ZeroSplines::ARM_ZeroSplines(ARM_Date& asOf, int numSplineCoeffs,
                                 double* splineCoeffs,
                                 double Int1, double Int2, 
                                 double Int3)
                                 :ARM_ZeroCurve(asOf)
{
    int    i;
    

    Init();

    SetName(ARM_ZERO_SPLINES);

    //    check input

    if ( numSplineCoeffs > K_MAX_NUM_SPLINE_COEFFS || numSplineCoeffs == 0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Number of spline coeff should be less than 20");
    }

    // Set variables
    
    itsNumSplineCoeffs = numSplineCoeffs;

    itsInt1=Int1;
    itsInt2=Int2;
    itsInt3=Int3;

    if ( splineCoeffs == NULL )
    {
       itsNumSplineCoeffs = 0;

       return;
    }

    ARM_Vector* parameters = new ARM_Vector(numSplineCoeffs);

    for (i = 0; i < numSplineCoeffs; i++) 
    {
        itsSplineCoeffs[i] = splineCoeffs[i];

        parameters->Elt(i) = splineCoeffs[i];
    }

    SetParameters(parameters);

    GenerateFields();
}



ARM_ZeroSplines::ARM_ZeroSplines(ARM_ZeroSplines& zeroSplines)
                :ARM_ZeroCurve(zeroSplines)
{
    Init();

    SetName(ARM_ZERO_SPLINES);

    BitwiseCopy(&zeroSplines);
}



ARM_ZeroSplines & ARM_ZeroSplines::operator = (ARM_ZeroSplines& zeroSplines)
{
    (*this).ARM_ZeroCurve::operator = (zeroSplines);


    BitwiseCopy(&zeroSplines);

    return(*this);
}


    
void ARM_ZeroSplines::EstimeBracket(ARM_Container* Assets, 
                                    ARM_Matrix* data,
                                    int* CST)
{
    int i;
    ARM_Vector*  splineCoeff = new ARM_Vector(9);
    ARM_Vector* splineCoeff2 = NULL;
    double* OPTION=NULL;



    if ( CST == NULL )
       CST = (int *) malloc (20*sizeof(int));

    if (GetParameters())
       splineCoeff2 = (ARM_Vector *) GetParameters()->Clone();
  
    for(i=0;i<itsNumSplineCoeffs;i++)
    {
        splineCoeff->Elt(i)=splineCoeff2->Elt(i);
    }
    splineCoeff->Elt(6)=itsInt1;
    splineCoeff->Elt(7)=itsInt2;
    splineCoeff->Elt(8)=itsInt3;

    itsNumSplineCoeffs=9;

    SetParameters(splineCoeff);

    OPTION=default_op(OPTION);
  
    OPTION[14]=6000.0;

    FitToMarketPrices(Assets, data,CST,OPTION);

    itsInt1=itsSplineCoeffs[6];
    itsInt2=itsSplineCoeffs[7];
    itsInt3=itsSplineCoeffs[8];

    itsNumSplineCoeffs=6;

    SetParameters(splineCoeff2);

    if (OPTION)
       free(OPTION);
}



int ARM_ZeroSplines::GetNumSplineCoeffs(void)
{
    return(itsNumSplineCoeffs);
}


 
int ARM_ZeroSplines::GetSplineCoeffs(double* splineCoeffs)
{
    MEMCPY(splineCoeffs, itsSplineCoeffs, itsNumSplineCoeffs*sizeof(double));

    return(1);
}



/*----------------------------------------------------------------------------*
    Returns the discount price with maturity yearTerm years from settlement, 
    computed from splines.
*----------------------------------------------------------------------------*/

double ARM_ZeroSplines::DiscountFunction(double yearTerm)
{
    double z = 0.0;
    double zeroShift = 0.0;



    if ( yearTerm < 0.0 )
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");

        return(-999999999);
    }

    if (GetBucketEndPeriod() > K_DOUBLE_TOL
        && yearTerm >= GetBucketStartPeriod() 
        && yearTerm < GetBucketEndPeriod())
    {
        zeroShift = 0.0001 * GetBPShift();
    }

    switch (itsNumSplineCoeffs) 
    {
        case (6) :
        {
            z = 1.0+itsSplineCoeffs[0]*yearTerm+itsSplineCoeffs[1]
                  *pow(yearTerm, 2.0) 
                  +itsSplineCoeffs[2]*pow(yearTerm, 3.0);
    
            if ( yearTerm > itsInt1)
               z += itsSplineCoeffs[3]*pow(yearTerm-itsInt1, 3.0);

            if ( yearTerm > itsInt2 ) 
               z += itsSplineCoeffs[4]*pow(yearTerm-itsInt2, 3.0);

            if ( yearTerm > itsInt3 )
               z += itsSplineCoeffs[5]*pow(yearTerm-itsInt3, 3.0);

            z *= exp(-zeroShift*yearTerm);
        };
        break;
        
        case (9) :
        {
            z = 1.0+itsSplineCoeffs[0]*yearTerm+itsSplineCoeffs[1]
                  *pow(yearTerm, 2.0) 
                  +itsSplineCoeffs[2]*pow(yearTerm, 3.0);
    
            if ( yearTerm > itsSplineCoeffs[6])
               z += itsSplineCoeffs[3]*pow(yearTerm-itsSplineCoeffs[6], 3.0);

            if ( yearTerm > itsSplineCoeffs[7] ) 
               z += itsSplineCoeffs[4]*pow(yearTerm-itsSplineCoeffs[7], 3.0);

            if ( yearTerm > itsSplineCoeffs[8] )
               z += itsSplineCoeffs[5]*pow(yearTerm-itsSplineCoeffs[8], 3.0);

            z *= exp(-zeroShift*yearTerm);
        };
        break;
            
        default: 
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Number of spline coeff should be equal to 6 or 9");
    }
    
    return(z);
}



/*----------------------------------------------------------------------------*
    Returns the d(discount price) / d(yearTerm) with maturity yearTerm 
         years from settlement, computed from splines.
*----------------------------------------------------------------------------*/
double ARM_ZeroSplines::D1DiscountFunction(double yearTerm)
{
    double    zp = 0.0;


    
    if ( yearTerm < 0.0 ) 
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");

       return(0.0);
    }

    switch (itsNumSplineCoeffs) 
    {
        case (6) :
        {
            zp = itsSplineCoeffs[0] + 2.0 * itsSplineCoeffs[1] * yearTerm 
              + 3.0 * itsSplineCoeffs[2] * pow(yearTerm, 2.0);
    
            if ( yearTerm > itsInt1 )
               zp += 3.0 * itsSplineCoeffs[3] * pow(yearTerm - itsInt1, 2.0);

            if ( yearTerm > itsInt2 )
               zp += 3.0 * itsSplineCoeffs[4] * pow(yearTerm - itsInt2, 2.0);

            if ( yearTerm > itsInt3 ) 
               zp += 3.0 * itsSplineCoeffs[5] * pow(yearTerm - itsInt3, 2.0);
        };
        break;
            
        default: 
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Number of spline coeff should be equal to 6");
    }
    
    return(zp);
}


    
/*----------------------------------------------------------------------------*
    Returns the d2(discount price) / d(yearTerm)2 with maturity 
    yearTerm years from settlement, computed from splines.
*----------------------------------------------------------------------------*/

double ARM_ZeroSplines::D2DiscountFunction(double yearTerm)
{
    double    zp;
    


    if ( yearTerm < 0.0 ) 
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");

       return(0.0);
    }

    switch (itsNumSplineCoeffs) 
    {
        case (6) :
        {
            zp = 2.0*itsSplineCoeffs[1]+6.0*itsSplineCoeffs[2]*yearTerm;
    
            if ( yearTerm > itsInt1 ) 
            {
                zp += 6.0*itsSplineCoeffs[3]*(yearTerm-itsInt1);
            }

            if (yearTerm > itsInt2) 
            {
                zp += 6.0*itsSplineCoeffs[4]*(yearTerm - itsInt2);
            }

            if (yearTerm > itsInt3) 
            {
                zp += 6.0*itsSplineCoeffs[5]*(yearTerm-itsInt3);
            }            
        };
        break;
            
        default: 
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Number of spline coeff should be equal to 6");
    }
    
    return(zp);
}


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
