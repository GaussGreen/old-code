/*----------------------------------------------------------------------------*
 
    zerosplisum.cpp
 
    This file implements the ARM_ZeroSpliSum class,     
    a class for computing a ARM_ZeroCurve using Summit's algorithm 
	to produced a new set of yearterms and cubic splined zero coupon yields
	from an initial set of yearterms and zero coupon yields 
	or a 'Raw' computed zero coupon yields from a Market Data input.

	Rmq : the core computing of the splined zero yields is identical with
	the way it is done in the zerospl3 class but Summit's algorithm differs
	in the following points :

	- the computing of discounts on a not calibrated maturity is not done with
	the splined coefficients; Summit first adds, as the 'Forward' methodology,
	extra points (distanced each 3 months inside two consecutive calibrated maturities)
	on which the yields are computed with the splined coefficients and then memorize
	this new set of year terms and zero yields; to compute a discount on a maturity
	not belonging this set Summit use the standard interpolation rules
	(i.e. 'Linear' or 'Continuous');

	- in the 'Fitted Zero' algorithm, Summit is taking account the set of year terms
	and zero coupon yields and adapt them to refit the initial Market Data
	(mainly from swap segment)
 

*----------------------------------------------------------------------------*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "linalg.h"
#include "util.h"
#include "zerosplsum.h"
#include "currency.h"
#include "expt.h"




/*---------------------------------------------------------------------------*/


ARM_ZeroSpliSum::ARM_ZeroSpliSum(ARM_Date& asOf, 
                                 ARM_Vector* yearTerms,
                                 ARM_Vector* zeroYields, 
                                 int lastBucketInt)
                                 :ARM_ZeroCurve(asOf)
{
    Init();


    if ( yearTerms->GetSize() != zeroYields->GetSize() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
        "The input year terms and yields must have the same size");
    }

    // Set variables

    if (lastBucketInt)
    {
       SetYearTerms(yearTerms);

       SetZeroRates(zeroYields);
    }
    else
    {
       int size = yearTerms->GetSize();

       ARM_Vector* terms = new ARM_Vector(yearTerms, size+2, 0, size-1);

       terms->Elt(size) = 1000.0;
       terms->Elt(size+1) = 1500.0;

       ARM_Vector* yields = new ARM_Vector(zeroYields, size+2, 0, size-1);

       yields->Elt(size) = zeroYields->Elt(size-1);
       yields->Elt(size+1) = zeroYields->Elt(size-1);

       SetYearTerms(terms);

       SetZeroRates(yields);
    }

    int Size = GetYearTerms()->GetSize();

    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);

    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
           BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    SumGenerateD2Rates();

    GenerateFields();
}



ARM_ZeroSpliSum::ARM_ZeroSpliSum(ARM_Date& asOf, 
                                 ARM_CRV_TERMS& terms, 
                                     ARM_Vector* mktData, 
                                     int MMVsFut, 
                                     int SwapVsFut, 
                                     int raw, 
                                     int Cont_Lin, 
                                     int lastBucketInt, 
                                     ARM_Currency* ccy) : ARM_ZeroCurve(asOf, 
                                                           terms, 
                                                           mktData, 
                                                           MMVsFut, 
                                                           SwapVsFut, 
                                                           raw, 
                                                           Cont_Lin, 
                                                           ccy)
{
    // Check variables

    Init();

    // Set variables

    if (!lastBucketInt)
    {
       int size = GetYearTerms()->GetSize();

       ARM_Vector* terms = new ARM_Vector(GetYearTerms(), size+1, 0, size-1);

       terms->Elt(size) = 1000.0;

       ARM_Vector* yields = new ARM_Vector(GetZeroRates(), size+1, 0, size-1);

       yields->Elt(size) = GetZeroRates()->Elt(size-1);        
        
       ARM_Vector* df = new ARM_Vector(GetDiscountFactors(), size+1, 0, size-1);

       df->Elt(size) = -log(GetZeroRates()->Elt(size-1)/100.0)/1000.0;        
        
       SetYearTerms(terms);

       SetZeroRates(yields);

       SetDiscountFactors(df);
    }

    int Size = GetYearTerms()->GetSize();
 
    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);
 
    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
            BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsLastBucketInt = lastBucketInt;

    SumGenerateD2Rates();
}



ARM_ZeroSpliSum::ARM_ZeroSpliSum(ARM_Date& asOf, 
                                 ARM_CRV_TERMS& terms, 
                                 ARM_Vector* mktData, 
                                 ARM_Container* bonds, 
                                 ARM_Vector* yields,
                                 int MMVsFut, 
                                 int lastBucketInt, 
                                 ARM_Currency* ccy):ARM_ZeroCurve(asOf, 
                                                           terms, 
                                                           mktData, 
                                                           bonds, 
                                                           yields, 
                                                           MMVsFut, 
                                                           ccy)
{
    // Check variables

    Init();

    // Set variables

    if (!lastBucketInt)
    {
       int size = GetYearTerms()->GetSize();

       ARM_Vector* terms = new ARM_Vector(GetYearTerms(), size+1, 0, size-1);

       terms->Elt(size) = 1000.0;

       ARM_Vector* yields = new ARM_Vector(GetZeroRates(), size+1, 0, size-1);

       yields->Elt(size) = GetZeroRates()->Elt(size-1);        
        
       ARM_Vector* df = new ARM_Vector(GetDiscountFactors(), size+1, 0, size-1);

       df->Elt(size) = -log(GetZeroRates()->Elt(size-1)/100.0)/1000.0;        
        
       SetYearTerms(terms);

       SetZeroRates(yields);

       SetDiscountFactors(df);
    }

    int Size = GetYearTerms()->GetSize();
 
    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);
 
    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
            BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsLastBucketInt = lastBucketInt;

    SumGenerateD2Rates();
}




ARM_ZeroSpliSum::ARM_ZeroSpliSum(ARM_ZeroSpliSum& zeroSpliSum)
                :ARM_ZeroCurve(zeroSpliSum)
{
    Init();

    BitwiseCopy(&zeroSpliSum);
}



ARM_ZeroSpliSum& ARM_ZeroSpliSum::operator= (ARM_ZeroSpliSum& zeroSpliSum)
{
    (*this).ARM_ZeroCurve::operator = (zeroSpliSum);

    BitwiseCopy(&zeroSpliSum);

    return(*this);
}



/*----------------------------------------------------------------------------*
   Returns the discount price with maturity yearTerm years from settlement, 
    computed from splines.
*----------------------------------------------------------------------------*/

double ARM_ZeroSpliSum::DiscountFunction(double yearTerm)
{
    double z, zeroShift=0.0, intYield;
    
    int    zrIsCloned = 0;

    ARM_Vector* zeroRates = NULL;


    if ( yearTerm < 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");
    }

    if (GetBucketEndPeriod() > K_DOUBLE_TOL
        && yearTerm >= GetBucketStartPeriod() 
        && yearTerm < GetBucketEndPeriod())
    {
       zeroShift = 0.0001 * GetBPShift();
    }
    

    if ( GetBPShift() > 0 )
    {
       zeroRates = (ARM_Vector *) GetZeroRates()->Clone();

       zrIsCloned = 1;
    }
    else
    {
       zeroRates = (ARM_Vector *) GetZeroRates();
    }

    intYield = ComputeSpline(yearTerm)/100.0;

    z = exp(-yearTerm*(intYield + zeroShift));

    return(z);
}



/*----------------------------------------------------------------------------*
    Returns the d(discount price) / d(yearTerm) with maturity yearTerm 
    years from settlement, computed from linear interpolation.
*----------------------------------------------------------------------------*/

double ARM_ZeroSpliSum::D1DiscountFunction(double yearTerm)
{
    double    zp;
    

    if ( yearTerm < 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");
    }

    
    zp = (DiscountFunction(yearTerm+0.000001)
          -DiscountFunction(yearTerm-0.000001))*500000.0; 
    
    return(zp);
}


    
/*----------------------------------------------------------------------------*
    Returns the d2(discount price) / d(yearTerm)2 with maturity yearTerm 
    years from settlement, computed from splines.
*----------------------------------------------------------------------------*/
double ARM_ZeroSpliSum::D2DiscountFunction(double yearTerm)
{
    double zp;


    if ( yearTerm < 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");
    }

    zp = (D1DiscountFunction(yearTerm+0.000001)
           -D1DiscountFunction(yearTerm-0.000001))*500000.0; 

    return(zp);
}



double ARM_ZeroSpliSum::ComputeSpline(double x)
{
    ARM_Vector* term = GetYearTerms();
    ARM_Vector* y = GetZeroRates();

    ARM_Vector *d2ydt = itsD2Rates;

    int size = y->GetSize();
    int n = size-1;

    double result = 0.0;

    int imin = 0; 
    int imax = n;
    int i;


    // Case out of array
/*
    if ( x >= term->Elt(n) )
    {
       result = (y->Elt(1)-y->Elt(0))/(term->Elt(1)- term->Elt(0));
       result *= x -  term->Elt(0);
       result += y->Elt(0);
    }

    if ( x <= term->Elt(0) )
    {
       result = (y->Elt(n)-y->Elt(n-1))/(term->Elt(n)-term->Elt(n-1));
       result *= x-term->Elt(n);
       result += y->Elt(n);
    }
*/
    while ( (imax-imin) > 1)
    {
       i = (imax+imin)/2;

       if (term->Elt(i) > x)
           imax = i;
       else
           imin = i;
           
    }

    double h = term->Elt(imax) - term->Elt(imin);

    if ( h == 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "All element of YearTerm must be different");
    }

    double A = (term->Elt(imax) - x)/h;
    double B = 1 - A;

    result = A*y->Elt(imin)+B*y->Elt(imax); // linear interpol

    result += ((A*A*A-A)*d2ydt->Elt(imin)+(B*B*B-B)*d2ydt->Elt(imax))*h*h/6.0;

    return(result);
}



void ARM_ZeroSpliSum::SumGenerateD2Rates(void)
{
    int i, Size;


    // x are dates, y are the zerorates and y2 the snd derivative of zerorates

    ARM_Vector* x = GetYearTerms();
    ARM_Vector* y = GetZeroRates();
    ARM_Vector *y2 = NULL;

    Size = x->GetSize();

//    ARM_Vector SndMb(Size);
//    ARM_TDiag  mat(Size);

    y2 = new ARM_Vector(Size);
    ARM_Vector u(Size-1);

	double sig;
	double p;

	y2->Elt(0) = 0.;
	u.Elt(0) = 0.;

	for (i=1;i<Size-1;i++)
	{
		sig = (x->Elt(i) - x->Elt(i-1))/(x->Elt(i+1) - x->Elt(i-1));
		p = sig * y2->Elt(i-1)+2.;
		y2->Elt(i) = (sig-1.)/p;
		u.Elt(i) = (y->Elt(i+1)-y->Elt(i))/(x->Elt(i+1)-x->Elt(i)) - (y->Elt(i)-y->Elt(i-1))/(x->Elt(i)-x->Elt(i-1));
		u.Elt(i) = (6.*u.Elt(i)/(x->Elt(i+1)-x->Elt(i-1))-sig*u.Elt(i-1))/p;
	}

	y2->Elt(Size-1) = 0.;

	for (i=Size-2;i>=0;i--)
	{
		y2->Elt(i) = y2->Elt(i)*y2->Elt(i+1) + u.Elt(i);
	}
//    y2 = mat.LinSolve(&SndMb);

    SetD2Rates(y2);
}



ARM_ZeroCurve* ARM_ZeroSpliSum::GenerateShiftCurve(ARM_CRV_TERMS& Term, ARM_Vector* epsilon)
{
    if (GetMktData() == NULL)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
          "The input curve cant be shifted");
    }

   ARM_ZeroSpliSum* zc;
   ARM_Date tmpdat;
   ARM_MarketData* MktData = (ARM_MarketData*) GetMktData()->Clone();


   // verification de l'existence des datamarkets
   if (!MktData)
	  return NULL;

   int szinp =0;
   szinp = epsilon->GetSize();

   int k=0;
   int l=0;

   for (l=0; l<szinp; l++)
   {
	  for (k=0; k<ARM_NB_TERMS; k++)
	  {
		 if (!strcmp(Term[l],MktData->itsMktTerms[k]))
		 {
			MktData->itsMktValue->Elt(k) = 
					 MktData->itsMktValue->Elt(k) + epsilon->Elt(l);

			k = ARM_NB_TERMS;
		 }
	  }
   }

   zc = new ARM_ZeroSpliSum(GetAsOfDate(),
                     MktData->itsMktTerms,
                     MktData->itsMktValue,
                     MktData->itsMMVsFut,
                     MktData->itsSwapVsFut,
                     MktData->itsraw,
                     MktData->itsCont_Lin,
                     this->itsLastBucketInt,
                     GetCurrencyUnit());

   return zc;
}


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/