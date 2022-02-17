/*
 * $Log: zerocbdf.cpp,v $
 * Revision 1.4  2004/01/21 14:31:08  jpriaudel
 * modif dans le bump
 *
 * Revision 1.3  2002/10/11 08:24:19  mab
 * Improvements
 *
 */


/*----------------------------------------------------------------------------*
 
    zerocbdf.cpp

*----------------------------------------------------------------------------*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "linalg.h"
#include "util.h"
#include "zerocbdf.h"
#include "currency.h"
#include "expt.h"




/*---------------------------------------------------------------------------*/


ARM_ZeroCubDiff::ARM_ZeroCubDiff(ARM_Date& asOf, 
                                 ARM_Vector* yearTerms,
                                 ARM_Vector* zeroYields, 
                                 int compoundMeth):ARM_ZeroCurve(asOf)
{
    Init();


    itsCompoundMeth  = compoundMeth;


    if ( yearTerms->GetSize() != zeroYields->GetSize() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
        "The input year terms and yields must have the same size");
    }


    //  On ajoute un point au debut et un point a la fin

    int size = yearTerms->GetSize();

    ARM_Vector* terms = new ARM_Vector(yearTerms, size+2, 0, size-1, 1);

    terms->Elt(0) = 0.0;
    terms->Elt(size+1) = 1000.0;

    ARM_Vector* yields = new ARM_Vector(zeroYields, size+2, 0, size-1, 1);

    yields->Elt(0) = linInterpol(yearTerms, 0.0, zeroYields);
    yields->Elt(size+1) = zeroYields->Elt(size-1);

    SetYearTerms(terms);

    SetZeroRates(yields);

    GenerateDateTerms();

    GenerateDiscountFactors();


    // Gestion du shift

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

    ComputePolynomCoef();
}



ARM_ZeroCubDiff::ARM_ZeroCubDiff(ARM_Date& asOf, 
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
    Init();


    // Set variables

    //  On ajoute un point au debut et un point a la fin
 
    ARM_Vector* yearTerms = GetYearTerms();
    ARM_Vector* zeroYields = GetZeroRates();

    int size = yearTerms->GetSize();
 
    ARM_Vector* yTerms = new ARM_Vector(yearTerms, size+2, 0, size-1, 1);
 
    yTerms->Elt(0) = 0.0;
    yTerms->Elt(size+1) = 1000.0;
 
    ARM_Vector* yields = new ARM_Vector(zeroYields, size+2, 0, size-1, 1);
 
    yields->Elt(0) = linInterpol(yearTerms, 0.0, zeroYields);
    yields->Elt(size+1) = zeroYields->Elt(size-1);
 
    SetYearTerms(yTerms);
 
    SetZeroRates(yields);

    GenerateDateTerms();

    GenerateDiscountFactors();


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

    itsCompoundMeth  = K_COMP_CONT;

    ComputePolynomCoef();
}



ARM_ZeroCubDiff::ARM_ZeroCubDiff(ARM_Date& asOf, 
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
    Init();

    SetName(ARM_ZERO_CUBDIFF);

    // Set variables

     //  On ajoute un point au debut et un point a la fin
 
    int size = GetYearTerms()->GetSize();
    ARM_Vector* zeroYields = GetZeroRates();
 
    ARM_Vector* yTerms = new ARM_Vector(GetYearTerms(), size+2, 0, size-1, 1);
 
    yTerms->Elt(0) = 0.0;
    yTerms->Elt(size+1) = 1000.0;
 
    ARM_Vector* Yields = new ARM_Vector(zeroYields, size+2, 0, size-1, 1);
 
    Yields->Elt(0) = linInterpol(GetYearTerms(), 0.0, zeroYields);
    Yields->Elt(size+1) = zeroYields->Elt(size-1);
 
    SetYearTerms(yTerms);
 
    SetZeroRates(Yields);


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

    itsCompoundMeth  = K_COMP_CONT;

    ComputePolynomCoef();
}



ARM_ZeroCubDiff::ARM_ZeroCubDiff(ARM_ZeroCubDiff& zeroSpliCub)
                   :ARM_ZeroCurve(zeroSpliCub)
{
    Init();

    SetName(ARM_ZERO_SPLICUB);

    BitwiseCopy(&zeroSpliCub);
}



ARM_ZeroCubDiff& ARM_ZeroCubDiff::operator= (ARM_ZeroCubDiff& zeroCubDiff)
{
    (*this).ARM_ZeroCurve::operator = (zeroCubDiff);

    BitwiseCopy(&zeroCubDiff);

    return(*this);
}



/*----------------------------------------------------------------------------*
   Returns the discount price with maturity yearTerm years from settlement, 
    computed from splines.
*----------------------------------------------------------------------------*/

double ARM_ZeroCubDiff::DiscountFunction(double yearTerm)
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

    intYield = ComputeInterpolValue(yearTerm)/100.0;

    if ( itsCompoundMeth == 0 )
    {
       z = exp(-yearTerm*(intYield + zeroShift));
    }
 
    if ( itsCompoundMeth == -1 )
    {
       z = 1.0/(1.0+yearTerm*(intYield + zeroShift));
    }
 
    if ( itsCompoundMeth > 0 )
    {
       z = pow(1.0+(zeroShift + intYield)/itsCompoundMeth,
                 -yearTerm*itsCompoundMeth);
    }

    return(z);
}



/*----------------------------------------------------------------------------*
    Returns the d(discount price) / d(yearTerm) with maturity yearTerm 
    years from settlement, computed from linear interpolation.
*----------------------------------------------------------------------------*/

double ARM_ZeroCubDiff::D1DiscountFunction(double yearTerm)
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
double ARM_ZeroCubDiff::D2DiscountFunction(double yearTerm)
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



double ARM_ZeroCubDiff::ComputeInterpolValue(double x)
{
    ARM_Vector *terms = GetYearTerms();
    ARM_Vector *yields = GetZeroRates();

    ARM_Matrix *polycoef = GetPolyCoef();

    int size = yields->GetSize();
    int n = size-1;

    double result = 0.0;



    // Case out of array

    if ( x >= terms->Elt(size-2) )
    {
       result = yields->Elt(size-2);

       return(result);
    }

    if ( x <= terms->Elt(1) )
    {
       result = linInterpol(terms, x, yields);

       return(result);
    }

    // recherche de l'index

    int index = n;
    int trouve = 0;

    while (index > 1 && !trouve)
    {
       if (terms->Elt(index) >= x)
          index--;
       else
          trouve = 1;
    }

    result = polycoef->Elt(index, 0);
    result += polycoef->Elt(index, 1)*x;
    result += polycoef->Elt(index, 2)*x*x;
    result += polycoef->Elt(index, 3)*x*x*x;

    return(result);
}



void ARM_ZeroCubDiff::ComputePolynomCoef(void)
{
    int i, size;
    double t1, t2, t3, t4;
    double r1, r2, r3, r4;
    double det = 0;


    ARM_Vector* terms = GetYearTerms();
    ARM_Vector* yields = GetZeroRates();

    size = terms->GetSize();

    ARM_Matrix* coef = new ARM_Matrix(size, 4, 0.0);

    ARM_Matrix mat4(4, 4, 0.0);
    ARM_Vector* SndMb = new ARM_Vector(4);

    for (i=1; i < size-2; i++)
    {
       t1 = terms->Elt(i-1);
       t2 = terms->Elt(i);
       t3 = terms->Elt(i+1);
       t4 = terms->Elt(i+2);
     
       r1 = yields->Elt(i-1);
       r2 = yields->Elt(i);
       r3 = yields->Elt(i+1);
       r4 = yields->Elt(i+2);

       mat4.Elt(0,0) = t2*t2*t2;
       mat4.Elt(0,1) = t2*t2;
       mat4.Elt(0,2) = t2;
       mat4.Elt(0,3) = 1;
       mat4.Elt(1,0) = t3*t3*t3;
       mat4.Elt(1,1) = t3*t3;
       mat4.Elt(1,2) = t3;
       mat4.Elt(1,3) = 1;
       mat4.Elt(2,0) = 3*t2*t2;
       mat4.Elt(2,1) = 2*t2;
       mat4.Elt(2,2) = 1;
       mat4.Elt(2,3) = 0;
       mat4.Elt(3,0) = 3*t3*t3;
       mat4.Elt(3,1) = 2*t3;
       mat4.Elt(3,2) = 1;
       mat4.Elt(3,3) = 0;

       SndMb->Elt(0) = r2;
       SndMb->Elt(1) = r3;
       SndMb->Elt(2) = ((r2-r1)/(t2-t1) + (r3-r2)/(t3-t2))/2;
       SndMb->Elt(3) = ((r3-r2)/(t3-t2) + (r4-r3)/(t4-t3))/2;

       mat4.LinSolve(SndMb, det);

       if ( det == 0.0 )
       {
          throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Error invert matrix");
       }

       coef->Elt(i, 0) = SndMb->Elt(3);
       coef->Elt(i, 1) = SndMb->Elt(2);
       coef->Elt(i, 2) = SndMb->Elt(1);
       coef->Elt(i, 3) = SndMb->Elt(0);
    }


    SetPolyCoef(coef);
}



ARM_ZeroCurve* ARM_ZeroCubDiff::GenerateShiftCurve(ARM_CRV_TERMS& Term, ARM_Vector* epsilon)
{
   if (GetMktData() == NULL)
   {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "The input curve cant be shifted");
   }

   ARM_ZeroCubDiff* zc;
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
      for (k = 0; k < ARM_NB_TERMS; k++)
      {
         if (!strcmp(Term[l],MktData->itsMktTerms[k]))
         {
            MktData->itsMktValue->Elt(k) = 
                    MktData->itsMktValue->Elt(k)+epsilon->Elt(l);

            k = ARM_NB_TERMS;
         }
      }
   }

   zc = new ARM_ZeroCubDiff(GetAsOfDate(),
                            MktData->itsMktTerms,
                            MktData->itsMktValue,
                            MktData->itsMMVsFut,
                            MktData->itsSwapVsFut,
                            MktData->itsraw,
                            MktData->itsCont_Lin,
                            0,
                            GetCurrencyUnit());

   return zc;
}



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/