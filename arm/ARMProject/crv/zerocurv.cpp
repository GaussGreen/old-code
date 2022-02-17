/*
 * $Log: zerocurv.cpp,v $
 * Revision 1.35  2004/05/03 09:26:17  jpriaudel
 * modif pour taux negatifs + calendrier JPY en ZGJ
 *
 * Revision 1.34  2004/01/20 16:40:02  mab
 * AddMonths replaced everywhere by : AddMonths(nb, OTO_END_OF_MONTH);
 *
 * Revision 1.33  2004/01/09 11:14:02  jpriaudel
 * correctif sur le TOY
 *
 * Revision 1.32  2003/12/15 15:00:36  arm
 * SetMktData now in .h : because of strange compilation PB
 *
 * Revision 1.31  2003/11/13 13:30:45  jpriaudel
 * modif pour le credit
 *
 * Revision 1.30  2003/07/17 06:38:09  ebenhamou
 * for CptDF, if not value, throw an exception
 * Previously shown as a warning
 *
 * Revision 1.29  2003/07/01 18:58:49  arm
 * int i added istead for (int i;
 *
 * Revision 1.28  2003/06/17 16:14:05  jpriaudel
 * correction d'un bug cans le TOY
 *
 * Revision 1.27  2003/03/27 17:03:34  mab
 * Modif du la version precedente
 *
 * Revision 1.26  2003/03/27 10:41:15  mab
 * Amelioration des performances sur la generation des daily DF
 * dans le TOY
 *
 * Revision 1.25  2003/02/14 09:33:19  mab
 * Correction in ForwardYield (CountYears matu1, matu2)
 *
 * Revision 1.24  2003/02/11 15:03:40  mab
 * Corrections
 * to take in account the flag K_SWAP in
 * SK2 curves
 *
 * Revision 1.23  2002/11/15 10:15:34  mab
 *  Change the format when USD in : ZCFromMarketRatesTOY
 *
 * Revision 1.21  2002/09/24 13:16:47  mab
 * Added : TOY management (See TOY section)
 *
 * Revision 1.20  2002/09/17 16:44:09  mab
 * improvements in : void ARM_ZeroCurve::CptFuturesZeroRates
 * (see new flag : int MMandSizeSeq0 = 0; )
 * void ARM_ZeroCurve::CptSwapZeroRates
 * (see int idxj = 0;)
 *
 * Revision 1.19  2002/08/08 10:23:29  mab
 * in the 18M case :
 * valTmp = Nb/12.0; replaced by:
 * valTmp = ((double) Nb)/((double) 12.0);
 * and :
 *                            newIndex = int (valTmp*fxPayFreq);
 * by: newIndex = (int) (floor(valTmp*fxPayFreq));
 *
 * Revision 1.18  2002/07/29 13:33:29  mab
 * View() improvement
 *
 * Revision 1.17  2002/07/18 13:43:24  mab
 * Introduction of 18M tenor
 *
 * Revision 1.16  2002/05/30 13:40:09  mab
 * in terms[ARM_NB_TERMS][6] : 6 replaced by 12
 * in order to deal correctly with dates instead of maturities only
 *
 * Revision 1.15  2001/10/08 13:33:07  mab
 * Ameliorations
 *
 * Revision 1.14  2001/04/03 11:59:39  nicolasm
 * Acceleration
 *
 * Revision 1.13  2001/03/12 13:15:33  abizid
 * Optim MC
 *
 * Revision 1.12  2000/12/11 16:09:07  mab
 * Ds la generation de courbes : Pour la partie monetaire
 * Rajout d'un parametre a AddMonth
 *
 * Revision 1.11  2000/05/02 10:04:19  mab
 * passage par reference de ARM_Matrix ds les generations de courbes
 * in : CptFuturesZeroRates , CptFuturesZeroRates , CptFuturesZeroRatesGen
 * CptSwapZeroRates
 *
 * Revision 1.10  1999/09/15 15:18:58  nicolasm
 * Modif du calcul du taux dans le cas du toy
 *
 * Revision 1.9  1999/09/10 14:27:59  nicolasm
 * Ajout du turn of the year dans la courbe swap
 *
 * Revision 1.8  1999/06/09 13:03:29  nicolasm
 * Erreur de logique && et || ds fct ZcFromMktRateGen
 *
 * Revision 1.7  1999/06/08 07:53:55  nicolasm
 * Ajout test 'XX' et non 'X' dans ZCFromMarketDataGen
 *
 * Revision 1.6  1999/04/07 17:22:00  ypilchen
 * Rajout de delete d' ARM_Vector
 * voir rubriques : MA : delete before !!
 *
 * Revision 1.5  1999/02/22 09:52:28  nicolasm
 * Ajout test sur la validiti des contrats dans ZCFromMarketDataGen
 *
 * Revision 1.4  1999/02/18 18:43:36  nicolasm
 * Ajout Constructeur pour contrat echeance mensuelle
 *
 * Revision 1.3  1999/01/28 17:15:16  nicolasm
 * Ajout paralleleshift
 *
 */

/*--------------------------------------------------------------------------*/
/*                                                                          */
/* FILE         : zerocurv.cpp                                              */
/*                                                                          */
/* DESCRIPTION  : This file implements the ARM_ZeroCurve class, a class for */
/*                managing zero curves                                      */
/* DATE         : Thu Aug  1 1996                                           */
/*                                                                          */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/*---- System Include ----*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>


/*---- Application Include ----*/

#include "firsttoinc.h"
#include "dates.h"
#include "fromto.h"
#include "interpol.h"
#include "zerocurv.h"
#include "expt.h"
//#include "ycmodel.h"
//#include "y2cmodel.h"
#include "containr.h"
#include "security.h"
//#include "bond.h"

//#include "swap.h"

#include "zeroint.h"


void intToStr(int input, char* str, int dummy)
{
    sprintf(str, "%d", input);
}

#ifndef unix 
#define ITOA _itoa
#else
#define ITOA intToStr
#endif


extern ARM_INDEX_TYPE GetDefaultIndexFromCurrency(char* ccy);



int GOTO_END_OF_MONTH = 1;



void ARM_ZeroCurve::Init(void)
{
    //  set default values of other variables

    itsDateTerms = NULL;

    itsYearTerms = NULL;

    itsBPShifts  = NULL;

    itsZeroRates = NULL;


    itsForwardDates = NULL;

    itsForwardRates = NULL;

    itsForwardDF    = NULL;

 
    itsDiscountFactors = NULL;

    itsParameters = NULL;

    itsBPShift = 0.0;

    itsBucketStartPeriod = 0.0;

    itsBucketEndPeriod = 0.0;

    itsCurrency = new ARM_Currency(ARM_DEFAULT_COUNTRY);

    itsDFSpotFwd = 0.0;

   // Initialisation des markets data

   itsMktData = NULL;

   itsSmoothingMethod = "NOTHING";

   itsFixDayCount = KNOBASE;
}


void ARM_ZeroCurve::GenerateYearTermsAndRates(void)
{
    int i, sz;


    sz = 4+9+15; // 1D, 1M, 3M, 6M, 1Y, 1.5Y, 2, ...,4.5Y, 5Y, 6Y,...,20Y

    if ( itsYearTerms == NULL )
    {
       itsYearTerms = new ARM_Vector(sz);

       if (itsZeroRates)
       {
          delete itsZeroRates;
          itsZeroRates = NULL;
       }

       itsZeroRates = new ARM_Vector(sz);

       if (itsDiscountFactors)
       {
          delete itsDiscountFactors;
          itsDiscountFactors = NULL;
       }

       itsDiscountFactors = new ARM_Vector(sz);

       (*itsYearTerms)[0] = 1.0/365.0;

       (*itsYearTerms)[1] = 30.0/365.0; 

       (*itsYearTerms)[2] = 90.0/365.0; 

       (*itsYearTerms)[3] = 120.0/365.0;


       (*itsYearTerms)[4] = 1;

       for ( i = 5; i < 13; i++)
       {
           (*itsYearTerms)[i] = (*itsYearTerms)[i-1]+0.5;
       }

       for ( i = 13; i < sz; i++)
       {
           (*itsYearTerms)[i] = (*itsYearTerms)[i-1]+1.0;
       }

       for ( i = 0; i < sz; i++)
       {
           (*itsZeroRates)[i] = this->DiscountYield((*itsYearTerms)[i]);

           (*itsDiscountFactors)[i] = this->DiscountPrice((*itsYearTerms)[i]);
       }
    }
}



void ARM_ZeroCurve::GenerateDateTerms(void)
{
     int i, sz;



    if (itsDateTerms)
    {
       delete itsDateTerms;

       itsDateTerms = NULL;
    }

    sz = itsYearTerms->GetSize();

    itsDateTerms = new ARM_Vector(sz);

    for ( i = 0; i < sz; i++)
    {
        (*itsDateTerms)[i] = (*itsYearTerms)[i]*365.0;
    }
}



void ARM_ZeroCurve::GenerateDiscountFactors(int compMeth)
{
    int i, sz;

    if (itsDiscountFactors)
    {
       delete itsDiscountFactors;

       itsDiscountFactors = NULL;
    }

    sz = itsYearTerms->GetSize();

    itsDiscountFactors = new ARM_Vector(sz);
    double* itsDiscountFactorsElt = itsDiscountFactors->GetElt();
    double* itsYearTermsElt       = itsYearTerms->GetElt();
    double* itsZeroRatesElt       = itsZeroRates->GetElt();

    if ( compMeth == 0 )
    {
       for ( i = 0; i < sz; i++)
       {
           itsDiscountFactorsElt[i] = 
                exp(-0.01*itsYearTermsElt[i]*itsZeroRatesElt[i]);
       }
    }

    if ( compMeth == -1 )
    {
       for ( i = 0; i < sz; i++)
       {
           itsDiscountFactorsElt[i] = 
                1.0/(1.0 + 0.01*itsYearTermsElt[i]*itsZeroRatesElt[i]);
       }
    }
    
    if ( compMeth > 0 )
    {
       for ( i = 0; i < sz; i++)
       {
           itsDiscountFactorsElt[i] = 
                pow(1.0 + 0.01*itsZeroRatesElt[i]/compMeth, 
                -itsYearTermsElt[i]*compMeth);
       }
    }
}



void ARM_ZeroCurve::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

	fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
    fprintf(fOut, "\n\n\t          >>>>>>>>>>>>> Zero Curve <<<<<<<<<<<<\n\n");

    ARM_AbstractMarketClass::View(id, fOut);

	fprintf(fOut, "\n\n \tSmoothing method : %s \n\n", itsSmoothingMethod.c_str());
 
	fprintf(fOut, "\n\n \tFix DayCount : %d \n\n", itsFixDayCount);

    if (itsMktData)
    {
       fprintf(fOut, "\n\t ------------> Market Data \n");
       int szMktData = itsMktData->itsMktValue->GetSize();

       fprintf(fOut, "\n \tPlot\tData \n\n");
       for (int i = 0; i < szMktData; i++)
       {
           fprintf(fOut, "\t%s\t %.5lf\n", 
                        itsMktData->itsMktTerms[i],
                        itsMktData->itsMktValue->Elt(i));
       }
    }

    int sz = itsYearTerms->GetSize();
    char d[20];
    int  i;
 
    fprintf(fOut, "\n\nDate\t\tDays\tRate\t\tDiscount \n\n");
   
    for (i = 0; i < sz; i++)
    {
        ARM_Date date, tmpDate;
        double   rate, dF;
        int      days;
       

        if ( (*itsYearTerms)[i] == 1000.0 )
        {
            break;
        }

        days = int(NEAREST_INT((*itsDateTerms)[i]));

        tmpDate = itsAsOfDate;

        date = tmpDate.AddDays(days);

        date.JulianToStrDate(d);
 
        rate = (*itsZeroRates)[i];

        if (itsDiscountFactors)
           dF = (*itsDiscountFactors)[i];
        else
           dF = 0.0;
 
        if ( (*itsYearTerms)[i] != 1000.0 )
        {
           fprintf(fOut, "%s\t%5d\t %.5lf \t %.5lf\n", d, days, rate, dF);
        }
    }

    if (itsForwardRates)
    {
       if (itsForwardRates)
       {
          ViewForward(id, fOut);
       }
    } 

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


 
/*--------------------------------------------------------------------------*/
/* Constructor                                                              */
/*--------------------------------------------------------------------------*/

ARM_ZeroCurve::ARM_ZeroCurve(ARM_Date& asOf, ARM_Currency* ccy)
{
    SetName(ARM_ZERO_CURVE);

    Init();


    // Set default values of other variables

    itsAsOfDate = asOf;

    SetCurrencyUnit(ccy);

    itsBPShift = 0.0;

    itsBucketStartPeriod = 0.0;

    itsBucketEndPeriod = 0.0;
}



ARM_ZeroCurve::ARM_ZeroCurve(ARM_Date& asOf, 
                             ARM_CRV_TERMS& terms, 
                             ARM_Vector* mktData, 
                             int MMVsFut, 
                             int SwapVsFut, 
                             int raw, 
                             int Cont_Lin, 
                             ARM_Currency* ccy,
							 int swapFrqId,
							 int fixDayCount)
{
    SetName(ARM_ZERO_CURVE);

    Init();


    itsAsOfDate = asOf;

    // Set variables Market Data

   itsMktData = new ARM_MarketData(terms, mktData, 
                                   MMVsFut, SwapVsFut, raw, Cont_Lin, 0);

    SetCurrencyUnit(ccy);

	SetFixDayCount(fixDayCount);

    ZCFromMarketRates(terms, mktData, MMVsFut, SwapVsFut, raw, Cont_Lin, ccy, swapFrqId, fixDayCount);
}



// Constructeur swap plus general
 
ARM_ZeroCurve::ARM_ZeroCurve(ARM_CRV_TERMS& terms,
                             ARM_Date& asOf,
                             ARM_Vector* mktData,
                             int MMVsFut,
                             int SwapVsFut,
                             int raw,
                             int Cont_Lin,
                             ARM_Currency* ccy)
{
    SetName(ARM_ZERO_CURVE);

    Init();
 
 
    itsAsOfDate = asOf;
 
    // Set variables
 
   itsMktData = new ARM_MarketData(terms, mktData, MMVsFut, 
                                   SwapVsFut, raw, Cont_Lin, 1);

    SetCurrencyUnit(ccy);
 
    ZCFromMarketRatesGen(terms, mktData, MMVsFut, 
                         SwapVsFut, raw, Cont_Lin, ccy);
}



ARM_ZeroCurve::ARM_ZeroCurve(ARM_Date& asOf, 
                             ARM_CRV_TERMS& terms, 
                             ARM_Vector* mktData, 
                             ARM_Container* bonds, 
                             ARM_Vector* bdYields,
                             int MMVsFut, 
                             ARM_Currency* ccy)
{
    SetName(ARM_ZERO_CURVE);

    Init();


    itsAsOfDate = asOf;

    // Set variables

    SetCurrencyUnit(ccy);

    ZCFromMarketRates(terms, mktData, bonds, bdYields, MMVsFut, ccy);
}



/*--------------------------------------------------------------------------*/
/* Constructor(Copy)                                                        */
/*--------------------------------------------------------------------------*/

ARM_ZeroCurve::ARM_ZeroCurve(const ARM_ZeroCurve& zeroCurve) 
              : ARM_AbstractMarketClass(zeroCurve)
{
    Init();

    BitwiseCopy(&zeroCurve);
}



ARM_ZeroCurve::~ARM_ZeroCurve(void)
{
    if (itsDateTerms)
       delete itsDateTerms;

    if (itsYearTerms)
       delete itsYearTerms;

    if (itsBPShifts)
       delete itsBPShifts;

    if (itsZeroRates)
       delete itsZeroRates;

    if (itsDiscountFactors)
       delete itsDiscountFactors;

    if (itsParameters)
       delete itsParameters;

    if (itsCurrency && itsCurrency != ARM_DEFAULT_CURRENCY)
        delete itsCurrency;

   if (itsMktData)
      delete itsMktData;
   itsMktData = NULL;
}



/*--------------------------------------------------------------------------*/
/* Assignment operator                                                      */
/*--------------------------------------------------------------------------*/

ARM_ZeroCurve& ARM_ZeroCurve::operator = (const ARM_ZeroCurve& zeroCurve)
{
    (*this).ARM_AbstractMarketClass::operator = (zeroCurve);

    BitwiseCopy(&zeroCurve);

    return(*this);
}



void ARM_ZeroCurve::SetCurrencyUnit(ARM_Currency* ccy) 
{
    if (itsCurrency && itsCurrency != ARM_DEFAULT_CURRENCY)
       delete itsCurrency;

    if (ccy)
       itsCurrency = (ARM_Currency *) ccy->Clone();
    else 
       itsCurrency = NULL;

    SetStrCurrency(itsCurrency->GetCcyName());
}


 
/*--------------------------------------------------------------------------*/
/* Returns discount price at settlement of zero coupon maturing at maturity */
/* Returns K_HUGE_DOUBLE if failed.                                         */
/*--------------------------------------------------------------------------*/
double ARM_ZeroCurve::DiscountPrice(ARM_Date& maturity)
{
    double yearTerm;

    double res = 0.0;
    
    yearTerm = (maturity.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;

    res = DiscountFunction(yearTerm);

    return(res);
}


 
/*--------------------------------------------------------------------------*/
/* Returns the discount price for yearTerm  maturity                        */
/* Returns K_HUGE_DOUBLE if failed.                                         */
/*--------------------------------------------------------------------------*/

double ARM_ZeroCurve::DiscountPrice(double yearTerm)
{
    return(DiscountFunction(yearTerm));
}



double ARM_ZeroCurve::DiscountYield(ARM_Date &maturity, int compMeth)
{
    double yearTerm, z, y;


    yearTerm = (maturity.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;
    
    if (yearTerm < K_DOUBLE_TOL) yearTerm = YC_INSTANTANEOUS_MAT;
    
    z = DiscountFunction(yearTerm);

    if ( compMeth == 0 ) 
    {
       y = -100.0*log(z)/yearTerm;
    }
    else if(compMeth == -1) 
    {
       y = 100.0*(1.0/z-1.0)/yearTerm;
    }   
    else if ( compMeth > 0 )
    {
       y = 100.0*compMeth*(pow(1.0/z, 1.0/(yearTerm*compMeth))-1.0);
    }   

    return(y);
}



double ARM_ZeroCurve::DiscountYield(double yearTerm, int compMeth)
{
    double z, y;


    if (yearTerm < K_DOUBLE_TOL) yearTerm = YC_INSTANTANEOUS_MAT;

    z = DiscountFunction(yearTerm);

    if ( compMeth == 0 ) 
    {
       y = -100.0*log(z)/yearTerm;
    }
    else if ( compMeth == -1 ) 
    {
       y = 100.0*(1.0/z-1.0)/yearTerm;
    }   
    else if ( compMeth > 0 )
    {
       y = 100.0*compMeth*( pow(1.0/z, 1.0/(yearTerm*compMeth))-1.0);
    }   

    return(y);
}



double ARM_ZeroCurve::ForwardPrice(ARM_Date& maturity1, ARM_Date& maturity2)
{
    double yearTerm1, yearTerm2, z1, z2, f;



    yearTerm1 = (maturity1.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;
    yearTerm2 = (maturity2.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;

    if ( yearTerm2 < yearTerm1)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Zero maturity date must be greater than forward date");
    }

    z1 = DiscountFunction(yearTerm1);
    z2 = DiscountFunction(yearTerm2);
    
    f = z2/z1;

    return(f);
}



double ARM_ZeroCurve::ForwardPrice(double yearTerm1, double yearTerm2)
{
    if ( yearTerm2 < yearTerm1 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Zero maturity date must be greater than forward date");
    }

    double z1, z2, f;

    z1 = DiscountFunction(yearTerm1);
    z2 = DiscountFunction(yearTerm2);
    
    f = z2/z1;

    return(f);
}



/*---------------------------------------------------------------------------*
      !!! NOTE : The ForwardYield methods are given for convenience :
      the basis used to compute the rate term is ACTUAL/365. 
      So take care and use the Yield Curve Model forward yield method
      which compute the rate term with the right basis of the security
*---------------------------------------------------------------------------*/
        
double ARM_ZeroCurve::ForwardYield(ARM_Date& maturity1, ARM_Date& maturity2, 
                                   int compMeth,
                                   int adjFwd)
{
    double yearTerm1, yearTerm2, z1, z2, f, term;


	if (maturity1 == maturity2)
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Yield maturity is equal to forward date");
	}

	if (maturity1 > maturity2)
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Yield maturity date must be greater than forward date");
	}

    yearTerm1 = (maturity1.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;
    yearTerm2 = (maturity2.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;

    z1 = DiscountFunction(yearTerm1);
    z2 = DiscountFunction(yearTerm2);

    if (!(adjFwd))
    {
       term = (maturity2.GetJulian()-maturity1.GetJulian())/K_YEAR_LEN;
    }
    else
    {
       term = CountYears(GetCurrencyUnit()->GetLiborIndexDayCount(),
                         maturity1.GetJulian(),
                         maturity2.GetJulian());
    }

    if ( compMeth == 0 ) 
    {
       f = -100.0*log(z2/z1)/term;
    }
    else if ( compMeth == -1 ) 
    {
       f = 100.0*(z1/z2-1.0)/term;
    }
    else if ( compMeth > 0 )
    {
	   f = 100.0*compMeth*(pow(z1/z2, 1.0/(compMeth*term))-1.0);
    }

    return(f);
}



double ARM_ZeroCurve::ForwardYieldWithDayCount(ARM_Date& maturity1,
                                               ARM_Date& maturity2, 
                                               int compMeth,
                                               int dayCount)
{
   double yearTerm1, yearTerm2, z1, z2, f, term;


	if ( maturity1 == maturity2 )
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Yield maturity is equal to forward date");
	}

	if ( maturity1 > maturity2 )
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Yield maturity date must be greater than forward date");
	}

    yearTerm1 = (maturity1.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;
    yearTerm2 = (maturity2.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;

    z1 = DiscountFunction(yearTerm1);
    z2 = DiscountFunction(yearTerm2);

   
    term = CountYears(dayCount,
                      maturity1.GetJulian(),
                      maturity2.GetJulian());

    if ( compMeth == 0 ) 
    {
       f = -100.0*log(z2/z1)/term;
    }
    else if ( compMeth == -1 ) 
    {
       f = 100.0*(z1/z2-1.0)/term;
    }
    else if ( compMeth > 0 )
    {
	   f = 100.0*compMeth*(pow(z1/z2, 1.0/(compMeth*term))-1.0);
    }

    return(f);
}



double ARM_ZeroCurve::ForwardYield(double yearTerm1, double yearTerm2,
                                   int compMeth)
{
    if ( yearTerm2 <= yearTerm1)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Zero maturity date must be greater than forward date");
    }

    double    z1, z2, f;

    z1 = DiscountFunction(yearTerm1);
    z2 = DiscountFunction(yearTerm2);
    
    if ( compMeth == 0 ) 
    {
       f = -100.0*log(z2/z1)/(yearTerm2 - yearTerm1);
    }
    else if ( compMeth == -1 ) 
    {
       f = 100.0*(z1/z2-1.0)/(yearTerm2-yearTerm1);
    }
    else if ( compMeth > 0 )
    {
       f = 100.0*compMeth
            *( pow(z1/z2, 1.0/(compMeth*(yearTerm2-yearTerm1)))-1.0);

    }

    return(f);
}



double ARM_ZeroCurve::CalcNumericalObjectSignature(void)
{
    int sz = itsYearTerms->GetSize();

    double signature = 0.0;

    for (int i = 0; i < sz; i++)
    {
        signature += (itsYearTerms->Elt(i)*itsZeroRates->Elt(i));
    }

    signature += GetAsOfDate().GetJulian();

    return(signature);
}





/*---------------------------------------------------------------------------*

    Return the continuous spot (or instantaneous) forward rate defined as 

    f(t) = -d(ln(Z(0,t))/dt, the limit of forward rate F(t, t+dt) of 

    maturity dt when dt->0

*---------------------------------------------------------------------------*/

double ARM_ZeroCurve::SpotForwardRate(ARM_Date &maturity)
{
    double yearTerm, z, zp, f;


    yearTerm = (maturity.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;

    if ( yearTerm < K_DOUBLE_TOL )
    {
       f = DiscountYield(YC_INSTANTANEOUS_MAT);

       return(f);
    }
    else if (itsDFSpotFwd == 0.0)
    {
       z = DiscountFunction(yearTerm);

       zp = D1DiscountFunction(yearTerm);

       f = -100.0*zp/z;

       return(f);
    }
    else
    {
       f = DiscountPrice(yearTerm)-
                  DiscountPrice(yearTerm+itsDFSpotFwd);
 
       f /= itsDFSpotFwd;
 
       f /= DiscountPrice(yearTerm);
 
       f  *= 100.0;
 
       return(f);
    }


    return(f);
}



double ARM_ZeroCurve::SpotForwardRate(double yearTerm)
{
    double z, zp, f;


    if ( yearTerm < K_DOUBLE_TOL )
    {
       f = DiscountYield(YC_INSTANTANEOUS_MAT);

       return(f);
    }
    else if (itsDFSpotFwd == 0.0)
    {
       z = DiscountFunction(yearTerm);

       zp = D1DiscountFunction(yearTerm);

       f = -100.0*zp/z;

       return(f);
    }
    else 
    {
       f = DiscountPrice(yearTerm)-
                  DiscountPrice(yearTerm+itsDFSpotFwd);
 
       f /= itsDFSpotFwd;
 
       f /= DiscountPrice(yearTerm);
 
       f  *= 100.0;
 
       return(f);
    }

}



/*--------------------------------------------------------------------------*
  Returns the BPV at settlement of zero coupon maturing at maturity 
*---------------------------------------------------------------------------*/

double ARM_ZeroCurve::DiscountPriceBPV(double yearTerm)
{
    double bpv = DiscountFunction(yearTerm);


    SetParallelShift(1.0);

    bpv -= DiscountFunction(yearTerm);

    GetRidOfShift();

    return(10000.0*bpv);
}



double ARM_ZeroCurve::DiscountPriceBPV(ARM_Date& maturity)
{
    double    yearTerm;
    
    yearTerm = (maturity.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;

    return(DiscountPriceBPV(yearTerm));
}



double ARM_ZeroCurve::ForwardPriceBPV(double fwdTerm, double zeroTerm)
{
    double bpv = ForwardPrice(fwdTerm, zeroTerm);

    SetParallelShift(1.0);

    bpv -= ForwardPrice(fwdTerm, zeroTerm);

    GetRidOfShift();

    return(10000.0*bpv);
}



double ARM_ZeroCurve::ForwardPriceBPV(ARM_Date& fwdDate, ARM_Date& zeroMat)
{
    double yearTerm1, yearTerm2;//, z1, z2, f;


    yearTerm1 = (fwdDate.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;
    yearTerm2 = (zeroMat.GetJulian()-itsAsOfDate.GetJulian())/K_YEAR_LEN;


    return(ForwardPriceBPV(yearTerm1, yearTerm2));
}



/******  THE FOLLOWING METHODS MUST BE OVERRIDEN BY SUBCLASSES  *******/



    
/*--------------------------------------------------------------------------*/
/* Returns the discount price with maturity yearTerm years from start date. */
/*   This method MUST be overriden by subclasses.                           */
/*--------------------------------------------------------------------------*/

double ARM_ZeroCurve::DiscountFunction(double yearTerm)
{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <DiscountFunction> method");
    return(K_HUGE_DOUBLE);
}


    
/*--------------------------------------------------------------------------*/
/* Returns the d(discount price) / d(yearTerm) with maturity yearTerm years */
/* from start date.                                                         */
/*   This method MUST be overriden by subclasses.                           */
/*--------------------------------------------------------------------------*/

double ARM_ZeroCurve::D1DiscountFunction(double yearTerm)
{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <D1DiscountFunction> method");
    return(K_HUGE_DOUBLE);
}



/*--------------------------------------------------------------------------*/
/* Returns the d(discount price) / d(yearTerm) with maturity yearTerm years */
/*  from start date.                                                        */
/*   This method MUST be overriden by subclasses.                           */
/*--------------------------------------------------------------------------*/

double ARM_ZeroCurve::D2DiscountFunction(double yearTerm)
{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <D2DiscountFunction> method");

    return(K_HUGE_DOUBLE);
}



double EvalEcart(ARM_Matrix* data, ARM_Vector* Params, void** parameters)
            
{ 
    ARM_Container* Assets;
    ARM_YCModel* ycModel;
    ARM_ZeroCurve* ZC;


    Assets = (ARM_Container *) parameters[1];
    ycModel = (ARM_YCModel * ) parameters[2];
    ZC = (ARM_ZeroCurve *) parameters[3];

    ARM_Security* sec;
    double resultat = 0.0;
    int i;
   
    ZC->SetParameters((ARM_Vector *) Params->Clone());
    
    ycModel->SetZeroCurve(ZC);

    sec = (ARM_Security *) Assets->Start();

    i = 0;

    double price;
    
    while (sec)
    {
        sec->SetModel(ycModel);

        price = sec->ComputePrice();

        data->Elt(i,2) = data->Elt(i,0)*(price-data->Elt(i,1));        

        resultat +=pow(price-data->Elt(i,2), 2);

        sec = (ARM_Security *) Assets->Next();

        i++;
    }

    return(resultat);
}



double ARM_ZeroCurve::FitToMarketPrices(ARM_Container* Assets,
                                        ARM_Matrix* data,
                                        int* CST,double* OP)

{    
    void *Divers[4];
    ARM_YCModel* YCMod = new ARM_YCModel(this);
    ARM_Matrix *Hessien;    
    ARM_Vector *p;
    int i, rc;

    int nbAssets=Assets->GetSize();
    
    
    Divers[0] = (void *) &nbAssets;
    Divers[1] = (void *) Assets;
    Divers[2] = (void *) YCMod;
    Divers[3] = (void *) this;

    
    p = itsParameters;
    
    Hessien = new ARM_Matrix (p->GetSize(), p->GetSize(), 0.0);

/*    rc = dfpmin(p,Hessien,&iter,&fret,parameters,data,EvalEcart);*/

    if ( CST == NULL )
    {
        CST= (int*) malloc(20*sizeof(int));

        for (i = 0; i < 20; i++) 
            CST[i] = 0;
    }

    for (i = 0; i < nbAssets; i++)
    {
        data->Elt(i,0) = sqrt(data->Elt(i,0));
    }

    rc = leastsq(p, CST,Divers, data, EvalEcart,OP);
    
    ARM_Security* sec;
    double resultat=0.0;
       
    SetParameters((ARM_Vector *) p->Clone());
        
    YCMod->SetZeroCurve(this);

    sec = (ARM_Security *) Assets->Start();
    
    i = 0;
        
    while (sec)
    {
        sec->SetModel(YCMod);
        
        data->Elt(i,2) = sec->ComputePrice() - data->Elt(i,1);

        resultat += SQR(data->Elt(i,2))*data->Elt(i,0);

        sec = (ARM_Security *) Assets->Next();

        i++;
    }
    
    resultat = pow(resultat,.5)/(i-1);

    delete YCMod;

    GenerateFields();

    return(resultat);

}        


void ARM_ZeroCurve::ZCFromMarketRates(ARM_CRV_TERMS& Terms, 
                                      ARM_Vector* data, 
                                      int MMVsFut, int SwapVsFut,
                                      int raw,
                                      int Cont_Lin,
                                      ARM_Currency* ccy,
									  int swapFrqId,
									  int fixDayCount)
{
    int SizeM=0, SizeF=1, SizeS=1, i=0, Nb=0;
    int Tx_Spot=0, Tx_D=0;
    ARM_Date matDate, spotDate, swapStartDate, indexSwapDate, startDate, prevDate, limite_AUD_DateQ;

	int firstSwapIndex = 1000;

	int compteur = 1;

    ARM_Date matuDateFormat;
 
    int CurrentMonth;

    ARM_Matrix MMRates(100, 3, 0.0);
    ARM_Matrix Futures(300, 3, 0.0); // regroupe Future, FRA, Serial
    ARM_Matrix SwapRates(1000, 2, MINDOUBLE);

    SetCurrencyUnit(ccy);

    char ccyName[10];
	ccy->CalcFloatPayCal(ccyName);

    int fwdRule = ccy->GetFwdRule();
    int spotDays = ccy->GetSpotDays();

    int fxPayFreq;
	if (swapFrqId == K_DEF_FREQ)
		fxPayFreq = ccy->GetFixedPayFreq();
	else
		fxPayFreq = swapFrqId;

    int fxPayFreqVal = fxPayFreq; //variable pour currency australie
	int lastAUDtriNb = 0; 
	int indexAUD = 0; 

	int swapFixDayCount;
	if (fixDayCount == KNOBASE)
		swapFixDayCount = ccy->GetFixedDayCount();
	else
		swapFixDayCount = fixDayCount;

    char matu, matu2;
    int year;
    int month = 0;

    int lastMMTerm;

    spotDate = itsAsOfDate;

	if (spotDays > 0)
		spotDate.NextBusinessDay(spotDays, ccyName);
	else
		Tx_Spot = 1;

    swapStartDate = itsAsOfDate;

	swapStartDate.NextBusinessDay(spotDays, ccyName);

	limite_AUD_DateQ = swapStartDate;

	limite_AUD_DateQ.AddYears(3); // Dans la courbe AUD les swaps <= 3Y sont trimestriels
    
    CurrentMonth = spotDate.GetMonth();
    
    matu = 'X';
    matu2 = 'X';
    i = 0;

	prevDate = itsAsOfDate;

    while ( Terms[i][0] != 'X' )
    {   
        if ( data->Elt(i) < MINDOUBLE ) 
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Non zero rates are required for building zero curve");
        }

        // si 1er caract alphanum -> Futures
        if (isalpha(Terms[i][0]))
        {
            GetMonthYearFromExpiryDate(Terms[i], &month, &year);
            matu='Z';
        }
        else // sinon MM, Serial, FRA ou Swap
        {
           if ((( Terms[i][2] == '.' ) || ( Terms[i][2] == '/' ))
               &&
               (( Terms[i][5] == '.' ) || ( Terms[i][5] == '/' ))
              ) // A date in a MM, Serial or FRA segment
           {
              if ( strlen(Terms[i]) == 10 ) // DD.MM.YYYY
              {
                  char c1, c2;

                  int d;
                  int m;
                  int y;

                  sscanf(Terms[i], "%d%c%d%c%4d", &d, &c1, &m, &c2, &y);

                  ARM_Date tmpDate(d, m, y);

                  matuDateFormat = tmpDate;

                  matu = 'S';
              }
              else if ( strlen(Terms[i]) == 8 ) // DD.MM.YY
              {
                  char c1, c2;

                  int d;
                  int m;
                  int y;

                  sscanf(Terms[i], "%d%c%d%c%2d", &d, &c1, &m, &c2, &y);

                  y = 2000+y;
    
                  ARM_Date tmpDate(d, m, y);

                  matuDateFormat = tmpDate;

                  matu = 'S';
              }
			  if (strlen (Terms[i]) > 10 ) // FRA's or Serial case: MM/DD/YY-MM/DD/YY
			  {
                  char c;

                  int d1,d2;
                  int m1,m2;
                  int y1,y2;

                  if ( strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0 )
				  {
                     sscanf(Terms[i], "%d%c%d%c%2d%c%d%c%d%c%2d", &m1, &c, &d1,
                            &c, &y1, &c, &m2, &c, &d2, &c, &y2 );
				  }
                  else
				  {
                     sscanf(Terms[i], "%d%c%d%c%2d%c%d%c%d%c%2d", &d1, &c, &m1,
                            &c, &y1, &c, &d2, &c, &m2, &c, &y2 );
				  }

                  y1 = 2000+y1;
                  ARM_Date tmpDate(d1, m1, y1);
                  startDate = tmpDate.GetJulian();// chge format from MM/JJ/YY to julian

                  y2 = 2000+y2;
                  ARM_Date tmpDate2(d2, m2, y2);
                  matDate = tmpDate2.GetJulian();// chge format from MM/JJ/YY to julian

                  matu = 'R';

			  }
              else
              {
                 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                                 "Invalid maturity in curve inputs");
              }
           }
           else
           {
              //sscanf(Terms[i], "%d%c", &Nb, &matu); 
              matu2 = 'X';
			  sscanf(Terms[i], "%d%c%c", &Nb, &matu,&matu2);

              matu = toupper(matu);
           }
        }

        // Conversion des Termes char -> Termes double

        if ( matu == 'D' ) // Ex : "1D"
        {    
            matDate = itsAsOfDate;
            
            if ( spotDays < 1 || Nb==1)
            {
               matDate.NextBusinessDay(Nb, ccyName);
               MMRates.Elt(SizeM, 1) = matDate.GetJulian();
            }
            
			matDate = itsAsOfDate;
			matDate.NextBusinessDay(Nb, ccyName);


			if (matDate <= spotDate)
			{
				startDate = prevDate;
				prevDate = matDate;
				if (matDate == spotDate)
					Tx_Spot = 1;
			}
			else
			{
				if (Tx_Spot == 1)
					startDate = spotDate;
				else
					startDate = itsAsOfDate;
			}

			
			MMRates.Elt(SizeM, 0) = data->Elt(i);

			MMRates.Elt(SizeM, 1) = matDate.GetJulian();

			MMRates.Elt(SizeM, 2) = startDate.GetJulian();

			Tx_D = 1;

            SizeM++;
        }
        else if ( matu == 'W' )  
        {   //  Ex : "1W"    

			if (Tx_D == 1) // supposed that MM segment sorted in input
			{
				startDate = spotDate;
			}
			else
			{
				startDate = itsAsOfDate;
			}
            
			matDate = startDate;

            matDate.AddDays(Nb*7);

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);


			MMRates.Elt(SizeM, 0) = data->Elt(i);
            
            MMRates.Elt(SizeM, 1) = matDate.GetJulian();

			MMRates.Elt(SizeM, 2) = startDate.GetJulian();

            SizeM++;
        }
        else if ((( matu == 'M' ) && ( Nb <= 12 )) 
                 ||
                 ( matu == 'S' )
                )
        {   //  Ex : "9M"    

            if ( matu == 'M' ) // cas Money market classique
            {
               matDate = spotDate;
            
               matDate.AddMonths(Nb, GOTO_END_OF_MONTH);
            
               lastMMTerm = Nb;
            }
            else // cas date
            {
               matDate = matuDateFormat; 


               if (((matuDateFormat-itsAsOfDate)/30.0) < 1 )
               {
                  lastMMTerm = 1;
               }
               else
               {
                  lastMMTerm = int((matuDateFormat-itsAsOfDate)/30.0);
               }
            }

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);


			if ( (Tx_D == 1) && (matu == 'M') ) // supposed that MM segment sorted in input
			{
				startDate = spotDate;
			}
			else
			{
				startDate = itsAsOfDate;
			}

			if ( matu == 'M' ) // cas Money market classique
            {
               matDate = startDate;
            
               matDate.AddMonths(Nb, GOTO_END_OF_MONTH);

			   matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);
            }

			MMRates.Elt(SizeM, 0) = data->Elt(i);
            
            MMRates.Elt(SizeM, 1) = matDate.GetJulian();

			MMRates.Elt(SizeM, 2) = startDate.GetJulian();

            SizeM++;
        }
        else if (( matu == 'Y')  // ->implicitement ce sont des taux de swap
                 ||
                 (( matu == 'M' ) && ( Nb > 12 ))
                )
        {   
            //  Ex : "15Y"

            matDate = swapStartDate;;
            
            //matDate.AddMonths(12*Nb/fxPayFreq);
            if (( matu == 'M' ) && ( Nb > 12 ))
            {
               matDate.AddMonths(Nb);
            }
            else
            {
               matDate.AddYears(Nb);
            }

			if (raw == K_FORWARD)
			{
				SwapRates.Elt(SizeS, 0) = data->Elt(i);

				SwapRates.Elt(SizeS, 1) = matDate.GetJulian(); // End Swap unadjusted

				firstSwapIndex = 1;

				if ( (matDate.GetJulian() <= limite_AUD_DateQ.GetJulian()) &&
					 (strcmp((const char*)ccyName,"AUD") == 0) )
					lastAUDtriNb = SizeS;

				SizeS++;
			}
			else
			{

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

			//|->
			// on construire dynamiquement la partie de swap de la courbe 
			int flagMonthOrYear;// flagMonthOrYear=0 if matu='M',Nb>12;  flagMonthOrYear=1 if matu='Y'
			int flagAUDTri;

			if (( matu == 'M' ) && ( Nb > 12 )) // 18M ...
				flagMonthOrYear = 0;
			else
				flagMonthOrYear = 1;
			
			if(matu==matu2 || Nb<=3) //par default, 1Y,2Y,3Y sont trimestielle
				flagAUDTri=1;
			else
				flagAUDTri=0;

			ConstructDynamSwap(	Nb, flagMonthOrYear, data->Elt(i), flagAUDTri, //flagAUDTri==1 si matu==matu2 eg.1YY
								matDate, 
								lastMMTerm, 
								Futures, SizeF,
								ccy, 
								swapStartDate,
								&SwapRates,
								fxPayFreq,SizeS, 
								compteur,
								firstSwapIndex,
								lastAUDtriNb);	//lastAUDtriNb: last nb of years in the part of the AUD trimestre
			//<-

			}
			
        }
        // traitement du segment Futures
        else if ( month > 0 )
        {  
            matDate.ChangeDate(1, month, year);

			matDate.PiborDelivery(ccy->GetCcyName());

            // tester si NextDate < = > itsAsOfDate
            if ( matDate.GetJulian() < itsAsOfDate.GetJulian()+ccy->GetSpotDays())
            {
               char BUF[100];

               sprintf(BUF, 
                      "Futures contract expiring : %s  %d",Terms[i], year);

               throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, BUF);
            }


			startDate = matDate;

			matDate.ChangeDate(1, month, year);

			matDate.AddMonths(3);

			matDate.PiborDelivery(ccy->GetCcyName());

			Futures.Elt(SizeF, 0) = 100.0-data->Elt(i);
            
            Futures.Elt(SizeF, 2) = startDate.GetJulian();

			Futures.Elt(SizeF, 1) = matDate.GetJulian();

            SizeF++;
        }
		else if (matu == 'R') // traitement des FRA
		{
			
		    Futures.Elt(SizeF, 0) = data->Elt(i);
            
            Futures.Elt(SizeF, 2) = startDate.GetJulian();

			Futures.Elt(SizeF, 1) = matDate.GetJulian();

            SizeF++;
		}
        else  
           break; // unrecognizable format

        i++;
    }

    
    if ( SizeF == 1 )  // pas de contrats futures ->Forcer les options 
    {
       SwapVsFut = K_SWAP;

       MMVsFut = K_MM;
    }

	if( SizeS == 1 )
	{
		SwapVsFut = K_FUT;
	}
	
    // il faut le taux swap 1Y sauf si ce sont les futures qui dominent 
    if (( SwapRates.Elt(fxPayFreq, 0) < MINDOUBLE+0.0001 ) && ( SwapVsFut == 1 ))
    {
       char BUF[100];

       sprintf(BUF, "1 year Swap rate not found");
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, BUF);
    }

	SizeS = 0;
    for (i = 0; i < SwapRates.GetNumLines(); i++) 
    {
       if ( SwapRates.Elt(i, 1) > MINDOUBLE )
          SizeS = i; // sizeS = maturite max de la courbe
    }

	SwapRates.Elt(0, 1) = swapStartDate.GetJulian();
	
	if (raw == K_FORWARD)
		indexAUD = lastAUDtriNb;
	else
		indexAUD = lastAUDtriNb*4;// last index of the AUD trimestrielle

    
	// calcul des taux ZC implicites dans chaque segment du marche
    
	
	CptMMZeroRatesFWD(&MMRates, SizeM, Cont_Lin);


	if (SizeF != 1)
       CptAllForwardsZeroRates(MMRates, SizeM, &Futures, SizeF, 
                           SizeS,
                           SwapRates, 
                           MMVsFut, Cont_Lin);

	if ( SizeS > 1 && SizeF > 1 )
	{
		if(SwapRates.Elt(firstSwapIndex, 1)-Futures.Elt(SizeF-1,1) > (366.00/fxPayFreq+10) ) // more than 12/fxPayFreq months 
																							 // (on suppose that it is more than (366.0/fxPayFreq+10) jours)
			SwapVsFut = K_FUT;
	}


    if (SizeS > 0) 
	{
		if (raw != K_FORWARD)
			CptSwapZeroRates(MMRates, SizeM, &SwapRates, SizeS, fxPayFreq,
                     Futures, SizeF, MMVsFut, SwapVsFut, raw, Cont_Lin,firstSwapIndex,indexAUD, swapFixDayCount); 
		else
			CptSwapZeroRatesFWD(MMRates, SizeM, &SwapRates, SizeS, fxPayFreq,
                     Futures, SizeF, MMVsFut, SwapVsFut, raw, Cont_Lin,firstSwapIndex,indexAUD, 3, K_UNADJUSTED, swapFixDayCount); //AR

	}
    // Preparation du tableau des taux ZC correspondant a generateCurve

    int zcSize = 0;

	if (strcmp((const char*)ccyName,"AUD") == 0)
	{
		fxPayFreqVal = 4;	
	}
	else
	{
		fxPayFreqVal = fxPayFreq;
	}

	if( SizeF>1 && SizeS>0 ) // manage the chevauchement between the Futures and Swaps
	{
		if( SwapVsFut==K_FUT || 
			( (raw == K_PAR) && ( SwapRates.Elt(firstSwapIndex, 1) - Futures.Elt(SizeF-1,1) > (366.00/fxPayFreqVal+10)) ))
		{
			while( SwapRates.Elt(firstSwapIndex, 1) - Futures.Elt(SizeF-1,1) > (366.00/fxPayFreqVal+10) )
			{
				firstSwapIndex--;
			}
		}
	}
	
	if( SizeF==1 && SizeS>0 && strcmp((const char*)ccyName,"AUD") == 0 && lastAUDtriNb!=0)
	{
		if((raw == K_PAR) && ( SwapRates.Elt(firstSwapIndex, 1) - MMRates.Elt(SizeM-1,1) > 366.00/fxPayFreqVal+10) ) 
		{
			while( SwapRates.Elt(firstSwapIndex, 1) - MMRates.Elt(SizeM-1,1) > 366.00/fxPayFreqVal+10 )
			{
				firstSwapIndex--;
			}			
		}
	}

	//manage the chevauchement between the MMRates and Futures
	int firstFutureIndex = 0;
	if(MMVsFut == K_MM)
	{
		firstFutureIndex = 1;
	}
	else
	{
		firstFutureIndex = 0;
	}

    for (i = 0; i < SizeM; i++)
    {
        if (MMVsFut != K_FUT) 
           zcSize++;
        else if ( MMRates.Elt(i, 1) < Futures.Elt(firstFutureIndex, 1) ) //else if ( MMRates.Elt(i, 1) < Futures.Elt(0, 1) ) 
           zcSize++;
    }

    for (i = firstFutureIndex; i < SizeF; i++)//for (i = 0; i < SizeF; i++)
    {
        if ( (MMVsFut != K_FUT) && 
             (MMRates.Elt(SizeM-1, 1) < Futures.Elt(i, 1)) )
        {
            if (SwapVsFut == K_FUT) 
               zcSize++;
            else if ( Futures.Elt(i, 1) < SwapRates.Elt(firstSwapIndex, 1) )  
			   zcSize++;
        }

        if  (MMVsFut == K_FUT) 
        {
            if (SwapVsFut == K_FUT) 
               zcSize++;

            else if ( Futures.Elt(i, 1) < SwapRates.Elt(firstSwapIndex, 1) )  
               zcSize++;
        }
    }

    for (i = firstSwapIndex; i <= SizeS; i++)
    {
        if (SwapVsFut != K_FUT) 
           zcSize++;

        else if ( SwapRates.Elt(i, 1) > Futures.Elt(SizeF-1, 1) ) 
           zcSize++;
    }

    // MA : delete before !!

    if (itsDateTerms)
       delete itsDateTerms;

    if (itsYearTerms)
       delete itsYearTerms;

    if (itsZeroRates)
       delete itsZeroRates;

    if (itsDiscountFactors)
       delete itsDiscountFactors;
    
    itsDateTerms = new ARM_Vector (zcSize, 0.0);
    itsYearTerms = new ARM_Vector (zcSize, 0.0);
    itsZeroRates = new ARM_Vector (zcSize, 0.0);
    itsDiscountFactors = new ARM_Vector(zcSize, 0.0);

    zcSize=0;

    for (i = 0; i < SizeM; i++)
    {
        if (MMVsFut != K_FUT) 
        {
            itsDateTerms->Elt(zcSize) = 
                MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();     
            
            itsDiscountFactors->Elt(zcSize) = MMRates.Elt(i, 0);     
            
            zcSize++;
        }
		else if ( MMRates.Elt(i, 1) < Futures.Elt(firstFutureIndex, 1) )//else if ( MMRates.Elt(i, 1) < Futures.Elt(0, 1) )
        {
            itsDateTerms->Elt(zcSize) = 
                MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();     
            
            itsDiscountFactors->Elt(zcSize) = MMRates.Elt(i, 0);     
            
            zcSize++;
        }
    }

    for (i = firstFutureIndex; i < SizeF; i++)//for (i = 0; i < SizeF; i++)
    {
        if  (MMVsFut == K_FUT) 
        {
            if (SwapVsFut == K_FUT) 
            {
                itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
                zcSize++;
            }
            else if ( Futures.Elt(i, 1) < SwapRates.Elt(firstSwapIndex, 1) )          
            {
                itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
                zcSize++;
            }
        }
        else if ( MMRates.Elt(SizeM-1, 1) < Futures.Elt(i, 1) )
        {
            if (SwapVsFut == K_FUT) 
            {
                itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
                zcSize++;
            }
            else if ( Futures.Elt(i, 1) < SwapRates.Elt(firstSwapIndex, 1) ) 
            {
                itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
                zcSize++;
            }
        }

    }

    for (i = firstSwapIndex; i <= SizeS; i++)
    {
        if (SwapVsFut != K_FUT)         
        {
            itsDateTerms->Elt(zcSize) = 
                SwapRates.Elt(i, 1) - itsAsOfDate.GetJulian();
            
            itsDiscountFactors->Elt(zcSize) = SwapRates.Elt(i, 0);     
            
            zcSize++;
        }
        else if ( SwapRates.Elt(i, 1) > Futures.Elt(SizeF-1, 1) )
        {
            itsDateTerms->Elt(zcSize) = 
                SwapRates.Elt(i, 1) - itsAsOfDate.GetJulian();
            
            itsDiscountFactors->Elt(zcSize) = SwapRates.Elt(i, 0);
            
            zcSize++;
        }
    }


    for (i=0; i < zcSize; i++)
    {
        itsYearTerms->Elt(i) = itsDateTerms->Elt(i)/365.0; 
        itsZeroRates->Elt(i) = -100.0*log(itsDiscountFactors->Elt(i))/
                               itsYearTerms->Elt(i);
    }

	/*	// Suppression des doublons
    ARM_Vector vNewDateTerms(zcSize, 0.0);
    ARM_Vector vNewYearTerms(zcSize, 0.0);
    ARM_Vector vNewZeroRates(zcSize, 0.0);
    ARM_Vector vNewDiscountFactors(zcSize, 0.0);
	compteur = 0;

    for (i=0; i < zcSize-1; i++)
	{
		if (itsDateTerms->Elt(i) != itsDateTerms->Elt(i+1))
		{
			vNewDateTerms.Elt(i-compteur) = itsDateTerms->Elt(i);
			vNewYearTerms.Elt(i-compteur) = itsYearTerms->Elt(i);
			vNewZeroRates.Elt(i-compteur) = itsZeroRates->Elt(i);
			vNewDiscountFactors.Elt(i-compteur) = itsDiscountFactors->Elt(i);

		}
		else
		{
			compteur ++;
		}
	}

	vNewDateTerms.Elt(zcSize-1-compteur) = itsDateTerms->Elt(zcSize-1);
	vNewYearTerms.Elt(zcSize-1-compteur) = itsYearTerms->Elt(zcSize-1);
	vNewZeroRates.Elt(zcSize-1-compteur) = itsZeroRates->Elt(zcSize-1);
	vNewDiscountFactors.Elt(zcSize-1-compteur) = itsDiscountFactors->Elt(zcSize-1);

	delete itsDateTerms;
	delete itsYearTerms;
	delete itsZeroRates;
	delete itsDiscountFactors;

	itsDateTerms = new ARM_Vector(&vNewDateTerms,0,zcSize-1-compteur);
	itsYearTerms = new ARM_Vector(&vNewYearTerms,0,zcSize-1-compteur);
	itsZeroRates = new ARM_Vector(&vNewZeroRates,0,zcSize-1-compteur);
	itsDiscountFactors = new ARM_Vector(&vNewDiscountFactors,0,zcSize-1-compteur);*/
}



void ARM_ZeroCurve::ZCFromMarketRates(ARM_CRV_TERMS& Terms, 
                                      ARM_Vector* data, 
                                      ARM_Container* bonds, 
                                      ARM_Vector* yields,
                                      int MMVsFut, 
                                      ARM_Currency* ccy)
{
    int SizeM=0, SizeF=1, SizeB=0, i, Nb;
    
    int Tx_Spot=0;
    
    ARM_Date matDate, spotDate, settleDate;

    int CurrentMonth, settleGap;

    ARM_Matrix MMRates(100, 3, 0.0);
    ARM_Matrix Futures(300, 3, 0.0);

    int fwdRule = ccy->GetFwdRule();
    int spotDays = ccy->GetSpotDays();

    char matu;
    int year;
    int month = 0;

    spotDate = itsAsOfDate;

 
    char ccyName[10];

    strcpy(ccyName, (const char*) GetCurrencyUnit()->GetCcyName());

	ccy->CalcFloatPayCal(ccyName);

    spotDate.NextBusinessDay(spotDays, ccyName);
    
    CurrentMonth = spotDate.GetMonth();
    

    matu = 'X';
    
    i = 0;

    while ( Terms[i][0] != 'X' )
    {   
       if ( data->Elt(i) < 1.0E-6 )
       {
          throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Non zero rates are required for building zero curve");
       }

        // si 1er caract alphanum -> Futures
        if (isalpha(Terms[i][0]))
        {
            GetMonthYearFromExpiryDate(Terms[i], &month, &year);

            matu='Z';
        }
        else // sinon MM ou Swap
        {
            sscanf(Terms[i], "%d%c", &Nb, &matu);

            matu = toupper(matu);
        }

        // Conversion des Termes char -> Termes double

        if ( matu == 'D' ) // Ex : "1D"
        {    
            matDate = itsAsOfDate;
            
            matDate.NextBusinessDay(Nb, ccyName);
    
            MMRates.Elt(SizeM, 0) = data->Elt(i);
            
            MMRates.Elt(SizeM, 1) = matDate.GetJulian();
            
            Tx_Spot = 1;

            SizeM++;
        }
        else if ( matu == 'W' )  
        {   //  Ex : "1W"    

            matDate = spotDate;
            
            matDate.AddDays(Nb*7);

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

            MMRates.Elt(SizeM, 0) = data->Elt(i);
            
            MMRates.Elt(SizeM, 1) = matDate.GetJulian();

            SizeM++;
        }
        else if ( matu == 'M' ) 
        {   //  Ex : "9M"    

            matDate = spotDate;
            
            matDate.AddMonths(Nb, GOTO_END_OF_MONTH);

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

            MMRates.Elt(SizeM, 0) = data->Elt(i);
            
            MMRates.Elt(SizeM, 1) = matDate.GetJulian();

            SizeM++;
        }
        // traitement du segment Futures
        else if ( month > 0)
        {  
            matDate.ChangeDate(1, month, year);

            matDate.PiborDelivery(); 

            // tester si NextDate < = > itsAsOfDate
            if ( matDate.GetJulian() <= itsAsOfDate.GetJulian() + 2)
            {
               char BUF[100];

               sprintf(BUF, "Futures contract expiring : %s  %d",
                       Terms[i], year);

               throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, BUF);
            }

            Futures.Elt(SizeF, 0) = 100.0 - data->Elt(i);
            
            Futures.Elt(SizeF, 1) = matDate.GetJulian();

            SizeF++;
        }
        else  
           break;

        i++;
    }

    if (Tx_Spot == 0)  // pas de taux spot
    {
        ARM_Matrix spotRate(1, 3, 0.0);
        
        spotRate.Elt(0, 1) = spotDate.GetJulian();

        MMRates = spotRate & MMRates;

        SizeM += 1;
    /*
        char BUF[100];
        sprintf(BUF,"Spot rate not found ");
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,BUF);
    */
    }
    
    if (SizeF == 1)  // pas de contrats futures ->Forcer les options 
    {
        MMVsFut = K_MM;
    }

    // Get bonds data

    SizeB = bonds->GetSize();
    
    ARM_Matrix BondYields(SizeB, 3, 0.0);

    ARM_Bond* bd =(ARM_Bond *) bonds->Start();
    
    i = 0;
        
    while (bd)
    {
        settleGap = bd->GetSettlementGap();

        settleDate = itsAsOfDate;

        bd->CptCFYearTerms(settleDate);

  //      settleDate.AddDays(settleGap);
        
  //      settleDate.GoodBusinessDay(K_FOLLOWING, ccyName);
        
        BondYields.Elt(i, 0) = bd->YieldToPrice(settleDate, yields->Elt(i));

        BondYields.Elt(i, 0) += bd->Accrued(settleDate);

        BondYields.Elt(i, 1) = bd->GetMaturity().GetJulian();

        BondYields.Elt(i, 2) = double(i);

        bd = (ARM_Bond *) bonds->Next();

        i++;
    }

    BondYields = BondYields.Sort(1);

    // calcul des taux ZC implicites dans chaque segment du marche
    
    CptMMZeroRates(&MMRates, SizeM);

    MMRates = ARM_Matrix(&MMRates, 0, SizeM-1);

    if (SizeF > 1)
       CptFuturesZeroRates(MMRates, SizeM, &Futures, SizeF, 
                           SizeB,
                           BondYields, 
                           MMVsFut, K_LINEAR);

    CptBondZeroRates(MMRates, bonds, &BondYields);

    // Preparation du tableau des taux ZC generes

    int zcSize = 0;

    for (i = 0; i < SizeM; i++)
    {
        if (MMVsFut != K_FUT) 
           zcSize++;
        else if ( MMRates.Elt(i, 1) < Futures.Elt(0, 1) ) 
           zcSize++;
    }

    for (i = 0; i < SizeF; i++)
    {
        if ( (MMVsFut != K_FUT) && 
             (MMRates.Elt(SizeM-1, 1) < Futures.Elt(i, 1)) )
        {
           if ( Futures.Elt(i, 1) < BondYields.Elt(1, 1) )  
               zcSize++;
        }

        if ( MMVsFut == K_FUT )
        {
           if ( Futures.Elt(i, 1) < BondYields.Elt(1, 1) )  
              zcSize++;
        }
    }

    for (i = 0; i < SizeB; i++)
    {
        if ( BondYields.Elt(i, 1) > Futures.Elt(SizeF-1, 1) ) 
           zcSize++;
    }

    // MA : delete before !

    if (itsDateTerms)
       delete itsDateTerms;

    if (itsYearTerms)
       delete itsYearTerms;

    if (itsZeroRates)
       delete itsZeroRates;

    if (itsDiscountFactors)
       delete itsDiscountFactors;


    itsDateTerms = new ARM_Vector (zcSize, 0.0);
    itsYearTerms = new ARM_Vector (zcSize, 0.0);
    itsZeroRates = new ARM_Vector (zcSize, 0.0);
    itsDiscountFactors = new ARM_Vector(zcSize, 0.0);

    zcSize = 0;

    for (i = 0; i < SizeM; i++)
    {
        if ( MMVsFut != K_FUT )
        {
           itsDateTerms->Elt(zcSize) = 
                MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();     
            
           itsDiscountFactors->Elt(zcSize) = MMRates.Elt(i, 0);     
            
           itsYearTerms->Elt(zcSize) = itsDateTerms->Elt(zcSize)/365.0;
            
           zcSize++;
        }
        else if ( MMRates.Elt(i, 1) < Futures.Elt(0, 1) )
        {
           itsDateTerms->Elt(zcSize) = 
                MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();     
            
           itsDiscountFactors->Elt(zcSize) = MMRates.Elt(i, 0);     
            
           itsYearTerms->Elt(zcSize) = itsDateTerms->Elt(zcSize)/365.0;
            
           zcSize++;
        }
    }

    for (i = 0; i < SizeF; i++)
    {
        if  (MMVsFut == K_FUT) 
        {
            if (Futures.Elt(i, 1) < BondYields.Elt(1, 1) )          
            {
               itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
               itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
               itsYearTerms->Elt(zcSize) = itsDateTerms->Elt(zcSize)/365.0;
            
               zcSize++;
            }
        }
        else if ( MMRates.Elt(SizeM-1, 1) < Futures.Elt(i, 1) )
        {
           if ( Futures.Elt(i, 1) < BondYields.Elt(1, 1) ) 
           {
              itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
              itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
              itsYearTerms->Elt(zcSize) = itsDateTerms->Elt(zcSize)/365.0;
            
              zcSize++;
           }
        }
    }

    for (i = 0; i < SizeB; i++)
    {
        if ( BondYields.Elt(i, 1) > Futures.Elt(SizeF-1, 1) )
        {
           itsDateTerms->Elt(zcSize) = 
                BondYields.Elt(i, 1) - itsAsOfDate.GetJulian();

           itsYearTerms->Elt(zcSize) = BondYields.Elt(i, 2);
            
           itsDiscountFactors->Elt(zcSize) = BondYields.Elt(i, 0);
            
           zcSize++;
        }
    }


    for (i = 0; i < zcSize; i++)
    {
        itsZeroRates->Elt(i) = -100.0*log(itsDiscountFactors->Elt(i))/
                                  itsYearTerms->Elt(i);
    }
}



// Calcul des DF a partir des futures

void ARM_ZeroCurve::CptFuturesZeroRates(ARM_Matrix& MMRates, int SizeM,
                                        ARM_Matrix* Future, int SizeF,
                                        int SizeS,
                                        ARM_Matrix& Swap,
                                        int MMVsFut, int Cont_Lin)
{
    ARM_Date D2_Cash,D1_Cash;
    ARM_Date Ctr_Ech1, Ctr_Ech2;
    int i, j, first;
    double N6;
    double Df2,Df1, Df_Ctr, TmP; //inter


    // char* ccyName = GetCurrencyUnit()->GetCcyName();
    
    char ccyName[10];

    strcpy(ccyName,(const char*)GetCurrencyUnit()->GetCcyName());

	GetCurrencyUnit()->CalcFloatPayCal(ccyName);


    // Procedure de tri des taux Futures

    for (i = 1; i < SizeF; i++)
    {
        for (j = i; j < SizeF ; j++)
        {  
            if ( Future->Elt(i, 1) > Future->Elt(j, 1) )
            {
               TmP = Future->Elt(j, 1);
               Future->Elt(j, 1) = Future->Elt(i, 1);
               Future->Elt(i, 1) = TmP;
                
               TmP = Future->Elt(j,0);
               Future->Elt(j, 0) = Future->Elt(i, 0);
               Future->Elt(i, 0) = TmP;
            }
        }
    }

    // Calcul des Df a partir des taux Futures

    for (i = 1; i < SizeF-1; i++)
    {
        if ( fabs(Future->Elt(i+1, 1) - Future->Elt(i, 1) - 90.0) > 10.0 )
        {
           char BUF[100];


           sprintf(BUF, "%s",
               "One or more future contract absent between two contracts ");

           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,BUF);
        }
    }

    first = 1;

    for (i = 1; (( i < SizeF ) && ( MMVsFut == K_MM )); i++)
    {
        if ( Future->Elt(i, 1) < MMRates.Elt(SizeM-1, 1) ) 
           first = i;
    }

    i = SizeM-1;

    while ( MMRates.Elt(i, 1) > Future->Elt(first, 1) && i>0) 
    {
        i--;
    }
   
    int inputDayCount = KACTUAL_360;

    if ( (strcmp(ccyName, "GBP") == 0) || (strcmp(ccyName, "CAD") == 0) )
    {
       inputDayCount = KACTUAL_365;
    }

    if ( i == SizeM-1 )
    {
		double dte_t0tF  = DaysBetweenDates(inputDayCount, 
                                            itsAsOfDate.GetJulian(), Future->Elt(first, 1))/365.0; 
		
		double dte_t0tMM = DaysBetweenDates(inputDayCount, 
                                            itsAsOfDate.GetJulian(), MMRates.Elt(SizeM-1, 1))/365.0; 

		double taux = -log(MMRates.Elt(SizeM-1,0))/(dte_t0tMM);

		Df_Ctr = exp(-taux * dte_t0tF);
    }
    else
    {
	   D1_Cash = (ARM_Date) MMRates.Elt(i, 1);
       Df1 = MMRates.Elt(i, 0);

       D2_Cash = (ARM_Date) MMRates.Elt(i+1, 1);
       Df2 = MMRates.Elt(i+1, 0);

	   Df_Ctr = Interpol(D1_Cash.GetJulian(),
                         D2_Cash.GetJulian(),
                         Future->Elt(first, 1),
                         Df1,
                         Df2,
                         Cont_Lin,
                         inputDayCount);
    }

    Ctr_Ech1 = (ARM_Date) Future->Elt(first, 1);

    Future->Elt(first-1, 0) = Df_Ctr;

    Future->Elt(first-1, 1) = Future->Elt(first, 1);;
    
    for (i = 0; i < first-1; i++) 
    {
        Future->Elt(i, 1) = 0.0;
    }

    i = first;

    while ( i < SizeF )
    {
        Ctr_Ech2 = Ctr_Ech1;

		Ctr_Ech2.PiborDelivery(ccyName);

        N6 = DaysBetweenDates(inputDayCount, Ctr_Ech1, Ctr_Ech2);

        if ( (strcmp(ccyName, "GBP") == 0) || (strcmp(ccyName, "CAD") == 0) )
        {
           Df_Ctr = Df_Ctr / ( 
                            ( 1 + (Future->Elt(i, 0)) / 100.0 
                            * ( N6  / 365.0 ))
                         );
        }
        else
        {
           Df_Ctr = Df_Ctr / (
                            ( 1 + (Future->Elt(i, 0)) / 100.0
                            * ( N6  / 360.0 ))
                         );
        }

        Future->Elt(i, 0) = Df_Ctr;
        Future->Elt(i, 1) = Ctr_Ech2.GetJulian();
        Ctr_Ech1 = Future->Elt(i, 1) ;

        i++;
    }
}


// Calcul des DF a partir des futures

void ARM_ZeroCurve::CptAllForwardsZeroRates(ARM_Matrix& MMRates, int SizeM,
                                        ARM_Matrix* AllFwd, int SizeF,
                                        int SizeS,
                                        ARM_Matrix& Swap,
                                        int MMVsFut, int Cont_Lin)
{
    ARM_Date D2_Cash,D1_Cash;
    ARM_Date Ctr_Ech1, Ctr_Ech2;
    int i, j, first;
    double N6;
    double Df2,Df1, Df_Ctr, TmP; //inter


    // char* ccyName = GetCurrencyUnit()->GetCcyName();
    
    char ccyName[10];

    strcpy(ccyName,(const char*)GetCurrencyUnit()->GetCcyName());

	GetCurrencyUnit()->CalcFloatPayCal(ccyName);


    // Procedure de tri des taux Futures

    for (i = 1; i < SizeF; i++)
    {
        for (j = i; j < SizeF ; j++)
        {  
            if ( AllFwd->Elt(i, 1) > AllFwd->Elt(j, 1) )
            {
               TmP = AllFwd->Elt(j, 1);
               AllFwd->Elt(j, 1) = AllFwd->Elt(i, 1);
               AllFwd->Elt(i, 1) = TmP;
                
               TmP = AllFwd->Elt(j,0);
               AllFwd->Elt(j, 0) = AllFwd->Elt(i, 0);
               AllFwd->Elt(i, 0) = TmP;

			   TmP = AllFwd->Elt(j,2);
               AllFwd->Elt(j, 2) = AllFwd->Elt(i, 2);
               AllFwd->Elt(i, 2) = TmP;
            }
        }
    }

    // Calcul des Df a partir des taux Futures

/*
    for (i = 1; i < SizeF-1; i++)
    {
        if ( fabs(Future->Elt(i+1, 1) - Future->Elt(i, 1) - 90.0) > 10.0 )
        {
           char BUF[100];


           sprintf(BUF, "%s",
               "One or more future contract absent between two contracts ");

           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,BUF);
        }
    }
*/

    // Recherche du 1er contrat à traiter :
	// le 1er si MMVsFUT == K_FUT
	// sinon celui dont la date de fin dépasse
    // la plus grande date d'échéance des taux MM

    first = 1;

    for (i = 1; (( i < SizeF ) && ( MMVsFut == K_MM )); i++)
    {
		
        if ( AllFwd->Elt(i, 1) < MMRates.Elt(SizeM-1, 1) ) 
           first = i;
    }

	// Recherche de la plus grande échéance MM inférieure
	// à la date de départ du 1er contrat à traiter

    i = SizeM-1;

    while ( MMRates.Elt(i, 1) > AllFwd->Elt(first, 2) && i>0) 
    {
        i--;
    }
   
    int inputDayCount = KACTUAL_360;

    if ( (strcmp(ccyName, "GBP") == 0) || (strcmp(ccyName, "CAD") == 0) )
    {
       inputDayCount = KACTUAL_365;
    }

    if ( i == SizeM-1 )
    {
		double dte_t0tF  = DaysBetweenDates(inputDayCount, 
                                            itsAsOfDate.GetJulian(), AllFwd->Elt(first, 2))/365.0; 
		
		double dte_t0tMM = DaysBetweenDates(inputDayCount, 
                                            itsAsOfDate.GetJulian(), MMRates.Elt(SizeM-1, 1))/365.0; 

		double taux = -log(MMRates.Elt(SizeM-1,0))/(dte_t0tMM);

		Df_Ctr = exp(-taux * dte_t0tF);
    }
    else
    {
	   D1_Cash = (ARM_Date) MMRates.Elt(i, 1);
       Df1 = MMRates.Elt(i, 0);

       D2_Cash = (ARM_Date) MMRates.Elt(i+1, 1);
       Df2 = MMRates.Elt(i+1, 0);

	   Df_Ctr = Interpol(D1_Cash.GetJulian(),
                         D2_Cash.GetJulian(),
                         AllFwd->Elt(first, 2),
                         Df1,
                         Df2,
                         Cont_Lin,
                         inputDayCount);
    }

    Ctr_Ech1 = (ARM_Date) AllFwd->Elt(first, 2); // Start Date du 1er forward

    AllFwd->Elt(first-1, 0) = Df_Ctr;

    AllFwd->Elt(first-1, 1) = Ctr_Ech1.GetJulian();

    
    for (i = 0; i < first-1; i++) 
    {
        AllFwd->Elt(i, 1) = 0.0;
    }

    i = first;

    while ( i < SizeF )
    {
        Ctr_Ech1 = AllFwd->Elt(i, 2) ; // Ctr_Ech2 = Ctr_Ech1;

		Ctr_Ech2 = AllFwd->Elt(i, 1) ; // Ctr_Ech2.PiborDelivery(ccyName);

        N6 = DaysBetweenDates(inputDayCount, Ctr_Ech1, Ctr_Ech2);

		// interpolate StartDisc Fut from already fitted MM's and/or future's disc 

		if (( MMVsFut == K_MM ) && (Ctr_Ech1 <= MMRates.Elt(SizeM-1, 1)))
			Df_Ctr = LocDiscInterp(&MMRates, SizeM-1,itsAsOfDate, Ctr_Ech1, Cont_Lin);
		else
			Df_Ctr = LocDiscInterp(AllFwd, i-1,itsAsOfDate, Ctr_Ech1, Cont_Lin);
		    
        if ( (strcmp(ccyName, "GBP") == 0) || (strcmp(ccyName, "CAD") == 0) )
        {
           Df_Ctr = Df_Ctr / ( 
                            ( 1 + (AllFwd->Elt(i, 0)) / 100.0 
                            * ( N6  / 365.0 ))
                         );
        }
        else
        {
           Df_Ctr = Df_Ctr / (
                            ( 1 + (AllFwd->Elt(i, 0)) / 100.0
                            * ( N6  / 360.0 ))
                         );
        }

        AllFwd->Elt(i, 0) = Df_Ctr;

        AllFwd->Elt(i, 1) = Ctr_Ech2.GetJulian();
        
        i++;
    }
}



void ARM_ZeroCurve::CptSwapZeroRates( ARM_Matrix& MMRates, int SizeM, 
									  ARM_Matrix* Swap, 
							          int& SizeS, int fxPayFreq, 
								      ARM_Matrix& Future, int SizeF, int MMVsFut, 
									  int SwapVsFut, int raw,
									  int Cont_Lin,
								      int indexAUD,
									  int fxDayCount)
{
	int firstIndexSwap = 1;
	
    CptSwapZeroRates(MMRates, SizeM, Swap, SizeS, fxPayFreq,
                                     Future, SizeF, 
                                     MMVsFut,
                                     SwapVsFut,
                                     raw,
                                     Cont_Lin,
									 firstIndexSwap,
									 indexAUD, 
									 fxDayCount);
}



void ARM_ZeroCurve::CptSwapZeroRates(ARM_Matrix& MMRates, int SizeM,
                                     ARM_Matrix* Swap, 
                                     int& SizeS, int fxPayFreq,
                                     ARM_Matrix& Future, int SizeF, 
                                     int MMVsFut,
                                     int SwapVsFut,
                                     int raw,
                                     int Cont_Lin,
									 int& firstIndexSwap, // for the method "raw"
									 int indexAUD, // indexAUD is the last index of the AUD with freq=4
									 int fxDayCount) 
{
    int i, j, k, iter = 0;
    double TmP, inter, tmp2, TmP2, startDF;

    ARM_Date startDate = itsAsOfDate;

	char ccyName[10];
	GetCurrencyUnit()->CalcFloatPayCal(ccyName);

    int spotDays   = GetCurrencyUnit()->GetSpotDays();

	if (fxDayCount == KNOBASE)
       fxDayCount = GetCurrencyUnit()->GetFixedDayCount();

    int fxPayFreqVal;

	if(strcmp(ccyName, "AUD")==0 && indexAUD>0)
		fxPayFreqVal = 4;
	else
	    fxPayFreqVal = fxPayFreq ;

    if ( spotDays == 0 ) 
       startDF = 1.0;
    else
    {
       startDate.NextBusinessDay(spotDays, ccyName);

       startDF = Interpol(MMRates.Elt(0, 1),
                          MMRates.Elt(1, 1),
                          startDate.GetJulian(),
                          MMRates.Elt(0, 0),
                          MMRates.Elt(1, 0),
                          Cont_Lin,
                          KACTUAL_360);
    }
    
    ARM_Matrix Df_interpol(200, 2, 0.0);
    double Df_Ctr;

    // Compute DF for maturities < 1 year when fixed leg payment frequency 
    // is > 1 

    if ( fxPayFreq > 1 )
    {
        if ( MMVsFut != K_FUT )
        {
           for (i = 1; i < fxPayFreqVal; i++) 
           {
               j = 0;

               while (( j < SizeM )
                      &&
                      ( Swap->Elt(i, 1) > MMRates.Elt(j, 1) )
                     )
               {
                   j++;
               }

               if ( j <= SizeM-1 ) // OK
               {
                  Swap->Elt(i, 0) = Interpol(MMRates.Elt(j-1, 1),
                                             MMRates.Elt(j, 1),
                                             Swap->Elt(i, 1),
                                             MMRates.Elt(j-1, 0),
                                             MMRates.Elt(j, 0),
                                             Cont_Lin,
                                             KACTUAL_360);
               }
               else
               {
                  int idxj = 0;

				  if(SizeF!=1)
				  {
					  while ( Swap->Elt(i, 1) > Future.Elt(idxj, 1) )
						  idxj++;

					  Swap->Elt(i, 0) = Interpol(Future.Elt(idxj-1, 1),
												 Future.Elt(idxj, 1),
												 Swap->Elt(i, 1),
												 Future.Elt(idxj-1, 0),
												 Future.Elt(idxj, 0),
												 Cont_Lin,
												 KACTUAL_360);
				  }
				  else
				  {
					  Swap->Elt(i, 0) = Interpol(MMRates.Elt(SizeM-2, 1),
												 MMRates.Elt(SizeM-1, 1),
												 Swap->Elt(i, 1),
												 MMRates.Elt(SizeM-2, 0),
												 MMRates.Elt(SizeM-1, 0),
												 Cont_Lin,
												 KACTUAL_360);
//					  ttesti = i;
//					  break;

				  }

               }

               Df_interpol.Elt(i, 0) = Swap->Elt(i, 0);
           }
        }
        else 
        {
           for (i = 1; i < fxPayFreqVal; i++) 
           {
               j = 0;

               while ( Swap->Elt(i, 1) > Future.Elt(j, 1) ) 
                     j++;

               Swap->Elt(i, 0) = Interpol(Future.Elt(j-1, 1),
                                    Future.Elt(j, 1),
                                    Swap->Elt(i, 1),
                                    Future.Elt(j-1, 0),
                                    Future.Elt(j, 0),
                                    Cont_Lin,
                                    KACTUAL_360);

               Df_interpol.Elt(i, 0) = Swap->Elt(i, 0);
            }
        }

    }

	int SizeSVal = SizeS ; 
  
	if ( SwapVsFut != K_FUT )
		SizeSVal = firstIndexSwap-1;


    for (i = fxPayFreqVal; i <= SizeSVal ; i++)//for (i = fxPayFreqVal; i <= SizeS ; i++)
    {
		for (j = 1; j < SizeF; j++)
        {
			// si FUT, on interpole les DF aux maturites 
            // des swaps < maturite du dernier future
			// si SWAP, on interpole les DF aux maturites 
            // des swaps < maturite du dernier future && swaps < first input market swap
            // en utilisant les DF calcules a partir des futures
            if ( Future.Elt(j, 1) > Swap->Elt(i, 1) )
			{
				Df_Ctr = Interpol(Future.Elt(j-1, 1),
                                     Future.Elt(j, 1),
                                     Swap->Elt(i, 1),
                                     Future.Elt(j-1, 0),
                                     Future.Elt(j, 0),
                                     Cont_Lin,
                                     KACTUAL_360);
                    
                Df_interpol.Elt(i, 0) = Df_Ctr;

                // calcul des taux de swap
                TmP = 0.0;
                    
                for (k = 1; k <= i; k++)
                {
					inter = CountYears(fxDayCount,
                                       ( k == 1 ? startDate :
                                           (ARM_Date) Swap->Elt(k-1, 1) ),
                                           (ARM_Date) Swap->Elt(k, 1));
    
                    TmP += Df_interpol.Elt(k, 0)*inter;
                }

                Swap->Elt(i, 0) = 100.0*(startDF-Df_interpol.Elt(i, 0))/TmP;

                break;
            }
        }
    }

    for (i = 1; i <= SizeS; i++)
    { 
        Df_interpol.Elt(i, 0) = 0.0;
    }

    // si methode PAR, on calcule les taux de swap manquant 
    // par interpolation lineaire
    // des taux de swap renseignes 


    if ( raw == K_PAR )
    {
        
		// Ajout : cas ou freq > 1; si le 12M existe en MM ou Fut, on prend le taux swap 1Y egal à ce taux
		
		int idx = -1;
		
		if (fxPayFreq > 1)
		{
		// interpolation of 1Y swap rate when fixed freq > 1

			if (SizeF>1)
			{
				// 1Y swap rate = 12M Fut if exists
				for (k = 0; k <= SizeF; k++)
				{
					if (Future.Elt(k,1) == Swap->Elt(fxPayFreqVal,1))
					{
						idx = k;
						break;
					}
				}
			}
			if (SizeF == 1 && SizeM > 1)
			{
				// 1Y swap rate = 12M MM rate if exists
				for (k = 0; k <= SizeM; k++)
				{
					if (MMRates.Elt(k,1) == Swap->Elt(fxPayFreqVal,1))
					{
						idx = k;
						break;
					}
				}
			}
		}

		// Fin ajout


		for (i = fxPayFreqVal ; i <= SizeS; i++) 
       {
           if ( Swap->Elt(i, 0) <= MINDOUBLE+0.0001 ) 
           {
              j = i+1;

              while ( Swap->Elt(j, 0) < MINDOUBLE+0.0001 )  
                    j++;
            
			  // Do the right interpolation
			  // Modif : le taux swap 1Y est egal au taux MM ou Fut s'ils existent

			  if (( (Swap->Elt(i-1, 0) <= MINDOUBLE+0.0001)|| (i-1 == 0) || ((i == fxPayFreqVal) && (idx != -1))) && ( SwapVsFut == K_SWAP ))
			  {
				
				 // 1Y swap rate

			     ARM_Date swapStartDate = itsAsOfDate;

                 swapStartDate.NextBusinessDay(spotDays, ccyName);

				 swapStartDate.AddYears(1);
 
                 swapStartDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

                 Swap->Elt(i, 1) = swapStartDate.GetJulian();

				 int idxj = 0;
				
				 double dfTemp;

				 if(SizeF>1)
				 {
					 idxj = 0;

					 while ( Swap->Elt(i, 1) > Future.Elt(idxj, 1) )
						   idxj++;

					 dfTemp = Interpol(Future.Elt(idxj-1, 1),
									  Future.Elt(idxj, 1),
									  Swap->Elt(i, 1),
									  Future.Elt(idxj-1, 0),
									  Future.Elt(idxj, 0),
									  Cont_Lin,
									  KACTUAL_360);
				 }
				 if(SizeF==1 && SizeM >1)
				 {
					 idxj = 0;

					 while ( Swap->Elt(i, 1) > MMRates.Elt(idxj, 1) )
						   idxj++;

					 dfTemp = Interpol(MMRates.Elt(idxj-1, 1),
									  MMRates.Elt(idxj, 1),
									  Swap->Elt(i, 1),
									  MMRates.Elt(idxj-1, 0),
									  MMRates.Elt(idxj, 0),
									  Cont_Lin,
									  KACTUAL_360);					
				 }
                // calcul des taux de swap 
                TmP = 0.0;
                
				// Ajout : Prend en compte le cas général où freq peut être > 1
				// On doit tenir compte des taux des freq intermediaires pour calculer le taux swap 1Y
				for (k = 1; k < i; k++)
				{

					inter = CountYears(fxDayCount,
										k==1 ? startDate : (ARM_Date) Swap->Elt(k-1, 1),
									   (ARM_Date) Swap->Elt(k, 1));

					TmP += Swap->Elt(k,0)*inter;  // les taux swap des freq precedentes sont en fait des ZC
					//TmP += (startDF - Swap->Elt(k, 0) * inter / 100 * TmP) / (1 + Swap->Elt(k, 0) * inter / 100);
				}
				// Fin ajout

				inter = CountYears(fxDayCount,
									i==1 ? startDate : (ARM_Date) Swap->Elt(i-1, 1),
									(ARM_Date) Swap->Elt(i, 1));
				TmP += dfTemp*inter;
				Swap->Elt(i, 0) = 100.0 * (startDF - dfTemp) / TmP;
				
			  }
			  else
			  {
				
				 double x = DaysBetweenDates(fxDayCount, 
                                        (ARM_Date) Swap->Elt(i-1, 1),
                                           (ARM_Date) Swap->Elt(i, 1));

			     double y = DaysBetweenDates(fxDayCount, (ARM_Date) Swap->Elt(i-1, 1),
                                              (ARM_Date) Swap->Elt(j, 1));

                 Swap->Elt(i, 0) = Swap->Elt(i-1, 0)+(Swap->Elt(j, 0)-Swap->Elt(i-1, 0))
                                   *(x/y);
			  }
           }
       }
    }


	ARM_Vector FlowFreq(200, MINDOUBLE);  // the flows whose payment frequency= 'fxPayFreq' (eg. 1Y,2Y...pour EUR; 6M,1Y,18M...pour USD)

	ARM_Vector FlowNoFreq(200, MINDOUBLE);// the flows whose payment frequency is more than 'fxPayFreq'

	int SizeFlowFreq=0, SizeNoFlowFreq=0, SizeAUDtri=0; // SizeAUDtri: (variable for the currency AUD),Nb of the dates with the freq 3 monthes
	
	PaymentFlowVector( Swap, SizeS, fxPayFreq, startDate, ccyName, indexAUD,
					   FlowFreq,   SizeFlowFreq,
					   FlowNoFreq, SizeNoFlowFreq,
					   SizeAUDtri);

	// calculer DF pour des points correspondant à freqence (variable : FlowFreq)
	double SwapEltAv,SwapEltDateAv;
	
	i=1;
	
	while(FlowFreq.Elt(i)<fxPayFreqVal && i<=SizeFlowFreq)
		i++;
    for (i; i <= SizeFlowFreq; i++)//for (i = fxPayFreq; i <= SizeS; i++) 
    {
        // calcul implicite de chaque taux non renseigne

        if (( Swap->Elt(FlowFreq.Elt(i), 0) < MINDOUBLE+0.0001 ) && ( raw == K_RAW )) 
        {
           // chercher les taux encadrant le taux a calculer 

           j = i+1;

           while ( Swap->Elt(FlowFreq.Elt(j), 0) < MINDOUBLE+0.0001 ) 
                 j++;
            
           //Df_Ctr = Swap->Elt(FlowFreq.Elt(i-1), 0);
		   Df_Ctr = i==1?Swap->Elt(FlowFreq.Elt(i-1), 0) : Future.Elt(SizeF-1, 0);
            
           do
           {
                for (k = i; k < j; k++)
                {
					if(Swap->Elt(FlowFreq.Elt(i-1),1)<=Future.Elt(SizeF-1, 1) && SwapVsFut != K_SWAP)
					{
						SwapEltDateAv = Future.Elt(SizeF-1, 1);
						SwapEltAv  = Future.Elt(SizeF-1, 0);
					}
					else
					{
						SwapEltDateAv = Swap->Elt(FlowFreq.Elt(i-1),1);
						SwapEltAv  = Swap->Elt(FlowFreq.Elt(i-1),0);
					}

					Swap->Elt(FlowFreq.Elt(k), 0) = Interpol(SwapEltDateAv,
															 Swap->Elt(FlowFreq.Elt(j), 1),
															 Swap->Elt(FlowFreq.Elt(k), 1),
															 SwapEltAv,
															 Df_Ctr,
															 Cont_Lin,
															 KACTUAL_360);

					// Derivee
                    Df_interpol.Elt(FlowFreq.Elt(k), 0) = (	-Swap->Elt(FlowFreq.Elt(k), 0)+
															Interpol(SwapEltDateAv,
																	 Swap->Elt(FlowFreq.Elt(j), 1),
																	 Swap->Elt(FlowFreq.Elt(k), 1),
																	 SwapEltAv,
																	 Df_Ctr + 1.0e-10,
																	 Cont_Lin,
																	 KACTUAL_360)
														  )/1.0e-10;
                    
                    Df_interpol.Elt(FlowFreq.Elt(k), 1) = 1.0;
                }


                TmP  = 0.0;
                TmP2 = 0.0;

                for (k = 1; k < j; k++)
                {
                    inter = CountYears(fxDayCount,
                           ( k==1 ? startDate :
                                (ARM_Date) Swap->Elt(FlowFreq.Elt(k-1), 1) ),
                                          (ARM_Date) Swap->Elt(FlowFreq.Elt(k), 1) );
                    
                    TmP  +=  Swap->Elt(FlowFreq.Elt(k), 0) * inter;
                    TmP2 +=  Df_interpol.Elt(FlowFreq.Elt(k), 0) * inter;
                }

                TmP  *= Swap->Elt(FlowFreq.Elt(j), 0)/100.0;
                TmP2 *= Swap->Elt(FlowFreq.Elt(j), 0)/100.0;

                TmP  = startDF-TmP; // Fonction d'iteration
                TmP2 = -TmP2;       // Derivee de la fonction d'iteration

                inter = CountYears(fxDayCount,
                          ( j==1 ? startDate :
                                   (ARM_Date) Swap->Elt(FlowFreq.Elt(j-1), 1) ),
                                              (ARM_Date) Swap->Elt(FlowFreq.Elt(j), 1));

                tmp2 = 1.0+Swap->Elt(FlowFreq.Elt(j), 0)*inter/100.0;


                Df_Ctr = Df_Ctr-(TmP-tmp2*Df_Ctr)/(TmP2-tmp2);
            }
            while ( fabs((TmP-tmp2*Df_Ctr)/(TmP2-tmp2)) > 1.0e-15 );
		   //while ( fabs(TmP-tmp2*Df_Ctr) > 1.0e-20 );
        
            Swap->Elt(FlowFreq.Elt(j),0) = Df_Ctr;

            i = j;
        }
        else 
        {

            // si meth = PAR, calculer les DF a partir 
            // des taux de swap interpoles

            //|->

            TmP = 0.0;

            for (j = 1; j < i; j++)
            {
                inter = CountYears(fxDayCount, 
                       ( j==1 ? startDate :
                            (ARM_Date) Swap->Elt(FlowFreq.Elt(j-1), 1) ),
                                      (ARM_Date) Swap->Elt(FlowFreq.Elt(j), 1) );

                TmP += Swap->Elt(FlowFreq.Elt(j), 0)*inter;  
            }

            TmP *= Swap->Elt(FlowFreq.Elt(i), 0)/100.0;

            TmP = startDF-TmP;

            inter = CountYears(fxDayCount, 
                      ( i==1 ? startDate :
                               (ARM_Date) Swap->Elt(FlowFreq.Elt(i-1), 1) ),
                                          (ARM_Date) Swap->Elt(FlowFreq.Elt(i), 1));

            Swap->Elt(FlowFreq.Elt(i), 0) = 1.0 + Swap->Elt(FlowFreq.Elt(i), 0)*inter/100.0;

            Swap->Elt(FlowFreq.Elt(i), 0) = TmP/Swap->Elt(FlowFreq.Elt(i), 0);
        }
	}

	// Calculer DF pour des points qui ne correspondant à freqence (variable : FlowNoFreq)
	i=1;
	while(FlowNoFreq.Elt(i)<=fxPayFreqVal && i<=SizeAUDtri )
		i++;

	for(i; i<=SizeAUDtri; i++) //Cas de AUD: il y a une partie trimestielle(freq=4) qui est différent que freqPay(2)
	{
		TmP = 0.0 ;
			
		for (j=1; j<i; j++)
		{
			inter = CountYears( fxDayCount, ( j==1 ? startDate :(ARM_Date) Swap->Elt(FlowNoFreq.Elt(j-1), 1) ),
								(ARM_Date) Swap->Elt(FlowNoFreq.Elt(j), 1) );
            
			TmP += Swap->Elt(FlowNoFreq.Elt(j), 0)*inter;  
		}
        
		TmP *= Swap->Elt(FlowNoFreq.Elt(i), 0)/100.0;

		TmP = startDF-TmP;

		inter = CountYears(fxDayCount, (ARM_Date) Swap->Elt(FlowNoFreq.Elt(i-1), 1),
							   (ARM_Date) Swap->Elt(FlowNoFreq.Elt(i), 1));

		Swap->Elt(FlowNoFreq.Elt(i), 0) = 1.0 + Swap->Elt(FlowNoFreq.Elt(i), 0)*inter/100.0;

		Swap->Elt(FlowNoFreq.Elt(i), 0) = TmP/Swap->Elt(FlowNoFreq.Elt(i), 0);

        Df_interpol.Elt(FlowNoFreq.Elt(i), 1) = 1.0;

	}

	for(i=SizeAUDtri+1; i<=SizeNoFlowFreq; i++) //for(i; i<=SizeNoFlowFreq; i++)
	{
		Swap->Elt(FlowNoFreq.Elt(i), 0) = CalculSwapFlowNoFreq(Swap->Elt(FlowNoFreq.Elt(i),1), 
																 Swap->Elt(FlowNoFreq.Elt(i),0), 
																 startDF,
																 MMRates,SizeM, 
																 Future, SizeF, 
																 Swap, SizeS,
																 fxPayFreq,
																 Cont_Lin,
																 fxDayCount,
																 spotDays,
																 ccyName);
	}

    // si meth = RAW, virer les taux manquant initialement 

    if ( raw == K_RAW )
    {
       iter = 0;

       for (i = 1; i <= SizeS; i++)
       {
           if ( Df_interpol.Elt(i, 1) > 0.0001 )
           {
              for (j = i; j < SizeS; j++)
              {
                  Df_interpol.Elt(j, 1) = Df_interpol.Elt(j+1, 1);

                  Swap->Elt(j, 0) = Swap->Elt(j+1, 0);
                  Swap->Elt(j, 1) = Swap->Elt(j+1, 1);
              }
			  
			  if(i<=firstIndexSwap)
				  firstIndexSwap--;

              SizeS--;

              i--;

           }
       }
    }
}



void ARM_ZeroCurve::CptBondZeroRates(ARM_Matrix& MMRates, 
                                     ARM_Container* bonds, 
                                     ARM_Matrix* bondPrices)
{
    int i, j, l, lastDFIndex, ij, kl, nbFlows, nbBonds;

    int fstCF, dayCount;
    
    double lastDF, lastDFDate;

    double bdPrice = 0.0, price = 0.0, price_h = 0.0, d_price = 0.0;

    ARM_Vector* cfValues = NULL;

    ARM_Vector* cfTerms = NULL;

    ARM_Vector* cfDates = NULL;

    ARM_Bond* bond;
    
    double DFkl, DFij, DFij_h, term1, term2;

    // Compute DF for maturities < 1 year when fixed leg payment frequency 
    // is > 1 

    nbBonds = bonds->GetSize();

    ARM_Matrix DF = MMRates;

    lastDFIndex = DF.GetNumLines() - 1;

    for (i = 0; i <= lastDFIndex; i++)
    {
        DF.Elt(i, 2) = (DF.Elt(i, 1) - itsAsOfDate.GetJulian())/365.0;
    }

    lastDF = DF.Elt(lastDFIndex, 0);

    lastDFDate = DF.Elt(lastDFIndex, 1);

    for (i = 0; i < nbBonds; i++)
    {
        l = int(bondPrices->Elt(i, 2));

        bond = (ARM_Bond *) bonds->GetNthItem(l);

        dayCount = bond->GetDayCount();

        bdPrice = bondPrices->Elt(i, 0);

        cfValues = bond->GetCashFlowValues();

        cfDates = bond->GetPaymentDates();

        cfTerms = bond->GetYearTerms();

        nbFlows = cfTerms->GetSize();

        fstCF = 0;

        while ( cfTerms->Elt(fstCF) < 0.0 ) 
              fstCF++;
            
        // Discount CF for which DF can be interpolated from known DF 

        for (kl = fstCF; kl < nbFlows; kl++) 
        {
            if ( cfDates->Elt(kl) <= lastDFDate )
            {
               j = 0;

               while (DF.Elt(j, 1) < cfDates->Elt(kl)) 
                     j++;

               term1 = DF.Elt(j-1, 2);
            
               term2 = DF.Elt(j, 2);
            
               DFkl = bondInterpol(term1,
                                term2,
                                cfTerms->Elt(kl),
                                DF.Elt(j-1, 0),
                                DF.Elt(j, 0),
                                K_LINEAR);
                
               bdPrice -= cfValues->Elt(kl)*DFkl;
            }
            else
            {
               lastDF = 0.69079;

               do
               {
                   price = price_h = bdPrice;

                   term1 = DF.Elt(lastDFIndex, 2);
                    
                   for (ij = kl; ij < nbFlows; ij++)
                   {
                       DFij = bondInterpol(term1,
                                    cfTerms->Elt(nbFlows-1),
                                    cfTerms->Elt(ij),
                                    DF.Elt(lastDFIndex, 0),
                                    lastDF,
                                    K_LINEAR);
                    
                       DFij_h = bondInterpol(term1,
                                    cfTerms->Elt(nbFlows-1),
                                    cfTerms->Elt(ij),
                                    DF.Elt(lastDFIndex, 0),
                                    lastDF - 1.0e-10,
                                    K_LINEAR);

                       price -= cfValues->Elt(ij)*DFij;

                       price_h -= cfValues->Elt(ij)*DFij_h;
                    }

                    d_price = (price - price_h)*1.0e10;

                    lastDF = lastDF - price/d_price;
                }
                while ( fabs(price/d_price) > 1.0e-10 );

                kl = nbFlows;

                ARM_Matrix last(1, 3, lastDF);

                lastDFDate = last.Elt(0, 1) = cfDates->Elt(nbFlows-1);

                last.Elt(0, 2) = cfTerms->Elt(nbFlows-1);

                DF = DF & last;

                lastDFIndex++;

                bondPrices->Elt(i, 0) = lastDF;

                bondPrices->Elt(i, 2) = cfTerms->Elt(nbFlows-1);
            }
        }
    }
}



double ARM_ZeroCurve::Interpol(double date1, double date2, 
                               double date_interpol,
                               double Df1, double Df2, 
                               int meth_inter, int base)
{
/*  
  base:

    KACTUAL_ACTUAL          1
    KACTUAL_365             2
    KACTUAL_360             3**
    K30_360                 4
    KACTUAL_REAL            5
*/
    double N1, N2, N3, N4, N5;
    double Df_Ctr;

    ARM_Date D1_Cash  = ( ARM_Date ) date1;
    ARM_Date D2_Cash  = ( ARM_Date ) date2;

    ARM_Date Ctr_Ech1 = ( ARM_Date ) date_interpol;

        
    N1= DaysBetweenDates (base, D1_Cash, Ctr_Ech1);
    N2= DaysBetweenDates (base, D1_Cash, D2_Cash);
    N3= DaysBetweenDates (base, itsAsOfDate, Ctr_Ech1);
    N4= DaysBetweenDates (base, itsAsOfDate, D1_Cash);
    N5= DaysBetweenDates (base, itsAsOfDate, D2_Cash);

    if (meth_inter == K_CONTINUOUS)
    {
       Df_Ctr = pow(Df1, (1.0-N1/N2))*pow(Df2,(N1/N2));
    }
    else // K_LINEAR
    {
       if ( N4 < 0.00001 ) 
          if (N2 < 0.00001)
			  Df_Ctr = 1.;
		  else
			  Df_Ctr = pow(Df2, N3*N3/(N2*N2));
       else 
	   {
		   if (N2 < 0.00001 || N3 > N5)
			   Df_Ctr = pow(Df2, (N3/N5)); // on cristallise le zéro à droite
		   else
			Df_Ctr = pow(Df1, (1.0-N1/N2)*(N3/N4))*
                    pow(Df2, (N1/N2)*(N3/N5));
	   }
    }

    return(Df_Ctr);
}


/*
double ARM_ZeroCurve::InterpolDF(ARM_Vector* date, ARM_Vector* DF, 
                                 double date_interpol, 
                                 double Df1, double Df2,
                                 int meth_inter, int base)
{
    double N1, N2, N3, N4, N5;
    double Df_Ctr;

    ARM_Date D1_Cash = ( ARM_Date ) date1;
    ARM_Date D2_Cash = ( ARM_Date ) date2;

    ARM_Date Ctr_Ech1=( ARM_Date ) date_interpol;

        
    N1 = DaysBetweenDates(base, D1_Cash, Ctr_Ech1);
    N2 = DaysBetweenDates(base, D1_Cash, D2_Cash);
    N3 = DaysBetweenDates(base, itsAsOfDate, Ctr_Ech1);
    N4 = DaysBetweenDates(base, itsAsOfDate, D1_Cash);
    N5 = DaysBetweenDates(base, itsAsOfDate, D2_Cash);

    if (meth_inter)
    {
       Df_Ctr = pow(Df1, (1.0-N1/N2))*pow(Df2, (N1/N2));
    }
    else
    {
       Df_Ctr = pow(Df1, (1.0-N1/N2)*(N3/N4))*pow(Df2, (N1/N2)*(N3/N5));
    }

    return(Df_Ctr);

}

*/



double ARM_ZeroCurve::bondInterpol(double date1, double date2, 
                                   double date_interpol,
                                   double Df1, double Df2, int meth_inter)
{
/*  
  base:

    KACTUAL_ACTUAL          1
    KACTUAL_365             2
    KACTUAL_360             3**
    K30_360                 4
    KACTUAL_REAL            5
*/
    double N1, N2, N3, N4, N5;
    double Df_Ctr;
        
    N1 = date_interpol-date1;
    N2 = date2-date1;
    N3 = date_interpol;
    N4 = date1;
    N5 = date2;

    if ( meth_inter == K_CONTINUOUS )
    {
       Df_Ctr = pow(Df1, (1.0-N1/N2))*pow(Df2, (N1/N2));
    }
    else // K_LINEAR
    {
       Df_Ctr = pow(Df1, (1.0-N1/N2)*(N3/N4))*pow(Df2, (N1/N2)*(N3/N5));
    }

    return(Df_Ctr);

}



void ARM_ZeroCurve::CptMMZeroRates(ARM_Matrix* MMRates, int SizeM)
{
    int i,j;

    double TmP = 1.0, interJJ, term;

    ARM_Date spotDate;

    ARM_Currency* ccy = GetCurrencyUnit();

    char ccyName[10];

    strcpy(ccyName,(const char*)GetCurrencyUnit()->GetCcyName());

	ccy->CalcFloatPayCal(ccyName);

/*
    if (strcmp((const char*)ccyName,"JPY") == 0)
	{
		if (GetCALYPSOVersion == 0)
			strcpy(ccyName,"ZGJ");
		else
			strcpy(ccyName,"GBPJPY");
	}

*/

    int spotDays = ccy->GetSpotDays();
    int dayCount = ccy->GetMMDayCount();
    //dayCount = KACTUAL_360;

    //  Procedure de tri des taux monetaires

    for (i = 0; i < SizeM ; i++)
    {
        for (j = i; j < SizeM; j++)
        {  
            if ( MMRates->Elt(i, 1) > MMRates->Elt(j, 1) )
            {
               TmP = MMRates->Elt(j, 1);
               MMRates->Elt(j, 1) = MMRates->Elt(i, 1);
               MMRates->Elt(i, 1) = TmP;
                
               TmP = MMRates->Elt(j, 0);
               MMRates->Elt(j, 0) = MMRates->Elt(i, 0);
               MMRates->Elt(i, 0) = TmP;
            }
        }
    }

	TmP=1.0;

    interJJ = MMRates->Elt(0, 1) - itsAsOfDate.GetJulian();

    spotDate = itsAsOfDate;

//	if(strcmp(ccyName, "CHF")==0)
//		spotDate.NextBusinessDay(1, ccyName);
//	else
		spotDate.NextBusinessDay(spotDays, ccyName);

    // Calcul des Df a partir des taux monetaires

/**** OLD CODE 
    if ( interJJ <= 7.0 )
    {
       spotDate = (ARM_Date) MMRates->Elt(0, 1);

       term = CountYears(dayCount, itsAsOfDate, ARM_Date(MMRates->Elt(0, 1)));

       TmP  = 1.0 + MMRates->Elt(0, 0)*term/100.0;
    }
*****/

    // MA : 23/05/2002 : In order to correct GBP case
    // we solve the pb in a general manner

    if ( interJJ <= 7.0 )
    {
       term = CountYears(dayCount, itsAsOfDate, spotDate);

       TmP  = 1.0 + MMRates->Elt(0, 0)*term/100.0;
    }
    
    for (i = 0; i < SizeM; i++)
    {    
	   if(spotDate.GetJulian()<MMRates->Elt(i, 1))
	   {
		   term = CountYears(dayCount, spotDate, ARM_Date(MMRates->Elt(i, 1)));
		   MMRates->Elt(i, 0) = 1.0/(TmP*(1+MMRates->Elt(i, 0)*term/100.0));
	   }
	   else
	   {
		   term = CountYears(dayCount,ARM_Date(MMRates->Elt(i, 1)), spotDate);
		   MMRates->Elt(i, 0) = (1+MMRates->Elt(i, 0)*term/100.0)/TmP;
	   }
    }
}


void ARM_ZeroCurve::CptMMZeroRatesFWD(ARM_Matrix* MMRates, int SizeM, int Cont_Lin)
{
    int i,j;

    double TmP = 1.0, term;

    ARM_Date spotDate, startDate, matDate, prevDate;

    ARM_Currency* ccy = GetCurrencyUnit();

    char ccyName[10];

    strcpy(ccyName,(const char*)GetCurrencyUnit()->GetCcyName());

	ccy->CalcFloatPayCal(ccyName);


    int spotDays = ccy->GetSpotDays();
    int dayCount = ccy->GetMMDayCount();
	int Tx_spot = 0;

	spotDate = itsAsOfDate;

	spotDate.NextBusinessDay(spotDays, ccyName);
	
	double discStart = 1;

    //  Procedure de tri des taux monetaires

    for (i = 0; i < SizeM ; i++)
    {
        for (j = i; j < SizeM; j++)
        {  
            if ( MMRates->Elt(i, 1) > MMRates->Elt(j, 1) )
            {
               TmP = MMRates->Elt(j, 1);
               MMRates->Elt(j, 1) = MMRates->Elt(i, 1);
               MMRates->Elt(i, 1) = TmP;
                
               TmP = MMRates->Elt(j, 0);
               MMRates->Elt(j, 0) = MMRates->Elt(i, 0);
               MMRates->Elt(i, 0) = TmP;

			   TmP = MMRates->Elt(j, 2);
               MMRates->Elt(j, 2) = MMRates->Elt(i, 2);
               MMRates->Elt(i, 2) = TmP;
            }
        }
    }

	TmP=1.0;

    // Calcul des Df a partir des taux monetaires

	// interpolate StartDisc from already fitted MM's disc 

   
    for (i = 0; i < SizeM; i++)
    {
       if (MMRates->Elt(i, 2) == itsAsOfDate.GetJulian())
		   discStart = 1.;
	   else
	   {
		   j = 0;
		   while ((ARM_Date(MMRates->Elt(j, 1)) <= ARM_Date(MMRates->Elt(i, 2))) && j < SizeM)
			   j++;

		   if (j > 0)
			   j--;

           TmP = LocDiscInterp(MMRates, j, itsAsOfDate, ARM_Date(MMRates->Elt(i, 2)), Cont_Lin, 1);
		   discStart = TmP; // exp(-TmP/100*(ARM_Date(MMRates->Elt(i, 2)) - itsAsOfDate)/365.);
	   }

	   term = CountYears(dayCount, ARM_Date(MMRates->Elt(i, 2)), ARM_Date(MMRates->Elt(i, 1)));
       TmP = 1+MMRates->Elt(i, 0)*term/100.0;
	   MMRates->Elt(i, 0) = discStart / TmP;
    }
}


void ARM_ZeroCurve::ParallelShift(double value)
{
    int size = itsZeroRates->GetSize();

    for (int i = 0; i < size; i++)
    {
         itsZeroRates->Elt(i) += value;
    }

    GenerateFields();
}



// Calcul de la courbe a partir de depots, futures et swaps

void ARM_ZeroCurve::ZCFromMarketRatesGen(ARM_CRV_TERMS& Terms, 
                                         ARM_Vector* data, 
                                         int MMVsFut, int SwapVsFut,
                                         int raw,
                                         int Cont_Lin,
                                         ARM_Currency* ccy)
{
    int SizeM=0, SizeF=1, SizeS=0, i, k, l, Nb;
    int SizeT = 0;
    int Tx_Spot=0;
    int isContrat = 0;
    int isToY = 0;
    ARM_Date matDate;
    ARM_Date dateContrat, dateUnder;


    int CurrentMonth;

    // Les matrices suivantes ont 2 colonnes, dans la premiere on met les taux
    // dans la seconde les dates
    ARM_Matrix MMRates(100, 2, 0.0);
    ARM_Matrix SwapRates(200, 2, 0.0);

    ARM_Matrix Futures3(100, 3, 0.0);
    ARM_Matrix Futures(100, 2, 0.0);

    ARM_Matrix toyTab(100, 3, 0.0);

    SetCurrencyUnit(ccy);

    char* ccyName = ccy->GetCcyName();
	ccy->CalcFloatPayCal(ccyName);

    int fwdRule   = ccy->GetFwdRule();
    int spotDays  = ccy->GetSpotDays();
    int fxPayFreq = ccy->GetFixedPayFreq();

    char matu;

    int month = 0;

    int lastMMTerm;

    ARM_Date spotDate = itsAsOfDate;

    ARM_Date swapStartDate = itsAsOfDate;

    swapStartDate.NextBusinessDay(spotDays, ccyName);
    
    CurrentMonth = spotDate.GetMonth();
    
    matu = 'X';

    i = 0;

    while ( Terms[i][0] != 'X' || Terms[i][1] != 'X' )
    {   
        isContrat = 0;
        isToY = 0;

        if (data->Elt(i) < 1.0E-6)
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Non zero rates are required for building zero curve");
        }

        // si 1er caract alphanum -> Futures ou "Turn of year"
        if (isalpha(Terms[i][0]))
        {
           if (Terms[i][0] == 'T') // cas Turn of the year
           {
              FromTOYCodeToDates(Terms[i], dateContrat, dateUnder);
              isToY = 1;
           }
           else
           {
              FromFutureCodeToDates(Terms[i], dateContrat, dateUnder);
              isContrat = 1;
           }

           matu='Z';
        }
        else // sinon MM ou Swap
        {
           sscanf(Terms[i], "%d%c", &Nb, &matu);

           matu = toupper(matu);
        }


        if ( matu == 'D' ) // Ex : "1D"
        {    
           matDate = itsAsOfDate;
            
           spotDate.NextBusinessDay(spotDays, ccyName);
    
    
           MMRates.Elt(SizeM, 0) = data->Elt(i);
            
            
           MMRates.Elt(SizeM, 1) = spotDate.GetJulian();

           if ( spotDays < 1 )
           {
              matDate.NextBusinessDay(Nb, ccyName);
              MMRates.Elt(SizeM, 1) = matDate.GetJulian();
           }
            
           Tx_Spot = 1;

           SizeM++;
        }
        else if ( matu == 'W' )  
        {   //  Ex : "1W"    

            matDate = spotDate;
            
            matDate.AddDays(Nb*7);

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

            MMRates.Elt(SizeM, 0) = data->Elt(i);
            
            MMRates.Elt(SizeM, 1) = matDate.GetJulian();

            SizeM++;
        }
        else if ( matu == 'M' ) 
        {   //  Ex : "9M"    

            matDate = spotDate;
            
            matDate.AddMonths(Nb, GOTO_END_OF_MONTH);

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

            MMRates.Elt(SizeM, 0) = data->Elt(i);
            
            MMRates.Elt(SizeM, 1) = matDate.GetJulian();

            SizeM++;

            lastMMTerm = Nb;
        }
        else if ( matu == 'Y')  // ->implicitement ce sont des taux de swap
        {   
           //  Ex : "15Y"

           matDate = swapStartDate;
            
           matDate.AddYears(Nb);

           matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

           if ( lastMMTerm < 12 )
           {
              SwapRates.Elt(Nb*fxPayFreq, 0) = data->Elt(i);
            
              SwapRates.Elt(Nb*fxPayFreq, 1) = matDate.GetJulian();
           }
           else if (lastMMTerm == 12 && Nb > 1)
           {
              SwapRates.Elt((Nb-1)*fxPayFreq, 0) = data->Elt(i);
            
              SwapRates.Elt((Nb-1)*fxPayFreq, 1) = matDate.GetJulian();
           }

           SizeS++;
        }
        // traitement du segment Futures
        else if ( isContrat )
        {  
            Futures3.Elt(SizeF, 0) = 100.0 - data->Elt(i);
            
            // on met la date du contrat dans la deuxieme colonne
            Futures3.Elt(SizeF, 1) = dateContrat.GetJulian();
            Futures3.Elt(SizeF, 2) = dateUnder.GetJulian();

            SizeF++;
        }
        else if ( isToY)
        {
            toyTab.Elt(SizeT, 0) = data->Elt(i);
 
 
            toyTab.Elt(SizeT, 1) = dateContrat.GetJulian();
            toyTab.Elt(SizeT, 2) = dateUnder.GetJulian();
 
            SizeT++;
        }
        else 
           break;

        i++;
    }

    // Control pas deux fois meme date de contrat

    for (k = 0; k < SizeF-1; k++)
    {
        for (l = k+1; l < SizeF; l++)
        {
            if (Futures3.Elt(k,2) == Futures3.Elt(l,2))
            {
               throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Contracts mustn't have same expiry dates ");

            }
        }
    }

    // on insert les toys dans le tableau des contrats future
 
    for (k=0; k < SizeT ; k++)
    {
        // on recherche le contrat qui contient la fin d'annee
        l=0; int indCtr;
        int trouve = 0;
 
        while (l <SizeF && !trouve)
        {
            if (Futures3.Elt(l,1) < toyTab.Elt(k,2) &&
                 toyTab.Elt(k,1) < Futures3.Elt(l,2) )
            {
                trouve = 1;
                indCtr = l;
            }
 
            l++;
        }
 
        if (!trouve)
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "turn of the year not included in future ! ");
 
        }
        // on ajoute 2 lignes au tableau des futures
        SizeF = SizeF+2;
 
        
        Futures3.Elt(SizeF-2,0) = toyTab.Elt(k,0);
        Futures3.Elt(SizeF-2,1) = toyTab.Elt(k,1);
        Futures3.Elt(SizeF-2,2) = toyTab.Elt(k,2);
 
        Futures3.Elt(SizeF-1,1) = toyTab.Elt(k,2);
        Futures3.Elt(SizeF-1,2) = Futures3.Elt(indCtr,2);
 
        Futures3.Elt(indCtr,2) = toyTab.Elt(k,1);

        double d1 = (Futures3.Elt(indCtr,2) - Futures3.Elt(indCtr,1))/360.0;
        double d2 = (Futures3.Elt(SizeF-2,2) - Futures3.Elt(SizeF-2,1))/360.0;
        double d3 = (Futures3.Elt(SizeF-1,2) - Futures3.Elt(SizeF-1,1))/360.0;

        double dtot = (Futures3.Elt(SizeF-1,2)-Futures3.Elt(indCtr,1))/360.0;

        double taux = 1 + Futures3.Elt(indCtr,0)*dtot/100.0;

        taux /= 1+Futures3.Elt(indCtr,0)*d1/100.0;

        taux /= 1+toyTab.Elt(k,0)*d2/100.0;

        taux -= 1;

        taux *= 100/d3;

        Futures3.Elt(SizeF-1,0) = taux;
    }


    if (Tx_Spot == 0)  // pas de taux spot
    {
        ARM_Matrix spotRate(1, 2, 0.0);
        
        spotRate.Elt(0, 1) = itsAsOfDate.GetJulian();

        MMRates = spotRate & MMRates;

        SizeM += 1;
    }
    
    if ( SizeF == 1 )  // pas de contrats futures ->Forcer les options 
    {
       SwapVsFut = K_SWAP;

       MMVsFut = K_MM;
    }

    // il faut le taux swap 1Y sauf si ce sont les futures qui dominent
    if (( SwapRates.Elt(fxPayFreq, 0) < 0.0001 ) && ( SwapVsFut == 1 ))
    {
       char BUF[100];

       sprintf(BUF, "1 year Swap rate not found");
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, BUF);
    }

    for (i = 0; i < SwapRates.GetNumLines(); i++) 
    {
       if ( SwapRates.Elt(i, 1) > 0.0 ) 
          SizeS = i; // sizeS = maturite max de la courbe
    }

    //XX On remplit les lignes des freq intermediaires (6m, 1.5Y 2.5Y ...) 
    //   ds le tableau SwapRates
    for (i = 1; i <= SizeS; i++)
    {
        if ( SwapRates.Elt(i, 1) < 0.0001 ) 
        {
           matDate = swapStartDate;

           matDate.AddMonths(i*12/fxPayFreq, GOTO_END_OF_MONTH);

           matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

           SwapRates.Elt(i, 1) = matDate.GetJulian();
        }
    }

    // calcul des taux ZC implicites dans chaque segment du marche
    
    CptMMZeroRates(&MMRates, SizeM);


    if (SizeF != 1)
       CptFuturesZeroRatesGen(MMRates, SizeM, Futures3, SizeF, SwapRates, 
            MMVsFut, Cont_Lin);

    // on recopie le tableau de 3 col ds celui 2 col
    for (k = 0; k < SizeF; k++)
    {
        Futures.Elt(k, 0) = Futures3.Elt(k, 0);
        Futures.Elt(k, 1) = Futures3.Elt(k, 1);
    }


    if (SizeS > 0)
        CptSwapZeroRates(MMRates, SizeM, &SwapRates, SizeS, fxPayFreq, 
                     Futures, SizeF, MMVsFut, SwapVsFut, raw, Cont_Lin);

    // Preparation du tableau des taux ZC correspondant a generateCurve

    int zcSize = 0;

    for (i = 0; i < SizeM; i++)
    {
        if (MMVsFut != K_FUT) 
           zcSize++;
        else if ( MMRates.Elt(i, 1) < Futures.Elt(0, 1) ) 
           zcSize++;
    }

    for (i = 0; i < SizeF; i++)
    {
        if ( (MMVsFut != K_FUT) && 
             (MMRates.Elt(SizeM-1, 1) < Futures.Elt(i, 1)) )
        {
            if (SwapVsFut == K_FUT) 
               zcSize++;
            else if ( Futures.Elt(i, 1) < SwapRates.Elt(fxPayFreq, 1) )  
               zcSize++;
        }

        if  (MMVsFut == K_FUT) 
        {
            if (SwapVsFut == K_FUT) 
               zcSize++;

            else if ( Futures.Elt(i, 1) < SwapRates.Elt(fxPayFreq, 1) )  
               zcSize++;
        }
    }

    for (i = fxPayFreq; i <= SizeS; i++)
    {
        if (SwapVsFut != K_FUT) 
           zcSize++;

        else if ( SwapRates.Elt(i, 1) > Futures.Elt(SizeF-1, 1) ) 
           zcSize++;
    }

    // MA : delete before !!

    if (itsDateTerms)
       delete itsDateTerms;

    if (itsYearTerms)
       delete itsYearTerms;

    if (itsZeroRates)
       delete itsZeroRates;

    if (itsDiscountFactors)
       delete itsDiscountFactors;

    itsDateTerms = new ARM_Vector (zcSize, 0.0);
    itsYearTerms = new ARM_Vector (zcSize, 0.0);
    itsZeroRates = new ARM_Vector (zcSize, 0.0);
    itsDiscountFactors = new ARM_Vector(zcSize, 0.0);

    zcSize=0;

    for (i = 0; i < SizeM; i++)
    {
        if (MMVsFut != K_FUT) 
        {
            itsDateTerms->Elt(zcSize) = 
                MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();     
            
            itsDiscountFactors->Elt(zcSize) = MMRates.Elt(i, 0);     
            
            zcSize++;
        }
        else if ( MMRates.Elt(i, 1) < Futures.Elt(0, 1) )
        {
            itsDateTerms->Elt(zcSize) = 
                MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();     
            
            itsDiscountFactors->Elt(zcSize) = MMRates.Elt(i, 0);     
            
            zcSize++;
        }
    }

    for (i = 0; i < SizeF; i++)
    {
        if  (MMVsFut == K_FUT) 
        {
            if (SwapVsFut == K_FUT) 
            {
                itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
                zcSize++;
            }
            else if (Futures.Elt(i, 1) < SwapRates.Elt(fxPayFreq, 1) )          
            {
                itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
                zcSize++;
            }
        }
        else if ( MMRates.Elt(SizeM-1, 1) < Futures.Elt(i, 1) )
        {
            if (SwapVsFut == K_FUT) 
            {
                itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
                zcSize++;
            }
            else if ( Futures.Elt(i, 1) < SwapRates.Elt(fxPayFreq, 1) ) 
            {
                itsDateTerms->Elt(zcSize) = 
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();     
                
                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);     
                
                zcSize++;
            }
        }

    }

    for (i = fxPayFreq; i <= SizeS; i++)
    {
        if (SwapVsFut != K_FUT)         
        {
            itsDateTerms->Elt(zcSize) = 
                SwapRates.Elt(i, 1) - itsAsOfDate.GetJulian();
            
            itsDiscountFactors->Elt(zcSize) = SwapRates.Elt(i, 0);     
            
            zcSize++;
        }
        else if ( SwapRates.Elt(i, 1) > Futures.Elt(SizeF-1, 1) )
        {
            itsDateTerms->Elt(zcSize) = 
                SwapRates.Elt(i, 1) - itsAsOfDate.GetJulian();
            
            itsDiscountFactors->Elt(zcSize) = SwapRates.Elt(i, 0);
            
            zcSize++;
        }
    }

    for (i = 0; i < zcSize; i++)
    {
        itsYearTerms->Elt(i) = itsDateTerms->Elt(i)/365.0; 
        itsZeroRates->Elt(i) = -100.0*log(itsDiscountFactors->Elt(i))/
                               itsYearTerms->Elt(i);
    }
}



void ARM_ZeroCurve::CptFuturesZeroRatesGen(ARM_Matrix& MMRates, int SizeM,
                                           ARM_Matrix& Future, int SizeF,
                                           ARM_Matrix& Swap,
                                           int MMVsFut, int Cont_Lin)
{
    ARM_Date D2_Cash;
    ARM_Date Ctr_Ech1, Ctr_Ech2;
    int i, j, first;
    double N6;
    double Df2, Df_Ctr, inter, TmP;



    // Procedure de tri des taux Futures sur la date d'echeance (2)

    for (i = 1; i < SizeF; i++)
    {
        for (j = i; j < SizeF ; j++)
        {  
            if ( Future.Elt(i, 2) > Future.Elt(j, 2) )
            {
               TmP = Future.Elt(j, 1);
               Future.Elt(j, 1) = Future.Elt(i, 1);
               Future.Elt(i, 1) = TmP;
                
               TmP = Future.Elt(j,0);
               Future.Elt(j, 0) = Future.Elt(i, 0);
               Future.Elt(i, 0) = TmP;
            
               TmP = Future.Elt(j,2);
               Future.Elt(j, 2) = Future.Elt(i, 2);
               Future.Elt(i, 2) = TmP;
            }
        }
    }

    // Calcul des Df a partir des taux Futures


    first = 1;

    for (i = 0; (( i < SizeF ) && ( MMVsFut == 1 )); i++)
    {
        if ( Future.Elt(i, 1) < MMRates.Elt(SizeM-1, 1) ) 
           first = i;
    }

    i = SizeM-1;


    while  ( MMRates.Elt(i, 1) > Future.Elt(first, 1) && i>0) 
    {
        i--;
    }
    

    if ( i == SizeM-1 )
    {
       D2_Cash = (ARM_Date) Swap.Elt(1, 1); // XX Contient le 1er tx de
                                            // swap 1Y (ou 2Y
                                            // XX s'il y a le 12M )
        
       TmP = MMRates.Elt(0, 0); // DF du 1 ou 2 Jours
       inter = DaysBetweenDates(K30_360, (ARM_Date) MMRates.Elt(0, 1) ,
                                         (ARM_Date) Swap.Elt(1, 1)  );
       Df2 = 1.0 + Swap.Elt(1, 0)*inter/360.0/100.0;
       Df2 = TmP/Df2;
    }
    else
    {
       D2_Cash = (ARM_Date) MMRates.Elt(i+1, 1);
       Df2 = MMRates.Elt(i+1, 0);
    }

    Ctr_Ech1= (ARM_Date) Future.Elt(first, 1);

    Df_Ctr=Interpol(MMRates.Elt(i, 1),
                    D2_Cash.GetJulian(),
                    Future.Elt(first, 1),
                    MMRates.Elt(i, 0),
                    Df2,
                    Cont_Lin,
                    KACTUAL_360);
    
    Future.Elt(first-1, 0) = Df_Ctr;

    Future.Elt(first-1, 1) = Future.Elt(first, 1);;
    
    for (i=0; i < first-1; i++) Future.Elt(i, 1) = 0.0;

    i = first;

    while ( i < SizeF )
    {
        Ctr_Ech1 = (ARM_Date) Future.Elt(i, 1);
        Ctr_Ech2 = (ARM_Date) Future.Elt(i, 2);

        double julEch1 = Future.Elt(i, 1);

    
        N6 = DaysBetweenDates (KACTUAL_360, Ctr_Ech1, Ctr_Ech2);

        double prevDF = 0.0;
        int found = 0;

        // recherche du DF a la date julEch1

        for (int k = first-1; k < i; k++)
        {
            if ( (Future.Elt(k,1) == julEch1) && (!found) ) 
            {
                found = 1;
                prevDF = Future.Elt(k,0);
            }
        }

        if (!found)
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Invalid Data");
        }

        Df_Ctr = prevDF / ( 
                            ( 1 + (Future.Elt(i, 0)) / 100.0 
                            * ( N6  / 360.0 ))
                         );

        Future.Elt(i, 0) = Df_Ctr;
        Future.Elt(i, 1) = Ctr_Ech2.GetJulian();

        i++;
    }
}


















            /**************************************/
            /*                                    */
            /*         TOY METHODS                */
            /*                                    */
            /**************************************/



//*******************************************************************************
//*** Author:    Herve FRANCOIS
//***
//*** Goal:  Constructor to Integrate the TOY on the american market
//***
//*** Detail:    ARM_ZeroCurve::ARM_ZeroCurve(AseOfDate,
//***                                         Term,
//***                                         Values,
//***                                         Priority to MM ou au Fut
//***                                         )
//***
//*** Date:        09/01/02
//***
//*******************************************************************************
           

ARM_ZeroCurve::ARM_ZeroCurve(ARM_Date& asOf,
                             ARM_CRV_TERMS& Terms,
                             ARM_Vector* mktData,
                             int MMVsFut,
                             int SwapVsFut,
                             int raw,
                             int Cont_Lin,
                             int Frequency,
                             ARM_Currency* ccy)
{
    SetName(ARM_ZERO_CURVE);

    Init();

    itsAsOfDate = asOf;

    // Set variables Market Data

   itsMktData = new ARM_MarketData(Terms, mktData, 
                                   MMVsFut, SwapVsFut, raw, Cont_Lin, 2);

    SetCurrencyUnit(ccy);

    ZCFromMarketRatesTOY(Terms, mktData, MMVsFut, SwapVsFut, 
                         raw, Cont_Lin, ccy, Frequency);
}



void ARM_ZeroCurve::ZCFromMarketRatesTOY(ARM_CRV_TERMS& Terms,
                                         ARM_Vector* data,
                                         int MMVsFut,
                                         int SwapVsFut,
                                         int raw,
                                         int Cont_Lin,
                                         ARM_Currency* ccy,
                                         int Frequency)
{
    // Matrix have 3 columns
    // In first column: rates next DF
    // In second column: dates
    // In third column: rates

    int SizeM=0, SizeFRA=0, SizeF=1, SizeS=0, i, Nb, cpt = 0;
    int Tx_Spot=0;
    ARM_Date matDate, spotDate, swapStartDate,date, dateContrat, dateUnder;

    ARM_Date matuDateFormat;

    int CurrentMonth;

    ARM_Matrix MMRates(100, 3, 0.0);
    ARM_Matrix Futures(500, 3, 0.0);
    ARM_Matrix SwapRates(500, 3, 0.0);
    ARM_Matrix FRARates(500, 3, 0.0);
    ARM_Matrix DF(500, 2, 0.0);
    ARM_Matrix Temp(2, 1, 0.0);


    SetCurrencyUnit(ccy);

    char* ccyName = ccy->GetCcyName();
    int fwdRule = ccy->GetFwdRule();
    int spotDays = ccy->GetSpotDays();
    int fxPayFreq = ccy->GetFixedPayFreq();

    ARM_INDEX_TYPE idxType = GetDefaultIndexFromCurrency(ccyName);

    ARM_Currency tmpCurrency(ccyName);

    ccyName = tmpCurrency.GetPayCalName(idxType);

    

    double StartNumber = -1;

    char matu;
    int year;
    int month = 0;

    int lastMMTerm;

    spotDate = itsAsOfDate;

    swapStartDate = itsAsOfDate;

    swapStartDate.NextBusinessDay(spotDays, ccyName);

    CurrentMonth = spotDate.GetMonth();

    matu = 'X';


    i = 0;

    while ( Terms[i][0] != 'X' )
    {
        if ( data->Elt(i) < 1.0E-6 )
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Non zero rates are required for building zero curve");
        }

        // si 1er caract alphanum -> Futures
        if (isalpha(Terms[i][0]))
        {
           GetMonthYearFromExpiryDate(Terms[i], &month, &year);
           matu='Z';
        }
        else // sinon MM ou Swap
        {
           if (strlen (Terms[i]) > 3 ) // FRA's case: MM/DD/YYtoMM/DD/YY
           {
              char c;

              int d1,d2;
              int m1,m2;
              int y1,y2;

              if ( strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0 )
              {
                 sscanf(Terms[i], "%d%c%d%c%2d%c%d%c%d%c%2d", &m1, &c, &d1,
                        &c, &y1, &c, &m2, &c, &d2, &c, &y2 );
              }
              else
              {
                 sscanf(Terms[i], "%d%c%d%c%2d%c%d%c%d%c%2d", &d1, &c, &m1,
                        &c, &y1, &c, &d2, &c, &m2, &c, &y2 );
              }

              y1 = 2000+y1;
              ARM_Date tmpDate(d1, m1, y1);
              dateContrat = tmpDate.GetJulian();// chge format from MM/JJ/YY to julian

              y2 = 2000+y2;
              ARM_Date tmpDate2(d2, m2, y2);
              dateUnder = tmpDate2.GetJulian();// chge format from MM/JJ/YY to julian

              matu = 'R';
           }
           else // case: MM or Swap
           {
              sscanf(Terms[i], "%d%c",&Nb, &matu);
              matu = toupper(matu);
           }
        }

        // Conversion  Termes char -> Termes double

        if ( matu == 'D' ) // Ex : "1D"
        {
            matDate = itsAsOfDate;

            spotDate.NextBusinessDay(spotDays, ccyName);


            MMRates.Elt(SizeM, 0) = data->Elt(i);
            MMRates.Elt(SizeM, 2) = data->Elt(i);


            MMRates.Elt(SizeM, 1) = spotDate.GetJulian();

            if ( spotDays < 1 )
            {
               matDate.NextBusinessDay(Nb, ccyName);
               MMRates.Elt(SizeM, 1) = matDate.GetJulian();
            }

            Tx_Spot = 1;

            SizeM++;
        }
        else if ( matu == 'W' )
        {   //  Ex : "1W"

            matDate = spotDate;

            matDate.AddDays(Nb*7);

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

            MMRates.Elt(SizeM, 0) = data->Elt(i);
            MMRates.Elt(SizeM, 2) = data->Elt(i);

            MMRates.Elt(SizeM, 1) = matDate.GetJulian();

            SizeM++;
        }
        else if (( matu == 'M' )
                 ||
                 ( matu == 'S' )
                )
        {   //  Ex : "9M"

            if ( matu == 'M' )
            {
               matDate = spotDate;

               matDate.AddMonths(Nb, GOTO_END_OF_MONTH);

               lastMMTerm = Nb;
            }
            else
            {
               matDate = matuDateFormat;


               if (((matuDateFormat-itsAsOfDate)/30.0) < 1 )
               {
                  lastMMTerm = 1;
               }
               else
               {
                  lastMMTerm = int((matuDateFormat-itsAsOfDate)/30.0);
               }
            }

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

            MMRates.Elt(SizeM, 0) = data->Elt(i);
            MMRates.Elt(SizeM, 2) = data->Elt(i);

            MMRates.Elt(SizeM, 1) = matDate.GetJulian();

            SizeM++;
        }
        else if ( matu == 'Y')  // ->implicitement ce sont des taux de swap
        {
            //  Ex : "15Y"

            matDate = swapStartDate;;

            //matDate.AddMonths(12*Nb/fxPayFreq);
            matDate.AddYears(Nb);

            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

           if ( lastMMTerm < 12 )
           {
              SwapRates.Elt(Nb*fxPayFreq, 0) = data->Elt(i);
              SwapRates.Elt(Nb*fxPayFreq, 2) = data->Elt(i);

              SwapRates.Elt(Nb*fxPayFreq, 1) = matDate.GetJulian();

              if ( StartNumber == -1 )
              {
                 StartNumber = Nb*fxPayFreq;
              }
           }
           else if (lastMMTerm == 12 && Nb>1)
           {
              SwapRates.Elt((Nb-1)*fxPayFreq, 0) = data->Elt(i);
              SwapRates.Elt((Nb-1)*fxPayFreq, 2) = data->Elt(i);

              SwapRates.Elt((Nb-1)*fxPayFreq, 1) = matDate.GetJulian();

              if ( StartNumber == -1 )
              {
                 StartNumber = (Nb-1)*fxPayFreq;
              }
           }

            SizeS++;
        }                      
        else if (matu == 'R' )  // FRA's case
        {
            FRARates.Elt(SizeFRA, 0) = data->Elt(i);
            FRARates.Elt(SizeFRA, 1) = dateContrat.GetJulian();
            FRARates.Elt(SizeFRA, 2) = dateUnder.GetJulian();

            SizeFRA++;
        }
        else if ( month > 0) // traitement du segment Futures
        {
            matDate.ChangeDate(1, month, year);

            matDate.PiborDelivery();

            // tester si NextDate < = > itsAsOfDate
            if ( matDate.GetJulian() <= itsAsOfDate.GetJulian() + 2)
            {
               char BUF[100];

               sprintf(BUF,
                      "Futures contract expiring : %s  %d",Terms[i], year);

               throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, BUF);
            }

            Futures.Elt(SizeF, 0) = 100.0 - data->Elt(i);
            Futures.Elt(SizeF, 2) = 100.0 - data->Elt(i);

            Futures.Elt(SizeF, 1) = matDate.GetJulian();

            SizeF++;
        }
        else
           break;

        i++;
    }
        
    if ( Tx_Spot == 0 )  // pas de taux spot
    {
        ARM_Matrix spotRate(1, 2, 0.0);

        spotRate.Elt(0, 1) = itsAsOfDate.GetJulian();

        MMRates = spotRate & MMRates;

        SizeM += 1;
    }

    if ( SizeF == 1 )  // pas de contrats futures ->Forcer les options
    {
       SwapVsFut = K_SWAP;

       MMVsFut = K_MM;
    }


/* OLD CODE
    if (( SwapRates.Elt(fxPayFreq, 0) < 0.0001 ) && ( SwapVsFut == 1 ))
    {
       char BUF[100];

       sprintf(BUF, "1 year Swap rate not found");
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, BUF);
    }
*/

    if (( SwapRates.Elt(fxPayFreq*2, 0) < 0.0001 ) && ( SwapVsFut == 1 ))
    {
       char BUF[100];

       sprintf(BUF, "2 years Swap rate not found");
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, BUF);
    }

    for (i = 0; i < SwapRates.GetNumLines(); i++)
    {
       if ( SwapRates.Elt(i, 1) > 0.0 )
          SizeS = i; // sizeS = maturite max de la courbe
    }

    for (i = 1; i <=SizeS; i++)
    {
        if ( SwapRates.Elt(i, 1) < 0.0001 )
        {
           matDate = swapStartDate;

           matDate.AddMonths((i)*12/fxPayFreq, GOTO_END_OF_MONTH);

           matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

           SwapRates.Elt(i, 1) = matDate.GetJulian();
        }
    }
  
    // calcul des taux ZC implicites dans chaque segment du marche
    CptMMZeroRates(&MMRates, SizeM);

    if ( SizeFRA >= 1 )
    {
       SizeFRA = CptDF(MMRates, SizeM, FRARates, SizeFRA, Cont_Lin);

       if (SizeF != 1)
          CptFuturesZeroRatesTOY(MMRates, SizeM, FRARates, SizeFRA, 
                                 &Futures, SizeF,
                                 SwapRates,
                                 MMVsFut, Cont_Lin);

        if ( SizeS > 0 )
           CptSwapZeroRatesTOY(MMRates, SizeM, FRARates, SizeFRA, 
                               &SwapRates, SizeS, fxPayFreq,
                               Futures, SizeF, MMVsFut,
                               SwapVsFut, raw, Cont_Lin);

        SizeFRA = CptTOYFRA(FRARates, SizeFRA, &Futures, SizeF,  Cont_Lin);
    }
    else
    {
       if ( SizeF != 1 )
          CptFuturesZeroRates(MMRates, SizeM, &Futures, SizeF, 
                              SizeS,
                              SwapRates,
                              MMVsFut, Cont_Lin);

        if ( SizeS > 0 )
           CptSwapZeroRates(MMRates, SizeM, &SwapRates, SizeS, fxPayFreq,
                            Futures, SizeF, MMVsFut, SwapVsFut, raw, Cont_Lin);

        if ( SizeF != 1 )
           SizeF = CptTOY( MMRates, SizeM, &Futures, SizeF,  Cont_Lin);
    }

    // Preparation du tableau des taux ZC correspondant a generateCurve

    while ( SwapRates.Elt(cpt,2) == 0 )
    {
        cpt++;
    }

    int zcSize = 0;

    if ( SizeFRA > 0 )
    {
       for (i = 0; i < SizeFRA; i++)
       {
           zcSize++;
       }
    }
    else
    {
       for (i = 0; i < SizeM; i++)
       {
           if (MMVsFut != K_FUT)
              zcSize++;
           else if ( MMRates.Elt(i, 1) < Futures.Elt(0, 1) )
              zcSize++;
       }
    }

    for (i = 0; i < SizeF; i++)
    {
        if ( (MMVsFut != K_FUT) &&
             (MMRates.Elt(SizeM-1, 1) < Futures.Elt(i, 1)) )
        {
            if (SwapVsFut == K_FUT)
               zcSize++;
            else if ( Futures.Elt(i, 1) < SwapRates.Elt(cpt, 1) )
               zcSize++;
        }

        if  (MMVsFut == K_FUT)
        {
            if (SwapVsFut == K_FUT)
               zcSize++;

            else if ( Futures.Elt(i, 1) < SwapRates.Elt(cpt, 1) )
               zcSize++;
        }
    }


    for (i = cpt; i <= SizeS; i++)
    {
        if (SwapVsFut != K_FUT)
           zcSize++;

        else if ( SwapRates.Elt(i, 1) > Futures.Elt(SizeF-1, 1) )
           zcSize++;
    }

    // MA : delete before !!

    if (itsDateTerms)
       delete itsDateTerms;

    if (itsYearTerms)
       delete itsYearTerms;

    if (itsZeroRates)
       delete itsZeroRates;

    if (itsForwardRates)
       delete itsForwardRates;

    if (itsForwardDF)
       delete itsForwardDF;

    if (itsForwardDates)
       delete itsForwardDates;

    if (itsDiscountFactors)
       delete itsDiscountFactors;


    itsDateTerms    = new ARM_Vector (zcSize, 0.0);
    itsYearTerms    = new ARM_Vector (zcSize, 0.0);
    itsZeroRates    = new ARM_Vector (zcSize, 0.0);
    itsForwardRates = new ARM_Vector (zcSize, 0.0);
    itsForwardDates  = new ARM_Vector (zcSize, 0.0);
    itsForwardDF    = new ARM_Vector (zcSize, 0.0);
    itsDiscountFactors = new ARM_Vector(zcSize, 0.0);


    zcSize = 0;

   
    if (SizeFRA > 0)
    {
        for (i = 0; i < SizeFRA; i++)
        {
            itsDateTerms->Elt(zcSize) =
            FRARates.Elt(i, 1) - itsAsOfDate.GetJulian();

            itsDiscountFactors->Elt(zcSize) = FRARates.Elt(i, 0);

            double value = FRARates.Elt(i, 0);
            value = itsDiscountFactors->Elt(zcSize);

            DF.Elt(zcSize, 0) = itsDiscountFactors->Elt(zcSize);
            DF.Elt(zcSize, 1) = FRARates.Elt(i, 1);

            zcSize++;

        }
    }
    else
    {
        for (i = 0; i < SizeM; i++)
        {
            if (MMVsFut != K_FUT)
            {
                itsDateTerms->Elt(zcSize) =
                MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();

                itsDiscountFactors->Elt(zcSize) = MMRates.Elt(i, 0);

                DF.Elt(zcSize, 0) = itsDiscountFactors->Elt(zcSize);
                DF.Elt(zcSize, 1) =MMRates.Elt(i, 1);

                zcSize++;
            }
            else if ( MMRates.Elt(i, 1) < Futures.Elt(0, 1) )
            {
                double tmp = MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();
                tmp = MMRates.Elt(i, 1) ;
                itsDateTerms->Elt(zcSize) =
                MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();

                itsDiscountFactors->Elt(zcSize) = MMRates.Elt(i, 0);

                DF.Elt(zcSize, 0) = itsDiscountFactors->Elt(zcSize);
                DF.Elt(zcSize, 1) = MMRates.Elt(i, 1);

                zcSize++;
            }
        }
    }
   
    for (i = 0; i < SizeF; i++)
    {
        if  (MMVsFut == K_FUT)
        {
            if (SwapVsFut == K_FUT)
            {
                itsDateTerms->Elt(zcSize) =
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();

                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);

                DF.Elt(zcSize, 0) = itsDiscountFactors->Elt(zcSize);
                DF.Elt(zcSize, 1) = Futures.Elt(i, 1);

                zcSize++;
            }
            else if (Futures.Elt(i, 1) < SwapRates.Elt(cpt, 1) )
            {
                itsDateTerms->Elt(zcSize) =
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();

                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);

                DF.Elt(zcSize, 0) = itsDiscountFactors->Elt(zcSize);
                DF.Elt(zcSize, 1) = Futures.Elt(i, 1);

                zcSize++;
            }
        }
        else if ( MMRates.Elt(SizeM-1, 1) < Futures.Elt(i, 1) )
        {
            if (SwapVsFut == K_FUT)
            {
                itsDateTerms->Elt(zcSize) =
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();

                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);

                DF.Elt(zcSize, 0) = itsDiscountFactors->Elt(zcSize);
                DF.Elt(zcSize, 1) = Futures.Elt(i, 1);

                zcSize++;
            }
            else if ( Futures.Elt(i, 1) < SwapRates.Elt(cpt, 1) )
            {
                itsDateTerms->Elt(zcSize) =
                    Futures.Elt(i, 1) - itsAsOfDate.GetJulian();

                itsDiscountFactors->Elt(zcSize) = Futures.Elt(i, 0);

                DF.Elt(zcSize, 0) = itsDiscountFactors->Elt(zcSize);
                DF.Elt(zcSize, 1) = Futures.Elt(i, 1);

                zcSize++;
            }
        }
    }

    for (i = cpt; i <= SizeS; i++)
    {
        if (SwapVsFut != K_FUT)
        {
            itsDateTerms->Elt(zcSize) =
                SwapRates.Elt(i, 1) - itsAsOfDate.GetJulian();

            double tmp = SwapRates.Elt(i, 0);
            itsDiscountFactors->Elt(zcSize) = SwapRates.Elt(i, 0);

            DF.Elt(zcSize, 0) = itsDiscountFactors->Elt(zcSize);
            DF.Elt(zcSize, 1) = SwapRates.Elt(i, 1);

            zcSize++;
        }
        else if ( SwapRates.Elt(i, 1) > Futures.Elt(SizeF-1, 1) )
        {
            itsDateTerms->Elt(zcSize) =
                SwapRates.Elt(i, 1) - itsAsOfDate.GetJulian();

            itsDiscountFactors->Elt(zcSize) = SwapRates.Elt(i, 0);

            DF.Elt(zcSize, 0) = itsDiscountFactors->Elt(zcSize);
            DF.Elt(zcSize, 1) = SwapRates.Elt(i, 1);

            zcSize++;
        }
    }

    for (i = 0;i < zcSize ; i++)
    {
        for (int j = i+1;j < zcSize; j++)
        {
            if ( DF.Elt(i, 1) > DF.Elt(j, 1) )
            {
               double TmP = DF.Elt(j, 1);
               DF.Elt(j, 1) = DF.Elt(i, 1);
               DF.Elt(i, 1) = TmP;

               TmP = DF.Elt(j, 0);
               DF.Elt(j, 0) = DF.Elt(i, 0);
               DF.Elt(i, 0) = TmP;
            }
        }
    }


    // Sort
    for (i = 0;i < zcSize ; i++)
    {
        for (int j = i;j < zcSize; j++)
        { 
            if ( itsDateTerms->Elt(i) > itsDateTerms->Elt(j))
            {
               double TmP = itsDateTerms->Elt(j);
               itsDateTerms->Elt(j) = itsDateTerms->Elt(i);
               itsDateTerms->Elt(i) = TmP;

               TmP = itsDiscountFactors->Elt(j);
               itsDiscountFactors->Elt(j) = itsDiscountFactors->Elt(i);
               itsDiscountFactors->Elt(i) = TmP;
            }
        }
    }

    // Add of forward rates
    i=1;

    date = spotDate;
    ARM_Date date3M;

    itsForwardDates->Elt(i) = date.GetJulian();

    Temp.Elt(0,0) = InterpolDF(date.GetJulian(), DF, zcSize,
                               MMRates, SizeM, Cont_Lin);

    date3M.AddMonths(3, GOTO_END_OF_MONTH);

    if (date3M < DF.Elt (zcSize-1, 1))
    {
       Temp.Elt(1,0) = InterpolDF(date3M.GetJulian(), DF, 
                                  zcSize, MMRates, SizeM, Cont_Lin);

       itsForwardRates->Elt(i) = 100*(Temp.Elt(0,0)
                                 /Temp.Elt(1,0)-1)*360/(date3M-date);
    }
   
    date = itsAsOfDate;
    date.AddMonths(3, GOTO_END_OF_MONTH);

    while (date.IsBusinessDay(ccyName) == 0)
    {
        date.NextBusinessDay(1, ccyName);
    }

    i=2;   

    while (( i < itsForwardDates->GetSize() )
           &&
           ( date.GetJulian() < DF.Elt (zcSize-1, 1)))
    {
        itsForwardDates->Elt(i) = date.GetJulian();

        Temp.Elt(0,0) = InterpolDF(date.GetJulian(), DF, 
                                   zcSize, MMRates, SizeM, Cont_Lin);

        date3M = date;
        date3M.AddMonths(3, GOTO_END_OF_MONTH);

        while (date3M.IsBusinessDay(ccyName) == 0)
        {
            date3M.NextBusinessDay(1, ccyName);
        }


        if (date3M > DF.Elt(zcSize-1, 1))
        {
            break;
        }
        else
        {
           Temp.Elt(1,0) = InterpolDF(date3M.GetJulian(),DF,zcSize,
                                      MMRates, SizeM, Cont_Lin);

           itsYearTerms->Elt(i) = itsForwardDates->Elt(i)/365.0;

           itsForwardRates->Elt(i) = 100*(Temp.Elt(0,0)
                                      /Temp.Elt(1,0)-1)*360.0/(date3M-date);

           itsForwardDF->Elt(i) = 1.0/exp(itsForwardRates->Elt(i)/100.0
                                          *itsYearTerms->Elt(i)/100.0);

           date = itsAsOfDate;

           date.AddMonths(3*i, GOTO_END_OF_MONTH);

           while (date.IsBusinessDay(ccyName) == 0)
           {
               date.NextBusinessDay(1, ccyName);
           }

           i++;
        }
    }

    itsForwardDates->Elt(0) = i;

    for (i=0; i < zcSize; i++)
    {
        itsYearTerms->Elt(i) = itsDateTerms->Elt(i)/365.0;
        itsZeroRates->Elt(i) = 100.0*log(1/itsDiscountFactors->Elt(i))/
                               itsYearTerms->Elt(i);

        double Tmp2 = itsZeroRates->Elt(i);
    }

    delete ccyName;
}



int ARM_ZeroCurve::CptDF(ARM_Matrix& MMRates, int SizeM,
                         ARM_Matrix& FRARates, int SizeFRA,
                         int Cont_Lin)
{
    ARM_Date matDate, DFbeforeEndDate, spotDate = itsAsOfDate, tmpDate, dateDF; 
    double DFendDate, DFstartDate, DBD360, Fwd;
    int SizeDF = 0;
    ARM_Matrix DF(500, 2, 0.0);
    int i,j;
    double TmP;

    // initialization
    ARM_Currency* ccy = GetCurrencyUnit();
    char* ccyName = ccy->GetCcyName();
    int dayCount = ccy->GetMMDayCount();


    ARM_INDEX_TYPE idxType = GetDefaultIndexFromCurrency(ccyName);

    ARM_Currency tmpCurrency(ccyName);

    ccyName = tmpCurrency.GetPayCalName(idxType);

    // Sort FRA

    for (i = 0;i < SizeFRA ; i++)
    {
        for (int j = i;j < SizeFRA; j++)
        {  
            if ( FRARates.Elt(i, 1) > FRARates.Elt(j, 1) )
            {

               double TmP = FRARates.Elt(j, 2);
               FRARates.Elt(j, 2) = FRARates.Elt(i, 2);
               FRARates.Elt(i, 2) = TmP;

               TmP = FRARates.Elt(j, 1);
               FRARates.Elt(j, 1) = FRARates.Elt(i, 1);
               FRARates.Elt(i, 1) = TmP;

               TmP = FRARates.Elt(j, 0);
               FRARates.Elt(j, 0) = FRARates.Elt(i, 0);
               FRARates.Elt(i, 0) = TmP;
            }
        }
    }

    // Pick up the 1 or 2 days DF            
    if ( (MMRates.Elt(0, 1) 
         - itsAsOfDate.GetJulian()) < 7 ) //Suppose to have 1D or 2D input
    {
        DF.Elt(SizeDF, 0) =MMRates.Elt(0, 0); //DF 1 or 2 Days 
        DF.Elt(SizeDF, 1) =MMRates.Elt(0, 1);
    
        SizeDF++;


        // DF 3 months
    
        // Look for the nearest input of 3 months

        int i = 0;
        matDate = MMRates.Elt(0, 1);

        while (((MMRates.Elt(i, 1) - matDate.GetJulian()) < 100 )&&(i<SizeM))
        {
            DF.Elt(SizeDF, 0) = MMRates.Elt(i, 2); 
            DF.Elt(SizeDF, 1) = MMRates.Elt(i, 1);
            i++;
        }


        // determination of the DF
        DBD360 = DaysBetweenDates(KACTUAL_360, (ARM_Date) DF.Elt(0, 1),
                                  (ARM_Date) DF.Elt(i-1, 1))/360.0;

        DF.Elt(SizeDF, 0) = DF.Elt(0, 0)/(1.0+DF.Elt(SizeDF, 0)/100.0*DBD360); 
    
        SizeDF++;


        // DF 6, 9, ... Months    

        matDate = DF.Elt(SizeDF-1, 1);
		
        for (i = 1; i < 8; i++)
        {
            matDate.AddMonths(3, GOTO_END_OF_MONTH);
            matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

            double temp  = DF.Elt(SizeDF, 1);

            if (matDate.GetJulian() < FRARates.Elt(SizeFRA-1,2))
            {
               DF.Elt(SizeDF, 1) = matDate.GetJulian();

               DBD360 = DaysBetweenDates(KACTUAL_360, 
                                         (ARM_Date) DF.Elt(SizeDF-1, 1), 
                                         (ARM_Date) DF.Elt(SizeDF, 1))/360.0;
        
               // determination of the forward 3 months thanks to 
               // a continuous interpolation
               Fwd = InterpolForward(DF.Elt(SizeDF-1,1),FRARates,SizeFRA,
                                     MMRates, SizeM,
                                     -1);
    
               DF.Elt(SizeDF, 0) = DF.Elt(SizeDF-1,0) / ( 1 + Fwd /100.0 * DBD360);

               SizeDF++;
            }
        }

        // DF at the end of the last FRA

        DF.Elt(SizeDF, 1) = FRARates.Elt(SizeFRA-1,2);

        double DBD365 = DaysBetweenDates(KACTUAL_360, 
                                         (ARM_Date) DF.Elt(SizeDF-1,1), 
                                         (ARM_Date) DF.Elt(SizeDF,1))/365.0;
        
        double DBDFRA360 = DaysBetweenDates(KACTUAL_360, 
                                    (ARM_Date) FRARates.Elt(SizeFRA-1,1),
                                    (ARM_Date) FRARates.Elt(SizeFRA-1,2))/360.0;
        
        double DBDFRA365 = DaysBetweenDates(KACTUAL_360, 
                              (ARM_Date) FRARates.Elt(SizeFRA-1,1), 
                              (ARM_Date) FRARates.Elt(SizeFRA-1,2))/365.0;
        
        DFendDate = DF.Elt(SizeDF-1,0)*exp(-log(1+FRARates.Elt(SizeFRA-1, 0)
                         /100*DBDFRA360)/DBDFRA365*DBD365);

        DF.Elt(SizeDF, 0) = DFendDate;
    
        SizeDF++;
    

        // DF at the start of the last FRA

        DF.Elt(SizeDF, 1) = FRARates.Elt(SizeFRA - 1, 1);

        DBD360 = DaysBetweenDates(KACTUAL_360, (ARM_Date) DF.Elt(SizeDF,1), 
                              (ARM_Date) DF.Elt(SizeDF-1,1))/360.0;
    
        DFstartDate = DFendDate * ( 1 + FRARates.Elt(SizeFRA-1, 0)/100 * DBD360);

        DF.Elt(SizeDF, 0) = DFstartDate;
    
        SizeDF++;

        dateDF = DF.Elt(SizeDF-1, 1);

        // 1st sort
        for (i = 0;i < SizeDF ; i++)
        {
            for (j = i+1;j < SizeDF; j++)
            {
                if ( DF.Elt(i, 1) > DF.Elt(j, 1) )
                {
                   TmP = DF.Elt(j, 1);
                   DF.Elt(j, 1) = DF.Elt(i, 1);
                   DF.Elt(i, 1) = TmP;

                   TmP = DF.Elt(j, 0);
                   DF.Elt(j, 0) = DF.Elt(i, 0);
                   DF.Elt(i, 0) = TmP;
                }
            }
        }

        int index = SizeDF - 1;

        // Daily DF
        while (dateDF > (DF.Elt(0, 1) + 1))
        {

            dateDF = dateDF.PreviousBusinessDay(ccyName);

            int k = index;
            while ( (k > 0) && (dateDF < DF.Elt(k, 1)) )// Looking for existing DF
            {
                k--;
            }

            if ( dateDF == DF.Elt(k, 1))
            {
               dateDF = dateDF.PreviousBusinessDay(ccyName);
            }

            if (dateDF.GetJulian() < MMRates.Elt(0,1))
            {
                break;
            }

            tmpDate = dateDF;

            tmpDate.AddMonths(3, GOTO_END_OF_MONTH);
            tmpDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

            if (tmpDate.IsBusinessDay(ccyName) == 0)
            {
               tmpDate = tmpDate.NextBusinessDay(1, ccyName);
            }

            DBD360 = DaysBetweenDates(KACTUAL_360, dateDF, tmpDate)/360.0;

            Fwd = InterpolForward(dateDF.GetJulian(),
                              FRARates, SizeFRA,
                              MMRates,  SizeM,
                              -1);// method basic linear interpolation

            DFendDate = InterpolDF(tmpDate.GetJulian(),DF,SizeDF,
                               MMRates, SizeM,
                               Cont_Lin);

            DFstartDate = DFendDate * (1 +  Fwd/100. * DBD360);

            index = k+1;

            k = SizeDF;

            while (k>index)
            {
                DF.Elt(k, 1) = DF.Elt(k-1, 1);
                DF.Elt(k, 0) = DF.Elt(k-1, 0);
                k--;
            }

            DF.Elt(index, 1) = dateDF.GetJulian();
            DF.Elt(index, 0) = DFstartDate;

            SizeDF++;
        }

        for (i=0; i < SizeDF; i++)
        {
            FRARates.Elt(i, 0) = DF.Elt(i, 0); 
            FRARates.Elt(i, 1) = DF.Elt(i, 1);
        }

        SizeFRA = SizeDF;
		
		/// the delete has to be done here
		/// since we exit
	    return SizeFRA;
    }

	/// other cases are errors
	delete ccyName;
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"MMRates are invalid");
}



double ARM_ZeroCurve::InterpolForward(double Date, ARM_Matrix& FRARates, 
                                      int SizeFRA,
                                      ARM_Matrix& MMRates, int SizeM,
                                      int Cont_Lin)
{
    int i,j;

    double check = 0, DF1, DF1Date, DF2, DF2Date, Forward;

    ARM_Date matDate;


    ARM_Currency* ccy = GetCurrencyUnit();
    char* ccyName = ccy->GetCcyName();
    int dayCount = ccy->GetMMDayCount();


    ARM_INDEX_TYPE idxType = GetDefaultIndexFromCurrency(ccyName);

    ARM_Currency tmpCurrency(ccyName);

    ccyName = tmpCurrency.GetPayCalName(idxType);


    // DF1

    for (i = 0; ( i < SizeFRA ); i++) // Looking for FRA with a maturity anterior at date
    {
        if ( Date > FRARates.Elt(i, 1)) 
        {
            DF1Date = FRARates.Elt(i, 1);
            DF1        = FRARates.Elt(i, 0);
            check = 1;
        }
    }

    if (check == 0 )  // no FRA before date, take the MM at 3M or interpolate it
    {
        // look for the 3M

        matDate = MMRates.Elt(0, 1);
        matDate.AddMonths(3, GOTO_END_OF_MONTH);
        matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);
        
        i= SizeM-1;
        j = -1;

        while ((MMRates.Elt(i, 1) != matDate.GetJulian()) && (i> -1))
        {
            if (MMRates.Elt(i, 1) > matDate.GetJulian())
            {
                j = i;
            }
            i--;
        }

        if (i > -1) //case 3M exists
        {
            DF1Date = MMRates.Elt(0, 1);
            DF1        = MMRates.Elt(i, 2);    
        }
        else // no 3M
        {
            DF1Date = MMRates.Elt(0, 1);
            
            if ( j != -1 )
            {
                DF1    = Interpol(MMRates.Elt(0, 1), MMRates.Elt(j, 1), 
                                  DF1Date, MMRates.Elt(0, 2), 
                                  MMRates.Elt(j, 2), Cont_Lin, dayCount);
            }
            else
            {
                double DBD360 = DaysBetweenDates(KACTUAL_360, 
                                       (ARM_Date) MMRates.Elt(0, 1), 
                                       (ARM_Date) MMRates.Elt(SizeM-1, 1))/360.0;
    
                double DBD365 = (DaysBetweenDates(KACTUAL_360, 
                                    (ARM_Date) MMRates.Elt(0, 1), 
                                    (ARM_Date) DF1Date)-2)/365.0;
    
                double z = -log(MMRates.Elt(SizeM-1, 0)
                            /MMRates.Elt(0, 0))/DBD360;


                DBD360 = DaysBetweenDates(KACTUAL_360, 
                                       (ARM_Date) MMRates.Elt(0, 1), 
                                       (ARM_Date) DF1Date)/360.0;

                // TO CHECK NOT SURE
                DF1 = MMRates.Elt(0, 0)*exp(-z*DBD360);
                DF1 = 1 / exp (DF1*DBD365);
            }
        }
    }
        
    //DF2

    j = SizeFRA-1;

    while ( FRARates.Elt(j, 1) >= Date ) // first FRA just after Date
    {
        DF2Date = FRARates.Elt(j, 1);
        DF2        = FRARates.Elt(j, 0);
        j--;

        if ( j < 0 )
           break;
    }


    // Interpol
    
        
    if (Cont_Lin == -1) // method basic linear interpolation
    {
        int N1= DaysBetweenDates (dayCount, DF1Date, Date);
        int N2= DaysBetweenDates (dayCount, DF1Date, DF2Date);        
        Forward = (DF2-DF1) * N1 / N2 + DF1;
    }
    else if (DF1Date == Date)
       Forward = DF1;
    else if (DF2Date == Date)
       Forward = DF2;
    else
    {
        Forward = Interpol(DF1Date,
                           DF2Date,
                           Date,
                           DF1,
                           DF2,
                           Cont_Lin,
                           dayCount);
    }

    delete ccyName;

    return(Forward);
}

    

double ARM_ZeroCurve::InterpolDF(double Date,ARM_Matrix& DF, int SizeDF,
                                 ARM_Matrix& MMRates, int SizeM,
                                 int Cont_Lin)
{
    ARM_Currency* ccy = GetCurrencyUnit();
    int dayCount = ccy->GetMMDayCount();
    int j;

    double check = 0, DF1, DF1Date, DF2, DF2Date, Forward;


    // DF1

    int i;
    for (i = 0; ( i < SizeDF ); i++) // Looking for FRA with a maturity anterior at date
    {
        if ( Date >= DF.Elt(i, 1)) 
        {
            DF1Date = DF.Elt(i, 1);
            DF1        = DF.Elt(i, 0);
            check = 1;
        }
    }

    if (check == 0 )  // no FRA before date, take the last MM or before
    {
        i = 1;
        DF1Date = MMRates.Elt(SizeM-1, 1);
        DF1        = MMRates.Elt(SizeM-1, 0);

        while ( Date <= MMRates.Elt(SizeM-i, 1)) 
        {
            DF1Date = MMRates.Elt(SizeM-i-1, 1);
            DF1     = MMRates.Elt(SizeM-i-1, 0);
            i++;
        }
    }
        

    //DF2

    j = SizeDF-1;

    if ((check == 0) && (DF1Date <MMRates.Elt(SizeM-1, 1))) // first MM just after 
    {
        i=SizeM-1;
        
        while ( MMRates.Elt(i, 1) >= Date) 
        {
            DF2Date = MMRates.Elt(SizeM-i, 1);
            DF2        = MMRates.Elt(SizeM-i, 0);
            i--;
        }
    }
    else
    {
        while ( DF.Elt(j, 1) >= Date ) // DF just after Date
        {
            DF2Date = DF.Elt(j, 1);
            DF2        = DF.Elt(j, 0);
            j--;

            if (j<0)
                break;
        }

    }
    
    // Interpol
    if (DF1Date == Date)
      Forward = DF1;
    else if (DF2Date == Date)
      Forward = DF2;
    else
    {
       Forward = Interpol(DF1Date,
                          DF2Date,
                          Date,
                          DF1,
                          DF2,
                          Cont_Lin,
                          dayCount);
    }

    return Forward;

}


int ARM_ZeroCurve::CptTOYFRA(ARM_Matrix& FRARates, int SizeFRA,
                             ARM_Matrix* Futures, int SizeF,
                             int Cont_Lin)
{
    ARM_Date matDate;
    int month1, year1;
    char *year = 0;
    char str1[10];
    int i ,j, k, l;
    double TmP;
    
    ARM_Currency* ccy = GetCurrencyUnit();
    char* ccyName = ccy->GetCcyName();
    int   spotDays = ccy->GetSpotDays();


    ARM_INDEX_TYPE idxType = GetDefaultIndexFromCurrency(ccyName);

    ARM_Currency tmpCurrency(ccyName);

    ccyName = tmpCurrency.GetPayCalName(idxType);

    
    // Procedure to sort rates  

    for (i = 0; i < SizeFRA ; i++)
    {
        for (j = i; j < SizeFRA; j++)
        {  
            if ( FRARates.Elt(i, 1) > FRARates.Elt(j, 1) )
            {
               TmP = FRARates.Elt(j, 1);
               FRARates.Elt(j, 1) = FRARates.Elt(i, 1);
               FRARates.Elt(i, 1) = TmP;

               TmP = FRARates.Elt(j, 2);
               FRARates.Elt(j, 2) = FRARates.Elt(i, 2);
               FRARates.Elt(i, 2) = TmP; 

               TmP = FRARates.Elt(j, 0);
               FRARates.Elt(j, 0) = FRARates.Elt(i, 0);
               FRARates.Elt(i, 0) = TmP;
            }
        }
    }


    // looking for year of 31/12 date
    int y = 0;
    char str2 [10];

    
    ARM_Date lastDEC2 (31, 12, 2000); 

    ARM_Date firstJAN2 (02, 01, 2001); 

    
    // Search the year of the TOY calculate
    while (FRARates.Elt(SizeFRA-1,1) - lastDEC2.GetJulian() > 0)
    {
        lastDEC2.AddYears(1);
        firstJAN2.AddYears(1);
        y++;
    }

    while (lastDEC2.IsBusinessDay(ccyName) == 0)
    {
        lastDEC2.PreviousBusinessDay(spotDays, ccyName);
    }
    
    while (firstJAN2.IsBusinessDay(ccyName) == 0)
    {
        firstJAN2.NextBusinessDay(spotDays, ccyName);
    }


    // determination of SEP's DF
    strcpy (str1,"SEP0");
    ITOA(y, str2,10);
    strcat(str1, str2);
    
    GetMonthYearFromExpiryDate(str1, &month1, &year1);
    matDate.ChangeDate(1, month1, year1);
    matDate.PiborDelivery();
    
    j=1;
    while (Futures->Elt(j, 1) < matDate.GetJulian())
        j++;

    // determination of DEC's DF
    strcpy (str1,"DEC0");
    ITOA(y, str2,10);
    strcat(str1, str2);

    GetMonthYearFromExpiryDate(str1, &month1, &year1);
    matDate.ChangeDate(1, month1, year1);
    matDate.PiborDelivery();

    k=1;
    while (Futures->Elt(k, 1) < matDate.GetJulian())
        k++;

    
    // determination of MAR's DF
    y++;
    strcpy (str1,"MAR0");
    ITOA(y, str2,10);
    strcat (str1, str2);
    GetMonthYearFromExpiryDate(str1, &month1, &year1);
    matDate.ChangeDate(1, month1, year1);
    matDate.PiborDelivery();

    l=1;
    while (Futures->Elt(l, 1) < matDate.GetJulian())
        l++;


    // DEC theorical

    double DECth = (Futures->Elt(j+1, 2)*(Futures->Elt(l, 1)
                   -Futures->Elt(k, 1)) 
                   +Futures->Elt(l+1, 2)*(Futures->Elt(k, 1)
                   -Futures->Elt(j, 1))) 
                   /(Futures->Elt(l, 1) - Futures->Elt(j, 1));

    // z
    double z = log(1+DECth/100.0*DaysBetweenDates(KACTUAL_360, 
                   (ARM_Date) Futures->Elt(k, 1),
                   (ARM_Date) Futures->Elt(l, 1))/360.0)
                /(DaysBetweenDates(KACTUAL_365, 
                  (ARM_Date)Futures->Elt(k, 1), 
                  (ARM_Date)Futures->Elt(l, 1))/365.0); 

    // TOYrate
    double TOYrate = log(1+Futures->Elt(k+1, 2)/100.0
                         *DaysBetweenDates(KACTUAL_360, 
                    (ARM_Date) Futures->Elt(k, 1),
                    (ARM_Date) Futures->Elt(l, 1))/360);

    TOYrate = TOYrate-z*DaysBetweenDates(KACTUAL_365, 
                    (ARM_Date) Futures->Elt(k, 1),
                    (ARM_Date) lastDEC2.GetJulian())/365.0;

    TOYrate = TOYrate-z*DaysBetweenDates(KACTUAL_365, 
                   (ARM_Date) firstJAN2.GetJulian(),
                   (ARM_Date)  Futures->Elt(l, 1))/365.0;

    TOYrate = TOYrate/(DaysBetweenDates(KACTUAL_365, 
              (ARM_Date) lastDEC2.GetJulian(),
               (ARM_Date) firstJAN2.GetJulian())/365.0);

    
    //insert 12/31
    double DBD365 = DaysBetweenDates(KACTUAL_365, 
                                     (ARM_Date) Futures->Elt(k, 1), 
                                     (ARM_Date) lastDEC2.GetJulian() );

    FRARates.Elt(SizeFRA, 0)= Futures->Elt(k, 0)*exp(-z*DBD365/365.0);
    FRARates.Elt(SizeFRA, 1)= lastDEC2.GetJulian();

    SizeFRA++;

    
    //insert 02/01
    DBD365 = DaysBetweenDates(KACTUAL_365, (ARM_Date) firstJAN2.GetJulian(),
                              (ARM_Date) Futures->Elt(l, 1)) ;

    FRARates.Elt(SizeFRA, 0)=Futures->Elt(l, 0)/exp(-z*DBD365/365.0);
    FRARates.Elt(SizeFRA, 1)= firstJAN2.GetJulian();

    SizeFRA++;

    delete ccyName;


    return(SizeFRA);
}



int ARM_ZeroCurve::CptTOY(ARM_Matrix& MMRates, int SizeM,
                          ARM_Matrix* Futures, int SizeF,
                          int Cont_Lin)
{
    ARM_Date matDate;
    int month1, year1;
    char *year = 0;
    char str1[10];
    int j, k, l;

    
    ARM_Currency* ccy = GetCurrencyUnit();
    char* ccyName = ccy->GetCcyName();
    int   spotDays = ccy->GetSpotDays();

    ARM_INDEX_TYPE idxType = GetDefaultIndexFromCurrency(ccyName);

    ARM_Currency tmpCurrency(ccyName);

    ccyName = tmpCurrency.GetPayCalName(idxType);

    // FIRST TOY

    // looking for year of 31/12 date
    int y = 0;
    char str2 [10];

    
    ARM_Date lastDEC (31, 12, 2000); 

    ARM_Date firstJAN (02, 01, 2001); 

    
    // find the current year

    while (Futures->Elt(1,1) - lastDEC.GetJulian() > 0)
    {
        lastDEC.AddYears(1);
        firstJAN.AddYears(1);
        y++;
    }

    while (lastDEC.IsBusinessDay(ccyName) == 0)
    {
        lastDEC.PreviousBusinessDay(spotDays, ccyName);
    }
    
    while (firstJAN.IsBusinessDay(ccyName) == 0)
    {
        firstJAN.NextBusinessDay(spotDays, ccyName);
    }


    // determination of SEP's DF
    strcpy (str1,"SEP0");
    ITOA(y, str2,10);
    strcat ( str1, str2);
    
    GetMonthYearFromExpiryDate(str1, &month1, &year1);
    matDate.ChangeDate(1, month1, year1);
    matDate.PiborDelivery();
    
    j=1;
    while (Futures->Elt(j, 1) < matDate.GetJulian())
        j++;

    // determination of DEC's DF
    strcpy (str1,"DEC0");
    ITOA(y, str2,10);
    strcat(str1, str2);
    GetMonthYearFromExpiryDate(str1, &month1, &year1);
    matDate.ChangeDate(1, month1, year1);
    matDate.PiborDelivery();

    k=1;
    while (Futures->Elt(k, 1) < matDate.GetJulian())
        k++;

    
    // determination of MAR's DF
    y++;
    strcpy (str1,"MAR0");
    ITOA(y, str2,10);
    strcat (str1, str2);
    GetMonthYearFromExpiryDate(str1, &month1, &year1);
    matDate.ChangeDate(1, month1, year1);
    matDate.PiborDelivery();

    l = 1;
    while (Futures->Elt(l, 1) < matDate.GetJulian())
        l++;


    // DEC theorical

    double DECth = (Futures->Elt(j+1, 2)*(Futures->Elt(l, 1)
                   -Futures->Elt(k, 1)) 
            +Futures->Elt(l+1, 2)*(Futures->Elt(k, 1)-Futures->Elt(j, 1))) 
            /(Futures->Elt(l, 1)-Futures->Elt(j, 1));

    // z
    double z = log(1+DECth/100.0*DaysBetweenDates(KACTUAL_360, 
               (ARM_Date) Futures->Elt(k, 1),
               (ARM_Date) Futures->Elt(l, 1))/360.0)
               /(DaysBetweenDates(KACTUAL_365, 
               (ARM_Date)Futures->Elt(k, 1), 
               (ARM_Date)Futures->Elt(l, 1))/365.0);

    // TOYrate
    double TOYrate = log(1+Futures->Elt(k+1, 2)/100.0
                         *DaysBetweenDates(KACTUAL_360, 
                          (ARM_Date) Futures->Elt(k, 1),
                          (ARM_Date)  Futures->Elt(l, 1))/360.0);

    TOYrate = TOYrate-z*DaysBetweenDates(KACTUAL_365, 
                (ARM_Date) Futures->Elt(k, 1),
                (ARM_Date)  lastDEC.GetJulian())/365.0;

    TOYrate = TOYrate-z*DaysBetweenDates(KACTUAL_365, 
                    (ARM_Date) firstJAN.GetJulian(),
                    (ARM_Date)  Futures->Elt(l, 1))/365.0;

    TOYrate = TOYrate/(DaysBetweenDates(KACTUAL_365, 
              (ARM_Date) lastDEC.GetJulian(),
              (ARM_Date) firstJAN.GetJulian())/365.0);

    //insert 12/31
    double DBD365 = DaysBetweenDates(KACTUAL_365, 
                            (ARM_Date) Futures->Elt(k, 1), 
                            (ARM_Date) lastDEC.GetJulian());

    Futures->Elt(SizeF, 0)= Futures->Elt(k, 0)*exp(-z*DBD365/365);
    Futures->Elt(SizeF, 1)= lastDEC.GetJulian();

    SizeF++;

    
    //insert 02/01
    DBD365 = DaysBetweenDates(KACTUAL_365, (ARM_Date) firstJAN.GetJulian(), 
                              (ARM_Date) Futures->Elt(l, 1)) ;

    Futures->Elt(SizeF, 0)=Futures->Elt(l, 0)/exp(-z*DBD365/365.0);
    Futures->Elt(SizeF, 1)= firstJAN.GetJulian();

    SizeF++;


    // SECOND TOY

    ARM_Date lastDEC2 (31, 12, 2000); 

    ARM_Date firstJAN2 (02, 01, 2001); 

    lastDEC2.AddYears(y);
    firstJAN2.AddYears(y);
    

    while (lastDEC2.IsBusinessDay(ccyName) == 0)
    {
        lastDEC2.PreviousBusinessDay(spotDays, ccyName);
    }
    
    while (firstJAN2.IsBusinessDay(ccyName) == 0)
    {
        firstJAN2.NextBusinessDay(spotDays, ccyName);
    }


    // determination of SEP's DF
    strcpy (str1,"SEP0");
    ITOA(y, str2,10);
    strcat ( str1, str2);
    
    GetMonthYearFromExpiryDate(str1, &month1, &year1);
    matDate.ChangeDate(1, month1, year1);
    matDate.PiborDelivery();
    
    j=1;
    while (Futures->Elt(j, 1) < matDate.GetJulian())
        j++;

    // determination of DEC's DF
    strcpy (str1,"DEC0");
    ITOA(y, str2,10);
    strcat(str1, str2);
    GetMonthYearFromExpiryDate(str1, &month1, &year1);
    matDate.ChangeDate(1, month1, year1);
    matDate.PiborDelivery();

    k=1;
    while (Futures->Elt(k, 1) < matDate.GetJulian())
        k++;

    
    // determination of MAR's DF
    y++;
    strcpy (str1,"MAR0");
    ITOA(y, str2,10);
    strcat(str1, str2);
    GetMonthYearFromExpiryDate(str1, &month1, &year1);
    matDate.ChangeDate(1, month1, year1);
    matDate.PiborDelivery();

    l=1;
    while (Futures->Elt(l, 1) < matDate.GetJulian())
        l++;


    // DEC theorical

    DECth = (Futures->Elt(j+1, 2)*(Futures->Elt(l, 1)-Futures->Elt(k, 1)) 
            +Futures->Elt(l+1, 2)*(Futures->Elt(k, 1)-Futures->Elt(j, 1))) 
            /(Futures->Elt(l, 1)-Futures->Elt(j, 1));

    // z
    z = log(1+DECth/100.0*DaysBetweenDates(KACTUAL_360, 
           (ARM_Date) Futures->Elt(k, 1),(ARM_Date) Futures->Elt(l, 1))/360.0)
            /(DaysBetweenDates(KACTUAL_365, (ARM_Date)Futures->Elt(k, 1), 
            (ARM_Date)Futures->Elt(l, 1))/365.0);

    // TOYrate
    TOYrate = log(1+Futures->Elt(k+1, 2)/100.0*
                  DaysBetweenDates(KACTUAL_360, 
               (ARM_Date) Futures->Elt(k, 1),
               (ARM_Date) Futures->Elt(l, 1))/360);

    TOYrate = TOYrate-z*DaysBetweenDates(KACTUAL_365, 
              (ARM_Date) Futures->Elt(k, 1),
                (ARM_Date) lastDEC2.GetJulian())/365.0;

    TOYrate = TOYrate-z*DaysBetweenDates(KACTUAL_365, 
             (ARM_Date) firstJAN2.GetJulian(), 
                (ARM_Date) Futures->Elt(l, 1))/365.0;

    TOYrate = TOYrate/(DaysBetweenDates(KACTUAL_365, 
               (ARM_Date) lastDEC2.GetJulian(),
               (ARM_Date) firstJAN2.GetJulian())/365.0);

    //insert 12/31
    DBD365 = DaysBetweenDates(KACTUAL_365, (ARM_Date) Futures->Elt(k, 1), 
                              (ARM_Date) lastDEC2.GetJulian());

    Futures->Elt(SizeF, 0)= Futures->Elt(k, 0)*exp(-z*DBD365/365);
    Futures->Elt(SizeF, 1)= lastDEC2.GetJulian();

    SizeF++;

    
    // insert 02/01
    DBD365 = DaysBetweenDates(KACTUAL_365, (ARM_Date) firstJAN2.GetJulian(), 
                              (ARM_Date) Futures->Elt(l, 1));
    Futures->Elt(SizeF, 0)=Futures->Elt(l, 0)/exp(-z*DBD365/365);
    Futures->Elt(SizeF, 1)= firstJAN2.GetJulian();

    SizeF++;

    delete ccyName;

    return(SizeF);
}



void ARM_ZeroCurve::ViewForward(char* id, FILE* fOut)
{
    char strDate[20];

 
    fprintf(fOut, "\n\n\t =====> 3 Months Forward Rates \n\n");

    itsAsOfDate.JulianToStrDate(strDate);
    fprintf(fOut, "\n\n\t AsOfDate : %s \n\n\n", strDate);

 
    int sz = itsForwardDates->GetSize();
    char d[20];
    int  i;
 
    fprintf(fOut, "\n\nDate\t\tDays\tForwardRate\t\t\n\n");
   
    //itsForwardDates->Elt(0) = number of Forward Rates
    for (i = 1; i < itsForwardDates->Elt(0); i++)
    {
        ARM_Date date, tmpDate;
        double   rate, dF;
        int      days;

//        itsYearTerms->Elt(i) = itsForwardDates->Elt(i)/365.0;

        if ( (*itsYearTerms)[i] == 1000.0 )
        {
           break;
        }

        days = int(NEAREST_INT(itsForwardDates->Elt(i)-itsAsOfDate.GetJulian()));

        tmpDate = itsAsOfDate;

        date = tmpDate.AddDays(days);

        date.JulianToStrDate(d);
 
        rate = (*itsForwardRates)[i];

        if (itsDiscountFactors)
           dF = (*itsForwardDF)[i];
        else
           dF = 0.0;

        if ( (*itsYearTerms)[i] != 1000.0 )
        {
           fprintf(fOut, "%s\t%5d\t %.5lf \t \n", d, days, rate);
        }
    }
}



void ARM_ZeroCurve::CptFuturesZeroRatesTOY(ARM_Matrix& MMRates, int SizeM,
                                           ARM_Matrix& FRARates, int SizeFRA,
                                           ARM_Matrix* Future, int SizeF,
                                           ARM_Matrix Swap,
                                           int MMVsFut, int Cont_Lin)
{
    ARM_Date D2_Cash;
    ARM_Date Ctr_Ech1, Ctr_Ech2;
    int i, j;
    double N6;
    double Df_Ctr, TmP;
    double DF1, DF1Date, DF2, DF2Date;


    char* ccyName = GetCurrencyUnit()->GetCcyName();
    
    // Procedure de tri des taux Futures

    for (i = 1; i < SizeF; i++)
    {
        for (j = i; j < SizeF ; j++)
        {  
            if ( Future->Elt(i, 1) > Future->Elt(j, 1) )
            {
               TmP = Future->Elt(j, 1);
               Future->Elt(j, 1) = Future->Elt(i, 1);
               Future->Elt(i, 1) = TmP;
                
               TmP = Future->Elt(j,0);
               Future->Elt(j, 0) = Future->Elt(i, 0);
               Future->Elt(i, 0) = TmP;
            }
        }
    }

    // Calcul des Df a partir des taux Futures

    for (i = 1; i < SizeF-1; i++)
    {
        if ( fabs(Future->Elt(i+1, 1) - Future->Elt(i, 1) - 90.0) > 10.0 )
        {
           char BUF[100];


           sprintf(BUF, "%s",
               "One or more future contract absent between two contracts ");

           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,BUF);
        }
    }


    // DF1
    
    for (i = 0; ( i < SizeFRA ); i++) // Looking for FRA with a maturity anterior at date
    {
        if ( Future->Elt(1, 1) > FRARates.Elt(i, 1)) 
        {
            DF1Date = FRARates.Elt(i, 1);
            DF1        = FRARates.Elt(i, 0);
        }
    }


    //DF2

    j = SizeFRA-1;

    while ( FRARates.Elt(j, 1) >= Future->Elt(1, 1) ) // first FRA just after Date
    {
        DF2Date = FRARates.Elt(j, 1);
        DF2        = FRARates.Elt(j, 0);
        j--;
        if (j<0)
           break;
    }
  
    // MA : 23/05/2002
    int inputDayCount = KACTUAL_360;

    if ( (strcmp(ccyName, "GBP") == 0) || (strcmp(ccyName, "CAD") == 0) )
    {
       inputDayCount = KACTUAL_365;
    }

    double tmp = Future->Elt(1, 1);


    Df_Ctr=Interpol(DF1Date,
                    DF2Date,
                    Future->Elt(1, 1),
                    DF1,
                    DF2,
                    Cont_Lin,
                    inputDayCount);
    
    Future->Elt(0, 0) = Df_Ctr;

    Future->Elt(0, 1) = Future->Elt(1, 1);
    
    Future->Elt(0, 2) = Future->Elt(1, 2);
    
    
    i = 1;

    tmp = Future->Elt(i, 0);
    Ctr_Ech1 = Future->Elt(1, 1);

    while ( i < SizeF )
    {
        Ctr_Ech2 = Ctr_Ech1;
        Ctr_Ech2.PiborDelivery();
    
        N6 = DaysBetweenDates(inputDayCount, Ctr_Ech1, Ctr_Ech2);

        if ( (strcmp(ccyName, "GBP") == 0) || (strcmp(ccyName, "CAD") == 0) )
        {
           Df_Ctr = Df_Ctr / ( 
                            ( 1 + (Future->Elt(i, 0)) / 100.0 
                            * ( N6  / 365.0 ))
                         );
        }
        else
        {
           Df_Ctr = Df_Ctr / (
                            ( 1 + (Future->Elt(i, 0)) / 100.0
                            * ( N6  / 360.0 ))
                         );
        }

        Future->Elt(i, 0) = Df_Ctr;
        Future->Elt(i, 1) = Ctr_Ech2.GetJulian();
        Ctr_Ech1 = Future->Elt(i, 1) ;

        i++;
    }
}



void ARM_ZeroCurve::CptSwapZeroRatesTOY(ARM_Matrix& MMRates, int SizeM,
                                        ARM_Matrix& FRARates, int SizeFRA,
                                        ARM_Matrix* Swap, 
                                        int& SizeS, int fxPayFreq,
                                        ARM_Matrix& Future, int SizeF, 
                                        int MMVsFut,
                                        int SwapVsFut,
                                        int raw,
                                        int Cont_Lin)
{
    int i, j, k, iter=0, cpt=0;
    double TmP, inter, tmp2, TmP2, startDF;

    ARM_Date startDate = itsAsOfDate;

                   
    char* ccyName = GetCurrencyUnit()->GetCcyName();
    int spotDays = GetCurrencyUnit()->GetSpotDays();
    int fxDayCount = GetCurrencyUnit()->GetFixedDayCount();

    ARM_INDEX_TYPE idxType = GetDefaultIndexFromCurrency(ccyName);

    ARM_Currency tmpCurrency(ccyName);

    ccyName = tmpCurrency.GetPayCalName(idxType);

    
    while ( Swap->Elt(cpt,2) == 0 )
    {
        cpt++;
    }


    double tempo = Swap->Elt(1, 0);

    if ( spotDays == 0 ) 
       startDF = 1.0;
    else
    {
       startDate.NextBusinessDay(spotDays, ccyName);

       startDF = Interpol(MMRates.Elt(0, 1),
                          MMRates.Elt(1, 1),
                          startDate.GetJulian(),
                          MMRates.Elt(0, 0),
                          MMRates.Elt(1, 0),
                          Cont_Lin,
                          KACTUAL_360);
    }
    
    ARM_Matrix Df_interpol(200, 2, 0.0);
    double Df_Ctr;

    // Compute DF for maturities < 1 year when fixed leg payment frequency 
    // is > 1 
    
    if ( fxPayFreq > 1 )
    {
        if ( MMVsFut != K_FUT )
        {
           for (i=1; i<fxPayFreq; i++)
           {
               j = 0;

               while (Swap->Elt(i, 1) > MMRates.Elt(j, 1)) 
                     j++;

               Swap->Elt(i, 0) = Interpol(MMRates.Elt(j-1, 1),
                                          MMRates.Elt(j, 1),
                                          Swap->Elt(i, 1),
                                          MMRates.Elt(j-1, 0),
                                          MMRates.Elt(j, 0),
                                          Cont_Lin,
                                          KACTUAL_360);

               Df_interpol.Elt(i, 0) = Swap->Elt(i, 0);
           }
        }
        else 
        {
           for (i=1; i< cpt; i++) 
           {
               j = 0;

           /*    while (Swap->Elt(i, 1) > Future.Elt(j, 1)) 
                     j++;

               Swap->Elt(i, 0) = Interpol(Future.Elt(j-1, 1),
                                    Future.Elt(j, 1),
                                    Swap->Elt(i, 1),
                                   Future.Elt(j-1, 0),
                                    Future.Elt(j, 0),
                                    Cont_Lin,
                                    KACTUAL_360);
             */
               if (Swap->Elt(i, 1) < FRARates.Elt(SizeFRA-1, 1))
               {
                  while (Swap->Elt(i, 1) > FRARates.Elt(j, 1)) 
                     j++;

                  Swap->Elt(i, 0) = Interpol(FRARates.Elt(j-1, 1),
                                             FRARates.Elt(j, 1),
                                             Swap->Elt(i, 1),
                                             FRARates.Elt(j-1, 0),
                                             FRARates.Elt(j, 0),
                                             Cont_Lin,
                                             KACTUAL_360);

               }
               else
               {
                  while (Swap->Elt(i, 1) > Future.Elt(j, 1)) 
                     j++;

                  Swap->Elt(i, 0) = Interpol(Future.Elt(j-1, 1),
                                             Future.Elt(j, 1),
                                             Swap->Elt(i, 1),
                                             Future.Elt(j-1, 0),
                                             Future.Elt(j, 0),
                                             Cont_Lin,
                                             KACTUAL_360);
               }

               Df_interpol.Elt(i, 0) = Swap->Elt(i, 0);
            }
        }
    }

    if ( SwapVsFut == K_FUT )
    {
        for (i=fxPayFreq; i<=SizeS ; i++)
        {
            for (j=1; j<SizeF; j++)
            {
                // si FUT, on interpole les DF aux maturites 
                // des swaps < maturite du dernier future
                // en utilisant les DF calcules a partir des futures

                if ( Future.Elt(j, 1) > Swap->Elt(i, 1) )
                {
                    Df_Ctr = Interpol(Future.Elt(j-1, 1),
                                      Future.Elt(j, 1),
                                      Swap->Elt(i, 1),
                                      Future.Elt(j-1, 0),
                                      Future.Elt(j, 0),
                                      Cont_Lin,
                                      KACTUAL_360);
                    
                    Df_interpol.Elt(i, 0) = Df_Ctr;

                    // calcul des taux de swap 
                    TmP = 0.0;
                    
                    for (k=1; k<=i; k++)
                    {
                        inter = CountYears(fxDayCount,
                        ( k==1 ? startDate :
                          (ARM_Date) Swap->Elt(k-1, 1) ),
                                           (ARM_Date) Swap->Elt(k, 1));
    
                        TmP += Df_interpol.Elt(k, 0)*inter;
                    }

                    Swap->Elt(i, 0) = 100.0*
                                     (startDF-Df_interpol.Elt(i, 0))/TmP;

                    break;
                }
            }
        }
    }

    for (i=1; i<=SizeS; i++)
    { 
        Df_interpol.Elt(i, 0) = 0.0;
    }

    // si methode PAR, on calcule les taux de swap manquant 
    // par interpolation lineaire
    // des taux de swap renseignes 


    if ( raw == K_PAR )
    {
       for (i = cpt; i <= SizeS; i++)
       {
           if ( Swap->Elt(i, 0) < 0.0001)
           {
              j = i+1;

              while ( Swap->Elt(j, 0) < 0.0001 ) 
                    j++;
            
              double temp = Swap->Elt(i, 1);

              temp = Swap->Elt(i-1, 0);
              temp = Swap->Elt(j, 0);
              temp = Swap->Elt(i, 0);

              Swap->Elt(i, 0) = Swap->Elt(i-1, 0) 
                                +(Swap->Elt(j, 0)-Swap->Elt(i-1, 0))
                         *DaysBetweenDates(fxDayCount, 
                                        (ARM_Date) Swap->Elt(i-1, 1),
                                           (ARM_Date) Swap->Elt(i, 1))
                 /DaysBetweenDates(fxDayCount, (ARM_Date) Swap->Elt(i-1, 1),
                                              (ARM_Date) Swap->Elt(j, 1));
           }
       }
    }
    
    
    for (i = cpt; i <= SizeS ; i++)
    {
        // calcul implicite de chaque taux non renseigne

        if ( (Swap->Elt(i, 0) < 0.0001) && (raw == K_RAW) )
        {
           // chercher les taux encadrant le taux a calculer 

           j = i+1;

           while ( Swap->Elt(j, 0) < 0.0001 ) 
                 j++;
            
           Df_Ctr = Swap->Elt(i-1, 0);
            
           do
           {
                for (k = i; k < j; k++)
                {
                    Swap->Elt(k, 0) = Interpol(Swap->Elt(i-1, 1),
                                            Swap->Elt(j, 1),
                                            Swap->Elt(k, 1),
                                            Swap->Elt(i-1, 0),
                                            Df_Ctr,
                                            Cont_Lin,
                                            KACTUAL_360);
                    
                    Df_interpol.Elt(k, 0) = (-Swap->Elt(k, 0)+
                                            Interpol(Swap->Elt(i-1, 1),
                                                    Swap->Elt(j, 1),
                                                    Swap->Elt(k, 1),
                                                    Swap->Elt(i-1, 0),
                                                    Df_Ctr + 1.0e-10,
                                                    Cont_Lin,
                                                    KACTUAL_360))/1.0e-10;
                    
                    Df_interpol.Elt(k, 1) = 1.0;
                }


                TmP = 0.0;
                TmP2 = 0.0;

                for (k = 1; k<j ;k++)
                {
                    inter = CountYears(fxDayCount,
                           ( k==1 ? startDate :
                                (ARM_Date) Swap->Elt(k-1, 1) ),
                                          (ARM_Date) Swap->Elt(k, 1) );
                    
                    TmP +=  Swap->Elt(k, 0) * inter;
                    TmP2 +=  Df_interpol.Elt(k, 0) * inter;
                }

                TmP *= Swap->Elt(j, 0)/100.0;
                TmP2 *= Swap->Elt(j, 0)/100.0;

                TmP = startDF-TmP;
                TmP2 = -TmP2;

                inter = CountYears(fxDayCount,
                          ( j==1 ? startDate :
                                   (ARM_Date) Swap->Elt(j-1, 1) ),
                                              (ARM_Date) Swap->Elt(j, 1));

                tmp2 = 1.0 + Swap->Elt(j, 0)*inter/100.0;

                //printf ("Df_Ctr1 :%10.10f ",Df_Ctr);

                Df_Ctr= Df_Ctr - (TmP - tmp2*Df_Ctr)/(TmP2 - tmp2);
            }
            while ( fabs((TmP - tmp2*Df_Ctr)/(TmP2 - tmp2)) > 1.0e-6 );
        
            Swap->Elt(j,0) = Df_Ctr;

            i = j;
        }
        else 
        {
            // si meth = PAR, calculer les DF a partir 
            // des taux de swap interpoles

            TmP = 0.0;
            
            for (j=1; j<i; j++)
            {
                inter = CountYears(fxDayCount, 
                       ( j==1 ? startDate :
                            (ARM_Date) Swap->Elt(j-1, 1) ),
                                      (ARM_Date) Swap->Elt(j, 1) );

                TmP += Swap->Elt(j, 0)*inter;
            }
        
            TmP *= Swap->Elt(i, 0)/100.0;

            TmP = startDF-TmP;

            inter = CountYears(fxDayCount, 
                      ( i==1 ? startDate :
                               (ARM_Date) Swap->Elt(i-1, 1) ),
                                          (ARM_Date) Swap->Elt(i, 1));
            ARM_Date d1, d2;

            if ( i > 1 )
            {
               d1 = (ARM_Date) Swap->Elt(i-1, 1);
               d2 = (ARM_Date) Swap->Elt(i, 1);
            }

            Swap->Elt(i, 0) = 1.0 + Swap->Elt(i, 0)*inter/100.0;

            Swap->Elt(i, 0) = TmP/Swap->Elt(i, 0);

        }
    }

    // si meth = RAW, virer les taux manquant initialement 

    if ( raw == K_RAW )
    {
       iter = 0;

       for (i = 1; i <= SizeS; i++)
       {
           if ( Df_interpol.Elt(i, 1) > 0.0001 )
           {
              for (j = i; j < SizeS; j++)
              {
                  Df_interpol.Elt(j, 1) = Df_interpol.Elt(j+1, 1);

                  Swap->Elt(j, 0) = Swap->Elt(j+1, 0);
                  Swap->Elt(j, 1) = Swap->Elt(j+1, 1);
              }

              SizeS--;

              i--;
           }
       }
    }
}


//--------------------------------------------------
// variable compteur: pour indiquer le prochain point correspondant à freq à établir
// variable indexDle: pour indiquer la position de actual date
// varialbe indexInt: le point correspondant à freq pour indexDle
// eg. si freq=1, actual donnée: 1.5Y; donc indexDle=1.5; indexInt = 1;
void ConstructDynamSwap(int NbMonthYear, int flagMonthOrYear, double dataIn,
						int flagAUDTri, //flagAUDTri==1 si matu==matu2 eg.1YY
						const ARM_Date& matDate,
						int lastMMTerm,
						const ARM_Matrix& Futures, int SizeF,
						ARM_Currency* ccy,						
						const ARM_Date& swapStartDate,
						ARM_Matrix* SwapRates,
						int fxPayFreq, int& SizeS, 
						int& compteurSwap,
						int& firstSwapIndex,
						int& lastAUDtriNb)	//lastAUDtriNb: last nb of years in the part of the AUD trimestre
						
{
	ARM_Date indexSwapDate;

    char ccyName[10];
	ccy->CalcFloatPayCal(ccyName);

//	char* ccyName = ccy->GetCcyName();
//    int fxPayFreq = ccy->GetFixedPayFreq();
	int fxPayFreqVal=fxPayFreq ;

	if (flagMonthOrYear==0) // 18M ...(matu =='M' && Nb>12)
	{
		double valTmp;
			
		valTmp = ((double) NbMonthYear)/((double) 12.0); //Here NbMonthYear is nb of monthes
		
		int indexInt = (int) (floor(valTmp*fxPayFreq));

		double indexDle = (valTmp*fxPayFreq);

		if ( compteurSwap > indexInt )  // compteurSwap > indexInt; on va établir des points avant compteur
		{
			if(matDate.GetJulian() > Futures.Elt(SizeF-1,1))
			{						
				if ( double(indexInt) == indexDle )
				{
					SwapRates->Elt(SizeS-1,0) = dataIn;
				}  //On peux l'enlever
				else
				{
					SwapRates->Elt(SizeS,1) = matDate.GetJulian();
					SwapRates->Elt(SizeS,0) = dataIn;
					SizeS++;
				}
			}
		}
		else  // compteurSwap <= indexInt; on va établir le point a la position de compteurSwap et des point
		      // de frequence jusqu'a actual position, si actual position n'est pas a frequence, on l'établit.
		{
			while(compteurSwap < indexInt)
			{
				indexSwapDate = swapStartDate;

				indexSwapDate.AddMonths(compteurSwap*12/fxPayFreq); 
 						
				indexSwapDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

				SwapRates->Elt(SizeS, 1) = indexSwapDate.GetJulian();

				compteurSwap++;

				SizeS++;
			}
			if ( (double)indexInt == indexDle )
			{
				SwapRates->Elt(SizeS, 1) = matDate.GetJulian();
				SwapRates->Elt(SizeS, 0) = dataIn;

				firstSwapIndex = MIN(SizeS,firstSwapIndex);

				compteurSwap++;
				SizeS++;
			}
			else
			{
				indexSwapDate = swapStartDate;
				indexSwapDate.AddMonths(compteurSwap*12/fxPayFreq); 
				indexSwapDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

				SwapRates->Elt(SizeS, 1) = indexSwapDate.GetJulian();

				compteurSwap++;
				SizeS++;
							
				firstSwapIndex = MIN(SizeS,firstSwapIndex);
						
				if(matDate.GetJulian() > Futures.Elt(SizeF-1,1))
				{
					SwapRates->Elt(SizeS, 1) = matDate.GetJulian();
					SwapRates->Elt(SizeS, 0) = dataIn;
					SizeS++;
				}
			}
		}
	}			
	else
    {
		int value;
		int addnbmois;
				
        //if ( lastMMTerm < 12 )
        {
  			if (strcmp(ccyName, "AUD")==0 && flagAUDTri==1 )// do with the currency "AUD"; two frequences,eg.1YY
			{
				fxPayFreqVal = 4;

				value = NbMonthYear*fxPayFreqVal;
			}
			else
			{
				fxPayFreqVal = fxPayFreq;

				value = lastAUDtriNb*4+(NbMonthYear-lastAUDtriNb)*fxPayFreqVal;	  
			}
				  
			while( compteurSwap< value )
			{
				if(strcmp(ccyName, "AUD")==0 && flagAUDTri==1)
				{
					addnbmois = compteurSwap*12/fxPayFreqVal;
							  
					lastAUDtriNb = NbMonthYear;
				}
				else
					addnbmois = (compteurSwap-lastAUDtriNb*4)*12/fxPayFreqVal+12*lastAUDtriNb;

				indexSwapDate = swapStartDate;
				indexSwapDate.AddMonths(addnbmois); 
				indexSwapDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

				SwapRates->Elt(SizeS, 1) = indexSwapDate.GetJulian();
					  
				SizeS++;
				compteurSwap++;
			}
			if(value < compteurSwap)
			{
				SwapRates->Elt(SizeS+NbMonthYear*fxPayFreqVal-compteurSwap, 0) = dataIn;
            }
			else
			{
				SwapRates->Elt(SizeS, 0) = dataIn;
            
				SwapRates->Elt(SizeS, 1) = matDate.GetJulian();

	 			firstSwapIndex = MIN(SizeS,firstSwapIndex);

				compteurSwap++;

				SizeS++;
			}

		}
	}
}
      
//---------------------------------------------------

//FlowAll:    //all the payment flows
//FlowFreq:   //variable to marque the index of the FlowAll when payment flows frequency =  payfreq
//FlowNoFreq: //variable to marque the index of the FlowAll when payment flows frequency is more than payfreq
			  //and the associated index FlowFreq
void PaymentFlowVector(ARM_Matrix* FlowAll, int SizeS,
					   int fxPayFreq,
					   const ARM_Date& startDate,
					   char* ccyName,
					   int indexAUD,
					   ARM_Vector& FlowFreq, int& SizeFlowFreq,     // index à partir 1
					   ARM_Vector& FlowNoFreq, int& SizeNoFlowFreq, // index à partir 1 
					   int& SizeAUDtri)
{
	ARM_Date indexSwapDate = startDate;

	indexSwapDate.AddMonths(12/fxPayFreq);
	indexSwapDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

	int i;
	for(i=1; i<=SizeS; i++)
	{
		if (FlowAll->Elt(i,1)==indexSwapDate.GetJulian())
		{
			SizeFlowFreq++;
			FlowFreq.Elt(SizeFlowFreq) = i;

			indexSwapDate = startDate;
			indexSwapDate.AddMonths((SizeFlowFreq+1)*12/fxPayFreq);
			indexSwapDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);
		}
		else
		{			
				SizeNoFlowFreq++;
			    
				FlowNoFreq.Elt(SizeNoFlowFreq) = i; //nb
			    //FlowNoFreq.Elt(SizeNoFlowFreq,1) = SizeFlowFreq;//index associe
				if(i<=indexAUD)
				{
					SizeAUDtri++;
				}
		}
	}
}



double ARM_ZeroCurve::InterpDfSwapMilieu(	double dateJul, 
											const ARM_Matrix& MMRates, int SizeM, 
											const ARM_Matrix& Future, int SizeF, 
											const ARM_Matrix& Swap, int SizeS,
											int Cont_Lin)
{
	int i=1;
	double Df_Ctr;
	int flag =0;
	
	if(dateJul<MMRates.Elt(0,1) || dateJul>Swap.Elt(SizeS,1))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "date is out of the range");

	while( i<SizeM )
	{
		if( dateJul>=MMRates.Elt(i-1,1) && dateJul<=MMRates.Elt(i,1) )	
		{
			flag = 1;		
			break;
		}
		else 
			i++;
	}
	
	if(flag==1) //dateIn est dans les MMrates
	{
		Df_Ctr = Interpol(MMRates.Elt(i-1,1),
						  MMRates.Elt(i,1),
						  dateJul,
						  MMRates.Elt(i-1,0),
						  MMRates.Elt(i,0),
						  Cont_Lin,
						  KACTUAL_360);
		return Df_Ctr;
	}		
	
	if (dateJul<Future.Elt(1,1))
		return Df_Ctr = Interpol(MMRates.Elt(SizeM-1,1),
								Future.Elt(1,1),
								dateJul,
								MMRates.Elt(SizeM-1,0),
								Future.Elt(1,0),
								Cont_Lin,
								KACTUAL_360);
	i=2;
	if(SizeF!=1)
	{
		while(i<=SizeF)
		{
			if( dateJul>=Future.Elt(i-1,1) && dateJul<Future.Elt(i,1) )	
			{
				flag = 1;		
				break;
			}
			else 
				i++;
		}
		
		if(flag==1) //dateIn est dans les Future
		{
			Df_Ctr = Interpol(Future.Elt(i-1,1),
							  Future.Elt(i,1),
							  dateJul,
							  Future.Elt(i-1,0),
							  Future.Elt(i,0),
							  Cont_Lin,
							  KACTUAL_360);
			return Df_Ctr;
		}
	}

	i=2;
	while(i<=SizeS)
	{
		if( dateJul>Swap.Elt(i-1,1) && dateJul<=Swap.Elt(i,1) )	
		{
			flag = 1;		
			break;
		}
		else 
			i++;
	}
	
	if(flag==1) //dateIn est dans les Swap
	{
		if(dateJul == Swap.Elt(i,1) ) 
			Df_Ctr = Interpol(Swap.Elt(i-1,1),
							  Swap.Elt(i+1,1),
							  dateJul,
							  Swap.Elt(i-1,0),
							  Swap.Elt(i+1,0),
							  Cont_Lin,
							  KACTUAL_360);
		else
			Df_Ctr = Interpol(Swap.Elt(i-1,1),
							  Swap.Elt(i,1),
							  dateJul,
							  Swap.Elt(i-1,0),
							  Swap.Elt(i,0),
							  Cont_Lin,
							  KACTUAL_360);
		return Df_Ctr;
	}

	return MINDOUBLE;
}



double ARM_ZeroCurve::CalculSwapFlowNoFreq(double dateJul, double txSwap, double startDF,
							const ARM_Matrix& MMRates, int SizeM, 
							const ARM_Matrix& Future, int SizeF, 
							const ARM_Matrix& Swap, int SizeS,
							int fxPayFreq,
							int Cont_Lin,
							int fxDayCount,
							int spotdays,
							char* ccyName)
{
	ARM_Date date(dateJul);
	
	ARM_Date startDate = itsAsOfDate;

	startDate.NextBusinessDay(spotdays,ccyName);

	double TmP=0.0, inter, result;
	
	ARM_Vector flowNoFreq(200,0.0);
	
	// construct the flows by constructing the startdates
	ARM_Vector* starts = CptStartDates(startDate, date,
									   fxPayFreq,
									   K_PREVIOUS, K_SHORTSTART, K_ADJUSTED,
									   ccyName);
	int sz = starts->GetSize();

	for(int j=0; j<sz; j++)
	{
		flowNoFreq.Elt(j) = starts->Elt(j);		
	}

	flowNoFreq.Elt(sz) = date.GetJulian();

	sz = sz+1;

	for(int i=0; i<sz-2; i++)
	{
		inter = CountYears(fxDayCount, (ARM_Date) flowNoFreq.Elt(i),
		                 (ARM_Date) flowNoFreq.Elt(i+1));
		TmP+=InterpDfSwapMilieu(flowNoFreq.Elt(i+1),
								MMRates, SizeM,
								Future, SizeF,
								Swap, SizeS,
								Cont_Lin)*inter; 
	}

	TmP *= txSwap/100.00;
	
	TmP = startDF - TmP;
	
	inter = CountYears(fxDayCount, (ARM_Date) flowNoFreq.Elt(sz-2),
	                   (ARM_Date) flowNoFreq.Elt(sz-1));

	result = 1.0 + txSwap*inter/100.0;
	
	result  = TmP/result;

	delete starts;
	return result;
}


		/*******************************************/
		/*                                         */
		/* Fonctions utilisées par Méthode Forward */
		/*                                         */
		/*       Andria Rajaona - 27/08/04         */
		/*                                         */
		/*******************************************/


void ARM_ZeroCurve::FillFutSwapDisc(ARM_Matrix* LocDiscData, int& LocDiscSize,
									  double FixRate, double FixDaycount, double FixedFreq,
									  double StartFwd, double ValSwap, double FwdDaycount, double FwdFreq,
									  ARM_Date &DtAsofFwd,ARM_Date &DtRollFwd, double DiscStartFwd,
									  ARM_Date &DtThFirstSw, ARM_Date &DtThEndSw, double DiscStartSw,
									  char* calName, int Cont_Lin, int FwdAdj)
{
	double PteFwd, ecart, px, pplush, pmoinsh, deriv;
	int nbiter;

	PteFwd = 0.;

	
	pplush = CptSwapFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  StartFwd, PteFwd + 0.0001, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);

	pmoinsh = CptSwapFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  StartFwd, PteFwd - 0.0001, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);

	px = CptSwapFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  StartFwd, PteFwd, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);


	deriv = (pplush - pmoinsh) / (2. * 0.0001);
 
	nbiter = 1;

	ecart = px - ValSwap;

	while (fabs(ecart) > 1.0e-12 && nbiter < 50)
	{
		nbiter ++;
		PteFwd = PteFwd - ecart / deriv;
		pplush = CptSwapFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  StartFwd, PteFwd + 0.0001, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);

		px = CptSwapFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  StartFwd, PteFwd, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);

		deriv = (pplush - px) /0.0001;

		ecart = px - ValSwap;
	}

}

double ARM_ZeroCurve::CptSwapFwdValue(ARM_Matrix *LocDiscData, int& LocDiscSize,
									  double FixRate, double FixDaycount, double FixedFreq,
									  double StartFwd, double PteFwd, double FwdDaycount, double FwdFreq,
									  ARM_Date &DtAsofFwd,ARM_Date &DtRollFwd, double DiscStartFwd,
									  ARM_Date &DtThFirstSw, ARM_Date &DtThEndSw, double DiscStartSw,
									  char* calName, int Cont_Lin, int FwdAdj)
{
	double res, CurFwd, inter, ThetaFwd;
	double DiscFact, PrevDisc, NextDisc, DiscEndSw, FltVal, FixVal, CurDisc, FixFlow;
	int i, j, nbmonths;
	ARM_Date DtThNext, DtAdjStart, DtAdjEnd, DtAdjStartFwd, DtAdjEndFwd, DtAdjInterp;
	ARM_Date DtAdjStartSw, DtAdjEndSw, DtPrevEnd, DtPrevMatu, DtNextMatu;

	DtAdjStart = DtThFirstSw;
	DtAdjStart.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjStartSw = DtAdjStart;
	DtAdjEnd = DtAdjStart;

	DtAdjEndSw = DtThEndSw;
	DtAdjEndSw.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjStartFwd = DtAsofFwd;
	DtAdjEndFwd = DtThEndSw;

	LocDiscData->Elt(0,1) = DtAdjStart.GetJulian();
	LocDiscData->Elt(0,0) = DiscStartSw;

	LocDiscData->Elt(1,1) = DtAsofFwd.GetJulian();
	LocDiscData->Elt(1,0) = DiscStartFwd;

	j = 1;

	DtAdjStartFwd = DtAsofFwd;

	if (FwdAdj == K_ADJUSTED)
		DtAdjStartFwd.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjEnd = DtAsofFwd;

	i = 0;

	while(DtAdjEnd < DtAdjEndSw)
	{
		i++;
		DtAdjEnd = DtRollFwd;
		DtAdjEnd.AddMonths(FwdFreq*(i-1), GOTO_END_OF_MONTH);
		if (FwdAdj == K_ADJUSTED)
			DtAdjEnd.GoodBusinessDay(K_MOD_FOLLOWING, calName);
	}

	DtAdjEndFwd = DtRollFwd;
	DtAdjEndFwd.AddMonths(FwdFreq*(i-1), GOTO_END_OF_MONTH);

	if (FwdAdj == K_ADJUSTED)
		DtAdjEndFwd.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DiscFact = DiscStartFwd;
	i = 1;
	DtAdjStart = DtAsofFwd;
	DtAdjEnd = DtAsofFwd;
	DtThNext = DtRollFwd;
	DtThNext.AddMonths(FwdFreq*i, GOTO_END_OF_MONTH);
	DtAdjEnd = DtThNext;
	DtPrevEnd = DtRollFwd;

	if (FwdAdj == K_ADJUSTED)
		DtAdjEnd.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjInterp = DtAdjEnd;

	while (DtAdjEnd < DtAdjEndSw)
	{
		DtAdjInterp = DtAdjEnd;

		inter = CountYears(K30_360, DtAdjStartFwd, DtAdjInterp);

		CurFwd = StartFwd + inter * PteFwd;
		ThetaFwd = CountYears(FwdDaycount, DtAdjStart, DtAdjEnd);
		DiscFact = DiscFact / (1. + CurFwd * ThetaFwd / 100.);
		j++;
		LocDiscData->Elt(j,1) = DtAdjEnd.GetJulian();
		LocDiscData->Elt(j,0) = DiscFact;
		
		DtAdjStart = DtAdjEnd;
		i++;
		DtThNext = DtRollFwd;
		DtThNext.AddMonths(FwdFreq*i, GOTO_END_OF_MONTH);
		DtPrevEnd = DtAdjEnd;
		DtAdjEnd = DtThNext;
		
		if (FwdAdj == K_ADJUSTED)
			DtAdjEnd.GoodBusinessDay(K_MOD_FOLLOWING, calName);
	}

	PrevDisc = DiscFact;
	DtPrevMatu = DtPrevEnd;

	DtAdjInterp = DtAdjEnd;
	inter = CountYears(K30_360, DtAdjStartFwd,DtAdjInterp);

	CurFwd = StartFwd + inter * PteFwd;
	ThetaFwd = CountYears(FwdDaycount, DtAdjStart,DtAdjEnd);
	NextDisc = DiscFact / (1. + CurFwd * ThetaFwd / 100.);
	DtNextMatu = DtAdjInterp;

	DiscFact = Interpol(DtPrevMatu.GetJulian(),
						  DtNextMatu.GetJulian(),
						  DtAdjEndSw.GetJulian(),
						  PrevDisc,
						  NextDisc,
						  Cont_Lin,
						  KACTUAL_360);

	j++;
	
	DiscEndSw = DiscFact;

	LocDiscData->Elt(j,1) = DtAdjEndSw.GetJulian();
	LocDiscData->Elt(j,0) = DiscEndSw;

	FltVal = 100.*(DiscStartSw-DiscEndSw);

	LocDiscSize = j;

	i = 1;
	FixVal = 0.;
	DtAdjEnd = DtAdjEndSw;

	nbmonths = 12/FixedFreq;

	DtThNext = DtThEndSw;
	DtThNext.AddMonths(nbmonths*-i, GOTO_END_OF_MONTH);
	DtAdjStart = DtThNext;
	DtAdjStart.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	while ( DtAdjStart > DtAdjStartSw )
	{
		CurDisc = LocDiscInterp(LocDiscData, LocDiscSize, itsAsOfDate, DtAdjEnd, Cont_Lin); // à corriger

		inter = CountYears(FixDaycount, DtAdjStart, DtAdjEnd);
		FixFlow = FixRate * inter;

		FixVal += FixFlow * CurDisc;
		DtAdjEnd = DtAdjStart;
		i++;
		DtThNext = DtThEndSw;
		DtThNext.AddMonths(nbmonths*-i, GOTO_END_OF_MONTH);
		DtAdjStart = DtThNext;
		DtAdjStart.GoodBusinessDay(K_MOD_FOLLOWING, calName);
	}

	DtAdjStart = DtAdjStartSw;

	CurDisc = LocDiscInterp(LocDiscData, LocDiscSize, itsAsOfDate, DtAdjEnd, Cont_Lin); // à corriger

	inter = CountYears(FixDaycount, DtAdjStart, DtAdjEnd);
	FixFlow = FixRate * inter;

	FixVal += FixFlow * CurDisc;

	res = FixVal - FltVal;

	

	return(res);
}



double ARM_ZeroCurve::LocDiscInterp(ARM_Matrix *TabDiscData, int TabDiscSize,
									ARM_Date &AsofDate, ARM_Date &CurDate,
									int Cont_Lin, int IntFlatLeft, double startDisc)
{
	double res;
	int i;
	double DtPrevMatu, DtNextMatu;
	double PrevDisc, NextDisc;

	i = 0;
	PrevDisc = startDisc; // 1. by default
	NextDisc = startDisc;
	DtPrevMatu = AsofDate.GetJulian();
	DtNextMatu = AsofDate.GetJulian();

	while (i<=TabDiscSize && DtNextMatu <= CurDate.GetJulian())
	{
		DtPrevMatu = DtNextMatu;
		PrevDisc = NextDisc;
		if (TabDiscData->Elt(i,1)> 0.)
		{
			DtNextMatu = TabDiscData->Elt(i,1);
			NextDisc = TabDiscData->Elt(i,0);
		}
		i++;
	}

	if ((DtNextMatu <= CurDate.GetJulian()) 
		|| ((IntFlatLeft) && (CurDate.GetJulian() < TabDiscData->Elt(0,1))) )
	{
		if (DtNextMatu == AsofDate.GetJulian())
			res = startDisc;
		else // Extrapole Flat on both side
		{
			res = (CurDate.GetJulian() - AsofDate.GetJulian()) / (DtNextMatu - AsofDate.GetJulian());
			res = pow(NextDisc,res);
		}
	}
	else if (CurDate <= AsofDate)
		res = startDisc;
	else
		res = Interpol(DtPrevMatu,
						DtNextMatu,
						CurDate.GetJulian(),
						PrevDisc,
						NextDisc,
						Cont_Lin,
						KACTUAL_360);
	return(res);
}



double ARM_ZeroCurve::LocSensi(ARM_Matrix *TabDiscData, int TabDiscSize,
							ARM_Date &AsofDate, ARM_Date &DtThStart, ARM_Date &DtThEnd,
							char* FreqUnit, int FreqQty, int FixDayCount,
							int Cont_Lin, char* calName)
{
	int i;
	double res, inter, CurDisc;

	ARM_Date DtCur, DtAdjPrev, DtAdjNext, DtAdjStart, DtAdjEnd;

	DtCur = DtThEnd;
	DtAdjEnd = DtCur;
	DtAdjEnd.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjNext = DtAdjEnd;

	i = 1;
	DtCur.AddMonths(FreqQty*-i, GOTO_END_OF_MONTH);
	DtAdjPrev = DtCur;
	DtAdjPrev.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjStart = DtThStart;
	DtAdjStart.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	res = 0.;

	while (DtAdjPrev > DtAdjStart)
	{
		CurDisc = LocDiscInterp(TabDiscData, TabDiscSize, AsofDate, DtAdjNext, Cont_Lin, 1);

		inter = CountYears(FixDayCount, DtAdjPrev, DtAdjNext);

		res += inter * CurDisc;

		DtAdjNext = DtAdjPrev;

		i++;
		DtCur = DtThEnd;
		DtCur.AddMonths(FreqQty*-i, GOTO_END_OF_MONTH);
		DtAdjPrev = DtCur;
		DtAdjPrev.GoodBusinessDay(K_MOD_FOLLOWING, calName);
	}

	DtAdjPrev = DtAdjStart;

	CurDisc = LocDiscInterp(TabDiscData, TabDiscSize, AsofDate, DtAdjNext, Cont_Lin, 1);

	inter = CountYears(FixDayCount, DtAdjPrev, DtAdjNext);

	res += inter * CurDisc;

	return(res);
}



double LocCptFwd(ARM_Date &DtStart, ARM_Date &DtEnd, double DiscStart, double DiscEnd, int DayCount)
{
	double res, inter;

	inter = CountYears(DayCount, DtStart, DtEnd);

	if (inter < MINDOUBLE)
		res = 0.;
	else
	{
		res = 100. * (DiscStart / DiscEnd - 1.);
		res = res / inter;
	}

	return(res);
}



double ARM_ZeroCurve::LocFwdInterp(ARM_Matrix *TabDiscData, int TabDiscSize,
									ARM_Date &AsofDate, ARM_Date &DtStartFwd, ARM_Date &DtEndFwd,
									int Cont_Lin, int FwdDayCount)
{
	double res, DiscStart, DiscEnd;

	DiscStart = LocDiscInterp(TabDiscData, TabDiscSize, AsofDate, DtStartFwd, Cont_Lin);

	DiscEnd = LocDiscInterp(TabDiscData, TabDiscSize, AsofDate, DtEndFwd, Cont_Lin);

	res = LocCptFwd(DtStartFwd, DtEndFwd, DiscStart, DiscEnd, FwdDayCount);

	return(res);
}



void FillTabFirstDisc(ARM_Matrix& MMRates, int SizeM,
					  ARM_Matrix& Future, int SizeF,
					  ARM_Matrix* Swap, int SizeS,
					  int MMVsFUT, int SwapVsFut,
					  ARM_Matrix* TabFirstDisc, int& SizeT)
{

	int i, j;
	SizeT = 0;

	if ( MMVsFUT != K_FUT ) // MM overlap FUT
	{
		i = 0;
		while (i < SizeM)
		{
			TabFirstDisc->Elt(i,0) = MMRates.Elt(i, 0);
			TabFirstDisc->Elt(i,1) = MMRates.Elt(i, 1);
			i++;
		}

		j = 0;
		while (j < SizeF)
		{
			if (Future.Elt(j, 1) > MMRates.Elt(SizeM, 1))
			{
				TabFirstDisc->Elt(i,0) = Future.Elt(j, 0);
				TabFirstDisc->Elt(i,1) = Future.Elt(j, 1);
				i++;
			}
			j++;
		}

		SizeT = i;
	}
	else // FUT overlap MM
	{
		i = 0;
		while (i < SizeM)
		{
			if (SizeF == 0 || MMRates.Elt(i, 1) < Future.Elt(1, 2))
			{
				TabFirstDisc->Elt(i,0) = MMRates.Elt(i, 0);
				TabFirstDisc->Elt(i,1) = MMRates.Elt(i, 1);
			}
			i++;
		}

		j = 0;
		while (j < SizeF)
		{
			TabFirstDisc->Elt(i,0) = Future.Elt(j, 0);
			TabFirstDisc->Elt(i,1) = Future.Elt(j, 1);
			j++;
			i++;
		}

		SizeT = i;
	}

	if (SwapVsFut == K_SWAP)
	{
		while(TabFirstDisc->Elt(i-1,1) > Swap->Elt(1,1) && i >0)
		{
			TabFirstDisc->Elt(i-1,0) = 0.;
			TabFirstDisc->Elt(i-1,1) = 0.;
			i--;
		}
		SizeT = i;
	}

}



void ARM_ZeroCurve::CptSwapZeroRatesFWD(ARM_Matrix& MMRates, int SizeM,
                                     ARM_Matrix* Swap, 
                                     int& SizeS, int fxPayFreq,
                                     ARM_Matrix& Future, int SizeF, 
                                     int MMVsFut,
                                     int SwapVsFut,
                                     int raw,
                                     int Cont_Lin,
									 int& firstIndexSwap, // for the calcul by "raw"
									 int indexAUD, // indexAUD is the last index of the AUD with freq=4
									 int FwdMatu, // Maturity of forward in months, 3 by default
									 int FwdAdj, // Adjust the Forward Interpolation date
									 int FixDayCount)
{
	int i, j, SizeT, SizeL, SizeI, FirstS;
	double sensi_prev, ValSwapPrev;
	double DiscSpot, DiscStartFutSwap, DiscPrevSwap, nextRate, prevRate, StartFwd;

	ARM_Date DtThNext, DtAdjStart, DtAdjEnd, DtThPrevMatu, DtThNextMatu, DtAdjPrevMatu, DtAdjNextMatu;
	ARM_Date DtAdjStartSw, DtAdjEndSw, DtPrevEnd, DtLastFirstDisc, DtRollFwd, DtStartFwd, DtEndFwd;

	char calName[10];
	GetCurrencyUnit()->CalcFloatPayCal(calName);

	int spotDays   = GetCurrencyUnit()->GetSpotDays();
	if (FixDayCount == KNOBASE)
       FixDayCount = GetCurrencyUnit()->GetFixedDayCount();
	int FwdDayCount = GetCurrencyUnit()->GetLiborIndexDayCount();

    int fxPayFreqVal, nbmonths;

	if(strcmp(calName, "AUD")==0 && indexAUD>0)
		fxPayFreqVal = 4;
	else
	    fxPayFreqVal = fxPayFreq ;


	// On merge les discount MM et FUT selon la règle MMVsFUT

	ARM_Matrix TabFirstDisc(700, 2, 0.0);


	FillTabFirstDisc(MMRates, SizeM, Future, SizeF, Swap, SizeS, MMVsFut, SwapVsFut, &TabFirstDisc, SizeI);

	SizeT = SizeI;
	
	// On détermine la date et la valeur du dernier discount précédent

	DtLastFirstDisc = (ARM_Date) TabFirstDisc.Elt(SizeT-1,1);
	DiscStartFutSwap = TabFirstDisc.Elt(SizeT-1,0);
	DtRollFwd = DtLastFirstDisc;


	// On cherche la maturité du 1er swap coté

	i = 1;
	
	nextRate = Swap->Elt(i,0);
	DtThNextMatu = (ARM_Date) Swap->Elt(1,1);
	DtAdjNextMatu = DtThNextMatu;
	DtAdjNextMatu.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	if (SwapVsFut == K_FUT && DtAdjNextMatu < DtLastFirstDisc) // FUT overlap SWAP
	{
		while (DtAdjNextMatu < DtLastFirstDisc)
		{
			i++;
			nextRate = Swap->Elt(i,0);
			DtThNextMatu = (ARM_Date) Swap->Elt(i,1);
			DtAdjNextMatu = DtThNextMatu;
			DtAdjNextMatu.GoodBusinessDay(K_MOD_FOLLOWING, calName);
		}
	}

	FirstS = i;
	


	// On cherche la maturité du swap précédant le dernier discount généré

	DtThPrevMatu = DtThNextMatu;

	nbmonths = 12/fxPayFreqVal;

	while (DtThPrevMatu > DtLastFirstDisc)
		DtThPrevMatu.AddMonths(-nbmonths, GOTO_END_OF_MONTH);

	DtAdjPrevMatu = DtThPrevMatu;
	DtAdjPrevMatu.GoodBusinessDay(K_MOD_FOLLOWING, calName);


	// calcul de la sensi et du discount à la matu du précédent swap (non coté)

	ARM_Date spotDate = itsAsOfDate;

	if (spotDays > 0)
		spotDate.NextBusinessDay(spotDays, calName);

	DiscSpot = LocDiscInterp(&TabFirstDisc, SizeT, itsAsOfDate, spotDate, Cont_Lin, 1); // extrapole Flat if needed

	DiscPrevSwap = LocDiscInterp(&TabFirstDisc, SizeT, itsAsOfDate, DtAdjPrevMatu, Cont_Lin);

	sensi_prev = LocSensi(&TabFirstDisc, SizeT,
							itsAsOfDate, spotDate, DtThPrevMatu,
							"M", nbmonths, FixDayCount,
							Cont_Lin, calName);


	// calcul du taux implicite du swap à la maturité précédente
	// calcul de la valeur des flux du swap courant antérieurs à la maturité du swap précédent

	ARM_Matrix LocDiscData(100, 2, 0.); // extend for future domain hedging

	if (sensi_prev > 0)
	{
		prevRate = 100.* (DiscSpot - DiscPrevSwap) / sensi_prev;
		ValSwapPrev = (nextRate - prevRate) * sensi_prev;
	

	
	// détermine les discounts de raccord entre la maturité du dernier généré 
	// dans la partie courte et celle du 1er swap coté

		DtStartFwd = DtLastFirstDisc;
		DtStartFwd.AddMonths(-FwdMatu, GOTO_END_OF_MONTH);
		if (DtStartFwd < spotDate)
			DtStartFwd = spotDate;

		DtEndFwd = DtLastFirstDisc;

		StartFwd = LocFwdInterp(&TabFirstDisc, SizeT,
							itsAsOfDate, DtStartFwd, DtEndFwd,Cont_Lin, FwdDayCount);
	


		
		FillFutSwapDisc(&LocDiscData, SizeL,
					nextRate, FixDayCount, fxPayFreqVal,
					StartFwd, -ValSwapPrev, FwdDayCount, FwdMatu,
					DtLastFirstDisc,DtRollFwd, DiscStartFutSwap,
					DtThPrevMatu, DtThNextMatu, DiscPrevSwap,
					calName, Cont_Lin, FwdAdj);
	}

	else
	{
		ValSwapPrev = 0;
		StartFwd = nextRate;

		FillFutSwapDisc(&LocDiscData, SizeL,
					nextRate, FixDayCount, fxPayFreqVal,
					StartFwd, -ValSwapPrev, FwdDayCount, 12,
					DtLastFirstDisc,DtRollFwd, DiscStartFutSwap,
					DtThPrevMatu, DtThNextMatu, DiscPrevSwap,
					calName, Cont_Lin, FwdAdj);
	}

	// transvase les nouveaux discounts

	for (i=2;i<=SizeL;i++)
	{
		TabFirstDisc.Elt(SizeT,0) = LocDiscData.Elt(i,0);
		TabFirstDisc.Elt(SizeT,1) = LocDiscData.Elt(i,1);
		SizeT++;
	}


	DtLastFirstDisc = (ARM_Date) LocDiscData.Elt(i-1,1);
	DiscStartFutSwap = LocDiscData.Elt(i-1,0);

			
	// if (FwdAdj == K_ADJUSTED)
		DtRollFwd = DtAdjNextMatu;
	// else
		// DtRollFwd = DtThNextMatu;

	// on boucle ensuite sur les prochains swaps cotés

	for (i=FirstS+1; i<=SizeS;i++)
	{
		nextRate = Swap->Elt(i,0);

		DtThNextMatu = (ARM_Date) Swap->Elt(i,1);
		DtAdjNextMatu = DtThNextMatu;
		DtAdjNextMatu.GoodBusinessDay(K_MOD_FOLLOWING, calName);

		DtThPrevMatu = DtThNextMatu;

		
		if(strcmp(calName, "AUD")==0 && i<=indexAUD)
			fxPayFreqVal = 4;
		else
			fxPayFreqVal = fxPayFreq ;

		nbmonths = 12/fxPayFreqVal;


		while (DtThPrevMatu > DtLastFirstDisc)
			DtThPrevMatu.AddMonths(-nbmonths, GOTO_END_OF_MONTH);

		DtAdjPrevMatu = DtThPrevMatu;
		DtAdjPrevMatu.GoodBusinessDay(K_MOD_FOLLOWING, calName);

		sensi_prev = LocSensi(&TabFirstDisc, SizeT,
							itsAsOfDate, spotDate, DtThPrevMatu,
							"M", nbmonths, FixDayCount,
							Cont_Lin, calName);


		if (sensi_prev > 0)
		{
			DiscPrevSwap = LocDiscInterp(&TabFirstDisc, SizeT, itsAsOfDate, DtAdjPrevMatu, Cont_Lin); 
			prevRate = 100.* (DiscSpot - DiscPrevSwap) / sensi_prev;
			ValSwapPrev = (nextRate - prevRate) * sensi_prev;
		}
		else
			ValSwapPrev = 0;

		// détermine les discounts de raccord entre la maturité du dernier généré 
		// dans la partie courte et celle du 1er swap coté

		DtStartFwd = DtAdjPrevMatu;
		DtStartFwd.AddMonths(-FwdMatu, GOTO_END_OF_MONTH);
		DtEndFwd = DtAdjPrevMatu;

		DtStartFwd = DtLastFirstDisc;
		DtStartFwd.AddMonths(-FwdMatu, GOTO_END_OF_MONTH);
		DtEndFwd = DtLastFirstDisc;


		StartFwd = LocFwdInterp(&TabFirstDisc, SizeT,
							itsAsOfDate, DtStartFwd, DtEndFwd,Cont_Lin, FwdDayCount);

		ARM_Matrix LocDiscData(200, 2, 0.);


		FillFutSwapDisc(&LocDiscData, SizeL,
					nextRate, FixDayCount, fxPayFreqVal,
					StartFwd, -ValSwapPrev, FwdDayCount, FwdMatu,
					DtLastFirstDisc,DtRollFwd, DiscStartFutSwap,
					DtThPrevMatu, DtThNextMatu, DiscPrevSwap,
					calName, Cont_Lin, FwdAdj);


		// transvase les nouveaux discounts

		for (j=2;j<=SizeL;j++)
		{
			
			TabFirstDisc.Elt(SizeT,0) = LocDiscData.Elt(j,0);
			TabFirstDisc.Elt(SizeT,1) = LocDiscData.Elt(j,1);
			SizeT++;
		}

		DtLastFirstDisc = (ARM_Date) LocDiscData.Elt(j-1,1);
		DiscStartFutSwap = LocDiscData.Elt(j-1,0);

		
		// if (FwdAdj == K_ADJUSTED)
			DtRollFwd = DtAdjNextMatu;
		// else
			// DtRollFwd = DtThNextMatu;
	}


	// on met enfin à jour le tableau des taux de swap par les discounts générés
	// par le fitting des swaps (en ne reprenant pas ceux du MM et FUT

	SizeS = 0;
	for (i=SizeI;i<=SizeT;i++)
	{
		if (TabFirstDisc.Elt(i,1) > 0.)
		{
			SizeS++;
			Swap->Elt(SizeS,0) = TabFirstDisc.Elt(i,0);
			Swap->Elt(SizeS,1) = TabFirstDisc.Elt(i,1);
		}	
	}


	firstIndexSwap = 1;

}


/*****************************************************/
/*                                                   */
/*			Bump Forwards Methods                    */
/*			   ( AR - 10/10/04)                      */
/*                                                   */
/*****************************************************/



ARM_ZeroCurve* ARM_ZeroCurve::GenerateShiftCurveFwd(
      char Term[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS] , ARM_Vector* epsilon,
	  int FreqFwd, int AdjFwd)
{
	if (GetMktData() == NULL)
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "The input curve can't be shifted");
	}

	int i, j, SizeT, SizeL;
	double sensi_prev, ValSwapPrev;
	double DiscSpot=1., DiscPrevSwap, nextRate, prevRate, thetaFwd, DiscStart;
	double cur_spread, bump_value, curRate, cur_disc, cur_fwd, prev_disc, prev_fwd, TxFut, DiscPrevFwd, cur_matu;
	char cur_seg, prev_seg, bump_seg;

	ARM_Date DtThNext, DtAdjStart, DtAdjEnd, DtThPrevMatu, DtThNextMatu, DtAdjPrevMatu, DtAdjNextMatu;
	ARM_Date DtAdjStartSw, DtAdjEndSw, DtPrevEnd, DtLastFirstDisc, DtRollFwd, DtStartFwd, DtEndFwd;
	ARM_Date asofDate, spotDate, prev_datefwd, DTCur, DtNext, DtPrev, startFut, endFut, cur_datefwd, next_datefwd;
	ARM_Date prev_dateplot, next_dateplot, limite_AUD_DateQ;

	char calName[10]; char* cur_plot=""; char* prev_plot="";
	GetCurrencyUnit()->CalcFloatPayCal(calName);

	int spotDays   = GetCurrencyUnit()->GetSpotDays();
    // int FixDayCount = GetCurrencyUnit()->GetFixedDayCount();
    int FixDayCount = GetFixDayCount();
	int FwdDayCount = GetCurrencyUnit()->GetLiborIndexDayCount();
	int fxPayFreq = GetCurrencyUnit()->GetFixedPayFreq();

	if (FixDayCount == KNOBASE)
		FixDayCount = GetCurrencyUnit()->GetFixedDayCount();

	int Cont_Lin = itsMktData->itsCont_Lin; // itsInterpolMeth;

	int MMVsFut = itsMktData->itsMMVsFut; // Priority MM/FUT

	int SwapVsFut = itsMktData->itsSwapVsFut; // Priority Fut/SWAP

	int priority = K_YES;


	int fxPayFreqVal, nbmonths;

	limite_AUD_DateQ = spotDate;

	limite_AUD_DateQ.AddYears(3); // Dans la courbe AUD les swaps <= 3Y sont trimestriels

	
	prev_seg='M';

	ARM_ZeroLInterpol* zc=NULL;
	ARM_Date tmpdat;
	ARM_MarketData* MktData = (ARM_MarketData*) GetMktData()->Clone();

	asofDate = GetAsOfDate();
	spotDate = asofDate;

	if (spotDays > 0)
	{
		spotDate.NextBusinessDay(spotDays, calName);
		DiscSpot = DiscountPrice(spotDate);
	}

	
	
	// verification de l'existence des datamarkets
	if (!MktData)
		return NULL;

	int szinp =0;
	szinp = epsilon->GetSize();

	int k=0;
	int l=0;

	DtThPrevMatu = GetAsOfDate();
	DtAdjPrevMatu = DtThPrevMatu;

	DtAdjStart = asofDate;
	DtThNextMatu = spotDate;
	
	ARM_Matrix TabFirstDisc(700, 2, 0.);

	cur_disc = 1.;
	SizeT = 0;
	prev_seg = 'M';
	prev_datefwd = asofDate;


	for (k=0; k<MktData->itsMktValue->GetSize(); k++)
	{
		// strcpy(cur_plot, Term[k]);
		cur_plot = MktData->itsMktTerms[k];

		curRate = MktData->itsMktValue->Elt(k);

		// determine nature du plot en MM, FUT ou SERIAL, SWAP

		if ((prev_seg == 'M') && (cur_seg == 'M') && (prev_datefwd <= spotDate))
			DtAdjStart = prev_datefwd;
		else
			DtAdjStart = spotDate;

		DecodePlot (cur_plot, asofDate, spotDate, calName, cur_seg, DtThNextMatu, DtAdjStart);

		if (cur_seg == 'F')
		{
			startFut = DtAdjStart;
			endFut = DtThNextMatu;
			startFut.GoodBusinessDay(K_MOD_FOLLOWING, calName);
			endFut.GoodBusinessDay(K_MOD_FOLLOWING, calName);
		}

		next_dateplot = DtThNextMatu;
		DtAdjNextMatu = DtThNextMatu;
		DtAdjNextMatu.GoodBusinessDay(K_MOD_FOLLOWING, calName);

		nextRate = MktData->itsMktValue->Elt(k);
		cur_spread = 0.;

		if( (strcmp(calName, "AUD")==0) && 
			(DtThNextMatu.GetJulian() <= limite_AUD_DateQ.GetJulian()) )
			fxPayFreqVal = 4;
		else
			fxPayFreqVal = fxPayFreq ;

		nbmonths = 12/fxPayFreqVal;


		// recuperation du bump relatif au plot
		bump_value = 0.;
		bump_seg = ' ';

		for (l=0; l<szinp; l++)
		{
			if (!strcmp(Term[l],MktData->itsMktTerms[k]))
			{
				bump_value = epsilon->Elt(l);
				bump_seg = cur_seg;

				// Cas Futur
				if ( MktData->itsMktValue->Elt(k) > 50.0 )
					MktData->itsMktValue->Elt(k) = 
                           MktData->itsMktValue->Elt(k) - epsilon->Elt(l);
				else
					MktData->itsMktValue->Elt(k) = 
                           MktData->itsMktValue->Elt(k) + epsilon->Elt(l);

				l=szinp;
			}
		}

		ARM_Matrix LocDiscData(300, 2, 0.);

		switch(prev_seg)
			{
				case 'Y' :  // Swap segment
					// On cherche la maturité du swap précédant le dernier discount généré

					DtThPrevMatu = DtThNextMatu;

					while (DtThPrevMatu > prev_dateplot)
						DtThPrevMatu.AddMonths(-nbmonths, GOTO_END_OF_MONTH);

					DtAdjPrevMatu = DtThPrevMatu;
					DtAdjPrevMatu.GoodBusinessDay(K_MOD_FOLLOWING, calName);

					if( (strcmp(calName, "AUD")==0) && 
						(DtThNextMatu.GetJulian() <= limite_AUD_DateQ.GetJulian()) )
						fxPayFreqVal = 4;
					else
						fxPayFreqVal = fxPayFreq ;

					nbmonths = 12/fxPayFreqVal;


					// calcul de la sensi et du discount à la matu du précédent swap (non coté)

					DiscPrevSwap = LocDiscInterp(&TabFirstDisc, SizeT, asofDate, DtAdjPrevMatu, Cont_Lin);

					sensi_prev = LocSensi(&TabFirstDisc, SizeT,
											asofDate, spotDate, DtThPrevMatu,
											"M", nbmonths, FixDayCount,
											Cont_Lin, calName);


					// calcul du taux implicite du swap à la maturité précédente
					// calcul de la valeur des flux du swap courant antérieurs à la maturité du swap précédent
	
					if (sensi_prev > 0)
					{
						prevRate = 100.* (DiscSpot - DiscPrevSwap) / sensi_prev;
						ValSwapPrev = ((nextRate+bump_value) - prevRate) * sensi_prev;
					}
					else
						ValSwapPrev = 0;

					if( (strcmp(calName, "AUD")==0) && 
						(DtThNextMatu.GetJulian() <= limite_AUD_DateQ.GetJulian()) )
						fxPayFreqVal = 4;
					else
						fxPayFreqVal = fxPayFreq ;

					nbmonths = 12/fxPayFreqVal;


					cur_spread = FillDiscFValToSpr(&LocDiscData, SizeL, -ValSwapPrev,
									  nextRate+bump_value, FixDayCount, fxPayFreqVal,
									  FwdDayCount, FreqFwd,
									  prev_datefwd, prev_datefwd, DiscPrevFwd,
									  DtThPrevMatu, DtThNextMatu, DiscPrevSwap,
									  calName, Cont_Lin, AdjFwd);

					// transvase les nouveaux discounts

					for (i=2;i<=SizeL;i++)
					{
						TabFirstDisc.Elt(SizeT,0) = LocDiscData.Elt(i,0);
						TabFirstDisc.Elt(SizeT,1) = LocDiscData.Elt(i,1);
						SizeT++;
					}

					prev_datefwd = LocDiscData.Elt(i-1,1);
					DiscPrevFwd = LocDiscData.Elt(i-1,0);

					break;

				default : // Other cases in MM, FUT or SERIAL
					switch(cur_seg)
					{
					case 'M' : // Money Market segment
						DtAdjNextMatu = DtThNextMatu;
						DtAdjNextMatu.GoodBusinessDay(K_MOD_FOLLOWING, calName);
						prev_disc = LocDiscInterp(&TabFirstDisc, SizeT, asofDate, DtAdjPrevMatu, Cont_Lin);
						cur_fwd = ForwardYield(DtAdjPrevMatu, DtAdjNextMatu, -1, 1); // Forward monétaire ajusté du daycount
						// cur_spread = spread_MM(cur_fwd, asofDate, spotDate, DtAdjPrevMatu, DtAdjNextMatu, nextRate+bump_value, prev_disc, DiscSpot, FwdDayCount);
						DiscStart = LocDiscInterp(&TabFirstDisc, SizeT, asofDate, DtAdjStart, Cont_Lin);
						cur_spread = spread_MM(cur_fwd, asofDate, DtAdjStart, DtAdjPrevMatu, DtAdjNextMatu, nextRate+bump_value, prev_disc, DiscStart, FwdDayCount);

						thetaFwd = CountYears(FwdDayCount, DtAdjPrevMatu, DtAdjNextMatu);
						cur_disc = prev_disc / (1. + (cur_fwd+(cur_spread/100.))/100.*thetaFwd);
						SizeT = InsereDisc(DtAdjNextMatu, cur_disc, &TabFirstDisc, SizeT, K_YES);

						if (DtAdjNextMatu == spotDate)
							DiscSpot = cur_disc;
						
						break;
					case 'Y' : // Swap segment
						DtThPrevMatu = DtThNextMatu;

						if( (strcmp(calName, "AUD")==0) && 
							(DtThNextMatu.GetJulian() <= limite_AUD_DateQ.GetJulian()) )
							fxPayFreqVal = 4;
						else
							fxPayFreqVal = fxPayFreq ;

						nbmonths = 12/fxPayFreqVal;

						while (DtThPrevMatu > prev_datefwd)
							DtThPrevMatu.AddMonths(-nbmonths, GOTO_END_OF_MONTH);

						DtAdjPrevMatu = DtThPrevMatu;
						DtAdjPrevMatu.GoodBusinessDay(K_MOD_FOLLOWING, calName);


						// calcul de la sensi et du discount à la matu du précédent swap (non coté)

						DiscPrevSwap = LocDiscInterp(&TabFirstDisc, SizeT, asofDate, DtAdjPrevMatu, Cont_Lin);

						sensi_prev = LocSensi(&TabFirstDisc, SizeT,
												asofDate, spotDate, DtThPrevMatu,
												"M", nbmonths, FixDayCount,
												Cont_Lin, calName);


						// calcul du taux implicite du swap à la maturité précédente
						// calcul de la valeur des flux du swap courant antérieurs à la maturité du swap précédent
	
						DiscPrevFwd = LocDiscInterp(&TabFirstDisc, SizeT, asofDate, prev_datefwd, Cont_Lin);

						if (sensi_prev > 0.)
						{
							prevRate = 100.* (DiscSpot - DiscPrevSwap) / sensi_prev;
							ValSwapPrev = ((nextRate+bump_value) - prevRate) * sensi_prev;
						

							// calcul de la sensi et du discount à la matu du précédent MM ou FUT

							
							cur_spread = FillDiscFValToSpr(&LocDiscData, SizeL, -ValSwapPrev,
									  nextRate+bump_value, FixDayCount, fxPayFreqVal,
									  FwdDayCount, FreqFwd,
									  prev_datefwd, prev_datefwd, DiscPrevFwd,
									  DtThPrevMatu, DtThNextMatu, DiscPrevSwap,
									  calName, Cont_Lin, AdjFwd);
						}
						else
						{
							ValSwapPrev = 0.;

							cur_spread = FillDiscFValToSpr(&LocDiscData, SizeL, -ValSwapPrev,
									  nextRate+bump_value, FixDayCount, fxPayFreqVal,
									  FwdDayCount, 12,
									  prev_datefwd, prev_datefwd, DiscPrevFwd,
									  DtThPrevMatu, DtThNextMatu, DiscPrevSwap,
									  calName, Cont_Lin, AdjFwd);
						}

						// transvase les nouveaux discounts

						for (i=2;i<=SizeL;i++)
						{
							TabFirstDisc.Elt(SizeT,0) = LocDiscData.Elt(i,0);
							TabFirstDisc.Elt(SizeT,1) = LocDiscData.Elt(i,1);
							SizeT++;
						}

						prev_datefwd = LocDiscData.Elt(i-1,1);
						DiscPrevFwd = LocDiscData.Elt(i-1,0);
						break;

					default : // case FUTURES or SERIAL
						// determine les dates de départ et d'échéance du taux ou du contrat

						if (curRate > 50.)
							TxFut = 100. - curRate;
						else
							TxFut = curRate;

						if ((cur_seg == 'F') || (cur_seg == 'S'))
						{
							prev_disc = LocDiscInterp(&TabFirstDisc, SizeT, asofDate, startFut, Cont_Lin);
							prev_datefwd = startFut;

							if (MMVsFut == K_FUT)
							{
								priority = K_YES;
								SizeT = InsereDisc(startFut, prev_disc, &TabFirstDisc, SizeT, priority);
							}
							else
								priority = K_NO;

							
							thetaFwd = CountYears(FwdDayCount, startFut, endFut);
								
							cur_disc = prev_disc / (1. + (TxFut+bump_value)/100.*thetaFwd);
							
							// cur_disc = DiscExtrapol_MFUT(asofDate, prev_datefwd, startFut, endFut,
							//								TxFut+bump_value, prev_disc, FwdDayCount);
							

							SizeT = InsereDisc(endFut, cur_disc, &TabFirstDisc, SizeT, priority);
						}

					} // end switch cur_seg
		} // end switch prev_seg

		/*
		DTCur = DtPrev;
		prev_datefwd = DTCur;
		prev_datefwd.GoodBusinessDay(K_MOD_FOLLOWING, calName);
		cur_datefwd = prev_datefwd;
		next_datefwd = DtNext;
		next_datefwd.GoodBusinessDay(K_MOD_FOLLOWING, calName);
		*/

		prev_plot = cur_plot;
		prev_fwd = cur_fwd;
		prev_datefwd = DtAdjNextMatu;
		prev_dateplot = next_dateplot;
		DtPrev = DtNext;
		DtAdjPrevMatu = DtAdjNextMatu;
		prev_disc = cur_disc;
		prev_seg = cur_seg;

	}

	ARM_Vector* zeroYields = new ARM_Vector(SizeT);
	ARM_Vector* yearTerms = new ARM_Vector(SizeT);

	for (j=0;j<SizeT;j++)
	{
		DTCur = (ARM_Date) TabFirstDisc.Elt(j,1);
		cur_matu = ((ARM_Date) TabFirstDisc.Elt(j,1) - asofDate)/K_YEAR_LEN;
		cur_disc = TabFirstDisc.Elt(j,0);

		yearTerms->Elt(j) = cur_matu;
		zeroYields->Elt(j) = -100. * log(cur_disc) / cur_matu;
	}

	int lastBucketInt = 0;

	zc = new ARM_ZeroLInterpol(asofDate, 
                     yearTerms,
                     zeroYields, 
                     K_COMP_CONT, lastBucketInt, 
                     Cont_Lin, 
                     GetCurrencyUnit());

    if (yearTerms)
       delete yearTerms;

    if (zeroYields)
       delete zeroYields;

	zc->SetMktData(MktData->itsMktTerms,
                        MktData->itsMktValue,
                        MktData->itsMMVsFut,
                        MktData->itsSwapVsFut,
                        MktData->itsraw,
                        Cont_Lin);

	delete MktData;

	return(zc);
}


void ARM_ZeroCurve::DecodePlot ( char* Terms, ARM_Date &asofDate, ARM_Date &spotDate, char* ccyName, char &seg, ARM_Date& matDate, ARM_Date& startDate)
{
	int Nb, month, year;
	ARM_Date matuDateFormat;
	char matu, matu2;
	
	// startDate = spotDate; // case of MM before spotDate to be checked
	matuDateFormat = spotDate; // for init purpose
	Nb=0; month=0; year=0;

	// si 1er caract alphanum -> Futures
	if (isalpha(Terms[0]))
	{
		GetMonthYearFromExpiryDate(Terms, &month, &year);
		matu='Z';
		seg = 'F';
		ARM_Date tmpDate(1, month, year);
		matDate = tmpDate;
	}
	else // sinon MM Serial ou Swap
	{
		if	(	(( Terms[2] == '.' ) || ( Terms[2] == '/' ))
			&&
			(( Terms[5] == '.' ) || ( Terms[5] == '/' ))
			) // A date in Money Market segment
		{
			if ( strlen(Terms) == 10 ) // DD.MM.YYYY
			{
				char c1, c2;
				
				int d;
				int m;
				int y;
				
				sscanf(Terms, "%d%c%d%c%4d", &d, &c1, &m, &c2, &y);
				
				ARM_Date tmpDate(d, m, y);
				
				matuDateFormat = tmpDate;
				
				matu = 'S';
				seg = 'S';
			}
			else if ( strlen(Terms) == 8 ) // DD.MM.YY
			{
				char c1, c2;
				
				int d;
				int m;
				int y;
				
				sscanf(Terms, "%d%c%d%c%2d", &d, &c1, &m, &c2, &y);
				
				y = 2000+y;
				
				ARM_Date tmpDate(d, m, y);
				
				matuDateFormat = tmpDate;
				
				matu = 'S';
				seg = 'S';
			}
			if (strlen (Terms) == 17 ) // FRA's or Serial case: MM/DD/YY-MM/DD/YY
			{
                char c='-';

                int d1,d2;
                int m1,m2;
                int y1,y2;
				d1=0; d2=0; m1=0; m2=0; y1=0; y2=0;

                if ( strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0 )
				{
                   sscanf(Terms, "%d%c%d%c%2d%c%d%c%d%c%2d", &m1, &c, &d1,
                          &c, &y1, &c, &m2, &c, &d2, &c, &y2 );
				}
                else
				{
                   sscanf(Terms, "%d%c%d%c%2d%c%d%c%d%c%2d", &d1, &c, &m1,
                          &c, &y1, &c, &d2, &c, &m2, &c, &y2 );
				}

				

                y1 = 2000+y1;
                ARM_Date tmpDate(d1, m1, y1);
                // startDate = tmpDate.GetJulian();// chge format from MM/JJ/YY to julian
				startDate = tmpDate;

                y2 = 2000+y2;
                ARM_Date tmpDate2(d2, m2, y2);
                // matDate = tmpDate2.GetJulian();// chge format from MM/JJ/YY to julian
				matDate = tmpDate2;

                matu = 'R';
				seg = 'F';
			}
			else
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					"Invalid maturity in curve inputs");
			}
		}
		else
		{
			//sscanf(Terms[i], "%d%c", &Nb, &matu); 
			matu2 = 'X';
			sscanf(Terms, "%d%c%c", &Nb, &matu,&matu2);
			
			matu = toupper(matu);
		}
	}
	
	// Conversion des Termes char -> Termes double
	
	if ( matu == 'D' ) // Ex : "1D"
	{    
		// startDate = asofDate;
		matDate = asofDate;
		
		// if ( (spotDate==asofDate) || Nb==1)
		// {
			matDate.NextBusinessDay(Nb, ccyName);
		// }
		seg = 'M';
	}
	else if ( matu == 'W' )  
	{   //  Ex : "1W"    
		
		matDate = spotDate;
		
		matDate.AddDays(Nb*7);
		
		matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

		seg = 'M';
	}
	else if ((( matu == 'M' ) && ( Nb <= 12 )) 
		||
		( matu == 'S' )
		)
	{   //  Ex : "9M"    
		
		if ( matu == 'M' ) // cas Money market classique
		{
			matDate = spotDate;
            
			matDate.AddMonths(Nb, GOTO_END_OF_MONTH);
		}
		else // cas date
		{
			matDate = matuDateFormat; 
		}
		
		matDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);
	}
	else if (( matu == 'Y')  // ->implicitement ce sont des taux de swap
		||
		(( matu == 'M' ) && ( Nb > 12 ))
		)
	{   
		//  Ex : "15Y"
		
		matDate = spotDate;;
		
		if (( matu == 'M' ) && ( Nb > 12 ))
		{
			matDate.AddMonths(Nb, GOTO_END_OF_MONTH);
		}
		else
		{
			matDate.AddYears(Nb);
		}

		seg = 'Y'; // Swap segment
	}
	// traitement du segment Futures
	else if ( seg == 'F' && month > 0 )
	{  
		startDate = matDate;

		startDate.ChangeDate(1, month, year);
		
		startDate.PiborDelivery();

		matDate.AddMonths(3);
		
		matDate.PiborDelivery();
	}
	
}


double ARM_ZeroCurve::FillDiscFValToSpr(ARM_Matrix* LocDiscData, int& LocDiscSize, double ValSwap,
									  double FixRate, double FixDaycount, double FixedFreq,
									  double FwdDaycount, double FwdFreq,
									  ARM_Date &DtAsofFwd,ARM_Date &DtRollFwd, double DiscStartFwd,
									  ARM_Date &DtThFirstSw, ARM_Date &DtThEndSw, double DiscStartSw,
									  char* calName, int Cont_Lin, int FwdAdj)
{
	double SprFwd, ecart, px, pplush, pmoinsh, deriv, cur_fwd;
	int nbiter, j;
	ARM_Date DtThNext, DtAdjStart, DtAdjEnd, DtAdjEndSw;

	ARM_Matrix TabFwd(200, 2, 0.);

	// Fill Tab of Fwd from ref curve

	DtAdjStart = DtRollFwd;
	DtAdjStart.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjEndSw = DtThEndSw;
	DtAdjEndSw.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	j=0;
	DtThNext = DtRollFwd;
	DtAdjEnd = DtRollFwd;

	while (DtAdjEnd < DtAdjEndSw)
	{
		j=j+1;
		DtThNext = DtRollFwd;
		DtThNext.AddMonths(j*FwdFreq,GOTO_END_OF_MONTH);
		if (DtThNext > DtAdjEndSw)
			DtAdjEnd = DtAdjEndSw;
		else
			DtAdjEnd = DtThNext;

		if(FwdAdj == K_ADJUSTED)
			DtAdjEnd.GoodBusinessDay(K_MOD_FOLLOWING, calName);

		cur_fwd = ForwardYield(DtAdjStart, DtAdjEnd, -1, 1); // Forward monétaire ajusté du daycount

		TabFwd.Elt(j-1,1) = DtAdjStart.GetJulian();
		TabFwd.Elt(j-1,0) = cur_fwd;
		DtAdjStart = DtAdjEnd;
	}

	SprFwd = 0.;


	pplush = CptSwapSprFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  TabFwd, SprFwd + 0.0001, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);

	pmoinsh = CptSwapSprFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  TabFwd, SprFwd - 0.0001, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);

	px = CptSwapSprFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  TabFwd, SprFwd, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);

	deriv = (pplush - pmoinsh) / (2. * 0.0001);
 
	nbiter = 1;

	ecart = px - ValSwap;

	while (fabs(ecart) > 1.0e-12 && nbiter < 50)
	{
		nbiter ++;
		SprFwd = SprFwd - ecart / deriv;
		pplush = CptSwapSprFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  TabFwd, SprFwd + 0.0001, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);

		px = CptSwapSprFwdValue(LocDiscData, LocDiscSize,
									  FixRate, FixDaycount, FixedFreq,
									  TabFwd, SprFwd, FwdDaycount, FwdFreq,
									  DtAsofFwd,DtRollFwd, DiscStartFwd,
									  DtThFirstSw, DtThEndSw, DiscStartSw,
									  calName, Cont_Lin, FwdAdj);

		deriv = (pplush - px) /0.0001;

		ecart = px - ValSwap;
	}

	return(SprFwd);
}


double ARM_ZeroCurve::CptSwapSprFwdValue(ARM_Matrix *LocDiscData, int& LocDiscSize,
									  double FixRate, double FixDaycount, double FixedFreq,
									  ARM_Matrix TabFwd, double SprFwd, double FwdDaycount, double FwdFreq,
									  ARM_Date &DtAsofFwd,ARM_Date &DtRollFwd, double DiscStartFwd,
									  ARM_Date &DtThFirstSw, ARM_Date &DtThEndSw, double DiscStartSw,
									  char* calName, int Cont_Lin, int FwdAdj)
{
	double res, CurFwd, inter, ThetaFwd;
	double DiscFact, PrevDisc, NextDisc, DiscEndSw, FltVal, FixVal, CurDisc, FixFlow;
	int i, j, l, nbmonths;
	ARM_Date DtThNext, DtAdjStart, DtAdjEnd, DtAdjStartFwd, DtAdjEndFwd, DtAdjInterp;
	ARM_Date DtAdjStartSw, DtAdjEndSw, DtPrevEnd, DtPrevMatu, DtNextMatu, asofDate;

	asofDate = GetAsOfDate();
	DtAdjStart = DtThFirstSw;
	DtAdjStart.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjStartSw = DtAdjStart;
	DtAdjEnd = DtAdjStart;

	DtAdjEndSw = DtThEndSw;
	DtAdjEndSw.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjStartFwd = DtAsofFwd;
	DtAdjEndFwd = DtThEndSw;

	LocDiscData->Elt(0,1) = DtAdjStart.GetJulian();
	LocDiscData->Elt(0,0) = DiscStartSw;

	LocDiscData->Elt(1,1) = DtAsofFwd.GetJulian();
	LocDiscData->Elt(1,0) = DiscStartFwd;

	j = 1; // indice discounts

	DtAdjStartFwd = DtAsofFwd;

	if (FwdAdj == K_ADJUSTED)
		DtAdjStartFwd.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjEnd = DtAsofFwd;

	i = 0; // indice dates forwards

	while(DtAdjEnd < DtAdjEndSw)
	{
		i++;
		DtAdjEnd = DtRollFwd;
		DtAdjEnd.AddMonths(FwdFreq*(i-1), GOTO_END_OF_MONTH);
		if (FwdAdj == K_ADJUSTED)
			DtAdjEnd.GoodBusinessDay(K_MOD_FOLLOWING, calName);
	}

	DtAdjEndFwd = DtRollFwd;
	DtAdjEndFwd.AddMonths(FwdFreq*(i-1), GOTO_END_OF_MONTH);

	if (FwdAdj == K_ADJUSTED)
		DtAdjEndFwd.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DiscFact = DiscStartFwd;
	i = 1;
	DtAdjStart = DtAsofFwd;
	DtAdjEnd = DtAsofFwd;
	DtThNext = DtRollFwd;
	DtThNext.AddMonths(FwdFreq*i, GOTO_END_OF_MONTH);
	DtAdjEnd = DtThNext;
	DtPrevEnd = DtRollFwd;

	if (FwdAdj == K_ADJUSTED)
		DtAdjEnd.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	DtAdjInterp = DtAdjEnd;

	l = 0; // indice forwards

	while (DtAdjEnd < DtAdjEndSw)
	{
		DtAdjInterp = DtAdjEnd;

		// CurFwd = StartFwd + inter * PteFwd;
		CurFwd = TabFwd.Elt(l,0) + SprFwd / 100.;
		ThetaFwd = CountYears(FwdDaycount, DtAdjStart, DtAdjEnd);
		DiscFact = DiscFact / (1. + CurFwd * ThetaFwd / 100.);
		l++;
		j++;
		LocDiscData->Elt(j,1) = DtAdjEnd.GetJulian();
		LocDiscData->Elt(j,0) = DiscFact;
		
		DtAdjStart = DtAdjEnd;
		i++;
		DtThNext = DtRollFwd;
		DtThNext.AddMonths(FwdFreq*i, GOTO_END_OF_MONTH);
		DtPrevEnd = DtAdjEnd;
		DtAdjEnd = DtThNext;
		
		if (FwdAdj == K_ADJUSTED)
			DtAdjEnd.GoodBusinessDay(K_MOD_FOLLOWING, calName);

		if (DtAdjEnd > DtAdjEndSw)
			DtAdjEnd = DtAdjEndSw;
	}

	PrevDisc = DiscFact;
	DtPrevMatu = DtPrevEnd;

	DtAdjInterp = DtAdjEnd;

	CurFwd = TabFwd.Elt(l,0) + SprFwd / 100.;
	ThetaFwd = CountYears(FwdDaycount, DtAdjStart,DtAdjEnd);
	NextDisc = DiscFact / (1. + CurFwd * ThetaFwd / 100.);
	DtNextMatu = DtAdjInterp;

	DiscFact = Interpol(DtPrevMatu.GetJulian(),
						  DtNextMatu.GetJulian(),
						  DtAdjEndSw.GetJulian(),
						  PrevDisc,
						  NextDisc,
						  Cont_Lin,
						  KACTUAL_360);

	j++;

	DiscEndSw = DiscFact;

	LocDiscData->Elt(j,1) = DtAdjEndSw.GetJulian();
	LocDiscData->Elt(j,0) = DiscEndSw;

	
	FltVal = 100.*(DiscStartSw-DiscEndSw);

	LocDiscSize = j;

	i = 1;
	FixVal = 0.;
	DtAdjEnd = DtAdjEndSw;

	nbmonths = 12/FixedFreq;

	DtThNext = DtThEndSw;
	DtThNext.AddMonths(nbmonths*-i, GOTO_END_OF_MONTH);
	DtAdjStart = DtThNext;
	DtAdjStart.GoodBusinessDay(K_MOD_FOLLOWING, calName);

	while ( DtAdjStart > DtAdjStartSw )
	{
		CurDisc = LocDiscInterp(LocDiscData, LocDiscSize, asofDate, DtAdjEnd, Cont_Lin); // à corriger

		inter = CountYears(FixDaycount, DtAdjStart, DtAdjEnd);
		FixFlow = FixRate * inter;

		FixVal += FixFlow * CurDisc;
		DtAdjEnd = DtAdjStart;
		i++;
		DtThNext = DtThEndSw;
		DtThNext.AddMonths(nbmonths*-i, GOTO_END_OF_MONTH);
		DtAdjStart = DtThNext;
		DtAdjStart.GoodBusinessDay(K_MOD_FOLLOWING, calName);
	}

	DtAdjStart = DtAdjStartSw;

	CurDisc = LocDiscInterp(LocDiscData, LocDiscSize, asofDate, DtAdjEnd, Cont_Lin); // à corriger

	inter = CountYears(FixDaycount, DtAdjStart, DtAdjEnd);
	FixFlow = FixRate * inter;

	FixVal += FixFlow * CurDisc;

	res = FixVal - FltVal;

	return(res);
}


double ARM_ZeroCurve::spread_MM(double fwd, ARM_Date &asofDate, ARM_Date &spotDate, ARM_Date &startDate, ARM_Date &endDate,
				 double TxEnd, double DiscStart, double DiscSpot, int FwdDaycount)
{
	double res, termStart, termEnd;

	if (endDate > spotDate)
	{
		termStart = CountYears(FwdDaycount, spotDate, endDate);
		termEnd = CountYears(FwdDaycount, startDate, endDate);
		res = 100.*((1. + (TxEnd/100.)*termStart)*(DiscStart/DiscSpot)-1.)/termEnd;
	}
	else
	{
		termStart = CountYears(FwdDaycount, asofDate, endDate);
		termEnd = CountYears(FwdDaycount, startDate, endDate);
		res = 100.*((1. + (TxEnd/100.)*termStart)*(DiscSpot/DiscSpot)-1.)/termEnd;
	}

	res = (res-fwd)*100.;

	return(res);
}


double ARM_ZeroCurve::spread_MFUT(double fwd, ARM_Date &asofDate, ARM_Date &prevDate, ARM_Date &startFut, ARM_Date &endFut,
				 double TxFut, double DiscPrev, double DiscStart, int FwdDaycount)
{
	double res, termFut, termFwd, DiscEnd;

	termFut = CountYears(FwdDaycount, startFut, endFut);
	termFwd = CountYears(FwdDaycount, prevDate, endFut);

	DiscEnd = DiscStart/(1.+(TxFut/100.)*termFut);

	res = 1.+(TxFut/100.)*termFut;
	res = 100.*(res*DiscPrev/DiscStart-1.)/termFwd;
	res = (res-fwd)*100.;

	return(res);
}


double ARM_ZeroCurve::DiscExtrapol_MFUT(ARM_Date &asofDate, ARM_Date &prevDate, ARM_Date &startFut, ARM_Date &endFut,
				 double TxFut, double DiscPrev, int FwdDaycount)
{
	double termFut, termPrev, termStartFut, termEndFut, DiscEnd, prorata, zeroPrev, DiscFwd, zeroEnd, zeroStart, DiscStart;

	int mode=1; // non passé en paramètre pour l'instant
	termFut = CountYears(FwdDaycount, startFut, endFut);
	termPrev = (prevDate.GetJulian() - asofDate.GetJulian())/365.;
	termStartFut = (startFut.GetJulian() - asofDate.GetJulian())/365.;
	termEndFut = (endFut.GetJulian() - asofDate.GetJulian())/365.;

	zeroPrev = -100.*log(DiscPrev)/termPrev; // zero en taux reel : 0,10 pour 10%
	DiscFwd = (1.+(TxFut/100.)*termFut);

	if (mode == 0) // on aligne zeroPrev, zeroStartFut et zeroEndFut
	{
		prorata = (startFut.GetJulian() - prevDate.GetJulian())/365.;
		prorata = prorata / ((endFut.GetJulian() - prevDate.GetJulian())/365.);
		zeroEnd = (zeroPrev/100.)*termStartFut*(1.-prorata) + log(DiscFwd);
		zeroEnd = zeroEnd/(termEndFut - termStartFut*prorata)*100.;
		DiscEnd = exp(-(zeroEnd/100.)*termEndFut);
	}
	else // zeroStartFut = zeroPrev
	{
		zeroStart = zeroPrev;
		DiscStart = exp(-(zeroStart/100.)*termStartFut);
		DiscEnd = DiscStart/DiscFwd;
	}

	return(DiscEnd);
}



// insere un discount dans un tableau trié par maturité

int ARM_ZeroCurve::InsereDisc(ARM_Date &curDate, double curDisc, ARM_Matrix* Tabdisc, int curSize, int priority)
{
	int i, j;
	
	double prev_matu, next_matu;

	i = 0; 
	j = 0;

	prev_matu = Tabdisc->Elt(0,1);
	next_matu = prev_matu;

	if ( curSize == 0 )
	{
		Tabdisc->Elt(0,1) = curDate.GetJulian();
		Tabdisc->Elt(0,0) = curDisc;
		curSize++;
	}
	else if ( curDate.GetJulian() > Tabdisc->Elt(curSize-1, 1) )
	{
		Tabdisc->Elt(curSize,1) = curDate.GetJulian();
		Tabdisc->Elt(curSize,0) = curDisc;
		curSize++;
	}
	else
	{
		ARM_Matrix NewTab(curSize+1, 2, 0.);

		while ( i < curSize )
		{
			prev_matu = next_matu;
			next_matu = Tabdisc->Elt(i,1);

			if (( next_matu > curDate.GetJulian() ) && ( prev_matu < curDate.GetJulian() ))
			{
				NewTab.Elt(j,1) = curDate.GetJulian();
				NewTab.Elt(j,0) = curDisc;
				j++;
				
				NewTab.Elt(j,1) = Tabdisc->Elt(i,1);
				NewTab.Elt(j,0) = Tabdisc->Elt(i,0);
			}
			else if ( next_matu == curDate.GetJulian() )
			{
				if (priority == K_YES)
				{
					NewTab.Elt(j,1) = curDate.GetJulian();
					NewTab.Elt(j,0) = curDisc;
				}
				else
				{
					NewTab.Elt(j,1) = Tabdisc->Elt(i,1);
					NewTab.Elt(j,0) = Tabdisc->Elt(i,0);
				}
			}
			else
			{
				NewTab.Elt(j,1) = Tabdisc->Elt(i,1);
				NewTab.Elt(j,0) = Tabdisc->Elt(i,0);
			}

			next_matu = Tabdisc->Elt(i,1);

			i++;
			j++;
		}

		curSize = curSize+(j-i);

		for (i = 0; i < curSize; i++)
		{
			Tabdisc->Elt(i,1) = NewTab.Elt(i,1);
			Tabdisc->Elt(i,0) = NewTab.Elt(i,0);
		}

	}

	return (curSize);
}



double ComputeBasisSpread(char* matu, // "2D", "1M", ..., "1Y", ..., "10Y", ...
                          ARM_ZeroCurve* crv1, 
                          ARM_ZeroCurve* adjBScrv1,                       
                          ARM_ZeroCurve* crv2, 
                          ARM_ZeroCurve* adjBScrv2,
                          int spotUSD)
{
    double spread = 0.0;

    ARM_Date AsOfDate = crv1->GetAsOfDate();

    ARM_Y2CModel model1(crv1, adjBScrv1);
    ARM_Y2CModel model2(crv2, adjBScrv2);

    int spotDays1 = ( spotUSD < 0 ) ? crv1->GetCurrencyUnit()->GetSpotDays(): spotUSD;
    int spotDays2 = ( spotUSD < 0 ) ? crv2->GetCurrencyUnit()->GetSpotDays(): spotUSD;

    ARM_SwapLeg leg1(matu, AsOfDate, crv1->GetCurrencyUnit(), 
                     K_RCV, K_NX_BOTH, spotDays1);

    ARM_SwapLeg leg2(matu, AsOfDate, crv2->GetCurrencyUnit(), 
                     K_PAY, K_NX_BOTH, spotDays2);

    ARM_Swap swap(&leg1, &leg2);

    double swapPrice = 0.0;
    int Leg1or2      = 2;


    spread = swap.ComputeImpliedSpreadWith2Models(swapPrice, 
                                             &model1,
                                             &model2,
                                             Leg1or2);
    
    return(-1.0*spread*100.0);
}



// The result is in JPY in this case
ARM_ZeroCurve* GenerateTwoCcyBSAdjusted(ARM_ZeroCurve* crv1, // Exp: EUR: The reference CCY
                                        ARM_ZeroCurve* adjBScrv1,
                                        ARM_ZeroCurve* crv2, // Exp: JPY
                                        ARM_ZeroCurve* adjBScrv2,
                                        int spreadMatuSize,
                                        ARM_CRV_TERMS& spreadsMatu,
										int inputAsSpread,
                                        int spreadsOnly)
{
	 if (strcmp(crv1->GetCurrencyUnit()->GetCcyName(), "USD") == 0 )
	 {
		 ARM_ZeroCurve* resCrv = (ARM_ZeroCurve*) adjBScrv2->Clone();
		 return(resCrv);
	 }

    ARM_ZeroLInterpol* spreadsCrv = NULL;
    ARM_ZeroLInterpol* spreadedZc = NULL;

    ARM_CRV_TERMS matus; 

    int matuSize;

    // Remember the absolute reference in Basis context is always USD so!:
    ARM_Currency USD("USD");
    double spotUSD(USD.GetSpotDays());

    static ARM_CRV_TERMS matusSummitAUD[] = {"2D", "1M", "2M", "3M", "6M", "9M",
                                     "1Y", "2Y", "3Y", "4Y", "5Y", "6Y",
                                     "7Y", "8Y", "9Y", "10Y", "12Y", "15Y",
                                     "20Y", "25Y", "30Y", "40Y"};

    static ARM_CRV_TERMS matusSummitCCY[] = {"2D", "1M", "2M", "3M", "6M", "9M",
                                     "1Y", "2Y", "3Y", "4Y", "5Y", "6Y",
                                     "7Y", "8Y", "9Y", "10Y", "12Y", "15Y",
                                     "20Y", "25Y", "30Y"};

    int MMFreq   = K_MONTHLY;
    int swapFreq1 = K_SEMIANNUAL; // for currency of curve1
    int swapFreq2 = K_SEMIANNUAL; // for currency of curve2

    ARM_Currency ccy1 = *(crv1->GetCurrencyUnit());
    ARM_Currency ccy2 = *(crv2->GetCurrencyUnit());

    // Spreads curve Generation

    if (( strcmp(ccy1.GetCcyName(), "EUR") == 0 )
        ||
        ( strcmp(ccy1.GetCcyName(), "SEK") == 0 )
       )
    {
       swapFreq1 = K_ANNUAL;
    }

    if (( strcmp(ccy2.GetCcyName(), "EUR") == 0 )
        ||
        ( strcmp(ccy2.GetCcyName(), "SEK") == 0 )
       )
    {
       swapFreq2 = K_ANNUAL;
    }

    if ( spreadMatuSize == 0 )
    {
       if ( strcmp(ccy2.GetCcyName(), "AUD") == 0 )
       {
          memcpy(matus, matusSummitAUD, sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS); 

          matuSize = 22;
       }
       else
       {
          memcpy(matus, matusSummitCCY, sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

          matuSize = 21;
       }
    }
    else
    {
       memcpy(matus, spreadsMatu, sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

       matuSize = spreadMatuSize;
    }

    ARM_Vector spreadOfMatu(matuSize);

	if ( inputAsSpread )
	{
		// adjBScrv are spread curves, not adjusted BS curves => create them
		ARM_Date asOf = crv1->GetAsOfDate();

		adjBScrv1 = new ARM_ZeroLInterpol(asOf, adjBScrv1, crv1, MMFreq, swapFreq1, &ccy1);
		adjBScrv2 = new ARM_ZeroLInterpol(asOf, adjBScrv2, crv2, MMFreq, swapFreq2, &ccy2);
	}

    for (int i = 0; i < matuSize; i++)
    {
        try
        {
            spreadOfMatu[i] = ComputeBasisSpread(matus[i],
                                                 crv1, adjBScrv1,
                                                 crv2, adjBScrv2,
                                                 spotUSD);
        }

        catch(Exception& e)
        {
            throw e;
        }

        catch(...)
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "GenerateTwoCcyBSAdjusted,ComputeBasisSpread(): Unrecognized Failure");
        }
    }

    spreadsCrv = new ARM_ZeroLInterpol(crv2->GetAsOfDate(), 
                                       &spreadOfMatu,                                
                                       matus, 
                                       K_COMP_PROP,
                                       0, 
                                       K_LINEAR, 
                                       &ccy2);

	if ( inputAsSpread )
	{
		delete adjBScrv1;
		delete adjBScrv2;
	}

    if (spreadsOnly) 
    {       
       return(spreadsCrv);
    }
    else
    {
       spreadedZc = new ARM_ZeroLInterpol(crv2->GetAsOfDate(),
                                          (ARM_ZeroCurve *) spreadsCrv, 
                                          (ARM_ZeroCurve *) crv2,
                                          MMFreq, 
                                          swapFreq2,
                                          &ccy2);

       delete spreadsCrv;

       return(spreadedZc);
    }
}


/*--------------------------------------------------------------------------*/
/*---- End of file ----*/
